#!/usr/bin/env python3
"""
Horizon Mapper for N.I.N.A.
===========================
1. PLAN    – generate an Advanced Sequencer sequence/template that images a strip of sky
             around your (estimated) horizon, optionally keeping frames level with a rotator.
2. PROCESS – read the resulting FITS frames, place them in an Alt/Az mosaic, detect the
             horizon, and write NINA .hrz files (solid obstructions / incl. thin branches).
3. COMPARE – compare with an existing horizon file (usable sky area, plots).

Requires: Python 3.9+, numpy, Pillow   (pip install numpy pillow)
Run:      python horizon_mapper.py      (GUI)
          python horizon_mapper.py --help (command line)
"""
import argparse, hashlib, json, math, os, queue, sys, threading, time, traceback
from concurrent.futures import ThreadPoolExecutor

import numpy as np
try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:  # GUI will complain later
    Image = None

APP_NAME = "Horizon Mapper for NINA"
APP_VERSION = "1.0.1"

# =============================================================================
#  Small numeric helpers (no scipy needed)
# =============================================================================

def box2d(a, k):
    """k x k mean filter (reflect padding)."""
    p = k // 2
    b = np.pad(a, p, mode="reflect")
    c = np.cumsum(np.cumsum(b, 0), 1)
    c = np.pad(c, ((1, 0), (1, 0)))
    return (c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]) / (k * k)


def box1d(v, k):
    p = k // 2
    b = np.pad(v, p, mode="edge")
    c = np.concatenate([[0.0], np.cumsum(b)])
    return (c[k:] - c[:-k]) / k


def median_wrap(v, k):
    p = k // 2
    st = np.stack([np.roll(v, s) for s in range(-p, p + 1)])
    return np.median(st, axis=0)


# =============================================================================
#  Astronomy / geometry
# =============================================================================

def parallactic_angle(alt, az, lat):
    """Position angle (deg, E of N) of the direction to the zenith at alt/az. Time independent."""
    A, h, p = map(math.radians, (az, alt, lat))
    dec = math.asin(math.sin(p) * math.sin(h) + math.cos(p) * math.cos(h) * math.cos(A))
    ha = math.atan2(-math.sin(A) * math.cos(h),
                    math.cos(p) * math.sin(h) - math.sin(p) * math.cos(h) * math.cos(A))
    return math.degrees(math.atan2(math.sin(ha),
                                   math.tan(p) * math.cos(dec) - math.sin(dec) * math.cos(ha)))


def sun_altaz(t_utc, lat, lon):
    jd = 2440587.5 + t_utc / 86400.0
    n = jd - 2451545.0
    L = (280.46 + 0.9856474 * n) % 360
    g = math.radians((357.528 + 0.9856003 * n) % 360)
    lam = math.radians(L + 1.915 * math.sin(g) + 0.02 * math.sin(2 * g))
    eps = math.radians(23.439 - 4e-7 * n)
    ra = math.degrees(math.atan2(math.cos(eps) * math.sin(lam), math.cos(lam))) % 360
    dec = math.asin(math.sin(eps) * math.sin(lam))
    gmst = (280.46061837 + 360.98564736629 * n) % 360
    H = math.radians((gmst + lon - ra) % 360)
    p = math.radians(lat)
    alt = math.degrees(math.asin(math.sin(p) * math.sin(dec) + math.cos(p) * math.cos(dec) * math.cos(H)))
    az = math.degrees(math.atan2(-math.cos(dec) * math.sin(H),
                                 math.sin(dec) * math.cos(p) - math.cos(dec) * math.sin(p) * math.cos(H))) % 360
    return alt, az


# =============================================================================
#  Horizon files
# =============================================================================

def load_hrz(path):
    pts = []
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            parts = line.replace(",", " ").split()
            if len(parts) >= 2:
                try:
                    pts.append((float(parts[0]) % 360, float(parts[1])))
                except ValueError:
                    pass
    if not pts:
        raise ValueError("No azimuth/altitude pairs found in " + path)
    pts.sort()
    return pts


def hrz_alt(pts, az):
    az %= 360
    ext = [(pts[-1][0] - 360, pts[-1][1])] + pts + [(pts[0][0] + 360, pts[0][1])]
    for (a0, h0), (a1, h1) in zip(ext, ext[1:]):
        if a0 <= az <= a1:
            return h0 if a1 == a0 else h0 + (h1 - h0) * (az - a0) / (a1 - a0)
    return pts[0][1]


def hrz_array(pts, az):
    return np.array([hrz_alt(pts, a) for a in az])


def write_hrz(path, az, alt, step=1.0):
    """Write NINA horizon file; each output azimuth uses the max within +-step/2 (keeps narrow peaks)."""
    az = np.asarray(az); alt = np.asarray(alt)
    with open(path, "w", newline="\r\n") as f:
        a = 0.0
        while a < 360 - 1e-9:
            d = (az - a + 180) % 360 - 180
            m = np.abs(d) <= step / 2
            v = float(np.max(alt[m])) if m.any() else float(np.interp(a, az, alt, period=360))
            f.write(f"{a:06.2f}      {v:5.2f}\n")
            a += step


def sky_area_sqdeg(az_deg, alt_deg):
    """Solid angle above the horizon line up to the zenith (square degrees)."""
    az = np.asarray(az_deg); alt = np.asarray(alt_deg)
    order = np.argsort(az); az = az[order]; alt = alt[order]
    d = np.diff(np.concatenate([az, [az[0] + 360]]))
    return float(np.sum((1 - np.sin(np.radians(alt))) * np.radians(d)) * (180 / math.pi) ** 2)


# =============================================================================
#  PLAN – NINA sequence generation
# =============================================================================

SEQ = "NINA.Sequencer"
ITEMS_T = ("System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer."
           "SequenceItem.ISequenceItem, NINA.Sequencer]], System.ObjectModel")
COND_T = ("System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer."
          "Conditions.ISequenceCondition, NINA.Sequencer]], System.ObjectModel")
TRIG_T = ("System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer."
          "Trigger.ISequenceTrigger, NINA.Sequencer]], System.ObjectModel")


class PlanParams:
    def __init__(self, **kw):
        # optics
        self.sensor_w_px = 6248; self.sensor_h_px = 4176; self.pixel_um = 3.76
        self.focal_mm = 408.0
        # site
        self.lat = 61.57
        # area
        self.mode = "band"            # "band" | "horizon"
        self.alt_min = 20.0; self.alt_max = 45.0
        self.az_start = 0.0; self.az_end = 360.0
        self.horizon_files = []       # 1 or 2 files for mode "horizon"
        self.margin_below = 4.0; self.margin_above = 4.0
        self.floor = 20.0; self.ceiling = 85.0
        self.edge_frames = 1.0        # extra frame widths beside steep horizon edges (0 = off)
        self.edge_slope = 3.0         # "steep" = more than this many deg altitude per deg azimuth
        self.overlap = 0.10
        # rotator
        self.rotator = True; self.rot_offset = -58.04; self.rot_sign = 1; self.rot_tolerance = 3.0
        self.orientation = "auto"     # auto | landscape | portrait
        # exposure
        self.exposure = 0.008; self.gain = 0; self.offset = -1; self.binning = 2
        # sequence
        self.unpark = True; self.park = True; self.meridian_split = True
        self.template = False; self.collapsed = True
        # time model
        self.sec_per_frame = 8.1; self.sec_per_rot = 5.0
        for k, v in kw.items():
            setattr(self, k, v)

    @property
    def fov(self):
        w = 2 * math.degrees(math.atan(self.sensor_w_px * self.pixel_um / 1000 / 2 / self.focal_mm))
        h = 2 * math.degrees(math.atan(self.sensor_h_px * self.pixel_um / 1000 / 2 / self.focal_mm))
        return w, h


class _Ids:
    def __init__(self): self.n = 0
    def next(self): self.n += 1; return str(self.n)


def _dms(v):
    v = max(v, 0.0); d = int(v); mf = (v - d) * 60; m = int(mf); s = round((mf - m) * 60, 2)
    if s >= 60: s, m = 0, m + 1
    if m >= 60: m, d = 0, d + 1
    return d, m, s


def _container(ids, type_name, name, parent_id, expanded=True):
    cid = ids.next()
    return {
        "$id": cid,
        "$type": f"{SEQ}.Container.{type_name}, {SEQ}",
        "Strategy": {"$type": f"{SEQ}.Container.ExecutionStrategy.SequentialStrategy, {SEQ}"},
        "Name": name,
        "Conditions": {"$id": ids.next(), "$type": COND_T, "$values": []},
        "IsExpanded": expanded,
        "Items": {"$id": ids.next(), "$type": ITEMS_T, "$values": []},
        "Triggers": {"$id": ids.next(), "$type": TRIG_T, "$values": []},
        "Parent": {"$ref": parent_id} if parent_id else None,
        "ErrorBehavior": 0, "Attempts": 1,
    }, cid


def _item(ids, type_name, parent_id, **extra):
    d = {"$id": ids.next(), "$type": f"{SEQ}.SequenceItem.{type_name}, {SEQ}"}
    d.update(extra)
    d.update({"Parent": {"$ref": parent_id}, "ErrorBehavior": 0, "Attempts": 1})
    return d


def _slew(ids, pid, az, alt):
    azd, azm, azs = _dms(az % 360); ald, alm, als = _dms(alt)
    c = {"$id": ids.next(), "$type": "NINA.Astrometry.InputTopocentricCoordinates, NINA.Astrometry",
         "AzDegrees": azd, "AzMinutes": azm, "AzSeconds": azs,
         "AltDegrees": ald, "AltMinutes": alm, "AltSeconds": als}
    return _item(ids, "Telescope.SlewScopeToAltAz", pid, Coordinates=c)


def _rotate(ids, pid, angle):
    return _item(ids, "Rotator.MoveRotatorMechanical", pid, MechanicalPosition=round(angle % 360, 2))


def _expose(ids, pid, p):
    b = {"$id": ids.next(), "$type": "NINA.Core.Model.Equipment.BinningMode, NINA.Core",
         "X": p.binning, "Y": p.binning}
    return _item(ids, "Imaging.TakeExposure", pid, ExposureTime=p.exposure, Gain=p.gain,
                 Offset=p.offset, Binning=b, ImageType="LIGHT", ExposureCount=0)


def band_limits(p):
    """Return function az -> (lo, hi) of the strip to image.

    In horizon mode the strip follows the current horizon (and an optional second file) with
    margins. Where a horizon rises or drops more steeply than `p.edge_slope`, the strip is
    widened sideways by `p.edge_frames` frame widths and spans the whole height of the jump,
    so both sides of tree and roof edges are imaged.
    """
    if not (p.mode == "horizon" and p.horizon_files):
        return lambda az: (max(p.alt_min, p.floor), min(p.alt_max, p.ceiling))
    hz = [load_hrz(f) for f in p.horizon_files]
    step = 0.25
    grid = np.arange(0, 360, step)
    vals = np.array([[hrz_alt(h, a) for a in grid] for h in hz])          # (files, n)
    lo = vals.min(0) - p.margin_below
    hi = vals.max(0) + p.margin_above
    if p.edge_frames > 0:
        fw = max(p.fov)            # wider side: safe whichever way the camera is turned
        n = len(grid)
        lo_new, hi_new = lo.copy(), hi.copy()
        for v in vals:
            slope = np.abs(np.roll(v, -1) - v) / step
            for i in np.where(slope >= p.edge_slope)[0]:
                j = (i + 1) % n
                top = max(v[i], v[j]) + p.margin_above
                bot = min(v[i], v[j]) - p.margin_below
                w = p.edge_frames * fw / max(math.cos(math.radians(min(top, 80.0))), 0.17)
                k = int(math.ceil(w / step))
                idx = np.arange(i - k, j + k + 1) % n
                lo_new[idx] = np.minimum(lo_new[idx], bot)
                hi_new[idx] = np.maximum(hi_new[idx], top)
        lo, hi = lo_new, hi_new
    lo = np.maximum(lo, p.floor); hi = np.minimum(hi, p.ceiling)

    def f(az):
        x = (az % 360) / step
        i0 = int(math.floor(x)) % len(grid); i1 = (i0 + 1) % len(grid); t = x - math.floor(x)
        return lo[i0] * (1 - t) + lo[i1] * t, hi[i0] * (1 - t) + hi[i1] * t
    return f


def _effective_fov(w, h, tilt):
    t = math.radians(tilt)
    return w * math.cos(t) - h * math.sin(t), h * math.cos(t) - w * math.sin(t)


def _grid(p, fw, fh, lim):
    span_az = (p.az_end - p.az_start) % 360 or 360
    full = span_az >= 360 - 1e-9
    samples = [p.az_start + span_az * k / 720 for k in range(721)]
    lows = [lim(a)[0] for a in samples]; highs = [lim(a)[1] for a in samples]
    amin, amax = min(lows), max(highs)
    if amax <= amin:
        return []
    span = amax - amin
    step = fh * (1 - p.overlap)
    n_rows = 1 if span <= fh else math.ceil((span - fh) / step) + 1
    if n_rows == 1:
        rows = [amin + span / 2]
    else:
        lo, hi = amin + fh / 2, amax - fh / 2
        rows = [lo + i * (hi - lo) / (n_rows - 1) for i in range(n_rows)]
    # Each row is responsible for its own altitude zone (half-way to the neighbouring rows; the
    # outer rows reach their frame edge). The zones together cover [amin, amax], and every frame
    # covers its row's whole zone, so a row only needs frames where the strip enters its zone.
    zones = []
    for i, alt in enumerate(rows):
        z0 = alt - fh / 2 if i == 0 else (rows[i - 1] + alt) / 2
        z1 = alt + fh / 2 if i == len(rows) - 1 else (alt + rows[i + 1]) / 2
        zones.append((z0, z1))
    ds = 0.25                                                    # azimuth sampling step
    ns = int(round(span_az / ds))
    e_lim = [lim(p.az_start + k * ds) for k in range(ns + 1)]   # cell edges
    c_lim = [lim(p.az_start + (k + 0.5) * ds) for k in range(ns)]  # cell centres
    s_lim = [(min(e_lim[k][0], c_lim[k][0], e_lim[k + 1][0]), max(e_lim[k][1], c_lim[k][1], e_lim[k + 1][1]))
             for k in range(ns)]                                 # strip extent within each cell
    grid = []
    for alt, (z0, z1) in zip(rows, zones):
        low_edge = max(alt - fh / 2, 0.0)
        cov = fw / max(math.cos(math.radians(low_edge)), 1e-3)   # narrowest azimuth width of the frame
        stp = cov * (1 - p.overlap)
        need = [lo_ < z1 and hi_ > z0 for lo_, hi_ in s_lim]
        if not any(need):
            continue
        if all(need) and full:
            n = math.ceil(360 / stp)
            grid.append((alt, sorted((p.az_start + k * 360 / n) % 360 for k in range(n)))); continue
        # contiguous needed stretches [start, end) in degrees from az_start (joined across 0/360)
        runs = []; k = 0
        while k < ns:
            if need[k]:
                j = k
                while j < ns and need[j]:
                    j += 1
                runs.append([k * ds, j * ds]); k = j
            else:
                k += 1
        if full and len(runs) > 1 and runs[0][0] == 0 and runs[-1][1] >= span_az - 1e-9:
            last = runs.pop(); runs[0][0] = last[0] - 360
        count = lambda L: 1 if L <= cov else math.ceil((L - cov) / stp) + 1
        merged = []
        for r in runs:      # join neighbouring stretches when that needs no more frames
            if merged and count(r[1] - merged[-1][0]) <= count(merged[-1][1] - merged[-1][0]) + count(r[1] - r[0]):
                merged[-1][1] = r[1]
            else:
                merged.append(list(r))
        azs = []
        for s0, s1 in merged:
            s0 -= 0.05; s1 += 0.05                               # small safety margin at the ends
            L = s1 - s0; n = count(L)
            if n == 1:
                azs.append(s0 + L / 2)
            else:
                azs += [s0 + cov / 2 + i * (L - cov) / (n - 1) for i in range(n)]
        grid.append((alt, sorted((p.az_start + a) % 360 for a in azs)))
    if p.meridian_split:
        east = [(a, [z for z in zs if z < 180]) for a, zs in grid]
        west = [(a, [z for z in zs if z >= 180]) for a, zs in reversed(grid)]
        rows_all = [r for r in east if r[1]] + [r for r in west if r[1]]
    else:
        rows_all = grid
    out = []
    for alt, zs in rows_all:
        out.append((alt, list(reversed(zs)) if len(out) % 2 else zs))
    return out


def make_plan(p):
    """Return dict with frames [(alt, az)], fov used, orientation, rotator angles."""
    w, h = p.fov
    lim = band_limits(p)
    options = []
    orients = ["landscape", "portrait"] if (p.rotator and p.orientation == "auto") else \
        [p.orientation if p.orientation != "auto" else "landscape"]
    for o in orients:
        fw, fh = (h, w) if o == "portrait" else (w, h)
        if p.rotator:
            fw, fh = _effective_fov(fw, fh, p.rot_tolerance)
        g = _grid(p, fw, fh, lim)
        options.append((sum(len(z) for _, z in g), o, g, fw, fh))
    n, orient, grid, fw, fh = min(options, key=lambda t: t[0])
    frames, rot_moves = assign_rotator(p, grid, orient)
    secs = len(frames) * (p.sec_per_frame + p.exposure) + rot_moves * p.sec_per_rot
    return dict(grid=grid, frames=frames, orientation=orient, fw=fw, fh=fh, rot_moves=rot_moves,
                seconds=secs, lim=lim)


def assign_rotator(p, grid, orient="landscape"):
    frames = []; current = None; rot_moves = 0
    for alt, azs in grid:
        for az in azs:
            rot = None
            if p.rotator:
                sky = parallactic_angle(alt, az, p.lat) + (90 if orient == "portrait" else 0)
                mech = p.rot_sign * (sky - p.rot_offset)
                cands = [(mech + k * 180) % 360 for k in range(2)]
                tgt = min(cands) if current is None else min(cands, key=lambda c: abs(c - current))
                if current is None or abs(tgt - current) > p.rot_tolerance:
                    rot = tgt; current = tgt; rot_moves += 1
            frames.append((alt, az, rot))
    return frames, rot_moves


def fillin_boxes(det, extra_below=8.0, extra_above=6.0, floor=0.0, ceiling=85.0, pad_az=1.0,
                 lines=None, tol=1.0):
    """Turn 'clear to bottom' / 'blocked to top' azimuths into (az0, az1, alt0, alt1, label) boxes.

    `lines` = (az, solid, branches): the final (possibly edited) horizon lines. When given, an
    azimuth only counts if the lines really reach the edge of the imaged area there: the higher
    line within `tol` degrees of the top for 'blocked to top', the lower line within `tol` of the
    bottom for 'clear to bottom'. The detector's raw flags also fire on dark frame edges high in
    the sky that the solid/branches lines have already rejected.
    """
    def keep(o):
        if lines is None:
            return True
        laz, ls, lb = (np.asarray(x, float) for x in lines)
        top_ = max(np.interp(o["az"], laz, ls, period=360), np.interp(o["az"], laz, lb, period=360))
        bot_ = min(np.interp(o["az"], laz, ls, period=360), np.interp(o["az"], laz, lb, period=360))
        return top_ >= o["hi"] - tol if o["flag"] == "blocked_to_top" else bot_ <= o["lo"] + tol
    items = sorted((o["az"], o["flag"], o["lo"], o["hi"]) for o in det["details"]
                   if o["flag"] in ("clear_to_bottom", "blocked_to_top") and keep(o))
    groups = []
    for az, flag, lo, hi in items:
        if groups and groups[-1]["flag"] == flag and az - groups[-1]["az1"] <= 1.01:
            g = groups[-1]; g["az1"] = az; g["lo"] = min(g["lo"], lo); g["hi"] = max(g["hi"], hi)
        else:
            groups.append(dict(flag=flag, az0=az, az1=az, lo=lo, hi=hi))
    # merge a group that wraps over 0/360
    if len(groups) > 1 and groups[0]["flag"] == groups[-1]["flag"] and groups[0]["az0"] <= 0.01 \
            and groups[-1]["az1"] >= 359.0:
        g = groups.pop(); groups[0]["az0"] = g["az0"] - 360
        groups[0]["lo"] = min(groups[0]["lo"], g["lo"]); groups[0]["hi"] = max(groups[0]["hi"], g["hi"])
    merged = []
    for g in groups:   # merge areas whose padded ranges touch
        if merged and merged[-1]["flag"] == g["flag"] and g["az0"] - merged[-1]["az1"] <= 2 * pad_az + 0.5:
            m = merged[-1]; m["az1"] = g["az1"]; m["lo"] = min(m["lo"], g["lo"]); m["hi"] = max(m["hi"], g["hi"])
        else:
            merged.append(dict(g))
    boxes = []
    for g in merged:
        a0, a1 = g["az0"] - pad_az, g["az1"] + 0.5 + pad_az
        if g["flag"] == "clear_to_bottom":
            h0, h1 = max(floor, g["lo"] - extra_below), g["lo"] + 1.0
            label = f"Az {a0 % 360:.0f}-{a1 % 360:.0f} below {g['lo']:.0f}°"
        else:
            h0, h1 = g["hi"] - 1.0, min(ceiling, g["hi"] + extra_above)
            label = f"Az {a0 % 360:.0f}-{a1 % 360:.0f} above {g['hi']:.0f}°"
        if h1 > h0:
            boxes.append((a0, a1, h0, h1, label))
    return boxes


def steep_edge_boxes(az, lines, frame_w, edge_frames=1.0, edge_slope=3.0, min_jump=2.0,
                     margin=1.0, floor=0.0, ceiling=85.0):
    """Fill-in areas around the steep rises and drops of the measured horizon line(s).

    For every run of the line steeper than `edge_slope` (deg altitude per deg azimuth) whose
    total jump is at least `min_jump` degrees, an area is made that spans the jump in altitude
    (plus `margin`) and reaches `edge_frames` frame widths (scaled by 1/cos(alt)) to both sides.
    Areas from the different lines that overlap are merged. Returns (az0, az1, alt0, alt1, label).
    """
    if edge_frames <= 0:
        return []
    az = np.asarray(az, float); n = len(az)
    step = float(np.median(np.diff(az))) if n > 1 else 0.5
    raw = []
    for v in lines:
        v = np.asarray(v, float)
        steep = np.abs(np.roll(v, -1) - v) / step >= edge_slope           # segment i -> i+1
        i = 0
        while i < n:
            if not steep[i]:
                i += 1; continue
            j = i
            while j + 1 < n and steep[j + 1] and np.sign(v[(j + 2) % n] - v[j + 1]) == np.sign(v[(i + 1) % n] - v[i]):
                j += 1                                          # same direction: one rise or one drop
            seg = v[np.arange(i, j + 2) % n]
            bot, top = float(seg.min()), float(seg.max())
            if top - bot >= min_jump:
                w = edge_frames * frame_w / max(math.cos(math.radians(min(top + margin, 80.0))), 0.17)
                raw.append([az[i] - w, az[(j + 1) % n] + (360 if j + 1 >= n else 0) + w,
                            max(floor, bot - margin), min(ceiling, top + margin)])
            i = j + 2
    raw.sort(key=lambda r: r[0])
    merged = []      # only identical edges (seen on both lines) are merged; frames shared by
    for r in raw:    # neighbouring areas are planned once by plan_fillin
        m = merged[-1] if merged else None
        if m and abs(r[0] - m[0]) < 1.0 and abs(r[1] - m[1]) < 1.0:
            m[0] = min(m[0], r[0]); m[1] = max(m[1], r[1]); m[2] = min(m[2], r[2]); m[3] = max(m[3], r[3])
        else:
            merged.append(list(r))
    return [(a0, a1, h0, h1, f"Steep edge Az {a0 % 360:.0f}-{a1 % 360:.0f} ({h0:.0f}-{h1:.0f}°)")
            for a0, a1, h0, h1 in merged if h1 > h0]


def frame_coverage(mos, alt, az, half_w, half_h):
    """Fraction of a frame footprint (alt/az box) that the mosaic already has data for."""
    M = mos["M"]; res = mos["res"]; a0 = mos["alt0"]
    NH, NA = M.shape
    hs = np.linspace(alt - half_h, alt + half_h, 9); zs = np.linspace(az - half_w, az + half_w, 13)
    r = np.floor((hs - a0) / res).astype(int); c = np.floor((zs % 360) / 360 * NA).astype(int) % NA
    ok = (r >= 0) & (r < NH)
    if not ok.any():
        return 0.0
    sub = M[np.ix_(r[ok], c)]
    return float(np.isfinite(sub).sum()) / (len(hs) * len(zs))


def order_color(t):
    """Colour for position t (0 = first frame, 1 = last): green → yellow → magenta."""
    stops = [(0.0, (0, 170, 80)), (0.5, (240, 200, 0)), (1.0, (220, 0, 160))]
    t = min(max(t, 0.0), 1.0)
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t <= t1:
            u = (t - t0) / (t1 - t0)
            return "#%02x%02x%02x" % tuple(int(a + (b - a) * u) for a, b in zip(c0, c1))
    return "#dc00a0"


def draw_order(canvas, plan, X, Y, tag="order", numbers="rows", flip=True, width=2):
    """Draw the shooting order of plan['frames'] on a Tk canvas: a colour-graded path with arrows,
    frame numbers (at row starts or on every frame), Start/End markers and the meridian flip."""
    fr = plan["frames"]; n = len(fr)
    if not n:
        return
    pts = [(X(a), Y(alt)) for alt, a, _ in fr]
    step = max(1, n // 150)
    for i in range(n - 1):
        (xa, ya), (xb, yb) = pts[i], pts[i + 1]
        col = order_color(i / max(n - 1, 1))
        if abs(xb - xa) > abs(X(180) - X(0)):        # wraps over 0°/360°: no line across the image
            continue
        canvas.create_line(xa, ya, xb, yb, fill=col, width=width, tags=tag,
                           arrow="last" if i % step == 0 else "none", arrowshape=(7, 8, 3))
    idx = range(n) if numbers == "all" else []
    if numbers == "rows":
        fi = 0; idx = []
        for _, azs in plan["grid"]:
            if azs:
                idx.append(fi)
            fi += len(azs)
    for i in idx:
        x, y = pts[i]
        canvas.create_text(x, y - 9, text=str(i + 1), font=("Segoe UI", 8, "bold"),
                           fill=order_color(i / max(n - 1, 1)), tags=tag)
    for i, txt in ((0, "Start"), (n - 1, "End")):
        x, y = pts[i]
        canvas.create_oval(x - 5, y - 5, x + 5, y + 5, outline="black", width=2,
                           fill=order_color(i / max(n - 1, 1)), tags=tag)
        canvas.create_text(x + 8, y + 10, text=txt, anchor="w", font=("Segoe UI", 9, "bold"), tags=tag)
    if flip:
        for i in range(1, n):
            if fr[i - 1][1] < 180 <= fr[i][1]:
                x, y = pts[i]
                canvas.create_text(x, y + 14, text="↻ meridian flip", font=("Segoe UI", 8),
                                   fill="#b00060", tags=tag)
                break


def plan_fillin(p, boxes, mos=None, max_cover=0.97):
    """Plan frames covering the given boxes (az0, az1, alt0, alt1, label).

    With `mos`, frames of 'Steep edge' areas whose footprint the mosaic already covers
    (more than `max_cover`) are left out, so only the gaps near the horizon line are shot.
    """
    w, h = p.fov
    fw, fh = _effective_fov(w, h, p.rot_tolerance) if p.rotator else (w, h)
    ov = p.overlap
    east = [b for b in boxes if ((b[0] + b[1]) / 2) % 360 < 180]
    west = [b for b in boxes if ((b[0] + b[1]) / 2) % 360 >= 180]
    order = sorted(east, key=lambda b: (b[0] + b[1]) / 2 % 360) + \
        sorted(west, key=lambda b: (b[0] + b[1]) / 2 % 360) if p.meridian_split else \
        sorted(boxes, key=lambda b: (b[0] + b[1]) / 2 % 360)
    grid = []; labels = []; done = []          # done: centres of frames already planned

    def covered(alt, az):
        for h, z in done:
            if abs(alt - h) < fh * (0.5 - ov / 2) and \
                    abs((az - z + 180) % 360 - 180) < fw / max(math.cos(math.radians(h)), 1e-3) * (0.5 - ov / 2):
                return True
        return False

    for a0, a1, h0, h1, label in order:
        span_h = h1 - h0
        nr = 1 if span_h <= fh else math.ceil((span_h - fh) / (fh * (1 - ov))) + 1
        rows = [(h0 + h1) / 2] if nr == 1 else [h0 + fh / 2 + i * (span_h - fh) / (nr - 1) for i in range(nr)]
        for ri, alt in enumerate(rows):
            cov = fw / max(math.cos(math.radians(max(alt - fh / 2, 0))), 1e-3)
            span_a = a1 - a0
            n = 1 if span_a <= cov else math.ceil((span_a - cov) / (cov * (1 - ov))) + 1
            azs = [(a0 + a1) / 2] if n == 1 else [a0 + cov / 2 + k * (span_a - cov) / (n - 1) for k in range(n)]
            azs = [z % 360 for z in azs if not covered(alt, z % 360)]   # overlapping areas: shoot once
            if mos is not None and label.startswith("Steep edge"):
                hw = fw / max(math.cos(math.radians(alt)), 1e-3) / 2
                azs = [z for z in azs if frame_coverage(mos, alt, z, hw, fh / 2) <= max_cover]
            if not azs:
                continue
            if ri % 2:
                azs.reverse()
            done.extend((alt, z) for z in azs)
            grid.append((alt, azs)); labels.append(label)
    frames, rot_moves = assign_rotator(p, grid, "landscape")
    secs = len(frames) * (p.sec_per_frame + p.exposure) + rot_moves * p.sec_per_rot
    return dict(grid=grid, frames=frames, orientation="landscape", fw=fw, fh=fh, rot_moves=rot_moves,
                seconds=secs, lim=lambda az: (0, 0), labels=labels)


def plan_to_nina(p, plan, name="Horizon mapping"):
    ids = _Ids()
    if p.template:
        root, rid = _container(ids, "SequentialContainer", name, None, True)
        start = body = end = root; sid = bid = eid = rid
    else:
        root, rid = _container(ids, "SequenceRootContainer", name, None)
        start, sid = _container(ids, "StartAreaContainer", "Start", rid)
        body, bid = _container(ids, "TargetAreaContainer", "Targets", rid)
        end, eid = _container(ids, "EndAreaContainer", "End", rid)
        root["Items"]["$values"] = [start, body, end]
    if p.unpark:
        start["Items"]["$values"].append(_item(ids, "Telescope.UnparkScope", sid))
    fi = 0
    for gi, (alt, azs) in enumerate(plan["grid"]):
        side = "E" if azs and azs[0] < 180 else "W"
        title = f"Alt {alt:.1f}° {side} ({len(azs)} frames)"
        if plan.get("labels"):
            title = f"{plan['labels'][gi]} – " + title
        row, rowid = _container(ids, "SequentialContainer", title, bid, not p.collapsed)
        for az in azs:
            a, z, rot = plan["frames"][fi]; fi += 1
            row["Items"]["$values"].append(_slew(ids, rowid, z, a))
            if rot is not None:
                row["Items"]["$values"].append(_rotate(ids, rowid, rot))
            row["Items"]["$values"].append(_expose(ids, rowid, p))
        body["Items"]["$values"].append(row)
    if p.park:
        end["Items"]["$values"].append(_item(ids, "Telescope.ParkScope", eid))
    return root


def save_plan(p, plan, path):
    if p.template and not path.lower().endswith(".template.json"):
        path = (path[:-5] if path.lower().endswith(".json") else path) + ".template.json"
    root = plan_to_nina(p, plan, os.path.splitext(os.path.basename(path))[0].replace(".template", ""))
    with open(path, "w", encoding="utf-8") as f:
        json.dump(root, f, indent=2, ensure_ascii=False)
    return path


# =============================================================================
#  PROCESS – FITS reading
# =============================================================================

def read_fits_header(f):
    h = {}; n = 0
    while True:
        block = f.read(2880)
        if len(block) < 2880:
            raise ValueError("truncated FITS header")
        n += 2880
        for i in range(0, 2880, 80):
            card = block[i:i + 80].decode("ascii", "replace")
            if card.startswith("END") and card[3:].strip() == "":
                return h, n
            if card[8:10] == "= ":
                v = card[10:].strip()
                if v.startswith("'"):
                    v = v[1:].split("'")[0].strip()
                else:
                    v = v.split("/")[0].strip()
                    try:
                        v = float(v)
                    except ValueError:
                        pass
                h[card[:8].strip()] = v


def read_fits(path, data=True):
    with open(path, "rb") as f:
        h, _ = read_fits_header(f)
        if not data:
            return h, None
        nx, ny, bp = int(h["NAXIS1"]), int(h["NAXIS2"]), int(h["BITPIX"])
        dt = {8: ">u1", 16: ">i2", 32: ">i4", -32: ">f4", -64: ">f8"}[bp]
        a = np.frombuffer(f.read(nx * ny * abs(bp) // 8), dtype=dt).reshape(ny, nx).astype(np.float32)
    a = a * float(h.get("BSCALE", 1.0)) + float(h.get("BZERO", 0.0))
    return h, a


def reduce_frame(a, bayer, target_w=400):
    """Return luminance image downsampled to ~target_w wide and the total binning factor."""
    if bayer:
        ny, nx = a.shape; a = a[:ny // 2 * 2, :nx // 2 * 2]
        a = (a[0::2, 0::2] + a[0::2, 1::2] + a[1::2, 0::2] + a[1::2, 1::2]) / 4
        f0 = 2
    else:
        f0 = 1
    k = max(1, int(round(a.shape[1] / target_w)))
    ny, nx = a.shape[0] // k * k, a.shape[1] // k * k
    a = a[:ny, :nx].reshape(ny // k, k, nx // k, k).mean((1, 3))
    return a.astype(np.float32), f0 * k


def find_fits(folders, recursive=True):
    out = []
    for fo in folders:
        if os.path.isfile(fo):
            out.append(fo); continue
        for root, dirs, files in os.walk(fo):
            for fn in files:
                if fn.lower().endswith((".fits", ".fit", ".fts")):
                    out.append(os.path.join(root, fn))
            if not recursive:
                break
    return sorted(set(out))


# =============================================================================
#  PROCESS – mosaic
# =============================================================================

class Frame:
    __slots__ = ("path", "alt", "az", "rot", "pier", "scale", "img", "t")


def frame_delta(fr, K, lat):
    """Angle (deg) of the frame's 'up' relative to the local vertical."""
    flip = 180.0 if (isinstance(fr.pier, str) and fr.pier.lower().startswith("east")) else 0.0
    return fr.rot + K - parallactic_angle(fr.alt, fr.az, lat) + flip


def sample_frame(img, fr, delta, scale, AZ, ALT):
    """Bilinear sample of a frame at canvas az/alt (deg arrays)."""
    ny, nx = img.shape
    h0 = math.radians(fr.alt)
    A = np.radians(AZ - fr.az); H = np.radians(ALT)
    cosc = math.sin(h0) * np.sin(H) + math.cos(h0) * np.cos(H) * np.cos(A)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.degrees(np.cos(H) * np.sin(A) / cosc)
        y = np.degrees((math.cos(h0) * np.sin(H) - math.sin(h0) * np.cos(H) * np.cos(A)) / cosc)
    d = math.radians(delta)
    u = x * math.cos(d) + y * math.sin(d); w = -x * math.sin(d) + y * math.cos(d)
    col = nx / 2 + u / scale; row = ny / 2 - w / scale
    ok = (cosc > 0) & (col >= 0) & (col < nx - 1) & (row >= 0) & (row < ny - 1)
    out = np.full(AZ.shape, np.nan, np.float32)
    c = col[ok]; r = row[ok]; c0 = c.astype(int); r0 = r.astype(int); fc = c - c0; fr_ = r - r0
    out[ok] = (img[r0, c0] * (1 - fc) * (1 - fr_) + img[r0, c0 + 1] * fc * (1 - fr_)
               + img[r0 + 1, c0] * (1 - fc) * fr_ + img[r0 + 1, c0 + 1] * fc * fr_)
    return out


def default_cache_dir():
    base = os.environ.get("LOCALAPPDATA") or os.path.join(os.path.expanduser("~"), ".cache")
    return os.path.join(base, "HorizonMapper", "cache")


def _cache_key(path):
    st = os.stat(path)
    return hashlib.sha1(f"{os.path.abspath(path)}|{st.st_size}|{st.st_mtime_ns}|v1".encode()).hexdigest()


def _load_one(path, cache_dir):
    """Return (meta dict, reduced image) or (None, reason). Uses/creates a cache entry."""
    cf = None
    if cache_dir:
        try:
            cf = os.path.join(cache_dir, _cache_key(path) + ".npz")
            if os.path.exists(cf):
                z = np.load(cf, allow_pickle=False)
                meta = json.loads(str(z["meta"]))
                return meta, z["img"].astype(np.float32)
        except Exception:
            cf = None
    h, a = read_fits(path)
    if "CENTALT" not in h or "CENTAZ" not in h:
        return None, "no CENTALT/CENTAZ (mount not connected?)"
    if str(h.get("IMAGETYP", "LIGHT")).upper() not in ("LIGHT", "SNAPSHOT"):
        return None, f"image type {h.get('IMAGETYP')}"
    pix = float(h.get("XPIXSZ", 0) or 0); foc = float(h.get("FOCALLEN", 0) or 0)
    if pix <= 0 or foc <= 0:
        return None, "XPIXSZ/FOCALLEN missing in header"
    bayer = bool(str(h.get("BAYERPAT", "")).strip())
    img, k = reduce_frame(a, bayer)
    meta = dict(alt=float(h["CENTALT"]), az=float(h["CENTAZ"]),
                rot=float(h.get("ROTATOR", h.get("ROTATANG", 0.0)) or 0.0),
                pier=str(h.get("PIERSIDE", "")), t=str(h.get("DATE-LOC", "")),
                scale=206.265 * pix / foc / 3600.0 * k,
                lat=float(h["SITELAT"]) if "SITELAT" in h else None)
    if cf:
        try:
            os.makedirs(cache_dir, exist_ok=True)
            tmp = cf + ".tmp.npz"
            np.savez(tmp, img=np.clip(img, -65000, 65000).astype(np.float16), meta=json.dumps(meta))
            os.replace(tmp, cf)
        except Exception:
            pass
    return meta, img


def load_frames(paths, log, progress, bias=None, cache_dir=None, workers=None):
    """Read and reduce all frames (in parallel, using the cache when possible)."""
    frames = []; skipped = 0; lat = None; done = 0; cached_before = 0
    if cache_dir:
        try:
            names = set(os.listdir(cache_dir))
            cached_before = sum(1 for p in paths if _cache_key(p) + ".npz" in names)
        except Exception:
            cached_before = 0
    workers = workers or min(8, (os.cpu_count() or 4))
    results = [None] * len(paths)

    def job(i):
        try:
            return i, _load_one(paths[i], cache_dir)
        except Exception as e:
            return i, (None, str(e))

    with ThreadPoolExecutor(max_workers=workers) as ex:
        for i, res in ex.map(job, range(len(paths))):
            results[i] = res; done += 1
            if done % 20 == 0 or done == len(paths):
                progress(done / max(len(paths), 1) * 0.35, f"Reading frames {done}/{len(paths)}")
    for pth, (meta, img) in zip(paths, results):
        if meta is None:
            skipped += 1
            log(f"  skipped {os.path.basename(pth)}: {img}")
            continue
        fr = Frame(); fr.path = pth
        fr.alt, fr.az, fr.rot, fr.pier, fr.t, fr.scale = (meta["alt"], meta["az"], meta["rot"],
                                                          meta["pier"], meta["t"], meta["scale"])
        fr.img = img
        if lat is None and meta.get("lat") is not None:
            lat = meta["lat"]
        frames.append(fr)
    if cache_dir and cached_before:
        log(f"{cached_before} of {len(paths)} frames came from the cache")
    # bias: lowest percentile over all frames (dark obstructions approach the bias)
    if frames:
        b = bias if bias is not None else min(float(np.percentile(f.img, 0.1)) for f in frames) - 5
        for f in frames:
            f.img = np.log(np.clip(f.img - b, 1, None)).astype(np.float32)
        log(f"Loaded {len(frames)} frames ({skipped} skipped); bias level used {b:.0f} ADU")
        if bias is None and b > 3000:
            log("  Warning: no dark frames found, so the automatic bias level may be too high. "
                "Enter the camera's bias level (e.g. ~500 ADU for offset 50) in Settings.")
    return frames, lat


def clear_cache(cache_dir):
    n = 0
    if cache_dir and os.path.isdir(cache_dir):
        for fn in os.listdir(cache_dir):
            if fn.endswith(".npz"):
                try:
                    os.remove(os.path.join(cache_dir, fn)); n += 1
                except OSError:
                    pass
    return n


def _highpass(img):
    h = box2d(img, 3) - box2d(img, 31 if min(img.shape) > 62 else 9)
    return (h / (h.std() + 1e-6)).astype(np.float32)


def calibrate_rotation(frames, lat, guess, log, search_full=False):
    """Find K (deg) so that frame up = ROTATOR + K (+180 on east pier) best matches overlaps."""
    tex = np.array([float(np.std(_highpass(f.img))) * float(np.std(f.img)) for f in frames])
    best_i = int(np.argmax(tex))
    c = frames[best_i]
    near = [i for i, f in enumerate(frames)
            if abs((f.az - c.az + 180) % 360 - 180) * math.cos(math.radians(f.alt)) < 7 and abs(f.alt - c.alt) < 5]
    near = sorted(near, key=lambda i: -tex[i])[:30]
    if len(near) < 3:
        log("Rotation calibration: not enough overlapping frames, using the given offset.")
        return guess
    H = {i: _highpass(frames[i].img) for i in near}
    azs = [frames[i].az for i in near]; alts = [frames[i].alt for i in near]
    a0 = c.az - 9 / max(math.cos(math.radians(c.alt)), 0.2); a1 = c.az + 9 / max(math.cos(math.radians(c.alt)), 0.2)
    AZ, ALT = np.meshgrid(np.arange(a0, a1, 0.05), np.arange(min(alts) - 2, max(alts) + 2, 0.05))

    def score(K):
        S = np.stack([sample_frame(H[i], frames[i], frame_delta(frames[i], K, lat), frames[i].scale, AZ, ALT)
                      for i in near])
        V = ~np.isnan(S); num = den = 0.0
        for a in range(len(S)):
            for b in range(a + 1, len(S)):
                m = V[a] & V[b]; k = int(m.sum())
                if k < 200:
                    continue
                r = np.corrcoef(S[a][m], S[b][m])[0, 1]
                if np.isfinite(r):
                    num += r * k; den += k
        return num / den if den else -1

    if search_full:
        cands = list(range(0, 360, 10))
    else:
        g = guess % 360
        cands = [g + d for d in (-8, -4, 0, 4, 8)] + [g + 180 + d for d in (-8, -4, 0, 4, 8)]
    res = sorted(((score(K), K) for K in cands), reverse=True)
    s0, K0 = res[0]
    fine = sorted(((score(K0 + d), K0 + d) for d in np.arange(-3, 3.01, 1.0)), reverse=True)
    s1, K1 = fine[0]
    log(f"Rotation calibration: best offset {K1 % 360:.1f}° (overlap correlation {s1:.3f}; "
        f"given {guess % 360:.1f}°)")
    return K1 % 360


def build_mosaic(frames, lat, K, res, log, progress):
    alts_all = [f.alt for f in frames]
    alt0 = math.floor(min(alts_all) - 2.5); alt1 = math.ceil(max(alts_all) + 2.5)
    alt0 = max(alt0, -5); alt1 = min(alt1, 90)
    NA = int(round(360 / res)); NH = int(round((alt1 - alt0) / res))
    alts = alt0 + np.arange(NH) * res
    patches = []
    for k, f in enumerate(frames):
        half = 0.5 * math.hypot(f.img.shape[1], f.img.shape[0]) * f.scale + 0.1
        i0 = max(int((f.alt - half - alt0) / res), 0); i1 = min(int((f.alt + half - alt0) / res) + 1, NH)
        top = min(f.alt + half, 89.0)
        ra = min(half / max(math.cos(math.radians(top)), 0.02), 180)
        j0 = int(math.floor((f.az - ra) / res)); j1 = int(math.ceil((f.az + ra) / res)) + 1
        if j1 - j0 >= NA:
            j0, j1 = 0, NA
        J = np.arange(j0, j1) % NA
        AZ, ALT = np.meshgrid(J * res, alts[i0:i1])
        v = sample_frame(f.img, f, frame_delta(f, K, lat), f.scale, AZ, ALT)
        patches.append((i0, i1, J, v))
        if k % 20 == 0:
            progress(0.45 + 0.25 * k / len(frames), f"Placing frames {k}/{len(frames)}")
    # brightness offsets between neighbours
    n = len(frames)
    cen = np.array([[f.az, f.alt] for f in frames])
    rows = []; rhs = []
    for a in range(n):
        dz = np.abs((cen[a + 1:, 0] - cen[a, 0] + 180) % 360 - 180) * math.cos(math.radians(cen[a, 1]))
        da = np.abs(cen[a + 1:, 1] - cen[a, 1])
        cand = np.where((dz < 6) & (da < 4))[0] + a + 1
        ia0, ia1, Ja, va = patches[a]
        posa = {j: x for x, j in enumerate(Ja)}
        for c in cand:
            ic0, ic1, Jc, vc = patches[c]
            i0, i1 = max(ia0, ic0), min(ia1, ic1)
            if i1 <= i0:
                continue
            common = [x for x, j in enumerate(Jc) if j in posa]
            if len(common) < 5:
                continue
            xa = [posa[Jc[x]] for x in common]
            P = va[i0 - ia0:i1 - ia0][:, xa]; Q = vc[i0 - ic0:i1 - ic0][:, common]
            m = ~np.isnan(P) & ~np.isnan(Q)
            if m.sum() < 200:
                continue
            w = math.sqrt(m.sum())
            r = np.zeros(n); r[a] = w; r[c] = -w
            rows.append(r); rhs.append(float(np.median(P[m] - Q[m])) * w)
        if a % 50 == 0:
            progress(0.70 + 0.15 * a / n, "Matching brightness")
    rows.append(np.ones(n) * 10); rhs.append(0.0)
    off = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0]
    log(f"Brightness matched using {len(rows) - 1} overlapping pairs")
    acc = np.zeros((NH, NA), np.float32); ws = np.zeros((NH, NA), np.float32)
    for k, (i0, i1, J, v) in enumerate(patches):
        m = ~np.isnan(v)
        sub = acc[i0:i1][:, J]; sw = ws[i0:i1][:, J]
        sub[m] += v[m] - off[k]; sw[m] += 1
        acc[i0:i1, J] = sub; ws[i0:i1, J] = sw
    M = np.where(ws > 0, acc / np.maximum(ws, 1), np.nan).astype(np.float32)
    return dict(M=M, alt0=alt0, res=res)


def detect_horizon(mos, solid_frac=0.6, thin_frac=0.85, step=0.5, max_branch_gap=6.0, log=print):
    M = mos["M"]; res = mos["res"]; alt0 = mos["alt0"]
    NH, NA = M.shape; alts = alt0 + np.arange(NH) * res
    cov = ~np.isnan(M)
    S = np.where(box2d(cov.astype(np.float32), 5) > 0.5,
                 box2d(np.nan_to_num(M, nan=0.0), 5) / np.maximum(box2d(cov.astype(np.float32), 5), 1e-6), np.nan)
    win = max(3, int(round(0.3 / res)))
    out = []
    for a in np.arange(0, 360, step):
        j0 = int(a / res); j1 = int((a + step) / res)
        col = S[:, j0:j1]; c = ~np.isnan(col)
        if c.sum() < 50:
            out.append(dict(az=float(a), solid=None, thin=None, lo=None, hi=None, flag="nodata")); continue
        rws = np.where(c.any(1))[0]; lo, hi = float(alts[rws.min()]), float(alts[rws.max()])
        sky = np.nanpercentile(col, 97)

        def top(mask, thr):
            f = box1d(mask.mean(1).astype(float), win); cf = box1d(c.mean(1).astype(float), win)
            frac = np.where(cf > 0.2, f / np.maximum(cf, 1e-6), 0)
            idx = np.where(frac >= thr)[0]
            return float(alts[idx.max()]) if len(idx) else None

        s = top((col < sky + math.log(solid_frac)) & c, 0.5)
        t = top((col < sky + math.log(thin_frac)) & c, 0.3)
        flag = ""
        if s is None and t is None:
            flag = "clear_to_bottom"
        elif max(x for x in (s, t) if x is not None) >= hi - 0.6:
            flag = "blocked_to_top"
        out.append(dict(az=float(a), solid=s, thin=t, lo=lo, hi=hi, flag=flag))
    az = np.array([o["az"] for o in out])

    def arr(key):
        v = np.array([np.nan if o[key] is None else o[key] for o in out])
        ok = ~np.isnan(v)
        if not ok.any():
            return np.full(len(az), np.nan)
        return np.interp(az, az[ok], v[ok], period=360)

    solid = arr("solid"); thin = arr("thin")
    for i, o in enumerate(out):
        if o["flag"] == "clear_to_bottom":
            solid[i] = thin[i] = o["lo"]
    solid_f = median_wrap(solid, 3)
    branches = np.maximum(solid_f, median_wrap(np.where(thin - solid < max_branch_gap, thin, solid), 3))
    nodata = [o["az"] for o in out if o["flag"] == "nodata"]
    if nodata:
        log(f"No data at {len(nodata)} azimuth samples (interpolated): "
            f"{nodata[0]:.1f}..{nodata[-1]:.1f}°")
    unsure = [(o["az"], o["flag"]) for o in out if o["flag"] in ("clear_to_bottom", "blocked_to_top")]
    measured = np.array([o["flag"] != "nodata" for o in out])
    return dict(az=az, solid=solid_f, branches=branches, details=out, uncertain=unsure, measured=measured)


def widen_steep_edges(az, alt, side_deg=0.0, steep_slope=3.0):
    """Extend the horizon sideways at steep edges.

    Where the horizon rises or falls faster than `steep_slope` degrees of altitude per degree of
    azimuth, the higher side is extended `side_deg` degrees sideways (a local maximum filter).
    Gentle slopes and flat stretches are left unchanged.
    """
    alt = np.asarray(alt, float)
    if side_deg <= 0 or len(alt) < 3:
        return alt.copy()
    az = np.asarray(az, float)
    step = float(np.median(np.diff(np.sort(az)))) or 0.5
    k = max(1, int(round(side_deg / step)))
    # steep segments between neighbouring samples (wrap-around included)
    slope = np.abs(np.roll(alt, -1) - alt) / step
    steep = slope >= steep_slope
    out = alt.copy()
    for i in np.where(steep)[0]:
        j = (i + 1) % len(alt)
        hi_idx, hi = (i, alt[i]) if alt[i] >= alt[j] else (j, alt[j])
        # spread the higher value k samples into the lower side
        direction = 1 if hi_idx == i else -1
        for d in range(1, k + 1):
            n = (hi_idx + direction * d) % len(alt)
            if out[n] < hi:
                out[n] = hi
    return out


def fill_unmeasured(det, old_pts):
    """Use the old horizon wherever nothing was imaged (partial-sky runs)."""
    if old_pts is None or det["measured"].all():
        return det
    d = dict(det)
    oldv = hrz_array(old_pts, det["az"])
    m = ~det["measured"]
    d["solid"] = np.where(m, oldv, det["solid"]); d["branches"] = np.where(m, oldv, det["branches"])
    return d


# =============================================================================
#  Images
# =============================================================================

def _font(size):
    for name in ("arial.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            pass
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        return ImageFont.load_default()


def render_mosaic(mos, lines=(), scale=2, legend=True):
    """lines: list of (label, az_array, alt_array, rgb)."""
    M = mos["M"]; res = mos["res"]; alt0 = mos["alt0"]; NH, NA = M.shape
    v = M[np.isfinite(M)]
    lo, hi = (np.percentile(v, [1, 99.5]) if v.size else (0, 1))
    g = np.clip((M - lo) / max(hi - lo, 1e-6), 0, 1)
    g8 = (np.nan_to_num(g, nan=0) * 255).astype(np.uint8)
    rgb = np.dstack([g8, g8, g8]); rgb[np.isnan(M)] = (28, 28, 36)
    im = Image.fromarray(rgb[::-1]).resize((max(NA // scale, 1), max(NH // scale, 1)), Image.LANCZOS)
    PL, PB, PT = 70, 90 if legend else 40, 34
    can = Image.new("RGB", (im.width + PL + 20, im.height + PB + PT), "white"); can.paste(im, (PL, PT))
    d = ImageDraw.Draw(can); f = _font(24)
    X = lambda a: PL + a / res / scale
    Y = lambda h: PT + im.height - 1 - (h - alt0) / res / scale
    for a in range(0, 361, 15):
        d.line([X(a), PT, X(a), PT + im.height], fill=(70, 130, 200))
        if a % 30 == 0:
            d.text((X(a) - 16, PT + im.height + 6), f"{a}°", fill="black", font=f)
    for hh in range(int(math.ceil(alt0 / 10) * 10), int(alt0 + NH * res) + 1, 10):
        d.line([PL, Y(hh), PL + im.width, Y(hh)], fill=(70, 130, 200))
        d.text((6, Y(hh) - 12), f"{hh}°", fill="black", font=f)
    for lab, a in (("N", 0), ("E", 90), ("S", 180), ("W", 270), ("N", 360)):
        d.text((X(a) - 8, 4), lab, fill="black", font=f)
    top = alt0 + NH * res
    for label, az, alt, col in lines:
        o = np.argsort(az); az = np.asarray(az)[o]; alt = np.asarray(alt)[o]
        az = np.concatenate([az, [az[0] + 360]]); alt = np.concatenate([alt, [alt[0]]])
        d.line([(X(a), Y(min(max(h, alt0), top))) for a, h in zip(az, alt)], fill=col, width=3)
    if legend:
        x = PL; y = PT + im.height + 46
        for label, _, _, col in lines:
            d.rectangle([x, y + 6, x + 40, y + 20], fill=col); d.text((x + 50, y), label, fill="black", font=f)
            x += int(d.textlength(label, font=f)) + 110
    return can


# =============================================================================
#  GUI
# =============================================================================

def run_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    try:
        from PIL import ImageTk
    except ImportError:
        ImageTk = None

    root = tk.Tk()
    root.title(f"{APP_NAME}  v{APP_VERSION}")
    root.geometry("1280x860")
    nb = ttk.Notebook(root); nb.pack(fill="both", expand=True)

    def labeled(parent, r, text, var, width=10, col=0, tip=None):
        ttk.Label(parent, text=text).grid(row=r, column=col, sticky="w", padx=4, pady=2)
        e = ttk.Entry(parent, textvariable=var, width=width); e.grid(row=r, column=col + 1, sticky="w", padx=4, pady=2)
        if tip:
            ttk.Label(parent, text=tip, foreground="#666").grid(row=r, column=col + 2, sticky="w")
        return e

    # ------------------------------------------------------------------ PLAN
    plan_tab = ttk.Frame(nb); nb.add(plan_tab, text="1. Plan sequence")
    left = ttk.Frame(plan_tab); left.pack(side="left", fill="y", padx=6, pady=6)
    right = ttk.Frame(plan_tab); right.pack(side="left", fill="both", expand=True, padx=6, pady=6)
    P = PlanParams()
    V = {}
    current_hz = tk.StringVar(value="")   # shared: Plan tab horizon file -> Build tab comparison

    def var(name, value, kind=tk.StringVar):
        V[name] = kind(value=value); return V[name]

    box = ttk.LabelFrame(left, text="Camera & telescope"); box.pack(fill="x", pady=3)
    labeled(box, 0, "Sensor width (px)", var("sensor_w_px", P.sensor_w_px))
    labeled(box, 1, "Sensor height (px)", var("sensor_h_px", P.sensor_h_px))
    labeled(box, 2, "Pixel size (µm)", var("pixel_um", P.pixel_um))
    labeled(box, 3, "Focal length (mm)", var("focal_mm", P.focal_mm))
    fov_lbl = ttk.Label(box, text=""); fov_lbl.grid(row=4, column=0, columnspan=3, sticky="w", padx=4)

    def from_fits():
        f = filedialog.askopenfilename(title="Any FITS frame from this setup",
                                       filetypes=[("FITS", "*.fits *.fit *.fts"), ("All", "*.*")])
        if not f:
            return
        h, _ = read_fits(f, data=False)
        b = float(h.get("XBINNING", 1) or 1)
        V["sensor_w_px"].set(int(float(h["NAXIS1"]) * b)); V["sensor_h_px"].set(int(float(h["NAXIS2"]) * b))
        V["pixel_um"].set(round(float(h.get("XPIXSZ", 0)) / b, 3)); V["focal_mm"].set(h.get("FOCALLEN", ""))
        if "SITELAT" in h:
            V["lat"].set(round(float(h["SITELAT"]), 4))
        if "OBJCTROT" in h and "ROTATOR" in h:
            V["rot_offset"].set(round((float(h["OBJCTROT"]) - float(h["ROTATOR"]) + 180) % 360 - 180, 2))
            messagebox.showinfo(APP_NAME, "Rotator offset taken from OBJCTROT - ROTATOR.\n"
                                          "This is only valid if that frame was plate-solved/rotated by NINA.")
        update_fov()
    ttk.Button(box, text="Read from a FITS file…", command=from_fits).grid(row=5, column=0, columnspan=2, sticky="w", padx=4, pady=3)

    box = ttk.LabelFrame(left, text="Site & area to image"); box.pack(fill="x", pady=3)
    labeled(box, 0, "Latitude (°)", var("lat", P.lat))
    V["mode"] = tk.StringVar(value="horizon")
    ttk.Radiobutton(box, text="Fixed altitude band", variable=V["mode"], value="band").grid(row=1, column=0, sticky="w")
    ttk.Radiobutton(box, text="Around current horizon", variable=V["mode"], value="horizon").grid(row=1, column=1, sticky="w")
    labeled(box, 2, "Band: alt min (°)", var("alt_min", P.alt_min))
    labeled(box, 3, "Band: alt max (°)", var("alt_max", P.alt_max))
    hzf = ttk.Frame(box); hzf.grid(row=4, column=0, columnspan=3, sticky="w", padx=4, pady=(4, 0))
    ttk.Label(hzf, text="Current horizon (.hrz)").pack(anchor="w")
    hz_list = tk.Listbox(hzf, height=2, width=40); hz_list.pack(anchor="w")

    def sync_current():
        files = list(hz_list.get(0, "end"))
        if files:
            current_hz.set(files[0])

    def add_hz():
        for f in filedialog.askopenfilenames(title="Current horizon file", filetypes=[("Horizon", "*.hrz *.txt *.csv"), ("All", "*.*")]):
            if hz_list.size() < 2:
                hz_list.insert("end", f)
        sync_current()

    def clear_hz():
        files = list(hz_list.get(0, "end"))
        hz_list.delete(0, "end")
        if files and current_hz.get() == files[0]:
            current_hz.set("")
    hzb = ttk.Frame(hzf); hzb.pack(anchor="w", pady=2)
    ttk.Button(hzb, text="Add current horizon…", command=add_hz).pack(side="left")
    ttk.Button(hzb, text="Clear", command=clear_hz).pack(side="left", padx=4)
    ttk.Label(hzf, text="The first file is your current horizon; it is also used on the Build tab.\n"
                        "Optional second file (e.g. an all-sky estimate): the strip spans both.",
              foreground="#666", justify="left").pack(anchor="w")
    labeled(box, 7, "Margin below (°)", var("margin_below", P.margin_below))
    labeled(box, 8, "Margin above (°)", var("margin_above", P.margin_above))
    labeled(box, 9, "Lowest altitude (°)", var("floor", P.floor))
    labeled(box, 13, "Extra frames at steep edges", var("edge_frames", P.edge_frames), tip="0 = off")
    labeled(box, 14, "Steep = more than (°/°)", var("edge_slope", P.edge_slope), tip="alt per az")
    labeled(box, 10, "Azimuth start (°)", var("az_start", P.az_start))
    labeled(box, 11, "Azimuth end (°)", var("az_end", P.az_end))
    labeled(box, 12, "Overlap (0-1)", var("overlap", P.overlap))

    box = ttk.LabelFrame(left, text="Rotator"); box.pack(fill="x", pady=3)
    V["rotator"] = tk.BooleanVar(value=True)
    ttk.Checkbutton(box, text="Keep frames level with rotator", variable=V["rotator"]).grid(row=0, column=0, columnspan=2, sticky="w")
    labeled(box, 1, "Offset (sky - mech, °)", var("rot_offset", P.rot_offset))
    labeled(box, 2, "Skip moves below (°)", var("rot_tolerance", P.rot_tolerance))
    labeled(box, 3, "Direction (+1/-1)", var("rot_sign", P.rot_sign))
    V["orientation"] = tk.StringVar(value="auto")
    ttk.Label(box, text="Orientation").grid(row=4, column=0, sticky="w", padx=4)
    ttk.Combobox(box, textvariable=V["orientation"], values=["auto", "landscape", "portrait"], width=10, state="readonly").grid(row=4, column=1, sticky="w")

    box = ttk.LabelFrame(left, text="Exposure & sequence"); box.pack(fill="x", pady=3)
    labeled(box, 0, "Exposure (s)", var("exposure", P.exposure))
    labeled(box, 1, "Gain (-1 = default)", var("gain", P.gain))
    labeled(box, 2, "Offset (-1 = default)", var("offset", P.offset))
    labeled(box, 3, "Binning", var("binning", P.binning))
    V["meridian_split"] = tk.BooleanVar(value=True); V["unpark"] = tk.BooleanVar(value=True)
    V["park"] = tk.BooleanVar(value=True); V["template"] = tk.BooleanVar(value=False); V["collapsed"] = tk.BooleanVar(value=True)
    ttk.Checkbutton(box, text="East half first, then west (one meridian flip)", variable=V["meridian_split"]).grid(row=4, column=0, columnspan=3, sticky="w")
    ttk.Checkbutton(box, text="Unpark at start", variable=V["unpark"]).grid(row=5, column=0, sticky="w")
    ttk.Checkbutton(box, text="Park at end", variable=V["park"]).grid(row=5, column=1, sticky="w")
    ttk.Checkbutton(box, text="Save as template (else full sequence)", variable=V["template"]).grid(row=6, column=0, columnspan=3, sticky="w")
    ttk.Checkbutton(box, text="Rows collapsed (use Load sequence)", variable=V["collapsed"]).grid(row=7, column=0, columnspan=3, sticky="w")
    labeled(box, 8, "Seconds per frame", var("sec_per_frame", P.sec_per_frame))
    labeled(box, 9, "Seconds per rotator move", var("sec_per_rot", P.sec_per_rot))

    btns = ttk.Frame(left); btns.pack(fill="x", pady=6)
    plan_info = ttk.Label(right, text="", font=("Segoe UI", 10, "bold")); plan_info.pack(anchor="w")
    canvas = tk.Canvas(right, bg="white", height=560); canvas.pack(fill="both", expand=True)
    state = {"plan": None, "params": None}

    def collect(require_area=True):
        p = PlanParams()
        num = ["sensor_w_px", "sensor_h_px", "pixel_um", "focal_mm", "lat", "alt_min", "alt_max", "margin_below",
               "margin_above", "floor", "az_start", "az_end", "overlap", "rot_offset", "rot_tolerance", "exposure",
               "sec_per_frame", "sec_per_rot", "edge_frames", "edge_slope"]
        for k in num:
            setattr(p, k, float(V[k].get()))
        for k in ("sensor_w_px", "sensor_h_px"):
            setattr(p, k, int(float(V[k].get())))
        for k in ("gain", "offset", "binning", "rot_sign"):
            setattr(p, k, int(float(V[k].get())))
        for k in ("rotator", "meridian_split", "unpark", "park", "template", "collapsed"):
            setattr(p, k, bool(V[k].get()))
        p.orientation = V["orientation"].get(); p.mode = V["mode"].get()
        p.horizon_files = list(hz_list.get(0, "end"))
        if p.mode == "horizon" and not p.horizon_files and require_area:
            raise ValueError("Add your current horizon file, or choose 'Fixed altitude band'.")
        return p

    def update_fov(*_):
        try:
            p = collect() if V["mode"].get() == "band" or hz_list.size() else PlanParams(
                sensor_w_px=int(float(V["sensor_w_px"].get())), sensor_h_px=int(float(V["sensor_h_px"].get())),
                pixel_um=float(V["pixel_um"].get()), focal_mm=float(V["focal_mm"].get()))
            w, h = p.fov
            fov_lbl.config(text=f"Field of view: {w:.2f}° × {h:.2f}°")
        except Exception:
            fov_lbl.config(text="")

    def draw_plan():
        canvas.delete("all")
        p, plan = state["params"], state["plan"]
        if not plan:
            return
        W = canvas.winfo_width() or 900; H = canvas.winfo_height() or 540
        x0, y0, x1, y1 = 50, 20, W - 15, H - 35
        amin, amax = 0, 90
        X = lambda a: x0 + (a % 360) / 360 * (x1 - x0)
        Y = lambda h: y1 - (h - amin) / (amax - amin) * (y1 - y0)
        for a in range(0, 361, 30):
            canvas.create_line(X(a) if a < 360 else x1, y0, X(a) if a < 360 else x1, y1, fill="#ddd")
            canvas.create_text(X(a) if a < 360 else x1, y1 + 14, text=f"{a}°")
        for h in range(0, 91, 10):
            canvas.create_line(x0, Y(h), x1, Y(h), fill="#ddd"); canvas.create_text(x0 - 22, Y(h), text=f"{h}°")
        az = np.arange(0, 360, 1.0)
        lo_hi = [plan["lim"](a) for a in az]
        poly = [(X(a), Y(hi)) for a, (lo, hi) in zip(az, lo_hi)] + [(X(a), Y(lo)) for a, (lo, hi) in reversed(list(zip(az, lo_hi)))]
        canvas.create_polygon(*[c for pt in poly for c in pt], fill="#d6e6ff", outline="")
        for alt, a, rot in plan["frames"]:
            cw = plan["fw"] / max(math.cos(math.radians(alt)), 0.05) / 2; ch = plan["fh"] / 2
            canvas.create_rectangle(X(a - cw), Y(alt + ch), X(a + cw), Y(alt - ch), outline="#4a7bd0")
        colors = ["#d00000", "#00a000"]
        for i, f in enumerate(p.horizon_files):
            pts = load_hrz(f)
            line = [c for a in az for c in (X(a), Y(hrz_alt(pts, a)))]
            canvas.create_line(*line, fill=colors[i % 2], width=2)
        canvas.create_line(x0, Y(p.floor), x1, Y(p.floor), fill="#e0b000", dash=(4, 3))
        if show_order.get():
            draw_order(canvas, plan, X, Y, flip=p.meridian_split)
            canvas.create_text(x1, y0 + 2, anchor="ne", font=("Segoe UI", 9),
                               text="Shooting order: green = first → yellow → magenta = last · numbers = first frame of each row")

    def preview():
        try:
            p = collect()
            plan = make_plan(p)
        except Exception as e:
            messagebox.showerror(APP_NAME, str(e)); return
        state["plan"], state["params"] = plan, p
        w, h = p.fov
        rows = len(plan["grid"])
        msg = (f"{len(plan['frames'])} frames in {rows} rows · {plan['orientation']} · "
               f"{plan['rot_moves']} rotator moves · ~{plan['seconds'] / 3600:.1f} h  "
               f"(FOV {w:.2f}°×{h:.2f}°)")
        plan_info.config(text=msg); draw_plan()

    def save():
        if not state["plan"]:
            preview()
        if not state["plan"]:
            return
        p = state["params"]
        ext = ".template.json" if p.template else ".json"
        f = filedialog.asksaveasfilename(defaultextension=".json", initialfile="Horizon mapping" + ext,
                                         filetypes=[("NINA sequence", "*.json")])
        if f:
            out = save_plan(p, state["plan"], f)
            messagebox.showinfo(APP_NAME, f"Saved:\n{out}\n\n" + (
                "Copy it to NINA's template folder and restart NINA." if p.template else
                "Open it in the Advanced Sequencer with 'Load sequence'."))

    def sun_check():
        try:
            lat = float(V["lat"].get()); lon = float(sun_lon.get())
        except ValueError:
            messagebox.showerror(APP_NAME, "Enter a longitude to check the sun."); return
        t0 = time.time(); lines = []
        for hrs in range(0, 13):
            alt, az = sun_altaz(t0 + hrs * 3600, lat, lon)
            lines.append(f"+{hrs:2d} h: sun alt {alt:5.1f}°, az {az:5.1f}°")
        messagebox.showinfo(APP_NAME, "Sun position (from now):\n" + "\n".join(lines) +
                            "\n\nKeep the lowest imaged altitude well above the sun when shooting by day.")

    ttk.Button(btns, text="Preview plan", command=preview).pack(side="left", padx=3)
    ttk.Button(btns, text="Save NINA sequence…", command=save).pack(side="left", padx=3)
    show_order = tk.BooleanVar(value=True)
    ttk.Checkbutton(btns, text="Show shooting order", variable=show_order,
                    command=lambda: draw_plan()).pack(side="left", padx=6)
    sun_lon = tk.StringVar(value="23.62")
    sb = ttk.Frame(left); sb.pack(fill="x")
    ttk.Label(sb, text="Longitude (°E)").pack(side="left", padx=4)
    ttk.Entry(sb, textvariable=sun_lon, width=8).pack(side="left")
    ttk.Button(sb, text="Sun check", command=sun_check).pack(side="left", padx=4)
    canvas.bind("<Configure>", lambda e: draw_plan())
    for k in ("sensor_w_px", "sensor_h_px", "pixel_um", "focal_mm"):
        V[k].trace_add("write", update_fov)
    update_fov()

    # ------------------------------------------------------------------ PROCESS
    proc = ttk.Frame(nb); nb.add(proc, text="2. Build horizon")
    pl = ttk.Frame(proc); pl.pack(side="left", fill="y", padx=6, pady=6)
    pr = ttk.Frame(proc); pr.pack(side="left", fill="both", expand=True, padx=6, pady=6)
    box = ttk.LabelFrame(pl, text="FITS folders (searched recursively)"); box.pack(fill="x", pady=3)
    folders = tk.Listbox(box, height=5, width=48); folders.pack(fill="x", padx=4, pady=2)
    fb = ttk.Frame(box); fb.pack(fill="x")
    ttk.Button(fb, text="Add folder…", command=lambda: (lambda d: d and folders.insert("end", d))(
        filedialog.askdirectory())).pack(side="left", padx=3)
    ttk.Button(fb, text="Remove", command=lambda: [folders.delete(i) for i in reversed(folders.curselection())]).pack(side="left")
    box = ttk.LabelFrame(pl, text="Settings"); box.pack(fill="x", pady=3)
    pv = {}
    def pvar(n, v): pv[n] = tk.StringVar(value=v); return pv[n]
    labeled(box, 0, "Rotator offset guess (°)", pvar("K", "-58.04"), tip="sky - mech")
    pv["auto"] = tk.BooleanVar(value=True)
    ttk.Checkbutton(box, text="Calibrate rotation from overlaps", variable=pv["auto"]).grid(row=1, column=0, columnspan=3, sticky="w")
    pv["full"] = tk.BooleanVar(value=False)
    ttk.Checkbutton(box, text="Full 360° search (offset unknown)", variable=pv["full"]).grid(row=2, column=0, columnspan=3, sticky="w")
    labeled(box, 3, "Latitude (°, blank = FITS)", pvar("lat", ""))
    labeled(box, 4, "Mosaic resolution (°)", pvar("res", "0.03"))
    labeled(box, 5, "Solid: darker than", pvar("solid", "0.60"), tip="× sky")
    labeled(box, 6, "Branches: darker than", pvar("thin", "0.85"), tip="× sky")
    labeled(box, 7, "Bias level (blank = auto)", pvar("bias", ""))
    pvar("side", "0"); pvar("steep", "3")   # shown in the toolbar above the image
    pv["cache"] = tk.BooleanVar(value=True)
    ttk.Checkbutton(box, text="Cache reduced frames (fast re-runs)", variable=pv["cache"]).grid(row=10, column=0, columnspan=2, sticky="w")
    ttk.Button(box, text="Clear cache", command=lambda: messagebox.showinfo(
        APP_NAME, f"Removed {clear_cache(default_cache_dir())} cached frames.")).grid(row=10, column=2, sticky="w")
    box = ttk.LabelFrame(pl, text="Current horizon (compare / fill gaps)"); box.pack(fill="x", pady=3)
    old_var = current_hz   # filled automatically from the Plan tab
    ttk.Entry(box, textvariable=old_var, width=40).pack(side="left", padx=4)
    ttk.Button(box, text="Browse…", command=lambda: old_var.set(filedialog.askopenfilename(
        filetypes=[("Horizon", "*.hrz *.txt"), ("All", "*.*")]) or old_var.get())).pack(side="left")
    bb = ttk.Frame(pl); bb.pack(fill="x", pady=4)
    pb = ttk.Progressbar(pl, length=340, maximum=1.0); pb.pack(fill="x", pady=4)
    status = ttk.Label(pl, text="Idle"); status.pack(anchor="w")
    logbox = tk.Text(pl, height=10, width=58, wrap="word"); logbox.pack(fill="both", expand=True, pady=3)

    # ---- editor (right side)
    tb = ttk.Frame(pr); tb.pack(fill="x")
    ed = {"target": tk.StringVar(value="none"), "brush": tk.StringVar(value="1.0"), "zoom": tk.StringVar(value="6")}
    ttk.Label(tb, text="Edit:").pack(side="left")
    for txt, val in (("off", "none"), ("solid line", "solid"), ("branches line", "branches")):
        ttk.Radiobutton(tb, text=txt, variable=ed["target"], value=val).pack(side="left", padx=2)
    ttk.Label(tb, text="   brush width (°)").pack(side="left")
    ttk.Entry(tb, textvariable=ed["brush"], width=5).pack(side="left")
    ttk.Label(tb, text="   zoom (px/°)").pack(side="left")
    zc = ttk.Combobox(tb, textvariable=ed["zoom"], values=["3", "6", "12", "20"], width=4, state="readonly")
    zc.pack(side="left")
    ttk.Button(tb, text="Undo", command=lambda: undo()).pack(side="left", padx=6)
    ttk.Button(tb, text="Reset edits", command=lambda: reset_edits()).pack(side="left")
    tb2 = ttk.Frame(pr); tb2.pack(fill="x", pady=(3, 0))
    ttk.Label(tb2, text="Side margin at steep edges (°)").pack(side="left")
    ttk.Spinbox(tb2, textvariable=pv["side"], from_=0, to=10, increment=0.5, width=5).pack(side="left", padx=(2, 8))
    ttk.Label(tb2, text="steep = more than").pack(side="left")
    ttk.Spinbox(tb2, textvariable=pv["steep"], from_=0.5, to=20, increment=0.5, width=5).pack(side="left", padx=2)
    ttk.Label(tb2, text="° altitude per ° azimuth   (0 = off; 1–2° recommended)", foreground="#666").pack(side="left")
    show_fill = tk.BooleanVar(value=False)
    ttk.Checkbutton(tb2, text="Show fill-in frames", variable=show_fill,
                    command=lambda: generate_fill(quiet=True) if show_fill.get() and fill["plan"] is None
                    else draw_fill()).pack(side="left", padx=(12, 2))
    ttk.Label(pr, text="Left-drag on the image to redraw the selected line · right-click: restore detected value · Ctrl+Z: undo",
              foreground="#666").pack(anchor="w")
    cf = ttk.Frame(pr); cf.pack(fill="both", expand=True)
    cv = tk.Canvas(cf, bg="#1c1c24", highlightthickness=0)
    hs = ttk.Scrollbar(cf, orient="horizontal", command=cv.xview)
    vs = ttk.Scrollbar(cf, orient="vertical", command=cv.yview)
    cv.configure(xscrollcommand=hs.set, yscrollcommand=vs.set)
    cv.grid(row=0, column=0, sticky="nsew"); vs.grid(row=0, column=1, sticky="ns"); hs.grid(row=1, column=0, sticky="ew")
    cf.rowconfigure(0, weight=1); cf.columnconfigure(0, weight=1)
    stats_lbl = ttk.Label(pr, text="", justify="left"); stats_lbl.pack(anchor="w", pady=3)
    q = queue.Queue(); result = {}
    E = {"solid": None, "branches": None, "undo": [], "bg": {}, "items": {}}
    MARG = 40

    def log(m): q.put(("log", m))
    def progress(f, m): q.put(("prog", (f, m)))

    def worker(opts):
        try:
            paths = find_fits(opts["folders"])
            log(f"Found {len(paths)} FITS files")
            frames, lat_f = load_frames(paths, log, progress, opts["bias"],
                                        default_cache_dir() if opts["cache"] else None)
            if not frames:
                raise RuntimeError("No frames with CENTALT/CENTAZ found.")
            lat = opts["lat"] if opts["lat"] is not None else lat_f
            if lat is None:
                raise RuntimeError("Latitude not in FITS headers; enter it in Settings.")
            K = opts["K"]
            if opts["auto"]:
                progress(0.37, "Calibrating rotation")
                K = calibrate_rotation(frames, lat, K, log, opts["full"])
            mos = build_mosaic(frames, lat, K, opts["res"], log, progress)
            progress(0.88, "Detecting horizon")
            det = detect_horizon(mos, opts["solid"], opts["thin"], log=log)
            q.put(("done", dict(mos=mos, det=det, K=K, lat=lat, n=len(frames))))
        except Exception as e:
            q.put(("error", f"{e}\n{traceback.format_exc()}"))

    def read_opts():
        return dict(folders=list(folders.get(0, "end")), K=float(pv["K"].get()), auto=pv["auto"].get(),
                    full=pv["full"].get(), lat=float(pv["lat"].get()) if pv["lat"].get().strip() else None,
                    res=float(pv["res"].get()), solid=float(pv["solid"].get()), thin=float(pv["thin"].get()),
                    bias=float(pv["bias"].get()) if pv["bias"].get().strip() else None, cache=pv["cache"].get())

    def start():
        if not folders.size():
            messagebox.showerror(APP_NAME, "Add at least one folder with FITS frames."); return
        try:
            opts = read_opts()
        except ValueError as e:
            messagebox.showerror(APP_NAME, f"Check settings: {e}"); return
        logbox.delete("1.0", "end"); result.clear(); E["bg"].clear()
        threading.Thread(target=worker, args=(opts,), daemon=True).start()

    def redetect():
        if "mos" not in result:
            messagebox.showerror(APP_NAME, "Build the mosaic first."); return
        try:
            o = read_opts()
        except ValueError as e:
            messagebox.showerror(APP_NAME, f"Check settings: {e}"); return
        result["det"] = detect_horizon(result["mos"], o["solid"], o["thin"], log=log)
        log(f"Re-detected with solid {o['solid']:.2f}, branches {o['thin']:.2f}")
        reset_edits(); report_uncertain()

    def old_pts():
        if old_var.get().strip():
            try:
                return load_hrz(old_var.get())
            except Exception as e:
                log(f"Current horizon file: {e}")
        return None

    def base_det():
        return fill_unmeasured(result["det"], old_pts())

    def reset_edits():
        if "det" not in result:
            return
        d = base_det()
        E["solid"] = np.array(d["solid"], float); E["branches"] = np.array(d["branches"], float)
        E["undo"].clear(); redraw_all(); refresh_fill()

    def undo():
        if E["undo"]:
            E["solid"], E["branches"] = E["undo"].pop(); redraw_lines(); update_stats(); refresh_fill()

    def side_settings():
        try:
            return max(float(pv["side"].get()), 0.0), max(float(pv["steep"].get()), 0.1)
        except ValueError:
            return 0.0, 3.0

    def final_lines():
        az = result["det"]["az"]; side, steep = side_settings()
        # The two lines are independent: an edited branches line may go below the solid line.
        s_ = widen_steep_edges(az, E["solid"], side, steep)
        b_ = widen_steep_edges(az, E["branches"], side, steep)
        return az, s_, b_

    def zoom():
        try:
            return float(ed["zoom"].get())
        except ValueError:
            return 6.0

    def geom():
        mos = result["mos"]; z = zoom()
        top = mos["alt0"] + mos["M"].shape[0] * mos["res"]
        return z, mos["alt0"], top

    def to_xy(az, alt):
        z, a0, top = geom()
        return MARG + az * z, MARG + (top - alt) * z

    def from_xy(x, y):
        z, a0, top = geom()
        return (x - MARG) / z, top - (y - MARG) / z

    def background():
        z = zoom()
        if z in E["bg"]:
            return E["bg"][z]
        mos = result["mos"]; M = mos["M"]
        v = M[np.isfinite(M)]
        lo, hi = np.percentile(v, [1, 99.5]) if v.size else (0, 1)
        g = (np.nan_to_num(np.clip((M - lo) / max(hi - lo, 1e-6), 0, 1), nan=0) * 255).astype(np.uint8)
        rgb = np.dstack([g, g, g]); rgb[np.isnan(M)] = (28, 28, 36)
        W = int(360 * z); H = int(M.shape[0] * mos["res"] * z)
        im = Image.fromarray(rgb[::-1]).resize((max(W, 1), max(H, 1)), Image.BILINEAR)
        E["bg"][z] = ImageTk.PhotoImage(im) if ImageTk else None
        return E["bg"][z]

    def redraw_all(*_):
        cv.delete("all")
        if "mos" not in result or E["solid"] is None:
            return
        z, a0, top = geom()
        bg = background()
        if bg is not None:
            cv.create_image(MARG, MARG, anchor="nw", image=bg)
        for a in range(0, 361, 10):
            x, _ = to_xy(a, top)
            cv.create_line(x, MARG, x, MARG + (top - a0) * z, fill="#3f6fb0" if a % 30 else "#6fa0e0")
            if a % 30 == 0 or z >= 12:
                cv.create_text(x, MARG - 12, text=f"{a}°", fill="#ddd")
        for hgt in range(int(math.ceil(a0 / 5) * 5), int(top) + 1, 5):
            _, y = to_xy(0, hgt)
            cv.create_line(MARG, y, MARG + 360 * z, y, fill="#3f6fb0")
            cv.create_text(MARG - 18, y, text=f"{hgt}°", fill="#ddd")
        cv.configure(scrollregion=(0, 0, MARG * 2 + 360 * z, MARG * 2 + (top - a0) * z))
        E["items"] = {}
        redraw_lines(); update_stats()

    def coords(az, alt):
        z, a0, top = geom()
        order = np.argsort(az); az = np.asarray(az)[order]; alt = np.clip(np.asarray(alt)[order], a0, top)
        az = np.concatenate([az, [az[0] + 360]]); alt = np.concatenate([alt, [alt[0]]])
        pts = []
        for a, h in zip(az, alt):
            pts += list(to_xy(a, h))
        return pts

    def redraw_lines():
        if E["solid"] is None:
            return
        az, s_, b_ = final_lines()
        specs = [("solid", s_, "#ff3c28", 2), ("branches", b_, "#00dcff", 2)]
        op = old_pts()
        if op is not None:
            specs.insert(0, ("old", hrz_array(op, az), "#f0c800", 2))
        elif "old" in E["items"]:
            cv.delete(E["items"].pop("old"))
        for key, vals, col, w in specs:
            c = coords(az, vals)
            if key in E["items"]:
                cv.coords(E["items"][key], *c)
            else:
                E["items"][key] = cv.create_line(*c, fill=col, width=w)
        draw_fill()

    def paint(ev, push):
        tgt = ed["target"].get()
        if tgt == "none" or E["solid"] is None:
            return
        if push:
            E["undo"].append((E["solid"].copy(), E["branches"].copy()))
            E["undo"][:] = E["undo"][-50:]
        x, y = cv.canvasx(ev.x), cv.canvasy(ev.y)
        a, h = from_xy(x, y)
        _, a0, top = geom()
        h = min(max(h, max(a0, 0.0)), min(top, 90.0))   # keep edits inside the shown altitude range
        try:
            bw = max(float(ed["brush"].get()), 0.1)
        except ValueError:
            bw = 1.0
        az = result["det"]["az"]
        m = np.abs((az - a + 180) % 360 - 180) <= bw / 2
        E[tgt][m] = h
        redraw_lines()

    def restore(ev):
        tgt = ed["target"].get()
        if tgt == "none" or E["solid"] is None:
            return
        E["undo"].append((E["solid"].copy(), E["branches"].copy()))
        a, _ = from_xy(cv.canvasx(ev.x), cv.canvasy(ev.y))
        try:
            bw = max(float(ed["brush"].get()), 0.1)
        except ValueError:
            bw = 1.0
        d = base_det(); az = result["det"]["az"]
        m = np.abs((az - a + 180) % 360 - 180) <= bw / 2
        E[tgt][m] = np.asarray(d[tgt])[m]
        redraw_lines(); update_stats(); refresh_fill()

    current_hz.trace_add("write", lambda *a: redraw_all())
    for k in ("side", "steep"):
        pv[k].trace_add("write", lambda *a: (redraw_lines(), update_stats(), refresh_fill()) if E.get("solid") is not None else None)
    cv.bind("<ButtonPress-1>", lambda e: paint(e, True))
    cv.bind("<B1-Motion>", lambda e: paint(e, False))
    cv.bind("<ButtonRelease-1>", lambda e: (update_stats(), refresh_fill()))
    cv.bind("<Button-3>", restore)
    cv.bind("<MouseWheel>", lambda e: cv.xview_scroll(int(-e.delta / 40), "units"))
    root.bind("<Control-z>", lambda e: undo())
    zc.bind("<<ComboboxSelected>>", redraw_all)

    def update_stats():
        if E["solid"] is None:
            return
        az, s_, b_ = final_lines()
        a_s = sky_area_sqdeg(az, s_); a_b = sky_area_sqdeg(az, b_)
        txt = (f"{result['n']} frames · rotation offset {result['K']:.1f}° · latitude {result['lat']:.3f}°"
               f"{'  · edited' if E['undo'] else ''}"
               f"{'  · side margin %.1f° at steep edges' % side_settings()[0] if side_settings()[0] > 0 else ''}\n"
               f"Usable sky above horizon:  solid {a_s:,.0f} sq.deg ({a_s / 206.26:.1f}% of hemisphere)   "
               f"with branches {a_b:,.0f} sq.deg ({a_b / 206.26:.1f}%)")
        op = old_pts()
        if op is not None:
            azf = np.arange(0, 360, 0.1); a_o = sky_area_sqdeg(azf, hrz_array(op, azf))
            txt += (f"\nCurrent horizon {a_o:,.0f} sq.deg → new with branches {a_b - a_o:+,.0f} "
                    f"({(a_b / a_o - 1) * 100:+.1f}%), solid {a_s - a_o:+,.0f} ({(a_s / a_o - 1) * 100:+.1f}%)")
        if not result["det"]["measured"].all():
            frac = result["det"]["measured"].mean() * 100
            txt += (f"\nOnly {frac:.0f}% of azimuths were imaged; the rest "
                    + ("uses the current horizon." if op is not None else "is interpolated – load the current horizon to fill it."))
        stats_lbl.config(text=txt)

    def report_uncertain():
        unc = result["det"]["uncertain"]
        fill.update(plan=None, boxes=[], p=None); fill_info.config(text="")
        if show_fill.get():
            generate_fill(quiet=True)       # new horizon: show its fill-in right away
        else:
            draw_fill()
        d = result["det"]
        boxes = fillin_boxes(d, lines=(d["az"], d["solid"], d["branches"])) if unc else []
        if boxes:
            log(f"{len(boxes)} areas where the horizon line reaches the edge of the imaged area "
                f"– use 'Generate fill-in frames' to plan them.")
        else:
            log("No uncertain areas – every azimuth has sky and obstruction inside the imaged area.")

    def save_outputs():
        if E["solid"] is None:
            messagebox.showerror(APP_NAME, "Build the horizon first."); return
        which = save_which.get()
        az, s_, b_ = final_lines()
        vals = b_ if which == "with branches" else s_
        name = time.strftime("horizon_%Y-%m-%d") + ("_with_branches" if which == "with branches" else "") + ".hrz"
        f = filedialog.asksaveasfilename(title="Save horizon file", defaultextension=".hrz", initialfile=name,
                                         filetypes=[("NINA horizon", "*.hrz"), ("All", "*.*")])
        if not f:
            return
        write_hrz(f, az, vals)
        messagebox.showinfo(APP_NAME, f"Saved the {which} horizon:\n{f}\n\n"
                                      "In NINA: Options → General → Astrometry → Horizon.")

    def compute_fill():
        """(params, boxes, plan) for the fill-in frames; raises ValueError with a readable message."""
        if "det" not in result or E.get("solid") is None:
            raise ValueError("Build the horizon first.")
        try:
            p = collect(require_area=False)
        except Exception as e:
            raise ValueError(f"Check the settings on the Plan tab: {e}")
        try:
            ef = max(float(fill_edge.get()), 0.0); es = max(float(fill_slope.get()), 0.1)
        except ValueError:
            raise ValueError("Check the fill-in steep-edge values.")
        az_, s_, b_ = final_lines()
        boxes = fillin_boxes(result["det"], floor=p.floor, ceiling=p.ceiling, lines=(az_, s_, b_))
        boxes = boxes + steep_edge_boxes(az_, (s_, b_), max(p.fov), ef, es, floor=p.floor, ceiling=p.ceiling)
        if not boxes:
            return p, [], None
        plan = plan_fillin(p, boxes, mos=result.get("mos"))
        used = set(plan["labels"])
        boxes = [b for b in boxes if b[4] in used]        # drop steep-edge areas that are fully imaged
        return p, boxes, (plan if boxes else None)

    fill = {"p": None, "boxes": [], "plan": None}       # the last generated fill-in plan

    def generate_fill(quiet=False):
        """Plan the fill-in frames with the current lines and settings, store and show them."""
        try:
            p, boxes, plan = compute_fill()
        except ValueError as e:
            if not quiet:
                messagebox.showerror(APP_NAME, str(e))
            log(f"Fill-in not generated: {e}"); return False
        except Exception as e:                          # never fail silently
            log(f"Fill-in failed: {type(e).__name__}: {e}")
            if not quiet:
                messagebox.showerror(APP_NAME, f"Fill-in failed: {e}")
            return False
        fill.update(p=p, boxes=boxes, plan=plan)
        if plan:
            n_edge = sum(1 for b in boxes if b[4].startswith("Steep edge"))
            txt = (f"Fill-in: {len(boxes)} areas ({n_edge} steep edges), {len(plan['frames'])} frames, "
                   f"~{plan['seconds'] / 60:.0f} min")
        else:
            txt = "Fill-in: nothing to shoot – every edge is imaged."
        fill_info.config(text=txt)
        if not quiet:
            log(txt + ("\n  " + "\n  ".join(b[4] for b in boxes) if boxes else ""))
            show_fill.set(True)
        draw_fill()
        return True

    def refresh_fill():
        """Re-generate after a change of lines or settings, if a fill-in plan is being shown."""
        if fill["plan"] is not None or (show_fill.get() and fill["p"] is not None):
            generate_fill(quiet=True)

    def draw_fill():
        cv.delete("fill")
        plan, boxes = fill["plan"], fill["boxes"]
        if not show_fill.get() or "mos" not in result or E.get("solid") is None or not plan:
            return
        X = lambda a: to_xy(a, 0)[0]; Y = lambda h: to_xy(0, h)[1]
        for a0, a1, h0, h1, label in boxes:
            for sh in (0, 360, -360):
                if a1 + sh < 0 or a0 + sh > 360:
                    continue
                cv.create_rectangle(X(a0 + sh), Y(h1), X(a1 + sh), Y(h0), outline="#ff60ff", dash=(5, 3),
                                    width=2, tags="fill")
            cv.create_text(X(max(a0, 0)) + 3, Y(h1) - 8, anchor="w", text=label, fill="#ff9cff",
                           font=("Segoe UI", 8), tags="fill")
        for alt, a, rot in plan["frames"]:
            cw = plan["fw"] / max(math.cos(math.radians(alt)), 0.05) / 2; ch = plan["fh"] / 2
            for sh in (0, 360, -360):
                if a + cw + sh < 0 or a - cw + sh > 360:
                    continue
                cv.create_rectangle(X(a - cw + sh), Y(alt + ch), X(a + cw + sh), Y(alt - ch),
                                    outline="#ffa000", width=1, tags="fill")
        draw_order(cv, plan, lambda a: X(a % 360), Y, tag="fill",
                   numbers="all" if zoom() >= 20 else "rows", flip=False, width=1)

    def save_fill():
        if fill["plan"] is None and not generate_fill():
            return
        if fill["plan"] is None:
            messagebox.showinfo(APP_NAME, "Nothing to fill in – every edge is imaged."); return
        p, plan = fill["p"], fill["plan"]
        f = filedialog.asksaveasfilename(defaultextension=".json", initialfile="Horizon fill-in.json",
                                         filetypes=[("NINA sequence", "*.json")])
        if f:
            out = save_plan(p, plan, f)
            messagebox.showinfo(APP_NAME, f"Saved {len(plan['frames'])} frames:\n{out}\n\n"
                                          "Camera, rotator and exposure settings come from the Plan tab.")

    ttk.Button(bb, text="Build mosaic & horizon", command=start).grid(row=0, column=0, padx=2, pady=2, sticky="ew")
    ttk.Button(bb, text="Re-detect (new thresholds)", command=redetect).grid(row=0, column=1, padx=2, pady=2, sticky="ew")
    fb = ttk.Frame(bb); fb.grid(row=3, column=0, columnspan=2, padx=2, pady=2, sticky="ew")
    ttk.Button(fb, text="Generate fill-in frames", command=generate_fill).pack(side="left")
    ttk.Button(fb, text="Save fill-in sequence…", command=save_fill).pack(side="left", padx=4)
    fill_info = ttk.Label(fb, text="", foreground="#555"); fill_info.pack(side="left", padx=4)
    fe = ttk.Frame(bb); fe.grid(row=2, column=0, columnspan=2, padx=2, pady=2, sticky="w")
    # Start from the Plan tab's values and follow them whenever they change there.
    fill_edge = tk.StringVar(value=V["edge_frames"].get()); fill_slope = tk.StringVar(value=V["edge_slope"].get())
    V["edge_frames"].trace_add("write", lambda *a: fill_edge.set(V["edge_frames"].get()))
    V["edge_slope"].trace_add("write", lambda *a: fill_slope.set(V["edge_slope"].get()))
    ttk.Label(fe, text="Fill-in: extra frames at steep edges").pack(side="left")
    ttk.Spinbox(fe, textvariable=fill_edge, from_=0, to=5, increment=0.5, width=4).pack(side="left", padx=(2, 8))
    ttk.Label(fe, text="steep = more than").pack(side="left")
    ttk.Spinbox(fe, textvariable=fill_slope, from_=0.5, to=20, increment=0.5, width=4).pack(side="left", padx=2)
    ttk.Label(fe, text="°/°  (0 = off; copied from the Plan tab)", foreground="#666").pack(side="left")
    for _v in (fill_edge, fill_slope):
        _v.trace_add("write", lambda *a: refresh_fill())
    sv = ttk.Frame(bb); sv.grid(row=1, column=0, columnspan=2, padx=2, pady=2, sticky="w")
    save_which = tk.StringVar(value="solid")
    ttk.Combobox(sv, textvariable=save_which, values=["solid", "with branches"], width=12,
                 state="readonly").pack(side="left")
    ttk.Button(sv, text="Save horizon file…", command=save_outputs).pack(side="left", padx=2)

    def poll():
        try:
            while True:
                kind, val = q.get_nowait()
                if kind == "log":
                    logbox.insert("end", val + "\n"); logbox.see("end")
                elif kind == "prog":
                    pb["value"] = val[0]; status.config(text=val[1])
                elif kind == "done":
                    result.update(val); pb["value"] = 1.0; status.config(text="Done")
                    reset_edits(); report_uncertain()
                elif kind == "error":
                    status.config(text="Error"); logbox.insert("end", val + "\n")
                    messagebox.showerror(APP_NAME, val.splitlines()[0])
        except queue.Empty:
            pass
        root.after(150, poll)
    poll()

    # ------------------------------------------------------------------ HELP
    help_tab = ttk.Frame(nb); nb.add(help_tab, text="Help")
    t = tk.Text(help_tab, wrap="word"); t.pack(fill="both", expand=True)
    t.insert("1.0", HELP_TEXT); t.config(state="disabled")
    root.mainloop()


HELP_TEXT = """HORIZON MAPPER FOR NINA   (the version is shown in the window title)

WORKFLOW
1. Plan sequence
   • Enter camera/telescope (or 'Read from a FITS file…'). The field of view is computed.
   • Choose the area: a fixed altitude band, or a strip around your current horizon file
     (optionally plus a second estimate, e.g. from an all-sky camera); the strip spans both,
     plus margins. The current horizon is also filled in on the Build tab for comparison.
   • 'Extra frames at steep edges' (horizon mode): where the current horizon rises or drops
     more steeply than 'Steep = more than' (altitude ° per azimuth °, 3 ≈ 72°), the strip is
     widened sideways by this many frame widths and covers the whole height of the jump, so
     the sides of trees and roofs are imaged. 1 is a good default; 0 = off.
     Fractions work: 0.5 usually gives two frames per row at a narrow edge.
   • Each row gets frames only where the strip enters it, as few as the overlap allows.
   • 'Lowest altitude' keeps frames (and the sun) out of the low sky.
   • 'Preview plan' draws the frames; with 'Show shooting order' (on by default) a colour-graded
     path shows the order (green = first → yellow → magenta = last), with the first frame number
     of each row, Start/End markers and the meridian flip.
   • Rotator: keeps every frame level with the horizon, so overlap can be ~10 %.
     Offset = sky position angle − mechanical angle (from NINA's rotator panel after a
     plate solve). A wrong 180° doesn't matter here.
   • East-half-first ordering means the mount flips over the meridian only once.
   • Preview, then 'Save NINA sequence…'. Open it with Load sequence in the Advanced
     Sequencer (or save as template and put it in NINA's template folder).
   • Daytime under even cloud gives the best silhouettes. Take a test frame first:
     open sky should read ~10,000–30,000 ADU. Use 'Sun check' and keep the lowest
     altitude well above the sun.

2. Build horizon
   • Add the folder(s) with the frames. NINA's headers provide pointing (CENTALT/CENTAZ),
     rotator angle, pier side, pixel size, focal length and latitude.
   • 'Calibrate rotation from overlaps' measures the true camera angle from how well
     neighbouring frames line up (use 'Full 360° search' if you have no idea of the offset).
   • The horizon is the highest altitude where the view is blocked:
       solid      – pixels darker than 0.60 × the local sky (trees, roofs, walls)
       branches   – darker than 0.85 × sky (also thin bare branches)
   • Fill-in areas are where the solid or branches line (with your edits) reaches the edge of
     the imaged area – 'clear to bottom' or 'blocked to top' – so the true horizon may lie
     outside it. Dark frame edges high in the sky that the lines rejected are ignored.
   • 'Generate fill-in frames' plans the frames and shows them on the image, with a summary in
     the log. 'Save fill-in sequence…' saves exactly the plan that is shown (camera, rotator and
     exposure settings come from the Plan tab). Shoot it, add its folder and build again.
   • 'Fill-in: extra frames at steep edges' (above the buttons): around every steep rise or
     drop of the new horizon line, an area ('Steep edge …') spanning the jump and this many
     frame widths to each side is checked; only frames over gaps in the mosaic are shot.
     0 = off. The values are copied from the Plan tab and follow changes made there; you can
     still change them here. Frames shared by overlapping areas are shot only once.
   • 'Show fill-in frames' (toolbar above the image) draws the fill-in areas (dashed magenta),
     their frames (orange) and the shooting order (green = first → magenta = last, numbers =
     first frame of each row; every frame is numbered at zoom 20). It follows your edits and
     the fill-in settings: after a change the plan is generated again and redrawn.
   • 'Re-detect' applies new thresholds instantly without re-reading the frames.
   • 'Side margin at steep edges' (toolbar above the image): where the horizon rises or drops steeply (the side of a
     tree, the end of a roof), the higher value is extended that many degrees sideways, so a
     target just beside the obstruction isn't scheduled too early. Gentle slopes are not
     changed. 'Steep' means more than the given degrees of altitude per degree of azimuth
     (3 ≈ a 72° slope). 1-2° is a sensible margin; 0 switches it off. The lines on screen and
     the saved files include the margin; your edits are made on the unwidened line.
   • The two lines are edited independently – the branches line may be drawn below the solid
     line (e.g. where the "solid" part is really a thin, see-through tree top).
   • Editing: choose 'solid line' or 'branches line' above the image, then left-drag to
     redraw the line (brush width in degrees), right-click to restore the detected value,
     Ctrl+Z to undo. Zoom and scroll to work precisely. Edits are used when saving.
   • Reduced frames are cached (in %LOCALAPPDATA%\\HorizonMapper\\cache), so building again
     – e.g. after adding fill-in frames – only reads the new files. 'Clear cache' frees the space.
   • 'Save horizon file…' saves one NINA .hrz file: choose 'solid' or 'with branches' next to
     the button. Edits and the side margin are included.
     In NINA: Options → General → Astrometry → Horizon file.

NOTES
   • Very close objects (observatory walls) appear as smooth blurred shapes; they still
     count as obstruction for the telescope.
   • Frames without CENTALT/CENTAZ (e.g. no mount connected) are skipped.
"""


# =============================================================================
#  Command line
# =============================================================================

def main():
    # A windowed .exe has no console: give print() somewhere harmless to write.
    if sys.stdout is None:
        sys.stdout = open(os.devnull, "w")
    if sys.stderr is None:
        sys.stderr = open(os.devnull, "w")
    ap = argparse.ArgumentParser(description=APP_NAME)
    sub = ap.add_subparsers(dest="cmd")
    pp = sub.add_parser("plan", help="generate a NINA sequence")
    pp.add_argument("--fov", nargs=2, type=float, metavar=("W", "H"))
    pp.add_argument("--sensor", nargs=3, type=float, metavar=("WPX", "HPX", "UM"), default=[6248, 4176, 3.76])
    pp.add_argument("--focal", type=float, default=408)
    pp.add_argument("--lat", type=float, default=61.57)
    pp.add_argument("--band", nargs=2, type=float, metavar=("MIN", "MAX"))
    pp.add_argument("--hrz", nargs="+", help="1-2 horizon files")
    pp.add_argument("--below", type=float, default=4); pp.add_argument("--above", type=float, default=4)
    pp.add_argument("--floor", type=float, default=20)
    pp.add_argument("--edge-frames", type=float, default=1.0,
                    help="extra frame widths beside steep edges of the horizon file(s) (0 = off)")
    pp.add_argument("--edge-slope", type=float, default=3.0,
                    help="an edge is steep above this many degrees altitude per degree azimuth")
    pp.add_argument("--az", nargs=2, type=float, default=[0, 360])
    pp.add_argument("--overlap", type=float, default=0.10)
    pp.add_argument("--no-rotator", action="store_true")
    pp.add_argument("--rot-offset", type=float, default=-58.04)
    pp.add_argument("--rot-tol", type=float, default=3)
    pp.add_argument("--exposure", type=float, default=0.008)
    pp.add_argument("--gain", type=int, default=0)
    pp.add_argument("--bin", type=int, default=2)
    pp.add_argument("--template", action="store_true")
    pp.add_argument("-o", "--out", default="Horizon mapping.json")
    bp = sub.add_parser("build", help="build mosaic + horizon from FITS folders")
    bp.add_argument("folders", nargs="+")
    bp.add_argument("--rot-offset", type=float, default=-58.04)
    bp.add_argument("--no-calibrate", action="store_true")
    bp.add_argument("--full-search", action="store_true")
    bp.add_argument("--lat", type=float)
    bp.add_argument("--res", type=float, default=0.03)
    bp.add_argument("--old", help="old .hrz for comparison")
    bp.add_argument("--side-margin", type=float, default=0.0,
                    help="extend the higher side of steep horizon edges this many degrees sideways")
    bp.add_argument("--steep-slope", type=float, default=3.0,
                    help="an edge is steep above this many degrees altitude per degree azimuth")
    bp.add_argument("--no-cache", action="store_true", help="don't use/create the frame cache")
    bp.add_argument("--fillin", metavar="FILE.json", help="also write a NINA fill-in sequence for uncertain areas")
    bp.add_argument("--fillin-exposure", type=float, default=0.008)
    bp.add_argument("--fillin-floor", type=float, default=15.0)
    bp.add_argument("--fillin-edge-frames", type=float, default=1.0,
                    help="for --fillin: extra frame widths beside steep horizon edges (0 = off)")
    bp.add_argument("--fillin-edge-slope", type=float, default=3.0,
                    help="for --fillin: steep = more than this many deg altitude per deg azimuth")
    bp.add_argument("--sensor", nargs=3, type=float, metavar=("WPX", "HPX", "UM"), default=[6248, 4176, 3.76],
                    help="for --fillin: unbinned sensor size and pixel size")
    bp.add_argument("--focal", type=float, default=408, help="for --fillin: focal length (mm)")
    bp.add_argument("-o", "--outbase", default="horizon")
    a = ap.parse_args()

    if a.cmd is None:
        run_gui(); return
    if a.cmd == "plan":
        p = PlanParams(sensor_w_px=int(a.sensor[0]), sensor_h_px=int(a.sensor[1]), pixel_um=a.sensor[2],
                       focal_mm=a.focal, lat=a.lat, floor=a.floor, az_start=a.az[0], az_end=a.az[1],
                       overlap=a.overlap, rotator=not a.no_rotator, rot_offset=a.rot_offset,
                       rot_tolerance=a.rot_tol, exposure=a.exposure, gain=a.gain, binning=a.bin,
                       template=a.template, margin_below=a.below, margin_above=a.above,
                       edge_frames=a.edge_frames, edge_slope=a.edge_slope)
        if a.hrz:
            p.mode = "horizon"; p.horizon_files = a.hrz
        else:
            p.mode = "band"; p.alt_min, p.alt_max = a.band or (20, 45)
        plan = make_plan(p)
        out = save_plan(p, plan, a.out)
        print(f"{out}: {len(plan['frames'])} frames, {len(plan['grid'])} rows, {plan['orientation']}, "
              f"{plan['rot_moves']} rotator moves, ~{plan['seconds'] / 3600:.1f} h")
    elif a.cmd == "build":
        t0 = time.time()
        log = lambda m: print(m, flush=True)
        prog = lambda f, m: None
        frames, lat_f = load_frames(find_fits(a.folders), log, prog,
                                    cache_dir=None if a.no_cache else default_cache_dir())
        lat = a.lat if a.lat is not None else lat_f
        K = a.rot_offset if a.no_calibrate else calibrate_rotation(frames, lat, a.rot_offset, log, a.full_search)
        mos = build_mosaic(frames, lat, K, a.res, log, prog)
        det = detect_horizon(mos, log=log)
        if a.old:
            det = fill_unmeasured(det, load_hrz(a.old))
        if a.side_margin > 0:
            det = dict(det)
            det["solid"] = widen_steep_edges(det["az"], det["solid"], a.side_margin, a.steep_slope)
            det["branches"] = np.maximum(widen_steep_edges(det["az"], det["branches"], a.side_margin, a.steep_slope),
                                         det["solid"])
        write_hrz(a.outbase + ".hrz", det["az"], det["solid"])
        write_hrz(a.outbase + "_with_branches.hrz", det["az"], det["branches"])
        lines = [("Solid obstructions", det["az"], det["solid"], (230, 40, 30)),
                 ("Incl. thin branches", det["az"], det["branches"], (0, 200, 240))]
        if a.old:
            az = np.arange(0, 360, 0.25); lines.append(("Old horizon", az, hrz_array(load_hrz(a.old), az), (240, 200, 0)))
        render_mosaic(mos, lines).save(a.outbase + "_mosaic.png")
        print(f"usable sky: solid {sky_area_sqdeg(det['az'], det['solid']):.0f} sq.deg, "
              f"branches {sky_area_sqdeg(det['az'], det['branches']):.0f} sq.deg")
        boxes = fillin_boxes(det, floor=a.fillin_floor,
                             lines=(det["az"], det["solid"], det["branches"])) if det["uncertain"] else []
        if boxes:
            print(f"{len(det['uncertain'])} uncertain azimuth samples in {len(boxes)} areas:",
                  "; ".join(b[4] for b in boxes))
        if a.fillin:
            p = PlanParams(sensor_w_px=int(a.sensor[0]), sensor_h_px=int(a.sensor[1]), pixel_um=a.sensor[2],
                           focal_mm=a.focal, lat=lat, rot_offset=a.rot_offset, exposure=a.fillin_exposure,
                           floor=a.fillin_floor)
            boxes = boxes + steep_edge_boxes(det["az"], (det["solid"], det["branches"]), max(p.fov),
                                             a.fillin_edge_frames, a.fillin_edge_slope,
                                             floor=p.floor, ceiling=p.ceiling)
            plan = plan_fillin(p, boxes, mos=mos) if boxes else None
            if plan and plan["frames"]:
                n_edge = len({l for l in plan["labels"] if l.startswith("Steep edge")})
                print(f"fill-in: {save_plan(p, plan, a.fillin)} ({len(plan['frames'])} frames, "
                      f"{n_edge} steep-edge areas with gaps, ~{plan['seconds'] / 60:.0f} min)")
            else:
                print("fill-in: nothing to shoot")
        print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
