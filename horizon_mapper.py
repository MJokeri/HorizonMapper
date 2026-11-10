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
APP_VERSION = "1.1.0"
# DEVELOPMENT DEFAULTS: debug log and log windows on. Set DEV_MODE = False before a release
# (README → "Release checklist"). Users can still switch both in the Debug menu.
DEV_MODE = True

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


def altaz_to_hadec(alt, az, lat):
    """Alt/az → hour angle and declination, degrees (HA positive to the west)."""
    A, h, p = math.radians(az), math.radians(alt), math.radians(lat)
    dec = math.asin(max(-1.0, min(1.0, math.sin(p) * math.sin(h) + math.cos(p) * math.cos(h) * math.cos(A))))
    ha = math.atan2(-math.sin(A) * math.cos(h), math.cos(p) * math.sin(h) - math.sin(p) * math.cos(h) * math.cos(A))
    return math.degrees(ha), math.degrees(dec)


def hadec_to_altaz(ha, dec, lat):
    """Hour angle / declination → alt/az, degrees."""
    H, d, p = math.radians(ha), math.radians(dec), math.radians(lat)
    alt = math.asin(max(-1.0, min(1.0, math.sin(p) * math.sin(d) + math.cos(p) * math.cos(d) * math.cos(H))))
    az = math.atan2(-math.sin(H) * math.cos(d), math.cos(p) * math.sin(d) - math.sin(p) * math.cos(d) * math.cos(H))
    return math.degrees(alt), math.degrees(az) % 360


def angular_sep(alt1, az1, alt2, az2):
    a1, a2 = math.radians(alt1), math.radians(alt2)
    c = math.sin(a1) * math.sin(a2) + math.cos(a1) * math.cos(a2) * math.cos(math.radians(az1 - az2))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def slew_path(p_from, p_to, lat, n=24):
    """Approximate path of an equatorial mount between two alt/az pointings: both axes move at once,
    so hour angle and declination change linearly. Returns alt/az points (ends included)."""
    h1, d1 = altaz_to_hadec(p_from[0], p_from[1], lat)
    h2, d2 = altaz_to_hadec(p_to[0], p_to[1], lat)
    dh = (h2 - h1 + 180) % 360 - 180
    return [hadec_to_altaz(h1 + dh * k / n, d1 + (d2 - d1) * k / n, lat) for k in range(n + 1)]


def path_sun_sep(points, sun):
    """Smallest distance (degrees) between a list of alt/az points and the Sun (alt, az)."""
    return min(angular_sep(a, z, sun[0], sun[1]) for a, z in points) if points else 180.0


def safe_waypoint(p_from, p_to, lat, sun, sun_min):
    """An intermediate pointing that keeps both legs of a slew at least `sun_min` from the Sun:
    first move only one axis, then the other; else via the meridian high in the sky. None if no
    such route is found."""
    h1, d1 = altaz_to_hadec(p_from[0], p_from[1], lat)
    h2, d2 = altaz_to_hadec(p_to[0], p_to[1], lat)
    cands = [hadec_to_altaz(h1, d2, lat), hadec_to_altaz(h2, d1, lat),
             hadec_to_altaz(0.0, (d1 + d2) / 2, lat), hadec_to_altaz(0.0, lat, lat)]
    for w in cands:
        if w[0] < 5:
            continue
        if min(path_sun_sep(slew_path(p_from, w, lat), sun), path_sun_sep(slew_path(w, p_to, lat), sun)) >= sun_min:
            return w
    return None


def parse_start(date_txt, time_txt, now=None):
    """Local start time from a date ('2026-11-12', '12.11.2026', '12.11.' or blank = today) and a time
    ('13:30', '13.30' or blank = the current time of day). Raises ValueError for anything else."""
    now = time.time() if now is None else now
    lt = time.localtime(now)
    d = (date_txt or "").strip()
    if not d:
        y, mo, da = lt.tm_year, lt.tm_mon, lt.tm_mday
    elif "-" in d:
        y, mo, da = (int(x) for x in d.split("-")[:3])
    else:
        parts = [x for x in d.split(".") if x != ""]
        if len(parts) == 2:
            da, mo, y = int(parts[0]), int(parts[1]), lt.tm_year
        elif len(parts) == 3:
            da, mo, y = int(parts[0]), int(parts[1]), int(parts[2])
            y += 2000 if y < 100 else 0
        else:
            raise ValueError("date")
    t = (time_txt or "").strip()
    if t:
        hh, mm = (int(x) for x in t.replace(".", ":").split(":")[:2])
        ss = 0
    else:
        hh, mm, ss = lt.tm_hour, lt.tm_min, lt.tm_sec
    if not (1 <= mo <= 12 and 1 <= da <= 31 and 0 <= hh <= 23 and 0 <= mm <= 59):
        raise ValueError("date or time out of range")
    return time.mktime((y, mo, da, hh, mm, ss, 0, 0, -1))


def frame_times(frames, order, t0, sec_per_frame, sec_per_rot=0.0):
    """Estimated time each frame in `order` is shot, starting at t0 (rotator moves add time)."""
    t, out, last_rot = t0, {}, None
    for i in order:
        r = frames[i][2]
        if r is not None and r != last_rot:
            t += sec_per_rot; last_rot = r
        out[i] = t
        t += sec_per_frame
    return out


def sun_risk(frames, order, t0, sec_per_frame, sec_per_rot, lat, lon, sun_min, start_pos=None):
    """Which frames of a run come closer to the Sun than `sun_min`, at the time they are shot, and
    which slews (from the previous pointing, or `start_pos` for the first) pass too close on the way.
    Returns dict(frames=[(i, t, sep)], paths=[(i, t, sep)], end=t_end, track=[(t, alt, az)])."""
    times = frame_times(frames, order, t0, sec_per_frame, sec_per_rot)
    fr, pa = [], []
    prev = start_pos
    for i in order:
        t = times[i]; sun = sun_altaz(t, lat, lon)
        alt, az = frames[i][0], frames[i][1]
        if sun[0] > -sun_min:                          # the Sun can matter at all
            sep = angular_sep(alt, az, sun[0], sun[1])
            if sep < sun_min:
                fr.append((i, t, sep))
            elif prev is not None:
                psep = path_sun_sep(slew_path(prev, (alt, az), lat), sun)
                if psep < sun_min:
                    pa.append((i, t, psep))
        prev = (alt, az)
    t_end = (max(times.values()) + sec_per_frame) if times else t0
    steps = max(int((t_end - t0) / 600), 1)
    track = [(t0 + (t_end - t0) * k / steps,) + sun_altaz(t0 + (t_end - t0) * k / steps, lat, lon)
             for k in range(steps + 1)]
    return dict(frames=fr, paths=pa, end=t_end, track=track)


def sun_circle(sun_alt, sun_az, radius, n=72):
    """Alt/az points at `radius` degrees around the Sun (for drawing the keep-away zone)."""
    a, r = math.radians(sun_alt), math.radians(radius)
    pts = []
    for k in range(n + 1):
        b = 2 * math.pi * k / n
        alt2 = math.asin(max(-1.0, min(1.0, math.sin(a) * math.cos(r) + math.cos(a) * math.sin(r) * math.cos(b))))
        daz = math.atan2(math.sin(b) * math.sin(r) * math.cos(a), math.cos(r) - math.sin(a) * math.sin(alt2))
        pts.append((math.degrees(alt2), (sun_az + math.degrees(daz)) % 360))
    return pts


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
        self.edge_frames = 0.5        # extra frame widths beside steep horizon edges (0 = off)
        self.edge_slope = 2.0         # "steep" = more than this many deg altitude per deg azimuth
        self.overlap = 0.10
        # rotator
        self.rotator = True; self.rot_offset = -58.04; self.rot_sign = 1; self.rot_tolerance = 3.0
        self.orientation = "auto"     # auto | landscape | portrait
        # exposure
        self.exposure = 0.008; self.gain = 0; self.offset = -1; self.binning = 2
        # sequence
        self.unpark = True; self.park = True; self.meridian_split = True
        self.bottom_up = False        # with the east/west split: west half bottom→top too
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
        west = [(a, [z for z in zs if z >= 180]) for a, zs in (grid if p.bottom_up else reversed(grid))]
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


def gap_boxes(mos, az, lines, near=1.0, min_gap=0.02, pad_az=0.5, floor=0.0, ceiling=85.0,
              names=("solid", "branches")):
    """Fill-in areas where a horizon line runs through, or within `near` degrees of, a spot of the
    mosaic without pictures (inside the mosaic or past its edge).

    For every azimuth sample and line, the box line ± `near` (altitude) × ± `near`/cos(alt)
    (azimuth) is checked; if more than `min_gap` of it has no data, the sample is marked.
    Neighbouring marked samples (within 1°) form one area around that stretch of the line.
    """
    if near <= 0 or mos is None:
        return []
    M = mos["M"]; res = mos["res"]; a0 = mos["alt0"]; NH, NA = M.shape
    az = np.asarray(az, float)
    boxes = []
    for v, name in zip(lines, names):
        v = np.asarray(v, float)
        bad = np.zeros(len(az), bool)
        for i, (a, h) in enumerate(zip(az, v)):
            hw = near / max(math.cos(math.radians(min(h, 80.0))), 0.17)
            rs = np.floor((np.linspace(h - near, h + near, 9) - a0) / res).astype(int)
            cs = np.floor((np.linspace(a - hw, a + hw, 9) % 360) / 360 * NA).astype(int) % NA
            inside = (rs >= 0) & (rs < NH)
            vals = M[np.ix_(rs[inside], cs)]
            missing = (~inside).sum() * len(cs) + (~np.isfinite(vals)).sum()
            bad[i] = missing / (len(rs) * len(cs)) > min_gap
        idx = np.where(bad)[0]
        if not len(idx):
            continue
        step = float(np.median(np.diff(az))) if len(az) > 1 else 0.5
        groups = [[idx[0], idx[0]]]
        for k in idx[1:]:
            if (az[k] - az[groups[-1][1]]) <= 1.0 + 1e-9:
                groups[-1][1] = k
            else:
                groups.append([k, k])
        if len(groups) > 1 and az[groups[0][0]] - step < 1e-9 and az[groups[-1][1]] + step >= 360 - 1e-9:
            last = groups.pop(); groups[0][0] = last[0] - len(az)       # join across 0/360
        for g0, g1 in groups:
            ii = np.arange(g0, g1 + 1) % len(az)
            seg = v[ii]
            h0 = max(floor, float(seg.min()) - near); h1 = min(ceiling, float(seg.max()) + near)
            za = az[g0 % len(az)] - (360 if g0 < 0 else 0)
            zb = az[g1] + (360 if az[g1] < za else 0)
            za -= pad_az + near; zb += step + pad_az + near
            if h1 > h0:
                boxes.append((za, zb, h0, h1, f"Gap Az {za % 360:.0f}-{zb % 360:.0f} ({name} line, {h0:.0f}-{h1:.0f}°)"))
    return boxes


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

    With `mos`, frames of 'Steep edge' and 'Gap' areas whose footprint the mosaic already covers
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
            if mos is not None and label.startswith(("Steep edge", "Gap")):
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
#  DEBUG LOG – every action, with timing and thread (menu: Debug → Debug log)
# =============================================================================

import logging

_DBG = {"on": False, "logger": None, "path": None}


def app_data_dir():
    base = os.environ.get("LOCALAPPDATA") or os.path.join(os.path.expanduser("~"), ".local", "share")
    d = os.path.join(base, "HorizonMapper"); os.makedirs(d, exist_ok=True)
    return d


def load_settings():
    try:
        with open(os.path.join(app_data_dir(), "settings.json"), encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return {}


def save_settings(**kw):
    st = load_settings(); st.update(kw)
    try:
        with open(os.path.join(app_data_dir(), "settings.json"), "w", encoding="utf-8") as fh:
            json.dump(st, fh, indent=1)
    except Exception:
        pass


def debug_enable(on=True):
    """Switch the debug log on/off. The file goes to <app data>/logs/debug_<date_time>.log."""
    _DBG["on"] = bool(on)
    if on and _DBG["logger"] is None:
        d = os.path.join(app_data_dir(), "logs"); os.makedirs(d, exist_ok=True)
        path = os.path.join(d, time.strftime("debug_%Y-%m-%d_%H%M%S.log"))
        lg = logging.getLogger("horizonmapper.debug"); lg.setLevel(logging.DEBUG); lg.propagate = False
        h = logging.FileHandler(path, encoding="utf-8")
        h.setFormatter(logging.Formatter("%(asctime)s.%(msecs)03d [%(threadName)s] %(message)s", "%H:%M:%S"))
        lg.addHandler(h)
        _DBG["logger"], _DBG["path"] = lg, path
        lg.debug(f"=== {APP_NAME} {APP_VERSION} debug log – Python {sys.version.split()[0]}, {sys.platform} ===")
    if _DBG["logger"]:
        _DBG["logger"].debug(f"debug log {'ON' if on else 'OFF'}")
    return _DBG["path"]


def dlog(msg):
    if _DBG["on"] and _DBG["logger"] is not None:
        try:
            _DBG["logger"].debug(msg)
        except Exception:
            pass


def dlog_exc(where):
    dlog(f"EXCEPTION in {where}:\n{traceback.format_exc()}")


# =============================================================================
#  CAPTURE – drive NINA through its "Advanced API" plugin
# =============================================================================

class NinaApiError(RuntimeError):
    pass


class NinaApi:
    """Minimal client for the NINA Advanced API plugin (v2, http://host:port/v2/api)."""

    def __init__(self, host="localhost", port=1888, timeout=30.0):
        self.base = f"http://{host}:{int(port)}/v2/api"
        self.timeout = timeout

    def get(self, path, timeout=None, **params):
        import urllib.request, urllib.parse
        q = {k: (str(v).lower() if isinstance(v, bool) else v) for k, v in params.items() if v is not None}
        url = self.base + path + ("?" + urllib.parse.urlencode(q) if q else "")
        t0 = time.time()
        dlog(f"API → {path} {q if q else ''}")
        try:
            with urllib.request.urlopen(url, timeout=timeout or self.timeout) as r:
                data = json.loads(r.read().decode("utf-8", "replace"))
        except Exception as e:                      # connection refused, timeout, HTTP error, bad JSON
            dlog(f"API ✗ {path} after {time.time() - t0:.2f} s: {type(e).__name__}: {e}")
            raise NinaApiError(f"NINA API {path}: {e}") from e
        if _DBG["on"]:
            resp = data.get("Response")
            if isinstance(resp, dict) and path == "/equipment/mount/info":
                resp = {k: resp.get(k) for k in ("Slewing", "TrackingEnabled", "AtPark", "RightAscension",
                                                 "Declination", "Altitude", "Azimuth", "SideOfPier")}
            dlog(f"API ← {path} {time.time() - t0:.2f} s Success={data.get('Success')} "
                 f"Error={data.get('Error')!r} Response={str(resp)[:200]}")
        if not data.get("Success", False):
            if data.get("Error"):
                raise NinaApiError(f"NINA API {path}: {data['Error']}")
            raise NinaApiError(f"NINA API {path}: NINA reported that the command failed"
                               + (f" (its message: '{data.get('Response')}')" if data.get("Response") else ""))
        return data.get("Response")

    # convenience wrappers
    def version(self):        return self.get("/version")
    def mount(self):          return self.get("/equipment/mount/info")
    def rotator(self):        return self.get("/equipment/rotator/info")
    def camera(self):         return self.get("/equipment/camera/info")
    def profile(self):        return self.get("/profile/show", active=True)

    def stop_all(self):
        """Stop every movement and exposure at once; errors are collected, not raised."""
        errs = []
        for path in ("/equipment/mount/slew/stop", "/equipment/rotator/stop-move",
                     "/equipment/camera/abort-exposure"):
            try:
                self.get(path, timeout=5)
            except NinaApiError as e:
                errs.append(str(e))
        return errs


def altaz_to_radec_date(alt, az, lat, lst_hours):
    """Alt/az → apparent RA/Dec of date (JNow), degrees (refraction ignored)."""
    A, h, p = math.radians(az), math.radians(alt), math.radians(lat)
    sd = math.sin(p) * math.sin(h) + math.cos(p) * math.cos(h) * math.cos(A)
    dec = math.asin(max(-1.0, min(1.0, sd)))
    ha = math.atan2(-math.sin(A) * math.cos(h),
                    math.cos(p) * math.sin(h) - math.sin(p) * math.cos(h) * math.cos(A))
    return math.degrees((math.radians(lst_hours * 15.0) - ha) % (2 * math.pi)), math.degrees(dec)


def altaz_to_radec_j2000(alt, az, lat, lst_hours, t_utc=None):
    """Convert a topocentric alt/az direction to J2000 RA/Dec in degrees (good to ~1 arcmin:
    precession applied; nutation, aberration and refraction ignored - fine for 2-3° frames)."""
    t_utc = time.time() if t_utc is None else t_utc
    A, h, p = math.radians(az), math.radians(alt), math.radians(lat)
    sd = math.sin(p) * math.sin(h) + math.cos(p) * math.cos(h) * math.cos(A)
    dec = math.asin(max(-1.0, min(1.0, sd)))
    ha = math.atan2(-math.sin(A) * math.cos(h),
                    math.cos(p) * math.sin(h) - math.sin(p) * math.cos(h) * math.cos(A))
    ra = (math.radians(lst_hours * 15.0) - ha) % (2 * math.pi)              # of date
    # precess date -> J2000 (IAU 1976), rotating back with the inverse matrix
    T = ((t_utc / 86400.0 + 2440587.5) - 2451545.0) / 36525.0
    as2r = math.pi / 180 / 3600
    zeta = (2306.2181 * T + 0.30188 * T * T + 0.017998 * T ** 3) * as2r
    z = (2306.2181 * T + 1.09468 * T * T + 0.018203 * T ** 3) * as2r
    th = (2004.3109 * T - 0.42665 * T * T - 0.041833 * T ** 3) * as2r
    cz, sz, cZ, sZ, ct, st = math.cos(zeta), math.sin(zeta), math.cos(z), math.sin(z), math.cos(th), math.sin(th)
    P = [[cz * ct * cZ - sz * sZ, -sz * ct * cZ - cz * sZ, -st * cZ],
         [cz * ct * sZ + sz * cZ, -sz * ct * sZ + cz * cZ, -st * sZ],
         [cz * st, -sz * st, ct]]                                            # J2000 -> date
    v = [math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)]
    w = [sum(P[r][c] * v[r] for r in range(3)) for c in range(3)]           # P^T v: date -> J2000
    return math.degrees(math.atan2(w[1], w[0])) % 360, math.degrees(math.asin(max(-1, min(1, w[2]))))


def sun_separation(alt, az, lat, lon, t_utc=None):
    """Angle between a pointing (alt/az) and the Sun, degrees."""
    s_alt, s_az = sun_altaz(time.time() if t_utc is None else t_utc, lat, lon)
    a1, a2 = math.radians(alt), math.radians(s_alt)
    c = math.sin(a1) * math.sin(a2) + math.cos(a1) * math.cos(a2) * math.cos(math.radians(az - s_az))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


class CaptureRunner:
    """Shoots a list of frames through NINA, one at a time: slew → rotator → capture → find file.

    Every wait polls NINA, so pause() takes effect at once: it stops the mount, the rotator and
    the exposure, and the run holds until resume(); the interrupted frame is then repeated from
    its start. Any error or timeout pauses the run (it is never skipped silently). stop() ends it.
    Events go to `emit(kind, **data)`: 'log', 'state', 'frame' (index, status, path), 'done'.
    """

    def __init__(self, api, frames, emit, exposure=1.0, gain=None, binning=None, target="Horizon mapping",
                 settle=1.0, sun_min=30.0, floor=0.0, unpark=True, park=False, image_folder="",
                 slew_timeout=180.0, rot_timeout=120.0, rot_tol=0.5, start_index=0, poll=0.5, save_wait=5.0,
                 order=None, status=None, paths=None, names=None):
        self.api, self.emit = api, emit
        self.frames = list(frames)
        # rotator target for every frame (a plan only lists the moves; carry the angle forward)
        self.rot = []; cur = None
        for _, _, r in self.frames:
            cur = r if r is not None else cur
            self.rot.append(cur)
        self.exposure, self.gain, self.binning, self.target = exposure, gain, binning, target
        self.settle, self.sun_min, self.floor = settle, sun_min, floor
        self.unpark, self.park, self.image_folder = unpark, park, image_folder
        self.slew_timeout, self.rot_timeout, self.rot_tol, self.poll = slew_timeout, rot_timeout, rot_tol, poll
        # frames to shoot, in this order (default: all from start_index); the rest keep their state
        self.order = list(order) if order is not None else list(range(start_index, len(self.frames)))
        self.pos = 0
        self.index = self.order[0] if self.order else len(self.frames)
        n = len(self.frames)
        self.status = list(status) if status else ["pending"] * n
        self.paths = list(paths) if paths else [None] * n
        self.names = list(names) if names else [None] * n   # file name reported by NINA, if any
        self.save_wait = save_wait
        self.t_start = time.time()
        self._pause = threading.Event(); self._stop = threading.Event(); self._skip = threading.Event()
        self.state = "idle"; self.thread = None; self.pause_reason = ""
        self._gen = 0; self._step_gen = 0            # pause counter: any pause interrupts the running step
        self._dir_hint = None

    # ---------------------------------------------------------------- control (any thread)
    def start(self):
        self.thread = threading.Thread(target=self._run, daemon=True); self.thread.start()

    def pause(self, reason="Paused by user"):
        """Stop all movement and exposure now and hold."""
        dlog(f"runner.pause({reason!r}) called, state {self.state}")
        if self.state in ("done", "stopped"):
            return
        self.pause_reason = reason
        self._gen += 1
        self._pause.set()
        errs = self.api.stop_all()
        self._set_state("paused")
        self.emit("log", msg=f"PAUSED: {reason}" + (f" (stop commands: {'; '.join(errs)})" if errs else
                                                    " – mount, rotator and camera stopped."))

    def resume(self):
        dlog(f"runner.resume() called, state {self.state}")
        if self._pause.is_set():
            self._pause.clear(); self._set_state("running")
            self.emit("log", msg=f"Resumed at frame {self.index + 1}.")

    def skip(self):
        """While paused: skip the current frame and continue with the next one when resumed."""
        dlog(f"runner.skip() called at frame {self.index + 1}")
        self._skip.set()

    def stop(self):
        dlog(f"runner.stop() called, state {self.state}")
        self._stop.set(); self._pause.clear()
        self.api.stop_all()
        self.emit("log", msg="Stopping: mount, rotator and camera stopped.")

    # ---------------------------------------------------------------- internals
    def _set_state(self, st):
        dlog(f"runner state {self.state} → {st} (frame {self.index + 1})")
        self.state = st; self.emit("state", state=st)

    class _Interrupted(Exception):
        pass

    def _interrupted(self):
        return self._stop.is_set() or self._pause.is_set() or self._gen != self._step_gen

    def _check(self):
        if self._interrupted():
            dlog("step interrupted (pause/stop)")
            raise self._Interrupted()

    def _wait(self, seconds):
        end = time.time() + seconds
        while time.time() < end:
            self._check(); time.sleep(min(self.poll, max(end - time.time(), 0.01)))

    def _wait_until(self, cond, timeout, what):
        t0 = time.time(); end = t0 + timeout
        dlog(f"waiting for {what} (timeout {timeout:.0f} s)")
        while True:
            self._check()
            if cond():
                dlog(f"… {what}: done after {time.time() - t0:.1f} s")
                return
            if time.time() > end:
                dlog(f"… {what}: TIMEOUT after {timeout:.0f} s")
                raise NinaApiError(f"timeout after {timeout:.0f} s waiting for {what}")
            time.sleep(self.poll)

    def _hold(self):
        """Block while paused. Returns False if the run was stopped."""
        while self._pause.is_set() and not self._stop.is_set():
            time.sleep(0.1)
        return not self._stop.is_set()

    def _locate(self, filename):
        """Find a file NINA saved, by name, under its image folder (newest folders first)."""
        if not filename:
            return None
        if self._dir_hint and os.path.isfile(os.path.join(self._dir_hint, filename)):
            return os.path.join(self._dir_hint, filename)
        root = self.image_folder
        if not root or not os.path.isdir(root):
            return None
        for d, _, files in os.walk(root):
            if filename in files:
                self._dir_hint = d
                return os.path.join(d, filename)
        return None

    def _prepare(self):
        m = self.api.mount()
        if not m.get("Connected"):
            raise NinaApiError("the mount is not connected in NINA")
        if m.get("AtPark"):
            if not self.unpark:
                raise NinaApiError("the mount is parked (enable 'Unpark at start' or unpark it in NINA)")
            self.emit("log", msg="Unparking the mount…")
            self.api.get("/equipment/mount/unpark", timeout=120)
            self._wait_until(lambda: not self.api.mount().get("AtPark"), 120, "unpark")
        # Tracking is not needed: exposures are milliseconds long and every slew is aimed at where
        # the target is at that moment. So the tracking command (which some drivers, e.g. the ZWO
        # ASI Mount, answer with "failed") is only used if a slew is refused (see _slew).
        dlog(f"tracking at start: {m.get('TrackingEnabled')}")
        self._tracking_info_shown = False
        if not self.api.camera().get("Connected"):
            raise NinaApiError("the camera is not connected in NINA")
        if any(r is not None for r in self.rot) and not self.api.rotator().get("Connected"):
            raise NinaApiError("the plan uses the rotator, but no rotator is connected in NINA")
        if self.binning:
            self.api.get("/equipment/camera/set-binning", binning=f"{int(self.binning)}x{int(self.binning)}")
        return m

    def _frame(self, i, site):
        alt, az, _ = self.frames[i]
        dlog(f"--- frame {i + 1}/{len(self.frames)}: alt {alt:.2f}°, az {az:.2f}°, rotator {self.rot[i]}")
        lat, lon = site
        if alt < self.floor - 1e-6:
            return "skipped", f"below the lowest altitude ({alt:.1f}° < {self.floor:.1f}°)"
        sep = sun_separation(alt, az, lat, lon)
        if sep < self.sun_min:
            return "skipped", f"too close to the Sun ({sep:.0f}° < {self.sun_min:.0f}°)"
        # 1 – slew (target computed now, from the mount's own sidereal time)
        m = self.api.mount()
        lst_at_slew = float(m["SiderealTime"])
        # the slew path must keep away from the Sun too: go via a safe intermediate position if needed
        sun = sun_altaz(time.time(), lat, lon)
        if sun[0] > -self.sun_min and m.get("RightAscension") is not None:
            ha_m = (lst_at_slew * 15.0 - float(m["RightAscension"]) * 15.0 + 180) % 360 - 180
            here = hadec_to_altaz(ha_m, float(m["Declination"]), lat)
            psep = path_sun_sep(slew_path(here, (alt, az), lat), sun)
            if psep < self.sun_min:
                w = safe_waypoint(here, (alt, az), lat, sun, self.sun_min)
                if w is None:
                    return "skipped", f"the slew there would pass {psep:.0f}° from the Sun (no safe route)"
                self.emit("log", msg=f"Frame {i + 1}: the direct slew would pass {psep:.0f}° from the Sun – "
                                     f"going via alt {w[0]:.0f}°, az {w[1]:.0f}° first.")
                wra, wdec = altaz_to_radec_j2000(w[0], w[1], lat, lst_at_slew)
                self.api.get("/equipment/mount/slew", ra=round(wra, 5), dec=round(wdec, 5),
                             waitForResult=True, timeout=self.slew_timeout + 15)
                self._check()
                m = self.api.mount(); lst_at_slew = float(m["SiderealTime"])
        ra, dec = altaz_to_radec_j2000(alt, az, lat, lst_at_slew)
        self.emit("frame", index=i, status="slewing")
        # NINA's own slew handling decides when the slew is over: the request waits for it, in a
        # background thread, while this thread keeps checking for PAUSE / Stop (which stop the slew).
        box = {}

        def do_slew():
            try:
                box["r"] = self.api.get("/equipment/mount/slew", ra=round(ra, 5), dec=round(dec, 5),
                                        waitForResult=True, timeout=self.slew_timeout + 15)
            except NinaApiError as e:
                # some drivers only slew while tracking: switch it on (its answer may say "failed"
                # even when it works) and try the slew once more
                try:
                    if not self.api.mount().get("TrackingEnabled", True) and not self._stop.is_set() \
                            and not self._pause.is_set():
                        dlog(f"slew refused ({e}) while not tracking: switching tracking on, retrying once")
                        try:
                            self.api.get("/equipment/mount/tracking", mode=0)
                        except NinaApiError as e2:
                            dlog(f"tracking command answered: {e2}")
                        time.sleep(1.0)
                        box["r"] = self.api.get("/equipment/mount/slew", ra=round(ra, 5), dec=round(dec, 5),
                                                waitForResult=True, timeout=self.slew_timeout + 15)
                        return
                except NinaApiError as e3:
                    e = e3
                box["e"] = e
        th = threading.Thread(target=do_slew, daemon=True); th.start()
        self._wait_until(lambda: not th.is_alive(), self.slew_timeout + 20, "the slew")
        self._check()                                   # a pause stops the slew: not an arrival
        if "e" in box:
            raise box["e"]
        # some mounts report the end a moment early: wait until it is not slewing and stands still
        last = [None]

        def settled():
            """Not slewing, and standing still: RA/Dec unchanged (tracking) or hour angle/Dec
            unchanged (not tracking) since the previous check."""
            mm = self.api.mount()
            ra = float(mm.get("RightAscension") or 0) * 15.0; dec = float(mm.get("Declination") or 0)
            ha = float(mm.get("SiderealTime") or 0) * 15.0 - ra
            prev = last[0]; last[0] = (ra, ha, dec)
            if prev is None or mm.get("Slewing"):
                return False
            k = math.cos(math.radians(dec))
            d_ra = abs((ra - prev[0] + 180) % 360 - 180) * k
            d_ha = abs((ha - prev[1] + 180) % 360 - 180) * k
            return abs(dec - prev[2]) < 0.003 and min(d_ra, d_ha) < 0.003
        self._wait_until(settled, 30, "the mount to stand still after the slew")
        m = self.api.mount()
        # arrival: compare RA/Dec of date (some drivers report alt/az poorly); recompute the target for
        # the current sidereal time, so a mount without tracking does not count as off target
        if m.get("RightAscension") is not None and m.get("Declination") is not None:
            ra_d, dec_d = altaz_to_radec_date(alt, az, lat, float(m["SiderealTime"]))
            ra_m, dec_m = float(m["RightAscension"]) * 15.0, float(m["Declination"])
            if m.get("TrackingEnabled"):          # tracking: the mount stays on the RA it slewed to
                ra_d = altaz_to_radec_date(alt, az, lat, float(lst_at_slew))[0]
            err = max(abs(dec_m - dec_d), abs((ra_m - ra_d + 180) % 360 - 180) * math.cos(math.radians(dec_d)))
            where = f"RA {ra_m / 15:.3f} h, Dec {dec_m:.2f}° instead of RA {ra_d / 15:.3f} h, Dec {dec_d:.2f}°"
        else:
            err = max(abs(float(m.get("Altitude", alt)) - alt),
                      abs((float(m.get("Azimuth", az)) - az + 180) % 360 - 180) * math.cos(math.radians(alt)))
            where = f"alt {m.get('Altitude')}°, az {m.get('Azimuth')}° instead of {alt:.1f}°, {az:.1f}°"
        if err > 2.0:
            raise NinaApiError(f"after the slew the mount is at {where} – slew stopped or limits hit?")
        if not m.get("TrackingEnabled") and not getattr(self, "_tracking_info_shown", False):
            self._tracking_info_shown = True
            self.emit("log", msg="The mount is not tracking – that is fine here: the exposures are very short "
                                 "and every slew is aimed at the current sky position.")
        # 2 – rotator
        tgt = self.rot[i]
        if tgt is not None:
            r = self.api.rotator()
            if abs((float(r.get("MechanicalPosition", 0)) - tgt + 180) % 360 - 180) > self.rot_tol:
                self.emit("frame", index=i, status="rotating")
                self.api.get("/equipment/rotator/move-mechanical", position=round(tgt, 2))
                time.sleep(self.poll)

                def rot_done():
                    rr = self.api.rotator()
                    return (not rr.get("IsMoving")) and \
                        abs((float(rr.get("MechanicalPosition", 0)) - tgt + 180) % 360 - 180) <= max(self.rot_tol, 1.0)
                self._wait_until(rot_done, self.rot_timeout, "the rotator")
        self._wait(self.settle)
        # 3 – capture (NINA saves the FITS with its usual headers: CENTALT, CENTAZ, ROTATOR, …)
        self.emit("frame", index=i, status="exposing")
        n0 = int(self.api.get("/image-history", count=True) or 0)
        self.api.get("/equipment/camera/capture", duration=self.exposure, gain=self.gain, save=True,
                     imageType="LIGHT", targetName=self.target, omitImage=True, skipAutoStretch=True)

        def captured():
            try:
                return self.api.get("/equipment/camera/capture", getResult=True, omitImage=True) \
                    != "Capture already in progress"
            except NinaApiError:
                return False
        self._wait_until(captured, self.exposure + 120, "the exposure and download")
        # 4 – NINA saves in the background and the file can take a while to appear: look briefly,
        #     then move on; the folder scan (match_late_files) picks the file up when it arrives
        name = None
        end = time.time() + self.save_wait
        while time.time() < end:
            self._check()
            if int(self.api.get("/image-history", count=True) or 0) > n0:
                h = self.api.get("/image-history", index=n0)
                entry = h[0] if isinstance(h, list) and h else (h if isinstance(h, dict) else {})
                name = entry.get("Filename"); break
            time.sleep(self.poll)
        self.names[i] = name
        path = self._locate(name)
        dlog(f"frame {i + 1}: NINA file name {name!r}, found in the image folder: {path!r}")
        return "done", path

    def _run(self):
        self._set_state("running")
        try:
            m = self._prepare()
        except self._Interrupted:
            m = None
        except NinaApiError as e:
            self.emit("log", msg=f"Cannot start: {e}"); self._set_state("stopped"); self.emit("done", ok=False); return
        site = (float(m["SiteLatitude"]), float(m["SiteLongitude"])) if m else (0.0, 0.0)
        while self.pos < len(self.order) and not self._stop.is_set():
            self.index = self.order[self.pos]
            if self._pause.is_set():
                if not self._hold():
                    break
                if m is None:                           # paused during preparation: prepare again
                    try:
                        m = self._prepare(); site = (float(m["SiteLatitude"]), float(m["SiteLongitude"]))
                    except self._Interrupted:
                        continue
                    except NinaApiError as e:
                        self.pause(f"error – {e}"); continue
            if self._skip.is_set():
                self._skip.clear()
                self.status[self.index] = "skipped"; self.emit("frame", index=self.index, status="skipped")
                self.emit("log", msg=f"Frame {self.index + 1} skipped by user."); self.pos += 1; continue
            i = self.index
            self._step_gen = self._gen
            try:
                st, info = self._frame(i, site)
            except self._Interrupted:
                self.status[i] = "pending"; self.emit("frame", index=i, status="pending"); continue
            except NinaApiError as e:
                dlog(f"frame {i + 1} error: {e}")
                if self._interrupted():                     # caused by our own stop commands
                    self.status[i] = "pending"; self.emit("frame", index=i, status="pending"); continue
                self.status[i] = "error"; self.emit("frame", index=i, status="error")
                self.pause(f"error at frame {i + 1} – {e}"); continue
            except Exception as e:                      # unexpected: pause, never run on blindly
                dlog_exc(f"frame {i + 1}")
                self.status[i] = "error"; self.emit("frame", index=i, status="error")
                self.pause(f"unexpected error at frame {i + 1} – {type(e).__name__}: {e}"); continue
            self.status[i] = st
            if st == "done":
                self.paths[i] = info
            self.emit("frame", index=i, status=st, path=info if st == "done" else None,
                      name=self.names[i] if st == "done" else None)
            if st == "skipped":
                self.emit("log", msg=f"Frame {i + 1} skipped: {info}")
            self.pos += 1
        stopped = self._stop.is_set()
        if self.park and not stopped:
            try:
                self.emit("log", msg="Parking the mount…"); self.api.get("/equipment/mount/park", timeout=10)
            except NinaApiError as e:
                self.emit("log", msg=f"Park failed: {e}")
        self._set_state("stopped" if stopped else "done")
        self.emit("done", ok=not stopped)


def read_nina_sequence(path):
    """Frames [(alt, az, rotator or None)] from a NINA sequence or template (as saved by this app,
    or any sequence using 'Slew to Alt/Az', 'Move rotator mechanical' and 'Take exposure').
    Also returns the first exposure's settings: dict(exposure, gain, binning)."""
    with open(path, encoding="utf-8-sig") as fh:
        root = json.load(fh)
    frames = []; settings = {}
    st = {"alt": None, "az": None, "rot": None}

    def dms(d, m, s_):
        d, m, s_ = float(d or 0), float(m or 0), float(s_ or 0)
        sign = -1 if (d < 0 or m < 0 or s_ < 0) else 1
        return sign * (abs(d) + abs(m) / 60 + abs(s_) / 3600)

    def walk(node):
        if isinstance(node, dict):
            t = str(node.get("$type", ""))
            if "SlewScopeToAltAz" in t:
                c = node.get("Coordinates") or {}
                st["alt"] = dms(c.get("AltDegrees"), c.get("AltMinutes"), c.get("AltSeconds"))
                st["az"] = dms(c.get("AzDegrees"), c.get("AzMinutes"), c.get("AzSeconds")) % 360
            elif "MoveRotatorMechanical" in t:
                st["rot"] = float(node.get("MechanicalPosition", 0))
            elif "TakeExposure" in t or "TakeManyExposures" in t:
                if st["alt"] is not None:
                    frames.append((st["alt"], st["az"], st["rot"])); st["rot"] = None
                if not settings:
                    b = node.get("Binning") or {}
                    settings.update(exposure=node.get("ExposureTime"), gain=node.get("Gain"), binning=b.get("X"))
            items = node.get("Items")
            if isinstance(items, dict):
                for it in items.get("$values", []):
                    walk(it)
            for k, v in node.items():                 # templates nest containers in other keys too
                if k not in ("Items", "Parent", "Coordinates") and isinstance(v, (dict, list)):
                    walk(v)
        elif isinstance(node, list):
            for it in node:
                walk(it)
    walk(root)
    return frames, settings


def rows_of(frames):
    """Group consecutive frames of equal altitude into rows: [(alt, [az, …])]."""
    grid = []
    for alt, az, _ in frames:
        if grid and abs(grid[-1][0] - alt) < 1e-3:
            grid[-1][1].append(az)
        else:
            grid.append((alt, [az]))
    return grid


def fits_complete(path):
    """True when a FITS file is fully written: its size covers the header and all the image data
    (NINA writes files in the background, so a file can exist but still be short or empty)."""
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as fh:
            h, n = read_fits_header(fh)
        naxis = int(h.get("NAXIS", 0))
        data = abs(int(h.get("BITPIX", 16))) // 8
        for k in range(1, naxis + 1):
            data *= int(h.get(f"NAXIS{k}", 0))
        return size >= n + data and data > 0
    except Exception:
        return False


def frame_thumbnail(path, K, lat, half_w, half_h, px=110, cache_dir=None):
    """Project one captured frame onto a small alt/az grid around its own centre (the same maths as
    the mosaic). Returns dict(az0, az1, alt0, alt1, img) with img rows from top (high) to bottom."""
    meta, img = _load_one(path, cache_dir)
    if meta is None:
        raise ValueError(img)
    fr = Frame()
    fr.path, fr.alt, fr.az, fr.rot, fr.pier, fr.scale, fr.img, fr.t = (
        path, meta["alt"], meta["az"], meta["rot"], meta["pier"], meta["scale"], img, meta.get("t"))
    lat = meta.get("lat") if meta.get("lat") is not None else lat
    cw = half_w / max(math.cos(math.radians(min(abs(fr.alt) + half_h, 85.0))), 0.05)
    n_az = px; n_alt = max(int(px * half_h / max(cw * math.cos(math.radians(fr.alt)), 1e-3)), 8)
    az = np.linspace(fr.az - cw, fr.az + cw, n_az); alt = np.linspace(fr.alt + half_h, fr.alt - half_h, n_alt)
    AZ, ALT = np.meshgrid(az, alt)
    v = sample_frame(img, fr, frame_delta(fr, K, lat), fr.scale, AZ, ALT)
    return dict(az0=fr.az - cw, az1=fr.az + cw, alt0=fr.alt - half_h, alt1=fr.alt + half_h, img=v)


def thumbnail_rgba(img):
    """Stretched greyscale RGBA PIL image of a thumbnail array (transparent where there is no data)."""
    from PIL import Image
    ok = np.isfinite(img)
    out = np.zeros(img.shape + (4,), np.uint8)
    if ok.any():
        v = np.log(np.maximum(img[ok] - np.nanmin(img) + 1.0, 1.0))
        lo, hi = np.percentile(v, [1, 99.5]); hi = hi if hi > lo else lo + 1
        g = (np.clip((v - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)
        out[ok, 0] = out[ok, 1] = out[ok, 2] = g; out[ok, 3] = 255
    return Image.fromarray(out, "RGBA")


def match_late_files(folder, frames, waiting, names=None, since=0.0, seen=None, tol=0.6):
    """Find FITS files that arrived in NINA's image folder after the frame was taken.

    `waiting`: indices of frames taken but still without a file. A file is matched by the name NINA
    reported, or else by the pointing in its header (CENTALT/CENTAZ within `tol` degrees, nearest
    frame). Only files modified after `since` are looked at; `seen` (dict path → header pointing or
    None) caches headers between calls. Returns {frame index: path}."""
    seen = {} if seen is None else seen
    waiting = set(waiting)
    if not waiting or not folder or not os.path.isdir(folder):
        return {}
    want = {os.path.basename(n): i for i, n in enumerate(names or []) if n and i in waiting}
    found = {}
    for d, _, files in os.walk(folder):
        for fn in files:
            if not fn.lower().endswith((".fits", ".fit", ".fts")):
                continue
            p = os.path.join(d, fn)
            if fn in want and want[fn] in waiting:
                found[want[fn]] = p; waiting.discard(want[fn]); continue
            if p in seen:
                continue
            try:
                if os.path.getmtime(p) < since:
                    seen[p] = None; continue
                with open(p, "rb") as fh:
                    h, _ = read_fits_header(fh)
                seen[p] = (float(h["CENTALT"]), float(h["CENTAZ"]))
            except Exception:
                continue          # still being written, or not one of ours: try again next time
    taken = set(found.values())
    for p, pos in seen.items():
        if not pos or p in taken or not waiting:
            continue
        alt, az = pos
        best = min(waiting, key=lambda i: max(abs(frames[i][0] - alt),
                                              abs((frames[i][1] - az + 180) % 360 - 180) * math.cos(math.radians(alt))))
        a2, z2 = frames[best][0], frames[best][1]
        if max(abs(a2 - alt), abs((z2 - az + 180) % 360 - 180) * math.cos(math.radians(alt))) <= tol:
            found[best] = p; waiting.discard(best); taken.add(p)
    return found


CAPTURE_FILE_VERSION = 1


def save_capture_file(path, frames, status, paths, names, fw, fh, source="", settings=None,
                      horizon_files=None):
    """Save a capture's progress (frames, each frame's state, file paths and names, settings) as JSON,
    so it can be resumed later with 'Load saved capture'. Written to a temporary file first."""
    data = dict(app=APP_NAME, version=APP_VERSION, file_version=CAPTURE_FILE_VERSION,
                saved=time.strftime("%Y-%m-%d %H:%M:%S"), source=source, fw=fw, fh=fh,
                horizon_files=list(horizon_files or []), settings=dict(settings or {}),
                frames=[dict(n=i + 1, alt=round(a, 5), az=round(z, 5), rot=r, status=st, path=p, name=nm)
                        for i, ((a, z, r), st, p, nm) in enumerate(zip(frames, status, paths, names))])
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh_:
        json.dump(data, fh_, indent=1)
    os.replace(tmp, path)
    return path


def load_capture_file(path):
    """Read a saved capture. Frames whose saved file no longer exists are set back to 'pending'.
    Returns dict(frames, status, paths, names, fw, fh, source, settings, horizon_files, missing)."""
    with open(path, encoding="utf-8") as fh_:
        d = json.load(fh_)
    if d.get("file_version", 0) > CAPTURE_FILE_VERSION:
        raise ValueError("this capture file was saved by a newer version of the app")
    frames, status, paths, names, missing = [], [], [], [], []
    for i, f_ in enumerate(d["frames"]):
        frames.append((float(f_["alt"]), float(f_["az"]), f_.get("rot")))
        st, p = f_.get("status", "pending"), f_.get("path")
        if st in ("slewing", "rotating", "exposing"):          # saved mid-frame: shoot it again
            st = "pending"
        if st == "done" and not (p and os.path.isfile(p)):
            missing.append(i); st = "pending"; p = None
        status.append(st); paths.append(p); names.append(f_.get("name"))
    return dict(frames=frames, status=status, paths=paths, names=names, fw=float(d["fw"]), fh=float(d["fh"]),
                source=d.get("source", ""), settings=d.get("settings", {}),
                horizon_files=d.get("horizon_files", []), missing=missing, saved=d.get("saved", ""))


def frames_to_shoot(status, first=1, last=None, reshoot=False):
    """Indices of the frames a run should shoot: within frames first..last (1-based, inclusive),
    skipping those already taken unless `reshoot`."""
    n = len(status)
    last = n if last is None else min(int(last), n)
    first = max(int(first), 1)
    return [i for i in range(first - 1, last) if reshoot or status[i] != "done"]


def write_capture_list(paths, folder=None):
    """Write the captured FITS paths to a list file that the Build tab accepts like a folder."""
    folder = folder or os.path.join(os.path.dirname(default_cache_dir()), "captures")
    os.makedirs(folder, exist_ok=True)
    fn = os.path.join(folder, time.strftime("capture_%Y-%m-%d_%H%M%S.txt"))
    with open(fn, "w", encoding="utf-8") as f:
        f.write("# Horizon Mapper capture list – one FITS file per line\n")
        for p in paths:
            if p and os.path.isabs(p):
                f.write(p + "\n")
    return fn


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
            if fo.lower().endswith(".txt"):          # capture list: one FITS path per line
                with open(fo, encoding="utf-8") as fh:
                    out += [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")
                            and os.path.isfile(ln.strip())]
            else:
                out.append(fo)
            continue
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


def widen_steep_edges(az, alt, side_deg=0.0, steep_slope=2.0):
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

def run_gui(selftest=False):
    import types
    import tkinter as _tk_mod
    from tkinter import ttk as _ttk_mod, filedialog, messagebox
    # Private copies of the tkinter / ttk namespaces: the hover and click-logging wrappers below
    # replace classes here only. Never patch the tkinter modules themselves – ttk builds Combobox
    # and Spinbox on ttk.Entry internally, and a patched Entry breaks them.
    tk = types.SimpleNamespace(**{k: getattr(_tk_mod, k) for k in dir(_tk_mod) if not k.startswith("__")})
    ttk = types.SimpleNamespace(**{k: getattr(_ttk_mod, k) for k in dir(_ttk_mod) if not k.startswith("__")})
    try:
        from PIL import ImageTk
    except ImportError:
        ImageTk = None

    debug_enable(load_settings().get("debug_log", DEV_MODE))   # default: on during development

    # hover explanations: shown after a short pause over any setting, choice or button
    _selftest_hovers = []

    class Hover:
        def __init__(self, w, text):
            self.w, self.text, self.tw, self.job = w, text, None, None
            _selftest_hovers.append(self)
            for ev, fn in (("<Enter>", self.sched), ("<Leave>", self.hide), ("<ButtonPress>", self.hide)):
                try:
                    w.bind(ev, fn, add="+")
                except Exception:
                    pass

        def sched(self, e=None):
            self.hide(); self.job = self.w.after(600, self.show)

        def show(self):
            try:
                x, y = self.w.winfo_pointerx() + 14, self.w.winfo_pointery() + 18
                self.tw = tk.Toplevel(self.w); self.tw.wm_overrideredirect(True); self.tw.wm_geometry(f"+{x}+{y}")
                tk.Label(self.tw, text=self.text, background="#ffffe8", relief="solid", borderwidth=1,
                         wraplength=360, justify="left", padx=6, pady=3).pack()
            except Exception:
                self.tw = None

        def hide(self, e=None):
            if self.job:
                try:
                    self.w.after_cancel(self.job)
                except Exception:
                    pass
                self.job = None
            if self.tw is not None:
                try:
                    self.tw.destroy()
                except Exception:
                    pass
                self.tw = None
    no_tip = set()
    last_label = {}                       # parent → hover text of the last label made in it

    def tip_for(text):
        t = TIPS.get(text) or TIPS.get(str(text).strip())
        return t

    def attach_tip(w, text, kind):
        t = tip_for(text) if text else None
        if t:
            Hover(w, t); return t
        if text and kind != "label":
            no_tip.add(f"{kind}: {str(text).strip()}")
        return None

    def _hovered(cls, kind):
        class H(cls):
            def __init__(self, *a, hover=None, **kw):
                super().__init__(*a, **kw)
                parent = a[0] if a else kw.get("master")
                if kind == "label":
                    t = attach_tip(self, kw.get("text"), kind)
                    last_label[str(parent)] = t
                elif kind in ("entry", "spinbox", "combobox", "listbox"):
                    t = hover or last_label.get(str(parent))  # an input takes its label's explanation
                    if t:
                        Hover(self, t)
                    else:
                        no_tip.add(f"{kind} in {parent}")
                else:
                    if hover:
                        Hover(self, hover)
                    else:
                        attach_tip(self, kw.get("text"), kind)
        return H
    ttk.Label = _hovered(ttk.Label, "label"); ttk.Entry = _hovered(ttk.Entry, "entry")
    ttk.Spinbox = _hovered(ttk.Spinbox, "spinbox"); ttk.Combobox = _hovered(ttk.Combobox, "combobox")

    def draw_sun_overlay(cv, X, Y, risk, frames, fw, fh, sun_min, x0, x1, tag="sun"):
        """Sun at the start (filled) and end (ring) of the run, its track, the keep-away zone around
        both, frames too close to the Sun (red) and frames whose slew passes near it (red cross)."""
        half = (x1 - x0) / 2

        def poly(pts, **kw):                       # split where the azimuth wraps around 0°/360°
            seg = []
            for a, z in pts:
                xy = (X(z), Y(a))
                if seg and abs(xy[0] - seg[-1][0]) > half:
                    if len(seg) > 1:
                        cv.create_line(*[c for p_ in seg for c in p_], tags=tag, **kw)
                    seg = []
                seg.append(xy)
            if len(seg) > 1:
                cv.create_line(*[c for p_ in seg for c in p_], tags=tag, **kw)
        tr = risk.get("track") or []
        if not tr:
            return
        poly([(a, z) for _, a, z in tr], fill="#ff9800", width=2)
        for k, (t_, a, z) in enumerate((tr[0], tr[-1])):
            poly(sun_circle(a, z, sun_min), fill="#ff9800", dash=(6, 4) if k else (), width=2)
            x, y = X(z), Y(a)
            cv.create_oval(x - 8, y - 8, x + 8, y + 8, outline="#e65100", width=2,
                           fill="#ffd54f" if k == 0 else "", tags=tag)
            cv.create_text(x + 11, y - 11, anchor="w", text=time.strftime("%H:%M", time.localtime(t_)),
                           fill="#e65100", font=("Segoe UI", 8, "bold"), tags=tag)
        for i, t_, sep in risk.get("frames", []):
            alt, az, _ = frames[i]; cw = fw / max(math.cos(math.radians(alt)), 0.05) / 2
            kx = abs(X(az + 1) - X(az)) if abs(X(az + 1) - X(az)) < half else abs(X(az) - X(az - 1))
            cv.create_rectangle(X(az) - cw * kx, Y(alt + fh / 2), X(az) + cw * kx, Y(alt - fh / 2),
                                outline="#d50000", width=3, fill="#ffcdd2", tags=tag)
        for i, t_, sep in risk.get("paths", []):
            alt, az, _ = frames[i]; x, y = X(az), Y(alt)
            cv.create_line(x - 6, y - 6, x + 6, y + 6, fill="#d50000", width=3, tags=tag)
            cv.create_line(x - 6, y + 6, x + 6, y - 6, fill="#d50000", width=3, tags=tag)

    def sun_summary(risk, sun_min):
        tr = risk.get("track") or []
        if not tr:
            return ""
        t0_, t1_ = tr[0][0], tr[-1][0]
        day = "" if time.strftime("%Y%m%d", time.localtime(t0_)) == time.strftime("%Y%m%d") else \
            time.strftime("%d.%m.%Y ", time.localtime(t0_))
        if max(a for _, a, _ in tr) < -sun_min:
            return f"Sun {day}{time.strftime('%H:%M', time.localtime(t0_))}–{time.strftime('%H:%M', time.localtime(t1_))}: below the horizon – no risk"
        nf, np_ = len(risk.get("frames", [])), len(risk.get("paths", []))
        return (f"Sun {day}{time.strftime('%H:%M', time.localtime(t0_))}–{time.strftime('%H:%M', time.localtime(t1_))}: "
                + (f"{nf} frames within {sun_min:g}° (red)" if nf else "no frame too close")
                + (f", {np_} slews pass near it (×, detoured)" if np_ else ""))

    # every button / check box / radio button click goes to the debug log
    def _logged(cls, kind):
        class Logged(cls):
            def __init__(self, *a, **kw):
                cmd = kw.get("command")
                if callable(cmd):
                    label = kw.get("text") or kind

                    def wrapped(*ca, _cmd=cmd, _label=label):
                        dlog(f"UI {kind}: {_label!r}")
                        return _cmd(*ca)
                    kw["command"] = wrapped
                super().__init__(*a, **kw)
        return Logged
    ttk.Button = _hovered(_logged(ttk.Button, "button"), "button")
    ttk.Checkbutton = _hovered(_logged(ttk.Checkbutton, "check"), "check")
    ttk.Radiobutton = _hovered(_logged(ttk.Radiobutton, "radio"), "radio")
    tk.Button = _hovered(_logged(tk.Button, "button"), "button")

    root = tk.Tk()
    root.title(f"{APP_NAME}  v{APP_VERSION}")

    def tk_error(exc, val, tb):
        """Errors inside window callbacks: log them (the windowed exe has no console) and show them."""
        msg = "".join(traceback.format_exception(exc, val, tb))
        dlog(f"EXCEPTION in a window callback:\n{msg}")
        if selftest:                                      # CI: fail instead of showing a dialog
            print(msg, file=sys.stderr); os._exit(1)
        try:
            messagebox.showerror(APP_NAME, f"Unexpected error: {val}\n\nDetails are in the debug log.")
        except Exception:
            pass
    root.report_callback_exception = tk_error
    threading.excepthook = lambda a: dlog(f"EXCEPTION in thread {a.thread.name if a.thread else '?'}:\n"
                                          + "".join(traceback.format_exception(a.exc_type, a.exc_value, a.exc_traceback)))

    # Debug menu
    menubar = tk.Menu(root); root.config(menu=menubar)
    dbg_menu = tk.Menu(menubar, tearoff=False); menubar.add_cascade(label="Debug", menu=dbg_menu)
    dbg_var = tk.BooleanVar(value=_DBG["on"])

    def toggle_debug():
        p = debug_enable(dbg_var.get()); save_settings(debug_log=dbg_var.get())
        if dbg_var.get() and p:
            messagebox.showinfo(APP_NAME, f"Debug log on:\n{p}")

    def open_log_folder():
        d = os.path.join(app_data_dir(), "logs"); os.makedirs(d, exist_ok=True)
        try:
            os.startfile(d)                                   # Windows
        except Exception:
            messagebox.showinfo(APP_NAME, d)
    dbg_menu.add_checkbutton(label="Debug log (all actions)", variable=dbg_var, command=toggle_debug)
    logs_var = tk.BooleanVar(value=load_settings().get("show_logs", DEV_MODE))   # default: DEV_MODE
    dbg_menu.add_checkbutton(label="Show log windows", variable=logs_var)
    dbg_menu.add_separator()
    dbg_menu.add_command(label="Open log folder", command=open_log_folder)
    dbg_menu.add_command(label="Show current log file", command=lambda: messagebox.showinfo(
        APP_NAME, _DBG["path"] or "The debug log is off."))

    # watchdog: notes in the debug log whenever the window was blocked for more than a second
    ui_beat = {"t": time.time(), "doing": "idle"}

    def heartbeat():
        ui_beat["t"] = time.time(); root.after(200, heartbeat)

    def watchdog():
        warned = False
        while True:
            time.sleep(0.5)
            gap = time.time() - ui_beat["t"]
            if gap > 1.0 and not warned:
                dlog(f"UI BLOCKED for {gap:.1f} s so far – busy with: {ui_beat['doing']}"); warned = True
            elif gap <= 1.0 and warned:
                dlog("UI responsive again"); warned = False
    root.after(200, heartbeat)
    threading.Thread(target=watchdog, daemon=True, name="watchdog").start()
    root.geometry("1600x900")
    nb = ttk.Notebook(root); nb.pack(fill="both", expand=True)
    scroll_panels = []          # (holder frame, canvas, scrollbar) of the settings bars
    bars = []                   # (groups frame, fixed frame) of each settings bar

    def topbar(parent):
        """Settings bar across the top of a tab. Returns (bar, fixed): `bar` holds the setting groups
        side by side and scrolls sideways when they don't fit; `fixed` is a strip on the right
        (actions, log) that always stays visible."""
        outer = ttk.Frame(parent); outer.pack(side="top", fill="x", padx=6, pady=(6, 2))
        fixed = ttk.Frame(outer); fixed.pack(side="right", fill="y", padx=(6, 0))
        holder = ttk.Frame(outer); holder.pack(side="left", fill="both", expand=True)
        c = tk.Canvas(holder, highlightthickness=0, borderwidth=0)
        sb = ttk.Scrollbar(holder, orient="horizontal", command=c.xview)
        inner = ttk.Frame(c)
        c.create_window(0, 0, anchor="nw", window=inner)
        c.configure(xscrollcommand=sb.set)
        c.pack(side="top", fill="x", expand=True)

        def update(_=None):
            c.configure(scrollregion=(0, 0, inner.winfo_reqwidth(), inner.winfo_reqheight()),
                        height=max(inner.winfo_reqheight(), fixed.winfo_reqheight()))
            if inner.winfo_reqwidth() > c.winfo_width() > 1:
                if not sb.winfo_ismapped():
                    sb.pack(side="bottom", fill="x")
            else:
                if sb.winfo_ismapped():
                    sb.pack_forget()
                c.xview_moveto(0)
        inner.bind("<Configure>", update); c.bind("<Configure>", update); fixed.bind("<Configure>", update)
        scroll_panels.append((holder, c, sb)); bars.append((inner, fixed))
        ttk.Separator(parent, orient="horizontal").pack(side="top", fill="x", padx=6)
        return inner, fixed

    def group(parent, text):
        """One group of settings in a top bar."""
        g = ttk.LabelFrame(parent, text=text); g.pack(side="left", fill="y", padx=3, pady=2, anchor="n")
        return g

    def panel_wheel(e):
        """Mouse wheel over a settings bar scrolls it sideways (lists and text boxes keep their own)."""
        w = root.winfo_containing(e.x_root, e.y_root)
        if w is None or w.winfo_class() in ("Listbox", "Text", "TCombobox", "TSpinbox"):
            return
        for outer, c, sb in scroll_panels:
            if str(w).startswith(str(outer)) and sb.winfo_ismapped():
                step = -1 if (getattr(e, "num", 0) == 4 or getattr(e, "delta", 0) > 0) else 1
                c.xview_scroll(step * 3, "units"); return
    root.bind_all("<MouseWheel>", panel_wheel, add="+")
    root.bind_all("<Button-4>", panel_wheel, add="+"); root.bind_all("<Button-5>", panel_wheel, add="+")

    def labeled(parent, r, text, var, width=10, col=0, tip=None):
        ttk.Label(parent, text=text).grid(row=r, column=col, sticky="w", padx=4, pady=2)
        e = ttk.Entry(parent, textvariable=var, width=width); e.grid(row=r, column=col + 1, sticky="w", padx=4, pady=2)
        if tip:
            ttk.Label(parent, text=tip, foreground="#666").grid(row=r, column=col + 2, sticky="w")
        return e

    # ------------------------------------------------------------------ PLAN
    plan_tab = ttk.Frame(nb); nb.add(plan_tab, text="1. Plan sequence")
    cap_tab = ttk.Frame(nb); nb.add(cap_tab, text="2. Capture (NINA)")      # filled in below
    left, plan_fixed = topbar(plan_tab)
    right = ttk.Frame(plan_tab); right.pack(side="top", fill="both", expand=True, padx=6, pady=6)
    P = PlanParams()
    V = {}
    current_hz = tk.StringVar(value="")   # shared: Plan tab horizon file -> Build tab comparison

    def var(name, value, kind=tk.StringVar):
        V[name] = kind(value=value); return V[name]

    box = group(left, "Camera & telescope")
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

    box = group(left, "Site & area to image")
    labeled(box, 0, "Latitude (°)", var("lat", P.lat))
    V["mode"] = tk.StringVar(value="horizon")
    ttk.Radiobutton(box, text="Fixed altitude band", variable=V["mode"], value="band").grid(row=1, column=0, sticky="w")
    ttk.Radiobutton(box, text="Around current horizon", variable=V["mode"], value="horizon").grid(row=1, column=1, sticky="w")
    labeled(box, 2, "Band: alt min (°)", var("alt_min", P.alt_min))
    labeled(box, 3, "Band: alt max (°)", var("alt_max", P.alt_max))
    hzf = group(left, "Current horizon (.hrz)")
    hz_list = tk.Listbox(hzf, height=2, width=40); hz_list.pack(anchor="w", padx=4)
    Hover(hz_list, "Your current horizon file(s). The first one is also used on the Build tab.")

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
    hzb = ttk.Frame(hzf); hzb.pack(anchor="w", pady=2, padx=4)
    ttk.Button(hzb, text="Add current horizon…", command=add_hz).pack(side="left")
    ttk.Button(hzb, text="Clear", command=clear_hz).pack(side="left", padx=4)
    ttk.Label(hzf, text="The first file is your current horizon; it is also used on the Build tab.\n"
                        "Optional second file (e.g. an all-sky estimate): the strip spans both.",
              foreground="#666", justify="left").pack(anchor="w", padx=4)
    box = group(left, "Strip to image")
    labeled(box, 0, "Margin below (°)", var("margin_below", P.margin_below), width=6)
    labeled(box, 1, "Margin above (°)", var("margin_above", P.margin_above), width=6)
    labeled(box, 2, "Lowest altitude (°)", var("floor", P.floor), width=6)
    labeled(box, 3, "Overlap (0-1)", var("overlap", P.overlap), width=6)
    labeled(box, 0, "Extra frames at steep edges", var("edge_frames", P.edge_frames), width=6, col=3, tip="0 = off")
    labeled(box, 1, "Steep = more than (°/°)", var("edge_slope", P.edge_slope), width=6, col=3, tip="alt per az")
    labeled(box, 2, "Azimuth start (°)", var("az_start", P.az_start), width=6, col=3)
    labeled(box, 3, "Azimuth end (°)", var("az_end", P.az_end), width=6, col=3)

    box = group(left, "Rotator")
    V["rotator"] = tk.BooleanVar(value=True)
    ttk.Checkbutton(box, text="Keep frames level with rotator", variable=V["rotator"]).grid(row=0, column=0, columnspan=2, sticky="w")
    labeled(box, 1, "Offset (sky - mech, °)", var("rot_offset", P.rot_offset))
    labeled(box, 2, "Skip moves below (°)", var("rot_tolerance", P.rot_tolerance))
    labeled(box, 3, "Direction (+1/-1)", var("rot_sign", P.rot_sign))
    V["orientation"] = tk.StringVar(value="auto")
    ttk.Label(box, text="Orientation").grid(row=4, column=0, sticky="w", padx=4)
    ttk.Combobox(box, textvariable=V["orientation"], values=["auto", "landscape", "portrait"], width=10, state="readonly").grid(row=4, column=1, sticky="w")

    box = group(left, "Exposure")
    labeled(box, 0, "Exposure (s)", var("exposure", P.exposure))
    labeled(box, 1, "Gain (-1 = default)", var("gain", P.gain))
    labeled(box, 2, "Offset (-1 = default)", var("offset", P.offset))
    labeled(box, 3, "Binning", var("binning", P.binning))
    labeled(box, 4, "Seconds per frame", var("sec_per_frame", P.sec_per_frame))
    labeled(box, 5, "Seconds per rotator move", var("sec_per_rot", P.sec_per_rot))
    box = group(left, "Sequence")
    V["meridian_split"] = tk.BooleanVar(value=True); V["unpark"] = tk.BooleanVar(value=True)
    V["park"] = tk.BooleanVar(value=True); V["template"] = tk.BooleanVar(value=False); V["collapsed"] = tk.BooleanVar(value=True)
    ttk.Checkbutton(box, text="East half first, then west (one meridian flip)", variable=V["meridian_split"]).grid(row=0, column=0, columnspan=3, sticky="w")
    V["bottom_up"] = tk.BooleanVar(value=False)
    ttk.Checkbutton(box, text="Rows bottom to top in both halves", variable=V["bottom_up"]).grid(row=4, column=0, columnspan=3, sticky="w")
    ttk.Checkbutton(box, text="Unpark at start", variable=V["unpark"]).grid(row=1, column=0, sticky="w")
    ttk.Checkbutton(box, text="Park at end", variable=V["park"]).grid(row=1, column=1, sticky="w")
    ttk.Checkbutton(box, text="Save as template (else full sequence)", variable=V["template"]).grid(row=2, column=0, columnspan=3, sticky="w")
    ttk.Checkbutton(box, text="Rows collapsed (use Load sequence)", variable=V["collapsed"]).grid(row=3, column=0, columnspan=3, sticky="w")

    btns = ttk.LabelFrame(plan_fixed, text="Plan"); btns.pack(side="top", fill="both", expand=True, pady=2)
    plan_info = ttk.Label(right, text="", font=("Segoe UI", 10, "bold")); plan_info.pack(anchor="w")
    ttk.Label(right, text="Wheel: zoom · drag: move · double-click: whole sky", foreground="#666").pack(anchor="w")
    canvas = tk.Canvas(right, bg="white", height=300); canvas.pack(fill="both", expand=True)
    state = {"plan": None, "params": None, "view": (0.0, 360.0, 0.0, 90.0), "box": (50, 20, 900, 540)}

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
        for k in ("rotator", "meridian_split", "bottom_up", "unpark", "park", "template", "collapsed"):
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
        az0, az1, amin, amax = state["view"]
        state["box"] = (x0, y0, x1, y1)
        X = lambda a: x0 + ((a % 360) - az0) / (az1 - az0) * (x1 - x0)
        Y = lambda h: y1 - (h - amin) / (amax - amin) * (y1 - y0)

        def ticks(lo, hi, steps):
            st = next((t for t in steps if (hi - lo) / t <= 12), steps[-1])
            return [v for v in np.arange(math.ceil(lo / st) * st, hi + 1e-9, st)]
        XA = lambda a: x0 + (a - az0) / (az1 - az0) * (x1 - x0)       # no wrap, for tick marks
        a_ticks = ticks(az0, az1, (1, 2, 5, 10, 15, 30)); h_ticks = ticks(amin, amax, (1, 2, 5, 10))
        for a in a_ticks:
            canvas.create_line(XA(a), y0, XA(a), y1, fill="#ddd")
        for h in h_ticks:
            canvas.create_line(x0, Y(h), x1, Y(h), fill="#ddd")
        az = np.arange(0, 360, 0.25)
        lo_hi = [plan["lim"](a) for a in az]
        poly = [(X(a), Y(hi)) for a, (lo, hi) in zip(az, lo_hi)] + [(X(a), Y(lo)) for a, (lo, hi) in reversed(list(zip(az, lo_hi)))]
        canvas.create_polygon(*[c for pt in poly for c in pt], fill="#d6e6ff", outline="")
        kx = (x1 - x0) / (az1 - az0)
        for alt, a, rot in plan["frames"]:
            cw = plan["fw"] / max(math.cos(math.radians(alt)), 0.05) / 2; ch = plan["fh"] / 2
            for sh in (0, 360, -360):              # a frame across 0°/360° is drawn on both sides
                xc = XA(a % 360 + sh)
                if xc + cw * kx < x0 - 50 or xc - cw * kx > x1 + 50:
                    continue
                canvas.create_rectangle(xc - cw * kx, Y(alt + ch), xc + cw * kx, Y(alt - ch), outline="#4a7bd0")
        colors = ["#d00000", "#00a000"]
        for i, f in enumerate(p.horizon_files):
            pts = load_hrz(f)
            line = [c for a in az for c in (X(a), Y(hrz_alt(pts, a)))]
            canvas.create_line(*line, fill=colors[i % 2], width=2)
        canvas.create_line(x0, Y(p.floor), x1, Y(p.floor), fill="#e0b000", dash=(4, 3))
        if show_order.get():
            draw_order(canvas, plan, X, Y, flip=p.meridian_split)
        # the Sun during the planned session (start time + plan duration)
        sun_txt = ""
        try:
            t0_ = plan_start_time(); keep = float(sun_keep.get()); lat_ = float(V["lat"].get()); lon_ = float(sun_lon.get())
            key = (id(plan), round(t0_ / 60), keep, lat_, lon_, V["sec_per_frame"].get(), V["sec_per_rot"].get())
            if state.get("sunrisk", (None,))[0] != key:
                rk = sun_risk(plan["frames"], list(range(len(plan["frames"]))), t0_, float(V["sec_per_frame"].get()),
                              float(V["sec_per_rot"].get()), lat_, lon_, keep)
                state["sunrisk"] = (key, rk)
            rk = state["sunrisk"][1]
            draw_sun_overlay(canvas, X, Y, rk, plan["frames"], plan["fw"], plan["fh"], keep, x0, x1)
            sun_txt = sun_summary(rk, keep)
        except (ValueError, KeyError):
            sun_txt = "Sun: enter latitude, longitude, start date / time (e.g. 12.11.2026, 13:30) and distance to show it"
        # Tk does not clip: cover what lies outside the plot area, then draw the axis labels on top
        for r in ((0, 0, W, y0), (0, y1, W, H), (0, 0, x0, H), (x1, 0, W, H)):
            canvas.create_rectangle(*r, fill="white", outline="")
        canvas.create_rectangle(x0, y0, x1, y1, outline="#bbb")
        for a in a_ticks:
            canvas.create_text(XA(a), y1 + 14, text=f"{a:g}°")
        for h in h_ticks:
            canvas.create_text(x0 - 22, Y(h), text=f"{h:g}°")
        if sun_txt:
            canvas.create_text(x0, y0 - 10, anchor="w", font=("Segoe UI", 9), fill="#e65100", text=sun_txt)
        if show_order.get():
            canvas.create_text(x1, y0 - 10, anchor="ne", font=("Segoe UI", 9),
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

    ttk.Button(btns, text="Preview plan", command=preview).pack(side="top", fill="x", padx=4, pady=2)
    ttk.Button(btns, text="Save NINA sequence…", command=save).pack(side="top", fill="x", padx=4, pady=2)
    show_order = tk.BooleanVar(value=True)
    ttk.Checkbutton(btns, text="Show shooting order", variable=show_order,
                    command=lambda: draw_plan()).pack(side="top", anchor="w", padx=4, pady=2)
    sun_lon = tk.StringVar(value="23.62")
    sb = ttk.Frame(btns); sb.pack(side="top", fill="x", pady=2)
    ttk.Label(sb, text="Longitude (°E)").pack(side="left", padx=4)
    ttk.Entry(sb, textvariable=sun_lon, width=8).pack(side="left")
    ttk.Button(sb, text="Sun check", command=sun_check).pack(side="left", padx=4)
    sb2 = ttk.Frame(btns); sb2.pack(side="top", fill="x", pady=2)
    sun_start = tk.StringVar(value=""); sun_keep = tk.StringVar(value="30"); sun_date = tk.StringVar(value="")
    ttk.Label(sb2, text="Start date").pack(side="left", padx=4)
    ttk.Entry(sb2, textvariable=sun_date, width=10).pack(side="left")
    ttk.Label(sb2, text="Start time (HH:MM)").pack(side="left", padx=(8, 4))
    ttk.Entry(sb2, textvariable=sun_start, width=6).pack(side="left")
    ttk.Label(sb2, text="Sun distance (°)").pack(side="left", padx=(8, 4))
    ttk.Entry(sb2, textvariable=sun_keep, width=4).pack(side="left")
    for _v in (sun_start, sun_keep, sun_lon, sun_date):
        _v.trace_add("write", lambda *a: (state.pop("sunrisk", None), draw_plan()))

    def plan_start_time():
        """Planned start in local time: 'Start date' (blank = today) at 'Start time' (blank = now)."""
        return parse_start(sun_date.get(), sun_start.get())
    canvas.bind("<Configure>", lambda e: draw_plan())

    def plan_point(x, y):
        x0, y0, x1, y1 = state["box"]; az0, az1, h0, h1 = state["view"]
        return az0 + (x - x0) / (x1 - x0) * (az1 - az0), h0 + (y1 - y) / (y1 - y0) * (h1 - h0)

    def set_view(az0, az1, h0, h1):
        w = min(max(az1 - az0, 5.0), 360.0); h = min(max(h1 - h0, 3.0), 90.0)
        az0 = min(max(az0, 0.0), 360.0 - w); h0 = min(max(h0, 0.0), 90.0 - h)
        state["view"] = (az0, az0 + w, h0, h0 + h); draw_plan()

    def plan_wheel(e):
        a, h = plan_point(e.x, e.y)
        up = getattr(e, "num", 0) == 4 or getattr(e, "delta", 0) > 0
        f = 0.8 if up else 1.25
        az0, az1, h0, h1 = state["view"]
        set_view(a - (a - az0) * f, a + (az1 - a) * f, h - (h - h0) * f, h + (h1 - h) * f)
        return "break"

    def plan_drag_start(e):
        state["drag"] = (e.x, e.y, state["view"])

    def plan_drag(e):
        if "drag" not in state:
            return
        sx, sy, (az0, az1, h0, h1) = state["drag"]
        x0, y0, x1, y1 = state["box"]
        da = -(e.x - sx) / (x1 - x0) * (az1 - az0); dh = (e.y - sy) / (y1 - y0) * (h1 - h0)
        set_view(az0 + da, az1 + da, h0 + dh, h1 + dh)

    canvas.bind("<MouseWheel>", plan_wheel)
    canvas.bind("<Button-4>", plan_wheel); canvas.bind("<Button-5>", plan_wheel)
    canvas.bind("<ButtonPress-1>", plan_drag_start); canvas.bind("<B1-Motion>", plan_drag)
    canvas.bind("<Double-Button-1>", lambda e: set_view(0, 360, 0, 90))
    for k in ("sensor_w_px", "sensor_h_px", "pixel_um", "focal_mm"):
        V[k].trace_add("write", update_fov)
    update_fov()

    # ------------------------------------------------------------------ PROCESS
    proc = ttk.Frame(nb); nb.add(proc, text="3. Build horizon")
    pl, build_fixed = topbar(proc)
    pr = ttk.Frame(proc); pr.pack(side="top", fill="both", expand=True, padx=6, pady=6)
    box = group(pl, "FITS folders (searched recursively)")
    folders = tk.Listbox(box, height=5, width=44); folders.pack(fill="x", padx=4, pady=2)
    Hover(folders, "Folders (and capture lists) with the FITS frames to build the horizon from.")
    fb = ttk.Frame(box); fb.pack(fill="x")
    ttk.Button(fb, text="Add folder…", command=lambda: (lambda d: d and folders.insert("end", d))(
        filedialog.askdirectory())).pack(side="left", padx=3)
    ttk.Button(fb, text="Remove", command=lambda: [folders.delete(i) for i in reversed(folders.curselection())]).pack(side="left")
    ttk.Button(fb, text="Import from Capture…", command=lambda: import_capture()).pack(side="left", padx=3)
    box = group(pl, "Settings")
    pv = {}
    def pvar(n, v): pv[n] = tk.StringVar(value=v); return pv[n]
    labeled(box, 0, "Rotator offset guess (°)", pvar("K", "-58.04"), tip="sky - mech")
    pv["auto"] = tk.BooleanVar(value=True)
    ttk.Checkbutton(box, text="Calibrate rotation from overlaps", variable=pv["auto"]).grid(row=1, column=0, columnspan=3, sticky="w")
    pv["full"] = tk.BooleanVar(value=False)
    ttk.Checkbutton(box, text="Full 360° search (offset unknown)", variable=pv["full"]).grid(row=2, column=0, columnspan=3, sticky="w")
    labeled(box, 3, "Latitude (°, blank = FITS)", pvar("lat", ""))
    labeled(box, 4, "Mosaic resolution (°)", pvar("res", "0.03"))
    labeled(box, 0, "Solid: darker than", pvar("solid", "0.60"), width=6, col=3, tip="× sky")
    labeled(box, 1, "Branches: darker than", pvar("thin", "0.85"), width=6, col=3, tip="× sky")
    labeled(box, 2, "Bias level (blank = auto)", pvar("bias", ""), width=6, col=3)
    pvar("side", "0"); pvar("steep", "2")   # shown in the toolbar above the image
    pv["cache"] = tk.BooleanVar(value=True)
    ttk.Checkbutton(box, text="Cache reduced frames", variable=pv["cache"]).grid(row=3, column=3, columnspan=2, sticky="w", padx=4)
    ttk.Button(box, text="Clear cache", command=lambda: messagebox.showinfo(
        APP_NAME, f"Removed {clear_cache(default_cache_dir())} cached frames.")).grid(row=4, column=3, sticky="w", padx=4)
    box = group(pl, "Current horizon (compare / fill gaps)")
    old_var = current_hz   # filled automatically from the Plan tab
    ttk.Entry(box, textvariable=old_var, width=34,
              hover=TIPS["Current horizon (compare / fill gaps)"]).pack(side="top", anchor="w", padx=4, pady=2)
    ttk.Button(box, text="Browse…", command=lambda: old_var.set(filedialog.askopenfilename(
        filetypes=[("Horizon", "*.hrz *.txt"), ("All", "*.*")]) or old_var.get())).pack(side="top", anchor="w", padx=4)
    run_box = group(pl, "Build")
    fill_box = group(pl, "Fill-in frames")
    hz_box = group(pl, "Horizon file")
    bb = ttk.Frame(run_box); bb.pack(side="top", fill="x", padx=4, pady=2)
    pb = ttk.Progressbar(run_box, length=200, maximum=1.0); pb.pack(side="top", fill="x", padx=4, pady=4)
    status = ttk.Label(run_box, text="Idle", wraplength=220); status.pack(side="top", anchor="w", padx=4)
    log_box = ttk.LabelFrame(build_fixed, text="Log"); log_box.pack(side="left", fill="y", padx=(6, 0), pady=2)
    logbox = tk.Text(log_box, height=8, width=56, wrap="word"); logbox.pack(fill="both", expand=True, padx=2, pady=2)

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
    ttk.Label(pr, text="Wheel: zoom · Shift+wheel / middle-drag / scrollbars: move · left-drag: redraw the selected line"
              " · right-click: restore detected value · Ctrl+Z: undo",
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
        while len(E["bg"]) > 4:                         # keep memory in check while zooming
            E["bg"].pop(next(iter(E["bg"])))
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
    def wheel_zoom(e):
        if "mos" not in result or E.get("solid") is None:
            return
        x, y = cv.canvasx(e.x), cv.canvasy(e.y)
        a, h = from_xy(x, y)                                  # keep this point under the mouse
        up = getattr(e, "num", 0) == 4 or getattr(e, "delta", 0) > 0
        z = zoom(); nz = min(max(z * (1.25 if up else 0.8), 2.0), 24.0)
        if abs(nz - z) < 1e-6:
            return
        ed["zoom"].set(f"{nz:.3g}"); redraw_all()
        nx, ny = to_xy(a, h); z2, a0, top = geom()
        W = MARG * 2 + 360 * z2; H = MARG * 2 + (top - a0) * z2
        cv.xview_moveto(max(nx - e.x, 0) / W); cv.yview_moveto(max(ny - e.y, 0) / H)
        return "break"

    cv.bind("<MouseWheel>", wheel_zoom)
    cv.bind("<Button-4>", wheel_zoom); cv.bind("<Button-5>", wheel_zoom)
    cv.bind("<Shift-MouseWheel>", lambda e: cv.xview_scroll(int(-e.delta / 40), "units"))
    cv.bind("<ButtonPress-2>", lambda e: cv.scan_mark(e.x, e.y))
    cv.bind("<B2-Motion>", lambda e: cv.scan_dragto(e.x, e.y, gain=1))
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

    def lines_fingerprint():
        """Identifies the horizon as it would be saved now (lines with edits and side margin)."""
        if E.get("solid") is None:
            return None
        az, s_, b_ = final_lines()
        return hash((np.round(s_, 3).tobytes(), np.round(b_, 3).tobytes()))

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
        E["saved_fp"] = lines_fingerprint()
        messagebox.showinfo(APP_NAME, f"Saved the {which} horizon:\n{f}\n\n"
                                      "In NINA: Options → General → Astrometry → Horizon.")

    def import_capture():
        """Bring a capture to this tab: its file list (as a folder entry), the rotator offset and the
        current horizon; then offer to build. Source: the Capture tab's capture, or a saved one."""
        st_, pa_, src_ = [], [], None
        if C.get("frames"):
            st_, pa_, _ = progress_lists()
        n_taken = sum(1 for x, p_ in zip(st_, pa_) if x == "done" and p_ and os.path.isfile(p_))
        cap = None
        if n_taken:
            ans = messagebox.askyesnocancel(APP_NAME, f"Import the capture on the Capture tab ({n_taken} frames taken)?\n\n"
                                                      "No = choose a saved capture file instead.")
            if ans is None:
                return
            if ans:
                src_ = C.get("source") or "Capture tab"
                hz_ = list(getattr(C.get("params"), "horizon_files", None) or [])
        if src_ is None:
            d_ = os.path.join(app_data_dir(), "captures"); os.makedirs(d_, exist_ok=True)
            fn = filedialog.askopenfilename(title="Import a saved capture", initialdir=d_,
                                            filetypes=[("Horizon Mapper capture", "*.hmcap.json"), ("JSON", "*.json")])
            if not fn:
                return
            try:
                cap = load_capture_file(fn)
            except Exception as e:
                dlog_exc("import capture"); messagebox.showerror(APP_NAME, f"Could not read the capture: {e}"); return
            st_, pa_, src_, hz_ = cap["status"], cap["paths"], os.path.basename(fn), cap["horizon_files"]
        files = [p_ for x, p_ in zip(st_, pa_) if x == "done" and p_ and os.path.isfile(p_)]
        if not files:
            messagebox.showinfo(APP_NAME, "That capture has no frames with files yet."); return
        lst = write_capture_list(files)
        if lst not in folders.get(0, "end"):
            folders.insert("end", lst)
        notes = [f"{len(files)} frames from {src_} (file list {os.path.basename(lst)})"]
        try:
            pv["K"].set(V["rot_offset"].get()); notes.append(f"rotator offset {V['rot_offset'].get()}° from the Plan tab")
        except Exception:
            pass
        if hz_ and os.path.isfile(hz_[0]):
            current_hz.set(hz_[0]); notes.append(f"current horizon {os.path.basename(hz_[0])}")
        pv["lat"].set("")                                  # read from the FITS headers
        log("Imported from capture: " + "; ".join(notes) + ".")
        dlog("import capture: " + "; ".join(notes))
        if messagebox.askyesno(APP_NAME, "Imported:\n• " + "\n• ".join(notes) + "\n\nBuild the mosaic and horizon now?"):
            start()

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
        try:
            gn = max(float(fill_near.get()), 0.0)
        except ValueError:
            raise ValueError("Check the fill-in 'near gaps' value.")
        boxes = boxes + gap_boxes(result.get("mos"), az_, (s_, b_), gn, floor=p.floor, ceiling=p.ceiling)
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
            n_gap = sum(1 for b in boxes if b[4].startswith("Gap"))
            txt = (f"Fill-in: {len(boxes)} areas ({n_edge} steep edges, {n_gap} gaps), {len(plan['frames'])} frames, "
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

    ttk.Button(bb, text="Build mosaic & horizon", command=start).pack(side="top", fill="x", pady=2)
    ttk.Button(bb, text="Re-detect (new thresholds)", command=redetect).pack(side="top", fill="x", pady=2)
    fb = ttk.Frame(fill_box); fb.grid(row=3, column=0, padx=4, pady=2, sticky="w")
    ttk.Button(fb, text="Generate fill-in frames", command=generate_fill).pack(side="left")
    ttk.Button(fb, text="Save fill-in sequence…", command=save_fill).pack(side="left", padx=4)
    fill_info = ttk.Label(fill_box, text="", foreground="#555", wraplength=380)
    fill_info.grid(row=4, column=0, padx=4, sticky="w")
    fe = ttk.Frame(fill_box); fe.grid(row=0, column=0, padx=4, pady=2, sticky="w")
    # Start from the Plan tab's values and follow them whenever they change there.
    fill_edge = tk.StringVar(value=V["edge_frames"].get()); fill_slope = tk.StringVar(value=V["edge_slope"].get())
    V["edge_frames"].trace_add("write", lambda *a: fill_edge.set(V["edge_frames"].get()))
    V["edge_slope"].trace_add("write", lambda *a: fill_slope.set(V["edge_slope"].get()))
    ttk.Label(fe, text="Extra frames at steep edges ").pack(side="left")
    ttk.Spinbox(fe, textvariable=fill_edge, from_=0, to=5, increment=0.5, width=4).pack(side="left", padx=(2, 8))
    ttk.Label(fe, text="steep = more than").pack(side="left")
    ttk.Spinbox(fe, textvariable=fill_slope, from_=0.5, to=20, increment=0.5, width=4).pack(side="left", padx=2)
    ttk.Label(fe, text="°/°", foreground="#666").pack(side="left")
    ttk.Label(fill_box, text="0 = off; steep-edge values are copied from the Plan tab", foreground="#666").grid(
        row=1, column=0, padx=4, sticky="w")
    fg = ttk.Frame(fill_box); fg.grid(row=2, column=0, padx=4, pady=2, sticky="w")
    fill_near = tk.StringVar(value="1.0")
    ttk.Label(fg, text="Horizon line near gaps within (°)").pack(side="left")
    ttk.Spinbox(fg, textvariable=fill_near, from_=0, to=5, increment=0.5, width=4).pack(side="left", padx=2)
    ttk.Label(fg, text="(0 = off)", foreground="#666").pack(side="left")
    for _v in (fill_edge, fill_slope, fill_near):
        _v.trace_add("write", lambda *a: refresh_fill())
    save_which = tk.StringVar(value="solid")
    ttk.Combobox(hz_box, textvariable=save_which, values=["solid", "with branches"], width=14,
                 hover="Which line to save: 'solid' (buildings, trunks, dense trees) or 'with branches' (also thin branches).",
                 state="readonly").pack(side="top", anchor="w", padx=4, pady=2)
    ttk.Button(hz_box, text="Save horizon file…", command=save_outputs).pack(side="top", fill="x", padx=4, pady=2)

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


    # ------------------------------------------------------------------ CAPTURE (NINA Advanced API)
    # picture + settings on the left, the log on the right; drag the divider to resize the log
    cap_pane = ttk.PanedWindow(cap_tab, orient="horizontal"); cap_pane.pack(fill="both", expand=True)
    cmain = ttk.Frame(cap_pane)
    clog_box = ttk.LabelFrame(cap_pane, text="Log")
    clog_sb = ttk.Scrollbar(clog_box, orient="vertical")
    clog = tk.Text(clog_box, width=52, wrap="word", yscrollcommand=clog_sb.set)
    clog_sb.config(command=clog.yview)
    clog_sb.pack(side="right", fill="y"); clog.pack(side="left", fill="both", expand=True, padx=(2, 0), pady=2)
    cap_pane.add(cmain, weight=1); cap_pane.add(clog_box, weight=0)
    cbar, cfixed = topbar(cmain)
    C = {"runner": None, "frames": [], "source": "", "api": None, "polling": False, "q": queue.Queue(),
         "view": (0.0, 360.0, 0.0, 90.0), "box": (50, 20, 900, 400), "cur": None,
         "thumbs": {}, "paths": {}, "plan": None, "params": None, "photos": []}
    cv_ = {}

    def cvar(n, v, kind=tk.StringVar):
        cv_[n] = kind(value=v); return cv_[n]

    g = group(cbar, "NINA – Advanced API plugin")
    labeled(g, 0, "Host", cvar("host", "192.168.1.55" if DEV_MODE else "localhost"), width=12)   # DEV_MODE: test rig
    labeled(g, 1, "Port", cvar("port", "1888"), width=12)
    ttk.Button(g, text="Test connection", command=lambda: cap_connect()).grid(row=2, column=0, columnspan=2,
                                                                          sticky="w", padx=4, pady=2)
    conn_lbl = ttk.Label(g, text="Not connected", foreground="#666", wraplength=260, justify="left")
    conn_lbl.grid(row=3, column=0, columnspan=3, sticky="w", padx=4)
    g = group(cbar, "Frames to shoot")
    cvar("src", "plan")
    src_widgets = []
    src_widgets.append(ttk.Radiobutton(g, text="Plan from the Plan tab", variable=cv_["src"], value="plan"))
    src_widgets[-1].grid(row=0, column=0, sticky="w", padx=4)
    src_widgets.append(ttk.Radiobutton(g, text="Fill-in frames from the Build tab", variable=cv_["src"], value="fill"))
    src_widgets[-1].grid(row=1, column=0, sticky="w", padx=4)
    sf = ttk.Frame(g); sf.grid(row=2, column=0, sticky="w")
    src_widgets.append(ttk.Radiobutton(sf, text="Saved sequence file", variable=cv_["src"], value="file"))
    src_widgets[-1].pack(side="left", padx=4)
    src_widgets.append(ttk.Button(sf, text="Browse…", command=lambda: (lambda fn: fn and (cv_["file"].set(fn), cv_["src"].set("file"), source_changed()))(
        filedialog.askopenfilename(title="Plan or fill-in sequence saved for NINA",
                                   filetypes=[("NINA sequence", "*.json"), ("All", "*.*")]))))
    src_widgets[-1].pack(side="left")
    cvar("file", "")
    frames_lbl = ttk.Label(g, text="No frames loaded", foreground="#666", wraplength=280)
    frames_lbl.grid(row=3, column=0, sticky="w", padx=4)
    g = group(cbar, "Capture")
    labeled(g, 0, "Exposure (s)", cvar("exposure", V["exposure"].get()), width=7)
    labeled(g, 1, "Gain (-1 = NINA default)", cvar("gain", V["gain"].get()), width=7)
    labeled(g, 2, "Binning", cvar("binning", V["binning"].get()), width=7)
    labeled(g, 3, "Settle after moves (s)", cvar("settle", "1.0"), width=7)
    for k in ("exposure", "gain", "binning"):          # follow the Plan tab
        V[k].trace_add("write", lambda *a, k=k: cv_[k].set(V[k].get()))
    g = group(cbar, "Safety")
    labeled(g, 0, "Keep away from the Sun (°)", cvar("sun", "30"), width=6)
    labeled(g, 1, "Never below altitude (°)", cvar("floor", "10"), width=6)
    cvar("unpark", True, tk.BooleanVar); cvar("park", V["park"].get(), tk.BooleanVar)
    ttk.Checkbutton(g, text="Unpark at start ", variable=cv_["unpark"]).grid(row=2, column=0, sticky="w", padx=4)
    ttk.Checkbutton(g, text="Park when finished", variable=cv_["park"]).grid(row=3, column=0, sticky="w", padx=4)
    g = group(cbar, "Range")
    labeled(g, 0, "From frame", cvar("first", "1"), width=6)
    labeled(g, 1, "To frame (blank = last)", cvar("last", ""), width=6)
    cvar("reshoot", False, tk.BooleanVar)
    ttk.Checkbutton(g, text="Re-shoot frames already taken", variable=cv_["reshoot"]).grid(
        row=2, column=0, columnspan=3, sticky="w", padx=4)
    range_lbl = ttk.Label(g, text="", foreground="#666", wraplength=230); range_lbl.grid(row=3, column=0, columnspan=3, sticky="w", padx=4)
    g = group(cbar, "Files")
    labeled(g, 0, "Target name", cvar("target", "Horizon mapping"), width=18)
    ttk.Label(g, text="NINA image folder (found automatically)").grid(row=1, column=0, columnspan=3, sticky="w", padx=4)
    ttk.Entry(g, textvariable=cvar("folder", ""), width=30, hover=TIPS["NINA image folder (found automatically)"]).grid(row=2, column=0, columnspan=2, sticky="w", padx=4)
    cvar("thumbs", True, tk.BooleanVar); cvar("flip", False, tk.BooleanVar)

    def folder_changed(*a):
        C["folder_snap"] = cv_["folder"].get(); C["seen"] = {}
        r_ = C.get("runner")
        if r_ is not None:                      # the run's own look-up uses the new folder from now on
            r_.image_folder = C["folder_snap"]; r_._dir_hint = None
        dlog(f"image folder set to {C['folder_snap']!r} (exists: {os.path.isdir(C['folder_snap'])})")
    cv_["folder"].trace_add("write", folder_changed)
    ttk.Checkbutton(g, text="Show finished frames in the picture", variable=cv_["thumbs"],
                    command=lambda: draw_cap()).grid(row=3, column=0, columnspan=3, sticky="w", padx=4)
    ttk.Checkbutton(g, text="Turn previews 180°", variable=cv_["flip"],
                    command=lambda: cap_rethumb()).grid(row=4, column=0, columnspan=3, sticky="w", padx=4)
    ttk.Button(g, text="…", width=3, hover="Choose NINA's image folder.", command=lambda: cv_["folder"].set(
        filedialog.askdirectory() or cv_["folder"].get())).grid(row=2, column=2, sticky="w")

    # fixed on the right: run controls (Pause always visible) and the log
    runb = ttk.LabelFrame(cfixed, text="Run"); runb.pack(side="left", fill="y", pady=2)
    pause_btn = tk.Button(runb, text="⏸  PAUSE  (Esc)", bg="#c62828", fg="white", activebackground="#8e0000",
                          activeforeground="white", font=("Segoe UI", 12, "bold"), relief="raised", bd=3,
                          state="disabled", command=lambda: cap_pause())
    pause_btn.pack(side="top", fill="x", padx=4, pady=(4, 6), ipady=6)
    rb = ttk.Frame(runb); rb.pack(side="top", fill="x", padx=4)
    start_btn = ttk.Button(rb, text="▶ Start", command=lambda: cap_start()); start_btn.grid(row=0, column=0, sticky="ew", padx=1, pady=1)
    resume_btn = ttk.Button(rb, text="▶ Resume", state="disabled", command=lambda: cap_resume()); resume_btn.grid(row=0, column=1, sticky="ew", padx=1, pady=1)
    skip_btn = ttk.Button(rb, text="Skip frame", state="disabled", command=lambda: cap_skip()); skip_btn.grid(row=1, column=0, sticky="ew", padx=1, pady=1)
    stop_btn = ttk.Button(rb, text="■ End run", state="disabled", command=lambda: cap_stop()); stop_btn.grid(row=1, column=1, sticky="ew", padx=1, pady=1)

    class PctBar:
        """Progress bar with the percentage written in its centre."""
        def __init__(self, parent):
            self.c = tk.Canvas(parent, height=20, highlightthickness=1, highlightbackground="#a0a0a0",
                               bg="#f4f4f4", width=220)
            self.v, self.note = 0.0, ""
            self.c.bind("<Configure>", lambda e: self.draw())

        def set(self, v, note=""):
            self.v, self.note = min(max(float(v), 0.0), 1.0), note; self.draw()

        def __setitem__(self, k, v):
            if k == "value":
                self.set(v, self.note)

        def draw(self):
            try:
                w = self.c.winfo_width() or 220; h = self.c.winfo_height() or 20
                self.c.delete("all")
                if self.v > 0:
                    self.c.create_rectangle(0, 0, w * self.v, h, fill="#4caf50", outline="")
                self.c.create_text(w / 2, h / 2, text=f"{self.v * 100:.0f} %" + (f"  ({self.note})" if self.note else ""),
                                   font=("Segoe UI", 9, "bold"), fill="#111")
            except Exception:
                pass

    cap_pb = PctBar(runb); cap_pb.c.pack(side="top", fill="x", padx=4, pady=4)
    Hover(cap_pb.c, "Progress of the current run: frames done or skipped out of the frames this run shoots. "
                    "Without a run: frames taken out of all frames.")
    cap_status = ttk.Label(runb, text="Idle", wraplength=240, justify="left"); cap_status.pack(side="top", anchor="w", padx=4)

    slb = ttk.LabelFrame(cfixed, text="Save / Load Capture"); slb.pack(side="left", fill="y", padx=(6, 0), pady=2)
    save_cap_btn = ttk.Button(slb, text="Save current capture…", command=lambda: cap_save())
    save_cap_btn.pack(side="top", fill="x", padx=4, pady=(4, 2))
    load_cap_btn = ttk.Button(slb, text="Load saved capture…", command=lambda: cap_open())
    load_cap_btn.pack(side="top", fill="x", padx=4, pady=2)

    cright = ttk.Frame(cmain); cright.pack(side="top", fill="both", expand=True, padx=6, pady=6)
    dev_lbl = ttk.Label(cright, text="Devices: –", font=("Segoe UI", 10)); dev_lbl.pack(anchor="w")
    lg = ttk.Frame(cright); lg.pack(anchor="w", fill="x")
    ttk.Label(lg, text="Blue = to do · yellow = in progress · green = done (picture in its slot; light green = file not"
              " arrived yet) · orange = skipped"
              " · red = error · wheel: zoom · drag: move · double-click: whole sky", foreground="#666").pack(side="left")
    cvar("order", True, tk.BooleanVar)
    ttk.Checkbutton(lg, text="Show shooting order", variable=cv_["order"], command=lambda: draw_cap()).pack(side="left", padx=8)
    ccan = tk.Canvas(cright, bg="white", height=300); ccan.pack(fill="both", expand=True)

    def clog_add(m):
        clog.insert("end", time.strftime("%H:%M:%S ") + m + "\n"); clog.see("end")

    def api_now():
        return NinaApi(cv_["host"].get().strip() or "localhost", int(cv_["port"].get() or 1888), timeout=15)

    def cap_connect():
        """Check the plugin, read device states and NINA's image folder (in the background)."""
        api = api_now()

        def job():
            try:
                ver = api.version()
                prof = api.profile() or {}
                folder = ((prof.get("ImageFileSettings") or {}).get("FilePath")) or ""
                ftype = ((prof.get("ImageFileSettings") or {}).get("FileType")) or ""
                dlog(f"NINA image settings: {prof.get('ImageFileSettings')}")
                parts = []
                for name, fn in (("mount", api.mount), ("camera", api.camera), ("rotator", api.rotator)):
                    try:
                        parts.append(f"{name} {'connected' if fn().get('Connected') else 'NOT connected'}")
                    except NinaApiError:
                        parts.append(f"{name} ?")
                C["q"].put(("connected", (ver, folder, ftype, parts)))
            except NinaApiError as e:
                C["q"].put(("conn_error", str(e)))
        threading.Thread(target=job, daemon=True).start()
        conn_lbl.config(text="Connecting…")

    def poll_devices():
        """Device status line, every 2 s while connected."""
        def job():
            while C["polling"]:
                api = C["api"]
                try:
                    m = api.mount()
                    C["mount_pos"] = (float(m.get("Altitude", 0)), float(m.get("Azimuth", 0)))
                    if m.get("SiteLatitude") is not None:
                        C["site"] = (float(m["SiteLatitude"]), float(m["SiteLongitude"]))
                    txt = (f"Mount: alt {m.get('Altitude', 0):.1f}°  az {m.get('Azimuth', 0):.1f}°  "
                                            f"{str(m.get('SideOfPier', '')).replace('pier', 'pier ')}"
                                            f"{'  PARKED' if m.get('AtPark') else ''}{'  slewing' if m.get('Slewing') else ''}")
                except NinaApiError:
                    txt = "Mount: ?"
                try:
                    r = api.rotator()
                    txt += f"   ·   Rotator: {r.get('MechanicalPosition', 0):.1f}° mech{'  moving' if r.get('IsMoving') else ''}" \
                        if r.get("Connected") else "   ·   Rotator: not connected"
                except NinaApiError:
                    pass
                try:
                    c = api.camera(); txt += f"   ·   Camera: {c.get('CameraState', '?')}" if c.get("Connected") else "   ·   Camera: not connected"
                except NinaApiError:
                    txt += "   ·   Camera: ?"
                C["q"].put(("dev", txt)); time.sleep(2.0)
        if not C["polling"]:
            C["polling"] = True; threading.Thread(target=job, daemon=True).start()

    def cap_load(after=None, quiet=False):
        """Frames from the Plan tab, the Build tab fill-in, or a saved sequence file. The window
        settings are read here; the slow part (planning, reading the file) runs in the background."""
        r = C.get("runner")
        if r and r.state == "running":
            dlog("Load frames refused: run is running")
            messagebox.showinfo(APP_NAME, "A run is in progress. Pause it first to load other frames.")
            return False
        if r and r.state == "paused":
            if not messagebox.askyesno(APP_NAME, "Loading other frames ends the paused run.\n"
                                                 "Frames taken so far are kept (and added to the Build tab).\n\n"
                                                 "Load the new frames?"):
                return False
            dlog("Load frames while paused: ending the paused run first")
            threading.Thread(target=r.stop, daemon=True, name="stop-for-load").start()
        if C.get("loading"):
            if time.time() - C.get("loading_since", 0) < 60:
                dlog("Load frames ignored: frames are already being loaded")
                return False
            dlog("Load frames: previous load seems stuck (>60 s) – starting again")
        dlog(f"Load frames: source {cv_['src'].get()!r}, run state {r.state if r else None}")
        try:
            p = collect(require_area=False)
            src = cv_["src"].get(); fn = cv_["file"].get()
            p_full = collect() if src == "plan" else None
            if src == "fill" and not fill.get("plan"):
                raise ValueError("Generate fill-in frames on the Build tab first.")
            if src == "file" and not fn:
                raise ValueError("Choose the sequence file (Browse…).")
            fill_plan = dict(fill["plan"]) if src == "fill" else None
        except Exception as e:
            dlog(f"Load frames refused: {e}")
            frames_lbl.config(text=f"No frames: {e}")
            if not quiet:
                messagebox.showerror(APP_NAME, str(e))
            return False
        C["loading_since"] = time.time()
        C["loading"] = True; cap_buttons(C["runner"].state if C.get("runner") else "idle")
        frames_lbl.config(text="Loading frames…")
        dlog(f"loading frames from {src}")

        def job():
            try:
                st = {}
                if src == "fill":
                    plan = fill_plan; label = "fill-in frames"
                elif src == "file":
                    frames, st = read_nina_sequence(fn)
                    if not frames:
                        raise ValueError("No 'Slew to Alt/Az' + 'Take exposure' steps found in that file.")
                    w, h = p.fov
                    fw, fh = _effective_fov(w, h, p.rot_tolerance) if p.rotator else (w, h)
                    plan = dict(frames=frames, grid=rows_of(frames), fw=fw, fh=fh); label = os.path.basename(fn)
                else:
                    plan = make_plan(p_full); label = "Plan tab plan"
                if "grid" not in plan:
                    plan["grid"] = rows_of(plan["frames"])
                C["q"].put(("loaded", (plan, label, p, st, after)))
            except Exception as e:
                dlog_exc("loading frames")
                C["q"].put(("load_error", (str(e), quiet)))
        threading.Thread(target=job, daemon=True, name="load-frames").start()
        return True

    def cap_loaded(plan, label, p, st, after):
        try:
            _cap_loaded(plan, label, p, st, after)
        finally:
            C["loading"] = False
            r_ = C.get("runner")
            cap_buttons(r_.state if r_ else "idle")

    def _cap_loaded(plan, label, p, st, after):
        C["loading"] = False
        for k in ("exposure", "gain", "binning"):
            if st.get(k) is not None:
                cv_[k].set(st[k])
        C["frames"] = list(plan["frames"]); C["plan"] = plan; C["source"] = label; C["params"] = p
        C["status"] = ["pending"] * len(C["frames"]); C["runner"] = None; C["cur"] = None
        n_ = len(C["frames"])
        C["prog"] = dict(paths=[None] * n_, names=[None] * n_); C["from_capture"] = None
        update_range_label()
        C["thumbs"] = {}; C["paths"] = {}; C["photo_cache"] = {}; C["thumb_retry"] = {}; C["thumb_attempts"] = {}
        try:
            ex = float(cv_["exposure"].get())
        except ValueError:
            ex = 1.0
        mins = len(C["frames"]) * (ex + 6.0) / 60
        frames_lbl.config(text=f"{len(C['frames'])} frames ({label}), roughly {mins:.0f} min")
        dlog(f"loaded {len(C['frames'])} frames from {label}"); update_progress()
        cap_buttons("idle"); request_draw()
        if after:
            after()

    def progress_lists():
        """(status, paths, names) of the frames: from the run if there is one, else from the loaded
        capture or the freshly loaded frames."""
        r = C.get("runner")
        if r is not None and r.frames == C["frames"]:
            return list(r.status), list(r.paths), list(r.names)
        pr = C.get("prog") or {}
        n_ = len(C["frames"])
        return (list(C.get("status") or ["pending"] * n_), list(pr.get("paths") or [None] * n_),
                list(pr.get("names") or [None] * n_))

    def range_values():
        n_ = len(C["frames"])
        try:
            first = int(float(cv_["first"].get() or 1))
            last = int(float(cv_["last"].get())) if cv_["last"].get().strip() else n_
        except ValueError:
            raise ValueError("From frame / To frame must be frame numbers.")
        if not (1 <= first <= max(n_, 1)) or last < first:
            raise ValueError(f"The frame range must lie within 1–{n_}, with From ≤ To.")
        return first, min(last, n_)

    def capture_sun_risk(order, t0=None):
        """Sun risk for the frames in `order`, shot from now (or t0) at the Plan tab's time estimates."""
        lat_, lon_ = C.get("site") or (float(V["lat"].get()), float(sun_lon.get()))
        return sun_risk(C["frames"], order, time.time() if t0 is None else t0, float(V["sec_per_frame"].get()),
                        float(V["sec_per_rot"].get()), lat_, lon_, float(cv_["sun"].get()), start_pos=C.get("mount_pos"))

    def update_progress():
        """Percentage bar: progress of the current run, or frames taken overall when there is none."""
        r = C.get("runner"); st_ = C.get("status") or []
        if r is not None and r.order and r.frames == C["frames"]:
            n_ = len(r.order); d_ = sum(1 for i in r.order if st_[i] in ("done", "skipped"))
            cap_pb.set(d_ / n_, f"{d_} / {n_} this run")
        elif C["frames"]:
            n_ = len(C["frames"]); d_ = sum(1 for x in st_ if x == "done")
            cap_pb.set(d_ / n_, f"{d_} / {n_} taken")
        else:
            cap_pb.set(0.0, "")

    def update_range_label(*a):
        if not C["frames"]:
            range_lbl.config(text=""); return
        st_ = progress_lists()[0]
        try:
            first, last = range_values()
            todo = frames_to_shoot(st_, first, last, cv_["reshoot"].get())
            range_lbl.config(text=f"{len(todo)} of frames {first}–{last} to shoot "
                                  f"({sum(1 for x in st_ if x == 'done')} of {len(st_)} taken)")
        except ValueError as e:
            range_lbl.config(text=str(e))
    for _k in ("first", "last", "reshoot"):
        cv_[_k].trace_add("write", update_range_label)

    def capture_settings():
        keys = ("host", "port", "exposure", "gain", "binning", "settle", "sun", "floor", "unpark", "park",
                "target", "folder", "first", "last", "reshoot", "src", "file")
        return {k: cv_[k].get() for k in keys}

    def write_capture(path):
        st_, pa_, nm_ = progress_lists()
        p_ = C.get("params")
        hz_files = list(getattr(p_, "horizon_files", None) or []) if p_ is not None else []
        return save_capture_file(path, C["frames"], st_, pa_, nm_, C["plan"]["fw"], C["plan"]["fh"],
                                 source=C.get("source", ""), settings=capture_settings(), horizon_files=hz_files)

    def autosave():
        """After every frame: keep the progress in last_capture.hmcap.json (for resuming after a crash)."""
        if not C["frames"] or not C.get("plan"):
            return
        try:
            d = os.path.join(app_data_dir(), "captures"); os.makedirs(d, exist_ok=True)
            write_capture(os.path.join(d, "last_capture.hmcap.json"))
        except Exception as e:
            dlog(f"autosave failed: {e}")

    def cap_save():
        if not C["frames"] or not C.get("plan"):
            messagebox.showinfo(APP_NAME, "There is no capture to save yet."); return
        d = os.path.join(app_data_dir(), "captures"); os.makedirs(d, exist_ok=True)
        fn = filedialog.asksaveasfilename(
            title="Save current capture", initialdir=d, defaultextension=".hmcap.json",
            initialfile=time.strftime(f"{cv_['target'].get() or 'capture'} %Y-%m-%d %H%M.hmcap.json"),
            filetypes=[("Horizon Mapper capture", "*.hmcap.json"), ("JSON", "*.json")])
        if not fn:
            return
        try:
            write_capture(fn)
        except Exception as e:
            dlog_exc("saving capture"); messagebox.showerror(APP_NAME, f"Could not save: {e}"); return
        st_ = progress_lists()[0]
        clog_add(f"Capture saved ({sum(1 for x in st_ if x == 'done')} of {len(st_)} frames taken): {fn}")
        dlog(f"capture saved to {fn}")

    def cap_open():
        r = C.get("runner")
        if r and r.state == "running":
            messagebox.showinfo(APP_NAME, "A run is in progress. Pause it first."); return
        if r and r.state == "paused":
            if not messagebox.askyesno(APP_NAME, "Loading a saved capture ends the paused run.\n"
                                                 "Its frames are kept – save it first if you need it.\n\nContinue?"):
                return
            threading.Thread(target=r.stop, daemon=True, name="stop-for-open").start()
        d = os.path.join(app_data_dir(), "captures"); os.makedirs(d, exist_ok=True)
        fn = filedialog.askopenfilename(title="Load saved capture", initialdir=d,
                                        filetypes=[("Horizon Mapper capture", "*.hmcap.json"), ("JSON", "*.json")])
        if not fn:
            return
        try:
            cap = load_capture_file(fn)
        except Exception as e:
            dlog_exc("loading capture"); messagebox.showerror(APP_NAME, f"Could not load: {e}"); return
        dlog(f"capture loaded from {fn}: {len(cap['frames'])} frames, "
             f"{sum(1 for x in cap['status'] if x == 'done')} taken, {len(cap['missing'])} files missing")
        C["opening"] = True                      # don't let restored settings trigger a source reload
        try:
            for k, v in (cap["settings"] or {}).items():
                if k in cv_ and k not in ("src", "file"):
                    try:
                        cv_[k].set(v)
                    except Exception:
                        pass
        finally:
            C["opening"] = False
        import types
        C["frames"] = cap["frames"]
        C["plan"] = dict(frames=cap["frames"], grid=rows_of(cap["frames"]), fw=cap["fw"], fh=cap["fh"])
        C["params"] = types.SimpleNamespace(mode="horizon", horizon_files=cap["horizon_files"])
        C["source"] = f"saved capture {os.path.basename(fn)}"; C["from_capture"] = fn
        C["status"] = list(cap["status"]); C["prog"] = dict(paths=cap["paths"], names=cap["names"])
        C["runner"] = None; C["cur"] = None
        C["thumbs"] = {}; C["paths"] = {}; C["photo_cache"] = {}; C["thumb_retry"] = {}; C["thumb_attempts"] = {}
        for i, (st_, pth) in enumerate(zip(cap["status"], cap["paths"])):
            if st_ == "done" and pth:
                C["paths"][i] = pth; cap_thumb(i, pth)
        taken = sum(1 for x in cap["status"] if x == "done")
        frames_lbl.config(text=f"{len(cap['frames'])} frames from the saved capture ({taken} taken)")
        clog_add(f"Loaded saved capture {os.path.basename(fn)} (saved {cap['saved']}): {taken} of "
                 f"{len(cap['frames'])} frames taken" + (f"; {len(cap['missing'])} files missing – those frames "
                 f"will be shot again: {', '.join(str(i + 1) for i in cap['missing'][:20])}"
                 f"{' …' if len(cap['missing']) > 20 else ''}" if cap["missing"] else ""))
        update_range_label(); update_progress(); cap_buttons("idle"); request_draw()

    def source_changed(*a):
        """Choosing a frame source shows its frames at once (asks first if that replaces progress)."""
        if C.get("opening") or C.get("loading"):
            return
        r = C.get("runner")
        if r and r.state in ("running", "paused"):
            return
        st_ = progress_lists()[0] if C["frames"] else []
        taken = sum(1 for x in st_ if x == "done")
        if taken and not messagebox.askyesno(
                APP_NAME, f"Show the frames from the new source?\n\nThe current frames and their progress "
                          f"({taken} taken) are replaced. They are also kept in the automatic save "
                          f"'last_capture.hmcap.json' – or use 'Save current capture…' first."):
            dlog("source change cancelled by user"); return
        cap_load(quiet=True)
    cv_["src"].trace_add("write", lambda *a: root.after_idle(source_changed))

    def thumb_K():
        """Rotator offset for the previews: the Build tab's calibrated value if there is one."""
        K = result.get("K")
        if K is None:
            try:
                K = float(V["rot_offset"].get()) + 180.0          # measured: K = offset + 180
            except ValueError:
                K = 0.0
        return K + (180.0 if cv_["flip"].get() else 0.0)

    def cap_thumb(i, path):
        """Project a finished frame into its slot (in the background)."""
        if not path or not os.path.isfile(path) or not C.get("plan"):
            return
        fw, fh = C["plan"]["fw"], C["plan"]["fh"]
        try:
            lat = float(V["lat"].get())
        except ValueError:
            lat = 60.0
        K = thumb_K(); q_ = C["q"]

        def job():
            if not fits_complete(path):                    # still being written: try again shortly
                dlog(f"preview frame {i + 1}: file not complete yet ({os.path.getsize(path) if os.path.exists(path) else 0} bytes)")
                q_.put(("thumb_retry", (i, path, "file not complete yet"))); return
            try:
                th = frame_thumbnail(path, K, lat, fw / 2, fh / 2, cache_dir=default_cache_dir())
                q_.put(("thumb", (i, th, thumbnail_rgba(th["img"]))))
            except Exception as e:
                dlog(f"preview frame {i + 1} failed: {type(e).__name__}: {e}")
                q_.put(("thumb_retry", (i, path, f"{type(e).__name__}: {e}")))
        threading.Thread(target=job, daemon=True, name=f"preview-{i + 1}").start()

    def cap_rethumb():
        C["thumbs"] = {}; C["photo_cache"] = {}; C["thumb_retry"] = {}; C["thumb_attempts"] = {}
        for i, pth in list(C.get("paths", {}).items()):
            cap_thumb(i, pth)
        draw_cap()

    def cap_buttons(st):
        running, paused = st == "running", st == "paused"
        pause_btn.config(state="normal" if running else "disabled")
        resume_btn.config(state="normal" if paused else "disabled")
        skip_btn.config(state="normal" if paused else "disabled")
        stop_btn.config(state="normal" if (running or paused) else "disabled")
        start_btn.config(state="disabled" if (running or paused) else "normal")
        busy = running or C.get("loading")
        load_cap_btn.config(state="disabled" if busy else "normal")
        for w_ in src_widgets:
            try:
                w_.config(state="disabled" if (running or paused) else "normal")
            except Exception:
                pass

    def cap_start():
        if C["runner"] and C["runner"].state in ("running", "paused"):
            return
        if C.get("loading"):
            return
        st0 = progress_lists()[0] if C["frames"] else []
        fresh = not any(x == "done" for x in st0) and not C.get("from_capture")
        if not C["frames"] or (fresh and not C.get("started_fresh_load")):
            # nothing taken yet: take the frames from the source as it is now (e.g. Plan tab changes)
            C["started_fresh_load"] = True
            if not cap_load(after=cap_start):
                C["started_fresh_load"] = False
            return
        C["started_fresh_load"] = False
        try:
            ex = float(cv_["exposure"].get()); gain = int(float(cv_["gain"].get()))
            binning = int(float(cv_["binning"].get())); settle = float(cv_["settle"].get())
            sun = float(cv_["sun"].get()); floor_ = float(cv_["floor"].get())
            first, last = range_values()
        except ValueError as e:
            messagebox.showerror(APP_NAME, f"Check the capture, safety and range values. {e}"); return
        status, paths, names = progress_lists()
        order = frames_to_shoot(status, first, last, cv_["reshoot"].get())
        if not order:
            messagebox.showinfo(APP_NAME, f"Frames {first}–{last} are all taken already.\n"
                                          "Tick 'Re-shoot frames already taken' to shoot them again."); return
        try:
            rk = capture_sun_risk(order)
        except Exception:
            dlog_exc("sun risk at start"); rk = {"frames": [], "paths": []}
        if rk["frames"] or rk["paths"]:
            fmt = lambda t_: time.strftime("%H:%M", time.localtime(t_))
            fr_ = rk["frames"]; pa_ = rk["paths"]
            lines = []
            if fr_:
                azs = [C["frames"][i][1] for i, _, _ in fr_]
                lines.append(f"• {len(fr_)} frames will be closer than {sun:g}° to the Sun when their turn comes "
                             f"(azimuth {min(azs):.0f}–{max(azs):.0f}°, around {fmt(fr_[0][1])}–{fmt(fr_[-1][1])}). "
                             "They are skipped automatically – never slewed to.")
            if pa_:
                lines.append(f"• {len(pa_)} slews would pass near the Sun on the way; they go via a safe "
                             "intermediate position (or are skipped if there is none).")
            dlog(f"sun risk at start: {len(fr_)} frames, {len(pa_)} slews")
            ans = messagebox.askyesnocancel(
                APP_NAME, "Sun check for this run:\n\n" + "\n".join(lines) +
                "\n\nYes = shoot as planned (the protections above apply)\n"
                "No = leave those frames out of this run (they stay 'to do' for later)\nCancel = don't start",
                icon="warning")
            if ans is None:
                return
            if ans is False:
                bad = {i for i, _, _ in fr_}
                order = [i for i in order if i not in bad]
                if not order:
                    messagebox.showinfo(APP_NAME, "No frames left to shoot in this range."); return
        taken_in = sum(1 for i in range(first - 1, last) if status[i] == "done")
        if not messagebox.askokcancel(APP_NAME,
                f"Shoot {len(order)} frames through NINA?\n"
                f"Range: frames {first}–{last} of {len(C['frames'])} ({C['source']})"
                + (f"; {taken_in} already taken are {'shot again' if cv_['reshoot'].get() else 'skipped'}" if taken_in else "")
                + ".\n\nThe mount, rotator and camera will move on their own.\n"
                "Stay nearby – PAUSE (or Esc) stops every movement at once.", icon="warning"):
            return
        for i in order:                                   # frames to shoot start from scratch
            status[i] = "pending"
        q_ = C["q"]
        runner = CaptureRunner(api_now(), C["frames"], lambda kind, **d: q_.put(("run", (kind, d, runner))),
                               exposure=ex, gain=(gain if gain >= 0 else None), binning=binning,
                               target=cv_["target"].get() or "Horizon mapping", settle=settle, sun_min=sun,
                               floor=floor_, unpark=cv_["unpark"].get(), park=cv_["park"].get(),
                               image_folder=cv_["folder"].get(), order=order, status=status, paths=paths,
                               names=names)
        C["runner"] = runner; C["status"] = list(runner.status); update_progress()
        clog_add(f"Start: {len(order)} frames ({first}–{last}), exposure {ex} s, sun ≥ {sun}°, altitude ≥ {floor_}°.")
        dlog(f"start: order {order[:10]}{' …' if len(order) > 10 else ''}")
        if not C["polling"]:
            C["api"] = api_now(); poll_devices()
        runner.start(); cap_buttons("running"); scan_folder(); request_draw()

    def cap_pause(reason="Paused by user"):
        r = C["runner"]
        if r and r.state == "running":
            threading.Thread(target=r.pause, args=(reason,), daemon=True).start()
            cap_status.config(text="Pausing – stopping mount, rotator and camera…")

    def cap_resume():
        r = C["runner"]
        if r and r.state == "paused":
            r.resume(); cap_buttons("running")

    def cap_skip():
        r = C["runner"]
        if r and r.state == "paused":
            r.skip(); clog_add(f"Frame {r.index + 1} will be skipped – press Resume to continue.")

    def cap_stop():
        r = C["runner"]
        if r and r.state in ("running", "paused"):
            if messagebox.askyesno(APP_NAME, "End the run? Frames taken so far are kept, and Start later "
                                             "shoots what is still missing.\n(PAUSE holds the run instead.)"):
                threading.Thread(target=r.stop, daemon=True).start()

    root.bind_all("<Escape>", lambda e: (dlog("UI key: Esc"), cap_pause("Esc pressed")), add="+")

    def waiting_frames(r=None):
        r = r or C["runner"]
        if not r:
            return []
        return [i for i, s_ in enumerate(r.status) if s_ == "done" and not (r.paths[i] and os.path.isfile(r.paths[i]))]

    def update_capture_list(r=None):
        """(Re)write the run's file list and make sure the Build tab has it."""
        r = r or C["runner"]
        paths = [p for p in r.paths if p and os.path.isabs(p) and os.path.isfile(p)]
        if not paths:
            return None
        lst = getattr(r, "list_file", None)
        if lst:
            with open(lst, "w", encoding="utf-8") as fh:
                fh.write("# Horizon Mapper capture list – one FITS file per line\n")
                fh.writelines(p + "\n" for p in paths)
        else:
            lst = r.list_file = write_capture_list(paths)
        if lst not in folders.get(0, "end"):
            folders.insert("end", lst)
        return lst

    def scan_folder():
        """Every few seconds: look in NINA's image folder for files of frames still without one
        (NINA saves in the background; files can arrive late). Keeps going for 10 min after a run."""
        if C.get("scanning"):
            return
        C["scanning"] = True; C["seen"] = {}
        C["folder_snap"] = cv_["folder"].get()          # kept up to date by a trace on the field

        def job():
            idle_until = None
            while C.get("scanning"):
                r = C["runner"]
                if r is None:
                    break
                wait = waiting_frames()
                active = r.state in ("running", "paused")
                if not active and not wait:
                    break
                if not active:
                    idle_until = idle_until or time.time() + 600
                    if time.time() > idle_until:
                        C["q"].put(("log", f"{len(wait)} frames still have no file after 10 min – stopped looking."))
                        break
                if wait:
                    try:
                        hits = match_late_files(C["folder_snap"], r.frames, wait, r.names,
                                                since=r.t_start - 120, seen=C["seen"])
                    except Exception as e:
                        hits = {}; C["q"].put(("log", f"Folder scan failed: {e}"))
                    fo = C["folder_snap"]
                    dlog(f"folder scan in {fo!r} (exists: {os.path.isdir(fo) if fo else False}): "
                         f"{len(wait)} frames waiting {[i + 1 for i in wait][:10]}, names {[r.names[i] for i in wait][:3]}, "
                         f"{len(C['seen'])} files looked at, {len(hits)} found")
                    for i, pth in hits.items():
                        C["q"].put(("late_file", (i, pth)))
                time.sleep(4.0)
            C["scanning"] = False
        threading.Thread(target=job, daemon=True).start()

    def cap_finished(ok, r=None):
        r = r or C["runner"]
        if r is C["runner"]:
            autosave(); update_range_label()
        done = sum(1 for s_ in r.status if s_ == "done")
        msg = f"Run {'finished' if ok else 'stopped'}: {done} frames taken"
        lst = update_capture_list(r)
        if lst:
            msg += f"; the list of files was added to the Build tab:\n{lst}"
        wait = waiting_frames(r)
        current = r is C["runner"]
        if wait and current:
            msg += (f"\n\n{len(wait)} files have not appeared in NINA's image folder yet. The app keeps looking "
                    "for 10 minutes and adds them to the picture and the list as they arrive.")
        elif wait:
            msg += f"\n\n{len(wait)} files of that run had not appeared in NINA's image folder yet."
        clog_add(msg.replace("\n", " "))
        if current:
            cap_buttons(r.state)
        messagebox.showinfo(APP_NAME, msg)

    def request_draw():
        C["dirty"] = True

    def cap_poll():
        """Handle events from the worker threads – for at most ~80 ms per tick so the window stays
        responsive (PAUSE must always work), and redraw the picture at most ~4 times a second."""
        t0 = time.time(); n = 0
        try:
            while time.time() - t0 < 0.08:
                kind, val = C["q"].get_nowait(); n += 1
                ui_beat["doing"] = f"capture event {kind}"
                if kind == "connected":
                    ver, folder, ftype, parts = val
                    conn_lbl.config(text=f"Advanced API {ver}: " + ", ".join(parts))
                    if folder and not cv_["folder"].get():
                        cv_["folder"].set(folder)
                    if ftype and str(ftype).upper() not in ("FITS", "0"):
                        clog_add(f"Note: NINA saves {ftype} files – set NINA to FITS for the Build tab.")
                    clog_add("Connected to NINA: " + ", ".join(parts))
                    C["api"] = api_now(); poll_devices()
                elif kind == "conn_error":
                    conn_lbl.config(text="Not connected – is NINA running with the Advanced API plugin enabled?")
                    clog_add(val)
                elif kind == "dev":
                    dev_lbl.config(text=val)
                elif kind == "thumb_retry":
                    i, pth, why = val
                    att = C.setdefault("thumb_attempts", {}).get(i, 0) + 1
                    C["thumb_attempts"][i] = att
                    if att <= 60:                          # every 3 s for up to 3 minutes
                        C.setdefault("thumb_retry", {})[i] = (pth, time.time() + 3.0)
                    else:
                        clog_add(f"Frame {i + 1}: no preview after 3 minutes ({why}).")
                elif kind == "thumb":
                    i, th, im = val
                    C.get("thumb_attempts", {}).pop(i, None)
                    if i < len(C["frames"]):
                        C.setdefault("thumbs", {})[i] = (th, im)
                        for key in [k_ for k_ in C.get("photo_cache", {}) if k_[0] == i]:
                            del C["photo_cache"][key]
                        request_draw()
                elif kind == "log":
                    clog_add(val)
                elif kind == "loaded":
                    cap_loaded(*val)
                elif kind == "load_error":
                    msg_, quiet_ = val
                    C["loading"] = False; cap_buttons("idle"); frames_lbl.config(text=f"No frames: {msg_}")
                    if not quiet_:
                        messagebox.showerror(APP_NAME, msg_)
                elif kind == "late_file":
                    i, pth = val
                    r = C["runner"]
                    if r and i < len(r.paths) and not (r.paths[i] and os.path.isfile(r.paths[i])):
                        r.paths[i] = pth; C.setdefault("paths", {})[i] = pth; C["autosave_due"] = True
                        cap_thumb(i, pth)
                        if r.state not in ("running", "paused"):
                            update_capture_list()
                        n_wait = len(waiting_frames())
                        clog_add(f"Frame {i + 1}: file arrived – {os.path.basename(pth)}"
                                 + (f" ({n_wait} still waiting)" if n_wait else ""))
                elif kind == "run":
                    k, d, r = val
                    if r is not C["runner"]:                  # a run that was replaced (Load frames while paused)
                        if k == "log":
                            clog_add(d["msg"])
                        elif k == "done":
                            cap_finished(d.get("ok", False), r)
                        continue
                    if k == "log":
                        clog_add(d["msg"])
                        if d["msg"].startswith("PAUSED"):
                            cap_status.config(text=d["msg"]); cap_buttons("paused")
                            try:
                                root.bell()
                            except Exception:
                                pass
                    elif k == "state":
                        cap_buttons(d["state"])
                        if d["state"] == "running":
                            cap_status.config(text="Running")
                    elif k == "frame":
                        i = d["index"]; C["status"][i] = d["status"]; C["cur"] = i
                        if d["status"] in ("done", "skipped", "error"):
                            C["autosave_due"] = True
                        if d["status"] == "done" and d.get("path"):
                            C.setdefault("paths", {})[i] = d["path"]; cap_thumb(i, d["path"])
                        elif d["status"] == "done":
                            scan_folder()                  # file not there yet: look for it later
                        nf = len(C["frames"]); update_progress()
                        if d["status"] in ("slewing", "rotating", "exposing"):
                            alt, az, _ = C["frames"][i]
                            cap_status.config(text=f"Frame {i + 1}/{nf}: {d['status']} – alt {alt:.1f}°, az {az:.1f}°")
                        request_draw()
                    elif k == "done":
                        cap_finished(d.get("ok", False))
        except queue.Empty:
            pass
        except Exception:
            dlog_exc("cap_poll")
        if C.get("autosave_due") and time.time() - C.get("last_autosave", 0) > 2.0:
            C["autosave_due"] = False; C["last_autosave"] = time.time(); autosave(); update_range_label()
        for i, (pth, due) in list(C.get("thumb_retry", {}).items()):
            if time.time() >= due:
                del C["thumb_retry"][i]; cap_thumb(i, pth)
        if n:
            dlog(f"UI handled {n} capture events in {time.time() - t0:.3f} s" + (" – more queued" if not C["q"].empty() else ""))
        if C.get("dirty") and time.time() - C.get("last_draw", 0) > 0.25:
            C["dirty"] = False; ui_beat["doing"] = "drawing the capture picture"
            t1 = time.time(); draw_cap(); C["last_draw"] = time.time()
            dlog(f"capture picture drawn in {C['last_draw'] - t1:.3f} s ({len(C.get('thumbs', {}))} previews)")
        ui_beat["doing"] = "idle"
        root.after(50 if not C["q"].empty() else 150, cap_poll)
    cap_poll()

    STATUS_COL = {"pending": "#4a7bd0", "slewing": "#f9a825", "rotating": "#f9a825", "exposing": "#f9a825",
                  "done": "#2e7d32", "skipped": "#ef6c00", "error": "#c62828"}

    def draw_cap():
        """Like the Plan preview (strip, horizon, frames, shooting order) plus each frame's state;
        finished frames are shown in their slots."""
        ccan.delete("all"); C["photos"] = []
        if len(C.setdefault("photo_cache", {})) > 3000:      # zoom levels pile up: start afresh
            C["photo_cache"] = {}
        W = ccan.winfo_width() or 900; H = ccan.winfo_height() or 400
        x0, y0, x1, y1 = 50, 20, W - 15, H - 30
        az0, az1, h0, h1 = C["view"]; C["box"] = (x0, y0, x1, y1)
        XA = lambda a: x0 + (a - az0) / (az1 - az0) * (x1 - x0)
        X = lambda a: XA(a % 360)
        Y = lambda h: y1 - (h - h0) / (h1 - h0) * (y1 - y0)
        kx = (x1 - x0) / (az1 - az0)
        st_ = next((t for t in (1, 2, 5, 10, 15, 30) if (az1 - az0) / t <= 12), 30)
        at = [v for v in np.arange(math.ceil(az0 / st_) * st_, az1 + 1e-9, st_)]
        sh_ = next((t for t in (1, 2, 5, 10) if (h1 - h0) / t <= 12), 10)
        ht = [v for v in np.arange(math.ceil(h0 / sh_) * sh_, h1 + 1e-9, sh_)]
        for a in at:
            ccan.create_line(XA(a), y0, XA(a), y1, fill="#e3e3e3")
        for h in ht:
            ccan.create_line(x0, Y(h), x1, Y(h), fill="#e3e3e3")
        plan = C.get("plan"); p = C.get("params")
        if plan and plan.get("lim"):                               # strip to image, as on the Plan tab
            azs = np.arange(0, 360.01, 0.25); lh = [plan["lim"](a % 360) for a in azs]
            poly = [c for a, (lo, hi) in zip(azs, lh) for c in (XA(a), Y(hi))] + \
                   [c for a, (lo, hi) in reversed(list(zip(azs, lh))) for c in (XA(a), Y(lo))]
            ccan.create_polygon(*poly, fill="#d6e6ff", outline="")
        try:
            floor_ = float(cv_["floor"].get())
            ccan.create_line(x0, Y(floor_), x1, Y(floor_), fill="#e0b000", dash=(4, 3))
        except ValueError:
            pass
        if plan and C["frames"]:
            fw, fh = plan["fw"], plan["fh"]
            for i, (alt, az, _) in enumerate(C["frames"]):
                cw = fw / max(math.cos(math.radians(alt)), 0.05) / 2
                st = C["status"][i] if i < len(C["status"]) else "pending"
                th = C.get("thumbs", {}).get(i) if cv_["thumbs"].get() else None
                for sh in (0, 360, -360):
                    xc = XA(az % 360 + sh)
                    if xc + cw * kx < x0 - 50 or xc - cw * kx > x1 + 50:
                        continue
                    if th is not None and ImageTk is not None:
                        t, im = th
                        tx0, tx1 = XA(t["az0"] - az + az % 360 + sh), XA(t["az1"] - az + az % 360 + sh)
                        ty0, ty1 = Y(t["alt1"]), Y(t["alt0"])
                        wpx, hpx = int(round(tx1 - tx0)), int(round(ty1 - ty0))
                        if 2 <= wpx <= 4000 and 2 <= hpx <= 4000:
                            key = (i, wpx, hpx)                     # re-use the resized preview
                            ph = C["photo_cache"].get(key)
                            if ph is None:
                                ph = C["photo_cache"][key] = ImageTk.PhotoImage(im.resize((wpx, hpx)))
                            C["photos"].append(ph)
                            ccan.create_image(tx0, ty0, image=ph, anchor="nw")
                    col = STATUS_COL.get(st, "#4a7bd0")
                    wid = 3 if st in ("slewing", "rotating", "exposing") else (2 if st in ("error", "skipped") else 1)
                    if st == "done" and th is not None:
                        col = "#43a047"
                    elif st == "done" and not C.get("paths", {}).get(i):
                        col = "#9ccc65"                      # taken, file not in the folder yet
                    ccan.create_rectangle(xc - cw * kx, Y(alt + fh / 2), xc + cw * kx, Y(alt - fh / 2),
                                          outline=col, width=wid)
            if p is not None:                                      # current horizon file(s)
                for hf in (p.horizon_files or []) if p.mode == "horizon" else []:
                    try:
                        pts = load_hrz(hf)
                        a_ = np.arange(0, 360.01, 0.5)
                        ccan.create_line(*[c for a in a_ for c in (XA(a), Y(hrz_alt(pts, a % 360)))],
                                         fill="#d00000", width=2)
                    except Exception:
                        pass
            try:
                st_now = C.get("status") or []
                todo = [i for i in range(len(C["frames"])) if (st_now[i] if i < len(st_now) else "pending") != "done"]
                r_ = C.get("runner")
                if r_ is not None and r_.state in ("running", "paused"):
                    todo = [i for i in r_.order[r_.pos:]]
                key = (len(todo), todo[:1], int(time.time() // 30), cv_["sun"].get())
                if C.get("sunrisk", (None,))[0] != key:
                    C["sunrisk"] = (key, capture_sun_risk(todo))
                draw_sun_overlay(ccan, X, Y, C["sunrisk"][1], C["frames"], plan["fw"], plan["fh"],
                                 float(cv_["sun"].get()), x0, x1)
                C["sun_txt"] = sun_summary(C["sunrisk"][1], float(cv_["sun"].get()))
            except Exception:
                C["sun_txt"] = ""
            if cv_["order"].get():
                draw_order(ccan, plan, X, Y, tag="order", numbers="rows",
                           flip=any(z < 180 for _, z, _ in C["frames"]) and any(z >= 180 for _, z, _ in C["frames"]))
            if C["cur"] is not None and C["cur"] < len(C["frames"]):
                alt, az, _ = C["frames"][C["cur"]]
                x, y = X(az), Y(alt)
                ccan.create_oval(x - 7, y - 7, x + 7, y + 7, outline="#000", width=3)
        for r_ in ((0, 0, W, y0), (0, y1, W, H), (0, 0, x0, H), (x1, 0, W, H)):
            ccan.create_rectangle(*r_, fill="white", outline="")
        ccan.create_rectangle(x0, y0, x1, y1, outline="#bbb")
        if C.get("sun_txt"):
            ccan.create_text(x0, y0 - 10, anchor="w", font=("Segoe UI", 9), fill="#e65100", text=C["sun_txt"])
        for a in at:
            ccan.create_text(XA(a), y1 + 12, text=f"{a:g}°")
        for h in ht:
            ccan.create_text(x0 - 22, Y(h), text=f"{h:g}°")

    def cap_set_view(a0, a1, b0, b1):
        w = min(max(a1 - a0, 5.0), 360.0); h = min(max(b1 - b0, 3.0), 90.0)
        a0 = min(max(a0, 0.0), 360.0 - w); b0 = min(max(b0, 0.0), 90.0 - h)
        C["view"] = (a0, a0 + w, b0, b0 + h); draw_cap()

    def cap_wheel(e):
        x0, y0, x1, y1 = C["box"]; az0, az1, h0, h1 = C["view"]
        a = az0 + (e.x - x0) / (x1 - x0) * (az1 - az0); h = h0 + (y1 - e.y) / (y1 - y0) * (h1 - h0)
        f_ = 0.8 if (getattr(e, "num", 0) == 4 or getattr(e, "delta", 0) > 0) else 1.25
        cap_set_view(a - (a - az0) * f_, a + (az1 - a) * f_, h - (h - h0) * f_, h + (h1 - h) * f_)
        return "break"

    def cap_drag(e):
        if "drag" not in C:
            return
        sx, sy, (az0, az1, h0, h1) = C["drag"]; x0, y0, x1, y1 = C["box"]
        da = -(e.x - sx) / (x1 - x0) * (az1 - az0); dh = (e.y - sy) / (y1 - y0) * (h1 - h0)
        cap_set_view(az0 + da, az1 + da, h0 + dh, h1 + dh)

    ccan.bind("<Configure>", lambda e: draw_cap())
    ccan.bind("<MouseWheel>", cap_wheel); ccan.bind("<Button-4>", cap_wheel); ccan.bind("<Button-5>", cap_wheel)
    ccan.bind("<ButtonPress-1>", lambda e: C.__setitem__("drag", (e.x, e.y, C["view"])))
    ccan.bind("<B1-Motion>", cap_drag)
    ccan.bind("<Double-Button-1>", lambda e: cap_set_view(0, 360, 0, 90))

    # ------------------------------------------------------------------ HELP
    help_tab = ttk.Frame(nb); nb.add(help_tab, text="Help")
    t = tk.Text(help_tab, wrap="word"); t.pack(fill="both", expand=True)
    t.insert("1.0", HELP_TEXT); t.config(state="disabled")

    def show_logs(on):
        """Debug → Show log windows: the Capture log pane and the Build tab log."""
        try:
            panes = [str(p_) for p_ in (cap_pane.panes() or ())]
            if on and str(clog_box) not in panes:
                cap_pane.add(clog_box, weight=0)
            elif not on and str(clog_box) in panes:
                cap_pane.forget(clog_box)
            if on and not log_box.winfo_ismapped():
                log_box.pack(side="left", fill="y", padx=(6, 0), pady=2)
            elif not on:
                log_box.pack_forget()
        except Exception:
            dlog_exc("show_logs")
        dlog(f"log windows {'shown' if on else 'hidden'}")
    logs_var.trace_add("write", lambda *a: (save_settings(show_logs=logs_var.get()), show_logs(logs_var.get())))
    show_logs(logs_var.get())

    if no_tip:
        dlog(f"controls without hover text ({len(no_tip)}): {sorted(no_tip)}")

    def on_close():
        r = C.get("runner")
        if r and r.state in ("running", "paused"):
            if not messagebox.askyesno(APP_NAME, "A capture run is active.\n\nStop it (mount, rotator and "
                                                 "camera are stopped) and exit?", icon="warning"):
                nb.select(cap_tab); return
            r.stop()
        fp = lines_fingerprint()
        if fp is not None and fp != E.get("saved_fp"):
            if not messagebox.askyesno(APP_NAME, "The horizon has not been saved since it was built or "
                                                 "last edited.\n\nExit without saving it?", icon="warning"):
                nb.select(proc); return
        root.destroy()
    root.protocol("WM_DELETE_WINDOW", on_close)

    # Start wide enough for every settings group (and the fixed boxes), if the screen allows it.
    try:
        root.update_idletasks()
        need = max(inner.winfo_reqwidth() + fixed.winfo_reqwidth() for inner, fixed in bars) + 40
        sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
        w = min(max(need, 1280), sw - 20); h = min(900, sh - 80)
        root.geometry(f"{w}x{h}+{max((sw - w) // 2, 0)}+{max((sh - h) // 2 - 20, 0)}")
    except Exception:
        pass
    if selftest:
        # build the real window, show a hover, switch through the tabs, then close (used by CI)
        def _selftest():
            try:
                for t_ in nb.tabs():
                    nb.select(t_); root.update()
                hov = next((o for o in _selftest_hovers if o.w.winfo_exists()), None)
                if hov:
                    hov.show(); root.update(); hov.hide()
                print(f"GUI self-test OK: {len(nb.tabs())} tabs, {len(_selftest_hovers)} hover texts, "
                      f"{len(no_tip)} controls without one")
            except Exception:
                traceback.print_exc(); os._exit(1)
            root.destroy()
        root.after(1500, _selftest)
    root.mainloop()


# Hover explanations for every setting, choice and button (key = the text shown in the window).
# Add an entry for every new control – the debug log lists controls that still have none.
TIPS = {
    # ---- Plan tab: camera & telescope
    "Sensor width (px)": "Width of the camera sensor in pixels (ASI2600MC Pro: 6248).",
    "Sensor height (px)": "Height of the camera sensor in pixels (ASI2600MC Pro: 4176).",
    "Pixel size (µm)": "Size of one sensor pixel in micrometres (ASI2600MC Pro: 3.76).",
    "Focal length (mm)": "Focal length of the telescope (Askar SQA85: 408 mm). With the sensor size it gives the field of view.",
    "Orientation": "How the camera is turned: landscape, portrait, or auto (whichever needs fewer frames).",
    "Read from a FITS file…": "Fill in sensor size, pixel size and focal length from the header of one of your FITS files.",
    # ---- Plan tab: site & area
    "Latitude (°)": "Your site's latitude in degrees (north positive). Needed for the rotator and the Sun check.",
    "Fixed altitude band": "Image a band of constant altitude all the way round (set the two band limits).",
    "Around current horizon": "Image a strip that follows your current horizon file(s), with the margins below and above it.",
    "Band: alt min (°)": "Lowest altitude of the band to image (fixed-band mode).",
    "Band: alt max (°)": "Highest altitude of the band to image (fixed-band mode).",
    "Add current horizon…": "Add a NINA .hrz horizon file. The first one is your current horizon; a second one (e.g. an all-sky estimate) widens the strip.",
    "Clear": "Remove the horizon files from the list.",
    "Margin below (°)": "How far below the current horizon the strip starts.",
    "Margin above (°)": "How far above the current horizon the strip ends.",
    "Lowest altitude (°)": "No frames are planned below this altitude.",
    "Overlap (0-1)": "Fraction by which neighbouring frames overlap (0.1 = 10 %). More overlap = more frames, better stitching.",
    "Extra frames at steep edges": "Where the horizon rises or drops steeply, widen the strip sideways by this many frame widths (0 = off).",
    "Steep = more than (°/°)": "A horizon edge counts as steep above this many degrees of altitude per degree of azimuth (2 ≈ 63°).",
    "Azimuth start (°)": "Azimuth where the planned area starts (0 = north, 90 = east).",
    "Azimuth end (°)": "Azimuth where the planned area ends. Same as start = all the way round.",
    # ---- Plan tab: rotator
    "Keep frames level with rotator": "Turn the rotator for every row so that each frame stays level with the horizon.",
    "Offset (sky - mech, °)": "Rotator offset: sky angle minus mechanical angle (measured for your setup: -58.04°).",
    "Direction (+1/-1)": "Whether the sky angle grows (+1) or shrinks (-1) when the mechanical angle grows.",
    "Skip moves below (°)": "Rotator moves smaller than this are skipped to save time (the frames then tilt slightly).",
    # ---- Plan tab: exposure & sequence
    "Exposure (s)": "Exposure time of each frame in seconds (daylight: a few milliseconds).",
    "Gain (-1 = default)": "Camera gain for each frame; -1 uses the camera's default in NINA.",
    "Offset (-1 = default)": "Camera offset for each frame; -1 uses the camera's default in NINA.",
    "Binning": "Pixel binning (2 = 2×2). Binning shortens downloads and is plenty for horizon mapping.",
    "Seconds per frame": "Time estimate per frame (slew, settle, exposure, download) – only used for the run-time estimate.",
    "Seconds per rotator move": "Time estimate per rotator move – only used for the run-time estimate.",
    "East half first, then west (one meridian flip)": "Shoot all frames east of the meridian first, then the west half, so the mount flips only once.",
    "Rows bottom to top in both halves": "Also shoot the west half from the bottom up (default: west half top to bottom, the shorter path).",
    "Unpark at start": "Unpark the mount before the first frame.",
    "Park at end": "Park the mount when the sequence has finished.",
    "Save as template (else full sequence)": "Save a NINA template (drag into a sequence) instead of a complete sequence file.",
    "Rows collapsed (use Load sequence)": "Save the rows collapsed in NINA's sequencer (shown collapsed only with Load sequence).",
    "Preview plan": "Work out the plan with the current settings and draw it.",
    "Save NINA sequence…": "Save the plan as a NINA Advanced Sequencer file.",
    "Show shooting order": "Draw the order in which the frames are taken (green = first → magenta = last).",
    "Sun check": "Show where the Sun is at the given time and whether any frame is too close to it.",
    "Longitude (°E)": "Your site's longitude in degrees east – used for the Sun check.",
    "Start date": "Planned start date: 2026-11-12, 12.11.2026 or 12.11. (this year); blank = today. With the start time it sets where the preview draws the Sun.",
    "Start time (HH:MM)": "Planned start time (local; blank = the current time of day). The preview shows the Sun's path from then until the plan's estimated end.",
    "Sun distance (°)": "Keep-away distance from the Sun for the preview: frames inside it are drawn red, slews passing it get a red cross.",
    # ---- Capture tab
    "Host": "Address of the computer running NINA (localhost = this computer, or its IP address, e.g. 192.168.1.55).",
    "Port": "Port of NINA's Advanced API plugin (default 1888).",
    "Test connection": "Check that NINA's Advanced API answers, show the devices and find NINA's image folder.",
    "Plan from the Plan tab": "Shoot the plan made with the Plan tab's settings.",
    "Fill-in frames from the Build tab": "Shoot the fill-in frames generated on the Build tab.",
    "Saved sequence file": "Shoot the frames of a sequence saved for NINA (plan or fill-in).",
    "Browse…": "Choose a file.",
    "Gain (-1 = NINA default)": "Camera gain for the frames; -1 keeps NINA's current camera setting.",
    "Settle after moves (s)": "Wait this long after each slew or rotator move before exposing.",
    "Keep away from the Sun (°)": "Frames closer to the Sun than this are skipped and never slewed to.",
    "Never below altitude (°)": "Frames below this altitude are skipped and never slewed to.",
    "Unpark at start ": "Unpark the mount before the first frame.",
    "Park when finished": "Park the mount when the run has finished (not after End run).",
    "From frame": "First frame to shoot (frame numbers follow the shooting order shown in the picture).",
    "To frame (blank = last)": "Last frame to shoot; leave blank for the last frame of the plan.",
    "Re-shoot frames already taken": "Also shoot frames in the range that are already taken (normally they are left out).",
    "Target name": "Target name NINA uses for the files (and their folder, depending on NINA's file pattern).",
    "NINA image folder (found automatically)": "Folder where NINA saves the images. Filled in by Test connection; can be changed at any time.",
    "Show finished frames in the picture": "Draw each taken frame's picture in its slot.",
    "Turn previews 180°": "Turn the frame pictures upside down if they appear the wrong way round.",
    "⏸  PAUSE  (Esc)": "Stop slew, rotator and exposure at once and hold. Esc does the same from anywhere in the app.",
    "▶ Start": "Shoot the frames in the range that are not taken yet (asks for confirmation first).",
    "▶ Resume": "Continue a paused run; the interrupted frame is shot again from its start.",
    "Skip frame": "While paused: leave out the current frame and continue with the next one on Resume.",
    "■ End run": "End the run. Frames taken so far are kept; Start later shoots what is still missing.",
    "Save current capture…": "Save the frames, the progress and the file paths, so the capture can be resumed later.",
    "Load saved capture…": "Load a saved capture: its progress, settings and the pictures of the frames already taken.",
    # ---- Build tab
    "Import from Capture…": "Bring a capture's frames (from the Capture tab or a saved capture), the rotator offset and the current horizon here, then build.",
    "Add folder…": "Add a folder with FITS frames (searched recursively). A capture list (.txt) works too.",
    "Remove": "Remove the selected folder from the list.",
    "Rotator offset guess (°)": "Starting value for the rotator offset; the calibration refines it.",
    "Calibrate rotation from overlaps": "Measure the rotator offset from how neighbouring frames overlap.",
    "Full 360° search (offset unknown)": "Search all rotator offsets (slower) – use when the offset is unknown.",
    "Latitude (°, blank = FITS)": "Site latitude; leave blank to read it from the FITS headers.",
    "Mosaic resolution (°)": "Size of one mosaic pixel in degrees (0.03 is fine; smaller = sharper but slower).",
    "Solid: darker than": "A pixel counts as solid obstruction when darker than this fraction of the sky brightness.",
    "Branches: darker than": "A pixel counts as branches when darker than this fraction of the sky brightness.",
    "Bias level (blank = auto)": "Camera bias (black level) subtracted from the frames; blank = estimated automatically.",
    "Cache reduced frames": "Keep reduced copies of the frames, so building again is much faster.",
    "Clear cache": "Delete the cached reduced frames.",
    "Build mosaic & horizon": "Read the frames, stitch the mosaic and detect the horizon lines.",
    "Re-detect (new thresholds)": "Detect the horizon again with new thresholds, without rebuilding the mosaic.",
    "Extra frames at steep edges ": "Fill-in: how many frame widths to add beside steep rises and drops of the new horizon (0 = off).",
    "steep = more than": "An edge counts as steep above this many degrees of altitude per degree of azimuth.",
    "Horizon line near gaps within (°)": "Fill-in: also cover places where the horizon line is this close to a spot without pictures (0 = off).",
    "Generate fill-in frames": "Plan extra frames where the horizon is uncertain, near gaps or at steep edges, and show them.",
    "Save fill-in sequence…": "Save the fill-in frames as a NINA sequence (or shoot them on the Capture tab).",
    "Save horizon file…": "Save the selected horizon line as a NINA .hrz file.",
    "Edit:": "Which line a left-drag on the image redraws.",
    "off": "No editing: dragging on the image does nothing.",
    "solid line": "Redraw the solid-horizon line by left-dragging on the image.",
    "branches line": "Redraw the branches line by left-dragging on the image.",
    "brush width (°)": "Width in azimuth of each brush stroke when redrawing a line.",
    "zoom (px/°)": "Image scale in pixels per degree (the mouse wheel zooms too).",
    "Undo": "Undo the last line edit (Ctrl+Z).",
    "Reset edits": "Remove all line edits and go back to the detected lines.",
    "Side margin at steep edges (°)": "Saved horizon: widen it sideways by this many degrees at steep edges (0 = off).",
    "Show fill-in frames": "Draw the fill-in areas, frames and their order on the image.",
    "Current horizon (compare / fill gaps)": "Your current horizon: drawn for comparison and used for azimuths without data.",
}


HELP_TEXT = """HORIZON MAPPER FOR NINA   (the version is shown in the window title)

WORKFLOW   1. Plan sequence → 2. Capture (NINA) → 3. Build horizon
1. Plan sequence
   • Enter camera/telescope (or 'Read from a FITS file…'). The field of view is computed.
   • Choose the area: a fixed altitude band, or a strip around your current horizon file
     (optionally plus a second estimate, e.g. from an all-sky camera); the strip spans both,
     plus margins. The current horizon is also filled in on the Build tab for comparison.
   • 'Extra frames at steep edges' (horizon mode): where the current horizon rises or drops
     more steeply than 'Steep = more than' (altitude ° per azimuth °, 2 ≈ 63°), the strip is
     widened sideways by this many frame widths and covers the whole height of the jump, so
     the sides of trees and roofs are imaged. Default 0.5; 0 = off.
     Fractions work: 0.5 usually gives two frames per row at a narrow edge.
   • Each row gets frames only where the strip enters it, as few as the overlap allows.
   • 'Lowest altitude' keeps frames (and the sun) out of the low sky.
   • 'Preview plan' draws the frames; with 'Show shooting order' (on by default) a colour-graded
     path shows the order (green = first → yellow → magenta = last), with the first frame number
     of each row, Start/End markers and the meridian flip.
     Mouse wheel zooms around the pointer, drag moves the view, double-click shows the whole sky.
   • The settings sit in a bar across the top, with the picture below using the full width.
     The bar scrolls sideways (wheel or scrollbar) when the window is too narrow; the buttons
     on its right ('Plan' here, the log on the Build tab) always stay visible. At start the window
     is made wide enough for all groups if the screen allows.
   • Closing the window asks first if the horizon has not been saved since it was built or edited.
   • Order: 'East half first, then west' = one meridian flip; east half bottom→top, west half
     top→bottom unless 'Rows bottom to top in both halves' is ticked.
   • Rotator: keeps every frame level with the horizon, so overlap can be ~10 %.
     Offset = sky position angle − mechanical angle (from NINA's rotator panel after a
     plate solve). A wrong 180° doesn't matter here.
   • East-half-first ordering means the mount flips over the meridian only once.
   • Preview, then 'Save NINA sequence…'. Open it with Load sequence in the Advanced
     Sequencer (or save as template and put it in NINA's template folder).
   • Daytime under even cloud gives the best silhouettes. Take a test frame first:
     open sky should read ~10,000–30,000 ADU. Keep the lowest altitude well above the sun.
   • The Sun in the preview: set 'Start date' (e.g. 12.11.2026; blank = today), 'Start time
     (HH:MM)' (blank = now) and 'Sun distance (°)'. The
     preview draws the Sun's track from the start to the plan's estimated end with the keep-away
     zone; frames too close at their time are red, slews passing near it get a red ×.

2. Capture (NINA)  – shoot the frames from here, through NINA's Advanced API plugin
   • Setup once: install 'Advanced API' in NINA (Plugins), enable its server (port 1888), let
     NINA save FITS, connect mount, camera and rotator in NINA. No NINA sequence at the same time.
   • 'Test connection' (in the NINA group) shows the devices, starts the live device line and finds NINA's image
     folder. 'Frames to shoot': the Plan tab plan, the Build tab fill-in frames, or a saved
     sequence file (Browse…; exposure, gain and binning come from it). Choosing a source shows its
     frames; a fresh Start takes them from the source as it is then.
   • Range: 'From frame' / 'To frame' (numbers in shooting order) limit the run; frames already
     taken are left out unless 'Re-shoot frames already taken' is ticked.
   • 'Save current capture…' stores frames, progress, file paths and settings (.hmcap.json);
     'Load saved capture…' restores them, shows the pictures from the original paths and sets
     frames whose file is gone back to 'to do'. Start then shoots only what is missing.
     Progress is also saved automatically after every frame (captures\\last_capture.hmcap.json).
   • The picture is drawn like the Plan preview (strip, horizon, frames, shooting order); blue =
     to do, yellow = in progress, green = done, orange = skipped, red = error. Finished frames
     are shown in their slots (projected like the mosaic; 'Turn previews 180°' if upside down).
   • Safety: frames nearer the Sun than the limit (default 30°) or below the altitude limit
     (default 10°) are skipped and never slewed to. Unpark at start / park when finished.
     Slew paths are checked as well: a slew that would pass near the Sun goes via a safe
     intermediate position (or the frame is skipped). At Start, frames that will be too close
     to the Sun during the run are listed: Yes = shoot as planned, No = leave them out.
   • 'Start' asks first. Per frame: slew → rotator → settle → expose; NINA saves the FITS.
   • ⏸ PAUSE (red button) or Esc anywhere: stops slew, rotator and exposure at once and holds.
     'Resume' repeats the interrupted frame; 'Skip frame' skips it; 'End run' ends the run.
   • Any error or timeout pauses the run automatically – nothing is skipped silently.
   • The NINA image folder can be changed during a run – file look-ups always use the current
     value.
   • The log runs down the right side of the tab (scrollable). The progress bar shows the
     percentage of the run. 'Save / Load Capture' sits right of 'Run'.
   • A preview is made only when the FITS file is complete (NINA may still be writing it);
     otherwise it is tried again every 3 s for up to 3 minutes.
   • Files that NINA saves late are picked up by a scan of its image folder every few seconds
     (matched by name or by the pointing in the header; light green = file not arrived yet), also
     for 10 minutes after the run.
   • At the end the list of files is added to the Build tab ('Start' again can continue a
     stopped run).

HOVER over any setting, choice or button for a short explanation.

DEBUG MENU
   • 'Debug log (all actions)' writes every click, NINA request (with timing and answer), capture
     step, pause, error and window stall to %LOCALAPPDATA%\\HorizonMapper\\logs. On by default
     for now. 'Open log folder' opens it – please attach the log when reporting a problem.
   • 'Show log windows' shows or hides the log panes (the Capture log can be made wider or
     narrower by dragging its divider).

3. Build horizon
   • 'Import from Capture…' adds a capture's files (Capture tab or saved .hmcap.json), the
     rotator offset and the current horizon here, then offers to build.
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
   • 'Fill-in: horizon line near gaps within (°)': wherever the solid or branches line runs
     through a spot without pictures, or comes within this many degrees of one, a 'Gap …' area
     is added – also where the horizon is not steep. Only frames over the gap are shot.
     Default 1°, 0 = off.
   • 'Show fill-in frames' (toolbar above the image) draws the fill-in areas (dashed magenta),
     their frames (orange) and the shooting order (green = first → magenta = last, numbers =
     first frame of each row; every frame is numbered at zoom 20). It follows your edits and
     the fill-in settings: after a change the plan is generated again and redrawn.
   • 'Re-detect' applies new thresholds instantly without re-reading the frames.
   • 'Side margin at steep edges' (toolbar above the image): where the horizon rises or drops steeply (the side of a
     tree, the end of a roof), the higher value is extended that many degrees sideways, so a
     target just beside the obstruction isn't scheduled too early. Gentle slopes are not
     changed. 'Steep' means more than the given degrees of altitude per degree of azimuth
     (2 ≈ a 63° slope; default 2). 1-2° is a sensible margin; 0 switches it off. The lines on screen and
     the saved files include the margin; your edits are made on the unwidened line.
   • The two lines are edited independently – the branches line may be drawn below the solid
     line (e.g. where the "solid" part is really a thin, see-through tree top).
   • Editing: choose 'solid line' or 'branches line' above the image, then left-drag to
     redraw the line (brush width in degrees), right-click to restore the detected value,
     Ctrl+Z to undo. Edits are used when saving.
   • Mouse wheel zooms the image around the pointer (2–24 px/°, or pick a level in 'zoom');
     Shift+wheel, middle-drag or the scrollbars move it.
   • Reduced frames are cached (in %LOCALAPPDATA%\\HorizonMapper\\cache), so building again
     – e.g. after adding fill-in frames – only reads the new files. 'Clear cache' frees the space.
   • 'Save horizon file…' (below the fill-in buttons) saves one NINA .hrz file: choose 'solid'
     or 'with branches' next to the button. Edits and the side margin are included.
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
    ap.add_argument("--selftest-gui", action="store_true",
                    help="open the window, show a hover, switch through the tabs and close (for CI)")
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
    pp.add_argument("--edge-frames", type=float, default=0.5,
                    help="extra frame widths beside steep edges of the horizon file(s) (0 = off)")
    pp.add_argument("--edge-slope", type=float, default=2.0,
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
    pp.add_argument("--bottom-up", action="store_true",
                    help="shoot the rows bottom to top in the west half too (default: west half top to bottom)")
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
    bp.add_argument("--steep-slope", type=float, default=2.0,
                    help="an edge is steep above this many degrees altitude per degree azimuth")
    bp.add_argument("--no-cache", action="store_true", help="don't use/create the frame cache")
    bp.add_argument("--fillin", metavar="FILE.json", help="also write a NINA fill-in sequence for uncertain areas")
    bp.add_argument("--fillin-exposure", type=float, default=0.008)
    bp.add_argument("--fillin-floor", type=float, default=15.0)
    bp.add_argument("--fillin-edge-frames", type=float, default=0.5,
                    help="for --fillin: extra frame widths beside steep horizon edges (0 = off)")
    bp.add_argument("--fillin-gap-near", type=float, default=1.0,
                    help="for --fillin: also cover places where a horizon line is within this many degrees "
                         "of a spot without pictures (0 = off)")
    bp.add_argument("--fillin-edge-slope", type=float, default=2.0,
                    help="for --fillin: steep = more than this many deg altitude per deg azimuth")
    bp.add_argument("--sensor", nargs=3, type=float, metavar=("WPX", "HPX", "UM"), default=[6248, 4176, 3.76],
                    help="for --fillin: unbinned sensor size and pixel size")
    bp.add_argument("--focal", type=float, default=408, help="for --fillin: focal length (mm)")
    bp.add_argument("-o", "--outbase", default="horizon")
    a = ap.parse_args()

    if a.cmd is None:
        run_gui(selftest=a.selftest_gui); return
    if a.cmd == "plan":
        p = PlanParams(sensor_w_px=int(a.sensor[0]), sensor_h_px=int(a.sensor[1]), pixel_um=a.sensor[2],
                       focal_mm=a.focal, lat=a.lat, floor=a.floor, az_start=a.az[0], az_end=a.az[1],
                       overlap=a.overlap, rotator=not a.no_rotator, rot_offset=a.rot_offset,
                       rot_tolerance=a.rot_tol, exposure=a.exposure, gain=a.gain, binning=a.bin,
                       template=a.template, margin_below=a.below, margin_above=a.above,
                       edge_frames=a.edge_frames, edge_slope=a.edge_slope, bottom_up=a.bottom_up)
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
            boxes = boxes + gap_boxes(mos, det["az"], (det["solid"], det["branches"]), a.fillin_gap_near,
                                      floor=p.floor, ceiling=p.ceiling)
            plan = plan_fillin(p, boxes, mos=mos) if boxes else None
            if plan and plan["frames"]:
                n_edge = len({l for l in plan["labels"] if l.startswith("Steep edge")})
                n_gap = len({l for l in plan["labels"] if l.startswith("Gap")})
                print(f"fill-in: {save_plan(p, plan, a.fillin)} ({len(plan['frames'])} frames, "
                      f"{n_edge} steep-edge and {n_gap} gap areas, ~{plan['seconds'] / 60:.0f} min)")
            else:
                print("fill-in: nothing to shoot")
        print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
