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
import argparse, json, math, os, queue, sys, threading, time, traceback

import numpy as np
try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:  # GUI will complain later
    Image = None

APP_NAME = "Horizon Mapper for NINA"

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
    """Return function az -> (lo, hi) of the strip to image."""
    if p.mode == "horizon" and p.horizon_files:
        hz = [load_hrz(f) for f in p.horizon_files]
        def f(az):
            v = [hrz_alt(h, az) for h in hz]
            return max(min(v) - p.margin_below, p.floor), min(max(v) + p.margin_above, p.ceiling)
        return f
    return lambda az: (max(p.alt_min, p.floor), min(p.alt_max, p.ceiling))


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
    grid = []
    for alt in rows:
        low_edge = max(alt - fh / 2, 0.0)
        cov = fw / max(math.cos(math.radians(low_edge)), 1e-3)
        stp = cov * (1 - p.overlap)
        if full:
            n = math.ceil(360 / stp); azs = [p.az_start + k * 360 / n for k in range(n)]
        else:
            n = 1 if span_az <= cov else math.ceil((span_az - cov) / stp) + 1
            a0 = p.az_start + cov / 2; a1 = p.az_start + span_az - cov / 2
            azs = [a0 + k * (a1 - a0) / max(n - 1, 1) for k in range(n)] if n > 1 else [p.az_start + span_az / 2]
        keep = []
        for az in azs:
            los, his = zip(*[lim(az + cov / 2 * k / 4) for k in range(-4, 5)])
            if alt + fh / 2 >= min(los) and alt - fh / 2 <= max(his):
                keep.append(az % 360)
        if keep:
            grid.append((alt, sorted(keep)))
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
    secs = len(frames) * (p.sec_per_frame + p.exposure) + rot_moves * p.sec_per_rot
    return dict(grid=grid, frames=frames, orientation=orient, fw=fw, fh=fh, rot_moves=rot_moves,
                seconds=secs, lim=lim)


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
    for alt, azs in plan["grid"]:
        side = "E" if azs and azs[0] < 180 else "W"
        row, rowid = _container(ids, "SequentialContainer",
                                f"Alt {alt:.1f}° {side} ({len(azs)} frames)", bid, not p.collapsed)
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


def load_frames(paths, log, progress, bias=None):
    frames = []; skipped = 0; lat = None
    for i, pth in enumerate(paths):
        try:
            h, a = read_fits(pth)
            if "CENTALT" not in h or "CENTAZ" not in h:
                skipped += 1; continue
            if str(h.get("IMAGETYP", "LIGHT")).upper() not in ("LIGHT", "SNAPSHOT"):
                skipped += 1; continue
            bayer = bool(str(h.get("BAYERPAT", "")).strip())
            img, k = reduce_frame(a, bayer)
            fr = Frame(); fr.path = pth
            fr.alt = float(h["CENTALT"]); fr.az = float(h["CENTAZ"])
            fr.rot = float(h.get("ROTATOR", h.get("ROTATANG", 0.0)) or 0.0)
            fr.pier = h.get("PIERSIDE", "")
            fr.t = h.get("DATE-LOC", "")
            pix = float(h.get("XPIXSZ", 0) or 0); foc = float(h.get("FOCALLEN", 0) or 0)
            if pix <= 0 or foc <= 0:
                raise ValueError("XPIXSZ/FOCALLEN missing in header")
            fr.scale = 206.265 * pix / foc / 3600.0 * k        # deg per reduced pixel
            fr.img = img
            if lat is None and "SITELAT" in h:
                lat = float(h["SITELAT"])
            frames.append(fr)
        except Exception as e:
            skipped += 1
            log(f"  skipped {os.path.basename(pth)}: {e}")
        if i % 20 == 0:
            progress(i / max(len(paths), 1) * 0.35, f"Reading frames {i}/{len(paths)}")
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
    root.title(APP_NAME)
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
    ttk.Radiobutton(box, text="Around horizon file(s)", variable=V["mode"], value="horizon").grid(row=1, column=1, sticky="w")
    labeled(box, 2, "Band: alt min (°)", var("alt_min", P.alt_min))
    labeled(box, 3, "Band: alt max (°)", var("alt_max", P.alt_max))
    hz_list = tk.Listbox(box, height=2, width=34); hz_list.grid(row=4, column=0, columnspan=2, sticky="w", padx=4)

    def add_hz():
        for f in filedialog.askopenfilenames(title="Horizon file(s)", filetypes=[("Horizon", "*.hrz *.txt *.csv"), ("All", "*.*")]):
            if hz_list.size() < 2:
                hz_list.insert("end", f)
    ttk.Button(box, text="Add .hrz…", command=add_hz).grid(row=5, column=0, sticky="w", padx=4)
    ttk.Button(box, text="Clear", command=lambda: hz_list.delete(0, "end")).grid(row=5, column=1, sticky="w")
    ttk.Label(box, text="Two files (e.g. old + all-sky estimate): the strip spans both.", foreground="#666").grid(row=6, column=0, columnspan=3, sticky="w", padx=4)
    labeled(box, 7, "Margin below (°)", var("margin_below", P.margin_below))
    labeled(box, 8, "Margin above (°)", var("margin_above", P.margin_above))
    labeled(box, 9, "Lowest altitude (°)", var("floor", P.floor))
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

    def collect():
        p = PlanParams()
        num = ["sensor_w_px", "sensor_h_px", "pixel_um", "focal_mm", "lat", "alt_min", "alt_max", "margin_below",
               "margin_above", "floor", "az_start", "az_end", "overlap", "rot_offset", "rot_tolerance", "exposure",
               "sec_per_frame", "sec_per_rot"]
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
        if p.mode == "horizon" and not p.horizon_files:
            raise ValueError("Add at least one horizon file, or choose 'Fixed altitude band'.")
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
    folders = tk.Listbox(box, height=6, width=48); folders.pack(fill="x", padx=4, pady=2)
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
    box = ttk.LabelFrame(pl, text="Compare"); box.pack(fill="x", pady=3)
    old_var = tk.StringVar(value="")
    ttk.Entry(box, textvariable=old_var, width=40).pack(side="left", padx=4)
    ttk.Button(box, text="Old .hrz…", command=lambda: old_var.set(filedialog.askopenfilename(
        filetypes=[("Horizon", "*.hrz *.txt"), ("All", "*.*")]) or old_var.get())).pack(side="left")
    pb = ttk.Progressbar(pl, length=340, maximum=1.0); pb.pack(fill="x", pady=4)
    status = ttk.Label(pl, text="Idle"); status.pack(anchor="w")
    logbox = tk.Text(pl, height=14, width=58, wrap="word"); logbox.pack(fill="both", expand=True, pady=3)
    img_lbl = ttk.Label(pr); img_lbl.pack(fill="both", expand=True)
    stats_lbl = ttk.Label(pr, text="", justify="left"); stats_lbl.pack(anchor="w")
    q = queue.Queue(); result = {}

    def log(m): q.put(("log", m))
    def progress(f, m): q.put(("prog", (f, m)))

    def worker(opts):
        try:
            paths = find_fits(opts["folders"])
            log(f"Found {len(paths)} FITS files")
            frames, lat_f = load_frames(paths, log, progress, opts["bias"])
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
            if det["uncertain"]:
                log("Check these azimuths (edge of imaged area): " +
                    ", ".join(f"{a:.1f}° {f.replace('_', ' ')}" for a, f in det["uncertain"][:30]))
            q.put(("done", dict(mos=mos, det=det, K=K, lat=lat, n=len(frames))))
        except Exception as e:
            q.put(("error", f"{e}\n{traceback.format_exc()}"))

    def start():
        if not folders.size():
            messagebox.showerror(APP_NAME, "Add at least one folder with FITS frames."); return
        try:
            opts = dict(folders=list(folders.get(0, "end")), K=float(pv["K"].get()), auto=pv["auto"].get(),
                        full=pv["full"].get(), lat=float(pv["lat"].get()) if pv["lat"].get().strip() else None,
                        res=float(pv["res"].get()), solid=float(pv["solid"].get()), thin=float(pv["thin"].get()),
                        bias=float(pv["bias"].get()) if pv["bias"].get().strip() else None)
        except ValueError as e:
            messagebox.showerror(APP_NAME, f"Check settings: {e}"); return
        logbox.delete("1.0", "end"); result.clear()
        threading.Thread(target=worker, args=(opts,), daemon=True).start()

    def current_det():
        det = result["det"]
        if old_var.get().strip():
            try:
                det = fill_unmeasured(det, load_hrz(old_var.get()))
            except Exception as e:
                log(f"Old horizon: {e}")
        return det

    def lines_for_display():
        det = current_det()
        L = [("Solid obstructions", det["az"], det["solid"], (230, 40, 30)),
             ("Incl. thin branches", det["az"], det["branches"], (0, 200, 240))]
        if old_var.get().strip():
            try:
                pts = load_hrz(old_var.get()); az = np.arange(0, 360, 0.25)
                L.append(("Old horizon", az, hrz_array(pts, az), (240, 200, 0)))
            except Exception as e:
                log(f"Old horizon: {e}")
        return L

    def show_result():
        if "det" not in result:
            return
        im = render_mosaic(result["mos"], lines_for_display(), scale=2)
        result["image"] = im
        if ImageTk:
            w = max(pr.winfo_width() - 10, 600)
            disp = im.resize((w, int(im.height * w / im.width)))
            result["tk"] = ImageTk.PhotoImage(disp); img_lbl.config(image=result["tk"])
        det = current_det()
        a_s = sky_area_sqdeg(det["az"], det["solid"]); a_b = sky_area_sqdeg(det["az"], det["branches"])
        txt = (f"{result['n']} frames · rotation offset {result['K']:.1f}° · latitude {result['lat']:.3f}°\n"
               f"Usable sky above horizon:  solid {a_s:,.0f} sq.deg ({a_s / 206.26:.1f}% of hemisphere)   "
               f"with branches {a_b:,.0f} sq.deg ({a_b / 206.26:.1f}%)")
        if old_var.get().strip():
            try:
                pts = load_hrz(old_var.get()); az = np.arange(0, 360, 0.1)
                a_o = sky_area_sqdeg(az, hrz_array(pts, az))
                txt += (f"\nOld horizon {a_o:,.0f} sq.deg → with branches {a_b - a_o:+,.0f} "
                        f"({(a_b / a_o - 1) * 100:+.1f}%), solid {a_s - a_o:+,.0f} ({(a_s / a_o - 1) * 100:+.1f}%)")
            except Exception:
                pass
        if not result["det"]["measured"].all():
            frac = result["det"]["measured"].mean() * 100
            txt += (f"\nOnly {frac:.0f}% of azimuths were imaged; the rest "
                    + ("uses the old horizon." if old_var.get().strip() else "is interpolated – load the old horizon to fill it."))
        stats_lbl.config(text=txt)

    def save_outputs():
        if "det" not in result:
            messagebox.showerror(APP_NAME, "Build the horizon first."); return
        d = filedialog.askdirectory(title="Folder for results")
        if not d:
            return
        base = os.path.join(d, time.strftime("horizon_%Y-%m-%d"))
        det = current_det()
        write_hrz(base + ".hrz", det["az"], det["solid"])
        write_hrz(base + "_with_branches.hrz", det["az"], det["branches"])
        render_mosaic(result["mos"], lines_for_display(), scale=2).save(base + "_mosaic.png")
        with open(base + "_details.json", "w") as f:
            json.dump(dict(rotation_offset=result["K"], latitude=result["lat"], details=det["details"]), f, indent=1)
        messagebox.showinfo(APP_NAME, f"Saved to {d}:\n{os.path.basename(base)}.hrz\n"
                                      f"{os.path.basename(base)}_with_branches.hrz\n"
                                      f"{os.path.basename(base)}_mosaic.png\n{os.path.basename(base)}_details.json")

    bb = ttk.Frame(pl); bb.pack(fill="x", pady=4)
    ttk.Button(bb, text="Build mosaic & horizon", command=start).pack(side="left", padx=3)
    ttk.Button(bb, text="Refresh view", command=show_result).pack(side="left", padx=3)
    ttk.Button(bb, text="Save results…", command=save_outputs).pack(side="left", padx=3)

    def poll():
        try:
            while True:
                kind, val = q.get_nowait()
                if kind == "log":
                    logbox.insert("end", val + "\n"); logbox.see("end")
                elif kind == "prog":
                    pb["value"] = val[0]; status.config(text=val[1])
                elif kind == "done":
                    result.update(val); pb["value"] = 1.0; status.config(text="Done"); show_result()
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


HELP_TEXT = """HORIZON MAPPER FOR NINA

WORKFLOW
1. Plan sequence
   • Enter camera/telescope (or 'Read from a FITS file…'). The field of view is computed.
   • Choose the area: a fixed altitude band, or a strip around one or two horizon files
     (e.g. your current horizon and a rough estimate). The strip spans both, plus margins.
   • 'Lowest altitude' keeps frames (and the sun) out of the low sky.
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
   • Azimuths listed as 'clear to bottom' or 'blocked to top' reached the edge of the
     imaged area – shoot a few extra frames there for a precise value.
   • 'Save results…' writes two .hrz files, a mosaic PNG and a details JSON.
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
    bp.add_argument("-o", "--outbase", default="horizon")
    a = ap.parse_args()

    if a.cmd is None:
        run_gui(); return
    if a.cmd == "plan":
        p = PlanParams(sensor_w_px=int(a.sensor[0]), sensor_h_px=int(a.sensor[1]), pixel_um=a.sensor[2],
                       focal_mm=a.focal, lat=a.lat, floor=a.floor, az_start=a.az[0], az_end=a.az[1],
                       overlap=a.overlap, rotator=not a.no_rotator, rot_offset=a.rot_offset,
                       rot_tolerance=a.rot_tol, exposure=a.exposure, gain=a.gain, binning=a.bin,
                       template=a.template, margin_below=a.below, margin_above=a.above)
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
        frames, lat_f = load_frames(find_fits(a.folders), log, prog)
        lat = a.lat if a.lat is not None else lat_f
        K = a.rot_offset if a.no_calibrate else calibrate_rotation(frames, lat, a.rot_offset, log, a.full_search)
        mos = build_mosaic(frames, lat, K, a.res, log, prog)
        det = detect_horizon(mos, log=log)
        if a.old:
            det = fill_unmeasured(det, load_hrz(a.old))
        write_hrz(a.outbase + ".hrz", det["az"], det["solid"])
        write_hrz(a.outbase + "_with_branches.hrz", det["az"], det["branches"])
        lines = [("Solid obstructions", det["az"], det["solid"], (230, 40, 30)),
                 ("Incl. thin branches", det["az"], det["branches"], (0, 200, 240))]
        if a.old:
            az = np.arange(0, 360, 0.25); lines.append(("Old horizon", az, hrz_array(load_hrz(a.old), az), (240, 200, 0)))
        render_mosaic(mos, lines).save(a.outbase + "_mosaic.png")
        print(f"usable sky: solid {sky_area_sqdeg(det['az'], det['solid']):.0f} sq.deg, "
              f"branches {sky_area_sqdeg(det['az'], det['branches']):.0f} sq.deg")
        if det["uncertain"]:
            print("check:", ", ".join(f"{x:.1f}°({f})" for x, f in det["uncertain"]))
        print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
