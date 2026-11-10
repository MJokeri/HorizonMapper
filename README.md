# Horizon Mapper for NINA

The version is shown in the window title (currently 1.0.1).

Plans a NINA Advanced Sequencer run that photographs the sky around your horizon, then turns the
FITS frames into a measured NINA horizon file (`.hrz`).

## Install (Windows)

**Easiest:** download `HorizonMapper-windows.zip` from the repository's
[Releases](../../releases) page. If there's no release yet, open the latest
[Actions → Build Windows exe](../../actions) run and download the build from there.
Unzip it and run `HorizonMapper.exe`. No Python is needed.

**From source:** install Python 3.10 or newer from https://www.python.org (tick **"Add python.exe to PATH"**).
Then double-click **`Start Horizon Mapper.bat`**. The first time, it installs `numpy` and `Pillow`.
`build_exe.bat` builds the `.exe` yourself.

## 1. Plan sequence

- **Camera & telescope:** use **Read from a FITS file…** with any frame from the setup. It fills
  in sensor, pixel size, focal length and latitude.
- **Area:**
  - **Fixed altitude band**, or
  - **Around current horizon:** under **Current horizon (.hrz)**, add the horizon file you use now. An
    optional second file (for example an all-sky estimate) widens the strip to cover both, plus the
    margins. The first file is also filled in automatically on the Build tab.
  - **Extra frames at steep edges** (horizon mode): where the current horizon rises or drops more
    steeply than **Steep = more than (°/°)** (default 3° altitude per 1° azimuth, about 72°), the strip
    is widened sideways by this many frame widths on both sides and covers the whole height of the
    jump. This images the sides of trees and roofs properly. Default 1; 0 = off; 2 for extra safety.
    Each extra frame width costs roughly 10–15% more frames with a typical tree line.
  - **Lowest altitude** sets a floor for the strip. Keep it well above the sun when you shoot by day
    (use **Sun check**).
- **Rotator:** keeps every frame level, so 10% overlap is enough.
  - Offset = sky angle − mechanical angle, as shown in NINA's rotator panel after a plate solve.
  - Being off by 180° doesn't matter.
- **Preview plan**, then **Save NINA sequence…**. Open the file in the Advanced Sequencer with
  **Load sequence**. Templates go into NINA's template folder.
  - **Show shooting order** (on by default) draws the order over the preview: a path from green
    (first frame) through yellow to magenta (last), with arrows, the first frame number of each row,
    Start/End markers and the meridian flip.
- **Conditions:** even overcast daylight gives the sharpest silhouettes. Take one test frame first;
  open sky should read about 10,000–30,000 ADU.

## 2. Build horizon

- Add the folder(s) with the frames. Subfolders are searched too.
- Leave **Calibrate rotation from overlaps** on. If you don't know the rotator offset at all, also tick
  **Full 360° search**.
- **Current horizon (compare / fill gaps):** filled in from the Plan tab, or chosen with **Browse…**. It
  is drawn in yellow, used for the usable-sky comparison, and used for any azimuths that weren't
  imaged.
- Click **Build mosaic & horizon**.
- **Save horizon file…** saves one NINA `.hrz` file. Choose the line next to the button:
  - **solid:** trees, roofs and walls
  - **with branches:** also avoids thin bare branches

  Edits and the side margin are included. Save twice if you want both files.
- Azimuths reported as **clear to bottom** or **blocked to top** reached the edge of the imaged area.
  **Plan fill-in frames…** turns them into a short NINA sequence. Shoot it, add its folder, and build
  again. Thanks to the cache, only the new frames are read.
  - **Fill-in: extra frames at steep edges** (below the button): every steep rise or drop of the new
    horizon line (steeper than **steep = more than**, default 3°/°, and at least 2° high) gets a
    **Steep edge** area. It spans the jump (±1°) and this many frame widths to each side of it. Only
    frames over gaps in the mosaic are shot, so edges that are already imaged cost nothing.
    Default 1, 0 = off. Both values are copied from the Plan tab's **Extra frames at steep edges**
    and **Steep = more than**, and follow changes made there; you can still change them here.
    Frames shared by overlapping areas are shot only once. CLI:
    `build … --fillin FILE.json --fillin-edge-frames N --fillin-edge-slope S`.
- **Show fill-in frames** (toolbar above the image) draws the fill-in areas (dashed magenta), their
  frames (orange) and the shooting order: a path from green (first) to magenta (last), with the first
  frame of each row numbered (every frame at zoom 20). It updates live as you edit the lines or change
  the fill-in settings, and switches on by itself when you click **Plan fill-in frames…**.
- **Re-detect** applies new thresholds instantly, without re-reading the frames.
- **Side margin at steep edges** (in the toolbar above the image): where the horizon rises or drops steeply (the side of a tree, the
  end of a roof), the higher value is extended this many degrees sideways. That stops NINA
  starting a target that is still right next to the obstruction. Gentle slopes stay unchanged.
  - "Steep" means more than the given degrees of altitude per degree of azimuth (default 3, about a
    72° slope).
  - 1–2° is a sensible margin; 0 switches it off.
  - The lines on screen and the saved files include the margin.
- **Editing the line:** above the image, choose **solid line** or **branches line**, then:
  - **left-drag** to redraw it (the brush width is in degrees),
  - **right-click** to restore the detected value,
  - **Ctrl+Z** to undo.

  Zoom (up to 20 px/°) and scroll to work precisely. Your edits are what gets saved. The two lines are
  independent: you can draw the branches line below the solid line, for example where the "solid" part
  is really a thin tree top you can see through.
- **Frame cache:** reduced frames are cached in `%LOCALAPPDATA%\HorizonMapper\cache` (about 200 kB per
  frame), so rebuilding takes seconds instead of minutes. **Clear cache** frees the space.
- Use the result in NINA: Options → General → Astrometry → Horizon.

## Command line (optional)

```
python horizon_mapper.py plan --hrz old.hrz --floor 20 --exposure 0.008 -o "Horizon mapping.json"
python horizon_mapper.py build "Y:\Photo Import\New Horizons" --old old.hrz -o horizon --fillin "Horizon fill-in.json" --side-margin 2
```

## Requirements of the frames

- NINA FITS files with `CENTALT`, `CENTAZ`, `XPIXSZ` and `FOCALLEN` in the header (NINA writes these
  when the mount is connected).
- `ROTATOR` and `PIERSIDE` are used when present. Colour (Bayer) and mono cameras both work.

## Changes

See [CHANGELOG.md](CHANGELOG.md).
