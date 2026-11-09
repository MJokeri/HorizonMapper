# Horizon Mapper for NINA

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
  - **Around horizon file(s):** add your current `.hrz`, or two files (for example your old horizon and a
    rough estimate). The strip covers both, plus the margins.
  - **Lowest altitude** sets a floor for the strip. Keep it well above the sun when you shoot by day
    (use **Sun check**).
- **Rotator:** keeps every frame level, so 10% overlap is enough.
  - Offset = sky angle − mechanical angle, as shown in NINA's rotator panel after a plate solve.
  - Being off by 180° doesn't matter.
- **Preview plan**, then **Save NINA sequence…**. Open the file in the Advanced Sequencer with
  **Load sequence**. Templates go into NINA's template folder.
- **Conditions:** even overcast daylight gives the sharpest silhouettes. Take one test frame first;
  open sky should read about 10,000–30,000 ADU.

## 2. Build horizon

- Add the folder(s) with the frames. Subfolders are searched too.
- Leave **Calibrate rotation from overlaps** on. If you don't know the rotator offset at all, also tick
  **Full 360° search**.
- Add your old `.hrz` under **Compare**. It is used for the comparison, and for any azimuths that
  weren't imaged.
- Click **Build mosaic & horizon**, then **Save results…**. You get:
  - `horizon_DATE.hrz`: solid obstructions (trees, roofs, walls)
  - `horizon_DATE_with_branches.hrz`: also avoids thin bare branches
  - `horizon_DATE_mosaic.png`: the 360° mosaic with the horizons drawn on it
  - `horizon_DATE_details.json`: per-azimuth details
- Azimuths reported as **clear to bottom** or **blocked to top** reached the edge of the imaged area.
  **Plan fill-in frames…** turns them into a short NINA sequence. Shoot it, add its folder, and build
  again. Thanks to the cache, only the new frames are read.
- **Re-detect** applies new thresholds instantly, without re-reading the frames.
- **Editing the line:** above the image, choose **solid line** or **branches line**, then:
  - **left-drag** to redraw it (the brush width is in degrees),
  - **right-click** to restore the detected value,
  - **Ctrl+Z** to undo.

  Zoom (up to 20 px/°) and scroll to work precisely. Your edits are what gets saved.
- **Frame cache:** reduced frames are cached in `%LOCALAPPDATA%\HorizonMapper\cache` (about 200 kB per
  frame), so rebuilding takes seconds instead of minutes. **Clear cache** frees the space.
- Use the result in NINA: Options → General → Astrometry → Horizon.

## Command line (optional)

```
python horizon_mapper.py plan --hrz old.hrz --floor 20 --exposure 0.008 -o "Horizon mapping.json"
python horizon_mapper.py build "Y:\Photo Import\New Horizons" --old old.hrz -o horizon --fillin "Horizon fill-in.json"
```

## Requirements of the frames

- NINA FITS files with `CENTALT`, `CENTAZ`, `XPIXSZ` and `FOCALLEN` in the header (NINA writes these
  when the mount is connected).
- `ROTATOR` and `PIERSIDE` are used when present. Colour (Bayer) and mono cameras both work.

## Changes

See [CHANGELOG.md](CHANGELOG.md).
