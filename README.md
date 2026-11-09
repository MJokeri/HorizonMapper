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
  - **Around current horizon:** under **Current horizon (.hrz)**, add the horizon file you use now. An
    optional second file (for example an all-sky estimate) widens the strip to cover both, plus the
    margins. The first file is also filled in automatically on the Build tab.
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
- **Re-detect** applies new thresholds instantly, without re-reading the frames.
- **Side margin at steep edges:** where the horizon rises or drops steeply (the side of a tree, the
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
