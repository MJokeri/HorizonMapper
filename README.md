# Horizon Mapper for NINA

The version is shown in the window title (currently 1.1.0).

Plans the frames that photograph the sky around your horizon, shoots them (through NINA's Advanced
API plugin, or by saving a NINA Advanced Sequencer sequence), then turns the FITS frames into a
measured NINA horizon file (`.hrz`).

Tabs: **1. Plan sequence** → **2. Capture (NINA)** → **3. Build horizon**.

On every tab the settings sit in a bar across the top and the picture fills the width below it.
The bar scrolls sideways (mouse wheel or scrollbar) when the window is too narrow; the Plan
buttons (and the log on the Build tab) on its right always stay visible. At start the window is made
wide enough for all groups if the screen allows. Closing the window asks first if the horizon has
not been saved since it was built or last edited.

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
    steeply than **Steep = more than (°/°)** (default 2° altitude per 1° azimuth, about 63°), the strip
    is widened sideways by this many frame widths on both sides and covers the whole height of the
    jump. This images the sides of trees and roofs properly. Default 0.5; 0 = off; 1–2 for extra safety.
    Fractions work too: with 0.5 a narrow edge strip is usually covered by two frames per row.
    Each extra frame width costs roughly 10–15% more frames with a typical tree line.
  - Frames are placed per row only over the stretches where the strip enters that row, using as
    few frames as the overlap allows. Every point of the strip is still covered.
  - **Lowest altitude** sets a floor for the strip. Keep it well above the sun when you shoot by day
    (use **Sun check**).
- **Order** (Sequence group): **East half first, then west** gives a single meridian flip. The east
  half is shot bottom to top and, by default, the west half top to bottom (the shortest path). Tick
  **Rows bottom to top in both halves** to start the west half at the bottom too.
  CLI: `plan --bottom-up`.
- **Rotator:** keeps every frame level, so 10% overlap is enough.
  - Offset = sky angle − mechanical angle, as shown in NINA's rotator panel after a plate solve.
  - Being off by 180° doesn't matter.
- **Preview plan**, then **Save NINA sequence…**. Open the file in the Advanced Sequencer with
  **Load sequence**. Templates go into NINA's template folder.
  - **Show shooting order** (on by default) draws the order over the preview: a path from green
    (first frame) through yellow to magenta (last), with arrows, the first frame number of each row,
    Start/End markers and the meridian flip.
  - **Mouse wheel** zooms the preview around the pointer, **drag** moves it, **double-click** shows
    the whole sky again.
- **Conditions:** even overcast daylight gives the sharpest silhouettes. Take one test frame first;
  open sky should read about 10,000–30,000 ADU.

## 2. Capture (NINA)

Shoots the frames directly from Horizon Mapper. NINA still talks to the equipment; Horizon Mapper
tells it what to do through the **Advanced API** plugin, one frame at a time. (Saving a sequence on
the Plan tab and running it in NINA still works as before.)

**Setup (once):** in NINA install **Advanced API** (Plugins tab), enable its server in the plugin
settings (default port 1888), and set NINA to save **FITS** files. Connect the mount, camera and
rotator in NINA as usual. Don't run a NINA sequence at the same time.

1. **Test connection** checks the plugin, shows which devices are connected, starts the live device
   line (mount alt/az, pier side, rotator, camera) and fills in NINA's image folder.
2. **Frames to shoot:** the plan from the Plan tab, the fill-in frames generated on the Build tab,
   or a **saved sequence file** (a plan or fill-in sequence saved for NINA; **Browse…**). Exposure,
   gain and binning are taken from the file. **Load frames** shows them.
   - The picture looks like the Plan preview: the strip to image, the current horizon, the lowest
     altitude, the frames and the shooting order (**Show shooting order**).
   - Each frame's outline shows its state: blue = to do, yellow (thick) = in progress, green = done,
     orange = skipped, red = error.
   - **When a frame is finished, its picture is drawn in its slot.** It is projected with the same
     maths as the mosaic, using the rotator offset (the Build tab's calibrated value if there is
     one). If the previews come out upside down, tick **Turn previews 180°**. **Show finished
     frames in the picture** turns them off.
   - Wheel zooms, drag moves, double-click shows the whole sky.
3. **Capture** (exposure, gain, binning) follows the Plan tab; offset is NINA's camera setting.
   **Settle** waits after each move.
4. **Safety:** frames closer to the Sun than **Keep away from the Sun** (default 30°) or below **Never
   below altitude** (default 10°) are skipped, never slewed to. **Unpark at start** / **Park when
   finished**.
5. **▶ Start** asks for confirmation, then for every frame: slew (alt/az converted to J2000 RA/Dec
   with the mount's own sidereal time) → rotator to its mechanical angle → settle → expose. NINA saves
   the FITS with its usual headers, and Horizon Mapper finds the file in NINA's image folder.

**Pause and errors**

- The big red **⏸ PAUSE** button, or **Esc** anywhere in the app, immediately stops the slew, the
  rotator and the exposure, and holds. **▶ Resume** repeats the interrupted frame from its start.
  While paused, **Skip frame** skips it.
- **Any error or timeout pauses the run automatically** (with a beep and the reason in the log).
  Examples: NINA not answering, a failed slew, the mount ending up somewhere else, or a stuck
  rotator. Nothing is skipped silently. Fix the cause and **Resume**, or **■ Stop run**.
- **■ Stop run** ends the run and keeps the frames taken so far. **Start** again offers to continue
  from where it stopped. Closing the window during a run asks first and stops everything.

When the run ends, the list of captured files is written to
`%LOCALAPPDATA%\HorizonMapper\captures\` and added to the Build tab's folder list. The Build tab
accepts such a list like a folder.

## 3. Build horizon

- Add the folder(s) with the frames. Subfolders are searched too.
- Leave **Calibrate rotation from overlaps** on. If you don't know the rotator offset at all, also tick
  **Full 360° search**.
- **Current horizon (compare / fill gaps):** filled in from the Plan tab, or chosen with **Browse…**. It
  is drawn in yellow, used for the usable-sky comparison, and used for any azimuths that weren't
  imaged.
- Click **Build mosaic & horizon**.
- **Save horizon file…** (below the fill-in buttons) saves one NINA `.hrz` file. Choose the line next
  to the button:
  - **solid:** trees, roofs and walls
  - **with branches:** also avoids thin bare branches

  Edits and the side margin are included. Save twice if you want both files.
- **Fill-in frames** cover the places where the solid or branches line (including your edits) reaches
  the edge of the imaged area (**clear to bottom** or **blocked to top**), plus the steep edges below.
  Dark frame edges high in the sky that the lines already rejected are ignored.
  - **Generate fill-in frames** plans them, draws them on the image and writes a summary to the log.
  - **Save fill-in sequence…** saves exactly the plan that is shown as a short NINA sequence.
  Shoot it, add its folder, and build again. Thanks to the cache, only the new frames are read.
  - **Fill-in: extra frames at steep edges** (above the buttons): every steep rise or drop of the new
    horizon line (steeper than **steep = more than**, default 2°/°, and at least 2° high) gets a
    **Steep edge** area. It spans the jump (±1°) and this many frame widths to each side of it. Only
    frames over gaps in the mosaic are shot, so edges that are already imaged cost nothing.
    Default 0.5, 0 = off. Both values are copied from the Plan tab's **Extra frames at steep edges**
    and **Steep = more than**, and follow changes made there; you can still change them here.
    Frames shared by overlapping areas are shot only once. CLI:
    `build … --fillin FILE.json --fillin-edge-frames N --fillin-edge-slope S`.
  - **Fill-in: horizon line near gaps within (°)**: wherever the solid or branches line runs through
    a spot without pictures, or comes within this many degrees of one, a **Gap** area is added there,
    whether or not the horizon is steep. Only frames over the gap are shot. Default 1°, 0 = off.
    CLI: `--fillin-gap-near D`.
- **Show fill-in frames** (toolbar above the image) draws the fill-in areas (dashed magenta), their
  frames (orange) and the shooting order: a path from green (first) to magenta (last), with the first
  frame of each row numbered (every frame at zoom 20). After you edit the lines or change the
  fill-in settings, the plan is generated again and redrawn, so the image always shows what
  **Save fill-in sequence…** will write.
- **Re-detect** applies new thresholds instantly, without re-reading the frames.
- **Side margin at steep edges** (in the toolbar above the image): where the horizon rises or drops steeply (the side of a tree, the
  end of a roof), the higher value is extended this many degrees sideways. That stops NINA
  starting a target that is still right next to the obstruction. Gentle slopes stay unchanged.
  - "Steep" means more than the given degrees of altitude per degree of azimuth (default 2, about a
    63° slope).
  - 1–2° is a sensible margin; 0 switches it off.
  - The lines on screen and the saved files include the margin.
- **Editing the line:** above the image, choose **solid line** or **branches line**, then:
  - **left-drag** to redraw it (the brush width is in degrees),
  - **right-click** to restore the detected value,
  - **Ctrl+Z** to undo.

  **Mouse wheel** zooms the image around the pointer (2–24 px/°; the **zoom** box offers fixed
  levels); **Shift+wheel**, **middle-drag** or the scrollbars move it. Your edits are what gets saved. The two lines are
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
