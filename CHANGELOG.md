# Changelog

## 1.1.0 – 2026-11-10

### Added
- **2. Capture (NINA) tab:** shoots the Plan tab plan or the Build tab fill-in frames directly,
  through NINA's **Advanced API** plugin. For each frame: slew (alt/az → J2000 RA/Dec from the mount's
  sidereal time), mechanical rotator move, settle, exposure. NINA saves the FITS with its usual
  headers; the files are found in NINA's image folder (read from the profile). The file list is added
  to the Build tab, which now accepts such list files like folders.
- **Pause:** a large red **⏸ PAUSE** button and **Esc** stop slew, rotator and exposure immediately and
  hold; **Resume** repeats the interrupted frame, **Skip frame** skips it, **End run** ends the run
  (a later Start can continue). Errors and timeouts pause the run automatically.
- **Safety limits:** Sun avoidance (default 30°) and a minimum altitude (default 10°); frames outside
  them are skipped, never slewed to. Confirmation before starting and when closing during a run.
- Live device line (mount position, pier side, rotator, camera state).
- **Capture picture like the Plan preview** (strip, current horizon, frames, shooting order) with each
  frame's state; **finished frames are drawn in their slots**, projected like the mosaic (calibrated
  rotator offset from the Build tab if available; **Turn previews 180°** option).
- **Frames from a saved sequence file:** besides the Plan tab plan and the Build tab fill-in, any plan
  or fill-in sequence saved for NINA can be loaded (exposure, gain and binning are read from it).
- **Plan tab – Rows bottom to top in both halves:** with the east/west split, the west half can also
  be shot bottom to top (default stays top to bottom, the shortest path). CLI: `plan --bottom-up`.

- **Save current capture… / Load saved capture…** (Capture tab, replacing Load frames):
  - A `.hmcap.json` file stores the frames, each frame's state, the file path and NINA file name of
    every frame taken, and the capture, safety, range and file settings.
  - Loading restores everything and shows the pictures from the original paths. A frame whose file
    is gone is set back to "to do".
  - The progress is also saved automatically after every frame (`captures\last_capture.hmcap.json`).
- **Range: From frame / To frame** (Capture tab): a run shoots only those frames (in shooting order),
  leaving out frames already taken unless **Re-shoot frames already taken** is ticked. A stopped or
  loaded capture is resumed simply by pressing Start.
- Choosing a frame source shows its frames at once; a fresh Start uses the source as it is then.
- **Hover explanations** for every setting, choice and button on all tabs. New controls get theirs
  from the `TIPS` table; the debug log lists any control without one.
- **Debug → Show log windows** (on during development): shows or hides the Capture and Build logs.
  The Capture log sits in its own pane; drag the divider to resize it.
- `DEV_MODE` switch for the development defaults (debug log and log windows on), plus a release
  checklist in the README.
- Capture tab: default host 192.168.1.55 (for testing).
- Capture tab: the log runs down the right side at full height, with a scrollbar.
- **Debug menu → Debug log (all actions):** every click, NINA request (timing and answer), capture
  step, wait, pause and state change, folder scan and redraw, plus errors from window callbacks and
  background threads, goes to `%LOCALAPPDATA%\HorizonMapper\logs\`. A watchdog logs whenever the
  window is blocked for more than a second. The setting is remembered and on by default for now.

### Fixed
- **The app crashed at start** ("super(type, obj): obj must be an instance or subtype of type"). The
  hover and click-logging wrappers replaced classes inside the tkinter library itself, and ttk
  builds its Combobox and Spinbox on ttk.Entry internally. The wrappers now live in the app's own
  copy of the tkinter names, and the library stays untouched.
- **The Windows build now opens the real window** (`--selftest-gui`), both as a Python script and as
  the exe. It switches through the tabs, shows a hover and closes, so a start-up crash fails the
  build instead of reaching users.
- **The window could freeze during a capture run, making PAUSE and Stop unresponsive.**
  - Every event from the run redrew the whole picture, and the handler only gave control back to
    the window once all queued events were done. With hundreds of frames and previews, redrawing
    became slower than the events arrived.
  - Now events are handled for at most ~80 ms at a time, the picture is redrawn at most four times
    a second, and resized previews are cached.
  - Background threads no longer read window settings, which on Windows can stall both sides.
- **Load frames during a paused run locked the app.** It recomputed the plan on the window's thread
  and dropped the paused run, so Resume, Stop and PAUSE had nothing to act on.
  - Frames are now prepared in the background.
  - Loading while paused first asks, then ends the paused run properly (its frames are kept and
    listed for the Build tab).
  - While a run is in progress, Load frames is disabled.
- **Some frames never got a preview:** the file was found while NINA was still writing it (empty or
  short), the preview failed and was not tried again. Now a file is used only when its size covers
  the header and all the image data. Incomplete files, and failed previews, are retried every 3 s
  for up to 3 minutes.
- Loading frames can no longer get stuck: the internal "loading" flag is always cleared, a stale one
  is reset after 60 s, and every decision is written to the debug log.
- **Changing the NINA image folder during a run had no effect** (the folder was read once at Start),
  so no files were found and no pictures appeared. The folder is now read whenever files are looked
  up, and the debug log shows the folder, whether it exists, and NINA's file names.
- Files that NINA saves late no longer hold up or pause the run. The app scans NINA's image folder
  every few seconds and matches new FITS files to their frames, by name or by the pointing in the
  header. They are added to the picture and the file list as they arrive, also for 10 minutes
  after the run. Waiting frames are shown in light green.
- Capture would not start with mounts whose driver reports switching tracking on as failed (seen
  with the ZWO "ASI Mount": NINA answered "Tracking mode changed" with failure). Tracking is now
  tried, checked and only noted in the log if it stays off; every slew is still verified.
- Slews are no longer judged too early: the app now waits for NINA's own slew to finish (the
  request waits, while PAUSE still works at once), then until the mount stands still, before
  checking where it is. Before, a mount slow to report "slewing" (ZWO ASI Mount) was taken as
  arrived at the old position, and the automatic pause then stopped the slew it had just begun.
- Slew arrival is checked against the mount's RA/Dec (recomputed for the current sidereal time),
  so it works whether or not the mount tracks.
- Error messages from NINA that carry no error text now say that NINA reported a failure.

### Changed
- **Stop run** is now **End run**.
- **Tracking is no longer switched on at the start of a capture.** The exposures are very short and
  every slew is aimed at the current sky position, so tracking isn't needed. This avoids a confusing
  "failed" reply from drivers such as the ZWO ASI Mount, which answer the tracking command that way
  even when it works. If a slew is refused while the mount isn't tracking, tracking is switched on
  and that slew is tried once more. A non-tracking mount gets one plain information line in the log.
- Tabs renumbered: 1. Plan sequence, 2. Capture (NINA), 3. Build horizon.

## 1.0.1 – 2026-11-09

### Added
- **Plan tab – Extra frames at steep edges:** where the current horizon rises or drops steeply
  (threshold **Steep = more than**, default 2° altitude per degree of azimuth), the planned strip is
  widened sideways by the chosen number of frame widths (default 0.5, 0 = off) and spans the full
  height of the jump. CLI: `plan --edge-frames N --edge-slope S`. With the Tikankuja horizon,
  0 / 0.5 / 1 / 2 extra frames give 569 / 638 / 713 / 873 frames (steep = 2°/°).
- **Plan tab – fewer frames:** each row now places its frames only over the stretches where the
  strip enters that row's altitude zone, as few as the overlap allows. Before, frames sat on a fixed
  grid and every grid frame that touched the strip was kept, so a narrow strip (for example with
  0.5 extra frames at a steep edge) often got three frames where two cover it. Horizon plans need
  10–13% fewer frames (Tikankuja: 647 → 569 without extra edge frames); full-circle band plans are
  unchanged. Coverage was checked on a 0.2° grid: no gaps.
- **Build tab – Generate fill-in frames / Save fill-in sequence…:** generating and saving are
  separate buttons. Generate draws the plan and logs a summary; Save writes exactly the plan shown.
- **Build tab – Fill-in: extra frames at steep edges:** the same idea for the fill-in frames.
  Every steep rise or drop of the new horizon line (at least 2° high) gets a **Steep edge** area
  spanning the jump and the chosen number of frame widths to each side; only frames over gaps in
  the mosaic are shot. The fields below the button take their values from the Plan tab and follow
  changes made there. CLI: `build --fillin-edge-frames N --fillin-edge-slope S`.
- **Build tab – Fill-in: horizon line near gaps within (°):** wherever the solid or branches line
  runs through a spot without pictures, or comes within this distance of one (default 1°), a **Gap**
  area is added, whether or not the horizon is steep there. Only frames over the gap are shot.
  CLI: `build --fillin-gap-near D`.
- **Build tab – Show fill-in frames:** draws the fill-in areas, their frames and the shooting order
  on the horizon image, updating live.
- **Plan tab – Show shooting order:** the preview shows the order of frames as a colour-graded path
  (green = first → magenta = last) with row numbers, Start/End and the meridian flip.

- **Mouse zoom:** the wheel zooms the Plan preview (drag moves it, double-click shows the whole sky)
  and the Build tab image (around the pointer, 2–24 px/°; Shift+wheel or middle-drag moves it)
  instead of scrolling sideways.
- **Settings bar at the top:** on both tabs the settings are grouped in a horizontal bar across the
  top and the picture fills the full width below it (it is always wider than high). The bar scrolls
  sideways when the window is too narrow; the action buttons (and the Build tab's log) on its right
  stay visible. Long groups were split into short columns; the default window is 1600×900.

- **Exit check:** closing the window asks for confirmation if the horizon has not been saved since
  it was built or last edited.
- **Window size:** at start the window is made wide enough for all settings groups, as far as the
  screen allows.

### Fixed
- Plan preview: a frame lying across 0°/360° was drawn as a box across the whole preview.
- Fill-in areas are checked against the solid and branches lines (with edits). Before, the detector's
  raw "blocked to top" flags also fired on dark frame edges high in the sky that the lines had
  rejected, which planned many useless frames at the top of the imaged area.
- The fill-in overlay now always follows the settings and line edits: the plan is generated again
  and redrawn after every change, and errors are written to the log instead of failing silently.
- Fill-in frames shared by overlapping areas are now shot only once.

### Changed
- Defaults: **Extra frames at steep edges** 0.5 (was 1) and **Steep = more than** 2°/° (was 3), on the
  Plan tab, for fill-in and in the CLI; the Build toolbar's side-margin steepness is also 2.
- **Build tab:** the **Build** group (Build mosaic & horizon, Re-detect, progress) sits in the settings
  bar between **Current horizon** and **Fill-in frames**; only the log stays fixed on the right.
- **Build tab:** **Save horizon file…** and its solid / with branches choice moved below the fill-in
  buttons, so the order follows the workflow: build, fill in, save.
- **Side margin at steep edges:** the setting and its steepness threshold moved from the Build tab's
  Settings box to the toolbar directly above the image, as spin boxes. The lines on screen update as
  you change them. On smaller screens the old place could be cut off and hard to find.
- The window title shows the version number, so you can tell which build you are running.

## 1.0.0 – 2026-11-09

### Added
- **Plan fill-in frames:** finds the azimuths whose horizon reached the edge of the imaged area
  ("clear to bottom" or "blocked to top"), groups neighbouring ones into areas, and writes a short
  NINA sequence covering them. Available in the GUI and as `build --fillin FILE.json`.
- **Horizon line editor** in the Build tab:
  - left-drag redraws the solid or branches line, with an adjustable brush width;
  - right-click restores the detected value, and Ctrl+Z undoes;
  - zoom up to 20 px/° and scroll.

  Saved files use the edited lines.
- **Re-detect** applies new detection thresholds without re-reading the frames.
- **Side margin at steep edges:** where the horizon rises or drops more steeply than a set slope
  (default 3° of altitude per degree of azimuth), the higher value is extended a chosen number of
  degrees sideways. Gentle slopes are unchanged. It's in the GUI Settings, and on the command line as
  `--side-margin` and `--steep-slope`. The value is recorded in the details JSON. Off by default.
- **Frame cache:** reduced frames are stored in `%LOCALAPPDATA%\HorizonMapper\cache` and reused, and
  frames are read in parallel. On a 129-frame set a rebuild took 6 s instead of 59 s. `--no-cache`
  and **Clear cache** control it.
- **Windows .exe:** a GitHub Actions workflow builds `HorizonMapper.exe` with PyInstaller,
  smoke-tests it, uploads it as an artifact and attaches it to `v*` releases. `build_exe.bat`
  builds it locally.

### Fixed
- Dragging outside the image while editing no longer sets impossible altitudes. Edits are clamped to
  the altitude range shown.

### Changed
- **Saving:** **Save horizon file…** now writes a single `.hrz` file, with a file dialog. A selector
  next to the button picks the solid line or the line with branches. The GUI no longer writes the
  mosaic PNG and details JSON. The command-line `build` still writes all outputs.
- **Current horizon:** the Plan tab's horizon field is now labelled **Current horizon (.hrz)**. The
  first file added there is filled in automatically as the Build tab's comparison and gap-filling
  horizon. Clearing it on the Plan tab clears it there too. The Build tab box is renamed **Current
  horizon (compare / fill gaps)**.
- When you edit, the branches line may now go below the solid line. The lines are no longer forced
  together on screen or when saving. Detection still starts the branches line at or above the solid
  line.
- The build workflow uses Node.js 24 action versions (`checkout` v6, `setup-python` v6,
  `upload-artifact` v6, `action-gh-release` v3).
- A windowed `.exe` started with command-line arguments no longer fails on the missing console.

## 0.1.0 – 2026-11-09
- First version: sequence planning (fixed band or around horizon files, rotator levelling,
  east-half-first ordering, NINA sequence/template output), mosaic building with automatic rotation
  calibration and brightness matching, two-level horizon detection, `.hrz` output, comparison with an
  old horizon, and a command-line mode.
