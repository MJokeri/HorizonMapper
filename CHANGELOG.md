# Changelog

## 1.0.1 – 2026-11-09

### Added
- **Plan tab – Extra frames at steep edges:** where the current horizon rises or drops steeply
  (threshold **Steep = more than**, default 3° altitude per degree of azimuth), the planned strip is
  widened sideways by the chosen number of frame widths (default 1, 0 = off) and spans the full
  height of the jump. CLI: `plan --edge-frames N --edge-slope S`. With the sample site horizon,
  0 / 1 / 2 extra frames give 647 / 733 / 819 frames.
- **Build tab – Fill-in: extra frames at steep edges:** the same idea for **Plan fill-in frames…**.
  Fields below the button set the number of extra frame widths (default 1, 0 = off) and the
  steepness threshold. A fill-in area whose side lies at a steep rise or drop of the measured
  (and edited) horizon is extended sideways, keeping its altitude range; it is labelled
  "+ steep edge". CLI: `build --fillin-edge-frames N --fillin-edge-slope S`.
- **Build tab – Show fill-in frames:** draws the fill-in areas, their frames and the shooting order
  on the horizon image, updating live.
- **Plan tab – Show shooting order:** the preview shows the order of frames as a colour-graded path
  (green = first → magenta = last) with row numbers, Start/End and the meridian flip.

### Fixed
- Fill-in frames shared by overlapping areas are now shot only once.

### Changed
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
