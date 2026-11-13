# Horizon Mapper for NINA

The version is shown in the window title (currently 1.1.0).

Plans the frames that photograph the sky around your horizon, shoots them (through NINA's Advanced
API plugin, directly through ASCOM Alpaca without NINA, or by saving a NINA Advanced Sequencer
sequence), then turns the FITS frames into a measured NINA horizon file (`.hrz`) and, if you wish, a **Stellarium** landscape.

Tabs: **1. Plan sequence** → **2. Capture** → **3. Build horizon**.

**Capture works two ways:** through **NINA** (Advanced API plugin), or through **ASCOM Alpaca without
NINA** (e.g. ASCOM Remote), also from another computer: Horizon Mapper on your own PC, the devices
on the observatory PC, and the images saved on your PC. The Alpaca method has been run on real
equipment (ZWO ASI mount, ASI2600MC, ZWO CAA rotator through ASCOM Remote). See
[2. Capture](#2-capture) for the setup, and note that Alpaca has no password: use it only on your
own network or a VPN.

> ⚠ **Experimental software – use entirely at your own risk.** It moves real equipment and
> may point it near the **Sun**, which can instantly blind and destroy equipment. Read the
> [Disclaimer and safety warning](#disclaimer-and-safety-warning) before use.

**What is new in 1.1.0** (details in [CHANGELOG.md](CHANGELOG.md)):
- **Capture tab:** shoots the frames itself, through NINA or through **ASCOM Alpaca without NINA**
  (also from another computer, images saved on your own computer).
- **Safety and Peripherals:** rotator range (full, 180°, 90°), sun and altitude limits, and **flat panel**
  open/close buttons with options to open it at the start of a run and close it at the end.
- **Run control:** pause/resume, **From/To frame** with a question when frames are already taken,
  **Optimise route** (fewer meridian crossings, shorter run), camera settings changeable during a run,
  saved captures named `_completed` or `_partial`.
- **Build tab:** **Save Stellarium landscape…** next to the NINA horizon file.

**Hover** over any setting, choice or button for a short explanation of what it does.

On every tab the settings sit in a bar across the top and the picture fills the width below it.
The bar scrolls sideways (mouse wheel or scrollbar) when the window is too narrow; the Plan
buttons (and the log on the Build tab) on its right always stay visible. At start the window is made
wide enough for all groups if the screen allows. Closing the window asks first if the horizon has
not been saved since it was built or last edited.

## Disclaimer and safety warning

**Please read this before you use the software. By installing, running or distributing Horizon
Mapper for NINA you accept the terms below.**

### 1. Experimental software – no warranty, no liability

Horizon Mapper for NINA is **experimental software**. It sends movement and exposure commands to
real equipment (mount, rotator, camera) and may contain errors that cause that equipment to move or
behave unexpectedly.

- **No warranty.** The software is provided **"AS IS"**, without warranty of any kind, express or
  implied, including but not limited to the implied warranties of merchantability, fitness for a
  particular purpose and non-infringement. The entire risk as to the quality, performance and
  results of the software is with you. Should the software prove defective, you assume the cost of
  all necessary servicing, repair or correction (GNU GPL v3, section 15).
- **No liability.** To the maximum extent permitted by applicable law, in no event shall the
  developer, the copyright holders, the distributors or any other party who modifies and/or conveys
  the software be liable to you or to any third party for any damages or losses of any kind –
  including, without limitation, damage to or loss of equipment, property or data, personal injury
  (including injury to the eyes), or any direct, indirect, special, incidental or consequential
  damages – arising out of or in connection with the use of, or the inability to use, the software,
  even if such party has been advised of the possibility of such damages (GNU GPL v3, section 16).
- **Use at your own risk and responsibility.** You alone are responsible for deciding whether and
  how to use the software, for checking every setting and plan before running it, for the safe set-up,
  balance, cabling and limits of your equipment, and for supervising the equipment while it is
  controlled by the software. Keep people, animals and objects clear of moving equipment and be ready
  to stop it at any time (the **⏸ PAUSE** button / **Esc**, and the controls of NINA, the mount and its
  power supply).
- **Applicable law.** Where applicable law does not allow a limitation or exclusion of liability
  (for example for intent, gross negligence or injury to life, body or health under mandatory consumer
  protection law), the limitation applies to the greatest extent that law permits; the remaining
  provisions stay in effect. As provided in GNU GPL v3, section 17, this disclaimer is to be
  interpreted in the way that most closely approximates an absolute waiver of all civil liability.
- This section summarises and supplements the terms of the [GNU General Public License v3](LICENSE),
  under which the software is licensed. In case of conflict, the licence prevails.

### 2. ⚠ Danger: the Sun

Horizon frames are often taken **in daylight**, so the Sun may be in the sky while the mount moves
and the camera exposes.

- **Pointing a telescope, finder, guide scope or camera lens at the Sun – or even slewing across or
  close to it – can instantly cause severe and permanent eye injury or blindness** to anyone looking
  through the optics, and can **damage or destroy cameras, sensors, filters, mirrors, focusers,
  guide cameras and other equipment** through concentrated light and heat; it can also create a
  **risk of fire**. NASA states that looking at even part of the bright Sun through a telescope,
  binoculars or camera lens that has no purpose-made solar filter fixed over its front end injures
  the eye severely and at once. The American Academy of Ophthalmology warns that
  viewing the Sun even briefly can permanently damage the retina, and that concentrated solar rays
  can melt or damage filters.
- **Never look through any optics during a daytime capture.** Cover or remove finder scopes and
  guide scopes that are not needed.
- **Only properly certified solar filters securely attached to the *front* of the optics** keep
  the Sun's light and heat out of the instrument (American Astronomical Society). Eyepiece filters,
  eclipse glasses, sunglasses, ND filters or improvised materials are **not** safe in or behind
  optics. Unless you deliberately do solar imaging with such filters fitted and checked, the optics
  must **never** point at or near the Sun.
- **The Sun avoidance in this software is an aid only and is no guarantee of safety.** The Sun
  position, frame timing and slew paths are estimates; the real path of a mount (meridian flips,
  firmware behaviour, axis order, clock, location, time zone or settings errors, communication
  failures or interrupted commands) can differ, and a mount can pass near the Sun even when no frame
  is planned there. Keep **Keep away from the Sun** generous, check the Sun overlay and the time and
  location settings before every run, prefer capturing when the Sun is low, behind an obstruction or
  below the horizon, and **stay with the equipment and supervise every daytime run**, ready to stop
  it immediately.
- If you are unsure, do not capture in daylight. Seek advice from an experienced astronomer before
  using any solar filter with a camera or telescope.

### Sources

- GNU General Public License v3, sections 15–17 (this repository's [LICENSE](LICENSE));
  <https://www.gnu.org/licenses/gpl-3.0.html>
- NASA – Eclipse eye safety: <https://science.nasa.gov/eclipses/safety/> and
  <https://science.nasa.gov/eclipses/future-eclipses/eclipse-2024/safety/>
- American Astronomical Society – Solar filters for optics (telescopes, binoculars & cameras):
  <https://eclipse.aas.org/eye-safety/optics-filters>
- American Academy of Ophthalmology – Photographing an eclipse without damaging your eyes:
  <https://www.aao.org/newsroom/news-releases/detail/learn-how-to-photograph-eclipse-without-damaging-y>

## Debug log

**Debug → Debug log (all actions)** writes everything the app does to a log file in
`%LOCALAPPDATA%\HorizonMapper\logs\`, with millisecond times and the thread that did it:
- every click and Esc press;
- the start-up time: lines `startup: N s – …` show how long each part of the start-up took (useful if the
  program starts slowly);
- every NINA request with its timing and answer;
- every step, wait, pause and state change of a capture run;
- folder scans and redraw times;
- errors, including ones the window would otherwise hide.

A watchdog adds a line whenever the window is blocked for more than a second.
**Debug → Show log windows** shows or hides the log panes. The Capture log sits beside the
picture, and you can drag the divider to make it wider or narrower. **Debug → Open log
folder** opens the folder. The setting is remembered; it is on by default during development.
Please attach the log when reporting a problem.

## Development notes

- Every new setting, choice or button needs a hover text: add an entry to `TIPS` in
  `horizon_mapper.py` (key = the text shown in the window). The debug log lists controls without
  one ("controls without hover text").
- Update the README, the in-app Help and `CHANGELOG.md` with every change.
- `python horizon_mapper.py --selftest-gui` opens the real window, switches through the tabs, shows a
  hover and closes. The Windows build runs it for the script and for the exe.
- Never patch the `tkinter` / `ttk` modules themselves: the wrappers in `run_gui` work on private
  copies of their names, because ttk builds some widgets on others internally.

### Release checklist

1. Set `DEV_MODE = False` in `horizon_mapper.py`. This switches the development defaults off:
   - the debug log;
   - the log windows.

   Users can still change all three.
2. Set `APP_VERSION` and finish the `CHANGELOG.md` entry.
3. Tag the release (`vX.Y.Z`); the GitHub Actions build attaches the exe to the release.

## Install (Windows)

**Easiest:** download `HorizonMapper-windows.zip` from the repository's
[Releases](../../releases) page. If there's no release yet, open the latest
[Actions → Build Windows exe](../../actions) run and download the build from there.
**Unzip it to a normal folder first** (right-click → Extract All) and run `HorizonMapper.exe` from there. Running it
directly from inside the ZIP makes Windows unpack it to a temporary folder and scan it, and the start can take
about 10 seconds. No Python is needed.

**If it starts slowly** (the single-file `.exe` unpacks itself at every start, and virus scanners can make that
take 10 seconds or more): use the **folder version**, `HorizonMapper-windows-folder.zip` (artifact
`HorizonMapper-windows-folder`). Unzip it anywhere and run `HorizonMapper\HorizonMapper-folder.exe`; it does
not unpack anything at start.

**From source:** install Python 3.10 or newer from https://www.python.org (tick **"Add python.exe to PATH"**).
Then double-click **`Start Horizon Mapper.bat`**. The first time, it installs `numpy` and `Pillow`.
`build_exe.bat` builds the `.exe` yourself.

## 1. Plan sequence

- **Camera & telescope:** use **Read from a FITS file…** with any frame from the setup. It fills
  in sensor, pixel size, focal length and latitude.
- **Site (latitude and longitude):** the latitude is in **Site & area to image**, the longitude next to
  **Sun check**. Both are **remembered in `settings.json`** (`site_lat`, `site_lon`) the moment you change
  them (also when read from a FITS file), and are the defaults the next time the program starts. The
  command line `plan --lat` uses the saved latitude as its default. `settings.json` is in
  `%LOCALAPPDATA%\HorizonMapper`. The first time, a neutral default (60°, 25°) is shown and saved: enter
  your own site once.
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
  - **Lowest altitude** sets a floor for the strip. Keep it well above the sun when you shoot by day.
- **The Sun in the preview:** enter **Start date** and **Start time (HH:MM)**, plus **Sun distance
  (°)** (default 30), under the Plan buttons. The date can be `2026-11-12`, `12.11.2026` or
  `12.11.` (this year); blank means today, and a blank time means the current time of day; the longitude is the one next to **Sun check**. The preview
  then draws:
  - the Sun at the start (filled) and at the plan's estimated end (ring), with its track between;
  - the keep-away zone around both;
  - **frames that would be too close to the Sun when their turn comes, in red** (each frame's time
    comes from the shooting order and the time estimates);
  - a **red ×** on frames whose slew from the previous frame would pass near the Sun.

  A line at the top sums it up. Try different start times or azimuth ranges until nothing is red.
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

## 2. Capture

Shoots the frames directly from Horizon Mapper, one frame at a time. Choose under **Connect
through** how the equipment is reached; everything else on the tab works the same for both:

| | **NINA** | **Alpaca (without NINA)** |
|---|---|---|
| Who drives the devices | NINA, told what to do through its **Advanced API** plugin | Horizon Mapper itself, through an **ASCOM Alpaca** server |
| What must run | NINA with the plugin | An Alpaca server, e.g. **ASCOM Remote** (part of ASCOM Platform 7), or devices that speak Alpaca themselves |
| Where the images are saved | NINA's image folder, on NINA's computer | **On the computer running Horizon Mapper**, in the folder you choose |
| Remote use | Horizon Mapper here, NINA on the observatory computer | Horizon Mapper here, ASCOM Remote on the observatory computer |

(Saving a sequence on the Plan tab and running it in NINA still works as before.)

**Setup for NINA (once):** in NINA install **Advanced API** (Plugins tab), enable its server in the
plugin settings (default port 1888), and set NINA to save **FITS** files. Connect the mount, camera
and rotator in NINA as usual. Don't run a NINA sequence at the same time.

**Setup for Alpaca (once), on the computer the equipment is connected to:**

1. Install **ASCOM Platform 7** (it includes **ASCOM Remote**) and the ASCOM drivers of your mount,
   camera and rotator. The camera needs its **ASCOM** driver: ASCOM Remote cannot use the native
   camera driver NINA may be using.
2. Start **ASCOM Remote Server**, press **Setup** and add the mount (Telescope), camera and rotator
   on the device tab.
   - **Server address:** on the *Server Configuration* tab, ASCOM Remote starts listening on
     `localhost (127.0.0.1)` only, which other computers cannot reach. Choose this computer's
     **network address** (or all addresses) from the drop-down list, accept the administrator
     prompt, and restart it. Without this, Horizon Mapper gets "Bad Request" (Windows' own web
     service refusing the address) or "connection refused".
   - **Firewall (Windows, on the observatory computer):** allow inbound **TCP 11111** (the Alpaca
     port, used for everything) and **UDP 32227** (only for **Find**). The discovery port is UDP,
     not TCP. If there was no firewall prompt when ASCOM Remote first started, add the rules by hand.
   - Test from the other computer: open `http://<address>:11111/management/v1/description` in a
     browser. It should show a short text starting with `{"Value"`.
3. Close NINA, or at least disconnect the devices there, if the drivers cannot be shared.
4. Alpaca has **no password and no encryption**: use it only inside your own network, or through a
   VPN – never open the port to the internet.

In Horizon Mapper, with **Alpaca (without NINA)** chosen:

- **Alpaca host** / **Alpaca port**: where the Alpaca server runs (`localhost` if on this computer). **Find** looks for Alpaca servers on the local network and
  fills these in; it can miss a server (it needs UDP 32227 open), so entering the address by hand
  always works. Problems are explained in the Capture log, and the debug log has the details.
- **Test connection** lists the server's devices in **Mount**, **Camera** and **Rotator** (choose
  `none` if there is no rotator). If the server offers several devices of one kind (e.g. a main and
  a guide camera), choose the right ones and press **Test connection** again; the chosen devices
  are then connected and the live device line starts. Start also connects them.
- **Save images in (on this computer)**: the frames are downloaded from the camera and saved here
  as FITS files, in a subfolder per date, named
  `<target>_<date>_<time>_f<frame number>_alt.._az...fits`. The default is
  `Documents\Horizon Mapper\Captures`. The files carry the headers the Build tab reads, the same
  ones NINA writes (`CENTALT`, `CENTAZ`, `ROTATOR`, `PIERSIDE`, `XPIXSZ`, `FOCALLEN`, `SITELAT`, …).
  Pixel size comes from the camera and **focal length from the Plan tab**, so fill it in there.
  Because the files are written by this program, they are complete the moment a frame is done:
  the pictures appear at once and nothing has to be looked for afterwards.
- **Gain**: -1 leaves the camera's gain as it is; otherwise it is set in the camera driver.
  Binning is set in the camera, and the whole sensor is read.
- The method, address, devices and folder are remembered for the next time.
- **Rotator drivers differ.** Some do something else with the standard "move to mechanical angle"
  call (the ZWO CAA driver moves *by* that angle instead). Horizon Mapper starts with MoveAbsolute
  (corrected by the driver's offset) and checks where the rotator ends up; if it is wrong it says so
  in the log and tries the next method (MoveMechanical, then a relative Move), and remembers the one
  that works for the next run. If none reaches the angle, the run pauses with the numbers.
- The live device line refreshes every 2 s during a run and every 5 s otherwise (as long as the tab
  is connected), so the Alpaca server's own request log keeps filling slowly while the app is open.

With **NINA** chosen:

1. **Host** / **Port**: where NINA runs (port 1888 by default). The host is `localhost` or the address of the NINA computer. **Test connection** (in the same group) checks the plugin, shows which devices are connected, starts the live device
   line (mount alt/az, pier side, rotator, camera) and fills in NINA's image folder.
2. **Frames to shoot:** the plan from the Plan tab, the fill-in frames generated on the Build tab,
   or a **saved sequence file** (a plan or fill-in sequence saved for NINA; **Browse…**). Exposure,
   gain and binning are taken from the file. Choosing a source shows its frames straight away
   (prepared in the background). If frames have already been taken, it asks first, because their
   progress is replaced. A fresh **Start** always takes the frames from the source as it is at that
   moment, so later changes on the Plan tab are used.
   - The picture looks like the Plan preview: the strip to image, the current horizon, the lowest
     altitude, the frames and the shooting order (**Show shooting order**).
   - Each frame's outline shows its state: blue = to do, yellow (thick) = in progress, green = done,
     orange = skipped, red = error.
   - **When a frame is finished, its picture is drawn in its slot.** It is projected with the same
     maths as the mosaic, using the rotator offset (the Build tab's calibrated value if there is
     one). If the previews come out upside down, tick **Turn previews 180°**. **Show finished
     frames in the picture** turns them off.
   - Wheel zooms, drag moves, double-click shows the whole sky.
3. **Capture** (exposure, gain, binning) follows the Plan tab; offset stays as set in NINA or in
   the camera driver. **Settle** waits after each move. **The values can be changed while a run is
   going** (also while paused): the new exposure, gain, binning and settle time are used from the next
   frame, and the log says so. A saved sequence file only suggests its exposure, gain and binning the
   first time it is loaded; it never overwrites values you have set, and Start no longer resets them.
4. **Range:** **From frame** and **To frame** (blank = the last frame) limit the run to those frames,
   numbered in shooting order (the numbers in the picture). Frames already taken are left out
   unless **Re-shoot frames already taken** is ticked. If **From frame** (above 1) is a frame that is
   already taken – for example to redo part of a run after pausing and ending it – Start asks whether to
   shoot the range again from that frame (frames before it stay as they are) or only the frames not
   taken yet. The line below says how many will be shot.
   **Optimise route** (on by default) puts the frames still to shoot in a shorter order when you press
   Start: the planned order follows the fill-in areas, which overlap, so the mount zig-zags, crosses
   the meridian again and again (every flip costs about a minute) and the rotator turns back and forth.
   The optimiser shoots the east half (azimuth < 180°) first and then the west half – one flip – and
   shortens the path inside each half, with the rotator angles recomputed for the new order (always within
   the Plan tab's tolerance of the ideal angle). The log tells the estimated saving. Frames already taken
   stay where they are; **the numbers of the frames still to shoot change** (also in the saved capture).
   Switch it off to shoot in the planned order. In a real 553-frame fill-in run it cut the meridian
   crossings from 11 to 1 and the estimated time by about 13 %.
5. **Safety and Peripherals:** frames closer to the Sun than **Keep away from the Sun** (default 30°), below **Never
   below altitude** (default 10°) or above **Never above altitude** (default 90° = no limit; lower
   it where the camera would hit the tripod) are skipped, never slewed to. **Unpark at start** /
   **Park when finished**.
   - **Rotator limit:** **Rotator range** defaults to *Half turn (180°)*: the rotator stays between
     **Range starts at** and that angle + 180° (mechanical). A planned angle outside is shot turned
     by 180°, which covers exactly the same sky with a rectangular sensor and spares the cables.
     *Quarter turn (90°)* keeps it within 90°; an angle that a 180° turn cannot bring inside goes to
     the nearer end, so some frames are tilted by up to 45° and the coverage can have gaps.
     *Full turn* goes to the planned angles. (This replaces NINA's own rotator range setting, which
     does not apply to the mechanical moves used here.)
   - **Flat panel:** **Open flat panel** / **Close flat panel** operate the panel's cover through NINA
     (the flat device connected in NINA, Advanced API `set-cover`) or Alpaca (choose the **Flat
     panel** device, an ASCOM CoverCalibrator; it is *none* by default). Closing asks first. The check
     boxes **Open flat panel at start** and **Close flat panel when finished** do it automatically around
     a run: the run waits until the cover is open (and stops with a message if it cannot be opened), and
     the cover is closed after parking when the run has finished (not after **End run**). Without the
     automatic opening, a run warns at the start if the cover is closed, because it blocks the sky.
   - **The mount is watched while it slews.** About once a second its real position is compared
     with the Sun. If it comes much closer than planned (within 60 % of the keep-away distance) –
     because the mount took another route than estimated – the movement is stopped and the run
     pauses. If the mount cannot be reached during a slew, the run pauses as well.
   - **The mount's own limits stay in charge.** Altitude, meridian and cable limits are set in the
     mount or its driver; neither NINA nor Alpaca can read them. If a limit stops a slew, the
     arrival check below notices that the mount is not where it should be and pauses the run.
   - **Slew paths are checked too.** Before each slew the app estimates the mount's path from where
     it actually is (both axes moving together). If the path would pass closer to the Sun than the
     limit, it first slews to a safe intermediate position (one axis first, or via the meridian high
     in the sky). If there is no safe route, the frame is skipped.
   - **Sun check at Start:** the frames to shoot are timed through the run. If any would be too
     close to the Sun when its turn comes, or its slew would pass near it, Start lists them. Then
     choose: **Yes** = shoot as planned (the protections above apply), **No** = leave those frames
     out of this run (they stay "to do"), **Cancel**.
   - The Capture picture shows the Sun now, its track over the rest of the run and the keep-away
     zone, with risky frames in red.
6. **▶ Start** asks for confirmation. Tracking is left as it is: the exposures are very short and
   every slew is aimed at the current sky position, so it isn't needed. If the mount isn't tracking,
   the log says so once. Only if a slew is refused while the mount isn't tracking does the app
   switch tracking on and try that slew once more. Every slew is checked against the mount's
   RA/Dec. Then, for every frame: slew (alt/az converted to RA/Dec with the
   mount's own sidereal time: J2000 for NINA, the mount's own epoch for Alpaca) → rotator to its
   mechanical angle → settle → expose.
   - **NINA** saves the FITS with its usual headers, and Horizon Mapper finds the file in NINA's
     image folder.
   - **Alpaca:** Horizon Mapper downloads the image and saves the FITS itself (see above). If
     tracking had to be switched on for the slews, it is switched off again at the end of the run.
     Before the first frame the app checks that the mount, camera and (if used) rotator are
     connected, that the mount reports its site, and that the image folder can be written.

**Pause and errors**

- The big red **⏸ PAUSE** button, or **Esc** anywhere in the app, immediately stops the slew, the
  rotator and the exposure, and holds. **▶ Resume** repeats the interrupted frame from its start.
  While paused, **Skip frame** skips it.
- **Any error or timeout pauses the run automatically** (with a beep and the reason in the log).
  Examples: NINA or the Alpaca server not answering, a failed slew, the mount ending up somewhere
  else, a stuck rotator, a camera error, or an image that cannot be saved. Nothing is skipped silently. Fix the cause and **Resume**, or **■ End run**.
- **■ End run** ends the run and keeps the frames taken so far. **Start** again shoots only the
  frames still missing in the range. Closing the window during a run asks first and stops
  everything.

The **Run** group's progress bar shows the percentage of the current run in its centre (frames done
or skipped out of the frames this run shoots). Without a run, it shows frames taken out of all.

**Saving and resuming a capture** (group **Save / Load Capture**, right of **Run**)

- **Save current capture…** writes a `.hmcap.json` file with all the frames, each frame's state
  (taken, skipped, to do), the file path and NINA file name of every frame taken, and the capture,
  safety, range and file settings. It works at any time, also during a run. The suggested file name ends
  in `_completed` when every frame has been taken and in `_partial` otherwise (frames still to shoot,
  skipped or with an error), for example `Horizon mapping 2026-11-11 2205_partial.hmcap.json`. You can
  change the name; the automatic `last_capture.hmcap.json` keeps its fixed name.
- **Load saved capture…** restores all of that: frames, progress, settings and the current horizon
  lines.
  - Frames already taken get their pictures back from the original file paths.
  - A frame whose file no longer exists is set back to "to do", and the log lists those frames.
  - **Start** then shoots only what is still missing in the chosen range.
- The progress is also **saved automatically** after every frame to
  `%LOCALAPPDATA%\HorizonMapper\captures\last_capture.hmcap.json`. After a crash or a closed window,
  load that file to continue.

A preview is only made once its FITS file is complete: its size covers the header and all the
image data. NINA writes files in the background, so a file can be there but still empty or short.
Incomplete files, and previews that fail, are tried again every 3 seconds for up to 3 minutes.

The image folder (**NINA image folder**, or **Save images in** for Alpaca) can be changed at any
time, also during a run. File look-ups and the
folder scan always use the current value. The Capture log runs down the right side of the tab,
over the full height, with a scrollbar.

With the NINA method, NINA saves files in the background, and they can take a while to appear in
its image folder (which must be reachable from this computer, e.g. as a shared folder). The
app doesn't wait for them: it scans the folder every few seconds and matches each new FITS file to
its frame, by the name NINA reported or else by the pointing in the file header. The picture shows
these frames in light green until their file arrives. Scanning continues for 10 minutes after the
run, and late files are added to the picture and the file list.

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
- **Save Stellarium landscape…** saves the same selected line (solid or with branches, edits and side
  margin included) as a **Stellarium** polygonal landscape, so the measured horizon is drawn in
  Stellarium and blocks the sky there. See [Stellarium landscape](#stellarium-landscape) below.
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
- **Import from Capture…** (FITS folders group) brings a capture to this tab: the capture on the
  Capture tab, or a saved `.hmcap.json`. It does three things, then offers to build right away:
  - adds the list of the capture's files to the folder list;
  - takes the rotator offset from the Plan tab;
  - takes the current horizon from the capture.
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
- Use the result in Stellarium: see below.

### Stellarium landscape

**Save Stellarium landscape…** (Build tab) writes **one ZIP** holding **one folder** with `landscape.ini`,
`horizon.txt` and, once the mosaic has been built, **`panorama.png`**: a picture of your horizon made from
the mosaic. Below the horizon line it shows the mosaic's grey tones (dark obstructions, lighter gaps;
flat dark ground below the imaged altitudes); above the line it is transparent, so the sky shows through.
The picture is 3600 × 1800 pixels (0.1° per pixel), left edge = east, and is a *spherical* (photo)
landscape in Stellarium. It is grey only (the mosaic has no colour). Stellarium darkens landscapes at night, so the file sets `minimal_brightness = 0.4` to keep the picture visible; you can change it with the **Minimal brightness** option on Stellarium's Landscape tab. The horizon polygon is included
as the horizon line, so the line and the altitude list still work. Without a built mosaic (for example
if you only loaded lines) the file is the plain polygon landscape without the picture. The command
line writes `<name>_stellarium.zip` and `<name>_with_branches_stellarium.zip` next to the `.hrz` files.

- **Install:** in Stellarium press **F4 → Landscape → Add/remove landscapes… → Install a new landscape
  from a ZIP archive**, pick the ZIP and select the new landscape in the list. Or unzip the folder into
  `%APPDATA%\Stellarium\landscapes` and restart Stellarium. One ZIP holds exactly one landscape.
- **Newer version:** remove the old landscape in Stellarium first (Add/remove landscapes → Remove), then
  install the new ZIP.
- **Contents:** `horizon.txt` has one line per degree, `azimuth altitude`, azimuth 0–359 counted from
  true north towards east (the same convention as Stellarium's own polygonal landscapes), altitude in
  degrees; negative altitudes become 0. `landscape.ini` is a polygonal landscape in the layout of
  Stellarium's sample landscapes (`polygonal_angle_rotatez = 0.00001`, ground colour, horizon line
  colour, `minimal_brightness`). Everything below the line is drawn as ground.
- **Check:** in Stellarium, with the azimuth grid on, the horizon line should be at the same altitude
  as in the picture on the Build tab (for example, look east and compare).
- The **Location from landscape** option is not used: the file has no location, so Stellarium keeps your
  own location.

## Command line (optional)

```
python horizon_mapper.py plan --hrz old.hrz --floor 20 --exposure 0.008 -o "Horizon mapping.json"
python horizon_mapper.py build "Y:\Photo Import\New Horizons" --old old.hrz -o horizon --fillin "Horizon fill-in.json" --side-margin 2
```

## Requirements of the frames

- NINA FITS files with `CENTALT`, `CENTAZ`, `XPIXSZ` and `FOCALLEN` in the header (NINA writes these
  when the mount is connected).
- Frames shot by Horizon Mapper itself (NINA or Alpaca) already have these headers; Alpaca frames need
  no NINA at all.
- `ROTATOR` and `PIERSIDE` are used when present. Colour (Bayer) and mono cameras both work.

## Changes

See [CHANGELOG.md](CHANGELOG.md).
