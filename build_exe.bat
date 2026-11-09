@echo off
cd /d "%~dp0"
echo Building HorizonMapper.exe ...
python -m pip install --user numpy pillow pyinstaller
python -m PyInstaller --noconfirm --onefile --windowed --name HorizonMapper horizon_mapper.py
echo.
echo Done: dist\HorizonMapper.exe
pause
