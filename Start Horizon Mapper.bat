@echo off
cd /d "%~dp0"
python -c "import numpy, PIL" 2>NUL || (
  echo Installing numpy and Pillow ...
  python -m pip install --user -r requirements.txt
)
start "" pythonw horizon_mapper.py
