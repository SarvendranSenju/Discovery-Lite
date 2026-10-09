@echo off
REM Builds Discovery Lite into dist\DiscoveryLite\ as a self-contained folder:
REM   DiscoveryLite.exe            <- double-click this one to launch the app
REM   redock_and_validate.exe      <- backend workers the GUI calls via QProcess
REM   dock_only.exe                    (see worker_command() in redock_gui.py)
REM   find_pockets.exe
REM   prepare_structures.exe
REM
REM Run this ON WINDOWS, from the project folder (with your venv activated),
REM after `pip install -r requirements.txt pyinstaller`.
REM
REM NOTE: this only bundles the Python/Qt side of the app. AutoDock Vina,
REM OpenBabel (obabel), ADFRsuite (prepare_receptor), and P2Rank (prank, needs
REM a Java runtime) are separate tools this app shells out to - they are NOT
REM bundled here and must be installed and on PATH on whatever Windows machine
REM runs the built app, same as today.

setlocal

echo === Checking prerequisites ===
python -m PyInstaller --version >nul 2>nul
if errorlevel 1 (
    echo.
    echo ERROR: could not run 'python -m PyInstaller'.
    echo Install it first with:   pip install pyinstaller
    echo ^(and make sure you're using the same Python/venv where PyQt6 is installed^)
    goto :error
)

for %%F in (redock_gui.py redock_and_validate.py dock_only.py find_pockets.py prepare_structures.py session_memory.py results_log.py viewer_tab.py) do (
    if not exist "%%F" (
        echo.
        echo ERROR: %%F not found in this folder.
        echo Run this script from inside the discovery_lite_docking_tools folder,
        echo with ALL of these files present ^(including dock_only.py and
        echo find_pockets.py, which must be copied in alongside the rest^).
        goto :error
    )
)

echo === Cleaning previous build ===
rmdir /s /q build 2>nul
rmdir /s /q dist 2>nul

echo === Freezing backend worker scripts (onefile) ===
python -m PyInstaller --onefile --console --name redock_and_validate redock_and_validate.py || goto :error
python -m PyInstaller --onefile --console --name dock_only dock_only.py || goto :error
python -m PyInstaller --onefile --console --name find_pockets find_pockets.py || goto :error
python -m PyInstaller --onefile --console --name prepare_structures prepare_structures.py || goto :error

echo === Freezing the GUI (onedir - recommended for QtWebEngine apps) ===
python -m PyInstaller --windowed --name DiscoveryLite redock_gui.py || goto :error

echo === Copying worker executables next to DiscoveryLite.exe ===
copy /y dist\redock_and_validate.exe dist\DiscoveryLite\ || goto :error
copy /y dist\dock_only.exe dist\DiscoveryLite\ || goto :error
copy /y dist\find_pockets.exe dist\DiscoveryLite\ || goto :error
copy /y dist\prepare_structures.exe dist\DiscoveryLite\ || goto :error

echo.
echo === Done. Launch dist\DiscoveryLite\DiscoveryLite.exe ===
pause
goto :eof

:error
echo.
echo === Build failed - see the error above ===
pause
exit /b 1
