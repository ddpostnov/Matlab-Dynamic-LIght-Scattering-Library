@echo off
setlocal enableextensions
REM ===========================================================================
REM  launchMyographWorkbench.bat
REM  Double-click to open the Myograph workbench GUI (guiMyograph) without
REM  starting MATLAB by hand.  MATLAB must be installed; this file finds the
REM  newest version automatically.  Keep it inside the GUIs folder - it locates
REM  the repository relative to its own location (repo root = parent of GUIs),
REM  so the whole repo can be moved freely.
REM  (For a console-free launch, use launchMyographWorkbench.vbs or the shortcut.)
REM ===========================================================================

REM --- folder that holds this script = the GUIs folder; repo root is its parent ---
set "GUIDIR=%~dp0"
if "%GUIDIR:~-1%"=="\" set "GUIDIR=%GUIDIR:~0,-1%"
for %%I in ("%GUIDIR%\..") do set "REPOROOT=%%~fI"

REM --- find the newest installed MATLAB (R20xx), highest version first ---
set "MLROOT="
for /f "delims=" %%D in ('dir /b /ad /o-n "C:\Program Files\MATLAB\R*" 2^>nul') do (
    if not defined MLROOT set "MLROOT=C:\Program Files\MATLAB\%%D"
)
if not defined MLROOT (
    echo Could not find MATLAB under "C:\Program Files\MATLAB".
    echo Edit this file and set MATLABEXE to the full path of your matlab.exe.
    pause
    exit /b 1
)
set "MATLABEXE=%MLROOT%\bin\matlab.exe"

REM --- launch: put the whole library on the path, open the GUI, wait, then quit ---
REM     guiMyograph (in GUIs\) depends on Core\Myograph\ etc., so the entire tree
REM     must be on the path - genpath(repo root) mirrors the launchers' STEP 0.
REM     MATLAB starts minimised and closes itself when the workbench window closes.
"%MATLABEXE%" -nosplash -minimize -sd "%REPOROOT%" -r "addpath(genpath('%REPOROOT%')); try, h=guiMyograph('Visible','on'); uiwait(h); catch e, disp(getReport(e)); disp('--- press any key to close ---'); pause; end; exit"

endlocal
