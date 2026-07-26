@echo off
setlocal enableextensions
REM ===========================================================================
REM  launchMyographWorkbench.bat
REM  Double-click to open the Myograph workbench GUI without starting MATLAB by
REM  hand.  MATLAB must be installed; this file finds the newest version
REM  automatically.  Keep it inside the Myograph folder - it locates the code
REM  relative to its own location, so the whole repo can be moved freely.
REM  (For a console-free launch, use launchMyographWorkbench.vbs or the shortcut.)
REM ===========================================================================

REM --- folder that holds this script = the Myograph code folder ---
set "MYODIR=%~dp0"
if "%MYODIR:~-1%"=="\" set "MYODIR=%MYODIR:~0,-1%"

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

REM --- launch: put Myograph on the path, open the GUI, wait for it, then quit ---
REM     MATLAB starts minimised and closes itself when the workbench window is closed.
"%MATLABEXE%" -nosplash -minimize -sd "%MYODIR%" -r "addpath('%MYODIR%'); try, h=myographWorkbench('Visible','on'); uiwait(h); catch e, disp(getReport(e)); disp('--- press any key to close ---'); pause; end; exit"

endlocal
