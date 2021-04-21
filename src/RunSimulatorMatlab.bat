@echo off
tasklist /FI "IMAGENAME eq matlab.exe" 2>NUL | find /I /N "matlab.exe">NUL
if not "%ERRORLEVEL%"=="0" (
   matlab.exe -r "matlab.engine.shareEngine; addpath('./Matlab')"
)

CircuitSimulator.exe
   