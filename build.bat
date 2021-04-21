@echo off

cd src
cmake -S "." -B "../bld" -GNinja || goto :error
cmake --build "../bld" || goto :error
goto EOF

:error
pause

