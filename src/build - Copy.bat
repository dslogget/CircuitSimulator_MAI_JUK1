@echo off

cmake -DCMAKE_BUILD_TYPE=Release -S "." -B "../releasebld" -GNinja || goto :error
cmake --build "../releasebld" || goto :error
goto EOF

:error
pause
:EOF
