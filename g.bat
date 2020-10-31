@echo off

set incs=-DNDEBUG -DSTRICT -DNDEBUG -DWIN32 -D_WIN32 -DWINDOWS -Izlib ^
-D__LARGE64_FILES -DFILE_OFFSET_BITS=64 -Wno-multichar

set opts=-fstrict-aliasing -fomit-frame-pointer -ffast-math ^
-fexceptions -fno-stack-protector -fno-stack-check -fno-check-new ^
-flto -ffat-lto-objects -Wl,-flto -fuse-linker-plugin -Wl,-O -Wl,--sort-common -Wl,--as-needed -ffunction-sections

rem -fprofile-use -fprofile-correction 

set gcc=C:\MinGW\bin\g++.exe -march=k8
set path=%gcc%\..\



"%gcc%/../gcc" -c -O9 -s %incs% %opts% -static @list_c

%gcc% -O9 -s %incs% %opts% -static -fpermissive fp8sk.cpp *.o -o fp8sk.exe

del *.o
pause