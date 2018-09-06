cmake -DCMAKE_INSTALL_PREFIX=%PREFIX% %SRC_DIR% -DCMAKE_INSTALL_BINDIR=bin

"C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\MSBuild\15.0\Bin\MSBuild.exe" %SRC_DIR%\ALL_BUILD.vcxproj
