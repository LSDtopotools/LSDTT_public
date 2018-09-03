#!/bin/sh
if [ -d "$build" ]; then
  # Control will enter here if $DIRECTORY exists.
  rm -r build
fi

mkdir build/
cp CMakeLists.txt build/
cd build/
cmake .
make
mv lsdtt-terraces ../
cd ..
