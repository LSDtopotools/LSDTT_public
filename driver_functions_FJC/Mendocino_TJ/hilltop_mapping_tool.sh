#!/bin/sh
cd build/
cmake .
make
mv hilltop_mapping_tool.out ../
cd ..
