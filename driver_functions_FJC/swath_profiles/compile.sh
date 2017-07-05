#!/bin/sh
cd build/
cmake .
make
mv get_swath_from_latlong.out ../
cd ..
