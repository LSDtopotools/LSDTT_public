#!/bin/sh
cd build/
cmake .
make
mv get_terraces_from_outlet.out ../
cd ..
