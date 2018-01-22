#!/bin/sh
echo "I am going to build the elevation profile swath mapper for you."
if [ ! -d "build" ]; then
	echo "You don't have a build directory. I will make one."
	mkdir build/
fi
echo "Let me copy the CMakeLists into the appropriate file in the build directory."
cp CMakeLists_elevation_profiles.txt ./build/CMakeLists.txt
echo "Moving into the build directory"
cd build/
echo "Calling cmake."
cmake .
echo "Now calling make."
make
echo "Moving the program back to the driver_function directory."
mv get_elevation_profiles.out ../
echo "Going into the driver_function directory."
cd ..
echo "All done! Have a nice day."
