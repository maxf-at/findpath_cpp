#! /bin/sh

# main function
#g++ ./src/main.cpp ./include/findpath.cpp -o3 -std=c++20 -o main -lRNA -lm -I./include/ -I/usr/include/python3.8 -I./usr/lib/ViennaRNA/ -DEBUG_FINDPATH -fno-lto -pthread -march=native

# pybind11 compile
# g++ ./include/findpath.cpp -o3 -std=c++20 -shared -fPIC $(python3 -m pybind11 --includes) -lRNA -lm -I./include/ -I/usr/include/python3.8 -I./usr/lib/ViennaRNA/ -DEBUG_FINDPATH -fno-lto -pthread -march=native -o findpath$(python3-config --extension-suffix) 
clang++ ./include/findpath.cpp -o3 -std=c++20 -shared -fPIC $(python3 -m pybind11 --includes) -lRNA -I./include/ -fno-lto -pthread -march=native -o findpath$(python3-config --extension-suffix) 

