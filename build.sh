#! /bin/sh

# main function
#g++ ./src/main.cpp ./include/findpath.cpp -o3 -std=c++20 -o main -lRNA -lm -I./include/ -I/usr/include/python3.8 -I./usr/lib/ViennaRNA/ -DEBUG_FINDPATH -fno-lto -pthread -march=native

# pybind11 compile 
# for gcc 10.2
# g++ ./include/findpath.cpp -o3 -std=c++20 -shared -fPIC $(python3 -m pybind11 --includes) -lRNA -lm -I./include/ -I/usr/include/python3.8 -I./usr/lib/ViennaRNA/ -DEBUG_FINDPATH -fno-lto -pthread -march=native -o findpath$(python3-config --extension-suffix) 

# for gcc 9.3
# g++ ./include/findpath.cpp -o3 -std=c++2a -fconcepts -shared -fPIC $(python3 -m pybind11 --includes) -lRNA -lm -I./include/ -I/usr/include/python3.8 -I./usr/lib/ViennaRNA/ -DEBUG_FINDPATH -pthread -march=native -o findpath$(python3-config --extension-suffix) 

# for clang 11.0
clang++ ./include/findpath.cpp -o3 -std=c++20 -shared -fPIC $(python3 -m pybind11 --includes) -lRNA -I./include/ -fno-lto -pthread -march=native -o findpath$(python3-config --extension-suffix) 
