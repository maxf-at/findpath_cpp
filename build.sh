#! /bin/sh

#g++ ./src/main.cpp ./include/findpath.cpp -o3 -std=c++20 -o main -lRNA -lm -I./include/ -I/usr/include/python3.8 -I./usr/lib/ViennaRNA/ -DEBUG_FINDPATH -fno-lto -pthread -march=native

# g++ ./src/example.cpp -o3 -std=c++20 -fPIC $(python3 -m pybind11 --includes) -o example$(python3-config --extension-suffix)

# g++ -O3 -Wall -shared -std=c++11 -I./include/ -fPIC $(python3 -m pybind11 --includes) ./src/example.cpp -o example$(python3-config --extension-suffix)

# g++ -O3 -Wall -std=c++20 -I./include/ -lRNA -fPIC $(python3 -m pybind11 --includes) ./src/example.cpp ./include/findpath.cpp -o example$(python3-config --extension-suffix) -fno-lto -pthread

# g++ ./include/findpath.cpp -o3 -std=c++20 -shared -fPIC -o findpath.o -lRNA -lm -I./include/ -I/usr/include/python3.8 -I./usr/lib/ViennaRNA/ -DEBUG_FINDPATH -fno-lto -pthread -march=native
# g++ -O3 -Wall -std=c++20 -I./include/ -shared -fPIC $(python3 -m pybind11 --includes) ./src/example.cpp -o part2.o -fno-lto -pthread
# g++ -O3 -Wall -shared -std=c++20 -fPIC $(python3 -m pybind11 --includes) *.o -o example$(python3-config --extension-suffix) 

# g++ ./include/findpath.cpp -o3 -std=c++20 -shared -fPIC $(python3 -m pybind11 --includes) -lRNA -lm -I./include/ -I/usr/include/python3.8 -I./usr/lib/ViennaRNA/ -DEBUG_FINDPATH -fno-lto -pthread -march=native -o findpath$(python3-config --extension-suffix) 

clang++ ./include/findpath.cpp -o3 -std=c++20 -shared -fPIC $(python3 -m pybind11 --includes) -lRNA -I./include/ -fno-lto -pthread -march=native -o findpath$(python3-config --extension-suffix) 
