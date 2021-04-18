cmake_minimum_required(VERSION 3.2)

project(findpath-project)

#Bring the headers, such as Student.h into the project
include_directories(include)


set (CMAKE_CXX_FLAGS "-fexceptions -lRNA -lm -I./usr/lib/ViennaRNA/")


#find_library(libRNA ClassLibrary /lib)
#target_link_libraries(findpath-project lRNA)

#g++ main.cpp -o3 -std=c++20 -o main findpath_B.cpp -lRNA -lm -I./usr/lib/ViennaRNA/ -DEBUG_FINDPATH -fno-lto -pthread -march=native

#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -o3 -std=c++20 -lRNA -lm -I./usr/lib/ViennaRNA/ -fno-lto -pthread -march=native")


set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

#set(DIVISIBLE_INSTALL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
#set(DIVISIBLE_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR})
#set(DIVISIBLE_INSTALL_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

#set(DIVISION_HEADERS_DIR ${PROJECT_SOURCE_DIR}/src/division)

#include_directories(${DIVISIBLE_INSTALL_INCLUDE_DIR})
#include_directories(${DIVISION_HEADERS_DIR})

add_subdirectory(src)
#add_subdirectory(test)

install(TARGETS findpath DESTINATION ./)

