#include <pybind11/pybind11.h>


#include <findpath.hpp>


int add(int i, int j) {
    return i + j;


    char* seq = "AAAA";
    char* s1  = "(())";
    char* s2  = "(())";

    float search_width_multiplier = 2;

    testfunc(seq, s1, s2, search_width_multiplier);


}

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    m.def("add", &add, "A function which adds two numbers");
}

