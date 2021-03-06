#include <findpath.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// all functions below are exported as Python functions

int init_single_findpath(std::string sequence, std::string s1, std::string s2,
                         float search_width_multiplier, bool mp = true, int en_limit = INT_MAX - 1)
{
    const char* c_sequence = sequence.c_str();
    const char* c_s1       = s1.c_str();
    const char* c_s2       = s2.c_str();

    vrna_md_t md;

    // md.temperature = 32; // args.temp FLOAT
    // md.dangles = 2; // args.dangles
    // md.special_hp = // not args.noTetra
    // md.noGU = args.noGU
    // md.noGUclosure = args.noClosingGU

//   8     vrna_md.noLP = 1
//   7     vrna_md.logML = 0
//   6     vrna_md.temperature = args.temp
//   5     vrna_md.dangles = args.dangles
//   4     vrna_md.special_hp = not args.noTetra
//   3     vrna_md.noGU = args.noGU
//   2     vrna_md.noGUclosure = args.noClosingGU


    set_model_details(&md);
    vrna_fold_compound_t* fc = vrna_fold_compound(c_sequence, &md, VRNA_OPTION_EVAL_ONLY);

    short* pt_1 = vrna_ptable(c_s1);
    short* pt_2 = vrna_ptable(c_s2);

    int bp_dist            = vrna_bp_distance(c_s1, c_s2);
    int final_search_width = bp_dist * search_width_multiplier;

    if (final_search_width < 2) { final_search_width = 2; }

    single_findpath test;
    auto            result = test.init(fc, pt_1, pt_2, final_search_width, mp, en_limit);

    // s_graph G_inner{fc, pt_1, pt_2, bp_dist, result};
    // G_inner.display_path();

    vrna_fold_compound_free(fc);
    free(pt_1);
    free(pt_2);

    if (result.size() > 0) {
        return result[0].max_en;
    } else {
        return INT_MAX - 1;
    }
}

#include <pybind11/iostream.h>

// indirect single
int init_single_findpath_i(std::string sequence, std::string s1, std::string s2,
                           float search_width_multiplier, std::vector<int> add_moves = {},
                           bool mp = true, int en_limit = INT_MAX - 1, bool Verbose = false)
{
    py::scoped_ostream_redirect stream(std::cout,                                 // std::ostream&
                                       py::module_::import("sys").attr("stdout")  // Python output
    );

    const char* c_sequence = sequence.c_str();
    const char* c_s1       = s1.c_str();
    const char* c_s2       = s2.c_str();

    vrna_md_t md;
    set_model_details(&md);
    vrna_fold_compound_t* fc = vrna_fold_compound(c_sequence, &md, VRNA_OPTION_EVAL_ONLY);

    short* pt_1 = vrna_ptable(c_s1);
    short* pt_2 = vrna_ptable(c_s2);

    int bp_dist            = vrna_bp_distance(c_s1, c_s2);
    int final_search_width = bp_dist * search_width_multiplier;

    if (final_search_width < 2) { final_search_width = 2; }

    single_findpath_i test;
    auto result = test.init(fc, pt_1, pt_2, add_moves, final_search_width, mp, en_limit, Verbose);

    // s_graph G_inner{fc, pt_1, pt_2, bp_dist, result};
    // G_inner.display_path();
    if (result.size() > 0) {
        return result[0].max_en;
    } else {
        return INT_MAX - 1;
    }


    // single_findpath_i(std::string sequence, std::string s1, std::string s2, std::vector<int> add_moves = {}, int search_width = 0, float search_width_multiplier=0,
    //                 bool mp = true, int en_limit = INT_MAX - 1, const py::dict& model_details = {})


}

int init_multi_findpath(std::string sequence, std::string s, std::vector<std::string> destinations,
                        float search_width_multiplier, int en_limit, bool mp = true)
{
    const char* c_sequence = sequence.c_str();
    vrna_md_t   md;

    set_model_details(&md);
    vrna_fold_compound_t* fc = vrna_fold_compound(c_sequence, &md, VRNA_OPTION_EVAL_ONLY);

    // single_findpath test;
    // auto            result = test.init(fc, pt_1, pt_2, final_search_width, mp);

    multi_findpath(fc, s, destinations, search_width_multiplier, en_limit);

    // free(c_sequence);
    free(fc);
    // free(md);

    return 0;
}

int init_merge_findpath(std::string sequence, std::string s1, std::string s2,
                        float search_width_multiplier, bool mp = true)
{
    auto test   = findpath(sequence, mp);
    auto result = test.init(s1, s2, search_width_multiplier);

    // if (display_path) { result.display_path(); }
    return result.max_en;
}

int init_merge_ext_findpath(std::string sequence, std::string s1, std::string s2,
                            float search_width_multiplier, bool mp = true)
{
    auto test   = findpath(sequence, mp);
    auto result = test.init_ext(s1, s2, search_width_multiplier);

    // if (display_path) { result.display_path(); }
    return result.max_en;
}

int init_mfe_findpath(std::string sequence, std::string s1, std::string s2,
                      float search_width_multiplier, bool mp = true)
{
    auto result = mfe_findpath(sequence, s1, s2, search_width_multiplier, mp);

    return result;
}

int init_vrna_findpath(std::string sequence, std::string s1, std::string s2,
                       float search_width_multiplier, bool mp = true)
{
    const char* c_sequence = sequence.c_str();
    const char* c_s1       = s1.c_str();
    const char* c_s2       = s2.c_str();

    vrna_md_t md;
    set_model_details(&md);
    vrna_fold_compound_t* fc = vrna_fold_compound(c_sequence, &md, VRNA_OPTION_EVAL_ONLY);

    int bp_dist            = vrna_bp_distance(c_s1, c_s2);
    int final_search_width = bp_dist * search_width_multiplier;

    auto result = vrna_path_findpath_saddle(fc, c_s1, c_s2, final_search_width);

    return result;
}

PYBIND11_MODULE(findpath, m)
{
    m.doc() = "pybind11 findpath";  // optional module docstring

    // findpath(std::string sequence, bool mp = true, const py::dict &model_dict = {});

    py::class_<findpath>(m, "findpath_class")
        .def(py::init<std::string, bool, py::dict>(), py::arg("sequence"), py::arg("mp"), py::arg("model_details") = py::dict())
        .def("get_en", &findpath::get_en)
        .def("get_path", &findpath::get_path)
        .def("get_sections", &findpath::get_sections)
        .def("init", &findpath::init_python);

    py::class_<single_findpath>(m, "findpath_single")
        .def(py::init<std::string, std::string, std::string, int, float, int, bool, py::dict>(), 
            py::arg("sequence"), py::arg("s1"), py::arg("s2"),             
            py::arg("search_width") = 0, 
            py::arg("search_width_multiplier") = 0, py::arg("mp") = true, py::arg("en_limit") = INT_MAX - 1,
            py::arg("model_details") = py::dict())
        .def("get_en", &single_findpath::get_en)
        .def("get_path", &single_findpath::get_path);

    py::class_<single_findpath_i>(m, "findpath_single_i")
        .def(py::init<std::string, std::string, std::string, std::vector<int>, int, float, int, bool, py::dict>(), 
            py::arg("sequence"), py::arg("s1"), py::arg("s2"),    
            py::arg("add_moves") = std::vector<int>{},         
            py::arg("search_width") = 0, 
            py::arg("search_width_multiplier") = 0, py::arg("mp") = true, py::arg("en_limit") = INT_MAX - 1,
            py::arg("model_details") = py::dict())
        .def("get_en", &single_findpath_i::get_en)
        .def("get_path", &single_findpath_i::get_path);

    m.def("init_single_findpath", &init_single_findpath, "single_findpath", py::arg("sequence"),
          py::arg("s1"), py::arg("s2"), py::arg("sw"), py::arg("mp") = true,
          py::arg("en_limit") = INT_MAX - 1);

    m.def("init_single_findpath_i", &init_single_findpath_i, "single_findpath_i",
          py::arg("sequence"), py::arg("s1"), py::arg("s2"), py::arg("sw"),
          py::arg("add_moves") = std::vector<int>{}, py::arg("mp") = true,
          py::arg("en_limit") = INT_MAX - 1, py::arg("Verbose") = false);

    m.def("init_multi_findpath", &init_multi_findpath, "multi_findpath");

    m.def("init_vrna_findpath", &init_vrna_findpath, "vrna_findpath");
    m.def("init_mfe_findpath", &init_mfe_findpath, "vrna_findpath");
    m.def("init_merge_findpath", &init_merge_findpath, "merge_findpath");
    m.def("init_merge_ext_findpath", &init_merge_ext_findpath, "merge_findpath");
}

// main func
int testfunc(std::string seq, std::string s1, std::string s2, float search_width_multiplier)
{
    // set model params
    // vrna_md_t md;
    // set_model_details(&md);
    // vrna_fold_compound_t* fc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

    // short* pt_1 = vrna_ptable(s1);
    // short* pt_2 = vrna_ptable(s2);

    // int bp_dist            = vrna_bp_distance(s1, s2);
    // int final_search_width = bp_dist * search_width_multiplier;

    // single_findpath test;
    // test.init(fc, pt_1, pt_2, final_search_width);

    return 0;
}