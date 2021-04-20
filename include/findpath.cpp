#include <findpath.hpp>

#include <pybind11/pybind11.h>
namespace py = pybind11;


// all functions below are exported as Python functions

int init_single_findpath(std::string sequence, std::string s1, std::string s2,
                         float search_width_multiplier, bool mp = true)
{

    const char* c_sequence = sequence.c_str();
    const char* c_s1       = s1.c_str();
    const char* c_s2       = s2.c_str();

    // testfunc(c_sequence, c_s1, c_s2, search_width_multiplier);

    vrna_md_t md;
    set_model_details(&md);
    vrna_fold_compound_t* fc = vrna_fold_compound(c_sequence, &md, VRNA_OPTION_EVAL_ONLY);

    short* pt_1 = vrna_ptable(c_s1);
    short* pt_2 = vrna_ptable(c_s2);



    int bp_dist            = vrna_bp_distance(c_s1, c_s2);
    int final_search_width = bp_dist * search_width_multiplier;

    // std::cout << "start" << final_search_width << "\n";

    if (final_search_width < 2) {
        final_search_width = 2;
    }

    single_findpath test;
    auto            result = test.init(fc, pt_1, pt_2, final_search_width, mp);

    // s_graph G_inner{fc, pt_1, pt_2, bp_dist, result};
    // G_inner.display_path();

    return result[0].max_en;
}

int init_merge_findpath(std::string sequence, std::string s1, std::string s2,
                        float search_width_multiplier, bool mp = true)
{
    const char* c_sequence = sequence.c_str();
    const char* c_s1       = s1.c_str();
    const char* c_s2       = s2.c_str();

    vrna_md_t md;
    set_model_details(&md);
    // vrna_fold_compound_t* fc = vrna_fold_compound(c_sequence, &md, VRNA_OPTION_EVAL_ONLY);
    vrna_fold_compound_t* fc = vrna_fold_compound(c_sequence, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    short* pt_1 = vrna_ptable(c_s1);
    short* pt_2 = vrna_ptable(c_s2);

    int bp_dist            = vrna_bp_distance(c_s1, c_s2);
    int final_search_width = bp_dist * search_width_multiplier;




    auto test   = findpath(fc, pt_1, pt_2, search_width_multiplier);
    auto result = test.init();

    // if (display_path) { result.display_path(); }
    return result.max_en;
}

int init_merge_ext_findpath(std::string sequence, std::string s1, std::string s2,
                            float search_width_multiplier, bool mp = true)
{
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

    auto test   = findpath(fc, pt_1, pt_2, search_width_multiplier);
    // auto result = test.init();
    auto result = test.init_ext();

    // if (display_path) { result.display_path(); }
    return result.max_en;
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

int init_mfe_findpath(std::string sequence, std::string s1, std::string s2,
                       float search_width_multiplier, bool mp = true)
{
    const char* c_sequence = sequence.c_str();
    const char* c_s1       = s1.c_str();
    const char* c_s2       = s2.c_str();

    vrna_md_t md;
    set_model_details(&md);
    vrna_fold_compound_t* fc = vrna_fold_compound(c_sequence, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    short* pt_1 = vrna_ptable(c_s1);
    short* pt_2 = vrna_ptable(c_s2);


    auto result = mfe_findpath(fc, pt_1, pt_2, search_width_multiplier);

    return result;
}



PYBIND11_MODULE(findpath, m)
{
    m.doc() = "pybind11 findpath";  // optional module docstring

    m.def("init_single_findpath", &init_single_findpath, "single_findpath");
    m.def("init_vrna_findpath", &init_vrna_findpath, "vrna_findpath");
    m.def("init_mfe_findpath", &init_mfe_findpath, "vrna_findpath");
    m.def("init_merge_findpath", &init_merge_findpath, "merge_findpath");
    m.def("init_merge_ext_findpath", &init_merge_ext_findpath, "merge_findpath");
}

// main func
void testfunc(const char* seq, const char* s1, const char* s2, float search_width_multiplier)
{
    // set model params
    vrna_md_t md;
    set_model_details(&md);
    vrna_fold_compound_t* fc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

    short* pt_1 = vrna_ptable(s1);
    short* pt_2 = vrna_ptable(s2);

    int bp_dist            = vrna_bp_distance(s1, s2);
    int final_search_width = bp_dist * search_width_multiplier;

    single_findpath test;
    test.init(fc, pt_1, pt_2, final_search_width);

    return;
}