
#include "main.h"



int testfunc(std::string sequence, std::string s1, std::string s2, float search_width_multiplier)
{

    bool mp = true;
    int en_limit = INT_MAX - 1;


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


    // free (fc);
    // free (md);

    // vrna_fold_compound_free(fc);
    // free(pt_1);
    // free(pt_2);

    // free(c_sequence);
    // free(c_s1);
    // free(c_s2);

    return 0;
}

int main(int argc, const char** argv)
{


    // cxxopts::Options options("test", "A brief description");

    // options.add_options()
    //     ("a,sequence", "Sequence", cxxopts::value<std::string>())
    //     ("b,s1", "s1", cxxopts::value<std::string>())
    //     ("c,s2", "s2", cxxopts::value<std::string>())
    //     ("d,debug", "Enable debugging", cxxopts::value<bool>()->default_value("false"))
    //     ("m,search_width_multiplier", "search width multiplier", cxxopts::value<float>()->default_value("2"))
    //     ("h,help", "Print usage");


    // // set 3 options as positional arguments
    // options.parse_positional({"sequence", "s1", "s2"});
    // auto arguments = options.parse(argc, argv);

    // if (arguments.count("help")) {
    //     std::cout << options.help() << std::endl;
    //     exit(0);
    // }

    // if (not arguments.count("sequence") or not arguments.count("s1") or not arguments.count("s2"))
    // {
    //     std::cout << options.help() << std::endl;
    //     exit(1);
    // }

    // const char* seq = arguments["sequence"].as<std::string>().c_str();
    // const char* s1 = arguments["s1"].as<std::string>().c_str();
    // const char* s2 = arguments["s2"].as<std::string>().c_str();

    // float search_width_multiplier = arguments["search_width_multiplier"].as<float>();


    // debug sample seqs

    std::string sequence = "AAAAUAAUGUACCGGACAUUCGCGCACGACCACCAUAUGGCAGAGCAUGUGUCUGUGGACCCACUAUAGCUGGGGCGCUUAACCCCAGAAAAGUAUCUUCGGUCUAUGCCUCACACGCAGCCUCCUAUUAGCAGCUCUCCUGGCCCACAAUUUUAUUAAAAGUCCAAGUUGGACUGACAAAACGCGUGCGGUGUCCUAGGGAUUGGUGGCAUAACCAGCGGUUUAAAAGCUGUGUAUAUCCGCAGCAAAUCACCGGAAAGCGGCGUUAUUAGCACCACAAAUUGAUGGUUGGUACGAGUACAAUUGCGCCGCAUAAAACCAGAGAUUCUACCCUCAAUCGGUUCUUAAGACGUACUGCGCGUUUCACCAGACCACAAUGCAGGGCGGCACCGUUAGGCAACACAACGAGACUACUCAUGCACAUAAGGAAGGUUAUCGCCAUAGACAUGGCGCGGCAGCGCAGAAUGUUUAAAUCUAAAUCUGGUAUGGGAGGCGUGCCCGUUGGUAUGAAGAAAUUUGCUGGGAGAAAAAGUCUAAGGCCUUGAAUCCGGCGGGUCUUAAUACUUACCUACAAAAUCAUCAGGCUGUACUUCCUGUAUC";
    std::string s1 =       "........(((((((((((((((((..((((.((.((((((((((...((.((((((((((....(((.((((((.......)))))).....)))....))))..(((((...(((..(((((((....(((......((((.................((((((...))))))........((((((((((((.(((((((((((...........((((......((((((......)))))).....))))....(((((((((((.(.((((((......))..)))).).)))).....))))))).........(((........))).))))))))))).)))..)))))))))....)))).......))))))).)))..))).))))).))))..)).))..))).))).)))).))..))))..((((((....)))))).....)))).)))))............))))))))((((((((.(((.(.((((.........((((..(((((.....((.((((((((((....))..)))))))))).))).))..)))))))).).))))))..))))).....";
    std::string s2 =       ".............((((((((((((.(((((.((.((((((((((...((.((((((((((...(((..((((((.......))))))....))).....))))..(((((...(((..(((((((....(((......((((.................((((((...))))))........((((((((((((.(((((((((((((.....))..((((......((((((......)))))).....))))....((((((((((((..((((........)))))))))...........))))))).........(((........))).))))))))))).)))..)))))))))....)))).......))))))).)))..))).))))).))))..)).))..))).))).)))).))..)))...((((((....))))))))...)))).))))))))............(((((((((((((.(((.((.(((((.((((..(((((((((.......(((...)))......))))))))))))).))))).))..............)))))..)))))))))))";

   
    float search_width_multiplier = 2;

    // auto en = init_single_findpath(sequence, s1, s2, search_width_multiplier);
    auto en = testfunc(sequence, s1, s2, search_width_multiplier);

    fmt::print ("result: {}\n", en);



    exit(0);
}
