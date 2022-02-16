

auto mfe_findpath(std::string sequence, std::string s1, std::string s2,
                      float search_width_multiplier, bool mp = true)
{
    // std::cout << "mfe calculation \n";


    const char* c_sequence = sequence.c_str();
    const char* c_s1       = s1.c_str();
    const char* c_s2       = s2.c_str();

    vrna_md_t md;
    set_model_details(&md);
    vrna_fold_compound_t* fc =
        vrna_fold_compound(c_sequence, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    short* pt1 = vrna_ptable(c_s1);
    short* pt2 = vrna_ptable(c_s2);


    auto start = std::chrono::system_clock::now();

    std::vector<short> bp_difference{};
    std::vector<short> bp_common{};

    // Ronnys code below

    for (int i = 1; i <= pt1[0]; i++) {
        // common and different base pairs
        if (pt1[i] == pt2[i] and i < pt1[i]) {
            bp_common.push_back(i);
        } else {
            bp_difference.push_back(i);
        }
        // forbid all base pairs
        for (int j = i + 4; j <= pt1[0]; j++) {
            vrna_hc_add_bp(fc, i, j,
                           VRNA_CONSTRAINT_CONTEXT_NONE | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
        }
    }

    // // (re-)allow base pairs in which s1 and s2 differ
    for (const int i : bp_difference) {
        if (i < pt1[i]) {
            vrna_hc_add_bp(fc, i, pt1[i],
                           VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
        }
        if (i < pt2[i]) {
            vrna_hc_add_bp(fc, i, pt2[i],
                           VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
        }
    }
    // force base pairs s1 and s2 have in common
    for (const int i : bp_common) {
        vrna_hc_add_bp(fc, i, pt1[i],
                       VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_ENFORCE);
    }

    char* mfe_structure = new char[pt1[0] + 1];
    vrna_mfe(fc, mfe_structure);

    short* pt_mfe = vrna_ptable(mfe_structure);
    std::string s_mfe = vrna_db_from_ptable(pt_mfe);

    int bp_dist_1 = bp_distance(pt1, pt_mfe);
    int bp_dist_2 = bp_distance(pt_mfe, pt2);



    if (std::min(bp_dist_1, bp_dist_2) == 0) {
        // launch regular merge findpath between s1 ans s2

        // std::cout << "launch regular merge fp \n";

        delete[] mfe_structure;
        vrna_fold_compound_free(fc);
        free(pt1);
        free(pt2);
        free(pt_mfe);

        auto fp_call = findpath(sequence, mp);
        auto result  = fp_call.init(s1, s2, search_width_multiplier);

        return result.max_en;
    }

    // std::cout << "launch part 1 \n";
    // std::cout << sequence << "\n";

    auto end     = std::chrono::system_clock::now();
    auto elapsed = end - start;
    start        = std::chrono::system_clock::now();

    auto part_1   = findpath(sequence, mp);
    auto result_1 = part_1.init(s1, s_mfe, search_width_multiplier);

    end           = std::chrono::system_clock::now();
    auto elapsed2 = end - start;
    start         = std::chrono::system_clock::now();

    // std::cout << "s1 -> mfe " << bp_dist_1 << "\n";
    // std::cout << vrna_db_from_ptable(pt1) << " / " << vrna_eval_structure_pt(fc, pt1) << "\n";
    // std::cout << mfe_structure << " / " << vrna_eval_structure_pt(fc, pt_mfe) << " " << result_1.max_en << "\n";

    auto part_2   = findpath(sequence, mp);
    auto result_2 = part_2.init(s_mfe, s2, search_width_multiplier);

    end           = std::chrono::system_clock::now();
    auto elapsed3 = end - start;
    start         = std::chrono::system_clock::now();

    // std::cout << "mfe -> s2 " << bp_dist_2 << "\n";
    // std::cout << mfe_structure << " / " << vrna_eval_structure_pt(fc, pt_mfe) << " " << result_2.max_en << "\n";
    // std::cout << vrna_db_from_ptable(pt2) << " / " << vrna_eval_structure_pt(fc, pt2) << "\n";


    // std::cout << "elapsed: " << elapsed.count() / 1000.0 << "/" << elapsed2.count() / 1000.0 <<
    // "/"
    //           << elapsed3.count() / 1000.0 << "\n";

    delete[] mfe_structure;
    vrna_fold_compound_free(fc);
    free(pt1);
    free(pt2);
    free(pt_mfe);

    return std::max(result_1.max_en, result_2.max_en);
}