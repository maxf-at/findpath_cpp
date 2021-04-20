

auto mfe_findpath(vrna_fold_compound_t* fc, short* pt1, short* pt2, float search_width_multiplier)
{
    // std::cout << "mfe calculation \n";

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

    int bp_dist_1 = bp_distance(pt1, pt_mfe);
    int bp_dist_2 = bp_distance(pt_mfe, pt2);

    if (std::min(bp_dist_1, bp_dist_2) == 0) {
        // launch regular merge findpath between s1 ans s2

        auto fp_call = findpath(fc, pt1, pt2, search_width_multiplier);
        auto result  = fp_call.init();

        return result.max_en;
    }

    auto end     = std::chrono::system_clock::now();
    auto elapsed = end - start;
    start        = std::chrono::system_clock::now();

    auto part_1   = findpath(fc, pt1, pt_mfe, search_width_multiplier);
    auto result_1 = part_1.init();

    end           = std::chrono::system_clock::now();
    auto elapsed2 = end - start;
    start         = std::chrono::system_clock::now();

    // std::cout << "s1 -> mfe\n";
    // std::cout << vrna_db_from_ptable(pt1) << "\n";
    // std::cout << mfe_structure << " " << result_1.max_en << "\n";
    // std::cout << vrna_db_from_ptable(pt2) << "\n";

    auto part_2   = findpath(fc, pt_mfe, pt2, search_width_multiplier);
    auto result_2 = part_2.init();

    end           = std::chrono::system_clock::now();
    auto elapsed3 = end - start;
    start         = std::chrono::system_clock::now();

    // std::cout << "mfe -> s2\n";
    // std::cout << vrna_db_from_ptable(pt1) << "\n";
    // std::cout << mfe_structure << " " << result_2.max_en << "\n";
    // std::cout << vrna_db_from_ptable(pt2) << "\n";
    // std::cout << "elapsed: " << elapsed.count() / 1000.0 << "/" << elapsed2.count() / 1000.0 <<
    // "/"
    //           << elapsed3.count() / 1000.0 << "\n";

    delete[] mfe_structure;

    return std::max(result_1.max_en, result_2.max_en);
}