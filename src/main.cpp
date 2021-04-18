
#include "main.h"

auto main(int argc, const char** argv) -> int
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




    int          buffer = 81;
    char         line[1000];
    char *       seq, *s1, *s2;
    int          E, maxkeep = 1000;
    int          verbose = 0, i;
    vrna_path_t *route, *r;

    float search_width_multiplier = 2;

    for (i = 1; i < argc; i++) {
        switch (argv[i][1]) {
            case 'm':
                if (strcmp(argv[i], "-m") == 0) sscanf(argv[++i], "%f", &search_width_multiplier);

                break;
            case 'v': verbose = !strcmp(argv[i], "-v"); break;
            case 'd':
                if (strcmp(argv[i], "-d") == 0) sscanf(argv[++i], "%d", &dangles);

                break;
                // default: usage();
        }
    }

    cut_point = -1;

    fgets(line, 1000, stdin);
    strtok(line, "\n");
    seq = vrna_cut_point_remove(line, &cut_point);

    fgets(line, 1000, stdin);
    strtok(line, "\n");
    s1 = vrna_cut_point_remove(line, &cut_point);

    fgets(line, 1000, stdin);
    strtok(line, "\n");
    s2 = vrna_cut_point_remove(line, &cut_point);


    testfunc(seq, s1, s2, search_width_multiplier);

    exit(0);
}
