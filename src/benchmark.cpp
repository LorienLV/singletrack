#include <getopt.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <memory>

#include "alignment.h"
#include "penalties.h"
#include "dp_aligner.h"

struct CMDArgs {
    std::string dataset = "";
    int match = 0;
    int mismatch = -1;
    int gapo = -1;
    bool gapo_set = false;
    int gape = -1;
    int gapo2 = -1;
    bool gapo2_set = false;
    int gape2 = -1;
    bool gape2_set = false;
};

void help() {
    std::cout << "Usage: benchmark [OPTIONS]\n"
                 "Options:\n"
                 "  -d, --dataset <file>  Dataset file\n"
                 "  --match <int>         The match score [default=0]\n"
                 "  --mismatch <int>      The mismatch score [default=-1]\n"
                 "  --gapo <int>          The gap open penalty (Optional)\n"
                 "  --gape <int>          The gap extension penalty [default=-1]\n"
                 "  --gapo2 <int>         The second gap open penalty for dual "
                                         "gap-affine (Optional)\n"
                 "  --gape2 <int>         The second gap extension penalty for "
                                         "dual gap-affine (Optional)\n";
}

CMDArgs parse_args(int argc, char *const *argv) {
    static const option long_options[] = {{"dataset", required_argument, 0, 'd'},
                                          {"match", required_argument, 0, 0},
                                          {"mismatch", required_argument, 0, 0},
                                          {"gapo", required_argument, 0, 0},
                                          {"gape", required_argument, 0, 0},
                                          {"gapo2", required_argument, 0, 0},
                                          {"gape2", required_argument, 0, 0},
                                          {0, 0, 0, 0}};

    CMDArgs args;

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "d:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'd':
                args.dataset = optarg;
                break;
            case 0:
                if (strcmp(long_options[option_index].name, "match") == 0)
                    args.match = std::stoi(optarg);
                else if (strcmp(long_options[option_index].name, "mismatch") == 0)
                    args.mismatch = std::stoi(optarg);
                else if (strcmp(long_options[option_index].name, "gapo") == 0) {
                    args.gapo = std::stoi(optarg);
                    args.gapo_set = true;
                }
                else if (strcmp(long_options[option_index].name, "gape") == 0)
                    args.gape = std::stoi(optarg);
                else if (strcmp(long_options[option_index].name, "gapo2") == 0) {
                    args.gapo2 = std::stoi(optarg);
                    args.gapo2_set = true;
                }
                else if (strcmp(long_options[option_index].name, "gape2") == 0) {
                    args.gape2 = std::stoi(optarg);
                    args.gape2_set = true;
                }
                break;
            default:
                std::cerr << "Invalid option" << std::endl;
                help();
                exit(1);
        }
    }

    return args;
}

int main(int argc, char *const *argv) {
    CMDArgs args = parse_args(argc, argv);

    if (args.dataset.empty()) {
        std::cerr << "Missing dataset file.\n";
        help();
        return 1;
    }

    // Open the file
    std::ifstream dataset(args.dataset);

    if (!dataset.is_open()) {
        std::cerr << "Could not open dataset file\n";
        help();
        return 1;
    }

    if ((args.gapo2_set || args.gape2_set) && !args.gapo_set) {
        std::cerr << "gapo must be set if gapo2 or gape2 is set\n";
        help();
        return 1;
    }

    if ((args.gapo2_set && !args.gape2_set) || (args.gape2_set && !args.gapo2_set)) {
        std::cerr << "gape2 must be set if gapo2 is set\n";
        help();
        return 1;
    }

    std::unique_ptr<Penalties> penalties;
    if (args.gapo_set) {
        if (args.gapo2_set) {
            penalties = std::make_unique<Penalties>(args.match,
                                                    args.mismatch,
                                                    args.gapo,
                                                    args.gape,
                                                    args.gapo2,
                                                    args.gape2);
        }
        else {
            penalties = std::make_unique<Penalties>(args.match,
                                                    args.mismatch,
                                                    args.gapo,
                                                    args.gape);
        }
    }
    else {
        penalties = std::make_unique<Penalties>(args.match, args.mismatch, args.gape);
    }

    if (penalties->type() == Penalties::Type::Linear) {
        std::cout << "--- Using Gap-Linear ---\n";
        std::cout << "  Match: " << penalties->match() << "\n";
        std::cout << "  Mismatch: " << penalties->mismatch() << "\n";
        std::cout << "  Gap extension: " << penalties->gape() << "\n";
    }
    else if (penalties->type() == Penalties::Type::Affine) {
        std::cout << "--- Using Gap-Affine ---\n";
        std::cout << "  Match: " << penalties->match() << "\n";
        std::cout << "  Mismatch: " << penalties->mismatch() << "\n";
        std::cout << "  Gap open: " << penalties->gapo() << "\n";
        std::cout << "  Gap extension: " << penalties->gape() << "\n";
    }
    else {
        std::cout << "--- Using Dual Gap-Affine ---\n";
        std::cout << "  Match: " << penalties->match() << "\n";
        std::cout << "  Mismatch: " << penalties->mismatch() << "\n";
        std::cout << "  Gap open 1: " << penalties->gapo() << "\n";
        std::cout << "  Gap extension 1: " << penalties->gape() << "\n";
        std::cout << "  Gap open 2: " << penalties->gapo2() << "\n";
        std::cout << "  Gap extension 2: " << penalties->gape2() << "\n";
    }
    std::cout << "\n";

    // Read the file to get the maximum sizes.
    int max_size_target = -1;
    int max_size_query = -1;
    while (!dataset.eof()) {
        std::string target;
        std::string query;

        std::getline(dataset, target);
        std::getline(dataset, query);

        if (target.empty() || query.empty()) {
            break;
        }

        max_size_target = std::max(max_size_target, static_cast<int>(target.size()));
        max_size_query = std::max(max_size_query, static_cast<int>(query.size()));
    }
    // Reset ifstream to read it again.
    dataset.clear();
    dataset.seekg(0, std::ios::beg);

    DPAligner dpa(*penalties, max_size_target, max_size_query);

    std::string target;
    target.reserve(max_size_target);

    std::string query;
    query.reserve(max_size_query);

    while (!dataset.eof()) {
        std::getline(dataset, target);
        std::getline(dataset, query);

        if (target.empty() || query.empty()) {
            break;
        }

        std::string_view target_view(target.data() + 1, target.size() - 1);
        std::string_view query_view(query.data() + 1, query.size() - 1);

        auto dpa_cigar = dpa.align(target_view, query_view);

        std::string_view dpa_cigar_view(dpa_cigar.data(), dpa_cigar.size());

        auto dpa_alignment = alignment::cigar_to_alignment(dpa_cigar_view,
                                                           target_view,
                                                           query_view);

        const int dpa_score = alignment::alignment_score(*penalties,
                                                         dpa_alignment.first,
                                                         dpa_alignment.second);

        if (dpa_score != alignment::cigar_score(*penalties, dpa_cigar_view)) {
            std::cerr << "Score mismatch\n";
            return EXIT_FAILURE;
        }

        std::cout << "Score: " << dpa_score << "\n";
        std::cout << "Alignment:\n";
        std::cout << dpa_alignment.first << "\n";
        std::cout << dpa_alignment.second << "\n";
    }

    return EXIT_SUCCESS;
}