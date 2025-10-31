//C
#include <getopt.h>
#include <stdlib.h>
//C++
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <string.h>
#include <utility>
#include <vector>
//KSW2
#include "ksw2/ksw2.h"
#include "ksw2/kseq.h"
#include "ksw2/kalloc.h"
//Single Track
#include "ksw2_singletrack.hpp"

typedef enum {
    LINEAR, 
    AFFINE,
    AFFINE2P
} aligner_t;

struct cmd_args_t {
    //input
    std::string dataset;
    //output
    std::string output_path;
    //Settings
    bool only_score   = false; 
    bool single_track = false; 
    int32_t threads   = 1;
    //penalties
    int32_t match   = 0;
    int32_t mismatch= 1;
    int32_t gapo    = 0;
    int32_t gape    = 1;
    int32_t gapo2   = 0;
    int32_t gape2   = 0;
    //detect (LINEAR, AFFINE, AFFINE2p)
    bool gapo_set  = false;
    bool gapo2_set = false;
    bool gape2_set = false;
    // constructors
    cmd_args_t() = default;
    cmd_args_t(const cmd_args_t& other) = default;
    cmd_args_t(cmd_args_t&& other) noexcept = default;
    cmd_args_t& operator=(const cmd_args_t& other) = default;
    cmd_args_t& operator=(cmd_args_t&& other) noexcept = default;
    ~cmd_args_t() = default;
};

void help() {
    std::cout <<
        "Usage: benchmark [OPTIONS]\n"
        "Options:\n"
        "  -h, --help            Show this help message\n"
        "  -d, --dataset <file>  Dataset file\n"
        "  -o, --output <file>   Output file\n"
        "  --match <int>         The match score [default=0]\n"
        "  --mismatch <int>      The mismatch score [default=1]\n"
        "  --gapo <int>          The gap open penalty (Optional)\n"
        "  --gape <int>          The gap extension penalty [default=1]\n"
        "  --gapo2 <int>         The second gap open penalty for dual gap-affine (Optional)\n"
        "  --gape2 <int>         The second gap extension penalty for dual gap-affine (Optional)\n"
        "  --single-track        Apply SingleTrack strategy\n"
        "  --only-score          Compute only alignment score (no CIGAR)\n";
}

cmd_args_t parse_arguments(int argc, char *const *argv) {
    static const struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"dataset", required_argument, 0, 'd'},
        {"only-score", no_argument, 0, 0},
        {"single-track", no_argument, 0, 0},
        {"output", required_argument, 0, 'o'},
        {"match", required_argument, 0, 0},
        {"mismatch", required_argument, 0, 0},
        {"gapo", required_argument, 0, 0},
        {"gape", required_argument, 0, 0},
        {"gapo2", required_argument, 0, 0},
        {"gape2", required_argument, 0, 0},
        {0, 0, 0, 0}
    };
    cmd_args_t args;
    int opt, option_index = 0;
    while ((opt = getopt_long(argc, argv, "hd:o:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                help();
                exit(0);
            case 'd':
                args.dataset = std::string(optarg);
                break;
            case 'o':
                args.output_path = std::string(optarg);
                break;
            case 0:
                if (strcmp(long_options[option_index].name, "match") == 0)
                    args.match = std::stoi(optarg);
                else if (strcmp(long_options[option_index].name, "mismatch") == 0)
                    args.mismatch = std::stoi(optarg);
                else if (strcmp(long_options[option_index].name, "gape") == 0)
                    args.gape = std::stoi(optarg);
                else if (strcmp(long_options[option_index].name, "only-score") == 0)
                    args.only_score = true;
                else if (strcmp(long_options[option_index].name, "single-track") == 0)
                    args.single_track = true;
                else if (strcmp(long_options[option_index].name, "gapo") == 0) {
                    args.gapo = std::stoi(optarg);
                    args.gapo_set = true;
                }
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
                std::cerr <<"Invalid option" << std::endl;
                help();
                exit(EXIT_FAILURE);
        }
    }
    if (args.dataset.empty()) {
        std::cerr <<  "Dataset file is required" << std::endl;
        help();
        exit(EXIT_FAILURE);
    }
    return args;
}

aligner_t get_aligner(const cmd_args_t& args)
{
    if (!args.gapo_set) {
        std::cout << "Running LINEAR\n";
        return LINEAR; 
    }
    if (!args.gapo2_set) {
        std::cout << "Running AFFINE\n";
        return AFFINE; 
    }
    if (!args.gape2_set && args.gapo2_set){
        std::cerr << "gape2 must be set if gapo2 is set" << std::endl;
        exit(EXIT_FAILURE);    
    }
    std::cout << "Running AFFINE2P\n";
    return AFFINE2P;
}

void write_result(std::optional<std::ofstream>& out, const cmd_args_t& args, ksw_extz_t* ez)
{
    if (!out) return;
    if (args.only_score) {
        *out << ez->score << "\n";
        return;
    }
    *out << ez->score << "\t";
    for (int i = 0; i < ez->n_cigar; ++i) {
        const uint32_t len = ez->cigar[i] >> 4;
        const uint32_t op  = ez->cigar[i] & 0xf;
        *out << len << "MID"[op];
    }
    *out << "\n";
}

void aligner_linear(const cmd_args_t& args) {
    std::ifstream dataset(args.dataset);
    if (!dataset.is_open()) {
        std::cerr << "Error: File could not be opened." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::optional<std::ofstream> out;
    if (!args.output_path.empty()) {
        out.emplace(args.output_path);
        if (!out->is_open()) {
            std::cerr << "Error: Output file could not be opened.\n";
            exit(EXIT_FAILURE);
        }
    }
    const int8_t match = static_cast<int8_t>(args.match);
    const int8_t mmatch= static_cast<int8_t>(-args.mismatch);
    const int8_t gape  = static_cast<int8_t>(args.gape); 
    // build the encoding table
    uint8_t c[256];
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; 
    void* km = km_init(); 
    std::chrono::duration<double> ksw2_timer(0);
    while (!dataset.eof()) {
        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));
        std::string target;
        std::string query;
        std::getline(dataset, target);
        std::getline(dataset, query);
        if (target.empty() || query.empty()) {
            break;
        }
        const size_t qlen = query.size()  - 1; 
        const size_t tlen = target.size() - 1;
        uint8_t *ts = new uint8_t[tlen];
        uint8_t *qs = new uint8_t[qlen];
        for (size_t i = 0; i < tlen; ++i) ts[i] = c[static_cast<uint8_t>(target[i+1])]; // encode to 0/1/2/3
        for (size_t i = 0; i < qlen; ++i) qs[i] = c[static_cast<uint8_t>(query[i+1])];
        const auto start_timer = std::chrono::high_resolution_clock::now();
        ksw_extf2_sse(km, static_cast<int32_t>(qlen), qs, static_cast<int32_t>(tlen), ts, match, mmatch, gape, -1, -1, &ez);
        const auto end_timer = std::chrono::high_resolution_clock::now();
        ksw2_timer += (end_timer - start_timer);
        write_result(out, args, &ez);
        delete[] ts; 
        delete[] qs; 
        kfree(km, ez.cigar);
    }
    km_destroy(km); 
    std::cout << "Time Alignment: " << ksw2_timer.count() << "s\n"; 
}

void aligner_affine(const cmd_args_t& args) {
    std::ifstream dataset(args.dataset);
    if (!dataset.is_open()) {
        std::cerr << "Error: File could not be opened." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::optional<std::ofstream> out;
    if (!args.output_path.empty()) {
        out.emplace(args.output_path);
        if (!out->is_open()) {
            std::cerr << "Error: Output file could not be opened.\n";
            exit(EXIT_FAILURE);
        }
    }
    const int flag = KSW_EZ_APPROX_MAX | (args.only_score ? KSW_EZ_SCORE_ONLY : 0); 
    const int8_t a = static_cast<int8_t>(args.match); 
    const int8_t b = static_cast<int8_t>(args.mismatch < 0? args.mismatch : -args.mismatch);
    const int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    const int8_t gapo    = static_cast<int8_t>(args.gapo); 
    const int8_t gape    = static_cast<int8_t>(args.gape);     
    // build the encoding table
    uint8_t c[256];
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; 
    void* km = km_init(); 
    std::chrono::duration<double> ksw2_timer(0);
    while (!dataset.eof()) {
        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));
        std::string target;
        std::string query;
        std::getline(dataset, target);
        std::getline(dataset, query);
        if (target.empty() || query.empty()) {
            break;
        }
        const size_t qlen = query.size()  - 1; 
        const size_t tlen = target.size() - 1;
        uint8_t *ts = new uint8_t[tlen];
        uint8_t *qs = new uint8_t[qlen];
        for (size_t i = 0; i < tlen; ++i) ts[i] = c[static_cast<uint8_t>(target[i+1])]; // encode to 0/1/2/3
        for (size_t i = 0; i < qlen; ++i) qs[i] = c[static_cast<uint8_t>(query[i+1])];
        const auto start_timer = std::chrono::high_resolution_clock::now();
        ksw_extz2_sse(km, static_cast<int32_t>(qlen), qs, static_cast<int32_t>(tlen), ts, 5, mat, gapo, gape, -1, -1, 0, flag, &ez);
        const auto end_timer = std::chrono::high_resolution_clock::now();
        ksw2_timer += (end_timer - start_timer);
        write_result(out, args, &ez);
        delete[] ts; 
        delete[] qs; 
        kfree(km, ez.cigar);
    }
    km_destroy(km); 
    std::cout << "Time Alignment: " << ksw2_timer.count() << "s\n"; 
}

void aligner_affine2p(const cmd_args_t& args) {
    std::ifstream dataset(args.dataset);
    if (!dataset.is_open()) {
        std::cerr << "Error: File could not be opened." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::optional<std::ofstream> out;
    if (!args.output_path.empty()) {
        out.emplace(args.output_path);
        if (!out->is_open()) {
            std::cerr << "Error: Output file could not be opened.\n";
            exit(EXIT_FAILURE);
        }
    }
    const int flag = KSW_EZ_APPROX_MAX | (args.only_score ? KSW_EZ_SCORE_ONLY : 0); 
    const int8_t a = static_cast<int8_t>(args.match);
    const int8_t b = static_cast<int8_t>(args.mismatch < 0? args.mismatch : -args.mismatch);
    const int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    const int8_t gapo    = static_cast<int8_t>(args.gapo); 
    const int8_t gape    = static_cast<int8_t>(args.gape); 
    const int8_t gapo2   = static_cast<int8_t>(args.gapo2); 
    const int8_t gape2   = static_cast<int8_t>(args.gape2);
    // build the encoding table
    uint8_t c[256];
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; 
    void* km = km_init(); 
    std::chrono::duration<double> ksw2_timer(0);
    while (!dataset.eof()) {
        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));
        std::string target;
        std::string query;
        std::getline(dataset, target);
        std::getline(dataset, query);
        if (target.empty() || query.empty()) {
            break;
        }
        const size_t qlen = query.size()  - 1; 
        const size_t tlen = target.size() - 1;
        uint8_t *ts = new uint8_t[tlen];
        uint8_t *qs = new uint8_t[qlen];
        for (size_t i = 0; i < tlen; ++i) ts[i] = c[static_cast<uint8_t>(target[i+1])]; // encode to 0/1/2/3
        for (size_t i = 0; i < qlen; ++i) qs[i] = c[static_cast<uint8_t>(query[i+1])]; 
        const auto start_timer = std::chrono::high_resolution_clock::now();
        ksw_extd2_sse(km, static_cast<int32_t>(qlen), qs, static_cast<int32_t>(tlen), ts, 5, mat, gapo, gape, gapo2, gape2, -1, -1, 0, flag, &ez); 
        const auto end_timer = std::chrono::high_resolution_clock::now();
        ksw2_timer += (end_timer - start_timer);
        write_result(out, args, &ez);
        delete[] ts; 
        delete[] qs; 
        kfree(km, ez.cigar);
    }
    km_destroy(km); 
    std::cout << "Time Alignment: " << ksw2_timer.count() << "s\n"; 
}

void aligner_affine_singletrack(const cmd_args_t& args) { 
    std::ifstream dataset(args.dataset);
    if (!dataset.is_open()) {
        std::cerr << "Error: File could not be opened." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::optional<std::ofstream> out;
    if (!args.output_path.empty()) {
        out.emplace(args.output_path);
        if (!out->is_open()) {
            std::cerr << "Error: Output file could not be opened.\n";
            exit(EXIT_FAILURE);
        }
    }
    const int flag = KSW_EZ_APPROX_MAX | (args.only_score ? KSW_EZ_SCORE_ONLY : 0); 
    const int8_t a = static_cast<int8_t>(args.match);
    const int8_t b = static_cast<int8_t>(args.mismatch < 0? args.mismatch : -args.mismatch);
    const int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    const int8_t gapo    = static_cast<int8_t>(args.gapo); 
    const int8_t gape    = static_cast<int8_t>(args.gape);   
    // build the encoding table
    uint8_t c[256];
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; 
    void* km = km_init(); 
    std::chrono::duration<double> ksw2_timer(0);
    while (!dataset.eof()) {
        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));
        std::string target;
        std::string query;
        std::getline(dataset, target);
        std::getline(dataset, query);
        if (target.empty() || query.empty()) {
            break;
        }
        const size_t qlen = query.size()  - 1; 
        const size_t tlen = target.size() - 1;
        uint8_t *ts = new uint8_t[tlen];
        uint8_t *qs = new uint8_t[qlen];
        for (size_t i = 0; i < tlen; ++i) ts[i] = c[static_cast<uint8_t>(target[i+1])]; // encode to 0/1/2/3
        for (size_t i = 0; i < qlen; ++i) qs[i] = c[static_cast<uint8_t>(query[i+1])];
        const auto start_timer = std::chrono::high_resolution_clock::now();
        ksw_extz2_singletrack_sse(km, static_cast<int32_t>(qlen), qs, static_cast<int32_t>(tlen), ts, 5, mat, gapo, gape, -1, -1, 0, flag, &ez);
        const auto end_timer = std::chrono::high_resolution_clock::now();
        ksw2_timer += (end_timer - start_timer);
        write_result(out, args, &ez);
        delete[] ts; 
        delete[] qs; 
        kfree(km, ez.cigar);
    }
    km_destroy(km); 
    std::cout << "Time Alignment: " << ksw2_timer.count() << "s\n"; 
}

void aligner_affine2p_singletrack(const cmd_args_t& args) { 
    std::ifstream dataset(args.dataset);
    if (!dataset.is_open()) {
        std::cerr << "Error: File could not be opened." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::optional<std::ofstream> out;
    if (!args.output_path.empty()) {
        out.emplace(args.output_path);
        if (!out->is_open()) {
            std::cerr << "Error: Output file could not be opened.\n";
            exit(EXIT_FAILURE);
        }
    }
    const int flag = KSW_EZ_APPROX_MAX | (args.only_score ? KSW_EZ_SCORE_ONLY : 0); 
    const int8_t a = static_cast<int8_t>(args.match);
    const int8_t b = static_cast<int8_t>(args.mismatch < 0? args.mismatch : -args.mismatch);
    const int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    const int8_t gapo    = static_cast<int8_t>(args.gapo); 
    const int8_t gape    = static_cast<int8_t>(args.gape); 
    const int8_t gapo2   = static_cast<int8_t>(args.gapo2); 
    const int8_t gape2   = static_cast<int8_t>(args.gape2);
    // build the encoding table
    uint8_t c[256];
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; 
    void* km = km_init(); 
    std::chrono::duration<double> ksw2_timer(0);
    while (!dataset.eof()) {
        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));
        std::string target;
        std::string query;
        std::getline(dataset, target);
        std::getline(dataset, query);
        if (target.empty() || query.empty()) {
            break;
        }
        const size_t qlen = query.size()  - 1; 
        const size_t tlen = target.size() - 1;
        uint8_t *ts = new uint8_t[tlen];
        uint8_t *qs = new uint8_t[qlen];
        for (size_t i = 0; i < tlen; ++i) ts[i] = c[static_cast<uint8_t>(target[i+1])]; // encode to 0/1/2/3
        for (size_t i = 0; i < qlen; ++i) qs[i] = c[static_cast<uint8_t>(query[i+1])]; 
        const auto start_timer = std::chrono::high_resolution_clock::now();
        ksw_extd2_singletrack_sse(km,static_cast<int32_t>(qlen), qs, static_cast<int32_t>(tlen), ts, 5, mat, gapo, gape, gapo2, gape2, -1, -1, 0, flag, &ez); 
        const auto end_timer = std::chrono::high_resolution_clock::now();
        ksw2_timer += (end_timer - start_timer);
        write_result(out, args, &ez);
        delete[] ts; 
        delete[] qs; 
        kfree(km, ez.cigar);
    }
    km_destroy(km); 
    std::cout << "Time Alignment: " << ksw2_timer.count() << "s\n"; 
}

int main(int argc, char *const *argv) {
    if (argc <= 1) {
        help();
        return EXIT_SUCCESS;
    }
    const cmd_args_t args = parse_arguments(argc, argv);
    aligner_t aligner_opt = get_aligner(args);
    switch(aligner_opt) {
        case LINEAR: 
            aligner_linear(args);
            break; 
        case AFFINE:
            if (args.single_track)
                aligner_affine_singletrack(args);
            else
                aligner_affine(args);
            break;
        case AFFINE2P:
            if (args.single_track)
                aligner_affine2p_singletrack(args);
            else
                aligner_affine2p(args);
            break;
        default: 
            break; 
    }
    return EXIT_SUCCESS;
}
