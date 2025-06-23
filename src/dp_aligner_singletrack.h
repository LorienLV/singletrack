#pragma once

#include <string>
#include <vector>

#include "penalties.h"

class DPAlignerSingletrack {
public:
    /**
     * Create a DPAlignerSingletrack instance with the given penalties that will
     * allocate space to align sequences one of them with a maximum size of @p
     * max_size1 and the other with a maximum size of @p max_size2. It
     * does not matter which sequence is the target and which is the query.
     *
     * @param penalties The penalties.
     * @param max_size1 The maximum size of one of the sequences.
     * @param max_size2 The maximum size of the other sequence.
     */
    DPAlignerSingletrack(const Penalties &penalties, int max_size1, int max_size2);

    /**
     * Align the target and query sequences using the Dynamic Programming.
     * If the target or query size exceeds the maximum size given at
     * construction time, an exception is thrown.
     *
     * @param target The target sequence.
     * @param query The query sequence.
     * @return The CIGAR of the alignment.
     */
    std::string align(std::string_view target, std::string_view query);

    /**
     * Get the size in bytes of all the matrices used by the DPAlignerSingletrack
     * object.
     *
     * @return The size in bytes of the DPAlignerSingletrack object.
     */
    int memory_usage();

private:
    /**
     * Alignment function for gap-linear penalties.
     *
     * @param target The target sequence.
     * @param query The query sequence.
     */
    void align_glinear(std::string_view target, std::string_view query);

    /**
     * Alignment function for gap-affine penalties.
     *
     * @param target The target sequence.
     * @param query The query sequence.
     */
    void align_gaffine(std::string_view target, std::string_view query);

    /**
     * Alignment function for dual gap-affine penalties.
     *
     * @param target The target sequence.
     * @param query The query sequence.
     */
    void align_dgaffine(std::string_view target, std::string_view query);

    /**
     * Traceback the alignment for gap-linear penalties.
     *
     * @tparam swapped If true, the target and query are swapped.
     * @param target The target sequence.
     * @param query The query sequence.
     * @return The CIGAR of the alignment.
     */
    template <bool swapped>
    std::string traceback_glinear(std::string_view target, std::string_view query);

    /**
     * Traceback the alignment for gap-affine penalties.
     *
     * @tparam swapped If true, the target and query are swapped.
     * @param target The target sequence.
     * @param query The query sequence.
     * @return The CIGAR of the alignment.
     */
    template <bool swapped>
    std::string traceback_gaffine(std::string_view target, std::string_view query);

    /**
     * Traceback the alignment for dual gap-affine penalties.
     *
     * @tparam swapped If true, the target and query are swapped.
     * @param target The target sequence.
     * @param query The query sequence.
     * @return The CIGAR of the alignment.
     */
    template <bool swapped>
    std::string traceback_dgaffine(std::string_view target, std::string_view query);

    /**
     * Add @p l characters @p c to the end of the CIGAR string @p cigar.
     *
     * @param cigar A pointer to the CIGAR string.
     * @param c The character to push.
     * @param l The number of times to push the character.
     */
    void push_to_cigar(char *cigar, char c, int l) {
        for (int i = 0; i < l; ++i) {
            cigar[i] = c;
        }
    }

    Penalties penalties_;

    int max_size_target_;
    int max_size_query_;

    std::vector<int> mmatrix_;
    int &mmatrix(int i, int j) { return mmatrix_[i * (max_size_target_ + 1) + j]; }

    std::vector<int> drow_;
    std::vector<int> drow2_;
};