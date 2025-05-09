#pragma once

#include <string>
#include <vector>

#include "alignment.h"
#include "penalties.h"

class DPAligner {
public:
    /**
     * Create a DPAligner instance with the given penalties that will allocate
     * space to align sequences one of them with a maximum size of @p
     * max_size1 and the other with a maximum size of @p max_size2. It
     * does not matter which sequence is the target and which is the query.
     *
     * @param penalties The penalties.
     * @param max_size1 The maximum size of one of the sequences.
     * @param max_size2 The maximum size of the other sequence.
     */
    DPAligner(const Penalties &penalties,
              int max_size1,
              int max_size2);

    /**
     * Align the target and query sequences using the Dynamic Programming.
     * If the target or query size exceeds the maximum size given at construction
     * time, an exception is thrown.
     *
     * @param target The target sequence.
     * @param query The query sequence.
     * @return The CIGAR of the alignment.
     */
    std::string align(std::string_view target, std::string_view query);

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

    Penalties _penalties;

    int _max_size_target;
    int _max_size_query;

    std::vector<int> _mmatrix;
    int &mmatrix(int i, int j) { return _mmatrix[i * (_max_size_target + 1) + j]; }

    std::vector<int> _imatrix;
    int &imatrix(int i, int j) { return _imatrix[i * (_max_size_target + 1) + j]; }

    std::vector<int> _dmatrix;
    int &dmatrix(int i, int j) { return _dmatrix[i * (_max_size_target + 1) + j]; }

    std::vector<int> _imatrix2;
    int &imatrix2(int i, int j) { return _imatrix2[i * (_max_size_target + 1) + j]; }

    std::vector<int> _dmatrix2;
    int &dmatrix2(int i, int j) { return _dmatrix2[i * (_max_size_target + 1) + j]; }
};