#include "dp_aligner_mod.h"

#include <algorithm>
#include <iostream>
#include <limits>

DPAlignerMod::DPAlignerMod(const Penalties& penalties,
                           const int max_size1,
                           const int max_size2)
    : _penalties(penalties),
      // We want the bigger sequence on the left.
      _max_size_target(std::min(max_size1, max_size2)),
      _max_size_query(std::max(max_size1, max_size2)),
      _mmatrix((_max_size_target + 1) * (_max_size_query + 1)) {

    // Initialize M.
    mmatrix(0, 0) = 0;

    if (_penalties.type() == Penalties::Type::Linear) {
        for (int j = 1; j < _max_size_target + 1; ++j) {
            mmatrix(0, j) = j * _penalties.gape();
        }
        for (int i = 1; i < _max_size_query + 1; ++i) {
            mmatrix(i, 0) = i * _penalties.gape();
        }
    }
    else if (_penalties.type() == Penalties::Type::Affine) {
        for (int j = 1; j < _max_size_target + 1; ++j) {
            mmatrix(0, j) = _penalties.gapo() + j * _penalties.gape();
        }
        for (int i = 1; i < _max_size_query + 1; ++i) {
            mmatrix(i, 0) = _penalties.gapo() + i * _penalties.gape();
        }
    }
    else if (_penalties.type() == Penalties::Type::DualAffine) {
        for (int j = 1; j < _max_size_target + 1; ++j) {
            mmatrix(0, j) = std::max(_penalties.gapo() + j * _penalties.gape(),
                                     _penalties.gapo2() + j * _penalties.gape2());
        }
        for (int i = 1; i < _max_size_query + 1; ++i) {
            mmatrix(i, 0) = std::max(_penalties.gapo() + i * _penalties.gape(),
                                     _penalties.gapo2() + i * _penalties.gape2());
        }
    }

    if (_penalties.type() == Penalties::Type::Affine ||
        _penalties.type() == Penalties::Type::DualAffine) {
        _drow.resize(_max_size_target + 1, 0);
    }

    if (_penalties.type() == Penalties::Type::DualAffine) {
        _drow2.resize(_max_size_target + 1, 0);
    }
}

void DPAlignerMod::align_glinear(std::string_view target, std::string_view query) {
    for (int i = 1; i < std::ssize(query) + 1; ++i) {
        for (int j = 1; j < std::ssize(target) + 1; ++j) {
            const int ins = mmatrix(i, j - 1) + _penalties.gape();
            const int del = mmatrix(i - 1, j) + _penalties.gape();
            const int sub = mmatrix(i - 1, j - 1) + subs(target[j - 1], query[i - 1]);

            mmatrix(i, j) = std::max({ins, del, sub});
        }
    }
}

void DPAlignerMod::align_gaffine(std::string_view target, std::string_view query) {
    // Initialize the D row.
    for (int j = 1; j < std::ssize(target) + 1; ++j) {
        _drow[j] = mmatrix(0, j) + _penalties.gapo();
    }

    for (int i = 1; i < std::ssize(query) + 1; ++i) {

        int _icell = mmatrix(i, 0) + _penalties.gapo();

        for (int j = 1; j < std::ssize(target) + 1; ++j) {

            const int ins = std::max({
                mmatrix(i, j - 1) + _penalties.gapo() + _penalties.gape(),
                _icell + _penalties.gape()
            });

            const int del = std::max({
                mmatrix(i - 1, j) + _penalties.gapo() + _penalties.gape(),
                _drow[j] + _penalties.gape()
            });

            const int sub = mmatrix(i - 1, j - 1) + subs(target[j - 1], query[i - 1]);

            _icell = ins;
            _drow[j] = del;
            mmatrix(i, j) = std::max({sub, ins, del});
        }
    }
}

void DPAlignerMod::align_dgaffine(std::string_view target, std::string_view query) {
    // Initialize the D1 and D2 rows.
    for (int j = 1; j < std::ssize(target) + 1; ++j) {
        _drow[j] = mmatrix(0, j) + _penalties.gapo();
        _drow2[j] = mmatrix(0, j) + _penalties.gapo2();
    }

    for (int i = 1; i < std::ssize(query) + 1; ++i) {

        int _icell = mmatrix(i, 0) + _penalties.gapo();
        int _icell2 = mmatrix(i, 0) + _penalties.gapo2();

        for (int j = 1; j < std::ssize(target) + 1; ++j) {

            const int ins1 = std::max({
                mmatrix(i, j - 1) + _penalties.gapo() + _penalties.gape(),
                _icell + _penalties.gape()
            });

            const int ins2 = std::max({
                mmatrix(i, j - 1) + _penalties.gapo2() + _penalties.gape2(),
                _icell2 + _penalties.gape2()
            });

            const int del1 = std::max({
                mmatrix(i - 1, j) + _penalties.gapo() + _penalties.gape(),
                _drow[j] + _penalties.gape()
            });

            const int del2 = std::max({
                mmatrix(i - 1, j) + _penalties.gapo2() + _penalties.gape2(),
                _drow2[j] + _penalties.gape2()
            });

            const int sub = mmatrix(i - 1, j - 1) + subs(target[j - 1], query[i - 1]);

            _icell = ins1;
            _icell2 = ins2;
            _drow[j] = del1;
            _drow2[j] = del2;
            mmatrix(i, j) = std::max({sub, ins1, ins2, del1, del2});
        }
    }
}

template <bool swapped>
std::string DPAlignerMod::traceback_glinear(std::string_view target, std::string_view query) {
    int i = static_cast<int>(std::ssize(query));
    int j = static_cast<int>(std::ssize(target));

    char ins_char = 'I';
    char del_char = 'D';

    if constexpr (swapped) {
        std::swap(ins_char, del_char);
    }

    std::string cigar;
    cigar.resize_and_overwrite(target.size() + query.size(), [&](char* buf, size_t size) {
        (void)size; // Unused.

        int idx = 0;
        while (i > 0 || j > 0) {
            if (j > 0 && mmatrix(i, j) == mmatrix(i, j - 1) + _penalties.gape()) {
                buf[idx] = ins_char;

                j -= 1;
            }
            else if (i > 0 && mmatrix(i, j) == mmatrix(i - 1, j) + _penalties.gape()) {
                buf[idx] = del_char;

                i -= 1;
            }
            else {
                buf[idx] = (target[j - 1] == query[i - 1]) ? 'M' : 'X';

                i -= 1;
                j -= 1;
            }


            idx += 1;
        }

        return idx;
    });

    std::reverse(cigar.begin(), cigar.end());

    return cigar;
}

template <bool swapped>
std::string DPAlignerMod::traceback_gaffine(std::string_view target, std::string_view query) {
    int i = static_cast<int>(std::ssize(query));
    int j = static_cast<int>(std::ssize(target));

    int nindels = 0;

    char ins_char = 'I';
    char del_char = 'D';

    if constexpr (swapped) {
        std::swap(ins_char, del_char);
    }

    bool in_mmatrix = true;

    std::string cigar;
    cigar.resize_and_overwrite(target.size() + query.size(), [&](char* buf, size_t size) {
        (void)size; // Unused.

        int idx = 0;

        while (i > 0 || j > 0) {
            if (in_mmatrix) {
                if (i > 0 && j > 0 &&
                    mmatrix(i, j) == mmatrix(i - 1, j - 1) +
                                     subs(target[j - 1], query[i - 1])) {

                    buf[idx] = (target[j - 1] == query[i - 1]) ? 'M' : 'X';

                    idx += 1;

                    i -= 1;
                    j -= 1;
                }
                else {
                    in_mmatrix = false;

                    // Freeze i and j. We are in either I or D. Find a coherent path
                    // back to M.
                    nindels = 0;
                }
            }
            else {
                auto find_path_back_to_m = [&](const int new_i,
                                               const int new_j,
                                               const int expected_score,
                                               const char indel_char) {

                    if (new_i < 0 || new_j < 0) {
                        return false;
                    }

                    if (expected_score != mmatrix(new_i, new_j)) {
                        return false;
                    }

                    i = new_i;
                    j = new_j;

                    in_mmatrix = true;

                    for (int k = 0; k < nindels; ++k) {
                        buf[idx] = indel_char;
                        idx += 1;
                    }

                    return true;
                };

                ++nindels;

                const int indel_i = i - nindels;
                const int indel_j = j - nindels;
                const int back_to_m_score = mmatrix(i, j) -
                                            nindels * _penalties.gape() -
                                            _penalties.gapo();

                if (find_path_back_to_m(i, indel_j, back_to_m_score, ins_char)) {
                    // Ins path.
                    continue;
                }
                else if (find_path_back_to_m(indel_i, j, back_to_m_score, del_char)) {
                    // Del path.
                    continue;
                }
            }
        }

        return idx;
    });

    std::reverse(cigar.begin(), cigar.end());

    return cigar;
}

template <bool swapped>
std::string DPAlignerMod::traceback_dgaffine(std::string_view target, std::string_view query) {
    int i = static_cast<int>(std::ssize(query));
    int j = static_cast<int>(std::ssize(target));

    int nindels = 0;

    char ins_char = 'I';
    char del_char = 'D';

    if constexpr (swapped) {
        std::swap(ins_char, del_char);
    }

    bool in_mmatrix = true;

    std::string cigar;
    cigar.resize_and_overwrite(target.size() + query.size(), [&](char* buf, size_t size) {
        (void)size; // Unused.

        int idx = 0;

        while (i > 0 || j > 0) {
            if (in_mmatrix) {
                if (i > 0 && j > 0 &&
                    mmatrix(i, j) == mmatrix(i - 1, j - 1) +
                                     subs(target[j - 1], query[i - 1])) {

                    buf[idx] = (target[j - 1] == query[i - 1]) ? 'M' : 'X';

                    idx += 1;

                    i -= 1;
                    j -= 1;
                }
                else {
                    in_mmatrix = false;

                    // Freeze i and j. We are in either I1, I2, D1 or D2.
                    // Find a coherent path back to M.
                    nindels = 0;
                }
            }
            else {
                auto find_path_back_to_m = [&](const int new_i,
                                               const int new_j,
                                               const int expected_score,
                                               const char indel_char) {
                    if (new_i < 0 || new_j < 0) {
                        return false;
                    }

                    if (expected_score != mmatrix(new_i, new_j)) {
                        return false;
                    }

                    i = new_i;
                    j = new_j;

                    in_mmatrix = true;

                    for (int k = 0; k < nindels; ++k) {
                        buf[idx] = indel_char;
                        idx += 1;
                    }

                    return true;
                };

                ++nindels;

                const int indel_i = i - nindels;
                const int indel_j = j - nindels;
                const int back_to_m_score1 = mmatrix(i, j) -
                                            nindels * _penalties.gape() -
                                            _penalties.gapo();
                const int back_to_m_score2 = mmatrix(i, j) -
                                            nindels * _penalties.gape2() -
                                            _penalties.gapo2();

                if (find_path_back_to_m(i, indel_j, back_to_m_score1, ins_char)) {
                    // Ins1 path.
                    continue;
                }
                else if (find_path_back_to_m(indel_i, j, back_to_m_score1, del_char)) {
                    // Del1 path.
                    continue;
                }
                else if (find_path_back_to_m(i, indel_j, back_to_m_score2, ins_char)) {
                    // Ins2 path.
                    continue;
                }
                else if (find_path_back_to_m(indel_i, j, back_to_m_score2, del_char)) {
                    // Del2 path.
                    continue;
                }
            }
        }

        return idx;
    });

    std::reverse(cigar.begin(), cigar.end());

    return cigar;
}

std::string DPAlignerMod::align(std::string_view target, std::string_view query) {

    // We always want the bigger sequence on the left.
    bool swapped = false;
    if (target.size() > query.size()) {
        std::swap(target, query);
        swapped = true;
    }

    if (std::ssize(target) > _max_size_target || std::ssize(query) > _max_size_query) {
        throw std::length_error("Target or query size exceeds the maximum size");
    }

    if (_penalties.type() == Penalties::Type::Linear) {
        align_glinear(target, query);
        if (swapped) {
            return traceback_glinear<true>(target, query);
        }
        else {
            return traceback_glinear<false>(target, query);
        }
    }
    else if (_penalties.type() == Penalties::Type::Affine) {
        align_gaffine(target, query);
        if (swapped) {
            return traceback_gaffine<true>(target, query);
        }
        else {
            return traceback_gaffine<false>(target, query);
        }
    }
    else if (_penalties.type() == Penalties::Type::DualAffine) {
        align_dgaffine(target, query);
        if (swapped) {
            return traceback_dgaffine<true>(target, query);
        }
        else {
            return traceback_dgaffine<false>(target, query);
        }
    }
    else {
        throw std::runtime_error("Unknown penalties type");
    }
}
