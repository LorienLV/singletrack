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

    // Initialize I and D.
    if (_penalties.type() == Penalties::Type::Affine ||
        _penalties.type() == Penalties::Type::DualAffine) {

        _imatrix.resize(_mmatrix.size());
        _dmatrix.resize(_mmatrix.size());

        imatrix(0, 0) = 0;
        imatrix(0, 0) = 0;

        for (int j = 1; j < _max_size_target + 1; ++j) {
            imatrix(0, j) = _penalties.gapo() + j * _penalties.gape();
            dmatrix(0, j) = mmatrix(0, j) + _penalties.gapo();
        }
        for (int i = 1; i < _max_size_query + 1; ++i) {
            dmatrix(i, 0) = _penalties.gapo() + i * _penalties.gape();
            imatrix(i, 0) = mmatrix(i, 0) + _penalties.gapo();
        }
    }

    // Initialize I2 and D2.
    if (_penalties.type() == Penalties::Type::DualAffine) {

        _imatrix2.resize(_mmatrix.size());
        _dmatrix2.resize(_mmatrix.size());

        imatrix2(0, 0) = 0;
        dmatrix2(0, 0) = 0;

        for (int j = 1; j < _max_size_target + 1; ++j) {
            imatrix2(0, j) = _penalties.gapo2() + j * _penalties.gape2();
            dmatrix2(0, j) = mmatrix(0, j) + _penalties.gapo2();
        }

        for (int i = 1; i < _max_size_query + 1; ++i) {
            dmatrix2(i, 0) = _penalties.gapo2() + i * _penalties.gape2();
            imatrix2(i, 0) = mmatrix(i, 0) + _penalties.gapo2();
        }
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
    for (int i = 1; i < std::ssize(query) + 1; ++i) {
        for (int j = 1; j < std::ssize(target) + 1; ++j) {

            const int ins = std::max({
                mmatrix(i, j - 1) + _penalties.gapo() + _penalties.gape(),
                imatrix(i, j - 1) + _penalties.gape(),
            });

            const int del = std::max({
                mmatrix(i - 1, j) + _penalties.gapo() + _penalties.gape(),
                dmatrix(i - 1, j) + _penalties.gape(),
            });

            const int sub = mmatrix(i - 1, j - 1) + subs(target[j - 1], query[i - 1]);

            imatrix(i, j) = ins;
            dmatrix(i, j) = del;
            mmatrix(i, j) = std::max({sub, ins, del});
        }
    }
}

void DPAlignerMod::align_dgaffine(std::string_view target, std::string_view query) {
    for (int i = 1; i < std::ssize(query) + 1; ++i) {
        for (int j = 1; j < std::ssize(target) + 1; ++j) {

            const int ins1 = std::max({
                mmatrix(i, j - 1) + _penalties.gapo() + _penalties.gape(),
                imatrix(i, j - 1) + _penalties.gape(),
            });

            const int ins2 = std::max({
                mmatrix(i, j - 1) + _penalties.gapo2() + _penalties.gape2(),
                imatrix2(i, j - 1) + _penalties.gape2(),
            });

            const int del1 = std::max({
                mmatrix(i - 1, j) + _penalties.gapo() + _penalties.gape(),
                dmatrix(i - 1, j) + _penalties.gape(),
            });

            const int del2 = std::max({
                mmatrix(i - 1, j) + _penalties.gapo2() + _penalties.gape2(),
                dmatrix2(i - 1, j) + _penalties.gape2(),
            });

            const int sub = mmatrix(i - 1, j - 1) + subs(target[j - 1], query[i - 1]);

            imatrix(i, j) = ins1;
            imatrix2(i, j) = ins2;
            dmatrix(i, j) = del1;
            dmatrix2(i, j) = del2;
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

    int j_ins = j;
    int i_del = i;

    int score_ins = 0;
    int score_del = 0;

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
                    j_ins = j;
                    i_del = i;

                    score_ins = mmatrix(i, j);
                    score_del = mmatrix(i, j);
                }
            }
            else {
                // --- Follow the ins path ---

                // We have found a path back to M.
                if (j_ins > 0 && score_ins == mmatrix(i, j_ins - 1) +
                                              _penalties.gapo() + _penalties.gape()) {
                    j_ins -= 1;

                    const int nins = j - j_ins;
                    for (int k = 0; k < nins; ++k) {
                        buf[idx] = ins_char;
                        idx += 1;
                    }

                    j = j_ins;
                    in_mmatrix = true;

                    continue;
                }
                // Keep following the ins path.
                else if (j_ins > 0) {
                    score_ins -= _penalties.gape();
                    j_ins -= 1;
                }

                // --- Follow the del path ---

                // We have found a path back to M.
                if (i_del > 0 && score_del == mmatrix(i_del - 1, j) +
                                              _penalties.gapo() + _penalties.gape()) {
                    i_del -= 1;

                    const int ndel = i - i_del;
                    for (int k = 0; k < ndel; ++k) {
                        buf[idx] = del_char;
                        idx += 1;
                    }

                    i = i_del;
                    in_mmatrix = true;

                    continue;
                }
                // Keep following the del path.
                else if (i_del > 0) {
                    score_del -= _penalties.gape();
                    i_del -= 1;
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

    int j_ins1 = j;
    int j_ins2 = j;
    int i_del1 = i;
    int i_del2 = i;

    int score_ins1 = 0;
    int score_ins2 = 0;
    int score_del1 = 0;
    int score_del2 = 0;

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
                    j_ins1 = j;
                    j_ins2 = j;
                    i_del1 = i;
                    i_del2 = i;

                    score_ins1 = mmatrix(i, j);
                    score_ins2 = mmatrix(i, j);
                    score_del1 = mmatrix(i, j);
                    score_del2 = mmatrix(i, j);
                }
            }
            else {
                // --- Follow the ins1 path ---

                // We have found a path back to M.
                if (j_ins1 > 0 && score_ins1 == mmatrix(i, j_ins1 - 1) +
                                                _penalties.gapo() + _penalties.gape()) {
                    j_ins1 -= 1;

                    const int nins = j - j_ins1;
                    for (int k = 0; k < nins; ++k) {
                        buf[idx] = ins_char;
                        idx += 1;
                    }

                    j = j_ins1;
                    in_mmatrix = true;

                    continue;
                }
                // Keep following the ins1 path.
                else if (j_ins1 > 0) {
                    score_ins1 -= _penalties.gape();
                    j_ins1 -= 1;
                }

                // --- Follow the ins2 path ---

                // We have found a path back to M.
                if (j_ins2 > 0 && score_ins2 == mmatrix(i, j_ins2 - 1) +
                                                _penalties.gapo2() + _penalties.gape2()) {
                    j_ins2 -= 1;

                    const int nins = j - j_ins2;
                    for (int k = 0; k < nins; ++k) {
                        buf[idx] = ins_char;
                        idx += 1;
                    }

                    j = j_ins2;
                    in_mmatrix = true;

                    continue;
                }
                // Keep following the ins2 path.
                else if (j_ins2 > 0) {
                    score_ins2 -= _penalties.gape2();
                    j_ins2 -= 1;
                }

                // --- Follow the del1 path ---

                // We have found a path back to M.
                if (i_del1 > 0 && score_del1 == mmatrix(i_del1 - 1, j) +
                                                _penalties.gapo() + _penalties.gape()) {
                    i_del1 -= 1;

                    const int ndel = i - i_del1;
                    for (int k = 0; k < ndel; ++k) {
                        buf[idx] = del_char;
                        idx += 1;
                    }

                    i = i_del1;
                    in_mmatrix = true;

                    continue;
                }
                // Keep following the del1 path.
                else if (i_del1 > 0) {
                    score_del1 -= _penalties.gape();
                    i_del1 -= 1;
                }

                // --- Follow the del2 path ---

                // We have found a path back to M.
                if (i_del2 > 0 && score_del2 == mmatrix(i_del2 - 1, j) +
                                                _penalties.gapo2() + _penalties.gape2()) {
                    i_del2 -= 1;

                    const int ndel = i - i_del2;
                    for (int k = 0; k < ndel; ++k) {
                        buf[idx] = del_char;
                        idx += 1;
                    }

                    i = i_del2;
                    in_mmatrix = true;

                    continue;
                }
                // Keep following the del2 path.
                else if (i_del2 > 0) {
                    score_del2 -= _penalties.gape2();
                    i_del2 -= 1;
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
