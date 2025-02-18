#include "dp_aligner.h"

#include <algorithm>
#include <iostream>
#include <limits>

DPAligner::DPAligner(const Penalties& penalties,
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

void DPAligner::align_glinear(std::string_view target, std::string_view query) {
    for (int i = 1; i < std::ssize(query) + 1; ++i) {
        for (int j = 1; j < std::ssize(target) + 1; ++j) {
            const int ins = mmatrix(i, j - 1) + _penalties.gape();
            const int del = mmatrix(i - 1, j) + _penalties.gape();
            const int sub = mmatrix(i - 1, j - 1) + subs(target[j - 1], query[i - 1]);

            mmatrix(i, j) = std::max({ins, del, sub});
        }
    }
}

void DPAligner::align_gaffine(std::string_view target, std::string_view query) {
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

void DPAligner::align_dgaffine(std::string_view target, std::string_view query) {
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
std::string DPAligner::traceback_glinear(std::string_view target, std::string_view query) {
    int i = static_cast<int>(std::ssize(query));
    int j = static_cast<int>(std::ssize(target));

    std::string cigar;
    cigar.resize_and_overwrite(target.size() + query.size(), [&](char* buf, size_t size) {
        (void)size; // Unused.

        int idx = 0;
        while (i > 0 || j > 0) {
            if (i > 0 && mmatrix(i, j) == mmatrix(i - 1, j) + _penalties.gape()) {
                if constexpr (swapped) {
                    buf[idx] = 'I';
                }
                else {
                    buf[idx] = 'D';
                }

                i -= 1;
            }
            else if (j > 0 && mmatrix(i, j) == mmatrix(i, j - 1) + _penalties.gape()) {
                if constexpr (swapped) {
                    buf[idx] = 'D';
                }
                else {
                    buf[idx] = 'I';
                }

                j -= 1;
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
std::string DPAligner::traceback_gaffine(std::string_view target, std::string_view query) {
    enum class Matrix {
        M,
        I,
        D
    };

    int i = static_cast<int>(std::ssize(query));
    int j = static_cast<int>(std::ssize(target));

    Matrix curr_matrix = Matrix::M;

    std::string cigar;
    cigar.resize_and_overwrite(target.size() + query.size(), [&](char* buf, size_t size) {
        (void)size; // Unused.

        int idx = 0;

        while (i > 0 || j > 0) {
            if (curr_matrix == Matrix::M) {
                if (j > 0 && mmatrix(i, j) == imatrix(i, j)) {
                    curr_matrix = Matrix::I;
                }
                else if (i > 0 && mmatrix(i, j) == dmatrix(i, j)) {
                    curr_matrix = Matrix::D;
                }
                else {
                    buf[idx] = (target[j - 1] == query[i - 1]) ? 'M' : 'X';

                    idx += 1;

                    i -= 1;
                    j -= 1;
                }
            }
            else if (j > 0 && curr_matrix == Matrix::I) {
                if (imatrix(i, j) != imatrix(i, j - 1) + _penalties.gape()) {
                    curr_matrix = Matrix::M;
                }

                if constexpr (swapped) {
                    buf[idx] = 'D';
                }
                else {
                    buf[idx] = 'I';
                }

                idx += 1;

                j -= 1;
            }
            else if (i > 0 && curr_matrix == Matrix::D) {
                if (dmatrix(i, j) != dmatrix(i - 1, j) + _penalties.gape()) {
                    curr_matrix = Matrix::M;
                }

                if constexpr (swapped) {
                    buf[idx] = 'I';
                }
                else {
                    buf[idx] = 'D';
                }

                idx += 1;

                i -= 1;
            }
            else {
                curr_matrix = Matrix::M;
            }
        }

        return idx;
    });

    std::reverse(cigar.begin(), cigar.end());

    return cigar;
}

template <bool swapped>
std::string DPAligner::traceback_dgaffine(std::string_view target, std::string_view query) {
    enum class Matrix {
        M,
        I1,
        D1,
        I2,
        D2
    };

    int i = static_cast<int>(std::ssize(query));
    int j = static_cast<int>(std::ssize(target));

    Matrix curr_matrix = Matrix::M;

    std::string cigar;
    cigar.resize_and_overwrite(target.size() + query.size(), [&](char* buf, size_t size) {
        (void)size; // Unused.

        int idx = 0;

        while (i > 0 || j > 0) {
            if (curr_matrix == Matrix::M) {
                if (j > 0 && mmatrix(i, j) == imatrix(i, j)) {
                    curr_matrix = Matrix::I1;
                }
                else if (j > 0 && mmatrix(i, j) == imatrix2(i, j)) {
                    curr_matrix = Matrix::I2;
                }
                else if (i > 0 && mmatrix(i, j) == dmatrix(i, j)) {
                    curr_matrix = Matrix::D1;
                }
                else if (i > 0 && mmatrix(i, j) == dmatrix2(i, j)) {
                    curr_matrix = Matrix::D2;
                }
                else {
                    buf[idx] = (target[j - 1] == query[i - 1]) ? 'M' : 'X';

                    idx += 1;

                    i -= 1;
                    j -= 1;
                }
            }
            else if (j > 0 && curr_matrix == Matrix::I1) {
                if (imatrix(i, j) != imatrix(i, j - 1) + _penalties.gape()) {
                    curr_matrix = Matrix::M;
                }

                if constexpr (swapped) {
                    buf[idx] = 'D';
                }
                else {
                    buf[idx] = 'I';
                }

                idx += 1;

                j -= 1;
            }
            else if (j > 0 && curr_matrix == Matrix::I2) {
                if (imatrix2(i, j) != imatrix2(i, j - 1) + _penalties.gape2()) {
                    curr_matrix = Matrix::M;
                }

                if constexpr (swapped) {
                    buf[idx] = 'D';
                }
                else {
                    buf[idx] = 'I';
                }

                idx += 1;

                j -= 1;
            }
            else if (i > 0 && curr_matrix == Matrix::D1) {
                if (dmatrix(i, j) != dmatrix(i - 1, j) + _penalties.gape()) {
                    curr_matrix = Matrix::M;
                }

                if constexpr (swapped) {
                    buf[idx] = 'I';
                }
                else {
                    buf[idx] = 'D';
                }

                idx += 1;

                i -= 1;
            }
            else if (i > 0 && curr_matrix == Matrix::D2) {
                if (dmatrix2(i, j) != dmatrix2(i - 1, j) + _penalties.gape2()) {
                    curr_matrix = Matrix::M;
                }

                if constexpr (swapped) {
                    buf[idx] = 'I';
                }
                else {
                    buf[idx] = 'D';
                }

                idx += 1;

                i -= 1;
            }
            else {
                curr_matrix = Matrix::M;
            }
        }

        return idx;
    });

    std::reverse(cigar.begin(), cigar.end());

    return cigar;
}

std::string DPAligner::align(std::string_view target, std::string_view query) {

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
