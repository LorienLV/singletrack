#include "dp_aligner_base.h"

#include <algorithm>
#include <iostream>
#include <limits>

DPAlignerBase::DPAlignerBase(const Penalties& penalties,
                             const int max_size1,
                             const int max_size2)
    : penalties_(penalties),
      // We want the bigger sequence on the left.
      max_size_target_(std::max(max_size1, max_size2)),
      max_size_query_(std::max(max_size1, max_size2)),
      mmatrix_((max_size_target_ + 1) * (max_size_query_ + 1)) {

    // Initialize M.
    mmatrix(0, 0) = 0;

    if (penalties_.type() == Penalties::Type::Linear) {
        for (int j = 1; j < max_size_target_ + 1; ++j) {
            mmatrix(0, j) = j * penalties_.gape();
        }
        for (int i = 1; i < max_size_query_ + 1; ++i) {
            mmatrix(i, 0) = i * penalties_.gape();
        }
    }
    else if (penalties_.type() == Penalties::Type::Affine) {
        for (int j = 1; j < max_size_target_ + 1; ++j) {
            mmatrix(0, j) = penalties_.gapo() + j * penalties_.gape();
        }
        for (int i = 1; i < max_size_query_ + 1; ++i) {
            mmatrix(i, 0) = penalties_.gapo() + i * penalties_.gape();
        }
    }
    else if (penalties_.type() == Penalties::Type::DualAffine) {
        for (int j = 1; j < max_size_target_ + 1; ++j) {
            mmatrix(0, j) = std::min(penalties_.gapo() + j * penalties_.gape(),
                                     penalties_.gapo2() + j * penalties_.gape2());
        }
        for (int i = 1; i < max_size_query_ + 1; ++i) {
            mmatrix(i, 0) = std::min(penalties_.gapo() + i * penalties_.gape(),
                                     penalties_.gapo2() + i * penalties_.gape2());
        }
    }

    // Initialize I and D.
    if (penalties_.type() == Penalties::Type::Affine ||
        penalties_.type() == Penalties::Type::DualAffine) {

        imatrix_.resize(mmatrix_.size());
        dmatrix_.resize(mmatrix_.size());

        imatrix(0, 0) = 0;
        imatrix(0, 0) = 0;

        for (int j = 1; j < max_size_target_ + 1; ++j) {
            imatrix(0, j) = penalties_.gapo() + j * penalties_.gape();
            dmatrix(0, j) = mmatrix(0, j) + penalties_.gapo();
        }
        for (int i = 1; i < max_size_query_ + 1; ++i) {
            dmatrix(i, 0) = penalties_.gapo() + i * penalties_.gape();
            imatrix(i, 0) = mmatrix(i, 0) + penalties_.gapo();
        }
    }

    // Initialize I2 and D2.
    if (penalties_.type() == Penalties::Type::DualAffine) {

        imatrix2_.resize(mmatrix_.size());
        dmatrix2_.resize(mmatrix_.size());

        imatrix2(0, 0) = 0;
        dmatrix2(0, 0) = 0;

        for (int j = 1; j < max_size_target_ + 1; ++j) {
            imatrix2(0, j) = penalties_.gapo2() + j * penalties_.gape2();
            dmatrix2(0, j) = mmatrix(0, j) + penalties_.gapo2();
        }

        for (int i = 1; i < max_size_query_ + 1; ++i) {
            dmatrix2(i, 0) = penalties_.gapo2() + i * penalties_.gape2();
            imatrix2(i, 0) = mmatrix(i, 0) + penalties_.gapo2();
        }
    }

    traceback_duration_ = std::chrono::duration<double>::zero();
}

size_t DPAlignerBase::memory_usage() {
    return mmatrix_.capacity() * sizeof(mmatrix_[0]) +
           imatrix_.capacity() * sizeof(imatrix_[0]) +
           dmatrix_.capacity() * sizeof(dmatrix_[0]) +
           imatrix2_.capacity() * sizeof(imatrix2_[0]) +
           dmatrix2_.capacity() * sizeof(dmatrix2_[0]);
}

void DPAlignerBase::align_glinear(std::string_view target, std::string_view query) {
    for (int i = 1; i < std::ssize(query) + 1; ++i) {
        for (int j = 1; j < std::ssize(target) + 1; ++j) {
            const int ins = mmatrix(i, j - 1) + penalties_.gape();
            const int del = mmatrix(i - 1, j) + penalties_.gape();
            const int sub = mmatrix(i - 1, j - 1) +
                            penalties_.subs(target[j - 1], query[i - 1]);

            mmatrix(i, j) = std::min({ins, del, sub});
        }
    }
}

void DPAlignerBase::align_gaffine(std::string_view target, std::string_view query) {
    for (int i = 1; i < std::ssize(query) + 1; ++i) {
        for (int j = 1; j < std::ssize(target) + 1; ++j) {

            const int ins = std::min({
                mmatrix(i, j - 1) + penalties_.gapo() + penalties_.gape(),
                imatrix(i, j - 1) + penalties_.gape(),
            });

            const int del = std::min({
                mmatrix(i - 1, j) + penalties_.gapo() + penalties_.gape(),
                dmatrix(i - 1, j) + penalties_.gape(),
            });

            const int sub = mmatrix(i - 1, j - 1) +
                            penalties_.subs(target[j - 1], query[i - 1]);

            imatrix(i, j) = ins;
            dmatrix(i, j) = del;
            mmatrix(i, j) = std::min({sub, ins, del});
        }
    }
}

void DPAlignerBase::align_dgaffine(std::string_view target, std::string_view query) {
    for (int i = 1; i < std::ssize(query) + 1; ++i) {
        for (int j = 1; j < std::ssize(target) + 1; ++j) {

            const int ins1 = std::min({
                mmatrix(i, j - 1) + penalties_.gapo() + penalties_.gape(),
                imatrix(i, j - 1) + penalties_.gape(),
            });

            const int ins2 = std::min({
                mmatrix(i, j - 1) + penalties_.gapo2() + penalties_.gape2(),
                imatrix2(i, j - 1) + penalties_.gape2(),
            });

            const int del1 = std::min({
                mmatrix(i - 1, j) + penalties_.gapo() + penalties_.gape(),
                dmatrix(i - 1, j) + penalties_.gape(),
            });

            const int del2 = std::min({
                mmatrix(i - 1, j) + penalties_.gapo2() + penalties_.gape2(),
                dmatrix2(i - 1, j) + penalties_.gape2(),
            });

            const int sub = mmatrix(i - 1, j - 1) +
                            penalties_.subs(target[j - 1], query[i - 1]);

            imatrix(i, j) = ins1;
            imatrix2(i, j) = ins2;
            dmatrix(i, j) = del1;
            dmatrix2(i, j) = del2;
            mmatrix(i, j) = std::min({sub, ins1, ins2, del1, del2});
        }
    }
}

template <bool swapped>
std::string DPAlignerBase::traceback_glinear(std::string_view target,
                                             std::string_view query) {
    int i = static_cast<int>(std::ssize(query));
    int j = static_cast<int>(std::ssize(target));

    std::string cigar;
    cigar.resize_and_overwrite(target.size() + query.size(), [&](char* buf, size_t size) {
        (void)size; // Unused.

        int idx = 0;
        while (i > 0 || j > 0) {
            if (i > 0 && mmatrix(i, j) == mmatrix(i - 1, j) + penalties_.gape()) {
                if constexpr (swapped) {
                    buf[idx] = 'I';
                }
                else {
                    buf[idx] = 'D';
                }

                i -= 1;
            }
            else if (j > 0 && mmatrix(i, j) == mmatrix(i, j - 1) + penalties_.gape()) {
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
std::string DPAlignerBase::traceback_gaffine(std::string_view target,
                                             std::string_view query) {
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
                if (imatrix(i, j) != imatrix(i, j - 1) + penalties_.gape()) {
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
                if (dmatrix(i, j) != dmatrix(i - 1, j) + penalties_.gape()) {
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
std::string DPAlignerBase::traceback_dgaffine(std::string_view target,
                                              std::string_view query) {
    auto start = std::chrono::high_resolution_clock::now();

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

        char ins_char = 'I';
        char del_char = 'D';

        if constexpr (swapped) {
            std::swap(ins_char, del_char);
        }

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
                if (imatrix(i, j) != imatrix(i, j - 1) + penalties_.gape()) {
                    curr_matrix = Matrix::M;
                }

                buf[idx] = ins_char;
                idx += 1;

                j -= 1;
            }
            else if (j > 0 && curr_matrix == Matrix::I2) {
                if (imatrix2(i, j) != imatrix2(i, j - 1) + penalties_.gape2()) {
                    curr_matrix = Matrix::M;
                }

                buf[idx] = ins_char;
                idx += 1;

                j -= 1;
            }
            else if (i > 0 && curr_matrix == Matrix::D1) {
                if (dmatrix(i, j) != dmatrix(i - 1, j) + penalties_.gape()) {
                    curr_matrix = Matrix::M;
                }

                buf[idx] = del_char;
                idx += 1;

                i -= 1;
            }
            else if (i > 0 && curr_matrix == Matrix::D2) {
                if (dmatrix2(i, j) != dmatrix2(i - 1, j) + penalties_.gape2()) {
                    curr_matrix = Matrix::M;
                }

                buf[idx] = del_char;
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

    auto end = std::chrono::high_resolution_clock::now();
    traceback_duration_ += end - start;

    return cigar;
}

std::string DPAlignerBase::align(std::string_view target, std::string_view query) {
    // We always want the bigger sequence on the left.
    bool swapped = false;
    if (target.size() > query.size()) {
        // std::swap(target, query);
        // swapped = true;
    }

    if (std::ssize(target) > max_size_target_ || std::ssize(query) > max_size_query_) {
        throw std::length_error("Target or query size exceeds the maximum size");
    }

    if (penalties_.type() == Penalties::Type::Linear) {
        align_glinear(target, query);
        if (swapped) {
            return traceback_glinear<true>(target, query);
        }
        else {
            return traceback_glinear<false>(target, query);
        }
    }
    else if (penalties_.type() == Penalties::Type::Affine) {
        align_gaffine(target, query);
        if (swapped) {
            return traceback_gaffine<true>(target, query);
        }
        else {
            return traceback_gaffine<false>(target, query);
        }
    }
    else if (penalties_.type() == Penalties::Type::DualAffine) {
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
