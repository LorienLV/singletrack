#include "dp_aligner_singletrack.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <chrono>

DPAlignerSingletrack::DPAlignerSingletrack(const Penalties& penalties,
                                           const int max_size1,
                                           const int max_size2)
    : penalties_(penalties),
      // We want the bigger sequence on the left.
      max_size_target_(std::min(max_size1, max_size2)),
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

    if (penalties_.type() == Penalties::Type::Affine ||
        penalties_.type() == Penalties::Type::DualAffine) {
        drow_.resize(max_size_target_ + 1, 0);
    }

    if (penalties_.type() == Penalties::Type::DualAffine) {
        drow2_.resize(max_size_target_ + 1, 0);
    }

    backtrace_duration_ = std::chrono::duration<double>::zero();
}

size_t DPAlignerSingletrack::memory_usage() {
    std::cerr << "TIME IN BACKTRACE: "
              << backtrace_duration_.count() << " seconds\n";

    return mmatrix_.capacity() * sizeof(mmatrix_[0]) +
           drow_.capacity() * sizeof(drow_[0]) +
           drow2_.capacity() * sizeof(drow2_[0]);
}

void DPAlignerSingletrack::align_glinear(std::string_view target,
                                         std::string_view query) {
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

void DPAlignerSingletrack::align_gaffine(std::string_view target,
                                         std::string_view query) {
    // Initialize the D row.
    for (int j = 1; j < std::ssize(target) + 1; ++j) {
        drow_[j] = mmatrix(0, j) + penalties_.gapo();
    }

    for (int i = 1; i < std::ssize(query) + 1; ++i) {

        int _icell = mmatrix(i, 0) + penalties_.gapo();

        for (int j = 1; j < std::ssize(target) + 1; ++j) {

            const int ins = std::min({
                mmatrix(i, j - 1) + penalties_.gapo() + penalties_.gape(),
                _icell + penalties_.gape()
            });

            const int del = std::min({
                mmatrix(i - 1, j) + penalties_.gapo() + penalties_.gape(),
                drow_[j] + penalties_.gape()
            });

            const int sub = mmatrix(i - 1, j - 1) +
                            penalties_.subs(target[j - 1], query[i - 1]);

            _icell = ins;
            drow_[j] = del;
            mmatrix(i, j) = std::min({sub, ins, del});
        }
    }
}

void DPAlignerSingletrack::align_dgaffine(std::string_view target,
                                          std::string_view query) {
    // Initialize the D1 and D2 rows.
    for (int j = 1; j < std::ssize(target) + 1; ++j) {
        drow_[j] = mmatrix(0, j) + penalties_.gapo();
        drow2_[j] = mmatrix(0, j) + penalties_.gapo2();
    }

    for (int i = 1; i < std::ssize(query) + 1; ++i) {

        int _icell = mmatrix(i, 0) + penalties_.gapo();
        int _icell2 = mmatrix(i, 0) + penalties_.gapo2();

        for (int j = 1; j < std::ssize(target) + 1; ++j) {

            const int ins1 = std::min({
                mmatrix(i, j - 1) + penalties_.gapo() + penalties_.gape(),
                _icell + penalties_.gape()
            });

            const int ins2 = std::min({
                mmatrix(i, j - 1) + penalties_.gapo2() + penalties_.gape2(),
                _icell2 + penalties_.gape2()
            });

            const int del1 = std::min({
                mmatrix(i - 1, j) + penalties_.gapo() + penalties_.gape(),
                drow_[j] + penalties_.gape()
            });

            const int del2 = std::min({
                mmatrix(i - 1, j) + penalties_.gapo2() + penalties_.gape2(),
                drow2_[j] + penalties_.gape2()
            });

            const int sub = mmatrix(i - 1, j - 1) +
                            penalties_.subs(target[j - 1], query[i - 1]);

            _icell = ins1;
            _icell2 = ins2;
            drow_[j] = del1;
            drow2_[j] = del2;
            mmatrix(i, j) = std::min({sub, ins1, ins2, del1, del2});
        }
    }
}

template <bool swapped>
std::string DPAlignerSingletrack::traceback_glinear(std::string_view target,
                                                    std::string_view query) {
    int i = static_cast<int>(std::ssize(query));
    int j = static_cast<int>(std::ssize(target));

    char ins_char = 'I';
    char del_char = 'D';

    if constexpr (swapped) {
        std::swap(ins_char, del_char);
    }

    int score = 0;
    std::string cigar;
    cigar.resize_and_overwrite(target.size() + query.size(), [&](char* buf, size_t size) {
        (void)size; // Unused.

        int idx = 0;
        while (i > 0 || j > 0) {
            if (j > 0 && mmatrix(i, j) == mmatrix(i, j - 1) + penalties_.gape()) {
                buf[idx] = ins_char;

                j -= 1;
            }
            else if (i > 0 && mmatrix(i, j) == mmatrix(i - 1, j) + penalties_.gape()) {
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
std::string DPAlignerSingletrack::traceback_gaffine(std::string_view target,
                                                    std::string_view query) {
    int i = static_cast<int>(std::ssize(query));
    int j = static_cast<int>(std::ssize(target));

    char ins_char = 'I';
    char del_char = 'D';

    if constexpr (swapped) {
        std::swap(ins_char, del_char);
    }

    bool in_mmatrix = true;

    int l = 0; // Length of the gap.

    int score = 0;
    std::string cigar;
    cigar.resize_and_overwrite(target.size() + query.size(), [&](char* buf, size_t size) {
        (void)size; // Unused.

        int idx = 0;

        while (i > 0 || j > 0) {
            if (in_mmatrix) {
                if (i > 0 && j > 0 &&
                    mmatrix(i, j) == mmatrix(i - 1, j - 1) +
                                     penalties_.subs(target[j - 1], query[i - 1])) {

                    buf[idx] = (target[j - 1] == query[i - 1]) ? 'M' : 'X';

                    idx += 1;

                    i -= 1;
                    j -= 1;
                }
                else {
                    in_mmatrix = false;

                    // Freeze i and j. We are in either I or D. Find a coherent path
                    // back to M.
                    l = 0;
                }
            }
            else {
                ++l;

                const int back_to_m_score = mmatrix(i, j) - penalties_.gapo() -
                                            l * penalties_.gape();

                if (j - l >= 0 && back_to_m_score == mmatrix(i, j - l)) {
                    push_to_cigar(buf + idx, ins_char, l);
                    idx += l;

                    j = j - l;

                    in_mmatrix = true;

                    continue;
                }
                else if (i - l >= 0 && back_to_m_score == mmatrix(i - l, j)) {
                    push_to_cigar(buf + idx, del_char, l);
                    idx += l;

                    i = i - l;

                    in_mmatrix = true;

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
std::string DPAlignerSingletrack::traceback_dgaffine(std::string_view target, std::string_view query) {
    auto start = std::chrono::high_resolution_clock::now();

    int i = static_cast<int>(std::ssize(query));
    int j = static_cast<int>(std::ssize(target));

    char ins_char = 'I';
    char del_char = 'D';

    if constexpr (swapped) {
        std::swap(ins_char, del_char);
    }

    bool in_mmatrix = true;

    int l = 0; // Length of the gap.

    std::string cigar;
    cigar.resize_and_overwrite(target.size() + query.size(), [&](char* buf, size_t size) {
        (void)size; // Unused.

        int idx = 0;

        while (i > 0 || j > 0) {
            if (in_mmatrix) {
                if (i > 0 && j > 0 &&
                    mmatrix(i, j) == mmatrix(i - 1, j - 1) +
                                     penalties_.subs(target[j - 1], query[i - 1])) {

                    buf[idx] = (target[j - 1] == query[i - 1]) ? 'M' : 'X';

                    idx += 1;

                    i -= 1;
                    j -= 1;
                }
                else {
                    in_mmatrix = false;

                    // Freeze i and j. We are in either I1, I2, D1 or D2.
                    // Find a coherent path back to M.
                    l = 0;
                }
            }
            else {
                ++l;

                const int back_to_m_score1 = mmatrix(i, j) - penalties_.gapo() -
                                             l * penalties_.gape();
                const int back_to_m_score2 = mmatrix(i, j) - penalties_.gapo2() -
                                             l * penalties_.gape2();

                if (j - l >= 0 && back_to_m_score1 == mmatrix(i, j - l)) {
                    push_to_cigar(buf + idx, ins_char, l);
                    idx += l;

                    j = j - l;

                    in_mmatrix = true;

                    continue;
                }
                else if (i - l >= 0 && back_to_m_score1 == mmatrix(i - l, j)) {
                    push_to_cigar(buf + idx, del_char, l);
                    idx += l;

                    i = i - l;

                    in_mmatrix = true;

                    continue;
                }
                else if (j - l >= 0 && back_to_m_score2 == mmatrix(i, j - l)) {
                    push_to_cigar(buf + idx, ins_char, l);
                    idx += l;

                    j = j - l;

                    in_mmatrix = true;

                    continue;
                }
                else if (i - l >= 0 && back_to_m_score2 == mmatrix(i - l, j)) {
                    push_to_cigar(buf + idx, del_char, l);
                    idx += l;

                    i = i - l;

                    in_mmatrix = true;

                    continue;
                }
            }
        }

        return idx;
    });

    std::reverse(cigar.begin(), cigar.end());

    auto end = std::chrono::high_resolution_clock::now();
    backtrace_duration_ += end - start;

    return cigar;
}

std::string DPAlignerSingletrack::align(std::string_view target, std::string_view query) {
    // We always want the bigger sequence on the left.
    bool swapped = false;
    if (target.size() > query.size()) {
        std::swap(target, query);
        swapped = true;
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
