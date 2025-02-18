#include "alignment.h"

#include <iostream>
#include <cassert>

namespace alignment {

int cigar_score(const Penalties &penalties, const std::string_view cigar) {

    int score = 0;
    int score2 = 0; // Score for dual gap-affine.

    int target_idx = 0;
    int query_idx = 0;

    char prev_op = ' ';

    for (int i = 0; i < std::ssize(cigar); i++) {
        char op = cigar[i];
        if ((prev_op == 'I' || prev_op == 'D') && op != prev_op &&
            penalties.type() == Penalties::Type::DualAffine) {
            // We have finished a chain of gaps, get the best score.
            score = std::max(score, score2);
        }

        if (op == 'M' || op == '=') {
            score += penalties.match();

            target_idx += 1;
            query_idx += 1;
        }
        else if (op == 'X') {
            score += penalties.mismatch();

            target_idx += 1;
            query_idx += 1;
        }
        else if (op == 'I') {
            if (prev_op != 'I') {
                score2 = score + penalties.gapo2();
                score += penalties.gapo();
            }

            score += penalties.gape();
            score2 += penalties.gape2();

            target_idx += 1;
        }
        else if (op == 'D') {
            if (prev_op != 'D') {
                score2 = score + penalties.gapo2();
                score += penalties.gapo();
            }

            score += penalties.gape();
            score2 += penalties.gape2();

            query_idx += 1;
        }

        prev_op = op;
    }

    // In case the sequence ends with a gap.
    if ((cigar.back() == 'I' || cigar.back() == 'D') &&
        penalties.type() == Penalties::Type::DualAffine) {
        score = std::max(score, score2);
    }

    return score;
}

int alignment_score(const Penalties &penalties,
                    const std::string_view target_alignment,
                    const std::string_view query_alignment) {

    assert(target_alignment.size() == query_alignment.size());

    char prev_op = ' ';

    int score = 0;
    int score2 = 0; // Score for dual gap-affine.

    for (int i = 0; i < std::ssize(target_alignment); i++) {
        if ((prev_op == 'I' || prev_op == 'D') && target_alignment[i] != '-' &&
            query_alignment[i] != '-' &&
            penalties.type() == Penalties::Type::DualAffine) {
            // We have finished a chain of gaps, get the best score.
            score = std::max(score, score2);
        }

        if (target_alignment[i] == query_alignment[i]) {
            score += penalties.match();
            prev_op = 'M';
        }
        else if (target_alignment[i] == '-') {
            if (prev_op != 'D') {
                score2 = score + penalties.gapo2();
                score += penalties.gapo();
            }

            score += penalties.gape();
            score2 += penalties.gape2();

            prev_op = 'D';
        }
        else if (query_alignment[i] == '-') {
            if (prev_op != 'I') {
                score2 = score + penalties.gapo2();
                score += penalties.gapo();
            }

            score += penalties.gape();
            score2 += penalties.gape2();

            prev_op = 'I';
        }
        else {
            score += penalties.mismatch();
            prev_op = 'X';
        }
    }

    if ((target_alignment.back() == '-' || query_alignment.back() == '-') &&
        penalties.type() == Penalties::Type::DualAffine) {
        // We have finished a chain of gaps, get the best score.
        score = std::max(score, score2);
    }

    return score;
}

namespace {
/**
 * Transform the CIGAR to an alignment, either the target or the query,
 * controlled by the template parameter and store the result in the buffer
 * provided by the user. The buffer must have at least the same size as the
 * CIGAR string.
 *
 * @tparam target If true, the target alignment is generated, otherwise the
 * query alignment.
 * @param cigar The CIGAR string.
 * @param seq_alignment The target or query sequence.
 * @param buffer The buffer to store the alignment. It must have at least the
 * same size as the CIGAR string.
 */
template <bool target>
void cigar_to_alignment_impl(std::string_view cigar,
                             std::string_view seq,
                             char *buffer) {

    int idx = 0;
    for (int i = 0; i < std::ssize(cigar); ++i) {
        char op = cigar[i];

        if (op == 'M' || op == 'X' || op == '=') {
            buffer[i] = seq[idx];
            idx += 1;
        }
        else if (op == 'I') {
            if constexpr (target) {
                buffer[i] = seq[idx];
                idx += 1;
            }
            else {
                buffer[i] = '-';
            }
        }
        else if (op == 'D') {
            if constexpr (target) {
                buffer[i] = '-';
            }
            else {
                buffer[i] = seq[idx];
                idx += 1;
            }
        }
    }
}
} // namespace

std::pair<std::string, std::string> cigar_to_alignment(std::string_view cigar,
                                                       std::string_view target,
                                                       std::string_view query) {

    std::string target_alignment;
    target_alignment.resize_and_overwrite(cigar.size(),
                                          [&](char *buf, [[maybe_unused]] size_t size) {
        cigar_to_alignment_impl<true>(cigar, target, buf);
        return cigar.size();
    });

    std::string query_alignment;
    query_alignment.resize_and_overwrite(cigar.size(),
                                         [&](char *buf, [[maybe_unused]] size_t size) {
        cigar_to_alignment_impl<false>(cigar, query, buf);
        return cigar.size();
    });

    return {target_alignment, query_alignment};
}

void cigar_to_alignment(std::string_view cigar,
                        std::string_view target,
                        std::string_view query,
                        char *target_alignment,
                        char *query_alignment) {

    cigar_to_alignment_impl<true>(cigar, target, target_alignment);
    cigar_to_alignment_impl<false>(cigar, query, query_alignment);
}

namespace {
/**
 * Transform an alignment into a CIGAR string and store the result in the buffer
 * provided by the user. The buffer must have at least the same size as the
 * alignment strings.
 *
 * @param target_alignment The target alignment.
 * @param query_alignment The query alignment.
 * @param buffer The buffer to store the CIGAR string. It must have at least the
 * same size as the alignment strings.
 */
void alignment_to_cigar_impl(const std::string_view target_alignment,
                             const std::string_view query_alignment,
                             char *buffer) {

    if (target_alignment.size() != query_alignment.size()) {
        throw std::runtime_error(
            "The target and query alignments must have the same size.");
    }

    for (int i = 0; i < std::ssize(target_alignment); i++) {
        if (target_alignment[i] == query_alignment[i]) {
            buffer[i] = 'M';
        }
        else if (target_alignment[i] == '-') {
            buffer[i] = 'D';
        }
        else if (query_alignment[i] == '-') {
            buffer[i] = 'I';
        }
        else {
            buffer[i] = 'X';
        }
    }
}
} // namespace

std::string alignment_to_cigar(const std::string_view target_alignment,
                               const std::string_view query_alignment) {

    std::string cigar;
    cigar.resize_and_overwrite(target_alignment.size(),
                               [&](char *buf, [[maybe_unused]] size_t size) {
        alignment_to_cigar_impl(target_alignment, query_alignment, buf);
        return target_alignment.size();
    });

    return cigar;
}

void alignment_to_cigar(std::string_view target_alignment,
                        std::string_view query_alignment,
                        char *cigar) {
    alignment_to_cigar_impl(target_alignment, query_alignment, cigar);
}

}   // namespace alignment
