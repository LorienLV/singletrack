#pragma once

#include <string>

#include "penalties.h"

/**
 * Helper functions to handle CIGARS and alignments.
 *
 * A CIGAR is a sequence of characters with the following possible values:
 * - 'M' or '=': Match or mismatch.
 * - 'I': Insertion.
 * - 'D': Deletion.
 * - 'X': Substitution.
 *
 * An alignment is a pair of strings with the same length where each entry is
 * either a letter of the alphabet, or a gap represented by the '-' character.
 *
 */

namespace alignment {

/**
 * Check if the CIGAR is coherent with the target and query, i.e., all the
 * characters in the sequences are present in the CIGAR.
 *
 * @param cigar The CIGAR string.
 * @param target The target sequence.
 * @param query The query sequence.
 * @return true if the CIGAR is coherent with the target and query, false
 * otherwise.
 */
bool cigar_coherent(std::string_view cigar,
                    std::string_view target,
                    std::string_view query);

/**
 * Compute the score of a CIGAR given the penalties and the target and query
 * that originated the CIGAR.
 *
 * @param penalties The penalties.
 * @param cigar The CIGAR.
 * @param target The target sequence.
 * @param query The query sequence.
 * @return The score of the CIGAR.
 */
int cigar_score(const Penalties &penalties, std::string_view cigar);

/**
 * Check if the alignment is coherent with the target and query, i.e., all the
 * characters in the sequences are present in the alignment.
 *
 * @param target_alignment The target alignment.
 * @param query_alignment The query alignment.
 * @param target The target sequence.
 * @param query The query sequence.
 * @return true if the alignment is coherent with the target and query, false
 * otherwise.
 */
bool alignment_coherent(std::string_view target_alignment,
                        std::string_view query_alignment,
                        std::string_view target,
                        std::string_view query);

/**
 * Compute the score of an alignment given the penalties. The target_alignment
 * and query_alignment strings must have the same size, otherwise the behavior
 * is undefined.
 *
 * @param penalties The penalties.
 * @param target_alignment The target alignment.
 * @param query_alignment The query alignment.
 * @return The score of the alignment.
 */
int alignment_score(const Penalties &penalties,
                    std::string_view target_alignment,
                    std::string_view query_alignment);

/**
 * Transform a CIGAR string into an alignment. Return the alignment as a pair
 * <target_alignment, query_alignment>.
 *
 * @param cigar The CIGAR string.
 * @param target The target sequence.
 * @param query The query sequence.
 * @return The alignment corresponding to the CIGAR as a pair <target_alignment,
 * query_alignment>.
 */
std::pair<std::string, std::string> cigar_to_alignment(std::string_view cigar,
                                                       std::string_view target,
                                                       std::string_view query);

/**
 * Transform an alignment into a CIGAR string. The target_alignment and
 * query_alignment strings must have the same size, otherwise the function will
 * throw an exception.
 *
 * @param target_alignment The target alignment.
 * @param query_alignment The query alignment.
 * @return The CIGAR string corresponding to the alignment.
 */
std::string alignment_to_cigar(std::string_view target_alignment,
                               std::string_view query_alignment);

} // namespace alignment
