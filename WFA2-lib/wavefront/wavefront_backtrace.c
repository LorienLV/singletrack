/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: WaveFront-Alignment module for backtracing alignments
 */

#include "utils/commons.h"
#include "wavefront_backtrace.h"

/*
 * Wavefront type
 */
#define BACKTRACE_TYPE_BITS                   4 // 4-bits for piggyback
#define BACKTRACE_TYPE_MASK 0x000000000000000Fl // Extract mask

#define BACKTRACE_PIGGYBACK_SET(offset,backtrace_type) \
  (( ((int64_t)(offset)) << BACKTRACE_TYPE_BITS) | backtrace_type)

#define BACKTRACE_PIGGYBACK_GET_TYPE(offset) \
  ((offset) & BACKTRACE_TYPE_MASK)
#define BACKTRACE_PIGGYBACK_GET_OFFSET(offset) \
  ((offset) >> BACKTRACE_TYPE_BITS)

typedef enum {
  backtrace_M       = 9,
  backtrace_D2_ext  = 8,
  backtrace_D2_open = 7,
  backtrace_D1_ext  = 6,
  backtrace_D1_open = 5,
  backtrace_I2_ext  = 4,
  backtrace_I2_open = 3,
  backtrace_I1_ext  = 2,
  backtrace_I1_open = 1,
} backtrace_type;

static const uint64_t matches_lut = 0x4D4D4D4D4D4D4D4Dul; // Matches LUT = "MMMMMMMM"
static const uint64_t mismatches_lut = 0x5858585858585858ul; // Mismatches LUT = "XXXXXXXX"
static const uint64_t insertions_lut = 0x4949494949494949ul; // Insertions LUT = "IIIIIIII"
static const uint64_t deletions_lut = 0x4444444444444444ul; // Deletions LUT = "DDDDDDDD"

/**
 * Add @p num operations of type @p op to the CIGAR. @p lut_op represents the
 * 64-bit representation of the operation (8 times the same operation). This is
 * used when @p num is greater than or equal to 8, in which case the operation
 * is added in blocks of 8 operations to the CIGAR.
 *
 * @param cigar The CIGAR.
 * @param op The character representing the operation to add (e.g., 'M', 'X',
 * 'I', 'D').
 * @param lut_op The 64-bit representation of @p op, i.e., the same operation
 * repeated 8 times (e.g., 0x4D4D4D4D4D4D4D4Dul for 'M').
 * @param num The number of instances of the operation to add to the CIGAR.
 */
void wavefront_backtrace_add_nop_to_cigar(cigar_t *const cigar,
                                          const char op,
                                          const uint64_t lut_op,
                                          int num) {
  char* operations = cigar->operations + cigar->begin_offset;
  // Update offset first
  cigar->begin_offset -= num;
  // Blocks of 8-operations
  while (num >= 8) {
    operations -= 8;
    *((uint64_t*)(operations+1)) = lut_op;
    num -= 8;
  }
  // Remaining operations
  for (int i = 0; i < num; ++i) {
    *operations = op;
    --operations;
  }
}

/**
 * Add @p num operations of type @p lut_op (8 repetitions of a character) to the
 * CIGAR. If @p num is not multiple of 8, then @p num + (8 - (num % 8))
 * operations will be added to the CIGAR, but the begin_offset will be updated
 * using @p num. In order to do this, the cigar needs a padding of at least 7
 * elements.
 *
 * @param cigar The CIGAR.
 * @param lut_op The 64-bit representation of an operation, i.e., the same
 * operation repeated 8 times (e.g., 0x4D4D4D4D4D4D4D4Dul for 'M').
 * @param num The number of instances of the operation to add to the CIGAR.
 */
void wavefront_backtrace_add_lut_to_cigar(cigar_t *const cigar,
                                          const uint64_t lut_op,
                                          int num) {
  uint64_t* operations = (uint64_t*)(cigar->operations + cigar->begin_offset + 1);
  // Update offset first
  cigar->begin_offset -= num;
  // Blocks of 8-operations. We may add more than num operations if num is
  // not a multiple of 8. We need enough space for this.
  while (num > 0) {
    --operations;
    *operations = lut_op;
    num -= 8;
  }
}

/*
 * Backtrace Trace Patch Match/Mismsmatch
 */
int64_t wavefront_backtrace_misms(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo <= k &&
      k <= mwavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k]+1,backtrace_M);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
void wavefront_backtrace_matches(
    wavefront_aligner_t* const wf_aligner,
    const int k,
    wf_offset_t offset,
    int num_matches,
    cigar_t* const cigar) {
  wavefront_backtrace_add_nop_to_cigar(cigar, 'M', matches_lut, num_matches);
}
/*
 * Backtrace Trace Patch Deletion
 */
int64_t wavefront_backtrace_del1_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo <= k+1 &&
      k+1 <= mwavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k+1],backtrace_D1_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_del2_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo <= k+1 &&
      k+1 <= mwavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k+1],backtrace_D2_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_del1_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const d1wavefront = wf_aligner->wf_components.d1wavefronts[score];
  if (d1wavefront != NULL &&
      d1wavefront->lo <= k+1 &&
      k+1 <= d1wavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(d1wavefront->offsets[k+1],backtrace_D1_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_del2_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const d2wavefront = wf_aligner->wf_components.d2wavefronts[score];
  if (d2wavefront != NULL &&
      d2wavefront->lo <= k+1 &&
      k+1 <= d2wavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(d2wavefront->offsets[k+1],backtrace_D2_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
/*
 * Backtrace Trace Patch Insertion
 */
int64_t wavefront_backtrace_ins1_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo <= k-1 &&
      k-1 <= mwavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k-1]+1,backtrace_I1_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_ins2_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo <= k-1 &&
      k-1 <= mwavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k-1]+1,backtrace_I2_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_ins1_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const i1wavefront = wf_aligner->wf_components.i1wavefronts[score];
  if (i1wavefront != NULL &&
      i1wavefront->lo <= k-1 &&
      k-1 <= i1wavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(i1wavefront->offsets[k-1]+1,backtrace_I1_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_ins2_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const i2wavefront = wf_aligner->wf_components.i2wavefronts[score];
  if (i2wavefront != NULL &&
      i2wavefront->lo <= k-1 &&
      k-1 <= i2wavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(i2wavefront->offsets[k-1]+1,backtrace_I2_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
/*
 * Backtrace wavefronts
 */
void wavefront_backtrace_linear(
    wavefront_aligner_t* const wf_aligner,
    const int alignment_score,
    const int alignment_k,
    const wf_offset_t alignment_offset) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  const wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  const distance_metric_t distance_metric = penalties->distance_metric;
  // Prepare cigar
  cigar_t* const cigar = wf_aligner->cigar;
  cigar_clear(cigar);
  cigar->end_offset = cigar->max_operations - 1;
  cigar->begin_offset = cigar->max_operations - 2;
  cigar->operations[cigar->end_offset] = '\0';
  // Compute starting location
  int score = alignment_score;
  int k = alignment_k;
  int h = WAVEFRONT_H(alignment_k,alignment_offset);
  int v = WAVEFRONT_V(alignment_k,alignment_offset);
  wf_offset_t offset = alignment_offset;
  // Account for ending insertions/deletions
  if (v < pattern_length) {
    int i = pattern_length - v;
    while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --i;};
  }
  if (h < text_length) {
    int i = text_length - h;
    while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --i;};
  }
  // Trace the alignment back
  while (v > 0 && h > 0 && score > 0) {
    // Compute scores
    const int mismatch = score - penalties->mismatch;
    const int gap_open1 = score - penalties->gap_opening1;
    // Compute source offsets
    const int64_t misms = (distance_metric != indel) ?
        wavefront_backtrace_misms(wf_aligner,mismatch,k) :
        WAVEFRONT_OFFSET_NULL;
    const int64_t ins = wavefront_backtrace_ins1_open(wf_aligner,gap_open1,k);
    const int64_t del = wavefront_backtrace_del1_open(wf_aligner,gap_open1,k);
    const int64_t max_all = MAX(misms,MAX(ins,del));
    // Check source score
    if (max_all < 0) break; // No source
    // Traceback Matches
    const int max_offset = BACKTRACE_PIGGYBACK_GET_OFFSET(max_all);
    const int num_matches = offset - max_offset;
    wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
    offset = max_offset;
    // Update coordinates
    v = WAVEFRONT_V(k,offset);
    h = WAVEFRONT_H(k,offset);
    if (v <= 0 || h <= 0) break;
    // Traceback Operation
    const backtrace_type backtrace_type = BACKTRACE_PIGGYBACK_GET_TYPE(max_all);
    switch (backtrace_type) {
      case backtrace_M:
        score = mismatch;
        cigar->operations[(cigar->begin_offset)--] = 'X';
        --offset;
        break;
      case backtrace_I1_open:
        score = gap_open1;
        cigar->operations[(cigar->begin_offset)--] = 'I';
        --k; --offset;
        break;
      case backtrace_D1_open:
        score = gap_open1;
        cigar->operations[(cigar->begin_offset)--] = 'D';
        ++k;
        break;
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.4\n");
        exit(1);
        break;
    }
    // Update coordinates
    v = WAVEFRONT_V(k,offset);
    h = WAVEFRONT_H(k,offset);
  }
  // Account for last operations
  if (v > 0 && h > 0) {
    // Account for beginning series of matches
    const int num_matches = MIN(v,h);
    wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
    v -= num_matches;
    h -= num_matches;
  }
  // Account for beginning insertions/deletions
  while (v > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --v;};
  while (h > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --h;};
  // Set CIGAR
  ++(cigar->begin_offset);
  cigar->score = alignment_score;
}
void wavefront_backtrace_affine(
    wavefront_aligner_t* const wf_aligner,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int alignment_score,
    const int alignment_k,
    const wf_offset_t alignment_offset) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  const wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  const distance_metric_t distance_metric = penalties->distance_metric;
  // Prepare cigar
  cigar_t* const cigar = wf_aligner->cigar;
  cigar_clear(cigar);
  cigar->end_offset = cigar->max_operations - 1;
  cigar->begin_offset = cigar->max_operations - 2;
  cigar->operations[cigar->end_offset] = '\0';
  // Compute starting location
  affine2p_matrix_type matrix_type = component_end;
  int score = alignment_score;
  int k = alignment_k;
  int h = WAVEFRONT_H(alignment_k,alignment_offset);
  int v = WAVEFRONT_V(alignment_k,alignment_offset);
  wf_offset_t offset = alignment_offset;
  // Account for ending insertions/deletions
  if (component_end == affine2p_matrix_M) { // ends-free
    if (v < pattern_length) {
      int i = pattern_length - v;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --i;};
    }
    if (h < text_length) {
      int i = text_length - h;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --i;};
    }
  }
  // Trace the alignment back
  while (v > 0 && h > 0 && score > 0) {
    // Compute scores
    const int mismatch = score - penalties->mismatch;
    const int gap_open1 = score - penalties->gap_opening1 - penalties->gap_extension1;
    const int gap_open2 = score - penalties->gap_opening2 - penalties->gap_extension2;
    const int gap_extend1 = score - penalties->gap_extension1;
    const int gap_extend2 = score - penalties->gap_extension2;
    // Compute source offsets
    int64_t max_all;
    switch (matrix_type) {
      case affine2p_matrix_M: {
        const int64_t misms = wavefront_backtrace_misms(wf_aligner,mismatch,k);
        const int64_t ins1_open = wavefront_backtrace_ins1_open(wf_aligner,gap_open1,k);
        const int64_t ins1_ext  = wavefront_backtrace_ins1_ext(wf_aligner,gap_extend1,k);
        const int64_t max_ins1 = MAX(ins1_open,ins1_ext);
        const int64_t del1_open = wavefront_backtrace_del1_open(wf_aligner,gap_open1,k);
        const int64_t del1_ext  = wavefront_backtrace_del1_ext(wf_aligner,gap_extend1,k);
        const int64_t max_del1 = MAX(del1_open,del1_ext);
        if (distance_metric == gap_affine) {
          max_all = MAX(misms,MAX(max_ins1,max_del1));
          break;
        }
        const int64_t ins2_open = wavefront_backtrace_ins2_open(wf_aligner,gap_open2,k);
        const int64_t ins2_ext  = wavefront_backtrace_ins2_ext(wf_aligner,gap_extend2,k);
        const int64_t max_ins2 = MAX(ins2_open,ins2_ext);
        const int64_t del2_open = wavefront_backtrace_del2_open(wf_aligner,gap_open2,k);
        const int64_t del2_ext  = wavefront_backtrace_del2_ext(wf_aligner,gap_extend2,k);
        const int64_t max_del2 = MAX(del2_open,del2_ext);
        const int64_t max_ins = MAX(max_ins1,max_ins2);
        const int64_t max_del = MAX(max_del1,max_del2);
        max_all = MAX(misms,MAX(max_ins,max_del));
        break;
      }
      case affine2p_matrix_I1: {
        const int64_t ins1_open = wavefront_backtrace_ins1_open(wf_aligner,gap_open1,k);
        const int64_t ins1_ext  = wavefront_backtrace_ins1_ext(wf_aligner,gap_extend1,k);
        max_all = MAX(ins1_open,ins1_ext);
        break;
      }
      case affine2p_matrix_I2: {
        const int64_t ins2_open = wavefront_backtrace_ins2_open(wf_aligner,gap_open2,k);
        const int64_t ins2_ext  = wavefront_backtrace_ins2_ext(wf_aligner,gap_extend2,k);
        max_all = MAX(ins2_open,ins2_ext);
        break;
      }
      case affine2p_matrix_D1: {
        const int64_t del1_open = wavefront_backtrace_del1_open(wf_aligner,gap_open1,k);
        const int64_t del1_ext  = wavefront_backtrace_del1_ext(wf_aligner,gap_extend1,k);
        max_all = MAX(del1_open,del1_ext);
        break;
      }
      case affine2p_matrix_D2: {
        const int64_t del2_open = wavefront_backtrace_del2_open(wf_aligner,gap_open2,k);
        const int64_t del2_ext  = wavefront_backtrace_del2_ext(wf_aligner,gap_extend2,k);
        max_all = MAX(del2_open,del2_ext);
        break;
      }
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.1\n");
        exit(1);
        break;
    }
    // Check source score
    if (max_all < 0) break; // No source
    // Traceback matches
    if (matrix_type == affine2p_matrix_M) {
      const int max_offset = BACKTRACE_PIGGYBACK_GET_OFFSET(max_all);
      const int num_matches = offset - max_offset;
      wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
      offset = max_offset;
      // Update coordinates
      v = WAVEFRONT_V(k,offset);
      h = WAVEFRONT_H(k,offset);
      if (v <= 0 || h <= 0) break;
    }
    // Traceback operation
    const backtrace_type backtrace_type = BACKTRACE_PIGGYBACK_GET_TYPE(max_all);
    switch (backtrace_type) {
      case backtrace_M:
        score = mismatch;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_I1_open:
        score = gap_open1;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_I1_ext:
        score = gap_extend1;
        matrix_type = affine2p_matrix_I1;
        break;
      case backtrace_I2_open:
        score = gap_open2;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_I2_ext:
        score = gap_extend2;
        matrix_type = affine2p_matrix_I2;
        break;
      case backtrace_D1_open:
        score = gap_open1;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_D1_ext:
        score = gap_extend1;
        matrix_type = affine2p_matrix_D1;
        break;
      case backtrace_D2_open:
        score = gap_open2;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_D2_ext:
        score = gap_extend2;
        matrix_type = affine2p_matrix_D2;
        break;
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.2\n");
        exit(1);
        break;
    }
    switch (backtrace_type) {
      case backtrace_M:
        cigar->operations[(cigar->begin_offset)--] = 'X';
        --offset;
        break;
      case backtrace_I1_open:
      case backtrace_I1_ext:
      case backtrace_I2_open:
      case backtrace_I2_ext:
        cigar->operations[(cigar->begin_offset)--] = 'I';
        --k; --offset;
        break;
      case backtrace_D1_open:
      case backtrace_D1_ext:
      case backtrace_D2_open:
      case backtrace_D2_ext:
        cigar->operations[(cigar->begin_offset)--] = 'D';
        ++k;
        break;
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.3\n");
        exit(1);
        break;
    }
    // Update coordinates
    v = WAVEFRONT_V(k,offset);
    h = WAVEFRONT_H(k,offset);
  }
  // Account for last operations
  if (matrix_type == affine2p_matrix_M) {
    if (v > 0 && h > 0) {
      // Account for beginning series of matches
      const int num_matches = MIN(v,h);
      wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
      v -= num_matches;
      h -= num_matches;
    }
    // Account for beginning insertions/deletions
    while (v > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --v;};
    while (h > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --h;};
  } else {
    // DEBUG
    if (v != 0 || h != 0 || (score != 0 && penalties->match == 0)) {
      fprintf(stderr,"[WFA::Backtrace] I?/D?-Beginning backtrace error\n");
      fprintf(stderr,">%.*s\n",pattern_length,sequences->pattern);
      fprintf(stderr,"<%.*s\n",text_length,sequences->text);
      exit(-1);
    }
  }
  // Set CIGAR
  ++(cigar->begin_offset);
  cigar->score = alignment_score;
}



#if 1
/**
 * Check that the cigar produced by the M-only backtrace is coherent, i.e.,
 * it has the expected score and matches and mismatches are consistent with
 * the pattern and text.
 * 
 * @param penalties The penalties.
 * @param cigar The cigar produced by the M-only backtrace.
 * @param sequences The sequences.
 * @param expected_score The expected score of the alignment.
 * @param affine2p True if the alignment is affine2p, false otherwise.
 * @return True if the cigar is coherent, false otherwise.
 */
static bool
check_cigar_backtrace_affine_m_only(const wavefront_penalties_t* const penalties,
                                    const cigar_t* const cigar,
                                    const wavefront_sequences_t* const sequences,
                                    const int expected_score,
                                    const bool affine2p) {

  const char* const pattern = sequences->pattern;
  const char* const text = sequences->text;

  int score = 0;
  int score2 = 0;   // For Dual affine.

  int v = 0;
  int h = 0;

  char prev_op = ' ';

  for (int i = cigar->begin_offset; i < cigar->end_offset; ++i) {
    const char op = cigar->operations[i];

    if ((prev_op == 'I' || prev_op == 'D') && op != prev_op && affine2p) {
      // We have finished a chain of gaps, get the best score.
      score = MIN(score, score2);
    }

    if (op == 'M') {
      if (pattern[v] != text[h]) {
        return false;
      }

      score += penalties->match;

      ++v;
      ++h;
    }
    else if (op == 'X') {
      if (pattern[v] == text[h]) {
        return false;
      }

      score += penalties->mismatch;

      ++v;
      ++h;
    }
    else if (op == 'I') {
      if (prev_op != 'I') {
        score2 = score + penalties->gap_opening2;
        score += penalties->gap_opening1;
      }

      score += penalties->gap_extension1;
      score2 += penalties->gap_extension2;

      ++h;
    }
    else if (op == 'D') {
      if (prev_op != 'D') {
        score2 = score + penalties->gap_opening2;
        score += penalties->gap_opening1;
      }

      score += penalties->gap_extension1;
      score2 += penalties->gap_extension2;

      ++v;
    }
    else {
      return false;
    }

    prev_op = op;
  }

  // In case the sequence ends with a gap.
  if ((prev_op == 'I' || prev_op == 'D') && affine2p) {
    score = MIN(score, score2);
  }

  if (score != expected_score) {
    return false;
  }

  return true;
}

#endif

/**
 * Retrieve the cigar of the alignment for gap-affine and dual gap-affine
 * using exclusively the information in the M matrix (M wavefronts). This is
 * slightly slower than the general backtracking, but it enables storing only
 * the M matrix when performing the alignment (I1, I2, D1 and D2 only need a
 * small scope).
 *
 * @param wf_aligner The wavefront aligner.
 * @param component_begin unused.
 * @param component_end The component (matrix) where the backtracking starts.
 * @param alignment_score The score of the alignment.
 * @param alignment_k The diagonal that contains the cell (N, M), where the
 * backtracking starts.
 * @param alignment_offset The offset of the cell (N, M), where the backtracking
 * starts.
 */
void wavefront_backtrace_affine_m_only(
    wavefront_aligner_t* const wf_aligner,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int alignment_score,
    const int alignment_k,
    const wf_offset_t alignment_offset) {

  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  const wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  const distance_metric_t distance_metric = penalties->distance_metric;

  // TODO: Do this elsewhere.
  // We need padding in order to perform the backwards extend.
  sequences->pattern[-1] = '!';
  sequences->text[-1] = '?';

  // Prepare cigar
  // WARNING: We want padding in the cigar so we can always use LUTs.
  cigar_t* const cigar = wf_aligner->cigar;
  cigar_clear(cigar);
  cigar->end_offset = cigar->max_operations - 1;
  cigar->begin_offset = cigar->max_operations - 2;
  cigar->operations[cigar->end_offset] = '\0';

  // Compute starting location
  affine2p_matrix_type matrix_type = component_end;
  int score = alignment_score;
  int k = alignment_k;
  int h = WAVEFRONT_H(alignment_k,alignment_offset);
  int v = WAVEFRONT_V(alignment_k,alignment_offset);
  wf_offset_t offset = alignment_offset;
  wf_offset_t offset_orig = offset;

  // Account for ending insertions/deletions
  if (component_end == affine2p_matrix_M) { // ends-free
    if (v < pattern_length) {
      int i = pattern_length - v;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --i;};
    }
    if (h < text_length) {
      int i = text_length - h;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --i;};
    }
  }

  bool in_mmatrix = component_end == affine2p_matrix_M;
  int l = 0; // Length of the current chain of insertions or deletions.

  // Trace the alignment back
  while (score != 0) {
    if (in_mmatrix) {
      v = WAVEFRONT_V(k, offset);
      h = WAVEFRONT_H(k, offset);

      int nmatches = 0;

      // TODO: Function to backwards extend?
#if __BYTE_ORDER == __LITTLE_ENDIAN
      // Blocked backwards extend.
      const uint64_t* pattern_blocks = (uint64_t*)(wf_aligner->sequences.pattern + v - 8);
      const uint64_t* text_blocks = (uint64_t*)(wf_aligner->sequences.text + h - 8);

      uint64_t cmp = *pattern_blocks ^ *text_blocks;
      while (__builtin_expect(cmp==0,0)) {
        // Next blocks
        --pattern_blocks;
        --text_blocks;

        nmatches += 8;

        // Compare
        cmp = *pattern_blocks ^ *text_blocks;
      }

      const int equal_left_bits = __builtin_clzl(cmp);
      const int equal_chars = DIV_FLOOR(equal_left_bits, 8);

      nmatches += equal_chars;
#else
      // Char-wise backwards extend.
      const char* pattern_ptr = wf_aligner->sequences.pattern + v - 1;
      const char* text_ptr = wf_aligner->sequences.text + h - 1;

      while (*pattern_ptr == *text_ptr) {
        // Next chars
        --pattern_ptr;
        --text_ptr;

        ++nmatches;
      }
#endif

      const int mismatch = score - penalties->mismatch;
      const wavefront_t* const mwavefront = (mismatch >= 0) ?
        wf_aligner->wf_components.mwavefronts[mismatch]
        : NULL;

      offset_orig = offset - nmatches;

      // If we come from a mismatch, then the backwards extend is correct.
      if (mwavefront != NULL &&
          mwavefront->lo <= k &&
          k <= mwavefront->hi &&
          mwavefront->offsets[k] + 1 == offset_orig) {

        offset = offset_orig - 1;

        score = mismatch;

        wavefront_backtrace_add_lut_to_cigar(cigar, matches_lut, nmatches);
        cigar->operations[(cigar->begin_offset)--] = 'X';
      }
      else {
        // Otherwise, we come from either I1, D2, I2 or D2.
        // Freeze v and h and start searching a path back to M.

        // In WFA, for a given diagonal and score, we only store the furthest
        // reaching offset. We do not know which was the offset prior to
        // extending it. To account for that, we must allow the indel to come
        // at any point between the offset previous performing the backwards
        // extension and the current offset.
        //
        // Example, in the following DP table, where there is a chain of 3
        // matches an insertion can come from any of the positions marked with
        // '>'.
        //
        //      A  A  A
        //  A > M
        //  A    > M
        //  A       > M
        //
        in_mmatrix = false;

        l = 0;
      }
    } else {
      // We are searching a path back to M in I1, D1, I2 and D2.
      // Any path that leads to M is valid.
      ++l;

      const int k_ins = k - l;
      const int k_del = k + l;

      const int indel1 = score - penalties->gap_opening1 -
                         l * penalties->gap_extension1;
      const wavefront_t* const mwavefront1 = (indel1 >= 0) ?
        wf_aligner->wf_components.mwavefronts[indel1] : NULL;

      // I1 path.
      if (mwavefront1 != NULL &&
          mwavefront1->lo <= k_ins &&
          k_ins <= mwavefront1->hi &&
          mwavefront1->offsets[k_ins] + l >= offset_orig &&
          mwavefront1->offsets[k_ins] + l <= offset) {

        const int nmatches = offset - (mwavefront1->offsets[k_ins] + l);

        wavefront_backtrace_add_lut_to_cigar(cigar, matches_lut, nmatches);
        wavefront_backtrace_add_lut_to_cigar(cigar, insertions_lut, l);

        k = k_ins;
        offset = mwavefront1->offsets[k_ins];
        score = indel1;

        in_mmatrix = true;

        continue;
      }

      // D1 path.
      if (mwavefront1 != NULL &&
          mwavefront1->lo <= k_del &&
          k_del <= mwavefront1->hi &&
          mwavefront1->offsets[k_del] >= offset_orig &&
          mwavefront1->offsets[k_del] <= offset) {

        const int nmatches = offset - mwavefront1->offsets[k_del];

        wavefront_backtrace_add_lut_to_cigar(cigar, matches_lut, nmatches);
        wavefront_backtrace_add_lut_to_cigar(cigar, deletions_lut, l);

        k = k_del;
        offset = mwavefront1->offsets[k_del];
        score = indel1;

        in_mmatrix = true;

        continue;
      }

      if (distance_metric != gap_affine_2p) {
        continue;
      }

      const int indel2 = score - penalties->gap_opening2 -
                         l * penalties->gap_extension2;
      const wavefront_t* const mwavefront2 = (indel2 >= 0) ?
        wf_aligner->wf_components.mwavefronts[indel2] : NULL;

      // I2 path.
      if (mwavefront2 != NULL &&
          mwavefront2->lo <= k_ins &&
          k_ins <= mwavefront2->hi &&
          mwavefront2->offsets[k_ins] + l >= offset_orig &&
          mwavefront2->offsets[k_ins] + l <= offset) {

        const int nmatches = offset - (mwavefront2->offsets[k_ins] + l);

        wavefront_backtrace_add_lut_to_cigar(cigar, matches_lut, nmatches);
        wavefront_backtrace_add_lut_to_cigar(cigar, insertions_lut, l);

        k = k_ins;
        offset = mwavefront2->offsets[k_ins];
        score = indel2;

        in_mmatrix = true;

        continue;
      }

      // D2 path.
      if (mwavefront2 != NULL &&
          mwavefront2->lo <= k_del &&
          k_del <= mwavefront2->hi &&
          mwavefront2->offsets[k_del] >= offset_orig &&
          mwavefront2->offsets[k_del] <= offset) {

        const int nmatches = offset - mwavefront2->offsets[k_del];

        wavefront_backtrace_add_lut_to_cigar(cigar, matches_lut, nmatches);
        wavefront_backtrace_add_lut_to_cigar(cigar, deletions_lut, l);

        k = k_del;
        offset = mwavefront2->offsets[k_del];
        score = indel2;

        in_mmatrix = true;

        continue;
      }
    }
  }

  // Account for last operations.
  // TODO: This only works if mismatch, and indels has > 0 score.
  // Check that we can assume this. Otherwise, insert boundary condition to the loop.
  assert(penalties->mismatch > 0);
  if (distance_metric == gap_affine || distance_metric == gap_affine_2p) {
    assert(penalties->gap_extension1 > 0);
  }
  if (distance_metric == gap_affine_2p) {
    assert(penalties->gap_extension2 > 0);
  }

  wavefront_backtrace_add_lut_to_cigar(cigar, matches_lut, offset);
  offset = 0;

  // DEBUG
  if (offset != 0 || k != 0 || (score != 0 && penalties->match == 0)) {
    fprintf(stderr,"[WFA::Backtrace] I?/D?-Beginning backtrace error\n");
    fprintf(stderr,">%.*s\n",pattern_length,sequences->pattern);
    fprintf(stderr,"<%.*s\n",text_length,sequences->text);
    exit(-1);
  }

  // Dirty way of freeing the slabs used by I1, I2, D1 and D2.
  // TODO: This should be done in a better way.
  for (int i = 0; i <= alignment_score; ++i) {
    if (distance_metric == gap_affine) {
      if (wf_aligner->wf_components.i1wavefronts[i]) {
        wf_aligner->wavefront_slab->memory_used -= wavefront_get_size(wf_aligner->wf_components.i1wavefronts[i]);
        wavefront_free(wf_aligner->wf_components.i1wavefronts[i],wf_aligner->wavefront_slab->mm_allocator);
        wf_aligner->wf_components.i1wavefronts[i]->status = wavefront_status_deallocated;
      }
      if (wf_aligner->wf_components.d1wavefronts[i]) {
        wf_aligner->wavefront_slab->memory_used -= wavefront_get_size(wf_aligner->wf_components.d1wavefronts[i]);
        wavefront_free(wf_aligner->wf_components.d1wavefronts[i],wf_aligner->wavefront_slab->mm_allocator);
        wf_aligner->wf_components.d1wavefronts[i]->status = wavefront_status_deallocated;
      }
    }
    if (distance_metric == gap_affine_2p) {
      if (wf_aligner->wf_components.i2wavefronts[i]) {
        wf_aligner->wavefront_slab->memory_used -= wavefront_get_size(wf_aligner->wf_components.i2wavefronts[i]);
        wavefront_free(wf_aligner->wf_components.i2wavefronts[i],wf_aligner->wavefront_slab->mm_allocator);
        wf_aligner->wf_components.i2wavefronts[i]->status = wavefront_status_deallocated;
      }
      if (wf_aligner->wf_components.d2wavefronts[i]) {
        wf_aligner->wavefront_slab->memory_used -= wavefront_get_size(wf_aligner->wf_components.d2wavefronts[i]);
        wavefront_free(wf_aligner->wf_components.d2wavefronts[i],wf_aligner->wavefront_slab->mm_allocator);
        wf_aligner->wf_components.d2wavefronts[i]->status = wavefront_status_deallocated;
      }
    }
  }

  // Set CIGAR
  ++(cigar->begin_offset);
  cigar->score = alignment_score;

#if 0
  const bool ok = check_cigar_backtrace_affine_m_only(penalties,
                                                      cigar,
                                                      sequences,
                                                      alignment_score,
                                                      distance_metric == gap_affine_2p);
  if (!ok) {
    fprintf(stderr, "[WFA::Backtrace] M-only backtrace not coherent\n");
    exit(-1);
  }
#endif
}
/*
 * Backtrace from BT-Buffer (pcigar)
 */
void wavefront_backtrace_pcigar(
    wavefront_aligner_t* const wf_aligner,
    const int alignment_k,
    const int alignment_offset,
    const pcigar_t pcigar_last,
    const bt_block_idx_t prev_idx_last) {
  // Parameters
  wf_backtrace_buffer_t* const bt_buffer =  wf_aligner->wf_components.bt_buffer;
  // Traceback pcigar-blocks
  bt_block_t bt_block_last = {
      .pcigar = pcigar_last,
      .prev_idx = prev_idx_last
  };
  bt_block_t* const init_block = wf_backtrace_buffer_traceback_pcigar(bt_buffer,&bt_block_last);
  // Fetch initial coordinate
  const int init_position_offset = init_block->pcigar;
  wf_backtrace_init_pos_t* const backtrace_init_pos =
      vector_get_elm(bt_buffer->alignment_init_pos,init_position_offset,wf_backtrace_init_pos_t);
  // Unpack pcigar blocks (packed alignment)
  const int begin_v = backtrace_init_pos->v;
  const int begin_h = backtrace_init_pos->h;
  const int end_v = WAVEFRONT_V(alignment_k,alignment_offset);
  const int end_h = WAVEFRONT_H(alignment_k,alignment_offset);
  if (wf_aligner->penalties.distance_metric <= gap_linear) {
    wf_backtrace_buffer_unpack_cigar_linear(
        bt_buffer,&wf_aligner->sequences,
        begin_v,begin_h,end_v,end_h,wf_aligner->cigar);
  } else {
    wf_backtrace_buffer_unpack_cigar_affine(
        bt_buffer,&wf_aligner->sequences,
        begin_v,begin_h,end_v,end_h,wf_aligner->cigar);
  }
}
