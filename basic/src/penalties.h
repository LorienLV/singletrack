#pragma once

#include <span>
#include <string>
#include <vector>

/**
 * A class containing the penalties for the alignment algorithm. The objective
 * is always to minimize the score.
 *
 */

class Penalties {
public:
    enum Type {
        Linear,
        Affine,
        DualAffine
    };

    /**
     * Create a gap-linear penalties object.
     *
     * @param match The match score.
     * @param mismatch The mismatch score.
     * @param gape The gap extension penalty.
     */
    Penalties(int match,
              int mismatch,
              int gape);

    /**
     * Create a gap-affine penalties object.
     *
     * @param match The match score.
     * @param mismatch The mismatch score.
     * @param gapo The gap open penalty.
     * @param gape The gap extension penalty.
     */
    Penalties(int match,
              int mismatch,
              int gapo,
              int gape);

    /**
     * Create a dual gap-affine penalties object.
     *
     * @param match The match score.
     * @param mismatch The mismatch score.
     * @param gapo The first gap open penalty.
     * @param gape The first gap extension penalty.
     * @param gapo2 The second gap open penalty.
     * @param gape2 The second gap extension penalty.
     */
    Penalties(int match,
              int mismatch,
              int gapo,
              int gape,
              int gapo2,
              int gape2);

    /**
     * Get the type of penalties.
     *
     * @return The type of penalties.
     */
    Type type() const { return type_; }

    /**
     * Get the match score.
     *
     * @return The match score.
     */
    int match() const { return match_; }

    /**
     * Get the mismatch score.
     *
     * @return The mismatch score.
     */
    int mismatch() const { return mismatch_; }

    /**
     * Get the gap open penalty.
     *
     * @return The gap open penalty.
     */
    int gape() const { return ins_; }

    /**
     * Get the gap extension penalty. If type is Linear, this function returns 0.
     *
     * @return The gap extension penalty.
     */
    int gapo() const { return gapo_; }

    /**
     * Get the gap extension penalty for dual gap-affine. If type is Linear or
     * Affine, this function returns 0.
     *
     * @return The gap extension penalty for dual gap-affine.
     */
    int gape2() const { return ins2_; }

    /**
     * Get the gap open penalty for dual gap-affine. If type is Linear or Affine,
     * this function returns 0.
     *
     * @return The gap open penalty for dual gap-affine.
     */
    int gapo2() const { return gapo2_; }

    /**
     * Get the score of the substitution between two characters.
     *
     * @param t The target character.
     * @param q The query character.
     */
    int subs(char t, char q) const {
        return (t == q) ? match_ : mismatch_;
    }

private:
    Type type_;

    int match_ = 0;
    int mismatch_ = 0;

    int gapo_ = 0; // Gap open.
    int ins_ = 0; // Insertion.
    int del_ = 0; // Deletion.

    int gapo2_ = 0; // Gap open for dual gap-affine.
    int ins2_ = 0; // Insertion for dual gap-affine.
    int del2_ = 0; // Deletion for dual gap-affine.
};
