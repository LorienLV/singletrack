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
    Type type() const { return _type; }

    /**
     * Get the match score.
     *
     * @return The match score.
     */
    int match() const { return _match; }

    /**
     * Get the mismatch score.
     *
     * @return The mismatch score.
     */
    int mismatch() const { return _mismatch; }

    /**
     * Get the gap open penalty.
     *
     * @return The gap open penalty.
     */
    int gape() const { return _ins; }

    /**
     * Get the gap extension penalty. If type is Linear, this function returns 0.
     *
     * @return The gap extension penalty.
     */
    int gapo() const { return _gapo; }

    /**
     * Get the gap extension penalty for dual gap-affine. If type is Linear or
     * Affine, this function returns 0.
     *
     * @return The gap extension penalty for dual gap-affine.
     */
    int gape2() const { return _ins2; }

    /**
     * Get the gap open penalty for dual gap-affine. If type is Linear or Affine,
     * this function returns 0.
     *
     * @return The gap open penalty for dual gap-affine.
     */
    int gapo2() const { return _gapo2; }

    /**
     * Get the score of the substitution between two characters.
     *
     * @param t The target character.
     * @param q The query character.
     */
    int subs(char t, char q) const {
        return (t == q) ? _match : _mismatch;
    }

private:
    Type _type;

    int _match = 0;
    int _mismatch = 0;

    int _gapo = 0; // Gap open.
    int _ins = 0; // Insertion.
    int _del = 0; // Deletion.

    int _gapo2 = 0; // Gap open for dual gap-affine.
    int _ins2 = 0; // Insertion for dual gap-affine.
    int _del2 = 0; // Deletion for dual gap-affine.
};
