#pragma once

#include <span>
#include <string>
#include <vector>

/**
 * A class containing the penalties for the alignment algorithm.
 *
 */

class Penalties {
public:
    enum Type {
        Linear,
        Affine,
        DualAffine
    };

    // TODO: Linear
    Penalties(int match,
              int mismatch,
              int gape);

    // TODO: Affine
    Penalties(int match,
              int mismatch,
              int gapo,
              int gape);

    // TODO: Dual affine
    Penalties(int match,
              int mismatch,
              int gapo,
              int gape,
              int gapo2,
              int gape2);

    // TODO:
    Type type() const { return _type; }

    // TODO:
    int match() const { return _match; }

    // TODO:
    int mismatch() const { return _mismatch; }

    // TODO:
    int gape() const { return _ins; }

    // TODO:
    int gapo() const { return _gapo; }

    // TODO:
    int gape2() const { return _ins2; }

    // TODO:
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
