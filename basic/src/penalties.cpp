#include "penalties.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

Penalties::Penalties(int match, int mismatch, int gape)
    : type_(Type::Linear), match_(match), mismatch_(mismatch), ins_(gape), del_(gape) {}

Penalties::Penalties(int match, int mismatch, int gapo, int gape)
    : type_(Type::Affine),
      match_(match),
      mismatch_(mismatch),
      gapo_(gapo),
      ins_(gape),
      del_(gape) {}

Penalties::Penalties(int match, int mismatch, int gapo, int gape, int gapo2, int gape2)
    : type_(Type::DualAffine),
      match_(match),
      mismatch_(mismatch),
      gapo_(gapo),
      ins_(gape),
      del_(gape),
      gapo2_(gapo2),
      ins2_(gape2),
      del2_(gape2) {}
