#include "penalties.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

Penalties::Penalties(int match, int mismatch, int gape)
    : _type(Type::Linear), _match(match), _mismatch(mismatch), _ins(gape), _del(gape) {}

Penalties::Penalties(int match, int mismatch, int gapo, int gape)
    : _type(Type::Affine),
      _match(match),
      _mismatch(mismatch),
      _gapo(gapo),
      _ins(gape),
      _del(gape) {}

Penalties::Penalties(int match, int mismatch, int gapo, int gape, int gapo2, int gape2)
    : _type(Type::DualAffine),
      _match(match),
      _mismatch(mismatch),
      _gapo(gapo),
      _ins(gape),
      _del(gape),
      _gapo2(gapo2),
      _ins2(gape2),
      _del2(gape2) {}
