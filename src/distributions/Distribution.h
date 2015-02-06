#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

// Copyright (c) 2014 The EDDA Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

// Distribution:
// Define generic functions, which can be overriden by specific template classes

#include <typeinfo>
#include "edda_export.h"

namespace edda { namespace dist {

// a dummy distribution for demostrantion
// class Real: this is only for data storage stype, which is float as default to reduce the size.  Otherwise in general the input/output types should be double

template <class Real = float>
class EDDA_EXPORT Distribution {
public:
    // get probability
    inline Real getProb(const Real x) const { return 0; }
    inline Real getMean() const { return 0; }
    inline Real getStd() const { return 0; }
    inline Real getVariance() const { return 0; }
    inline Real getSample() const { return 0; }

    // random variable +=
    inline Distribution& operator+=(const Distribution& rhs) { return *this; }
    // random variable -=
    inline Distribution& operator-=(const Distribution& rhs) { return *this; }
    // // random variable +=
    inline Distribution& operator+=(const double r) { return *this;  }
    // random variable *=
    inline Distribution& operator*=(const double r) { return *this; }
};


// random variable +
template<class T>
inline T operator+(const T& lhs, const T& rhs) {
    T h(lhs);
    return h += rhs;
}

// random variable -
template<class T>
inline T operator-(const T& lhs, const T& rhs) {
    T h(lhs);
    return h -= rhs;
}

// cdf
template <class T>
inline double cdf(const T& t, double x)
{
    std::cout << "Generic computation of cdf: Not implemented" << std::endl;
    return 0;
}


} } // namespace dist, edda


#endif // DISTRIBUTION_H
