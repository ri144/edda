#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

// Copyright (c) 2014 The EDDA Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.


#include "edda_export.h"

namespace edda {

/// Defines a generic distribution class
template<class Child>
class EDDA_EXPORT Distribution {
protected:
    Child child;
public:
    // construction
    explicit Distribution() {}
    explicit Distribution(Child child_) : child( child_ ) {}
    Distribution(const Distribution<Child>& x) { child=x.child; }

    // assign
    Distribution& operator=(const Distribution<Child>& x) { child=x.child; return *this; }

    ~Distribution() {}

    // get probability
    inline double getProb(const double x) const { return child.getProb; }
    inline double getMean() const { return child.getMean(); }
    inline double getVariance() const { return child.getVariance(); }
    inline double getStd() const { return child.getStd(); }

    inline Distribution<Child>& operator+=(const Distribution<Child>& rhs) { child += rhs.child; return *this; }
    inline Distribution<Child>& operator-=(const Distribution<Child>& rhs) { child -= rhs.child; return *this; }
    inline Distribution<Child>& operator+=(const double r) { child += r; return *this; }
    inline Distribution<Child>& operator*=(const double r) { child *= r; return *this; }
};

template<class T>
inline Distribution<T> operator+(const Distribution<T>& lhs, const Distribution<T>& rhs) {
    Distribution<T> h(lhs);
    return h += rhs;
}
template<class T>
inline Distribution<T> operator-(const Distribution<T>& lhs, const Distribution<T>& rhs) {
    Distribution<T> h(lhs);
    return h -= rhs;
}


} // namespace itl


#endif // DISTRIBUTION_H
