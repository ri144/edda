// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

/// Distributions:
/// Define generic functions, which can be overriden by specific template classes
///

#include <iostream>
#include <typeinfo>
#include <iostream>

#include <common.h>
#include <edda_export.h>
#include <core/vector_matrix.h>

namespace edda {
namespace dist {

///
/// \brief The Distribution class is a root class for all distribution-type classes.
///
/// This is useful for applications to identify whether a class is a distribution-type class, by using ENABLE_IF_BASE_OF()
///
class EDDA_EXPORT Distribution {
    void NotImpError() const {std::cerr << "Functions in this class should be overloaded."; throw NotImplementedException(); }
};

class EDDA_EXPORT ContinuousDistribution : public Distribution{
};

class EDDA_EXPORT DiscreteDistribution : public Distribution{
};

// Here are generic functions that can be reused by new distributions if not
// implemented.

///
/// \brief random variable -=
///
template<class T, ENABLE_IF_BASE_OF(T, Distribution) >
inline T operator-=(const T& lhs, const T& rhs) {
    T h(lhs);
    return h += (-rhs);
}

///
/// \brief random variable +
///
template<class T, ENABLE_IF_BASE_OF(T, Distribution) >
inline T operator+(const T& lhs, const T& rhs) {
    T h(lhs);
    return h += rhs;
}

///
/// \brief random variable -
///
template<class T, ENABLE_IF_BASE_OF(T, Distribution) >
inline T operator-(const T& lhs, const T& rhs) {
    T h(lhs);
    return h -= rhs;
}

///
/// \brief random variable *
///
template<class T, ENABLE_IF_BASE_OF(T, Distribution) >
inline T operator*(const T& lhs, const double x) {
    T h(lhs);
    return h *= x;
}

///
/// \brief Compute CDF of a distribution
///
template <class T, ENABLE_IF_BASE_OF(T, Distribution) >
inline double getCdf(const T& dist, double x)
{
    std::cerr << "Generic computation not implemented" << std::endl;
    throw NotImplementedException();
    return 0;
}

///
/// \brief Compute the mean of the distribution
///
template <class T, ENABLE_IF_BASE_OF(T, Distribution) >
double getMean(const T &dist)
{
    std::cerr << "Generic computation not implemented" << std::endl;
    throw NotImplementedException();
    return 0;
}

///
/// \brief Get a Monte-Carlo sample of the distribution. We rely on the specific implementation from each distribution.
///
#if 0
template <class T, ENABLE_IF_BASE_OF(T, Distribution)>
double getSample(const T &dist)
{
    std::cerr << "Generic computation not implemented" << std::endl;
    throw NotImplementedException();
    return 0;
}
#endif

///
/// \brief Get vector mean
///
template <class T, int N, ENABLE_IF_BASE_OF(T, Distribution) >
inline Vector<double, N> getMean(const Tuple<T, N> &dist)
{
    Vector<double,N> result;
    for (int i=0; i<N; i++)
        result[i] = getMean(dist[i]);
    return result;
}

}  // namespace dist
}  // namespace edda

#endif // DISTRIBUTION_H_
