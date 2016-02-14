// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_GAUSSIAN_H_
#define DIST_GAUSSIAN_H_

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <boost/math/distributions.hpp>

#include "distribution.h"
#include "core/statistics.h"

namespace edda {
namespace dist {

// ------------------------------------------------------------------------------
///
/// \brief Defines a Gaussian class
///
struct EDDA_EXPORT Gaussian: public ContinuousDistribution {
  Real mean, var;
  // constructor
  Gaussian(): Gaussian(0, (Real)1.) {}
  Gaussian(Real m, Real var) { this->mean=m; this->var=var; }

};

// ------------------------------------------------------------------------------
// Below defines Gaussian related generic functions

///
/// \brief Return mean
///
inline double getMean(const Gaussian &dist)
{
    return dist.mean;
}

///
/// \brief Return variance
///
inline double getVar(const Gaussian &dist)
{
    return (double) dist.var;
}

///
/// \brief Return PDF of x
///
inline double getPdf(const Gaussian &dist, const double x)
{
    if (dist.var==0) {
        return ( fabs(x-dist.mean) < EPS )? 1.: 0;
    }
    return exp( -0.5 * pow(x-dist.mean, 2) / dist.var ) / sqrt(2. * dist.var * M_PI);
}

///
/// \brief Return a sample
///
inline double getSample(const Gaussian &dist)
{
    return box_muller((double)dist.mean, sqrt(dist.var) );
}

///
/// \brief Return CDF of x
///
inline double getCdf(const Gaussian &dist, double x)
{
  if (dist.var==0) {
    return x > dist.mean ? 1 : 0;
  }
  // TODO: need to implement on our own
  boost::math::normal_distribution<double> normal (dist.mean, sqrt(dist.var) );
  return boost::math::cdf<>(normal, x);
}

///
/// \brief Print itself
///
inline std::ostream& operator<<(std::ostream& os, const Gaussian &dist)
{
    os <<  "<Gaussian: mean=" << getMean(dist) << ", variance=" << getVar(dist) << ">" ;
    return os;
}

// ------------------------------------------------------------------------------
// Below defines Gaussian related arithmetics

///
/// \brief random variable with unary -
///
inline Gaussian& operator-(Gaussian &x)
{
    x.mean = -x.mean;
    return x;
}

///
/// \brief random variable +=
///
inline Gaussian& operator+=(Gaussian &x, const Gaussian& rhs) {
    x.mean += rhs.mean;
    x.var += rhs.var;
    return x;
}

///
/// \brief random variable += with scalar
///
inline Gaussian& operator+=(Gaussian &x, const double r) {
    x.mean += r;
    return x;
}

///
/// \brief random variable *= with scalar
///
inline Gaussian& operator*=(Gaussian &x, const double r) {
    x.mean *= r;
    x.var *= r*r;
    return x;
}

}  // namespace dist
}  // namespace edda

#endif  // DIST_GAUSSIAN_H_
