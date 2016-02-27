// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_GAUSSIAN_H_
#define DIST_GAUSSIAN_H_


#include <cmath>
#include <cstdlib>
#include <iostream>

#include <boost/math/distributions.hpp>
#include <thrust/random.h>

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
  __host__ __device__
  Gaussian(): Gaussian(0, (Real)1.) {}
  __host__ __device__
  Gaussian(Real m, Real var): mean(m), var(var) { }

};

// ------------------------------------------------------------------------------
// Below defines Gaussian related generic functions

///
/// \brief Return mean
///
__host__ __device__
inline double getMean(const Gaussian &dist)
{
    return dist.mean;
}

///
/// \brief Return variance
///
__host__ __device__
inline double getVar(const Gaussian &dist)
{
    return (double) dist.var;
}

///
/// \brief Return PDF of x
///
__host__ __device__
inline double getPdf(const Gaussian &dist, const double x)
{
    if (dist.var==0) {
        return ( fabs(x-dist.mean) < EPS )? 1.: 0;
    }
    return exp( -0.5 * pow(x-dist.mean, 2) / dist.var ) / sqrt(2. * dist.var * M_PI);
}

///
/// \brief Return a random sample
///
__host__
inline double getSample(const Gaussian &dist)
{
    return box_muller((double)dist.mean, (double)sqrt(dist.var) );
}

///
/// \brief Return a random sample using random engine
///
__host__ __device__
inline double getSample(const Gaussian &dist, thrust::default_random_engine &rng)
{
  thrust::random::normal_distribution<double> ndist(dist.mean, sqrt(dist.var) );
  return ndist(rng);
}

#if 0
namespace detail {
  // Returns the erf() of a value (not super precice, but ok)
  // ref: http://www.cplusplus.com/forum/beginner/62864/
  // ref: http://math.stackexchange.com/questions/97/how-to-accurately-calculate-the-error-function-erfx-with-a-computer
  // seems wrong. Don't use.
  __host__ __device__
  inline double erf(double x)
  {
   double y = 1.0 / ( 1.0 + 0.3275911 * x);
   return 1 - (((((
          + 1.061405429  * y
          - 1.453152027) * y
          + 1.421413741) * y
          - 0.284496736) * y
          + 0.254829592) * y)
          * exp (-x * x);
  }
}
#endif

///
/// \brief Return CDF of x
///
/// The underlining erf() function uses Cuda implementation
///
__host__ __device__
inline double getCdf(const Gaussian &dist, double x)
{
  if (dist.var==0) {
    return x >= dist.mean ? 1 : 0;
  }
  //return 0.5 * (1 + boost::math::erf((x - dist.mean) / (sqrt(2.*dist.var))));
  return 0.5 * (1 + erf((x - dist.mean) / (sqrt(2.*dist.var))));
}

///
/// \brief Return CDF of x
///
__host__
inline double getCdfPrecise(const Gaussian &dist, double x)
{
  if (dist.var==0) {
    return x >= dist.mean ? 1 : 0;
  }
  // TODO: need to implement on our own for Cuda to work
  boost::math::normal_distribution<double> normal (dist.mean, sqrt(dist.var) );
  return boost::math::cdf<>(normal, x);
}

///
/// \brief Print itself
///
__host__
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
__host__ __device__
inline Gaussian& operator-(Gaussian &x)
{
    x.mean = -x.mean;
    return x;
}

///
/// \brief random variable +=
///
__host__ __device__
inline Gaussian& operator+=(Gaussian &x, const Gaussian& rhs) {
    x.mean += rhs.mean;
    x.var += rhs.var;
    return x;
}

///
/// \brief random variable += with scalar
///
__host__ __device__
inline Gaussian& operator+=(Gaussian &x, const double r) {
    x.mean += r;
    return x;
}

///
/// \brief random variable *= with scalar
///
__host__ __device__
inline Gaussian& operator*=(Gaussian &x, const double r) {
    x.mean *= r;
    x.var *= r*r;
    return x;
}

}  // namespace dist
}  // namespace edda

#endif  // DIST_GAUSSIAN_H_
