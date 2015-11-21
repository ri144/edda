// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_GAUSSIAN_H_
#define DIST_GAUSSIAN_H_

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <boost/math/distributions.hpp>

#include "edda_export.h"
#include "distributions/distribution.h"
#include "math/statistics.h"

namespace edda {
namespace dist {

/// Defines a simple Gaussian distribution storage
/// Advanced storage allows storing mean and variance in different ways
template <class Real = float>
class EDDA_EXPORT SimpleGaussianStorage {
  Real mean_, var_;
public:
  inline Real &mean() { return mean_; }
  inline Real &var() { return var_; }
  inline const Real &mean() const { return mean_; }
  inline const Real &var() const { return var_; }
};

/// Defines a Gaussian class
template<typename Real = float, class Storage=SimpleGaussianStorage<Real>>
class EDDA_EXPORT Gaussian: public ContinuousDistribution {
  Storage storage;
public:
  // constructor
  explicit Gaussian(): Gaussian(0, 1.) {}
  explicit Gaussian(Real m, Real var) { storage.mean()=m; storage.var()=var; }

  // The implementation of a Gaussian class should have two functions:
  inline Real &mean() { return storage.mean(); }
  inline Real &var() { return storage.var(); }
  inline const Real &mean() const { return storage.mean(); }
  inline const Real &var() const { return storage.var(); }
};

/// Below defines Gaussian related generic functions
template<typename Real, class Storage>
inline double getMean(const Gaussian<Real, Storage> &dist)
{
    return dist.mean();
}

// get variance
template<typename Real, class Storage>
inline double getVar(const Gaussian<Real, Storage> &dist)
{
    return (double) dist.var();
}

// get pdf of x
template<typename Real, class Storage>
inline double getPdf(const Gaussian<Real, Storage> &dist, const double x)
{
    if (dist.var()==0) {
        return ( fabs(x-dist.mean()) < EPS )? 1.: 0;
    }
    return exp( -0.5 * pow(x-dist.mean(), 2) / dist.var() ) / sqrt(2. * dist.var() * M_PI);
}

template<typename Real, class Storage>
inline double getSample(const Gaussian<Real, Storage> &dist)
{
    return box_muller((double)dist.mean(), sqrt(dist.var()) );
}

// random variable with unary -
template<typename Real, class Storage>
inline Gaussian<Real, Storage>& operator-(Gaussian<Real, Storage> &x)
{
    x.mean() = -x.mean();
    return x;
}
// random variable +=
template<typename Real, class Storage>
inline Gaussian<Real, Storage>& operator+=(Gaussian<Real, Storage> &x, const Gaussian<Real, Storage>& rhs) {
    x.mean() += rhs.mean();
    x.var() += rhs.var();
    return x;
}
// random variable += with scalar
template<typename Real, class Storage>
inline Gaussian<Real, Storage>& operator+=(Gaussian<Real, Storage> &x, const double r) {
    x.mean() += r;
    return x;
}
// random variable *= with scalar
template<typename Real, class Storage>
inline Gaussian<Real, Storage>& operator*=(Gaussian<Real, Storage> &x, const double r) {
    x.mean() *= r;
    x.var() *= r*r;
    return x;
}
// get CDF
template<typename Real, class Storage>
inline double getCdf(const Gaussian<Real, Storage> &dist, double x)
{
  // TODO: need to implement on our own
  boost::math::normal_distribution<double> normal (dist.mean(), sqrt(dist.var()) );
  return boost::math::cdf<>(normal, x);
}
// output
template<typename Real, class Storage>
std::ostream& operator<<(std::ostream& os, const Gaussian<Real, Storage>& dist)
{
    os <<  "<Gaussian: mean=" << getMean(dist) << ", variance=" << getVar(dist) << ">" ;
    return os;
}

}  // namespace dist
}  // namespace edda

#endif  // DIST_GAUSSIAN_H_
