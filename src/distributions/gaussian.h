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

///
/// \brief Defines a simple Gaussian distribution storage, where mean and variance are stored together.
///
template <class Real = float>
class EDDA_EXPORT SimpleGaussianStorage {
  Real mean_, var_;
public:
  inline Real &mean() { return mean_; }
  inline Real &var() { return var_; }
  inline const Real &mean() const { return mean_; }
  inline const Real &var() const { return var_; }
};

// ------------------------------------------------------------------------------
///
/// \brief Defines a Gaussian class
///
template<typename Real = float, class Storage=SimpleGaussianStorage<Real> >
class EDDA_EXPORT Gaussian: public ContinuousDistribution, public Storage {
public:
  // constructor
  explicit Gaussian(): Gaussian(0, (Real)1.) {}
  explicit Gaussian(Real m, Real var) { this->mean()=m; this->var()=var; }
};

// ------------------------------------------------------------------------------
// Below defines Gaussian related generic functions

///
/// \brief Return mean
///
template<typename Real, class Storage>
inline double getMean(const Gaussian<Real, Storage> &dist)
{
    return dist.mean();
}

///
/// \brief Return variance
///
template<typename Real, class Storage>
inline double getVar(const Gaussian<Real, Storage> &dist)
{
    return (double) dist.var();
}

///
/// \brief Return PDF of x
///
template<typename Real, class Storage>
inline double getPdf(const Gaussian<Real, Storage> &dist, const double x)
{
    if (dist.var()==0) {
        return ( fabs(x-dist.mean()) < EPS )? 1.: 0;
    }
    return exp( -0.5 * pow(x-dist.mean(), 2) / dist.var() ) / sqrt(2. * dist.var() * M_PI);
}

///
/// \brief Return a sample
///
template<typename Real, class Storage>
inline double getSample(const Gaussian<Real, Storage> &dist)
{
    return box_muller((double)dist.mean(), sqrt(dist.var()) );
}

///
/// \brief Return CDF of x
///
template<typename Real, class Storage>
inline double getCdf(const Gaussian<Real, Storage> &dist, double x)
{
  if (dist.var()==0) {
    return x > dist.mean() ? 1 : 0;
  }
  // TODO: need to implement on our own
  boost::math::normal_distribution<double> normal (dist.mean(), sqrt(dist.var()) );
  return boost::math::cdf<>(normal, x);
}

///
/// \brief Print itself
///
template<typename Real, class Storage>
std::ostream& operator<<(std::ostream& os, const Gaussian<Real, Storage>& dist)
{
    os <<  "<Gaussian: mean=" << getMean(dist) << ", variance=" << getVar(dist) << ">" ;
    return os;
}

// ------------------------------------------------------------------------------
// Below defines Gaussian related arithmetics

///
/// \brief random variable with unary -
///
template<typename Real, class Storage>
inline Gaussian<Real, Storage>& operator-(Gaussian<Real, Storage> &x)
{
    x.mean() = -x.mean();
    return x;
}

///
/// \brief random variable +=
///
template<typename Real, class Storage>
inline Gaussian<Real, Storage>& operator+=(Gaussian<Real, Storage> &x, const Gaussian<Real, Storage>& rhs) {
    x.mean() += rhs.mean();
    x.var() += rhs.var();
    return x;
}

///
/// \brief random variable += with scalar
///
template<typename Real, class Storage>
inline Gaussian<Real, Storage>& operator+=(Gaussian<Real, Storage> &x, const double r) {
    x.mean() += r;
    return x;
}

///
/// \brief random variable *= with scalar
///
template<typename Real, class Storage>
inline Gaussian<Real, Storage>& operator*=(Gaussian<Real, Storage> &x, const double r) {
    x.mean() *= r;
    x.var() *= r*r;
    return x;
}

}  // namespace dist
}  // namespace edda

#endif  // DIST_GAUSSIAN_H_
