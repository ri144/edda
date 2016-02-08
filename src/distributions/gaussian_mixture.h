// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef GAUSSIAN_MIXTURE_H
#define GAUSSIAN_MIXTURE_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <boost/math/distributions.hpp>
#include "common.h"
#include "distribution.h"
#include "gaussian.h"
#include "core/statistics.h"

namespace edda {
namespace dist {

struct GMMTuple {
  union{
    struct{Real m,v,w;};  // mean, variance, weight
    Real p[3];
  };
};

// ------------------------------------------------------------------------------
///
/// \brief Defines a Gaussian Mixture class
///
struct EDDA_EXPORT GaussianMixture: public ContinuousDistribution {
  std::vector<GMMTuple> models;

  // constructor
  explicit GaussianMixture() {}
  explicit GaussianMixture(std::vector<GMMTuple> models_) {
    models = models_;
    normalizeWeights();
  }
  void normalizeWeights() {
    Real sum = 0;
    size_t i;
    for (i=0; i<models.size(); i++) {
      sum += models[i].w;
    }
    if (sum == 0) return;
    for (i=0; i<models.size(); i++) {
      models[i].w /= sum;
    }
  }
};

// ------------------------------------------------------------------------------
// Below defines GaussianMixture related generic functions

///
/// \brief Return mean
///
inline double getMean(const GaussianMixture &dist)
{
  double mean = 0;
  for (int i=0; i<dist.models.size(); i++)
  {
    mean += dist.models[i].m * dist.models[i].w;
  }
  return mean;
}

///
/// \brief Return variance
///
inline double getVar(const GaussianMixture &dist)
{
    throw NotImplementedException();
}

///
/// \brief Return PDF of x
///
inline double getPdf(const GaussianMixture &dist, const double x)
{
  double p=0;
  for (int i=0; i<dist.models.size(); i++)
  {
    p += getPdf( Gaussian(dist.models[i].m, dist.models[i].v), x ) * dist.models[i].w;
  }
  return p;
}

///
/// \brief Return a sample
///
inline double getSample(const GaussianMixture &dist)
{
  float ratio = rand() / (float)RAND_MAX;
  float accumulated = 0;
  for (int i=0; i<dist.models.size(); i++)
  {
    accumulated += dist.models[i].w;
    if (ratio < accumulated) {
      return getSample( Gaussian(dist.models[i].m, dist.models[i].v) );
    }
  }
  // return sample from the last model
  return getSample( Gaussian(dist.models[dist.models.size()-1].m, dist.models[dist.models.size()-1].v) );
}

///
/// \brief Return CDF of x
///
inline double getCdf(const GaussianMixture &dist, double x)
{
  double cdf=0;
  for (int i=0; i<dist.models.size(); i++)
  {
    cdf += getCdf(Gaussian(dist.models[i].m, dist.models[i].v), x) * dist.models[i].w;
  }
  return cdf;
}

///
/// \brief Print itself
///
inline std::ostream& operator<<(std::ostream& os, const GaussianMixture &dist)
{
  os <<  "<GaussianMixture(m,v,w):";
  for (int i=0; i<dist.models.size(); i++)
    os << " (" << dist.models[i].m << "," << dist.models[i].v << "," << dist.models[i].w << ")";
  os << ">";
  return os;
}

// ------------------------------------------------------------------------------
// Below defines Gaussian related arithmetics

///
/// \brief random variable with unary -
///
inline GaussianMixture& operator-(GaussianMixture &x)
{
  for (int i=0; i<x.models.size(); i++)
    x.models[i].m = -x.models[i].m;
  return x;
}

///
/// \brief random variable +=
///
inline GaussianMixture& operator+=(GaussianMixture &x, const GaussianMixture& rhs) {
  throw NotImplementedException();
}

///
/// \brief random variable += with scalar
///
inline GaussianMixture& operator+=(GaussianMixture &x, const double r) {
  for (int i=0; i<x.models.size(); i++)
    x.models[i].m += r;
  return x;
}

///
/// \brief random variable *= with scalar
///
inline GaussianMixture& operator*=(GaussianMixture &x, const double r) {
  for (int i=0; i<x.models.size(); i++) {
    x.models[i].m *= r;
    x.models[i].v *= r;
  }
  return x;
}

}  // namespace dist
}  // namespace edda


#endif  // GAUSSIAN_MIXTURE_H
