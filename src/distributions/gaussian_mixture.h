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
template <int GMMs>
class EDDA_EXPORT GaussianMixture: public ContinuousDistribution {

  void modelReduction(const std::vector<GMMTuple> &models_) {
    if (models_.size() <=  GMMs) {
      reset();
      for (size_t i=0; i<models_.size(); i++)
        models[i] = models_[i];

    } else {

      throw NotImplementedException();
    }
  }

  ///
  /// \brief Reset all model weights to 0
  ///
  __host__ __device__
  inline void reset () {
    for (int i=0; i<GMMs; i++)
      models[i].w = 0;
  }

public:
  Tuple<GMMTuple, GMMs> models;

  // constructor
  __host__ __device__
  GaussianMixture() { reset(); }

  __host__ __device__
  template <int GMMs_>
  void assign (const Tuple<GMMTuple, GMMs_> &models) {
    std::vector<GMMTuple> vmodels_;
    for (int i=0; i<GMMs_; i++)
    {
      GMMTuple t = models[i];
      vmodels_.push_back(t);
    }
    modelReduction(vmodels_);
  }

  GaussianMixture(const std::vector<GMMTuple> &models_) {
    modelReduction(models_);
  }

  ///
  /// \brief Scale sum of weights to 1.
  ///
  __host__ __device__
  void normalizeWeights() {
    double sum = 0;
    int i;
    for (i=0; i<GMMs; i++) {
      sum += models[i].w;
    }
    if (sum == 0) return;
    for (i=0; i<GMMs; i++) {
      models[i].w /= sum;
    }
  }
};

// ------------------------------------------------------------------------------
// Below defines GaussianMixture related generic functions

///
/// \brief Return mean
///
__host__ __device__
template <int GMMs>
inline double getMean(const GaussianMixture<GMMs> &dist)
{
  double mean = 0;
  for (int i=0; i<GMMs; i++)
  {
    mean += dist.models[i].m * dist.models[i].w;
  }
  return mean;
}

///
/// \brief Return variance.
///
/// if f(x) = sum_i( w_i * f_i(x) ), v_i = variance of f_i, and m_i = mean of f_i, then
/// var(f) = sum_i( w_i * v_i + w_i * m_i^2 ) - (sum_i( w_i * m_i) )^2.
///
/// ref: http://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians
/// (code not verified)
///
__host__ __device__
template <int GMMs>
inline double getVar(const GaussianMixture<GMMs> &dist)
{
  // Let the first summation as term1 and second as term2
  double term1=0, term2=0;
  for (int i=0; i<GMMs; i++)
  {
    const GMMTuple & model = dist.models[i];
    term1 += (double)model.w * model.v + (double)model.w * model.m * model.m ;
    term2 += (double)model.w * model.m;
  }
  return term1 - term2 * term2;
}

///
/// \brief Return PDF of x
///
__host__ __device__
template <int GMMs>
inline double getPdf(const GaussianMixture<GMMs> &dist, const double x)
{
  double p=0;
  for (int i=0; i<GMMs; i++)
  {
    p += getPdf( Gaussian(dist.models[i].m, dist.models[i].v), x ) * dist.models[i].w;
  }
  return p;
}

///
/// \brief Return a sample
///
/// Note: To ensure correct sampling distribution, the sum of weights
/// should be normalized to 1 before calling this function.
///
__host__ __device__
template <int GMMs>
inline double getSample(const GaussianMixture<GMMs> &dist)
{
  float ratio = rand() / (float)RAND_MAX;
  float accumulated = 0;
  for (int i=0; i<GMMs; i++)
  {
    accumulated += dist.models[i].w;
    if (ratio < accumulated) {
      return getSample( Gaussian(dist.models[i].m, dist.models[i].v) );
    }
  }
  // return sample from the last model
  return getSample( Gaussian(dist.models[GMMs-1].m, dist.models[GMMs-1].v) );
}

///
/// \brief Return CDF of x
///
__host__ __device__
template <int GMMs>
inline double getCdf(const GaussianMixture<GMMs> &dist, double x)
{
  double cdf=0;
  for (int i=0; i<GMMs; i++)
  {
    cdf += getCdf(Gaussian(dist.models[i].m, dist.models[i].v), x) * dist.models[i].w;
  }
  return cdf;
}

///
/// \brief Print itself
///
__host__ __device__
template <int GMMs>
inline std::ostream& operator<<(std::ostream& os, const GaussianMixture<GMMs> &dist)
{
  os <<  "<GaussianMixture(m,v,w):";
  for (int i=0; i<GMMs; i++)
    os << " (" << dist.models[i].m << "," << dist.models[i].v << "," << dist.models[i].w << ")";
  os << ">";
  return os;
}

// ------------------------------------------------------------------------------
// Below defines Gaussian related arithmetics

///
/// \brief random variable with unary -
///
__host__ __device__
template <int GMMs>
inline GaussianMixture<GMMs>& operator-(GaussianMixture<GMMs> &x)
{
  for (int i=0; i<GMMs; i++)
    x.models[i].m = -x.models[i].m;
  return x;
}

///
/// \brief random variable +=
///
__host__ __device__
template <int GMMs>
inline GaussianMixture<GMMs>& operator+=(GaussianMixture<GMMs> &x, const GaussianMixture<GMMs>& rhs) {
  throw NotImplementedException();
}

///
/// \brief random variable += with scalar
///
__host__ __device__
template <int GMMs>
inline GaussianMixture<GMMs>& operator+=(GaussianMixture<GMMs> &x, const double r) {
  for (int i=0; i<GMMs; i++)
    x.models[i].m += r;
  return x;
}

///
/// \brief random variable *= with scalar
///
__host__ __device__
template <int GMMs>
inline GaussianMixture<GMMs>& operator*=(GaussianMixture<GMMs> &x, const double r) {
  for (int i=0; i<GMMs; i++) {
    x.models[i].m *= r;
    x.models[i].v *= r;
  }
  return x;
}

}  // namespace dist
}  // namespace edda


#endif  // GAUSSIAN_MIXTURE_H
