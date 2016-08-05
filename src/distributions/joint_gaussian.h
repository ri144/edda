// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_JOINT_GAUSSIAN_H_
#define DIST_JOINT_GAUSSIAN_H_

#include <cstdlib>
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES  // For Visual Studio
#include <math.h>

#include <boost/limits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include "common.h"
#include "distribution_tag.h"
#include "core/statistics.h"
#include "invert_matrix.h"


namespace edda {
namespace dist {

// ------------------------------------------------------------------------------
///
/// \brief Defines a Gaussian class
///
struct EDDA_EXPORT JointGaussian: public ContinuousDistributionTag, public JointDistributionTag {
  ublas_vector mean; // boost's vector

  // constructor
  __host__ __device__
  JointGaussian() {
    mean = ublas::zero_vector<Real>(3);
    setCovariance( ublas::identity_matrix<Real>(3) );
  }
  __host__ __device__
  JointGaussian(const ublas_vector &mean, const ublas_matrix &cov) {
    assert(mean.size() == cov.size1() && mean.size() == cov.size2());
    this->mean = mean;
    setCovariance(cov);
  }

  __host__ __device__
  void setCovariance(const ublas_matrix &cov) {
    this->cov = cov;
    bool r = invert_matrix(this->cov, inv_cov);
    if (!r)
      std::cout << "Error: inverse covariance cannot be computed for matrix: " << cov << std::endl;
    // for debugging
    std::cout << "(JointGaussian debug) Inverse func: " << inv_cov << std::endl;
    std::cout << "(JointGaussian debug) " << ublas::prod(cov, inv_cov) << std::endl;
    this->det = determinant(this->cov);
    std::cout << "(JointGaussian debug) Det: " << det << std::endl;
  }

  __host__ __device__
  const ublas_matrix &getCovariance() const {  return this->cov ; }

  __host__ __device__
  const ublas_matrix &getInvCovariance() const {  return this->inv_cov ; }

  __host__ __device__
  double getDet() const {  return this->det ; }

private:
  ublas_matrix cov;
  ublas_matrix inv_cov;
  double det;
};

// ------------------------------------------------------------------------------
// Below defines JointGaussian related generic functions
__host__ __device__
inline std::vector<Real> getJointMean(const JointGaussian &dist)
{
  std::vector<Real> m(dist.mean.size());
  std::copy(dist.mean.begin(), dist.mean.end(), m.begin());
  return m;
}

///
/// \brief Return PDF of x
///
__host__ __device__
inline double getJointPdf(const JointGaussian &dist, const std::vector<Real> x_)
{
  int k = dist.mean.size();
  assert(x_.size() == k);
  ublas_vector x;  std::copy(x_.begin(), x_.end(), x.begin());
  x = x-dist.mean;
  ublas_vector tmp;
  tmp = ublas::prod(ublas::trans(x), dist.getInvCovariance());
  double exp = -0.5 * ublas::inner_prod(tmp, x);
  return exp / sqrt( powf(2*M_PI, k) * dist.getDet() );
}

///
/// \brief Return a random sample
///
__host__
inline std::vector<Real> getJointSample(const JointGaussian &dist)
{
  // TODO: Cholesky Decomposition
  return getJointMean(dist);
}

///
/// \brief Return a random sample using random engine
///
__host__ __device__
inline std::vector<Real> getJointSample(const JointGaussian &dist, thrust::default_random_engine &rng)
{
  // TODO: Cholesky Decomposition
  return getJointMean(dist);
}

///
/// \brief Print itself
///
__host__
inline std::ostream& operator<<(std::ostream& os, const JointGaussian &dist)
{
    os <<  "<JointGaussian: mean=" << dist.mean << ", covariance=" << dist.getCovariance() << ">" ;
    return os;
}

__host__ __device__
inline std::string getName(const JointGaussian &dist) {
    return "JointGaussian";
}


}  // namespace dist
}  // namespace edda

#endif  // DIST_GAUSSIAN_H_
