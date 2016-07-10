// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

/// Histogram: TODO
///

#include <iostream>
#include <typeinfo>
#include <iostream>

#include "common.h"
#include "edda_export.h"
#include "core/vector_matrix.h"
#include "core/tuple.h"
#include "core/thrust_common.h"
#include "distribution_tag.h"

namespace edda {
namespace dist {

///
/// \brief The Distribution class is a root class for all distribution-type classes.
///
/// This is useful for applications to identify whether a class is a distribution-type class, by using ENABLE_IF_BASE_OF()
///
class EDDA_EXPORT Histogram: public DiscreteDistributionTag {

};

///
/// \brief Compute the mean of the distribution
///
inline double getMean(const Histogram &dist)
{
    std::cerr << "Generic computation not implemented" << std::endl;
    throw NotImplementedException();
    return 0;
}

///
/// \brief Compute Variance
///
inline double getVar(const Histogram &dist)
{
  std::cerr << "Generic computation not implemented" << std::endl;
  throw NotImplementedException();
  return 0;
}

///
/// \brief Return PDF of x
///
inline double getPdf(const Histogram &dist, const double x)
{
  std::cerr << "Generic computation not implemented" << std::endl;
  throw NotImplementedException();
  return 0;
}

///
/// \brief Get a Monte-Carlo sample of the distribution. We rely on the specific implementation from each distribution.
///
inline double getSample(const Histogram &dist)
{
  std::cerr << "Generic computation not implemented" << std::endl;
  throw NotImplementedException();
  return 0;
}


///
/// \brief Return CDF of x
///
__host__ __device__
inline double getCdf(const Histogram &dist, double x)
{
  std::cerr << "Generic computation not implemented" << std::endl;
  throw NotImplementedException();
  return 0;
}

///
/// \brief Print itself
///
inline std::ostream& operator<<(std::ostream& os, const Histogram &dist)
{
  return os;
}

__host__ __device__
inline const char *getName(const Histogram &x) {
    return "Histogram";
}

}  // namespace dist
}  // namespace edda

#endif // HISTOGRAM_H
