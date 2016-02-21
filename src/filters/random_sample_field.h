#ifndef RANDOM_SAMPLE_FIELD
#define RANDOM_SAMPLE_FIELD

#include "common.h"
#include "core/thrust_common.h"
#include "distributions/variant.h"
#include "distributions/gaussian_mixture.h"
#include "core/thrust_common_functors.h"

namespace edda {


template<class Dist>
inline ReturnStatus randomSampleField(const thrust::device_vector<Dist> &field, thrust::device_vector<Real> &out )
{
  out.clear();
  out.resize(field.size());

  thrust::transform(field.begin(), field.end(), out.begin(), GetSample_functor() );

  return ReturnStatus::SUCCESS;
}


} // edda

#endif // RANDOM_SAMPLE_FIELD
