#ifndef RANDOM_SAMPLE_FIELD
#define RANDOM_SAMPLE_FIELD

#include "common.h"
#include "core/thrust_common.h"
#include "distributions/variant.h"
#include "distributions/gaussian_mixture.h"

namespace edda {

struct GetSample_functor {
  __host__ __device__
  template<class Dist>
  double operator() (const Dist &x ) const { return dist::getSample(x); }
};

template <class InputIterator, class OutputIterator>
void randomSampleField(InputIterator begin, InputIterator end, OutputIterator out)
{
  thrust::transform(begin, end, out, GetSample_functor());
}

} // edda

#endif // RANDOM_SAMPLE_FIELD
