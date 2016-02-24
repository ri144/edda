#ifndef RANDOM_SAMPLE_FIELD
#define RANDOM_SAMPLE_FIELD

#include <ctime>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include "common.h"
#include "core/thrust_common.h"
#include "distributions/variant.h"
#include "distributions/gaussian_mixture.h"

namespace edda {



struct GenRand: public thrust::unary_function<int, thrust::default_random_engine>
{
  __device__
  thrust::default_random_engine operator () (int idx)
  {
      thrust::default_random_engine randEng;
      randEng.discard(idx*2);
      return randEng;
  }
};

inline thrust::transform_iterator<GenRand, thrust::counting_iterator<int> > randomEngineIterator(int seed) {
  return thrust::make_transform_iterator(thrust::make_counting_iterator(seed), GenRand()) ;
}



struct GetSample_functor {
  template <class Dist> //, ENABLE_IF_BASE_OF(Dist, dist::Distribution)>
  __host__ __device__
  double operator() (thrust::tuple<Dist, thrust::default_random_engine> &tuple)
  {
    Dist &dist = thrust::get<0>(tuple);
    return dist::getSample(dist, thrust::get<1>(tuple));
  }
};

///
/// Output iterator must have elements in type 'Real'
///
template <class InputIterator, class OutputIterator>
void randomSampleField(InputIterator begin, InputIterator end, OutputIterator out)
{
  static int seed ;
  seed += time(NULL) % 10000;
  int n = end-begin;
  thrust::transform( thrust::make_zip_iterator(thrust::make_tuple(begin, randomEngineIterator(seed)) ) ,
                     thrust::make_zip_iterator(thrust::make_tuple(end, randomEngineIterator(seed+n)) ),
                     out, GetSample_functor()) ;
  seed += n;
}

} // edda

#endif // RANDOM_SAMPLE_FIELD
