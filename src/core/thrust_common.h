#ifndef THRUST_COMMON_H
#define THRUST_COMMON_H

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random.h>

#include <cmath>
#include <ctime>

namespace edda{

// cross-platform math
// singleton class
class Math
{
private:
  thrust::default_random_engine rng;
  thrust::uniform_real_distribution<float> dist;

  Math() {
    rng = thrust::default_random_engine(time(NULL));
    dist = thrust::uniform_real_distribution<float> (0, RAND_MAX);
  }

  Math(Math const&) ;
  void operator=(Math const&); // Don't implement

public:
  __host__ __device__
  static Math &getInstance() {
    static Math instance;
    return instance;
  }

  __host__ __device__
  static float rand()  {
    return getInstance().dist(getInstance().rng);
  }

};

} // edda

#endif // THRUST_COMMON_H
