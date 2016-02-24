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
namespace Math
{

  ///
  /// \brief Return random value between 0..1
  ///
  /// Need to pass a unique index value
  ///
  struct Rand{
    __host__ __device__
    float operator() (int index)
    {
      thrust::default_random_engine rng(index<<1);
      return rng();
    }
  };

#if 0
  namespace detail{
    struct BoxMullerGPU{
      __device__
      void operator() (float& u1, float& u2){
        float   r = sqrtf(-2.0f * logf(u1));
        float phi = 2 * M_PI * u2;
        u1 = r * __cosf(phi);
        u2 = r * __sinf(phi);
      }
    };
  }

  struct BoxMuller{
    __device__
    void operator() (int index){
      float p = Rand(index);
      float   r = sqrtf(-2.0f * logf(u1));
      float phi = 2 * M_PI * u2;
      u1 = r * __cosf(phi);
      u2 = r * __sinf(phi);
    }
  };
#endif
} // Math

} // edda

#endif // THRUST_COMMON_H
