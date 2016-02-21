#ifndef THRUST_COMMON_FUNCTIONS_H
#define THRUST_COMMON_FUNCTIONS_H

#include "thrust_common.h"

namespace edda{


struct GetSample_functor {
  __host__ __device__
  template<class Dist>
  double operator() (const Dist &x ) const { return dist::getSample(x); }
};


} // edda

#endif // THRUST_COMMON_FUNCTIONS_H
