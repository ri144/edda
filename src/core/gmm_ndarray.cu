#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include "gmm_ndarray.h"

namespace edda{


namespace detail{

  __host__ __device__
  dist::GaussianMixture<MAX_GMMs> MakeStridedGmm::operator() (int idx) const
  {
    dist::GaussianMixture<MAX_GMMs> gmm;
    int narray = narrays;
    for (int i=0; i<narray; i++) {
      gmm.models[i/3].p[i%3] = dDataPtrArray[i]->get_val(idx);
    }
    if (narray==2)
      gmm.models[0].w = 1.;
    return gmm;
  }
} // detail


} // edda
