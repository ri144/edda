#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include "gmm_ndarray.h"

namespace edda{

namespace detail{


  __host__ __device__
  dist::GaussianMixture<MAX_GMMs> MakeStridedGmm::operator() (int idx)
  {
    dist::GaussianMixture<MAX_GMMs> gmm;
    int narray = (int)dataArray.size();
    for (int i=0; i<narray; i++) {
      gmm.models[i/3].p[i%3] = dataArray[i].data()[idx];
    }
    if (narray==2)
      gmm.models[0].w = 1.;
    return gmm;
  }
} // detail

thrust::transform_iterator<detail::MakeStridedGmm, thrust::counting_iterator<int> >
GmmNdArray::begin() {
  return thrust::make_transform_iterator( thrust::make_counting_iterator(0), detail::MakeStridedGmm(this->dataArray) );
}

thrust::transform_iterator<detail::MakeStridedGmm, thrust::counting_iterator<int> >
GmmNdArray::end() {
  return thrust::make_transform_iterator( thrust::make_counting_iterator((int)dataArray[0].num_of_elems()), detail::MakeStridedGmm(this->dataArray) );
}

__host__ __device__
dist::GaussianMixture<MAX_GMMs> GmmNdArray::get_val(int idx) const {
  //detail::MakeStridedGmm f(this->dataArray);
  //return f(idx);

  dist::GaussianMixture<MAX_GMMs> gmm;
  int narray = (int)dataArray.size();
  for (int i=0; i<narray; i++) {
    gmm.models[i/3].p[i%3] = dataArray[i].get_val(idx);
  }
  if (narray==2)
    gmm.models[0].w = 1.;
  return gmm;
}

} // edda
