#ifndef GMM_ND_ARRAY_H
#define GMM_ND_ARRAY_H

#include <memory>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/host_vector.h>

#include <core/ndarray.h>
#include <core/thrust_common.h>
#include <distributions/gaussian_mixture.h>

namespace edda{

namespace detail {

  struct MakeStridedGmm: public thrust::unary_function<int, dist::GaussianMixture<MAX_GMMs> > {
    // [GMMs][n] array
    thrust::host_vector<NdArray<Real> > &dataArray;

    __host__ __device__
    MakeStridedGmm(thrust::host_vector<NdArray<Real> > &dataArray_) : dataArray(dataArray_) {}

    __host__ __device__
    dist::GaussianMixture<MAX_GMMs> operator() (int idx);
  };

} // detail

class GmmNdArray{
public:
  // [GMMs][n] array
  thrust::host_vector<NdArray<Real> > dataArray;

  GmmNdArray(const thrust::host_vector<NdArray<Real> > &data_) : dataArray(data_) {  }

  __host__ __device__
  dist::GaussianMixture<MAX_GMMs> get_val(int idx) const;

  thrust::transform_iterator<detail::MakeStridedGmm, thrust::counting_iterator<int> > begin() ;

  thrust::transform_iterator<detail::MakeStridedGmm, thrust::counting_iterator<int> > end() ;

};


} // edda

#endif // #ifndef GMM_ND_ARRAY_H
