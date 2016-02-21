#ifndef GMM_ND_ARRAY_H
#define GMM_ND_ARRAY_H

#include <memory>
#include <thrust/iterator/constant_iterator.h>
#include <vtkDataArray.h>

#include <core/ndarray.h>
#include <core/thrust_common.h>
#include <distributions/gaussian_mixture.h>

namespace edda{

#if 0
namespace detail {
  struct MakeStridedGmm : public thrust::unary_function<int,dist::Variant > {
    // [GMMs][n] array
    const std::vector<NdArray<Real> > &dataArray;
    MakeStridedGmm(const std::vector<NdArray<Real> > &dataArray_) : dataArray(dataArray_) {}

    __host__ __device__
    dist::GaussianMixture<MAX_GMMs> operator() (int idx) {
      dist::GaussianMixture<MAX_GMMs> gmm;
      int narray = (int)dataArray.size();
      for (int i=0; i<narray; i++) {
        gmm.models[i/3].p[i%3] = dataArray[i].data()[idx];
      }
      if (narray==2)
        gmm.models[0].w = 1.;
      return gmm;
    }

#if 0
    typedef thrust::tuple<Real, Real, Real> Tuple3;

    __host__ __device__
    inline dist::GMMTuple makeModel(Tuple3 t) {
      dist::GMMTuple gt;
      gt.p[0] = thrust::get<0>(t);
      gt.p[1] = thrust::get<1>(t);
      gt.p[2] = thrust::get<2>(t);
      return gt;
    }

    dist::Variant operator() (const thrust::tuple<Tuple3> &tuple) {
      dist::GaussianMixture<1> gmm;
      gmm.models[0] = makeModel(thrust::get<0> (tuple));
      return dist::Variant(gmm);
    }
#endif
  };

} // detail
#endif

class GmmNdArray{
public:
  // [GMMs][n] array
  std::vector<NdArray<Real> > dataArray;

  GmmNdArray(const std::vector<NdArray<Real> > &data_) : dataArray(data_) {}

  thrust::transform_iterator<GmmNdArray, thrust::counting_iterator<int> > begin() {
    return thrust::make_transform_iterator( thrust::make_counting_iterator(0), *this );
  }

  thrust::transform_iterator<GmmNdArray, thrust::counting_iterator<int> > end() {
    return thrust::make_transform_iterator( thrust::make_counting_iterator((int)dataArray[0].num_of_elems()), *this );
  }

  __host__ __device__
  dist::GaussianMixture<MAX_GMMs> operator() (int idx) {
    return get_val(idx);
  }

  __host__ __device__
  dist::GaussianMixture<MAX_GMMs> get_val(int idx) {

    dist::GaussianMixture<MAX_GMMs> gmm;
    int narray = (int)dataArray.size();
    for (int i=0; i<narray; i++) {
      gmm.models[i/3].p[i%3] = dataArray[i].data()[idx];
    }
    if (narray==2)
      gmm.models[0].w = 1.;
    return gmm;
  }
};


} // edda

#endif // #ifndef GMM_ND_ARRAY_H
