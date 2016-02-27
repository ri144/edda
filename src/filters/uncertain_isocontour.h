#ifndef UNCERTAIN_ISOCONTOUR_H_
#define UNCERTAIN_ISOCONTOUR_H_

#include <memory>
#include <common.h>
#include <dataset/dataset.h>
#include <distributions/distribution.h>
#include <core/shared_ary.h>
#include <dataset/abstract_data_array.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "core/ndarray.h"
#include "core/gmm_ndarray.h"
#include "distributions/gaussian_mixture.h"


namespace edda{

/// \breif Compute level crossing probablility.
///
/// \param gmmArray  input array, in size dim[0]*dim[1]*dim[2]
/// \param isov   isovalue to query
/// \param[out] probField   output array.
///
///  Compute cell-wise probability of level crossing given distributions on the grid points.
/// The output will be allocated with size (dim[0]-1)*(dim[1]-1)*(dim[2]-1).
ReturnStatus levelCrossing(std::shared_ptr<GmmNdArray> gmmArray_, int dim[3], double isov, std::shared_ptr<NdArray<float> > &probField);

/// \breif Compute level crossing probablility (Not using Thrust).  Compute cell-wise probability of level crossing given distributions on the grid points.
///
/// \param input  input array, in size dim[0]*dim[1]*dim[2]
/// \param isov   isovalue to query
/// \param[out] probField   output array.
/// You can pass an empty array of probField.  The output should be allocated with size (dim[0]-1)*(dim[1]-1)*(dim[2]-1).
ReturnStatus levelCrossingSerial(AbstractDataArray *input, int dim[3], double isov, float *probField);

// serial
///
/// \brief Compute level crossing probabilities given a Dataset class.
///
/// Currently the input dataset is assumed to be in regular cartesian grids.
/// For a given dataset in <w,h,d> dimensions, the output is in <w-1, h-1, d-1> dimensions,
/// to describe the level crossing probabilities for each cell.
///
template <typename Dist, ENABLE_IF_BASE_OF(Dist, dist::Distribution)>
std::shared_ptr<Dataset<Real> > uncertainIsocontour(std::shared_ptr<Dataset<Dist> > dataset, double isov)
{
  int *dim = dataset->getDimension();
  int new_dim[3]={dim[0]-1, dim[1]-1, dim[2]-1};
  shared_ary<float> output(new_dim[0]*new_dim[1]*new_dim[2]);

  levelCrossingSerial(dataset->getArray(), dataset->getDimension(), isov, output.get());

  return make_Dataset<float>(new RegularCartesianGrid(new_dim[0], new_dim[1], new_dim[2]),
                             new ScalarArray<Real>(output) // TODO
      );
}

} // namespace edda
#endif // UNCERTAIN_ISOCONTOUR_H_
