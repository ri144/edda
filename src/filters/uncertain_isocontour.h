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
ReturnStatus levelCrossing(std::shared_ptr<GmmNdArray> gmmArray, int dim[3], double isov, std::shared_ptr<NdArray<float> > &probField);

/// \breif Compute level crossing probablility (Not using Thrust).  Compute cell-wise probability of level crossing given distributions on the grid points.
///
/// \param input  input array, in size dim[0]*dim[1]*dim[2]
/// \param isov   isovalue to query
/// \param[out] probField   output array.
/// You can pass an empty array of probField.  The output should be allocated with size (dim[0]-1)*(dim[1]-1)*(dim[2]-1).
inline ReturnStatus levelCrossingSerial(AbstractDataArray *input, int dim[3], double isov, float *probField)
{
  int i,j,k;
  int count=0;
  int new_dim[3]={dim[0]-1, dim[1]-1, dim[2]-1};

  if (input->getLength() < dim[0]*dim[1]*dim[2]) {
    return ReturnStatus::FAIL;
  }

  //precompute cdf for speedup
  shared_ary<float> cdfField (dim[0]*dim[1]*dim[2]);
  count = 0;
  for (k=0; k<dim[2]; k++)
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++) {
        // compute level crossing
        cdfField[count] = dist::getCdf( input->getScalar(count) , isov);
        count++;
      }

  count = 0;
  for (k=0; k<new_dim[2]; k++)
    for (j=0; j<new_dim[1]; j++)
      for (i=0; i<new_dim[0]; i++) {
        // compute level crossing
        // = 1 - prob. of all larger than isovalue - prob. of all less than isovalue
        double cdf[8];
#define IJK_TO_IDX(i,j,k)  (i+dim[0]*(j+dim[1]*(k)))
        cdf[0] = cdfField[IJK_TO_IDX(i  ,j  ,k)];
        cdf[1] = cdfField[IJK_TO_IDX(i+1,j  ,k  )];
        cdf[2] = cdfField[IJK_TO_IDX(i  ,j+1,k  )];
        cdf[3] = cdfField[IJK_TO_IDX(i+1,j+1,k  )];
        cdf[4] = cdfField[IJK_TO_IDX(i  ,j  ,k+1)];
        cdf[5] = cdfField[IJK_TO_IDX(i+1,j  ,k+1)];
        cdf[6] = cdfField[IJK_TO_IDX(i  ,j+1,k+1)];
        cdf[7] = cdfField[IJK_TO_IDX(i+1,j+1,k+1)];
#undef IJK_TO_IDX
        double prob1=1., prob2=1.;
        for (int l=0; l<8; l++) {
          prob1 *= cdf[l];
          prob2 *= 1.-cdf[l];
        }

        probField[count] = 1.-prob1-prob2;
        count++;
      }

  return ReturnStatus::SUCCESS;
}

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
