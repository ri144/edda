#ifndef UNCERTAIN_ISOCONTOUR_H_
#define UNCERTAIN_ISOCONTOUR_H_

#include "dataset.h"
#include "distributions/distribution.h"

namespace edda{

// Currently the dataset is assumed in regular cartesian grids
template <typename Dist, ENABLE_IF_BASE_OF(Dist, dist::Distribution)>
std::shared_ptr<Dataset<float> > uncertainIsocontour(std::shared_ptr<Dataset<Dist> > dataset, float isov)
{
  ReturnStatus r;
  int *dim = dataset->getDimension();

  int i,j,k;
  int count=0;

  shared_ary<float> cdfField (dim[0]*dim[1]*dim[2]); //precompute cdf for speedup
  for (k=0; k<dim[2]; k++)
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++) {
        // compute level crossing
        cdfField[count] = dist::getCdf(dataset->at_comp(i,j,k ,r), isov);
        count++;
      }

  int new_dim[3]={dim[0]-1, dim[1]-1, dim[2]-1};
  shared_ary<float> probField (new_dim[0]*new_dim[1]*new_dim[2]);
  count = 0;
  for (k=0; k<new_dim[2]; k++)
    for (j=0; j<new_dim[1]; j++)
      for (i=0; i<new_dim[0]; i++) {
        // compute level crossing
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
        for (int l=0; l<8; l++)
        {
          prob1 *= cdf[l];
          prob2 *= 1-cdf[l];
        }

        probField[count] = 1.-prob1-prob2;
        count++;
      }

  return std::shared_ptr<Dataset<float> > (
        new Dataset<float> (
          new RegularCartesianGrid(new_dim[0], new_dim[1], new_dim[2]),
          new DataArray<float>(probField)
          ));
}

} // namespace edda
#endif // UNCERTAIN_ISOCONTOUR_H_
