#include <cassert>
#include "file_reader.h"
#include "dataset.h"
#include "distributions/gaussian.h"
#include "math/vector_matrix.h"

using namespace std;

namespace edda{

/// array loader

shared_ary<Gaussianf> loadGaussianRawArray(string meanfile, string stdfile, size_t len)
{
  shared_ary<float> pMean = loadRawFile<float>(meanfile, len);
  shared_ary<float> pStd = loadRawFile<float>(stdfile, len);

  // Create Gaussian array
  Gaussianf *pData = new Gaussianf[len];
  size_t i;
  for (i=0; i<len; i++)
  {
    pData[i] = Gaussianf(pMean[i], pStd[i]*pStd[i]);
  }
  // return smart pointer of the array
  return shared_ary<Gaussianf> (pData, len);
}

shared_ary<Gaussianf3> loadVec3GaussianRawArray(string meanfile, string stdfile, size_t len) {
  shared_ary<Tuple3<float> > pMean3 = loadRawFile<Tuple3<float> >(meanfile, len);
  shared_ary<Tuple3<float> > pStd3 = loadRawFile<Tuple3<float> >(stdfile, len);

  // Create Gaussian array
  Gaussianf3 x;
  shared_ary<Gaussianf3> data(new Gaussianf3[len], len);
  size_t i;
  for (i=0; i<len; i++)
  {
    data[i] = Gaussianf3(
                Gaussianf(pMean3[i][0], pStd3[i][0]*pStd3[i][0]),
                Gaussianf(pMean3[i][1], pStd3[i][1]*pStd3[i][1]),
                Gaussianf(pMean3[i][2], pStd3[i][2]*pStd3[i][2]) );
  }
  return data;
}

///////////////////////////////////////
/// dataset creator

/// create a regular grid dataset of gaussian distribution
shared_ptr<Dataset<Gaussianf> > loadGaussianRegularGrids(string &meanfile, string &stdfile, int dim[3])
{
  shared_ary<Gaussianf> array = loadGaussianRawArray(meanfile, stdfile, dim[0]*dim[1]*dim[2]);
  Dataset<Gaussianf> *dataset = new Dataset<Gaussianf> (
        new RegularCartesianGrid (dim[0], dim[1], dim[2]),
        new DataArray< Gaussianf >( array )
    );
  return shared_ptr<Dataset<Gaussianf> > (dataset);
}

/// create a regular grid dataset of random values from gaussian distributions
shared_ptr<Dataset<double> > loadGaussianSamplingRegularGrids(string &meanfile, string &stdfile, int dim[3])
{
  shared_ary<Gaussianf> array = loadGaussianRawArray(meanfile, stdfile, dim[0]*dim[1]*dim[2]);
  Dataset<double> *dataset = new Dataset<double> (
        new RegularCartesianGrid (dim[0], dim[1], dim[2]),
        new DataArray< Gaussianf, GetItemSampled >( array )
    );
  return shared_ptr<Dataset<double> > (dataset);
}

/// create a regular grid dataset of 3d gaussian distribution
shared_ptr<Dataset<Gaussianf3> > loadVec3GaussianRegularGrids(string &meanfile, string &stdfile, int dim[3])
{
  shared_ary<Gaussianf3> array = loadVec3GaussianRawArray(meanfile, stdfile, dim[0]*dim[1]*dim[2]);
  Dataset<Gaussianf3> *dataset = new Dataset<Gaussianf3> (
            new RegularCartesianGrid (dim[0], dim[1], dim[2]),
            new DataArray<Gaussianf3> (array) );
  return shared_ptr<Dataset<Gaussianf3> > (dataset);
}

/// create a regular grid dataset of random samples from 3d gaussian distribution
shared_ptr<Dataset<VECTOR3> > loadVec3GaussianSamplingRegularGrids(string &meanfile, string &stdfile, int dim[3])
{
  shared_ary<Gaussianf3> array = loadVec3GaussianRawArray(meanfile, stdfile, dim[0]*dim[1]*dim[2]);
  Dataset<VECTOR3> *dataset = new Dataset<VECTOR3> (
              new RegularCartesianGrid (dim[0], dim[1], dim[2]),
              new DataArray<Gaussianf3,GetItemSampledVector> (array)
  );
  return shared_ptr<Dataset<VECTOR3> > (dataset);
}

} // edda
