#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <vector>
#include <string>
#include <cstring>
#include "dataset.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "core/shared_ary.h"

namespace edda{

typedef dist::Gaussian<float> Gaussianf;
typedef Vector3<Gaussianf> Gaussianf3;

template<typename T>
shared_ary<T> loadRawFile(const std::string &fname, size_t len) {
  FILE * fIn;
  T *pData;

  fIn = fopen(fname.c_str(), "rb");
  assert(fIn != NULL);
  pData = new T[len];
  size_t read_len = fread(pData, sizeof(T), len, fIn);
  fclose(fIn);
  assert(len == read_len);
  return shared_ary<T>(pData, len);
}

shared_ary<Gaussianf> loadGaussianRawArray(std::string meanfile, std::string stdfile, size_t len) ;
shared_ary<Gaussianf3> loadVec3GaussianRawArray(std::string meanfile, std::string stdfile, size_t len);

//----------------------------------------------------------------------------------------------------------
// dataset creators

/// create a regular grid dataset of gaussian distribution
std::shared_ptr<Dataset<Gaussianf> > loadGaussianRegularGrids(std::string &meanfile, std::string &stdfile, int dim[3]);

/// create a regular grid dataset of random values from gaussian distributions
std::shared_ptr<Dataset<double> > loadGaussianSamplingRegularGrids(std::string &meanfile, std::string &stdfile, int dim[3]);

/// create a regular grid dataset of 3d gaussian distribution
std::shared_ptr<Dataset<Gaussianf3> > loadVec3GaussianRegularGrids(std::string &meanfile, std::string &stdfile, int dim[3]);

/// create a regular grid dataset of random samples from 3d gaussian distribution
std::shared_ptr<Dataset<VECTOR3> > loadVec3GaussianSamplingRegularGrids(std::string &meanfile, std::string &stdfile, int dim[3]);

} // edda


#endif // FILE_READER_H_
