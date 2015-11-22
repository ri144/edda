#ifndef FILE_WRITER_H_
#define FILE_WRITER_H_

#include <vector>
#include <string>
#include <cstring>
#include "dataset.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "core/shared_ary.h"
#include "io/nrrd.h"
#include "path.h"

namespace edda{

// Input: dataset, filename prefix
// Output: raw file (.raw) and nrrd file (.nrrd)
// Assume the input dataset is regular grids and the array uses raw array
template<typename T>
void writeRawNrrdFile(std::shared_ptr<Dataset<T> > dataset, std::string filename)
{
  std::string rawfilename = filename+".raw";
  // output raw file
  FILE *fp = fopen(rawfilename.c_str(), "wb");
  DataArray<T> *array = dynamic_cast<DataArray<T> *>( dataset->getArray() );
  assert(array!=NULL);
  fwrite( array->getRawArray().get() , dataset->getArray()->getLength(), sizeof(float), fp);
  fclose(fp);
  // output nrrd file
  int *new_dim = dynamic_cast<RegularCartesianGrid*>( dataset->getGrid() )->getDimension();
  write_nrrd_3d((filename+".nrrd").c_str(), rawfilename.c_str(), new_dim[0], new_dim[1], new_dim[2], "float");
}

} // namespace edda

#endif // FILE_WRITER_H_
