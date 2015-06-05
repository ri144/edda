#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <vector>
#include <string>
#include <cstring>
#include <math.h>

#include "core/shared_ary.h"

namespace edda{

template<typename T>
shared_ary<T> read_raw(const std::string &fname, size_t len) {
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


} // edda


#endif // FILE_READER_H_
