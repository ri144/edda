// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef NDARRAY_H_
#define NDARRAY_H_

#include <cassert>
#include <iostream>

#include "thrust_common.h"

namespace edda {

const int kNdArrayMaxDims = 8;

template <typename Type>
class NdArray {

 public:

  NdArray(Type* data, const std::initializer_list<int>& dims) {
    assert(dims.size() <= kNdArrayMaxDims);

    num_of_dims_ = dims.size();

    auto it = dims.begin();
    for (int i = 0; i < num_of_dims_; ++i)
      dims_[i] = *(it + i);

    UpdateStrides();

    int num_of_elems = 1;
    for (int i = 0; i < num_of_dims_; ++i)
      num_of_elems *= dims_[i];

    data_.swap( thrust::device_vector<Type> (num_of_elems) );
    thrust::copy(data, data+num_of_elems, data_.begin() );
  }

  NdArray(Type* data, int num_of_dims, int* dims) {
    assert(num_of_dims <= kNdArrayMaxDims);

    num_of_dims_ = num_of_dims;

    for (int i = 0; i < num_of_dims_; ++i)
      dims_[i] = dims[i];

    UpdateStrides();

    int num_of_elems = 1;
    for (int i = 0; i < num_of_dims_; ++i)
      num_of_elems *= dims_[i];

    data_.resize( num_of_elems) ;
    thrust::copy(data, data+num_of_elems, data_.begin() );

  }

  NdArray(int num_of_dims, int* dims) {
    assert(num_of_dims <= kNdArrayMaxDims);

    num_of_dims_ = num_of_dims;

    for (int i = 0; i < num_of_dims_; ++i)
      dims_[i] = dims[i];

    UpdateStrides();

    int num_of_elems = 1;
    for (int i = 0; i < num_of_dims_; ++i)
      num_of_elems *= dims_[i];

    data_.resize( num_of_elems) ;
  }

#if 0
  NdArray(Type* data, int num_of_dims, int* dims, int* strides) {
    assert(num_of_dims <= kNdArrayMaxDims);

    num_of_dims_ = num_of_dims;

    for (int i = 0; i < num_of_dims_; ++i) {
      dims_[i] = dims[i];
      strides_[i] = strides[i];
    }

    if (memory_type_ == kHost) {
      data_ = data;
      ownership_ = false;
    }

#ifdef WITH_CUDA
    if (memory_type_ == kDevice) {
      int num_of_elems = 1;
      for (int i = 0; i < num_of_dims_; ++i)
        num_of_elems *= dims_[i];

      CudaAllocate(data_, num_of_elems);
      CudaMemCopyHostToDevice(data_, data, num_of_elems);

      ownership_ = true;
    }
#endif  // WITH_CUDA
  }
#endif

  ~NdArray() {
  }

  __host__ __device__ int num_of_dims() const {
    return num_of_dims_;
  }

  __host__ __device__ const int* dims() const {
    return dims_;
  }

  __host__ __device__ thrust::device_vector<Type> & data() {
    return data_;
  }

  __host__ __device__ const thrust::device_vector<Type> & data() const {
    return data_;
  }

  __host__ __device__ Type* get_ptr(const std::initializer_list<int>& ind) {
    int nd = num_of_dims_;
    int* strides = strides_;
    Type* dptr = data_;
    auto it = ind.begin();
    while (nd--) {
      dptr += (*strides++) * (*it++);
    }
    return dptr;
  }

  __host__ __device__ Type get_val(const std::initializer_list<int>& ind) {
    return *(get_ptr(ind));
  }

  __host__ __device__ void set_val(
      const std::initializer_list<int>& ind, const Type& val) {
    *(get_ptr(ind)) = val;
  }

  __host__ __device__ void Reshape(
      const std::initializer_list<int>& newshape) {

    // total size check, todo: slice check
    int newsize=1, oldsize=1;
    auto it = newshape.begin();
    for (; it != newshape.end(); ++it)
      newsize *= *it;
    int i;
    for (i = 0; i < dims_[i]; ++i)
      oldsize *= dims_[i];
    assert (newsize == oldsize);

    num_of_dims_ = newshape.size();

    it = newshape.begin();
    for (; it != newshape.end(); ++it)
      dims_[i] = *it;

    UpdateStrides();
  }

  // NdArray* Slice(const std::initializer_list<int>& from,
  //                const std::initializer_list<int>& to) {
  //   int* dims = new int[num_of_dims_];
  //   auto itf = from.begin(), itt = to.begin();
  //   for (int i = 0; i < num_of_dims_; ++i) {
  //     dims[i] = *(itt + i) - *(itf + i) + 1;
  //   }

  //   int* strides = new int[num_of_dims_];
  //   for (int i = 0; i < num_of_dims_; ++i) {
  //     strides[i] = strides_[i];
  //   }

  //   return new NdArray(get_ptr(from), num_of_dims_, dims, strides);
  // }

 private:
  __host__ __device__ void UpdateStrides() {
    int stride = 1;
    for (int i = num_of_dims_ - 1; i >= 0; --i) {
      strides_[i] = stride;
      stride *= dims_[i] ? dims_[i] : 1;
    }
  }

  thrust::device_vector<Type> data_;
  int num_of_dims_ = 0;
  int dims_[kNdArrayMaxDims];
  int strides_[kNdArrayMaxDims];
  //bool ownership_ = false;
};

}  // namespace dv

#endif  // NDARRAY_H_
