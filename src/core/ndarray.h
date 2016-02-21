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

  thrust::device_ptr<Type> dev_ptr;
  int num_of_dims_ = 0;
  int num_of_elems_ = 0;
  int dims_[kNdArrayMaxDims];
  int strides_[kNdArrayMaxDims];
  bool ownership_ = true;


public:
  ///
  /// \brief Create an empty array
  ///
  NdArray(int num_of_dims, int* dims) {

    init_shape(num_of_dims, dims);

    dev_ptr = thrust::device_malloc<Type>(num_of_elems_);

  }

  ///
  /// \brief Create an empty array
  ///
  NdArray(const std::initializer_list<int>& dims) {

    init_shape(dims);

    dev_ptr = thrust::device_malloc<Type>(num_of_elems_);

  }

  ///
  /// \brief Create a device array and copy host content to it.
  ///
  NdArray(Type* data, int num_of_dims, int* dims)
    : NdArray(num_of_dims, dims)
  {
    thrust::copy(data, data+num_of_elems_, dev_ptr );
  }

  ///
  /// \brief Create a device array and copy host content to it.
  ///
  NdArray(Type* data, const std::initializer_list<int>& dims)
    : NdArray( dims )
  {
    thrust::copy(data, data+num_of_elems_, dev_ptr );
  }

  ///
  /// \brief pass a device pointer.
  ///
  /// Note: NdArray now does not own the array
  ///
  NdArray(thrust::device_ptr<Type> dev_ptr, int num_of_dims, int* dims)
    : NdArray( num_of_dims, dims )
  {
    this->dev_ptr  = dev_ptr;
    this->ownership_ = false;
  }


  ~NdArray() {
    if (this->ownership_) {
      thrust::device_free( dev_ptr );
    }
  }

  template <typename OutputIterator>
  void copy_to_host(OutputIterator out) {

    thrust::copy(begin(), end(), out);
  }

  __host__ __device__ int num_of_dims() const {
    return num_of_dims_;
  }

  __host__ __device__ int num_of_elems() const {
    return num_of_elems_;
  }

  __host__ __device__ const int* dims() const {
    return dims_;
  }

  __host__ __device__ void set_ownership(bool own) {
    this->ownership_ = own;
  }

  __host__ __device__ thrust::device_ptr<Type> & data() {
    return dev_ptr;
  }

  __host__ __device__ const thrust::device_ptr<Type> begin() const {
    return dev_ptr;
  }

  __host__ __device__ const thrust::device_ptr<Type> end() const {
    return dev_ptr + num_of_elems_;
  }

  __host__ __device__ thrust::device_ptr<Type> get_ptr(const std::initializer_list<int>& ind) {
    int nd = num_of_dims_;
    int* strides = strides_;
    thrust::device_ptr<Type> dptr = dev_ptr;
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

    // total size check
    int newsize=1;
    auto it = newshape.begin();
    for (; it != newshape.end(); ++it)
      newsize *= *it;
    assert (newsize == num_of_elems_);

    num_of_dims_ = newshape.size();

    it = newshape.begin();
    for (int i = 0; i < num_of_dims_; ++i, ++it)
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
  __host__ __device__ void init_shape(const std::initializer_list<int>& dims) {

    assert(dims.size() <= kNdArrayMaxDims);

    num_of_dims_ = dims.size();

    auto it = dims.begin();
    for (int i = 0; i < num_of_dims_; ++i)
      dims_[i] = *(it + i);

    UpdateStrides();

    num_of_elems_ = 1;
    for (int i = 0; i < num_of_dims_; ++i)
      num_of_elems_ *= dims_[i];
  }

  __host__ __device__ void init_shape(int num_of_dims, int* dims) {

    assert(num_of_dims <= kNdArrayMaxDims);

    num_of_dims_ = num_of_dims;

    for (int i = 0; i < num_of_dims_; ++i)
      dims_[i] = dims[i];

    UpdateStrides();

    num_of_elems_ = 1;
    for (int i = 0; i < num_of_dims_; ++i)
      num_of_elems_ *= dims_[i];

  }

  __host__ __device__ void UpdateStrides() {
    int stride = 1;
    for (int i = num_of_dims_ - 1; i >= 0; --i) {
      strides_[i] = stride;
      stride *= dims_[i] ? dims_[i] : 1;
    }
  }

};

}  // namespace dv

#endif  // NDARRAY_H_
