// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DATA_ARRAY_H_
#define DATA_ARRAY_H_

#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
#include "abstract_data_array.h"
#include "vector_matrix.h"
#include "shared_ary.h"
#include "distributions/distribution.h"

namespace edda {

//---------------------------------------------------------------------------------------
/// \brief A simple implementation of DataArray.
///
/// This is the class that holds the actual array, in a smart array.
/// \param T The element type of an array.
///
template<typename T>
class ScalarArray: public AbstractDataArray
{
protected:
  shared_ary<T> array;
public:
  ScalarArray(shared_ary<T> array) { this->array = array; }

  virtual ~ScalarArray() { }

  virtual size_t getLength() { return array.getLength(); }

  virtual int getNumComponents() { return 1; }

  virtual dist::Variant getScalar(size_t idx) { return array[idx]; }

  virtual std::vector<dist::Variant> getVector(size_t idx) { return std::vector<dist::Variant> (1, array[idx] ); }

  virtual boost::any getItem(size_t idx) { return boost::any( array[idx] );  }

  virtual void setItem(size_t idx, int component, const boost::any &item) { array[idx] = boost::any_cast<T>( item );  }

  virtual boost::any getRawArray() { return boost::any(array); }
};


//---------------------------------------------------------------------------------------

template<typename T, int N>
class VectorArray: public AbstractDataArray
{
protected:
  shared_ary<Vector<T,N> > array;
public:
  VectorArray(shared_ary<Vector<T,N> > array) { this->array = array; }

  virtual ~VectorArray() { }

  virtual size_t getLength() { return array.getLength(); }

  virtual int getNumComponents() { return 1; }

  virtual dist::Variant getScalar(size_t idx) {
      throw std::runtime_error("Requesting scalar in a VectorArray");
    }

  virtual std::vector<dist::Variant> getVector(size_t idx) {
    std::vector<dist::Variant> v(N);
    for (int i=0; i<N; i++)
      v[i] = array[idx][i];
    return v; }

  virtual boost::any getItem(size_t idx) { return boost::any( array[idx] );  }

  virtual void setItem(size_t idx, int component, const boost::any &item) { array[idx] = boost::any_cast<Vector<T,N> >( item );  }

  virtual boost::any getRawArray() { return boost::any(array); }
};



} // namespace edda

#endif // DATA_ARRAY_H_
