// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DATA_ARRAY_H_
#define DATA_ARRAY_H_

#include <cassert>
#include <iostream>
#include <memory>
#include <boost/any.hpp>
#include "math/vector_matrix.h"
#include "core/shared_ary.h"
#include "distributions/distribution.h"

namespace edda {

///
/// \brief The AbstractDataArray class used in class Dataset
///
class AbstractDataArray
{
public:
    AbstractDataArray() {}

    virtual ~AbstractDataArray() {}

    ///
    /// Get the number of elements
    ///
    virtual size_t getLength() =0;

    ///
    /// Get data at the given index
    /// Note: We don't return reference because the return data can be derived from the original data.  Use setItem() to change the data content.
    ///
    virtual boost::any getItem(size_t idx) =0;

    ///
    /// Set data at the given index
    ///
    virtual void setItem(size_t idx, const boost::any &item) =0;
};

//---------------------------------------------------------------------------------------
///
/// \brief A simple implementation of DataArray.
/// This is the class that holds the actual array, in a smart pointer.
///
template<typename T>
class DataArray: public AbstractDataArray
{
protected:
  shared_ary<T> array;
public:
  DataArray(shared_ary<T> array) { this->array = array; }

  virtual ~DataArray() { }

  virtual boost::any getItem(size_t idx) { return boost::any( array[idx] );  }

  virtual void setItem(size_t idx, const boost::any &item) { array[idx] = boost::any_cast<T>( item );  }

  virtual size_t getLength() { return array.getLength(); }

  virtual shared_ary<T> getRawArray() { return array; }
};

//---------------------------------------------------------------------------------------
///
/// \brief Returns a sampling of the distribution.
/// Derived from DataArray
///
template<typename Dist, ENABLE_IF_BASE_OF(Dist, dist::Distribution)>
class SampledDataArray: public DataArray<Dist>
{
public:
  SampledDataArray(shared_ary<Dist> array) : DataArray<Dist>(array) {}

  virtual ~SampledDataArray() { }

  virtual boost::any getItem(size_t idx) { return boost::any ( dist::getSample(this->array[idx]) );    }

};

//---------------------------------------------------------------------------------------
///
/// \brief For each element of the tuple, sample a value from the distribution.
/// getItem() Return a Tuple of float types (Can be VECTOR3 or VECTOR4)
/// Derived from DataArray
///
template<typename TupleType, typename OutputType = Vector<float, TupleType::LENGTH> >
class SampledIndepTupleArray: public DataArray<TupleType>
{
public:
  SampledIndepTupleArray(shared_ary<TupleType> array) : DataArray<TupleType>(array) {}

  virtual boost::any getItem(size_t idx) {
    TupleType &data_dist = this->array[idx];
    OutputType data_sampled;

    // sample from each dimension
    for (int i=0; i < TupleType::LENGTH; i++)
      data_sampled[i] = dist::getSample(data_dist[i]);

    return boost::any ( data_sampled );
  }
};


} // namespace edda

#endif // DATA_ARRAY_H_
