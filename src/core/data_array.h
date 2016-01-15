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
/// \brief The AbstractDataArray class used in the Dataset class
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

    ///
    /// Get the array
    ///
    virtual boost::any getRawArray() = 0;
};

//---------------------------------------------------------------------------------------
struct GetItemAsIs ;
struct GetItemSampled ;
struct GetItemSampledVector ;
///
/// \brief A simple implementation of DataArray.
///
/// This is the class that holds the actual array, in a smart pointer.
/// \param T The element type of an array.
/// \param GetItemPolicy Allows to output in different representations of the distribution.  Possible choices: GetItemAsIs, GetItemSampled, GetItemSampledVector.
///
template<typename T, class GetItemPolicy = GetItemAsIs >
class DataArray: public AbstractDataArray
{
protected:
  shared_ary<T> array;
  GetItemPolicy GetItem;
public:
  DataArray(shared_ary<T> array) { this->array = array; }

  virtual ~DataArray() { }

  virtual boost::any getItem(size_t idx) { return boost::any( GetItem(array[idx]) );  }

  virtual void setItem(size_t idx, const boost::any &item) { array[idx] = boost::any_cast<T>( item );  }

  virtual size_t getLength() { return array.getLength(); }

  virtual boost::any getRawArray() { return boost::any(array); }
};


//---------------------------------------------------------------------------------------
// GetItemPolicy implementations for DataArray
///
/// \brief Describes GetItemPolicy for DataArray.  Returns as is the input.
///
struct GetItemAsIs {
  template<class T>
  inline T &operator() (T &x) { return x; }
};

///
/// \brief Describes GetItemPolicy for DataArray.  It assumes the input is a distribution and outputs a sample of it.
///
struct GetItemSampled {
  template<class T>
  inline double operator() (T &x) { return dist::getSample(x); }
};

///
/// \brief Describes GetItemPolicy for DataArray.  It assumes the input is a tuple of distributions and outputs a vector of sampled values in floats.
///
struct GetItemSampledVector {
  template<typename VectorType, typename OutputType = Vector<float, VectorType::LENGTH> >
  OutputType operator() (const VectorType &tuple)
  {
    OutputType tuple_sampled;

    // sample from each element
    for (int i=0; i < VectorType::LENGTH; i++)
      tuple_sampled[i] = dist::getSample(tuple[i]);

    return tuple_sampled;
  }
};

} // namespace edda

#endif // DATA_ARRAY_H_
