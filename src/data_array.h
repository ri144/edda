// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DATA_ARRAY_H_
#define DATA_ARRAY_H_

#include <cassert>
#include <iostream>
#include <memory>
#include <boost/any.hpp>
#include "vector_matrix.h"
#include "core/shared_ary.h"
#include "distributions/distribution.h"

namespace edda {


class AbstractDataArray
{
public:
    AbstractDataArray() {}

    virtual ~AbstractDataArray() {}

    // get number of elements
    virtual size_t getLength() =0;

    // We don't return reference because the return data can be derived from raw data
    virtual boost::any getItem(size_t idx) =0;

};

#if 0
class AbstractVECTOR3Array: public AbstractDataArray
{
public:
    // We don't return reference because the return data can be derived from raw data
    virtual VECTOR3 getVector3(size_t idx) =0;

    virtual boost::any getItem(size_t idx) {
        return boost::any ( getVector3(idx) );
    }
};
#endif

/////////////////////////////////////////////////////////
template<typename T>
class DataArray: public AbstractDataArray
{
protected:
  shared_ary<T> array;
public:
  DataArray(shared_ary<T> array) { this->array = array; }

  virtual ~DataArray() { }

  virtual boost::any getItem(size_t idx) { return boost::any( array[idx] );  }

  virtual size_t getLength() { return array.getLength(); }

  virtual shared_ary<T> getRawArray() { return array; }
};

template<typename Dist, ENABLE_IF_BASE_OF(Dist, dist::Distribution)>
class SampledDataArray: public DataArray<Dist>
{
public:
  SampledDataArray(shared_ary<Dist> array) : DataArray<Dist>(array) {}

  virtual ~SampledDataArray() { }

  virtual boost::any getItem(size_t idx) { return boost::any ( dist::getSample(this->array[idx]) );    }

};

// For each element of the tuple, sample a value from the distribution.  Return a Tuple of float types (Can be VECTOR3 or VECTOR4)
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
#if 0
/////////////////////////////////////////////////////////
class VECTOR3Array: public AbstractVECTOR3Array
{
protected:
    VECTOR3 *array;
public:
    VECTOR3Array(VECTOR3 *array, size_t len) { this->array = array; this->len = len; }

    ~VECTOR3Array() { delete[] array; }

    virtual VECTOR3 getVector (size_t idx) {
        assert( idx >=0 && idx < len );
        return array[idx];
    }
};

template<typename Dist>
class SampledVECTOR3Array: public AbstractVECTOR3Array
{
protected:
    Dist *array;
public:
    SampledVECTOR3Array(Dist *array, size_t len) { this->array = array; this->len = len; }

    ~SampledVECTOR3Array() { delete[] array; }

    VECTOR3 &operator[] (size_t idx) {
        assert( idx >=0 && idx < len );
        return array[idx].getSample();
    }
};

#endif


} // namespace edda

#endif // DATA_ARRAY_H_
