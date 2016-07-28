// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DISTR_ARRAY_H_
#define DISTR_ARRAY_H_

#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "abstract_distr_array.h"
#include <core/vector_matrix.h>
#include <core/shared_ary.h>
#include "distributions/distribution.h"

namespace edda {

//---------------------------------------------------------------------------------------
/// \brief A simple implementation of DataArray.
///
/// This is the class that holds the actual array, in a smart array.
/// \param T The element type of an array.
///
template<typename Distr, ENABLE_IF_BASE_OF(Distr, dist::DistributionTag)>
class DistrArray: public AbstractDistrArray
{
protected:
  shared_ary<Distr> array;
public:
  DistrArray(shared_ary<Distr> array) { this->array = array; }

  virtual ~DistrArray() { }

  virtual size_t getLength() { return array.getLength(); }

  virtual int getNumComponents() { return 1; }

  virtual dist::Variant getDistr(size_t idx) { return array[idx]; }

  virtual std::vector<dist::Variant> getDistrVector(size_t idx) { return std::vector<dist::Variant> (1, array[idx] ); }

  virtual Real getScalar(size_t idx) { return getSample(array[idx]); }

  virtual std::vector<Real> getVector(size_t idx) { return std::vector<Real> (1, Real(getSample(array[idx])) ); }

  //virtual boost::any getItem(size_t idx) { return boost::any( array[idx] );  }

  //virtual void setItem(size_t idx, int component, const boost::any &item) { array[idx] = boost::any_cast<T>( item );  }

  virtual boost::any getRawArray() { return boost::any(array); }

  virtual std::string getDistrName() { return getName(Distr()); }
};


//---------------------------------------------------------------------------------------

template<typename Distr, int N, ENABLE_IF_BASE_OF(Distr, dist::DistributionTag)>
class DistrVectorArray: public AbstractDistrArray
{
protected:
  shared_ary<Vector<Distr,N> > array;
public:
  DistrVectorArray(shared_ary<Vector<Distr,N> > array) { this->array = array; }

  virtual ~DistrVectorArray() { }

  virtual size_t getLength() { return array.getLength(); }

  virtual int getNumComponents() { return N; }

  virtual dist::Variant getDistr(size_t idx) {
      throw std::runtime_error("Requesting scalar in a VectorArray");
    }

  virtual std::vector<dist::Variant> getDistrVector(size_t idx) {
    std::vector<dist::Variant> v(N);
    for (int i=0; i<N; i++) {
      v[i] = array[idx][i];
    }
    return v;
  }

  virtual Real getScalar(size_t idx) {
    throw std::runtime_error("Requesting scalar in a VectorArray");
  }

  virtual std::vector<Real> getVector(size_t idx) {
    std::vector<Real> v(N);
    for (int i=0; i<N; i++) {
      v[i] = getSample(array[idx][i]);
    }
    return v;
  }

  //virtual boost::any getItem(size_t idx) { return boost::any( array[idx] );  }

  //virtual void setItem(size_t idx, int component, const boost::any &item) { array[idx] = boost::any_cast<Vector<T,N> >( item );  }

  virtual boost::any getRawArray() { return boost::any(array); }

  virtual std::string getDistrName() { return getName(Distr()); }
};



} // namespace edda

#endif // DISTR_ARRAY_H_
