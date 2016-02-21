// Copyright 2015 The Edda Authors.s All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef ABSTRACT_SAMPLING_ARRAY_H
#define ABSTRACT_SAMPLING_ARRAY_H

#include <cassert>
#include <iostream>
#include <memory>
#include "abstract_data_array.h"

namespace edda {

//---------------------------------------------------------------------------------------
/// \brief Take a distribution data array and output a sample
/// \param Dist The distribution type of the input data array.
///
template <class T>
class DataSamplingArray: public AbstractDataArray
{
protected:
  AbstractDataArray * array;
public:
  DataSamplingArray(AbstractDataArray *array_): array(array_) { }

  virtual ~DataSamplingArray() { if (array) delete array; }

  virtual size_t getLength() { return array->getLength(); }

  virtual int getNumComponents() { return array->getNumComponents(); }

  virtual dist::Variant getScalar(size_t idx) { return getSample( array->getScalar(idx) );  }

  virtual std::vector<dist::Variant> getVector(size_t idx) {
    std::vector<dist::Variant> v = array->getVector(idx);
    std::vector<dist::Variant> out;
    for (size_t i=0; i<v.size(); i++)
      out[i] = getSample(v[i]);
    return out;
  }

  ///
  /// getSample returns double or tuple types
  ///
  virtual boost::any getItem(size_t idx) { return dist::getSample (
            boost::any_cast<T>( array->getItem(idx) )
          ); }

  virtual void setItem(size_t idx, int component, const boost::any &item) { array->setItem( idx, component, item );  }

  virtual boost::any getRawArray() { return array->getRawArray(); }
};


} // namespace edda

#endif // ABSTRACT_SAMPLING_ARRAY_H
