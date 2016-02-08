// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef SAMPLING_DATA_ARRAY_H_
#define SAMPLING_DATA_ARRAY_H_

#include <cassert>
#include <iostream>
#include <memory>
#include "abstract_data_array.h"
#include "vector_matrix.h"
#include "shared_ary.h"
#include "distributions/distribution.h"

namespace edda {

//---------------------------------------------------------------------------------------
/// \brief Take a distribution data array and output a sample
/// \param Dist The distribution type of the input data array.
///
template<typename Dist, ENABLE_IF_BASE_OF(Dist, dist::Distribution)>
class SamplingDataArray: public AbstractDataArray
{
protected:
  std::shared_ptr<AbstractDataArray> array;
public:
  SamplingDataArray(std::shared_ptr<AbstractDataArray> array_): array(array_) { }

  virtual ~SamplingDataArray() { }

  virtual boost::any getItem(size_t idx) {
    return boost::any(
          getSample( boost::any_cast<Dist>(array->getItem(idx) ) )
          );
  }

  virtual void setItem(size_t idx, const boost::any &item) { array->setItem( idx, boost::any_cast<Dist>( item ) );  }

  virtual size_t getLength() { return array->getLength(); }

  virtual boost::any getRawArray() { return array->getRawArray(); }
};


} // namespace edda

#endif // SAMPLING_DATA_ARRAY_H_
