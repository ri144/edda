#ifndef ABSTRACT_DISTR_ARRAY_H
#define ABSTRACT_DISTR_ARRAY_H

#include <string>
#include <vector>
#include <boost/any.hpp>
#include "distributions/variant.h"

namespace edda {

///
/// \brief The AbstractDataArray class used in the Dataset class
///
class AbstractDistrArray
{
public:

  virtual ~AbstractDistrArray() {}

  ///
  /// Get the number of elements
  ///
  virtual size_t getLength() =0;

  ///
  /// Get the number of element components.  This is used for elements in variable-length vector or matrix.
  ///
  virtual int getNumComponents() = 0;

  ///
  /// Set the target component index of vector data when calling getScalar()
  ///
  virtual void SetTargetComponent(int idx) = 0;

  ///
  /// Get the target component index of vector data
  ///
  virtual int GetTargetComponent() = 0;

  ///
  /// Get random sampling of scalar distribution
  ///
  virtual Real getScalar(size_t idx) =0;

  ///
  /// Get random sampling of vector distribution
  ///
  virtual std::vector<Real> getVector(size_t idx)=0;

  virtual dist::Variant getDistr(size_t idx) =0;

  virtual std::vector<dist::Variant> getDistrVector(size_t idx)=0;

  //virtual void setDistr(size_t idx, dist::Variant) =0;

  //virtual void setDistrVector(size_t idx, std::vector<dist::Variant>)=0;

  ///
  /// Get the array
  ///
  virtual boost::any getRawArray() = 0;

  ///
  /// Get gistribution name for data writer
  ///
  virtual std::string getDistrName() = 0;

};

} // edda
#endif
