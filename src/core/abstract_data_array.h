#ifndef ABSTRACT_DATA_ARRAY_H
#define ABSTRACT_DATA_ARRAY_H

#include <boost/any.hpp>
#include <boost/variant/variant_fwd.hpp>

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
  /// Get the number of element components.  This is used for elements in variable-length vector or matrix.
  ///
  virtual int getComponents() = 0;

  ///
  /// \brief Get data at the given index
  /// \param idx Index of array
  /// \param component Select the component of the array item (vector component).
  ///
  /// Note: We don't return reference because the return data can be derived from the original data.  Use setItem() to change the data content.
  ///
  virtual boost::any getItem(size_t idx, int component=0) =0;

  ///
  /// Set data at the given index
  ///
  virtual void setItem(size_t idx, int component, const boost::any &item) =0;

  ///
  /// Get the array
  ///
  virtual boost::any getRawArray() = 0;

};

} // edda
#endif
