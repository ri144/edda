#ifndef ABSTRACT_DATA_ARRAY_H
#define ABSTRACT_DATA_ARRAY_H

#include <boost/any.hpp>

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

} // edda
#endif
