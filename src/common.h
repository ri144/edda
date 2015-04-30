/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Streamlines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef COMMON_H
#define COMMON_H

#include <limits>

namespace edda {

const double	DEG_TO_RAD = 0.0174532925199432957692;	// PI / 180
const double	RAD_TO_DEG = 57.2957795130823208768;	// 180 / PI
const double	PIBY2 = 1.57079632679489661923;			// PI / 2

// Do not use this in the future:
const double	EPS = 1.0E-6;
// Instead call this function:  std::numeric_limits<T>::epsilon()
// or use this C++11 alias
template <typename T> using limits = std::numeric_limits<T>;

enum ReturnStatus { SUCCESS=0, FAIL };

#if defined( WIN32 )
typedef long long int64_t;
#endif

// May not be needed in the release version but just keep here for now
class NotImplementedException
{};

} // namespace edda

#endif  // COMMON_H
