/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Streamlines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef _HEADER_H_
#define _HEADER_H_

namespace edda {

const double	DEG_TO_RAD = 0.0174532925199432957692;	// PI / 180
const double	RAD_TO_DEG = 57.2957795130823208768;	// 180 / PI
const double	PIBY2 = 1.57079632679489661923;			// PI / 2
const double	EPS = 1.0E-6;
const int		OCT = 8;

enum ReturnStatus { SUCCESS=0, FAIL };

#if defined( WIN32 )
typedef long long int64_t;
#endif

} // namespace edda

#endif
