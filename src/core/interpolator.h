// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef INTERPOLATOR_H_
#define INTERPOLATOR_H_

#include "core/tuple.h"

namespace edda {

// barycentric interpolation
template <class T>
inline T baryInterp(const T &v1, const T &v2, const T &v3, const T &v4, double coeff[3])
{
    return v1 +
            (v2-v1)*coeff[0] +
            (v3-v1)*coeff[1] +
            (v4-v1)*coeff[2];
}

// linear interpolation
template <class T, typename coeffT=float>
inline T lerp(const T &x, const T &y, const coeffT ratio)
{
    return x * (1. - ratio) + y * ratio;
}

// bilinear interpolation
template <class T, typename coeffT=float>
inline T biLerp(const T &ll, const T &hl, const T &lh, const T &hh, const coeffT coeff[2])
{
    return lerp(lerp(ll, hl, coeff[0]), lerp(lh, hh, coeff[0]), coeff[1]);
}

// trilinear interpolation
template <class T, typename coeffT=float>
inline T triLerp(const T &lll, const T &hll, const T &lhl, const T &hhl,
                 const T &llh, const T &hlh, const T &lhh, const T &hhh, const coeffT coeff[3])
{
    return lerp(biLerp(lll, hll, lhl, hhl, coeff),
                 biLerp(llh, hlh, lhh, hhh, coeff),
                 coeff[2]);
}

///////////////////////////////////////////////
/// The classes to be used in data model
////////////////////////////////////////////////


template<class T>
inline T cubeLerp(const Tuple8<T> points, const double coeff[3])
{
    return TriLerp(points[0], points[1], points[2], points[3],
                  points[4], points[5], points[6], points[7],
                  coeff);
}

template<class T>
inline T tetraLerp(const Tuple4<T> points, const double coeff[3])
{
    return baryInterp(points[0], points[1], points[2], points[3], coeff);
}



#if 0
// Used when T in edda:dist types
template<class T>
class StochasticCubeInterpolator
{
    static double operator() (
            const std::vector<T> &vData,
            const double coeff[3]
            )
    {
       return triLerp(vData[0].getSample(), vData[1][i].getSample(), vData[2][i].getSample(), vData[3][i].getSample(),
                     vData[4][i].getSample(), vData[5][i].getSample(), vData[6][i].getSample(), vData[7][i].getSample(),
                     coeff);
    }
};

template<class T>
class StochasticTetraInterpolator
{
    static double operator () (
            const std::vector<T> &vData,
            const double coeff[3]
        )
    {
        // not tested
        return baryInterp(vData[0].getSample(), vData[1].getSample(), vData[2].getSample(), vData[3].getSample(), coeff);
    }
};
#endif

}  // namespace edda

#endif  // INTERPOLATOR_H_
