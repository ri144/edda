/*************************************************************************
*************************************************************************/

#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include "Tuple.h"

namespace edda {

// barycentric interpolation
template<class T>
inline T baryInterp(const T &v1, const T &v2, const T &v3, const T &v4, double coeff[3])
{
    return v1 +
            (v2-v1)*coeff[0] +
            (v3-v1)*coeff[1] +
            (v4-v1)*coeff[2];
}

// linear interpolation
template <class T>
inline T lerp(const T &x, const T &y, const double ratio)
{
    return x * (1. - ratio) + y * ratio;
}

// bilinear interpolation
template <class T>
inline T biLerp(const T &ll, const T &hl, const T &lh, const T &hh, const double coeff[2])
{
    return Lerp(Lerp(ll, hl, coeff[0]), lerp(lh, hh, coeff[0]), coeff[1]);
}

// trilinear interpolation
template <class T>
inline T triLerp(const T &lll, const T &hll, const T &lhl, const T &hhl,
                 const T &llh, const T &hlh, const T &lhh, const T &hhh, const double coeff[3])
{
    return Lerp(BiLerp(lll, hll, lhl, hhl, coeff),
                 BiLerp(llh, hlh, lhh, hhh, coeff),
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

} // namespace edda

#endif
