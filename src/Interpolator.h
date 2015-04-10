/*************************************************************************
*************************************************************************/

#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include <vector>

namespace edda {

// barycentric interpolation
template<class T>
inline T BaryInterp(const T &v1, const T &v2, const T &v3, const T &v4, double coeff[3])
{
    return v1 +
            (v2-v1)*coeff[0] +
            (v3-v1)*coeff[1] +
            (v4-v1)*coeff[2];
}

// linear interpolation
template <class T>
inline T Lerp(const T &x, const T &y, const double ratio)
{
    return x * (1. - ratio) + y * ratio;
}

// bilinear interpolation
template <class T>
inline T BiLerp(const T &ll, const T &hl, const T &lh, const T &hh, const double coeff[2])
{
    return Lerp(Lerp(ll, hl, coeff[0]), lerp(lh, hh, coeff[0]), coeff[1]);
}

// trilinear interpolation
template <class T>
inline T TriLerp(const T &lll, const T &hll, const T &lhl, const T &hhl,
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
class CubeInterpolator
{
    T operator() (
            const std::vector<T> &vData,
            const double coeff[3]
            )
    {
       return TriLerp(vData[0], vData[1], vData[2], vData[3],
                     vData[4], vData[5], vData[6], vData[7],
                     coeff);
    }
};

template<class T>
class TetraInterpolator
{
    T operator () (
            const std::vector<T> &vData,
            const double coeff[3]
        )
    {
        // not tested
        return BaryInterp(&vData[0], coeff);
    }
};

// Use when T in edda:dist types
template<class T>
class StochasticCubeInterpolator
{
    double operator() (
            const std::vector<T> &vData,
            const double coeff[3]
            )
    {
       return TriLerp(vData[0].getSample(), vData[1].getSample(), vData[2].getSample(), vData[3].getSample(),
                     vData[4].getSample(), vData[5].getSample(), vData[6].getSample(), vData[7].getSample(),
                     coeff);
    }
};

template<class T>
class StochasticTetraInterpolator
{
    double operator () (
            const std::vector<T> &vData,
            const double coeff[3]
        )
    {
        // not tested
        return BaryInterp(vData[0].getSample(), vData[1].getSample(), vData[2].getSample(), vData[3].getSample(), coeff);
    }
};

} // namespace edda

#endif
