// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_GAUSSIAN_H_
#define DIST_GAUSSIAN_H_

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <boost/math/distributions.hpp>

#include "edda_export.h"
#include "distributions/distribution.h"
#include "math/statistics.h"

namespace edda {
namespace dist {

/// Defines a Gaussian distribution class
template <class Real = float>
class EDDA_EXPORT Gaussian: public ContinuousDistribution {
    Real mean, std;
public:
    typedef Real real_type;
    // construction
    explicit Gaussian(): mean(0), std(1) {}
    explicit Gaussian(Real m, Real s): mean(m), std(s) {}

    ~Gaussian() {}

    // get probability (pdf)
    inline double getProb(const double x) const {
        if (std==0) {
            return (fabs(x-mean)< std::numeric_limits<Real>::epsilon() )? 1: 0;
        }
        return exp( -0.5 * (x-mean)*(x-mean) / std / std ) /
            (std* sqrt(2.*M_PI));
    }
    inline double getSample() const {
        return box_muller(mean, std);
    }


    // random variable +=
    inline Gaussian& operator+=(const Gaussian& rhs) {
        mean += rhs.mean;
        std = sqrt(std*std + rhs.std*rhs.std);
        return *this;
    }
    // random variable -=
    inline Gaussian& operator-=(const Gaussian& rhs) {
        mean -= rhs.mean;
        std = sqrt(std*std + rhs.std*rhs.std);
        return *this;
    }
    // // random variable +=
    inline Gaussian& operator+=(const double r) {
        mean += r;
        return *this;
    }
    // random variable *=
    inline Gaussian& operator*=(const double r) {
        mean *= r;
        std *= r;
        return *this;
    }

    // additional functions
    inline Real getMean() const {
        return mean;
    }
    inline Real getStd() const {
        return std;
    }
    inline Real getVar() const {
        return std*std;
    }

};

template <typename Real>
std::ostream& operator<<(std::ostream& os, const Gaussian<Real>& dist)
{
    os <<  "<Gaussian: mean=" << dist.getMean() << ", std=" << dist.getStd() << ">" ;
    return os;
}

// get CDF
template <class Real>
double getCdf(const Gaussian<Real> &dist, double x)
{
    boost::math::normal_distribution<Real> normal (dist.getMean(), dist.getStd());
    return boost::math::cdf<>(normal, x);
}

template <class Real>
double getMean(const Gaussian<Real> &dist)
{
    return dist.getMean();
}


}  // namespace dist
}  // namespace edda

#endif  // DIST_GAUSSIAN_H_
