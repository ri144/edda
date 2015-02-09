// Copyright (c) 2014 The EDDA Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_GAUSSIAN_H
#define DIST_GAUSSIAN_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <boost/math/distributions.hpp>
#include "edda_export.h"
#include "Distribution.h"
#include "statistics.h"
#include "header.h"

namespace edda { namespace dist {

/// Defines a Gaussian distribution class
template <class Real = float>
class EDDA_EXPORT Gaussian {
    Real mean, std;
public:
    typedef Real real_type;
    // construction
    explicit Gaussian(): mean(0), std(1) {}
    explicit Gaussian(Real m, Real s): mean(m), std(s) {}

    ~Gaussian() {}

    // get probability
    inline Real getProb(const double x) const {
        if (std==0) {
            return (fabs(x-mean)<EPS)? 1: 0;
        }
        return exp( -0.5 * (x-mean)*(x-mean) / std / std ) /
            (std* sqrt(2.*M_PI));
    }
    inline Real getMean() const {
        return mean;
    }
    inline Real getStd() const {
        return std;
    }
    inline Real getVariance() const {
        return std*std;
    }
    inline Real getSample() const {
        return box_muller(mean, std);
    }
    // for debugging
    inline void print() {
        std::cout <<  "<Gaussian: mean=" << mean << ", std=" << std << ">" ;
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
};


// get CDF
template <class Real>
double cdf(const Gaussian<Real> &dist, double x)
{
    boost::math::normal_distribution<Real> normal (dist.getMean(), dist.getStd());
    return boost::math::cdf<>(normal, x);
}


} } // namespace dist, edda

#endif // DIST_GAUSSIAN_H
