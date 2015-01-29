// Copyright (c) 2014 The EDDA Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_GAUSSIAN_H
#define DIST_GAUSSIAN_H

#include <cmath>
#include <cstdlib>
#include "edda_export.h"
#include "statistics.h"

namespace edda {

/// Defines a Gaussian distribution class
class EDDA_EXPORT Gaussian {
protected:
    double mean;
    double std;
public:
    // construction
    explicit Gaussian(): mean(0), std(1) {}
    explicit Gaussian(double m, double s): mean(m), std(s) {}

    // assign
    Gaussian& operator=(const Gaussian& x) {
        mean = x.mean;  std = x.std;
        return *this;
    }

    // assign
    template<class Distribution>
    Gaussian& operator=(const Distribution& x) {
        mean = x.getMean(); std = x.getStd();
        return *this;
    }

    ~Gaussian() {}

    // get probability
    inline double getProb(const double x) const {
        return exp( -0.5 * (x-mean)*(x-mean) / std / std ) /
                (std* sqrt(2.*M_PI));
    }
    inline double getMean() const {
        return mean;
    }
    inline double getVariance() const {
        return std*std;
    }
    inline double getStd() const {
        return std;
    }

    inline Gaussian& operator+=(const Gaussian& rhs) {
        mean += rhs.mean;
        std = sqrt(std*std + rhs.std*rhs.std);
        return *this;
    }

    inline Gaussian& operator-=(const Gaussian& rhs) {
        mean -= rhs.mean;
        std = sqrt(std*std + rhs.std*rhs.std);
        return *this;
    }
    inline Gaussian& operator+=(const double r) {
        mean += r;
        return *this;
    }
    inline Gaussian& operator*=(const double r) {
        mean *= r;
        std *= r;
        return *this;
    }
};

} // namespace itl

#endif // DIST_GAUSSIAN_H
