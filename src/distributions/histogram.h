// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

/// Histogram: TODO
///

#include <iostream>
#include <typeinfo>
#include <iostream>

#include "common.h"
#include "edda_export.h"
#include "core/vector_matrix.h"
#include "core/tuple.h"
#include "core/thrust_common.h"
#include "distribution_tag.h"
#include "core/shared_ary.h"

namespace edda {
namespace dist {

///
/// \brief The Distribution class is a root class for all distribution-type classes.
///
/// This is useful for applications to identify whether a class is a distribution-type class, by using ENABLE_IF_BASE_OF()
///
class EDDA_EXPORT Histogram : public DiscreteDistributionTag {
public:

  Histogram(){}

  Histogram(double* histData){
	  m_nBins = histData[0];
	  m_minValue = histData[1];
	  m_maxValue = histData[2];
	  m_binWidth = (m_maxValue - m_minValue) / (Real)(m_nBins);

	  m_cdf.resize(m_nBins);
	  for (int b = 0; b < m_nBins; b++){
		  m_cdf[b] = histData[b+3];
	  }
  }

  Histogram(float *dataPoints, int points, const Real _minValue, const Real _maxValue, const int _nBins){
    m_nBins = _nBins;

    m_minValue = _minValue;
    m_maxValue = _maxValue;

    m_binWidth = (m_maxValue - m_minValue) / (Real)(m_nBins);

    m_cdf.resize(m_nBins);

    //modeling and convert to cdf
    for (int i = 0; i < points; i++){
      int b = valueToBinsIndex(dataPoints[i]);
      m_cdf[b] ++;
    }

    for (int b = 1; b < m_nBins; b++)
      m_cdf[b] += m_cdf[b - 1];

    for (int b = 0; b < m_nBins; b++){
      m_cdf[b] /= (Real)points;
    }

  }

  Real getMean() const{
    Real mean = 0;
    Real cValue = m_minValue + m_binWidth / 2.0;
    for (int i = 0; i < m_nBins; i++){
      if (i == 0) mean += m_cdf[0] * cValue;
      else mean += (m_cdf[i] - m_cdf[i - 1]) * cValue;
      cValue += m_binWidth;
    }
    return mean;
  }

  Real getVar() const{
    //using Var(X) = E[X^2] - (E[X])^2  to compute Var
    Real mean = 0;
    Real mSqrt = 0;
    Real cValue = m_minValue + m_binWidth / 2.0;
    for (int i = 0; i < m_nBins; i++){
      if (i == 0) mean += m_cdf[0] * cValue;
      else mean += (m_cdf[i] - m_cdf[i - 1]) * cValue;

      if (i == 0) mSqrt += m_cdf[0] * cValue * cValue;
      else mSqrt += (m_cdf[i] - m_cdf[i - 1]) * cValue * cValue;

      cValue += m_binWidth;
    }
    return mSqrt - mean*mean;
  }

  Real getPdf(const double x) const{
    int b = valueToBinsIndex(x);

    if (b == 0) return m_cdf[0];
    else return (m_cdf[b] - m_cdf[b - 1]);
  }

  Real getCdf(const double x) const{
    int b = valueToBinsIndex(x);

    return m_cdf[b];
  }

  Real getSample() const{
    Real sample;
    Real r = rand() / (float)RAND_MAX;

    int low = 0;
    int high = m_nBins - 1;

    //bineary searh in cdf
    do{
      int mid = floor((low + high) / 2.0);

      Real prevCdf;
      if (mid == 0) prevCdf = -0.1;
      else prevCdf = m_cdf[mid - 1];

      if (prevCdf < r && r <= m_cdf[mid]){
        //sample here
        Real subr = rand() / (float)RAND_MAX;
        Real s, t;
        rangeOfBin(s, t, mid);
        sample = s + m_binWidth * subr;
        break;
      }
      else if (m_cdf[mid] > r){
        high = mid - 1;
      }
      else{
        low = mid + 1;
      }
    } while (1);

    return sample;
  }

  void output(std::ostream& os) const{
    for (int b = 0; b < m_nBins; b++)
      os << "Bin: " << b << ": " << m_cdf[b] << std::endl;
  }

  int getBins(){
	  return m_nBins;
  }

  float getMaxValue(){
	  return m_maxValue;
  }

  float getMinValue(){
	  return m_minValue;
  }

  
  float getBinValue(int b){
	  //This usually return accumlative prob
	  return m_cdf[b];
  }

private:
  int valueToBinsIndex(Real v) const {
    if (v == m_maxValue)
      return (m_nBins - 1);
    else
      return floor((v - m_minValue) / m_binWidth);
  }

  void rangeOfBin(Real& s, Real &t, int b) const{
    s = m_minValue + b * m_binWidth;
    t = m_minValue + (b + 1) * m_binWidth;
  }

  int m_nBins;
  Real m_minValue;
  Real m_maxValue;
  Real m_binWidth;
  std::vector<Real> m_cdf;
};

///
/// \brief Compute the mean of the distribution
///
inline double getMean(const Histogram &dist)
{
  return (double)dist.getMean();
}

///
/// \brief Compute Variance
///
inline double getVar(const Histogram &dist)
{
  return (double)dist.getVar();
}

///
/// \brief Return PDF of x
///
inline double getPdf(const Histogram &dist, const double x)
{
  return (double)dist.getPdf(x);
}

///
/// \brief Get a Monte-Carlo sample of the distribution. We rely on the specific implementation from each distribution.
///
inline double getSample(const Histogram &dist)
{
  return (double)dist.getSample();
}


///
/// \brief Return CDF of x
///
//__host__ __device__
inline double getCdf(const Histogram &dist, double x)
{
  return (double)dist.getCdf(x);
}

///
/// \brief Print itself
///
inline std::ostream& operator<<(std::ostream& os, const Histogram &dist)
{
  dist.output(os);
  return os;
}

__host__ __device__
inline const char *getName(const Histogram &x) {
  return "Histogram";
}

}  // namespace dist
}  // namespace edda

#endif // HISTOGRAM_H
