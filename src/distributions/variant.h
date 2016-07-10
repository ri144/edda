#ifndef VARIANT_H
#define VARIANT_H

/// Experimental functionality
//#define BOOST_VARIANT_USE_RELAXED_GET_BY_DEFAULT

#include <boost/variant.hpp>

#include "distribution_tag.h"
#include "gaussian.h"
#include "gaussian_mixture.h"
//#include "histogram.h"

namespace edda{
namespace dist{

  typedef boost::variant<Real, Gaussian, DefaultGaussianMixture  > _Variant;

  struct Variant : public _Variant, public DistributionTag {
    Variant() : _Variant() {}
    Variant(const Real &obj) : _Variant (obj) {}
    Variant(const Gaussian &obj) : _Variant (obj) {}
    Variant(const DefaultGaussianMixture &obj) : _Variant (obj) {}
    //Variant(const Histogram &obj) : _Variant (obj) {}
  };

  namespace detail{
    struct _getPdf : public boost::static_visitor<double> {
      double x;
      template <class T> inline double operator() (const T& dist) { return getPdf(dist, x); }
    };
    struct _getCdf : public boost::static_visitor<double>  {
      double x;
      template <class T> inline double operator() (const T& dist) { return getCdf(dist, x); }
    };
    struct _getMean : public boost::static_visitor<double> {
      template <class T> inline double operator() (const T& dist) { return getMean(dist); }
    };
    struct _getVar: public boost::static_visitor<double> {
      template <class T> inline double operator() (const T& dist) { return getVar(dist); }
    };
    struct _getSample : public boost::static_visitor<double> {
      template <class T> inline double operator() (const T& dist) { return getSample(dist); }
    };
    struct _getName : public boost::static_visitor<const char *> {
      template <class T> inline const char *operator() (const T& dist) { return getName(dist); }
    };
  } // namespace detail

  inline double getPdf(const Variant &dist, double x) {
    detail::_getPdf f; f.x = x;
    return boost::apply_visitor( f, dist );
  }
  inline double getCdf(const Variant &dist, double x) {
    detail::_getCdf f; f.x = x;
    return boost::apply_visitor( f, dist );
  }
  inline double getMean(const Variant &dist)  {
    detail::_getMean f;
    return boost::apply_visitor( f, dist );
  }
  inline double getVar(const Variant &dist)  {
    detail::_getVar f;
    return boost::apply_visitor( f, dist );
  }
  inline double getSample(const Variant &dist) {
    detail::_getSample f;
    return boost::apply_visitor( f, dist );
  }  
  inline const char *getName(const Variant &dist) {
    detail::_getName f;
    return boost::apply_visitor( f, dist );
  }


  ///
  /// \brief Return a vector sample
  ///
  template <class Dist, int N, ENABLE_IF_BASE_OF(Dist, DistributionTag) >
  inline Vector<Real, N> getSample(const Vector<Dist, N> &v)
  {
    Vector<Real, N> out;
    for (int i=0; i<N; i++)
      out[i] = getSample(v[i]);
    return out;
  }

} // namespace dist
} // namespace edda

#endif // #define VARIANT_H
