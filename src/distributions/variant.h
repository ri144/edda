#include <boost/variant.hpp>

#include "distribution.h"
#include "gaussian.h"
#include "gaussian_mixture.h"

namespace edda{
namespace dist{

  typedef boost::variant<Real, Gaussian, GaussianMixture> _Variant;

  struct Variant : public _Variant, public Distribution {
    Variant() : _Variant() {}
    Variant(const Real &obj) : _Variant (obj) {}
    Variant(const Gaussian &obj) : _Variant (obj) {}
    Variant(const GaussianMixture &obj) : _Variant (obj) {}
  };

  namespace detail{
    struct _getPdf : public boost::static_visitor<double> {
      double x;
      _getPdf(double _x) : x(_x) {}
      template <class T> inline double operator() (const T& dist) { return getPdf(dist, x); }
      inline double operator() (const Real& value) { return (x==value)?1.:0; }
    };
    struct _getCdf : public boost::static_visitor<double>  {
      double x;
      _getCdf(double _x) : x(_x) {}
      template <class T> inline double operator() (const T& dist) { return getCdf(dist, x); }
      inline double operator() (const Real& value) { return (x<=value)?1.:0; }
    };
    struct _getMean : public boost::static_visitor<double> {
      template <class T> inline double operator() (const T& dist) { return getMean(dist); }
      inline double operator() (const Real& value) { return value; }
    };
    struct _getVar: public boost::static_visitor<double> {
      template <class T> inline double operator() (const T& dist) { return getVar(dist); }
      inline double operator() (const Real& value) { return 0; }
    };
    struct _getSample : public boost::static_visitor<double> {
      template <class T> inline double operator() (const T& dist) { return getSample(dist); }
      inline double operator() (const Real& value) { return value; }
    };
#if 0
    struct _plus_assign: public boost::static_visitor<double> {
      template <class T> inline double operator() (const T& lhs, const T& rhs) { lhs += rhs; return lhs; }
    };
    struct _subtract_assign: public boost::static_visitor<Variant> {
      template <class T> inline Variant operator() (const T& lhs, const T& rhs) { lhs -= rhs; return lhs; }
    };
#endif
  } // namespace detail

  inline double getPdf(const Variant &dist, double x) {
    detail::_getPdf f(x);
    return boost::apply_visitor( f, dist);
  }
  inline double getCdf(const Variant &dist, double x) {
    detail::_getCdf f(x);
    return boost::apply_visitor( f, dist );
  }
  inline double getMean(const Variant &dist)  {
    detail::_getMean f;
    return boost::apply_visitor( f, dist);
  }
  inline double getVar(const Variant &dist)  {
    detail::_getVar f;
    return boost::apply_visitor( f, dist);
  }
  inline double getSample(const Variant &dist) {
    detail::_getSample f;
    return boost::apply_visitor( f, dist);
  }
#if 0
  inline Variant operator-=(const Variant& lhs, const Variant& rhs) {
    detail::_subtract_assign f;
    return boost::apply_visitor( f, lhs, rhs);
  }

  ///
  /// \brief random variable *
  ///
  template<class T, ENABLE_IF_BASE_OF(T, Distribution) >
  inline T operator*(const T& lhs, const double x) {
      T h(lhs);
      return h *= x;
  }
#endif


  ///
  /// \brief Return a vector sample
  ///
  template <class Dist, int N, ENABLE_IF_BASE_OF(Dist, Distribution) >
  inline Vector<Real, N> getSample(const Vector<Dist, N> &v)
  {
    Vector<Real, N> out;
    for (int i=0; i<N; i++)
      out[i] = getSample(v[i]);
    return out;
  }

} // namespace dist
} // namespace edda
