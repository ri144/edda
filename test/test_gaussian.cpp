#include <iostream>

#include "distributions/distribution.h"
#include "distributions/gaussian.h"
#include "core/interpolator.h"
#include "test_common.h"

using namespace std;

int main(int argc, char* argv[]) {
  edda::dist::Gaussian<double> g1(1.,1.);
  edda::dist::Gaussian<double> g2(2.,2.);
  edda::dist::Gaussian<double> g3;

  cout << "g1= "<< g1 << endl;
  cout << "g1= "<< g2 << endl;

  // adding two random variables
  g3=g1+g2;

  cout << "g1+g2= " << g3 << endl;
  TEST( g1.getMean() == edda::dist::getMean(g1) );
  TEST( approx_equal(g3.getMean(), 3. ) );
  TEST( approx_equal(g3.getStd(), 2.23606797749979 ) );

  // linear interpolation
  g3 =  edda::lerp(g1, g2, .1);

  cout << "lerp(g1, g2, .1) = " << g3 << endl;
  TEST( approx_equal(g3.getMean(), 1.1));
  TEST( approx_equal(g3.getStd(), 0.921954));


  // 1.96 std from the mean (one side) will cover 97.5% of the distribution
  double cdf = edda::dist::getCdf(g1, g1.getMean()+1.96);
  cout << "cdf(g1,2.96)=" << cdf << endl;
  TEST( approx_equal(cdf, 0.975002) );

  cout << "A random sample of g3: " << g3.getSample() << endl;

  cout << "size of Gaussian: " << sizeof(g1) << endl;
  TEST( sizeof(g1) == sizeof(double)*2 );

  //edda::dist::Distribution d1;
  //edda::dist::cdf(d1, 2.); // generic function is currently not implmented

  // edda::dist::cdf(3, 2.); // This will not pass compilation

  return 0;
}
