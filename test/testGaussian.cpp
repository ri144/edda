#include <iostream>
#include "distributions/Gaussian.h"
#include "distributions/Distribution.h"
#include "Interpolator.h"

using namespace std;

int main()
{
    edda::dist::Gaussian<double> g1(1.,1.);
    edda::dist::Gaussian<double> g2(2.,2.);
    edda::dist::Gaussian<double> g3;

    cout << "g1= "; g1.print(); cout << endl;
    cout << "g1= "; g2.print(); cout << endl;

    // adding two random variables
    g3=g1+g2;

    cout << "g1+g2= "; g3.print(); cout << endl;
    
    // linear interpolation
    g3 =  edda::lerp(g1, g2, .1);
    cout << "lerp(g1, g2, .1) = "; g3.print();  cout << endl;

    // 1.96 std from the mean (one side) will cover 97.5% of the distribution
    cout << "cdf(g1,2.96)=" << edda::dist::cdf(g1, 2.96) << endl;

    cout << "A random sample of g3: " << g3.getSample() << endl;

    cout << "size of Gaussian: " << sizeof(g1) << endl;

    edda::dist::Distribution<double> d1;
    //cdf(d1, 2); // generic function is currently not implmented

    return 0;
}
