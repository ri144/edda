#include <iostream>
#include "distributions/Gaussian.h"
#include "distributions/Distribution.h"

using namespace std;

int main()
{
    edda::dist::Gaussian<> g1(1.,1.),
            g2(2.,2.);
    edda::dist::Gaussian<> g3;
    g3=g1+g2;
    cout << "g1+g2= <" << g3.getMean() << "," << g3.getStd() << ">" << endl;

    cout << "cdf(g1,2)=" << edda::dist::cdf(g1, 2) << endl;

    cout << "A random sample of g3: " << g3.getSample() << endl;

    cout << "size of Gaussian: " << sizeof(g1) << endl;

    edda::dist::Distribution<double> d1;
    cdf(d1, 2);

    return 0;
}
