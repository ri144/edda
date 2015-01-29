#include <iostream>
#include "distributions/Gaussian.h"
#include "distributions/Distribution.h"

using namespace edda;
using namespace std;

int main()
{
    Gaussian g1(1.,1.),
            g2(2.,2.);
    Gaussian g3 = g1;
    g3+=g2;
    cout << "g1+g2= <" << g3.getMean() << "," << g3.getStd() << ">" << endl;

    Distribution<Gaussian> d1(g1), d2(g2);
    Distribution<Gaussian> d3;
    d3 = d1+d2;
    cout << "d1+d2= <" << d3.getMean() << "," << d3.getStd() << ">" << endl;




    return 0;
}
