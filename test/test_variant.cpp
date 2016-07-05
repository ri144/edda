#include <cstdlib>
#include <ctime>
#include <iostream>

#include <boost/variant/variant.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>

#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>

#include <boost/format.hpp>
#include "distributions/distribution.h"
#include "distributions/variant.h"
using namespace edda;
using namespace std;
const int iterations_count = 10; //0000000;

#define TEST_NEW

double use_virtual() {

    dist::DefaultGaussianMixture gmm;
    double out;
    dist::Distribution *dist = new dist::GaussianMixtureWrapper(gmm);

    for (int i = 0; i < iterations_count; i++) {
#ifdef TEST_NEW
        dist = new dist::GaussianMixtureWrapper(gmm);
#endif
        //dist->operator<<( std::cout ) << std::endl;
        cout << sizeof(dist::GaussianMixtureWrapper) << endl;
        out = dist->getSample();
    }
    return out;
}

double use_variant() {
  dist::DefaultGaussianMixture gmm;
  double out;
  dist::Variant dist = gmm;

    for (int i = 0; i < iterations_count; i++) {
#ifdef TEST_NEW
        dist = gmm;
#endif
        cout << sizeof(dist) << endl;
        out = getSample(dist);
    }
    return out;
}

int main() {
    using namespace boost::posix_time;
    srand(time(NULL));

    ptime start, end;
    time_duration d1, d2;

    double out;

    // virtual
    start = microsec_clock::universal_time();
    out = use_virtual();
    end = microsec_clock::universal_time();
    std::cout << out << std::endl;

    // store result
    d1 = end - start;

    // variant
    start = microsec_clock::universal_time();
    out = use_variant();
    end = microsec_clock::universal_time();
    std::cout << out << std::endl;

    // store result
    d2 = end - start;

    // output
    std::cout <<
        boost::format(
            "Virtual: %1%\n"
            "Variant: %2%\n"
        ) % d1 % d2;
}
