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
const int iterations_count = 100000000;

double use_virtual() {

    dist::DefaultGaussianMixture gmm;
    double out;

    for (int i = 0; i < iterations_count; i++) {
        dist::Distribution *dist = new dist::GaussianMixtureWrapper(gmm);
        //dist->operator<<( std::cout ) << std::endl;
        out = dist->getSample();
    }
    return out;
}

double use_variant() {
  dist::DefaultGaussianMixture gmm;
  double out;

    for (int i = 0; i < iterations_count; i++) {

        dist::Variant var = gmm;
        out = getSample(gmm);
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
