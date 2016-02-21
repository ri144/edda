#include <iostream>

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <iostream>

#include "core/shared_ary.h"
#include "distributions/distribution.h"
#include "distributions/gaussian.h"
#include "core/interpolator.h"
#include "test_common.h"
#include "filters/random_sample_field.h"


using namespace std;
using namespace edda;

int main(int argc, char* argv[]) {
  dist::Gaussian g1(1.,1.);
  dist::Gaussian g2(2.,2.);
  dist::Gaussian g3=g1+g2;



  thrust::host_vector<dist::Gaussian> h_array(3);
  h_array[0] = g1;
  h_array[1] = g2;
  h_array[2] = g3;
  thrust::device_vector<dist::Gaussian> d_array = h_array;

  thrust::device_vector<Real> d_out;

  randomSampleField(d_array, d_out);

  thrust::host_vector<Real> h_out = d_out;

  for (int i=0; i<3; i++)
    cout << "Output: " << h_out[i] << endl;


  return 0;
}
