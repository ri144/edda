#include <iostream>

#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "dataset.h"
#include "io/file_reader.h"
#include "core/interpolator.h"
#include "core/shared_ary.h"
#include "core/data_array.h"

using namespace std;
using namespace edda;
using namespace edda::dist;

typedef Gaussian<float> Gaussianf;
typedef Vector3<Gaussianf> Gaussianf3;

shared_ary<Gaussianf3> load_gaussian_dataset(string meanfile, string stdfile, int *dim) {
  cout << "Loading..." << endl;
  size_t len = dim[0]*dim[1]*dim[2];
  shared_ary<Tuple3<float> > pMean3 = loadRawFile<Tuple3<float> >(meanfile, len);
  shared_ary<Tuple3<float> > pStd3 = loadRawFile<Tuple3<float> >(stdfile, len);

  // Create Gaussian array
  Gaussianf3 x;
  shared_ary<Gaussianf3> data(new Gaussianf3[len], len);
  int i;
  for (i=0; i<len; i++)
  {
    data[i] = Gaussianf3(
                Gaussianf(pMean3[i][0], pStd3[i][0]),
                Gaussianf(pMean3[i][1], pStd3[i][1]),
                Gaussianf(pMean3[i][2], pStd3[i][2]) );
  }
  return data;
}


int main(int argc, char **argv) {
    ReturnStatus r;

    cout << "Input arguments: <mean file> <std file> <w> <h> <d>" << endl;
    if (argc<6)
        return -1;
    string meanfile, stdfile;
    meanfile = argv[1];
    stdfile = argv[2];
    int dim[3];
    dim[0] = atoi(argv[3]);
    dim[1] = atoi(argv[4]);
    dim[2] = atoi(argv[5]);
    cout << "dim: " << dim[0] << "," << dim[1] << "," << dim[2] << endl;

    shared_ary<Gaussianf3> distData = load_gaussian_dataset(meanfile, stdfile, dim);

    /////////////////////////////////////
    cout << "value interpolation after sampling:" << endl;

    Dataset<VECTOR3> dataset1 (
                new RegularCartesianGrid (dim[0], dim[1], dim[2]),
                new SampledIndepTupleArray<Gaussianf3> (distData )
    );

    cout << "value interpolation after sampling:" << endl;

    float x;
    for (x=0; x<10; x+=.5)
    {
        VECTOR3 sampled_val;
        r = dataset1.at_phys( VECTOR3(x, 10, 10), sampled_val );
        cout << sampled_val << " ";
    }
    cout << endl;

    /////////////////////////////////////
    cout << "value sampling after distribution interpolation :" << endl;

    Dataset<Gaussianf3> dataset2 (
                new RegularCartesianGrid (dim[0], dim[1], dim[2]),
                new DataArray<Gaussianf3> (distData) );

    for (x=0; x<10; x+=.5)
    {
        Gaussianf3 sampled_vec_dist;
        r = dataset2.at_phys( VECTOR3(x, 10, 10), sampled_vec_dist );
        cout << sampled_vec_dist << " ";
        cout << "\t" << dist::getSample(sampled_vec_dist[0]) << " "
                     << dist::getSample(sampled_vec_dist[1]) << " "
                     << dist::getSample(sampled_vec_dist[2]) <<  endl;
    }
    cout << endl;



    return 0;
}
