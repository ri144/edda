#include <iostream>

#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "dataset/dataset.h"
#include "io/file_reader.h"
#include "core/interpolator.h"
#include "core/shared_ary.h"

using namespace std;
using namespace edda;
using namespace edda::dist;

int main(int argc, char **argv) {
    edda::ReturnStatus r;

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

    shared_ary<Gaussian> distData = edda::loadGaussianRawArray(meanfile, stdfile, dim[0]*dim[1]*dim[2]);

    /////////////////////////////////////
    cout << "value interpolation after sampling:" << endl;

    Dataset<float> dataset1 (
                new RegularCartesianGrid (dim[0], dim[1], dim[2]),
                new DataArray<Gaussian, GetItemSampled>(distData )
            );
    float x;
    for (x=0; x<10; x+=.5)
    {
        float sampled_val;
        r = dataset1.at_phys( VECTOR3(x, 10, 10), sampled_val );
        cout << sampled_val << " ";
    }
    cout << endl;

    /////////////////////////////////////
    cout << "value sampling after distribution interpolation :" << endl;

    Dataset<Gaussian> dataset2 (
                new RegularCartesianGrid (dim[0], dim[1], dim[2]),
                new DataArray<Gaussian> (distData)
            );

    for (x=0; x<10; x+=.5)
    {
        Gaussian sampled_dist;
        r = dataset2.at_phys( VECTOR3(x, 10, 10), sampled_dist );
        cout << sampled_dist <<  "\t" << dist::getSample(sampled_dist) <<  endl;
    }
    cout << endl;



    return 0;
}
