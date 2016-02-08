#include <iostream>

#include "edda.h"
#include "distributions/gaussian.h"
#include "io/file_reader.h"

using namespace std;
using namespace edda;
using namespace edda::dist;

typedef Vector3<Gaussian> Gaussian3;


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

    shared_ary<Gaussian3> distData = edda::loadVec3GaussianRawArray(meanfile, stdfile, dim[0]*dim[1]*dim[2]);

    /////////////////////////////////////
    cout << "value interpolation after sampling:" << endl;

    Dataset<VECTOR3> dataset1 (
                new RegularCartesianGrid (dim[0], dim[1], dim[2]),
                new DataArray<Gaussian3, GetItemSampledVector> (distData )
    );

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

    Dataset<Gaussian3> dataset2 (
                new RegularCartesianGrid (dim[0], dim[1], dim[2]),
                new DataArray<Gaussian3> (distData) );

    for (x=0; x<10; x+=.5)
    {
        Gaussian3 sampled_vec_dist;
        r = dataset2.at_phys( VECTOR3(x, 10, 10), sampled_vec_dist );
        cout << sampled_vec_dist << " ";
        cout << "\t" << dist::getSample(sampled_vec_dist[0]) << " "
                     << dist::getSample(sampled_vec_dist[1]) << " "
                     << dist::getSample(sampled_vec_dist[2]) <<  endl;
    }
    cout << endl;



    return 0;
}
