#include <iostream>
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "structured_grid_dataset.h"
#include "core/interpolator.h"
#include "core/shared_ary.h"

using namespace std;
using namespace edda;
using namespace edda::dist;

typedef Gaussian<float> fGaussian;

shared_ary<float> read_raw(string fname, int* dimension) {
  FILE * fIn;
  int totalNum;
  float *pData;

  cout << fname << endl;
  fIn = fopen(fname.c_str(), "rb");
  assert(fIn != NULL);
  totalNum = dimension[0] * dimension[1] * dimension[2];
  pData = new float[totalNum * 3];
  fread(pData, sizeof(float), totalNum*3, fIn);
  fclose(fIn);
  return shared_ary<float>(pData, totalNum*3);
}

shared_ary<fGaussian> load_gaussian_dataset(string meanfile, string stdfile, int *dim) {
  shared_ary<float> pMean = read_raw(meanfile, dim);
  shared_ary<float> pStd = read_raw(stdfile, dim);
  // build gaussian for x
  int i,j,k,count=0;
  shared_ary<fGaussian> data(new fGaussian[dim[0]*dim[1]*dim[2]], dim[0]*dim[1]*dim[2]);  // better to use managed array
  for (k=0; k<dim[2]; k++)
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++)
      {
        data[count] = fGaussian(pMean[count], pStd[count]);
        count++;
      }
  return data;
}


int main(int argc, char **argv) {
    edda::ReturnStatus r;

    cout << "isoProbField <mean file> <std file> <w> <h> <d>" << endl;
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

    shared_ary<fGaussian> distData = load_gaussian_dataset(meanfile, stdfile, dim);

    /////////////////////////////////////
    cout << "value interpolation after sampling:" << endl;

    StructuredGridDataset<float> dataset1 (
                new RegularCartesianGrid (dim[0], dim[1], dim[2]),
                new SampledDataArray<fGaussian> (distData )
    );

    cout << "value interpolation after sampling:" << endl;

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

    StructuredGridDataset<fGaussian> dataset2 (
                new RegularCartesianGrid (dim[0], dim[1], dim[2]),
                new GeneralDataArray<fGaussian> (distData) );

    for (x=0; x<10; x+=.5)
    {
        fGaussian sampled_dist;
        r = dataset2.at_phys( VECTOR3(x, 10, 10), sampled_dist );
        sampled_dist.print();
        cout << "\t" << sampled_dist.getSample() <<  endl;
    }
    cout << endl;



    return 0;
}
