#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "distributions/Gaussian.h"
#include "distributions/Distribution.h"
#include "VectorMatrix.h"
#include "nrrd.h"
#include "Field.h"

using namespace std;
using namespace edda;

typedef dist::Gaussian<> Gaussian ;


float* read(string fname, int* dimension)
{

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
  return(pData);
}


Gaussian *loadData(string meanfile, string stdfile, int *dim)
{
    float *pMean = (float *)read(meanfile, dim);
    float *pStd = (float *)read(stdfile, dim);
    // build gaussian for x
    int i,j,k,count=0;
    Gaussian *pData = new Gaussian[dim[0]*dim[1]*dim[2]];
    for (k=0; k<dim[2]; k++)
        for (j=0; j<dim[1]; j++)
            for (i=0; i<dim[0]; i++)
            {
                pData[count] = Gaussian(pMean[count], pStd[count]);
                count++;
            }
    return pData;
}

int main(int argc, char **argv)
{
    cout << "isoProbField <mean file> <std file> <w> <h> <d> <iso-value>" << endl;
    string meanfile, stdfile;
    meanfile = argv[1];
    stdfile = argv[2];
    int dim[3];
    dim[0] = atoi(argv[3]);
    dim[1] = atoi(argv[4]);
    dim[2] = atoi(argv[5]);
    cout << "dim: " << dim[0] << "," << dim[1] << "," << dim[2] << endl;
    Gaussian *pDistData = loadData(meanfile, stdfile, dim);
    float iso = atoi(argv[6]);

    //Solution<Gaussian> *solution = new Solution<Gaussian>(&pDistData, dim[0]*dim[1]*dim[2], 0);
    //RegularCartesianGrid *grid = new RegularCartesianGrid(dim[0], dim[1], dim[2]);
    //StructuredField<Gaussian> field(solution, grid);

    float *probField = new float[dim[0]*dim[1]*dim[2]];
    int i,j,k;
    int count=0;
    for (k=0; k<dim[2]; k++)
        for (j=0; j<dim[1]; j++)
            for (i=0; i<dim[0]; i++)
            {
                probField[count] = pDistData[count].getProb(iso);
                count++;
            }
    FILE *fp = fopen("result.raw", "wb");
    fwrite(probField, count, sizeof(float), fp);
    fclose(fp);
    write_nrrd_3d("result.nrrd", "result.raw", dim[0], dim[1], dim[2], "float");
}
