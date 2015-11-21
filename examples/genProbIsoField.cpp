// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

// A probably wrong method for showing isosurface of uncertain data using pdf
// But for histogram data it can be fine

#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/nrrd.h"

using namespace std;
using namespace edda;

typedef dist::Gaussian<> Gaussian;

float* readraw(string fname, int* dimension) {
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

vector<Gaussian> loadData(string meanfile, string stdfile, int *dim) {
  float *pMean = (float *)readraw(meanfile, dim);
  float *pStd = (float *)readraw(stdfile, dim);
  // build gaussian for x
  int i,j,k,count=0;
  vector<Gaussian> pData(dim[0]*dim[1]*dim[2]);  // better to use managed array
  for (k=0; k<dim[2]; k++)
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++)
      {
        pData[count] = Gaussian(pMean[count], pStd[count]);
        count++;
      }
  return pData;
}

int main(int argc, char **argv) {
  cout << "isoProbField <mean file> <std file> <w> <h> <d> <iso-value>" << endl;
  if (argc<7)
      return -1;
  string meanfile, stdfile;
  meanfile = argv[1];
  stdfile = argv[2];
  int dim[3];
  dim[0] = atoi(argv[3]);
  dim[1] = atoi(argv[4]);
  dim[2] = atoi(argv[5]);
  cout << "dim: " << dim[0] << "," << dim[1] << "," << dim[2] << endl;
  vector<Gaussian> pDistData = loadData(meanfile, stdfile, dim);
  float iso = atoi(argv[6]);

  vector<float> probField (dim[0]*dim[1]*dim[2]);
  int i,j,k;
  int count=0;
  for (k=0; k<dim[2]; k++)
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++) {
        probField[count] = dist::getPdf(pDistData[count], iso);
        count++;
      }
  FILE *fp = fopen("result.raw", "wb");
  fwrite(&probField[0], count, sizeof(float), fp);
  fclose(fp);
  write_nrrd_3d("result.nrrd", "result.raw", dim[0], dim[1], dim[2], "float");
}
