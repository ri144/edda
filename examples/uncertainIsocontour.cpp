// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/file_reader.h"
#include "io/file_writer.h"
#include "filters/uncertain_isocontour.h"

using namespace std;
using namespace edda;

typedef dist::Gaussian<> Gaussian;

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
  float isov = atof(argv[6]);

  shared_ptr<Dataset<Gaussianf> > dataset = loadGaussianRegularGrids(meanfile, stdfile, dim);
  shared_ptr<Dataset<float> > output = uncertainIsocontour(dataset, isov);

  cout << "Output: probField.raw, probField.nrrd" << endl;
  writeRawNrrdFile(output, "probField");
}
