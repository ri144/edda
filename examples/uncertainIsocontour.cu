// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "edda.h"
#include "io/file_reader.h"
#include "io/file_writer.h"

using namespace std;
using namespace edda;
using namespace edda::dist;

int main(int argc, char **argv) {
  cout << "isoProbField <info file> <iso-value>" << endl;
  if (argc<=2)
    return -1;
  string info_file = argv[1];
  float isov = atof(argv[2]);

  // load data
  shared_ptr<Dataset<Gaussian> > dataset = loadData<Gaussian>(info_file);

  // compute isocontour
  shared_ptr<Dataset<Real> > output = uncertainIsocontour(dataset, isov);

  // save into file
  cout << "Output: probField.raw, probField.nrrd" << endl;
  writeRawNrrdFile(output, "probField");
}
