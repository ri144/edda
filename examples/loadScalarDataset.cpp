#include <iostream>
#include <string>
#include <cstdio>
#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/edda_reader.h"
#include "dataset/dataset.h"
#include "io/edda_reader.h"

using namespace std;
using namespace edda;

int main(int argc, char **argv)
{
  srand(time(NULL));  // random seeding

  cout << "Loading sample file" << endl;
  //string filename = string(SAMPLE_DATA_PATH) + "/isabel_pressure_small.vti";
  //string filename = string(SAMPLE_DATA_PATH) + "/out_92651_0.vts";
  string filename = "D:/Data/Edda/single_passage_turbine/vts_files_generated/out_0_9602.vts";

  // load data with random sampling
  shared_ptr<Dataset<Real> > dataset = loadEddaScalarDataset(filename, "");

  VECTOR3 pos;
  Real value;
  int i;

  //pos = VECTOR3(10,10,10);
  pos = VECTOR3(0.0,0.01,0.01);
  dataset->at_phys(pos, value);
  cout << pos << ": " << value << endl;

  pos = VECTOR3(2.1,2.1,2.1);
  dataset->at_phys(pos, value);
  cout << pos << ": " << value << endl;

  value = dataset->at_comp(5, 5, 5);
  cout << "at_comp(5,5,5) : " << value << endl;

  return 0;
}
