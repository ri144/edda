#include <iostream>
#include <string>
#include <cstdio>
#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/edda_vtk_reader.h"
#include "dataset/dataset.h"

using namespace std;
using namespace edda;

int main(int argc, char **argv)
{
  srand(time(NULL));  // random seeding

  cout << "Loading sample file" << endl;
  //string filename = string(SAMPLE_DATA_PATH) + "/isabel_pressure_small.vti";
  string filename = string(SAMPLE_DATA_PATH) + "/out_92651_0.vts";
  //string filename = "D:/Data/Edda/single_passage_turbine/vts_files_generated/out_0_9602.vts";

  // load data with random sampling
  shared_ptr<Dataset<Real> > dataset = loadEddaScalarDataset(filename, "");

  VECTOR3 pos;
  Real value;
  int i;

  pos = VECTOR3(5, 5, 5);
  value = dataset->at_comp((int)pos[0], (int)pos[1], (int)pos[2]);
  cout << "at_comp " << pos << " : " << value << endl;

  pos = VECTOR3(-0.02, -0.4312, 0.006); //for out_92651_0.vts //vert id 1474 (4,7,3)
  //pos = VECTOR3(-0.03, 0.416, 0.208); //for out_0_9602.vts //vert id 831 (3,8,6)
  //pos = VECTOR3(-0.022, 0.426, 0.186); //another for out_0_9602.vts //vert id 543 (3,4,4)

  if (dataset->at_phys(pos, value) == SUCCESS)
	  cout << pos << ": " << value << endl;
  else
	  cout << pos << ": " << "fail to get value at given position" << endl;

  pos = VECTOR3(0.0, 0.0, 0.0);
  if (dataset->at_phys(pos, value) == SUCCESS)
	  cout << pos << ": " << value << endl;
  else
	  cout << pos << ": " << "fail to get value at given position" << endl;

 


  return 0;
}
