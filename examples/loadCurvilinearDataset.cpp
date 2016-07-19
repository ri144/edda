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

  VECTOR3 pos, pos2, pos3;
  pos = VECTOR3(5, 5, 5);
  pos2 = VECTOR3(2.1, 2.1, 2.1);
  pos3 = VECTOR3(0.0, 0.0, 0.0);

  cout << "Loading sample file" << endl;

  ////for testing regular grid
  //string filename = string(SAMPLE_DATA_PATH) + "/isabel_pressure_small.vti";
  
  //for testing curvilinear grid
  string filename = string(SAMPLE_DATA_PATH) + "/out_92651_0.vts";
  pos2 = VECTOR3(-0.02, -0.4312, 0.006); //for out_92651_0.vts //vert id 1474 (4,7,3)
  //however, error occurs when testing out_0_9602.vts with pos2 = VECTOR3(-0.03, 0.416, 0.208); //for out_0_9602.vts //vert id 831 (3,8,6)

  // load data with random sampling
  shared_ptr<Dataset<Real> > dataset = loadEddaScalarDataset(filename, "");

  
  Real value;
  int i,j,k;
  i = (int)pos[0];
  j = (int)pos[1];
  k = (int)pos[2];
  
  value = dataset->at_comp(i,j,k);
  cout << "at_comp " << i << "," << j << "," << k << " : " << value << endl;

  if (dataset->at_phys(pos2, value) == SUCCESS)
	  cout << pos2 << ": " << value << endl;
  else
	  cout << pos2 << ": " << "fail to get value at given position" << endl;

  
  if (dataset->at_phys(pos3, value) == SUCCESS)
	  cout << pos3 << ": " << value << endl;
  else
	  cout << pos3 << ": " << "fail to get value at given position" << endl;

  //dist::Variant distr = dataset->at_comp_distr(i,j,k);
  //cout << "at_comp_distr " << i << "," << j << "," << k << " : " << distr << endl;


 
  return 0;
}
