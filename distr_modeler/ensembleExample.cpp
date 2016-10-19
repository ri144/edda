#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include <netcdfcpp.h>

#include "distributionModeler.h"

using namespace std;
using namespace edda;

const int NC_ERR = 2;

int main(int argc, char* argv[])
{
  if(argc < 7)
  {
    cout << "USAGE: check the parameter list!!\n";
    exit(11);
  }
  string filename = argv[1];
  string xDimName = argv[2];
  string yDimName = argv[3];
  string zDimName = argv[4];
  string ensDimName = argv[5];
  string varName = argv[6];

  //try loading the netcdf file.
  NcFile dataFile(filename.c_str(), NcFile::ReadOnly);
  if(!dataFile.is_valid())
  {
    cerr << "ERROR[" << NC_ERR << "]" << "reading file" << endl;
    exit(NC_ERR);
  }

  //update the data dimensions by checking the nc file.
  int xDim = dataFile.get_dim(xDimName.c_str())->size();
  int yDim = dataFile.get_dim(yDimName.c_str())->size();
  int zDim = dataFile.get_dim(zDimName.c_str())->size();
  int ensDim = dataFile.get_dim(ensDimName.c_str())->size();
  ensDim = 10;

  NcVar *var = dataFile.get_var(varName.c_str());

  cout << "Dimensions : (" << xDim << "," << yDim << "," << zDim << ")" << endl;
  cout << "Ensemble Count : " << ensDim << endl;

  float *data;
  data = new float[ensDim];

  //edda modeling:
  DistributionModeler<dist::GaussianMixture<2>> dm;
  dm.assignGrid(xDim, yDim, zDim);
  //dm.initHistogram(64);

    

  for(size_t z=0; z<zDim; z++)
  {
    for(size_t y=0; y<yDim; y++)
    {
      for(size_t x=0; x<xDim; x++)
      {
        if(!var->set_cur(0, x, y, z))
        {
          cerr << "ERROR[" << NC_ERR << "]" << "setting current pointer" << endl;
          exit(NC_ERR); 
        }
        if(!var->get(data, ensDim, 1, 1, 1))
        {
          cerr << "ERROR[" << NC_ERR << "]" << "getting data" << endl;
          exit(NC_ERR);
        }

        dm.computeGMM(data, ensDim, z*xDim*yDim + y*xDim + x);       
      }
    }
  }
  dm.model();
  dm.writeToVTK("newTempEnsGMM.vti","newTempEnsGMM_");

  
    
    return 0;
}
