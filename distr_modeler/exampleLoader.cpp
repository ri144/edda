#include <iostream>
#include <stdlib.h>

#include <cmath>
#include <vector>
#include <sstream>

#include "distr_modeler.h"

using namespace std;
using namespace edda;


int main(int argc, char* argv[])
{
  DistrModeler dm(GMM);
  cout << dm.getType() << endl;

  //loading a netCDF ensemble file with dimension and variable names.
  dm.loader("/media/WD3/oceanEnsemble/pe_dif_sep2_98.nc", "tlat", "tlon", "outlev", "time", "temp");

  //cout << "distribution at location {45,26,10}: " <<  dm.getDistributionAt(45,26,10) << endl;
  dm.writeToVTK("tempEns.vti","tempEns_");//TODO: write to .vti files.

    
  return 0;
}

