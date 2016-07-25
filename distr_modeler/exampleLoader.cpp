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
  DistrModeler dm1(GMM);
  
  //loading a netCDF ensemble file with dimension and variable names.

  dm1.loader("/media/WD3/oceanEnsemble/pe_dif_sep2_98.nc", "tlat", "tlon", "outlev", "time", "temp");
  dm1.writeToVTK("tempEnsGmm.vti","tempEnsGmm_");


  DistrModeler dm2(HIST);
  dm2.initHistogram(64);
  dm2.loader("/media/WD3/oceanEnsemble/pe_dif_sep2_98.nc", "tlat", "tlon", "outlev", "time", "temp");
  dm2.writeToVTK("tempEnsHist.vti","tempEnsHist_");

    
  return 0;
}

