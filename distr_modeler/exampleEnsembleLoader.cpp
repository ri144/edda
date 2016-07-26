
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
	if(argc < 7)
	{
		cout << "USAGE: check the parameter list!!\n";
		exit(11);
	}
	string filename = argv[1];
	string xdimName = argv[2];
	string ydimName = argv[3];
	string zdimName = argv[4];
	string ensdimName = argv[5];
	string varName = argv[6];

  	DistrModeler dm1(GMM);
   	//loading a netCDF ensemble file with dimension and variable names.
  	dm1.loader(filename, xdimName, ydimName, zdimName, ensdimName, varName);
  	dm1.writeToVTK("tempEnsGmm.vti","tempEnsGmm_");

  	DistrModeler dm2(HIST);
  	dm2.initHistogram(64);
   	//loading a netCDF ensemble file with dimension and variable names.
  	dm2.loader(filename, xdimName, ydimName, zdimName, ensdimName, varName);
  	dm2.writeToVTK("tempEnsHist.vti","tempEnsHist_");

  
    
  	return 0;
}