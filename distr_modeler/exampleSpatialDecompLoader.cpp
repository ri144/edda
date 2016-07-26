
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
	if(argc < 8)
	{
		cout << "USAGE: check the parameter list!!\n";
		exit(11);
	}
	string filename = argv[1];
	int xdim = atoi(argv[2]);
	int ydim = atoi(argv[3]);
	int zdim = atoi(argv[4]);
	int blockXdim = atoi(argv[5]);
	int blockYdim = atoi(argv[6]);
	int blockZdim = atoi(argv[7]);

  	DistrModeler dm1(GMM2); 
  	dm1.loader(filename, xdim, ydim, zdim, blockXdim, blockYdim, blockZdim);
  	dm1.writeToVTK("isabelGmm.vti","isabelGmm_");

  	DistrModeler dm2(HIST); 
  	dm2.initHistogram(64);
  	dm2.loader(filename, xdim, ydim, zdim, blockXdim, blockYdim, blockZdim);
  	dm2.writeToVTK("isabelHist.vti","isabelHist_");

  
    
  	return 0;
}