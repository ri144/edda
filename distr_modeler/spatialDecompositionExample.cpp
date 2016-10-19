#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include "distributionModeler.h"

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

	float *inData;
	inData = new float[xdim*ydim*zdim];	
    
    FILE *fIn = fopen(filename.c_str(),"rb");
    if(!fIn)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }
    size_t readStatus = fread(inData, sizeof(float), xdim*ydim*zdim, fIn);
    fclose(fIn);

    int newW(0),newH(0),newD(0);

    if(xdim % blockXdim == 0)
    	newW = xdim / blockXdim;
    else
    {
    	fprintf(stderr, "Wrong dimension combination\n");
        exit(14);
    }
    if(ydim % blockYdim == 0)
    	newH = ydim / blockYdim;
    else
    {
    	fprintf(stderr, "Wrong dimension combination\n");
        exit(14);
    }
    if(zdim % blockZdim == 0)
    	newD = zdim / blockZdim;
    else
    {
    	fprintf(stderr, "Wrong dimension combination\n");
        exit(14);
    }

    //edda modeling:
	DistributionModeler<dist::GaussianMixture<2>> dm;
	dm.assignGrid(newW, newH, newD);
	//dm.initHistogram(64);

	int counter =0;
    for(size_t z=0; z<zdim; z += blockZdim)
    {
    	for(size_t y=0; y<ydim; y += blockYdim)
    	{
    		for(size_t x=0; x<xdim; x += blockXdim)
    		{

				float *data;
				data = new float[blockXdim*blockYdim*blockZdim];
    			int i = 0;
    			for(size_t zz=z; zz<z+blockZdim; zz++)
    			{
    				for(size_t yy=y; yy<y+blockYdim; yy++)
    				{
    					for(size_t xx=x; xx<x+blockXdim; xx++)
    					{
    						data[i] = inData[zz*xdim*ydim + yy*xdim + xx];
    						i++;
    					}
    				}
    			}
    			dm.computeGMM(data, blockXdim*blockYdim*blockZdim, counter);
    			counter++;
    		}
    	}
    }
    dm.model();
  	dm.writeToVTK("newIsabelGMM.vti","newIsabelGMM_"); 
    
  	return 0;
}