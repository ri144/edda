#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include <distributions/gaussian_mixture.h>
#include "distributions/histogram.h"
#include <distributions/estimate_gmm.h>
#include "distributions/variant.h"
#include "dataset/dataset.h"
#include "dataset/distr_array.h"

#include <sys/time.h>


#include "VecDistrModeler.h"

using namespace std;
using namespace edda;
typedef unsigned long long timestamp_t;

static timestamp_t
get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}


int main(int argc, char* argv[])
{
	if(argc < 8)
	{
		cout << "USAGE: check the parameter list!!\n";
		exit(11);
	}
	string filename1 = argv[1];
    string filename2 = argv[2];
    string filename3 = argv[3];
	int xdim = atoi(argv[4]);
	int ydim = atoi(argv[5]);
	int zdim = atoi(argv[6]);
	int blockXdim = atoi(argv[7]);
	int blockYdim = atoi(argv[8]);
	int blockZdim = atoi(argv[9]);

	float *inData1;
	inData1 = new float[xdim*ydim*zdim];	
    
    FILE *fIn1 = fopen(filename1.c_str(),"rb");
    if(!fIn1)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }
    size_t readStatus1 = fread(inData1, sizeof(float), xdim*ydim*zdim, fIn1);
    fclose(fIn1);

    float *inData2;
    inData2 = new float[xdim*ydim*zdim];    
    
    FILE *fIn2 = fopen(filename2.c_str(),"rb");
    if(!fIn2)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }
    size_t readStatus2 = fread(inData2, sizeof(float), xdim*ydim*zdim, fIn2);
    fclose(fIn2);

    float *inData3;
    inData3 = new float[xdim*ydim*zdim];    
    
    FILE *fIn3 = fopen(filename3.c_str(),"rb");
    if(!fIn3)
    {
        fprintf(stderr, "Error opening file\n");
        exit(13);
    }
    size_t readStatus3 = fread(inData3, sizeof(float), xdim*ydim*zdim, fIn3);
    fclose(fIn3);

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

    //edda ensemble data modeling
	VecDistrModeler<dist::Histogram> dm;   //create a distr modeler class for Histogram.
	dm.initVectorComponent(3);
    dm.assignGrid(newW, newH, newD);    //define the final distribution grid dimensions.
    dm.initHistogram(64);       //initialize histogram bin count.

    cout << "starting partitioning\n";
	int counter =0;
    timestamp_t t0 = get_timestamp();

    for(size_t z=0; z<zdim; z += blockZdim)
    {
    	for(size_t y=0; y<ydim; y += blockYdim)
    	{
    		for(size_t x=0; x<xdim; x += blockXdim)
    		{

				float *data;
				data = new float[blockXdim*blockYdim*blockZdim*3];
    			int i = 0;
    			for(size_t zz=z; zz<z+blockZdim; zz++)
    			{
    				for(size_t yy=y; yy<y+blockYdim; yy++)
    				{
    					for(size_t xx=x; xx<x+blockXdim; xx++)
    					{
    						data[i] = inData1[zz*xdim*ydim + yy*xdim + xx];
    						i++;
                            data[i] = inData2[zz*xdim*ydim + yy*xdim + xx];
                            i++;
                            data[i] = inData3[zz*xdim*ydim + yy*xdim + xx];
                            i++;
    					}
    				}
    			}
    			//dm.computeGMM(data, blockXdim*blockYdim*blockZdim, counter);
                dm.computeHistogram(data, blockXdim*blockYdim*blockZdim, counter);
    			counter++;
                cout << "completed " << counter-1 << "[" << z << "][" << y << "][" << x << "]" << "block\n";
    		}
    	}
    }
    dm.model();    //combine the distribution array with the assigned grid.
  	//dm.writeToVTK("newIsabelHist.vti","newIsabelHist_");    //store in vtk format. 

    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;

    cout << "time:" << secs << endl;

    shared_ptr<Dataset<VECTOR3>> eddaDS = dm.getDataset();

    DistrArray * array1 = eddaDS->getArray();
    int *dim = eddaDS->getDimension();

    cout << "DIM:" << *dim << "," << *(dim+1) << "," << *(dim+2) << endl;

    
  	return 0;
}