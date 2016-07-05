#include <iostream>
#include <stdlib.h>

#include "distributions/distribution.h"
#include "distributions/gaussian.h"
#include "core/interpolator.h"

//#include "vtkDistributionWriter.h"
//#include "testClass.h"

#include <netcdfcpp.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>

using namespace std;
using namespace edda;
using namespace dist;

static const int NC_ERR = 2;

//shared_ary<Gaussian> loadEnsembleNetCDF(string filename, string xDimName, string yDimName, string zDimName, string ensDimName, string varName);

int main(int argc, char* argv[])
{
	NcError err(NcError::verbose_nonfatal);

	int xDim(0), yDim(0), zDim(0), ensDim(0);

  if(argc > 7){
    cerr << "Too few arguments provide. Check Usage." << endl;
    exit(13);
  }

  string fileName = argv[1];
  
	/*string fileName = "/media/WD3/oceanEnsemble/pe_dif_sep2_98.nc";
	string xDimName = "tlat";
	string yDimName = "tlon";
	string zDimName = "outlev";
	string ensDimName = "time";
	string varName = "temp";*/

  NcFile dataFile(fileName.c_str(), NcFile::ReadOnly);
  if(!dataFile.is_valid())
  {
    cerr << "ERROR[" << NC_ERR << "]" << "reading file" << endl;
    exit(NC_ERR);
  }
  
  string xDimName = argv[2];
  string yDimName = argv[3];
  string zDimName = argv[4];
  string ensDimName = argv[5];
  string varName = argv[6];

  xDim = dataFile.get_dim(xDimName.c_str())->size();
  yDim = dataFile.get_dim(yDimName.c_str())->size();
  zDim = dataFile.get_dim(zDimName.c_str())->size();
  ensDim = dataFile.get_dim(ensDimName.c_str())->size();
  //ensDim = 10;
  NcVar *var = dataFile.get_var(varName.c_str());

  cout << "Dimensions : (" << xDim << "," << yDim << "," << zDim << ")" << endl;
  cout << "Ensemble Count : " << ensDim << endl;
  	
  float *data;

  data = new float[ensDim];
  // Create Gaussian array
  Gaussian *pData = new Gaussian[xDim*yDim*zDim];

 	
 	for(int z=0; z<zDim; z++)
 	{
 		for(int y=0; y<yDim; y++)
 		{
 			for(int x=0; x<xDim; x++)
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
        //TODO::Will call the GMM API from here with data[] as the input array of points.
 				float sum = 0;
 				for(int i=0; i<ensDim; i++)
 					sum += data[i];
 				sum = sum / ensDim;
 				float var = 0;
 				for(int i=0; i<ensDim; i++)
 					var += (sum - data[i])*(sum - data[i]);
 				var = var / ensDim;
 				
        pData[z*xDim*yDim + y*xDim + x] = Gaussian(sum, var);
 			}
 		}
 	}

  cout << "DISTRIBUTION TEST:" << endl;
  cout << pData[750] << endl;
  cout << dist::getMean(pData[750]) << ":::" << dist::getVar(pData[750]) << endl;

 	//vtkDistibutionWriter distWriter;
  //distWriter.write("mean1.vti", 10, 12, 12);

 	string outfilename = "meanData1.vti";
  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
	imageData->SetDimensions(xDim,yDim,zDim);
 	vtkSmartPointer<vtkDoubleArray> mu =  vtkSmartPointer<vtkDoubleArray>::New();
	mu->SetNumberOfComponents(1);
  mu->SetNumberOfTuples(xDim*yDim*zDim);

  vtkSmartPointer<vtkDoubleArray> sd =  vtkSmartPointer<vtkDoubleArray>::New();
 	sd->SetNumberOfComponents(1);
  sd->SetNumberOfTuples(xDim*yDim*zDim);


 	int* dims = imageData->GetDimensions();

  for(int i=0; i<dims[0]*dims[1]*dims[2]; i++)
  {
    float mu1 = dist::getMean(pData[i]);
    float sd1 = sqrt(dist::getVar(pData[i]));
  	mu->SetTuple(i, &mu1);
   	sd->SetTuple(i, &sd1);
  }

  imageData->GetPointData()->AddArray(mu);
  mu->SetName("mean0");

  imageData->GetPointData()->AddArray(sd);
  sd->SetName("stdev0");

  vtkSmartPointer<vtkXMLImageDataWriter> writer =  vtkSmartPointer<vtkXMLImageDataWriter>::New();
 	writer->SetFileName(outfilename.c_str());
#if VTK_MAJOR_VERSION <= 5
 		writer->SetInputConnection(imageData->GetProducerPort());
#else
 		writer->SetInputData(imageData);
#endif
 	writer->Write();
  cout << "Distribution File Created!" << endl;
  /*cout << "num of dimensions = " << dataFile.num_dims() << endl;
    cout << "num of variables = " << dataFile.num_vars() << endl;
    cout << "num of attributes = " << dataFile.num_atts() << endl;
    
    NcDim *dim = dataFile.get_dim("time");
    cout << "dim size = " << dim->size() << endl;*/
    
  return 0;
  	
}

