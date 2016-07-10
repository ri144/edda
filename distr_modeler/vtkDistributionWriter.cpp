#include "vtkDistributionWriter.h"

using namespace std;
using namespace edda;
using namespace dist;

vtkDistibutionWriter::vtkDistibutionWriter(){
	//filename = fname;
	//cout << distData[749] << endl;
}
void vtkDistibutionWriter::write(string fname, int xdim, int ydim, int zdim){
	filename = fname;
	cout << filename << endl;
}

string vtkDistibutionWriter::getFilename(){
	return filename;
}