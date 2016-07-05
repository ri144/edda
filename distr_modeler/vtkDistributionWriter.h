#ifndef VTK_DISTRIBUTION_WRITER_H_
#define VTK_DISTRIBUTION_WRITER_H_

#include <iostream>
#include <stdlib.h>
#include <string>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>

#include "distributions/distribution.h"
#include "distributions/gaussian.h"

using namespace std;
using namespace edda;
using namespace dist;

class vtkDistibutionWriter{
private:
	string filename;
public:
	vtkDistibutionWriter();
	~vtkDistibutionWriter();
	void write(string fname, int xdim, int ydim, int zdim);
	string getFilename();
};

#endif // VTK_DISTRIBUTION_WRITER_H_
