#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkStringArray.h>
#include <vtkFieldData.h>
#include "distributions/gaussian_mixture.h"
#include "distributions/histogram.h"

#include "edda_vtk_writer.h"

using namespace std;

namespace edda{

// add a GMM array for vtkPointData
void addVtkGmmArrays(vtkPointData *vtk_point_data, AbstractDistrArray *array, const string &array_name)
{
  vector<string> filenames;
  char name[256];
  int i;
  int n = array->getLength();
  int nc = array->getNumComponents();

  printf("%d %d\n", n, nc);
  getchar();

  for (i=0; i<MAX_GMs*3; i++)
  {
    vtkFloatArray *vtk_array = vtkFloatArray::New();
    vtk_array->SetNumberOfComponents( nc );
    vtk_array->SetNumberOfTuples( n );

    if (i%3==0) {
      sprintf(name, "%smean%d", array_name.c_str(), i/3);
      vtk_array->SetName(name);
      for (int j=0; j<n; j++)
      {
        vector<dist::Variant> vdist = array->getDistrVector(j);
        for (int c=0; c<nc; c++)
        {
			((float *)vtk_array->GetVoidPointer(j))[c] = boost::get<dist::DefaultGaussianMixture>(vdist[c]).models[i / 3].m;
        }
      }

    } else if (i%3==1) {
      sprintf(name, "%svar%d", array_name.c_str(), i/3);
      vtk_array->SetName(name);

      for (int j=0; j<n; j++)
      {
        vector<dist::Variant> vdist = array->getDistrVector(j);
        for (int c=0; c<nc; c++)
        {
			((float *)vtk_array->GetVoidPointer(j))[c] = boost::get<dist::DefaultGaussianMixture>(vdist[c]).models[i / 3].v;
        }
      }

    } else {
      sprintf(name, "%sweight%d", array_name.c_str(), i/3);
      vtk_array->SetName(name);

      for (int j=0; j<n; j++)
      {
        vector<dist::Variant> vdist = array->getDistrVector(j);
        for (int c=0; c<nc; c++)
        {
          ((float *)vtk_array->GetVoidPointer(j))[c] = boost::get<dist::DefaultGaussianMixture>( vdist[c] ).models[i/3].w;
        }
      }

    }

    vtk_point_data->AddArray(vtk_array);
  }
}

void addVtkHistoArrays(vtkPointData *vtk_point_data, AbstractDistrArray *array, const string &array_name)
{
	vector<string> filenames;
	char name[256];
	vector<dist::Variant> vdist = array->getDistrVector(0);
	int bins = boost::get<dist::Histogram>(vdist[0]).getBins();
	int n = array->getLength();
	int nc = bins + 3; //3 is number of bins, value_min, value_max

	vtkFloatArray *vtk_array = vtkFloatArray::New();
	vtk_array->SetNumberOfComponents(nc);
	vtk_array->SetNumberOfTuples(n);

	sprintf(name, "Variable_0" );
	vtk_array->SetName(name);

	float * tuple = (float*)malloc(sizeof(float)*nc);
	for (int j = 0; j < n; j++)
	{
		vector<dist::Variant> vdist = array->getDistrVector(j);

		tuple[0] = boost::get<dist::Histogram>(vdist[0]).getBins();
		tuple[1] = boost::get<dist::Histogram>(vdist[0]).getMinValue();
		tuple[2] = boost::get<dist::Histogram>(vdist[0]).getMaxValue();
		for (int b = 0; b < bins; b++)
		{
			tuple[b+3] = boost::get<dist::Histogram>(vdist[0]).getBinValue(b);
		}

		vtk_array->SetTuple(j, tuple);
	}
	free(tuple);
	vtk_point_data->AddArray(vtk_array);
}

void setDistrType(vtkFieldData* vtk_data, const char *distrName, const string &array_name_prefix)
{
  vtkStringArray *vtk_array = vtkStringArray::New();
  vtk_array->InsertNextValue(distrName);
  string array_name = array_name_prefix + "distr_type";
  vtk_array->SetName( array_name.c_str() );
  vtk_data->AddArray(vtk_array);
}

// edda exported function
void writeEddaVtkDataset(shared_ptr<Dataset<Real> > dataset, const string &edda_file, const string &array_name_prefix)
{
  CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(dataset->getGrid() ) ;

  if (cartesianGrid) {
    // write vti
    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
    int *dims = dataset->getDimension();
    image->SetDimensions(dims);
    image->SetExtent(0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1);

    AbstractDistrArray *array = dataset->getArray();
    const char *dName = array->getDistrName();
    if (strcmp(dName, "GaussianMixture") == 0) {
      addVtkGmmArrays(image->GetPointData(), array, array_name_prefix);

    } else if (strcmp(dName, "Histogram") == 0) {

      addVtkHistoArrays(image->GetPointData(), array, array_name_prefix);
    } else {
      cout << "Edda VTK Writer: Unsupported array type" << endl;
      throw NotImplementedException();
    }
    setDistrType(image->GetFieldData(), dName, array_name_prefix);

    printf("Saving converted file to %s.\n", edda_file.c_str());
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(edda_file.c_str());
    writer->SetInputData(image);
    writer->SetEncodeAppendedData(0); // no base-64 encoding
    writer->SetCompressorTypeToZLib();
    writer->Write();

  } else {

    // TODO for other grid types
    throw NotImplementedException();
  }
}

}
