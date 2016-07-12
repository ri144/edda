#include <memory>
#include <vtkDataSet.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkImageData.h>
#include <vtkNew.h>

#include "io/path.h"
#include "dataset/abstract_sampling_array.h"

#include "io/gmm_vtk_data_array.h"
#include "io/edda_reader.h"

using namespace std;

namespace edda {

template <typename T>
shared_ptr<Dataset<T> > loadEddaRandomSamplingDataset(const string &edda_file, const string &array_name)
{
  string ext = getFileExtension(edda_file);
  if (ext.compare("vti")==0) {
     vtkNew<vtkXMLImageDataReader> reader;
     reader->SetFileName(edda_file.c_str());
     reader->Update();
     vtkImageData *vtkdata = reader->GetOutput();

     int *dim = vtkdata->GetDimensions();
     printf("dim: %d %d %d\n", dim[0], dim[1], dim[2]);

     // TODO : check if gmm dataset
     shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
                                   new RegularCartesianGrid(dim[0], dim[1], dim[2]),
                                   new AbstractSamplingArray( new GmmVtkDataArray( vtkdata->GetPointData(), array_name.c_str() ) )
         );
     // or is histogram dataset

     return dataset;

  } else if (ext.compare("vts")==0){ // structured grids

    vtkNew<vtkXMLStructuredGridReader> reader;
    reader->SetFileName(edda_file.c_str());
    reader->Update();
    vtkStructuredGrid *vtkdata = reader->GetOutput();

    int *dim = vtkdata->GetDimensions();
    printf("dim: %d %d %d\n", dim[0], dim[1], dim[2]);

    // TODO : check if gmm dataset
		// TODO: check if is Curvilinear grids
	int numPoints = dim[0] * dim[1] * dim[2];
	CVertex* pVertexGeom = new CVertex[numPoints];
	double tempd[3];
	for (int k = 0; k < dim[2]; k++){
		for (int j = 0; j < dim[1]; j++){
			for (int i = 0; i < dim[0]; i++){
				int ind = k*dim[1] * dim[0] + j*dim[0] + i;// !!! need to confirm !!!
				vtkdata->GetPoint(i, j, k, tempd);
				pVertexGeom[ind].position[0] = tempd[0];
				pVertexGeom[ind].position[1] = tempd[1];
				pVertexGeom[ind].position[2] = tempd[2];

			}
		}
	}

    shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
			new CurvilinearGrid(dim, pVertexGeom),
            new AbstractSamplingArray( new GmmVtkDataArray( vtkdata->GetPointData(), array_name.c_str() ) )
        );
		// or is other structured grids
    // or is histogram dataset

    return dataset;
  } else {
    printf("File format of %s not supported\n", edda_file.c_str());
    exit(1);
  }
}

shared_ptr<Dataset<Real> > loadEddaScalarDataset(const string &edda_file, const string &array_name)
{
  return loadEddaRandomSamplingDataset<Real>(edda_file, array_name);
}

shared_ptr<Dataset<VECTOR3> > loadEddaVector3Dataset(const string &edda_file, const string &array_name)
{
  return loadEddaRandomSamplingDataset<VECTOR3>(edda_file, array_name);
}


} // edda
