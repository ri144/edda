#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkStringArray.h>
#include <vtkFieldData.h>
#include "distributions/gaussian_mixture.h"
#include "distributions/histogram.h"

#include "edda_writer.h"

using namespace std;

namespace edda{

	const dist::GMMTuple getGmmModels(dist::Variant &distr, int GMs, int model)
	{
		switch (GMs)
		{
		case 2:
			return boost::get<dist::GaussianMixture<2> >(distr).models[model];
		case 3:
			return boost::get<dist::GaussianMixture<3> >(distr).models[model];
		case 4:
			return boost::get<dist::GaussianMixture<4> >(distr).models[model];
		case 5:
			return boost::get<dist::GaussianMixture<5> >(distr).models[model];
		default:
			throw runtime_error("The Gaussian mixture models exceeds default size");
		}

	}

	void writeGmmArrays(ofstream & myFile, DistrArray *array, const string &array_name, int GMs)
//		vtkPointData *vtk_point_data, DistrArray *array, const string &array_name, int GMs)
	{
		
		//1. number of Gaussian components
		//2. number of array components
		//3. length of array //perhaps check if it fits the dimension size

		myFile.write((char*)(&GMs), sizeof (int));
		int nc = array->getNumComponents();
		myFile.write((char*)(&nc), sizeof (int));
		int n = array->getLength();
		myFile.write((char*)(&n), sizeof (int));


		float* gmData = new float[GMs * 3 * nc*n];
		for (int i = 0; i<GMs * 3; i++){
			for (int j = 0; j < n; j++){
				vector<dist::Variant> vdist = array->getDistrVector(j);
				for (int c = 0; c < nc; c++){
					int index = i*n*nc + j*nc + c;

					if (i % 3 == 0)
						gmData[index] = getGmmModels(vdist[c], GMs, i / 3).m;
					else if (i % 3 == 1)
						gmData[index] = getGmmModels(vdist[c], GMs, i / 3).v;
					else
						gmData[index] = getGmmModels(vdist[c], GMs, i / 3).w;
				}
			}

		}

		myFile.write((char*)(gmData), sizeof (float)*GMs * 3 * nc*n);

		delete gmData;

		/*
		printf("Gaussian Models in GaussianMixture=%d\n", GMs);
		char name[256];
		int i;

		for (i = 0; i<GMs * 3; i++)
		{
			vtkFloatArray *vtk_array = vtkFloatArray::New();
			vtk_array->SetNumberOfComponents(nc);
			vtk_array->SetNumberOfTuples(n);

			if (i % 3 == 0) {
				sprintf(name, "%smean%d", array_name.c_str(), i / 3);
				vtk_array->SetName(name);
				for (int j = 0; j<n; j++)
				{
					vector<dist::Variant> vdist = array->getDistrVector(j);
					for (int c = 0; c<nc; c++)
					{
						((float *)vtk_array->GetVoidPointer(j))[c] = getGmmModels(vdist[c], GMs, i / 3).m;
					}
				}

			}
			else if (i % 3 == 1) {
				sprintf(name, "%svar%d", array_name.c_str(), i / 3);
				vtk_array->SetName(name);

				for (int j = 0; j<n; j++)
				{
					vector<dist::Variant> vdist = array->getDistrVector(j);
					for (int c = 0; c<nc; c++)
					{
						((float *)vtk_array->GetVoidPointer(j))[c] = getGmmModels(vdist[c], GMs, i / 3).v;
					}
				}

			}
			else {
				sprintf(name, "%sweight%d", array_name.c_str(), i / 3);
				vtk_array->SetName(name);

				for (int j = 0; j<n; j++)
				{
					vector<dist::Variant> vdist = array->getDistrVector(j);
					for (int c = 0; c<nc; c++)
					{
						((float *)vtk_array->GetVoidPointer(j))[c] = getGmmModels(vdist[c], GMs, i / 3).w;
					}
				}

			}

			vtk_point_data->AddArray(vtk_array);
		}

		*/
	}


	// edda exported function
	void writeEddaDataset(shared_ptr<Dataset<Real> > dataset, const string &edda_file, const string &array_name_prefix)
	{
		ofstream myFile(edda_file.c_str(), ios::out | ios::binary);

		//1. "EDDA" as a marker
		//2. version number
		//3. gridType. 1: Regular CartesianGrid.
		//if Regular CartesianGrid:
		//	4. size
		//	5. spacing
		//	6. distr type. 1: GaussianMixture

		char eddaFileMark[4] = { 'E', 'D', 'D', 'A' };
		myFile.write(eddaFileMark, sizeof (char)*4);
		char majorVersion = 0, minorVersion = 1;
		myFile.write(&majorVersion, sizeof (char));
		myFile.write(&minorVersion, sizeof (char));

		
		CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(dataset->getGrid());

		if (cartesianGrid) {
			int gridTypeNumber = 1;
			myFile.write((char*)(&gridTypeNumber), sizeof (int));
			int *dims = dataset->getDimension();
			myFile.write((char*)(dims), sizeof (int)* 3);
			float spacing[3];
			dataset->getSpacing(spacing[0], spacing[1], spacing[2]);
			myFile.write((char*)(spacing), sizeof (float)*3);



			
			// write vti

			DistrArray *array = dataset->getArray();
			string dName = array->getDistrName();

			if (dName.compare(0, 15, "GaussianMixture") == 0) {
				// Only compare the first 15 chars because this string ends with the number of Gaussian models
				// Specified in edda::dist::GaussianMixture
				

				int distrTypeNumber = 1;
				myFile.write((char*)(&distrTypeNumber), sizeof (int));

				writeGmmArrays(myFile, array, array_name_prefix, stoi(dName.substr(15)));

				//addVtkGmmArrays(image->GetPointData(), array, array_name_prefix, stoi(dName.substr(15)));

			}
			/*
			
			else if (dName.compare("Histogram") == 0) {
				addVtkHistoArrays(image->GetPointData(), array, array_name_prefix);

			}
			else {
				cout << "Edda VTK Writer: Unsupported array type" << endl;
				throw NotImplementedException();
			}
			setDistrType(image->GetFieldData(), dName, array_name_prefix);

			*/
			printf("Saving converted file to %s.\n", edda_file.c_str());

			
			myFile.close();

		}
		else {

			// TODO for other grid types
			throw NotImplementedException();
		}

		
	}


}