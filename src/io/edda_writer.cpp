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

	void writeGmmArrays(ofstream & myFile, DistrArray *array, int GMs)
	{
		//1. number of array components
		//2. Is # of Gaussian components uniform (1-yes, 0-no)
		int nc = array->getNumComponents();
		myFile.write((char*)(&nc), sizeof (int));
		int isNumGaussianUniform = 1;
		myFile.write((char*)(&isNumGaussianUniform), sizeof (int));

		if (isNumGaussianUniform == 1){
			//3. # of Gaussian components
			myFile.write((char*)(&GMs), sizeof (int));

			//write data
			int n = array->getLength();
			float* gmData = new float[GMs * 3 * nc*n];		
			for (int c = 0; c < nc; c++){
				for (int j = 0; j < n; j++){
					for (int i = 0; i < GMs * 3; i++){
						int index = c*n*GMs * 3 + j*GMs * 3 + i;
						vector<dist::Variant> vdist = array->getDistrVector(j);
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
		}
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

	void writeHistoArrays(ofstream & myFile, DistrArray *array)
	{
		//1. number of array components
		//2. Is # of bins components uniform (1-yes, 0-no)
		int nc = array->getNumComponents();
		myFile.write((char*)(&nc), sizeof (int));
				
		int isNumBinsUniform = 0;
		myFile.write((char*)(&isNumBinsUniform), sizeof (int));

		if (isNumBinsUniform == 1){
			//3. # of Gaussian components
			vector<dist::Variant> vdist = array->getDistrVector(0);
			int nbins = boost::get<dist::Histogram>(vdist[0]).getBins();
			myFile.write((char*)(&nbins), sizeof (int));
		}

		//4. Is min and max value uniform (1-yes, 0-no)
		int isMinMaxValueUniform = 0;
		myFile.write((char*)(&isMinMaxValueUniform), sizeof (int));

		if (isMinMaxValueUniform == 1){
			//5. Min value & Max value
			vector<dist::Variant> vdist = array->getDistrVector(0);
			float minv = boost::get<dist::Histogram>(vdist[0]).getMinValue();
			float maxv = boost::get<dist::Histogram>(vdist[0]).getMaxValue();
			myFile.write((char*)(&minv), sizeof (float));
			myFile.write((char*)(&maxv), sizeof (float));
		}

		if (isNumBinsUniform == 1 && isMinMaxValueUniform == 1){
			//not tested yet
			vector<dist::Variant> vdist = array->getDistrVector(0);
			int nbins = boost::get<dist::Histogram>(vdist[0]).getBins();
			float* tuple = (float*)malloc(sizeof(float)*nbins);

			int n = array->getLength(); //may need to process fake dataset such as n==0?
			for (int j = 0; j < n; j++)
			{
				vector<dist::Variant> vdist = array->getDistrVector(j);
				for (int b = 0; b < nbins; b++)
				{
					tuple[b] = boost::get<dist::Histogram>(vdist[0]).getBinValue(b);
				}
				myFile.write((char*)(tuple), sizeof(float)*nbins);
			}

			free(tuple);
		}
		else if (isNumBinsUniform == 1 && isMinMaxValueUniform != 1){
			throw NotImplementedException();		
		}
		else if (isNumBinsUniform != 1 && isMinMaxValueUniform == 1){
			throw NotImplementedException();
		}
		else{ //(isNumBinsUniform == 0 && isMinMaxValueUniform == 0)

			int n = array->getLength(); //may need to process fake dataset such as n==0?
			
			int* headerBinary_nbins = (int*)malloc(sizeof(int)*n);
			float* headerBinary_minMaxV = (float*)malloc(sizeof(float)*2*n);

			for (int j = 0; j < n; j++)
			{
				vector<dist::Variant> vdist = array->getDistrVector(j);

				dist::Histogram curHist = boost::get<dist::Histogram>(vdist[0]);
				int nbins = curHist.getBins();
				headerBinary_nbins[j] = nbins;
				float minv = curHist.getMinValue();
				float maxv = curHist.getMaxValue();
				headerBinary_minMaxV[2 * j] = minv;
				headerBinary_minMaxV[2 * j + 1] = maxv;
			}
			myFile.write((char*)(headerBinary_nbins), sizeof(float)*n);
			myFile.write((char*)(headerBinary_minMaxV), sizeof(float)*2*n);
			free(headerBinary_nbins);
			free(headerBinary_minMaxV);

			for (int j = 0; j < n; j++)
			{
				vector<dist::Variant> vdist = array->getDistrVector(j);
				dist::Histogram curHist = boost::get<dist::Histogram>(vdist[0]);

				int nbins = curHist.getBins();
				float * tuple = (float*)malloc(sizeof(float)*nbins);

				for (int b = 0; b < nbins; b++)
				{
					tuple[b] = curHist.getBinValue(b);
				}

				myFile.write((char*)(tuple), sizeof(float)*nbins);

				free(tuple);
			}
		}
	}

	// edda exported function

	template <typename T>
	void writeEddaDatasetTemplate(shared_ptr<Dataset<T> > dataset, const string &edda_file)
	//void writeEddaDatasetTemplate(shared_ptr<Dataset<Real> > dataset, const string &edda_file)
	{
		const string array_name_prefix = "";

		ofstream myFile(edda_file.c_str(), ios::out | ios::binary);

		//1. "EDDA" as a marker
		//2. version number
		//3. gridType. 1: Regular CartesianGrid.
		//if Regular CartesianGrid:
		//	4. dimension
		//	5. spacing
		//	6. distr type. 1: GaussianMixture. 2: Histogram
		//	if GaussianMixture:
		//		use the function writeGmmArrays()

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

			DistrArray *array = dataset->getArray();
			string dName = array->getDistrName();

			if (dName.compare(0, 15, "GaussianMixture") == 0) {
				// Only compare the first 15 chars because this string ends with the number of Gaussian models
				// Specified in edda::dist::GaussianMixture			

				int distrTypeNumber = 1;
				myFile.write((char*)(&distrTypeNumber), sizeof (int));

				writeGmmArrays(myFile, array, stoi(dName.substr(15)));
			}
			
			
			else if (dName.compare("Histogram") == 0) {
				int distrTypeNumber = 2;
				myFile.write((char*)(&distrTypeNumber), sizeof (int));

				writeHistoArrays(myFile, array);
			}
			else {
				cout << "Edda Writer: Unsupported array type" << endl;
				throw NotImplementedException();
			}

			printf("Saving converted file to %s.\n", edda_file.c_str());		
			myFile.close();
		}
		else {

			// TODO for other grid types
			throw NotImplementedException();
		}		
	}

	void writeEddaDataset(shared_ptr<Dataset<VECTOR3> > dataset, const string &edda_file)
	{
		writeEddaDatasetTemplate(dataset, edda_file);
	}

	void writeEddaDataset(shared_ptr<Dataset<Real> > dataset, const string &edda_file)
	{
		writeEddaDatasetTemplate(dataset, edda_file);
	}
}