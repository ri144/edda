#include <memory>
#include <vtkDataSet.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkFieldData.h>
#include <vtkStringArray.h>
#include <vtkStdString.h>
#include <vtkPoints.h>

#include "io/path.h"

#include "io/gmm_vtk_data_array.h"
#include "io/edda_reader.h"
#include "dataset/distr_array.h"
#include "distributions/histogram.h"

using namespace std;
using namespace edda;
using namespace dist;

namespace edda {



	template <typename T>
	shared_ptr<Dataset<T> > loadEddaDatasetNew(const string &edda_file, const string &array_name_prefix)
	{
		ifstream myfile(edda_file.c_str(), ios::binary);

		streampos begin = myfile.tellg();
		myfile.seekg(0, ios::end);
		streampos end = myfile.tellg();
		int fileByteSize = end - begin;

		myfile.seekg(0, ios::beg);

		char eddaFileMark[4];
		myfile.read(eddaFileMark, sizeof(char)*4);
		if (eddaFileMark[0] != 'E' || eddaFileMark[1] != 'D' || eddaFileMark[2] != 'D' || eddaFileMark[3] != 'A'){
			myfile.close();
			//TODO: throw a proper exception
			exit(0);
		}

		char majorVersion, minorVersion;
		myfile.read(&majorVersion, sizeof(char));
		myfile.read(&minorVersion, sizeof(char));

		if (majorVersion == 0 && minorVersion == 1){
			int gridTypeNumber;
			myfile.read((char*)(&gridTypeNumber), sizeof(int));
			if (gridTypeNumber == 1){
				int dims[3];
				myfile.read((char*)(&dims), sizeof (int)* 3);
				float spacing[3];
				myfile.read((char*)(&spacing), sizeof (float)* 3);

				int distrTypeNumber = 1;
				myfile.read((char*)(&distrTypeNumber), sizeof (int));
				if (distrTypeNumber == 1){
					int GMs;
					myfile.read((char*)(&GMs), sizeof (int));
					int nc;
					myfile.read((char*)(&nc), sizeof (int));
					int n;
					myfile.read((char*)(&n), sizeof (int));
					float* gmData = new float[GMs * 3 * nc*n];
					myfile.read((char*)(gmData), sizeof (float)*GMs * 3 * nc*n);



					if (nc == 1){
						shared_ary<GaussianMixture<5>> s_gmAry(new GaussianMixture<5>[n], n);
						//TODO:: based on GMs may use other than GMM5

						for (int j = 0; j < n; j++){

							std::vector<GMMTuple> models;
							for (int i = 0; i<GMs * 3; i+=3){
								for (int c = 0; c < nc; c++){
									int index = i*n*nc + j*nc + c;
									float p[3] = { gmData[i*n*nc + j*nc + c], gmData[(i + 1)*n*nc + j*nc + c], gmData[(i + 2)*n*nc + j*nc + c] };
									GMMTuple curG;
									curG.m = p[0];
									curG.v = p[1];
									curG.w = p[2];
									curG.p[0] = p[0];
									curG.p[1] = p[1];
									curG.p[2] = p[2];
									models.push_back(curG);
								}
							}

							GaussianMixture<5> curGM(models);
							s_gmAry[j] = models;
						}
						DistrArray * abs_array = new ScalarDistrArray<GaussianMixture<5>>(s_gmAry);

						shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
							new RegularCartesianGrid(dims[0], dims[1], dims[2]),
							abs_array);

						return dataset;

					}
					else{
						delete gmData;
						myfile.close();
						throw NotImplementedException();
					}
					delete gmData;
				}
				else{
					myfile.close();
					throw NotImplementedException();
				}
			}
			else{
				myfile.close();
				throw NotImplementedException();
			}
		}
		else{
			myfile.close();
			throw NotImplementedException();
		}

		myfile.close();
		/*
		string ext = getFileExtension(edda_file);

		if (ext.compare("vti") == 0) {
			vtkNew<vtkXMLImageDataReader> reader;
			reader->SetFileName(edda_file.c_str());
			reader->Update();
			vtkImageData *vtkdata = reader->GetOutput();

			int *dim = vtkdata->GetDimensions();

			// check dataset type
			string type = getDistrType(vtkdata->GetFieldData(), array_name_prefix);

			if (type.compare(0, 15, "GaussianMixture") == 0) {
				shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
					new RegularCartesianGrid(dim[0], dim[1], dim[2]),
					new GmmVtkDataArray(vtkdata->GetPointData(), array_name_prefix.c_str())
					);
				return dataset;


			}
			else if (type.compare("Histogram") == 0)
			{
				shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
					new RegularCartesianGrid(dim[0], dim[1], dim[2]),
					genHistoArray(vtkdata->GetPointData())
					);
				return dataset;
			}
			else {
				cout << "Unknown distribution type: " << type << endl;
				exit(1);
			}


		} // end of if (ext.compare("vti")==0)
		else if (ext.compare("vts") == 0){ // structured grids

			vtkNew<vtkXMLStructuredGridReader> reader;
			reader->SetFileName(edda_file.c_str());
			reader->Update();
			vtkStructuredGrid *vtkdata = reader->GetOutput();

			int *dim = vtkdata->GetDimensions();

			float *point_ary = (float *)vtkdata->GetPoints()->GetVoidPointer(0);

			// check dataset type
			string type = getDistrType(vtkdata->GetFieldData(), array_name_prefix);

			if (type.compare(0, 15, "GaussianMixture") == 0) {
				shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
					new CurvilinearGrid(dim, point_ary),
					new GmmVtkDataArray(vtkdata->GetPointData(), array_name_prefix.c_str())
					);
				return dataset;

			}
			else if (type.compare("Histogram") == 0)
			{
				shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
					new CurvilinearGrid(dim, point_ary),
					genHistoArray(vtkdata->GetPointData())
					);
				return dataset;
			}
			else {
				cout << "Unknown distribution type: " << type << endl;
				exit(1);
			}


		}
		else {
			printf("File format of %s not supported\n", edda_file.c_str());
			exit(1);
		}
		*/
	}



	shared_ptr<Dataset<Real> > loadEddaDataset(const string &edda_file, const string &array_name_prefix)
	{
		return loadEddaDatasetNew<Real>(edda_file, array_name_prefix);
	}
}