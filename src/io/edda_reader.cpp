#include <memory>
#include <vector>
#include <stdexcept>
#include <string>
#include <fstream>

#include "edda_reader.h"
#include "dataset/distr_array.h"
#include "dataset/dataset.h"
#include "core/interpolator.h"
#include "../test/bmp_image.h"

#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "distributions/variant.h"

#include "distributions/gmm.h"
#include "distributions/gaussian_mixture.h""
#include "distributions/joint_gaussian.h"
#include "distributions/joint_GMM.h"
#include "distributions/histogram.h"
#include "distributions/joint_histogram.h"

#include "distributions/distribution_modeler.h"


using namespace std;
using namespace edda;
using namespace dist;

namespace edda {

	DistrArray * readMixArray(ifstream & myfile, int n)
	{
		//only used to test jointGMM
		bool useModeler = false;
		DistributionModeler mv_dm(n);
		int blockSize = 20;
		string filename = string(SAMPLE_DATA_PATH) + "/test_img.bmp";
		BMPImage image(filename.c_str());
		//number of row and col after down-sampleing
		int dsVs = image.height / blockSize;
		int dsUs = image.width / blockSize;
		dsUs = 2, dsVs = 2;
		//Joint GMM array
		shared_ary<JointGMM> array(new JointGMM[dsVs*dsUs], dsVs*dsUs);
		thrust::default_random_engine rng;//random engine for getJointSample()
		Real* varR = (Real*)malloc(sizeof(Real)*blockSize*blockSize);
		Real* varG = (Real*)malloc(sizeof(Real)*blockSize*blockSize);
		Real* varB = (Real*)malloc(sizeof(Real)*blockSize*blockSize);



		dist::Variant* distArray;
		distArray = new dist::Variant[n];

		float* distData = new float[333 * 3]; // !!! this array is designed to avoid the time cost by multiple memory allocation and deletion. better check if the size is big enough before each time using it !!!

		for (int index = 0; index < n; index++){
			int distrTypeNumber = 1;
			myfile.read((char*)(&distrTypeNumber), sizeof (int));
			if (distrTypeNumber == 1){ //gaussian

				int nGM;
				myfile.read((char*)(&nGM), sizeof (int));

				myfile.read((char*)(distData), sizeof (float)*nGM * 3);

				dist::GMM new_gmm = dist::GMM(nGM);
				for (int m = 0; m<nGM; m++)
				{
					new_gmm.models[m].m = distData[3 * m];
					new_gmm.models[m].v = distData[3 * m + 1];
					new_gmm.models[m].w = distData[3 * m + 2];
				}
				distArray[index] = new_gmm;
			}
			else if (distrTypeNumber == 2){ //histogram
				int nbins;
				myfile.read((char*)(&nbins), sizeof (int));

				myfile.read((char*)(distData), sizeof (float));
				myfile.read((char*)(distData + 1), sizeof (float));
				myfile.read((char*)(distData + 2), sizeof (float)*nbins);

				distArray[index] = Histogram(distData, nbins);
			}
			else if (distrTypeNumber == 3){ //joint gmm
				int nVar, nComp;
				myfile.read((char*)(&nVar), sizeof (int));
				myfile.read((char*)(&nComp), sizeof (int));

				int sizeEachComp = 1 + nVar + (nVar + 1)*nVar / 2;
				myfile.read((char*)(distData), sizeof (float)* sizeEachComp * nComp);
				
				/*
				// original way. should not need to use distributionModeler
				std::vector<Real> weights(nComp);
				std::vector<JointGaussian> gaus(nComp);
				for (int c = 0; c < nComp; c++)
				{
					weights[c] = distData[sizeEachComp * c];

					ublas_vector mean = ublas_vector(nVar, 0);
					ublas_matrix cov = ublas::zero_matrix<Real>(nVar, nVar);;
					for (int jj = 0; jj < nVar; jj++){
						mean(jj) = distData[sizeEachComp * c + 1 + jj];
					}
					int count = 0;
					for (int jj = 0; jj < nVar; jj++){
						for (int i = jj; i < nVar; i++){
							cov(jj, i) = distData[c*sizeEachComp + 1 + nVar + count];//cov row major
							count++;
						}
					}
					for (int jj = 0; jj < nVar; jj++){
						for (int i = 0; i < jj; i++){
							cov(jj, i) = cov(i, jj);
						}
					}
					gaus[c] = JointGaussian(mean, cov);
				}
				dist::JointGMM new_distr(weights, gaus, nVar, nComp);
				distArray[index] = new_distr;
				*/

				
				//test methods
				useModeler = true;
				
				ublas_vector w = ublas_vector(nComp, 0);
				ublas_matrix m(nComp, nVar);
				ublas_matrix covs(nComp* nVar, nVar);
				for (int c = 0; c < nComp; c++)
				{
					w[c] = distData[sizeEachComp * c];
					ublas_vector mean = ublas_vector(nVar, 0);
					ublas_matrix cov = ublas::zero_matrix<Real>(nVar, nVar);;
					for (int jj = 0; jj < nVar; jj++){
						mean(jj) = distData[sizeEachComp * c + 1 + jj];
					}
					int count = 0;
					for (int jj = 0; jj < nVar; jj++){
						for (int i = jj; i < nVar; i++){
							cov(jj, i) = distData[c*sizeEachComp + 1 + nVar + count];//cov row major
							count++;
						}
					}
					for (int jj = 0; jj < nVar; jj++){
						for (int i = 0; i < jj; i++){
							cov(jj, i) = cov(i, jj);
						}
					}
					//copy single GM to full GMM
					for (int jj = 0; jj < nVar; jj++){
						m(c, jj) = mean(jj);
					}
					subrange(covs, c*nVar, (c + 1)*nVar, 0, nVar) = cov;
				}
				
				//test methods 1
				//dist::JointGMM gmm;
				//gmm.setGMM(nVar, nComp, w, m, covs);
				//mv_dm.setJointGMM(gmm, index);

				//test methods 2
				mv_dm.setJointGMM(nVar, nComp, w, m, covs, index);

				/*
				//test methods 3
				int dsV = index / dsUs, dsU = index % dsUs;

				printf("processing index %d %d\n", dsV, dsU);
				//prepare training vectors to OpenCV EM
				int cnt = 0;
				for (int v = dsV*blockSize; v<(dsV + 1)*blockSize; v++) {//row
					for (int u = dsU*blockSize; u<(dsU + 1)*blockSize; u++) {//col
						varR[cnt] = (Real)(image.bitmapImage[(v*image.width + u) * 3 + 0]);
						varG[cnt] = (Real)(image.bitmapImage[(v*image.width + u) * 3 + 1]);
						varB[cnt] = (Real)(image.bitmapImage[(v*image.width + u) * 3 + 2]);
						cnt++;
					}
				}
				//EM in Edda
				std::vector<Real*> trainSamples;
				trainSamples.push_back(varR);
				trainSamples.push_back(varG);
				trainSamples.push_back(varB);
				mv_dm.computeJointGMM(trainSamples, blockSize*blockSize * 1, 2, index);
				*/
			}
			else{
				myfile.close();
				throw NotImplementedException();
			}
		}

		delete[] distData;

		if (useModeler){
			return mv_dm.getDistrArray();
		}
		else{
			shared_ary<dist::Variant> pArray(new dist::Variant[n], n);
			for (int i = 0; i < n; i++)
			{
				pArray[i] = distArray[i];
			}
			DistrArray * abs_array = new ScalarDistrArray<dist::Variant>(pArray);
			return abs_array;
		}
	}


	template <typename T>
	shared_ptr<Dataset<T> > loadEddaDatasetTemplate(const string &edda_file)
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
			/*
			//not supported any more
			*/
		}
		else if (majorVersion == 0 && minorVersion == 2){
			int gridTypeNumber;
			myfile.read((char*)(&gridTypeNumber), sizeof(int));
			if (gridTypeNumber == 1){
				int dims[3];
				myfile.read((char*)(&dims), sizeof (int)* 3);
				float spacing[3];
				myfile.read((char*)(&spacing), sizeof (float)* 3);

				
				int numDistrArray = 1;
				myfile.read((char*)(&numDistrArray), sizeof (int));
				
				std::vector<DistrArray *> dVec(numDistrArray);
				for (int i = 0; i < numDistrArray; i++){
					dVec[i] = readMixArray(myfile, dims[0] * dims[1] * dims[2]);
				}
				
				myfile.close();

				return make_shared<Dataset<T>>(new RegularCartesianGrid(dims[0], dims[1], dims[2]), dVec);
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

		printf("Read file from %s.\n", edda_file.c_str());
	}



	shared_ptr<Dataset<Real> > loadEddaScalarDataset_noneVTK(const string &edda_file)
	{
		return loadEddaDatasetTemplate<Real>(edda_file);
	}
}