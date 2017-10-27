#include "distributions/joint_GMM.h"
#include "distributions/joint_gaussian.h"
#include "vector"
#include "../test/bmp_image.h"
#include "dataset/distr_array.h"
#include "io/edda_vtk_reader.h"
#include "io/edda_vtk_writer.h"
#include "io/edda_reader.h"
#include "io/edda_writer.h"
#include "distributions/histogram.h"
#include "dataset/distr_array.h"

#include <string>

#include "distributions/distribution_modeler.h"

//this file uses components that are currently under construction, so it is temporarily turned off. Will be back soon.


using namespace edda;
using namespace std;
using namespace edda::dist;

//using image to test joint Gaussian GMM
int main()
{
	//load testing image from disk

	string filename = string(SAMPLE_DATA_PATH) + "/test_img.bmp";

	BMPImage image(filename.c_str());

	//down-sampled block size
	int blockSize = 20;
	int nVar = 3;//number of variables, rgb is 3 variables
	int nGmmComp = 2;//number of Gaussian compoenents


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


	DistributionModeler mv_dm(dsVs*dsUs * 1);
	int counter = 0;


	//loop: go through each block
	for (int dsV = 0; dsV < dsVs; dsV++){
		for (int dsU = 0; dsU < dsUs; dsU++){
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
			//array[dsV*dsUs + dsU] = eddaComputeJointGMM(trainSamples, blockSize*blockSize, nGmmComp);

			//JointGMM test = eddaComputeJointGMM(trainSamples, blockSize*blockSize, nGmmComp);
			//array[dsV*dsUs + dsU] = test;

			//int s = 9;

			mv_dm.computeJointGMM(trainSamples, blockSize*blockSize * 1, nGmmComp, counter);
			counter++;
		}
	}



	std::vector<DistrArray *> dVec;
	dVec.push_back(mv_dm.getDistrArray());
	Dataset<Real> *ds = new Dataset<Real>(new RegularCartesianGrid(dsUs, dsVs, 1), dVec);
	shared_ptr<Dataset<Real>> shr_ds(ds);

	DistrArray *arr = shr_ds->getArray(0);
	for (int j = 0; j < 4; j++){
		dist::Variant curDist = arr->getDistr(j);
		dist::JointGMM curJGMM = boost::get<dist::JointGMM>(curDist);
		//int nVar = curJGMM.getNumVariables();
		//int nComp = curJGMM.getNumComponents();
		//cout << "nVar " << nVar << " nComp " << nComp << endl;
		cout << "joint GMM No. " << j << ":" << endl;
		cout << curJGMM << endl;
	}



	//write the dataset using the writer
	writeEddaDataset(shr_ds, "testDataImage.edda");

	//read the dataset using the reader
	shared_ptr<Dataset<Real>> shr_ds2 = loadEddaScalarDataset_noneVTK("testDataImage.edda");


	// safe to free data, after constructing the distribution
	free(varR);
	free(varG);
	free(varB);





	int numDistrArray1 = shr_ds->getNumDistrArray();
	int numDistrArray2 = shr_ds2->getNumDistrArray();
	if (numDistrArray1 != numDistrArray2){
		cout << "Joint GMM IO failed! number of arrays changed! " << endl;
		return 0;
	}
	DistrArray *array1 = shr_ds->getArray(0); //we know there is only 1 array
	DistrArray *array2 = shr_ds2->getArray(0); //we know there is only 1 array
	int n1 = array1->getLength();
	int n2 = array2->getLength();
	if (n1 != n2){
		cout << "Joint GMM IO failed! length of arrays changed! " << endl;
		return 0;
	}

	cout << "Compare one single GMM from original dataset, and the dataset after IO: " << endl << endl;

	dist::Variant curDist1 = array1->getDistr(n1 / 2);
	dist::JointGMM curJGMM1 = boost::get<dist::JointGMM>(curDist1);
	cout << "joint GMM No. " << n1 / 2 << " from the original dataset:" << endl;
	cout << curJGMM1 << endl << endl;
	dist::Variant curDist2 = array2->getDistr(n1 / 2);
	dist::JointGMM curJGMM2 = boost::get<dist::JointGMM>(curDist2);
	cout << "joint GMM No. " << n1 / 2 << " from the dataset after IO:" << endl;
	cout << curJGMM2 << endl;

	double dif = 0.0;

	for (int j = 0; j < n1; j++){
		dist::Variant curDist1 = array1->getDistr(j);
		dist::Variant curDist2 = array2->getDistr(j);

		string s = getName(curDist2);
		if (s.compare(0, 15, "JointGMM") != 0) {
			cout << "Joint GMM IO failed! arrays element not joint GMM! " << endl;
			return 0;
		}

		dist::JointGMM curJGMM1 = boost::get<dist::JointGMM>(curDist1);
		dist::JointGMM curJGMM2 = boost::get<dist::JointGMM>(curDist2);
		cout << "joint GMM No. " << j << " from IO result:" << endl;
		cout << curJGMM1 << endl << endl;

		int nVar1 = curJGMM1.getNumVariables();
		int nComp1 = curJGMM1.getNumComponents();
		int nVar2 = curJGMM2.getNumVariables();
		int nComp2 = curJGMM2.getNumComponents();
		if (nVar1 != nVar2 || nComp1 != nComp2){
			cout << "Joint GMM IO failed! nVar or nComp of at least one joint GMM changed! " << endl;
			return 0;
		}

		for (int c = 0; c < nComp1; c++)
		{
			dif += abs(curJGMM1.getWeight(c) - curJGMM2.getWeight(c));

			dist::JointGaussian jg1 = curJGMM1.getJointGaussian(c);
			const ublas_vector curMean1 = jg1.getMean();
			const ublas_matrix curCov1 = jg1.getCovariance();
			dist::JointGaussian jg2 = curJGMM2.getJointGaussian(c);
			const ublas_vector curMean2 = jg2.getMean();
			const ublas_matrix curCov2 = jg2.getCovariance();
			for (int i = 0; i < nVar; i++){
				dif += abs(curMean1(i) - curMean2(i));
			}
			for (int j = 0; j < nVar; j++){
				for (int i = j; i < nVar; i++){
					dif += abs(curCov1(j, i) - curCov2(j, i));
				}
			}
		}
	}

	cout << "the total difference between the modeler result and the IO result is: " << dif << endl;

	return 1;
}