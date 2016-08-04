#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iostream>

#include "io/edda_vtk_reader.h"
#include "io/edda_vtk_writer.h"
#include "distributions/histogram.h"
#include "dataset/distr_array.h"
using namespace edda;
using namespace std;
using namespace edda::dist;

void loadFloatRawData(char* path, int D0, int D1, int D2, float* data){
	FILE* fp = fopen(path, "rb");
	fread(data, sizeof(float), D0*D1*D2, fp);
	fclose(fp);
}

AbstractDistrArray* computeDistributionArrayFromBlocks(float* data, int blockSize, int dsD0, int dsD1, int dsD2, int D0, int D1, int D2){
	int nblock = dsD0*dsD1*dsD2;
	float *dataPoint = (float*)malloc(sizeof(float)*blockSize*blockSize*blockSize);
	shared_ary<Histogram> array(new Histogram[nblock], nblock);
	int blkCnt = 0;
	for (int blkD0 = 0; blkD0 < dsD0; blkD0++){
		int d0s = blkD0 * blockSize;
		int d0t = d0s + blockSize;
		if (d0t > D0) d0t = D0;
		for (int blkD1 = 0; blkD1 < dsD1; blkD1++){
			int d1s = blkD1 * blockSize;
			int d1t = d1s + blockSize;
			if (d1t > D1) d1t = D1;
			for (int blkD2 = 0; blkD2 < dsD2; blkD2++){
				int d2s = blkD2 * blockSize;
				int d2t = d2s + blockSize;
				if (d2t > D2) d2t = D2;
				int cnt = 0;

				for (int d0 = d0s; d0 < d0t; d0++){
					for (int d1 = d1s; d1 < d1t; d1++){
						for (int d2 = d2s; d2 < d2t; d2++){
							dataPoint[cnt++] = data[d0*D1 * D2 + d1*D2 + d2];
						}
					}
				}
				array[blkCnt++] = eddaComputeHistogram(dataPoint, cnt, 128, -1000, 1000);
				//array[blkCnt++] = eddaComputeHistogram(dataPoint, cnt, 128);
			}
		}
	}
	return new DistrArray<Histogram>(array);
}

void writeResampleData(shared_ptr<Dataset<Real> > dataset2, int dsD0, int dsD1, int dsD2, int upsample)
{
	int* dim;
	dim = dataset2->getDimension();

	FILE* ofp = fopen("histOutput.raw", "wb");
	float tmp[1];
	float step = 1.0 / (float)upsample;
	for (float k = 0; k < dim[2]; k += step){
		printf("%f\n", k);
		for (float j = 0; j < dim[1]; j += step){
			for (float i = 0; i < dim[0]; i += step){
				VECTOR3 pos = VECTOR3(i, j, k);
				float v;
				dataset2->at_phys(pos, v);
				tmp[0] = v;
				fwrite(tmp, sizeof(float), 1, ofp);
			}
		}
	}
	fclose(ofp);
	printf("Dimenstion: %d %d %d\n", dim[0] * upsample, dim[1] * upsample, dim[2] * upsample);
}

int main()
{
	int blockSize = 8;	//local block size
	int D0 = 100, D1 = 500, D2 = 500;	//raw data dimenstion
	char path[500] = "C:\\GravityLabDataSet\\IsabelMultiVariable\\Pf07.raw";	//raw data path
	int dsD0 = int(ceil(D0 / float(blockSize)));	//resolution of blocks
	int dsD1 = int(ceil(D1 / float(blockSize)));
	int dsD2 = int(ceil(D2 / float(blockSize)));

	//load raw data set
	float* gtDetData = (float*)malloc(sizeof(float)*D0*D1*D2);
	loadFloatRawData(path, D0, D1, D2, gtDetData);

	//create distribution array from blocks
	AbstractDistrArray * abstract_array = computeDistributionArrayFromBlocks(gtDetData, blockSize, dsD0, dsD1, dsD2, D0, D1, D2);

	//writer: write distribution array to disk
	shared_ptr<Dataset<Real> > dataset = make_Dataset<Real>(
		new RegularCartesianGrid(dsD2, dsD1, dsD0),
		abstract_array
		);
	writeEddaVtkDataset(dataset, "testHist.vti", "test_");

	//reader: read distribution array from disk
	shared_ptr<Dataset<Real> > dataset2 = loadEddaScalarDataset("testHist.vti", "test_");

	//sample and output raw data to disk
	writeResampleData(dataset2, dsD0, dsD1, dsD2, 8);

	cout << "Done, press any key to finish!" << endl;
	getchar();
}