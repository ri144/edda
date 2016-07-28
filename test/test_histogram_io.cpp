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


int main()
{
	shared_ary<Histogram> array (new Histogram[8], 8);
	for (int i = 0; i < 8; i++){
		float data[1];
		data[0] = 50 + i;
		//array[i] = Histogram(data, 1, 50, 60, 20);
		//array[i] = eddaComputeHistogram(data, 1, 50, 60, 20);
		array[i] = eddaComputeHistogram(data, 1, 20);
	}
	AbstractDistrArray * abstract_array = new DistrArray<Histogram>(array);

	for (int i = 0; i<8; i++)
	{
		printf("%d\n", i);
		cout << i << ": " << abstract_array->getScalar(i) <<
			": sample = " << getSample(abstract_array->getScalar(i)) << endl;
	}

	shared_ptr<Dataset<Real> > dataset = make_Dataset<Real>(
		new RegularCartesianGrid(2,2,2),
		abstract_array
	);
	writeEddaVtkDataset(dataset, "testHist.vti", "test_");
	cout << "finish histogram modeling and write into file, press any key to continue!" << endl;
	//getchar();
																															
																																
	shared_ptr<Dataset<Real> > dataset2 = loadEddaScalarDataset("testHist.vti", "test_");
	cout << dataset2->getArray()->getDistrName() << endl;
	int* dim;
	dim = dataset2->getDimension();
	cout << dim[0] << dim[1] << dim[2]  << endl;

	for (int i = 0; i < dim[0]; i++)
		for (int j = 0; j < dim[1]; j++)
			for (int k = 0; k < dim[2]; k++)
				cout << "at_comp(" <<i << "," << j << "," << k << ") : " << dataset2->at_comp(i, j, k) << endl;
  cout << "load histogram from file and sampling, press any key to finish!" << endl;
  //getchar();
  return 0;
}

