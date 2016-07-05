#include <cstdlib>
#include <ctime>
#include <iostream>

#include "dataset/dataset.h"
#include "distributions/distribution.h"
#include "distributions/variant.h"
#include "dataset/abstract_data_array.h"
using namespace edda;
using namespace std;
using namespace edda::dist;

DefaultGaussianMixture gen_gmm()
{
  DefaultGaussianMixture gmm;
  for (int i=0; i<MAX_GMs; i++)
  {
    gmm.models[i].m = i;  gmm.models[i].v = 1; gmm.models[i].w=1;
  }
  gmm.normalizeWeights();

  return gmm;
}

AbstractDataArray *make_Gaussian_array() {
  shared_ary<Gaussian> array (new Gaussian[10], 10);
  AbstractDataArray * abstract_array = new ScalarArray<Gaussian>(array);
  return abstract_array;
}


AbstractDataArray * make_GMM_array() {
  shared_ary<DefaultGaussianMixture> array (new DefaultGaussianMixture[10], 10);
  for (int i=0; i<10; i++)
    array[i] = gen_gmm();
  AbstractDataArray * abstract_array = new ScalarArray<DefaultGaussianMixture>(array);
  return abstract_array;
}

AbstractDataArray * make_mixed_array() {
  shared_ary<Variant> array (new Variant[10], 10);
  // add Gaussian
  array[0] = Gaussian(1,2);
  // assign GaussianMixture
  array[1] = gen_gmm();

  AbstractDataArray * abstract_array = new ScalarArray<Variant>(array);
  return abstract_array;
}


int main()
{

  AbstractDataArray * array1 = make_Gaussian_array();
  AbstractDataArray * array2 = make_GMM_array();
  AbstractDataArray * array3 = make_mixed_array();

  int i;
  cout << "array1: " << endl;
  for (i=0; i<10; i++)
  {
    cout << i << ": " << array1->getScalar(i) <<
            ": sample = " << getSample(array1->getScalar(i)) << endl;
  }
  cout << endl ;

  cout << "array2: "<< endl;
  for (i=0; i<10; i++)
  {
    cout << i << ": " << array2->getScalar(i) <<
            ": sample = " << getSample(array2->getScalar(i)) << endl;
  }
  cout << endl;

  cout << "array3: "<< endl;
  for (i=0; i<2; i++)
  {
    cout << i << ": " << array3->getScalar(i) <<
            ": sample = " << getSample(array3->getScalar(i)) << endl;
  }
  cout << endl;

}
