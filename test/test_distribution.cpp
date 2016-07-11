#include <cstdlib>
#include <ctime>
#include <iostream>

#include "distributions/distribution.h"
#include "distributions/variant.h"
#include "dataset/distr_array.h"
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

GaussianMixture<2> gen_gmm2()
{
  GaussianMixture<2> gmm;
  for (int i=0; i<2; i++)
  {
    gmm.models[i].m = i;  gmm.models[i].v = 1; gmm.models[i].w=1;
  }
  gmm.normalizeWeights();

  return gmm;
}


AbstractDistrArray *make_Gaussian_array() {
  shared_ary<Gaussian> array (new Gaussian[10], 10);
  AbstractDistrArray * abstract_array = new DistrArray<Gaussian>(array);
  return abstract_array;
}


AbstractDistrArray * make_GMM_array() {
  shared_ary<DefaultGaussianMixture> array (new DefaultGaussianMixture[10], 10);
  for (int i=0; i<10; i++)
    array[i] = gen_gmm();
  AbstractDistrArray * abstract_array = new DistrArray<DefaultGaussianMixture>(array);
  return abstract_array;
}

AbstractDistrArray * make_GMM2_array() {
  shared_ary<GaussianMixture<2> > array (new GaussianMixture<2>[10], 10);
  for (int i=0; i<10; i++)
    array[i] = gen_gmm2();
  AbstractDistrArray * abstract_array = new DistrArray<GaussianMixture<2> >(array);
  return abstract_array;
}

AbstractDistrArray * make_hybrid_array() {
  shared_ary<Variant> array (new Variant[10], 10);
  // add Gaussian
  array[0] = Gaussian(1,2);
  // assign GaussianMixture
  array[1] = gen_gmm();

  AbstractDistrArray * abstract_array = new DistrArray<Variant>(array);
  return abstract_array;
}


int main()
{

  AbstractDistrArray * array1 = make_Gaussian_array();
  AbstractDistrArray * array2 = make_GMM_array();
  AbstractDistrArray * array22 = make_GMM2_array();
  AbstractDistrArray * array3 = make_hybrid_array();

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
