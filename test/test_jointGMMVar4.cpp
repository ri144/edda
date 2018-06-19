#include "distributions/joint_GMM.h"
#include "distributions/joint_gaussian.h"
#include "vector"
#include "bmp_image.h"
#include "dataset/distr_array.h"
#include "io/edda_reader.h"
#include "io/edda_writer.h"
#include "distributions/histogram.h"
#include "dataset/distr_array.h"

using namespace edda;
using namespace std;
using namespace edda::dist;

int main()
{
    //This test data set has four variable and very small and large data range on diffrent varaible
    //check whether the training results have inf or nan value
    //load data
    FILE* fp = fopen("../edda/sample_data/wideRangeJointGMMTest.bin", "rb");
    Real* var0 = (Real*)malloc(sizeof(Real)*1000);
    Real* var1 = (Real*)malloc(sizeof(Real)*1000);
    Real* var2 = (Real*)malloc(sizeof(Real)*1000);
    Real* var3 = (Real*)malloc(sizeof(Real)*1000);
    for( int i=0; i< 1000; i++ ){
        float tmp[4];
        fread(tmp, sizeof(float), 4, fp );
        var0[i] = tmp[0];
        var1[i] = tmp[1];
        var2[i] = tmp[2];
        var3[i] = tmp[3];
    }
    fclose(fp);

    std::vector<Real*> trainSamples;
    trainSamples.clear();
    trainSamples.push_back(var0);
    trainSamples.push_back(var1);
    trainSamples.push_back(var2);
    trainSamples.push_back(var3);

    JointGMM gmm = eddaComputeJointGMM(trainSamples, 1000, 8, 2);
    std::cout << "Check whether any Nan or Inf value in the training result" << std::endl;
}
