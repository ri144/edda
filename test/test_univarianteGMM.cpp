#include "distributions/joint_GMM.h"
#include "distributions/joint_gaussian.h"
#include "vector"
#include "bmp_image.h"
#include "dataset/distr_array.h"
#include "io/edda_reader.h"
#include "io/edda_writer.h"
#include "distributions/gmm.h"
#include "dataset/distr_array.h"
#include "distributions/estimate_gmm.h"

using namespace edda;
using namespace std;
using namespace edda::dist;

float getVoxel(float* data, int d0, int d1, int d2, int D0, int D1, int D2)
{
	int idx = d2*D1*D0 + d1*D0 + d0;
	return data[idx];
}

void setVoxel(float* data, int d0, int d1, int d2, int D0, int D1, int D2, float value)
{
	int idx = d2*D1*D0 + d1*D0 + d0;
	data[idx] = value;
}

//using image to test joint Gaussian GMM
int main()
{
	Real data[100] = {  0.10405,  0.528,     -0.31456,   -1.34501,   -1.29526,    0.07432,
						-0.19956,   -0.6546,     0.31801,   -0.89027,   50.11134,   49.98048,
						49.16001,   47.70179,   51.45653,   50.31664,   47.33587,   49.57357,
						50.39379,   49.77186,   50.58033,   49.02673,   50.17517,   49.94652,
						49.81694,   49.77897,   50.19976,   50.93272,   49.46988,   49.59276,
						50.16056,   49.87985,   50.3856 ,   50.71829,   51.29119,   49.88356,
						47.7227 ,   49.93038,   50.35387,   49.81304,   99.84676,   97.56749,
						100.50798,   99.67597,   98.48892,   99.12858,   99.13517,  100.60875,
						100.56164,  101.51475,  100.64792,   98.64835,   98.59079,  101.13073,
						101.56669,   99.76225,  100.5588 ,   98.49511,   98.05608,   98.82598,
						99.64281,   99.47862,   99.76989,   99.50899,  100.6793 ,  101.42755,
						100.0362 ,  102.03   ,   99.3656 ,   99.4749 ,  100.38773,   99.6452,
						101.17705,   99.35889,  101.32269,  100.19418,  102.56545,   99.53589,
						99.79731,  100.14565,   97.81897,  100.60227,  100.48085,  100.10932,
						98.4556 ,   98.45344,  100.58662,  101.17518,  101.59446,   99.10456,
						98.9692 ,   99.72806,   98.02427,   99.41107,  100.85179,  101.6346,
						100.27916,  101.64055,  100.41087,  100.19136};

	GMM gmm = eddaComputeGMM(data, 100, 3);
	std::cout << gmm << std::endl;
}

