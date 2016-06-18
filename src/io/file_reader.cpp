#include <cassert>
#include <distributions/gaussian.h>
#include <io/file_reader.h>
#include <dataset/dataset.h>
#include <core/vector_matrix.h>

using namespace std;

namespace edda{

using namespace dist;
/// array loader

shared_ary<Gaussian> loadGaussianRawArray(string meanfile, string stdfile, size_t len)
{
  shared_ary<float> pMean = loadRawFile<float>(meanfile, len);
  shared_ary<float> pStd = loadRawFile<float>(stdfile, len);

  // Create Gaussian array
  Gaussian *pData = new Gaussian[len];
  size_t i;
  for (i=0; i<len; i++)
  {
    pData[i] = Gaussian(pMean[i], pStd[i]*pStd[i]);
  }
  // return smart pointer of the array
  return shared_ary<Gaussian> (pData, len);
}

shared_ary<Gaussian3> loadVec3GaussianRawArray(string meanfile, string stdfile, size_t len) {
  shared_ary<Tuple3<float> > pMean3 = loadRawFile<Tuple3<float> >(meanfile, len);
  shared_ary<Tuple3<float> > pStd3 = loadRawFile<Tuple3<float> >(stdfile, len);

  // Create Gaussian array
  Gaussian3 x;
  shared_ary<Gaussian3> data(new Gaussian3[len], len);
  size_t i;
  for (i=0; i<len; i++)
  {
    data[i] = Gaussian3(
                Gaussian(pMean3[i][0], pStd3[i][0]*pStd3[i][0]),
                Gaussian(pMean3[i][1], pStd3[i][1]*pStd3[i][1]),
                Gaussian(pMean3[i][2], pStd3[i][2]*pStd3[i][2]) );
  }
  return data;
}

///////////////////////////////////////
/// dataset creator

/// create a regular grid dataset of gaussian distribution
shared_ptr<Dataset<Gaussian> > loadGaussianRegularGrids(string &meanfile, string &stdfile, int dim[3])
{
  shared_ary<Gaussian> array = loadGaussianRawArray(meanfile, stdfile, dim[0]*dim[1]*dim[2]);
  Dataset<Gaussian> *dataset = new Dataset<Gaussian> (
        new RegularCartesianGrid (dim[0], dim[1], dim[2]),
        new ScalarArray< Gaussian >( array )
    );
  return shared_ptr<Dataset<Gaussian> > (dataset);
}

/// create a regular grid dataset of random values from gaussian distributions
shared_ptr<Dataset<double> > loadGaussianSamplingRegularGrids(string &meanfile, string &stdfile, int dim[3])
{
  shared_ary<Gaussian> array = loadGaussianRawArray(meanfile, stdfile, dim[0]*dim[1]*dim[2]);
  Dataset<double> *dataset = new Dataset<double> (
        new RegularCartesianGrid (dim[0], dim[1], dim[2]),
        new AbstractSamplingArray( new ScalarArray< Gaussian >( array ) )
    );
  return shared_ptr<Dataset<double> > (dataset);
}

/// create a regular grid dataset of 3d gaussian distribution
shared_ptr<Dataset<Gaussian3> > loadVec3GaussianRegularGrids(string &meanfile, string &stdfile, int dim[3])
{
  shared_ary<Gaussian3> array = loadVec3GaussianRawArray(meanfile, stdfile, dim[0]*dim[1]*dim[2]);
  Dataset<Gaussian3> *dataset = new Dataset<Gaussian3> (
            new RegularCartesianGrid (dim[0], dim[1], dim[2]),
            new VectorArray<Gaussian, 3> (array) );
  return shared_ptr<Dataset<Gaussian3> > (dataset);
}

/// create a regular grid dataset of random samples from 3d gaussian distribution
shared_ptr<Dataset<VECTOR3> > loadVec3GaussianSamplingRegularGrids(string &meanfile, string &stdfile, int dim[3])
{
  shared_ary<Gaussian3> array = loadVec3GaussianRawArray(meanfile, stdfile, dim[0]*dim[1]*dim[2]);
  Dataset<VECTOR3> *dataset = new Dataset<VECTOR3> (
              new RegularCartesianGrid (dim[0], dim[1], dim[2]),
              new AbstractSamplingArray( new VectorArray<Gaussian, 3> (array) )
  );
  return shared_ptr<Dataset<VECTOR3> > (dataset);
}

namespace detail {
void print(boost::property_tree::ptree const& pt)
{
    boost::property_tree::ptree::const_iterator end = pt.end();
    std::cout << "{" << std::endl;
    for (boost::property_tree::ptree::const_iterator it = pt.begin(); it != end; ++it) {
        std::cout << it->first << ": " << it->second.get_value<std::string>() << std::endl;
        print(it->second);
    }
    std::cout << "}" << std::endl;
}
} // detail


} // edda
