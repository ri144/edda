#include <iostream>

#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/file_reader.h"
#include "filters/stream_tracer.h"

using namespace std;
using namespace edda;
using namespace edda::dist;

typedef Gaussian<float> Gaussianf;
typedef Vector3<Gaussianf> Gaussianf3;

shared_ary<Gaussianf3> load_gaussian_dataset(string meanfile, string stdfile, int *dim) {
  cout << "Loading..." << endl;
  size_t len = dim[0]*dim[1]*dim[2];
  shared_ary<Tuple3<float> > pMean3 = loadRawFile<Tuple3<float> >(meanfile, len);
  shared_ary<Tuple3<float> > pStd3 = loadRawFile<Tuple3<float> >(stdfile, len);

  // Create Gaussian array
  Gaussianf3 x;
  shared_ary<Gaussianf3> data(new Gaussianf3[len], len);
  int i;
  for (i=0; i<len; i++)
  {
    data[i] = Gaussianf3(
                Gaussianf(pMean3[i][0], pStd3[i][0]),
                Gaussianf(pMean3[i][1], pStd3[i][1]),
                Gaussianf(pMean3[i][2], pStd3[i][2]) );
  }
  return data;
}

template<typename DataType, class GetPositionPolicy>
void run(StreamTracer<DataType, GetPositionPolicy> &streamTracer, const list<DataType> &seeds)
{
    list<list<DataType> >traces;
    streamTracer.compute(seeds, traces);

    for (auto traceItr = traces.begin(); traceItr != traces.end(); ++traceItr)
    {
        auto &singleTrace = *traceItr;
        cout << "trace: ";
        for (auto posItr = singleTrace.begin(); posItr!=singleTrace.end(); ++posItr)
        {
            cout << *posItr << " ";
        }
        cout << endl;
    }
    cout << endl;
}


int main(int argc, char **argv) {
    cout << "Input arguments: <mean file> <std file> <w> <h> <d>" << endl;
    if (argc<6)
        return -1;
    string meanfile, stdfile;
    meanfile = argv[1];
    stdfile = argv[2];
    int dim[3];
    dim[0] = atoi(argv[3]);
    dim[1] = atoi(argv[4]);
    dim[2] = atoi(argv[5]);
    cout << "dim: " << dim[0] << "," << dim[1] << "," << dim[2] << endl;

    shared_ary<Gaussianf3> distData = load_gaussian_dataset(meanfile, stdfile, dim);

    /////////////////////////////////////
    cout << "Streamline computed by Monte-Carlo sampling from the uncertain velocity field:" << endl;
    {
        shared_ptr<Dataset<VECTOR3> > dataset  = make_Dataset<VECTOR3>(
                    new RegularCartesianGrid (dim[0], dim[1], dim[2]),
                    new DataArray<Gaussianf3, GetItemSampledVector> (distData )
        );

        list<VECTOR3> seeds;
        seeds.push_back(VECTOR3(10.f, 10.f, 10.f));

        StreamTracer<VECTOR3> streamTracer(dataset);

        run(streamTracer, seeds);
    }
    /////////////////////////////////////
    cout << "Streamline computed with distribution:" << endl;
    {
        shared_ptr<Dataset<Gaussianf3> > dataset = make_Dataset<Gaussianf3>(
                    new RegularCartesianGrid (dim[0], dim[1], dim[2]),
                    new DataArray<Gaussianf3> (distData)
                );


        list<Gaussianf3> seeds;
        seeds.push_back(Gaussianf3(Gaussianf(10.f, 0), Gaussianf(10.f, 0), Gaussianf(10.f, 0)));

        StreamTracer<Gaussianf3, GetPositionFromDistributionMean>
                streamTracer(dataset, GetPositionFromDistributionMean());

        run(streamTracer, seeds);
    }
    return 0;
}
