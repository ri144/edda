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
  cout << "Input arguments: <info file>" << endl;
  if (argc<2)
      return -1;
  string filename;
  filename = argv[1];

  /////////////////////////////////////
  cout << "Streamline computed by Monte-Carlo sampling from the uncertain velocity field:" << endl;
  {
    shared_ptr<Dataset<VECTOR3> > dataset = loadVectorData<VECTOR3, GetItemSampledVector>(filename);

    list<VECTOR3> seeds;
    seeds.push_back(VECTOR3(10.f, 10.f, 10.f));

    StreamTracer<VECTOR3> streamTracer(dataset);

    run(streamTracer, seeds);
  }
  /////////////////////////////////////
  cout << "Streamline computed with distribution:" << endl;
  {
    shared_ptr<Dataset<Gaussianf3> > dataset = loadVectorData<Gaussianf3>(filename);

    list<Gaussianf3> seeds;
    seeds.push_back(Gaussianf3(Gaussianf(10.f, 0), Gaussianf(10.f, 0), Gaussianf(10.f, 0)));

    StreamTracer<Gaussianf3, GetPositionFromDistributionMean> streamTracer(dataset);

    run(streamTracer, seeds);
  }
  return 0;
}
