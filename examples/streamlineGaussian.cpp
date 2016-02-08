#include <iostream>

#include "edda.h"
#include "io/file_reader.h"

using namespace std;
using namespace edda;
using namespace edda::dist;

typedef Vector3<Gaussian> Gaussian3;

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
    shared_ptr<Dataset<Gaussian3> > dataset = loadVectorData<Gaussian3>(filename);

    list<Gaussian3> seeds;
    seeds.push_back(Gaussian3(Gaussian(10.f, 0), Gaussian(10.f, 0), Gaussian(10.f, 0)));

    StreamTracer<Gaussian3, GetPositionFromDistributionMean> streamTracer(dataset);

    run(streamTracer, seeds);
  }
  return 0;
}
