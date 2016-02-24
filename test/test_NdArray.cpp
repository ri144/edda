#include "test_common.h"
#include <core/ndarray.h>

using namespace edda;

int main ()
{
  int dim=1;
  NdArray<Real> p(1, &dim);
  p.set_val({0}, 3);
  {
          NdArray<Real> q = p;
          TEST(q.get_val({0})==3);
  }
  return 0;

}
