#include <iostream>
#include <vtkSmartPointer.h>

#include "vtkEddaReader.h"

#include "edda.h"

using namespace std;
//using namespace edda;

#define vsp_new(cls,x) vtkSmartPointer<cls> x = vtkSmartPointer<cls>::New()

int main(int argc, char *argv[])
{
  if (argc<=1) {
    cout << "Please provide file name." << endl;
    return 1;
  }
  vsp_new(vtkEddaReader, reader);
  reader->SetFileName(argv[1]);
  reader->Update();

  vtkImageData *image;
  image = reader->GetOutput();

  cout << image->GetPointData()->GetNumberOfArrays() << endl;

  vtkDataArray * array;
  array = vtkDataArray::SafeDownCast(image->GetPointData()->GetArray("VectorSampling"));
  cout << "points=" << image->GetNumberOfPoints() << endl;
  for (int i=0; i<image->GetNumberOfPoints(); i++)
  {
    double *v = array->GetTuple3(i);
    printf("%f %f %f\n", v[0], v[1], v[2]);
  }


  return 0;
}
