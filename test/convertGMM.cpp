#include <vector>
#include <cstdio>
#include <string>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include "vtk/vtk_common.h"
#include <io/file_reader.h>
#include <cmath>

using namespace std;
using namespace edda;

int main(int argc, char **argv)
{
  printf("This program converts raw arrays into gmm in vtk formats. The dataset is supposed to be regular grids.\n");
  printf("Input: components w h d mean1 std1 [weight1 mean2 std2 weight2 ...]\n");
  printf("mean1, std1, weight1... are file names\n");
  if (argc <= 6)
    return 1;
  int components = atoi(argv[1]);
  int w = atoi(argv[2]);
  int h = atoi(argv[3]);
  int d = atoi(argv[4]);
  int n = w*h*d;

  vsp_new(vtkImageData, image);
  image->SetExtent(0, w-1, 0, h-1, 0, d-1);
  vector<string> filenames;
  int i;
  for (i=5; i<argc; i++)
    filenames.push_back(argv[i]);
  for (i=0; i<filenames.size(); i++)
  {
    vsp_new(vtkFloatArray , array);
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(n);
    printf("Reading %s\n", filenames[i].c_str());
    FILE *fp = fopen(filenames[i].c_str(), "rb");
    fread(array->GetVoidPointer(0), n, sizeof(float)*components, fp);
    char name[256];
    if (i%3==0) {
      sprintf(name, "mean%d", i/3);
      array->SetName(name);
    } else if (i%3==1) {
      sprintf(name, "var%d", i/3);
      array->SetName(name);
      // convert from std to variance.
      for (int j=0; j<n; j++)
        *(float *)array->GetVoidPointer(j) = powf(*(float *)array->GetVoidPointer(j) , 2.);
    } else {
      sprintf(name, "weight%d", i/3+1);
      array->SetName(name);
    }

    image->GetPointData()->AddArray(array);
  }
  vsp_new(vtkXMLImageDataWriter, writer);
  writer->SetFileName("merged.vti");
  writer->SetInputData(image);
  writer->Write();
  printf("Saved to %s.\n", "merged.vti");
  return 0;
}
