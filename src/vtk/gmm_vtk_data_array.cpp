#include <algorithm>
#include <cstring>
#include "gmm_vtk_data_array.h"
using namespace std;
namespace edda{
using namespace dist;

GmmVtkDataArray::GmmVtkDataArray(vtkFieldData *fieldData, const char *arrayNamePrefix)  {
  char meanArrayName[1024];
  char stdevArrayName[1024];
  char varArrayName[1024];
  char weightArrayName[1024];
  for (int i=0; i<fieldData->GetNumberOfArrays(); i++)
  {
    sprintf(meanArrayName, "%smean%d", arrayNamePrefix, (int)i/3 );
    sprintf(stdevArrayName, "%sstdev%d", arrayNamePrefix, (int)i/3 );
    sprintf(varArrayName, "%svar%d", arrayNamePrefix, (int)i/3 );
    sprintf(weightArrayName, "%sweight%d", arrayNamePrefix, (int)i/3 );

    vtkSmartPointer<vtkDataArray> meanArray = fieldData->GetArray(meanArrayName);
    vtkSmartPointer<vtkDataArray> stdevArray = fieldData->GetArray(stdevArrayName);
    vtkSmartPointer<vtkDataArray> varArray = fieldData->GetArray(varArrayName);
    vtkSmartPointer<vtkDataArray> weightArray = fieldData->GetArray(weightArrayName);
    if (i==0 && meanArray && weightArray==0) {
      // allows when only one mean and one variance are provided
      // create weight Array
      weightArray = vtkSmartPointer<vtkFloatArray>::New();
      int n = meanArray->GetNumberOfTuples();
      weightArray->SetNumberOfComponents(1);
      weightArray->SetNumberOfTuples(n);
      for (int j=0; j<n; j++)
        *(float *)weightArray->GetVoidPointer(j) = 1.f ;
      n = weightArray->GetNumberOfTuples();
    }
    if (meanArray && varArray && weightArray) {
      arrays.push_back(meanArray);
      arrays.push_back(varArray);
      arrays.push_back(weightArray);
    } else if (meanArray && stdevArray && weightArray) {
      // convert stdev to variance
      arrays.push_back(meanArray);
      vsp_new(vtkFloatArray, varArray);
      int n = stdevArray->GetNumberOfTuples();
      varArray->SetNumberOfComponents(1);
      varArray->SetNumberOfTuples(n);
      for (int j=0; j<n; j++)
        *(float *)varArray->GetVoidPointer(j) = pow(stdevArray->GetTuple1(j), 2) ;
      arrays.push_back(varArray);
      arrays.push_back(weightArray);
    }
  }
  if (arrays.size() == 0) {
    printf("Warning: no array assigned to GmmVtkArray\n");
    return;
  }
  length = arrays[0]->GetNumberOfTuples();
  for (size_t i=1; i<arrays.size(); i++)
  {
    length = std::min(length, (size_t)arrays[i]->GetNumberOfTuples());
  }
}

GmmVtkDataArray::GmmVtkDataArray(std::vector<vtkSmartPointer<vtkDataArray> > arrays_) {
  if (arrays_.size() == 0) {
    printf("Warning: no array assigned to GmmVtkArray\n");
    return;
  }
  if (arrays_.size() % 3 != 0) {
    printf("Warning: GmmVtkArray: some input arrays are truncated\n");
  }
  for (size_t i=0; i<arrays_.size()/3; i++) {
    this->arrays.push_back(arrays_[i]);
    this->arrays.push_back(arrays_[i+1]);
    this->arrays.push_back(arrays_[i+2]);
  }

  length = this->arrays[0]->GetNumberOfTuples();
  for (size_t i=0; i<this->arrays.size(); i++)
  {
    length = min(length, (size_t)this->arrays[i]->GetNumberOfTuples());
  }
}

boost::any GmmVtkDataArray::getItem(size_t idx) {
  std::vector<GMMTuple> models ( arrays.size()/3 );
  for (size_t i=0; i<arrays.size(); i++) {
    models[i/3].p[i%3] = arrays[i]->GetTuple1(idx);
  }
  return boost::any( GaussianMixture(models) );
}

void GmmVtkDataArray::setItem(size_t idx, const boost::any &item) {
  // not tested
  GaussianMixture gmm = boost::any_cast<GaussianMixture>( item );
  for (size_t i=0; i<arrays.size(); i++) {
    arrays[i]->SetTuple1(idx, gmm.models[i/3].p[i%3]);
  }
}
}; //edda
