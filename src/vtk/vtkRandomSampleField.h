#ifndef VTK_RANDOM_SAMPLE_FIELD_H
#define VTK_RANDOM_SAMPLE_FIELD_H

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkDataSetAlgorithm.h"
#include "vtkDataSetAttributes.h" // needed for vtkDataSetAttributes::FieldList

#include "gmm_vtk_data_array.h"

class vtkIdTypeArray;
class vtkCharArray;
class vtkMaskPoints;

class VTKFILTERSCORE_EXPORT vtkRandomSampleField : public vtkDataSetAlgorithm
{
public:
  static vtkRandomSampleField *New();
  vtkTypeMacro(vtkRandomSampleField,vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  double Isov;

  vtkRandomSampleField();
  ~vtkRandomSampleField();

  virtual int RequestData(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);
  virtual int RequestInformation(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);

  // non-vtk functions
  virtual void InitializeData(vtkDataSet* input,
                              vtkDataSet* output);
  virtual void SampleDataArray(std::shared_ptr<edda::GmmVtkDataArray> dataArray, vtkSmartPointer<vtkFieldData> output);

private:
  vtkRandomSampleField(const vtkRandomSampleField&);  // Not implemented.
  void operator=(const vtkRandomSampleField&);  // Not implemented.

};

#endif // VTK_RANDOM_SAMPLE_FIELD_H
