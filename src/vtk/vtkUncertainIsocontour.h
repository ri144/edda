#ifndef VTK_UNCERTAIN_ISOCONTOUR_H
#define VTK_UNCERTAIN_ISOCONTOUR_H

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkDataSetAlgorithm.h"
#include "vtkDataSetAttributes.h" // needed for vtkDataSetAttributes::FieldList

class vtkIdTypeArray;
class vtkCharArray;
class vtkMaskPoints;

class VTKFILTERSCORE_EXPORT vtkUncertainIsocontour : public vtkDataSetAlgorithm
{
public:
  static vtkUncertainIsocontour *New();
  vtkTypeMacro(vtkUncertainIsocontour,vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  ///
  /// \brief Get isovalue
  ///
  vtkGetMacro(Isov, double);

  ///
  /// \brief Set isovalue
  ///
  vtkSetMacro(Isov, double);

protected:
  double Isov;

  vtkUncertainIsocontour();
  ~vtkUncertainIsocontour();

  virtual int RequestData(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);
  virtual int RequestInformation(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);

  virtual void InitializeData(vtkDataSet* input,
                              vtkDataSet* output);

  virtual void Compute(vtkDataSet* input, int *dim,
                              vtkDataSet* output);

private:
  vtkUncertainIsocontour(const vtkUncertainIsocontour&);  // Not implemented.
  void operator=(const vtkUncertainIsocontour&);  // Not implemented.

};

#endif // VTK_UNCERTAIN_ISOCONTOUR_H
