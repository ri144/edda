#ifndef VTK_EDDA_READER_H
#define VTK_EDDA_READER_H


#include <vtkIOXMLModule.h> // For export macro
#include <vtkAlgorithm.h>
#include <vtkImageAlgorithm.h>

#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>

class vtkImageData;

///
/// \brief The vtkEddaReader class reads an internal .info file format and creates a dataset with vtkSamplingArray.  This class is experimental.
///
class vtkEddaReader : public vtkImageAlgorithm
{
public:
  vtkTypeMacro(vtkEddaReader,vtkImageAlgorithm)
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkEddaReader *New();

  // Description:
  // Get the reader's output.
  vtkImageData *GetOutput();
  vtkImageData *GetOutput(int idx);

  // Description:
  // For the specified port, copy the information this reader sets up in
  // SetupOutputInformation to outInfo
  //virtual void CopyOutputInformation(vtkInformation *outInfo, int port);

  // ----------

  // Description:
  // Get/Set the name of the input file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  vtkEddaReader();
  ~vtkEddaReader();

  virtual int FillOutputPortInformation(int, vtkInformation*);

  //--------

  // The input file's name.
  char* FileName;

  virtual int RequestInformation(
    vtkInformation*, vtkInformationVector**, vtkInformationVector* );
  virtual int RequestData(
    vtkInformation*, vtkInformationVector**, vtkInformationVector* );

private:
  vtkEddaReader(const vtkEddaReader&);  // Not implemented.
  void operator=(const vtkEddaReader&);  // Not implemented.
};

#endif
