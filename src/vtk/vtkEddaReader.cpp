
#include <cstdlib>
#include <ctime>

#include "vtkEddaReader.h"

#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"

#include "io/file_reader.h"
#include "vtkSamplingArray.h"

vtkStandardNewMacro(vtkEddaReader);

using namespace edda;

//----------------------------------------------------------------------------
vtkEddaReader::vtkEddaReader()
{
  this->FileName = NULL;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);

  // initialize randomness
  srand(time(NULL));
}

//----------------------------------------------------------------------------
vtkEddaReader::~vtkEddaReader()
{
}

//----------------------------------------------------------------------------
void vtkEddaReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
vtkImageData* vtkEddaReader::GetOutput()
{
  return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkImageData* vtkEddaReader::GetOutput(int idx)
{
  return vtkImageData::SafeDownCast(this->GetOutputDataObject(idx));
}


//----------------------------------------------------------------------------
int vtkEddaReader::FillOutputPortInformation(int, vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkEddaReader::RequestInformation(
  vtkInformation*, vtkInformationVector**, vtkInformationVector* outVec)
{
  vtkInformation* outInfo = outVec->GetInformationObject(0);
  int extent[6];

  // ... read file to find available extent ...
  printf("Loading file...\n");
  extent[0] = 0;
  extent[1] = 3;
  extent[2] = 0;
  extent[3] = 3;
  extent[4] = 0;
  extent[5] = 3;

  //store that in the pipeline
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

  // ... store other information ...

  return 1;
}

  //----------------------------------------------------------------------------
int vtkEddaReader:: RequestData(
  vtkInformation*, vtkInformationVector**, vtkInformationVector* outVec)
{
  vtkInformation* outInfo = outVec->GetInformationObject(0);
  vtkImageData* outData = vtkImageData::SafeDownCast
    (outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int extent[6] = {0,-1,0,-1,0,-1};
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);

  outData->SetExtent(extent);

  // ... read data for this extent from the file ...
  printf("Reading data...\n");

  std::shared_ptr<Dataset<VECTOR3> > dataset = loadVectorData<VECTOR3>(this->FileName, true);

  // Specify the size of the image data
#if VTK_MAJOR_VERSION <= 5
  //outData->SetNumberOfScalarComponents(3);
  //outData->SetScalarTypeToFloat();
#else
  //outData->AllocateScalars(VTK_FLOAT,3);
#endif
  const int *dim = dataset->getDimension();
  outData->SetExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);

  vtkSamplingArray *array = vtkSamplingArray::New();
  array->SetEddaArray(dataset->getArray());
  array->SetName("VectorSampling");

  outData->GetPointData()->AddArray(array);
  outData->GetPointData()->SetActiveVectors("VectorSampling");
#if 0
  vtkFloatArray *array = vtkFloatArray::New();
  array->SetNumberOfComponents(3);
  array->SetNumberOfTuples(dataset->getArray()->getLength());
  for (int i=0; i<dataset->getArray()->getLength(); i++) {
    VECTOR3 v = boost::any_cast<VECTOR3>( dataset->getArray()->getItem(i) );
    array->SetTuple3(i, v[0], v[1], v[2]);
  }
#endif

  //for (int i=0; i<64; i++)
  //  dynamic_cast<vtkFloatArray *>(outData->GetPointData()->GetScalars())->SetTuple1(i, (float)i);

  return 1;
}

#if 0

//----------------------------------------------------------------------------
// Note that any changes (add or removing information) made to this method
// should be replicated in CopyOutputInformation
void vtkEddaReader::SetupOutputInformation(vtkInformation *outInfo)
{
  this->Superclass::SetupOutputInformation(outInfo);

  outInfo->Set(vtkDataObject::ORIGIN(), this->Origin, 3);
  outInfo->Set(vtkDataObject::SPACING(), this->Spacing, 3);
}


//----------------------------------------------------------------------------
void vtkEddaReader::CopyOutputInformation(vtkInformation *outInfo, int port)
{
  this->Superclass::CopyOutputInformation(outInfo, port);
  vtkInformation *localInfo =
    this->GetExecutive()->GetOutputInformation(port);
  if (localInfo->Has(vtkDataObject::ORIGIN()))
    {
    outInfo->CopyEntry(localInfo, vtkDataObject::ORIGIN());
    }
  if (localInfo->Has(vtkDataObject::SPACING()))
    {
    outInfo->CopyEntry(localInfo, vtkDataObject::SPACING());
    }
}


#endif
