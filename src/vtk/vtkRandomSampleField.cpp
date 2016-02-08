#include <iostream>

#include "vtkRandomSampleField.h"

#include "vtk_common.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "gmm_vtk_data_array.h"

using namespace std;
using namespace edda;

vtkStandardNewMacro(vtkRandomSampleField);

//----------------------------------------------------------------------------
vtkRandomSampleField::vtkRandomSampleField()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkRandomSampleField::~vtkRandomSampleField()
{
}


void vtkRandomSampleField::InitializeData(vtkDataSet* input,
  vtkDataSet* output)
{
  // First, copy the input geometry to the output as a starting point
  output->CopyStructure(input);
}

//----------------------------------------------------------------------------
int vtkRandomSampleField::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  this->InitializeData(input, output);

  // process point data
  shared_ptr<GmmVtkDataArray> dataArray(new GmmVtkDataArray(input->GetPointData()) );

  // has point data?
  if (dataArray->getLength() > 0) {
    // create output array
    vsp_new(vtkFloatArray, out_vtkArray);
    out_vtkArray->SetNumberOfComponents(1);
    out_vtkArray->SetNumberOfTuples(dataArray->getLength());
    out_vtkArray->SetName("RandomSample");

    for (size_t i=0; i<dataArray->getLength(); i++)
    {
      *(float *)out_vtkArray->GetVoidPointer(i) = dist::getSample(boost::any_cast<dist::GaussianMixture>( dataArray->getItem(i) ) );
    }
    output->GetPointData()->AddArray(out_vtkArray);
  }

  // process cell data
  dataArray = shared_ptr<GmmVtkDataArray>(new GmmVtkDataArray(input->GetCellData()) );

  if (dataArray->getLength())
  {
    // create output array
    vsp_new(vtkFloatArray, out_vtkArray);
    out_vtkArray->SetNumberOfComponents(1);
    out_vtkArray->SetNumberOfTuples(dataArray->getLength());
    out_vtkArray->SetName("RandomSample");

    for (size_t i=0; i<dataArray->getLength(); i++)
    {
      *(float *)out_vtkArray->GetVoidPointer(i) = dist::getSample( boost::any_cast<dist::GaussianMixture>( dataArray->getItem(i) ) );
    }
    output->GetCellData()->AddArray(out_vtkArray);
  }

  return 1;
}

//----------------------------------------------------------------------------
int vtkRandomSampleField::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
               6);

  return 1;
}

//----------------------------------------------------------------------------
int vtkRandomSampleField::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int usePiece = 0;

  // What ever happened to CopyUpdateExtent in vtkDataObject?
  // Copying both piece and extent could be bad.  Setting the piece
  // of a structured data set will affect the extent.
  vtkDataObject* output = outInfo->Get(vtkDataObject::DATA_OBJECT());
  if (output &&
      (!strcmp(output->GetClassName(), "vtkUnstructuredGrid") ||
       !strcmp(output->GetClassName(), "vtkPolyData")))
    {
    usePiece = 1;
    }

  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

    inInfo->Set(
      vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT()), 6);

  return 1;
}

//----------------------------------------------------------------------------
void vtkRandomSampleField::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
