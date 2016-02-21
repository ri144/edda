#include <iostream>

#include "vtkUncertainIsocontour.h"

#include "vtk_common.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCharArray.h"
#include "vtkIdTypeArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "gmm_vtk_data_array.h"
#include <filters/uncertain_isocontour.h>

using namespace std;
using namespace edda;

vtkStandardNewMacro(vtkUncertainIsocontour);

//----------------------------------------------------------------------------
vtkUncertainIsocontour::vtkUncertainIsocontour()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkUncertainIsocontour::~vtkUncertainIsocontour()
{
}


void vtkUncertainIsocontour::InitializeData(vtkDataSet* input,
  vtkDataSet* output)
{
  // First, copy the input to the output as a starting point
  output->CopyStructure(input);
}

//----------------------------------------------------------------------------
int vtkUncertainIsocontour::RequestData(
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

  // check whether input is structured grids or image data
  int *dim = NULL;
  vtkImageData *image = vtkImageData::SafeDownCast(input);
  vtkStructuredGrid *sgrid = vtkStructuredGrid::SafeDownCast(input);
  if (image) {
    dim = image->GetDimensions();
  }
  else if (sgrid) {
    dim = sgrid->GetDimensions();
  }
  else {
    return 0;
  }

  this->InitializeData(input, output);

  // process point data
  shared_ptr<GmmVtkDataArray> dataArray(new GmmVtkDataArray(input->GetPointData()) );
  shared_ptr<GmmNdArray> gmmNdArray = dataArray->genNdArray();

  // has point data?
  if (dataArray->getLength() > 0) {
    int out_length = (dim[0]-1)*(dim[1]-1)*(dim[2]-1);

    shared_ptr<NdArray<float> > out_ndarray;

    //ReturnStatus r = levelCrossingSerial(dataArray.get(), dim, this->Isov, (float *)out_vtkArray->GetVoidPointer(0));
    levelCrossing(gmmNdArray, dim, this->Isov, out_ndarray);

    // create output array
    vsp_new(vtkFloatArray, out_vtkArray);
    out_vtkArray->SetNumberOfComponents(1);
    out_vtkArray->SetNumberOfTuples(out_length);
    out_vtkArray->SetName("ProbField");
    // copy from device to host
    thrust::copy( out_ndarray->data().begin(), out_ndarray->data().end(), (float *)out_vtkArray->GetVoidPointer(0) );

    // release ownership of shared_ary
    //if (r!=ReturnStatus::SUCCESS) {
    //  return 0;
    //}

    output->GetCellData()->AddArray(out_vtkArray);
  } else {
    // has cell data?
    //dataArray = shared_ptr<GmmVtkDataArray>(new GmmVtkDataArray(input->GetCellData()) );
  }

  //output->PrintSelf(cout, vtkIndent(0));

  return 1;
}

//----------------------------------------------------------------------------
int vtkUncertainIsocontour::RequestInformation(
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
int vtkUncertainIsocontour::RequestUpdateExtent(
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
void vtkUncertainIsocontour::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
