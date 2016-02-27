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


void vtkUncertainIsocontour::Compute(vtkDataSet* input, int *dim,
  vtkDataSet* output)
{
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
    out_ndarray->copy_to_host((float *)out_vtkArray->GetVoidPointer(0));

    // release ownership of shared_ary
    //if (r!=ReturnStatus::SUCCESS) {
    //  return 0;
    //}

    output->GetCellData()->AddArray(out_vtkArray);
  } else {
    // has cell data?
    //dataArray = shared_ptr<GmmVtkDataArray>(new GmmVtkDataArray(input->GetCellData()) );
  }

}
