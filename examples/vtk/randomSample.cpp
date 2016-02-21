// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

/// Read a distribution data set and generate a new field of random samples.

#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <vtkPiecewiseFunction.h>
#include <vtkSmartPointer.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkVolume.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkColorTransferFunction.h>
#include <vtkOutlineFilter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkCellDataToPointData.h>
#include <vtkXMLStructuredGridWriter.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include "io/file_reader.h"
#include "io/file_writer.h"
#include "io/path.h"
#include "vtk/vtk_common.h"
#include "vtk/vtkRandomSampleField.h"

using namespace std;
using namespace edda;

vtkSmartPointer<vtkDataSet> process_vtk_file(string vtk_file)
{
  vtkSmartPointer <vtkDataSet> dataset;
  if (getFileExtension(vtk_file).compare("vts")==0) {
    vsp_new(vtkXMLStructuredGridReader, reader);
    reader->SetFileName(vtk_file.c_str());
    reader->Update();
    dataset = reader->GetOutput();

  } else if (getFileExtension(vtk_file).compare("vti")==0) {
    vsp_new(vtkXMLImageDataReader, reader);
    reader->SetFileName(vtk_file.c_str());
    reader->Update();
    dataset = reader->GetOutput();
  } else {
    printf("File format not supported\n");
    exit(1);
  }

  // Edda filter: random sample field
  vsp_new(vtkRandomSampleField , randomSample);
  randomSample->SetInputData(dataset);
  randomSample->Update();

  // cell data to point data
  vsp_new(vtkCellDataToPointData, cell2point);
  cell2point->SetInputData(randomSample->GetOutput());
  cell2point->Update();

  cell2point->GetOutput()->PrintSelf(cout, vtkIndent(0));

  printf("Saving output file\n");
  if (getFileExtension(vtk_file).compare("vti")==0) {
    vsp_new(vtkXMLImageDataWriter, writer);
    writer->SetFileName("ProbField.vti");
    writer->SetInputData(randomSample->GetOutput());
    writer->Write();
  } else {
    vsp_new(vtkXMLStructuredGridWriter, writer);
    writer->SetFileName("ProbField.vts");
    writer->SetInputData(randomSample->GetOutput());
    writer->Write();

  }

  return cell2point->GetOutput();
}


int main(int argc, char **argv) {
  cout << "randomSample <info file> " << endl;
  if (argc<=1)
    return -1;
  string input_file = argv[1];

  vtkSmartPointer<vtkDataSet> probField;
  probField = process_vtk_file(input_file);


  // Volume render
  vsp_new(vtkPiecewiseFunction, alphaChannelFunc);
  alphaChannelFunc->AddPoint(0, 0);
  alphaChannelFunc->AddPoint(1., 1.);

  vsp_new(vtkColorTransferFunction, colorFunc);
  colorFunc->AddRGBPoint(0, 0.0, 0.0, 1.0);
  colorFunc->AddRGBPoint(.5, 1.0, 1.0, 1.0);
  colorFunc->AddRGBPoint(1., 1.0, 0.0, 0.0);

  vsp_new(vtkVolumeProperty, volumeProperty);
  volumeProperty->ShadeOff();
  volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);
  volumeProperty->SetColor(colorFunc);
  volumeProperty->SetScalarOpacity(alphaChannelFunc);

  vsp_new(vtkSmartVolumeMapper, volumeMapper);
  volumeMapper->SetBlendModeToComposite(); // composite first
  volumeMapper->SetInputData(probField);

  vsp_new(vtkVolume, volume);
  volume->SetMapper(volumeMapper);
  volume->SetProperty(volumeProperty);

  // outline
  vsp_new(vtkOutlineFilter, outline);
  outline->SetInputData(probField);
  vsp_new(vtkPolyDataMapper, outlineMapper);
  outlineMapper->SetInputConnection(outline->GetOutputPort());
  vsp_new(vtkActor, outlineActor);
  outlineActor->SetMapper(outlineMapper);
  outlineActor->GetProperty()->SetColor(.5,.5,.5);

  // Render window
  vsp_new(vtkRenderer, renderer);
  renderer->AddViewProp(volume);
  renderer->AddActor(outlineActor);
  renderer->SetBackground(0, 0, 0);
  renderer->ResetCamera();

  vsp_new(vtkRenderWindow, renderWin);
  renderWin->AddRenderer(renderer);
  vsp_new(vtkRenderWindowInteractor, renderInteractor);
  renderInteractor->SetRenderWindow(renderWin);

  vsp_new(vtkInteractorStyleTrackballCamera, style);
  renderInteractor->SetInteractorStyle(style);

  renderWin->SetSize(400, 400);
  renderInteractor->Initialize();

  renderWin->Render();

  // 3D texture mode. For coverage.
  //volumeMapper->SetRequestedRenderModeToRayCastAndTexture();
  //renderWin->Render();

  // Software mode, for coverage. It also makes sure we will get the same
  // regression image on all platforms.
  //volumeMapper->SetRequestedRenderModeToRayCast();
  //renderWin->Render();

  renderInteractor->Start();

}
