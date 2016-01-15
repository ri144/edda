// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <vtkPiecewiseFunction.h>
#include <vtkSmartPointer.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkImageData.h>
#include <vtkVolume.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkColorTransferFunction.h>
#include <vtkOutlineFilter.h>

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
#include "filters/uncertain_isocontour.h"

using namespace std;
using namespace edda;

#define vsp_new(cls,x) vtkSmartPointer<cls> x = vtkSmartPointer<cls>::New()

typedef dist::Gaussian<> Gaussian;

int main(int argc, char **argv) {
  cout << "isoProbField <info file> <iso-value>" << endl;
  if (argc<=2)
    return -1;
  string info_file = argv[1];
  float isov = atof(argv[2]);

  // load data
  shared_ptr<Dataset<Gaussianf> > dataset = loadData<Gaussianf>(info_file);

  // uncertain isocontour
  shared_ptr<Dataset<float> > output = uncertainIsocontour(dataset, isov);

  // Convert to vtk image data
  vsp_new(vtkImageData, image);
  image->SetDimensions(output->getDimension());
  image->AllocateScalars(VTK_FLOAT,1);

  shared_ary<float> array = boost::any_cast<shared_ary<float> >( output->getArray()->getRawArray() );
  std::copy( &array[0], &array[0]+output->getArray()->getLength(),
            (float *)image->GetScalarPointer(0,0,0));

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
  volumeMapper->SetInputData(image);

  vsp_new(vtkVolume, volume);
  volume->SetMapper(volumeMapper);
  volume->SetProperty(volumeProperty);

  // outline
  vsp_new(vtkOutlineFilter, outline);
  outline->SetInputData(image);
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
