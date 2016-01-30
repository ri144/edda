// Reference: VTK/Examples/VisualizationAlgorithms/Python/StreamlinesWithLineWidget.py

// Compatibility with Paraview
//#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
//#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL) 
// Compatibility with Paraview

#include <stdio.h>
#include <stdlib.h>

#include <list>
#include <iterator>

#include "vtkDataSet.h"
#include "vtkSmartPointer.h"
#include "vtkLineSource.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkOutlineFilter.h"
#include "vtkProperty.h"
#include <vtkXMLImageDataReader.h>
#include <vtkInterpolatedVelocityField.h>
#include <vtkImageData.h>
#include <vtkLineWidget.h>
#include <vtkCallbackCommand.h>
// streamline
#include "vtkStreamLine.h"

#include "vtkEddaReader.h"
#include "vtkSamplingArray.h"

using namespace std;


vtkLineWidget *lineWidget;
vtkRenderWindow *renWin;
vtkPolyData *seeds ;

vtkSmartPointer<vtkDataSet> getData(char *filename)
{
  vtkEddaReader *reader = vtkEddaReader::New();
  reader->SetFileName(filename);
  reader->Update();

  return reader->GetOutput();
}

void computeStreamlines(vtkObject* caller, unsigned long eventId, void *clientdata, void *calldata)
{
	double *point1 = lineWidget->GetPoint1();
	double *point2 = lineWidget->GetPoint2();
	printf("LineWidget Point1: %lf %lf %lf, Point2: %lf %lf %lf\n", point1[0], point1[1], point1[2], point2[0], point2[1], point2[2]);

	lineWidget->GetPolyData(seeds);
	renWin->Render();
}

int main(int argc, char **argv)
{
  // choose one of the following:
  //vtkSmartPointer<vtkDataSet> data = getData_plot3d();
  //vtkSmartPointer<vtkDataSet> data = getData_vts();
  vtkSmartPointer<vtkDataSet> data = getData(argv[1]);

  // rake
  lineWidget = vtkLineWidget::New();
  lineWidget->SetInputData(data);
  lineWidget->SetResolution(21); // 22 seeds along the line
  lineWidget->SetAlignToYAxis();
  lineWidget->PlaceWidget();
  lineWidget->ClampToBoundsOn();
  seeds = vtkPolyData::New();
  lineWidget->GetPolyData(seeds);


  // use vtkStreamLine
  vtkStreamLine *streamer = vtkStreamLine::New();
  streamer->SetInputData(data);
  streamer->SetSourceData(seeds);
  streamer->SetMaximumPropagationTime(10);
  streamer->SetIntegrationStepLength(.5);
  streamer->SetStepLength(.05);
  streamer->SetIntegrationDirectionToIntegrateBothDirections();
  streamer->VorticityOff();

  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
  mapper->SetInputConnection(streamer->GetOutputPort());
  mapper->SetScalarRange(data->GetScalarRange());
  vtkActor *actor = vtkActor::New();
  actor->SetMapper(mapper);


  //
  // outline
  //
  vtkOutlineFilter *outline = vtkOutlineFilter::New();
  outline->SetInputData(data);

  vtkPolyDataMapper *outlineMapper = vtkPolyDataMapper::New();
  outlineMapper->SetInputConnection(outline->GetOutputPort());

  vtkActor *outlineActor = vtkActor::New();
  outlineActor->SetMapper(outlineMapper);
  outlineActor->GetProperty()->SetColor(0,0,0);

  //
  // renderer
  //
  vtkRenderer *ren = vtkRenderer::New();
  renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren);
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
  iren->SetInteractorStyle(style);

  // line widget interactor
  lineWidget->SetInteractor(iren);
  lineWidget->SetDefaultRenderer(ren);
  vtkCallbackCommand *callback = vtkCallbackCommand::New();
  callback->SetCallback(computeStreamlines);
  lineWidget->AddObserver(vtkCommand::EndInteractionEvent, callback);

  ren->AddActor(actor);
  ren->AddActor(outlineActor);
  ren->SetBackground(.5,.5,.5);

  renWin->SetSize(500,500);

  iren->Initialize();
  renWin->Render();
  iren->Start();

  return 0;
}



