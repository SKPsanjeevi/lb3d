#!/usr/bin/env python

# This program gives you how to visualize streanlines with Geometry.
# Please read the manual which is given into the subversion before running this program from your system.

# Importing visualization toolkit from you your system.

import vtk,sys
from math import cos, sin,pi
from vtk.util.colors import *

# Rendering stuff for visualization, If you want to change the screen size, please just change the parameters into ren.Setsize(_, _)

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()	
renWin.AddRenderer(ren)
#renWin.SetSize(210, 190)
renWin.SetSize(210, 200)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# If you want to change the background of the screen please do it here.
ren.SetBackground(155, 133, 122)
#ren.GetActiveCamera().Zoom(3)
#ren.SetActiveCamera(cam1)

#Reading the Geometry File from the command line. Please see manual for details.  

def isoFile(filename):
 lines=0
 geomfile=open(filename,"rb")
 m=1
 spherepoints= vtk.vtkPoints()
 radii=vtk.vtkFloatArray()
 cubepoints=vtk.vtkPoints()
 while geomfile:
   line = geomfile.readline()
   if line == "":
    break
   line = str.rstrip(line)
   cols =line.split('	')
   
   x=cols[0]    
   y=cols[1]
   z=cols[2]
   r=cols[3]
   if m==1:
     sr=r
    
     spherepoints.InsertNextPoint(float(x),float(y), float(z))
   radii.InsertNextValue(float(r))
   if m>1:
     cubepoints.InsertNextPoint(float(x),float(y),float(z))
    
   m+=1  

# Putting the data into cellarray
 
 cubecell=vtk.vtkCellArray()
 cubecell.InsertNextCell(m-1)

 for j in range(0,m-1):
  cubecell.InsertCellPoint(j)
# Getting Points for spherevisualization

 spheredata=vtk.vtkPolyData()
 spheredata.SetPoints(spherepoints)

 vector=vtk.vtkFloatArray() 
 vector.SetNumberOfComponents(3)
 vector.SetNumberOfTuples(m)
 
 pline=vtk.vtkCellArray()
   
 pline.InsertNextCell(1)
 for i in range(0,1):
  pline.InsertCellPoint(i)
 spheredata.SetPolys(pline)
 spheredata.GetPointData().SetVectors(vector)
 

# Cube Rendering from Geometry file

 CubeData = vtk.vtkPolyData()
 CubeData.SetPoints(cubepoints)
 CubeData.GetPointData().SetVectors(vector)
# CubeSource is using as a source of vtkGlyph3D to visualize cube 

 Cube=vtk.vtkCubeSource()
 Cubesurface=vtk.vtkGlyph3D()
 Cubesurface.SetSource(Cube.GetOutput())
 Cubesurface.SetInput(CubeData)

# Data mapping for cube visualization

 CubeMapper = vtk.vtkPolyDataMapper()
 CubeMapper.SetInput(Cubesurface.GetOutput())
 CubeMapper.GlobalImmediateModeRenderingOn()
 CubeActor = vtk.vtkLODActor()
 CubeActor.SetMapper(CubeMapper)
 #CubeActor.GetProperty().SetColor(1, 0, 0)
 CubeActor.GetProperty().SetDiffuseColor(1, 0, 0);
 #CubeActor.GetProperty().SetSpecular(.3);
 #CubeActor.GetProperty().SetSpecularPower(20);
 #CubeActor.GetProperty().SetOpacity(0.5)
 #Sphere Rendering

# Sphere visualization. If you to want chage the radius, resolution, please  do changes in the numeric parameter only. 
 
 spheres=vtk.vtkSphereSource()
 spheres.SetRadius(16.0)
 spheres.SetThetaResolution(32)
 spheres.SetPhiResolution(32)

 Glyph=vtk.vtkGlyph3D()
 Glyph.SetSource(spheres.GetOutput())
 Glyph.SetInput(spheredata)

# Mapping of sphere data and also adding actors for rendering the spheres.Here you can also chnage the color of shere changing the parameters into SetColor. You can put any values and observe changes but (0,0,0) is used for black and (1,1,1) is for white.

 sphereMapper = vtk.vtkPolyDataMapper()
 sphereMapper.SetInput(Glyph.GetOutput())
 sphereMapper.GlobalImmediateModeRenderingOn()
 sphereActor = vtk.vtkLODActor()
 sphereActor.SetMapper(sphereMapper)
 sphereActor.GetProperty().SetColor(1, 0, 1) 

 ren.AddActor(CubeActor)
 ren.AddActor(sphereActor)
 

# Selecting the starting point for intergration of the vtkRungeKutta 4. Multiples sources are used here. if you want you can also use the single source for the starting point of integration if needed. In this case just remove the vtkLinesource keeping only onle of them. You can chage the number of Streamlines by changing the resolution into Setresolution.

line0=vtk.vtkLineSource()
line0.SetPoint1(0,0,20)
line0.SetPoint2(128,128,20)
line0.SetResolution(50)

line1=vtk.vtkLineSource()
line1.SetPoint1(50,0,5)
line1.SetPoint2(50,191,15)
line1.SetResolution(35)

line2=vtk.vtkLineSource()
line2.SetPoint1(191,0,15)
line2.SetPoint2(0,191,5)
line2.SetResolution(50)
 

#Reading StructuredPoints for Streamline visualization. For details of header file structure of vtk file please see the manual

def readpoints(filename):
 reader = vtk.vtkStructuredPointsReader()
 reader.SetFileName(filename)
 reader.Update()

 
 outline=vtk.vtkOutlineFilter()
 outline.SetInputConnection(reader.GetOutputPort())

 lineMapper = vtk.vtkDataSetMapper()
 lineMapper.SetInputConnection(outline.GetOutputPort())
 lineActor = vtk.vtkActor()
 lineActor.SetMapper(lineMapper)


 pdata=vtk.vtkPolyData()
 points=vtk.vtkPoints()   
 #Selecting source for Streamtracer 
 ns=200
 r=100
 Z=25
 for i in range(0, ns):
  a=(2*i*pi)/ns
  X=r*cos(a)
  Y=r*sin(a)
 # for Z in range(0,256):
  points.InsertNextPoint(float(X),float(Y),float(Z))

 pdata.SetPoints(points)
 
 integ0 = vtk.vtkRungeKutta4()
 integ1 = vtk.vtkRungeKutta4()
 integ2 = vtk.vtkRungeKutta4()

# Defining controller for vtkDistributedStreamTracer

 controller=vtk.vtkMPIController()
 Stream0 = vtk.vtkDistributedStreamTracer() 
 #Stream0 = vtk.vtkStreamTracer() 
 Stream0.SetInputConnection(reader.GetOutputPort())
 Stream0.SetSource(line0.GetOutput())
 #Stream0.SetSource(pdata)

# Setting the parameters for Integration. Here you can change the integration parameters according to your desire.
  
 Stream0.SetMaximumPropagation(255)
 Stream0.SetController(controller)
 Stream0.SetInitialIntegrationStepUnitToCellLengthUnit()
 Stream0.SetMaximumIntegrationStep(2000)
 Stream0.SetInitialIntegrationStep(0.5)
 Stream0.SetIntegrator(integ0)
 Stream0.SetIntegrationDirectionToBoth()

 Stream0.Update()

# Visualizing streamline as a tube.Here you can change the radius of the streamtube by changing the values of radius.
 
 streamTube0 = vtk.vtkTubeFilter()
 streamTube0.SetInputConnection(Stream0.GetOutputPort())
 streamTube0.SetRadius(0.25)
 streamTube0.SetNumberOfSides(16)
 #streamTube0.SetVaryRadiusToVaryRadiusByVector()
 mapStreamTube0 = vtk.vtkPolyDataMapper()
 mapStreamTube0.SetInputConnection(streamTube0.GetOutputPort())
 #mapStreamTube.SetScalarRange(reader.GetOutput())
 streamTubeActor0 = vtk.vtkActor()
 streamTubeActor0.SetMapper(mapStreamTube0)
 streamTubeActor0.GetProperty().BackfaceCullingOn()
 streamTubeActor0.GetProperty().SetColor(0, 1, 1) 

 Stream1 = vtk.vtkDistributedStreamTracer() 
 #Stream0 = vtk.vtkStreamTracer() 
 Stream1.SetInputConnection(reader.GetOutputPort())
 Stream1.SetSource(line1.GetOutput())
 #Stream0.SetSource(pdata)
 
 Stream1.SetMaximumPropagation(255)
 Stream1.SetController(controller)
 Stream1.SetInitialIntegrationStepUnitToCellLengthUnit()
 Stream1.SetMaximumIntegrationStep(2000)
 Stream1.SetInitialIntegrationStep(0.5)
 Stream1.SetIntegrator(integ1)
 Stream1.SetIntegrationDirectionToBoth()

 Stream1.Update()


 streamTube1= vtk.vtkTubeFilter()
 streamTube1.SetInputConnection(Stream0.GetOutputPort())
 streamTube1.SetRadius(0.25)
 streamTube1.SetNumberOfSides(12)
 #streamTube1.SetVaryRadiusToVaryRadiusByVector()
 mapStreamTube1 = vtk.vtkPolyDataMapper()
 mapStreamTube1.SetInputConnection(streamTube1.GetOutputPort())
 #mapStreamTube.SetScalarRange(reader.GetOutput())
 streamTubeActor1 = vtk.vtkActor()
 streamTubeActor1.SetMapper(mapStreamTube1)
 streamTubeActor1.GetProperty().BackfaceCullingOn()
 streamTubeActor1.GetProperty().SetColor(0.5, 0.25, 1) 
 #ren.AddActor(lineActor)

 Stream2 = vtk.vtkDistributedStreamTracer() 
 #Stream0 = vtk.vtkStreamTracer() 
 Stream2.SetInputConnection(reader.GetOutputPort())
 Stream2.SetSource(line2.GetOutput())
 #Stream2.SetSource(line1.GetOutput())
 #Stream2.SetSource(line0.GetOutput())
 
#Stream0.SetSource(pdata)
 
 Stream2.SetMaximumPropagation(255)
 Stream2.SetController(controller)
 Stream2.SetInitialIntegrationStepUnitToCellLengthUnit()
 Stream2.SetMaximumIntegrationStep(2000)
 Stream2.SetInitialIntegrationStep(0.5)
 Stream2.SetIntegrator(integ2)
 Stream2.SetIntegrationDirectionToBoth()

 Stream2.Update()


 streamTube2= vtk.vtkTubeFilter()
 streamTube2.SetInputConnection(Stream0.GetOutputPort())
 streamTube2.SetRadius(0.25)
 streamTube2.SetNumberOfSides(12)
 #streamTube2.SetVaryRadiusToVaryRadiusByVector()
 mapStreamTube2 = vtk.vtkPolyDataMapper()
 mapStreamTube2.SetInputConnection(streamTube2.GetOutputPort())
 #mapStreamTube.SetScalarRange(reader.GetOutput())
 streamTubeActor2 = vtk.vtkActor()
 streamTubeActor2.SetMapper(mapStreamTube2)
 streamTubeActor2.GetProperty().BackfaceCullingOn()
 streamTubeActor2.GetProperty().SetColor(0, .025, 0.125) 

 # ren.AddActor(lineActor)
# Adding the streanline to the actor for rendering

 ren.AddActor(streamTubeActor0)
 ren.AddActor(streamTubeActor1)
 #ren.AddActor(streamTubeActor2)

 renWin.Render()
 ren.ResetCamera()

# Selecting renderer good position and focal point. If you want  you can change the camera positin for viewing a good image that you need. 

 ren.GetActiveCamera().SetPosition(0, 1,0)
 ren.GetActiveCamera().SetFocalPoint(0,0,0)
 #ren.GetActiveCamera().SetViewUp(0,0 , 1)

 

 ren.GetActiveCamera().Dolly(1.4)
 ren.ResetCameraClippingRange()
	
 largestreamline=vtk.vtkRenderLargeImage()
 largestreamline.SetInput(ren)
 #largestreamline.SetMagnification(9)
 ren.ResetCamera()
 #cam1 = ren.GetActiveCamera().Zoom(1.5)
 #cam1 = ren.GetActiveCamera().Elevation(2)
 #cam1 = ren.GetActiveCamera().Azimuth(-5)
 #ren.SetActiveCamera(cam1)

 # Writing .PNG Images of the streamlines with images.
 
 writer = vtk.vtkPNGWriter()
 writer.SetInputConnection(largestreamline.GetOutputPort())
 writer.SetFileName("Streamline.png")
 writer.Modified()
 writer.Write()




#Main: Reading Files from Command line. Please see the manual for details.

i=1
while sys.argv[i].find('-') ==0:

  print "Geometry" 
lines=isoFile(sys.argv[i])


i=2
while i<len(sys.argv):
 if sys.argv[i].find('-') ==0:
  print "Commandline option: ",sys.argv[i]
 else: 
  readpoints(sys.argv[i])
    
 i=i+1

#iren.Start()	

