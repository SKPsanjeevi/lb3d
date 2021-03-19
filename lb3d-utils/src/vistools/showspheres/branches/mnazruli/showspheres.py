#!/usr/bin/env python
# -*- coding: <encoding name> -*-

# This program gives you how to visualize streanlines and the Geometry.
# Please read the manual which is given into the subversion before running this program from your system.

# Importing visualization toolkit from your system.

import vtk
import sys
import os
import math
import string
import time
from math import cos, sin,pi
from vtk.util.colors import tomato
import os.path
import spheredefaults as cfg
from array import array

# This function will help you to find the options for this program.
def printHelp():
 print """
 Usage: showspheres FilesAndOptions
 Files and options can be mixed. Options apply to all following files.
 E.g. you can change the text during the movie.
 Options (including the minus sign) can also appear in the data files.

 Options are:
 -color, -nocolor    Read color value between 0 and 1 from 5th column of input
  -write -nowrite     Do or do not write jpegs of the images displayed
  -fixedradius
  -nofixedradius      If nofixedradius is set, the radius will be read from
                      the 4th column of the data file. Otherwie use -radius=
  -radius=            Set radius if fixedradius is set.
  -anisotropy=        ratio of parallel and orthogonal half axes
  -elevation=         Observation angle,can be varied from 0 to 360 degrees
  -azimuth=           Obsevation angle, can be varied from 0 to +90/-90 degress
  -zoom=              To increase or decrease the image size.   
  -text=              Set the text that is displayed in the image
  -delay=             Delay in seconds between images
  
 """
 exit(0)




# This function is defined to print the spheredefaults and other process options in the command window

def printSettings():
 print "WriteImages:",writeImages
 print "color:",color
 print "fixedradius",fixedRadius
 print "radius:",radius
 print "anisotropy:",anisotropy
 print "delay:",delay
 print "elevation:", elevation
 print "azimuth:", azimuth
 print "zoom:", zoom

# This funsction is defined to change the elevation of the visualized image from the commandline

def setElevation(ea):
 ren.GetActiveCamera().Elevation(ea) 
 if first == 0:
  ren.GetActiveCamera().Modified()

# This funsction is defined to change the azimuth of the visualized image from the command line

def setAzimuth(aa):
 ren.GetActiveCamera().Azimuth(aa)
 if first == 0:
  ren.GetActiveCamera().Modified()

# This funsction is defined to increase or decrease of the visualized image from the command line

def setZoom(z):
 ren.GetActiveCamera().Zoom(z)
 if first == 0:
  ren.GetActiveCamera().Modified()



def setImagewrite(filename):
   toImage =vtk.vtkWindowToImageFilter()
   toImage.SetInput(renWin)
   toImage.Modified()
   writer = vtk.vtkJPEGWriter()
   writer.SetQuality(100)
   writer.ProgressiveOn()
   writer.SetInput(toImage.GetOutput())
   writer.SetFileName(filename)
   writer.Modified()
   print "writing to ",filename
   writer.Write()












  
# Rendering stuff for visualization, If you want to change the screen size, please just change the parameters into ren.Setsize(_, _)

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()	
renWin.AddRenderer(ren)
renWin.SetSize(600, 500)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# If you want to change the background of the screen please do it here by changing the numerical values of SetBackground (  ).

ren.SetBackground(155, 133, 122)

# Global function defined to transform the Sphere to Ellipsoid

sphereToEllipsoid = vtk.vtkTransform()

# Default value for sphere. To transform the Sphere to Ellipsid, the value of anisotropy should be less then 1 from the command line 

anisotropy = 1.0

# Function defined to prepare the text
def setText(txt):
 text.SetInput(txt)

#Taking Input fron Command line to change the elevation, azimuth, Writing Images and for ellipse visualization.

def processOption(s):
 global color
 global elevation
 global azimuth
 global zoom
 global writeImages
 global delay
 global radius
 global anisotropy
 global fixedRadius
 global geometryfile
 global sphereToEllipsoid
 global inputmode
 
 
 if s.find('-color') == 0:
  print "Setting color mode."
  color = 1
 
 if s.find('-nocolor') == 0:
  print "Unsetting color mode."
  color =0

 if s.find('-write') == 0:
  print "Writing images."
  writeImages = 1

 if s.find('-nowrite') == 0:
  print "Not writing images"
  writeImages =0

 if s.find('-fixedradius') == 0:
  print "Using fixed radius."
  fixedRadius = 1

 if s.find('-nofixedradius') == 0:
  print "Don't use fixed radius"
  fixedRadius =0
 if s.find('-radius') == 0:
  radius = s.split('=')[1]

 if s.find('-anisotropy') == 0:
  anisotropy = float(s.split('=')[1])
  sphereToEllipsoid.Scale(anisotropy,1,1)
  
 if s.find('-text') == 0:
  setText(s.split('=')[1])

 #if s.find('-geometry') == 0:
  #geometryfile = s.split('=')[1]

 if s.find('-elevation') == 0:
  setElevation(float(s.split('=')[1]))
   

 if s.find('-azimuth') == 0:
  setAzimuth(float(s.split('=')[1]))
  
 if s.find('-zoom') == 0:
  setZoom(float(s.split('=')[1]))

 if s.find('-delay') == 0:
  delay = float(s.split('=')[1])
 
 if s.find('-h') ==0:
  printHelp()

 if s.find('-help') ==0:
  printHelp()

    

#Reading Geometry File from Commamd Line in VTK format(Structured Points) to Generate Iso Surface and see the manual for details . 																																																																																																																																																						

data = vtk.vtkFloatArray()


#Reading position vector file either in .text or .asc format from the System for Streamlines or sphere visualization.
def readpoints(filename):
  
  file=open(filename,"r")
  xdim=32
  #int(raw_input("Please Enter the X-Dimension: "))
  ydim=64
  #int(raw_input("Please Enter the Y-Dimension: "))
  zdim=608
  #int(raw_input("Please Enter the Z-Dimension: "))
  lines=0
  a=0
  while file:
   line = file.readline()
   if line == "":
    break
   line = str.rstrip(line)
   cols=line.split('	') 
   if len(cols)  <3:
    if line.find(' ') == 0:
     print "Treating as options: ",line
     processOption(line)
    else: 
     print "Omitting line: ", line
     continue
   if len(cols)<11: 
    x=cols[0]
    y=cols[1]
    z=cols[2]
  
   else:
    x=cols[0]
    y=cols[1]
    z=cols[2]
    ox=cols[9]
    oy=cols[10]
    oz=cols[11]
    orientations.InsertNextTuple3(float(ox), float(oy), float(oz))
    if fixedRadius == 0:
      r=2*float(cols[3])
    else:
     if len(cols) <3:                     # was 4 (fj)
      print "Too few columns! Omitting: ",line
      continue

    r =2*float(radius)
   
    if color==1:
     c=cols[4]
     
    else:
     c=0
    orientations.InsertNextTuple3(float(ox), float(oy), float(oz))
    radii.InsertNextValue(float(r))
    tags.InsertNextValue(float(c))
   
   
    # some radii
    
  # p=float(z)
  # if m==0:
  #  zmin=p  
   
  # if m>=1:
  #   zmax=p
  #   p=(round(zmin)>=(zdim)) and (round(zmax)<=1)
  #   if p==1:
  #     d=m-1
       
  #     a=a+1
   #    if a==1: 
    #    t=d

    # zmin=zmax
    
   #orientations.InsertNextTuple3(float(ox), float(oy), float(oz)) 
   points.InsertNextPoint(float(x),float(y), float(z))
   
   lines=lines+1
  print filename, lines
  return lines	
# Putting the points from position vector file into Array
radii = vtk.vtkFloatArray()
#radii.SetName("radius")

# define the colours for the spheres
tags = vtk.vtkFloatArray()
#tags.SetName("tag")

data = vtk.vtkFloatArray()
data.SetNumberOfComponents(2)
data.SetName("data")

orientations = vtk.vtkFloatArray()
orientations.SetNumberOfComponents(3)
   
points=vtk.vtkPoints()


#Define Scalar and vectors
scalars=vtk.vtkFloatArray()
#scalars.SetNumberOfTuples(m)
 
vector=vtk.vtkFloatArray()
#vector.SetNumberOfComponents(3)
#vector.SetNumberOfTuples(m)

# Getting data for Unstructure Grid 
pdata=vtk.vtkUnstructuredGrid()     
pdata.SetPoints(points)
#pdata.GetPointData().AddArray(data)
pdata.GetPointData().SetActiveScalars("data")
#pdata.GetPointData().SetNormals(orientations)

# Visualizatizing Streamline as Streamtube
profileTubes = vtk.vtkTubeFilter()
profileTubes.SetNumberOfSides(4)
profileTubes.SetInput(pdata)
profileTubes.SetRadius(0.2)

 
 
#Mapping the data for Streamtube
map = vtk.vtkDataSetMapper()
map.SetInputConnection(profileTubes.GetOutputPort())

# Define actor for Streamtube rendering
triangulation = vtk.vtkActor()
triangulation.SetMapper(map)
triangulation.GetProperty().SetColor(1, 0, 0)

  

# Sphere visualization. If you to want chage the radius, resolution, please  do changes in the numeric parameter only. 
sphere = vtk.vtkSphereSource()
sphere.SetRadius(0.5)
sphere.SetPhiResolution(10)
sphere.SetThetaResolution(10)

# Transforming Sphere into ellipse  
sphereToEllipsoid = vtk.vtkTransform()
sphereToEllipsoid.Scale(anisotropy,1,1)

filter = vtk.vtkTransformPolyDataFilter()
filter.SetInputConnection(sphere.GetOutputPort())
filter.SetTransform(sphereToEllipsoid)


# Making the glyphs from UnstruncturedGrid
glyph = vtk.vtkGlyph3D()
glyph.SetInput(pdata)
glyph.SetSource(filter.GetOutput())
# glyph.SetSource(sphere.GetOutput())
glyph.ClampingOff()
glyph.SetScaleModeToScaleByScalar()
glyph.SetScaleFactor(1.0)
glyph.SetColorModeToColorByScalar()

glyph.OrientOn()
glyph.SetVectorModeToUseNormal()

# set up the mapper
mapper = vtk.vtkPolyDataMapper()
mapper.SetInput(glyph.GetOutput())
mapper.ScalarVisibilityOn()
mapper.ColorByArrayComponent("data", 1)


# set up the actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetColor(1,0,0)
 

ren.AddActor(actor)


# Prepare text
text = vtk.vtkTextMapper()
t =""
text.SetInput(t)
text.GetTextProperty().SetColor(0,0,1)
textActor = vtk.vtkScaledTextActor()
textActor.SetMapper(text)
textActor.SetDisplayPosition(cfg.textX,cfg.textY)
ren.AddActor2D(textActor)



#Reading the Geometry File from the command line. Please see the manual from subversion for details. 

def isoFile(filename):
 
 reader = vtk.vtkStructuredPointsReader()
 reader.SetFileName(filename)
 reader.Update()

 iso = vtk.vtkMarchingCubes()
 iso.SetInput(reader.GetOutput())
 iso.SetValue(0,0.4)

# Set up mapper 
 poly_map = vtk.vtkPolyDataMapper()
 poly_map.SetInput(iso.GetOutput())

# Set up actor for iso surface rendering
 isoActor = vtk.vtkActor()
 isoActor.SetMapper(poly_map) 

# Define outline of the isosurface
 outline = vtk.vtkOutlineFilter()
 outline.SetInputConnection(reader.GetOutputPort())
 mapOutline = vtk.vtkPolyDataMapper()
 mapOutline.SetInputConnection(outline.GetOutputPort())
 outlineActor = vtk.vtkActor()
 outlineActor.SetMapper(mapOutline)
 outlineActor.GetProperty().SetColor(0, 0, 0)

# Set up the Back ground of iso surface 
 ren.SetBackground(1,1,1)

# Adding actors for rendering. If you want to see the outline of the isosurface then  please activate outlineactor here.
 ren.AddActor(outlineActor)  
 ren.AddActor(isoActor)
 
# Selecting renderer good position and focal point. If you want  you can change the camera positin for viewing your expected image. 
 
 ren.GetActiveCamera().SetPosition(0, 2,20)
 ren.GetActiveCamera().SetFocalPoint(-15,12,90)
 ren.GetActiveCamera().SetViewUp(180,10 ,-90 )

 ren.ResetCamera()
 ren.GetActiveCamera().Dolly(2.25)
 ren.ResetCameraClippingRange()
	
 largestreamline=vtk.vtkRenderLargeImage()
 largestreamline.SetInput(ren)
 
 largestreamline.SetMagnification(2)
 
# Writing .PNG Images of the streamlines with geometry
 writer = vtk.vtkPNGWriter()
 writer.SetInputConnection(largestreamline.GetOutputPort())
 writer.SetFileName("Streamline.png")
 writer.Modified()
 
 writer.Write()

#Reading the Geometry File from the command line. Please see the manual from subversion for details.

def geoFile(filename):
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
# Getting Points for sphere visualization

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
 

# Cube rendering from Geometry file

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

# Mapping of sphere data and also adding actors for rendering the spheres.Here you can also change the color of shere changing the parameters into SetColor. You can put any values and observe changes but (0,0,0) is used for black and (1,1,1) is for white.

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
 

#Reading StructuredPoints for Streamline visualization. For details of header file for structure points in vtk file format,please see the manual

def streampoints(filename):
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

# Selecting renderer good position and focal point. Here You can change the camera positin for viewing a good image that you need. 

 ren.GetActiveCamera().SetPosition(0, 1,0)
 ren.GetActiveCamera().SetFocalPoint(0,0,0)
 ren.GetActiveCamera().SetViewUp(0,0 , 1)

 

 ren.GetActiveCamera().Dolly(1.4)
 ren.ResetCameraClippingRange()
	
 largestreamline=vtk.vtkRenderLargeImage()
 largestreamline.SetInput(ren)
 largestreamline.SetMagnification(9)
 
 #cam1 = ren.GetActiveCamera().Zoom(z)
 #cam1 = ren.GetActiveCamera().Elevation(2)
 #cam1 = ren.GetActiveCamera().Azimuth(-5)
 #ren.SetActiveCamera(cam1)
 ren.ResetCamera()
 # Writing .PNG Images of the streamlines with images.
 
 writer = vtk.vtkPNGWriter()
 writer.SetInputConnection(largestreamline.GetOutputPort())
 writer.SetFileName("Streamline.png")
 writer.Modified()
 writer.Write()


# Define the Global value of elevation, azimuth, zoom, writing image
azimuth =cfg.azimuth
elevation =cfg.elevation
zoom =cfg.zoom 

color =cfg.color
writeImages =cfg.writeImages
delay=cfg.delay
radius =cfg.radius
fixedRadius = cfg.fixedRadius
geometryfile = ""
printSettings()
first=1

	

#Main: Taking Input from the system as  input file

i=1
while sys.argv[i].find('-') ==0:
 processOption(sys.argv[i])
 i=i+1 
 print "Number of arguments:",str(len(sys.argv[i]))

a= sys.argv[i]

b=[".txt", ".asc",".dat"]

if os.path.splitext(a)[1] in b:
  
  lines=readpoints(sys.argv[i])
  data.SetNumberOfTuples(lines)
  data.CopyComponent(0, radii, 0)
  data.CopyComponent(1, tags, 0)
  
else:
  
  lines=geoFile(sys.argv[i])
 	
iren.Initialize()
renWin.Render()

setElevation(elevation)
setAzimuth(azimuth)
setZoom(zoom)

i=2
while i<len(sys.argv):
 if sys.argv[i].find('-') ==0:
  print "Commandline option: ",sys.argv[i]
  processOption(sys.argv[i])

 else: 
  d=sys.argv[i] 
  e=[".vtk"]
  if os.path.splitext(d)[1] in e:
 
   isoFile(sys.argv[i])
  
  else: 
   
   streampoints(sys.argv[i])
  if writeImages == 1:
   setImagewrite(sys.argv[i].replace('.vtk','')+'.jpg')
 time.sleep(delay)
 i=i+1
 
#iren.Start()	
