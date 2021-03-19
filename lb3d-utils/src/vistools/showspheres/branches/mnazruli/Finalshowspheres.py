#!/usr/bin/env python
# -*- coding: <encoding name> -*-
import vtk
import sys
import os
import math
import string
from vtk.util.colors import tomato
import spheredefaults as cfg
from array import array

def processOption(s):
 global elevation
 global writeImages

 if s.find('-elevation') == 0:
  elevation = s.split('=')[1]
 if s.find('-write') == 0:
  print "Writing images."
  writeImages = 1
 if s.find('-nowrite') == 0:
  print "Not writing images"
  writeImages =0
#Writing Image and Video File
def dumpImage(filename):
   toImage =vtk.vtkWindowToImageFilter()
   toImage.SetInput(renWin)
   toImage.Modified()
   writer = vtk.vtkPNGWriter()
   writer.SetInput(toImage.GetOutput())
   writer.SetFileName(filename)
   writer.Modified()
   print "writing to ",filename	
   writer.Write()
   

   Fractal0 = vtk.vtkImageMandelbrotSource()
   Fractal0.SetWholeExtent( 0, 247, 0, 247, 0, 0 )
   Fractal0.SetProjectionAxes( 0, 1, 2 )
   Fractal0.SetOriginCX( -1.75, -1.25, 0, 0 )
   Fractal0.SetSizeCX( 2.5, 2.5, 2, 1.5 )
   Fractal0.SetMaximumNumberOfIterations(1000) 

   cast = vtk.vtkImageCast()
   cast.SetInputConnection(Fractal0.GetOutputPort())
   cast.SetOutputScalarTypeToUnsignedChar()

   table = vtk.vtkLookupTable()
   table.SetTableRange(0, 100)
   table.SetNumberOfColors(100)
   table.Build()
   table.SetTableValue(99, 0, 0, 0,0)

   colorize = vtk.vtkImageMapToColors()
   colorize.SetOutputFormatToRGB()
   colorize.SetLookupTable(table)
   colorize.SetInputConnection(cast.GetOutputPort())


   mpeg=vtk.vtkMPEG2Writer()
   #mpeg.SetInputConnection(toImage.GetOutputPort())
   #mpeg.SetFileName("TestMovie.mpg")
   #print "writing File TestMovie.mpg "
   #mpeg.Start()

   itr=1
   for itr in range(1,10):
     Fractal0.SetMaximumNumberOfIterations(itr)
     table.SetTableRange(0,itr)   
     table.SetNumberOfColors(itr) 
     table.ForceBuild()
     table.SetTableValue(itr-1,0,0,0,0)   
     mpeg.Write()
     itr+=1
   err=0
   exists=0
 
     

#ISO Surface

def isoFile(filename):
 
 reader = vtk.vtkStructuredPointsReader()
 reader.SetFileName(filename)
 reader.Update()

 iso = vtk.vtkMarchingCubes()
 iso.SetInput(reader.GetOutput())
 iso.SetValue(0,0.4)
 poly_map = vtk.vtkPolyDataMapper()
 poly_map.SetInput(iso.GetOutput())

 isoActor = vtk.vtkActor()
 isoActor.SetMapper(poly_map) 
 
 ren.SetBackground(1,1,1)
 
 elevation=vtk.vtkElevationFilter()
 elevation.SetInputConnection(iso.GetOutputPort())
 elevation.ReleaseDataFlagOn()

 normals = vtk.vtkPolyDataNormals()
 normals.SetInput(elevation.GetPolyDataOutput())
 normals.SetFeatureAngle(60)
 normals.ConsistencyOff()
 normals.SplittingOff()
 normals.ReleaseDataFlagOn()

 demMapper = vtk.vtkPolyDataMapper()
 demMapper.SetInputConnection(normals.GetOutputPort())
 demMapper.ImmediateModeRenderingOn()

 demActor = vtk.vtkLODActor()
 demActor.SetMapper(demMapper)
 #ren.AddActor(demActor)
 
  
 ren.AddActor(isoActor)
 
 
 cam1 = ren.GetActiveCamera()
 cam1.SetViewUp(0, 0, 1)
 cam1.SetFocalPoint(reader.GetOutput().GetCenter())
 cam1.SetPosition(0, 0, 0)

 ren.ResetCamera()
#Changing Elevation, Azimuth,and Zoom 
 elevation=input("Please Enter the Value of Elevation:")
 azimuth=input("Please Enter the Value of Azimuth:")
 zoom =input("Please Enter the Value of Zoom:")
 
 cam1.Elevation(elevation)
 cam1.Azimuth(azimuth)
 cam1.Zoom(zoom) 


 #Reading Text file from the System for plotting Streamlines
def readpoints(filename):
 
  points=vtk.vtkPoints()
  file=open(filename,"r")
  xdim=int(raw_input("Please Enter the X-Dimension: "))
  ydim=int(raw_input("Please Enter the Y-Dimension: "))
  zdim=int(raw_input("Please Enter the Z-Dimension: "))
  m=0
  a=0
  while file:
   line = file.readline()
   if line == "":
    break
   line = str.rstrip(line)
   cols =line.split('	')
  
   if len(cols)  <3:
    if line.find('-') == 0:
     print "Treating as options: ",line
     processOption(line)
    else: 
     print "Omitting line: ", line
     continue
    
   x=cols[1]
   y=cols[2]
   z=cols[3]
  
   p=float(z)
   if m==0:
    zmin=p  
   
   if m>=1:
     zmax=p
     p=(round(zmin)>=(zdim-2)) and (round(zmax)<=1)
     if p==1:
       d=m-1
       a=a+1
       if a==1: 
        t=d

     zmin=zmax
    
 
   points.InsertNextPoint(float(x),float(y), float(z))
   #print m
   m=m+1


 #Define Scalar and vectors
  scalars=vtk.vtkFloatArray()
  scalars.SetNumberOfTuples(m)
 
  vector=vtk.vtkFloatArray()
  vector.SetNumberOfComponents(3)
  vector.SetNumberOfTuples(m)
 
  pdata=vtk.vtkPolyData()
  pline=vtk.vtkCellArray()
   
  pline.InsertNextCell(t)
  for i in range(0,t):
   pline.InsertCellPoint(i)

    
  pline.InsertNextCell(d-t-1)
  for i in range(t+1,d):
   pline.InsertCellPoint(i)
   
  pline.InsertNextCell(m-d-1)
  for i in range(d+1,m):
   pline.InsertCellPoint(i)
   
    
   
   
  pdata.SetPoints(points)
  pdata.SetLines(pline)
  pdata.GetPointData().SetVectors(vector)
 
  
 
  profileTubes = vtk.vtkTubeFilter()
  profileTubes.SetNumberOfSides(4)
  profileTubes.SetInput(pdata)
  profileTubes.SetRadius(0.75)

 
  slineGrid=vtk.vtkStructuredGrid()
  #slineGrid.Allocate(10,10)
  #slineGrid.InsertNextCell(12,Idlist)
  slineGrid.GetPointData().SetScalars(scalars)
  slineGrid.GetPointData().SetVectors(vector)
  slineGrid.GetPointData().AddArray
  slineGrid.SetPoints(points)
 
 
 
 
  map = vtk.vtkDataSetMapper()
  map.SetInputConnection(profileTubes.GetOutputPort())

  triangulation = vtk.vtkActor()
  triangulation.SetMapper(map)
  triangulation.GetProperty().SetColor(1, 0, 0)
  
  ren.AddActor(triangulation)
  ren.SetBackground(1, 1, 1)
 
  
 

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()	
renWin.AddRenderer(ren)
renWin.SetSize(800, 300)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)


 
writeImages =cfg.writeImages
 	

#Main: Taking Input from the system as a file

i=1
while sys.argv[i].find('-') ==0:
 processOption(sys.argv[i])
lines = readpoints(sys.argv[i])



i=2
while i<len(sys.argv):
 if sys.argv[i].find('-') ==0:
  print "Commandline option: ",sys.argv[i]
  processOption(sys.argv[i])
 else: 
  isoFile(sys.argv[i])
  
  if writeImages == 1:
   dumpImage(sys.argv[i].replace('.vtk','')+'.PNG')
  
 i=i+1
 
iren.Start()	
