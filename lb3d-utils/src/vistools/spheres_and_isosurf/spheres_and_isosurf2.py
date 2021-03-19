#!/usr/bin/env python
#Beware: this is badly written, not every parameter has any effect at all, some may also do something different than what they are suposed to.

import vtk
import sys
import math


class spheres_and_isosurf():
	def __init__(self,fid,new_params):
		self.params=self.ParamStorage()
		for param in new_params:
			self.params.processOption(param)
		self.params.set_id(fid)
		self.ren = vtk.vtkRenderer()
		self.renWin = vtk.vtkRenderWindow()
		self.renWin.AddRenderer(self.ren)
		self.renWin.SetSize(self.params.hsize, self.params.vsize)
		if self.params.spheres ==1:
			self.init_spheres()
		self.isoFile()
		if self.params.interactive ==1:
			self.iren = vtk.vtkRenderWindowInteractor()
			self.iren.SetRenderWindow(self.renWin)
		if self.params.volume ==1:
			self.volren()
		self.ren.SetBackground(1, 1, 1)
		self.ren.ResetCameraClippingRange()
		# Interface transparent:
		if self.params.writeImages ==1:
			self.ren.Render()
			self.dumpImage('left')
		self.isoActor.GetProperty().SetOpacity(1)
		if self.params.spheres ==1:
			self.sphereActor.GetProperty().SetOpacity(self.params.sphereopacity)
		self.ren.Render()
		# Interface solid
		if self.params.writeImages ==1:
			self.dumpImage('right')
		if self.params.interactive ==1:
			self.iren.Start()
	def update(self,fid):
		self.params.set_id(fid)
		self.isoreader.SetFileName(self.params.colour_file)
		self.isoreader.Update()
		self.volreader.SetFileName(self.params.colour_file)
		self.volreader.Update()
		#self.isoreader.GetOutput().SetOrigin(1,1,1)
		if self.params.spheres ==1:
			self.sphereGlyph.SetInput(self.readpoints())
		self.ren.ResetCameraClippingRange()
		self.isoActor.GetProperty().SetOpacity(self.params.isoopacity)
		if self.params.spheres ==1:
			self.sphereActor.GetProperty().SetOpacity(1)
		self.ren.Render()
		self.dumpImage('left')
		self.isoActor.GetProperty().SetOpacity(1)
		if self.params.spheres ==1:
			self.sphereActor.GetProperty().SetOpacity(self.params.sphereopacity)
		self.ren.Render()
		self.dumpImage('right')

	#Writing Image and Video File
	def dumpImage(self,name):
		filename='rendered_'+name+'_'+self.params.id+'.png'
		toImage =vtk.vtkWindowToImageFilter()
		toImage.SetInput(self.renWin)
		toImage.Modified()
		writer = vtk.vtkPNGWriter()
		writer.SetInput(toImage.GetOutput())
		writer.SetFileName(filename)
		writer.Modified()
		print "writing to ",filename	
		writer.Write()

	def volren(self):
		filename=self.params.colour_file
		self.volreader = vtk.vtkStructuredPointsReader()
		self.volreader.SetFileName(filename)
		self.volreader.Update()
		#self.volreader.GetOutput().SetOrigin(1,1,1)
		#volrange=self.volreader.GetOutput().GetScalarRange()
		volrange=(-1.0,1.0)
		imageShift=vtk.vtkImageShiftScale()
		imageShift.SetShift(-1*volrange[0])
		imageShift.SetScale(255.0/(volrange[1]-volrange[0]))
		imageShift.SetOutputScalarTypeToUnsignedChar()
		imageShift.SetInput(self.volreader.GetOutput())

		volume_mapper = vtk.vtkVolumeRayCastMapper()
		volume_mapper.SetInput(imageShift.GetOutput())

		composite_function = vtk.vtkVolumeRayCastCompositeFunction()
		volume_mapper.SetVolumeRayCastFunction(composite_function)

		color_transfer_func = vtk.vtkColorTransferFunction()
		color_transfer_func.AddRGBPoint( 0, 0.0, 0.0, 1.0 )
		color_transfer_func.AddRGBPoint( (255/2)-10, 0.0, 0.0, 1.0 )
		color_transfer_func.AddRGBPoint( (255/2)+10, 1.0, 0.0, 0.0 )
		color_transfer_func.AddRGBPoint( 255, 1.0, 0.0, 0.0 )

		opacity_transfer_func = vtk.vtkPiecewiseFunction()
		opacity_transfer_func.AddPoint( 0, self.params.volumeopacitylow )
		opacity_transfer_func.AddPoint( 118, self.params.volumeopacityhigh )
		opacity_transfer_func.AddPoint( 138, self.params.volumeopacityhigh )
		opacity_transfer_func.AddPoint( 255, self.params.volumeopacitylow )

		volume_properties = vtk.vtkVolumeProperty()
		volume_properties.SetColor( color_transfer_func )
		volume_properties.SetScalarOpacity( opacity_transfer_func )

		volume = vtk.vtkVolume()
		volume.SetMapper( volume_mapper )
		volume.SetProperty( volume_properties )

		self.ren.AddVolume( volume )

	def isoFile(self):
 		filename=self.params.colour_file
		self.isoreader = vtk.vtkStructuredPointsReader()
		self.isoreader.SetFileName(filename)
		self.isoreader.Update()
		#self.isoreader.GetOutput().SetOrigin(1,1,1)
		iso = vtk.vtkMarchingCubes()
		iso.SetInput(self.isoreader.GetOutput())
		iso.SetValue(0,self.params.isovalue)
		poly_map = vtk.vtkPolyDataMapper()
		poly_map.SetInput(iso.GetOutput())
		poly_map.ScalarVisibilityOff()
		self.isoActor = vtk.vtkActor()
		self.isoActor.SetMapper(poly_map)
		self.isoActor.GetProperty().SetColor(self.params.isocolour[0],self.params.isocolour[1],self.params.isocolour[2]) 
		self.isoActor.GetProperty().SetOpacity(self.params.isoopacity)

		if self.params.isosurf ==1:
			self.ren.AddActor(self.isoActor)
 
# Create outline around data
#
		if self.params.outline ==1:
			outline = vtk.vtkOutlineFilter()
			outline.SetInput(self.isoreader.GetOutput())
			outlineMapper = vtk.vtkPolyDataMapper()
			outlineMapper.SetInput( outline.GetOutput() )
 	
			outlineActor = vtk.vtkActor()
			outlineActor.SetMapper( outlineMapper )
			outlineActor.GetProperty().SetColor( 0, 0, 0 )

			self.ren.AddActor(outlineActor)
 
		cam1 = self.ren.GetActiveCamera()
		cam1.SetViewUp(1, 0, 0)
		cam1.SetFocalPoint(self.isoreader.GetOutput().GetCenter())
		size=[2*coord+1 for coord in self.isoreader.GetOutput().GetCenter()]
		cam1.SetPosition(size[0]/2,size[1]/2,8*math.sqrt(size[2]**2+size[1]**2+size[0]**2))
		#cam1.SetPosition(0, 0, 0)
		#cam1.Dolly(0.25)
 
		cam1.Elevation(self.params.elevation)
		cam1.Azimuth(self.params.azimuth)
		cam1.Zoom(self.params.zoom) 




	#Reading Text file from the System for plotting Spheres
	def readpoints(self):
		r=self.params.r
		filename=self.params.md_file
		points=vtk.vtkPoints()
		file=open(filename,"r")
		#xdim=int(raw_input("Please Enter the X-Dimension: "))
		#ydim=int(raw_input("Please Enter the Y-Dimension: "))
		#zdim=int(raw_input("Please Enter the Z-Dimension: "))
		m=0
		#a=0
		while file:
			line = file.readline()
			if line == "":
				break
			line = str.rstrip(line)
			cols =line.split(' ')
  
			if len(cols)  <3:
				print "Omitting line: ", line
				continue
    
			x=cols[2]
			y=cols[1]
			z=cols[0]
  

 						
			points.InsertNextPoint(float(x)-1,float(y)-1, float(z)-1)
			#points.InsertNextPoint(float(x),float(y), float(z))
			m+=1

 #Copied Sphere code from Finalstreamtracermod  
		spheredata=vtk.vtkPolyData()
		spheredata.SetPoints(points)

		vector=vtk.vtkFloatArray() 
		vector.SetNumberOfComponents(3)
		vector.SetNumberOfTuples(m)

		pline=vtk.vtkCellArray()
   
		pline.InsertNextCell(1)
		for i in range(0,1):
			pline.InsertCellPoint(i)
		spheredata.SetPolys(pline)
		spheredata.GetPointData().SetVectors(vector)

		return spheredata
	# Sphere visualization. 
 
	def init_spheres(self):
		spheres=vtk.vtkSphereSource()
		spheres.SetRadius(self.params.r)
		spheres.SetThetaResolution(32)
		spheres.SetPhiResolution(32)

		spheredata=self.readpoints()

		self.sphereGlyph=vtk.vtkGlyph3D()
		self.sphereGlyph.SetSource(spheres.GetOutput())
		self.sphereGlyph.SetInput(spheredata)

# Mapping of sphere data and also adding actors for rendering the spheres.Here you can also chnage the color of shere changing the parameters into SetColor. You can put any values and observe changes but (0,0,0) is used for black and (1,1,1) is for white.

		sphereMapper = vtk.vtkPolyDataMapper()
		sphereMapper.SetInput(self.sphereGlyph.GetOutput())
		sphereMapper.GlobalImmediateModeRenderingOn()
		self.sphereActor = vtk.vtkLODActor()
		self.sphereActor.SetMapper(sphereMapper)
		self.sphereActor.GetProperty().SetColor(self.params.spherecolour[0],self.params.spherecolour[1] ,self.params.spherecolour[2] ) 
		self.sphereActor.GetProperty().SetOpacity(1)

		self.ren.AddActor(self.sphereActor)

	class ParamStorage():
		def __init__(self):
			self.elevation=10.0
			self.zoom=4.0
			self.azimuth=-10.0
			self.writeImages=0
			self.r=2.0
			self.hsize=800
			self.vsize=600
			self.md_file=""
			self.colour_file=""
			self.id=""
			self.interactive=0
			self.spheres=1
			self.outline=1
			self.isosurf=1
			self.isovalue=0.0
			self.sphereopacity=0.1
			self.isoopacity=0.4
			self.volume=1
			self.volumeopacitylow=0.01
			self.volumeopacityhigh=0.01
                        self.spherecolour=[0,1,0]
                        self.isocolour=[0.8, 0.2, 0.9]
			self.interactive=0
		def processOption(self,s):
			if s.startswith('-elevation='):
				self.elevation = float(s.split('=')[1])
				print "Elevation:", self.elevation
			if s.startswith('-zoom='):
				self.zoom = float(s.split('=')[1])
				print "Zoom:", self.zoom
			if s.startswith('-azimuth='):
				self.azimuth = float(s.split('=')[1])
				print "Azimuth:", self.azimuth
			if s.startswith('-radius='):
				self.r = float(s.split('=')[1])
				print "Radius:", self.r
			if s.startswith('-novolume'):
				print "No volume rendering"
				self.volume=0
			if s.startswith('-volumeopacitylow='):
				self.volumeopacitylow = float(s.split('=')[1])
				print "Volume Opacity Low:", self.volumeopacitylow
			if s.startswith('-volumeopacityhigh='):
				self.volumeopacityhigh = float(s.split('=')[1])
				print "Volume Opacity High:", self.volumeopacityhigh
			if s.startswith('-sphereopacity='):
				self.sphereopacity = float(s.split('=')[1])
				print "Sphere Opacity:", self.sphereopacity
			if s.startswith('-isoopacity='):
				self.isoopacity = float(s.split('=')[1])
				print "Surface Opacity:", self.isoopacity
			if s.startswith('-isosurfacevalue='):
				self.isovalue = float(s.split('=')[1])
				print "Isosurface Value:", self.isovalue		
			if s.startswith('-size='):
				self.hsize,self.vsize = [int(size) for size in s.split('=')[1].split('x')]
				print "Size:", self.hsize, '*', self.vsize
			if s.startswith('-id='):
				if self.set_id(s.split('=')[1]):
					print "Colourfile:",self.colour_file,"MD-File:",self.md_file
			if s.startswith('-write'):
				print "Writing images."
				self.writeImages = 1
			if s.startswith('-interactive'):
				print "Interactive session"
				self.interactive=1
			if s.startswith('-nospheres'):
				print "No spheres"
				self.spheres=0
			if s.startswith('-nooutline'):
				print "No outline"
				self.outline=0
			if s.startswith('-nosurface'):
				print "No IsoSurface"
				self.isosurf=0
			if s.startswith('-isocolour='):
				self.isocolour = [float(c) for c in s.split('=')[1].split(',')]
				print 'Iso colour:',self.isocolour
			if s.startswith('-spherecolour='):
				self.spherecolour = [float(c) for c in s.split('=')[1].split(',')]
				print 'Sphere colour:',self.spherecolour

		def set_id(self, fid):
			cfn="colour_"+fid+".vtk"
			mdfn="md-cfg_"+fid+".asc"

			#print cfn,mdfn
			try:
				if self.spheres ==1:
					mdf=open(mdfn,"r")
					mdf.close()
				cf=open(cfn,"r")
				cf.close()
			except IOError:
				print "IOError while checking ID"
				return False
			else:
				self.id=fid
				self.colour_file=cfn
				self.md_file=mdfn
				return True

if __name__ == "__main__":
	import sys
	spheres_and_isosurf(sys.argv[1],sys.argv[2:])

