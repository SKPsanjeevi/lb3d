#!/usr/bin/env python
# This is based on spheres_and_isosurf2.py, but heavily restructured

import vtk
import sys
import math
import array

class omnirender():
	def __init__(self,fid,new_params):
		# Initialize parameters
		self.params = self.ParamStorage()
		print "OMNIRENDER initializing parameters"
		for param in new_params:
			self.params.processOption(param)
		if ( self.params.interactive == 0 and self.params.writeimage == 0 ):
			print "Program will do nothing; specify either -interactive or -writeimage"
			sys.exit()
		self.params.set_id(fid)

		# Set up renderer and window
		self.ren = vtk.vtkRenderer()
		self.ren.SetBackground(1, 1, 1)

		# Set up isosurfaces
		if (self.params.colour == 1):
			self.addIsoSurface()
		if (self.params.rock == 1):
			self.addRockIsoSurface()

		# Set up volume renderings
		if (self.params.colour ==1 and self.params.volume == 1):
			self.addVolume()
		if (self.params.rock == 1 and self.params.rockvolume == 1):
			self.addRockVolume()

		# Set up outline
		if (self.params.outline == 1):
			self.addOutline()

		# Set up spheres
		if (self.params.spheres == 1):
			self.addSpheres()

		# Set up ellipsoids
		if (self.params.ellipsoids == 1):
			self.addEllipsoids()

		# Set up ellipsoids
		#if (self.params.ellipsoidellipsoid == 1):
			#self.addEllipsoid()

		# Set up camera angles
		self.setCamera()

		print "Creating render window..."

		self.renwin = vtk.vtkRenderWindow()
		print "  Render library: " + self.renwin.GetRenderLibrary()

		if (self.params.offscreenrendering == 1):
			if (self.params.writeimage == 1):
				print "  Setting off-screen rendering"
				self.renwin.OffScreenRenderingOn()

		self.renwin.AddRenderer(self.ren)
		self.renwin.SetSize(self.params.hsize, self.params.vsize)
		self.ren.ResetCameraClippingRange()

		# Writing image to file
		if (self.params.writeimage == 1):
			print "RENDERING..."
			self.ren.Render()
			print " Done rendering."
			self.dumpImage()

		# Starting interactive session
		if (self.params.interactive == 1):
			print "Starting interactive session..."
			self.iren = vtk.vtkRenderWindowInteractor()
			self.iren.SetRenderWindow(self.renwin)
			self.iren.Start()

	def addIsoSurface(self):
		print "Adding isosurface actor to renderer..."
 		filename = self.params.colour_file
		self.isoreader = vtk.vtkStructuredPointsReader()
		self.isoreader.SetFileName(filename)
		self.isoreader.Update()

		self.extractiso = vtk.vtkExtractVOI()
		if (self.params.xmax > 0 or self.params.ymax > 0 or self.params.zmax > 0):
			self.extractiso.SetVOI(self.params.zmin,self.params.zmax,self.params.ymin,self.params.ymax,self.params.xmin,self.params.xmax)
		self.extractiso.SetInput(self.isoreader.GetOutput())

		iso = vtk.vtkMarchingCubes()
		iso.SetInput(self.extractiso.GetOutput())
		iso.SetValue(0,self.params.isovalue)
		poly_map = vtk.vtkPolyDataMapper()
		poly_map.SetInput(iso.GetOutput())
		poly_map.ScalarVisibilityOff()
		self.isoactor = vtk.vtkActor()
		self.isoactor.SetMapper(poly_map)

		if (self.params.tuecolours == 1 ):
			#self.isoactor.GetProperty().SetColor(0.678431373,0.125490196,0.678431373)
			self.isoactor.GetProperty().SetColor(173.0/255.0,32.0/255.0,173.0/255.0)
			self.isoactor.GetProperty().SetOpacity(0.4)
		else:
			self.isoactor.GetProperty().SetColor(self.params.isocolour[0],self.params.isocolour[1],self.params.isocolour[2])
			self.isoactor.GetProperty().SetOpacity(self.params.isoopacity)

		if self.params.isosurface == 1:
			self.ren.AddActor(self.isoactor)

	def addRockIsoSurface(self):
		print "Adding rock isosurface actor to renderer..."
 		filename=self.params.rock_file
		self.rockisoreader = vtk.vtkStructuredPointsReader()
		self.rockisoreader.SetFileName(filename)
		self.rockisoreader.Update()

		self.extractrockiso = vtk.vtkExtractVOI()
		if (self.params.xmax > 0 or self.params.ymax > 0 or self.params.zmax > 0):
			self.extractrockiso.SetVOI(self.params.zmin,self.params.zmax,self.params.ymin,self.params.ymax,self.params.xmin,self.params.xmax)
		self.extractrockiso.SetInput(self.rockisoreader.GetOutput())

		rockiso = vtk.vtkMarchingCubes()
		rockiso.SetInput(self.extractrockiso.GetOutput())
		rockiso.SetValue(0,self.params.rockisovalue)
		rockpoly_map = vtk.vtkPolyDataMapper()
		rockpoly_map.SetInput(rockiso.GetOutput())
		rockpoly_map.ScalarVisibilityOff()
		self.rockactor = vtk.vtkActor()
		self.rockactor.SetMapper(rockpoly_map)

		if (self.params.tuecolours == 1 ):
			self.rockactor.GetProperty().SetColor(self.params.rockisocolour[0],self.params.rockisocolour[1],self.params.rockisocolour[2])
			self.rockactor.GetProperty().SetOpacity(self.params.rockisoopacity)
		else:
			self.rockactor.GetProperty().SetColor(self.params.rockisocolour[0],self.params.rockisocolour[1],self.params.rockisocolour[2])
			self.rockactor.GetProperty().SetOpacity(self.params.rockisoopacity)

		if self.params.rockisosurface == 1:
			self.ren.AddActor(self.rockactor)


	def addVolume(self):
		print "Adding volume actor to renderer..."
		filename=self.params.colour_file
		self.volreader = vtk.vtkStructuredPointsReader()
		self.volreader.SetFileName(filename)
		self.volreader.Update()
		#self.volreader.GetOutput().SetOrigin(1,1,1)
		#volrange=self.volreader.GetOutput().GetScalarRange()
		volrange=(-1.1,1.1)
		imageShift=vtk.vtkImageShiftScale()
		imageShift.SetShift(-1*volrange[0])
		imageShift.SetScale(255.0/(volrange[1]-volrange[0]))
		imageShift.SetOutputScalarTypeToUnsignedChar()
		imageShift.SetInput(self.volreader.GetOutput())

		extractVol = vtk.vtkExtractVOI()
		if (self.params.xmax > 0 or self.params.ymax > 0 or self.params.zmax > 0):
			extractVol.SetVOI(self.params.zmin,self.params.zmax,self.params.ymin,self.params.ymax,self.params.xmin,self.params.xmax)
		extractVol.SetInput(imageShift.GetOutput())

		volume_mapper = vtk.vtkVolumeRayCastMapper()
		volume_mapper.SetInput(extractVol.GetOutput())

		composite_function = vtk.vtkVolumeRayCastCompositeFunction()
		volume_mapper.SetVolumeRayCastFunction(composite_function)

		color_transfer_func = vtk.vtkColorTransferFunction()

		if (self.params.tuecolours == 1 ):
			color_transfer_func.AddRGBPoint( 0, 0.0625, 0.0625, 0.44921875 )
			color_transfer_func.AddRGBPoint( (255/2)-10, 0.0625, 0.0625, 0.44921875 )
			color_transfer_func.AddRGBPoint( (255/2)+10, 0.8359375, 0.0, 0.2890625 )
			color_transfer_func.AddRGBPoint( 255, 0.8359375, 0.0, 0.2890625 )
		else:
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

	def addRockVolume(self):
		print "Adding rock volume actor to renderer..."
		filename=self.params.rock_file
		self.rockvolreader = vtk.vtkStructuredPointsReader()
		self.rockvolreader.SetFileName(filename)
		self.rockvolreader.Update()
		#self.rockvolreader.GetOutput().SetOrigin(1,1,1)
		#volrange=self.rockvolreader.GetOutput().GetScalarRange()
		volrange=(0.0,1.0)
		imageShift=vtk.vtkImageShiftScale()
		imageShift.SetShift(-1*volrange[0])
		imageShift.SetScale(255.0/(volrange[1]-volrange[0]))
		imageShift.SetOutputScalarTypeToUnsignedChar()
		imageShift.SetInput(self.rockvolreader.GetOutput())

		extractVol = vtk.vtkExtractVOI()
		if (self.params.xmax > 0 or self.params.ymax > 0 or self.params.zmax > 0):
			extractVol.SetVOI(self.params.zmin,self.params.zmax,self.params.ymin,self.params.ymax,self.params.xmin,self.params.xmax)
		extractVol.SetInput(imageShift.GetOutput())

		volume_mapper = vtk.vtkVolumeRayCastMapper()
		volume_mapper.SetInput(extractVol.GetOutput())

		composite_function = vtk.vtkVolumeRayCastCompositeFunction()
		volume_mapper.SetVolumeRayCastFunction(composite_function)

		color_transfer_func = vtk.vtkColorTransferFunction()
		color_transfer_func.AddRGBPoint( 0,self.params.rockvolumecolourlow[0], self.params.rockvolumecolourlow[1], self.params.rockvolumecolourlow[2] )
		color_transfer_func.AddRGBPoint( 255, self.params.rockvolumecolourhigh[0], self.params.rockvolumecolourhigh[1], self.params.rockvolumecolourhigh[2] )

		opacity_transfer_func = vtk.vtkPiecewiseFunction()
		opacity_transfer_func.AddPoint( 0, self.params.rockvolumeopacitylow )
		opacity_transfer_func.AddPoint( 118, self.params.rockvolumeopacitylow )
		opacity_transfer_func.AddPoint( 138, self.params.rockvolumeopacityhigh )
		opacity_transfer_func.AddPoint( 255, self.params.rockvolumeopacityhigh )

		volume_properties = vtk.vtkVolumeProperty()
		volume_properties.SetColor( color_transfer_func )
		volume_properties.SetScalarOpacity( opacity_transfer_func )

		volume = vtk.vtkVolume()
		volume.SetMapper( volume_mapper )
		volume.SetProperty( volume_properties )

		self.ren.AddVolume( volume )

	def addOutline(self):
		print "Adding outline actor to renderer..."
		# Create outline around data
		outline = vtk.vtkOutlineFilter()
		if (self.params.colour == 1):
			outline.SetInput(self.isoreader.GetOutput())
		elif (self.params.rock == 1):
			outline.SetInput(self.rockisoreader.GetOutput())
		outlineMapper = vtk.vtkPolyDataMapper()
		outlineMapper.SetInput( outline.GetOutput() )

		outlineActor = vtk.vtkActor()
		outlineActor.SetMapper( outlineMapper )
		outlineActor.GetProperty().SetColor( 0, 0, 0 )
		outlineActor.GetProperty().SetLineWidth( self.params.outlinethickness )

		self.ren.AddActor(outlineActor)

		# If we are looking at just a part of the system, create an extra outline to make it easier to see which part is rendered
		if ( ( self.params.xmax > 0 or self.params.ymax > 0 or self.params.zmax > 0 ) and self.params.sliceoutline == 1):
			extractoutlineiso = vtk.vtkExtractVOI()
			extractoutlineiso.SetVOI(self.params.zmin,self.params.zmax,self.params.ymin,self.params.ymax,self.params.xmin,self.params.xmax)
			if (self.params.colour == 1):
				extractoutlineiso.SetInput(self.isoreader.GetOutput())
			else:
				extractoutlineiso.SetInput(self.rockisoreader.GetOutput())

			outline = vtk.vtkOutlineFilter()
			outline.SetInput(extractoutlineiso.GetOutput())
			outlineMapper = vtk.vtkPolyDataMapper()
			outlineMapper.SetInput( outline.GetOutput() )

			outlineActor = vtk.vtkActor()
			outlineActor.SetMapper( outlineMapper )
			outlineActor.GetProperty().SetColor( 0, 0, 0 )
			outlineActor.GetProperty().SetLineWidth( self.params.outlinethickness )

			self.ren.AddActor(outlineActor)

	def setCamera(self):
		cam1 = self.ren.GetActiveCamera()
		if (self.params.parallelprojection == 1):
			cam1.ParallelProjectionOn()
		cam1.SetViewUp(1, 0, 0)
		if (self.params.colour == 1):
			cam1.SetFocalPoint(self.isoreader.GetOutput().GetCenter())
			size=[2*coord+1 for coord in self.isoreader.GetOutput().GetCenter()]
		elif (self.params.rock == 1):
			cam1.SetFocalPoint(self.rockisoreader.GetOutput().GetCenter())
			size=[2*coord+1 for coord in self.rockisoreader.GetOutput().GetCenter()]
		else:
			size= [ self.params.systemsize[2], self.params.systemsize[1], self.params.systemsize[0] ]
			centre = [0.5*coord for coord in size]
			cam1.SetFocalPoint(centre)

		cam1.SetPosition(size[0]/2,size[1]/2,8*math.sqrt(size[2]**2+size[1]**2+size[0]**2))

		cam1.Elevation(self.params.elevation)
		cam1.Azimuth(self.params.azimuth)
		cam1.Roll(self.params.roll)
		cam1.Zoom(self.params.zoom)

	def dumpImage(self):
		filename = 'rendered-'+self.params.id+'.png'
		toimage = vtk.vtkWindowToImageFilter()
		toimage.SetInput(self.renwin)
		toimage.Modified()
		writer = vtk.vtkPNGWriter()
		writer.SetInput(toimage.GetOutput())
		writer.SetFileName(filename)
		writer.Modified()
		print "OMNIRENDER writing image",filename
		writer.Write()

	#Reading Text file from the System for plotting Spheres
	def readSpherePoints(self):
		print "  Reading sphere points..."
		filename = self.params.md_file
		points = vtk.vtkPoints()
		trpoints = vtk.vtkPoints()
		file = open(filename,"r")
		#xdim = int(raw_input("Please Enter the X-Dimension: "))
		#ydim = int(raw_input("Please Enter the Y-Dimension: "))
		#zdim = int(raw_input("Please Enter the Z-Dimension: "))
		m = 0
		tr = 0

		ntl = len(self.params.spheretracerlist)

		while file:
			line = file.readline()
			if line == "":
				break
			line = str.rstrip(line)
			cols =line.split(' ')

			if len(cols)  <3:
				print "Omitting line: ", line
				continue

			x=cols[0]
			y=cols[1]
			z=cols[2]

			if ( ( ( self.params.xmax == 0 ) or ( float(x) <= self.params.xmax and float(x) >= self.params.xmin ) ) and
			     ( ( self.params.ymax == 0 ) or ( float(y) <= self.params.ymax and float(y) >= self.params.ymin ) ) and
			     ( ( self.params.zmax == 0 ) or ( float(z) <= self.params.zmax and float(z) >= self.params.zmin ) ) ):

				if ( ntl > 0):
					if ( self.params.spheretracerlist.count(int(cols[self.params.sphereidcolumn])) == 0):
						points.InsertNextPoint(float(z)-1,float(y)-1, float(x)-1)
						#points.InsertNextPoint(float(z),float(y), float(x))
						m+=1
					else:
						trpoints.InsertNextPoint(float(z)-1,float(y)-1, float(x)-1)
						tr+=1
				elif ( self.params.spheretracerodd ):
					if ( int(cols[self.params.sphereidcolumn]) % 2 == 0 ):
						points.InsertNextPoint(float(z)-1,float(y)-1, float(x)-1)
						m+=1
					else:
						trpoints.InsertNextPoint(float(z)-1,float(y)-1, float(x)-1)
						tr+=1
				else:
					points.InsertNextPoint(float(z)-1,float(y)-1, float(x)-1)
					m+=1

		# Copied Sphere code from Finalstreamtracermod
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

		# Copied Sphere code from Finalstreamtracermod
		trspheredata=vtk.vtkPolyData()
		trspheredata.SetPoints(trpoints)

		vector=vtk.vtkFloatArray()
		vector.SetNumberOfComponents(3)
		vector.SetNumberOfTuples(tr)

		pline=vtk.vtkCellArray()

		pline.InsertNextCell(1)
		for i in range(0,1):
			pline.InsertCellPoint(i)
		trspheredata.SetPolys(pline)
		trspheredata.GetPointData().SetVectors(vector)
		print "  Succesfully read ",tr+m," sphere points."
		return [spheredata,trspheredata]

	# Reading Text file from the System for plotting ellipsoids
	def readPointsandOrientations(self):
		print "  Reading sphere points..."
		filename = self.params.md_file
		points = vtk.vtkPoints()
		trpoints = vtk.vtkPoints()
		orientations = vtk.vtkFloatArray()
		orientations.SetNumberOfComponents(3)
		trorientations = vtk.vtkFloatArray()
		trorientations.SetNumberOfComponents(3)
		file = open(filename,"r")
		#xdim = int(raw_input("Please Enter the X-Dimension: "))
		#ydim = int(raw_input("Please Enter the Y-Dimension: "))
		#zdim = int(raw_input("Please Enter the Z-Dimension: "))
		m = 0
		tr = 0

		ntl = len(self.params.ellipsoidtracerlist)

		while file:
			line = file.readline()
			if line == "":
				break
			line = str.rstrip(line)
			cols =line.split(' ')

			if len(cols)  <3:
				print "Omitting line: ", line
				continue

			x=cols[0]
			y=cols[1]
			z=cols[2]
			ox=cols[9]
			oy=cols[10]
			oz=cols[11]

			if ( ( ( self.params.xmax == 0 ) or ( float(x) <= self.params.xmax and float(x) >= self.params.xmin ) ) and
			     ( ( self.params.ymax == 0 ) or ( float(y) <= self.params.ymax and float(y) >= self.params.ymin ) ) and
			     ( ( self.params.zmax == 0 ) or ( float(z) <= self.params.zmax and float(z) >= self.params.zmin ) ) ):

				if ( ntl > 0):
					if (self.params.ellipsoidtracerlist.count(int(cols[self.params.ellipsoididcolumn])) == 0):
						points.InsertNextPoint(float(z)-1,float(y)-1, float(x)-1)
						orientations.InsertNextTuple3(float(oz), float(oy), float(ox))
						m+=1
					else:
						trpoints.InsertNextPoint(float(z)-1,float(y)-1, float(x)-1)
						trorientations.InsertNextTuple3(float(oz), float(oy), float(ox))
						tr+=1
				else:
					points.InsertNextPoint(float(z)-1,float(y)-1, float(x)-1)
					orientations.InsertNextTuple3(float(oz), float(oy), float(ox))
					m+=1

		ellipsoiddata = vtk.vtkUnstructuredGrid()
		ellipsoiddata.SetPoints(points)
		ellipsoiddata.GetPointData().SetNormals(orientations)

		trellipsoiddata = vtk.vtkUnstructuredGrid()
		trellipsoiddata.SetPoints(trpoints)
		trellipsoiddata.GetPointData().SetNormals(trorientations)

		print "  Succesfully read sphere points."
		return [ellipsoiddata,trellipsoiddata]


	# Sphere visualization.
	def addSpheres(self):
		spheres=vtk.vtkSphereSource()
		spheres.SetRadius(self.params.sphereradius)
		spheres.SetThetaResolution(32)
		spheres.SetPhiResolution(32)

		allspheredata = self.readSpherePoints()
		spheredata = allspheredata[0]
		trspheredata = allspheredata[1]

		print "Adding normal spheres actor to renderer..."
		self.sphereGlyph = vtk.vtkGlyph3D()
		self.sphereGlyph.SetSource(spheres.GetOutput())
		self.sphereGlyph.SetInput(spheredata)

		# Mapping of sphere data and also adding actors for rendering the spheres.
		# Here you can also change the color of sphere changing the parameters into SetColor.
		# You can put any values and observe changes but (0,0,0) is used for black and (1,1,1) is for white.

		sphereMapper = vtk.vtkPolyDataMapper()
		sphereMapper.SetInput(self.sphereGlyph.GetOutput())
		sphereMapper.GlobalImmediateModeRenderingOn()
		self.sphereactor = vtk.vtkLODActor()
		self.sphereactor.SetMapper(sphereMapper)
		self.sphereactor.GetProperty().SetColor(self.params.spherecolour[0],self.params.spherecolour[1] ,self.params.spherecolour[2] )
		self.sphereactor.GetProperty().SetOpacity(self.params.sphereopacity)

		self.ren.AddActor(self.sphereactor)

		# TRACERS
		print "Adding tracer spheres actor to renderer..."
		self.trsphereGlyph = vtk.vtkGlyph3D()
		self.trsphereGlyph.SetSource(spheres.GetOutput())
		self.trsphereGlyph.SetInput(trspheredata)

		# Mapping of sphere data and also adding actors for rendering the spheres.
		# Here you can also change the color of sphere changing the parameters into SetColor.
		# You can put any values and observe changes but (0,0,0) is used for black and (1,1,1) is for white.

		trsphereMapper = vtk.vtkPolyDataMapper()
		trsphereMapper.SetInput(self.trsphereGlyph.GetOutput())
		trsphereMapper.GlobalImmediateModeRenderingOn()
		self.trsphereactor = vtk.vtkLODActor()
		self.trsphereactor.SetMapper(trsphereMapper)
		self.trsphereactor.GetProperty().SetColor(self.params.spheretracercolour[0],self.params.spheretracercolour[1],self.params.spheretracercolour[2])
		self.trsphereactor.GetProperty().SetOpacity(self.params.spheretraceropacity)

		self.ren.AddActor(self.trsphereactor)


	# Ellipsoids visualization.
	def addEllipsoids(self):
		ellipsoids=vtk.vtkSphereSource()
		ellipsoids.SetRadius(self.params.sphereradius)
		ellipsoids.SetThetaResolution(self.params.ellipsoidresolution)
		ellipsoids.SetPhiResolution(self.params.ellipsoidresolution)

		sphereToEllipsoid = vtk.vtkTransform()
		sphereToEllipsoid.Scale(self.params.ellipsoidanisotropy,1,1)

		sphereToEllipsoidFilter = vtk.vtkTransformPolyDataFilter()
		sphereToEllipsoidFilter.SetInputConnection(ellipsoids.GetOutputPort())
		sphereToEllipsoidFilter.SetTransform(sphereToEllipsoid)

		allellipsoiddata = self.readPointsandOrientations()

		print "Adding normal ellipsoids actor to renderer..."
		self.ellipsoidGlyph = vtk.vtkGlyph3D()
		self.ellipsoidGlyph.SetSource(sphereToEllipsoidFilter.GetOutput())
		self.ellipsoidGlyph.SetInput(allellipsoiddata[0])
		self.ellipsoidGlyph.OrientOn()
		self.ellipsoidGlyph.SetVectorModeToUseNormal()

		# Mapping of ellipsoid data and also adding actors for rendering the ellipsoids.
		# Here you can also change the color of ellipsoid changing the parameters into SetColor.
		# You can put any values and observe changes but (0,0,0) is used for black and (1,1,1) is for white.

		ellipsoidMapper = vtk.vtkPolyDataMapper()
		ellipsoidMapper.SetInput(self.ellipsoidGlyph.GetOutput())
		ellipsoidMapper.GlobalImmediateModeRenderingOn()
		self.ellipsoidactor = vtk.vtkLODActor()
		self.ellipsoidactor.SetMapper(ellipsoidMapper)
		self.ellipsoidactor.GetProperty().SetColor(self.params.ellipsoidcolour[0],self.params.ellipsoidcolour[1],self.params.ellipsoidcolour[2])
		self.ellipsoidactor.GetProperty().SetOpacity(self.params.ellipsoidopacity)

		self.ren.AddActor(self.ellipsoidactor)

                # TRACERS

		print "Adding tracer ellipsoids actor to renderer..."
		self.trellipsoidGlyph = vtk.vtkGlyph3D()
		self.trellipsoidGlyph.SetSource(sphereToEllipsoidFilter.GetOutput())
		self.trellipsoidGlyph.SetInput(allellipsoiddata[1])
		self.trellipsoidGlyph.OrientOn()
		self.trellipsoidGlyph.SetVectorModeToUseNormal()

		# Mapping of ellipsoid data and also adding actors for rendering the ellipsoids.
		# Here you can also change the color of ellipsoid changing the parameters into SetColor.
		# You can put any values and observe changes but (0,0,0) is used for black and (1,1,1) is for white.

		trellipsoidMapper = vtk.vtkPolyDataMapper()
		trellipsoidMapper.SetInput(self.trellipsoidGlyph.GetOutput())
		trellipsoidMapper.GlobalImmediateModeRenderingOn()
		self.trellipsoidactor = vtk.vtkLODActor()
		self.trellipsoidactor.SetMapper(trellipsoidMapper)
		self.trellipsoidactor.GetProperty().SetColor(self.params.ellipsoidtracercolour[0],self.params.ellipsoidtracercolour[1],self.params.ellipsoidtracercolour[2])
		self.trellipsoidactor.GetProperty().SetOpacity(self.params.ellipsoidtraceropacity)

		self.ren.AddActor(self.trellipsoidactor)


	def update(self,fid):
		self.params.set_id(fid)
		if (self.params.colour == 1):
			if (self.params.isosurface == 1):
				self.isoreader.SetFileName(self.params.colour_file)
				self.isoreader.Update()
			if (self.params.volume == 1):
				self.volreader.SetFileName(self.params.colour_file)
				self.volreader.Update()

		if (self.params.rock == 1):
			if (self.params.rockisosurface == 1):
				self.rockisoreader.SetFileName(self.params.rock_file)
				self.rockisoreader.Update()
			if (self.params.rockvolume == 1):
				self.rockvolreader.SetFileName(self.params.rock_file)
				self.rockvolreader.Update()

		if (self.params.spheres == 1):
			self.sphereGlyph.SetInput(self.readSpherePoints()[0])
			self.trsphereGlyph.SetInput(self.readSpherePoints()[1])

		if (self.params.ellipsoids == 1):
			self.ellipsoidGlyph.SetInput(self.readPointsandOrientations()[0])
			self.trellipsoidGlyph.SetInput(self.readPointsandOrientations()[1])

		self.ren.ResetCameraClippingRange()

		# if (self.params.spheres == 1):
		# 	self.sphereactor.GetProperty().SetOpacity(1.0)
		# if (self.params.colour == 1):
		# 	self.isoactor.GetProperty().SetOpacity(self.params.isoopacity)
		# self.ren.Render()
		# self.dumpImage('left')

		# if (self.params.spheres == 1):
		# 	self.sphereactor.GetProperty().SetOpacity(self.params.sphereopacity)
		# if (self.params.colour == 1):
		# 	self.isoactor.GetProperty().SetOpacity(1.0)
		# self.ren.Render()
		# self.dumpImage('right')

		print "RENDERING..."
		self.ren.Render()
		print " Done rendering."
		self.dumpImage()

	class ParamStorage():
		def __init__(self):
			# Camera properties
			self.elevation = 10.0
			self.zoom = 4.0
			self.azimuth = -10.0
			self.roll = 0.0
			self.systemsize = [0.0, 0.0, 0.0]
			self.parallelprojection = 0

			# Fluid properties
			self.colour = 1

			# ISOSURFACES

			self.isosurface = 1
			self.isoopacity = 0.4
			self.isovalue = 0.0
			self.isocolour = [0.8, 0.2, 0.9]

			# VOLUME RENDERING
			self.volume = 1
			self.volumeopacitylow = 0.01
			self.volumeopacityhigh = 0.01

			# Rock properties
			self.rock = 1

			# ISOSURFACES
			self.rockisosurface = 1
			self.rockisoopacity = 0.2
			self.rockisovalue = 1.0
			self.rockisocolour = [0.7, 0.7, 0.7]

			# VOLUME RENDERING
			self.rockvolume = 1
			self.rockvolumecolourlow = [0.5, 0.5, 0.5]
			self.rockvolumecolourhigh = [0.5, 0.5, 0.5]
			self.rockvolumeopacitylow = 0.01
			self.rockvolumeopacityhigh = 0.01

			# Spheres
			self.spheres = 1
			self.sphereradius = 2.0
			self.sphereopacity = 0.1
			self.spherecolour = [1.0, 1.0 , 1.0]
			self.spheretracerodd = 0
			self.sphereidcolumn = 18
			self.spheretracerlist = [ ]
			self.spheretraceropacity = 1.0
			self.spheretracercolour = [0.0, 0.0 , 0.0]

			# ellipsoids
			self.ellipsoids = 0
			self.ellipsoidopacity = 0.1
			self.ellipsoidcolour = [0.0, 1.0 , 0.0]
			self.ellipsoididcolumn = 16
			self.ellipsoidtracerlist = [ ]
			self.ellipsoidtraceropacity = 1.0
			self.ellipsoidtracercolour = [0.0, 0.0 , 0.0]
			self.ellipsoidanisotropy=2
			self.ellipsoidresolution=32

			# Miscellaneous
			self.hsize = 800
			self.vsize = 600

			self.id=""
			self.md_file = ""
			self.colour_file = ""
			self.rock_file = ""

			self.writeimage = 0
			self.interactive = 0
			self.outline = 1
			self.outlinethickness = 1.0
			self.sliceoutline = 0
			self.offscreenrendering = 0
			self.tuecolours = 0

			self.xmin = 0
			self.xmax = 0
			self.ymin = 0
			self.ymax = 0
			self.zmin = 0
			self.zmax = 0

		def processOption(self,s):
			# Camera properties
			if s.startswith('-roll='):
				self.roll = float(s.split('=')[1])
				print "  Roll:", self.roll
			if s.startswith('-elevation='):
				self.elevation = float(s.split('=')[1])
				print "  Elevation:", self.elevation
			if s.startswith('-zoom='):
				self.zoom = float(s.split('=')[1])
				print "  Zoom:", self.zoom
			if s.startswith('-azimuth='):
				self.azimuth = float(s.split('=')[1])
				print "  Azimuth:", self.azimuth
			if s.startswith('-systemsize='):
				self.systemsize = [float(c) for c in s.split('=')[1].split(',')]
				print '  System size:',self.systemsize
			if s.startswith('-parallelprojection'):
				print "  With parallel projection"
				self.parallelprojection = 1

			# Fluid properties
			if s.startswith('-nocolourdata'):
				print "  No colour data"
				self.colour = 0

			# ISOSURFACES
			if s.startswith('-noisosurface'):
				print "  No isosurface"
				self.isosurface=0
			if s.startswith('-isoopacity='):
				self.isoopacity = float(s.split('=')[1])
				print "  Surface opacity:", self.isoopacity
			if s.startswith('-isovalue='):
				self.isovalue = float(s.split('=')[1])
				print "  Isosurface value:", self.isovalue
			if s.startswith('-isocolour='):
				self.isocolour = [float(c) for c in s.split('=')[1].split(',')]
				print '  Isosurface colour:',self.isocolour

			# VOLUME RENDERING
			if s.startswith('-novolume'):
				print "  No volume rendering"
				self.volume = 0
			if s.startswith('-volumeopacitylow='):
				self.volumeopacitylow = float(s.split('=')[1])
				print "  Volume opacity low:", self.volumeopacitylow
			if s.startswith('-volumeopacityhigh='):
				self.volumeopacityhigh = float(s.split('=')[1])
				print "  Volume opacity high:", self.volumeopacityhigh

			# Rock properties
			if s.startswith('-norockdata'):
				print "  No rock data"
				self.rock = 0

			# ISOSURFACES
			if s.startswith('-norockisosurface'):
				print "  No rock surface"
				self.rockisosurface=0
			if s.startswith('-rockisoopacity='):
				self.rockisoopacity = float(s.split('=')[1])
				print "  Rock isosurface opacity:", self.rockisoopacity
			if s.startswith('-rockisovalue='):
				self.rockisovalue = float(s.split('=')[1])
				print "  Rock isosurface value:", self.rockisovalue
			if s.startswith('-rockisocolour='):
				self.rockisocolour = [float(c) for c in s.split('=')[1].split(',')]
				print '  Rock isosurface colour:',self.rockisocolour

			# VOLUME RENDERING
			if s.startswith('-norockvolume'):
				print "  No rock volume rendering"
				self.rockvolume = 0
			if s.startswith('-rockvolumecolourlow='):
				self.rockvolumecolourlow = [float(c) for c in s.split('=')[1].split(',')]
				print "  Rock volume colour low:", self.rockvolumecolourlow
			if s.startswith('-rockvolumecolourhigh='):
				self.rockvolumecolourhigh = [float(c) for c in s.split('=')[1].split(',')]
				print "  Rock volume colour low:", self.rockvolumecolourhigh
			if s.startswith('-rockvolumeopacitylow='):
				self.rockvolumeopacitylow = float(s.split('=')[1])
				print "  Rock volume opacity low:", self.rockvolumeopacitylow
			if s.startswith('-rockvolumeopacityhigh='):
				self.rockvolumeopacityhigh = float(s.split('=')[1])
				print "  Rock volume opacity high:", self.rockvolumeopacityhigh


			# Spheres properties
			if s.startswith('-nospheres'):
				print "  No spheres"
				self.spheres = 0
			if s.startswith('-sphereradius='):
				self.sphereradius = float(s.split('=')[1])
				print "  Sphere radius:", self.sphereradius
			if s.startswith('-sphereopacity='):
				self.sphereopacity = float(s.split('=')[1])
				print "  Sphere opacity:", self.sphereopacity
			if s.startswith('-spherecolour='):
				self.spherecolour = [float(c) for c in s.split('=')[1].split(',')]
				print "  Sphere colour:", self.spherecolour
			if s.startswith('-sphereidcolumn='):
				self.sphereidcolumn = int(s.split('=')[1])
				print "  Sphere ID column:", self.sphereidcolumn
			if s.startswith('-spheretracerlist='):
				self.spheretracerlist = [int(c) for c in s.split('=')[1].split(',')]
				print "  Sphere tracer list:", self.spheretracerlist
			if s.startswith('-spheretraceropacity='):
				self.spheretraceropacity = float(s.split('=')[1])
				print "  Sphere tracer opacity:", self.spheretraceropacity
			if s.startswith('-spheretracercolour='):
				self.spheretracercolour = [float(c) for c in s.split('=')[1].split(',')]
				print "  Sphere tracer colour:", self.spheretracercolour


			# Ellipsoid properties
			if s.startswith('-ellipsoids'):
				print "  Using ellipsoids, no spheres"
				self.ellipsoids = 1
				self.spheres = 0
			if s.startswith('-ellipsoidanisotropy='):
				self.ellipsoidanisotropy = float(s.split('=')[1])
				print "  Anisotropy:", self.ellipsoidanisotropy
			if s.startswith('-ellipsoidresolution='):
				self.ellipsoidresolution = int(s.split('=')[1])
				print "  Ellipsoid resolution:", self.ellipsoidresolution
			if s.startswith('-ellipsoidopacity='):
				self.ellipsoidopacity = float(s.split('=')[1])
				print "  Ellipsoidopacity opacity:", self.ellipsoidopacity
			if s.startswith('-ellipsoidcolour='):
				self.ellipsoidcolour = [float(c) for c in s.split('=')[1].split(',')]
				print "  Ellipsoid colour:", self.ellipsoidcolour
			if s.startswith('-ellipsoididcolumn='):
				self.ellipsoididcolumn = int(s.split('=')[1])
				print "  Ellipsoid ID column:", self.ellipsoididcolumn
			if s.startswith('-ellipsoidtracerlist='):
				self.ellipsoidtracerlist = [int(c) for c in s.split('=')[1].split(',')]
				print "  Ellipsoid tracer list:", self.ellipsoidtracerlist
			if s.startswith('-ellipsoidtraceropacity='):
				self.ellipsoidtraceropacity = float(s.split('=')[1])
				print "  Ellipsoid tracer opacity:", self.ellipsoidtraceropacity
			if s.startswith('-ellipsoidtracercolour='):
				self.ellipsoidtracercolour = [float(c) for c in s.split('=')[1].split(',')]
				print "  Ellipsoid tracer colour:", self.ellipsoidtracercolour


			# Miscellaneous properties
			if s.startswith('-size='):
				self.hsize,self.vsize = [int(size) for size in s.split('=')[1].split('x')]
				print "  Size:", self.hsize, '*', self.vsize
			if s.startswith('-id='):
				if self.set_id(s.split('=')[1]):
					print "  Colourfile:", self.colour_file, "MD-File:", self.md_file
			if s.startswith('-writeimage'):
				print "  Writing image"
				self.writeimage = 1
			if s.startswith('-interactive'):
				print "  Interactive session"
				self.interactive = 1
			if s.startswith('-sliceoutline'):
				print "  With slice outline"
				self.sliceoutline = 1
			if s.startswith('-nooutline'):
				print "  No outline"
				self.outline = 0
			if s.startswith('-outlinethickness='):
				self.outlinethickness = float(s.split('=')[1])
				print "  Outline thickness:", self.outlinethickness
			if s.startswith('-offscreenrendering'):
				print "  With offscreen rendering"
				self.offscreenrendering = 1
			if s.startswith('-tuecolours'):
				print "  With TU/e colours"
				self.tuecolours = 1
				self.spherecolour = [ 0.0, 0.671875, 0.5078125 ]
				self.ellipsoidcolour = [ 0.0, 0.671875, 0.5078125 ]

			if s.startswith('-xmin'):
				self.xmin = int(s.split('=')[1])
				print "  X min:", self.xmin
			if s.startswith('-xmax'):
				self.xmax = int(s.split('=')[1])
				print "  X max:", self.xmax
			if s.startswith('-ymin'):
				self.ymin = int(s.split('=')[1])
				print "  Y min:", self.ymin
			if s.startswith('-ymax'):
				self.ymax = int(s.split('=')[1])
				print "  Y max:", self.ymax
			if s.startswith('-zmin'):
				self.zmin = int(s.split('=')[1])
				print "  Z min:", self.zmin
			if s.startswith('-zmax'):
				self.zmax = int(s.split('=')[1])
				print "  Z max:", self.zmax

		def set_id(self, fid):
			cfn="colour_"+fid+".vtk"
			mdfn="md-cfg_"+fid+".asc"
			rfn="rock_"+fid+".vtk"

			print "OMNIRENDER opening files for ID:", fid

			if (self.colour == 1):
				try:
					cf=open(cfn,"r")
					cf.close()
				except IOError:
					print "  IOError while trying to open colour file. Assuming -nocolourdata."
					self.colour = 0
				else:
					print "  Opened colour file."
					self.colour_file=cfn
			if (self.spheres == 1):
				try:
					mdf=open(mdfn,"r")
					mdf.close()
				except IOError:
					print "  IOError while trying to open spheres file. Assuming -nospheres."
					self.spheres = 0
				else:
					print "  Opened spheres file."
					self.md_file=mdfn
			if (self.ellipsoids == 1):
				try:
					mdf=open(mdfn,"r")
					mdf.close()
				except IOError:
					print "  IOError while trying to open ellipsoids file. Assuming -nospheres."
					self.spheres = 0
				else:
					print "  Opened spheres file."
					self.md_file=mdfn
			if (self.rock == 1):
				try:
					rf=open(rfn,"r")
					rf.close()
				except IOError:
					print "  IOError while trying to open rock file. Assuming -norockdata."
					self.rock = 0
				else:
					print "  Opened rock file."
					self.rock_file=rfn

			self.id=fid
			return True

if __name__ == "__main__":
	import sys
	if ( len (sys.argv) < 2 ):
		print "./omnirender.py <file id> [<params/flags>]"
		print "See README.omnirender for details"
		sys.exit()
	omnirender(sys.argv[1],sys.argv[2:])
