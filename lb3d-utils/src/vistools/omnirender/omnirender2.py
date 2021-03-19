#!/usr/bin/env python
# This is based on spheres_and_isosurf2.py / omnirender.py but is probably unrecognizable by now

import colorsys # For rainbow
import math
import re
import sys

from optparse import OptionParser, OptionGroup

# Showing versions of python and external module files
def print_versions():
    print( "\nPython       " + sys.version)

    try:
        import numpy
        print( "\nNumPy        " + numpy.__version__)
        print( "             " + numpy.__file__)
    except ImportError:
        print( "\nNumPy        not found.")

    try:
        import vtk
        print( "\nVTK          " + vtk.vtkVersion.GetVTKSourceVersion() ) # __version__ symbol does not exist
        print( "             " + vtk.__file__)
    except ImportError:
        print( "\nVTK          not found.")

    import os
    print( "\nomnirender2 " )
    print( "             " + os.path.realpath(__file__) )
    print( "" )
    exit()

def unitVector(v):
    return v / sum(v * v)**0.5

# Reading LB3D .asc file
def readParticleData(fname, poscol, idcol, ocol, rcol, tlist = [], inverse = False):
    ex = np.array([1.0, 0.0, 0.0])
    ey = np.array([0.0, 1.0, 0.0])

    points = vtk.vtkPoints()

    # Include tensor for scaling in different directions for anisotropic particles
    scale = vtk.vtkFloatArray()
    scale.SetNumberOfComponents(9)

    # Scalar attribute to set different colour per particle
    colours = vtk.vtkFloatArray()

    # Open file for reading
    file_md = open(fname,"r")

    # Count particles that are actually placed
    n_ptcl = 0

    while file_md:
        # End reading on empty lines, skip commented lines, then split lines into columns by whitespace
        line = file_md.readline()
        if (line == ""): break
        if (line[0] == "#"): continue
	line = str.rstrip(line)
        cols = line.split(' ')

        # If not explicitly specified, assume the id column is the last column in the file
        if (idcol < 0): idcol = len(cols) - 1

        # Get particle ID
        pid = int(cols[idcol])

        # pidcheck is a function which should return a boolean value based on the supplied particle id
        # this can be used to make multiple lists of particles to have different colours etc

        if ( ( not inverse and ( pid in tlist or len(tlist) == 0 ) ) or ( inverse and ( not pid in tlist ) ) ):
            # Positions are found in 3 adjacent columns starting at poscol
            x = float(cols[poscol+0])
            y = float(cols[poscol+1])
            z = float(cols[poscol+2])

            # Orientations are found in 3 adjacent columns starting at poscol
            if (ellipsoids):
                # particle's axis of rotational symmetry
                o = np.array([float(cols[ocol+0]),
                              float(cols[ocol+1]),
                              float(cols[ocol+2])])

                # two arbitrary axes perpendicular to o
                if (abs(np.dot(o,ex)) < abs(np.dot(o,ey))):
                    o1 = unitVector(np.cross(o,ex))
                else:
                    o1 = unitVector(np.cross(o,ey))
                o2 = np.cross(o,o1)

            # Insert into list iff particle is in the region to be rendered
            if ( ( options.xmax < 0 or ( x <= options.xmax and x >= options.xmin ) ) and
                 ( options.ymax < 0 or ( y <= options.ymax and y >= options.ymin ) ) and
                 ( options.zmax < 0 or ( z <= options.zmax and z >= options.zmin ) ) ):

                points.InsertNextPoint(z-1, y-1, x-1)
                if (ellipsoids):
                    if (options.particlepolydispersity):
                        r_orth = float(cols[rcol+0])
                        r_para = float(cols[rcol+1])
                    else:
                        r_orth = options.particleradius
                        r_para = options.particleradius * options.particleanisotropy

                    scale.InsertNextTuple9(o[2]*r_para, o[1]*r_para, o[0]*r_para,
                                           o1[2]*r_orth, o1[1]*r_orth, o1[0]*r_orth,
                                           o2[2]*r_orth, o2[1]*r_orth, o2[0]*r_orth)
                    if (options.particlecolourbyaspectratio):
                        colours.InsertNextValue(
                            (r_para/r_orth-options.particlecolourbyaspectratiomin)
                            /options.particlecolourbyaspectratiomax)
                n_ptcl += 1

    file_md.close()

    # Make a VTK dataset
    ptcldata=vtk.vtkPolyData()
    ptcldata.SetPoints(points)

    # If anistropic particles, save half-axes and orientations as well
    if (ellipsoids):
        ptcldata.GetPointData().SetTensors(scale)

    if (particlecolourbydata):
        ptcldata.GetPointData().SetScalars(colours)

    if ( options.verbose): print "  Succesfully read", n_ptcl, "particles."
    return ptcldata

# Make actor for particles
def getParticleActor(ptcldata, colour, opacity):

    # Set particle size and angular resolution
    particles = vtk.vtkSphereSource()
    particles.SetThetaResolution(options.particleresolution)
    particles.SetPhiResolution(options.particleresolution)
    if (ellipsoids):
        particles.SetRadius(1.0)
    else:
        particles.SetRadius(options.particleradius)

    # Different glyphs for spheres and ellipsoids
    if (ellipsoids):
        ptclGlyph = vtk.vtkTensorGlyph()
        ptclGlyph.ExtractEigenvaluesOff()
        if (not particlecolourbydata): ptclGlyph.ColorGlyphsOff()
    else:
        ptclGlyph = vtk.vtkGlyph3D()
    ptclGlyph.SetInput(ptcldata)
    ptclGlyph.SetSource(particles.GetOutput())

    # Mapper
    ptclMapper = vtk.vtkPolyDataMapper()
    ptclMapper.SetInput(ptclGlyph.GetOutput())
    ptclMapper.GlobalImmediateModeRenderingOn()

    # Actor
    ptclActor = vtk.vtkLODActor()
    ptclActor.SetMapper(ptclMapper)
    ptclActor.GetProperty().SetColor(colour)
    ptclActor.GetProperty().SetOpacity(opacity)

    # Add to global renderer
    return ptclActor

def getParticleActorJanus(ptcldata, colour, opacity, StartTheta, EndTheta):

    # Set particle size and angular resolution
    particles = vtk.vtkSphereSource()
    particles.SetThetaResolution(options.particleresolution)
    particles.SetPhiResolution(options.particleresolution)
    particles.SetStartTheta(StartTheta)
    particles.SetEndTheta(EndTheta)
    particles.SetStartPhi(0)
    particles.SetEndPhi(180)
    if (ellipsoids):
        particles.SetRadius(1.0)
    else:
        particles.SetRadius(options.particleradius)

    # Different glyphs for spheres and ellipsoids
    if (ellipsoids):
        ptclGlyph = vtk.vtkTensorGlyph()
        ptclGlyph.ExtractEigenvaluesOff()
        if (not particlecolourbydata): ptclGlyph.ColorGlyphsOff()
    else:
        ptclGlyph = vtk.vtkGlyph3D()
    ptclGlyph.SetInput(ptcldata)
    ptclGlyph.SetSource(particles.GetOutput())

    # Mapper
    ptclMapper = vtk.vtkPolyDataMapper()
    ptclMapper.SetInput(ptclGlyph.GetOutput())
    ptclMapper.GlobalImmediateModeRenderingOn()

    # Actor
    ptclActor = vtk.vtkLODActor()
    ptclActor.SetMapper(ptclMapper)
    ptclActor.GetProperty().SetColor(colour)
    ptclActor.GetProperty().SetOpacity(opacity)

    # Add to global renderer
    return ptclActor

# Write renwin to PNG file
def writePngImage():
    fname = 'rendered-' + fid + '.png'
    toimage = vtk.vtkWindowToImageFilter()
    toimage.SetInput(renwin)
    toimage.Modified()
    writer = vtk.vtkPNGWriter()
    writer.SetInput(toimage.GetOutput())
    writer.SetFileName(fname)
    writer.Modified()
    if ( options.verbose): print "Writing image", fname
    writer.Write()

    if ( options.addmetadata): 
        if (options.verbose): print "Adding metadata"
        from PIL import Image, PngImagePlugin
        im = Image.open(fname)
        im.info["Comment"] = "Generated by omnirender2.py with options: " + str(options) + "\nCommand line options: <argv>" + cmdstr + "</argv>"
        reserved = ('interlace', 'gamma', 'dpi', 'transparency', 'aspect')
        meta = PngImagePlugin.PngInfo()
        for k,v in im.info.iteritems():
            if k in reserved: continue
            meta.add_text(k, v, 0)
        im.save(fname, "PNG", pnginfo=meta)

# Make actor for text
def getTextActor():
    textActor = vtk.vtkTextActor()
    textActor.SetTextScaleModeToViewport()
    textActor.SetDisplayPosition( int(options.hsize*options.textposition[0]), int(options.vsize*options.textposition[1]))
    textActor.SetInput(options.text)

    tprop = textActor.GetTextProperty()
    tprop.SetColor(options.textcolour[0], options.textcolour[1], options.textcolour[2])
    tprop.SetFontSize(options.textsize)
    ff = tprop.GetFontFamilyFromString(options.textfont)
    tprop.SetFontFamily(ff)
    tprop.SetJustificationToCentered()
    # tprop.SetVerticalJustificationToCentered()
    if (options.textbold):
        tprop.BoldOn()
    if (options.textitalic):
        tprop.ItalicOn()
    #tprop.ShadowOn()
    return textActor

# Make outlines
def getOutlineActors():
    # Create outline around full data
    outline = vtk.vtkOutlineFilter()
    if (options.colouriso):
        outline.SetInput(colourisoreader.GetOutput())
    elif (options.colourvolume):
        outline.SetInput(colourvolumereader.GetOutput())
    elif (options.rockiso):
        outline.SetInput(rockisoreader.GetOutput())
    elif (options.rockvolume):
        outline.SetInput(rockvolumereader.GetOutput())
    else:
        return []

    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInput( outline.GetOutput() )
	
    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper( outlineMapper )
    outlineActor.GetProperty().SetColor( options.outlinecolour )
    outlineActor.GetProperty().SetLineWidth( options.outlinethickness )

    # If we are looking at just a part of the system, create an extra outline to make it easier to see which part is rendered
    if ( ( options.xmax > 0 or options.ymax > 0 or options.zmax > 0 ) and options.sliceoutline):
	extractoutlineiso = vtk.vtkExtractVOI()
	extractoutlineiso.SetVOI(options.zmin,options.zmax, options.ymin,options.ymax, options.xmin, options.xmax)
	if (options.colouriso or options.colourvolume):
            extractoutlineiso.SetInput(colourisoreader.GetOutput())
	else:
            extractoutlineiso.SetInput(rockisoreader.GetOutput())

	outline = vtk.vtkOutlineFilter()
	outline.SetInput(extractoutlineiso.GetOutput())
	outlineMapper = vtk.vtkPolyDataMapper()
	outlineMapper.SetInput( outline.GetOutput() )
		
	outlineActorVOI = vtk.vtkActor()
	outlineActorVOI.SetMapper( outlineMapper )
	outlineActorVOI.GetProperty().SetColor( options.outlinecolour )
	outlineActorVOI.GetProperty().SetLineWidth( options.outlinethickness )
        return [ outlineActor, outlineActorVOI ]
    else:
	return [ outlineActor ]

# Position the camera and set viewing direction
def setCamera():
    cam = ren.GetActiveCamera()
    if (options.parallelprojection):
        cam.ParallelProjectionOn()
    cam.SetViewUp(1, 0, 0)

    # Get position from system size as gotten from any vtk input file
    if (options.colouriso):
	cam.SetFocalPoint(colourisoreader.GetOutput().GetCenter())
	size=[2*coord+1 for coord in colourisoreader.GetOutput().GetCenter()]
    elif (options.colourvolume):
	cam.SetFocalPoint(colourvolumereader.GetOutput().GetCenter())
	size=[2*coord+1 for coord in colourvolumereader.GetOutput().GetCenter()]
    elif (options.rockiso):
        cam.SetFocalPoint(rockisoreader.GetOutput().GetCenter())
	size=[2*coord+1 for coord in rockisoreader.GetOutput().GetCenter()]
    elif (options.rockvolume):
        cam.SetFocalPoint(rockvolumereader.GetOutput().GetCenter())
	size=[2*coord+1 for coord in rockvolumereader.GetOutput().GetCenter()]
    # If this fails, get it from the systemsize parameter
    else:
        size= [ options.systemsize[2], options.systemsize[1], options.systemsize[0] ]
	centre = [0.5*coord for coord in size]
	cam.SetFocalPoint(centre)

    # Set camera position
    cam.SetPosition(size[0]/2,size[1]/2,8*math.sqrt(size[2]**2+size[1]**2+size[0]**2))

    # Rotate the camera position around the center point of the system
    cam.Elevation(options.elevation)
    cam.Azimuth(options.azimuth)
    cam.Roll(options.roll)
    cam.Zoom(options.zoom) 

# Make a fid from input
# Strip leading and trailing tags, so one can call the script with a complete filename
def getFid(fid):
    fid = re.sub(r"^colour_","",fid)
    fid = re.sub(r"^rock_","",fid)
    fid = re.sub(r"^md_cfg_","",fid)
    fid = re.sub(r"\.vtk$","",fid)
    fid = re.sub(r"\.asc$","",fid)
    fid = re.sub(r"\.h5$","",fid) # Also strip .h5 for convenience, but one still needs a VTK file!
    return fid

# Check if files exist and unset flags if needed
def checkFiles(fid):
    fname_colour = "colour_" + fid + ".vtk"
    fname_md = "md-cfg_" + fid + ".asc"
    fname_rock = "rock_" + fid + ".vtk"

    print "OMNIRENDER opening files for ID:", fid

    if (options.colouriso or options.colourvolume):
        try:
            file_colour = open(fname_colour,"r")
	    file_colour.close()
	except IOError:
            print "  IOError while trying to open colour file '" + fname_colour + "'. Not rendering any colour data."
	    options.colouriso = False
	    options.colourvolume = False
	else:
            if ( options.verbose): print "  Opened colour file '" + fname_colour + "'."

    if (options.particles):
        try:
            file_md = open(fname_md,"r")
	    file_md.close()
	except IOError:
            print "  IOError while trying to open MD file '" + fname_md + "'. Not rendering any particle data."
	    options.particles = False
	else:
            if ( options.verbose): print "  Opened MD file '" + fname_md + "'."

    if (options.rockiso or options.rockvolume):
        try:
            file_rock = open(fname_rock,"r")
	    file_rock.close()
	except IOError:
            print "  IOError while trying to open rock file '" + fname_rock + "'. Not rendering any rock data."
	    options.rockiso = False
	    options.rockvolume = False
        else:
            if ( options.verbose): print "  Opened rock file '" + fname_rock + "'."

    return (fname_colour, fname_md, fname_rock)

# Make an isosurface actor from the data in fname
def getIsoActor(fname, colour, opacity, value):
    isoreader = vtk.vtkStructuredPointsReader()
    isoreader.SetFileName(fname)
    isoreader.Update()

    extractiso = vtk.vtkExtractVOI()
    if (options.xmax > -1 or options.ymax > -1 or options.zmax > -1):
        extractiso.SetVOI(options.zmin,options.zmax, options.ymin,options.ymax, options.xmin,options.xmax)
    extractiso.SetInput(isoreader.GetOutput())

    iso = vtk.vtkMarchingCubes()
    iso.SetInput(extractiso.GetOutput())
    iso.SetValue(0,value)

    polymap = vtk.vtkPolyDataMapper()
    polymap.SetInput(iso.GetOutput())
    polymap.ScalarVisibilityOff()
    
    isoactor = vtk.vtkActor()
    isoactor.SetMapper(polymap)

    isoactor.GetProperty().SetColor(colour) 
    isoactor.GetProperty().SetOpacity(opacity)

    return (isoreader, isoactor)

# Make a volume volume from the data in fname
def getVolumeVolume(fname, rangelow, rangehigh, opacitylow, opacityhigh, color_transfer_func):
    volreader = vtk.vtkStructuredPointsReader()
    volreader.SetFileName(fname)
    volreader.Update()

    volrange = (rangelow, rangehigh)
    imageShift = vtk.vtkImageShiftScale()
    imageShift.SetShift(-1*volrange[0])
    imageShift.SetScale(255.0/(volrange[1]-volrange[0]))
    imageShift.SetOutputScalarTypeToUnsignedChar()
    imageShift.SetInput(volreader.GetOutput())

    extractVol = vtk.vtkExtractVOI()
    if (options.xmax > -1 or options.ymax > -1 or options.zmax > -1):
        if ( options.verbose): print "  Extracting subvolume ..."
        extractVol.SetVOI( options.zmin, options.zmax, options.ymin, options.ymax, options.xmin, options.xmax )
    extractVol.SetInput(imageShift.GetOutput())

    volume_mapper = vtk.vtkVolumeRayCastMapper()
    volume_mapper.SetInput(extractVol.GetOutput())

    composite_function = vtk.vtkVolumeRayCastCompositeFunction()
    volume_mapper.SetVolumeRayCastFunction(composite_function)

    # Setting opacity functions
    opacity_transfer_func = vtk.vtkPiecewiseFunction()
    opacity_transfer_func.AddPoint( 0, opacitylow )
    opacity_transfer_func.AddPoint( 118, opacityhigh )
    opacity_transfer_func.AddPoint( 138, opacityhigh )
    opacity_transfer_func.AddPoint( 255, opacitylow )

    # Adding colour and opacity to properties
    volume_properties = vtk.vtkVolumeProperty()
    volume_properties.SetColor( color_transfer_func )
    volume_properties.SetScalarOpacity( opacity_transfer_func )

    # Adding properties to volume
    volume = vtk.vtkVolume()
    volume.SetMapper( volume_mapper )
    volume.SetProperty( volume_properties )

    return ( volreader, volume )

# Return a VTK colorTransferFunc for the colour volume rendering
def getColourColourTransferFunc():
    colourTransferFunc = vtk.vtkColorTransferFunction()
    if (options.tuecolours == 1 ):
	colourTransferFunc.AddRGBPoint( 0, 0.0625, 0.0625, 0.44921875 )
	colourTransferFunc.AddRGBPoint( (255/2)-10, 0.0625, 0.0625, 0.44921875 )
	colourTransferFunc.AddRGBPoint( (255/2)+10, 0.8359375, 0.0, 0.2890625 )
	colourTransferFunc.AddRGBPoint( 255, 0.8359375, 0.0, 0.2890625 )
    else:
        colourTransferFunc.AddRGBPoint( 0, options.colourvolumecolourlow[0], options.colourvolumecolourlow[1], options.colourvolumecolourlow[2])
	colourTransferFunc.AddRGBPoint( (255/2)-10, options.colourvolumecolourlow[0], options.colourvolumecolourlow[1], options.colourvolumecolourlow[2])
	colourTransferFunc.AddRGBPoint( (255/2)+10, options.colourvolumecolourhigh[0], options.colourvolumecolourhigh[1], options.colourvolumecolourhigh[2])
	colourTransferFunc.AddRGBPoint( 255, options.colourvolumecolourhigh[0], options.colourvolumecolourhigh[1], options.colourvolumecolourhigh[2])
    return colourTransferFunc

# Return a VTK colorTransferFunc for the rock volume rendering
def getRockColourTransferFunc():
    colourTransferFunc = vtk.vtkColorTransferFunction()
    colourTransferFunc.AddRGBPoint( 0, options.rockvolumecolourlow[0], options.rockvolumecolourlow[1], options.rockvolumecolourlow[2] )
    colourTransferFunc.AddRGBPoint( 255, options.rockvolumecolourhigh[0], options.rockvolumecolourhigh[0], options.rockvolumecolourhigh[0] )
    return colourTransferFunc

# Split comma separated list in the option parser
def split_opt(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

# PROGRAM START

# Setting up option parser and all the options

parser = OptionParser(description='Script for rendering LB3D datasets.')

parser.add_option("-i", "--interactive", action="store_true", help="interactive mode", default=False)
parser.add_option("-m", "--addmetadata", action="store_true", help="add metadata to file (image writing mode only)", default=False)
parser.add_option("-o", "--offscreenrendering", action="store_true", help="enable offscreen rendering (image writing mode only)", default=False)
parser.add_option("-v", "--verbose", action="store_true", help="more verbose output", default=False)
parser.add_option("-V", "--version", action="store_true", help="show version information", default=False)
parser.add_option("-w", "--writeimage", action="store_true", help="write image to PNG only", default=False)

group = OptionGroup(parser, "Inclusion options")
group.add_option("--colouriso", action="store_true", help="render colour isosurface", default=False)
group.add_option("--colourvolume", action="store_true", help="render colour volume", default=False)
group.add_option("--no-outline", action="store_true", help="do not render outline", default=False)
group.add_option("--sliceoutline", action="store_true", help="render slice outline", default=False)
group.add_option("--particles", action="store_true", help="render particles", default=False)
group.add_option("--rockiso", action="store_true", help="render rock isosurface", default=False)
group.add_option("--rockvolume", action="store_true", help="render rock volume", default=False)
parser.add_option_group(group)

group = OptionGroup(parser, "Colour isosurface options")
group.add_option("--colourisocolour", help="colour of the colour isosurface", type="float", nargs=3, default=[1.0,1.0,1.0])
group.add_option("--colourisoopacity", help="opacity of the colour isosurface", type="float", default=1.0)
group.add_option("--colourisovalue", help="value of the colour isosurface", type="float", default=0.0)
parser.add_option_group(group)

group = OptionGroup(parser, "Colour volume options")
group.add_option("--colourvolumecolourhigh", help="colour parameter of the colour volume", type="float", nargs=3, default=[1.0, 0.0, 0.0])
group.add_option("--colourvolumecolourlow", help="colour parameter of the colour volume", type="float", nargs=3, default=[0.0, 0.0, 1.0])
group.add_option("--colourvolumeopacityhigh", help="opacity parameter of colour volume", type="float",default=0.01)
group.add_option("--colourvolumeopacitylow", help="opacity parameter of colour volume", type="float",default=0.01)
parser.add_option_group(group)

group = OptionGroup(parser, "Particle options")
group.add_option("--particleanisotropy", help="constant particle anisotropy for all particles",
                 type="float",default=1.0)
group.add_option("--particlepolydispersity", action="store_true",
                 help="particle half axes read from file", default=False)
group.add_option("--particlecolour", help="colour of the particles", type="float", nargs=3, default=[1.0,1.0,1.0])
group.add_option("--particlecolour_janus1", help="colour of the particles: Janus area 1", type="float", nargs=3, default=[0.0,0.0,1.0])
group.add_option("--particlecolour_janus2", help="colour of the particles: Janus area 2", type="float", nargs=3, default=[1.0,0.0,0.0])
group.add_option("--particle_janus_para", action="store_true", help="show parallel janus particle", default=False)
group.add_option("--particlecolourbyaspectratio", action="store_true",
                 help="color particles according to aspect ratio", default=False)
group.add_option("--particlecolourbyaspectratiomin",
                 help="aspect ratio associated with minimum colour value", type="float",default=0.0)
group.add_option("--particlecolourbyaspectratiomax",
                 help="aspect ratio associated with maximum colour value", type="float",default=10.0)
group.add_option("--particleidcol", help="particle id column",type="int",default=-1)
group.add_option("--particleopacity", help="particle opacity", type="float",default=1.0)
group.add_option("--particleocol", help="particle first orientation column",type="int",default=9)
group.add_option("--particleposcol", help="particle first position column",type="int",default=0)
group.add_option("--particlercol", help="particle first radius (half axes) column",type="int",
                 default=18)
group.add_option("--particleradius", help="particle radius", type="float",default=5.0)
group.add_option("--particleresolution", help="particle mesh resolution", type="int",default=32)
parser.add_option_group(group)

group = OptionGroup(parser, "Tracer options")
group.add_option("--tracercolour", help="colour of the tracers (for tracerstyle ='solid')", type="float", nargs=3, default=[0.0,0.0,0.0])
group.add_option("--tracerids", help="ids of tracer particles", type="string", action="callback", callback = split_opt)
group.add_option("--traceropacity", help="tracer opacity (for tracerstyle = 'solid')", type="float",default=1.0)
group.add_option("--tracerstyle", help="set tracer style", type="string",default="solid")
parser.add_option_group(group)

group = OptionGroup(parser, "Rock isosurface options")
group.add_option("--rockisocolour", help="colour of the rock isosurface", type="float", nargs=3, default=[1,1,1])
group.add_option("--rockisoopacity", help="opacity of the rock isosurface", type="float", default=1.0)
group.add_option("--rockisovalue", help="value of the rock isosurface", type="float", default=0.0)
parser.add_option_group(group)

group = OptionGroup(parser, "Rock volume options")
group.add_option("--rockvolumecolourhigh", help="colour parameter of the rock volume", type="float", nargs=3, default=[0.5, 0.5, 0.5])
group.add_option("--rockvolumecolourlow", help="colour parameter of the rock volume", type="float", nargs=3, default=[0.5, 0.5, 0.5])
group.add_option("--rockvolumeopacityhigh", help="opacity parameter of rock volume", type="float",default=0.01)
group.add_option("--rockvolumeopacitylow", help="opacity parameter of rock volume", type="float",default=0.01)
parser.add_option_group(group)

group = OptionGroup(parser, "Outline options")
group.add_option("--outlinecolour", help="colour of the outline", type="float", nargs=3, default=[0.0, 0.0, 0.0])
group.add_option("--outlinethickness", help="thickness of the outline", type="float", default=1.0)

group = OptionGroup(parser, "Camera options")
group.add_option('--azimuth',help="camera azimuthal rotation", type="float", default=0.0)
group.add_option('--elevation',help="camera elevation", type="float", default=0.0)
group.add_option('--roll',help="camera roll", type="float", default=0.0)
group.add_option("--systemsize", help="system size (when neither colour nor rock file is supplied)", type="int", nargs=3, default=[1,1,1])
group.add_option('--zoom',help="zoom factor", type="float", default=1.0)
group.add_option('--parallelprojection',action="store_true",help="enable parallel projection", default=False)
parser.add_option_group(group)

group = OptionGroup(parser, "Text options")
group.add_option('--text',help="text to be added", type="string", default="")
group.add_option("--textbold", action="store_true", help="render bold text", default=False)
group.add_option("--textcolour", help="colour of the text", type="float", nargs=3, default=[0.0, 0.0, 0.0])
group.add_option('--textfont',help="text font", type="string", default="Courier")
group.add_option("--textitalic", action="store_true", help="render italicized text", default=False)
group.add_option("--textposition", help="relative position of the text", type="float", nargs=2, default=[0.5, 0.9])
group.add_option('--textsize',help="font size of the text", type="int", default=24)
group.add_option("--timestamp", action="store_true", help="render a timestamp", default=False)
parser.add_option_group(group)

group = OptionGroup(parser, "Restriction options")
group.add_option('--xmin', help="render only x > xmin",type="int", default=-1)
group.add_option('--xmax', help="render only x < xmax",type="int", default=-1)
group.add_option('--ymin', help="render only y > ymin",type="int", default=-1)
group.add_option('--ymax', help="render only y < ymax",type="int", default=-1)
group.add_option('--zmin', help="render only z > zmin",type="int", default=-1)
group.add_option('--zmax', help="render only z < zmax",type="int", default=-1)
parser.add_option_group(group)

group = OptionGroup(parser, "Output window options")
group.add_option("--bgcolour", help="colour of the background", type="float", nargs=3, default=[1.0, 1.0, 1.0])
group.add_option('--hsize', help="horizontal size",type="int", default=1600)
group.add_option('--vsize', help="vertical size",type="int", default=900)
parser.add_option_group(group)

group = OptionGroup(parser, "Style options")
group.add_option("--frontview", action="store_true", help="set camera to view from the front", default=False)
group.add_option("--pickering", action="store_true", help="set options for a cubic Pickering emulsion / bijel", default=False)
group.add_option("--sideview", action="store_true", help="set camera to view from the side", default=False)
group.add_option("--topview", action="store_true", help="set camera to view from the top", default=False)
group.add_option("--tuecolours", action="store_true", help="use TU/e colour scheme", default=False)
parser.add_option_group(group)

# Parse!
(options, args) = parser.parse_args()

# Show version info and exit
if ( options.version ):
    print_versions()

# Attempt to load external modules
try:
    import vtk
except ImportError:
    print( "\nERROR importing module vtk.")
    print_versions()

try:
    import numpy as np
except ImportError:
    print( "\nERROR importing module numpy.")
    print_versions()

# Remember this for metadata
cmdlist = sys.argv[1:]
for k in args:
    cmdlist.remove(k)
cmdstr = ' '.join(cmdlist)

# Die if nothing will be done
if ( not ( options.interactive or options.writeimage ) ):
    print "Program will do nothing; specify either -i for interactive mode or -w to write a PNG file."
    exit(-1)

# Set derived values
if ( len(args) == 0 ):
    print "Please specify a filename or identifier."
    exit(-1)

fid = getFid(args[0])
ellipsoids = False
if ( options.particleanisotropy != 1.0 or options.particlepolydispersity): ellipsoids = True

particlecolourbydata = False
if ( options.particlecolourbyaspectratio): particlecolourbydata = True

# Process style options

if ( options.pickering ):
    options.colouriso = True
    options.colourvolume = True
    options.particles = True
    options.tuecolours = True
    options.colourvolumeopacitylow = 0.01
    options.colourvolumeopacityhigh = 0.01
    options.azimuth = -15.0
    options.elevation = 15.0
    options.zoom = 5.0
    options.hsize = 1920
    options.vsize = 1200

if ( options.frontview ):
    options.azimuth = 0.0
    options.elevation = 90.0
    options.roll = 0.0

if ( options.sideview ):
    options.azimuth = 90.0
    options.elevation = 0.0
    options.roll = 90.0

if ( options.topview ):
    options.azimuth = 0.0
    options.elevation = 0.0
    options.roll = 90.0

if ( options.tuecolours ):
    options.colourisocolour = [ 173.0/255.0, 32.0/255.0, 173.0/255.0 ]
    options.colourisoopacity = 0.4
    options.particlecolour =  [ 0.0, 0.671875, 0.5078125 ]
    options.particlecolour_janus1 =  [ 1.0, 0.0, 0.0 ]
    options.particlecolour_janus2 =  [ 0.0, 1.0, 0.0 ]
    options.textcolour = [ 0.0625, 0.0625, 0.44921875 ]

if ( options.timestamp ):
    # This assumes the default lb3d filename structure - grepping for 8 decimals after '_t'
    options.text = "t = " + re.search('(?<=_t)\d{8}',fid).group(0)

# Check for file existence and unset options if needed
(fname_colour, fname_md, fname_rock) = checkFiles(fid)

if ( options.verbose):
    print "Command line arguments: "
    print cmdstr
    print "Rendering with options: "
    print options

# Set up global renderer
ren = vtk.vtkRenderer()
ren.SetBackground(options.bgcolour)
    
# Set up isosurfaces
if ( options.colouriso ):
    if ( options.verbose): print "Adding colour isosurface ..."
    colourisoreader, colourisoactor = getIsoActor(fname_colour, options.colourisocolour, options.colourisoopacity, options.colourisovalue)
    ren.AddActor(colourisoactor)

if ( options.rockiso ):
    if ( options.verbose): print "Adding rock isosurface ..."
    rockisoreader, rockisoactor = getIsoActor(fname_rock, options.rockisocolour, options.rockisoopacity, options.rockisovalue)
    ren.AddActor(rockisoactor)

# Set up volume renderings
if ( options.colourvolume ): 
    if ( options.verbose): print "Adding colour volume ..."
    colourvolumereader, colourvolume = getVolumeVolume(fname_colour, -1.1, 1.1, options.colourvolumeopacitylow, options.colourvolumeopacityhigh, getColourColourTransferFunc())
    ren.AddVolume(colourvolume)

if ( options.rockvolume ): 
    if ( options.verbose): print "Adding rock volume ..."
    rockvolumereader, rockvolume = getVolumeVolume(fname_rock, 0.0, 1.0, options.rockvolumeopacitylow, options.rockvolumeopacityhigh, getRockColourTransferFunc())
    ren.AddVolume(rockvolume)

# Set up outline
if ( not options.no_outline ): 
    if ( options.verbose): print "Adding outlines ..."
    outlines = getOutlineActors()
    for outline in outlines: ren.AddActor(outline)

# Set up particles
if ( options.particles ):

    if ( options.tracerids is not None):
        # If tracers are defined...
        tids = map(int,options.tracerids)

        if ( options.tracerstyle == "solid" ):
            # Solid colour - just make one batch of the tracers and give them tracer colour and opacity
            if ( options.verbose): print "Reading tracer points (solid) ..."
            ptcldata = readParticleData(fname_md, options.particleposcol, options.particleidcol, options.particleocol, options.particlercol, tlist = tids, inverse = False )
            if ( options.verbose): print "Adding particles ..."
            ptclActor = getParticleActor(ptcldata, options.tracercolour, options.traceropacity)
            ren.AddActor(ptclActor)

            # Then render the rest of the particles with normal colour and opacity
            if ( options.verbose): print "Reading basic points (solid) ..."
            ptcldata = readParticleData(fname_md, options.particleposcol, options.particleidcol, options.particleocol, options.particlercol, tlist = tids, inverse = True )
            if ( options.verbose): print "Adding particles ..."
            ptclActor = getParticleActor(ptcldata, options.particlecolour, options.particleopacity)
            ren.AddActor(ptclActor)

        elif ( options.tracerstyle == "rainbow" ):
            # Rainbow - list of tracers getting equally spaced colours in hsv colourspace
            colours = []
            sat = 1
            value = 1
            n_tr = len(tids)
            
            # Calculate colours
            for h in range(0, n_tr + 1):
                hue = h / float(n_tr)
                colour = list(colorsys.hsv_to_rgb(hue, sat, value))
		colours.append(colour)

            # Add the tracers one by one
            i = 0
            for t in tids:
                if ( options.verbose): print "Reading tracer point (rainbow) ..."
                ptcldata = readParticleData(fname_md, options.particleposcol, options.particleidcol, options.particleocol, options.particlercol, tlist = [ t ], inverse = False )
                if ( options.verbose): print "Adding particle (rainbow) ..."
                ptclActor = getParticleActor(ptcldata, colours[i] , options.traceropacity)
                ren.AddActor(ptclActor)
                i += 1
            
            # Then render the rest of the particles with normal colour and opacity
            if ( options.verbose): print "Reading basic points (solid) ..."
            ptcldata = readParticleData(fname_md, options.particleposcol, options.particleidcol, options.particleocol, options.particlercol, tlist = tids, inverse = True )
            if ( options.verbose): print "Adding particles ..."
            ptclActor = getParticleActor(ptcldata, options.particlecolour, options.particleopacity)
            ren.AddActor(ptclActor)

        else:
            # Error!
            print "Unknown tracerstyle, aborting ..."
            exit(-1)

    elif ( options.particle_janus_para):
        # No tracers, just render Janus particles
        if ( options.verbose): print "Reading particle points ..."
        ptcldata = readParticleData(fname_md, options.particleposcol, options.particleidcol, options.particleocol, options.particlercol)
        if ( options.verbose): print "Adding particles ..."
        ptclActor = getParticleActorJanus(ptcldata, options.particlecolour_janus1, options.particleopacity, 90, 270)
        ren.AddActor(ptclActor)
        ptclActor = getParticleActorJanus(ptcldata, options.particlecolour_janus2, options.particleopacity, 0, 90)
        ren.AddActor(ptclActor)
        ptclActor = getParticleActorJanus(ptcldata, options.particlecolour_janus2, options.particleopacity, 270, 360)
        ren.AddActor(ptclActor)

    else:
        # No tracers, just render everything
        if ( options.verbose): print "Reading particle points ..."
        ptcldata = readParticleData(fname_md, options.particleposcol, options.particleidcol, options.particleocol, options.particlercol)
        if ( options.verbose): print "Adding particles ..."
        ptclActor = getParticleActor(ptcldata, options.particlecolour, options.particleopacity)
        ren.AddActor(ptclActor)
        
# Set up camera angles
if ( options.verbose): print "Setting up camera ..."
setCamera()

# Adding text
if ( options.text != ""):
    if ( options.verbose): print "Adding text ..."
    textActor = getTextActor()
    ren.AddActor(textActor)

# Start the rendering
if ( options.verbose): print "Creating render window ..."
renwin = vtk.vtkRenderWindow()
if ( options.verbose): print "  Render library: " + renwin.GetRenderLibrary()

if (options.offscreenrendering and options.writeimage):
    if ( options.verbose): print "  Setting off-screen rendering."
    renwin.OffScreenRenderingOn()

if ( options.verbose ): print "Adding renderer ..."
renwin.AddRenderer(ren)
if ( options.verbose ): print "Setting size ..."
renwin.SetSize(options.hsize, options.vsize)
if ( options.verbose ): print "Resetting clipping range ..."
ren.ResetCameraClippingRange()

# Writing image to file
if (options.writeimage):
    print "Rendering ..."
    ren.Render()
    print "Done rendering."
    writePngImage()

# Starting interactive session
if (options.interactive):
    print "Starting interactive session ..."
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renwin)
    iren.Start()

print "Done!"

