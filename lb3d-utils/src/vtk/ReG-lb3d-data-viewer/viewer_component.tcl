# Volume renderer mainly designed for 'colour' data produced by 
# lb3d in that it imitates an isosurface at a value of zero.
# Expects scalar data as input and converts it to byte data for 
# volume rendering.
#
# A Porter, 29.10.2002
catch {load vtktcl}

package require vtk 

# Set this to zero for 'viewer' behaviour
set is_component 0

#---------------------------------------------------------------------

proc make_vol {num_px data} {

    global ren1
    global data_min
    global data_max
    global red_stop
    global old_red_stop
    global blue_start
    global old_blue_start
    global opacity_step_start
    global opacity_step_stop

    # Set the range for the slider - I think GetScalarRange
    # returns a list.
    set range [$data GetScalarRange]
    set data_min [lindex $range 0]
    set data_max [lindex $range 1]

    puts "Data min. = $data_min, data max. = $data_max"

    vtkImageShiftScale imageShift
      imageShift SetShift [expr -1*$data_min]
      imageShift SetScale [expr 255.0/($data_max - $data_min)]
      imageShift SetOutputScalarTypeToUnsignedChar
      imageShift SetInput $data


    # Opacity transfer function
#    vtkPiecewiseFunction opacityTFuncBAK
#      opacityTFunc AddPoint 0 0.0
#      opacityTFunc AddPoint $opacity_step_start 0.0
#      opacityTFunc AddPoint $opacity_step_stop  0.99
#      opacityTFunc AddPoint 255 0.99

    vtkPiecewiseFunction opacityTFunc
      opacityTFunc AddPoint 0 0.0
      opacityTFunc AddPoint $opacity_step_start 0.0
      opacityTFunc AddPoint $opacity_step_stop  1.0
      opacityTFunc AddPoint 255 1.0

    vtkColorTransferFunction colorTFunc
      colorTFunc AddRGBPoint 0.0         1 0 0
      colorTFunc AddRGBPoint $red_stop   1 0 0
      colorTFunc AddRGBPoint $blue_start 0 0 1
      colorTFunc AddRGBPoint 255         0 0 1

    # Create the properties of the volume
    vtkVolumeProperty volumeProperty
      volumeProperty SetColor colorTFunc 
      volumeProperty SetScalarOpacity opacityTFunc
      volumeProperty ShadeOn

    # Use tri-linear rather than NN interpolation to improve image quality
    volumeProperty SetInterpolationTypeToLinear

    # Create objects specific to volumetric ray casting
    vtkVolumeRayCastCompositeFunction compositeFunction
  
    vtkVolumeRayCastMapper volumeMapper
      volumeMapper SetVolumeRayCastFunction compositeFunction
  
    # Parameters to improve resolution
    #  set dist volumeMapper GetSampleDistance
    volumeMapper SetSampleDistance 0.1

    # Safety check on number of px's requested
    if { $num_px < 1 } { 
	set num_px 1
    } elseif { $num_px > 16} {
	set num_px 16
    }
    volumeMapper SetNumberOfThreads $num_px
    set nThreads [volumeMapper GetNumberOfThreads]
    puts "VolumeRayCastMapper is using $nThreads threads"

    # Set input of volume mapper
    volumeMapper SetInput [imageShift GetOutput]

    # Create the volume
    vtkVolume volume
    volume SetMapper volumeMapper
    volume SetProperty volumeProperty

    # Add the volume mapper to the scene
    ren1 AddProp volume
}

#---------------------------------------------------------------------

proc make_iso { data } {

    global ren1
    global iso
    global iso_val
    global data_min
    global data_max

    # Set the range for the isosurface slider - I think GetScalarRange
    # returns a list.
    set range [$data GetScalarRange]
    set data_min [lindex $range 0]
    set data_max [lindex $range 1]

    puts "make_iso: iso_min = $data_min and iso_max = $data_max"

    # Set the initial value used to generate the isosurface
    set iso_val [expr 0.5*($data_max + $data_min)]
 #   set iso_val 0.40

    vtkMarchingCubes iso
      iso SetInput $data
      iso SetValue 0 $iso_val

    vtkPolyDataMapper poly_map
      poly_map SetInput [iso GetOutput]
      poly_map ScalarVisibilityOff

    vtkActor isoActor
      isoActor SetMapper poly_map
      eval [isoActor GetProperty] SetColor 0 0 1

    ren1 AddActor isoActor
}

#---------------------------------------------------------------------

proc make_glyphs {data} {

    global ren1

    vtkConeSource cone
      cone SetResolution 4
      cone SetAngle 30
      cone SetHeight 6

    vtkTransform transform
      transform Translate 1.5 0.0 0.0

    vtkTransformPolyDataFilter transformF
      transformF SetInput [cone GetOutput]
      transformF SetTransform transform

    vtkGlyph3D coneGlyph
      coneGlyph SetInput $data
      coneGlyph SetSource [transformF GetOutput]
      coneGlyph SetScaleFactor 0.25
      coneGlyph SetScaleModeToScaleByVector
      coneGlyph SetColorModeToColorByVector

    vtkPolyDataMapper coneMapper
      coneMapper SetInput [coneGlyph GetOutput]
      eval coneMapper SetScalarRange [[coneGlyph GetOutput] GetScalarRange]

    vtkActor coneActor
      coneActor SetMapper coneMapper
      [coneActor GetProperty] SetOpacity 0.99

    ren1 AddActor coneActor
}

#---------------------------------------------------------------------

proc make_hedgehog {data} {

    global ren1

    vtkHedgeHog hedgehog
      hedgehog SetInput $data
      hedgehog SetScaleFactor 2.0
      hedgehog SetVectorModeToUseVector

    vtkPolyDataMapper sgridMapper
      sgridMapper SetInput [hedgehog GetOutput]

    vtkActor sgridActor
      sgridActor SetMapper sgridMapper
      [sgridActor GetProperty] SetColor 1 0 0

    ren1 AddActor sgridActor
}

#---------------------------------------------------------------------

proc make_cut_plane {data} {

    global ren1
    global cutter
    global xnormal
    global ynormal
    global znormal
    global xorigin
    global yorigin
    global zorigin
    global xmax
    global xmin
    global ymax
    global ymin
    global zmax
    global zmin

    set range [$data GetScalarRange]
    set data_min [lindex $range 0]
    set data_max [lindex $range 1]

    # The plane that defines the cutting plane
    set xnormal 1
    set ynormal 0
    set znormal 0

    set xorigin [expr 0.5*($xmax+$xmin)]
    set yorigin [expr 0.5*($ymax+$ymin)]
    set zorigin [expr 0.5*($zmax+$zmin)]

    vtkPlane plane
      plane SetOrigin $xorigin $yorigin $zorigin
      plane SetNormal $xnormal $ynormal $znormal

    vtkCutter cutter
      cutter SetInput $data
      cutter SetCutFunction plane

    vtkLookupTable clut
      clut SetHueRange 0 0.67
      clut Build

    vtkPolyDataMapper cutterMapper
      cutterMapper SetInput [cutter GetOutput]
      cutterMapper SetScalarRange $data_min $data_max
      cutterMapper SetLookupTable clut

    vtkActor cut
      cut SetMapper cutterMapper

    ren1 AddActor cut
}

#---------------------------------------------------------------------

proc make_outline {data} {

    global ren1

    vtkOutlineFilter outline
      outline SetInput $data

    vtkPolyDataMapper outlineMapper
      outlineMapper SetInput [outline GetOutput]

    vtkActor outlineActor
      outlineActor SetMapper outlineMapper
      set outlineProp [outlineActor GetProperty]
      eval $outlineProp SetColor 0 0 0

    # Add the actors to the renderer, set the background and size
    #
    ren1 AddActor outlineActor
}

#---------------------------------------------------------------------

proc make_axes {data} {

    global ren1
    global xcenter ycenter zcenter xmax
    set xAxis x 
    set yAxis y
    set zAxis z
 
    set fscale [expr 8]
    set fdist [expr -$xmax/12]
    # Create axes.
    reader Update
    set bounds [[reader GetOutput] GetBounds]
    vtkAxes axes

    axes SetOrigin $fdist $fdist $fdist
    axes SetScaleFactor [expr [[reader GetOutput] GetLength]/($xmax / 5)]

	set xcenter [expr [[reader GetOutput] GetLength]/($xmax / 5)]
	set ycenter [expr [[reader GetOutput] GetLength]/($xmax / 5)]
	set zcenter [expr [[reader GetOutput] GetLength]/($xmax / 5)]

    vtkTubeFilter axesTubes
    axesTubes SetInput [axes GetOutput]
    axesTubes SetRadius [expr [axes GetScaleFactor]/35.0]
    axesTubes SetNumberOfSides 6
    
    vtkPolyDataMapper axesMapper
    axesMapper SetInput [axesTubes GetOutput]

    vtkActor axesActor
    axesActor SetMapper axesMapper

    # Label the axes.
    vtkVectorText XText
    XText SetText $xAxis 
#$xAxis

    vtkPolyDataMapper XTextMapper
    XTextMapper SetInput [XText GetOutput]

    vtkFollower XActor
    XActor SetMapper XTextMapper
    XActor SetScale $fscale $fscale $fscale
    XActor SetPosition $xcenter [expr 3*$fdist] [expr 3*$fdist]
    [XActor GetProperty] SetColor 0 0 0

    XActor SetCamera [ren1 GetActiveCamera]

    vtkVectorText YText
    YText SetText $yAxis

    vtkPolyDataMapper YTextMapper
    YTextMapper SetInput [YText GetOutput]

    vtkFollower YActor
    YActor SetMapper YTextMapper
    YActor SetScale $fscale $fscale $fscale
    YActor SetPosition [expr 3*$fdist] $ycenter [expr 3*$fdist]
    [YActor GetProperty] SetColor 0 0 0
    YActor SetCamera [ren1 GetActiveCamera]
    
    vtkVectorText ZText
    ZText SetText $zAxis 

    vtkPolyDataMapper ZTextMapper
    ZTextMapper SetInput [ZText GetOutput]
    
    vtkFollower ZActor
    ZActor SetMapper ZTextMapper
    ZActor SetScale $fscale $fscale $fscale
    ZActor SetPosition [expr 3*$fdist] [expr 3*$fdist] $zcenter
    [ZActor GetProperty] SetColor 0 0 0
    ZActor SetCamera [ren1 GetActiveCamera]

    # Add the actors to the renderer
    #
    ren1 AddActor XActor
    ren1 AddActor YActor
    ren1 AddActor ZActor
    ren1 AddActor axesActor 
}

#---------------------------------------------------------------------
proc make_initrot {} {
    global azimuth
    global elevation
    [ren1 GetActiveCamera] Azimuth  $azimuth
    [ren1 GetActiveCamera] Elevation $elevation
    [ren1 GetActiveCamera] SetParallelProjection 0 
    [ren1 GetActiveCamera] Modified

}
#---------------------------------------------------------------------

# To avoid default conversion of strings to lowercase
# (Cf. http://dbwww.essc.psu.edu/lasdoc/tae/programmer/tclp_3_2.html):

#DISABLE-FORCE_LOWER
#disable-force_lower

set do_vol 0
set do_iso 0
set do_hedgehog 0
set do_cut 0
set do_initrot 1

#Following 2 are ideal for CONIC PROJECTION
 set azimuth 90
 set elevation 270
# Following 2 values are ideal for the shearing study:
#set azimuth 90
#set elevation -90

set opacity_step_start 120
set old_opacity_step_start 120
set opacity_step_stop 130
set old_opacity_step_stop 130

set red_stop   120
set old_red_stop   120
set blue_start 130
set old_blue_start 130

set iso_val 0.0
set data_min -1.0
set data_max +1.0

set xrange 0
set yrange 0
set zrange 0

if {$argc < 1} {
  puts "Usage: vtk viewer_component.tcl  <gui_file> <num_px> <method: iso, vol, cut or hh> <structured_points_file> <vector_data_file>"
  exit
}

# Location of the file defining the gui
set gui_file [lindex $argv 0]

if {![file isfile $gui_file]} {

  puts "Cannot access gui file <$gui_file>"
  exit
}

# Number of processors/processes to use
set num_px [lindex $argv 1]

# Type of visualisation to do
set vis_type [lindex $argv 2]

if {$vis_type == "iso"} {

    set do_iso 1

    if {$argc != 4} {

	puts "No data-file name specified"
	exit
    }

    set input_file [lindex $argv 3]
    set input_lock_file $input_file
    append input_lock_file ".lock"

    puts "Making isosurface..."

} elseif {$vis_type == "vol"} {

    set do_vol 1

    if {$argc != 4} {

	puts "No data-file name specified"
	exit
    }

    set input_file [lindex $argv 3]
    set input_lock_file $input_file
    append input_lock_file ".lock"

    puts "Volume rendering..."

} elseif {$vis_type == "cut"} {

    set do_cut 1

    if {$argc != 4} {

	puts "No data-file name specified"
	exit
    }

    set input_file [lindex $argv 3]
    set input_lock_file $input_file
    append input_lock_file ".lock"

    puts "Generating cut-place..."

} elseif {$vis_type == "hh"} {

    set do_hedgehog 1

    if {$argc == 5} {

	# Have both density and vector data...
	puts "Generating isosurface & hedgehog..."
	set do_iso 1

	set input_file [lindex $argv 3]
	set input_lock_file $input_file
	append input_lock_file ".lock"

	set vector_file [lindex $argv 4]
	set vector_lock_file $vector_file
	append vector_lock_file ".lock"

    } elseif {$argc == 4} {

	# Only have vector data...
	puts "Generating hedgehog..."
	set vector_file [lindex $argv 3]
	puts "Vector file: $vector_file"
	set vector_lock_file $vector_file
	append vector_lock_file ".lock"

    } else {

	puts "Hedgehog method needs either either a file containing vector data or one file with density and another with vectors..."
	exit
    }
} else {

  puts "Valid options are: iso, vol, cut or hh (for hedgehog)"
  exit
}

# Create the RenderWindow, Renderer and interactor
vtkRenderer ren1
  # Off-white
  ren1 SetBackground 1 1 1 
  # Black
  #ren1 SetBackground 0 0 0 

vtkRenderWindow renWin
  renWin AddRenderer ren1
 # renWin SetSize 800 300
  renWin SetSize 480 480

# NGS:
if { $do_iso == 1 || $do_vol == 1 || $do_cut == 1} {
    renWin SetWindowName $input_file
} elseif { $do_hedgehog == 1 } {
    renWin SetWindowName $vector_file    
}

vtkRenderWindowInteractor iren
  iren SetRenderWindow renWin

# Necessary bits to save image to file
vtkWindowToImageFilter toImage
  toImage SetInput renWin

vtkJPEGWriter imageWriter
  imageWriter SetInput [toImage GetOutput] 
  imageWriter SetQuality 80

vtkPostScriptWriter psWriter
  psWriter SetInput [toImage GetOutput]

if { $do_iso == 1 || $do_vol == 1 || $do_cut == 1} {
    renWin SetWindowName $input_file
} elseif { $do_hedgehog == 1 } {
    renWin SetWindowName $vector_file    
}
#set save_image_filename_prefix $input_file
if { $do_iso == 1 || $do_vol == 1 || $do_cut == 1} {
    set save_image_filename_prefix $input_file
} elseif { $do_hedgehog == 1 } {
    set save_image_filename_prefix $vector_file    
}

set save_image_sequence_number 0

#set save_ps_filename_prefix $input_file
if { $do_iso == 1 || $do_vol == 1 || $do_cut == 1} {
    set save_ps_filename_prefix $input_file
} elseif { $do_hedgehog == 1 } {
    set save_ps_filename_prefix $vector_file
}

set save_ps_sequence_number 0

# Reader for structured-points data
vtkStructuredPoints structPts

if { $do_iso == 1 || $do_vol == 1 || $do_cut == 1} {

    if {![file isfile $input_file]} {

	puts "Cannot access data file <$input_file>"
	exit
    }

    vtkStructuredPointsReader reader
      reader SetFileName $input_file

    set valid [reader IsFileStructuredPoints]
    if { $valid == 0 } {

	puts "Data file >$input_file< invalid...quitting"
	exit
    }

    reader Update
    reader CloseVTKFile

    # Delete file once read - file is 'consumed'
    if {$is_component == 1} {
      file delete $input_lock_file
      file delete $input_file
    }

    set bounds [[reader GetOutput] GetBounds]
    set xmin [expr {round([lindex $bounds 0])}]
    set xmax [expr {round([lindex $bounds 1])}]
    set ymin [expr {round([lindex $bounds 2])}]
    set ymax [expr {round([lindex $bounds 3])}]
  #  set xmin 0
  #  set xmax 32
  #  set ymin 0
  #  set ymax 32
    set zmin [expr {round([lindex $bounds 4])}]
    set zmax [expr {round([lindex $bounds 5])}]

    puts "xmin = $xmin, xmax = $xmax"
    puts "ymin = $ymin, ymax = $ymax"
    puts "zmin = $zmin, zmax = $zmax"

    set xcrop_min $xmin
    set xcrop_max $xmax
    set ycrop_min $ymin
    set ycrop_max $ymax
    set zcrop_min $zmin
    set zcrop_max $zmax

    set xrange [expr $xcrop_max - $xcrop_min]
    set yrange [expr $ycrop_max - $ycrop_min]
    set zrange [expr $zcrop_max - $zcrop_min]
    set xcenter [expr ($xcrop_max + $xcrop_min)/2]
    set ycenter [expr ($ycrop_max + $ycrop_min)/2]
    set zcenter [expr ($zcrop_max + $zcrop_min)/2]

    vtkExtractVOI extractVol
      extractVol SetVOI $xcrop_min $xcrop_max $ycrop_min $ycrop_max \
                        $zcrop_min $zcrop_max
      extractVol SetInput [reader GetOutput]
      extractVol Update

    set structPts [extractVol GetOutput]
}

if { $do_vol == 1 } {
    make_vol $num_px $structPts
}

if { $do_iso == 1 } {
    make_iso $structPts
}

if { $do_cut == 1 } {
    make_cut_plane $structPts
}

if { $do_hedgehog == 1 } {

    if {![file isfile $vector_file]} {

	puts "Cannot access data file <$vector_file>"
	exit
    }

    # Reader for structured-points vector data
    vtkAVSStructuredPointsReader vectReader
    vectReader Treat3ComptAsVector
    vectReader SetFileName $vector_file

    set valid [vectReader IsFileStructuredPoints]
    if { $valid == 0 } {

	puts "Data file >$vector_file< invalid...quitting"
	exit
    }
    vectReader Update
    vectReader CloseVTKFile

    # Delete file once read - file is 'consumed'
    if {$is_component == 1} {
      file delete $vector_lock_file
      file delete $vector_file
    }

    if { $do_iso != 1 && $do_vol != 1 } {
	# Have no density file so need to get bounds from the 
	# vector data file
	set bounds [[vectReader GetOutput] GetBounds]
	set [expr {round(xmin [lindex $bounds 0])}]
	set [expr {round(xmax [lindex $bounds 1])}]
	set [expr {round(ymin [lindex $bounds 2])}]
	set [expr {round(ymax [lindex $bounds 3])}]
	set [expr {round(zmin [lindex $bounds 4])}]
	set [expr {round(zmax [lindex $bounds 5])}]

	puts "xmin = $xmin, xmax = $xmax"
	puts "ymin = $ymin, ymax = $ymax"
	puts "zmin = $zmin, zmax = $zmax"

	set xcrop_min $xmin
	set xcrop_max $xmax
	set ycrop_min $ymin
	set ycrop_max $ymax
	set zcrop_min $zmin
	set zcrop_max $zmax

	set xrange [expr $xcrop_max - $xcrop_min]
	set yrange [expr $ycrop_max - $ycrop_min]
	set zrange [expr $zcrop_max - $zcrop_min]
    }

    set sampleRate 4

    vtkExtractVOI extractVectVol
    extractVectVol SetVOI $xcrop_min $xcrop_max $ycrop_min $ycrop_max \
                          $zcrop_min $zcrop_max
    extractVectVol SetInput [vectReader GetOutput]
    # Sub-sample as well as crop
    extractVectVol SetSampleRate $sampleRate $sampleRate $sampleRate

    vtkStructuredPoints vectPts
    set vectPts [extractVectVol GetOutput]

    #make_hedgehog $vectPts
    make_glyphs $vectPts

    if { $do_iso != 1 && $do_vol != 1} {

	make_outline $vectPts
	make_axes $vectPts
    }
}

if { $do_iso == 1 || $do_vol == 1 || $do_cut == 1 } {
        make_outline $structPts
	#####make_axes $structPts
}

if { $do_initrot == 1 } {
   make_initrot 
}
iren Initialize

# render the image
#
#iren SetUserMethod {wm deiconify .vtkInteract}

# source the user interface
source $gui_file
