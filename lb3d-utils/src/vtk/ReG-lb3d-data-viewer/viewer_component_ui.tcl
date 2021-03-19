# GUI for the vtk-based viewer designed for data produced by lb3d.

# SAVE button
button .button_save -text "Save jpeg" -command "button_save" \
-bd 4 -relief groove

# SAVE postscript button
button .button_save_ps -text "Save postscript"  -command "button_save_ps" \
-bd 4 -relief groove

# QUIT button
button .button_quit -text "Exit" -command "exit 0" -bd 4 -relief groove

# cropper sliders
scale .scale_xcrop_min -label "Crop x min" -from $xmin -to $xmax -orient horizontal -resolution 1 -length 7c -variable xcrop_min
scale .scale_xcrop_max -label "Crop x max" -from $xmin -to $xmax -orient horizontal -resolution 1 -length 7c -variable xcrop_max
scale .scale_ycrop_min -label "Crop y min" -from $ymin -to $ymax -orient horizontal -resolution 1 -length 7c -variable ycrop_min
scale .scale_ycrop_max -label "Crop y max" -from $ymin -to $ymax -orient horizontal -resolution 1 -length 7c -variable ycrop_max
scale .scale_zcrop_min -label "Crop z min" -from $zmin -to $zmax -orient horizontal -resolution 1 -length 7c -variable zcrop_min
scale .scale_zcrop_max -label "Crop z max" -from $zmin -to $zmax -orient horizontal -resolution 1 -length 7c -variable zcrop_max

# Sample-rate slider
scale .scale_sample -label "Hedgehog sample rate" -from 1 -to 10 -orient horizontal -resolution 1 -length 7c -variable sampleRate

# resolution of opacity/colour transfer function sliders and initial
# values

set res2 [expr ($data_max - $data_min)/256]
# set red_stop_tmp [expr ($data_max + $data_min)/2]
# Chaged to 3/4 of slider width:
set xx [expr (3*($data_max - $data_min)/4 + $data_min)]
set red_stop_tmp [expr $data_max * 0.95]
set blue_start_tmp $data_max
set opacity_start_tmp [expr $xx * 0.95]
set opacity_stop_tmp $xx
# opacity transfer function step-start slider
scale .scale_opacity_start -label "Opacity step start" \
-from $data_min -to $data_max -orient horizontal -resolution $res2 \
-length 7c -variable opacity_start_tmp

# opacity transfer function step-stop slider
scale .scale_opacity_stop -label "Opacity step stop" \
-from $data_min -to $data_max -orient horizontal -resolution $res2 \
-length 7c -variable opacity_stop_tmp

scale .scale_red_stop -label "Red stop" -from $data_min -to $data_max \
-orient horizontal -resolution $res2 -length 7c -variable red_stop_tmp

# colour transfer function blue-start slider
scale .scale_blue_start -label "Blue start" -from $data_min -to $data_max \
-orient horizontal -resolution $res2 -length 7c -variable blue_start_tmp

# iso slider
set res [expr 0.05*($data_max - $data_min)]
puts "Setting iso slider resolution to $res"

scale .scale_iso -label "Iso value" -from $data_min -to $data_max \
                 -orient horizontal -resolution $res -length 7c \
                 -variable iso_val

# Sliders for cut plane
scale .scale_xorigin -label "Origin - x" -from $xmin -to $xmax -orient horizontal -resolution 1 -length 7c -variable xorigin

scale .scale_yorigin -label "Origin - y" -from $ymin -to $ymax -orient horizontal -resolution 1 -length 7c -variable yorigin

scale .scale_zorigin -label "Origin - z" -from $zmin -to $zmax -orient horizontal -resolution 1 -length 7c -variable zorigin

scale .scale_xnormal -label "Normal - x" -from -1 -to 1 -orient horizontal -resolution 0.1 -length 7c -variable xnormal

scale .scale_ynormal -label "Normal - y" -from -1 -to 1 -orient horizontal -resolution 0.1 -length 7c -variable ynormal

scale .scale_znormal -label "Normal - z" -from -1 -to 1 -orient horizontal -resolution 0.1 -length 7c -variable znormal

grid .button_save -row 0 -column 0
grid .button_save_ps -row 1 -column 0
grid .button_quit -row 2 -column 0

grid .scale_xcrop_min -row 0 -column 1
grid .scale_xcrop_max -row 0 -column 2
grid .scale_ycrop_min -row 1 -column 1
grid .scale_ycrop_max -row 1 -column 2
grid .scale_zcrop_min -row 2 -column 1
grid .scale_zcrop_max -row 2 -column 2

set row 3
set col 1

#pack .button_save .button_quit -side left
#pack .scale_xcrop_min .scale_xcrop_max .scale_ycrop_min -side top
#pack .scale_ycrop_max .scale_zcrop_min .scale_zcrop_max -side top

if {$do_hedgehog} {
    #pack .scale_sample
    grid .scale_sample -row $row -column $col
    set row [expr $row + 1]
}

if {$do_vol == 1} {
    #pack .scale_opacity_start -side top
    #pack .scale_opacity_stop .scale_red_stop .scale_blue_start -side top
    grid .scale_opacity_start -row $row -column $col
    grid .scale_opacity_stop -row $row -column [expr $col + 1]
    set row [expr $row + 1]
    grid .scale_red_stop -row $row -column $col
    grid .scale_blue_start -row $row -column [expr $col + 1]
    set row [expr $row + 1]
}

if {$do_iso == 1} {
    #pack .scale_iso
    grid .scale_iso -row $row -column $col
    set row [expr $row + 1]
}

if {$do_cut == 1} {
    #pack .scale_xorigin -side top
    #pack .scale_yorigin .scale_zorigin -side top
    #pack .scale_xnormal .scale_ynormal .scale_znormal -side top
    grid .scale_xorigin -row $row -column $col
    grid .scale_xnormal -row $row -column [expr $col + 1]
    set row [expr $row + 1]
    grid .scale_yorigin -row $row -column $col
    grid .scale_ynormal -row $row -column [expr $col + 1]
    set row [expr $row + 1]
    grid .scale_zorigin -row $row -column $col
    grid .scale_znormal -row $row -column [expr $col + 1]
    set row [expr $row + 1]
}

# Start the procedure that repeatedly checks for the appearance of
# a new data file
if {$is_component == 1} {

  # Check for new data once every second
  set time_interval 1000
  after $time_interval [list check_file]
}

#----------------- isosurface value changer -----------------------

bind .scale_iso <ButtonRelease> {

    global renWin
    global iso
    global do_iso

    set isoVal [.scale_iso get]

    iso SetValue 0 $isoVal
    iso Modified
    iso Update

    # update the display
    puts -nonewline "Rendering..."
    renWin Render

    puts "done"
}

#----------------- cropper -----------------------------------------

proc update_cropping {} {

    global renWin
    global extractVol
    global extractVectVol
    global structPts
    global do_iso
    global do_vol
    global do_hedgehog
    global do_cut
    global xcrop_min
    global xcrop_max
    global ycrop_min
    global ycrop_max
    global zcrop_min
    global zcrop_max
    global xrange
    global yrange
    global zrange
    global sampleRate

    set xrange [expr $xcrop_max - $xcrop_min]
    set yrange [expr $ycrop_max - $ycrop_min]
    set zrange [expr $zcrop_max - $zcrop_min]

    if { $do_iso == 1 || $do_vol == 1 || $do_cut == 1 } {

	extractVol SetVOI $xcrop_min $xcrop_max $ycrop_min $ycrop_max \
	                  $zcrop_min $zcrop_max

	if { $do_iso == 1 } {

	    [iso GetInput] SetUpdateExtent 0 $xrange 0 $yrange 0 $zrange
	}
	if { $do_vol == 1 } {

	    [imageShift GetInput] SetUpdateExtent 0 $xrange 0 $yrange 0 $zrange
	}
	if { $do_cut == 1 } {

	    [cutter GetInput] SetUpdateExtent 0 $xrange 0 $yrange 0 $zrange
	}

	extractVol Modified
	extractVol Update
    }

    if { $do_hedgehog == 1 } {

	set deltax [expr $xrange/$sampleRate]
	set deltay [expr $yrange/$sampleRate]
	set deltaz [expr $zrange/$sampleRate]

	[coneGlyph GetInput] SetUpdateExtent 0 $deltax 0 $deltay 0 $deltaz
	#[hedgehog GetInput] SetUpdateExtent 0 $deltax 0 $deltay 0 $deltaz

	extractVectVol SetVOI $xcrop_min $xcrop_max $ycrop_min $ycrop_max \
                              $zcrop_min $zcrop_max

	extractVectVol Modified
	extractVectVol Update
    }

    # update the display
    puts -nonewline "Rendering..."
    renWin Render

    puts "done"
}

#----------------------

bind .scale_xcrop_min <ButtonRelease> {

    global xcrop_min
    global xcrop_max

    if {$xcrop_min > $xcrop_max} {

	set xcrop_min [expr $xcrop_max - 1]
    }

    update_cropping
}

#----------------------

bind .scale_xcrop_max <ButtonRelease> {

    global xcrop_min
    global xcrop_max

    if {$xcrop_min > $xcrop_max} {

	set xcrop_max [expr $xcrop_min + 1]
    }

    update_cropping
}

#----------------------

bind .scale_ycrop_min <ButtonRelease> {

    global ycrop_min
    global ycrop_max

    if {$ycrop_min > $ycrop_max} {

	set ycrop_min [expr $ycrop_max - 1]
    }

    update_cropping
}

#----------------------

bind .scale_ycrop_max <ButtonRelease> {

    global ycrop_min
    global ycrop_max

    if {$ycrop_min > $ycrop_max} {

	set ycrop_max [expr $ycrop_min + 1]
    }

    update_cropping
}

#----------------------

bind .scale_zcrop_min <ButtonRelease> {

    global zcrop_min
    global zcrop_max

    if {$zcrop_min > $zcrop_max} {

	set zcrop_min [expr $zcrop_max - 1]
    }

    update_cropping
}

#----------------------

bind .scale_zcrop_max <ButtonRelease> {

    global zcrop_min
    global zcrop_max

    if {$zcrop_min > $zcrop_max} {

	set zcrop_max [expr $zcrop_min + 1]
    }

    update_cropping
}

#---------------- Sampling-rate changer ---------------------------

bind .scale_sample <ButtonRelease> {

    global extractVectVol
    global coneGlyph
    global sampleRate
    global xmin xmax
    global renWin
    global xrange
    global yrange
    global zrange

    # Work out what the update extent will be so we can set it
    # for the vtkGlyph3D module downstream of the sub-sampler
    set deltax [expr $xrange/$sampleRate]
    set deltay [expr $yrange/$sampleRate]
    set deltaz [expr $zrange/$sampleRate]

    [coneGlyph GetInput] SetUpdateExtent 0 $deltax 0 $deltay 0 $deltaz

    # Only now do we actually change the sample rate
    extractVectVol SetSampleRate $sampleRate $sampleRate $sampleRate
    extractVectVol Modified
    extractVectVol Update

    renWin Render
}

#-------------- opacity step-start value changer -------------------

bind .scale_opacity_start <ButtonRelease> {

  global opacityTFunc
  global renWin
  global old_opacity_step_start
  global opacity_step_start
  global opacity_step_stop
  global opacity_start_tmp
  global opacity_stop_tmp

  set opacity_step_stop [expr (($opacity_stop_tmp - $data_min)*255)/ \
      ($data_max-$data_min)]
  set opacity_step_start [expr (($opacity_start_tmp - $data_min)*255)/ \
      ($data_max-$data_min)]

  if {$opacity_step_stop < $opacity_step_start} {

    set opacity_step_start [expr $opacity_step_stop - 1]
    set opacity_start_tmp [expr $opacity_stop_tmp - $res]

      if {$opacity_step_start < 0} {

        set opacity_step_start 0
        set opacity_step_stop  1
      }
  }

  opacityTFunc RemovePoint $old_opacity_step_start
  opacityTFunc AddPoint $opacity_step_start 0.0
#  opacityTFunc AddPoint $opacity_step_start 1.0
  set old_opacity_step_start $opacity_step_start
  
  opacityTFunc Modified
  opacityTFunc Update
  
  # update the display
  puts -nonewline "Rendering..."
  renWin Render
  puts "done"
}

#--------------- opacity step-stop value changer --------------------

bind .scale_opacity_stop <ButtonRelease> {

  global opacityTFunc
  global renWin
  global old_opacity_step_start
  global opacity_step_start
  global opacity_step_stop
  global opacity_start_tmp
  global opacity_stop_tmp
  
  set opacity_step_stop [expr (($opacity_stop_tmp - $data_min)*255)/ \
      ($data_max-$data_min)]
  set opacity_step_start [expr (($opacity_start_tmp - $data_min)*255)/ \
      ($data_max-$data_min)]

  if {$opacity_step_stop < $opacity_step_start} {

      set opacity_step_stop [expr $opacity_step_start + 1]
      set opacity_stop_tmp [expr $opacity_start_tmp + $res]

      if {$opacity_step_stop > 255} {
        set opacity_step_stop 255
        set opacity_step_start 254
      }
  }

  opacityTFunc RemovePoint $old_opacity_step_stop
  opacityTFunc AddPoint $opacity_step_stop 1.0
#  opacityTFunc AddPoint $opacity_step_stop 0.0
  set old_opacity_step_stop $opacity_step_stop
  
  opacityTFunc Modified
  opacityTFunc Update
  
  # update the display
  puts -nonewline "Rendering..."
  renWin Render
  
  puts "done"
}

#------------- colour red-stop value changer --------------------

bind .scale_red_stop <ButtonRelease> {

  global colorTFunc
  global volume
  global renWin
  global old_red_stop
  global red_stop
  global red_stop_tmp
  
  set red_stop [expr (($red_stop_tmp-$data_min)*255)/($data_max-$data_min)]

  if { $red_stop > 255 } { set red_stop 255 }
  if { $red_stop < 0 } { set red_stop 0 }
   
  colorTFunc RemovePoint $old_red_stop
  colorTFunc AddRGBPoint $red_stop 1 0 0
  set old_red_stop $red_stop

  volume Modified
  volume Update
  
  # update the display
  puts -nonewline "Rendering..."
  renWin Render
  puts "done"
}

#------------- colour blue-start value changer --------------------

bind .scale_blue_start <ButtonRelease> {

  global colorTFunc
  global volume
  global renWin
  global old_blue_stop
  global blue_stop
  global blue_stop_tmp
 
  set blue_start [expr (($blue_start_tmp-$data_min)*255)/($data_max-$data_min)] 

  if { $blue_start > 255 } { set blue_start 255 }
  if { $blue_start < 0 } { set blue_start 0 }

  colorTFunc RemovePoint $old_blue_start
  colorTFunc AddRGBPoint $blue_start 0 0 1
  set old_blue_start $blue_start

  volume Modified
  volume Update
  
  # update the display
  puts -nonewline "Rendering..."
  renWin Render
  puts "done"
}

#------------- Event handlers for crop-plane sliders -------------

bind .scale_xorigin <ButtonRelease> {


  update_slicer
}

bind .scale_yorigin <ButtonRelease> {


  update_slicer
}

bind .scale_zorigin <ButtonRelease> {


  update_slicer
}

proc update_slicer {} {
  global renWin
  global plane
  global cutter
  global xorigin
  global yorigin
  global zorigin
  global xnormal
  global ynormal
  global znormal

  plane SetOrigin $xorigin $yorigin $zorigin
  plane SetNormal $xnormal $ynormal $znormal

  cutter Modified
  cutter Update

  # update the display
  puts -nonewline "Rendering..."
  renWin Render
  puts "done"
}

#------------- Plane normal editor --------------------------------

bind .scale_xnormal <ButtonRelease> {

  update_slicer
}

bind .scale_ynormal <ButtonRelease> {

  update_slicer
}

bind .scale_znormal <ButtonRelease> {

  update_slicer
}

#------------- SAVE BUTTON procedure -----------------------------

proc button_save {} {

  global save_image_filename_prefix
  global save_image_sequence_number

  toImage Modified

  imageWriter SetFileName $save_image_filename_prefix.$save_image_sequence_number.jpg
  imageWriter Write

  set save_image_sequence_number [expr $save_image_sequence_number + 1 ]
	
  .button_save configure -text "Save jpeg $save_image_sequence_number"
}

#------------- Save postscript procedure --------------------------

proc button_save_ps {} {

  global save_ps_filename_prefix
  global save_ps_sequence_number

  toImage Modified

  psWriter SetFileName $save_ps_filename_prefix.$save_ps_sequence_number.ps
  psWriter Write

  set save_ps_sequence_number [expr $save_ps_sequence_number + 1 ]
	
  .button_save_ps configure -text "Save postscript $save_ps_sequence_number"
}

#------------- Update gui and re-render ----------------------------

proc gui_update {} {

  global structPts
  global renWin
  global data_min
  global data_max
  global imageShift
  global do_iso
  global do_vol
  global do_hedgehog
  global do_cut
  global iso_val

# Update range of slider
  if { $do_iso == 1 || $do_vol == 1 || $do_cut} {

      set range [$structPts GetScalarRange]
      set data_min [lindex $range 0]
      set data_max [lindex $range 1]

      puts "Data min. = $data_min, data max. = $data_max"
  }

  if { $do_iso == 1 } {

      .scale_iso configure -from $data_min
      .scale_iso configure -to $data_max
      set res [expr 0.05*($data_max - $data_min)]
      .scale_iso configure -resolution $res

      # Update iso value if necessary
      if { $iso_val < $data_min || $iso_val > $data_max } {
	  set iso_val [expr 0.5*($data_max - $data_min)]
      }
  }

# Update the shift and scale for the ImageShiftScale filter
  if { $do_vol == 1 } {
      set scale [expr (255.0/($data_max - $data_min))]
      imageShift SetShift [expr -1*$data_min]
      imageShift SetScale $scale
      imageShift Modified
      imageShift Update

      set opacity_step_start 120
      set old_opacity_step_start 120
      set opacity_step_stop  130
      set old_opacity_step_stop  130

      set red_stop 120
      set old_red_stop 120
      set blue_start 130
      set old_blue_start 130
  }

  renWin Render
}

#------------- Procedure to check for new data file --------------

if {$is_component == 1} {

proc check_file {} {

  global reader
  global vectReader
  global input_file
  global input_lock_file
  global vector_file
  global vector_lock_file
  global do_iso
  global do_vol
  global do_hedgehog
  global do_cut
  global time_interval
  global extractVol
  global extractVectVol

  if { $do_iso == 1 || $do_vol == 1 || $do_cut == 1 } {

      # We look for presence of lock file rather than data file
      # in order to avoid (hopefully) a race condition
      if {[file isfile $input_lock_file]} {

	  reader SetFileName $input_file
	  reader OpenVTKFile
	  reader Modified
	  reader Update
	  reader CloseVTKFile

	  # 'consume' this file...
	  file delete $input_lock_file
	  file delete $input_file

          extractVol Modified
          extractVol Update

	  gui_update
      }
  }

  if { $do_hedgehog == 1 } {

      # We look for presence of lock file rather than data file
      # in order to avoid (hopefully) a race condition
      if {[file isfile $vector_lock_file]} {

	  vectReader SetFileName $vector_file
	  vectReader OpenVTKFile
	  vectReader Modified
	  vectReader Update
	  vectReader CloseVTKFile
	  # 'consume' this file...
	  file delete $vector_lock_file
	  file delete $vector_file

          extractVectVol Modified
          extractVectVol Update

	  gui_update
      }
  }

# ...and round we go again...
  after $time_interval [list check_file]
}

}
