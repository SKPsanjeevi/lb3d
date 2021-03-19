#!/bin/bash
export LD_LIBRARY_PATH=/usr/local64/stow/VTK-5.0.3/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local64/stow/VTK-5.0.3/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/usr/local64/stow/VTK-5.0.3/lib64/python2.5/site-packages:$PYTHONPATH




i=100;
while [ $i -lt 1100 ];do
echo spheres-00000$i.vtk.noheadtab;
echo flow.$i.vtk.newheader; 
echo Running Python Script:;
./showspheres.py /data/data1/mnazruli/nobackup/spheres-00000$i.vtk.noheadtab /data/data1/mnazruli/nobackup/flow.$i.vtk.newheader;
mv Streamline.$ii.png /data/data1/mnazruli/nobackup/Streamline.$ii.png; 
let i=$i+100; 
done

