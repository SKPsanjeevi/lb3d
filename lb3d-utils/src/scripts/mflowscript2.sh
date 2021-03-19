#!/usr/bin/perl -w


#######################
#Skript for z averaging with program $sumcut
# Frank Raischel 2/2007
#######################

 $ftype="1";          # 1=Double, 2=Float
 $dims="32,32,64";  # Lattice Size
 $scalarType=1;           # 1=Scalar, 2=2scalar
 $vdir="Y";         # Variable Directions
# Ab="2,256"       # Fixed Points In Other Directions
 $samplecoord="17,81"; #z coord range of sample
 $tau="1.0";	  # tau (~viscosity) of liquid
 $lattice="31.25" ; # lattice constant = real size of sample / z coordinates          



$sumcut="/data/data0/raischel/Porous/version4/lbe/utils/ftools/zvelocity";

$searchstring = "find . -name \"vel_v_*.h5\" -type f";

@vlist = `$searchstring`;

foreach $vfile (@vlist){
  chomp($vfile);
  $timestamp =  (split /v_/, $vfile)[1];
  #print $timestamp , "\n";
  $vfile = "vel_v_".$timestamp;
	$odfile = "od_v_".$timestamp;
  print $vfile , "\n";
  print $odfile , "\n";
  open(OUTFILE, ">/tmp/sumcut.$$");
  print OUTFILE $vfile,  "\n";
  print OUTFILE $odfile,  "\n";
  print OUTFILE $ftype , "\n";
  print OUTFILE $scalarType , "\n";
  print OUTFILE $samplecoord , "\n";
  print OUTFILE $tau , "\n";
  print OUTFILE $lattice , "\n";
  close(OUTFILE);
  #$execstring = "cat /tmp/sumcut.$$ | $sumcut >>/dev/null \n";
	$execstring = "cat /tmp/sumcut.$$ | $sumcut  \n";
	print "$execstring";
	`$execstring`;
  `rm /tmp/sumcut.$$`;
  $copystring = "mv flow.sum flow_$timestamp.dat\n";
	`$copystring`;
  print $copystring;
  }


exit(0);





