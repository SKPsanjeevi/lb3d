#!/usr/bin/perl -w

#######################
#Skript for z averaging with program $sumcut
# Frank Raischel 2/2007
#
# Permeability perm_tot in mu^2 
# check max v_z for sanity (speed of sound)
# check tau (->vis)
# lattice constant (lattice)  in mu
# samplecoord is the range that is beeing evaluated in flow direction (z presumed)
# dims FULL sample with walls and accn. reservoir if used 
###########################################

 $ftype="1";              # 1=Double, 2=Float
 $dims="70,80,90";        # Lattice Size
 $scalarType=1;           # 1=Scalar, 2=2scalar
 $samplecoord="14,77";    # z coord range of sample
 $channelwidth_xy="3,8";  # wall in pixel on each side format x_wall,y_wall
 $tau="1.";	          # tau (~viscosity) of liquid
 $lattice="20" ;          # lattice constant            


 $sumcut="./zvelocity";   # path to  zvelocity file

 #get all velocity files into the array 
 $searchstring = "find . -name \"vel_*.h5*\" -not -type d"; 
 
 #Decide for which steps flow files should be generated
 #@vlist = `$searchstring |grep 100000`;    # just or step 100000
 @vlist = `$searchstring `;                 # all steps



#----------------------------source do not change------------------------
foreach $vfile (@vlist){
  chomp($vfile);
  $timestamp =  (split /vel_/, $vfile)[1];
  #print $timestamp , "\n";
  $vfile = "vel_".$timestamp;
  $odfile = "od_".$timestamp;
  print $vfile , "\n";
  print $odfile , "\n";
  open(OUTFILE, ">/tmp/sumcut.$$");
  if($vfile=~"\.gz"){
    $cmd = " gunzip -c $vfile >vf.tmp ";
    `$cmd`;
    print OUTFILE   "vf.tmp \n";
    }
  else{
    print OUTFILE $vfile,  "\n";
    }

  if($odfile=~"\.gz"){
    $cmd = " gunzip -c $odfile >od.tmp ";
    `$cmd`;
    print OUTFILE   "od.tmp \n";
     }
  else{
    print OUTFILE $odfile,  "\n";
  }

  if($vfile=~"_100u"){
    $dims="128,128,128";
    $samplecoord="15,114";
    $channelwidth_xy="14,14";
    $lattice="40";
    }
  elsif($vfile=~"_200u"){
    $dims="256,256,256";
    $samplecoord="29,228";
    $channelwidth_xy="28,28";
    $lattice="20";
    }
  elsif($vfile=~"_400u"){
    $dims="512,512,512";
    $samplecoord="57,456";
    $channelwidth_xy="56,56";
    $lattice="10";
    }
  elsif($vfile=~"_800u"){
    $dims="800,800,800";
    $samplecoord="9,792";
    $channelwidth_xy="1,1";
    $lattice="5";
    }
  else{}


  print OUTFILE $ftype , "\n";
  print OUTFILE $scalarType , "\n";
  print OUTFILE $channelwidth_xy, "\n";
  print OUTFILE $samplecoord , "\n";
  print OUTFILE $tau , "\n";
  print OUTFILE $lattice , "\n";
  close(OUTFILE);
	$execstring = "cat /tmp/sumcut.$$ | $sumcut  \n";
	print "$execstring";
	`$execstring`;
  $copystring = "mv flow.sum flow_$timestamp.dat\n";
	`$copystring`;
  print $copystring;
  if($vfile=~"\.gz"){
    $cmd = " rm vf.tmp ";
    `$cmd`;
    }
  else{}
  if($odfile=~"\.gz"){
    $cmd = " rm od.tmp ";
    `$cmd`;
    }
  else{}
  }


exit(0);





