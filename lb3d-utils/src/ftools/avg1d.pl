#!/usr/bin/perl
#
# VERSION=1 :
# This script reads in a file (e.g. as supplied by avg1d.sh)
# containing "x f err c1..cn" columns, and outputs the avg of f for:
# MIN + LAG < x < MAX - LAG .
#
# This is useful to x-crop and average pxz*.X files.
#
# VERSION=0 :
# This is meant to be called directly by the user (i.e. not from
# a shell script) and direct the stdout to a file. It does same 
# as VERSION=1 but the averaging sample is:
# MAX - LAG < x < MAX
#
# This is useful to t-crop and average pxz*.Stats files.
#
# NGS Apr 2004

$VERBOSE =0;

$VERSION =0;
$infile  =$ARGV[0];
$nx      =$ARGV[1];
$lag     =$ARGV[2];
$res     =$ARGV[3];    # Resolution of first column
$avg     =0.0; $sum2 =0.0; $count = 0; 
$min     = +1e+90; $max = -1e+90;
# @f   = 0.0;  velid??
# @err = 0.0;

if($VERBOSE) {
    printf("\nReading $infile");
}

&ReadDataFile;
if($count!=0) {
    &ComputeStderr;
# &GetOutputFilename;  # This is done in the shell script
    &Output;
} else {
    print "\nNo suitable lines read in! Quitting...\n";
}


sub ReadDataFile {
    open INFILE,$infile or die "Cannot open $infile for read";
    while(<INFILE>) {
	@data = split(/\s+/, $_);
#       Careful: the data have as a first row a whitespace,
#       so discard first element, [0].
	if($data[0]=="") {
	    shift(@data);
	}
	$x  = $data[0];
	if((($VERSION==1) &&
	    (($x >= $res+$lag) && ($x <= $nx-$lag+$res)))
	    ||
	   (($VERSION==0) &&
	    (($x >= $nx-$lag+$res) && ($x <= $nx))))
	{
	    $f[$count]    = $data[1];  # array will be used later
	    $err[$count]  = $data[2];  # stderr of each f
	    $avg         += $f[$count];
	    if($f[$count] > $max) {
		$max = $f[$count];
	    }
	    if($f[$count] < $min) {
		$min=$f[$count];
	    }
	    $sum2 += ($f[$count])*($f[$count]);
	    $count++;
	    if($VERBOSE) {
		print "\navg = ",$avg,"\ncount = ",$count;
	    }
#	    print "\nCOND==";
#	    print "\nRead-in (x,f) =($x,$f)\n";
	}
    } #endwhile
    close INFILE;
    if($VERBOSE) {
	printf("\navg = %f\ncount = %d", $avg, $count);
    }
}

sub ComputeStderr {
    $avg    /=$count;
    $stderr =$sum2-($count)*($avg)*($avg);
    if($stderr <= 0) {
	print "\nStraight stderr <= 0 ! Quitting...\n"; 
	exit;
    }
    $perturb = 0.0;
    for($i=0; $i<($count--); $i++) {
	# perturbation to straight stderr
	$perturb += (sqrt(abs($f[$i] - $avg)/$stderr))*
	    ($err[$i]);             
    }
#    print "\nA1: f  = ", @f;
#    print "\nA2: avg= ", $avg;
#    print "\nA3: sqrt(diff)= ", @f;
#    print "\nB: stderr    = $stderr";    
#    print "\nC: err       = @err";        
#    print "\nperturb = sqrt(abs(A1-A2))/B*C = $perturb\n";
    $stderr /=($count-1)*($count);
    $stderr =sqrt($stderr);    
}

sub GetOutputFilename {
    $_=$infile;
    s/_t[0-9]{6}/_Ch$lag/;
    s/.X/.X.Stats/;
    $outfile=$_;
#    print $outfile;
}

sub Output {
    if($VERBOSE) {
	printf("\nSample size = %i",$count);
	printf("\n%10.8f < data < %10.8f",$min,$max);
	printf("\nAvg = %10.8f +- %10.8f (one sigma)\n",
	       $avg,$stderr);
    }
    else {
	&GetTimestep;
#	open OUTFILE,">$outfile" 
#	    or die "Cannot open $outfile for write";
	printf("%i %10.8e %10.8e $10.8e\n",
	       $ts,
	       $avg,
	       $stderr,
	       ($stderr + $perturb));
    }
}

sub GetTimestep {
    if(($VERSION==1)&&($infile=~/_t(.*)_.*.X/))  # If matching succeeds
    {
	$ts = $1;
    } elsif($VERSION==0) {
	$ts = $nx;
    }
}

