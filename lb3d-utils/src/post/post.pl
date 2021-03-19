#!/usr/sbin/perl

# Acme LBE output processor

$nprocs=2;

# Read input file.

open (INFILE,"input-file") or die ("Can't open input-file.");

# The inputs hash contains all of the input-file's key-value pairs.
# Note that this parser only works for one pair on each line. FIXME

%inputs=();

while(<INFILE>) {
	# Ignore the namelist specifiers and blank lines.
	if (!(/\&/) && !(/\//) && (/\S/) )
	{
		s/\s//g;	# Get rid of any whitespace.
		split(/=/);
		$inputs{@_[0]} = @_[1];
	}
}

close INFILE;

$inputs{'gr_out_file'}=~s/[\"\']//g;
$inputs{'folder'}=~s/[\"\']//g;

$gr=$inputs{'gr_out_file'};
$dir='../output/'.$inputs{'folder'};

for ($t=0;$t<=$inputs{'n_iteration'};$t+=$inputs{'n_sci'}) {
	for ($p=0;$p<$nprocs;$p++) {
		$filename=sprintf("$dir/sites_%s_t%04d_p%04d.all",$gr,$t,$p);
		print "Reading from \"$filename\".\n";
		open(FOO,$filename) or die("Can't open $filename");
		while(<FOO>){
			chomp;
			s/^\s*//; # Strip leading whitespace.
			s/\s+/ /g; # Strip multiple whitespace.
			($x,$y,$z,$a)=split(/\s+/,$_,4);
			$array[$x][$y][$z]=$a;
		}
		close FOO;
		`rm $filename`;

	}

	# Now write the array back out in AVS-friendly format.
	
	$outfilename=sprintf("$dir/sites_%s_t%04d.all",$gr,$t);
	open(OUTFILE,">$outfilename") or die("Can't open $outfilename");
	for $z (1..$inputs{'nz'}) {
	 for $y (1..$inputs{'ny'}) {
	  for $x (1..$inputs{'nx'}) {
		print OUTFILE "$x $y $z ",$array[$x][$y][$z],"\n";
	  }
	 }
	}
	close(OUTFILE);

	# Write a field file too.
	
	$ffilename=sprintf("$dir/sites_%s_t%04d.fld",$gr,$t);
	$outfilename=sprintf("sites_%s_t%04d.all",$gr,$t);
	open(OUTFILE,">$ffilename") or die("Can't open $ffilename");
	print OUTFILE "# AVS\n";
	print OUTFILE "ndim = 3\n";
	print OUTFILE "dim1 = ",$inputs{'nx'},"\n";
	print OUTFILE "dim2 = ",$inputs{'ny'},"\n";
	print OUTFILE "dim3 = ",$inputs{'nz'},"\n";
	print OUTFILE "nspace = 3\n";
	print OUTFILE "veclen = 7\n";
	print OUTFILE "data = float\n";
	print OUTFILE "field = uniform\n";
	for $i (1..7) {
		print OUTFILE "variable $i file=$outfilename filetype=ascii ";
		print OUTFILE "skip=0 offset=",$i+2," stride=10\n";
	}
	close(OUTFILE);
}

