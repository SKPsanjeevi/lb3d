#!/usr/bin/perl

#=DESCRIPTION===========================================================
#
# fortran-namelist-checker.pl
#
# Checks Fortran namelists against namelist definitions read from
# Fortran source code.
#
# version:  2011-06-21
# website:  http://fortran-namelist-checker.noschek.de
# author:   Florian Janoschek (fortran-namelist-checker@noschek.de)
# requires: Fortran::F90Namelist::Group (tested with v0.5.1 from CPAN),
#           perl (tested with v5.12.1)
#
#=USAGE=================================================================
#
# Uses Fortran::F90Namelist::Group to read files consisting of Fortran
# namelists and compares their content with 

use strict;
use warnings;

use Fortran::F90Namelist::Group;

if (@ARGV < 2) {
    print "usage $0 <path-to-code> <input-file[s] ...>\n";
    exit(-1);
}

my $cp = shift(@ARGV);

my %namelists = ();
while (<$cp/*.[fF]{03,95,90,77,[oO][rR],}>) {
    open(my $cf, '<', $_);
    my $n = 1;
    my $cnt_line = 0;
    my $nl;
    while (<$cf>) {
	next if (/^\s*[!#]/);

	if ($cnt_line) {
	    s/^\s*&//;
	    $nl .= $_;
	} else {
	    next unless (/^\s*namelist/i);
	    $nl = $_;
	}
	$cnt_line = 0;

	if ($nl =~ /&\s*$/) {
	    $cnt_line = 1;
	    $nl =~ s/&\s*$//;
	} else {
	    $nl =~ s/\s+//g;
	    $nl =~ s/^namelist\/([^\/]+)\///i;
	    $nl = "\L$nl";
	    my %e = map {$_ => 1} split(/,/,$nl);
	    $namelists{"\L$1"} = {%e};
	    # my @e = split(/,/,$nl);
	    # $namelists{"\L$1"} = [ @e ];
	}
    }
    close ($cf);
}

while (my $fn = shift(@ARGV)) {
    print "$fn:\n";
    my $g = Fortran::F90Namelist::Group->new()
	or die "Creation of namelist group object failed.\n";
    $g->parse(file => $fn);

    while (my $n = $g->pop()) {
	$| = 1;
	print ' '.$n->name().' ('.$n->nslots().' definitions)... ';
	die('namelist '.$n->name()." seems unknown to the code.\n")
	    unless exists($namelists{$n->name()});

	my $nl = $namelists{$n->name()};

#	print $nl;

	while (my ($name, $entry) = each(%{$n->hash()})) {
	    die('definition "'.$name.'=<'.join('> <',@{$entry->{'value'}})
		.">\" seems unknown to the code.\n")
		unless exists($nl->{$name});
	}

	print "ok\n";
    }
    print " => [file ok]\n";
}
print "No errors were found!\n";
