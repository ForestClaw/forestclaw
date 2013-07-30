#! /usr/bin/perl

use strict;
use Getopt::Std;

my (%opt, $average, $startproc, $found, $filestr);
my ($procs, $numad, $numex, $numrg, $tad, $tex, $trg);
my (%proclist, %proclevels, %procad, %procex, %procrg);
my ($level, $adarr, $exarr, $rgarr, $step);

# Process command line options

$opt{a} = 0;
$opt{s} = 1;
getopts ("as:", \%opt);
$average = $opt{a};
$startproc = $opt{s};

if ($average) {
	$filestr = "avg";
}
else {
	$filestr = "max";
}

# Grab data from output files

while (<>) {

	# Try the latest output convention first
	$found = 0;
	if ($average) {
if (m:Procs (\d+) advance (\d+) ([0-9.e+-]+) exchange (\d+) ([0-9.e+-]+) regrid (\d+) ([0-9.e+-]+):) {
		$found = 1;
		$procs = $1;
		$numad = $2;
		$tad = $3;
		$numex = $4;
		$tex = $5;
		$numrg = $6;
		$trg = $7;
}
	}
	else {
if (m:Max/P (\d+) advance (\d+) ([0-9.e+-]+) exchange (\d+) ([0-9.e+-]+) regrid (\d+) ([0-9.e+-]+):) {
		$found = 1;
		$procs = $1;
		$numad = $2;
		$tad = $3;
		$numex = $4;
		$tex = $5;
		$numrg = $6;
		$trg = $7;
}
	}

	# If the latest convention didn't match, try the older one
	if (!$found) {
	if ($average) {
if (m:Procs (\d+) advance (\d+) ([0-9.e+-]+) exchange ([0-9.e+-]+) regrid ([0-9.e+-]+):) {
		$found = 1;
		$procs = $1;
		$numad = $2;
		$tad = $3;
		$tex = $4;
		$trg = $5;
}
	}
	else {
if (m:Max/P (\d+) advance (\d+) ([0-9.e+-]+) exchange ([0-9.e+-]+) regrid ([0-9.e+-]+):) {
		$found = 1;
		$procs = $1;
		$numad = $2;
		$tad = $3;
		$tex = $4;
		$trg = $5;
}
	}
	}

	# Suppose we have identified relevant output
	if ($found) {

		if (!$proclist{$procs}) {
			my (%adhash, %exhash, %rghash);
			$proclist{$procs} = 1;
			$proclevels{$procs} = 1;
			$procad{$procs} = \%adhash;
			$procex{$procs} = \%exhash;
			$procrg{$procs} = \%rghash;
			#print "New proc $procs at $proclevels{$procs}\n";

		}
		else {
			$proclevels{$procs}++;
			#print "Exi proc $procs at $proclevels{$procs}\n";
		}

		$level = $proclevels{$procs};
		$adarr = $procad{$procs};
		$exarr = $procex{$procs};
		$rgarr = $procrg{$procs};
		$adarr->{$level} = $tad / $numad;
		$exarr->{$level} = $tex / $numad;
		$rgarr->{$level} = $trg / $numad;
	}
}

# Create strong scaling analysis

open ADFILE, ">advance_strong.txt";
open EXFILE, ">exchange_strong.txt";
open RGFILE, ">regrid_strong.txt";

foreach $procs (sort { $a <=> $b } keys %proclist) {
	print ADFILE "$procs";
	print EXFILE "$procs";
	print RGFILE "$procs";

	$adarr = $procad{$procs};
	$exarr = $procex{$procs};
	$rgarr = $procrg{$procs};
	foreach $level (sort { $a <=> $b } keys %$adarr) {
		print ADFILE "\t$adarr->{$level}";
		print EXFILE "\t$exarr->{$level}";
		print RGFILE "\t$rgarr->{$level}";
	}

	print ADFILE "\n";
	print EXFILE "\n";
	print RGFILE "\n";
}

close ADFILE;
close EXFILE;
close RGFILE;

# Create weak scaling analysis

open ADFILE, ">advance_weak.txt";
open EXFILE, ">exchange_weak.txt";
open RGFILE, ">regrid_weak.txt";

foreach $step (0, 1, 2) {
	$procs = $startproc * 2 ** (2 * $step);
	print ADFILE "$procs";
	print EXFILE "$procs";
	print RGFILE "$procs";

	$adarr = $procad{$procs};
	$exarr = $procex{$procs};
	$rgarr = $procrg{$procs};
	for ($level = 3; $level <= 8; ++$level) {
		print ADFILE "\t$adarr->{$level - (2 - $step)}";
		print EXFILE "\t$exarr->{$level - (2 - $step)}";
		print RGFILE "\t$rgarr->{$level - (2 - $step)}";
	}

	print ADFILE "\n";
	print EXFILE "\n";
	print RGFILE "\n";
}

close ADFILE;
close EXFILE;
close RGFILE;
