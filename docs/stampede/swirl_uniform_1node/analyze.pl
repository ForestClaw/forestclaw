#! /usr/bin/perl

use strict;

my ($procs, $numad, $tad, $tex);
my (%proclist, %proclevels, %procad, %procex);
my ($level, $adarr, $exarr, $step);

# Grab data from output files

while (<>) {

	if (m/Procs (\d+) advance (\d+) ([0-9.]+) exchange ([0-9.]+)/) {

		$procs = $1;
		$numad = $2;
		$tad = $3;
		$tex = $4;

		if (!$proclist{$procs}) {
			my (%adhash, %exhash);
			$proclist{$procs} = 1;
			$proclevels{$procs} = 1;
			$procad{$procs} = \%adhash;
			$procex{$procs} = \%exhash;
			#print "New proc $procs at $proclevels{$procs}\n";

		}
		else {
			$proclevels{$procs}++;
			#print "Exi proc $procs at $proclevels{$procs}\n";
		}

		$level = $proclevels{$procs};
		$adarr = $procad{$procs};
		$exarr = $procex{$procs};
		$adarr->{$level} = $tad / $numad;
		$exarr->{$level} = $tex / $numad;
	}
}

# Create strong scaling analysis

open ADFILE, ">advance_strong.txt";
open EXFILE, ">exchange_strong.txt";

foreach $procs (sort { $a <=> $b } keys %proclist) {
	print ADFILE "$procs";
	print EXFILE "$procs";

	$adarr = $procad{$procs};
	$exarr = $procex{$procs};
	foreach $level (sort { $a <=> $b } keys %$adarr) {
		print ADFILE "\t$adarr->{$level}";
		print EXFILE "\t$exarr->{$level}";
	}

	print ADFILE "\n";
	print EXFILE "\n";
}

close ADFILE;
close EXFILE;

# Create weak scaling analysis

open ADFILE, ">advance_weak.txt";
open EXFILE, ">exchange_weak.txt";

foreach $step (0, 1, 2) {
	$procs = 2 ** (2 * $step);
	print ADFILE "$procs";
	print EXFILE "$procs";

	$adarr = $procad{$procs};
	$exarr = $procex{$procs};
	for ($level = 3; $level <= 8; ++$level) {
		print ADFILE "\t$adarr->{$level - (2 - $step)}";
		print EXFILE "\t$exarr->{$level - (2 - $step)}";
	}

	print ADFILE "\n";
	print EXFILE "\n";
}

close ADFILE;
close EXFILE;
