#! /usr/bin/perl

use strict;

my ($procs, $numad, $tad, $tex);
my (%proclist, %proclevels, %procad, %procex);
my ($level, $adarr, $exarr);

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
			print "New proc $procs at $proclevels{$procs}\n";

		}
		else {
			$proclevels{$procs}++;
			print "Exi proc $procs at $proclevels{$procs}\n";
		}

		$level = $proclevels{$procs};
		$adarr = $procad{$procs};
		$exarr = $procex{$procs};
		$adarr->{$level} = $tad / $numad;
		$exarr->{$level} = $tex / $numad;
	}
}

foreach $procs (sort { $a <=> $b } keys %proclist) {
	print "Procs for $procs\n";

	$adarr = $procad{$procs};
	foreach $level (sort { $a <=> $b } keys %$adarr) {
		print "$level\n";
	}

}
