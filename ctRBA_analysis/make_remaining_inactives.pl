#!usr/bin/perl -w
#Written by: Wheaton Schroeder
#Written to take the list of minimal essential inactive calculated by get_essen_v_inactive.gms and the list of inactive reactions
#and spit out the list of 

use strict; 

#read the inactive list
open(INACTIVES, "<model_base/inactive_rxns.txt") or die "could not open/read inactive_rnxs.txt, reason: $!\n";
chomp(my @inactives = <INACTIVES>);

#remove "/" formatting
shift(@inactives);
pop(@inactives);

#clear out the quotation marks, newline characters
for(my $b = 0; $b <= $#inactives; $b++) {

	$inactives[$b] =~ s/\'//g;
	$inactives[$b] =~ s/\'//g;
	$inactives[$b] =~ s/\r//g;
	$inactives[$b] =~ s/\n//g;

}

#get the list of essential inactive reactions
open(ESSENTIALS, "<minimal_essential_inactives_av.txt") or die "could not open/read minimal_essential_inactives_av.txt, reason: $!\n";
chomp(my @essentials = <ESSENTIALS>);

#remove header text
shift(@essentials);

#get just the id of the essential reactions 
for(my $a = 0; $a <= $#essentials; $a++) {

	#get just the reaction id		
	(my $rxn_id, my $lb, my $flux, my $ub) = split(/\t/,$essentials[$a]);

	#store just the reaction ID
	$essentials[$a] = $rxn_id;
	
}

#create the output file
open(REMAINING, ">model_base/remaining_inactive_rxns_av.txt") or die "could not create/write remaining_inactive_rxns_av.txt, reason: $!\n";

#print the starting "/"

printf(REMAINING "\/\n");

#now we have a list of all inactive reactions, add a list of essential inactie reactions. write those in the former but not in the latter to the output
for(my $c = 0; $c <= $#inactives; $c++) {

	if (!($inactives[$c] ~~ @essentials)) {

		printf(REMAINING "\'%s\'\n", $inactives[$c]);

	}

}

#print the trailing "/"'
printf(REMAINING "\/");







