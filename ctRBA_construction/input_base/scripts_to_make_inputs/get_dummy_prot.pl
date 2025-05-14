#Written by: Wheaton Schroeder
#Latest version: 08/25/2022
#Written to read the Amino acid-based genome annotation of C. therm and calculate average protien length

use strict;
use Statistics::Basic qw(:all);

#array to store protien length data
my @prot_lens = ();

#read in the annotated genome
open(ANNGEN, "<rt.gbff") or die "could not read annoted genome file, reason: $!\n";
chomp(my @data = <ANNGEN>);

#write a result file
open(RESULTS, ">prod_dummy_res.txt") or die "could not read annoted genome file, reason: $!\n";

#combine the lines to get a data string
my $data_string = join(" ", @data);

#list of amino acids
my @aas = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');

#save counts of each amino acid
my %aa_counts = { };

#save total amino acid count
my $total_aas = 0;

#initialize the counts
for(my $a = 0; $a <= $#aas; $a++) {

  $aa_counts{$aas[$a]} = 0;

}

#search the data string for translations, get each match
#captured text is the amino acid sequence
while ($data_string =~ m/ \/translation=\"(.+?)\"/g) {

  my $aa_seq = $1;

  #clean up the amino acid sequence by removing all spaces
  $aa_seq =~ s/\s+//g;

  #get the length of the aa sequence
  my $aa_len = length($aa_seq);

  #add to amino acid length to the array
  push @prot_lens, $aa_len;

  printf("count of ");

  #create an amino acid count that we will then use to calculate the average composition
  for(my $b = 0; $b <= $#aas; $b++) {

    #initialize a temporary count
    my $temp_count = 0;

    while($aa_seq =~ /$aas[$b]/gi) {

      $temp_count++;
      $total_aas++;

    }

    printf("%s: %d\t", $aas[$b], $temp_count);

    #update the count
    $aa_counts{$aas[$b]} += $temp_count;

  } 

  printf("\n");

}

my $median = median(@prot_lens);

printf("\n\nmedian protein length: %d", $median);
printf(RESULTS "median protein length: %d", $median);

my $mean = mean(@prot_lens);

printf("\nmean protein length: %d\n\n", $mean);
printf(RESULTS "\nmean protein length: %d\n\n", $mean);

#get fractions of different AAs in total proteome
for(my $c = 0; $c <= $#aas; $c++) {

  my $temp_frac = $aa_counts{$aas[$c]} / $total_aas;

  printf("%s\t%.8f\n", $aas[$c], $temp_frac);
  printf(RESULTS "%s\t%.8f\n", $aas[$c], $temp_frac);

}