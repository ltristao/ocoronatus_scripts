#!/usr/bin/perl -w

#loop
open IN, "gene_alignments.txt";
while($line = <IN>){
	chomp($line);
	push(@files, $line);
}
close IN;

##IQTREE
foreach my $var (@files){
$path_out = '/home/luca/rebuild/species_trees/best_trees/';
	system("iqtree2 -m MFP -T 27 -s $var");
	}
