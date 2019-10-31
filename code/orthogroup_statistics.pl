#! usr/bin/env perl

use warnings;
use strict;

#This will take in a file or hash containing information on which gene in each genome corresponds to which bac120 gene
#as well as the red directory, which should be populated with multiple directories such as node_##_fasta
#in each node_##_fasta directory there should be a folder called OrthoFinder
#This script will go calculate purity and fragmentation statistics over all the of the node directories

my $bac_file = shift;
my %bac_genes;
my $red_dir = shift;

my $first = 1;
open(my $bac, "<", $bac_file);

#for now will create hash every time
foreach my $line (<$bac>){
	if($first == 1){ $first =0; next; }
	chomp $line;
	my @vals = split /\t/, $line;
	${ $bac_genes{$vals[0]} }{$vals[1]} = $vals[2];
}
close($bac);

my ($min_groups, $total_groups, $pure, $total) = (0,0,0,0);

opendir(my $RED, $red_dir);
while(readdir $RED){ #Find each node_##_fasta directory, which should contain orthofinder output
	if( $_ !~ /_fasta/ ){ next; }
	my $results_dir = `ls -td $red_dir/$_/OrthoFinder/*/ | grep -v \$/ | head -1`; #this is done to get the most recently created OrthoFinder output!
											#just in case multiple outputs exist. Better than *.
	chomp $results_dir;
	my $ortho_dir = "$results_dir/Orthogroups/";
	if(! -d $ortho_dir){ print "Expected orthogroup directory not found: $ortho_dir\n"; next; } #This should pass as long as orthofinder finished
	open(my $tsv, "<", "$ortho_dir/Orthogroups.tsv");
	my $header_line = readline($tsv);
	close($tsv);
	chomp $header_line;
	$header_line =~ s/\r//; #orthofinder adds this but chomp doesn't get it. Messes with matching later.
	my @genomes = split "\t", $header_line;
	shift @genomes; #first column name is Orthogroup header, so skip it
	$min_groups += 120;
	open(my $ortho_file, "<", "$ortho_dir/Orthogroups.txt");
	my @orthogroups = <$ortho_file>;
	close($ortho_file);
	foreach my $key (keys %bac_genes){#go through each bac120 gene
		my %bac_in_og = ();
		my %total_in_og = ();
		foreach my $genome (@genomes){
			if( exists ${ $bac_genes{$key} }{"$genome"}){#if genome has gene, find OG
				foreach my $line (@orthogroups){
					if($line =~ ${ $bac_genes{$key} }{$genome}){
						my ($OG) = $line =~ /^(.*?):/;
						$bac_in_og{$OG} += 1;#hash of OGs that contain this bac120 gene. value = num of bac120 genes in group
						my $count = $line =~ tr/ //;# number of spaces on that line equals number of genes in OG
						$total_in_og{$OG} = $count;
						last;
					}
				}
			}
		}
		my $num = keys %bac_in_og;
		$total_groups += $num;#total number of groups to get all genes that correspond to this bac120 gene
		foreach my $orthogroup (keys %bac_in_og){
			$pure += $bac_in_og{$orthogroup};#number of bac120 genes over all OGs with at least one
			$total += $total_in_og{$orthogroup};#number of total genes over all OGs with at least one bac120 gene
		}
	}
}
my $frag = $min_groups/$total_groups;
my $purity = $pure / $total;
open(my $output, ">", "$red_dir/orthogroup_statistics.txt");
print $output ("Fragmentation: $frag\n");
print $output ("Purity: $purity\n");

close($output);
