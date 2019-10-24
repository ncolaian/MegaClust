#! /usr/bin/env perl

# This script will be used to perform orthologue identification across different red scores
# We will be using the phylaamphora genes to determine when these algorithms start to "fall off"
# in the identification of orthologues
# our goal is to identify a RED score value that uses the most (if not all) number of genomes while
# Being robust in identifying 


###TO DO:
#ADD LOGGER STATEMENTS!!!!


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;
use List::Util qw(sum min max);

#My Variables
my $help = 0;
my $man = 0;
my $RED_code_full_path; #this is a path that will be filled in based on the GIT repository
my $phylo_tree;
my $percent_id = .95;
my $out_dir;
my $fasta_dir;
my $fasta_ext = "faa";
my $num_genomes = 20;
my $bac120_tsv = "bac120_gene_info.tsv";
my $hmm_dir = "";

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'RED_code|r=s' => \$RED_code_full_path,
	    'phylo_tree|tree=s' => \$phylo_tree,
	    'perc_id=f' => \$percent_id,
	    'out_dir|out=s' => \$out_dir,
	    'fasta_dir|db=s' => \$fasta_dir,
	    'fasta_ext|ext:s' => \$fasta_ext,
	    'num_genomes|num=s' => \$num_genomes,
		 'bac120_tsv|bt:s'	=>	\$bac120_tsv,
            'help|h' => \$help
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manual and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

## Main ##
check_input();
#step that gets red groups and creates directory structure. Will also create statistics for the RED groups
create_directories_based_on_RED();

#Step 2 is to run ortho finder on all the groups
run_orthofinder();

#Step 3: Get completeness statistics


## Subroutines ##

#will be finished when we know exactly what we need to work on
sub check_input {
	$logger->info("Checking for all necessary inputs\n");
	if(! defined $RED_code_full_path){
		pod2usage(-message => "ERROR- Required parameter not found: --RED_code\n", -exitval => 2);
	}
	if(! defined $phylo_tree){
		pod2usage(-message => "ERROR- Required parameter not found: --phylo_tree\n", -exitval => 2);
	}
	if(! defined $out_dir){
		pod2usage(-message => "ERROR- Required parameter not found: --out_dir\n", -exitval => 2);
	}
	if(! defined $fasta_dir){
		pod2usage(-message => "ERROR- Required parameter not found: --fasta_dir\n", -exitval => 2);
	}
	if(! defined $num_genomes){
		pod2usage(-message => "ERROR- Required parameter not found: --num_genomes\n", -exitval => 2);
	}
	if( ! -e $bac120_tsv ) {
		logger->info("Bac120 tsv file not found. A new one will be created in the out directory\n");
		get_bac120_genes();
	}
   $logger->info("All needed inputs were passed\n");
	return();
}

#This script will be divided into multiple parts
#I will describe what each will do and provide a methodology to perform this analysis
#We will first need to create a RED score matrix.
#I currently have some code for this in R and will ultimately use that script to create RED matrix
####IMPORTANT####
#this has been changed, based on the need to perform statistics on the ORTHOGROUP calling at each froup
sub get_defined_RED_groups {
   my ($id) = @_;
   #run the r-script to calculate RED score
   #pass tree (-t), percent_id (-p), and out_dir  (-o)
   `Rscript --no-save --no-restore $RED_code_full_path -t $phylo_tree -p $id -o $out_dir/$id.dir`;
   #get stastics along with creating the directory structure needed
   calculate_RED_group_statistics($id);
   return();
}

# create directories for each run based on RED
# Will then need to separete the data within each of these directories and create a metadata file that keeps track
# of where each genome is placed
###IMPORTANT### -> We will have to make an arbirtrary decision of what to do with RED clusters containing small number of genomes
# The first analysis should like at cluster analysis (size, average, variation, etc) at each RED value.
sub create_directories_based_on_RED {
   if ( ! -d $out_dir ) {
      mkdir $out_dir;
   }
   #I will create an array to store all of the directories needed
   my @perc_ids;
   #first create the test databases if need be
   if ( $percent_id == 0 ) {
      for ( my $i = .95; $i >= .50; $i = $i - 0.05) {
         push @perc_ids, $i;
      }
   }
   else {
      @perc_ids = ($percent_id);
   }
   
   #I will now create a directory for each of the elements in the array
   foreach my $ids ( @perc_ids ) {
      mkdir "$out_dir/$ids.dir";
      get_defined_RED_groups($ids);
   }
   
   return();
}

#this will also create the individual directories along with grabbing the fasta files
###DONE
sub calculate_RED_group_statistics {
   #####statistics i'd like to keep track of#####
   #Median group size
   #max group size
   #min group size
   #Singleton count
   #genomes used
   
   my ($id) = @_;
   
   open my $IN, "<", "$out_dir/$id.dir/red_groups.tsv"; 
   open my $OUT, ">", "$out_dir/$id.dir/group_stats.txt"; 
   
   my @info;
   while ( <$IN> ) {
      chomp $_;
      my @split_line = split(/\t/, $_); #0 -> node id and #1 -> genomes
      my @genomes = split(/,/, $split_line[1]);
      push @info, scalar(@genomes);
      mkdir "$out_dir/$id.dir/node_$split_line[0]_fasta/";
      #cp the fasta files to the correct directory
      foreach my $g ( @genomes ) {
         `cp $fasta_dir/$g.$fasta_ext $out_dir/$id.dir/node_$split_line[0]_fasta/`;
      }
   }
   close $IN;
   #print out the stats
   my $time = localtime();
   print $OUT "These stats are prepared at $id\% on $time\n";
   print $OUT "Total Genomes Used:\t" . sum(@info) . "\n";
   print $OUT "Median Group Size:\t" . median(@info) . "\n";
   print $OUT "Max Group Size:\t" . max(@info) . "\n";
   print $OUT "Min Group Size:\t" . min(@info) . "\n";
   print $OUT "Number of Singletons:\t" . ($num_genomes - sum(@info)) . "\n";
   close $OUT;
   return();
}

#We will then have to run Orthofinder on each cluster within each red value. This will obviously have to be parallelized.
#We may try to create a run sheet so that the submissions can be shared across users
sub run_orthofinder {
	#run orthofinder
	
   
}

#Here should go the script to create a tsv file that identifies the bac120 genes in all of the genomes
sub get_bac120_genes {
	#the goal here is to have the code that creates a nice metadata file that contains the identity of all the bac120 genomes in the database provided
	#There is an optioon to pass this file so that it only needs to be run once.
	
}

#should maybe provide a way to run this without running the rest of the script
sub submit_orthofinder_jobs {
   
}

#Once orthofinder is complete, we will need to get the statistics from orthofinder. The idea will be to use the genomes from each RED group and get 2 overall statistics
#1 Purity -> (total genes in OG clusters that are bac120 genes)/(total # of genes within those orthogroups)
#2 Fragmentation -> (120 * # of RED Groups)/(OGs that contain all of the bac120 genes in all of the RED groups)

#these two stats are from Isai's Nature Genetics story
sub get_ortholog_stats {
	
}

######we should also look into using other metrics to score RED scores#####

#We will then need to take this output and create an R-code to create the graphs we would like with this information
sub collapse_completeness_info {
   
}

sub create_graphs {
   
}

#this will get the median within an array without having to load any extra modules
sub median {
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2)
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

__END__
=head1 TITLE



=head1 VERSION



=head1 INCLUDED MODULES

Getopt::Long;
Pod::Usage;
Carp;
Readonly;
Path::Class;
Data::Dumper;
Log::Log4perl qw(:easy);
Log::Log4perl::CommandLine qw(:all);

=head1 INHERIT

=head1 SYNOPSIS

	cluster_group_analysis.pl
		--RED_code (--r)
		--phylo_tree (--tree)
		--out_dir (--out)
		--fasta_dir (--db)
		[--fasta_ext (--ext)]
		--num_genomes (--num)
		[--perc_id]

		[--help]
		[--man]

	--RED_code = Path to R script that calculates RED scores.
	--phylo_tree = Path to tree of all genomes in Newick format.
	--out_dir = Path to output directory.
	--fasta_dir = Path to directory containing a fasta for each genome to be included.
	--fasta_ext = Extension for fasta files (DEFAULT: "faa").
	--num_genomes = Number of genomes involved in the analysis.
	--perc_id = The percent identity analysis should be done at. Set to 0 to test over the interval from 0.5 to 0.95 (DEFAULT: 0.95).
	--help = Prints USAGE.
	--man = Prints man page.


=head1 PARAMETERS

=head1 CONFIGURATION AND ENVIRONMENT

    

=head1 DEPENDENCIES


    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 AUTHOR

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2019, Nicholas Colaianni
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut
