#! /usr/bin/env perl

# This script will be used to perform orthologue identification across different red scores
# We will be using the bac120 genes to determine when these algorithms start to "fall off"
# in the identification of orthologues
# our goal is to identify a RED score value that uses the most (if not all) number of genomes while
# Being robust in identifying orthologous gene clusters


###TO DO:
#ADD LOGGER STATEMENTS!!!!


use strict;
use warnings;

BEGIN { our $start_time = time(); } #this will keep track of the amount of time the script has been running

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;
use List::Util qw(sum min max);
use FindBin qw($Bin);

#My Variables
my $help = 0;
my $man = 0;
my $RED_code_full_path = "$Bin/get_RED_groups.R"; #this is a path that will be filled in based on the GIT repository
my $path_to_orthostats = "$Bin/orthogroup_statistics.pl";
my $phylo_tree;
my $percent_id = .95;
my $out_dir;
my $fasta_dir;
my $fasta_ext = "faa";
my $num_genomes; #this will be calculated to check the fasta extension and directory parameters
my $bac120_tsv = "bac120_gene_info.tsv";
my $hmm_dir = "$Bin/../data/bac120_hmms/";
my $orthofinder_dir = "";
my $max_queue = 500; #how many jobs can be handled by the job server at once
my $config_file;
my $re_run = 0;
my $runtime_max = 518400; #about 6 days in seconds

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
				'config_file|c:s'	=>	\$config_file,
            'RED_code|r:s' => \$RED_code_full_path,
	    'phylo_tree|tree:s' => \$phylo_tree,
	    'perc_id:f' => \$percent_id,
	    'out_dir|out:s' => \$out_dir,
	    'fasta_dir|db:s' => \$fasta_dir,
	    'fasta_ext|ext:s' => \$fasta_ext,
	    'bac120_tsv|bt:s' => \$bac120_tsv,
	    'hmm_dir|hd:s' => \$hmm_dir,
		 'ortho_dir|od:s'	=>	\$orthofinder_dir,
		 'max_queue|mq:i'	=>	\$max_queue,
		 'num_genomes|ng:i'	=> \$num_genomes,
		 'max_runtime|run:i'	=>	\$runtime_max, 
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
if ( $re_run == 0 ) { #skips this step in the case of reruns
	create_directories_based_on_RED();
}

#Step 2 is to run ortho finder on all the groups
#need to handle a re-run situation
run_orthofinder($re_run);

#Step 3: Get completeness statistics
get_ortholog_stats();

## Subroutines ##

#will be finished when we know exactly what we need to work on
sub check_input {
	$logger->info("Checking for all necessary inputs\n");
	if ( $config_file ) {
		get_input();
		$re_run = 1;
	}
	
	if(! defined $RED_code_full_path){
		pod2usage(-message => "ERROR- Required parameter not found: --RED_code\n", -exitval => 2);
	}
	if(! -e $RED_code_full_path){ #can also check if R is installed here. Could add a -h option to rscript so we can test if it specifically is the script we want.
		pod2usage(-message => "ERROR- RED code file $RED_code_full_path deos not exist. Make sure the path is correct.", exitval => 2);
	}
	if(! defined $phylo_tree){
		pod2usage(-message => "ERROR- Required parameter not found: --phylo_tree\n", -exitval => 2);
	}
	if(! -e $phylo_tree){
		pod2usage(-message => "ERROR- Tree file $phylo_tree does not exist. Make sure the path is correct.", exitval => 2);
	}
	if(! defined $out_dir){
		pod2usage(-message => "ERROR- Required parameter not found: --out_dir\n", -exitval => 2);
	}
	if(! defined $fasta_dir){
		pod2usage(-message => "ERROR- Required parameter not found: --fasta_dir\n", -exitval => 2);
	}
	if(! -d $fasta_dir){
		pod2usage(-message => "ERROR- The directory $fasta_dir does not exist. Make sure the path is correct.", exitval => 2);
	}
	if(! defined $num_genomes){
		$logger->info("Number of genomes was not passed. We will calculate it using the fasta directory and extension\n");
		$num_genomes = (`ls -1 $fasta_dir/*$fasta_ext | wc -l`) + 0;
		if ( $num_genomes < 2 ) {
			pod2usage(-message => "ERROR- No genomes were found using fasta extension and fasta directory options\nMake sure they are passed correctly: --fasta_dir and/or --fasta_ext\n", -exitval => 2);
		}
		$logger->info("$num_genomes genome files were found in $fasta_dir");
	}
	else{
		my $test = (`ls -lh $fasta_dir/*$fasta_ext | wc -l`) + 0;
		if ( $num_genomes != $test ) {
			pod2usage(-message => "ERROR- Genome number does not match the count obtained from the fasta extension and fasta directory options\nMake sure they are passed correctly: --fasta_dir and/or --fasta_ext\n", -exitval => 2);
		}
	}
	if( ! -e $bac120_tsv ) {
		if(! defined $hmm_dir){
			pod2usage(-message => "ERROR- Bac120 tsv file not found. Must provide either this file or a directory of marker hmms.", -exitval => 2);
		}
		else{
			$logger->info("Bac120 tsv file not found. A new one will be created in the out directory\n");
			$bac120_tsv = get_bac120_genes();
		}
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
	my $R_installed = system("Rscript --help >/dev/null 2>&1");
	if($R_installed != 0){
		die "R command not found. Is R installed/loaded?\n";
	}
   #run the r-script to calculate RED score
   #pass tree (-t), percent_id (-p), and out_dir  (-o)
   `Rscript --no-save --no-restore $Bin/find_special_nodes.R -t $phylo_tree -o $out_dir/$id.dir`;
   `Rscript --no-save --no-restore $RED_code_full_path -t $phylo_tree -p $id -n $out_dir/$id.dir/nodes_to_root.tsv -o $out_dir/$id.dir`;
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
   if ( ! -d "$out_dir/tmp" ) {
      mkdir "$out_dir/tmp";
   }
   #I will create an array to store all of the directories needed
   my @perc_ids;
   #first create the test databases if need be
   if ( $percent_id == 0 ) {
      for ( my $i = .80; $i >= .49; $i = $i - 0.05) {
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

#remember that some of the higher RED vals will have a lot of genomic clusters but a small number of genomes inside them
#while the lower number RED values will have a smaller amount of clusters, but more genomes inside each one
sub run_orthofinder {
	my ($run_situation) = @_;
	
	#run orthofinder
	my $run_val = system($orthofinder_dir . "orthofinder -h >/dev/null 2>&1");
	#first check on orthofinders installation/load
	if ( $run_val != 0 ) {
		die "Orthofinder is not found. Make sure to pass in the path with -od or make sure it is installed/loaded\n";
	}
	$run_val = system("blastn -h >/dev/null 2>&1");
	if ( $run_val != 0 ) {
		die "Blast cannot be run. Make sure to put the code in your PATH and/or make sure it is installed/loaded\n";
	}
	my @commands;
	my @names;

	if ( $run_situation == 0 ) { #this checks if the script has been re-reun or not
		$logger->info("Attempting to run orthofinder within each cluster at each possible RED value\n");
		$logger->info("Building orthofinder commands...\n");
		
		
		#I will go through each directory in the out_dir and create an orthofinder command for it
		opendir( my $ORIG_OUT, $out_dir );
		while ( readdir $ORIG_OUT ) {
			if( $_ !~ /\.dir/ ){ next; }
			my $reddir = $_;
			opendir( my $RED_DIR, "$out_dir/$_" );
			while ( readdir $RED_DIR ) {
				if ( $_ !~ /_fasta/ ) { next; }
				push @commands, "orthofinder -S blast -f $out_dir/$reddir/$_";
				push @names, "$reddir.$_";
				#$commands{"$reddir.$_"} = "orthofinder -S blast -f $out_dir/$reddir/$_\n"; #THIS NEEDS TO BE THREADED WITH A SPECIFIC VALUE
			}
			closedir($RED_DIR);
		}
		closedir($ORIG_OUT);
	}
	
	
	
	#keep re-running orthofinder untill all the commands have finished
	submit_orthofinder_jobs(\@commands, \@names); #will have to add a feature to check for initial or new orthofinder runs
	my ($command_aref, $name_aref) = check_orthofinder_jobs();
	@commands = @{$command_aref};
	@names = @{$name_aref};
	while ( scalar @commands > 0 ) {
		submit_orthofinder_jobs(\@commands, \@names);
		($command_aref, $name_aref) = check_orthofinder_jobs();
	}
   return();
}

#Here should go the script to create a tsv file that identifies the bac120 genes in all of the genomes
sub get_bac120_genes {
	#the goal here is to have the code that creates a nice metadata file that contains the identity of all the bac120 genomes in the database provided
	#There is an optioon to pass this file so that it only needs to be run once.
	my $hmmer_installed = system("hmmsearch -h >/dev/null 2>&1"); #this redirects the stderror and the output to /dev/null (Nice piece of code Griffen)
	if($hmmer_installed != 0){
		die "hmmsearch command not found. Is HMMER installed/loaded?\n";
	}
	my $outpath = "bac120_gene_info.tsv";
	$logger->info("Attempting to create bac120 gene id file at $out_dir/$outpath\n");
	
	#For each genome, find the gene that is most similar to each marker.
	#Write to a file with three columns: "bac120_gene	genome_id	gene_id"
	my @genomes;
	opendir(my $FASTAS, $fasta_dir) or die "Cannot open fasta directory: $fasta_dir";
	while(readdir $FASTAS){
		if($_ =~ /^\./){ next; } #skip . and .. or any hidden files
		my ($genome_file, $genome_id) = ($_, $_);
		$genome_id =~ s/\.$fasta_ext//; #Get id without extension to write to output
		push @genomes, $genome_id;
		$logger->debug("About to run hmmsearch for each marker hmm on genome $genome_id");
		while(`squeue -h -n hmmsearch | wc -l` >= $max_queue){
			sleep 30;
		} 
		`sbatch $Bin/run_hmmsearch.sh $hmm_dir $fasta_dir/$genome_file $out_dir/tmp/hmmsearch_out/$genome_id.txt`;
	}
	close $FASTAS;
	while(`squeue -h -n hmmsearch | wc -l` > 0){
		sleep 120;
		$logger->debug("Waiting for hmmsearch jobs to finish");
	}
	open(my $outfile, ">", "$out_dir/$outpath");
	print $outfile "marker_gene\tgenome_id\tgene_id\n";
	foreach my $genome (@genomes){
		open(my $results, "<", "$out_dir/tmp/hmmsearch_out/$genome.txt");
		my ($hmm, $gene);
		foreach my $line (<$results>){
			chomp $line;
			if($line =~ /^Accession/){
				my @vals = split ' ', $line;
				$hmm = $vals[1];
			}
			if($line =~ /^>>/){
				my @vals = split ' ', $line;
				print $outfile "$hmm\t$genome\t$vals[1]\n";
			}
		}
	}
	close($outfile);
	return("$out_dir/$outpath");
}

#should maybe provide a way to run this without running the rest of the script
sub submit_orthofinder_jobs {
   my ($jobs_aref, $names_aref) = @_;
	
	$logger->info("Submitting orthofinder jobs to sbatch");
	
	my @jobs_to_submit = @{$jobs_aref};

	#all the necessary parameters for sbtach
	my $time = "-t 168:00:00"; #currently set for 1 week
	my $memory = "--mem=10g"; #set to 4gb until we know we need something more
	my $name = "ortho";
	my $threads = "-n 16"; #-n
	my $nodes = "-N 1"; #-N -> want to keep everything together
	my $job_output_dir = "-o $out_dir/tmp"; #this will be the directory where the sbatch output will be stored
	my $partition = "-p general";
	
	#create a temp directory to handle output
	if ( ! -d  $job_output_dir ) {
		mkdir $job_output_dir;
	}
	
	#build job commands for sbatch and then run them. Keep track of any jobs that are over the maximum allowed in the lsf queue
	my @jobs_remaining;
	my $job_count = 0;

	if($re_run == 1){
		$logger->info("Getting the Jobs that were not submitted yet\n");
		open my $JOBS, "<", "$out_dir/.restart.jobs";
		while( <$JOBS> ) {
			chomp $_;
			push @jobs_to_submit, $_;
		}
		close $JOBS;
		`rm $out_dir/.restart.jobs`; #get rid of the old jobs folder
	}

	#this will handle resubmitted scripts
	if ( -s "$out_dir/ACTIVE_JOBS" ) {
		my @wcout = split ' ', `wc -l $out_dir/ACTIVE_JOBS`;
		$job_count = $wcout[0];
	}
	my $i = 0;
	foreach my $command ( @jobs_to_submit ) {
		if($re_run == 0){
			my $job_out_name = "$job_output_dir/${$names_aref}[$i].out";
			$command = "sbatch $partition $time $memory -J $name $threads $nodes $job_out_name --wrap=\"$command\"";
		}
		if ( $job_count > $max_queue ) {
			push @jobs_remaining, $command;
		}
		else{
			system($command);
		}
		$job_count++;
		$i++;
	}
	
	#I will now try to get the code to wait until orthofinder is finished running
	stall_for_orthofinder(\@jobs_remaining, $name);
	return();
}

sub stall_for_orthofinder {
	my ($jobs_remaining_aref, $name) = @_;
	$logger->info("Stalling script for orthofinder to complete running\n");
	
	#dereference the href
	my @jobs_remain = @{$jobs_remaining_aref};
	
	#Create Active job file
	`echo start > $out_dir/ACTIVE_JOBS`;
	
	#How long we want the script to run before restarting itself
	
	#We need to check if the jobs are still running
	while ( -s "$out_dir/ACTIVE_JOBS" ) {
		#this will keep running unntil there are no more jobs
		sleep 120; #sleep for 2 minutes
		`squeue --name=$name -h > $out_dir/ACTIVE_JOBS`;
		
		#check to see if we can start running any more jobs
		if ( scalar(@jobs_remain > 0) ) {
		 my @wcout = split ' ', `wc -l $out_dir/ACTIVE_JOBS`;
		 my $lines = $wcout[0];
		 if ( $lines < $max_queue ) {
			 my $command = shift @jobs_remain;
			 system($command);
		 }
		}
		
		#check how long we have been running for
		if ( (time() - our $start_time) > $runtime_max ) {
			restart(\@jobs_remain);
		}
    }
	
	return();
}

#need to check if orthofinder jobs finished correctly or if they should be re-run starting from a particular position
sub check_orthofinder_jobs {
	#EASIER TO KEEP SBATCH OUTPUT AS IT'S EASIER TO PARSE. STILL ISSUE TO RESUME BEACUSE YOU CANNOT STATE DESIRED OUTPUT PATH, SO IT WILL CREATE NEW ORTHOFINDER FOLDER INSIDE WORKINGDIRECTORY...
	my @commands;
	my @names;
	opendir( my $ORIG_OUT, $out_dir);
	while( readdir $ORIG_OUT ){
		if( $_ !~ /\.dir/ ){ next; }
		my $reddir = $_;
		opendir( my $RED_DIR, "$out_dir/$reddir" );
		while( readdir $RED_DIR ){
			if( $_ !~ /_fasta/ ){ next; }
			my $groupdir = $_;
#I will be making sure to get the directory where working directory is
			opendir( my $GROUP_DIR, "$out_dir/$reddir/$groupdir" );
			my $dir_containing_ortho_working_dir = "";
			while ( readdir $GROUP_DIR ) {
if ( $_ !~ /\.faa/ ) {
    $dir_containing_ortho_working_dir = $_;
}
			}
			close($GROUP_DIR);
			my $name = "$reddir.$groupdir";
			if(! -e "$out_dir/tmp/$name.out"){
				push @commands, "orthofinder -S blast -f $out_dir/$reddir/$groupdir";
				push @names, $name;
				next;
			}

			my $stage;
			open(my $logfile, "<", "$out_dir/tmp/$name.out");
			foreach my $line (<$logfile>){
				chomp $line;
				if( $line =~ m/Done all-versus-all sequence search/i ){ #have to move blast output and various other files to node_##_fasta dir to use -b and get output in riht place
					$stage = "blast";
				}
				if($line =~ m/Done orthogroups/i ){
					$stage = "done";
					last;
				}
			}
#I think this is checking before things can start running. This will allow the error to stop
			if ( $stage eq "" ) {
			    next;
			}
			close($logfile);
			if( $stage eq "done" ){
				if( -e "$out_dir/$reddir/$groupdir/$dir_containing_ortho_working_dir/WorkingDirectory/SpeciesIDs.txt" ){
					system("rm $out_dir/$reddir/$groupdir/$dir_containing_ortho_working_dir/WorkingDirectory/Blast*");
                                	system("rm $out_dir/$reddir/$groupdir/$dir_containing_ortho_working_dir/WorkingDirectory/Species*.fa");
                                	system("rm $out_dir/$reddir/$groupdir/$dir_containing_ortho_working_dir/WorkingDirectory/*IDs.txt");
				}
				next;
			}
			if( $stage eq "blast" ){
				#move the files to the node fasta directory
				
				push @commands, "orthofinder -b $out_dir/$reddir/$groupdir/$dir_containing_ortho_working_dir/WorkingDirectory";
				push @names, $name;
			}
		}
		closedir($RED_DIR);
	}
	closedir($ORIG_OUT);
	return(\@commands, \@names);
}

#Once orthofinder is complete, we will need to get the statistics from orthofinder. The idea will be to use the genomes from each RED group and get 2 overall statistics
#1 Purity -> (total genes in OG clusters that are bac120 genes)/(total # of genes within those orthogroups)
#2 Fragmentation -> (120 * # of RED Groups)/(OGs that contain all of the bac120 genes in all of the RED groups)

#these two stats are from Isai's Nature Genetics story
sub get_ortholog_stats {
	opendir( my $ORIG_OUT, "$out_dir" );
	while(readdir $ORIG_OUT){
		if( $_ !~ /\.dir/ ){ next; }
		my $reddir = $_;
		system("perl $path_to_orthostats $bac120_tsv $out_dir/$reddir");
	}
}

sub restart {
	my ( $remaining_jobs ) = @_;
	$logger->info("Restarting");
    
    my $config_file = save_config();
    $logger->debug("Restart config file: $config_file");
    
    my $jobs_file = save_jobs($remaining_jobs);
    $logger->debug("Restart jobs file: $jobs_file");
    
    # resubmit the master job (ie this script)
    my $command = "sbatch -p general -o $out_dir/tmp/reset_log.txt -e $out_dir/tmp/reset_log.txt -J restart -t 168:00:00 --wrap=\"";
    $command .= "perl $Bin/cluster_group_analysis.pl ";
    $command .= "--config_file $config_file ";
    #$command .= "--jobs_file $jobs_file ";
    $command .= "--debug \"";
    $logger->debug("Restart command: $command");
    
    `$command`;
    
    exit 0;  # die but with success.  :)
}

#be sure to update this as everything changes
sub save_config {
	#create a file that contains the restarted config parameters
	my $file = "$out_dir/.restart.conf";
	open my $CON, ">", $file or
		$logger->fatal("Cannot open file: $file");
	
	#print out all the parameters in this file
	print $CON "RED_code=$RED_code_full_path\n";
	print $CON "phylo_tree=$phylo_tree\n";
	print $CON "perc_id=$percent_id\n";
	print $CON "out_dir=$out_dir\n";
	print $CON "fasta_dir=$fasta_dir\n";
	print $CON "fasta_ext=$fasta_ext\n";
	print $CON "num_genomes=$num_genomes\n"; #this will be calculated to check the fasta extension and directory parameters
	print $CON "bac120_tsv=$bac120_tsv\n";
	if(defined $hmm_dir){ print $CON "hmm_dir=$hmm_dir\n"; }
	print $CON "ortho_dir=$orthofinder_dir\n"; #how many jobs can be handled by the job server at once
	print $CON "max_queue=$max_queue\n";
	print $CON "max_runtime=$runtime_max\n";
	
	close $CON;
	return($file);
}

sub save_jobs {
	my ($remaining_jobs_aref) = @_;
	#Create a file that contains all of the jobs that still remain in the queue
	my $file = "$out_dir/.restart.jobs";
    open my $JOBS, ">", $file or
        $logger->fatal("Cannot open file: $file");
		  
	my $print_string = join("\n", @$remaining_jobs_aref);
	print $JOBS $print_string;
	close $JOBS;
	return($file);
}

sub get_input {
	open my $CONFIG, "<", $config_file;
	my %options;
	while ( <$CONFIG> ) {
		chomp $_;
		my $op = (split /=/, $_)[0];
		my $value = (split /=/, $_)[1];
		$options{$op} = $value;
	}
	close $CONFIG;
	#`rm $config_file`;
	
	#go through and re apply all of the parameters
	$RED_code_full_path = $options{"RED_code"}; #this is a path that will be filled in based on the GIT repository
	$phylo_tree = $options{"phylo_tree"};
	$percent_id = $options{"perc_id"};
	$out_dir = $options{"out_dir"};
	$fasta_dir = $options{"fasta_dir"};
	$fasta_ext = $options{"fasta_ext"};
	$num_genomes = $options{"num_genomes"}; #this will be calculated to check the fasta extension and directory parameters
	$bac120_tsv = $options{"bac120_tsv"};
	$hmm_dir = $options{"hmm_dir"};
	$orthofinder_dir = $options{'ortho_dir'};
	$max_queue = $options{"max_queue"}; #how many jobs can be handled by the job server at once
	$runtime_max = $options{"max_runtime"};
	
	return();
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
		[--bac120_tsv (--bt)]
		[--hmm_dir (--hd)]

		[--help]
		[--man]

	--RED_code = Path to R script that calculates RED scores.
	--phylo_tree = Path to tree of all genomes in Newick format.
	--out_dir = Path to output directory.
	--fasta_dir = Path to directory containing a fasta for each genome to be included.
	--fasta_ext = Extension for fasta files (DEFAULT: "faa").
	--num_genomes = Number of genomes involved in the analysis.
	--perc_id = The percent identity analysis should be done at. Set to 0 to test over the interval from 0.5 to 0.95 (DEFAULT: 0.95).
	--bac120_tsv = A file containing gene identifications for all bac120 genes in each genome.
	--hmm_dir = Directory containing hmms of marker genes. REQUIRED IF (and only if) --bac120_tsv NOT PROVIDED.
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
	
=head1 AUTHOR(S)

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

Griffen Kingkinner

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
