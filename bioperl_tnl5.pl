#!/usr/bin/perl -w   
# bioperl_tnl5 by TNL5
#
# A BioPerl pipeline for sequence analysis.
#
# This script uses bioperl modules to run a protein BLAST [1] on a single sequence in fasta format. The top 10 results of the 
# BLAST search are outputted into file saved as "hits.txt".
# 
# The fasta sequences of the first top 10 hits are retreived from the GenPept database and saved as a file called "hits.fasta".
#
# A multiple sequence alignment is then carried out on 10 sequences and an output saved in clustalw format as "hits.align". A choice of either
# a clustalw [2] or muscle [3] alignment can be set the command line.
#
# Finally, the program clustalx is opened and the alignment is viewed
# 
# This script requires clustalx to be installed. Also requires either clustalw or muscle to be installed.
# 
# The script should be run as follows:
#	perl bioperl_tnl5 <sequence file name> <alignment type>
#
# Examples:
#								
#	perl bioperl_tnl5 P11802.fasta m        		- runs muscle alignment
# 	perl bioperl_tnl5 P11802.fasta c			- runs clustal alignment
#
#
# K&R style formatting used.
# 
# References:
#
# The script found at http://doc.bioperl.org/releases/bioperl-1.0.1/Bio/Tools/Run/RemoteBlast.html#CODE10 was modified for this assignment to
# prepare the "Multiple nested while loop to get BLAST and FASTA file" section
# 
# [1] Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
# [2] Larkin MA, Blackshields G, Brown NP, Chenna R, McGettigan PA, McWilliam H, Valentin F, Wallace IM, Wilm A, Lopez R, Thompson JD, Gibson TJ and Higgins DG 
# (2007) "ClustalW and ClustalX version 2"  Bioinformatics 23(21): 2947-2948. doi:10.1093/bioinformatics/btm404 
# [3] Edgar, Robert C. (2004), "MUSCLE: multiple sequence alignment with high accuracy and high throughput", Nucleic Acids Research 32(5), 1792-97.
#
# Thomas Lawson (TNL5) May 2013
#
use strict; 

#-----------------------------------------------
# External modules
#-----------------------------------------------

use Bio::SeqIO;
use Bio::Tools::Run::RemoteBlast;
use Bio::DB::GenPept;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::AlignIO;


#-----------------------------------------------
# Initilization of the code. Variables declared.
#-----------------------------------------------
my $file  = $ARGV[0];
my $alignChoice  = $ARGV[1];
my $blast;
my $prog;
my $db;
my $str;
my $input;
my $count = 1;
my $countinput = 0;
my $hit;
my $name;
my $hsp;
my $seq_object;
my $rid;
my $rc;	
my $muscle;
my $aln;
my $out;
my $result;
my $clustalw;
my @rids;
my @BLASTparams;
my @clustalparams;

#-----------------------------------------------
# Check of console inputs and file input
#-----------------------------------------------

# Check user inputed file name
if(!defined $ARGV[0]){
	die "No file defined \n";
}

# Check file exists
open (INPUT, "$file") or die "$file does not exist in your current directory\n";
close INPUT;

#Check user has inputed align info
if(!defined $alignChoice or $alignChoice !~ m/^[mMcC]+$/){ 
	die "Please select if Muscle alignment (m), or ClustalW alignment (c) required.
	e.g.
	perl bioperl_tnl5 P11802.fasta m        		- runs muscle alignment
	perl bioperl_tnl5 P11802.fasta c			- runs clustalw alignment\n";
}

# Create a input object ($str) from the SeqIO class using the input file identified from the command line (in fasta format)
$str = Bio::SeqIO->new(-file=>$file , -format => 'fasta' );

# Check to see how many sequences there are
while(my $seq = $str->next_seq) {
  $countinput++;
}

# if count greater than 2 print out a warning
if ($countinput++ >= 2){
	print  "WARNING: Input file has $countinput sequences. This script will only blast the first fasta sequence;\n";
}

#bioperl prints an exception to the console when the format is not in fasta so a check is not required here.


#-----------------------------------------------
# preparing parameters for the BLAST search
#-----------------------------------------------

# Parameters for BLAST 
$db   = 'swissprot';

#Assign parameters to HASH for BLAST
@BLASTparams = ( '-data' => $db );

# Create an object ($blast) from the RemoteBlast class with the defined paramters
$blast = Bio::Tools::Run::RemoteBlast->new(@BLASTparams);


#------------------------------------------------------------------------------------------------------------------------------------
# Preparing files from the BLAST search and sequences retrival from the database
#------------------------------------------------------------------------------------------------------------------------------------

# Create a new file called "hits.txt" this will contain the results of the BLAST
open (HITS, ">hits.txt") or die "can't create hits file'";

#Create an object ($db_object) from the database GenPept class. Used to retrieve sequences
my $db_object = Bio::DB::GenPept->new();

# $str reverted back to the original file after checking the number of sequences
$str = Bio::SeqIO->new(-file=>$file , -format => 'fasta' );

# create new output object ($out_object) from the SeqIO (in fasta format). This will create a new file called >hits.fasta
# containg the sequences of the top 10 alignments
my $out_object = Bio::SeqIO->new( '-file' => ">hits.fasta",'-format' => 'fasta');


#------------------------------------------------------------------------------------------------------------------------------------
# Multiple nested while loop to get BLAST and FASTA file
#
# Multiple nested while and for loops are carried out to BLAST the query sequence, save selected results as a text file and 
# get fasta sequences out of a database. 
# The first section is involved with connecting to the NCBI website and BLAST server and waiting for the result. This is attempted
# repeatedly as can sometimes take a while. If the internet connection is lost the script is stopped and an exception 
# message is sent to the console.
#
# if a result form the NCBI can be retrieved the BLAST result is obtained and custom report of the top 10 hits are saved to a file.
# The fasta sequences for these top 10 hits are taken from GenPept database and saved into a new file.
#
#------------------------------------------------------------------------------------------------------------------------------------

# Get the sequence from the query file
$input = $str->next_seq();
#Submit blast job to ncbi blast queue
$blast->submit_blast($input);
# Print out to console waiting message while retrieval ID is 
print STDERR "waiting...";

# While loop which gets the request ID of the blast search. allows the script to wait until BLAST report can be made as well
# NOTE: Only one request ID actually used here as only one BLAST being performed
while ( @rids = $blast->each_rid ) {
	# Get the one and only request ID
	$rid = $rids[0];
	# Retrieve the blast result info
    	$rc = $blast->retrieve_blast($rid);
	# if statment to see if BLAST output ready	
	if( !ref($rc) ) {
		# I think as removing the request ID stops the while loop this will stop 
		# while loop on this iteration or the next
		# and this is to be done when the blast result info is <0.
		# Note: I have tested with and without and it seems it is required
		if( $rc < 0 ) {
			$blast->remove_rid($rid);
		}
		
	# wait 5 seconds for result from NCBI then checks then repeats		
	print STDERR ".";
	sleep 5;
	} else {
		#Result is ready, so BLAST result is saved
		$result = $rc->next_result();
		#Removes the request ID so while loop can stop
		$blast->remove_rid($rid);
		#Print to hits file with the query name
		print HITS "Query Name: ", $result->query_name(), "\n";
		# While loop to get the top 10 hits					
		while ($count <= 10){
			# Count for loop			
			$count++;	
			# Get the result of the hit	
			$hit = $result->next_hit; 
			# Get the name of the hit
			$name = $hit->name;
			# Print the name to the hits file
			print HITS "$name \n";
			# Print sequence identifier, score, Expect-value and percent sequence identity	
			# NOTE: this is from the high scoring pair results 		
			$hsp = $hit->next_hsp;					
			print HITS "\tscore is ", $hsp->score, "\n";
			print HITS "\tsignificance is ", $hsp->significance, "\n";
			printf HITS ("\tPercent_id is %1.1f %% \n",$hsp->percent_identity),
			
			# download complete fasta sequence of the hit
			$seq_object = $db_object->get_Seq_by_id($name);
			# create a file of the fasta sequences
			$out_object->write_seq($seq_object);
			}
		
		}
}


#------------------------------------------------------------------------------------------------------------------------------------
# Multiple sequence alignment
#
# Clustal and muscle methods available. clustalw or muscle can be chosen by the user at the console
#
#------------------------------------------------------------------------------------------------------------------------------------

# Both clustalw and muscle use the same input file
my $inputfilename = 'hits.fasta';

if($alignChoice eq "m" or $alignChoice eq "M"){   

	#Build a muscle alignment factory
	$muscle = Bio::Tools::Run::Alignment::Muscle->new();

	# $aln is a SimpleAlign object.
	$aln = $muscle->align($inputfilename);
	
	# prepare a alignment file
	$out = Bio::AlignIO->new(-file   => ">hits.aln", -format => 'clustalw');
	
	# Create the file
	$out->write_aln($aln);
}else{
	#  Build a clustalw alignment factory
	@clustalparams = ('output' => "clustalw");
	$clustalw = Bio::Tools::Run::Alignment::Clustalw->new(@clustalparams);
	
	#  Pass the factory a list of sequences to be aligned.	
	$aln = $clustalw->align($inputfilename); # $aln is a SimpleAlign object.
	
	# prepare a alignment file
	$out = Bio::AlignIO->new(-file   => ">hits.aln", -format => 'clustalw');
	
	# Create the file
	$out->write_aln($aln);
}

# Open up clustalx with the alignment file
system "clustalx hits.aln";

close HITS; 
