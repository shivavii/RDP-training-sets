#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
fixSilvaDb.pl

PURPOSE:

INPUT:
--infile <string> : Fasta sequence file
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
    'infile=s'  => \$infile,
    'verbose'   => \$verbose,
    'help'      => \$help
);
if ($help) { print $usage; exit; }
die "File is empty of does not exists $infile\n" if((!-e $infile) and (!-s $infile));
die "--infile missing\n" unless($infile);

## MAIN
my %hash;

my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	my $header = $curr->header();
	$header =~ s/>//;

	my @header = split(/\s/, $header);
	my $id = shift(@header);
	my @lineage = split(/;/, $header);
  my $n = 6;

  #print($_."\n") foreach(@lineage);

  @lineage = @lineage[-$n..-1];
	
	my $newLineage = "";
	
	my $i=0;
	foreach(@lineage){
    $_ =~ s/ //g;
		$newLineage .= "k__Bacteria;" if($i == 0);
		$newLineage .= "p__Cyanobacteria;" if($i == 1);
		$newLineage .= "c__Chloroplast;" if($i == 2);
		$newLineage .= $_.";" if($i == 3);
		$newLineage .= $_.";" if($i == 4);
		$newLineage .= $_.";" if($i == 5);
		$i++;
	}

  $newLineage = "Root;".$newLineage;

	my $seq = $curr->seq();
	$seq =~ s/\s//g;
	$seq =~ s/U/T/gi;
		
	print STDOUT ">".$id." ".$newLineage."\n".$seq."\n";
}

