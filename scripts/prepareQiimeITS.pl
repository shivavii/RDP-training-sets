#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
prepareQiimeITS.pl

PURPOSE:

INPUT:
--fasta    <string> : Sequence file
--taxonomy <string> : Taxonomy file
				
OUTPUT:
STDOUT - fasta file with taxonomy in headers.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $fasta, $taxonomy);
my $verbose = 0;

GetOptions(
    'fasta=s' 	  => \$fasta,
    'taxonomy=s' 	=> \$taxonomy,
    'verbose' 	  => \$verbose,
    'help' 		    => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $counter = 0;
my %hash;
open(TAX, "<".$taxonomy) or die "Can't open $taxonomy\n";
while(<TAX>){
  chomp;
  my @row = split(/\t/, $_);
  my $id = shift(@row);
  my $taxonomy = join(" ", @row);

  #print STDERR "$id\n$taxonomy\n";
  $hash{$id} = $taxonomy;

}
close(TAX);

my $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
  my $header = $curr->header();
  my $seq    = $curr->seq();
  $header =~ s/>//;
  if(exists $hash{$header}){
    print STDOUT ">".$header." ".$hash{$header}."\n".$seq."\n";
  } 
  $counter++;
}

