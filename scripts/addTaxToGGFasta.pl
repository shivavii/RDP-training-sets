#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
addTaxToGGFasta.pl

PURPOSE:

INPUT:
--infile_fasta <string> : Fasta sequence file
--infile_tax <string    : Taxonomy file
        
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile_tax, $infile_fasta);
my $verbose = 0;

GetOptions(
    'infile_tax=s'    => \$infile_tax,
    'infile_fasta=s'  => \$infile_fasta,
    'verbose'         => \$verbose,
    'help'            => \$help
);
if ($help) { print $usage; exit; }
die "File is empty of does not exists $infile_tax\n" if((!-e $infile_tax) and (!-s $infile_tax));
die "File is empty of does not exists $infile_fasta\n" if((!-e $infile_fasta) and (!-s $infile_fasta));
die "--infile_tax missing\n" unless($infile_tax);
die "--infile_fasta missing\n" unless($infile_fasta);

## MAIN
my %hash;

# First loop through tax, store taxonomy in a hash
open(IN, '<'.$infile_tax) or die "Can't open file $infile_tax\n";
while(<IN>){
  chomp;
  my @row = split(/\t/, $_);
  my @tax = split (/;/, $row[1]);

  my $ID;
  if($row[0] =~ m/(\S+)/){
    $ID = $1;
  }else{
    die "Someting went wrong with matching seq ID.\n";
  }  

  my $tax_string;
  my $kingdom;
  my $phylum;
  my $class;
  my $order;
  my $family;
  my $genus;
  my $specie;

  foreach my $field (@tax){
    $field =~ s/\s//;
    if($field =~ m/k__/){
      $kingdom = $field;
    }elsif($field =~ m/p__/){
      $phylum = $field;
    }elsif($field =~ m/c__/){
      $class = $field;
    }elsif($field =~ m/o__/){
      $order = $field;
    }elsif($field =~ m/f__/){
      $family = $field;
    }elsif($field =~ m/g__/){
      $genus = $field;
    }elsif($field =~ m/s__/){
      $specie = $field;
    }
  }
  $tax_string = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$specie;
  $hash{$ID} = $tax_string;  
}
close(IN);
print STDERR "Parsed taxonomy file\n";

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
  my $header = $curr->header();
  $header =~ s/>//;
  my $ID;  

  #print "header:\t".$header."\n";
  if($header =~ m/(\S+)/){
    $ID = $1;
  }else{
    die "Someting went wrong with matching seq ID.\n";
  }  

  if(exists $hash{$ID}){
    print STDOUT ">".$ID." ".$hash{$ID}."\n".$curr->seq()."\n";
    delete $hash{$ID};
  }
}

