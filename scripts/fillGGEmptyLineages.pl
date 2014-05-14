#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
count_seq.pl

PURPOSE:

INPUT:
--infile <string>   : Sequence file
				
OUTPUT:
--outfile <string>  : Sequence file
--outfilef <string> : Sequence file that failed

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay@lbl.gov

ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $outfilef);
my $verbose = 0;

GetOptions(
  'infile=s'    => \$infile,
  'outfile=s'   => \$outfile,
  'outfilef=s'  => \$outfilef,
  'verbose'     => \$verbose,
  'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT, ">".$outfile) or die "Can't open file ".$outfile."\n";
open(OUTF, ">".$outfilef) or die "Can't open file ".$outfilef."\n";

my $counter = 0;
my $envir_samples = 0;
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	my $header = $curr->header;
	#$header =~ s/>//;
	my @row = split(/;/, $header);

	my $id = $1 if($header =~ m/(>.*)k__/);
	my $k; 
	my $p;
	my $c;
	my $o;
	my $f;
	my $g;
  my $s;

	foreach(@row){
		if($_ =~ m/(k__.*)/){
			$k = $1; 
		}elsif($_ =~ m/(p__.*)/){
			$p = $1;
		}elsif($_ =~ m/(c__.*)/){
			$c = $1; 
			#$c =~ s/ \(Class\)/_Class/;
		}elsif($_ =~ m/(o__.*)/){
			$o = $1; 
		}elsif($_ =~ m/(f__.*)/){
			$f = $1; 
		}elsif($_ =~ m/(g__.*)/){
			$g = $1; 
		}elsif($_ =~ m/(s__.*)/){
			$s = $1; 
		}
	}

	my $new_header = "";
	
	# Special care for genus entries labeled g__environmental samples
	if(defined($g)){
		if($g eq "g__environmental samples"){
			$g = "g__environmental samples".$envir_samples;
			$envir_samples++;
		}
	}
	if(defined($f)){
		if($f eq "f__environmental samples"){
			$f = "f__environmental samples".$envir_samples;
			$envir_samples++;
		}
	}
	if(defined($o)){
		if($o eq "f__environmental samples"){
			$o = "f__environmental samples".$envir_samples;
			$envir_samples++;
		}
	}
	if(defined($c)){
		if($c eq "c__environmental samples"){
			$c = "c__environmental samples".$envir_samples;
			$envir_samples++;
		}
	}

	# TODO Add support for species.

	if( (defined($k) && defined($p) && defined($c) && defined($o) && defined($f) && defined($g) && defined($s) )
		&& ($k ne "k__" && $p ne "p__" && $c ne "c__" && $o ne "o__" && $f ne "f__" && $g ne "g__" && $s ne "s__") ){
		$new_header .= $id." ".$k.";".$p.";".$c.";".$o.";".$f.";".$g.";".$s.";";
		print OUT $new_header."\n".$curr->seq."\n"; 
	
	}elsif( (defined($k) && defined($p) && defined($c) && defined($o) && defined($f) && defined($g))
		&& ($k ne "k__" && $p ne "p__" && $c ne "c__" && $o ne "o__" && $f ne "f__" && $g ne "g__") 
		&& (lc($k) ne "k__environmental samples" && lc($p) ne "p__environmental samples" && lc($c) ne "c__environmental samples" && lc($o) ne "o__environmental samples" && lc($f) ne "f__environmental samples" && lc($g) ne "g__environmental samples") ) {
		my $newg = $g;
    $newg =~ s/g__//;
    $new_header .= $id." ".$k.";".$p.";".$c.";".$o.";".$f.";".$g.";s__".$newg."SP;";
		print OUT $new_header."\n".$curr->seq."\n";

	}elsif( (defined($k) && defined($p) && defined($c) && defined($o) && defined($f))
		&& ($k ne "k__" && $p ne "p__" && $c ne "c__" && $o ne "o__" && $f ne "f__")
		&& (lc($k) ne "k__environmental samples" && lc($p) ne "p__environmental samples" && lc($c) ne "c__environmental samples" && lc($o) ne "o__environmental samples" && lc($f) ne "f__environmental samples")) {
		my $newf = $f;
		$newf =~ s/f__//;
		$new_header .= $id." ".$k.";".$p.";".$c.";".$o.";".$f.";g__".$newf."FA;s__".$newf."FA;";
		print OUT $new_header."\n".$curr->seq."\n";
	
	}elsif( (defined($k) && defined($p) && defined($c) && defined($o) )
		&& ($k ne "k__" && $p ne "p__" && $c ne "c__" && $o ne "o__")
		&& (lc($k) ne "k__environmental samples" && lc($p) ne "p__environmental samples" && lc($c) ne "c__environmental samples" && lc($o) ne "o__environmental samples")) {
		my $newo = $o;
		$newo =~ s/o__//;
		$new_header .= $id." ".$k.";".$p.";".$c.";".$o.";f__".$newo."OR;g__".$newo."OR;s__".$newo."OR;";
		print OUT $new_header."\n".$curr->seq."\n";
	
	}elsif( (defined($k) && defined($p) && defined($c))
		&& ($k ne "k__" && $p ne "p__" && $c ne "c__")
		&& (lc($k) ne "k__environmental samples" && lc($p) ne "p__environmental samples" && lc($c) ne "c__environmental samples")) {
		my $newc = $c;
		$newc =~ s/c__//;
		$new_header .= $id." ".$k.";".$p.";".$c.";o__".$newc."CL;o__".$newc."CL;f__".$newc."CL;g__".$newc."CL;s__".$newc."CL;";
		print OUT $new_header."\n".$curr->seq."\n";

	}elsif( (defined($k) && defined($p))
		&& ($k ne "k__" && $p ne "p__")
		&& (lc($k) ne "k__environmental samples" && lc($p) ne "p__environmental samples")) {
		my $newp = $p;
		$newp =~ s/p__//;
		$new_header .= $id." ".$k.";".$p.";c__".$newp."PH;o__".$newp."PH;o__".$newp."PH;f__".$newp."PH;g__".$newp."PH;s__".$newp."PH;";
		print OUT $new_header."\n".$curr->seq."\n";

	}elsif( (defined($k))
		&& ($k ne "k__")
		&& (lc($k) ne "k__environmental samples")) {
		my $newk = $k;
		$newk =~ s/k__//;
		$new_header .= $id." ".$k.";p__".$newk."KI;c__".$newk."KI;o__".$newk."KI;o__".$newk."KI;f__".$newk."KI;g__".$newk."KI;s__".$newk."KI;";
		print OUT $new_header."\n".$curr->seq."\n";

	}else{
		print OUTF $curr->header."\n".$curr->seq."\n";
	}
}
exit;
