#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use List::Util qw(sum);

my $usage=<<'ENDHERE';
NAME:
generateRDPHierTree.pl

PURPOSE:
Generates a hierarchy taxonomy and fasta corresponding file to use as RDP training model.

INPUT:
--fasta <fasta_infile>           : Fasta header has to formatted like this: >ID<white_space>Taxonomy_level_1;Taxonomy_level_2;Taxonomy_level_n
--tax_level <int>                : Can be one of the following values: 1-kingdom, 2-phylum, 3-class, 4-order, 5-family, 6-genus, 7-specie
--greengenes_correction          : boolean 0: False and 1: True (Default value = True). If using Greengenes version of May 9th 2011, should be true.
--fungi_ITS_correction           : boolean 0: False and 1: True (Default value = True). For ITS fungi.


OUTPUT:
--outfile_model <txt_outfile>    : txt file representing the taxonomy.
--outfile_fasta <fasta_outfile>  : fasta file representing the taxonomy.
--outfile_lineages <tab_outfile> : file in which are represented all unique lineages.

NOTES:
This script was written to generate an RDP hierarchical model compatible with the rdp_classifier v2.3. It is intended to receive a Greengenes 
formatted fasta  file in input. Here is an example of a greengenes fasta header:
  
>46 M32222.1 Methanothermus fervidus k__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; Unclassified; otu_120 
ACGTACGTACGT...
- accepted taxonomic values in headers are k__, p__ c__, o__, f__, g__ and s__
- no tabular (i.e. \t) spaces is used to separate fields, please use a single whitespace character.
- In the current version of Greengenes DB (May 9th 2011), some redundant classification can be encountered, Ex, two Clostridium genus having different lineages.
     p__Tenericutes  c__Erysipelotrichi  o__Erysipelotrichales   f__Clostridiaceae   g__Clostridium
     p__Firmicutes   c__Clostridia   o__Clostridiales    f__Clostridiaceae   g__Clostridium

BUGS/LIMITATIONS:
Please take extra care to format the header of fasta headers to respect the Greengenes header style.

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $fasta, $tax_level, $outfile_model, $outfile_fasta, $outfile_lineages, $log, $greengenes_correction, $fungi_ITS_correction);
my $verbose = 0;

## SCRIPTS
GetOptions(
  'fasta=s'                   => \$fasta,
  'outfile_model=s'           => \$outfile_model,
  'outfile_fasta=s'           => \$outfile_fasta,
  'log=s'                     => \$log,
  'tax_level=i'               => \$tax_level,
  'greengenes_correction=i'   => \$greengenes_correction,
  'fungi_ITS_correction=i'    => \$fungi_ITS_correction,
  'outfile_lineages=s'        => \$outfile_lineages,
  'verbose'                   => \$verbose,
  'help'                      => \$help
);
if ($help) { print $usage; exit; }

#VALIDATE
die("--fasta fasta file required\n")                         unless $fasta;
die("--outfile_model outfile_model required\n")              unless $outfile_model;
die("--outfile_fasta outfile_fasta required\n")              unless $outfile_fasta;
die("--tax_level tax_level required\n")                      unless $tax_level;
die("--outfile_lineages outfile_lineages file required\n")   unless $outfile_lineages;
die("--tax_level tax level must be a value between 1 and 7") if($tax_level > 7 || $tax_level == 0);
$greengenes_correction = 0 unless $greengenes_correction;
$fungi_ITS_correction = 0 unless $fungi_ITS_correction;
## MAIN

print STDERR "Generating taxonomy hierarchical file to use in RDP...\n";

my $counter = 0;

my $fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
open(OUT_FASTA, ">".$outfile_fasta) or die "Can't find file ".$outfile_fasta."\n";
open(OUT_MODEL, ">".$outfile_model) or die "Can't open file ".$outfile_model."\n";
open(LOG, ">".$log) or die "Can't open file ".$log."\n";

my @taxonomy;
my %hash = ();
my @array = ();
my %hash_id;
my %seen_genus = ();
my %seen_family = ();
my @complete_array = ();

my $j=0;
my $incertae_count = 0;

my $k__;
my $p__;
my $c__;
my $o__;
my $f__;
my $g__;
my $s__;

while( my $seq = $fasta_db->next_seq() ) {
  my @header = split(";", $seq->header());
  my @header_part_1 =  split(" ", $header[0]);
  my $id = $1 if($header_part_1[0] =~ m/>(\S+)/);  
  #print STDERR $id."\n" if($verbose);
    
  @taxonomy = ();
  
  #Determine the taxonomy dept
  $k__=0;
  $p__=0;
  $c__=0;
  $o__=0;
  $f__=0;
  $g__=0;
  $s__=0;

  $taxonomy[0] = "";
  $taxonomy[1] = "";
  $taxonomy[2] = "";
  $taxonomy[3] = "";
  $taxonomy[4] = "";
  $taxonomy[5] = "";
  $taxonomy[6] = "";

  foreach(@header){
    if($_ =~ m/(k__.*)/) {
      $taxonomy[0] = $1;
      $k__++;
    }elsif($_ =~ m/(p__.*)/){
      $taxonomy[1] = $1;
      $p__++;
    }elsif($_ =~ m/(c__.*)/){
      $taxonomy[2] = $1;
      $c__++;
    }elsif($_ =~ m/(o__.*)/){
      $taxonomy[3] = $1;
      $o__++;
    }elsif($_ =~ m/(f__.*)/){
      $taxonomy[4] = $1;
      $f__++;
    }elsif($_ =~ m/(g__.*)/){
      $taxonomy[5] = $1;
      $g__++;
    }elsif($_ =~ m/(s__.*)/){
      $taxonomy[6] = $1;
      $s__++;
    }
  }
  
  if($k__>1){
    $k__=1;
    print LOG "Duplicate Kingdom: ".$id."\n";
  }    
  if($p__>1){
    $p__=1;
    print LOG "Duplicate phylum: ".$id."\n";
  }    
  if($c__>1){
    $c__=1;
    print LOG "Duplicate class: ".$id."\n";
  }    
  if($o__>1){
    $o__=1;
    print LOG "Duplicate order: ".$id."\n";
  }    
  if($f__>1){
    $f__=1;
    print LOG "Duplicate family: ".$id."\n";
  }    
  if($g__>1){
    $g__=1;
    print LOG "Duplicate genus: ".$id."\n";
  }    
  if($s__>1){
    $s__=1;
    print LOG "Duplicate species: ".$id."\n";
  }
  
  #Manually correct some inconsistencies in the Greengenes DB in date of 9 May 2011.
  if($greengenes_correction == 1){
    if($taxonomy[4] eq "f__Rhodobacteraceae"){
      $taxonomy[3] = "o__Rhodobacterales";
    }
    if($taxonomy[4] eq "f__Chromatiaceae"){
      $taxonomy[3] = "o__Chromatiales";
    }
    if($taxonomy[4] eq "f__Alteromonadaceae"){
      $taxonomy[3] = "o__Alteromonadales";
    }
    if($taxonomy[4] eq "f__Piscirickettsiaceae"){
      $taxonomy[3] = "o__Thiotrichales";
    }
    if($taxonomy[5] eq "g__Streptomyces"){
      $taxonomy[4] = "f__Streptomycetaceae";
    }
    if($taxonomy[5] eq "g__Ruminococcus"){
      $taxonomy[4] = "f__Ruminococcaceae";
    }
    if($taxonomy[4] eq "f__Thiotrichaceae"){
      $taxonomy[3] = "o__Thiotrichales";
    }
    if($taxonomy[5] eq "g__Eubacterium"){
      $taxonomy[4] = "f__Eubacteriaceae";
    }
    if($taxonomy[5] eq "g__Bacillus"){
      $taxonomy[4] = "f__Bacillaceae";
    }
    if($taxonomy[5] eq "g__Bacteroides"){
      $taxonomy[4] = "f__Bacteroidaceae";
      $taxonomy[3] = "o__Bacteroidales";
      $taxonomy[2] = "c__Bacteroidia";
      $taxonomy[1] = "p__Bacteroidetes";
    }
    if($taxonomy[5] eq "g__Mycoplasma"){
      $taxonomy[4] = "f__Mycoplasmataceae";
      $taxonomy[3] = "o__Mycoplasmatales";
      $taxonomy[2] = "c__Mollicutes";
      $taxonomy[1] = "p__Tenericutes";
    }
    if($taxonomy[5] eq "g__Clostridium"){
      $taxonomy[4] = "f__Clostridiaceae";
      $taxonomy[3] = "o__Clostridiales";
      $taxonomy[2] = "c__Clostridia";
      $taxonomy[1] = "p__Firmicutes";
    }
    if($taxonomy[4] eq "f__Clostridiaceae"){
      $taxonomy[3] = "o__Clostridiales";
      $taxonomy[2] = "c__Clostridia";
      $taxonomy[1] = "p__Firmicutes";
    }
    if($taxonomy[6] =~ m/s__Clostridium phytofermentans/){
      $taxonomy[5] = "g__Clostridium";
      $taxonomy[4] = "f__Clostridiaceae";
      $taxonomy[3] = "o__Clostridiales";
      $taxonomy[2] = "c__Clostridia";
      $taxonomy[1] = "p__Firmicutes";
    }
    if($taxonomy[5] eq "g__Oscillatoria"){
      $taxonomy[4] = "f__Oscillatoriaceae";
    }
    if($taxonomy[4] eq "f__Thermoanaerobacterales Family III. Incertae Sedis"){
      $taxonomy[3] = "o__Thermoanaerobacterales";
    }
    if($taxonomy[5] eq "g__Mycoplana"){
      $taxonomy[4] = "f__Rhizobiaceae";
      $taxonomy[3] = "o__Rhizobiales";
    }
    if($taxonomy[4] eq "f__Sinobacteraceae"){
      $taxonomy[3] = "o__Xanthomonadales";
    }


    ## EUKS
    if($taxonomy[0] =~ /k__eukaryota/i){
    
      # ==============================================
      # These fields are problematic in silva euks r117
      # when combined with greengenes-13 reference.
      if($taxonomy[3] =~ m/(o__chlorophyta)/i){
        $taxonomy[3] =  $1."_euk";
        $taxonomy[4] =  $taxonomy[4]."_euk";
        $taxonomy[5] =  $taxonomy[5]."_euk)";
      }
      if($taxonomy[4] =~ m/(f__trebouxiophyceae family)/i){
        $taxonomy[4] =  $1."_euk";
        $taxonomy[5] =  $taxonomy[5]."_euk";
      }
      if($taxonomy[5] =~ m/(g__trebouxiophyceae)/i){
        $taxonomy[5] =  $1."_euk_(".$taxonomy[4].")_euk";
        #$taxonomy[5] =  $1."_euk";
      }
      if($taxonomy[4] =~ m/(f__ulvophyceae\s+family)/i){
        $taxonomy[4] =  $1."_euk";
        $taxonomy[5] =  $taxonomy[5]."_euk";
      }
      if($taxonomy[3] =~ m/(o__cercozoa)/i){
        $taxonomy[3] =  $1."_euk";
        $taxonomy[4] =  $taxonomy[4]."_euk";
        $taxonomy[5] =  $taxonomy[5]."_euk)";
      }
      # ==============================================

      if($taxonomy[5] =~ m/(g__massisteria)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__rhizamoeba)/i){
      }
      if($taxonomy[5] =~ m/(g__unidentified)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__cercozoa)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__lagenidium)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__cercomonadida)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")"."_(".$taxonomy[3].")";
      }
      if($taxonomy[5] =~ m/(g__saprolegnia)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__oikomonas)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__chrysochromulina)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__graesiella)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
        print STDERR "[DEBUG] ".$1."_(".$taxonomy[4].")\n";
      }
      if($taxonomy[5] =~ m/(g__pythium)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__spumella)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__gloeotilopsis)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__muriella)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__chrysophyceae)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__dinophyceae)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__monodus)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__nor26)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__oblongichytrium)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[4] =~ m/(f__uncultured)/i){
        $taxonomy[4] =  $1."_(".$taxonomy[3].")";
      }
      if($taxonomy[5] =~ m/(g__uncultured)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__incertae)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__isochrysis)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__cf\.)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__tribonema)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__characium)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__hartmannella)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__halophytophthora)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__reticulitermes)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__bangia)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__clade)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__botrydiopsis)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__novel)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__pentapharsodinium)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__flabellulidae)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__porphyra)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__thraustochytriidae)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__eukaryote)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__pseudulvella)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__botryococcus)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__ochromonas)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__bicosoecida)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__rhizamoeba)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
      if($taxonomy[5] =~ m/(g__porotermes)/i){
        $taxonomy[5] =  $1."_(".$taxonomy[4].")";
      }
    }
  }elsif($fungi_ITS_correction == 1){
    #print STDERR "Doing ITS taxonomy correction\n";
 
    next if($taxonomy[5] =~ m/unidentified/);
    
		if($taxonomy[4] eq "f__Venturiaceae"){
			$taxonomy[3] = "o__Venturiales";
		}
		if($taxonomy[5] eq "g__Cryptococcus"){
			$taxonomy[5] = "g__Cryptococcus(".$taxonomy[4].")";
		}
		if($taxonomy[5] eq "g__Perenniporia"){
			$taxonomy[5] = "g__Perenniporia(".$taxonomy[4].")";
		}
		if($taxonomy[4] eq "f__Incertae sedis"){
			#$taxonomy[2] = $taxonomy[2]."(".$taxonomy[1]."-".$taxonomy[2].")";
			#$taxonomy[3] = $taxonomy[3]."(".$taxonomy[1]."-".$taxonomy[2].")";
			#$taxonomy[4] = $taxonomy[4]."(".$taxonomy[1]."-".$taxonomy[2].")";
      next;
		}
		if($taxonomy[4] eq "f__Amorphothecaceae" && $taxonomy[5] eq "g__Amorphotheca"){
	    $taxonomy[5] = "g__Amorphotheca(".$taxonomy[4].")";  
    }
		if($taxonomy[2] eq "c__Lecanoromycetes" && $taxonomy[3] eq "o__Incertae sedis"){
	    $taxonomy[3] = "o__Incertae sedis(".$taxonomy[2].")";  
    }
		if($taxonomy[2] eq "c__Sordariomycetes" && $taxonomy[3] eq "o__Incertae sedis"){
	    $taxonomy[3] = "o__Incertae sedis(".$taxonomy[2].")";  
    }
		if($taxonomy[2] eq "c__Dothideomycetes" && $taxonomy[3] eq "o__Incertae sedis"){
	    $taxonomy[3] = "o__Incertae sedis(".$taxonomy[2].")";  
    }
  
    # for plants and metazoa
		if($taxonomy[5] eq "g__Sphenopus"){
			$taxonomy[5] = "g__Sphenopus(".$taxonomy[4].")";
		}
  }

  $incertae_count++ if($taxonomy[1] =~ m/Incertae sedis/i);
  $incertae_count++ if($taxonomy[2] =~ m/Incertae sedis/i);
  $incertae_count++ if($taxonomy[3] =~ m/Incertae sedis/i);
  $incertae_count++ if($taxonomy[4] =~ m/Incertae sedis/i);
  $incertae_count++ if($taxonomy[5] =~ m/Incertae sedis/i);
  
  my $sum = ($k__ + $p__ + $c__ + $o__ + $f__ + $g__ + $s__);  

  if($sum >= $tax_level){
    $complete_array[$j] = $taxonomy[0]."\t".$taxonomy[1]."\t".$taxonomy[2]."\t".$taxonomy[3]."\t".$taxonomy[4]."\t".$taxonomy[5]."\t".$taxonomy[6]."\n" if $tax_level == 7;
    $complete_array[$j] = $taxonomy[0]."\t".$taxonomy[1]."\t".$taxonomy[2]."\t".$taxonomy[3]."\t".$taxonomy[4]."\t".$taxonomy[5]."\n" if $tax_level == 6;
    $complete_array[$j] = $taxonomy[0]."\t".$taxonomy[1]."\t".$taxonomy[2]."\t".$taxonomy[3]."\t".$taxonomy[4]."\n" if $tax_level == 5;
    $complete_array[$j] = $taxonomy[0]."\t".$taxonomy[1]."\t".$taxonomy[2]."\t".$taxonomy[3]."\n" if $tax_level == 4;
    $complete_array[$j] = $taxonomy[0]."\t".$taxonomy[1]."\t".$taxonomy[2]."\n" if $tax_level == 3;
    $complete_array[$j] = $taxonomy[0]."\t".$taxonomy[1]."\n" if $tax_level == 2;
    $complete_array[$j] = $taxonomy[0]."\n" if $tax_level == 1;
  
    print OUT_FASTA ">".$id." Root;".$taxonomy[0].";".$taxonomy[1].";".$taxonomy[2].";".$taxonomy[3].";".$taxonomy[4].";".$taxonomy[5].";".$taxonomy[6]."\n".$seq->seq()."\n" if $tax_level == 7;  
    print OUT_FASTA ">".$id." Root;".$taxonomy[0].";".$taxonomy[1].";".$taxonomy[2].";".$taxonomy[3].";".$taxonomy[4].";".$taxonomy[5]."\n".$seq->seq()."\n" if $tax_level == 6;  
    print OUT_FASTA ">".$id." Root;".$taxonomy[0].";".$taxonomy[1].";".$taxonomy[2].";".$taxonomy[3].";".$taxonomy[4]."\n".$seq->seq()."\n" if $tax_level == 5;  
    print OUT_FASTA ">".$id." Root;".$taxonomy[0].";".$taxonomy[1].";".$taxonomy[2].";".$taxonomy[3]."\n".$seq->seq()."\n" if $tax_level == 4;  
    print OUT_FASTA ">".$id." Root;".$taxonomy[0].";".$taxonomy[1].";".$taxonomy[2]."\n".$seq->seq()."\n" if $tax_level == 3;  
    print OUT_FASTA ">".$id." Root;".$taxonomy[0].";".$taxonomy[1]."\n".$seq->seq()."\n" if $tax_level == 2;  
    print OUT_FASTA ">".$id." Root;".$taxonomy[0]."\n" if $tax_level == 1;  
  
    $j++;
  }
}

# Dereplicate
$j=0;
foreach my $line (@complete_array){
  chomp($line);
  if(exists $hash{$line}){
  }else{
    $hash{$line} = $line; 
  }
  $j++;
}

@array = ();
$j=0;
my $i=0;
while ( my ($key, $value) = each(%hash) ) {
  $array[$j] =  $key;
  my @row = split("\t", $key);
  foreach(@row){
    if(exists $hash_id{$_}){
      #do nothin'
    }else{
      $hash_id{$_} = ($i+2);
      $i++;
    }
    #$i++;
  }
  $j++;
}

# Generate hierarcichal model
@array = sort(@array);
open(LIN, ">".$outfile_lineages) or die "Can't open file ".$outfile_lineages."\n";
print LIN $_."\n" foreach(@array);

my $prev_domain = "";
my $prev_phylum = "";
my $prev_class = "";
my $prev_order = "";
my $prev_family = "";
my $prev_genus = "";
my $prev_species = "";

my $domain_id;
my $phylum_id;
my $class_id;
my $order_id;
my $family_id;
my $genus_id;
my $species_id;

$j=0;
my %seen = ();

print OUT_MODEL "1*Root*0*0*Domain\n";

foreach(@array){
  
  my @row = split("\t", $_);
  #foreach my $t (@row){
  #  print $t."\t";
  #}
  #print "\n";
  if($row[0] ne $prev_domain){
    print OUT_MODEL $hash_id{$row[0]}."*".$row[0]."*1*1*Kingdom\n";
    $domain_id = $hash_id{$row[0]};
    $prev_domain = $row[0];
    $seen{$hash_id{$row[0]}}++;

  }
  next if(@row == 1);

  if($row[1] ne $prev_phylum){
    print OUT_MODEL $hash_id{$row[1]}."*".$row[1]."*".$domain_id."*2*Phylum\n";
    $phylum_id = $hash_id{$row[1]};
    $prev_phylum = $row[1];
    $seen{$hash_id{$row[1]}}++;

  }else{
  
  }
  next if(@row == 2);

  if($row[2] ne $prev_class){
    print OUT_MODEL $hash_id{$row[2]}."*".$row[2]."*".$phylum_id."*3*Class\n";
    $class_id = $hash_id{$row[2]};
    $prev_class = $row[2];
    $seen{$hash_id{$row[2]}}++;
  
  }else{
  
  }
  next if(@row == 3);

  if($row[3] ne $prev_order){
    print OUT_MODEL $hash_id{$row[3]}."*".$row[3]."*".$class_id."*4*Order\n";
    $order_id = $hash_id{$row[3]};
    $prev_order = $row[3];
    $seen{$hash_id{$row[3]}}++;

  }else{
  
  }
  next if(@row == 4);

  if($row[4] ne $prev_family){
    print OUT_MODEL $hash_id{$row[4]}."*".$row[4]."*".$order_id."*5*Family\n";
    $family_id = $hash_id{$row[4]};
    $prev_family = $row[4];
    $seen{$hash_id{$row[4]}}++;

  }else{
  
  }
  next if(@row == 5);

  if($row[5] ne $prev_genus){
    print OUT_MODEL $hash_id{$row[5]}."*".$row[5]."*".$family_id."*6*Genus\n";
    $genus_id = $hash_id{$row[5]};
    $prev_genus = $row[5];
    $seen{$hash_id{$row[5]}}++;
  
  }else{
  
  }
  next if(@row == 6);
  
  if($row[6] ne $prev_species){
    print OUT_MODEL $hash_id{$row[6]}."*".$row[6]."*".$genus_id."*7*Species\n";
    $species_id = $hash_id{$row[6]};
    $prev_species = $row[6];
    $seen{$hash_id{$row[6]}}++;
  
  }else{
  
  }
}

print STDERR "Incertae sedis: ".$incertae_count."\n" if($verbose);

close(OUT_MODEL);
exit;
