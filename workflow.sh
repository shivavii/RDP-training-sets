################################
## CHLOROPLAST AND MITOCHONDRIA
##chloroplasts were fetched from http://chloroplast.cbio.psu.edu/organism.cgi

## We modified the lineages of that db and further modified it with the script below
## to make it compatible with the lineages from the greengenes database.
## ~/build/RDP-training-sets/modifyChloroDb.pl --infile chloro_raw.fasta > chloro.fasta
## chloro.fasta


## Mitochondria were fetched from http://megasun.bch.umontreal.ca/ogmp/projects/other/mt_list.html
## only a subset of 22 sequences were kept. We then manually constructed taxonomy.
## mito.fasta

################################
## SILVA EUKS rRNA SSU
## Download Silva SSU euks db and fix lineages.
~/build/RDP-training-sets/fixSilvaDb.pl --infile_fasta silva_r117.fasta > silva_r117_fixed.fasta #(euks only). 

################################
## GREENGENES.

## Download greengenes 13-5 release from Second Genome site, add taxonomy and "fix" lineages.
~/build/RDP-training-sets/addTaxToGGFasta.pl --infile_fasta ./gg_13_5.fasta --infile_tax ./gg_13_5_taxonomy > ./gg_13_5_withTaxonomy.fasta
~/build/RDP-training-sets/fillGGEmptyLineagesSpecie.pl --infile ./gg_13_5_withTaxonomy.fasta --outfile ./gg_13_5_withTaxonomy_corrected.fasta --outfilef ./gg_13_5_withTaxonomy_corrected.fasta.failed

################################
## Concatenate all 3 fasta files with their "fixed" lineages.
cat ./chloro.fasta ./mito.fasta ./silva_r117_fixed.fasta ./gg_13_5_withTaxonomy_corrected.fasta > gg_13_5_withTaxonomy_corrected_silva_r117_euk_mito_chloro.fasta

################################
# TRAIN SET FOR SPECIE LEVEL. For long reads, Sanger, PacBio.

mkdir -p specie
~/build/RDP-training-sets/generateRDPHierTree.pl \
 --fasta ./gg_13_5_withTaxonomy_corrected_specie_silva_r117_euk_mito_chloro.fasta \
 --tax_level 7 \
 --outfile_model ./specie.txt \
 --outfile_fasta ./specie.fasta \
 --log ./specie.log \
 --outfile_lineages ./lineages.tab \
 --greengenes_correction 1 \
 --verbose \

## Generate RDP training model using the files generated in the previous step
module load mugqic/RDP_classifier/2.5; 
RDPBIN=`which rdp_classifier.jar`
java -Xmx3g -cp $RDPBIN edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker \
 ./specie.txt \
 ./specie.fasta 1 version1 test \
 ./specie/ \

################################
# TRAIN SET FOR GENUS LEVEL. (short amplicons should never be classifier lower than the genus level.

mkdir -p genus
~/build/RDP-training-sets/generateRDPHierTree.pl \
 --fasta $DIR/gg_13_5_withTaxonomy_corrected_silva_r117_euk_mito_chloro.fasta \
 --tax_level 6 \
 --outfile_model $DIR/genus.txt \
 --outfile_fasta $DIR/genus.fasta \
 --log $DIR/genus.log \
 --outfile_lineages $DIR/lineages.tab   \
 --greengenes_correction 1 \
 --verbose \

##Generate RDP training model using the files generated in the previous step
module load mugqic/RDP_classifier/2.5; 
RDPBIN=`which rdp_classifier.jar`
java -Xmx3g -cp $RDPBIN edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker \
 $DIR/genus.txt \
 $DIR/genus.fasta 1 version1 test \
 $DIR/genus/ \




