#mkdir -p ~/DBs/unite/v2/RDP
#mkdir -p ~/DBs/unite/v2/RDP/genus

# Fetch sequences form 
#wget http://unite.ut.ee/sh_files/sh_qiime_release_13.05.2014.zip
#unzip sh_qiime_release_13.05.2014.zip

# Merge taxonomy with fasta:
#~/build/tools/prepareQiimeITS.pl --fasta --taxonomy > seqsWithTax.fasta



#cat seqsWithTax.fasta ITS1-ITS2_Metazoa_Viriplantae.fasta > seqsWithTaxAndPlantsMetazoa.fasta
# replace all ë with e.
FASTA=./seqsWithTaxAndPlantsMetazoa.fasta

mkdir -p ./RDP
mkdir -p ./RDP/genus

#Generate hierarchical .txt files and corresponding fasta files
~/build/RDP-training-sets/scripts/generateRDPHierTree.pl \
 --fasta $FASTA \
 --tax_level 6 \
 --outfile_model ./RDP/genus/genus.txt \
 --outfile_fasta ./RDP/genus/genus.fasta \
 --log ./RDP/genus/genus.log \
 --outfile_lineages ./RDP/genus/lineages.tab \
 --fungi_ITS_correction 1 \


module load mugqic/RDP_classifier/2.5;
RDPBIN=`which rdp_classifier.jar`
java -Xmx2g -cp $RDPBIN edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker \
 ./RDP/genus/genus.txt \
 ./RDP/genus/genus.fasta 1 version1 test \
 ./RDP/genus \
