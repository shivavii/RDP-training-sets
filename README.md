RDP-training-sets
=================

This repo hosts scripts related to RDP training set generation and formatting

The RDP classifier software is a widespread tool to classify clusters (aka OTUs) obtained from short rRNA amplicons sequencing data. Although OTUs can be classified using the web page of the RDP classifier consortium (http://rdp.cme.msu.edu/classifier/classifier.jsp), anyone who gets serious on analyzing rRNA amplicon data will, at some point, use that software from the command line (http://sourceforge.net/projects/rdp-classifier/). The RDP classifier package contains a function to train your own data sets (i.e. curated rRNA fasta file), but lacks scipts to generate intermediate files needed to generate these training sets. This repo's mission is to provide scripts facilitating RDP classifier training sets generation.

contacts:
Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
