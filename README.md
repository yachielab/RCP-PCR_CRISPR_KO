# RCP-PCR_CRISPR_KO: Analysis tool for RCP-PCR experiments on CRISPR/Cas9 clonal knock out experiments.
Current Version: 0.20
Release Date: December 27, 2017
Platform: Linux x64 system

Please contact dan.yamamoto.evans@gmail.com for quick response to resolve any bug or feature update.

## Installation

Please follow the following steps to install NanoBLASTer from source:

Clone NanoBLASTer source code: git clone https://github.com/ /RCP-PCR_CRISPR_KO.git
Go to the NanoBLASTer source directory: cd RCP-PCR_CRISPR_KO/src 
Build the NanoBLASTer project: make

## Input specifications

Use the following options to run NanoBLASTer:
-C: To specify one of the Parameters: -C10, -C25, or -C50
-r: To specify the name of Reference file (FASTA format)
-i: To specify the name of Reads file (FASTA format)
-o: To specify the prefix of Output file
-k: To specify the size of KMER
-a: To specify the size of ANCHOR
-l: To specify the min number of Clusters
-s: To run the program at higher sensitivity
-n: To specify the Number of reads to be aligned
-g: To specify the interval (or Gap) length between KMERs
-X: To configure NanoBLASTer for less memory using Single index
-h, or -?: To print this Help information.

## Usage examples


## Contact information

Please send your comments or bug reports to dan.yamamoto.evans@gmail.com
