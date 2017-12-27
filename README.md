## RCP-PCR_CRISPR_KO
Analysis tool for RCP-PCR experiments on CRISPR/Cas9 clonal knock out experiments.
This tool is for analyzing on local computers. Contact the developer for codes for cluster servers.

Current Version: 1.1

Release Date: December 27, 2017

Platform: Linux x64 system

Please contact dan.yamamoto.evans@gmail.com for quick response to resolve any bug or feature update.

## Installation

Please follow the following steps to install RCP-PCR_KO from source:

Clone RCP-PCR_CRISPR_KO source code: git clone https://github.com/DanYamamotoEvans/RCP-PCR_CRISPR_KO.git
Go to the source directory: cd RCP-PCR_CRISPR_KO/rcppcr_ko
Setup the codes: 

## Input specifications

Use the following options to run rcppcr_ko:

parser.add_argument('-in','--input_file', action=\
'store_true', default=False)
parser.add_argument('-out','--output_file', actio\
n='store_true', default=False)

parser.add_argument('-r','--ratio', type=int, hel\
p='Minimum threashold (0 < ratio < 0.5 ) to call \
mutation profile')
parser.add_argument('-c','--core_num', type=int, \
help='Number of cores for multi-processing.')




-in[--input_file]: To specify input fastq files (unprocessed).
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
