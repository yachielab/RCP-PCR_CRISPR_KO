## RCP-PCR_CRISPR_KO
Analysis tool for RCP-PCR experiments on CRISPR/Cas9 clonal knock out experiments.
This tool is for analyzing on local/cluster computers. 

Current Version: 1.2

Release Date: January 7, 2017

Platform: Tested on Linux x64 / MacOSX system

Please contact dan.yamamoto.evans [at] gmail.com for quick response to resolve any bug or feature update.

## Installation

Please follow the following steps to install RCP-PCR_KO from source:

>Clone RCP-PCR_CRISPR_KO source code: git clone https://github.com/DanYamamotoEvans/RCP-PCR_CRISPR_KO.git


### Requirements
Python version 2.7+ (2.7 reccomended)   
Perl version 5  
R version 3+   
BLAST+ (blastn version 2.6.0+)  

## Input specifications

Use the following options to run rcppcr_ko:

usage:  


    PULLPATH/rcppcr_ko/rcppcr_ko.py  
                    [-h] [-R1 INPUT_FILE_R1] [-R2 INPUT_FILE_R2] [-t TARGETS]  
                    [-out OUTPUT_NAME] [-r RATIO] [-c CORE_NUM]  
                    [-sge SGE_COMPUTING]    


optional arguments:  

	-R1 INPUT_FILE_R1, --input_file_R1 INPUT_FILE_R1    
			Input fastq file of read1 (eg. R1.fastq)  
    
	-R2 INPUT_FILE_R2, --input_file_R2 INPUT_FILE_R2   
        	Input fastq file of read2 (eg. R2.fastq)  
    
	-t TARGETS, --targets TARGETS   
		Input target informtion in csv format. (see wiki for detail)  
    
	-out OUTPUT_NAME, --output_name OUTPUT_NAME   
   
	-r RATIO, --ratio RATIO    
		Minimum threashold (0 < ratio < 0.5 ) to call mutation profile (Default = 0.1).   
      
	-c CORE_NUM, --core_num CORE_NUM   
		Number of cores for multi-processing.  
        
	-sge SGE_COMPUTING, --sge_computing SGE_COMPUTING   
		1 if computing on SGE computers.  
        
	-h, --help    
		Show this help message.  
        
## Usage examples
    python ~/GitHub/RCP-PCR_CRISPR_KO/rcppcr_ko/rcppcr_ko.py -R1 PULLPATH/test/test_R1.fastq -R2 PULLPATH/test/test_R2.fastq  -t test/test_target.csv -c 2   

## Contact information

Please send your comments or bug reports to dan.yamamoto.evans [at] gmail.com  
