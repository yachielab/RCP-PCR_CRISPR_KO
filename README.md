## RCP-PCR_CRISPR_KO
This program suit is implemented for genotyping isolated single clones from a CRISPR-treated population using RCP-PCR. 
You may find detials designs of the experiments and analysis scripts below.

### Experiment design
We have treated human cell culture samples with CRISPR/Cas9, targeting specific genes of interest. The cells were isolated in 96-well culture plates to obtain clonal isolates. After growth of the cells, we performed serial PCRs to identify the genotype.

RCP(Row Column Plate)-PCR is capable of multiplexing up to 16 (Rows) x 24 (Columns) x 35 (Plates) = ~13,000 samples to subject on a single high throughput short read sequencing run ([Yachie et al (2016) Molecular Systems Biology](http://msb.embopress.org/content/12/4/863)). Here, we use this stragety to genotype clone derived cell samples with a minor modification in the design (; we performed a PCR prior to RCP-PCR in order to prevent ordering RCP-PCR primer sets for each target).   

For each of the targets, we designed primers which would amplify the region with length less than 150 bp. We performed the first PCR with those primers which flanks a common RC-PCR primer landing site. The samples were purified in a 96-well format, and products were used as template for subsequent RCP-PCR procedures. 
The second PCR (Row Column-PCR) was performed with Frd primers having Row index, and Rvs primer having Column index, respectively. These RC-PCR primers were flanked with binding sites for the third PCR (Plate-PCR). Products from the second (Row Column-PCR) PCR were grouped together and purified to use for the third PCR (Plate-PCR), and sequenced in a high-throughput sequencing run. 
   
![RCP-PCR](https://www.embopress.org/cms/asset/fd513902-3d16-43ea-b723-fd3e602b8f59/msb156660-fig-0003ev-m.jpg)
Figure EV3 in [Yachie et al (2016) Molecular Systems Biology](http://msb.embopress.org/content/12/4/863).




### Analysis script design
We 


### Reference
- Yachie et al
- Suzuki et al (In preparation)
    
   




## Version

Current Version: 1.2

Release Date: January 7, 2017

Platform: Tested on Linux x64 / MacOSX system

Please contact dan.yamamoto.evans [at] gmail.com for quick response to resolve any bug or feature update.

## Installation

Clone RCP-PCR_CRISPR_KO source code: 

    git clone https://github.com/DanYamamotoEvans/RCP-PCR_CRISPR_KO.git


### Requirements
Python version 2.7+ (2.7 reccomended)   
Perl version 5  
R version 3+   
BLAST+ (blastn version 2.6.0+)  

Before you use;  
Check that your BLAST+ program path is through. Go to NCBI website and download BLAST+ locally if not present. Reffer the BLAST® Command Line Applications User Manul published by NCBI (https://www.ncbi.nlm.nih.gov/books/NBK279690/).  

## Input specifications

Use the following options to run rcppcr_ko:

usage:  


    PULLPATH/rcppcr_ko/rcppcr_ko.py  
                    [-h] [-R1 INPUT_FILE_R1] [-R2 INPUT_FILE_R2] [-t TARGETS]  
                    [-out OUTPUT_NAME] [-r RATIO] [-c CORE_NUM]  
                    [-sge SGE_COMPUTING]    


Required arguments:  

	-R1 INPUT_FILE_R1, --input_file_R1 INPUT_FILE_R1    
			Input fastq file of read1 (eg. R1.fastq)  
    
	-R2 INPUT_FILE_R2, --input_file_R2 INPUT_FILE_R2   
        	Input fastq file of read2 (eg. R2.fastq)  
    
	-t TARGETS, --targets TARGETS   
		Input target informtion in csv format. (see wiki for detail)  
    
optional arguments:  

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
    python FULLPATH/RCP-PCR_CRISPR_KO/rcppcr_ko/rcppcr_ko.py -R1 FULLPATH/test/test_R1.fastq -R2 FULLPATH/test/test_R2.fastq  -t PULLPATH/test/test_target.csv -c 2    
>Change 'FULLPATH' accordingly to your path to the directory.


## Example of reference file
Target reference file: Available in test/test_target.csv

    Target,Target_seq,gRNA_s,gRNA_e
    MED4_sg2,ATTAAGTGCCAATTTCACAGTC..AATAAATCAGACAATAGACT,38,57
    MED4_sg3,GGCTAAAGGATCTTGTGAATAG..NNAAGGAGAAAGGTTAGTAT,90,109


Target referece: ID of the target. Please put a unique name.

Target seq     : Target (amplified) sequence. 

gRNA_s         : Position where gRNA starts.

gRNA_e         : Position where gRNA ends.



## Example of output files








## Contact information

Please send your comments or bug reports to dan.yamamoto.evans [at] gmail.com  
