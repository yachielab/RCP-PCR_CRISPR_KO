## RCP-PCR CRISPR KO
This program suit is implemented for genotyping isolated single clones from a CRISPR-treated population using RCP-PCR.
Information regarding detailed experimental design and analysis scripts are available below.

### Experiment design
We treated human cell culture samples with CRISPR/Cas9, targeting specific genes of interest. The cells were isolated in 96-well culture plates to obtain clonal isolates. After growth of the cells, we performed serial PCRs to identify the genotype.

RCP(Row Column Plate)-PCR is capable of multiplexing up to 16 (Rows) x 24 (Columns) x 35 (Plates) = ~13,000 samples to be sequenced from a single high-throughput short read sequencing run ([Yachie et al](http://msb.embopress.org/content/12/4/863)). Here, we use this strategy to genotype clone derived cell samples with a minor modification in the design (; we performed a PCR prior to RCP-PCR to avoid ordering RCP-PCR primer sets for each target).   

For each of the targets, we designed primers which would amplify a target region with a length less than 150 bp. We performed the first PCR with primers that flank a common RC-PCR primer landing site. The samples were purified in a 96-well format, and products were used as template for subsequent RCP-PCR procedures.
The second PCR (Row Column-PCR) was performed with FWD primers with a Row index, and REV primers with a Column index, respectively. These RC-PCR primers were flanked with binding sites for the third PCR (Plate-PCR). Products from the second (Row Column-PCR) PCR were grouped together and purified to use for the third PCR (Plate-PCR), and sequenced in a high-throughput sequencing run.



### Analysis script design
We first identified the RCP-PCR indices and target region from each sequencing read. Simultaneously, the target region of the read is aligned to the wild-type sequence using the BLASTn software. The alignment result for de-multiplexed reads were aggregated to compute the allele frequency within each of the clone-derived sample. Any allele appearing less than 10% within each de-multiplexed sample was eliminated for further analysis to cancel out sequencing error. Finally, the remaining count data of each allele was used to identify the genotype of the isolated clones. Clones with frame-shift mutation on both alleles were used for further analysis.

### Reference
- [Yachie et al (2016) Molecular Systems Biology](http://msb.embopress.org/content/12/4/863)
- Suzuki et al (In preparation)



## Version

Current Version: 1.2

Release Date: January 7, 2017

Platform: Tested on Linux x64 / MacOSX system

Please contact dan.yamamoto.evans [at] gmail.com for quick response to resolve any bugs or for feature updates.

## Installation

Clone RCP-PCR_CRISPR_KO source code:

    git clone https://github.com/yachielab/RCP-PCR_CRISPR_KO.git


### Requirements
Python version 2.7+ (2.7 reccomended)   
Perl version 5  
R version 3+   
BLAST+ (blastn version 2.6.0+)  

Before you use;  
Check that your BLAST+ program path is through. Go to NCBI website and download BLAST+ locally if not present. Refer to the BLAST® Command Line Applications User Manual published by NCBI (https://www.ncbi.nlm.nih.gov/books/NBK279690/).  

## Software usage

Use the following options to run rcppcr_ko:

usage:  


    FULLPATH/rcppcr_ko/rcppcr_ko.py  
                    [-h] [-R1 INPUT_FILE_R1] [-R2 INPUT_FILE_R2] [-t TARGETS]  
                    [-out OUTPUT_NAME] [-r RATIO] [-c CORE_NUM]  
                    [-sge SGE_COMPUTING]    


Required arguments:  

	-in INPUT_FILE_R1, --input
			Input fastq files  (eg. test_R*.fastq)  

	-t TARGETS, --targets TARGETS   
		Input target informtion in csv format. (see wiki for detail)  

optional arguments:  

	-out OUTPUT_NAME, --output_name OUTPUT_NAME   

	-r RATIO, --ratio RATIO    
		Minimum threshold (0 < ratio < 0.5 ) to call mutation profile (Default = 0.1).   

	-c CORE_NUM, --core_num CORE_NUM   
		Number of cores for multi-processing.  

	-sge SGE_COMPUTING, --sge_computing SGE_COMPUTING   
		1 if computing on SGE computers.  

	-h, --help    
		Show this help message.  


## Command example for test data
    python FULLPATH/RCP-PCR_CRISPR_KO/rcppcr_ko/rcppcr_ko.py -in FULLPATH/test/"*.fastq"  -t FULLPATH/test/test_target.csv -c 2    
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



## Example output summary file

#### Output file (.csv)

    Column	Name					Description
    1	Target					Target 	name
    2	Plate					Plate 	index
    3	Row					Row 	index
    4	Col					Column 	index
    5	Total_reads				Total reads after de-multiplexing
    6	#Mutation_profiles_above_threashold	Number of mutation profiles (allelles) observed above threshold after demultiplexing
    7	Well_KO_stat				Genotype summary of clonal sample. This is based on each of the profiles (columns 9-16). Criteria for specific statement is shown below.
    8	Reads_support_stat(%)			Percentage of demultiplexed reads assigned to allelles above threshold.
    9	Profile#1				Indel summary of the most dominant allelles.
    10	#Reads#1_rate				Ratio of Profile#1 reads within total reads in demultiplexed sample.
    11	Standard_Error#1			Standard error of Profile#1 reads based on read number.
    12	Profile#1_deatl				BLASTn btop string of the most dominant allele.
    13	Profile#2				Indel summary of the second most dominant allelles.
    14	#Reads#2_rate				Ratio of Profile#2 reads within total reads in demultiplexed sample.
    15	Standard_Error#2			Standard error of Profile#2 reads based on read number.
    16	Profile#2_detail			BLASTn btop string of the second most dominant allele.


#### Well KO stat detail

    Well_KO_stat			Detail
    -				No allelles were above threshold.
    WT				All alleles were same as wild type.
    Indel(1profile)			One allele was above threshold, and it was an in-frame indel.
    Frameshift(1profile)		One allele was above threshold, and it was a frameshift mutation.
    Hetero				Two alleles were above threshold, and one was wildtype, another was a frameshift mutation.
    Homo-indel(non-frameshift)	Two alleles were above threshold, and both alleles were in-frame indels
    Homo-frameshift			Two alleles were above threshold, and both alleles were frameshift mutations.
    Too many profiles		More than two alleles was above threshold.




## Contact information

Please send your comments or bug reports to dan.yamamoto.evans [at] gmail.com  
