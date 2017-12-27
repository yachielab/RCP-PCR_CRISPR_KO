#!/usr/bin/env python

import sys
import glob, os
import time

import argparse

parser = argparse.ArgumentParser(description='argparse sample.')
parser.add_argument('-in','--input_file', action='store_true', default=False)
parser.add_argument('-out','--output_file', action='store_true', default=False)

parser.add_argument('-r','--ratio', type=int, help='Minimum threashold (0 < ratio < 0.5 ) to call mutation profile')
parser.add_argument('-c','--core_num', type=int, help='Number of cores for multi-processing.')
args = parser.parse_args()

print args













def main(seq_dir):
    if not os.path.isdir('Log_%s'%(seq_dir)):
        os.makedirs('Log_%s'%(seq_dir))
        os.makedirs('Log_%s/csv'%(seq_dir))
        os.makedirs('Log_%s/pdf'%(seq_dir))

    PATH = os.path.abspath(".")
    Miseqrunname = '12192017-yachielab'
    #print "Split and generate fasta files...."
    #print "perl codes/fastq2fasta.pl Data/%s Data/%s/miseq_data/%s_S1_L001_R1_001.fastq Data/%s/miseq_data/%s_S1_L001_R2_001.fastq\n"% (seq_dir,seq_dir,Miseqrunname,seq_dir,Miseqrunname)
    #os.system("perl codes/fastq2fasta.pl Data/%s Data/%s/miseq_data/%s_S1_L001_R1_001.fastq Data/%s/miseq_data/%s_S1_L001_R2_001.fastq"%(seq_dir,seq_dir,Miseqrunname,seq_dir,Miseqrunname))
    #print '....... Finished\n\n'
    

    #print "Generate sh.blast files...."
    #print "perl codes/primers_blast_wrapper09212017DY.pl %s/Data/%s %s/Data/%s/fragmented_fasta/*\n"%(PATH,seq_dir,PATH,seq_dir)
    #os.system("perl codes/primers_blast_wrapper09212017DY.pl %s/Data/%s %s/Data/%s/fragmented_fasta/*"%(PATH,seq_dir,PATH,seq_dir))
    #print '....... Finished\n'

    #print "Combine sh.blast files....."
    #print "perl codes/list_qsub.pl %s/Data/%s/blast/sh.primers_blast/* > Log_%s/sgeBLAST.sh\n"%(PATH,seq_dir,seq_dir)
    #os.system("perl codes/list_qsub.pl %s/Data/%s/blast/sh.primers_blast/* > Log_%s/sgeBLAST.sh"%(PATH,seq_dir,seq_dir))
    #print '....... Finished\n'
  
    #print "Run sgeBLAST.sh files......"
    #print "sh Log_%s/sgeBLAST.sh\n"%(seq_dir)
    #os.system("sh Log_%s/sgeBLAST.sh"%(seq_dir))
    #print '....... Finished\n'

    #print "Generate sh.identification files......."
    #print "perl codes/target-identification_wrapperDY.pl %s/Data/%s %s/Data/%s/fragmented_fasta/*\n"%(PATH,seq_dir,PATH,seq_dir)
    #os.system("perl codes/target-identification_wrapperDY.pl %s/Data/%s %s/Data/%s/fragmented_fasta/*"%(PATH,seq_dir,PATH,seq_dir))
    #print '....... Finished\n'

    #print "Combine sh.dentification files......"
    #print "perl codes/list_qsub.pl %s/Data/%s/QC/sh.identification/* > Log_%s/QC.sh\n"%(PATH,seq_dir,seq_dir)
    #os.system("perl codes/list_qsub.pl %s/Data/%s/QC/sh.identification/* > Log_%s/QC.sh"%(PATH,seq_dir,seq_dir))
    #print '....... Finished\n'

    #print "Running sgeIdentification.sh files......"
    #print "sh Log_%s/QC.sh\n"%(seq_dir)
    #os.system("sh Log_%s/QC.sh"%(seq_dir))
    #print '....... Finished\n'

    #print "Generate sh.readcounting ......"
    #print "perl codes/read-counting_wrapper09212017DY.pl %s/Data/%s %s/Data/%s/QC/out.identification/*\n"%(PATH,seq_dir,PATH,seq_dir)
    #os.system("perl codes/read-counting_wrapper09212017DY.pl %s/Data/%s %s/Data/%s/QC/out.identification/*"%(PATH,seq_dir,PATH,seq_dir))
    #print '....... Finished\n'

    #print "Combining sh.readcounting files......"
    #print "perl codes/list_qsub.pl %s/Data/%s/QC/sh.count/* > Log_%s/ReadCount.sh\n"%(PATH,seq_dir,seq_dir)
    #os.system("perl codes/list_qsub.pl %s/Data/%s/QC/sh.count/* > Log_%s/ReadCount.sh"%(PATH,seq_dir,seq_dir))
    #print '....... Finished\n'

    #print "Running readcount.sh files......"
    #print "sh Log_%s/ReadCount.sh\n"%(seq_dir)
    #os.system("qsub Log_%s/ReadCount.sh"%(seq_dir))
    #print '....... Finished\n'

    print "Merging out.count files......"
    print "perl codes/merge_data.pl %s/Data/%s/QC/out.count/* > Log_%s/merged_count.dmp\n"%(PATH,seq_dir,seq_dir)
    os.system("perl codes/merge_data.pl %s/Data/%s/QC/out.count/* > Log_%s/merged_count.dmp"%(PATH,seq_dir,seq_dir))
    print '....... Finished\n'
    
    with open("Log_%s/merged_count.dmp"%(seq_dir), "rt") as fin:
        with open("Log_%s/merged_count_forPy_temp.txt"%(seq_dir), "wt") as fout:
            c=0
            for line in fin:
                c+=1
                if (c == 1):
                    fout.write("{")
                else:
                    fout.write(line.replace('=>', ':'))
    fin.close()
    fout.close()
    
    with open("Log_%s/merged_count_forPy_temp.txt"%(seq_dir), "rt") as fin:
        with open("Log_%s/merged_count_forPy.txt"%(seq_dir), "wt") as fout:
            c=0
            for line in fin:
                fout.write(line.replace(';', ''))
    fin.close()
    fout.close()
    os.system("rm Log_%s/merged_count_forPy_temp.txt"%(seq_dir))
    #os.system("rm Log_%s/merged_count.dmp"%(seq_dir))
    

    #print "Calling KO ......"
    #print "python Call_KO.py  Log_%s/merged_count_forPy.txt plate_tag_assignment/target_info_11262017.csv  Plus 0 Log_%s"%(seq_dir,seq_dir)
    #os.system("python codes/Call_KO.py  Log_%s/merged_count_forPy.txt plate_tag_assignment/target_info_11262017.csv  Plus 0 Log_%s"%(seq_dir,seq_dir))





    print 'Finished running RCP-PCR program suite.'
    

    


if __name__ == '__main__':
    main(sys.argv[1]) #Input directory name eg. template.MiSeq
