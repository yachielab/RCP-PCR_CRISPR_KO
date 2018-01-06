#!/usr/bin/env python

import glob, os
import time
import argparse
import sys

def main(args,script_path):
    out_dir = args.output_name

    if not os.path.isdir('%s'%(out_dir)):
        if not os.path.isdir('%s/Log_%s'%(out_dir,out_dir)):
            os.makedirs('%s/Log_%s'%(out_dir,out_dir))
            os.makedirs('%s/Log_%s/csv'%(out_dir,out_dir))
            os.makedirs('%s/Log_%s/pdf'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/fragmented_fasta'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/blast'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/blast/sh.blast'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/blast/out.blast'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/QC'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/QC/sh.identification'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/QC/out.identification'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/db'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/db/fasta'%(out_dir,out_dir))


    PATH = os.path.abspath(".")
    script_paths = script_path.split("rcppcr_ko.py")[0]

    target_regions = args.targets
    common = [["DBU1-primer","CCATACGAGCACATTACGGG"],["DBD2-primer","CTTGACTGAGCGACTGAGG"],["PS1.0-primer","TAACTTACGGAGTCGCTCTACG"],["PS2.0-primer","GGATGGGATTCTTTAGGTCCTG"]]
    with open(target_regions,"r") as F:
        for line  in F:
            c +=1
            cols = line.split(",")
            if c >1:
                common.append(["%s_Target"%(cols[0]),cols[1]])
                common.append(["%s_Frd"%(cols[0]),cols[1][:25]])
                common.append(["%s_Rvs"%(cols[0]),cols[1][-25:]])
    LL2csv()



    print sjsj
    print "Split and generate fasta files...."
    print "perl %sfastq2fasta.pl %s/workdir_%s %s %s\n"% (script_paths,out_dir,out_dir,args.input_file_R1, args.input_file_R2)
    os.system("perl %sfastq2fasta.pl %s/workdir_%s %s %s\n"% (script_paths,out_dir,out_dir,args.input_file_R1, args.input_file_R2))
    print '....... Finished\n\n'

    print "Generate sh.blast files...."
    print "perl %sprimers_blast_wrapper09212017DY.pl %s %s/%s %s/%s/fragmented_fasta/*\n"%(db,script_paths,PATH,out_dir,PATH,out_dir)
    os.system("perl %sprimers_blast_wrapper09212017DY.pl %s %s/%s %s/%s/fragmented_fasta/*"%(db,script_paths,PATH,out_dir,PATH,out_dir))
    print '....... Finished\n'

    print args.sge_computing

    """
    print "Combine sh.blast files....."
    print "perl codes/list_qsub.pl %s/Data/%s/blast/sh.primers_blast/* > Log_%s/sgeBLAST.sh\n"%(PATH,out_dir,out_dir)
    os.system("perl codes/list_qsub.pl %s/Data/%s/blast/sh.primers_blast/* > Log_%s/sgeBLAST.sh"%(PATH,out_dir,out_dir))
    print '....... Finished\n'

    print "Run sgeBLAST.sh files......"
    print "sh Log_%s/sgeBLAST.sh\n"%(out_dir)
    os.system("sh Log_%s/sgeBLAST.sh"%(out_dir))
    print '....... Finished\n'

    print "Generate sh.identification files......."
    print "perl codes/target-identification_wrapperDY.pl %s %s %s %s/Data/%s %s/Data/%s/fragmented_fasta/*\n"%(PATH,out_dir,PATH,out_dir)
    os.system("perl codes/target-identification_wrapperDY.pl %s %s %s %s/Data/%s %s/Data/%s/fragmented_fasta/*"%(PATH,out_dir,PATH,out_dir))
    print '....... Finished\n'
    my $primers_fasta  = '/home/t14905dy/projects/RCP-PCR/KO_clone/Data/db/fasta/const-seq.fna';
    my $bar2num_file   = '/home/t14905dy/projects/RCP-PCR/KO_clone/Data/bar2num.txt';
    my $targets        = '/home/t14905dy/projects/RCP-PCR/KO_clone/plate_tag_assignment/targets12192017.csv';

















    print "Combine sh.dentification files......"
    print "perl codes/list_qsub.pl %s/Data/%s/QC/sh.identification/* > Log_%s/QC.sh\n"%(PATH,out_dir,out_dir)
    os.system("perl codes/list_qsub.pl %s/Data/%s/QC/sh.identification/* > Log_%s/QC.sh"%(PATH,out_dir,out_dir))
    print '....... Finished\n'

    print "Running sgeIdentification.sh files......"
    print "sh Log_%s/QC.sh\n"%(out_dir)
    os.system("sh Log_%s/QC.sh"%(out_dir))
    print '....... Finished\n'

    print "Generate sh.readcounting ......"
    print "perl codes/read-counting_wrapper09212017DY.pl %s/Data/%s %s/Data/%s/QC/out.identification/*\n"%(PATH,out_dir,PATH,out_dir)
    os.system("perl codes/read-counting_wrapper09212017DY.pl %s/Data/%s %s/Data/%s/QC/out.identification/*"%(PATH,out_dir,PATH,out_dir))
    print '....... Finished\n'

    print "Combining sh.readcounting files......"
    print "perl codes/list_qsub.pl %s/Data/%s/QC/sh.count/* > Log_%s/ReadCount.sh\n"%(PATH,out_dir,out_dir)
    os.system("perl codes/list_qsub.pl %s/Data/%s/QC/sh.count/* > Log_%s/ReadCount.sh"%(PATH,out_dir,out_dir))
    print '....... Finished\n'

    print "Running readcount.sh files......"
    print "sh Log_%s/ReadCount.sh\n"%(out_dir)
    os.system("qsub Log_%s/ReadCount.sh"%(out_dir))
    print '....... Finished\n'

    print "Merging out.count files......"
    print "perl codes/merge_data.pl %s/Data/%s/QC/out.count/* > Log_%s/merged_count.dmp\n"%(PATH,out_dir,out_dir)
    os.system("perl codes/merge_data.pl %s/Data/%s/QC/out.count/* > Log_%s/merged_count.dmp"%(PATH,out_dir,out_dir))
    print '....... Finished\n'

    with open("Log_%s/merged_count.dmp"%(out_dir), "rt") as fin:
        with open("Log_%s/merged_count_forPy_temp.txt"%(out_dir), "wt") as fout:
            c=0
            for line in fin:
                c+=1
                if (c == 1):
                    fout.write("{")
                else:
                    fout.write(line.replace('=>', ':'))
    fin.close()
    fout.close()

    with open("Log_%s/merged_count_forPy_temp.txt"%(out_dir), "rt") as fin:
        with open("Log_%s/merged_count_forPy.txt"%(out_dir), "wt") as fout:
            c=0
            for line in fin:
                fout.write(line.replace(';', ''))
    fin.close()
    fout.close()
    os.system("rm Log_%s/merged_count_forPy_temp.txt"%(out_dir))
    #os.system("rm Log_%s/merged_count.dmp"%(out_dir))


    #print "Calling KO ......"
    #print "python Call_KO.py  Log_%s/merged_count_forPy.txt plate_tag_assignment/target_info_11262017.csv  Plus 0 Log_%s"%(out_dir,out_dir)
    #os.system("python codes/Call_KO.py  Log_%s/merged_count_forPy.txt plate_tag_assignment/target_info_11262017.csv  Plus 0 Log_%s"%(out_dir,out_dir))




    """
    print 'Finished running RCP-PCR program suite.'

def LL2csv(LL,name):
    with open(name,"w") as F:
        for L in LL:
            F.write("%s\n" %  (",").join( [str(i) for i in L]))
    F.close()




if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='rcppcr_ko')
    parser.add_argument('-R1','--input_file_R1', default=False,help='Input file of R1.fastq')
    parser.add_argument('-R2','--input_file_R2', default=False,help='Input file of R2.fastq')
    parser.add_argument('-t','--targets', default=False,help='Input target informtion in csv format. (see wiki for detail)')
    parser.add_argument('-out','--output_name', default="output" )
    parser.add_argument('-r','--ratio', type=int, help='Minimum threashold (0 < ratio < 0.5 ) to call mutation profile',default=0.1)
    parser.add_argument('-c','--core_num', type=int, help='Number of cores for multi-processing on local computer.',default=1)
    parser.add_argument('-sge','--sge_computing', type=int, help='1 if computing on SGE computers.',default=0)
    args = parser.parse_args()
    print args
    error = 0
    if args.input_file_R1 is False:
        print "Input file error. Please input R1.fastq files as input. Use option; [-R1 example_R1.fastq]."
        error +=1
    if args.input_file_R2 is False:
        print "Input file error. Please input R2.fastq files as input. Use option; [-R1 example_R2.fastq]."
        error +=1

    if args.targets is False:
        print "Input file error. Please input .csv files as target input. Use option; [-t example.csv]."
        error +=1
    if error > 0:
        quit()


    script_path = os.path.realpath(__file__)
    main(args,script_path)
