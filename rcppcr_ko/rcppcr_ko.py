#!/usr/bin/env python
import glob, os
import time
import argparse
import sys
import multiprocessing as mp
import operator


def main(args,script_path):
    PATH = os.path.abspath(".")
    script_paths = script_path.split("rcppcr_ko.py")[0]


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
            os.makedirs('%s/workdir_%s/QC/sh.count'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/QC/out.count'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/db'%(out_dir,out_dir))
            os.makedirs('%s/workdir_%s/db/fasta'%(out_dir,out_dir))



    target_regions = args.targets
    common = [["DBU1-primer","CCATACGAGCACATTACGGG"],["DBD2-primer","CTTGACTGAGCGACTGAGG"],["PS1.0-primer","TAACTTACGGAGTCGCTCTACG"],["PS2.0-primer","GGATGGGATTCTTTAGGTCCTG"]]
    with open(target_regions,"r") as F:
        c = 0
        for line  in F:
            c +=1
            cols = line.split(",")
            if c >1:
                common.append(["%s_Target"%(cols[0]),cols[1]])
                common.append(["%s_Frd"%(cols[0]),cols[1][:25]])
                common.append(["%s_Rvs"%(cols[0]),rev_comp(cols[1][-25:])])
    LL2 = []
    for i in common:
        LL2.append(["c%s"%(i[0]),rev_comp(i[1])])
    common += LL2
    LL2fna(common,'%s/workdir_%s/db/fasta/const-seq.fna'%(out_dir,out_dir))

    print "Making BLAST+ database....."
    os.system("makeblastdb -in %s/workdir_%s/db/fasta/const-seq.fna -dbtype nucl"%(out_dir,out_dir))
    db = '%s/workdir_%s/db/fasta/const-seq.fna'%(out_dir,out_dir)


    print "\n\n\n\nSplit and generating fasta files...."
    print "perl %sfastq2fasta.pl %s/workdir_%s %s %s\n"% (script_paths,out_dir,out_dir,args.input_file_R1, args.input_file_R2)
    os.system("perl %sfastq2fasta.pl %s/workdir_%s %s %s\n"% (script_paths,out_dir,out_dir,args.input_file_R1, args.input_file_R2))
    print '....... Finished\n\n'

    print "Generate sh.blast files...."
    print "perl %sprimers_blast_wrapper09212017DY.pl %s %s/%s/workdir_%s %s/%s/workdir_%s/fragmented_fasta/*\n"%(script_paths,db,PATH,out_dir,out_dir,PATH,out_dir,out_dir)
    os.system("perl %sprimers_blast_wrapper09212017DY.pl %s %s/%s/workdir_%s %s/%s/workdir_%s/fragmented_fasta/*"%(script_paths,db,PATH,out_dir,out_dir,PATH,out_dir,out_dir))
    print '....... Finished\n'

    if args.sge_computing < 1:
        #print "%s/%s/workdir_%s/blast/sh.blast/"%(PATH,out_dir,out_dir)
        subprocces_sh(args.core_num,"%s/%s/workdir_%s/blast/sh.blast/"%(PATH,out_dir,out_dir))
    else:
        print "Combine sh.blast files....."
        print "perl %slist_qsub.pl %s/%s/workdir_%s/blast/sh.blast/* > %s/%s/Log_%s/sgeBLAST.sh\n"%(script_paths,PATH,out_dir,out_dir,PATH,out_dir,out_dir)
        os.system("perl %slist_qsub.pl %s/%s/workdir_%s/blast/sh.blast/* > %s/%s/Log_%s/sgeBLAST.sh"%(script_paths,PATH,out_dir,out_dir,PATH,out_dir,out_dir))
        print '....... Finished\n'

        print "Run sgeBLAST.sh files......"
        print "sh %s/%s/Log_%s/sgeBLAST.sh\n"%(PATH,out_dir,out_dir)
        os.system("sh %s/%s/Log_%s/sgeBLAST.sh"%(PATH,out_dir,out_dir))
        print '....... Finished\n'

    print "Generate sh.identification files......."
    print "perl %starget-identification_wrapperDY.pl %s %s/%s/workdir_%s/db/fasta/const-seq.fna %sbar2num.txt %s %s/%s/workdir_%s %s/%s/workdir_%s/fragmented_fasta/*\n"%(script_paths,script_paths,PATH,out_dir,out_dir,script_paths,args.targets,  PATH,out_dir,out_dir,PATH,out_dir,out_dir)
    os.system("perl %starget-identification_wrapperDY.pl %s %s/%s/workdir_%s/db/fasta/const-seq.fna %sbar2num.txt %s %s/%s/workdir_%s %s/%s/workdir_%s/fragmented_fasta/*"%(script_paths,script_paths,PATH,out_dir,out_dir,script_paths,args.targets,PATH,out_dir,out_dir,PATH,out_dir,out_dir))
    print '....... Finished\n'

    if args.sge_computing < 1:
        #print "%s/%s/workdir_%s/blast/sh.blast/"%(PATH,out_dir,out_dir)
        subprocces_sh(args.core_num,"%s/%s/workdir_%s/QC/sh.identification/"%(PATH,out_dir,out_dir))
    else:
        print "Combine sh.identification files....."
        print "perl %slist_qsub.pl %s/%s/workdir_%s/QC/sh.dentification/* > %s/%s/Log_%s/sgeQC.sh\n"%(script_paths,PATH,out_dir,out_dir,PATH,out_dir,out_dir)
        os.system("perl %slist_qsub.pl %s/%s/workdir_%s/QC/sh.identification/* > %s/%s/Log_%s/sgeQC.sh"%(script_paths,PATH,out_dir,out_dir,PATH,out_dir,out_dir))
        print '....... Finished\n'

        print "Run sgeQC.sh files......"
        print "sh %s/%s/Log_%s/sgeQC.sh\n"%(PATH,out_dir,out_dir)
        os.system("sh %s/%s/Log_%s/sgeQC.sh"%(PATH,out_dir,out_dir))
        print '....... Finished\n'


    print "Generate sh.readcounting ......"
    print "perl %sread-counting_wrapper09212017DY.pl %s %s/%s/workdir_%s %s/%s/workdir_%s/QC/out.identification/*\n"%(script_paths,script_paths,PATH,out_dir,out_dir,PATH,out_dir,out_dir)
    os.system("perl %sread-counting_wrapper09212017DY.pl %s %s/%s/workdir_%s %s/%s/workdir_%s/QC/out.identification/*\n"%(script_paths,script_paths,PATH,out_dir,out_dir,PATH,out_dir,out_dir))
    print '....... Finished\n'

    if args.sge_computing < 1:
        #print "%s/%s/workdir_%s/blast/sh.blast/"%(PATH,out_dir,out_dir)
        subprocces_sh(args.core_num,"%s/%s/workdir_%s/QC/sh.count/"%(PATH,out_dir,out_dir))
    else:
        print "Combine sh.identification files....."
        print "perl %slist_qsub.pl %s/%s/workdir_%s/QC/sh.count/* > %s/%s/Log_%s/sgeCount.sh\n"%(script_paths,PATH,out_dir,out_dir,PATH,out_dir,out_dir)
        os.system("perl %slist_qsub.pl %s/%s/workdir_%s/QC/sh.count/* > %s/%s/Log_%s/sgeCount.sh"%(script_paths,PATH,out_dir,out_dir,PATH,out_dir,out_dir))
        print '....... Finished\n'

        print "Run sgeCount.sh files......"
        print "sh %s/%s/Log_%s/sgeCount.sh\n"%(PATH,out_dir,out_dir)
        os.system("sh %s/%s/Log_%s/sgeCount.sh"%(PATH,out_dir,out_dir))
        print '....... Finished\n'


    print "Merging out.count files......"
    print "perl %smerge_data.pl  %s/%s/workdir_%s/QC/out.count/* > %s/%s/Log_%s/merged_count.dmp\n"%(script_paths,PATH,out_dir,out_dir,PATH,out_dir,out_dir)
    os.system("perl %smerge_data.pl  %s/%s/workdir_%s/QC/out.count/* > %s/%s/Log_%s/merged_count.dmp\n"%(script_paths,PATH,out_dir,out_dir,PATH,out_dir,out_dir))
    print '....... Finished\n'

    with open("%s/%s/Log_%s/merged_count.dmp"%(PATH,out_dir,out_dir), "rt") as fin:
        with open("%s/%s/Log_%s/merged_count_forPy_temp.txt"%(PATH,out_dir,out_dir), "wt") as fout:
            c=0
            for line in fin:
                c+=1
                if (c == 1):
                    fout.write("{")
                else:
                    fout.write(line.replace('=>', ':'))
    fin.close()
    fout.close()

    with open("%s/%s/Log_%s/merged_count_forPy_temp.txt"%(PATH,out_dir,out_dir), "rt") as fin:
        with open("%s/%s/Log_%s/merged_count_forPy.txt"%(PATH,out_dir,out_dir), "wt") as fout:
            c=0
            for line in fin:
                fout.write(line.replace(';', ''))
    fin.close()
    fout.close()
    os.system("rm %s/%s/Log_%s/merged_count_forPy_temp.txt"%(PATH,out_dir,out_dir))
    #os.system("rm Log_%s/merged_count.dmp"%(out_dir))


    print "Calling KO ......"
    print "python %sCall_KO.py  %s/%s/Log_%s/merged_count_forPy.txt %s Plus %f %s/%s/Log_%s %s"%(script_paths,PATH,out_dir,out_dir,args.targets,args.ratio,PATH,out_dir,out_dir,script_paths)
    os.system("python %sCall_KO.py  %s/%s/Log_%s/merged_count_forPy.txt %s Plus %f %s/%s/Log_%s %s"%(script_paths,PATH,out_dir,out_dir,args.targets,args.ratio,PATH,out_dir,out_dir,script_paths))



    print 'Finished running RCP-PCR program suite.'




def get_sh(x):
    files = []
    lis = os.listdir("%s" %x)
    for i in lis:
        if i[-3:] == ".sh" :
            files.append(i)
    return files



def subprocces_sh(core,sh_dir):
    command_L = get_sh(sh_dir)
    command_L = ["sh %s%s"%(sh_dir,i) for i in command_L]
    if core > len(command_L):
        core = len(command_L)
    n = len(command_L)/core
    pool = mp.Pool(core)
    feed = {}
    for i in range(core):
        feed.update({i:[]})
    lp = len(command_L)
    ky = range(0,lp)
    v = 0
    for i in xrange(lp):
        feed[v].append(command_L[0])
        del command_L[0]
        v +=1
        if (v==core):
            v = 0
    fed = []
    for i in range(0,core):
        fed.append(feed[i])
    print "Dividing and runnig jobs.....",
    results = pool.map(run_sh,fed)
    #print "Done"
    pass

def run_sh(sh_L):
    for sh in sh_L:
        print "Running %s......."%(sh)
        os.system(sh)
        #print "Done"


def LL2fna(LL,name):
    with open(name,"w") as F:
        for L in LL:
            F.write(">%s\n" %  ("\n").join( [str(i) for i in L]))
    F.close()


def rev_comp(seq):
    complement_dict = {'A':'T','T':'A','G':'C','C':'G'}
    return "".join([complement_dict[base] for base in reversed(seq)])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='rcppcr_ko')
    parser.add_argument('-R1','--input_file_R1', default=False,help='Input file of R1.fastq')
    parser.add_argument('-R2','--input_file_R2', default=False,help='Input file of R2.fastq')
    parser.add_argument('-t','--targets', default=False,help='Input target informtion in csv format. (see wiki for detail)')
    parser.add_argument('-out','--output_name', default="output" )
    parser.add_argument('-r','--ratio', type=int, help='Minimum threashold (0 < ratio < 0.5 ) to call mutation profile',default=0.1)
    parser.add_argument('-c','--core_num', type=int, help='Number of cores for multi-processing on local computer.',default=1)
    parser.add_argument('-sge','--sge_computing',default=0, type=int, help='1 if computing on SGE computers.')
    args = parser.parse_args()
    #print args
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
