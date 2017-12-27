#!/usr/bin/env python
"""
Python script to analyse Clonal KO RCP-PCR result.
"""
import sys
import re
import os
import pickle
import ast
from tqdm import tqdm
import json
import pprint

def main(dat,info,strand,abundance_threashold,di):
    big_d = reading(dat)
    info_d= csv2dict(info)
    #pprint.pprint(info_d)
    #pprint.pprint(big_d)
    LL_for_R = get_top_profiles(big_d,info_d,di)

    for target in LL_for_R:
        for plate in LL_for_R[target]:
            Fname = "%s/csv/%s_%s_for_Rplot.csv"%(di,target,plate)
            LL2csv(LL_for_R[target][plate],Fname)

            PDF_name = ("%s/pdf/%s_%s_plot.pdf"%(di,target,plate))

            os.system("Rscript codes/KO_heatmap.r %s %s"%(Fname,PDF_name))
            


def get_top_profiles(dat,info,di):
    out_dat = {}
    LL_for_csv = [["Target","Plate","Row","Col","Total_reads","#Mutation_profiles_above_threashold","KO_stat","KO_detail(OnlyFrameshiftMiutations)","Reads_support_stat(%)","Btop"  ]  ]
    for target_site in dat:
        #print target_site,
        out_dat[target_site] = {}

        for plate in dat[target_site]:
            #print plate,
            out_dat[target_site][plate] = []
            out_dat[target_site][plate].append( ["Row","Col","Profile","base_index","Char_stat"])
            for row in dat[target_site][plate]:
                #print row,
                for col in dat[target_site][plate][row]:
                    target_loci = dat[target_site][plate][row][col][strand]
                    tot_reads = sum(dat[target_site][plate][row][col][strand].values())
                    #print tot_reads
                    count = 0
                    gRNA_s = int(info[target_site]["gRNA_s"])
                    gRNA_e = int(info[target_site]["gRNA_e"])
                    for char in range(len(info[target_site]["Target_seq"])):
                        char_stat = "-"
                        PAM_frd_G = 0
                        PAM_frd_C = 0
                        PAM_rvs_G = 0
                        PAM_rvs_C = 0
                        if (gRNA_s <= (char+1) <=gRNA_e):
                            char_stat = "gRNA"
                        if (gRNA_s -3 <= (char+1) < gRNA_s-1):
                            if (info[target_site]["Target_seq"][char] == "G"):
                                PAM_frd_G +=1
                            if (info[target_site]["Target_seq"][char] == "C"):
                                PAM_frd_C +=1
                        if (gRNA_e  < (char+1) <+ gRNA_e+3):
                            if (info[target_site]["Target_seq"][char] == "G"):
                                PAM_rvs_G +=1
                            if (info[target_site]["Target_seq"][char] == "C"):
                                PAM_rvs_C +=1
                        
                                
                        for char in range(len(info[target_site]["Target_seq"])):
                            if  (PAM_rvs_G ==2):
                                if (gRNA_e< (char+1) <= gRNA_e+3):
                                    char_stat = "PAM"
                            elif (PAM_rvs_C ==2):
                                if (gRNA_e< (char+1) <= gRNA_e+3):
                                    char_stat = "PAM"

                            elif(PAM_frd_G ==2):
                                if (gRNA_s-3 <= (char+1) <= gRNA_e-1):
                                    char_stat = "PAM"
                            elif (PAM_rvs_C ==2):
                                if (gRNA_s-3<= (char+1) <= gRNA_s-1):
                                    char_stat = "PAM"
                            out_dat[target_site][plate].append( [row,col,"WT",char +1, char_stat])


                    sum_of_hits = 0
                    btops       = []
                    KO_stat = []
                    for mutation_profile in dat[target_site][plate][row][col][strand]:
                        hits = dat[target_site][plate][row][col][strand][mutation_profile]
                        hit_rat , error = ErrorProp_div(hits,tot_reads)
                        imputed_hitrate = hit_rat-error



                        #print mutation_profile, imputed_hitrate
                        if (imputed_hitrate >= abundance_threashold) :
                            count -=1
                            #print target_site,plate
                            #print target_site,plate,row,col,tot_reads,hits,imputed_hitrate, mutation_profile
                            btops.append(mutation_profile)
                            sum_of_hits += hits

                            mut_l = btop_deconvolute(info[target_site]["Target_seq"],mutation_profile)

                            for char in range(len(mut_l)):
                                out_dat[target_site][plate].append( [row,col, imputed_hitrate ,char +1,  mut_l[char] ])
                            if mut_l.count("Del") >0:
                                if (mut_l.count("Del")%3 != 0):
                                    KO_stat.append("Del(Frameshift)")
                                else:
                                     KO_stat.append("Del")
                            if mut_l.count("Ins") >0:
                                if (mut_l.count("Ins")%3 != 0):
                                    KO_stat.append("Ins(Frameshift)")
                                else:
                                    KO_stat.append("Ins")
                    if len(KO_stat) >= 2:
                        KO_detail = (":").join(KO_stat)
                        KO_stat = "Homo"
                    if len(KO_stat) ==1:
                        KO_detail = KO_stat[0]
                        KO_stat = "Hetero"
                    if len(KO_stat) ==0:
                        KO_detail = "-"
                        KO_stat = "WT"
                    LL_for_csv.append([target_site,plate,row,col,tot_reads,abs(count),KO_stat,KO_detail,round((float(sum_of_hits)/tot_reads)*100,2),(":").join(btops)  ]   )
                    na_num = 0
                    if count > -2:
                        for remaining in range(count+2  ):
                            na_num +=1
                            for char in range(len(info[target_site]["Target_seq"])):
                                out_dat[target_site][plate].append( [row,col,"N.A.%d"%(na_num) ,char+1,"-" ])


    LL2csv(LL_for_csv,"%s/%s_sumary.csv"%(di,di.split("Log_")[1]))
    return out_dat

################################
#### Small functions
################################

def csv2dict(f_name):
    d = {}
    with open(f_name,"r") as F:
        c =0
        for line in F:
            c+=1
            cols = line.split(",")
            #print cols
            if (c ==1):
                header = cols
            else:
                d[cols[0]] = {}
                for i in range(len(header[2:])+1):
                    d[cols[0]][header[i+1].split("\n")[0]] =cols[i+1].split("\n")[0]

    return d

def LL2csv(LL,name):
    with open(name,"w") as F:
        for L in LL:
            F.write("%s\n" %  (",").join( [str(i) for i in L]))
    F.close()



def ErrorProp_div(Hit,All):
    val = 0
    error = 0

    Hit = float(Hit)
    All = float(All)


    val = Hit/All
    d_All = All **(0.5)
    #print All,d_All
    d_Hit = Hit **(0.5)
    #print target, d_target
    error = ((d_All/All) ** 2) + ((d_Hit/Hit) ** 2)
    error **= 0.5
    error  *= val

    return (val,error);

def reading(name):
    with open(name, 'r') as F:
        ob_from_file = eval(F.read())
    F.close()
    return ob_from_file

def btop_deconvolute(target,btop):
    btop_l = [a for a in  re.split('(\d+)',btop) if (len(a) > 0)]

    mut_mem = [] #for showing % in heatmap
    #print btop_l
    for i in btop_l:
        try:
            mut_mem += ["-" for i in range(0,int(i))]
            #print len(mut_mem),mut_mem
        except ValueError:
            Mut = [i[j:j+2] for j in range(0, len(i), 2)]
            #print Mut
            for m in Mut:
                if m[0] == "-":
                    mut_mem.append("Del")
                elif m[1] == "-":
                    mut_mem[-1] = ("Ins")
                else:
                    mut_mem.append("Mut")
    #print len(target),len(mut_mem)
    return mut_mem




if __name__ == '__main__':
    dat  = sys.argv[1]
    info = sys.argv[2]
    strand = sys.argv[3] #"Plus" or "Minus"
    abundance_threashold = int(sys.argv[4])#Default = 0.1
    di = sys.argv[5]#Default = 0.1
    main(dat,info,strand,abundance_threashold,di)
