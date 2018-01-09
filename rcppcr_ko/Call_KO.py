#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
Python script to analyse Clonal KO RCP-PCR result.
"""
import sys
import re
import os
import pickle
import ast
import json

def main(dat,info,strand,abundance_threashold,di,p):
    big_d = reading(dat)
    info_d= csv2dict(info)
    #pprint.pprint(info_d)
    #pprint.pprint(big_d)
    heatmap_d, bar_d= get_top_profiles(big_d,info_d,di,abundance_threashold)
    #pprint.pprint(LL_for_R)
    for target in heatmap_d:
        if not os.path.isdir('%s/csv/%s'%(di,target)):
            os.makedirs('%s/csv/%s'%(di,target))
        Fname_h = "%s/csv/%s/heatmap.csv"%(di,target)
        Fname_b = "%s/csv/%s/barplot.csv"%(di,target)
        LL2csv(heatmap_d[target],Fname_h)
        LL2csv(bar_d[target],Fname_b)
        PDF_name = ("%s/pdf/%s_plot.pdf"%(di,target))
        os.system("Rscript %srcppcr_ko_heatmap.r %s %s %s"%(p,Fname_b,Fname_h,PDF_name))


def get_top_profiles(dat,info,di,abundance_threashold):
    out_dat = {}
    out_bar = {}
    LL_for_csv = [["Target","Plate","Row","Col","Total_reads","#Mutation_profiles_above_threashold","Well_KO_stat","Reads_support_stat(%)", "Profile#1","#Reads#1_rate","Standard_Error#1", "Profile#2","#Reads#2_rate","Standard_Error#2"  ]  ]

    for target_site in dat:
        #print target_site,
        out_dat[target_site] = []
        out_bar[target_site] = []

        out_dat[target_site].append( ["SampleID","base_index","Char_stat","ORDER"])
        out_bar[target_site].append( ["SampleID","Ratio","SE","ORDER"])


        #print tot_reads
        count = 0
        COUNTS = 0
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
            out_dat[target_site].append( ["Reference",char +1, char_stat,COUNTS])

        #print info[target_site]["Target_seq"][gRNA_s-4 : gRNA_s-2], info[target_site]["Target_seq"][gRNA_e+2 : gRNA_e+4]
        if (info[target_site]["Target_seq"][gRNA_s-3 -1 : gRNA_s-2 ] == "CC"):
            out_dat[target_site][gRNA_s-3] = ["Reference", gRNA_s-3, "PAM",COUNTS]
            out_dat[target_site][gRNA_s-2] = ["Reference", gRNA_s-2, "PAM",COUNTS]
            out_dat[target_site][gRNA_s-1] = ["Reference", gRNA_s-1, "PAM",COUNTS]

        elif (info[target_site]["Target_seq"][gRNA_e+2-1 : gRNA_e+3] == "GG"):
            out_dat[target_site][gRNA_e+1] = ["Reference", gRNA_e+1, "PAM",COUNTS]
            out_dat[target_site][gRNA_e+2] = ["Reference", gRNA_e+2, "PAM",COUNTS]
            out_dat[target_site][gRNA_e+3] = ["Reference", gRNA_e+3, "PAM",COUNTS]
        else:
            print "Error: NO PAM for %s" %target_site
            #print info[target_site]["Target_seq"][gRNA_s-3-1: gRNA_e+3-1]
            #print info[target_site]["Target_seq"][gRNA_s-3-1: gRNA_s-2-1]
            #print info[target_site]["Target_seq"][gRNA_e+2-1: gRNA_e+3-1]
        out_dat[target_site].append( ["Reference (%s)"%(target_site),0, "Ref", COUNTS])
        out_bar[target_site].append( ["Reference",0, 0, COUNTS]  )
        COUNTS -=1

        for plate in dat[target_site]:
            #print plate,
            for row in dat[target_site][plate]:
                #print row
                for col in dat[target_site][plate][row]:
                    L = []

                    Well_KO_stat = "-"
                    target_loci = dat[target_site][plate][row][col][strand]
                    tot_reads = sum(dat[target_site][plate][row][col][strand].values())

                    L += [target_site,plate,row,col,tot_reads]


                    prof = 0
                    sum_of_hits = 0
                    profiles       = []
                    #Frameshift = "- (WT)"
                    for mutation_profile in dat[target_site][plate][row][col][strand]:
                        Frameshift = "- (WT)"
                        hits = dat[target_site][plate][row][col][strand][mutation_profile]
                        if hits > 1:
                            hit_rat , error = ErrorProp_div(hits,tot_reads)
                            imputed_hitrate = hit_rat-error


                            #print mutation_profile, imputed_hitrate
                            if (imputed_hitrate >= abundance_threashold) :
                                prof +=1
                                count +=1
                                #print target_site,plate
                                #print target_site,plate,row,col,tot_reads,hits,imputed_hitrate, mutation_profile
                                sum_of_hits += hits

                                if prof == 1:
                                    out_dat[target_site].append( ["%s_%s_%s______________________________________________________"%(plate,row,col), 0, "Ref",COUNTS ])
                                    out_bar[target_site].append( ["%s_%s_%s______________________________________________________"%(plate,row,col), 0, 0 ,COUNTS])
                                    COUNTS -=1

                                mut_l = btop_deconvolute(info[target_site]["Target_seq"],mutation_profile)
                                for char in range(len(mut_l)):
                                    out_dat[target_site].append( ["%s_%s_%s_Profile#%d"%(plate,row,col,prof), char +1,  mut_l[char],COUNTS ])
                                if mut_l.count("Del") > 0:
                                    if (mut_l.count("Del")%3 != 0):
                                        Frameshift = "'++ (Frameshift)"
                                    else:
                                        Frameshift = "'+ (Indel)"
                                if mut_l.count("Ins") >0:
                                    if (mut_l.count("Ins")%3 != 0):
                                        Frameshift = "'++ (Frameshift)"
                                    else:
                                        Frameshift = "'+ (Indel)"
                                profiles += [Frameshift, round(hit_rat,2)*100,round(error,2)*100 ]


                                out_dat[target_site].append( ["%s_%s_%s_Profile#%d"%(plate,row,col,prof), 0,  Frameshift,COUNTS ])
                                out_bar[target_site].append( ["%s_%s_%s_Profile#%d"%(plate,row,col,prof), hit_rat,  error,COUNTS ])
                                COUNTS -=1


                    wt = profiles.count("'- (WT)")
                    frameshift = profiles.count("'++ (Frameshift)")
                    indel = profiles.count("'+ (Indel)")

                    total_profs = wt+frameshift+ indel

                    if total_profs == 2:
                        if frameshift ==2:
                            Well_KO_stat = "'+++ (Homo-frameshift)"
                        if indel ==2:
                            Well_KO_stat = "'++ (Homo-indel)"
                        if wt ==1:
                            Well_KO_stat = "'+ (Hetero)"
                    if total_profs==1:
                        if wt ==1:
                            Well_KO_stat = "'T"
                        if frameshift ==1:
                            Well_KO_stat = "'+++? (Frameshift; 1 profile)"
                        if indel ==1:
                            Well_KO_stat = "'++? (Indel; 1 profile)"
                    if total_profs > 2:
                        Well_KO_stat = "'?? (Too many profiles)"



                    L += [ total_profs,Well_KO_stat,float(sum_of_hits)/tot_reads*100 ]
                    L += profiles[:6]
                    LL_for_csv.append(L)


    LL2csv(LL_for_csv,"%s/%s_sumary.csv"%(di,di.split("Log_")[1]))
    return out_dat, out_bar

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

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r



if __name__ == '__main__':
    dat  = sys.argv[1]
    info = sys.argv[2]
    strand = sys.argv[3] #"Plus" or "Minus"
    abundance_threashold = float(sys.argv[4])#Default = 0.1
    di = sys.argv[5]#Directory name. Defalt "."
    p = sys.argv[6]
    main(dat,info,strand,abundance_threashold,di,p)
