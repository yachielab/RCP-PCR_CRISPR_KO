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
import numpy as np
import math
import xlrd
import pprint

def main(dat,info,strand,abundance_threashold,di,cellv,coord):
    big_d = reading(dat)
    info_d= csv2dict(info)
    target_list = extract_targets(info)
    rows,cols = rc_pcr_coordinate(coord)
    wells     = list2combinatoral(rows,cols)
    D         = {}
    plate_tag = ["P01-P01","P04-P04","P05-P05","P06-P06","P07-P07"]

    for target in target_list:
        D[target] = {}
        for plate in plate_tag:
            D[target][plate] = {}
            for well in wells:
                D[target][plate][well] = {}

    D2 = read_viability_files(cellv,coord,D)

    heatmap_d, bar_d,A,LL_for_csv = get_top_profiles(big_d,info_d,di,D2,wells,plate_tag,strand)



    #LL2csv(LL_for_csv,"%s/%s_sumary_10022018DEY.csv"%(di,di.split("Log_")[1]))
    LL2csv(LL_for_csv,"sumary_10022018DEY.csv")
    #LL2csv(LL_profile,"%s/%s_prenormed.csv"%(di,di.split("Log_")[1]))
    #LL2csv(LL_normed,"%s/%s_normed.csv"%(di,di.split("Log_")[1]))
    #LL2csv(LL_profile,"%s/%s_profile.csv"%(di,di.split("Log_")[1]))

    #pprint.pprint(LL_for_R)

    """for target in heatmap_d:
        if not os.path.isdir('%s/csv/%s'%(di,target)):
            os.makedirs('%s/csv/%s'%(di,target))
        Fname_h = "%s/csv/%s/heatmap.csv"%(di,target)
        Fname_b = "%s/csv/%s/barplot.csv"%(di,target)
        LL2csv(heatmap_d[target],Fname_h)
        LL2csv(bar_d[target],Fname_b)
        PDF_name = ("%s/pdf/%s_plot.pdf"%(di,target))
        os.system("Rscript %srcppcr_ko_heatmap.r %s %s %s"%(p,Fname_b,Fname_h,PDF_name))"""

def get_top_profiles(dat,info,di,d,wells,plate_tag,strand):
    A_RECORD = []
    out_dat = {}
    out_bar = {}
    LL_for_csv   = [["Target","Well","Cell_viability","Well_stat","#Mutation_profiles_above_threashold","KO_profile","Confidence_level","Frameshift_allels"]]
    LL_for_csv[0] += [ "#Total_reads_%s"%(i) for i in plate_tag]
    LL_for_csv[0] += [ "#Total_ratio_%s"%(i) for i in plate_tag]
    LL_for_csv[0]  +=["Plate_replicates_support_stat","#Plate_replicates_support_stat","Mean_Reads_support_stat(%)","Std_Reads_support_stat(%)", "Profile#1","Profile#2","Profile#3","Profile#4"  ]




    normed_dat,LL_for_error,Detail_LL_error,Detail_LL_normed = normalize_reads(dat,d,wells)
    #pprint.pprint(normed_dat)
    #print normed_dat
    dat = normed_dat
    SUM = {}
    F = {}
    for target_site in normed_dat:
        F[target_site] ={}
        #print target_site,
        out_dat[target_site] = []
        out_bar[target_site] = []

        out_dat[target_site].append( ["SampleID","base_index","Char_stat","ORDER"])
        out_bar[target_site].append( ["SampleID","Ratio","SE","ORDER"])

        SUM[target_site] = {}
        #print tot_reads
        count = 0
        COUNTS = 0
        gRNA_s = int(info[target_site]["gRNA_s"])
        gRNA_e = int(info[target_site]["gRNA_e"])

        temp_dat = []
        temp_bar = []
        for char in range(len(info[target_site]["Target_seq"])):
            char_stat = "-"
            PAM_frd_G = 0
            PAM_frd_C = 0
            PAM_rvs_G = 0
            PAM_rvs_C = 0
            if (gRNA_s <= (char+1) <=gRNA_e):
                char_stat = "gRNA"
            temp_dat.append( ["Reference (%s)"%(target_site),char, char_stat,COUNTS])
        p = 0
        #print target_site, info[target_site]["Target_seq"][gRNA_s-4 : gRNA_e+4]
        if (info[target_site]["Target_seq"][gRNA_s-3 -1 : gRNA_s-2 ] == "CC"):
            temp_dat[gRNA_s-3] = ["Reference (%s)"%(target_site), gRNA_s-3, "PAM",COUNTS]
            temp_dat[gRNA_s-2] = ["Reference (%s)"%(target_site), gRNA_s-2, "PAM",COUNTS]
            temp_dat[gRNA_s-1] = ["Reference (%s)"%(target_site), gRNA_s-1, "PAM",COUNTS]
            p +=1

        elif (info[target_site]["Target_seq"][gRNA_e+2-1 : gRNA_e+3] == "GG"):
            temp_dat[gRNA_e+1] = ["Reference (%s)"%(target_site), gRNA_e+1, "PAM",COUNTS]
            temp_dat[gRNA_e+2] = ["Reference (%s)"%(target_site), gRNA_e+2, "PAM",COUNTS]
            temp_dat[gRNA_e+3] = ["Reference (%s)"%(target_site), gRNA_e+3, "PAM",COUNTS]
        else:
            print "Error: NO PAM for %s" %target_site
            #print info[target_site]["Target_seq"][gRNA_s-3-1: gRNA_e+3-1]
            #print info[target_site]["Target_seq"][gRNA_s-3-1: gRNA_s-2-1]
            #print info[target_site]["Target_seq"][gRNA_e+2-1: gRNA_e+3-1]
        out_dat[target_site].append( ["Reference (%s)"%(target_site), -30, "Ref", COUNTS])
        out_bar[target_site].append( ["Reference (%s)"%(target_site),0, 0, COUNTS]  )
        COUNTS -=1

        if p ==1:
            MOVE = gRNA_s
            F[target_site] = gRNA_s
        else:
            MOVE = gRNA_e
            F[target_site] = gRNA_e

        for l in temp_dat:
            if abs(l[1]-MOVE) <= 30:
                 out_dat[target_site].append( [l[0],l[1]-MOVE,l[2],l[3]] )
        for l in temp_bar:
            if abs(l[1]-MOVE) <= 30:
                out_bar[target_site].append([l[0],l[1],l[2]])
        PROFILES = {}
        PROF_ID = 1
        HIT_RATE_MEM = {}
        HIT_MEM = {}
        for pos in wells:
            HIT_RATE_MEM[pos] = {}
            HIT_MEM[pos] = {}

            Well_KO_stat = [] #Add [%support_stat,set(profiles)]
            well = TR96_convert(pos.split("-")[0],pos.split("-")[1])
            #print well
            level = "Level3"
            L = []
            L += [target_site,well]
            PS = {}
            for plate in plate_tag:
                HIT_MEM[pos][plate] = {}
                HIT_RATE_MEM[pos][plate] = {}
                prof = 0
                sum_of_hits = 0
                profiles       = []
                try:
                    tot_reads = sum(normed_dat[target_site][plate][well][strand].values())
                    #Frameshift = "- (WT)"
                    for mutation_profile in [ i[0] for i in sorted(normed_dat[target_site][plate][well][strand].iteritems(), key=lambda (k,v):(v,k),reverse=True)]:
                        hits = normed_dat[target_site][plate][well][strand][mutation_profile]
                        if hits > 1:
                            hit_rat , error = ErrorProp_div(hits,tot_reads)
                            imputed_hitrate = hit_rat-error
                            #print mutation_profile, imputed_hitrate
                            if (imputed_hitrate >= abundance_threashold) :
                                prof +=1
                                count +=1
                                profiles.append(mutation_profile)
                                HIT_MEM[pos][plate][mutation_profile] = hit_rat
                                try:
                                    DUMMY = PROFILES[mutation_profile] -1
                                except KeyError:
                                    PROFILES[mutation_profile] =  PROF_ID
                                    PROF_ID +=1
                                #print target_site,plate
                                #print target_site,plate,row,col,tot_reads,hits,imputed_hitrate, mutation_profile
                                sum_of_hits += hits
                    for i in HIT_MEM[pos][plate]:
                        HIT_RATE_MEM[pos][plate][i] = round(float(HIT_MEM[pos][plate][i])/float(tot_reads),3)*100
                    Well_KO_stat.append([round(float(sum_of_hits)/float(tot_reads),4)*100,set(profiles),plate])

                    PS[plate] = set(profiles)
                except KeyError:
                    pass



            replicates = {}
            plates_used = {}
            for replicate in Well_KO_stat:
                try:
                    replicates[frozenset(replicate[1])].append(replicate[0])
                    plates_used[frozenset(replicate[1])].append(replicate[2])
                except KeyError:
                    replicates[frozenset(replicate[1])] = [replicate[0]]
                    plates_used[frozenset(replicate[1])] = [replicate[2]]
            #print target_site,pos,replicates
            #print plates_used
            #print replicates
            if len([ i for i in sorted(replicates.iteritems(), key=lambda kv: (len(kv[1]), kv[0]),reverse=True)]) > 0:
                most_confident_profiles = [ i for i in sorted(replicates.iteritems(), key=lambda kv: (len(kv[1]), kv[0]),reverse=True)][0][0]

                if len( replicates[most_confident_profiles]) > 1:
                    c85 = 0
                    c75 = 0
                    for i in replicates[most_confident_profiles]:
                        if i >= 0.85:
                            c85 +=1
                        if i >= 0.75:
                            c75 +=1
                    if c85 >= 2:
                        if len(most_confident_profiles) ==1:
                            level = "Level2"
                        elif len(most_confident_profiles) > 1:
                            level = "Level1"
                    elif c75 >= 2:
                        level = "Level2"
                NUM_frame = 0
                genotype = []
                G = []
                for i in most_confident_profiles:
                    #print i
                    G.append(i)
                    Frameshift = "WT"
                    mut_l = btop_deconvolute(info[target_site]["Target_seq"],i)
                    if mut_l.count("Del") > 0:
                        if (mut_l.count("Del")%3 != 0):
                            Frameshift = "Frameshift"
                            NUM_frame +=1
                        else:
                            Frameshift = "Indel"
                    elif mut_l.count("Ins") >0:
                        if (mut_l.count("Ins")%3 != 0):
                            Frameshift = "Frameshift"
                            NUM_frame +=1
                        else:
                            Frameshift = "Indel"
                    genotype.append(Frameshift)
                wt = genotype.count("WT")
                frameshift = genotype.count("Frameshift")
                indel = genotype.count("Indel")
                KO_stat= "WT"
                total_profs = wt+frameshift+ indel
                if total_profs > 0:
                    if total_profs == 2:
                        if frameshift ==2:
                            KO_stat = "Homo-Frameshift"
                        elif indel ==2:
                            KO_stat = "Homo-Indel"
                        elif wt ==1:
                            if indel ==1:
                                KO_stat = "Hetero-WT/Indel"
                            if frameshift ==1:
                                KO_stat = "Hetero-WT/Frameshift"
                        elif frameshift + indel ==2:
                            KO_stat = "Homo-Frameshift/Indel"
                    if total_profs==1:
                        if wt ==1:
                            KO_stat = "WT"
                        elif frameshift ==1:
                            KO_stat = "Homo-Frameshift;1profile"
                        elif indel ==1:
                            KO_stat = "Homo-Indel;1profile"
                    if total_profs > 2:
                        KO_stat = "Excess_profiles"
                        if frameshift == total_profs:
                            KO_stat = "Homo-Frameshift;>2profiles"

                final_stat = ("_").join([level,KO_stat])

                L += [NUM_frame,d[target_site][plate][pos]["Cell_viability"],final_stat,str(total_profs),KO_stat,level]
                val = 0
                plttg = ["P01-P01","P04-P04","P05-P05","P06-P06","P07-P07"]
                for plateee in plttg:
                    for STAAAT in Well_KO_stat:
                        if (plateee == STAAAT[2]):
                            val = STAAAT[0]
                    if val==0:
                        L += [ 0 ]
                    else:
                        L += [val]

                #L += [frameshift,len() (";").join(sorted(plates_used[most_confident_profiles])),len(plates_used[most_confident_profiles]),str(np.mean(replicates[most_confident_profiles])),str(np.std(replicates[most_confident_profiles]))]
                #L += G

                #print most_confident_profiles,replicates[most_confident_profiles]
                LL_for_csv.append(L)

                out_dat[target_site].append( ["%s_______________________________________________________"%(pos), -30, "Ref",COUNTS ])
                out_bar[target_site].append( ["%s_______________________________________________________"%(pos), 0, 0 ,COUNTS])
                COUNTS -=1

                PL = []
                for i in plate_tag:
                    try:
                        if (PS[i] == most_confident_profiles):
                            PL.append(i)
                    except KeyError:
                        pass


                for profile in  most_confident_profiles:
                    rat = np.mean([ HIT_RATE_MEM[pos][a][profile] for a in PL])
                    error = np.std([ HIT_RATE_MEM[pos][a][profile] for a in PL])
                    out_bar[target_site].append( ["%s_Profile#%d"%(pos,PROFILES[profile]),rat ,  error,COUNTS ])

                    mut_l = btop_deconvolute(info[target_site]["Target_seq"],profile)
                    for char in range(len(mut_l)):
                        if abs(char-MOVE) <= 30:
                            out_dat[target_site].append( ["%s_Profile#%d"%(pos, PROFILES[profile]), char+1+MOVE,mut_l[char],COUNTS ])


                    COUNTS -=1



        #for profile in most_confident_profiles:
        #L += most_confident_profiles


    return out_dat, out_bar,A_RECORD,LL_for_csv

################################
#### Small functions
################################

def normalize_reads(D,d,wells):
    #total_reads[target_site][profile] = count
    total_reads      =  {}
    LL_for_error     =  [["Target","Plate","Profile","#Total_reads","#Total_wells","#Mean_reads","#Wells(w/cells)","#Mean_reads(w/cells)","#Wells(w/o_cells)","#Mean_reads(w/o_cells)"]]
    Detail_LL_error  =  [["Target","Plate","Well","Row","Col","Cellv","Profile","#Reads"]]
    Detail_LL_normed =  [["Target","Plate","Well","Row","Col","Cellv","Profile","#Reads"]]
    TOT_reads_per_well = [["Target","Plate","Well","Cellv","Reads"]]
    plate_tag       =  ["P01-P01","P04-P04","P05-P05","P06-P06","P07-P07"]
    #normed[target_site][plate][row][col]["Plus"][profile]= normedcount
    normed = {}
    D_for_Error = {}

    for target_site in D:
        D_for_Error[target_site] = {}
        total_reads[target_site] = {}
        for plate in plate_tag:
            D_for_Error[target_site][plate] = {}
            total_reads[target_site][plate] = {}
            for row in D[target_site][plate]:
                for col in D[target_site][plate][row]:
                    if "%s-%s"%(row,col) in wells:
                        #print target_site,plate,row,col, D[target_site][plate][row][col].keys()
                        for profile in D[target_site][plate][row][col]["Plus"]:
                            total_reads[target_site][plate][profile] = {}
                            total_reads[target_site][plate][profile][0] = []
                            total_reads[target_site][plate][profile][1] = []
                            total_reads[target_site][plate][profile]["All"] = []


            for row in D[target_site][plate]:
                for col in D[target_site][plate][row]:
                    if "%s-%s"%(row,col) in wells:
                        c = 0
                        #print target_site,plate,row,col, D[target_site][plate][row][col].keys()
                        for profile in D[target_site][plate][row][col]["Plus"]:
                            count   = D[target_site][plate][row][col]["Plus"][profile]
                            Detail_LL_error.append([target_site,plate,"%s-%s"%(row,col),row,col,d[target_site][plate]["%s-%s"%(row,col)]["Cell_viability"],profile,count])
                            #print target_site,plate,row,col,"%s-%s"%(row,col),d[target_site][plate]["%s-%s"%(row,col)]["Cell_viability"]
                            total_reads[target_site][plate][profile][int(d[target_site][plate]["%s-%s"%(row,col)]["Cell_viability"])].append(float(count))
                            total_reads[target_site][plate][profile]["All"].append(float(count))
                            c += count
                        TOT_reads_per_well.append([target_site,plate,"%s-%s"%(row,col),d[target_site][plate]["%s-%s"%(row,col)]["Cell_viability"],c])
            for row in D[target_site][plate]:
                for col in D[target_site][plate][row]:
                    if "%s-%s"%(row,col) in wells:
                        #print target_site,plate,row,col, D[target_site][plate][row][col].keys()
                        for profile in D[target_site][plate][row][col]["Plus"]:
                            for i in [0,1]:
                                if len(total_reads[target_site][plate][profile][i]) == 0:
                                    total_reads[target_site][plate][profile][i] = [0]


    for target_site in D:
        normed[target_site] = {}
        for plate in D[target_site]:
            if plate in plate_tag:
                normed[target_site][plate] = {}
                #print target_site,plate,len([i for i in [plate] if i in D[target_site].keys()]),normed[target_site].keys()
                for row in D[target_site][plate]:
                    for col in D[target_site][plate][row]:
                        if "%s-%s"%(row,col) in wells:
                            normed[target_site][plate][TR96_convert(row,col)] = {}
                            normed[target_site][plate][TR96_convert(row,col)]["Plus"] = {}
                            for profile in D[target_site][plate][row][col]["Plus"]:
                                #print total_reads[target_site][plate][profile]
                                #print target_site,profile, (np.mean(total_reads[target_site][plate][profile][1])),(np.mean(total_reads[target_site][plate][profile][0]))
                                R = D[target_site][plate][row][col]["Plus"][profile] - (np.mean(total_reads[target_site][plate][profile][0]))
                                if (R > ((np.mean(total_reads[target_site][plate][profile][1]) - np.mean(total_reads[target_site][plate][profile][0]) )*0.1)):
                                    if R > 1:
                                        normed[target_site][plate][TR96_convert(row,col)]["Plus"][profile] = R
                                        Detail_LL_normed.append([target_site,plate,"%s-%s"%(row,col),row,col,d[target_site][plate]["%s-%s"%(row,col)]["Cell_viability"],profile,R])
                                else:
                                    normed[target_site][plate][TR96_convert(row,col)]["Plus"][profile] = 0
                                    normed[target_site][plate][TR96_convert(row,col)]["Plus"] = removekey(normed[target_site][plate][TR96_convert(row,col)]["Plus"],profile)
                                    Detail_LL_normed.append([target_site,plate,"%s-%s"%(row,col),row,col,d[target_site][plate]["%s-%s"%(row,col)]["Cell_viability"],profile,0])
                            if len(normed[target_site][plate][TR96_convert(row,col)]["Plus"]) ==0:
                                normed[target_site][plate][TR96_convert(row,col)] = removekey(normed[target_site][plate][TR96_convert(row,col)],"Plus")
                        if "%s-%s"%(row,col) in wells:
                            if len(normed[target_site][plate][TR96_convert(row,col)]) ==0:
                                normed[target_site][plate] = removekey(normed[target_site][plate],TR96_convert(row,col))
                if len(normed[target_site][plate]) ==0:
                    normed[target_site] = removekey(normed[target_site],plate)
        if len(normed[target_site]) ==0:
            normed = removekey(D,target_site)
    LL2csv(Detail_LL_error,"MiSeq01-20-2018_detail.csv")
    LL2csv(TOT_reads_per_well,"MiSeq01-20-2018_read_num.csv")

    #print dkdjk
    return normed,LL_for_error,Detail_LL_error,Detail_LL_normed


def fx(x):
    return x



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
                    mut_mem.append("Ins")
                else:
                    mut_mem.append("Mut")
    #print len(target),len(mut_mem)
    return mut_mem

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

def list2combinatoral(l1,l2):
    l =[]
    for i1 in l1:
         for i2 in l2:
             i = ("-").join([i1,i2])
             l.append(i)
    return l
def rc_pcr_coordinate(POS):
    if POS[2:] == "96":
        if POS[0] == "T":
            ROW = ["R01","R03","R05","R07","R09","R11","R13","R15"]
        if POS[0] == "B":
            ROW = ["R02","R04","R06","R08","R10","R12","R14","R16"]
        if POS[1] == "L":
            COL = ["C01","C03","C05","C07","C09","C11","C13","C15","C17","C19","C21","C23"]
        if POS[1] == "R":
            COL = ["C02","C04","C06","C08","C10","C12","C14","C16","C18","C20","C22","C24"]
    return ROW, COL
def summary2dict(summary,D):
    c=0
    with open(summary,"r") as F:
        for line in F:
            c +=1
            cols = line.split("\n")[0].split(",")
            #print cols
            if c ==1:
                keys = cols[4:]
            else:
                try:
                    for key in range(len(keys)):
                        try:
                            D[cols[0]][cols[1]][("-").join([cols[2],cols[3]])][keys[key]] = cols[key+4]
                        except IndexError:
                            D[cols[0]][cols[1]][("-").join([cols[2],cols[3]])][keys[key]] = ""

                    #D[cols[0]][cols[1]][("-").join([cols[2],cols[3]])]["Cell_viability"]="N.A."
                except KeyError:
                    pass
        F.close()
    return D

def extract_targets(f):
    targets = []
    with open(f,"r") as F:
        c = 1
        for line in F:
            if c > 1:
                cols = line.split(",")
                targets.append(cols[0])
            c+=1
        F.close()
    return targets


def read_viability_files(directory,corrd,D):

    l = get_xlsx(directory)
    for target in D:
        for plate in D[target]:
            xlsx = [ i for i in l if target in i]
            #print xlsx,target,plate
            if len(xlsx) ==0:
                try:
                    D = removekey(D,target)
                except KeyError:
                    pass
            elif len(xlsx) ==1:
                book = xlrd.open_workbook("%s/%s"%(directory,xlsx[0]))
                sh = book.sheet_by_index(0)
                for well in D[target][plate]:
                    r,c,pos = xlsx_pos_TR96_convert(well)
                    D[target][plate][well]["Cell_viability"] = sh.cell_value(rowx=r-1, colx=c-1)
            elif len(xlsx) > 1:
                print "Too many .xlsx files for Target: %s || %s"%(target,str(xlsx))
    return D


def get_xlsx(x):
    files = []
    lis = os.listdir("%s" %x)
    for i in lis:
        if i[-5:] == ".xlsx" :
            if i[0] != "~":
                files.append(i)
    return files

def xlsx_pos_TR96_convert(well):
    r = well.split("-")[0]
    c = well.split("-")[1]

    R = {"R01":4,"R03":5,"R05":6,"R07":7,"R09":8,"R11":9,"R13":10,"R15":11}
    C = {"C02":2,"C04":3,"C06":4,"C08":5,"C10":6,"C12":7,"C14":8,"C16":9,"C18":10,"C20":11,"C22":12,"C24":13}

    posR = {"R01":"A","R03":"B","R05":"C","R07":"D","R09":"E","R11":"F","R13":"G","R15":"H"}
    posC = {"C02":1,"C04":2,"C06":3,"C08":4,"C10":5,"C12":6,"C14":7,"C16":8,"C18":9,"C20":10,"C22":11,"C24":12}
    pos  = ("").join([posR[r],str(posC[c])])
    return R[r],C[c],pos
def TR96_convert(PosR,PosC):
    posR = {"R01":"A","R03":"B","R05":"C","R07":"D","R09":"E","R11":"F","R13":"G","R15":"H"}
    posC = {"C02":1,"C04":2,"C06":3,"C08":4,"C10":5,"C12":6,"C14":7,"C16":8,"C18":9,"C20":10,"C22":11,"C24":12}
    pos  = ("").join([posR[PosR],str(posC[PosC])])
    return pos





if __name__ == '__main__':
    dat  = sys.argv[1]
    info = sys.argv[2]
    strand = sys.argv[3] #"Plus" or "Minus"
    abundance_threashold = float(sys.argv[4])#Default = 0.1
    di = sys.argv[5]#Directory name. Defalt "."
    #p = sys.argv[6]
    cellv = sys.argv[6]
    coord = sys.argv[7]

    main(dat,info,strand,abundance_threashold,di,cellv,coord)
    # python Call_KO_merge_replicates.py
