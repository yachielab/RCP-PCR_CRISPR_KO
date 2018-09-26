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
import operator
big_d = reading(dat)
    info_d= csv2dict(info)
    #pprint.pprint(info_d)
    #pprint.pprint(big_d)
    a = get_top_profiles(big_d,info_d,di,abundance_threashold)
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


def get_top_profiles(dat,info,di,abundance_threashold):

    sample_index = ["P01-P01","P02-P02","P03-P03","P04-P04","P18-P18","P19-P19","P21-P21","P22-P22"]


    for target_site in dat:
        for plate in dat[target_site]:
            #print plate
            if plate in sample_index:
                if not os.path.isdir('%s/csv/%s/%s'%(di,target_site,plate)):
                    os.makedirs('%s/csv/%s/%s'%(di,target_site,plate))
                mutation_profiles      = {}#btop of profile : #of wells
                mutation_profiles_mem  = {}

                pos_profile_mem = {}


                mutation_profiles_plot = [["Rank","Position","Char_stat"]]
                mutation_pos_plot = [["Rank","Position","Frequency","Well"]]
                Rank = 0 # 0 is always for reference (PAM, sgRNA)

                strand = info[target_site]["Strand"]
                gRNA_s = int(info[target_site]["gRNA_s"])
                gRNA_e = int(info[target_site]["gRNA_e"])
                target_seq = info[target_site]["Target_seq"].upper()

                #print gRNA_s,gRNA_e
                #plot range will be positions -20 - +20 of the PAM/gRNA sequences.

                for char in range(-40,-20):
                    l = [Rank, char, "-"]
                    mutation_profiles_plot.append(l)
                for char in range(-20,0):
                    l = [Rank, char, "gRNA"]
                    mutation_profiles_plot.append(l)
                for char in range(0,3):
                    l = [Rank, char, "PAM"]
                    mutation_profiles_plot.append(l)
                for char in range(3,21):
                    l = [Rank, char, "-"]
                    mutation_profiles_plot.append(l)
                Rank -=1
                c = 0
                wt = [Rank,-41,"WT"]
                mutation_profiles_plot.append(l)

                if strand == "1":
                    for i in range(-40,21):
                        l = [Rank,i,target_seq[gRNA_s-20+c-1]]
                        c +=1
                        #print l
                        mutation_profiles_plot.append(l)
                    Rank -=1

                if strand == "-1":
                    gRNA_pos = range(gRNA_s,gRNA_e+1)
                    PAM      = range(gRNA_s-3,gRNA_s)
                    for i in range(-40,21):
                        seq = rev_comp( target_seq[gRNA_s-25:gRNA_e+20])
                        l = [Rank,i,seq[c] ]
                        c +=1
                        #print l
                        mutation_profiles_plot.append(l)
                    Rank -=1



            #Pos -41 is for frameshift info of profile




                for row in dat[target_site][plate]:
                    for col in dat[target_site][plate][row]:
                        position = TR96_convert(row,col)
                        if position != "N.A.":
                            pos_profile_mem[position] = {}

                            target_loci = dat[target_site][plate][row][col]["Plus"]
                            tot_reads = sum(dat[target_site][plate][row][col]["Plus"].values())
                            for mutation_profile in dat[target_site][plate][row][col]["Plus"]:
                                Frameshift = "WT"
                                hits = dat[target_site][plate][row][col]["Plus"][mutation_profile]
                                if hits > 1:
                                    hit_rat , error = ErrorProp_div(hits,tot_reads)
                                    imputed_hitrate = hit_rat-error
                                    if (imputed_hitrate >= abundance_threashold) :
                                        #print mutation_profile
                                        mut_l = btop_deconvolute(info[target_site]["Target_seq"],mutation_profile,strand)
                                        target_mut = []
                                        c = 0

                                        if strand == "1":
                                            for i in range(-40,21):
                                                target_mut.append(mut_l[gRNA_s-20+c-1])
                                                c +=1
                                        if strand == "-1":
                                            seq = rev_comp( target_seq[gRNA_s-25:gRNA_e+20])
                                            rev_mut_l = mut_l[::-1]
                                            #print target_site,strand, len(info[target_site]["Target_seq"])
                                            #print mut_l
                                            #print rev_mut_l
                                            for i in range(-40,21):
                                                target_mut.append(rev_mut_l[gRNA_s-25+c])
                                                c +=1

                                        try:
                                            mutation_profiles[("").join(target_mut)]  +=1
                                        except KeyError:
                                            mutation_profiles[("").join(target_mut)]   = 1
                                            mutation_profiles_mem[("").join(target_mut)] = target_mut

                                        pos_profile_mem[position][("").join(target_mut)] = imputed_hitrate


                                        if (target_mut.count("Mut_A")+target_mut.count("Mut_T")+target_mut.count("Mut_G")+target_mut.count("Mut_C")) > 0:
                                            Frameshift = "Mutation"
                                        if (target_mut.count("Ins_A")+target_mut.count("Ins_T")+target_mut.count("Ins_G")+target_mut.count("Ins_C")) >0:
                                            if ((target_mut.count("Ins_A")+target_mut.count("Ins_T")+target_mut.count("Ins_G")+target_mut.count("Ins_C"))%3 != 0):
                                                Frameshift = "Indel(frameshift)"
                                            else:
                                                Frameshift = "Indel(non-frameshift)"
                                        if target_mut.count("Del") > 0:
                                            if (target_mut.count("Del")%3 != 0):
                                                Frameshift = "Indel(frameshift)"
                                            else:
                                                Frameshift = "Indel(non-frameshift)"
                                        if  ((target_mut.count("Del")>0) and ((target_mut.count("Ins_A")+target_mut.count("Ins_T")+target_mut.count("Ins_G")+target_mut.count("Ins_C"))>0)):
                                            if ((target_mut.count("Ins_A")+target_mut.count("Ins_T")+target_mut.count("Ins_G")+target_mut.count("Ins_C")) +(target_mut.count("Del"))%3 != 0):
                                                Frameshift = "Indel(frameshift)"
                                            else:
                                                Frameshift = "Indel(non-frameshift)"

                positions96 = [
                'A1',	'A2',	'A3',	'A4',	'A5',	'A6',	'A7',	'A8',	'A9',	'A10',	'A11',	'A12',
                'B1',	'B2',	'B3',	'B4',	'B5',	'B6',	'B7',	'B8',	'B9',	'B10',	'B11',	'B12',
                'C1',	'C2',	'C3',	'C4',	'C5',	'C6',	'C7',	'C8',	'C9',	'C10',	'C11',	'C12',
                'D1',	'D2',	'D3',	'D4',	'D5',	'D6',	'D7',	'D8',	'D9',	'D10',	'D11',	'D12',
                'E1',	'E2',	'E3',	'E4',	'E5',	'E6',	'E7',	'E8',	'E9',	'E10',	'E11',	'E12',
                'F1',	'F2',	'F3',	'F4',	'F5',	'F6',	'F7',	'F8',	'F9',	'F10',	'F11',	'F12',
                'G1',	'G2',	'G3',	'G4',	'G5',	'G6',	'G7',	'G8',	'G9',	'G10',	'G11',	'G12',
                'H1',	'H2',	'H3',	'H4',	'H5',	'H6',	'H7',	'H8',	'H9',	'H10',	'H11',	'H12']

                positions96_d = {
                'A1':1,	'A2':2      ,	'A3': 3     ,	'A4':  4    ,	'A5': 5     ,	'A6': 6     ,	'A7':  7    ,	'A8': 8     ,	'A9':9      ,	'A10':10      ,	'A11':11      ,	'A12': 12     ,
                'B1': 13     ,	'B2':14      ,	'B3':15      ,	'B4':  16    ,	'B5': 17     ,	'B6': 18     ,	'B7':19      ,	'B8':20      ,	'B9': 21     ,	'B10': 22     ,	'B11': 23    ,	'B12':  24    ,
                'C1':25      ,	'C2':26      ,	'C3': 27     ,	'C4':28      ,	'C5': 29     ,	'C6': 30     ,	'C7':31      ,	'C8':  32    ,	'C9':   33   ,	'C10': 34     ,	'C11':35      ,	'C12':36      ,
                'D1':37      ,	'D2':38      ,	'D3':39      ,	'D4': 40     ,	'D5':41      ,	'D6':42      ,	'D7':43      ,	'D8':44     ,	'D9':45      ,	'D10':46      ,	'D11':47      ,	'D12': 48     ,
                'E1':49      ,	'E2': 50     ,	'E3':51      ,	'E4':52      ,	'E5':53      ,	'E6':54      ,	'E7':55      ,	'E8':56      ,	'E9':57      ,	'E10':58      ,	'E11':59      ,	'E12':60      ,
                'F1':61      ,	'F2':62      ,	'F3': 63     ,	'F4':64      ,	'F5': 65     ,	'F6':66      ,	'F7':67      ,	'F8':68      ,	'F9':69      ,	'F10':70      ,	'F11':71      ,	'F12':72      ,
                'G1':73      ,	'G2':74      ,	'G3':75      ,	'G4':76      ,	'G5':77      ,	'G6':78      ,	'G7':79      ,	'G8':80      ,	'G9':81      ,	'G10':82      ,	'G11':83      ,	'G12':84      ,
                'H1':85      ,	'H2':86      ,	'H3':87      ,	'H4':88      ,	'H5':89      ,	'H6':90      ,	'H7':91      ,	'H8':92      ,	'H9':93     ,	'H10':94      ,	'H11':95      ,	'H12': 96}

                sorted_profiles = sorted(mutation_profiles.items(), key=operator.itemgetter(1),reverse =True)
                ps              = [i[0] for i in sorted_profiles if i[0] != '-------------------------------------------------------------']

                for positions in positions96:
                    l = [0,positions96_d[positions],0,positions]
                    mutation_pos_plot.append(l)
                    called = 0
                    try:
                        test = pos_profile_mem[positions]
                    except KeyError:
                        #print "N.A.", positions,p
                        #print target_site,Frameshift,positions,0
                        l = [-1, positions96_d[positions],-1,positions]
                        mutation_pos_plot.append(l)
                        called +=1

                    if called == 0:
                        try:
                            l = [-1,positions96_d[positions], pos_profile_mem[positions]['-------------------------------------------------------------'],positions]
                            mutation_pos_plot.append(l)

                        except KeyError:
                            #print "Zero", positions,p
                            l = [-1,positions96_d[positions],0,positions]
                            mutation_pos_plot.append(l)


                        #print target_site,i,"WT","None"
                for p in ps:
                    pos = -40
                    mut_profile  = mutation_profiles_mem[p]

                    if (mut_profile.count("Mut_A")+mut_profile.count("Mut_T")+mut_profile.count("Mut_G")+mut_profile.count("Mut_C")) > 0:
                        Frameshift = "Mutation"
                    if (mut_profile.count("Ins_A")+mut_profile.count("Ins_T")+mut_profile.count("Ins_G")+mut_profile.count("Ins_C")) >0:
                        if ((mut_profile.count("Ins_A")+mut_profile.count("Ins_T")+mut_profile.count("Ins_G")+mut_profile.count("Ins_C"))%3 != 0):
                            Frameshift = "Indel(frameshift)"
                        else:
                            Frameshift = "Indel(non-frameshift)"
                    if mut_profile.count("Del") > 0:
                        if (mut_profile.count("Del")%3 != 0):
                            Frameshift = "Indel(frameshift)"
                        else:
                            Frameshift = "Indel(non-frameshift)"
                    if  ((mut_profile.count("Del")>0) and ((mut_profile.count("Ins_A")+mut_profile.count("Ins_T")+mut_profile.count("Ins_G")+mut_profile.count("Ins_C"))>0)):
                        if ((mut_profile.count("Ins_A")+mut_profile.count("Ins_T")+mut_profile.count("Ins_G")+mut_profile.count("Ins_C")) +(mut_profile.count("Del"))%3 != 0):
                            Frameshift = "Indel(frameshift)"
                        else:
                            Frameshift = "Indel(non-frameshift)"

                    l = [Rank,-41,Frameshift]
                    mutation_profiles_plot.append(l)

                    for i in mut_profile:
                        l = [Rank,pos,i]
                        mutation_profiles_plot.append(l)
                        pos +=1

                    if p != '-------------------------------------------------------------':
                        for positions in positions96:
                            called = 0
                            try:
                                test = pos_profile_mem[positions]
                            except KeyError:
                                #print "N.A.", positions,p
                                #print target_site,Frameshift,positions,0
                                l = [Rank, positions96_d[positions],-1,positions]
                                mutation_pos_plot.append(l)
                                called +=1

                            if called == 0:
                                try:
                                    l = [Rank,positions96_d[positions], pos_profile_mem[positions][p],positions]
                                    mutation_pos_plot.append(l)

                                except KeyError:
                                    #print "Zero", positions,p
                                    l = [Rank,positions96_d[positions],0,positions]
                                    mutation_pos_plot.append(l)








                    #for posi in positions96:
                    Rank -=1
                LL2csv(mutation_profiles_plot,'%s/csv/%s/%s/observed_profiles.csv'%(di,target_site,plate))
                LL2csv(mutation_pos_plot,'%s/csv/%s/%s/mutation_pos.csv'%(di,target_site,plate))
    return ""

################################
#### Small functions
################################

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

def btop_deconvolute(target,btop,strand):
    complement_dict = {'A':'T','T':'A','G':'C','C':'G',"N":"N"}
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
                    if strand =="1":
                        mut_mem.append("Ins_%s"%(m[0]))
                    if strand == "-1":
                        mut_mem.append("Ins_%s"%(complement_dict[m[0]]))
                else:
                    if strand =="1":
                        mut_mem.append("Mut_%s"%(m[0]))
                    if strand == "-1":
                        mut_mem.append("Mut_%s"%(complement_dict[m[0]]))


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
def rev_comp(seq):
    complement_dict = {'A':'T','T':'A','G':'C','C':'G',"N":"N"}
    return "".join([complement_dict[base] for base in reversed(seq)])


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
    try:
        posR = {"R01":"A","R03":"B","R05":"C","R07":"D","R09":"E","R11":"F","R13":"G","R15":"H"}
        posC = {"C02":1,"C04":2,"C06":3,"C08":4,"C10":5,"C12":6,"C14":7,"C16":8,"C18":9,"C20":10,"C22":11,"C24":12}
        pos  = ("").join([posR[PosR],str(posC[PosC])])
        return pos
    except KeyError:
        return "N.A."





if __name__ == '__main__':
    dat  = sys.argv[1]
    info = sys.argv[2]
    strand = sys.argv[3] #"Plus" or "Minus"
    abundance_threashold = float(sys.argv[4])#Default = 0.1
    di = sys.argv[5]#Directory name. Defalt "."
    p = sys.argv[6]
    #cellv = sys.argv[7]
    #coord = sys.argv[7]

    main(dat,info,strand,abundance_threashold,di,p)
