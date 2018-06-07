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

def main(dat,info,di,p):
    big_d = reading(dat)
    info_d= csv2dict(info)
    #pprint.pprint(info_d)
    #pprint.pprint(big_d)
    mutation_profiles(big_d,info_d,di)
    #pprint.pprint(LL_for_R)
    """
    for target in heatmap_d:
        if not os.path.isdir('%s/csv/%s'%(di,target)):
            os.makedirs('%s/csv/%s'%(di,target))
        Fname_h = "%s/csv/%s/heatmap.csv"%(di,target)
        Fname_b = "%s/csv/%s/barplot.csv"%(di,target)
        LL2csv(heatmap_d[target],Fname_h)
        LL2csv(bar_d[target],Fname_b)
        PDF_name = ("%s/pdf/%s_plot.pdf"%(di,target))
        os.system("Rscript %srcppcr_ko_heatmap.r %s %s %s"%(p,Fname_b,Fname_h,PDF_name))
    """

def mutation_profiles(dat,info,di):
    """
    Dataformat for each target-Plate-Row-Col:
    [
    #Seq    ["A","T","C",...],
    #->A    ["0.00","0.00","0.23","0.22",...],
    #->T    ["0.00","0.00","0.23","0.23",...],
    #->G    ["0.00","0.00","0.23","0.01",...],
    #->C    ["0.00","0.00","0.23","0.02",...],
    #->del  ["0.001","0.001","0.01","0.021",...],
    #->Ins  ["0.002","0.001","0.00","0.021",...]
    ]

    Final data format;

    [
    [Target, Plate, Row, Col, Mutation_type,Substitution_type, Position, Value],
    [],...
    ]

    """

    out_dat = {}
    out_bar = {}
    LL_for_csv = [["Target","Plate","Row","Col","Mutation_type","Substitution_type","Position","Value"]  ]
    raw_ratio  = [[ "Target","Plate","Row","Col","Mutation_type","Substitution_type"]]

    SUM = {}
    F = {}
    for target_site in dat:
        for plate in dat[target_site]:
            if ((plate == "P11-P11") or  (plate == "P12-P12")):
                for row in dat[target_site][plate]:
                    for col in dat[target_site][plate][row]:
                        #print target_site,plate,row,col
                        tot = 0
                        target_region_s = int(info[target_site]["gRNA_s"]) - 23
                        target_region_e = int(info[target_site]["gRNA_e"]) + 18
                        #print int(info[target_site]["gRNA_s"]),int(info[target_site]["gRNA_e"])
                        #print target_site,target_region_s,target_region_e,len(info[target_site]["Target_seq"][target_region_s:target_region_e]),info[target_site]["Target_seq"][target_region_s:target_region_e]

                        seq_array = list(info[target_site]["Target_seq"][target_region_s:target_region_e])
                        L = [target_site,plate,row,col,"Target","Target"]
                        L += seq_array
                        raw_ratio.append(L)

                        mut_A     = np.zeros(len(seq_array))
                        mut_T     = np.zeros(len(seq_array))
                        mut_G     = np.zeros(len(seq_array))
                        mut_C     = np.zeros(len(seq_array))
                        mut_del   = np.zeros(len(seq_array))
                        mut_ins   = np.zeros(len(seq_array))

                        for profiles in dat[target_site][plate][row][col]["Plus"]:
                            tot += dat[target_site][plate][row][col]["Plus"][profiles]
                        for profiles in dat[target_site][plate][row][col]["Plus"]:
                            mut_info = btop_deconvolute(info[target_site]["Target_seq"],profiles)

                            print mut_info
                            A_indices = [i for i, x in enumerate(mut_info[target_region_s:target_region_e]) if x == "A"]
                            T_indices = [i for i, x in enumerate(mut_info[target_region_s:target_region_e]) if x == "T"]
                            G_indices = [i for i, x in enumerate(mut_info[target_region_s:target_region_e]) if x == "G"]
                            C_indices = [i for i, x in enumerate(mut_info[target_region_s:target_region_e]) if x == "C"]
                            del_indices = [i for i, x in enumerate(mut_info[target_region_s:target_region_e]) if x == "Del"]
                            ins_indices = [i for i, x in enumerate(mut_info[target_region_s:target_region_e]) if x == "Ins"]



                            for A in A_indices:
                                mut_A[A] += float(dat[target_site][plate][row][col]["Plus"][profiles])
                            for T in T_indices:
                                mut_T[T] += float(dat[target_site][plate][row][col]["Plus"][profiles])
                            for G in G_indices:
                                mut_G[G] += float(dat[target_site][plate][row][col]["Plus"][profiles])
                            for C in C_indices:
                                mut_C[C] += float(dat[target_site][plate][row][col]["Plus"][profiles])
                            for Del in del_indices:
                                mut_del[Del] += float(dat[target_site][plate][row][col]["Plus"][profiles])
                            for ins in ins_indices:
                                mut_ins[ins] += float(dat[target_site][plate][row][col]["Plus"][profiles])


                        mut_A = mut_A/tot
                        mut_T = mut_T/tot
                        mut_G = mut_G/tot
                        mut_C = mut_C/tot
                        mut_del = mut_del/tot
                        mut_ins = mut_ins/tot
                        mut_tot = mut_A+mut_T+mut_G+mut_C+mut_del+mut_ins

                        for i in range(len(seq_array)):
                            LL_for_csv.append( [target_site,plate,row,col,"Deletion","-",i,mut_del[i]] )
                            LL_for_csv.append( [target_site,plate,row,col,"Insertion","-",i,mut_ins[i]] )
                            LL_for_csv.append( [target_site,plate,row,col,"Mutation","A",i,mut_A[i]] )
                            LL_for_csv.append( [target_site,plate,row,col,"Mutation","T",i,mut_T[i]] )
                            LL_for_csv.append( [target_site,plate,row,col,"Mutation","G",i,mut_G[i]] )
                            LL_for_csv.append( [target_site,plate,row,col,"Mutation","C",i,mut_C[i]] )


                        L = [target_site,plate,row,col,"Deletion","-"]
                        #print mut_del
                        L += [str(i) for i in  mut_del]
                        raw_ratio.append(L)
                        L = [target_site,plate,row,col,"Insertion","-"]
                        L += [str(i) for i in mut_ins]
                        raw_ratio.append(L)
                        L = [target_site,plate,row,col,"Mutation","A"]
                        L +=[str(i) for i in  mut_A]
                        raw_ratio.append(L)
                        L = [target_site,plate,row,col,"Mutation","T"]
                        L += [str(i) for i in mut_T]
                        raw_ratio.append(L)
                        L = [target_site,plate,row,col,"Mutation","G"]
                        L += [str(i) for i in mut_G]
                        raw_ratio.append(L)
                        L = [target_site,plate,row,col,"Mutation","C"]
                        L += [str(i) for i in mut_C]
                        raw_ratio.append(L)
                        L = [target_site,plate,row,col,"Mutation","Total"]
                        L += [str(i) for i in mut_tot]
                        raw_ratio.append(L)


    LL2csv(LL_for_csv,"mutation_sumary_for_plot.csv" )
    LL2csv(raw_ratio,"mutation_sumary.csv" )
    return out_dat, out_bar

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
                    mut_mem.append(m[0])
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
    di = sys.argv[3]#Directory name. Defalt "."
    p = sys.argv[4]

    main(dat,info,di,p)
