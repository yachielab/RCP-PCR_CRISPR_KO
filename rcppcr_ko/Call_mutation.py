#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
Python script to analyse mutation frequency on target positions from RCP-PCR samples.
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

if __name__ == '__main__':
    dat  = sys.argv[1]
    info = sys.argv[2]
    di = sys.argv[3]#Directory name. Defalt "."
    p = sys.argv[4]

    main(dat,info,di,p)
