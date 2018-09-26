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
import numpy as np
import math
import pprint
import operator

def main(dat,info,strand,abundance_threashold,di,p):
    big_d = reading(dat)
    info_d= csv2dict(info)
    #pprint.pprint(info_d)
    #pprint.pprint(big_d)
    a = get_top_profiles(big_d,info_d,di,abundance_threashold)


def get_top_profiles(dat,info,di,abundance_threashold):

    sample_index = ["P01-P01"]#,"P02-P02","P03-P03","P04-P04","P18-P18","P19-P19","P21-P21","P22-P22"]


    for target_site in dat:

        for plate in dat[target_site]:
            #print plate
            CLONE_NUM = 0
            if plate in sample_index:
                if not os.path.isdir('%s/csv/%s/%s'%(di,target_site,plate)):
                    os.makedirs('%s/csv/%s/%s'%(di,target_site,plate))
                mutation_profiles      = {}#btop of profile : #of wells
                mutation_profiles_mem  = {}

                pos_profile_mem = {}


                mutation_profiles_plot = [["Rank","Position","Char_stat"]]
                pam_ref_profiles_plot = [["Rank","Position","Char_stat","String"]]
                ref_profiles_plot = [["Rank","Position","Char_stat"]]
                mutation_pos_plot = [["Rank","Position","Frequency"]]
                readnum = [["Clone","Readnum"]]
                genotype_stat = [["Rank","Position","Stat"]]
                genotype_detail = [["Rank","Position","Stat"]]
                Rank = 0 # 0 is always for reference (PAM, sgRNA)

                strand = info[target_site]["Strand"]
                gRNA_s = int(info[target_site]["gRNA_s"])
                gRNA_e = int(info[target_site]["gRNA_e"])
                target_seq = info[target_site]["Target_seq"].upper()

                #print gRNA_s,gRNA_e
                #plot range will be positions -20 - +20 of the PAM/gRNA sequences.
                protospacer = "    Protospacer     "
                pam = "PAM"

                for char in range(-40,-20):
                    l = [1, char, "-",""]
                    pam_ref_profiles_plot.append(l)
                for char in range(-20,0):

                    l = [1, char, "gRNA",protospacer[char+20]]
                    pam_ref_profiles_plot.append(l)
                for char in range(0,3):
                    l = [1, char, "PAM",pam[char]]
                    pam_ref_profiles_plot.append(l)
                for char in range(3,21):
                    l = [1, char, "-",""]
                    pam_ref_profiles_plot.append(l)

                c = 0
                #wt = [0,-41,"WT"]
                #ref_profiles_plot.append(l)

                if strand == "1":
                    for i in range(-40,21):
                        l = [0,i,target_seq[gRNA_s-20+c-1]]
                        c +=1
                        #print l
                        ref_profiles_plot.append(l)

                if strand == "-1":
                    gRNA_pos = range(gRNA_s,gRNA_e+1)
                    PAM      = range(gRNA_s-3,gRNA_s)
                    for i in range(-40,21):
                        seq = rev_comp( target_seq[gRNA_s-25:gRNA_e+20])
                        l = [0,i,seq[c] ]
                        c +=1
                        #print l
                        ref_profiles_plot.append(l)




            #Pos -41 is for frameshift info of profile




                for row in dat[target_site][plate]:
                    for col in dat[target_site][plate][row]:
                        position = TR96_convert(row,col)
                        if position != "N.A.":
                            CLONE_NUM += 1
                            pos_profile_mem[CLONE_NUM] = {}

                            target_loci = dat[target_site][plate][row][col]["Plus"]
                            tot_reads = sum(dat[target_site][plate][row][col]["Plus"].values())
                            readnum.append([CLONE_NUM,tot_reads])
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

                                        pos_profile_mem[CLONE_NUM][("").join(target_mut)] = imputed_hitrate


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


                sorted_profiles = sorted(mutation_profiles.items(), key=operator.itemgetter(1),reverse =True)
                ps              = [i[0] for i in sorted_profiles ]#if i[0] != '-------------------------------------------------------------']


                        #print target_site,i,"WT","None"
                for p in ps:
                    pos = -40
                    mut_profile  = mutation_profiles_mem[p]

                    if (mut_profile.count("Mut_A")+mut_profile.count("Mut_T")+mut_profile.count("Mut_G")+mut_profile.count("Mut_C")) > 0:
                        Frameshift = "Mutation"
                        Frame_shift_detail = "AA -> AA"
                    elif (mut_profile.count("Ins_A")+mut_profile.count("Ins_T")+mut_profile.count("Ins_G")+mut_profile.count("Ins_C")) >0:
                        if ((mut_profile.count("Ins_A")+mut_profile.count("Ins_T")+mut_profile.count("Ins_G")+mut_profile.count("Ins_C"))%3 != 0):
                            Frameshift = "Ins"
                            Frame_shift_detail = "Out"
                        else:
                            Frameshift = "Ins"
                            Frame_shift_detail = "In"
                    elif mut_profile.count("Del") > 0:
                        if (mut_profile.count("Del")%3 != 0):
                            Frameshift = "Del"
                            Frame_shift_detail = "Out"
                        else:
                            Frameshift = "Del"
                            Frame_shift_detail = "In"
                    elif  ((mut_profile.count("Del")>0) and ((mut_profile.count("Ins_A")+mut_profile.count("Ins_T")+mut_profile.count("Ins_G")+mut_profile.count("Ins_C"))>0)):
                        if ((mut_profile.count("Ins_A")+mut_profile.count("Ins_T")+mut_profile.count("Ins_G")+mut_profile.count("Ins_C")) +(mut_profile.count("Del"))%3 != 0):
                            Frameshift = "Indel"
                            Frame_shift_detail = "Out"
                        else:
                            Frameshift = "Indel"
                            Frame_shift_detail = "In"
                    else:
                        Frameshift = "WT"
                        Frame_shift_detail = "-"

                    l = [Rank,0,Frameshift]
                    ll = [Rank,0,Frame_shift_detail]
                    #mutation_profiles_plot.append(l)

                    genotype_stat.append(l)
                    genotype_detail.append(ll)



                    for i in mut_profile:
                        l = [Rank,pos,i]
                        mutation_profiles_plot.append(l)
                        pos +=1

                    for clones in range(CLONE_NUM):
                        called = 0
                        try:
                            test = pos_profile_mem[clones]
                        except KeyError:
                            #print "N.A.", positions,p
                            #print target_site,Frameshift,positions,0
                            l = [Rank, clones+1,-1]
                            mutation_pos_plot.append(l)
                            called +=1

                        if called == 0:
                            try:
                                l = [Rank,clones, pos_profile_mem[clones][p]]
                                mutation_pos_plot.append(l)

                            except KeyError:
                                #print "Zero", positions,p
                                l = [Rank,clones,0]
                                mutation_pos_plot.append(l)







                    #for posi in positions96:
                    Rank -=1
                LL2csv(mutation_profiles_plot,'%s/csv/%s/%s/observed_profiles.csv'%(di,target_site,plate))
                LL2csv(mutation_pos_plot,'%s/csv/%s/%s/mutation_pos.csv'%(di,target_site,plate))
                LL2csv(ref_profiles_plot,'%s/csv/%s/%s/ref_profiles.csv'%(di,target_site,plate))
                LL2csv(pam_ref_profiles_plot,'%s/csv/%s/%s/pam_ref_profiles.csv'%(di,target_site,plate))
                LL2csv(genotype_stat,'%s/csv/%s/%s/genotype_stat.csv'%(di,target_site,plate))
                LL2csv(genotype_detail,'%s/csv/%s/%s/genotype_detail.csv'%(di,target_site,plate))
                LL2csv(readnum,'%s/csv/%s/%s/readnum.csv'%(di,target_site,plate))
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
