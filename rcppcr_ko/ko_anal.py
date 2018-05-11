#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
Short code to analyze outputs from Call_KO.py
"""
import sys
import re
import os
import pickle
import ast
import json
import numpy as np
import argparse
import pprintc_pcr_coordinate
import xlrd
def main(args):
    #print args
    plate_tag = ["P01-P01","P04-P04","P05-P05","P06-P06","P07-P07"]
    target_list = extract_targets(t)

    D         = {}
    rows,cols = rc_pcr_coordinate(args.rc_pcr_coordinate)
    wells     = list2combinatoral(rows,cols)
    for target in target_list:
        D[target] = {}
        for plate in plate_tag:
            D[target][plate] = {}
            for well in wells:
                D[target][plate][well] = {}

    dataD  = summary2dict(summary,D)
    #pprint.pprint(dataD)
    dataD2 = read_viability_files(args.cell_viability,args.rc_pcr_coordinate,dataD)
    pprint.pprint(dataD2)








#####################
#    Functions
#####################

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
                    r,c = xlsx_pos_TR96_convert(well)
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

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

def xlsx_pos_TR96_convert(well):
    r = well.split("-")[0]
    c = well.split("-")[1]

    R = {"R01":4,"R03":5,"R05":6,"R07":7,"R09":8,"R11":9,"R13":10,"R15":11}
    C = {"C02":2,"C04":3,"C06":4,"C08":5,"C10":6,"C12":7,"C14":8,"C16":9,"C18":10,"C20":11,"C22":12,"C24":13}

    posR = {"R01":"A","R03":"B","R05":"C","R07":"D","R09":"E","R11":"F","R13":"G","R15":"H"}
    posC = {"C02":1,"C04":2,"C06":3,"C08":4,"C10":5,"C12":6,"C14":7,"C16":8,"C18":9,"C20":10,"C22":11,"C24":12}
    pos  = ("").join([posR[r],posC[c]])
    return R[r],C[c],pos



if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='ko_anal')
    parser.add_argument('-s','--summary_file', default=False,help='Output of RCP-PCR_KO analysis script in csv format.')
    parser.add_argument('-t','--targets', default=False,help='Input target informtion in csv format.')
    parser.add_argument('-out','--output_name', default="output" )
    parser.add_argument('-rc','--rc_pcr_coordinate', default="TR96" )
    parser.add_argument('-cellv','--cell_viability', default="." )


    args = parser.parse_args()
    main(args)
