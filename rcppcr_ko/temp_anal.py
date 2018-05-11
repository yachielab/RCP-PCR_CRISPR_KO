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

def main(t,summary):
    target_list = extract_targets(t)
    pass









def extract_targets(f):
    targets = []
    with open(f,"r") as F:
        for line in F:
            cols = line.split(",")
            targets.append(cols[0])
        F.close()
    return targets


def read_viability_files(f):
    with open(f,"r") as F:
        print line






if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='ko_anal')
    parser.add_argument('-s','--summary_file', default=False,help='Output of RCP-PCR_KO analysis script in csv format.')
    parser.add_argument('-t','--targets', default=False,help='Input target informtion in csv format.')
    parser.add_argument('-out','--output_name', default="output" )
    args = parser.parse_args()
    main()
