#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 00:14:32 2021

@author: coyote
"""

import json
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Process WDL input and SraRunTable.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('inputs', metavar='I', nargs='+',
                    help='json inputs file with an empty ".samples" array')
parser.add_argument('num_reads', metavar='N', type=int, nargs='+',
                    help='number of reads per sample', default=2)
parser.add_argument('-s', '--sra-table', dest='sra', 
                    help='the Sra run table for project')
parser.add_argument('-r', '--read-group', type=str, 
                    help='map of RG tags and run table columns in tag/column\n \
                    separated by comma',
                    default="ID/Library Name,SM/Run,LB/Library Name,PL/Platform")
parser.add_argument('-f', '--folder', type=str, default="reads",
                    help='the folder which the raw reads are kept')

args = parser.parse_args()

run_table = pd.read_csv(args.sra)
run_table.columns

with open(args.inputs[0]) as f:
    inputs = json.load(f)

rg = [i.split("/") for i in args.read_group.split(",")]
samples = []

for idx, i in run_table.iterrows():
    if args.num_reads[0] == 1:
        s = {"outputPrefix" : i["Run"],
             "read1" : "reads/" + i["Run"] + ".fastq",
             "readgroup" : "@RG" + ''.join(["\\t" + r[0] + ":" + i[r[1]] for r in rg])
             }
        samples.append(s)
    else:
        s = {"outPutprefix" : i["Run"],
             "read1" : i["Run"] + "_1.fastq",
             "read2" : i["Run"] + "_2.fastq",
             "readgroup" : "@RG" + ''.join(["\\t" + r[0] + ":" + i[r[1]] for r in rg])
             }
        samples.append(s)
            

inputs[list({key for key in inputs.keys() if key.endswith("samples")})[0]] = samples    
print(json.dumps(inputs, indent=2))
