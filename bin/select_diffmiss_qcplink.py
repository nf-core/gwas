#!/usr/bin/env python3

import sys
import pandas as pd
import argparse

EOL = chr(10)

parser = argparse.ArgumentParser()

parser.add_argument("--missing",action='store', help="..",required = True)
parser.add_argument("--probcol",action='store', help="..",required = True)
parser.add_argument("--cut_diff_miss",action='store', help="..",required = True)
parser.add_argument("--failed",action='store', help="..",required = True)

args = parser.parse_args()


mfr     = pd.read_csv(args.missing,delim_whitespace=True)
probcol = args.probcol
cut_diff_miss = float(args.cut_diff_miss)

wanted = mfr[mfr[probcol]<cut_diff_miss]["SNP"].values
out = open(args.failed,"w")
out.write(EOL.join(map(str,wanted)))
out.close()
