#!/usr/bin/env python3

# This script selects those who fail heterozygosity constraints -- previous versions also looked at missingness too
# hence the name


from pandas import read_csv
import sys
import argparse

TAB = chr(9)

parser.add_argument("--het",action='store', help="..",required = True)
parser.add_argument("--cut_het_low",action='store', help="..",required = True)
parser.add_argument("--cut_het_high",action='store', help="..",required = True)
parser.add_argument("--outfname",action='store', help="..",required = True)

args = parser.parse_args()

if len(sys.argv)<=1:
  sys.argv = ["select_miss_het.qcplink.py","$het",float("${params.cut_het_low}"),float("${params.cut_het_high}"),"$outfname"]

hetf = args.het
cut_het_high=float(args.cut_het_high)
cut_het_low=float(args.cut_het_low);
outfname = args.outfname


het      = read_csv(hetf,delim_whitespace=True)
mean_het = (het["N(NM)"]-het["O(HOM)"])/het["N(NM)"]
failed   = het[(mean_het<cut_het_low) | (mean_het>cut_het_high)]

failed.to_csv(outfname,header=False,sep=TAB)
