#!/usr/bin/env python3
from __future__ import print_function
import sys
import os
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("--inpfname",action='store', help="Input filename",required = True)
parser.add_argument("--outfname",action='store', help="Output filename",required = True)
parser.add_argument("--remove_on_bp",action='store', help="??",required = True)


args = parser.parse_args()




def getChrom(chrom):
    try:
        result = int(chrom)
    except ValueError:
        result = chrom
    return result

def removeOnBP(fname, out):
    f = open(fname)
    line = f.readline().strip().split()
    old_chrom = getChrom(line[0])
    old_snp   = line[1]
    old_bp    = int(line[3])
    for line in f:
        data = line.strip().split()
        try:
            chrom = getChrom(data[0])
        except IndexError:
            print("****",line)
            print(data)
            print(old_snp)
            sys.exit(11)
        snp   = data[1]
        bp    = int(data[3])
        if (chrom, bp) == (old_chrom, old_bp):
            out.write("%s\\n"%snp)
            (old_chrom,old_snp,old_bp) = (chrom,snp,bp)
        elif (chrom, bp) > (old_chrom, old_bp):
            (old_chrom,old_snp,old_bp) = (chrom,snp,bp)
        else:
            print(old_chrom, old_bp, old_snp)
            print(chrom,bp,snp)
            print(""" 

              The BIM file <%s> is not sorted see <%s> and <%s>"
            """ % (fname, old_snp, snp))
            sys.exit(10)


os.system("hostname > hostname")
f=open(args.inpfname)
if not f:
    sys.exit("File <%s> not opened"%args.inpfname)
s_name = set()

out=open(args.outfname,"w")

for line in f:
    data=line.strip().split()
    snp_name = data[1]
    if snp_name in s_name:
        out.write( "%s\\n"%snp_name)
    else:
        s_name.add(data[1])
f.close()

if args.remove_on_bp in ["1",1,True,"True","true"]:
    removeOnBP(args.inpfname,out)

out.close()

