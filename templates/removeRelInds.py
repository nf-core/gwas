#!/usr/bin/env python3


import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--missing",action='store', help="..",required = True)
parser.add_argument("--ibd_genome",action='store', help="..",required = True)
parser.add_argument("--outfname",action='store', help="..",required = True)
parser.add_argument("--super_pi_hat",action='store', help="..",required = True)

args = parser.parse_args()

EOL=chr(10)

str_type = {"FID":str, "IID":str, "FID1":str, "IID1":str, "FID2":str, 'IID2':str}

imissf = pd.read_csv(args.missing,delim_whitespace=True,dtype=str_type)
imissf.set_index(["FID","IID"],inplace=True)
genomef = pd.read_csv(args.ibd_genome,delim_whitespace=True,usecols=["FID1","IID1","FID2","IID2","PI_HAT"],dtype=str_type)

outf   =open(args.outfname,"w")
super_pi_hat = float(args.super_pi_hat)

def getDegrees(remove):
   elts = set(imissf.index.values)
   degd = {}
   rel  = {}
   for elt in elts: 
       degd[elt]=0
       rel[elt]=[]
   deg = pd.Series(degd)
   for i,row in genomef.iterrows():
       x=tuple(row[["FID1","IID1"]].tolist())
       y=tuple(row[["FID2","IID2"]].tolist())
       if x in remove or y in remove : continue
       try:
           deg[x]=deg[x]+1
           deg[y]=deg[y]+1
           rel[x].append(y)
           rel[y].append(x)
       except:
           print("saw the problem")
   return rel, deg





remove = set(map (tuple,genomef[genomef['PI_HAT']>super_pi_hat][["FID1","IID1"]].to_records(index=False)))\
         | set(map(tuple,genomef[genomef['PI_HAT']>super_pi_hat][["FID2","IID2"]].to_records(index=False)))


rel, deg = getDegrees(remove) 

candidates = deg[deg>=1].sort_values(ascending=False)

for i,c in candidates.iteritems():
    if deg[i]>0:
        remove.add(i)
        deg[i]=deg[i]-1
        for other in rel[i]:
            deg[other]=deg[other]-1

                

remove = sorted(list(remove))
outf.write(EOL.join(map (lambda x: "%s %s"%(x[0],x[1]),remove)))
outf.close()
        
