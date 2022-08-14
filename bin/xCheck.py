#!/usr/bin/env python3


import argparse
import sys
import re
import os
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument("--x",action='store', help="..",required = True)
parser.add_argument("--f_hi_female",action='store', help="..",required = True)
parser.add_argument("--f_lo_male",action='store', help="..",required = True)
parser.add_argument("--out",action='store', help="..",required = True)
parser.add_argument("--missingness",action='store', help="..",required = True)

args = parser.parse_args()


# default    missingness = [0.01,0.03,0.05]

idtypes = dict(map(lambda x: (x,object),["FID","IID"]))

def getCsvI(fn,cols,names=None):
    ''' read a CSV file with IDs as index '''
    frm = pd.read_csv(fn,delim_whitespace=True,usecols=cols,dtype=idtypes)
    frm = frm.set_index(['FID','IID'])
    # nb bug/feature in pandas, if a column is an index then dtype is not applied to it
    # can have numbers as IDs but want them as string
    return frm


def getResForM(base,m):
    out = "%s-%s"%(base,m)
    os.system("plink --bfile %s --mind %s --check-sex --out %s"%(base,m,out))
    cols=["FID","IID","STATUS","F"]
    sf = getCsvI("%s.sexcheck"%out,cols)
    sf[m] = np.where(sf['STATUS']=='OK',"OK", np.where((sf['F']>f_hi_female) & (sf['F']<f_lo_male),"S","H"))
    return sf

def checkNoX(base,outfn):
    line = open("%s.fam"%base).readline()
    if not line:
        pd.to_pickle("NO",outfn)
        sys.exit(0)
    return


base   = args.x
f_hi_female = float(args.f_hi_female)
f_lo_male   = float(args.f_lo_male)
outfn  = args.out

checkNoX(base,outfn)

result = getResForM(base,1)
for m in args.missingness:
    try:
        curr   = getResForM(base,m)
    except IOError:
        result[m]=np.nan
        continue
    result = result.join(curr[m],how='outer')
result=result.drop(['STATUS','F'],axis=1)
pd.to_pickle(result,outfn)


