#!/usr/bin/env python
import numpy  as np
import argparse
import sys
import evfiles as ev

if sys.version_info[0] < 3:
    raise SystemExit('ERROR: you must use Python 3')

description="""
Use this script cut out certain columns from a phantom ev file, and output to stdout
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('-f','--file',   required=True, nargs=1,   metavar='evfile', help='phantom ev file')
parser.add_argument('-c','--columns',required=True, nargs='+', metavar='column', help='column names to extract')
args = parser.parse_args()

data   = ev.load(args.file[0])  # Load evfile into a pandas data frame
labels = args.columns           # List of columns to cut from file

# Get columns
column_names = list(data)
ncols        = len(column_names)

# Replace integer columns with corresponding column name
for i in range(len(labels)):
    l = labels[i]
    if l not in column_names:
        try: # to convert to integer
            icol = int(l)
        except:
            raise SystemExit('Could not find column "'+l+'"')
        else:
            if icol-1>ncols:
                raise SystemExit('Could not find column "'+l+'", out of range 1-'+str(ncols))
            labels[i] = column_names[icol-1]

# Construct numpy array with only the columns desired
d = np.array(data[labels])

ev.printev(data=d,labels=labels)
