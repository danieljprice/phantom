#!/usr/bin/env python
import pandas as pd
import numpy  as np
import argparse
import os
import sys
import re

if sys.version_info[0] < 3:
    raise SystemExit('ERROR: you must use Python 3')

# Printing options
np.set_printoptions(threshold=np.inf,floatmode='fixed',linewidth=np.inf)
fmt = '%18.10E'

description="""
Use this script cut out certain columns from a phantom ev file, and output to stdout
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('-f','--file',   required=True, nargs=1,   metavar='evfile', help='phantom ev file')
parser.add_argument('-c','--columns',required=True, nargs='+', metavar='column', help='column names to extract')
args = parser.parse_args()
file = args.file[0]

# Check that file exist
if not os.path.isfile(file):
    raise SystemExit('Error: the file "'+file+'" does not exist or is not a file')

# Read the columns names in the file using regex
with open(file,'r') as f:
    header       = f.readline().strip('\n')
    column_names = re.findall('\[\d+\s+(.*?)\]',header)

# Number of columns in file
ncols = len(column_names)

data = pd.read_csv(file,skiprows           = 0,
                        dtype              = np.float64,
                        delim_whitespace   = True,
                        skip_blank_lines   = True,
                        comment            = '#',
                        header             = None,
                        names              = column_names)

# List of columns to cut from file
labels = args.columns

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

# Construct the header string
header='# '
for i in range(len(labels)):
    label  = labels[i].rjust(12)
    no     = str(i+1).zfill(2)
    header = header + '[' + no + label + ']' + '   '

np.savetxt(sys.stdout.buffer,d,header=header,fmt=fmt,comments='')
