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
Use this script to take the time derivative of selected columns in a phantom ev file.
If no columns are given, the derivative of all columns is computed.

If time column cannot be determined automatically, it is assumed to be the first column.
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('-f', '--file',   required=True, nargs=1,   metavar='evfile', help='phantom ev file')
parser.add_argument('-c', '--columns', default=None, nargs='+', metavar='column', help='columns to take time derivative of')
parser.add_argument('-dt','--dtmin',   default=None, nargs=1,   metavar='dtmin',  help='minimum dt to take derivative over')
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

if labels is not None:
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

else:
    # Take derivative of all columns if none are requested
    labels = column_names[1:]

# Construct numpy array with only the columns desired
d = np.array(data[labels])
nderivs = len(labels)


# Determine which column is time.
tcol = set(column_names).intersection(set(['t','time','Time','TIME','T']))

if len(tcol)==1:
    time = np.array(data[list(tcol)[0]])
else:
    # If no unique column can be determined, assume time is the first column
    time = np.array(data)[:,0]

# Set dtmin from argparse. By default assume dtmin = 0
if args.dtmin is None:
    dtmin = 0.
else:
    try:
        dtmin = float(args.dtmin[0])
    except:
        raise SystemExit('Error: bad value for dtmin. Could not convert "'+args.dtmin[0]+'" to float.')

nrows = len(time)
evdot = np.zeros((nrows,1+nderivs))

# Take the time derivative of all the requested columns using the desired dtmin
irowprev   = 0
irow_deriv = 0
for irow in range(nrows):
    deltat = time[irow]-time[irowprev]
    if deltat > dtmin:
        evdot[irow_deriv,0]  = time[irow]
        evdot[irow_deriv,1:] = (d[irow,:]-d[irowprev,:])/deltat
        irowprev    = irow
        irow_deriv += 1

# Construct the header string
header='# [01        time]   '
for i in range(len(labels)):
    label  = labels[i]+' dot'
    label  = label.rjust(12)
    no     = str(i+2).zfill(2)
    header = header + '[' + no + label + ']' + '   '

np.savetxt(sys.stdout.buffer,evdot[:irow_deriv,:],header=header,fmt=fmt,comments='')
