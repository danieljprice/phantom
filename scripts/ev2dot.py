#!/usr/bin/env python
import numpy  as np
import argparse
import sys
import evfiles as ev

if sys.version_info[0] < 3:
    raise SystemExit('ERROR: you must use Python 3')

description="""
Use this script to take the time derivative of selected columns in a phantom ev file.
If no columns are given, the derivative of all columns is computed.

If time column cannot be determined automatically, it is assumed to be the first column.

If multiple ev files are given, they are assumed to be in order, and are first combined
with any overlap removed, before taking the derivative.
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('-f', '--file',   required=True, nargs='+', metavar='evfile', help='phantom ev file')
parser.add_argument('-c', '--columns', default=None, nargs='+', metavar='column', help='columns to take time derivative of')
parser.add_argument('-dt','--dtmin',   default=None, nargs=1,   metavar='dtmin',  help='minimum dt to take derivative over')
args = parser.parse_args()

data   = ev.load_files(args.file)  # Load evfiles into a single, combined, pandas data frame
labels = args.columns              # List of columns to take time derivative of

# Get columns
column_names = list(data)
ncols        = len(column_names)

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

column_labels = ['time'] + [l + ' dot' for l in labels]
ev.printev(data=evdot[:irow_deriv,:],labels=column_labels)
