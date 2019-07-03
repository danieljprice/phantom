#!/usr/bin/env python
import numpy  as np
import argparse
import sys
import evfiles as ev

if sys.version_info[0] < 3:
    raise SystemExit('ERROR: you must use Python 3')

description="""
Use this script to smooth out data in a phantom ev file, given a delta t
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('evfiles',           nargs='+',metavar='file',help='phantom ev file')
parser.add_argument('-dt',    default=0.,nargs=1,  metavar='dt',  help='dt to smooth over',type=float)

args = parser.parse_args()

data = ev.load_files(args.evfiles)
column_labels = list(data)

# Determine which column is time.
tcol = set(column_labels).intersection(set(['t','time','Time','TIME','T']))

if len(tcol)==1:
    time = np.array(data[list(tcol)[0]])
else:
    # If no unique column can be determined, assume time is the first column
    time = np.array(data)[:,0]

nrows    = len(time)
evsmooth = np.zeros((nrows,len(column_labels)))

d = np.array(data)

# Smooth out all the columns by given dt
irowprev   = 0
irow_smooth = 0
for irow in range(nrows):
    deltat = time[irow]-time[irowprev]
    if deltat > args.dt:
        evsmooth[irow_smooth,:]  = 0.5*(d[irow,:]+d[irowprev,:])
        irowprev     = irow
        irow_smooth += 1

ev.printev(evsmooth[:irow_smooth,:],column_labels)
