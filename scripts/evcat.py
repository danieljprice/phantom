#!/usr/bin/env python
import numpy  as np
import argparse
import sys
import evfiles as ev

if sys.version_info[0] < 3:
    raise SystemExit('ERROR: you must use Python 3')

description="""
Use this script to combine multiple phantom ev files together, removing overlap, and output to stdout
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('evfiles', nargs='+', metavar='file', help='phantom ev file')
args = parser.parse_args()

# Check that all files have identical headers (i.e. same number of columns)
with open(args.evfiles[0],'r') as file:
    prev_header = file.readline()
    column_labels = ev.get_column_names(args.evfiles[0])
for filename in args.evfiles[1:]:
    with open(filename,'r') as file:
        header = file.readline()
        if header == prev_header:
            prev_header = header
        else:
            raise SystemExit('Error: Column labels of file: "'+filename+'" are different to the rest')

nfiles = len(args.evfiles)

for i in range(nfiles-1):
    data      = np.array(ev.load(args.evfiles[i]  ))
    data_next = np.array(ev.load(args.evfiles[i+1]))

    # Find the index at which the current file begins to overlap with the next file
    for ii in range(len(data[:,0])):
        t         = data[ii,0]
        tnextfile = data_next[0,0]
        if t > tnextfile:
            break

    # Only write the header once
    if i!=0: column_labels = None

    # Print the unique parts of each file
    ev.printev(data[:ii-1,:],column_labels)

if nfiles>1:
    # Print all of the last file
    ev.printev(data_next[:,:],labels=None)
else:
    # Just print the contents of the 1 file present
    data = np.array(ev.load(args.evfiles[0]))
    ev.printev(data[:,:],column_labels)
