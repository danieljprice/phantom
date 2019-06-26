#!/usr/bin/env python
import pandas as pd
import numpy  as np
import argparse
import os
import sys

if sys.version_info[0] < 3:
    raise SystemExit('ERROR: you must use Python 3')

# Printing options
np.set_printoptions(threshold=np.inf,floatmode='fixed',linewidth=np.inf)
fmt = '%18.10E'

description="""
Use this script to combine multiple phantom ev files together, removing overlap, and output to stdout
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('evfiles', nargs='+', metavar='file', help='phantom ev file')
args = parser.parse_args()

# Check that all files exist
for file in args.evfiles:
    if not os.path.isfile(file):
        raise SystemExit('Error: the file "'+file+'" does not exist or is not a file')

# Check that all files have identical headers (i.e. same number of columns)
with open(args.evfiles[0],'r') as file:
    prev_header = file.readline()
for filename in args.evfiles:
    with open(filename,'r') as file:
        header = file.readline()
        if header == prev_header:
            prev_header = header
        else:
            raise SystemExit('Error: Column labels of file: "'+filename+'" are different to the rest')

header = header.strip('\n')
nfiles = len(args.evfiles)

for i in range(nfiles-1):
    data = np.array(pd.read_csv(args.evfiles[i],skiprows           = 0,
                                                dtype              = np.float64,
                                                delim_whitespace   = True,
                                                skip_blank_lines   = True,
                                                comment            = '#',
                                                header             = None))

    data_next = np.array(pd.read_csv(args.evfiles[i+1],skiprows           = 0,
                                                       dtype              = np.float64,
                                                       delim_whitespace   = True,
                                                       skip_blank_lines   = True,
                                                       comment            = '#',
                                                       header             = None))

    # Find the index at which the current file begins to overlap with the next file
    for ii in range(len(data[:,0])):
        t         = data[ii,0]
        tnextfile = data_next[0,0]
        if t > tnextfile:
            break

    # Only write the header once
    if i!=0: header=''

    # Print the unique parts of each file
    np.savetxt(sys.stdout.buffer,data[:ii-1,:],header=header,fmt=fmt,comments='')

if nfiles>1:
    # Print all of the last file
    np.savetxt(sys.stdout.buffer,data_next[:,:],header=header,fmt=fmt,comments='')
else:
    # Just print the contents of the 1 file present
    data = np.array(pd.read_csv(args.evfiles[0],skiprows           = 0,
                                                dtype              = np.float64,
                                                delim_whitespace   = True,
                                                skip_blank_lines   = True,
                                                comment            = '#',
                                                header             = None))
    np.savetxt(sys.stdout.buffer,data[:,:],header=header,fmt=fmt,comments='')
