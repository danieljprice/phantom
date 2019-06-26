#!/usr/bin/env python
import pandas as pd
import numpy  as np
import argparse
import os
import sys
import re

if sys.version_info[0] < 3:
    raise SystemExit('ERROR: you must use Python 3')

description="""
Use this script to print the column names from an ev file
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('evfile', nargs=1, metavar='file', help='phantom ev file')
args = parser.parse_args()
file = args.evfile[0]

# Check that file exist
if not os.path.isfile(file):
    raise SystemExit('Error: the file "'+file+'" does not exist or is not a file')

# Read the columns names in the file using regex
with open(file,'r') as f:
    header       = f.readline().strip('\n')
    column_names = re.findall('\[\d+\s+(.*?)\]',header)

for column in column_names:
    print(column)
