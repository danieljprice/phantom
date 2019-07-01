#!/usr/bin/env python
import argparse
import os
import sys
import re
import evfiles as ev

if sys.version_info[0] < 3:
    raise SystemExit('ERROR: you must use Python 3')

description="""
Use this script to print the column names from an ev file
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('evfile', nargs=1, metavar='file', help='phantom ev file')
args = parser.parse_args()
file = args.evfile[0]

column_names = ev.get_column_names(file)

for column in column_names:
    print(column)
