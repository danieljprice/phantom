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

data = ev.load_files(args.evfiles)
column_labels = list(data)

ev.printev(np.array(data),column_labels)
