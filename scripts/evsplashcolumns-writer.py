"""
Written by David Liptai, 2017.

This script takes a PHANTOM .ev file as an input and uses
the column labels in the header to create 'evsplash.columns',
which will be read by SPLASH when plotting the ev file.

This should save you the trouble of checking to see which
columns correspond to which variable when plotting, or having
to write your own columns file manually.

NOTE: Please use a verison of Python3.

"""

import os
import sys
import subprocess
from numpy import savetxt

def BASH(command):
    return subprocess.check_output(command,shell=True).decode().strip()

error_message = '\n Please provide a .ev file with the column labels you need as an input argument \n \nFAIL'

print('START')

try:
    filename = sys.argv[1]
except:
    quit("\n No valid input file given" + error_message)

if not os.path.isfile(filename): quit("\n Could not find '" + filename + "'" + error_message)

command = 'head -1 ' + filename
string = BASH(command)
string = string.lstrip('# [').rstrip(']').split(']   [')
string = [i[2:].strip() for i in string]

out_fname = 'evsplash.columns'
savetxt(out_fname,string,fmt='%s')
print("\n Column labels from: '" + command + "'")
print("         Written to: '" + out_fname + "'")
print('\nCOMPLETE')
