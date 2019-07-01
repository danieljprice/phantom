import os
import re
import pandas as pd
import numpy as np
import sys

def load(filename):
    # Check that file exist
    if not os.path.isfile(filename):
        raise SystemExit('Error: the file "'+filename+'" does not exist or is not a file')

    column_names = get_column_names(filename)

    return pd.read_csv(filename,skiprows           = 0,
                                dtype              = np.float64,
                                delim_whitespace   = True,
                                skip_blank_lines   = True,
                                comment            = '#',
                                header             = None,
                                names              = column_names)


def get_column_names(filename):
    # Read the columns names in the file using regex
    with open(filename,'r') as f:
        header       = f.readline().strip('\n')
        column_names = re.findall('\[\d+\s+(.*?)\]',header)
        return column_names

def printev(data,labels):

    if labels is not None:
        # Construct the header string
        header = '# '
        for i in range(len(labels)):
            label  = labels[i].rjust(12)
            no     = str(i+1).zfill(2)
            header = header + '[' + no + label + ']' + '   '
        header = header[:-1]
    else:
        header = ''

    np.set_printoptions(threshold=np.inf,floatmode='fixed',linewidth=np.inf)
    np.savetxt(sys.stdout.buffer,data,header=header,fmt='%18.10E',comments='')
