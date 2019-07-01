import os
import re
import pandas as pd
import numpy as np
import sys

def load_file(filename):
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


def load_files(filenames):

    # Check that all files have identical headers (i.e. same number of columns)
    nfiles = len(filenames)
    for i in range(nfiles):
        filename = filenames[i]
        with open(filename,'r') as file:
            if i==0:
                prev_header = file.readline()
                column_labels = get_column_names(filename)
            else:
                header = file.readline()
                if header == prev_header:
                    prev_header = header
                else:
                    raise SystemExit('Error: Column labels of file: "'+filename+'" are different to the rest')

    for i in range(nfiles-1):
        data      = np.array(load_file(filenames[i]  ))
        data_next = np.array(load_file(filenames[i+1]))

        # Find the index at which the current file begins to overlap with the next file
        for ii in range(len(data[:,0])):
            t         = data[ii,0]
            tnextfile = data_next[0,0]
            if t > tnextfile:
                break

        if i==0:
            data_combined = np.copy(data[:ii-1,:])
        else:
            data_combined = np.concatenate((data_combined[:,:],data[:ii-1,:]))

    if nfiles>1:
        # Add all of the last file
        data_combined = np.concatenate((data_combined[:,:],data_next[:,:]))

    else:
        # Just load the contents of the 1 file
        data_combined = np.array(load_file(filenames[0]))

    return pd.DataFrame(data_combined,columns=column_labels)


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
