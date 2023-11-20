# -*- coding: utf-8 -*-

"""
Module to read inputs/outputs from Daniel Price's Phantom
Highly experimental, only tested with recent dumps on a 64 bit machine with no MHD

Author: Lionel Siess & Stéven Toupin
Last tested: 07/07/2021
"""

from numpy import fromfile, frombuffer, memmap, dtype, asscalar

### The following routines are used to extract data from the binary file

# Read data from the binary file
def read_binary(f, count, type):
  try:
    dt = dtype(type)
    size = dt.itemsize*count
    buffer = f.read(size)
    data = frombuffer(buffer, count=count, dtype=type)
  except IOError:
    data = fromfile(f, count=count, dtype=type)
  return data

# Read a fortran record (size, data, size)
def read_fortran_record(f, type, memorymap=False):
  try:
    size1 = int(read_binary(f, 1, 'i4'))
  except IOError:
    size1 = int(read_binary(f, 1, 'i4'))
  #print (size1)
  itemsize = dtype(type).itemsize
  n = int(size1/itemsize)
  if memorymap:
      try:
        pstart = int(f.tell())
        data = memmap(f, dtype=type, mode='r', shape=(n), offset=f.tell())
        f.seek(pstart+size1)
      except OverflowError:
        f.seek(pstart)
        data = read_binary(f, n, type)
        print('Warning: unable to mmap (overflow) (chunk of size %d from offset %d to %d)' % (size1, pstart, pstart+size1))
  else:
      data = read_binary(f, n, type)
  size2 = int(read_binary(f, 1, 'i4'))
  if size1 != size2:
    print('Warning: inconsistent record sizes: ', size1, size2)
  return data

# Read string
def read_string(f):
  line=''.join(map(bytes.decode,read_fortran_record(f, 'S1')))
  return line

# Read a list of names (record with the number of names, followed by a string containing the space-tabbed names
def read_names(f):
  n, = read_fortran_record(f, 'i4')
  if n>0:
    s = read_string(f)
    l = int(len(s)/n)
    names = [s[i:i+l].strip() for i in range(0,len(s),l)] # Split the string
  else:
    names = []
  return names

# Read a list of names with the associated values
# If several values have the same name, they are stored in an array
def read_fortran_record_with_names(f, type):
  from collections import Counter
  names = read_names(f)
  count = Counter(names)
  out = dict()
  if len(names)>0:
    numbers = read_fortran_record(f, type)
    i = 0
    while i<len(names):
      name = names[i]
      c = count[name]
      l = numbers[i:i+c]
      if c==1: l = asscalar(l)
      out[name] = l
      i += c
  return out

def read_dump(filename, memorymap=False):
  """ Read a phantom dump

  Input: the (full-)dump file to read
  Output: a dictionnary containing all the data found in the file
  """

  # Open the file
  if filename.endswith('.bz2'):
    from bz2 import BZ2File
    f = BZ2File(filename, 'rb')
  else:
    f = open(filename, 'rb')

  # Initialize the output dictionnary
  dump = {'filename': filename}

  # Weird record
  read_fortran_record(f, 'i4') # ?

  # Auto description of the file
  dump['desc'] = read_string(f).strip()
  if dump['desc'][0] == 'F':
    variable_type = 'f8'
  else:
    variable_type = 'f4'

  # Number of particles, twice
  dump['part_numbers'] = read_fortran_record_with_names(f, 'i4')
  read_fortran_record_with_names(f, 'i4') # ?
  read_fortran_record_with_names(f, 'i4') # ?
  read_fortran_record_with_names(f, 'i4') # ?
  dump['part_numbers_long'] = read_fortran_record_with_names(f, 'i8')

  # Various quantities (time, gamma, hfact etc.)
  dump['quantities'] = read_fortran_record_with_names(f, variable_type)

  read_fortran_record_with_names(f, variable_type) # ?

  # Units
  dump['units'] = read_fortran_record_with_names(f, 'f8')

  # Number of blocks
  nblocks = int(read_fortran_record(f, 'i4'))

  # Read block sizes
  blocks = []
  for iblock in range(nblocks):
    block = dict()
    n = read_fortran_record(f, 'i4')
    #print('reading block',iblock,' n=',n)
    block['dim'] = n[0]
    block['nint'] = n[1:5]
    block['nintdb'] = n[6]
    block['ndouble'] = n[7]
    block['nsingle'] = n[8]
    blocks.append(block)

  # Read blocks
  for b in blocks:
    b['data'] = dict()
    #print(b['nintdb'],'x',b['nint'],b['ndouble'],b['nsingle'])
    for i in range(b['nintdb']):
      column_name = read_string(f).strip()
      column_data = read_fortran_record(f, 'i8', memorymap)
      #print('[nintdb ] ',i,column_name,type(column_data), variable_type)
      b['data'][column_name] = column_data
    for i in range(b['ndouble']):
      column_name = read_string(f).strip()
      column_data = read_fortran_record(f, 'f8', memorymap)
      #print('[ndouble] ',i,column_name,type(column_data), variable_type)
      b['data'][column_name] = column_data
    for i in range(b['nsingle']):
      column_name = read_string(f).strip()
      column_data = read_fortran_record(f, 'f4', memorymap)
      #print('[nsingle] ',i,column_name,type(column_data), variable_type)
      b['data'][column_name] = column_data
  dump['blocks'] = blocks

  # Close the file
  f.close()
  return dump


def write_hdf5(dump, outfile, compression="gzip", compression_opts=9):
  import h5py
  from numpy import float32, float64, ndarray
  """ Write an hdf5 file from a dump dictionnary
  """
  def write_dict_to_group(d, group, compression, compression_opts):
    """ Recursive function to write a tree branch of the dictionnary into the hdf5 file """
    for k in d.keys():
      #print('%s: %s' % (k, type(d[k])))
      if type(d[k]) is dict:
        subgroup = group.create_group(k)
        write_dict_to_group(d[k], subgroup, compression, compression_opts)
      if type(d[k]) is list:
        subgroup = group.create_group(k)
        for i,e in enumerate(d[k]):
          if type(d[k][i]) is dict:
            subsubgroup = subgroup.create_group(str(i))
            write_dict_to_group(d[k][i], subsubgroup, compression, compression_opts)
      if type(d[k]) in [bool,float,float64,float32,str,int]:
        group.attrs[k] = d[k]
      if isinstance(d[k], ndarray):
        group.create_dataset(k, data=d[k], compression=compression, compression_opts=compression_opts)

  f = h5py.File(outfile,'w')
  write_dict_to_group(dump, f, compression, compression_opts)
  f.close()


def read_hdf5(hdf5file):
  import h5py
  """ Read a dump stored in an hdf5 file
  """
  def read_group(f):
    if len(f.keys())>0 and f.keys()[0]=='0':
      dump = []
      for k in f.keys():
        d = read_group(f[k])
        dump.append(d)
    else:
      dump = dict()
      for k in f.keys():
        if type(f[k])==h5py._hl.group.Group:
          dump[k] = read_group(f[k])
        if type(f[k])==h5py._hl.dataset.Dataset:
          dump[k] = f[k]
      for k in f.attrs.keys():
        dump[k] = f.attrs[k]
    return dump
  f = h5py.File(hdf5file,'r')
  dump = read_group(f)
  return dump

def read_infile(infile):
  """ Read the .in  and .setup files
  """
  data = dict()
  suffix=['.in','.setup']
  for suff in suffix:
      for line in open(infile+suff, 'r'):
        line = line.strip()
        if line.startswith('#') or len(line)==0:
          continue
        line = line.split('!')[0]
        namevalue = line.split('=')
        if len(namevalue)==2:
          name, value = namevalue
          name = name.strip()
          value = value.strip()
          try:
            if '.' in value:
              value = float(value)
            else:
              value = int(value)
          except ValueError:
            pass
          if name in data.keys():
            if not isinstance(data[name], list):
              data[name] = [data[name]]
            data[name].append(value)
          else:
            data[name] = value
  return data

def read_ev_file(filename):
  """ Read a .ev file
  """
  from re import findall, sub
  from numpy import loadtxt
  if filename.endswith('.bz2'):
    from bz2 import BZ2File
    f = BZ2File(filename, 'r')
  else:
    f = open(filename, 'r')
  header = f.readline()
  column_names = findall(r'\[\d+\s+([^\]]+)\]', header)
  lines = f.readlines()
  f.close()
  try:
    data = loadtxt(lines,unpack=True)
  except ValueError:
    # 1.7976931349+308 floats
    lines = [sub(r'(\d)\+(\d+)\s', r'\1E+\2 ', l) for l in lines] # replace with 1.7976931349E+308
    data = loadtxt(lines,unpack=True)
  output = dict()
  for column_number, column_name in enumerate(column_names):
    output[column_name] = data[column_number]
  return output

def concatdicts(A, B):
  """ Merge two python dictionnaries
  """
  for k in B.keys():
    if k in A.keys():
      A[k] = concatenate([A[k], B[k]])
    else:
      A[k] = B[k]

def sourcemorerecent(sourcefile, targetfile):
  """ Returns True if sourcefile is more recent than targetfile
  """
  from os.path import isfile, getmtime
  if isfile(targetfile):
    time_source = getmtime(sourcefile)
    time_target = getmtime(targetfile)
    return time_source > time_target
  else:
    return True

def read_ev_files(location):
  """ Read all ev files in 'location' directory and merge their contents
Stops if one file is older than the previous one
  """
  from os.path import getmtime, join, isdir
  from os import listdir
  if isdir(location):
    files = [f for f in listdir(location) if f.endswith('.ev')]
    f = join(location,files[0])
    output = read_ev_file(f)
    time = getmtime(f)
    for filename in files[1:]:
      timeprev = time
      f = join(location,filename)
      time = getmtime(f)
      if (time > timeprev):
        data = read_ev_file(f)
        concatdicts(output, data)
  return output

def find_prefix(location, infile=None):
  """ Find the prefix of .in and phantom dump files
  """
  class NoInfile(Exception): pass
  class SeveralInfiles(Exception): pass
  class InfileNotFound(Exception): pass
  from os import listdir
  infiles = [f for f in listdir(location) if f.endswith('.in')]
  if len(infiles)==0:
    raise NoInfile
  if infile:
    if not infile in infiles:
      raise InfileNotFound
  else:
    if len(infiles) > 1:
      print('Infiles found:', infiles)
      raise SeveralInfiles
    else:
      infile = infiles[0]
  prefix = infile.rstrip('.in')
  return prefix

def dump_file_list(location, fulldumps=False, infile=None):
  """ Get the list of dump file names
  """
  from os import listdir
  prefix = find_prefix(location, infile)
  dumpfiles = sorted([f for f in listdir(location) if f.startswith(prefix+'_') and not f.endswith('.ascii')])
  if fulldumps:
    numbers = [int(s[len(prefix)+1:len(prefix)+6]) for s in dumpfiles]
    nfulldump = read_infile(prefix+'.in')['nfulldump']
    fulldumps = [dumpfiles[i] for i in range(len(numbers)) if numbers[i]%nfulldump==0]
    dumpfiles = fulldumps
  return dumpfiles

def get_first_dump(location, fulldumps=False, infile=None):
  """ Get the name of the first dump
  """
  return dump_file_list(location, fulldumps, infile)[0]

def get_last_dump(location, fulldumps=False, infile=None):
  """ Get the name of the last dump
  """
  return dump_file_list(location, fulldumps, infile)[-1]

def get_units(location, infile=None):
  """ Get the units from the first dump file
  """
  data = read_dump(get_first_dump(location, False, infile))
  return data['units']

def read_custom_dump(dumpfile, npart, columns):
  """ Read a custom dump (used in gailwind)
  """
  if dumpfile.endswith('.bz2'):
    from bz2 import BZ2File
    f = BZ2File(dumpfile, 'rb')
  else:
    f = open(dumpfile, 'rb')
  ncolumns = len(columns)
  time = read_fortran_record(f, 'f8')[0]
  m = read_fortran_record(f, 'f8').reshape([-1,ncolumns])
  data = dict()
  for i,column in enumerate(columns):
    data[column] = m[0:npart,i]
  f.close()
  return data
