#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt


def extract_headers_bet_spaces(filename):
    f = open(filename,"r")
    raw_data = f.read().split('\n')
    f.close()
    Row1, Headers = raw_data[0].strip(), []
    while len(Row1)>1:
        h = Row1.find(' ')
        Headers.append(Row1[:h])
        Row1 = Row1[h:].strip() + ' '
    return Headers

def extract_data_columns_bet_spaces(filename):
    f = open(filename,"r")
    raw_data = f.read().split('\n')
    f.close()
    Columns, Headers = [], extract_headers_bet_spaces(filename)
    for i in Headers:
        Columns.append([])
    for row in raw_data[1:-1]:
        for col in Columns:
         h = row.find(' ')
         col.append(float(row[:h]))
         row = row[h:].strip() + ' '
    return Columns

def extract_data_columns(R_data,ncols,sep):
    Columns = []
    for i in range(ncols):
        Columns.append([])  
    for col in R_data[1:-1]:
        split_col = col.split(sep)
        for j in enumerate(Columns):
            try:
                j[1].append(float(split_col[j[0]]))
            except:
                j[1].append(split_col[j[0]])
    return Columns


##### For general datafiles with normal separation (variable 'sep') #####




def extract_data(filename,sep=' ',pheaders=True):
    f = open(filename,"r")
    raw_data = f.read().split('\n')
    if sep == ' ':
        headers = extract_headers_bet_spaces(filename)
        columns = extract_data_columns_bet_spaces(filename)
 #############################################################       
    else:
        headers = raw_data[0].split(sep)
        columns = extract_data_columns(raw_data,len(headers),sep)
    if pheaders == True:
        print(headers) 
    formatted_columns = []
    for col in enumerate(columns):
        formatted_columns.append(np.array(col[1]))    

    Data = {}

    for h,c in zip(headers, formatted_columns):
        Data.update({h:c})
    return Data


##### For phantom .ev files, headers type ===>  [1 XX]  [2 YY]  [3 ZZ]   #####

def phantom_evdata(filename,pheaders=True):
    f = open(filename,"r")
    raw_data = f.read().split('\n')
    f.close()
    Row1 = raw_data[0]
    ncols = len(Row1.split(']'))-1
    headers, columns,  = [], []
    for i in Row1.split("]")[:-1]:
        columns.append([])
        headers.append(i.strip('#').strip().strip("[").strip().lstrip('1234567890').strip())



#    l_side, r_side = Row1.find('[')+3, Row1.find(']')-1 
#    width_header = Row1.find('[',Row1.find(']')) - Row1.find('[')
#    headers, columns,  = [], []
#    for i in range(ncols):
#        headers.append(Row1[l_side+i*width_header:r_side+i*width_header+1].strip())
#        columns.append([])
    if pheaders==True:
        print(headers)

    for i in raw_data[1:]:
        if i.strip()=='':
            continue
        S = i
        for j in range(ncols):
            S = S.lstrip()
            if j < ncols-1:
                try:
                    columns[j].append(float(S[:S.find(' ')]))
                except:
                    columns[j].append(S[:S.find(' ')])
            else:
                try:
                    columns[j].append(float(S))
                except:
                    columns[j].append(S)
            S = S[S.find(' '):]

    formatted_columns = []
    for col in enumerate(columns):
        formatted_columns.append(np.array(col[1])) 

    Data = {}

    for h,c in zip(headers, formatted_columns):
        Data.update({h:c})
    return Data

    

##### This will be a new module for a personalized default plot style  #####

def plot_format(xlab,ylab, labeling=False):
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    plt.xlabel(r'$' + xlab + '$', fontsize='x-large')
    plt.ylabel(r'$' + ylab + '$', fontsize='x-large')
    plt.tick_params(labelsize='15')
    if labeling == True:
        plt.legend(fontsize=15)
    elif labeling == False:
        pass
    else: 
        print('Option not valid for labeling. Set as True for show.')


#Conversion of units from Phantom to cgs, day, year... 

class constants:
    mass = 1.989E33 
    time = 1.594E3 
    dist = 6.96E10
    vel = dist/time
    dens = mass/dist**3
    spangmom = dist**2/time
    spener = (dist/time)**2
    ener = mass*spener
    angmom = mass*spangmom
    pressure = ener/dist**3 
    yr = time/(24*3600*365)
    day = time/(24*3600)
    hr = time/(3600)

    def __init__(self,mass=mass,time=time,dist=dist,yr=yr,day=day,
                    spangmom=spangmom,ener=ener,spener=spener,vel=vel,
                    angmom=angmom,dens=dens, pressure=pressure):
        """Phantom units in cgs"""
        self.G = G 
        self.mass = mass 
        self.time = time 
        self.dist = dist
        self.vel = vel
        self.dens = dens
        self.spangmom = spangmom
        self.spener = spener
        self.angmom = angmom
        self.ener = ener
        self.pressure = pressure
        self.yr = yr 
        self.day = day 
        self.hr = hr
    
        
