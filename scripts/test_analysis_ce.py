#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# unit tests for the analysis_common_envelope module
# written by Miguel Gonzalez-Bolivar, Feb 2023
#
import os
import dataread as dr

sep_file = 'separation_vs_time.ev'
bound_file = 'boundunbound_vs_time.ev'
energy_file = 'energy.ev'
evfiles = [sep_file, bound_file, energy_file]

#Delete all files from previous tests
#os.system('rm *.ev *txt')

#Compile, and run phantomanalysis
#os.system('sh ~/phantom/scripts/test_analysis_ce.sh')

#Check that .ev files were created
def check_file_existence(evfile):
    os.system('ls ' + evfile +' > check_file.txt')
    check_file = open('check_file.txt')

    if len(check_file.readline().strip('\n')) != len(evfile):
        print('Failed: ' + evfile + ' not found')
    else:
        print('Passed: ' + evfile + ' created')
    check_file.close()

#Check that .ev for NaN values in .ev files
def check_file_nan(evfile):
    os.system('grep -nir "nan" ' + evfile + ' > check_nan.txt')
    check_nan = open('check_nan.txt')
    if len(check_nan.readline().strip('\n')) != 0:
        print('Failed: NaN values found in ' + evfile)
    else:
        print('Passed: No NaN values found in ' + evfile)

#Check that values of a given column "col" in evfile are all positive or negative
#(pot energy has to be always negative, kin energy always positive, ...)
def check_sign_values(evfile,col,crit='pos'):
    values = dr.phantom_evdata(evfile,pheaders=False)[col]
    if (crit == 'pos') or (crit == 'p'):
        if any(values>0) == True:
            print('Passed: ' + col + ' in ' + evfile + ' has all values positive')
        else:
            print('Failed: ' + col + ' in ' + evfile + ' has negative values')
    if (crit == 'ned') or (crit == 'n'):
        if any(values<0) == True:
            print('Passed: ' + col + ' in ' + evfile + ' has all values negative')
        else:
            print('Failed: ' + col + ' in ' + evfile + ' has positive values')
    print(' ')


#Run existence test for all evfiles
for evf in evfiles:
    check_file_existence(evf)
    print(' ')

#Check that looks for NaN values in evfiles
for evf in evfiles:
    check_file_nan(evf)
    print(' ')

#Checks obvious signs in certain quantities
check_sign_values(energy_file,'kin energy')
check_sign_values(energy_file,'pot energy',crit='n')
check_sign_values(energy_file,'therm energy')
