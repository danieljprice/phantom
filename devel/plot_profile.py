# -*- coding: utf-8 -*-
# python3 script
from pylab import *
from readPhantomDump import *
from PhysicalConstantsCGS import k, amu, au
from sys import argv
from os import path

'''
use : python3 ./plot_profile.py wind_00012 wind.in
'''
R_rochelobe = 0. # draw a vertical line at the location of the Roche radius

s=sys.argv
if (len(s) != 3) :
    s=input('specify dump file and input file (e.g.: wind_000100 wind.in) ')
    fdump=s.split()[0]
    fparm=s.split()[1]
else:
    fdump=s[1]
    fparm=s[2]


if not 'dataloaded' in locals(): # Avoid reloading everything in spyder
    # Read the 1D data
    print('Reading 1D data')
    f = open('windprofile1D.dat','r')
    headers = f.readline()
    headers = headers.split()[1:]
    data1D = dict()
    m = loadtxt(f)
    for i,column in enumerate(headers):
      data1D[column] = m[:,i]
    del headers, m, i, column

    # Read the 3D data
    print ("Reading 3D data (Phantom dump)",fdump)
    data3D = read_dump(fdump)
    npart = data3D['part_numbers']['npartoftype'][0]
    #print('Reading 3D data (custom dump).')
    #data3D_custom = read_custom_dump('injectcustomdump.bz2', npart, ['Jstar', 'K0', 'K1', 'K2', 'K3', 'mu'])

    dataloaded = True

print ("Read input parameters",fparm)
wind_param = read_infile(fparm)

Rstar = wind_param['wind_inject_radius']*au
Rstar = au
gamma = 5./3.

# Calculate specific energy for the 1D model
data1D['e'] = k*data1D['T']/((gamma-1.)*data1D['mu']*amu)

# Extract 3D data
udist, umass, utime = [data3D['units'][u] for u in ('udist', 'umass', 'utime')]
hfact = data3D['quantities']['hfact']
mpart = data3D['quantities']['massoftype'][0]*umass
x, y, z, vx, vy, vz, h, u = [array(data3D['blocks'][0]['data'][c],dtype=float)*conv for c, conv in [('x', udist), ('y', udist), ('z', udist), ('vx', udist/utime), ('vy', udist/utime), ('vz', udist/utime), ('h', udist), ('u', udist**2/utime**2)]]
r = sqrt(x**2+y**2+z**2)/au
v = sqrt(vx**2+vy**2+vz**2)
rho = mpart * (hfact/h)**3
del x, y, z, vx, vy, vz, h, hfact, mpart

# Parameters for integration with A&A style
rcParams['font.serif'] = 'Times New Roman'
rcParams['font.sans-serif'] = 'Times New Roman'
rcParams['font.size'] = 8.966
rcParams['axes.titlesize'] = 9.962
rcParams['axes.linewidth'] = 0.398
rcParams['lines.linewidth'] = 0.398
rcParams['xtick.major.width'] = 0.398
rcParams['xtick.minor.width'] = 0.198
rcParams['ytick.major.width'] = 0.398
rcParams['ytick.minor.width'] = 0.198
rcParams['patch.linewidth'] = 0.398
rcParams['legend.fontsize'] = 8.966

# Plot velocity
figure(1) #, figsize=(3.5,2.62), dpi=360)
clf()
plot(data1D['r']/Rstar, data1D['v']/1.e5, color='r', label='1D')
plot(r/Rstar, v/1.e5, '.', color='b', markersize=1.592, rasterized=True, label='3D')
if R_rochelobe > 0:
    axvline(R_rochelobe, ls='dashed', label='Roche lobe', dashes=(3,3))
axis([0., 26., 0., 35.])
xlabel(r'$r$ [R$_*$]')
ylabel(r'$v$ [km s$^{-1}$]')
legend(loc='best')
tight_layout()
savefig('plot_v.pdf', dpi=360)

# Plot density
figure(2) #, figsize=(3.5,2.62), dpi=360)
clf()
semilogy(data1D['r']/Rstar, data1D['rho'], color='r', label='1D')
semilogy(r/Rstar, rho, '.', color='b', markersize=1.592, rasterized=True, label='3D')
if R_rochelobe > 0:
    axvline(R_rochelobe, ls='dashed', label='Roche lobe', dashes=(3,3))
axis([0., 26., 1.e-17, 1.e-12])
xlabel(r'$r$ [R$_*$]')
ylabel(r'$\rho$ [g cm$^{-3}$]')
legend(loc='best')
tight_layout()
savefig('plot_rho.pdf', dpi=360)

# Plot specific energy
figure(3) #, figsize=(3.5,2.62), dpi=360)
clf()
semilogy(data1D['r']/Rstar, data1D['e'], color='r', label='1D')
semilogy(r/Rstar, u, '.', color='b', markersize=1.592, rasterized=True, label='3D')
if R_rochelobe > 0:
    axvline(R_rochelobe, ls='dashed', label='Roche lobe', dashes=(3,3))
axis([0., 26., 10.**7.5, 10.**11.5])
xlabel(r'$r$ [R$_*$]')
ylabel(r'$e$ [erg g$^{-1}$]')
legend(loc='best')
tight_layout()
savefig('plot_e.pdf', dpi=360)
