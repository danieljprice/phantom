import re
import struct
import numpy as np
import matplotlib.pyplot as plt 

#---------------INITIALISATIONS---------------
#Use the same mass units as defaults (1e9 vs 1e10)
kpc_cgs= 3.08567758e21
G_cgs  = 6.67e-8
Mo_cgs = 1.99e33
umass_GizToGas = 1.  #1e9Mo
umass = 1e10 #* umass_GizToGas
udist = 1.0  #kpc
uvel  = np.sqrt( G_cgs * umass * Mo_cgs / (udist * kpc_cgs) )/1e5
udens = umass * Mo_cgs / (udist * kpc_cgs)**3.
utime = np.sqrt(1./(udens * G_cgs))
print '<<<<<<<<<UNITS>>>>>>>>>'
print 'udist =',udist,'[kpc]'
print 'umass =',umass,'[Mo]'
print 'uvel  =',uvel,'[km/s]'
print 'udens =',udens,'[g/cm3]'
print 'utime =',utime,'[s]'
sec2myr = 60.*60.*24.*365.*1e6




def read_magalie(filenom,gengas,gasfrac):


  #CONVERSION UNITS FOR MAGALIE:
  kms   = 1e3
  kpc   = 3.086e19
  Mo    = 1.99e30      
  #Convert to SI cos' I keep borking this up, sorry folks.

  #######IMPORTANT: these values are used to re-sale the disc
  uvel_n  = 200.   #uvel for GIZMO is in km/s
  udist_n = 3.5    #udist for GIZMO is in kpc, was 2.45     
  print 'WARNING:----------------------------------------'
  print 'The code requires some distance and vel scales'
  print 'I advise you look at the plotted rotation curves'
  print 'to ensure this looks like what you want.'
  print 'uvel  = ',uvel_n
  print 'udist = ',udist_n
  print '            <HIT RETURN TO COMPLY>'
  print 'WARNING:----------------------------------------'
  raw_input()


  umass_n = ((uvel_n*kms*uvel_n*kms) * (udist_n*kpc)) /Mo/ 6.67e-11
  umass_n = umass_n/1e10  #umass for GIZMO in 1e10Mo
  #Rountine to read in what nemo-magalie dumps out
  print 'Reading in data file:',filenom
  print 'WARNING: nbulge,ngas,ndisc,nbulge are hard-wired for now.'
  print '         These can be read in, but need the magalie param file.'

  print 'magalie ascii dumps are in the format:'
  print '  header: ntot, ntypes ,0'
  print '  data  : m, x, y, z, vx, vy, vz'
  print 'where particles are ordered by: disc->bulge->halo.'

  print 'First open up the magalie.in file for the IC parameters'
  ReadPar    = open('magalie.in')
  for i, line in enumerate(ReadPar):
    if i == 2:
      ndisc  = int(line.split()[0])
    if i == 24:
      nbulge = int(line.split()[0])
    if i == 36:
      nhalo  = int(line.split()[0])
    elif i>36:
      break
      ReadPar.close

  print 'Particle counts:'
  print '  ndisc =',ndisc
  print '  nbulge=',nbulge
  print '  nhalo =',nhalo
  ntot = ndisc+nbulge+nhalo

  pd=np.zeros((ndisc,3))
  pb=np.zeros((nbulge,3))
  ph=np.zeros((nhalo,3))
  vd=np.zeros((ndisc,3))
  vb=np.zeros((nbulge,3))
  vh=np.zeros((nhalo,3))

  Readin    = open(filenom)
  partCount = 0
  for i, line in enumerate(Readin):
    if i == 0:
        n_dbh=int(line.split()[0])
        ntype=int(line.split()[1])
        print 'No. of stars+bulge+DM partciles:',n_dbh
        if n_dbh!=(ndisc+nbulge+nhalo):
          print 'Particle numbers dont add up...'
          exit()
        if ntype!=3:
          print 'You have are missing a particle type.'
          exit()
    else:
        partCount += 1
        mi = float(line.split()[0])
        xi = float(line.split()[1])
        yi = float(line.split()[2])
        zi = float(line.split()[3])
        vxi= - float(line.split()[4])
        vyi= - float(line.split()[5])
        vzi= float(line.split()[6])
        if partCount<ndisc:
          #disc
          md  =  mi
          pd[partCount,0]=xi
          pd[partCount,1]=yi
          pd[partCount,2]=zi
          vd[partCount,0]=vxi
          vd[partCount,1]=vyi
          vd[partCount,2]=vzi
        elif partCount<ndisc+nbulge:
          #bulge
          mb  =  mi
          pb[partCount-ndisc,0]=xi
          pb[partCount-ndisc,1]=yi
          pb[partCount-ndisc,2]=zi
          vb[partCount-ndisc,0]=vxi
          vb[partCount-ndisc,1]=vyi
          vb[partCount-ndisc,2]=vzi
        elif partCount<ndisc+nbulge+nhalo:
          #halo
          mh  =  mi
          ph[partCount-ndisc-nbulge,0]=xi
          ph[partCount-ndisc-nbulge,1]=yi
          ph[partCount-ndisc-nbulge,2]=zi
          vh[partCount-ndisc-nbulge,0]=vxi
          vh[partCount-ndisc-nbulge,1]=vyi
          vh[partCount-ndisc-nbulge,2]=vzi
  Readin.close()
  print 'NOTE: we have flipped the velocities so the disc rotates clockwise.'
  print 'Masses:',md,mb,mh

  print 'Convert units form magalie->GADGET/GIZMO/GASOLINE/PHANTOM:'
  print 'udist =',udist_n
  print 'uvel  =',uvel_n
  print 'umass =',umass_n
  pd = pd*udist_n
  pb = pb*udist_n
  ph = ph*udist_n
  vd = vd*uvel_n
  vb = vb*uvel_n
  vh = vh*uvel_n
  md = md*umass_n
  mb = mb*umass_n
  mh = mh*umass_n

  gasoff  = 90.
  print 'Making a gas array by rotating the disc stars by',gasoff,'[deg]'
  print ' *ngas must be ndisc for now'
  print 'Frist assigning mass of gas as',gasfrac,'that of disc stars'
  mg   = md * float(gasfrac)
  ngas = ndisc 
  gasoff = gasoff  *  3.4159/180.
  pg =np.zeros((ngas,3))
  vg =np.zeros((ngas,3))
  rho=np.zeros(ngas)
  EN =np.zeros(ngas)
  i=0
  while i<ngas:
    pg[i,0] =  pd[i,0]*np.cos(gasoff) - pd[i,1]*np.sin(gasoff)
    pg[i,1] =  pd[i,0]*np.sin(gasoff) + pd[i,1]*np.cos(gasoff)
    pg[i,2] =  pd[i,2]
    vg[i,0] =  vd[i,0]*np.cos(gasoff) - vd[i,1]*np.sin(gasoff)
    vg[i,1] =  vd[i,0]*np.sin(gasoff) + vd[i,1]*np.cos(gasoff)
    vg[i,2] =  vd[i,2]
    i+=1
  print 'Gas array built'
  
  return ngas,nhalo,ndisc,nbulge,mg,md,mb,mh,pg,pd,pb,ph,vg,vd,vb,vh




print '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
print '<<<Here to write your dumps for you >>>'
print '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'

wantgas = raw_input('Feeling gassy? [y/y]')

if wantgas == 'y':
  print 'Setting our own gas disc from disc stars.'
  gasfrac = raw_input('Enter gas to stellar disc mass ratio [0.05-0.5 advised]')
  if gasfrac<=0.:
    print 'Get out.'

ngas,nhalo,ndisc,nbulge,mg,md,mb,mh,pg,pd,pb,ph,vg,vd,vb,vh = read_magalie('m.dat',wantgas,gasfrac)
print 'NEMO read-in done'

######PLOT SOME THINGS....
plt.figure(1)
plt.plot(pg[:,0],pg[:,1],'r.',markersize=0.2,lw=0)
plt.axes().set_aspect('equal')
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.savefig('gas')

plt.figure(2)
plt.plot(pb[:,0],pb[:,1],'g.',markersize=0.5,lw=0)
plt.axes().set_aspect('equal')
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.savefig('bulge')

plt.figure(3)
plt.plot(ph[:,0],ph[:,1],'b.',markersize=0.5,lw=0)
plt.axes().set_aspect('equal')
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.savefig('halo')

plt.figure(4)
plt.plot(pd[:,0],pd[:,1],'k.',markersize=0.5,lw=0)
plt.axes().set_aspect('equal')
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.savefig('disc')

plt.figure(5)
plt.plot(ph[:,0],ph[:,2],'b.',markersize=0.5,lw=0)
plt.plot(pg[:,0],pg[:,2],'r.',markersize=0.5,lw=0)
plt.plot(pd[:,0],pd[:,2],'k.',markersize=0.5,lw=0)
plt.plot(pb[:,0],pb[:,2],'g.',markersize=0.5,lw=0)
plt.axes().set_aspect('equal')
plt.savefig('side')

plt.figure(6)
plt.plot(ph[:,0],ph[:,1],'b.',markersize=0.5,lw=0)
plt.plot(pg[:,0],pg[:,1],'r.',markersize=0.5,lw=0)
plt.plot(pd[:,0],pd[:,1],'k.',markersize=0.5,lw=0)
plt.plot(pb[:,0],pb[:,1],'g.',markersize=0.5,lw=0)
plt.axes().set_aspect('equal')
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.savefig('top')

plt.figure(7)
plt.plot(np.sqrt(pd[:,0]**2+pd[:,1]**2),np.sqrt(vd[:,0]**2+vd[:,1]**2),'k.',markersize=0.5,lw=0)
plt.plot(np.sqrt(ph[:,0]**2+ph[:,1]**2),np.sqrt(vh[:,0]**2+vh[:,1]**2),'b.',markersize=0.5,lw=0)
plt.plot(np.sqrt(pb[:,0]**2+pb[:,1]**2),np.sqrt(vb[:,0]**2+vb[:,1]**2),'g.',markersize=0.5,lw=0)
plt.savefig('RC_starshalo.png')
plt.clf()
plt.plot(np.sqrt(pd[:,0]**2+pd[:,1]**2),np.sqrt(vd[:,0]**2+vd[:,1]**2),'k.',markersize=0.5,lw=0)
plt.plot(np.sqrt(pb[:,0]**2+pb[:,1]**2),np.sqrt(vb[:,0]**2+vb[:,1]**2),'g.',markersize=0.5,lw=0)
plt.savefig('RC_stars.png')
plt.xlim(0,20)
plt.ylim(0,300)
plt.clf()


print '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
print '<<<Here to write this sheet for you>>>'
print '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='

ntot  = ndisc+nbulge+ngas+nhalo

print 'Converting arrays'
pg = pg /udist
pd = pd /udist
pb = pb /udist
ph = ph /udist
m_g  = mg/umass_GizToGas
m_d  = md/umass_GizToGas
m_b  = mb/umass_GizToGas
m_h  = mh/umass_GizToGas
vg = vg /uvel
vd = vd /uvel
vb = vb /uvel
vh = vh /uvel
print 'Converting to SI units (yes this is redundant but cba).'
udSItips = 3.0857e19 # (kpc)
uvSItips = uvel * 1.e3  # (207km/s)
umSItips = 1.e10 * 1.99e30 # (1e10Mo)
pg = pg*udSItips
pd = pd*udSItips
pb = pb*udSItips
ph = ph*udSItips
m_g  = m_g*umSItips
m_d  = m_d*umSItips
m_b  = m_b*umSItips
m_h  = m_h*umSItips
vg = vg*uvSItips
vd = vd*uvSItips
vb = vb*uvSItips
vh = vh*uvSItips

print 'Writing out x4 asciifiles'
phased = 10.
phaseg = 0.
phaseb = 10.
phaseh = 10.

outg=open('asciifile_G','w')
iout=0
while iout<ngas:
  outg.write(str(pg[iout,0])+' '+str(pg[iout,1])+' '+str(pg[iout,2])+' '+str(m_g)+' '+str(vg[iout,0])+' '+str(vg[iout,1])+' '+str(vg[iout,2])+' '+str(phaseg)+'\n')
  iout+=1
outg.close()

outd=open('asciifile_D','w')
iout=0
while iout<ndisc:
  outd.write(str(pd[iout,0])+' '+str(pd[iout,1])+' '+str(pd[iout,2])+' '+str(m_d)+' '+str(vd[iout,0])+' '+str(vd[iout,1])+' '+str(vd[iout,2])+' '+str(phased)+'\n')
  iout+=1
outd.close()

outb=open('asciifile_B','w')
iout=0
while iout<nbulge:
  outb.write(str(pb[iout,0])+' '+str(pb[iout,1])+' '+str(pb[iout,2])+' '+str(m_b)+' '+str(vb[iout,0])+' '+str(vb[iout,1])+' '+str(vb[iout,2])+' '+str(phaseb)+'\n')
  iout+=1
outb.close()

outh=open('asciifile_H','w')
iout=0
while iout<nhalo:
  outh.write(str(ph[iout,0])+' '+str(ph[iout,1])+' '+str(ph[iout,2])+' '+str(m_h)+' '+str(vh[iout,0])+' '+str(vh[iout,1])+' '+str(vh[iout,2])+' '+str(phaseh)+'\n')
  iout+=1
outh.close()


print 'Last step, so write out the file needed for phantom to setup.'
outnom = 'galsetic.txt'
outo=open(outnom,'w')
iout=0
outo.write('Ngas'+' '+str(ngas)+'\n')
outo.write('Nstar'+' '+str(ndisc)+'\n')
outo.write('Nbulge'+' '+str(nbulge)+'\n')
outo.write('Ndark'+' '+str(nhalo)+'\n')
outo.close()
print '<<<<<<<<<<<<<<DONE>>>>>>>>>>>>>>>'
print 'Now place',outnom,'and the asciifiles in the directory of phantomsetup.'

print 'DONE'

