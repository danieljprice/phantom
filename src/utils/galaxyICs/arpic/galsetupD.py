import matplotlib.pyplot as plt
import numpy as np
import random
import math
from scipy.special import i0,i1,iv,k0,k1,kn
#import asciitable
from scipy.integrate import quad as spiq
from scipy.stats import maxwell as mbd
import sys



def readinoptions():
  print "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  print " __  __ _  _______ _____  _____  _____ "
  print "|  \/  | |/ /  __ \_   _|/ ____|/ ____|"
  print "| \  / | ' /| |  | || | | (___ | |     "
  print "| |\/| |  < | |  | || |  \___ \| |     "
  print "| |  | | . \| |__| || |_ ____) | |____ "
  print "|_|  |_|_|\_\_____/_____|_____/ \_____|"
  print "                                       "
  print "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"



  inputfile = (sys.argv)[1]
  ans1      = (sys.argv)[2]
  answers   = open(inputfile)
  ans=[]
  for line in answers:
    ans.append((line.split())[0])
  ansP      = ans[0]
  Npart_P   = int(ans[1])
  Npart_V   = int(ans[2])
  startogas = float(ans[3])
  mswing    = float(ans[4])
  GasDisc   = float(ans[5])
  
  print '-=-=-=-=-=-=-=-=-=-=-=-=-=-'
  if ans1=='o':
    print 'Plotting RC curve alone'
  elif ans1=='p':
    print 'Drawing positions'
  elif ans1=='v':
    print 'Drawing velocities (must be done after positions)'
  elif ans1=='f':
    print 'Finishing up.'
  else:
    print 'Must pass argument for [o/p/v] at run time.'
    exit()
    
  if ansP=='p':
    print 'Plotting as we go'
  else:
    print 'Not plotting (likely too many paritcles for matplotlib!)'    
  
  print 'Setting:',Npart_P,'disc particle positions'
  print 'Setting:',Npart_V,'disc particle velocities (Npart_P<=Npart_V)'

  if Npart_V>Npart_V:
    print '**ERROR**'
    print 'Npart_V must be <= Npart_P or positions will not be available!'
    exit()

  print 'With a star to gas ratio of:',startogas  
  if startogas>0.5:
    print '**WARNING**'
    print 'Are you sure you want more stars than gas? Useful mainly for testing.'
  print '-=-=-=-=-=-=-=-=-=-=-=-=-=-'

  print 'Using a gas disc mass of:',GasDisc,'x10^9 Mo'

  return ans1,ansP,Npart_P,Npart_V,startogas,mswing,GasDisc
 

"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
"In SI units:"
Mo  = 1.99e30
kpc = 1.0e3 * 3.0856e16
G   = 6.67e-11
km  = 1.0e3
"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
"Variables (set to Baba et al. 2009 rather than Hernquist 1993 values)"
r1  = 0.0 * kpc
r2  = 13. * kpc     #The full disk radius we will use
z1  = 0. 

rscatter = 0.001*kpc  #radius to shuffle particles in Gaussian to take into account finite r array

zfr = 0.40
z2  = +zfr*r2
"halo:"
r200= 122 * kpc
Cnfw= 5.
rh  = r200/Cnfw
Mh  = 6.3e11 * Mo
"disk:"
r0  = 3.0 * kpc
z0  = 0.3 * kpc
Md  = 1.3 * 3.2e10 * Mo
"Orbital Params:"
Qref= 1.0
Rref= 1.5 * r0



def Sofue12():
  f=open('/home/alex/SPHNG_seeding/Bart_codes/rotcurve_data/SofueWebsite/Sofue12/Table_GrandRC.dat.txt')
  r =np.zeros(129)
  v =np.zeros(129)
  ev=np.zeros(129)
  i=0
  for line in f:
    r[i] =line.split()[0]
    v[i] =line.split()[1]
    ev[i]=line.split()[2]
    i+=1
  plt.errorbar(r,v,yerr=ev,xerr=None,fmt='bo',label='Sofue12 (GRC)',markersize=2)

def sech(x):
  s = 1./np.cosh(x)
  return s

def RhoDisk(r,z):
  C = Md/(4.*math.pi*r0*r0*z0)
  rho = C*np.exp(-r/r0)*np.power(sech(z/z0),2.)
  return rho
  
def SigDisk(r):
  C = Md/(2.*math.pi*r0*r0)
  sigma = C*np.exp(-r/r0)
  return sigma

def SigDiskz(z):
  C = Md/(4.*math.pi*r0*r0)   #Not exact constant but it will be integrated out in pdf.
  sigma = C*np.power(sech(z/z0),2.)
  return sigma

def pdfDisk(r):
  NormConst = Normaliser(1)
  N_r   = r*np.exp(-r/r0)
  pdf   = NormConst*N_r
  return pdf

def pdfDiskz(z):
  NormConst = Normaliser(2)
  N_z   = np.power(sech(z/z0),2.)
  pdf   = NormConst*N_z
  return pdf
  
def cdfDisk(r):
  NormConst = Normaliser(1)
  N_rdr   = -r0*np.exp(-r/r0)*(r0+r)
  r   = 0.
  N_rdr0  = -r0*np.exp(-r/r0)*(r0+r)
  cdf = NormConst*(N_rdr-N_rdr0)
  return cdf

def cdfDiskz(z):
  "NOTE: the z-dist doesnt need a *z as r does, because it is already pre one units dist (Sigma was 1/r^2)"
  NormConst = Normaliser(2)
  N_zdz   = z0*np.tanh(z/z0)  #z0*z*np.tanh(z/z0) - z0*z0*np.log(np.cosh(z/z0))
  z   = 0.
  N_zdz0  = z0*np.tanh(z/z0)  #z0*z*np.tanh(z/z0) - z0*z0*np.log(np.cosh(z/z0))
  cdf = NormConst*(N_zdz-N_zdz0)
  return cdf

def Normaliser(flag):
  if flag==1:
    "Radial"
    r=r2
    c2 = -r0*np.exp(-r/r0)*(r0+r)
    r=r1
    c1 = -r0*np.exp(-r/r0)*(r0+r)
    return 1./(c2-c1)
  if flag==2:
    "Vertical"
    z=z2
    c2 = z0*np.tanh(z/z0)  #z0*z*np.tanh(z/z0) - z0*z0*np.log(np.cosh(z/z0))
    z=z1
    c1 = z0*np.tanh(z/z0)  #z0*z*np.tanh(z/z0) - z0*z0*np.log(np.cosh(z/z0))
    return 1./(c2-c1)
  
def randomdist(r,cd,N):
  Rpred=np.zeros(N)
  "Pull N random particles to place"
  j=0
  while j<N:
	z=np.random.rand()
	k=0
	while  (k<len(r)) and (z>=cd[k]):
	  k+=1
	if k==len(r):
	  Rpred[j] = r[-1]
	else:
	  Rpred[j] = r[k]
	j+=1
  return Rpred

def gausfunc(x,mu,sigma):
  y=np.exp(-(x/sigma)**2)/np.sqrt(2.*math.pi*sigma**2)
  return y

def PotDisk(n,r):
  "n: order of differentiation"
  rratio  =0.5*r/r0
  if n==1:
    dPhidr_disk  =(G*Md*r/(2.*r0*r0*r0)) * (i0(rratio)*k0(rratio)-i1(rratio)*k1(rratio))
    return dPhidr_disk
  elif n==2:
    d2Phidr2_disk=(G*Md  /(4.*r0*r0*r0)) * (-rratio*iv(2,rratio)*k1(rratio) + i0(rratio)*(2*k0(rratio)-3.*rratio*k1(rratio)) + i1(rratio)*(3.*rratio*k0(rratio)-2.*k1(rratio)+rratio*kn(2,rratio)) )
    return d2Phidr2_disk
  elif n=='v':
    dPhidr_disk  =(G*Md*r/(2.*r0*r0*r0)) * (i0(rratio)*k0(rratio)-i1(rratio)*k1(rratio))
    Vc  =  np.sqrt(r*dPhidr_disk)
    return Vc
  else:
    Phi_disk     =(-G*Md*r/(2.*r0*r0)) * (i0(rratio)*k1(rratio)-i1(rratio)*k0(rratio))
    return Phi_disk
  
def PotHalo(n,r):
  "n: order of differentiation"
  if n==1:
    dPhidr_halo  =( -np.log(1.+r/rh)/(r*r) + 1./((rh*r)*(1.+r/rh)) ) * -G*Mh/(np.log(1.+Cnfw)-Cnfw/(1.+Cnfw))
    return dPhidr_halo
  elif n==2:
    d2Phidr2_halo=( +2*np.log(1.+r/rh)/(r*r*r) - 1./((rh*rh*r)*(1.+r/rh)*(1.+r/rh)) - 2./((rh*r*r)*(1.+r/rh)) ) * -G*Mh/(np.log(1.+Cnfw)-Cnfw/(1.+Cnfw))
    return d2Phidr2_halo
  elif n=='v':
    dPhidr_halo  =( -np.log(1.+r/rh)/(r*r) + 1./((rh*r)*(1.+r/rh)) ) * -G*Mh/(np.log(1.+Cnfw)-Cnfw/(1.+Cnfw))
    Vc  =  np.sqrt(r*dPhidr_halo)
    return Vc
  else:
    Phi_halo     =(  np.log(1.+r/rh)/r ) * -G*Mh/(np.log(1.+Cnfw)-Cnfw/(1.+Cnfw))
    return Phi_halo

def kappa_freq(r):
  "Epicycle frequency"
  dPhidr_disk  =PotDisk(1,r)
  d2Phidr2_disk=PotDisk(2,r)
  dPhidr_halo  =PotHalo(1,r)
  d2Phidr2_halo=PotHalo(2,r)
  kappa = np.sqrt( (dPhidr_disk+dPhidr_halo)*3./r  +  (d2Phidr2_disk+d2Phidr2_halo)  )
  return kappa

def omega_freq(r):
  "Circular frequency"
  Vcirc = Vc(r)
  omega = Vcirc/r
  return omega

def Vc(r):
  "Circular velocity"
  Vcirc = np.sqrt(r*(PotDisk(1,r) + PotHalo(1,r)))
  return Vcirc
    
def ToomreConstFn(r):
  kappa_Rref = kappa_freq(Rref)
  Vref  = Qref * 3.36 * G * SigDisk(Rref)/kappa_Rref
  Const = (Vref*Vref) * np.exp(Rref/r0)
  return Const



"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
"RUN TIME"
"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
ans1,ansP,Npart_P,Npart_V,StarToGas,mswing,GassDisc=readinoptions()


"How to scale masses for different swing amplifications"
if mswing ==2. :
  print 'Setting mass fractions to give m=2 spiral'
  Md_scale = 1.5   #2.0 was a bit strong.
  Mh_scale = 0.5
  print ' ' 
  print '--WARNING--WARNING--WARNING--WARNING--'
  print '     REMEMBER THE HALO POTENTIAL      '
  print ' WILL NEED TO BE CHANGED ACCORDINGLY  '
  print '--WARNING--WARNING--WARNING--WARNING--'
  print ' ' 
elif mswing ==4. :
  print 'Setting mass fractions to give m=4 spiral'
  Md_scale = 0.65
  Mh_scale = 1.4
  print ' ' 
  print '--WARNING--WARNING--WARNING--WARNING--'
  print '     REMEMBER THE HALO POTENTIAL      '
  print ' WILL NEED TO BE CHANGED ACCORDINGLY  '
  print '--WARNING--WARNING--WARNING--WARNING--'
  print ' ' 
else:
  print 'Setting mass fractions to default (roughly m=3 spiral)'
  Md_scale = 1.0
  Mh_scale = 1.0
Md = Md_scale * Md
Mh = Mh_scale * Mh
Mg = GassDisc * 1.0e9 * Mo

print 'The bar stability parameter;'
print 'eta_b = Vmax / SQRT( G*Md/Rd ), is' 
Vmax = np.max( Vc(np.arange(0.1,13.0,0.1)*kpc) )
print 'Vmax=',Vmax/1000.,'km/s'
etab = Vmax/np.sqrt( G*Md/r0 )
print etab
if (etab<1.1):
  print 'Unstable to bar formation (etabar<1.1)'
else:
  print 'Stable to bar formation (etabar>1.1)'

if ans1=='o':
  print '------------------------------------'
  print 'Plotting Rc curve, no initialisation'
  print '------------------------------------'
  r=np.arange(0.01,13.0,0.01)*kpc
  p1 = PotDisk('v',r)
  p2 = PotHalo('v',r)
  pt = np.sqrt(np.power(p1,2)+np.power(p2,2))
  plt.figure(1)
  plt.plot(r/kpc,p1/km,'b--')
  plt.plot(r/kpc,p2/km,'r--')
  plt.plot(r/kpc,pt/km,'k-')
  #Sofue12()
  plt.xlabel('R [kpc]')
  plt.ylabel('Vc [km/s]')
  plt.xlim(0,14)
  plt.ylim(0,300)
  plt.savefig('rc_curveD.png')

  print 'Plotting swing amplitude as a fn of R'
  plt.figure(7)
  swingM = r * kappa_freq(r)**2 / (4. * 3.14159 * G * SigDisk(r))
  meanM  = np.mean(swingM)
  solarM = np.mean( swingM[(r>6.9*kpc) & (r<7.1*kpc)]  )
  plt.plot(r/kpc,swingM,'ko',markersize=1.0,label='$ m(R) $')
  plt.plot(((0.,15.)),meanM*np.ones(2),'r--',label='$<m> = $'+str(np.round(meanM,3)))
  plt.plot(((0.,15.)),solarM*np.ones(2),'g--',label='$m(R=7\\rm kpc) = $'+str(np.round(solarM,3)))
  plt.ylim(0,10)
  plt.xlim(0,13)
  plt.xlabel('R[kpc]')
  plt.ylabel('m')
  plt.legend()
  print 'Mean Swing =',meanM
  print 'Solar Swing=',solarM
  
  plt.savefig('swingAmp_unifrom.png')
    

if ans1=='p':
  Npart=Npart_P    #Work with just a set for velocity setting (Npart_V<=Npart_P).
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #DRAW THE RADII:
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  r  = np.arange(r1,r2,0.001*(r2-r1))
  
  cd1=np.zeros(len(r))
  pd1=np.zeros(len(r))
  s1 =np.zeros(len(r))
  i  =0
  for line in cd1:
      cd1[i] = cdfDisk(r[i])
      pd1[i] = pdfDisk(r[i])
      s1[i]  = SigDisk(r[i])
      i+=1
  
  if ansP=='p':
    Am =plt.figure(1,figsize=(14, 6))
    A1=Am.add_subplot(1,2,1)
    A1.plot(r/kpc,cd1,label='cdf')
    A1.plot(r/kpc,pd1*kpc,label='pdf')
    A1.plot(r/kpc,s1, label='Sigma')
    plt.xlabel('R [kpc]')
    plt.xlim(r1/kpc,r2/kpc)
    plt.ylim(0,1.5)
    plt.legend(loc=2)
    
  rpred=randomdist(r,cd1,Npart)
  if ansP=='p':
    A1.hist(rpred/kpc,20,normed=True)
  
  phipred=np.ones(len(rpred))
  i=0
  for line in phipred:
    phipred[i]=2.*math.pi*random.random()
    i+=1
  
  dr   =0.05*(r2-r1)
  rlow =r1
  rhigh=r2
  rb1  =np.arange(rlow,rhigh,dr)
  rb2  =np.arange(rlow+dr,rhigh+dr,dr)
  s    =np.zeros(len(rb1))
  
  i=0
  for thing in s:
    s[i]= len(rpred[(rpred>rb1[i]) & (rpred<rb2[i])]) / (3.14159 * ((rb2[i]/kpc)**2-(rb1[i]/kpc)**2))
    i+=1
  scaler=60.0

  if ansP=='p':
    A1.errorbar(0.5*(rb1+rb2)/kpc,s/scaler,yerr=np.sqrt(s)/scaler,xerr=0.,fmt='k-o')
    A2=Am.add_subplot(1,2,2)
    A2.plot(rpred*np.cos(phipred)/kpc,rpred*np.sin(phipred)/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.savefig('DistsRad_D.png')
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #DRAW THE VERTICAL DISTANCE:
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  z  = np.arange(z1,z2,0.001*(z2-z1))

  cd2=np.zeros(len(z))
  pd2=np.zeros(len(z))
  s2 =np.zeros(len(z))

  i  =0
  for line in cd2:
      cd2[i] = cdfDiskz(z[i])
      pd2[i] = pdfDiskz(z[i])
      s2[i]  = SigDiskz(z[i])
      i+=1
  if ansP=='p':
    Bm =plt.figure(2,figsize=(14, 6))
    B1=Bm.add_subplot(1,2,1)
    B1.plot(z/kpc,cd2,label='cdf')
    B1.plot(z/kpc,0.5*pd2*kpc,label='pdf')   #Must be half to count for only using half dist.
    B1.plot(z/kpc,s2, label='Sigma')
    plt.xlabel('z [kpc]')
    plt.xlim(-z2/kpc,z2/kpc)
    plt.ylim(0,2.2)
    plt.legend(loc=2)
  
  "The z-dist is symmetric, so now randomly assign above/below z=0."  
  zpred=randomdist(z,cd2,Npart)
  i=0
  while i<len(zpred):
    if (random.random()>=0.5):
      zpred[i]=-zpred[i]
    i+=1
  if ansP=='p':
    B1.hist(zpred/kpc,20,normed=True)
    B1.plot(-z/kpc,0.5*pd2*kpc,c='green')

  dz   =0.05*(2*z2)
  zlow =-z2
  zhigh=+z2
  zb1  =np.arange(zlow,zhigh,dz)
  zb2  =np.arange(zlow+dz,zhigh+dz,dz)
  s    =np.zeros(len(zb1))
  
  i=0
  for thing in s:
    s[i]= len(zpred[(zpred>zb1[i]) & (zpred<zb2[i])]) / (3.14159 * ((zb2[i]/kpc)-(zb1[i]/kpc)))
    i+=1
  scaler=3800.0


  print 'Add some random offset to account for finite r/z array'
  print 'Scattering by',rscatter/kpc,'[kpc]'
  i=0
  while i<len(rpred):
    sig_posr = random.gauss(0.,rscatter)
    rpred[i] +=sig_posr
    #Ensure that this doesn't make stuff go r<0kpc:
    rpred[i] = abs(rpred[i])
    sig_posz = random.gauss(0.,rscatter)
    zpred[i] +=sig_posz
    i+=1
  

  if ansP=='p':
    B1.errorbar(0.5*(zb1+zb2)/kpc,s/scaler,yerr=np.sqrt(s)/scaler,xerr=0.,fmt='k-o')
    B2=Bm.add_subplot(1,2,2)
    B2.plot(rpred*np.cos(phipred)/kpc,zpred/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    plt.savefig('DistZ_D.png')

  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #Save Outputs
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  np.savez('positionsD_'+str(Npart),r=rpred,phi=phipred,z=zpred)

  exit()

elif ans1=='v':

  Npart=Npart_V    #Work with just a set for velocity setting (Npart_V<=Npart_P).
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #DRAW THE RADII:
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  inputz=np.load('positionsD_'+str(Npart_P)+'.npz')
  r   = inputz['r'][0:Npart]
  phi = inputz['phi'][0:Npart]
  z   = inputz['z'][0:Npart]
  x   = r*np.cos(phi)
  y   = r*np.sin(phi)
  if ansP=='p':
    plt.figure(1)
    plt.plot(x/kpc,y/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.savefig('xyPlotD.png')
  
  '''
  First we need the radom components of the radial, azimuthal, and verical velocities.
  These a taken from Gaussians centred on 0, 0 and v_phi respectively.
  '''
  sig_z  =np.zeros(len(r))
  sig_r  =np.zeros(len(r))
  sig_phi=np.zeros(len(r))
  v_phi  =np.zeros(len(r))
  
  "Begin with assigning sig_z"
  sig_z2 = math.pi*G*z0*SigDisk(r)
  i=0
  while i<len(sig_z):
    sig_z[i] = random.gauss(0.,np.sqrt(sig_z2[i]))
    i+=1
  if ansP=='p':
    plt.figure(2,figsize=(14, 6))
    plt.subplot(1,2,1)
    plt.scatter(x/kpc,y/kpc,c=sig_z/km)
    cbar=plt.colorbar()
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.subplot(1,2,2)
    plt.plot(r/kpc,sig_z/km,'ko',alpha=0.2)
    plt.xlabel('r [kpc]')
    plt.ylabel('$V_z \\rm [km/s]$',fontsize=16)
    plt.savefig('sigZ_D.png')
  
  "Then assign sig_r"
  SigVConst =ToomreConstFn(r)
  print 'The Sig_v constant is:',np.sqrt(SigVConst)/km,'km/s'
  sig_r2 = SigVConst * np.exp(-r/r0)
  i=0
  while i<len(sig_r):
    sig_r[i] = random.gauss(0.,np.sqrt(sig_r2[i]))
    i+=1  
  if ansP=='p':
    plt.figure(3,figsize=(14, 6))
    plt.subplot(1,2,1)
    plt.scatter(x/kpc,y/kpc,c=sig_r/km)
    cbar=plt.colorbar()
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.subplot(1,2,2)
    plt.plot(r/kpc,sig_r/km,'ko',alpha=0.2)
    plt.xlabel('r [kpc]')
    plt.ylabel('$V_r \\rm [km/s]$',fontsize=16)
    plt.savefig('sigR_D.png')
  
  "Then sig_phi, and finally the azimuthal streaming velocity, v_phi"
  sig_phi2 = sig_r2 * kappa_freq(r)**2 / (4.*omega_freq(r)**2)
  v_phi2   = sig_r2 - sig_phi2 - sig_r2*r/r0 + (Vc(r))**2
  i=0
  while i<len(sig_phi):
    sig_phi[i] = random.gauss(np.sqrt(v_phi2[i]),np.sqrt(sig_phi2[i]))
    i+=1  
  if ansP=='p':
    plt.figure(4,figsize=(14, 6))
    plt.subplot(1,2,1)
    plt.scatter(x/kpc,y/kpc,c=sig_phi/km)
    cbar=plt.colorbar()
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.subplot(1,2,2)
    plt.plot(r/kpc,sig_phi/km,'ko',alpha=0.2)
    plt.xlabel('r [kpc]')
    plt.ylabel('$\sigma_\phi \\rm [km/s]$',fontsize=16)
    plt.savefig('sigT_D.png')
  
  sig_phi = -sig_phi
  vx = sig_r*np.cos(phi) - sig_phi*np.sin(phi)
  vy = sig_r*np.sin(phi) + sig_phi*np.cos(phi)
  vz = sig_z
  
  MWgas = Mg
  MWstr = Md
  Nstr  = StarToGas*Npart
  Ngas  = Npart-Nstr

  if Ngas>0:
    mi_g  = MWgas/Ngas
  else:
    mi_g  = 0.

  if Nstr>0:
    mi_s  = MWstr/Nstr
  else:
    mi_s  = 0.
  "Phase flags"
  p_g   =0
  p_s   =10
  

  "Set masses, in kg:"
  mass  =np.zeros(len(r))
  phase =np.zeros(len(mass))          #0 for gas, 10 for stars
  i=0
  while i<len(mass):  #Auto fill as stars
     mass[i] = mi_s
     phase[i]= p_s
     if i>Nstr:
       mass[i] =mi_g
       phase[i]=p_g
     i+=1
  

  swingM = r * kappa_freq(r)**2 / (4. * 3.14159 * G * SigDisk(r))
  meanM  = np.mean(swingM)
  solarM = np.mean(swingM[(r>6.9*kpc) & (r<7.1*kpc)])
  if ansP=='p':
    plt.figure(6)
    vr = np.sqrt(vx*vx+vy*vy+vz*vz)
    plt.plot(r/kpc,vr/1000.,'ko',markersize=4.0)
    plt.xlabel('R[kpc]')
    plt.ylabel('Vr[km/s]')
    plt.savefig('rc_setD.png')
   
    print 'Plotting swing amplitude as a fn of R'
    plt.figure(7)
    plt.plot(r/kpc,swingM,'ko',markersize=1.0,label='$ m(R) $')
    plt.plot(((0.,15.)),meanM*np.ones(2),'r--',label='$m_{\\rm avg,i} = $'+str(np.round(meanM,3)))
    plt.xlabel('R[kpc]')
    plt.ylabel('m')
    plt.ylim(0,10)
    plt.xlim(0,13)
    plt.legend()
    plt.savefig('swingAmpPost.png')
  print 'Mean Swing =',meanM
  print 'Solar Swing=',solarM
    
  outd=open('asciifile_D','w')
  outg=open('asciifile_G','w')
  iout=0
  while iout<Nstr:
    outd.write(str(x[iout])+' '+str(y[iout])+' '+str(z[iout])+' '+str(mass[iout])+' '+str(vx[iout])+' '+str(vy[iout])+' '+str(vz[iout])+' '+str(phase[iout])+'\n')
    iout+=1
  outd.close()
  while iout<Nstr+Ngas:
    outg.write(str(x[iout])+' '+str(y[iout])+' '+str(z[iout])+' '+str(mass[iout])+' '+str(vx[iout])+' '+str(vy[iout])+' '+str(vz[iout])+' '+str(phase[iout])+'\n')
    iout+=1
  outg.close()
  exit()


elif ans1=='f':
  print 'Last step, so write out the file needed for phantom to setup.'
  Nstr  = StarToGas*Npart_P
  Ngas  = Npart_P-Nstr
  Nbulge = 0
  Nhalo  = 0
  outnom = 'galsetic.txt'
  outo=open(outnom,'w')
  iout=0
  outo.write('Ngas'+' '+str(Ngas)+'\n')
  outo.write('Nstar'+' '+str(Nstr)+'\n')
  outo.write('Nbulge'+' '+str(Nbulge)+'\n')
  outo.write('Ndark'+' '+str(Nhalo)+'\n')
  outo.close()
  print '<<<<<<<<<<<<<<DONE>>>>>>>>>>>>>>>'
  print 'Now place',outnom,'and the asciifiles in the directory of phantomsetup.'

  
else:
  print "No velocities or positions set. Erm, what's the point?"
  exit()
  
