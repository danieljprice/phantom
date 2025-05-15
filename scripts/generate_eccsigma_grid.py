import numpy as np

######################################################
# Generate sigma and eccentricity profile files      #
# NB: xmin and xmax scake is not relevant as         #
# rescales radii, but make sure xmin/xmax is correct #
# if you want to conserve geometry of the profile    #
######################################################

def sigma_prof(x):
    p=1
    Rc=6
    return x**(-p)#*np.exp(-x/Rc)

def ecc_prof(x):
    e0=0.1
    q=1
    return e0*x**(-q) 

x=np.linspace(1,10,100)
ecc=ecc_prof(x)
sigma=sigma_prof(x)

data_ecc=np.zeros((x.size,2))
data_sig=np.zeros((x.size,2))

data_ecc[:,0]=x
data_ecc[:,1]=ecc
data_sig[:,0]=x
data_sig[:,1]=sigma

np.savetxt('sigma_grid.dat',data_sig,delimiter=' ',header='a,sigma')
np.savetxt('ecc_grid.dat',data_ecc,delimiter=' ',header='a,ecc')
