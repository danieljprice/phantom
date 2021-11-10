import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy import interpolate, optimize, constants
import time
import readPhantomDump as RD

start = time.time()
  
# CONSTANTS [SI] ========================================

au    = 1.496e11    
Gr    = 6.67408e-11 
Msol  = 1.9891e30   
kb    = 1.38064852e-23
mH    = 1.6735575e-27 
mu    = 2.38

# VARIABLES ========================================

# Physical parameters [SI]
M_star       = 1.2*Msol
R_star      = 1.25*au
T_star      = 5.0e4
T_exp       = 0.5
gamma_adia  = 1.4

# Numerical parameters [SI]
r_inject       = 2.0*au
N_bound_sphere = 5

# Plotting options [SI]
r_outer_plot = 10*R_star
r_inner_plot = 0.5*au

# FUNCTIONS ========================================
  
def cs_array(gamma,T):   
    c     = np.sqrt(gamma*kb*T/(mu*mH))
    return c
  
def Temp(r,Tmax,eps):
    return Tmax*(R_star/r)**eps

# INITIALISATION ========================================

which_plot = 'velocity'  # velocity, density, energy

#maindir = "/lhome/ward/Desktop/SPH_numerics/Phantom_res_study/data/"
#maindir = "/lhome/ward/Desktop/SPH_numerics/Phantom_res_study/data_quinticKernel_iter1/"
maindir = "/lhome/ward/Desktop/SPH_numerics/Phantom_res_study/PaperData/quinticKernel/standard/"

#RES = np.linspace(1,10,10)
#WSS = np.linspace(0.2,2,10)
RES = [10]
WSS = [0.8]

dispersion    = np.zeros((len(RES),len(WSS)))
h_to_rc_ratio = np.zeros((len(RES),len(WSS)))

setups = ['noDust_adia_transsonic','noDust_adia_supersonic']
#setups = os.listdir(maindir)


# LOOP THROUGH MODELS ========================================
  
for res in RES:
  
  print('===============')
  print('RES = '+str(res))
  print('===============')
  
  for wss in WSS:
    
    print('----------')
    print('WSS = '+str(wss))
    print('----------')
    
    ## Make Figure canvas
    fig, axs = plt.subplots(1, 1)

    for setup in setups:
      
      print(setup)
      
      setupdir       = maindir + setup + '/' 
      thermodynamics = setup.split('_')[1].strip()
      dynamics       = setup.split('_')[2].strip()
    
      model      = 'RES_'+str(int(res))+'_WSS_'+str(round(wss,1))
      modeldir   = setupdir + model +'/'
      dump_in    = modeldir + 'wind_00030.ascii'
      profile_in = modeldir + 'windprofile1D.dat'
        
      #data_dump = RD.read_dump(dump_in)
      
      ## Extract 3D data
      #udist, umass, utime = [data_dump['units'][u] for u in ('udist', 'umass', 'utime')]
      #hfact = data_dump['quantities']['hfact']
      #mpart = data_dump['quantities']['massoftype'][0]*umass
      #x, y, z, vx, vy, vz, h, u = [np.array(data_dump['blocks'][0]['data'][c],dtype=float)*conv for c, conv in [('x', udist), ('y', udist), ('z', udist), ('vx', udist/utime), ('vy', udist/utime), ('vz', udist/utime), ('h', udist), ('u', udist**2/utime**2)]]
      #r = np.sqrt(x**2+y**2+z**2)/au
      #v = np.sqrt(vx**2+vy**2+vz**2)
      #rho = mpart * (hfact/h)**3
      #del x, y, z, vx, vy, vz, h, hfact, mpart

      
      data_dump    = np.genfromtxt(dump_in)
      data_dump    = [*zip(*data_dump[:-1])]
      data_profile = np.genfromtxt(profile_in)
      data_profile = [*zip(*data_profile[:-1])]
      
      # Construct arrays from SPH data (to SI units)
      r_sph   = np.sqrt(np.array(data_dump[0])**2+np.array(data_dump[1])**2+np.array(data_dump[2])**2)*au
      v_sph   = np.sqrt(np.array(data_dump[6])**2+np.array(data_dump[7])**2+np.array(data_dump[8])**2)*1e-2
      rho_sph = np.array(data_dump[5])*1e+3
      u_sph   = np.array(data_dump[9])*1e-4
      
      # Only select SPH particles between R_inner and R_outer and sort according to R
      index   = np.where((r_sph>r_inner_plot) & (r_sph<r_outer_plot))
      r_sph   = r_sph[index[0]]
      v_sph   = v_sph[index[0]]
      rho_sph = rho_sph[index[0]]
      u_sph   = u_sph[index[0]]
      
      # Construct arrays from 1d profile
      r_solution   = np.array(data_profile[1])*1e-2
      v_solution   = np.array(data_profile[2])*1e-2
      rho_solution = np.array(data_profile[6])*1e+3
      P_solution   = np.array(data_profile[5])*1e-1
      if thermodynamics == 'adia':
        u_solution   = P_solution/((gamma_adia-1)*rho_solution)
      
      # Calculate position of last boundary shell
      if dynamics == 'transsonic':
        d_tang       = r_inject*2*(constants.golden*np.sqrt(5))**(-1/2)/(2*res-1)
        d_rad        = wss*d_tang
        r_last_shell = r_inject+N_bound_sphere*d_rad
      
      # Create sorted r array for interpolation of solution
      r_sol_sort     = r_sph[np.argsort(r_sph)]
        
      if which_plot == 'energy':
        
        # Interpolate 1d profile and resample to radii of SPH data        
        u_interp       = interpolate.interp1d(r_solution, u_solution, kind='cubic', bounds_error=False, fill_value='extrapolate')
        u_solution_new = u_interp(r_sph)
        u_sol_sort     = u_solution_new[np.argsort(r_sph)]
        
        # Plot quantities
        R_ref = R_star
        u_ref = 1.
      
        if dynamics == 'transsonic': axs.plot(r_sph[0]/R_ref,u_sph[0]/u_ref,'r.',label='simulation')
        axs.plot(r_sph/R_ref,u_sph/u_ref,'r.',markersize=3)
        if dynamics == 'transsonic': axs.plot(r_sph[0]/R_ref,u_solution_new[0]/u_ref,'k-',label='exact')
        axs.plot(r_sol_sort/R_ref,u_sol_sort/u_ref,'k-',markersize=3)
        if dynamics == 'transsonic': axs.axvline(x=r_inject/R_ref,ymin=0,ymax=3,color='magenta',linestyle='-',label=r'R$_{\rm inject}$')
        if dynamics == 'transsonic': axs.axvline(x=r_last_shell/R_ref,ymin=0,ymax=3,color='orange',linestyle='-',label='last boundary shell')
        axs.set_ylabel('u [J/kg]')
        axs.legend(loc='upper right')
        axs.set_xlim([1,r_outer_plot/R_star])
        axs.set_xlabel('distance [Rstar]')
        axs.set_yscale('log')
        axs.set_xscale('log')
        
      if which_plot == 'density':
        
        # Interpolate 1d profile and resample to radii of SPH data        
        rho_interp       = interpolate.interp1d(r_solution, rho_solution, kind='cubic', bounds_error=False, fill_value='extrapolate')
        rho_solution_new = rho_interp(r_sph)
        rho_sol_sort     = rho_solution_new[np.argsort(r_sph)]
        
        # Plot quantities
        R_ref = R_star
        rho_ref = 1.
      
        if dynamics == 'transsonic': axs.plot(r_sph[0]/R_ref,rho_sph[0]/rho_ref,'r.',label='simulation')
        axs.plot(r_sph/R_ref,rho_sph/rho_ref,'r.',markersize=3)
        if dynamics == 'transsonic': axs.plot(r_sph[0]/R_ref,rho_solution_new[0]/rho_ref,'k-',label='exact')
        axs.plot(r_sol_sort/R_ref,rho_sol_sort/rho_ref,'k-',markersize=3)
        if dynamics == 'transsonic': axs.axvline(x=r_inject/R_ref,ymin=0,ymax=3,color='magenta',linestyle='-',label=r'R$_{\rm inject}$')
        if dynamics == 'transsonic': axs.axvline(x=r_last_shell/R_ref,ymin=0,ymax=3,color='orange',linestyle='-',label='last boundary shell')
        axs.set_ylabel('rho [kg/m^3]')
        axs.legend(loc='upper right')
        axs.set_xlim([1,r_outer_plot/R_star])
        axs.set_xlabel('distance [Rstar]')
        axs.set_yscale('log')
        axs.set_xscale('log')
      
      elif which_plot == 'velocity':
        
        # Interpolate 1d profile and resample to radii of SPH data        
        v_interp       = interpolate.interp1d(r_solution, v_solution, kind='cubic', bounds_error=False, fill_value='extrapolate')
        v_solution_new = v_interp(r_sph)
        v_sol_sort     = v_solution_new[np.argsort(r_sph)]
        
        # Calculate the sound speed
        if thermodynamics == 'iso':
          gamma = 1.0
          T     = T_star*np.ones(r_solution.shape)
          cs    = cs_array(gamma,T)
        elif thermodynamics == 'Tr':
          gamma = 1.0
          T     = T_star*(R_star/r_solution)**T_exp
          cs    = cs_array(gamma,T)
        elif thermodynamics == 'adia':
          gamma = gamma_adia
          u_sph = np.array(data_dump[9])*1e-4
          T_sph = mu*mH*(gamma-1)*u_sph[index[0]]/kb

          sph_index_sort = np.argsort(r_sph)
          r_sph_sort     = r_sph[sph_index_sort]
          T_sph_sort     = T_sph[sph_index_sort]

          uniques    = np.unique(r_sph_sort,return_index=True)
          r_sph_sort = r_sph_sort[uniques[1]]
          T_sph_sort = T_sph_sort[uniques[1]]
          
          popt, pcov = optimize.curve_fit(Temp, r_sph_sort, T_sph_sort)
          T          = Temp(r_solution,*popt)          
          cs         = cs_array(gamma,T)
                  
        # Interpolate sounds speed profile        
        if dynamics == 'transsonic':
          cs_interp = interpolate.interp1d(r_solution, cs, kind='cubic', bounds_error=False, fill_value='extrapolate')
          
        # Calculate escape velocity array
        r_esc = np.linspace(r_inner_plot,r_outer_plot,200)
        v_esc = np.sqrt(2*Gr*M_star/r_esc)
        
        # Plot quantities
        R_ref = R_star
        V_ref = 1000.
        
        if dynamics == 'transsonic': axs.plot(r_sph[0]/R_ref,v_sph[0]/V_ref,'r.',label='simulation')
        axs.plot(r_sph/R_ref,v_sph/V_ref,'r.',markersize=3)
        if dynamics == 'transsonic': axs.plot(r_sph[0]/R_ref,v_solution_new[0]/V_ref,'k-',label='exact')
        axs.plot(r_sol_sort/R_ref,v_sol_sort/V_ref,'k-',markersize=3)
        if dynamics == 'transsonic': axs.axvline(x=r_inject/R_ref,ymin=0,ymax=3,color='magenta',linestyle='-',label=r'R$_{\rm inject}$')
        if dynamics == 'transsonic': axs.axvline(x=r_last_shell/R_ref,ymin=0,ymax=3,color='orange',linestyle='-',label='last boundary shell')
        if dynamics == 'transsonic': axs.plot(r_sol_sort/R_ref,cs_interp(r_sol_sort)/V_ref,color='green',linestyle='--',label='sound speed')
        if dynamics == 'transsonic': axs.plot(r_esc/R_ref,v_esc/V_ref,color='blue',linestyle='-',label=r'v$_{\rm esc}$')
        axs.set_ylabel('v [km/s]')
        if dynamics == 'supersonic': axs.set_ylim([0,1.1*max(v_sph/V_ref)])
        axs.legend(loc='lower right')
        axs.set_xlim([1,r_outer_plot/R_star])
        axs.set_xlabel('distance [Rstar]')
        
    plt.show()
    #if not os.path.exists(setupdir+'0_PLOTS/'):
      #os.makedirs(setupdir+'0_PLOTS/')
    #plt.savefig('/lhome/ward/Desktop/SPH_numerics/Paper_plots/'+setups[0]+'/'+model+'.png',dpi=200)
    #plt.savefig(setupdir+'0_PLOTS/'+model+'.png',dpi=200)
    plt.close()

end = time.time()
hours, rem = divmod(end-start, 3600)
minutes, seconds = divmod(rem, 60)
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))   
