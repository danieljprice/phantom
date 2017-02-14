"""
This code fits a normalised gaussian to the density PDF.

When run in a directory containing the density PDF files at each timestep,
it will search for files named 'cluster_?????_rho.pdf', and attempt to fit
a normalised gaussian to the data.
(Assuming first column is density, and second column is the PDF).

Once it has fit all the files it can find, it will return a file named
'sigma_vs_time.data', containing the standard deviation of the PDF (sigma)
as a function of time.

You also have the option to plot the data and its fit as it runs, so you can
see how well the fitting is doing.

NOTE: code is written as a function, so cannot be run on its own.
      It either needs to be called from another script, or run interactively
      through something such as iPython, or iPython Notebook.

Written by:
David Liptai, Monash University.
2015-2016
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def gaus(x,x0,sigma):
    return (1.0/np.sqrt(2.0*np.pi*sigma**2))*np.exp(-(x-x0)**2/(2.0*sigma**2))

def gaussian_fit_pdf():
    plt.close('all')

    SIGMA = []
    PLOT = False

    #def gaus(x,a,x0,sigma):
    #    return a*np.exp(-(x-x0)**2/(2*sigma**2))

    RMSE = []

    TIME = []

    if PLOT:
        plt.figure()

    nfiles_max = 1000
    for i in range(0,nfiles_max):
        #print(i)
        fname = 'cluster_{:0>5}_rho.pdf'.format(i)
        try:
            data = np.loadtxt(fname,skiprows=2)
            rho  = data[:,0]
            pdf  = data[:,1]
            time = float(np.genfromtxt(fname,skip_footer=len(rho)+1))
            TIME += [time]
        except:
            print('failed to find file: ',fname)
            break

        x = np.log(rho)

        n     = len(x)
        x0    = sum(x*pdf)/n
        sigma = sum((x-x0)**2)/n
        sigma = np.sqrt(sigma)
        #print(x0,sigma)

    #    params_opt, params_covar = curve_fit(gaus,x,pdf,p0=[1,x0,sigma])
        try:
            params_opt, params_covar = curve_fit(gaus,x,pdf,p0=[x0,sigma])
            print(time,abs(params_opt[1]))
            #SIGMA += [abs(params_opt[2])]
            SIGMA += [abs(params_opt[1])]
        #    print(params_opt[1])

            #y=gaus(x,params_opt[0],params_opt[1],params_opt[2])
            y=gaus(x,params_opt[0],params_opt[1])

            rmse = np.sqrt(np.mean((y-pdf)**2)/len(y))/(max(y)-min(y))
            RMSE += [rmse]
        #    print('a = ',params_opt[0],' x0 = ',params_opt[1],' sigma = ',params_opt[2],' rms error = ',rmse)
        #    print('rmse =',rmse)
            if PLOT:
                plt.cla()
                plt.title(i)
                plt.plot(x,pdf,label='data')
                plt.plot(x,y,label='fit')
                plt.plot(x,pdf-y,label='error')
                plt.legend()
                plt.show()
                plt.pause(1e-10)

        except:
            TIME=TIME[:-1]
            print('Could not fit: ',fname)
    '''
    plt.figure()
    plt.plot(TIME,RMSE)
    plt.xlabel('time')
    plt.ylabel('RMS error (normalised)')
    plt.show()

    plt.figure()
    plt.plot(TIME,abs(np.array(SIGMA)))
    plt.ylabel('sigma')
    plt.xlabel('time')
    plt.show()
    '''
    np.savetxt('sigma_vs_time.data',np.column_stack((TIME,SIGMA)))#,RMSE)))#,header="Time, Sigma, Normalised RMS error")
