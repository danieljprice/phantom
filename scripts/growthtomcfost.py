import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def pimp_my_sim(gdump_name,
                outdump_name,
                path_to_phantom,
                bins_per_dex=5,
                force_smax=False,
                smax_user=0.2,
                maxdustlarge=25,
                compil_logfile="compil.log",
                logfile_moddump="moddump.log",
                logfile_phantom="phantom.log",
                save_plot=False,
                scale="log",
                color="black",
                show=True):
    """
    pimp_my_sim transforms a dustgrowth dump into a multi large grains dump.
    it runs phantom on that dump for 1 time step to recompute the densities.
    the output dump is ready to be processed by mcfost through command line or pymcfost.
    
    input parameters are:
    
    gdump_name      : (str)   - name of dustgrowth dump (input)
    outdump_name    : (str)   - name of desired multi large grains dump (output)
    path_to_phantom : (str)   - path to phantom's directory
    bins_per_dex    : (int)   - number of bins per magnitude of size
    force_smax      : (bool)  - wether or not to force a maximum size for binning, else find it automatically
    smax_user       : (float) - size of forced maximum in cm
    maxdustlarge    : (int)   - maximum number of large dust species for memory allocation in phantom
    compil_logfile  : (str)   - name of logfile to store compilation output
    logfile_moddump : (str)   - name of logfile to store moddump run output
    logfile_phantom : (str)   - name of logfile to store phantom run output
    save_plot       : (bool)  - wether or not to save plot of size distribution in pdf file, else shows it interactively
    scale           : (str)   - npart axis scale, options are "linear" or "log"
    color           : (str)   - histogram color
    show            : (bool)  - wether or not to show the distribution, only applies if save_plot=True
    """
    
    ############### First part : moddump ###############

    # compiling phantommoddump using designated setup
    print("--> Compiling Phantommoddump using growthtomulti setup...")
    if os.path.exists("phantommoddump"):
        os.system("rm phantommoddump")
    make_exec(path_to_phantom, setup="growthtomulti", exec="moddump", compil_logfile=compil_logfile, options=f"MAXDUSTLARGE={maxdustlarge}")
    
    # check for errors during the compilation
    if os.path.exists("phantommoddump"):
        print(f"--> Compilation Successful.")
    else:
        print(f"ERROR: check {compil_logfile} for more informations")

    # create input param file for moddump
    with open("bin_param.txt", "w") as paramfile:
        line = f"{force_smax} {smax_user} {bins_per_dex}"
        paramfile.write(line)

    # run moddump on input dustgrowth dump
    print(f"--> Running Phantommoddump on input file {gdump_name}...")
    if os.path.exists(outdump_name+".in"):
        os.system(f"rm {outdump_name}*")
    os.system(f"./phantommoddump {gdump_name} {outdump_name}_00000.tmp &> {logfile_moddump}")

    # check for errors in the execution of phantommoddump
    if os.path.exists(outdump_name+"_00000.tmp"):
        print("--> Operation Successful.")
    else:
        print(f"ERROR: check {logfile_moddump} for more informations")

    # plot and open the size distribution to check for potential caveats or empty bins, force the user to be aware of it
    plot_distrib("bin_distrib.dat", dump=gdump_name, save_plot=save_plot, scale=scale, color=color, show=show)

    ############### Second part : run phantom for 1 time step on new dump ###############
    
    # compiling phantom for multi large grains using dustydisc setup
    print("--> Compiling Phantom using dustydisc setup...")
    if os.path.exists("phantom"):
        os.system("rm phantom")
    make_exec(path_to_phantom, setup="dustydisc", compil_logfile=compil_logfile, options=f"MAXDUSTLARGE={maxdustlarge}")
    
    # check for errors during the compilation
    if os.path.exists("phantom"):
        print(f"--> Compilation Successful.")
    else:
        print(f"ERROR: check {compil_logfile} for more informations")
    
    # run phantom on tmp dump
    print(f"--> Running Phantom with input file {outdump_name}.in...")
    os.system(f"./phantom {outdump_name}.in &> {logfile_phantom}")
    
    # check for errors during the run
    if os.path.exists(outdump_name+"_00000"):
        print("--> Operation Successful.\n    Thank you for using pimp_my_sim!")
        print("   --> You can now run mcfost using: mcfost <paramfile> -phantom -<options>")
    else:
        print(f"ERROR: check {logfile_phantom} for more informations")
    
def make_exec(path_to_phantom, setup, compil_logfile, exec="", options=""):
    os.system(f"{path_to_phantom}/scripts/writemake.sh {setup} > Makefile")
    os.system(f"make {exec} {options} &> {compil_logfile}")
    
def plot_distrib(file, dump, show=True, save_plot=False, scale="log", color=None):
    # read data
    names  = ("bin", "smin", "s", "smax", "npart")
    data   = pd.read_csv(file, delim_whitespace=True, names=names)
    
    # compute width of each bin for histogram
    ntypes = len(data.s)-1
    width_s = np.diff(data.s)
    width_s = np.append(width_s, data.smax[ntypes]-data.smin[ntypes])
    
    # plot size distribution
    plt.bar(data.s/1.e4, data.npart, width=width_s/1.e4, align="edge", edgecolor="white", label=f"{len(data.s)} bins", color=color)
    plt.xscale("log")
    plt.yscale(scale)
    plt.xlabel("s [cm]")
    plt.ylabel(r"$n_\mathrm{part}$")
    plt.legend()
    plt.title(dump)
    
    # save and show plot (or not!)
    print("--> Please check dust size binning")
    if save_plot:
        plt.savefig("distrib.pdf", bbox_inches="tight")
        if show:
            os.system("open distrib.pdf")
    elif show:
        plt.show()
    
    