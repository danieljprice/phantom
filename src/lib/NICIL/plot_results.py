#----------------------------------------------------------------------!
#                               N I C I L                              !
#           Non-Ideal mhd Coefficient and Ionisation Library           !
#         A Control script to allow the user to generate graphs        !
#                         of the relevant data.                        !
#                                                                      !
#                 Copyright (c) 2015-2017 James Wurster                !
#        See LICENCE file for usage and distribution conditions        !
#----------------------------------------------------------------------!
import os
import sys
#
# Defaults (these are useful values that may be temporarily altered for more useful graphing)
open_eps = "open"  # Program used to open the .eps files; leave blank to prevent auto-opening (e.g. open_eps="open")
setxr    = True    # Use the updated default x-range; else, gnuplot will automatically select the range
setyr    = True    # Use the updated default y-range; else, gnuplot will automatically select the range
incl_H12 = False   # Include H2 on plots
#
# Verify that Gnuplot exists.  If it does not, ask the user if they wish to contine
does_gnuplot_exist = os.system("gnuplot --version")
if (does_gnuplot_exist<>0):
  print " "
  print "Gnuplot does not exist on this computer"
  print "Would you like to continue and make the Gnuplot script?"
  print "This scipt can be used as a template to make the graphs in an exising graphing programme."
  makescript = raw_input("Press Y then Enter to make the script; otherwise press Enter to exit: ")
  if (makescript<>"Y" and makescript<>"y"): sys.exit() 
#
# Determine if an old or new version of Gnuplot is being used, where 5.0 and newer is defined as new
gnuplot5    = True
bashcommand = "gnuplot --version"
line        = os.popen(bashcommand).read().strip()
if (" 0." in line or " 1." in line or " 2." in line or " 3." in line or " 4." in line): gnuplot5 = False
#
# Print instructions if no default viewer is set
if (open_eps==""): 
  print " "
  print "!---------------------------------------------------------------------------------------------------!"
  print "! Graphs will not be automatically opened.                                                          !"
  print "! Please set 'open_eps' to the programme to which you wish to use to open the graphs automatically. !"
  print "!---------------------------------------------------------------------------------------------------!"
  print " "
#
# Determine if eta or sph comparison
escomp = ""
ierr     = 0
print "What would you like to graph?  The options are"
print "1: Graph the results from nicil_ex_eta"
print "2: Graph the results from nicil_ex_sph"
while (escomp<>"1" and escomp<>"2"):
  if (ierr==1): print "That is an invalid entry. Please try again"
  escomp = raw_input("Enter the number for the desired graph: ")
  if (escomp=="q" or escomp=="Q"): sys.exit()
  ierr     = 1
#
# Determine if plotting against number density or temperature
if (escomp=="1"):
  xaxis = ""
  ierr     = 0
  print " "
  print "What would you like on the horizontal axis?  The options are"
  print "1: number density with T = constant"
  print "2: number density with T = T(n) {barotropic EOS}"
  print "3: temperature with n = constant"
  print "4: temperature with n = n(T) {barotropic EOS}"
  print "5: all of the above"
  while (xaxis<>"1" and xaxis<>"2" and xaxis<>"3" and xaxis<>"4" and xaxis<>"5"):
    if (ierr==1): print "That is an invalid entry. Please try again"
    xaxis = raw_input("Enter the number for the desired axis: ")
    if (xaxis=="q" or xaxis=="Q"): sys.exit()
    ierr     = 1
else:
  xaxis = "1"
#
# Set file (character) IDs & define the plots
if (xaxis=="5"):
  imin = 1
  imax = 5
else:
  imin = int(xaxis)
  imax = imin + 1
graphscript = "Graphs/plot.gnuplot"
graphfiles = ""
a=open(graphscript,"w")
for i in range(imin,imax):
  gtype = i
  if (escomp=="1"):
    if (gtype==1):
      datafile = "data/eta_density.dat"
      graphfile = "Graphs/eta_density_Tcnst.eps"
    elif (gtype==2):
      datafile = "data/eta_barotropic.dat"
      graphfile = "Graphs/eta_density_barotropic.eps"
    elif (gtype==3):
      datafile = "data/eta_temperature.dat"
      graphfile = "Graphs/eta_temperature_ncnst.eps"
    elif (gtype==4):
      datafile = "data/eta_barotropic.dat"
      graphfile = "Graphs/eta_temperature_barotropic.eps"
  else:
    datafile = "data/sph_density.dat"
    graphfile = "Graphs/sph_density.eps"
  graphfiles = graphfiles+" "+graphfile
  if (not os.path.isfile(datafile) ): 
    print "Data file, ",datafile,", does not exist.  Exiting graphing script."
    sys.exit()
  #
  #
  # Make the graphing file
  a.write("reset \n")
  a.write("# Manually modifying the style types \n")
  if (gnuplot5):
    a.write("set style line  1 lc 1 lw 3 ps 2 \n")
    a.write("set style line  2 lc 2 lw 3 ps 2 \n")
    a.write("set style line  3 lc 3 lw 3 ps 2 \n")
    a.write("set style line  4 lc 4 lw 3 ps 2 \n")
    a.write("set style line  5 lc 6 lw 3 ps 2 \n")
    a.write("set style line  6 lc 5 lw 3 ps 2 \n")
    a.write("set style line  7 lc 7 lw 3 ps 2 \n")
    a.write("set style line  8 lc 8 lw 3 ps 2 \n")
    a.write("set style line  9 lc 9 lw 3 ps 2 \n")
    a.write("set style line 11 lc 1 lw 1 ps 2 dt 2 \n")
    a.write("set style line 12 lc 2 lw 1 ps 2 dt 2 \n")
    a.write("set style line 13 lc 3 lw 1 ps 2 dt 2 \n")
    a.write("set style line 14 lc 4 lw 1 ps 2 dt 2 \n")
    a.write("set style line 15 lc 6 lw 1 ps 2 dt 2 \n")
    a.write("set style line 17 lc 7 lw 1 ps 2 dt 2 \n")
    a.write("set style line 18 lc 8 lw 1 ps 2 dt 2 \n")
    a.write("set style line 19 lc 9 lw 1 ps 2 dt 2 \n")
    a.write("set style line 21 lc 1 lw 6 ps 2 \n")
    a.write("set style line 27 lc 7 lw 6 ps 2 \n")
    a.write("set style line 29 lc 9 lw 6 ps 2 \n")
  else:
    a.write("set style line  1 lt 1 lc 1 lw 3 ps 2 \n")
    a.write("set style line  2 lt 1 lc 2 lw 3 ps 2 \n")
    a.write("set style line  3 lt 1 lc 3 lw 3 ps 2 \n")
    a.write("set style line  4 lt 1 lc 4 lw 3 ps 2 \n")
    a.write("set style line  5 lt 1 lc 6 lw 3 ps 2 \n")
    a.write("set style line  6 lt 1 lc 5 lw 3 ps 2 \n")
    a.write("set style line  7 lt 1 lc 7 lw 3 ps 2 \n")
    a.write("set style line  8 lt 1 lc 8 lw 3 ps 2 \n")
    a.write("set style line  9 lt 1 lc 9 lw 3 ps 2 \n")
    a.write("set style line 11 lt 3 lc 1 lw 1 ps 2 \n")
    a.write("set style line 12 lt 3 lc 2 lw 1 ps 2 \n")
    a.write("set style line 13 lt 3 lc 3 lw 1 ps 2 \n")
    a.write("set style line 14 lt 3 lc 4 lw 1 ps 2 \n")
    a.write("set style line 15 lt 3 lc 6 lw 1 ps 2 \n")
    a.write("set style line 17 lt 3 lc 7 lw 1 ps 2 \n")
    a.write("set style line 18 lt 3 lc 8 lw 1 ps 2 \n")
    a.write("set style line 19 lt 3 lc 9 lw 1 ps 2 \n")
    a.write("set style line 21 lt 1 lc 1 lw 6 ps 2 \n")
    a.write("set style line 27 lt 1 lc 7 lw 6 ps 2 \n")
    a.write("set style line 29 lt 1 lc 9 lw 6 ps 2 \n")
  a.write("  \n")
  a.write("m = 2.310*1.6726219e-24 \n")
  if (escomp=="1"):
    a.write("set terminal postscript eps enhanced colour size 10,10 font 'Times-Roman,20' \n")
  else:
    a.write("set terminal postscript eps enhanced colour size 14,10 font 'Times-Roman,20' \n")
  a.write("set output '"+graphfile+"' \n")
  a.write("  \n")
  if (escomp=="1"):
    # Spacing of points for "every"
    a.write("set multiplot layout 2,2 \n")
    a.write("set log x \n")
    if (gtype==1 or gtype==2):
      a.write("set xl  'log n_n(cm^{-3})' \n")
      if (setxr):
        a.write("set xr  [1e2  :1e22  ] \n")
        a.write("set x2r [1e2*m:1e22*m] \n")
        a.write("set x2l 'log {/Symbol r}_n (g cm^{-3})' \n")
        a.write("set log x2 \n")
        a.write("set format x2 '%L' \n")
        a.write("set x2tics \n")
        a.write("set xtics nomirror \n")
      xmin = "1.0e2"
      x = "13"
    else:
      if (setxr): a.write("set xr [1e1:1e5] \n")
      a.write("set xl 'log T (K)' \n")
      xmin = "1.0e1"
      x = "2"
    a.write("set format x '%L' \n")
    a.write("set key top left \n")
    a.write(" \n")
    a.write("#Plot number densities \n")
    a.write("set log y \n")
    a.write("set yl 'log n (cm^{-3})' \n")
    if (setyr): a.write("set yr [1e-11:1e24] \n")
    a.write("set format y ' %3L' \n")
    a.write("plot '"+datafile+"' u "+x+":12 ti 'n_e'            w l ls 29 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":14 ti 'n_{i,Rh}'       w l ls 13 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":15 ti 'n_{i,Rm}'       w l ls  3 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":16 ti 'n_{i,T1}'       w l ls  4 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":17 ti 'n_{i,T2}'       w l ls 14 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":18 ti 'n_{g,R}(Z=-1)'  w l ls 27 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":19 ti 'n_{g,R}(Z= 0)'  w l ls  6 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":20 ti 'n_{g,R}(Z=+1)'  w l ls  5 ,\\\n")
    if (incl_H12):
      a.write("   '"+datafile+"' u "+x+":21 ti 'n_{H_2}'        w l ls  8 ,\\\n")
      a.write("   '"+datafile+"' u "+x+":22 ti 'n_H'            w l ls 18 ,\\\n")
    a.write("     '"+datafile+"' u ("+xmin+"):(1) ti ''         w d          \n") # To avoid crashing if no data
    a.write(" \n")
    a.write("#Plot ionisation densities \n")
    a.write("set log y \n")
    a.write("set yl 'log n (cm^{-3})' \n")
    if (setyr): a.write("set yr [1e-11:1e24] \n")
    a.write("set format y ' %3L' \n")
    a.write("plot '"+datafile+"' u "+x+":12 ti 'n_e'    w l ls 29 ,\\\n")
    if (incl_H12):
      a.write("   '"+datafile+"' u "+x+":23 ti 'H_2+'   w l ls  8 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":24 ti 'H+'     w l ls  5 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":25 ti 'He+'    w l ls  3 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":26 ti 'Na+'    w l ls  4 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":27 ti 'Mg+'    w l ls  1 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":28 ti 'K+'     w l ls  7 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":29 ti 'He++'   w l ls 13 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":30 ti 'Na++'   w l ls 14 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":31 ti 'Mg++'   w l ls 11 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":32 ti 'K++'    w l ls 17 ,\\\n")
    a.write("     '"+datafile+"' u ("+xmin+"):(1) ti '' w d          \n") # To avoid crashing if no data
    a.write(" \n")
    a.write("#Plot conductivities \n")
    a.write("set log y \n")
    a.write("set yl 'log {/Symbol s} (s^{-1})' \n")
    if (setyr): a.write("set yr [5e-8:1e22] \n")
    a.write("set key top left \n")
    a.write("set format y ' %3L' \n")
    a.write("plot '"+datafile+"' u "+x+":   7  ti '{/Symbol s}_O'     w l ls 21 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":($8>0.0? $8:0/0) ti '{/Symbol s}_H > 0' w l ls  3 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":($8<0.0?-$8:0/0) ti '{/Symbol s}_H < 0' w l ls  4 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":   9  ti '{/Symbol s}_P'     w l ls  5 ,\\\n")
    a.write("     '"+datafile+"' u ("+xmin+"):(1) ti ''               w d         \n") # To avoid crashing if no data
    a.write(" \n")
    a.write("#Plot resistivities \n")
    a.write("set key top left \n")
    a.write("set yl 'log {/Symbol h} (cm^2 s^{-1})' \n")
    if (setyr): a.write("set yr [1e0:1e24] \n")
    a.write("plot '"+datafile+"' u "+x+":   4  ti '{/Symbol h}_{OR}'     w l ls 21 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":($5>0.0? $5:0/0) ti '{/Symbol h}_{HE} > 0' w l ls  3 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":($5<0.0?-$5:0/0) ti '{/Symbol h}_{HE} < 0' w l ls  4 ,\\\n")
    a.write("     '"+datafile+"' u "+x+":   6  ti '{/Symbol h}_{AD}'     w l ls  5 ,\\\n")
    a.write("     '"+datafile+"' u ("+xmin+"):(1) ti ''                  w d         \n") # To avoid crashing if no data
    if ( (gtype==2 or gtype==4) and False ):
      a.write(" \n")
      a.write("#Plot Temperature-density relationship \n")
      if (gtype==2):
        if (setyr): a.write("set yr [1e0:1e5] \n")
        a.write("set yl 'log T (K)' \n")
        y = "2"
      else:
        if (setyr): a.write("set yr [1e5:1e25] \n")
        a.write("set yl 'log n_n(cm^{-3})' \n")
        y = "14"
      a.write("plot '"+datafile+"' u   "+x+"   :"+y+" ti ''     w l ls 9 ,\\\n")
      a.write("     '"+datafile+"' u ("+xmin+"):(1)   ti ''     w d         \n") # To avoid crashing if no data
  else:
    intval = "25"  # Spacing of points for "every"
    xmin   = "1.8"
    a.write("set multiplot layout 2,3 \n")
    if (setxr): a.write("set xr [1.8:15.5] \n")
    a.write("set xl '{/Symbol r} (10^{-19} g cm^{-3})' \n")
    a.write(" \n")
    a.write("#Plot magnetic density \n")
    a.write("set yl 'J_y (10^{-52} G cm^{-3})' \n")
    if (setyr): a.write("set yr [-9:0] \n")
    a.write("plot '"+datafile+"' u ($1/1.0e-19):($3/1.0e-52) ti '' w l ls 9 ,\\\n")
    a.write("     '"+datafile+"' u ("+xmin+"):(1)            ti '' w d         \n") # To avoid crashing if no data
    a.write(" \n")
    a.write("#Plot du/dt_{non-ideal} \n")
    a.write("set yl 'dudt_{non-ideal} (erg g^{-1} s^{-1})' \n")
    a.write("set log y \n")
    a.write("set format y '10^{%3L}' \n")
    if (setyr): a.write("set yr [1e-7:1e-4] \n")
    a.write("plot '"+datafile+"' u ($1/1.0e-19):8 ti '' w l ls 9 ,\\\n")
    a.write("     '"+datafile+"' u ("+xmin+"):(1) ti '' w d         \n") # To avoid crashing if no data
    a.write(" \n")
    a.write("#Plot dB/dt_{non-ideal} \n")
    a.write("set yl 'dBdt_{non-ideal,z} (G s^{-1})' \n")
    if (setyr): a.write("set yr [1e-15:1e-12] \n")
    a.write("plot '"+datafile+"' u ($1/1.0e-19):4 ti '' w l ls 9 ,\\\n")
    a.write("     '"+datafile+"' u ("+xmin+"):(1) ti '' w d         \n") # To avoid crashing if no data
    a.write(" \n")
    a.write("#Plot resistivities \n")
    a.write("set log y \n")
    a.write("set yl '{/Symbol h} (cm^2 s^{-1})' \n")
    a.write("set key top left \n")
    if (setyr): a.write("set yr [1e7:1e20] \n")
    a.write("plot '"+datafile+"' u ($1/1.0e-19):5     ti '{/Symbol h}_{OR}'     w l ls  1 ,\\\n")
    a.write("     '"+datafile+"' u ($1/1.0e-19):6     ti '{/Symbol h}_{HE} > 0' w l ls  3 ,\\\n")
    a.write("     '"+datafile+"' u ($1/1.0e-19):(-$6) ti '{/Symbol h}_{HE} < 0' w l ls  4 ,\\\n")
    a.write("     '"+datafile+"' u ($1/1.0e-19):7     ti '{/Symbol h}_{AD}'     w l ls  5 ,\\\n")
    a.write("     '"+datafile+"' u ("+xmin+"):(1)     ti ''                     w d          \n") # To avoid crashing if no data
    a.write(" \n")
    a.write("#Plot Hall & ion drift velocities (2 plots) \n")
    a.write("unset log y \n")
    a.write("set format y '%g' \n")
    a.write("set yl 'v_{drift} (cm s^{-1})' \n")
    a.write("set key top left \n")
    a.write("set yr [0:3] \n")
    a.write("plot '"+datafile+"' u ($1/1.0e-19):10  ti 'v_{hall,y} ' w l ls 3 ,\\\n")
    a.write("     '"+datafile+"' u ("+xmin+")  :(1) ti ''            w d         \n") # To avoid crashing if no data
    a.write(" \n")
    a.write("set yr [-800:0] \n")
    a.write("plot '"+datafile+"' u ($1/1.0e-19): 9  ti 'v_{ion,x}'   w l ls 5 ,\\\n")
    a.write("     '"+datafile+"' u ("+xmin+")  :(1) ti ''            w d         \n") # To avoid crashing if no data
    a.write(" \n")
  a.write("unset multiplot \n")
a.close()
#
#Make file
print "Graph script made: ",graphscript
if (does_gnuplot_exist==0): 
  os.system("gnuplot "+graphscript)
  print "Graph made: ",graphfiles
  if (open_eps<>""):  os.system(open_eps+" "+graphfiles)
