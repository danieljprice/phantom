#----------------------------------------------------------------------!
#                               N I C I L                              !
#           Non-Ideal mhd Coefficient and Ionisation Library           !
#         A Control script to allow the user to generate graphs        !
#                         of the relevant data.                        !
#                                                                      !
#                 Copyright (c) 2015-2019 James Wurster                !
#        See LICENCE file for usage and distribution conditions        !
#----------------------------------------------------------------------!
import os
import sys
#
# Defaults (these are useful values that may be temporarily altered for more useful graphing)
open_eps = ""  # Program used to open the .eps files; leave blank to prevent auto-opening (e.g. open_eps="open")
setxr    = True    # Use the updated default x-range; else, gnuplot will automatically select the range
setyr    = True    # Use the updated default y-range; else, gnuplot will automatically select the range
plotrhon = True    # Plot n_n & rho_n (else plot rho)
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
escomp = "0"
ierr     = 0
print "What would you like to graph?  The options are"
print "1: Graph the results from nicil_ex_eta"
print "2: Graph the results from nicil_ex_sph"
while (escomp<>"1" and escomp<>"2" and escomp<>""):
  if (ierr==1): print "That is an invalid entry. Please try again"
  escomp = raw_input("Enter the number for the desired graph [default=1]: ")
  if (escomp=="q" or escomp=="Q"): sys.exit()
  if (escomp==""): escomp="1"
  ierr     = 1
#
# Determine if plotting against number density or temperature
if (escomp=="1"):
  xaxis = ""
  ierr     = 0
  print " "
  print "What would you like on the horizontal axis?  The options are"
  print " 1: against number density with T = constant"
  print " 2: against number density with T = T(n) {barotropic EOS}"
  print " 3: against number density with T = T(n) {from ideal MHD core collapse}"
  print " 4: against number density with T = T(n) {from non-ideal MHD core collapse}"
  print " 5: against temperature with n = constant"
  print " 6: against temperature with n = n(T) {barotropic EOS}"
  print " 7: against temperature with n = n(T) {from ideal MHD core collapse}"
  print " 8: against temperature with n = n(T) {from non-ideal MHD core collapse}"
  print " 9: against zeta_cr"
  print "10: all of the above (assuming data exists)"
  print "11: B and T for options all available data (less zeta_cr and phase space)"
  print "12: rho-B phase space diagram"
  keep_trying = True
  while (keep_trying):
    cxaxis = raw_input("Enter the number for the desired axis [default=10]: ")
    if (cxaxis=="q" or cxaxis=="Q"): sys.exit()
    if (cxaxis==""):
      xaxis = "10"
      keep_trying = False
    else:
      try:
        isinstance(int(cxaxis),int)
        iinput = int(cxaxis)
        if (iinput < 13 and iinput > 0):
           keep_trying = False
           xaxis = str(iinput)
      except:
        keep_trying = True
else:
  xaxis = "1"
#
# Set file (character) IDs & define the plots
if (xaxis=="10"):
  imin = 1
  imax = 10
elif (xaxis=="11"):
  imin = 1
  imax = 5
elif (xaxis=="12"):
  imin = 0
  imax = 0
else:
  imin = int(xaxis)
  imax = imin + 1
graphscript = "Graphs/plot.gnuplot"
graphfiles  = ""
a=open(graphscript,"w")
for i in range(imin,imax):
  gtype = i
  xtype = 1
  if (escomp=="1"):
    if (gtype==1):
      datafile  = "data/eta_density.dat"
      graphfile = "Graphs/eta_density_Tcnst.eps"
    elif (gtype==2):
      datafile  = "data/eta_barotropic.dat"
      graphfile = "Graphs/eta_density_barotropic.eps"
    elif (gtype==3):
      datafile  = "data/eta_collapseIdeal.dat"
      graphfile = "Graphs/eta_density_collapseIdeal.eps"
    elif (gtype==4):
      datafile  = "data/eta_collapseNonideal.dat"
      graphfile = "Graphs/eta_density_collapseNonideal.eps"
    elif (gtype==5):
      xtype     = 2
      datafile  = "data/eta_temperature.dat"
      graphfile = "Graphs/eta_temperature_ncnst.eps"
    elif (gtype==6):
      xtype     = 2
      datafile  = "data/eta_barotropic.dat"
      graphfile = "Graphs/eta_temperature_barotropic.eps"
    elif (gtype==7):
      xtype     = 2
      datafile  = "data/eta_collapseIdeal.dat"
      graphfile = "Graphs/eta_temperature_collapseIdeal.eps"
    elif (gtype==8):
      xtype     = 2
      datafile  = "data/eta_collapseNonideal.dat"
      graphfile = "Graphs/eta_temperature_collapseNonideal.eps"
    elif (gtype==9):
      xtype     = 3
      datafile  = "data/eta_zeta.dat"
      graphfile = "Graphs/eta_zeta.eps"
  else:
    datafile  = "data/sph_density.dat"
    graphfile = "Graphs/sph_density.eps"
  if (os.path.isfile(datafile) and xaxis<>"11"):
    print "Plotting < "+graphfile+" > from file < "+datafile+" >."
    graphfiles = graphfiles+" "+graphfile
    #
    #
    # Make the graphing file
    a.write("reset \n")
    a.write("# Manually modifying the style types \n")
    if (gnuplot5):
      for i in range( 1,10): a.write("set style line  "+str(i)+" lc "+str(i   )+" lw 3 ps 2 \n")
      for i in range(11,20): a.write("set style line " +str(i)+" lc "+str(i-10)+" lw 1 ps 2 dt 2\n")
      for i in range(21,30): a.write("set style line " +str(i)+" lc "+str(i-20)+" lw 6 ps 2 \n")
    else:
      for i in range( 1,10): a.write("set style line  "+str(i)+" lt 1 lc "+str(i   )+" lw 3 ps 2 \n")
      for i in range(11,20): a.write("set style line " +str(i)+" lt 3 lc "+str(i-10)+" lw 1 ps 2 \n")
      for i in range(21,30): a.write("set style line " +str(i)+" lt 1 lc "+str(i-20)+" lw 6 ps 2 \n")
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
      if (xtype == 1):
        if (plotrhon):
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
          a.write("set xl  'log n(cm^{-3})' \n")
          if (setxr):
            a.write("set xr [1e2*m:1e22*m] \n")
            a.write("set xl 'log {/Symbol r} (g cm^{-3})' \n")
          xmin = "1.0e2*m"
          x = "1"
      elif(xtype==2):
        if (setxr): a.write("set xr [1e1:1e5] \n")
        a.write("set xl 'log T (K)' \n")
        xmin = "1.0e1"
        x = "2"
      elif(xtype==3):
        if (setxr): a.write("set xr [1e-32:1e-8] \n")
        a.write("set xl 'log {/Symbol z} (s^{-1})' \n")
        xmin = "1.0e-32"
        x = "33"
      a.write("set format x '%L' \n")
      a.write("set key top left \n")
      a.write(" \n")
      a.write("#Plot number densities \n")
      a.write("set log y \n")
      a.write("set yl 'log n (cm^{-3})' \n")
      if (setyr):
         if (xtype==3):
            a.write("set yr [1e-18:1e6] \n")
         else:
            a.write("set yr [1e-11:1e24] \n")
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
      if (setyr):
        if (xtype==3):
          a.write("set yr [1e-18:1e6] \n")
        else:
          a.write("set yr [1e-11:1e24] \n")
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
      if (setyr):
        if (xtype==3):
          a.write("set yr [1e-15:1e20] \n")
        else:
          a.write("set yr [5e-8:1e22] \n")
      a.write("set key top left \n")
      a.write("set format y ' %3L' \n")
      a.write("plot '"+datafile+"' u "+x+":   7  ti '{/Symbol s}_O'     w l ls 21 ,\\\n")
      a.write("     '"+datafile+"' u "+x+":($8>0.0? $8:0/0) ti '{/Symbol s}_H > 0' w l ls  3 ,\\\n")
      a.write("     '"+datafile+"' u "+x+":($8<0.0?-$8:0/0) ti '{/Symbol s}_H < 0' w l ls  4 ,\\\n")
      a.write("     '"+datafile+"' u "+x+":   9  ti '{/Symbol s}_P'     w l ls  2 ,\\\n")
      a.write("     '"+datafile+"' u ("+xmin+"):(1) ti ''               w d         \n") # To avoid crashing if no data
      a.write(" \n")
      a.write("#Plot resistivities \n")
      a.write("set key top left \n")
      a.write("set yl 'log {/Symbol h} (cm^2 s^{-1})' \n")
      if (setyr):
        if (xtype==3):
          a.write("set yr [1e6:1e25] \n")
        else:
          a.write("set yr [1e0:1e24] \n")
      a.write("plot '"+datafile+"' u "+x+":   4  ti '{/Symbol h}_{OR}'     w l ls 21 ,\\\n")
      a.write("     '"+datafile+"' u "+x+":($5>0.0? $5:0/0) ti '{/Symbol h}_{HE} > 0' w l ls  3 ,\\\n")
      a.write("     '"+datafile+"' u "+x+":($5<0.0?-$5:0/0) ti '{/Symbol h}_{HE} < 0' w l ls  4 ,\\\n")
      a.write("     '"+datafile+"' u "+x+":   6  ti '{/Symbol h}_{AD}'     w l ls  2 ,\\\n")
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
  elif (os.path.isfile(datafile) and xaxis=="11"):
    print "Plotting < "+graphfile+" > from file < "+datafile+" >."
    graphfiles = graphfiles+" "+graphfile
    #
    # Make the graphing file
    a.write("reset \n")
    a.write("# Manually modifying the style types \n")
    if (gnuplot5):
      a.write("set style line  1 lc 1 lw 3 ps 2 \n")
      a.write("set style line  3 lc 3 lw 3 ps 2 \n")
    else:
      a.write("set style line 1 lt 1 lc 1 lw 3 ps 2 \n")
      a.write("set style line 3 lt 1 lc 3 lw 3 ps 2 \n")
    a.write("  \n")
    a.write("m = 2.310*1.6726219e-24 \n")
    a.write("set terminal postscript eps enhanced colour size 10,5 font 'Times-Roman,20' \n")
    a.write("set output '"+graphfile+"' \n")
    a.write("  \n")
    a.write("set multiplot layout 1,2 \n")
    a.write("unset key \n")
    a.write("set log xy \n")
    a.write("set xr [1e2*m:1e22*m] \n")
    a.write("set xl 'log {/Symbol r} (g cm^{-3})' \n")
    a.write("set format x ' %3L' \n")
    xmin = "1.0e2*m"
    x = "13"
    a.write("#Plot Magnetic field \n")
    a.write("set yl 'log B (G)' \n")
    a.write("set yr [1e-8:1e6] \n")
    a.write("set format y ' %3L' \n")
    a.write("plot '"+datafile+"' u 1:3 w l ls 1 \n")
    a.write("#Plot Temperature \n")
    a.write("set yl 'log T (K)' \n")
    a.write("set yr [9:1e5] \n")
    a.write("set format y ' %3L' \n")
    a.write("plot '"+datafile+"' u 1:2 w l ls 1 \n")
    a.write("unset multiplot \n")
datafile = "data/eta_phase.dat"
if (os.path.isfile(datafile) and xaxis=="12"):
  graphfile  = "Graphs/eta_phase.eps"
  graphfiles = graphfiles+" "+graphfile
  print "Plotting < "+graphfile+" > from file < "+datafile+" >."
  #
  # Make the graphing file
  a.write("reset \n")
  a.write("# Manually modifying the style types \n")
  a.write("set terminal postscript eps enhanced colour size 7,5 font 'Times-Roman,35'\n")
  if (gnuplot5):
    a.write("set style line  1 lc 1 lw 3 lt 7 ps 2 \n") # remove the lt 7?
    a.write("set style line  2 lc 2 lw 3 lt 7 ps 2 \n")
    a.write("set style line  3 lc 3 lw 3 lt 7 ps 2 \n")
    a.write("set style line  5 lc 9 lw 3 lt 7 ps 2 \n")
  else:
    a.write("set style line  1 lc 1 lw 3 lt 7 ps 2 \n")
    a.write("set style line  2 lc 2 lw 3 lt 7 ps 2 \n")
    a.write("set style line  3 lc 3 lw 3 lt 7 ps 2 \n")
    a.write("set style line  5 lc 5 lw 3 lt 7 ps 2 \n")
  a.write("set terminal postscript eps enhanced colour size 7,5 font 'Times-Roman,35'\n")
  a.write("set output '"+graphfile+"' \n")
  a.write("set xr [1e-22:1] \n")
  a.write("set yr [1e-7:1e5] \n")
  a.write("set log xy  \n")
  a.write("set yl 'log(B) [G]' \n")
  a.write("set xl 'log({/Symbol r}) [g cm^{-3}]' \n")
  a.write("set format x '%L' \n")
  a.write("set format y '%L' \n")
  a.write("set label 'Ambipolar'     at 1e-19,1e3 front centre textcolor rgb 'white' \n")
  a.write("set label '(dissipative)' at 1e-19,1e2 front centre textcolor rgb 'white' \n")
  a.write("set label 'Hall > 0'      at 5e-15,1e-5 front centre  \n")
  a.write("set label '(dispersive)'  at 5e-15,1e-6 front centre  \n")
  a.write("set label 'Hall < 0'      at 5e-13,1e1 front centre  \n")
  a.write("set label '(dispersive)'  at 5e-13,1e0 front centre  \n")
  a.write("set label 'Ohmic'         at 1e-5,1e-4 front centre textcolor rgb 'white' \n")
  a.write("set label '(dissipative)' at 1e-5,1e-5 front centre textcolor rgb 'white' \n")
  a.write("unset key \n")
  a.write("plot '"+datafile+"' u 1:($7== 1?$2:0/0) ls 1 ti 'Ohmic',      ")
  a.write("     '"+datafile+"' u 1:($7== 2?$2:0/0) ls 2 ti 'Hall > 0',   ")
  a.write("     '"+datafile+"' u 1:($7==-2?$2:0/0) ls 5 ti 'Hall < 0',   ")
  a.write("     '"+datafile+"' u 1:($7== 3?$2:0/0) ls 3 ti 'Ambipolar' \n")
a.close()
#
#Make file
print "Graph script made: ",graphscript
if (does_gnuplot_exist==0): 
  os.system("gnuplot "+graphscript)
  if (graphfiles <>""):
    print "Graph made: ",graphfiles
    if (open_eps<>""):  os.system(open_eps+" "+graphfiles)
