#!/Users/dprice/anaconda3/bin/python3
#---------------------------------------------------------------
#
#  Daniel Price's automatic smoothing kernel generation script
#  Prints kernels and all related functions in suitable
#  manner for direct incorporation into SPH codes
#
#---------------------------------------------------------------
from __future__ import division
from sympy import *
q, x, y = symbols('q x y')

##############################################
#                                            #
# various functions to actually perform the  #
# symbolic calculations                      #
#                                            #
##############################################

#---------------------------------------------
# function to get normalisation constants for
# kernel in 1, 2 and 3 dimensions
#---------------------------------------------
def getnorm(w,R):
    c1D = sympify(1)/(2*integrate(w,(q,0,R)))
    c2D = sympify(1)/(integrate(2*pi*q*w,(q,0,R)))
    c3D = sympify(1)/(integrate(4*pi*q*q*w,(q,0,R)))
    return (c1D, c2D, c3D)

#-----------------------------------------------
# work out the integration constant by matching
# the different parts of the piecewise function
#-----------------------------------------------
def intconst(g):
    if isinstance(g, Piecewise):
       garg = list(g.args)
       for i, (e, c) in reversed(list(enumerate(g.args))):
           if i < len(g.args) - 1:
              (ep, cp) = garg[i+1]
              s = "%s" %c
              qval = sympify(s.split()[2])
              ge = simplify(e + (ep.subs(q,qval) - e.subs(q,qval)))
              garg[i] = (ge,c)
       tuple(garg)
       g = Piecewise(*garg)
    return g

#----------------------------------------------
# function to get force softening kernel and
# related function from the density kernel
#----------------------------------------------
def getkernelfuncs(w,R):
    dw  = diff(w,q)
    d2w = diff(dw,q)
    c1D, c2D, c3D = getnorm(w,R)
    #
    #--force softening function
    #
    fsoft = piecewise_fold(4*pi*c3D*integrate(f*q*q,q)/q**2)
    farg = list(fsoft.args)
    lastarg = len(fsoft.args) - 1
    farg[lastarg] = (sympify(1)/q**2,fsoft.args[lastarg].cond)
    #
    #--work out the integration constant for the force softening
    #  by matching the different parts of the piecewise function
    #
    if isinstance(fsoft, Piecewise):
       for i, (e, c) in reversed(list(enumerate(fsoft.args))):
           if i < lastarg:
              (ep, cp) = farg[i+1]
              s = "%s" %c
              qval = sympify(s.split()[2])
              fe = simplify(e + qval**2*(ep.subs(q,qval) - e.subs(q,qval))/q**2)
              farg[i] = (fe,c)
    tuple(farg)
    fsoft = Piecewise(*farg)
    #
    #--potential function
    #
    pot = integrate(fsoft,q)
    #
    #--work out the integration constant for the potential
    #
    parg = list(pot.args)
    if isinstance(pot, Piecewise):
       for i, (e, c) in reversed(list(enumerate(pot.args))):
           if i < len(pot.args) - 1:
              (ep, cp) = parg[i+1]
              s = "%s" %c
              qval = s.split()[2]
              pote = simplify(e + (ep.subs(q,qval) - e.subs(q,qval)))
              parg[i] = (pote,c)
    tuple(parg)
    pot = Piecewise(*parg)
    #
    #--derivative of potential with respect to h
    #
    dpotdh = pot
    parg = list(pot.args)
    if isinstance(pot, Piecewise):
       for i, (e, c) in enumerate(pot.args):
           ep = simplify(-e - q*diff(e,q))
           parg[i] = (ep, c)
       tuple(parg)
       dpotdh = Piecewise(*parg)

    return (dw, d2w, c1D, c2D, c3D, fsoft, pot, dpotdh)

#---------------------------------------------
# function to get the variance of the kernel
#---------------------------------------------
def getvar(w,R):
    c1D, c2D, c3D = getnorm(w,R)
    var = (integrate(c1D*q*q*w,(q,0,R)),
           integrate(c2D*q*q*2*pi*q*w,(q,0,R)),
           integrate(c3D*q*q*4*pi*q*q*w,(q,0,R)))
    varcubic = (sympify(1)/6, sympify(31)/49, sympify(9)/10)
    relvar = (var[0]/varcubic[0], var[1]/varcubic[1], var[2]/varcubic[2])
    reldev = (sqrt(1.0*relvar[0]),sqrt(1.0*relvar[1]),sqrt(1.0*relvar[2]))
    return (var,relvar,reldev)

#-------------------------------------------------------
# function to get the standard deviation of the kernel
# scaled relative to the cubic spline
#-------------------------------------------------------
def getreldev(w,R):
    var, relvar, reldev = getvar(w,R)
    return (reldev)

#--------------------------------------------------------
# Functions to return kernels that are constructed from
# other kernels in some way
#--------------------------------------------------------
def intkernel(wref,R):
    f, name = wref(R)
    g = piecewise_fold(integrate(-q*f,q))
    g = intconst(g)
    name = "Integrated %s" %(name)
    return(g,name)

def intkernel2(wref,R):
    f, name = wref(R)
    g = piecewise_fold(integrate(-q*f,q))
    g = intconst(g)
    g = piecewise_fold(integrate(-q*g,q))
    g = intconst(g)
    name = "Twice-integrated %s" %(name)
    return(g,name)

def intkernel3(wref,R):
    f, name = wref(R)
    g = piecewise_fold(integrate(-q*f,q))
    g = intconst(g)
    g = piecewise_fold(integrate(-q*g,q))
    g = intconst(g)
    g = piecewise_fold(integrate(-q*g,q))
    g = intconst(g)
    name = "Triple-integrated %s" %(name)
    return(g,name)

def doublehump(wref,R):
    f, name = wref(R)
    g = piecewise_fold(f*q*q)
    name = "Double-hump %s" %(name)
    return(g,name)

def doublehump3(wref,R):
    f, name = wref(R)
    g = piecewise_fold(f*q*q*q)
    name = "Double-hump-on-steroids %s" %(name)
    return(g,name)

def doublehump5(wref,R):
    f, name = wref(R)
    g = piecewise_fold(f*q*q*q*q*q)
    name = "Double-hump-on-overdrive %s" %(name)
    return(g,name)

##############################################
#                                            #
#  various output functions to print kernel  #
#  information in different ways             #
#                                            #
##############################################

#-------------------------------------------------------
# function to print the variance and standard deviation
# of a kernel, and to print these relative to the cubic
#-------------------------------------------------------
def printvariances(w,R):
    var, relvar, reldev = getvar(w,R)
    print ("\nVariance of kernel in 1, 2, 3D:")
    print (var[0],var[1],var[2])
    print ("\nVariance and standard dev relative to cubic:")
    print (relvar[0],relvar[1],relvar[2])
    print (reldev[0],reldev[1],reldev[2])
    print ("\nKernel radius required to get same std. dev as cubic:")
    print (2/reldev[0],2/reldev[1],2/reldev[2])
    print ("\neta = 1.2 is equivalent to:")
    print (1.2/reldev[0],1.2/reldev[1],1.2/reldev[2])
    return

#-----------------------------------------------------------
# function to print basic kernel information to the screen
#-----------------------------------------------------------
def printkernel(w,R):
    dw, d2w, c1D, c2D, c3D, fsoft, pot, dpotdh = getkernelfuncs(w,R)
    print ("\n%s W:" %name)
    print (w)
    print ("\nFirst derivative:")
    print (dw)
    print ("\n2nd derivative:")
    print (d2w)
    print ("\nnormalisation:")
    print ("[ %s, %s, %s ]" %(c1D,c2D,c3D))
    print ("\n3D normalisation of artificial viscosity term:")
    avnorm = -sympify(2)*pi/15*c3D*integrate(q*q*q*dw,(q,0,R))
    print (avnorm)
    print ("\n2D normalisation of artificial viscosity term:")
    avnorm = -pi/8*c2D*integrate(q*q*dw,(q,0,R))
    print (avnorm)
    printvariances(w,R)
    return

#-------------------------------------------------------------
# print start of a LaTeX table containing kernel information
#-------------------------------------------------------------
def printheader_latex():
    print ("\\begin{tabular}{|l|l|l|l|l|l|l|l|}\n")
    print ("\\hline\nName & Functional form & C$_{1D}$ & C$_{2D}$ & C$_{3D}$ & $\sigma^2_{1D}$ & $\sigma^2_{2D}$ & $\sigma^2_{3D}$\\\\ \n")

#-----------------------------------------------------------
# print end of a LaTeX table containing kernel information
#-----------------------------------------------------------
def printfooter_latex():
    print ("\\hline\\end{tabular}\n")

#---------------------------------------------------------------
# print contents of a LaTeX table containing kernel information
#---------------------------------------------------------------
def printkernel_latex(w,R):
    c1D, c2D, c3D = getnorm(w,R)
    var, relvar, reldev = getvar(w,R)
    print ("\\hline\n%s & $" %fmttex(name))
    print (latex(w))
    print ("$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\\\" %(latex(c1D),latex(c2D),latex(c3D),latex(var[0]),latex(var[1]),latex(var[2])))
    return

#--------------------------------
# format names for LaTeX output
#--------------------------------
def fmttex(s):
    import re
    s = re.sub("\^\S+","$\g<0>$", s)
    s = re.sub("\_\S+","$\g<0>$", s)
    return s

#-------------------------------------------------------------------------------
# utility to format output of real numbers correctly for Fortran floating point
#-------------------------------------------------------------------------------
def fmt(e):
    import re
    s = "%s" %e
    # add decimal points to numbers, but not if powers like q**2 (or *2 with one digit)
    # note that \g<0> gives first matching argument in the regex
    # rules are: (not **n)(not 0.123)(match ab0123) or (not *n with one digit)
    s = re.sub("((?!\*\*\d+)(?!\D\D\d+\.)\D\D\d+)|((!?\*\d+)\D\d+)|(/\d+)|((?!^\.\d+)^\d+)|((?!^-\d+\.)^-\d+)","\g<0>.", s)

    # replace 15*x with 15.*x as long as it is not **15*x
    s = re.sub("(?!\*\d+)(\D\d+)\*","\g<1>.*", s)

    f = sympify(s)
    #
    # expand if it makes it shorter
    #
    h = "%s" %(expand(f))
    #f = h
    if (len(h) <= len(s)):
       f = h
    g = "%s" %simplify(f)

    # replace 1.4000000 with 1.4
    g = re.sub("(\.[1-9]*)(0+)(\D|$)","\g<1>\g<3>", g)

    # only return simplify-ed strings if no fully expanded floats 0.345242545..
    if re.search("(\.\d\d\d\d\d+)",g):
       return s
    else:
       return g

#------------------------------------------------------------------------
# extended version of above, replacing q**2 with q2, q**3 with q2*q etc.
# and getting rid of excess zeros after decimal point, e.g. 7.0000->7
#------------------------------------------------------------------------
def fmte(e,useqsub,useodd):
    import re
    s = "%s" %fmt(e)
    #fs = ""
    #for arg in (split(s,' ')):
    #    fs = fs+arg
    f = sympify(s)
    g = "%s" %simplify(f)
    if len(g) <= len(s) + 1:
       s = g
    if (useqsub):
       s = re.sub("q\*\*12","q6*q6", s)
       s = re.sub("q\*\*11","q6*q4*q", s)
       s = re.sub("q\*\*10","q6*q4", s)
       s = re.sub("q\*\*8","q8", s)
       if (useodd):
          s = re.sub("q\*\*9","q9", s)
          s = re.sub("q\*\*7","q7", s)
          s = re.sub("q\*\*5","q5", s)
          s = re.sub("q\*\*3","q3", s)
       else:
          s = re.sub("q\*\*9","q8*q", s)
          s = re.sub("q\*\*7","q6*q", s)
          s = re.sub("q\*\*5","q4*q", s)
          s = re.sub("q\*\*3","q2*q", s)
       s = re.sub("q\*\*6","q6", s)
       s = re.sub("q\*\*4","q4", s)
    s = re.sub("q\*\*2","q2", s)
    s = re.sub("q\*\*\(-2\.\)","1./q2",s)
    s = re.sub("q\*\*-2\.0","1./q2",s)
    s = re.sub("q\*\*\(-2\.0\)","1./q2",s)
    # remove excess zeros after decimal place
    # handles case of 3.0*q4 -> 3.*q4, zeros must be followed by non-digit or end of line
    s = re.sub("(\.0+)(\D+|$)",".\g<2>", s)
    return s

#-----------------------------------------------------
# wrap output to 72 characters for Fortran 77 output
#-----------------------------------------------------
def wrapit(s,indent):
    if len(s) > 72:
       pos = s.rfind(" ", 6, 72)
       if pos == -1:
           pos = 72
       hunk = s[:pos]
       rest = s[pos:].lstrip()
       s = "%s &\n" %(hunk) + " "*indent
       while len(rest) > 0:
           pos = rest.rfind(" ", 0, 66)
           if pos == -1 or len(rest) < 66:
               pos = 66
           hunk = rest[:pos]
           rest = rest[pos:].lstrip()
           if len(rest) > 0:
              s = "%s%s &\n" %(s,hunk) + " "*indent
           else:
              s = "%s%s" % (s,hunk)
    return s

#------------------------------------
# wrap and indent Fortran 77 output
#------------------------------------
def wrapf77(s,indent):
    maxl = (72 - indent)
    if len(s) > maxl:
       pos = s.rfind(" ", 6, maxl)
       if pos == -1:
           pos = maxl
       hunk = s[:pos]
       rest = s[pos:].lstrip()
       s = "%s \n" %(hunk) + "     &"+ " "*(indent-6)
       while len(rest) > 0:
           pos = rest.rfind(" ", 0, (maxl-6))
           if pos == -1 or len(rest) < (maxl-6):
               pos = maxl - 6
           hunk = rest[:pos]
           rest = rest[pos:].lstrip()
           if len(rest) > 0:
              s = "%s%s\n" %(s,hunk) + "     &"+" "*(indent-6)
           else:
              s = "%s%s" % (s,hunk)
    return s

#---------------------------------------------------------
# wrappers for above routines, specific to output formats
# these define the line length and the indent
#---------------------------------------------------------

def fmtp(e):
    s = "%s" %fmte(e,True,False)
    s = wrapit(s,17)
    return s

def fmtn(e):
    s = "%s" %fmte(e,True,False)
    s = wrapit(s,25)
    return s

def fmts(e):
    s = "%s" %fmte(e,True,True)
    s = wrapf77(s,18)
    return s

def stripcond(e):
    import re
    s = "%s" %fmt(e,True,True)
    s = re.sub("q|<|>|\s","",s)
    return s

#-------------------------------
# print FORTRAN77 comment line
#-------------------------------
def printc(s):
    print ("c\nc--%s\nc" %(s))
    return s

#---------------------------------
# print kernel code for ndspmhd
#---------------------------------
def printkernel_ndspmhd(w,R,name):
    useoddq = False
    dw, d2w, c1D, c2D, c3D, fsoft, pot, dpotdh = getkernelfuncs(w,R)
    print ("!")
    print ("!--%s (auto-generated by kernels.py)" %name)
    print ("!")
    print ("    kernellabel = '%s' \n" %name)
    print ("    radkern = %.1f" %(R))
    print ("    radkern2 = radkern*radkern")
    print ("    dq2table = radkern2/real(ikern)")
    print ("    select case(ndim)")
    print ("      case(1)")
    print ("         cnormk = %s" %fmt(c1D))
    print ("      case(2)")
    print ("         cnormk = %s" %fmt(c2D))
    print ("      case(3)")
    print ("         cnormk = %s" %fmt(c3D))
    print ("    end select")
    print ("    do i=0,ikern")
    print ("       q2 = i*dq2table")
    print ("       q4 = q2*q2")
    print ("       q6 = q4*q2")
    print ("       q8 = q4*q4")
    print ("       q = sqrt(q2)")
    if isinstance(w, Piecewise):
       for i, (e, c) in enumerate(w.args):
           (de, dc) = dw.args[i]
           (d2e, d2c) = d2w.args[i]
           (fe,fc) = fsoft.args[i]
           (pe,pc) = pot.args[i]
           (pdhe,pdhc) = dpotdh.args[i]
           if i == 0:
              print ("       if (%s) then" %fmt(c))
           elif i == len(w.args)-1 and c == True:
              print ("       else")
           else:
              print ("       elseif (%s) then" %fmt(c))
           print ("          wkern(i)     = %s " %fmtn(e))
           print ("          grwkern(i)   = %s " %fmtn(de))
           print ("          grgrwkern(i) = %s " %fmtn(d2e))
           print ("          fsoft(i)     = %s " %fmtn(fe))
           print ("          potensoft(i) = %s " %fmtn(pe))
           print ("          dphidh(i)    = %s " %fmtn(pdhe))
       print ("       endif")
    else:
       print (w)
    print ("    enddo\n")

#---------------------------------
# print kernel code for sphNG
#---------------------------------
def printkernel_sphNG(w,R,name):
    import datetime
    dw, d2w, c1D, c2D, c3D, fsoft, pot, dpotdh = getkernelfuncs(w,R)
    print ("      SUBROUTINE ktable")
    print ("c*********************************************************")
    print ("c  This subroutine builds a table for the kernel,")
    print ("c  the gradient of the kernel, the mass fraction,")
    print ("c  and the potential energy.")
    print ("c  The entry is v**2.")
    print ("c")
    print ("c  DO NOT EDIT: AUTO-GENERATED by kernels.py")
    print ("c  KERNEL NAME: %s " %name)
    print ("c  AUTHOR: kernels.py, by Daniel Price")
    print ("c  GENERATED:",datetime.datetime.now())
    print ("c")
    print ("c*********************************************************")
    print ("      IMPLICIT NONE ! because life is worth living")
    print ("      INCLUDE 'idim'\n")
    print ("      REAL*8 sum, v2max, q, q2, q3, q4, q5, q6, q7, q8, q9")
    print ("      INTEGER i")
    print ("\n      INCLUDE 'COMMONS/physcon'")
    print ("      INCLUDE 'COMMONS/kerne'")
    print ("      INCLUDE 'COMMONS/table'")
    print ("      INCLUDE 'COMMONS/logun'")
    print ("      INCLUDE 'COMMONS/debug'")
    printc("Allow for tracing flow")
    print ("      IF (itrace.EQ.'all') WRITE(iprint, 99001)")
    print ("99001 FORMAT (' entry subroutine ktable')")
    printc("Maximum interaction length and step size")
    print ("      radkernel = %.1f" %(R))
    if isinstance(w, Piecewise):
       for i, (e, c) in enumerate(w.args):
           if (c != True and i < 2):
              print ("      part%ikernel = %s" %(i+1,stripcond(c)))
    print ("      v2max = radkernel*radkernel")
    print ("      dvtable = v2max/itable")
    print ("      ddvtable = itable/v2max")
    printc("Build tables")
    print ("      DO i=0,itable")
    print ("         q2 = i*dvtable")
    print ("         q = sqrt(q2)")
    print ("         q3 = q*q2")
    print ("         q4 = q*q3")
    print ("         q5 = q*q4")
    print ("         q6 = q*q5")
    print ("         q7 = q*q6")
    print ("         q8 = q*q7")
    print ("         q9 = q*q8")
    if isinstance(w, Piecewise):
       for i, (e, c) in enumerate(w.args):
           (de, dc) = dw.args[i]
           (d2e, d2c) = d2w.args[i]
           (fe,fc) = fsoft.args[i]
           (pe,pc) = pot.args[i]
           (pdhe,pdhc) = dpotdh.args[i]
           if i == 0:
              print ("         IF (%s) THEN" %fmt(c))
           elif i == len(w.args)-1 and c == True:
              print ("         ELSE")
           else:
              print ("         ELSEIF (%s) THEN" %fmt(c))
           print ("            sum = %s" %fmts(e))
           print ("            wij(i) = sum")
           print ("            sum = %s" %fmts(de))
           print ("            grwij(i) = sum")
           print ("            sum = %s" %fmts(q*q*fe))
           print ("            fmass(i) = sum")
           print ("            sum = %s" %fmts(pe))
           print ("            fpoten(i) = sum")
           print ("            sum = %s" %fmts(-pdhe))
           print ("            dphidh(i) = sum")
       print ("         ENDIF")
    print ("      ENDDO")
    printc("Normalisation constant")
    print ("      cnormk = %s" %fmt(c3D))
    print ("      selfnormkernel = %s" %fmt(w.subs(q,0)))
    print ("      part1potenkernel = 0.0 ! term already included in fpoten above")
    print ("      part2potenkernel = 0.0 ! see above")
    #--double hump normalisation
    wdrag = piecewise_fold(w*q*q)
    c3Ddrag = sympify(1)/(integrate(4*pi*q*q*wdrag,(q,0,R)))
    printc("For dust/gas drag, need double humped kernel")
    print ("      doublehumpnormk = %s" %fmt(c3Ddrag))
    print ("\n      RETURN")
    print ("      END")

#-------------------------------------------------
# print q4 = q2*q2 and similar definitions
# e.g. if the string q4 is found in the function
#-------------------------------------------------
def print_defs(indent,*args):
    import re
    gotq4 = False
    # look for q4, q6, q8 in function string
    for i in (4,6,8):
       str = "q%i" %i
       doPrint = False
       # match in any functions about to be printed
       for arg in args:
          if (re.search(str,arg)):
             doPrint = True
       # print definition
       if (doPrint):
          print (" "*indent+"q%i = q%i*q2" %(i,i-2))

#----------------------------------------
# print code declarations of q4, q6 etc.
#----------------------------------------
def print_decl(w):
    import re
    str = ""
    for qstr in ("q4","q6","q8"):
        printQ = False
        for i, (e, c) in enumerate(w.args):
            if (re.search(qstr,fmtp(e))):
               printQ = True
        if (printQ):
           if (len(str) > 0):
              str = str+", "+qstr
           else:
              str = qstr
    if (len(str) > 0):
       print (" real :: %s\n" %(str))
    else:
       print ("")

#---------------------------------
# print kernel code for Phantom
#---------------------------------
def printkernel_phantom(w,R,name):
    import datetime
    dw, d2w, c1D, c2D, c3D, fsoft, pot, dpotdh = getkernelfuncs(w,R)
    w0 = w.subs(q,0)
    dpotdh0 = dpotdh.subs(q,0)
    #
    #--double-hump kernel used in drag routines, with normalisation
    #
    wdrag = piecewise_fold(w*q*q)
    c3Ddrag = sympify(1)/(integrate(4*pi*q*q*wdrag,(q,0,R)))
    lb = "!"+"-"*62
    print ("!--------------------------------------------------------------------------!")
    print ("! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !")
    print ("! Copyright (c) 2007-2014 The Authors (see AUTHORS)                        !")
    print ("! See LICENCE file for usage and distribution conditions                   !")
    print ("! http://users.monash.edu.au/~dprice/phantom                               !")
    print ("!--------------------------------------------------------------------------!")
    print ("!+")
    print ("!  MODULE: kernel")
    print ("!")
    print ("!  DESCRIPTION:")
    print ("!   This module implements the %s kernel" %name)
    print ("!   DO NOT EDIT - auto-generated by kernels.py")
    print ("!")
    print ("!  REFERENCES: None")
    print ("!")
    print ("!  OWNER: Daniel Price")
    print ("!")
    print ("!  $Id:$")
    print ("!")
    print ("!  RUNTIME PARAMETERS: None")
    print ("!")
    print ("!  DEPENDENCIES: physcon")
    print ("!")
    print ("!  GENERATED:",datetime.datetime.now())
    print ("!+")
    print ("!--------------------------------------------------------------------------")
    print ("module kernel")
    print (" use physcon, only:pi")
    print (" implicit none")
    print (" character(len=%i), public :: kernelname = '%s'" %(len(name),name))
    print (" real, parameter, public  :: radkern  = %s" %fmt(R))
    print (" real, parameter, public  :: radkern2 = %s" %fmt(R*R))
    print (" real, parameter, public  :: cnormk = %s" %fmt(c3D))
    print (" real, parameter, public  :: wab0 = %s, gradh0 = -3.*wab0" %fmt(w0))
    print (" real, parameter, public  :: dphidh0 = %s" %fmtp(dpotdh0))
    print (" real, parameter, public  :: cnormk_drag = %s " %fmt(c3Ddrag))
    var, relvar, reldev = getvar(w,R)
    print (" real, parameter, public  :: hfact_default = %.1f " %(1.2/reldev[2]))
    #print " real, parameter, public  :: hfact_default = %s " %fmt(reldev[2])
    print ("\ncontains\n")
    print ("pure subroutine get_kernel(q2,q,wkern,grkern)")
    print (" real, intent(in)  :: q2,q")
    print (" real, intent(out) :: wkern,grkern")
    print_decl(w)
    print (" !--%s" %name)
    if isinstance(w, Piecewise):
       for i, (e, c) in enumerate(w.args):
           (de, dc) = dw.args[i]
           if i == 0:
              print (" if (%s) then" %fmt(c))
           elif i == len(w.args)-1 and c == True:
              print (" else")
           else:
              print (" elseif (%s) then" %fmt(c))
           print_defs(4,fmtp(e),fmtp(de))
           print ("    wkern  = %s" %fmtp(e))
           print ("    grkern = %s" %fmtp(de))
       print (" endif")
    else:
       print (w)
    print ("\nend subroutine get_kernel\n")
    print ("pure elemental real function wkern(q2,q)")
    print (" real, intent(in) :: q2,q")
    print_decl(w)
    if isinstance(w, Piecewise):
       for i, (e, c) in enumerate(w.args):
           if i == 0:
              print (" if (%s) then" %fmt(c))
           elif i == len(w.args)-1 and c == True:
              print (" else")
           else:
              print (" elseif (%s) then" %fmt(c))
           print_defs(4,fmtp(e))
           print ("    wkern = %s" %fmtp(e))
       print (" endif")
    else:
       print_defs(4,w)
       print ("    wkern = %s" %w)
    print ("\nend function wkern\n")
    print ("pure elemental real function grkern(q2,q)")
    print (" real, intent(in) :: q2,q")
    print_decl(dw)
    if isinstance(dw, Piecewise):
       for i, (e, c) in enumerate(dw.args):
           if i == 0:
              print (" if (%s) then" %fmt(c))
           elif i == len(w.args)-1 and c == True:
              print (" else")
           else:
              print (" elseif (%s) then" %fmt(c))
           print_defs(4,fmtp(e))
           print ("    grkern = %s" %fmtp(e))
       print (" endif")
    else:
       print_defs(4,fmtp(dw))
       print ("    grkern = %s " %fmtp(dw))
    print ("\nend function grkern\n")
    print ("pure subroutine get_kernel_grav1(q2,q,wkern,grkern,dphidh)")
    print (" real, intent(in)  :: q2,q")
    print (" real, intent(out) :: wkern,grkern,dphidh")
    print_decl(dpotdh)
    if isinstance(w, Piecewise):
       for i, (e, c) in enumerate(w.args):
           (de, dc) = dw.args[i]
           (dphie, dphic) = dpotdh.args[i]
           if i == 0:
              print (" if (%s) then" %fmt(c))
           elif i == len(w.args)-1 and c == True:
              print (" else")
           else:
              print (" elseif (%s) then" %fmt(c))
           print_defs(4,fmtp(e),fmtp(de),fmtp(dphie))
           print ("    wkern  = %s" %fmtp(e))
           print ("    grkern = %s" %fmtp(de))
           print ("    dphidh = %s" %fmtp(dphie))
       print (" endif")
    else:
       print_defs(4,fmtp(w),fmtp(dw),fmtp(dpotdh))
       print ("    wkern  = %s" %fmtp(w))
       print ("    grkern = %s" %fmtp(dw))
       print ("    dphidh = %s" %fmtp(dpotdh))
    print ("\nend subroutine get_kernel_grav1\n")
#    print "pure subroutine get_kernel_grav2(q2,q,grkern,potensoft,fsoft)"
#    print " real, intent(in)  :: q2,q"
#    print " real, intent(out) :: grkern,potensoft,fsoft\n"
#    if isinstance(dw, Piecewise):
#       for i, (de, c) in enumerate(dw.args):
#           (pote, potc) = pot.args[i]
#          (fe, fc) = fsoft.args[i]
#           if i == 0:
#              print " if (%s) then" %fmt(c)
#           elif i == len(dw.args)-1 and c == True:
#              print " else"
#           else:
#              print " elseif (%s) then" %fmt(c)
#           print "    grkern    = %s" %fmtp(de)
#           print "    potensoft = %s" %fmtp(pote)
#           print "    fsoft     = %s" %fmtp(fe)
#       print " endif"
#    else:
#       print "    wkern     = %s" %fmtp(w)
#       print "    grkern    = %s" %fmtp(dw)
#       print "    potensoft = %s" %fmtp(pot)
#       print "    fsoft     = %s" %fmtp(fsoft)
#    print "\nend subroutine get_kernel_grav2\n"
    print ("pure subroutine kernel_softening(q2,q,potensoft,fsoft)")
    print (" real, intent(in)  :: q2,q")
    print (" real, intent(out) :: potensoft,fsoft")
    print_decl(pot)
    if isinstance(dw, Piecewise):
       for i, (de, c) in enumerate(dw.args):
           (pote, potc) = pot.args[i]
           (fe, fc) = fsoft.args[i]
           if i == 0:
              print (" if (%s) then" %fmt(c))
           elif i == len(dw.args)-1 and c == True:
              print (" else")
           else:
              print (" elseif (%s) then" %fmt(c))
           print_defs(4,fmtp(pote),fmtp(fe))
           print ("    potensoft = %s" %fmtp(pote))
           print ("    fsoft     = %s" %fmtp(fe))
       print (" endif")
    else:
       print ("    potensoft = %s" %fmtp(pot))
       print ("    fsoft     = %s" %fmtp(fsoft))
    print ("\nend subroutine kernel_softening\n")
    print ("!------------------------------------------")
    print ("! double-humped version of the kernel for")
    print ("! use in drag force calculations")
    print ("!------------------------------------------")
    print ("pure elemental real function wkern_drag(q2,q)")
    print (" real, intent(in) :: q2,q")
    print_decl(wdrag)
    print (" !--double hump %s kernel" %name)
    if isinstance(wdrag, Piecewise):
       for i, (e, c) in enumerate(wdrag.args):
           if i == 0:
              print (" if (%s) then" %fmt(c))
           elif i == len(wdrag.args)-1 and c == True:
              print (" else")
           else:
              print (" elseif (%s) then" %fmt(c))
           print_defs(4,fmtp(e))
           print ("    wkern_drag = %s" %fmtp(e))
       print (" endif")
    else:
       print_defs(4,fmtp(wdrag))
       print ("    wkern_drag = %s" %fmtp(wdrag))
    print ("\nend function wkern_drag\n")
    print ("end module kernel")

def printalltex():
    R = sympify(2)
    printheader_latex()
    for x in m4, m5, m6, w2_1D, w4_1D, w6_1D, w2, w4, w6, intm4, intm5, intm6, int2m4, int3m4, f6:
        f, name = x(R)
        printkernel_latex(f,R)
    printfooter_latex()

def print_stddevs():
    for x in m4, m5, m6, w2_1D, w4_1D, w6_1D, w2, w4, w6, intm4, intm5, intm6, int2m4, int3m4, f6:
       f, name = x(R)
       reldev = getreldev(f,R)
       print (x.__name__,1.0/reldev[0],1.0/reldev[1],1.0/reldev[2])


#####################################
#  KERNEL DEFINITIONS START BELOW   #
#####################################

#-------------------------
# B-spline kernels
#-------------------------
def m4(R):
    f = Piecewise((sympify(1)/4*(R-q)**3 - (R/2 - q)**3,q < R/2), (sympify(1)/4*(R-q)**3, q < R), (0, True))
    return(f,'M_4 cubic')

def m5(R):
    term1 = sympify((R-q)**4)
    term2 = -5*(sympify(3)/5*R - q)**4
    term3 = 10*(sympify(1)/5*R - q)**4
    f = Piecewise((term1 + term2 + term3,q < sympify(1)/5*R), (term1 + term2, q < sympify(3)/5*R), (term1, q < R), (0, True))
    return(f,'M_5 quartic')

def m6(R):
    f = symbols('f',cls=Function)
    term1 = sympify((R-q)**5)
    term2 = -6*(sympify(2)/3*R - q)**5
    term3 = 15*(sympify(1)/3*R - q)**5
    f = Piecewise((term1 + term2 + term3,q < sympify(1)/3*R), (term1 + term2, q < sympify(2)/3*R), (term1, q < R), (0, True))
    return(f,'M_6 quintic')

#-------------------------
# Wendland kernels in 1D
#-------------------------
def w2_1D(R):
    f = Piecewise(((1 - q/R)**3*(1 + 3*q/R),q < R), (0, True))
    return(f,'Wendland 1D C^2')

def w4_1D(R):
    f = Piecewise(((1 - q/R)**5*(1 + 5*q/R + 8*(q/R)**2),q < R), (0, True))
    return(f,'Wendland 1D C^4')

def w6_1D(R):
    f = Piecewise(((1 - q/R)**7*(1 + 7*q/R + 19*(q/R)**2 + 21*(q/R)**3),q < R), (0, True))
    return(f,'Wendland 1D C^6')

#--------------------------
# Wendland kernels in 2/3D
#--------------------------
def w2(R):
    f = Piecewise(((1 - q/R)**4*(1 + 4*q/R),q < R), (0, True))
    return(f,'Wendland 2/3D C^2')

def w4(R):
    f = Piecewise(((1 - q/R)**6*(1 + 6*q/R + sympify(35)/3*(q/R)**2),q < R), (0, True))
    return(f,'Wendland 2/3D C^4')

def w6(R):
    f = Piecewise(((1 - q/R)**8*(1 + 8*q/R + 25*(q/R)**2 + 32*(q/R)**3),q < R), (0, True))
    return(f,'Wendland 2/3D C^6')

def sinq(R,n):
    f = Piecewise(((sin(pi*q/R)/q)**n,q < R), (0, True))
    name = "[sin(q)/q]**%i" %n
    return(f,name)

def pcubic(R):
    f = Piecewise((q**3 - 6*q + 6,q < 1),((2-q)**3,q < 2), (0, True))
    return(f,'A peaked cubic')

def bcubic(R):
    f = Piecewise(((sympify(10) - sympify(13)*q**2 + sympify(6)*q**3)/16,q < 1), ((2 - q)**2*(5 - sympify(2)*q)/16, q < 2),(0, True))
    return(f,'Better cubic')

#-----------------------------
# integrated B-spline kernels
#-----------------------------
def intm4(R):
    return(intkernel(m4,R))

def int2m4(R):
    return(intkernel2(m4,R))

def int3m4(R):
    return(intkernel3(m4,R))

def intm5(R):
    return(intkernel(m5,R))

def intm6(R):
    return(intkernel(m6,R))

#------------------
# Ferrers kernels
#------------------
def f2(R):
    f = Piecewise(((1 - (q/R)**2)**2,q < R), (0, True))
    return(f,'Ferrers n=4')

def f3(R):
    f = Piecewise(((1 - (q/R)**2)**3,q < R), (0, True))
    return(f,'Ferrers n=4')

def f4(R):
    f = Piecewise(((1 - (q/R)**2)**4,q < R), (0, True))
    return(f,'Ferrers n=4')

def f5(R):
    f = Piecewise(((1 - (q/R)**2)**5,q < R), (0, True))
    return(f,'Ferrers n=5')

def f6(R):
    f = Piecewise(((1 - (q/R)**2)**6,q < R), (0, True))
    return(f,'Ferrers n=6')

########################################################
#  The actual program
#  Change the lines below to print the kernel you want
########################################################

# set kernel range
#R = sympify(3)
#R = sympify(5)/2
R = sympify(2)
#R = symbols('R')

# define which kernel to use
#f, name = sinq(R,3)
f, name = m4(R)
#f, name = w2(R)

#printvariances(f,R)
#f, name = doublehump5(m4,R)

# print the desired output
#printkernel(f,R)
#printkernel_ndspmhd(f,R,name)
printkernel_phantom(f,R,name)
#printkernel_sphNG(f,R,name)
#printall_tex
