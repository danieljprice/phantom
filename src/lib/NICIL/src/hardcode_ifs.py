#----------------------------------------------------------------------!
#                               N I C I L                              !                                  
#           Non-Ideal mhd Coefficient and Ionisation Library           !
#            A script to generate new source code with most            !
#                     logical if-statements removed                    !
#                                                                      !
#                 Copyright (c) 2015-2019 James Wurster                !
#        See LICENCE file for usage and distribution conditions        !
#----------------------------------------------------------------------!
# The primary function of this script is to make a new copy of
# nicil.F90, where the if-statements of selected input parameters have
# been removed (except in the initialisation subroutines). The interior
# commands will be included or deleted based upon the choices of the
# logicals.  The original file will be saved as nicil_source.F90.
# In production runs, it is not necessary to be continually calling the
# if-statements for values that are always true or false.
# Secondary functions are to copy nicil_source.F90 back onto nicil.F90,
# to list the hardcoded changes in nicil.F90, and to diff the files.
#----------------------------------------------------------------------!
import os
import sys
#
print "Welcome to Nicil's hardcode_ifs.py script"
#--Filenames; make nicil_source.F90 if it does not yet exist
srcname = "nicil_source.F90"
outname = "nicil.F90"
tmpname = "nicil_tmp.F90"
difname = "nicilsrc_nicil.diff"
if (not os.path.isfile(srcname)):
  print "Created  "+srcname+" by copying "+outname
  os.system("cp "+outname+" "+srcname)
  # sanity check
  if (not os.path.isfile(srcname)):
    print "Failed to create "+srcname+".  Aborting"
    sys.exit()
#
#--List the options and ask the user which task they would like performed
print "Please choose from one of the following options:"
print "  1) Remove hardcoded if's from nicil.F90"
print "  2) List the hard coded options and their status in nicil.F90"
print "  3) Replace nicil.F90 with nicil_source.F90"
print "  4) diff nicil.F90 and nicil_source.F90"
print "  5) Developer: remove 'pure' statements"
print "  6) Developer: re-add 'pure' statements"
opt = str(raw_input("Please enter option now: "))
ask = True
while ( ask ):
  if (opt=="1" or opt=="2" or opt=="3" or opt=="4" or opt=="5" or opt=="6" or opt[0:1]=="q" or opt[0:1]=="Q"):
    ask = False
  else:
    opt = str(raw_input("That is not a valid input.  Please enter option now:"))
#
#----------------------------------------------------------------------!
#+
# Make a new copy if nicil.F90, where the parameter if's are hardcoded
#+
#----------------------------------------------------------------------!
if (opt=="1"):
  #
  #--The list of logicals to remove
  clogical = [];                       llogical = []
  clogical.append('use_ohm');          llogical.append('')
  clogical.append('use_hall');         llogical.append('')
  clogical.append('use_ambi');         llogical.append('')
  clogical.append('ion_rays');         llogical.append('')
  clogical.append('ion_thermal');      llogical.append('')
  clogical.append('zeta_of_rho');      llogical.append('')
  clogical.append('use_fdg_in');       llogical.append('')
  clogical.append('rho_is_rhogas');    llogical.append('')
  clogical.append('eta_constant');     llogical.append('')
  clogical.append('mod_beta');         llogical.append('')
  clogical.append('warn_verbose');     llogical.append('')
  clogical.append('reorder_Jacobian'); llogical.append('')
  #
  #--Determine the value of the required logicals
  a=open(srcname,'r')
  g_cnst       = 0
  use_massfrac = 0
  for line in a:
    if ("logical" in line and ("public" in line or "private" in line)):
      for i in range(0,len(clogical)):
        if (clogical[i] in line):
          if ("True" in line or "true" in line or "TRUE" in line):
            llogical[i]=True
          else:
            llogical[i] = False
    if ("g_cnst" in line and "= .true."  in line):
       g_cnst =  1
    elif ("g_cnst" in line and "= .false." in line):
       g_cnst = -1
    elif ("use_massfrac" in line and "= .true."  in line):
       use_massfrac =  1
    elif ("use_massfrac" in line and "= .false." in line):
       use_massfrac = -1
    elif ("END OF INPUT PARAMETERS" in line): break
  a.close
  if (g_cnst==0):
    print "There has been a problem determining the value of < g_cnst >.  Aborting."
    sys.exit()
  if (use_massfrac==0):
    print "There has been a problem determining the value of < use_massfrac >.  Aborting."
    sys.exit()
  print " "
  print "The following logicals will be hard-coded as defined:"
  for i in range(0,len(clogical)):
    print " "+clogical[i]+"="+str(llogical[i])
  print " "
  print "Warning: some variables may now be defined but never used"
  print "Warning: indentation may no longer be logical"
  print "Warning: orphaned comments may exist"
  print " "
  #
  #--Open files and write the modified output file
  a = open(srcname,"r")
  b = open(outname,"w")
  #
  if_level  = 0
  keep_text = True
  allow_mod = False
  for line in a:
    write_line = True
    # Include "hardcoded" comment in the definition lines
    if ("logical" in line and "public" in line):
      for i in range(0,len(clogical)):
        if (clogical[i] in line):
          write_line = False
          comment = line.find("!")
          if (comment==-1):
            spaces = "          "
          else:
            spaces = ""
          b.write(line[:comment]+spaces+"! HARDCODED="+str(llogical[i])+" "+line[comment:])
    #
    # Add statements such that hard coded variable will be defined in the output log
    elif ("end subroutine nicil_print_summary" in line):
      b.write(" write(iprint,'(a)') 'NICIL: This version has been modified by hardcodeifs.py' \n")
      for i in range(0,len(clogical)):
        b.write(" write(iprint,'(a)') 'NICIL: HARDCODED PARAMETER: "+clogical[i]+"="+str(llogical[i])+"'\n")
      b.write(" !\n")
    #
    # If-statement exists.  Remove as required
    elif (allow_mod  and " if" in line and "end" not in line):
      comment = line.find("!")
      cif     = line.find(" if")
      if (comment==-1 or cif < comment):
        ilog     = -1
        for i in range(0,len(clogical)):
          if ((clogical[i]) in line): ilog = i
        if (ilog > -1):
          write_line = False
          # an if-then statement
          if ("then" in line):
            if_level = if_level + 1
            if ((llogical[ilog] and "not" not in line) or (not llogical[ilog] and "not" in line)):
              keep_text = True
            else:
              keep_text = False
          else:
            # an if statement on a single line
            bracket = line.find(")")
            and_    = line.find(".and.")
            if ((llogical[ilog] and "not" not in line) or (not llogical[ilog] and "not" in line)):
              if (and_==-1):
                # a single if statement
                b.write(line[bracket+1:])
              else:
                # a double if statement
                b.write(" if ("+line[and_+5:])
        elif ("then" in line):
          if (if_level > 0): if_level = if_level + 1
    # In if-statement and encountered "else"
    elif (allow_mod and if_level==1 and ("else" in line or "elif" in line)):
      keep_text  = not keep_text
      write_line = False
    # In if-statement and encountered "end if"
    elif (allow_mod and if_level > 0 and ("end if" in line or "endif" in line)):
      if (if_level==1):
        write_line = False
        keep_text  = True
      if_level = if_level - 1
    # Replace loops over nelements with the number since nelements is a variable and not a parameter
    elif ("do" in line and "1,nelements" in line):
      write_line = False
      where = line.find("1,nelements")
      if (use_massfrac==1):
        b.write(line[:where]+"1,3\n")
      elif(use_massfrac==-1):
        b.write(line[:where]+"1,6\n")
      else:
        b.write(line+"\n")
    elif ("do" in line and "3,nelements" in line):
      write_line = False
      where = line.find("3,nelements")
      if (use_massfrac==1):
        b.write(line[:where]+"3,3\n")
      elif (use_massfrac==-1):
        b.write(line[:where]+"3,6\n")
      else:
        b.write(line+"\n")
    # Replace loops over na with the number since na is a variable and not a parameter
    elif ("do" in line and "1,na" in line):
      write_line = False
      where = line.find("1,na")
      if (g_cnst==1):
        b.write(line[:where]+"1,1\n")
      elif (g_cnst==-1):
        b.write(line[:where]+"1,na_max\n")
      else:
        b.write(line+"\n")
    elif ("end subroutine nicil_print_summary" in line):
      allow_mod = True
    #
    # No action required: copy line verbatim
    if (write_line and keep_text): b.write(line)
  a.close()
  b.close()
#----------------------------------------------------------------------!
#+
# List the hardcoded ifs currently in nicil.F90
#+
#----------------------------------------------------------------------!
elif (opt=="2"):
  print "The hardcoded logicals current in nicil.F90 are"
  hardcodes = False
  a=open(outname,"r")
  for line in a:
    if ("NICIL: HARDCODED PARAMETER" in line):
      print line[49:-2]
      hardcodes = True
  a.close()
  if (not hardcodes): print "none"
#----------------------------------------------------------------------!
#+
# Replace nicil.F90 with nicil_source.F90
#+
#----------------------------------------------------------------------!
elif (opt=="3"):
  print "Replaced  "+outname+" with "+srcname
  os.system("cp "+srcname+" "+outname)
#----------------------------------------------------------------------!
#+
# diff nicil.F90 and nicil_source.F90
#+
#----------------------------------------------------------------------!
elif (opt=="4"):
  #--Additional promtps
  sf  = str(raw_input("Print results to screen (s) or to file (f) [default=s]: "))
  sbs = str(raw_input("Diff side-by-side (y or n)? [default=n]: "))
  if (sf[0:1]=="f" or sf[0:1]=="F"):
    tofile  = " > "+difname
    comment = " with results printed to "+difname
  else:
    tofile  = ""
    comment = ""
  if (sbs[0:1]=="y" or sbs[0:1]=="Y"):
    sidebyside = " --side-by-side"
  else:
    sidebyside = ""
  #--The command
  os.system("diff "+srcname+" "+outname+sidebyside+tofile)
  print "diffed  "+srcname+" and "+outname+comment
#----------------------------------------------------------------------!
#+
# Re-write nicil.F90 without pure subroutines/functions
#+
#----------------------------------------------------------------------!
elif (opt=="5"):
  ictr = 0
  a = open(outname,"r")
  b = open(tmpname,"w")
  #
  for line in a:
    if ("pure " in line):
       b.write(line[5:-1]+"  !pure!\n")
       ictr = ictr + 1
    else:
       b.write(line)
  a.close()
  b.close()
  os.system("mv "+tmpname+" "+outname)
  print "'Pure' was removed "+str(ictr)+" times."
#----------------------------------------------------------------------!
#+
# Re-write nicil.F90 replacing pure subroutines/functions
# This will only work to undo opt=5
#+
#----------------------------------------------------------------------!
elif (opt=="6"):
  ictr = 0
  a = open(outname,"r")
  b = open(tmpname,"w")
  #
  for line in a:
    if ("!pure!" in line):
       b.write("pure "+line[:-9]+"\n")
       ictr = ictr + 1
    else:
       b.write(line)
  a.close()
  b.close()
  os.system("mv "+tmpname+" "+outname)
  print "'Pure' was re-added "+str(ictr)+" times."
#----------------------------------------------------------------------!
#+
# none of the above: exit
#+
#----------------------------------------------------------------------!
elif (opt[0:1]=="q" or opt[0:1]=="Q"):
  sys.exit()
#----------------------------------------------------------------------!