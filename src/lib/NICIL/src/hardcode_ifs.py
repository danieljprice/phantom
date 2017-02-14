#----------------------------------------------------------------------!
#                               N I C I L                              !                                  
#           Non-Ideal mhd Coefficient and Ionisation Library           !
#            A script to generate new source code with most            !
#                     logical if-statements removed                    !
#                                                                      !
#                 Copyright (c) 2015-2016 James Wurster                !
#        See LICENCE file for usage and distribution conditions        !
#----------------------------------------------------------------------!
# The primary function of this script is to make a new copy of
# nicil.F90 from nicil_source.F90, where the if-statements of selected
# input parameters have been removed (except in the initialisation
# subroutines).
# The interior commands will be included or deleted based upon the 
# logical in nicil_source.F90.
# In production runs, it is not necessary to be continually calling the
# if-statements for values that are always true or false.
# Secondary functions are to copy nicil_source.F90 back onto nicil.F90,
# to list the hardcoded changes in nicil.F90, and to diff the files.
#----------------------------------------------------------------------!
import os
import sys
#
#--List the options and ask the user which task they would like performed
print "Welcome to Nicil's hardcode_ifs.py.  Please choose from one of the following options:"
print "  1) Remove hardcoded if's from nicil.F90"
print "  2) List the hard coded options and their status in nicil.F90"
print "  3) Replace nicil.F90 with nicil_source.F90"
print "  4) diff nicil.F90 and nicil_source.F90"
opt = str(raw_input("Please enter option now: "))
ask = True
while ( ask ):
  if (opt=="1" or opt=="2" or opt=="3" or opt=="4" or opt[0:1]=="q" or opt[0:1]=="Q"):
    ask = False
  else:
    opt = str(raw_input("That is not a valid input.  Please enter option now:"))
#
#--Filenames
srcname = "nicil_source.F90"
outname = "nicil.F90"
difname = "nicilsrc_nicil.diff"
#
#----------------------------------------------------------------------!
#+
# Make a new copy if nicil.F90, where the parameter if's are hardcoded
#+
#----------------------------------------------------------------------!
if (opt=="1"):
  #
  #--The list of logicals to remove
  clogical = [];                     llogical = []
  clogical.append('use_ohm');        llogical.append('')
  clogical.append('use_hall');       llogical.append('')
  clogical.append('use_ambi');       llogical.append('')
  clogical.append('ion_rays');       llogical.append('')
  clogical.append('ion_thermal');    llogical.append('')
  clogical.append('eta_constant');   llogical.append('')
  clogical.append('mod_beta');       llogical.append('')
  clogical.append('eta_const_calc'); llogical.append('')
  clogical.append('warn_verbose');   llogical.append('')
  #
  #--Determine the value of the required logicals
  a=open(srcname,'r')
  for line in a:
    if ("logical" in line and "public" in line):
      for i in range(0,len(clogical)):
        if (clogical[i] in line):
          if ("True" in line or "true" in line or "TRUE" in line):
            llogical[i]=True
          else:
            llogical[i] = False
    if ("END OF INPUT PARAMETERS" in line): break
  a.close
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
  a = open(srcname,'r')
  b = open(outname,'w')
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
    elif ("Version 1.0: 8 Dec 2015:" in line):
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
# none of the above: exit
#+
#----------------------------------------------------------------------!
elif (opt[0:1]=="q" or opt[0:1]=="Q"):
  sys.exit()
#----------------------------------------------------------------------!