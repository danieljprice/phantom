#
# Analyse a Phantom dump file with python
#
# Written by David Liptai 2018
#

import pyphantom
from ctypes import *

class PhantomAnalysis(pyphantom.Simulation):
    def __init__(self, dumpfile):
      self.libph = CDLL('./libphantom.so')

      # Read dumpfile
      time,hfact,massofgas = self.read_dump(dumpfile)

      # Units
      udist, umass, utime, udens, umagfd = self.get_units()
      self.units = {'udist':udist,
                    'umass':umass,
                    'utime':utime,
                    'udens':udens,
                    'umagfd':umagfd}

      # Note that we divide by code units so that quantities are in code units by default

      # General stuff
      self.time   = time
      self.hfact  = hfact
      self.massofgas  = massofgas/umass
      ##########################
      # Gas quantities
      ########################
      npart = self.get_npart()
      self.npart       = npart
      self.xyzh        = self.get_part_xyzh(npart)/udist
      self.vxyz        = self.get_part_vxyz(npart)/(udist/utime)
      ####################################
      # Gas quantities: Try loading internal energy and magnetic field
      ####################################
      try:
          self.utherm      = self.get_part_u(npart)/(udist**2/utime**2)
      except: 
          print "Unable to load uterm, likely the quantity is not stored."
      try:
          self.temperature = self.get_part_temp(npart)
      except: 
       	  print	"Unable to load temperature, likely the quantity is not stored."
      try:
          self.bxyz        = self.get_part_bxyz(npart)/umagfd
      except: 
       	  print	"Unable to load magnetic field, likely the quantity is not stored."
      ####################

      # Point masses
      nptmass             = self.get_nptmass()
      self.nptmass        = nptmass
      self.ptmass_xyzmh   = self.get_ptmass_xyzmh(nptmass)
      self.ptmass_xyzmh[0,:] = self.ptmass_xyzmh[0,:]/udist
      self.ptmass_xyzmh[1,:] = self.ptmass_xyzmh[1,:]/udist
      self.ptmass_xyzmh[2,:] = self.ptmass_xyzmh[2,:]/udist
      self.ptmass_xyzmh[3,:] = self.ptmass_xyzmh[3,:]/umass
      self.ptmass_xyzmh[4,:] = self.ptmass_xyzmh[4,:]/udist
      self.ptmass_vxyz    = self.get_ptmass_vxyz(nptmass)/(udist/utime)
      self.ptmass_spinxyz = self.get_ptmass_spinxyz(nptmass)/(udist**2*umass/utime)
