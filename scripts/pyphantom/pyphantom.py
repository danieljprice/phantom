#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ctypes import *
from numpy import zeros, shape, loadtxt, concatenate
from os import listdir, path
import re

# Exceptions
class UnableToReadDumpFile(Exception):
  pass
class IncorrectNumberOfParticles(Exception):
  pass
class SpecificEnergyNotStored(Exception):
  pass
class TemperatureNotStored(Exception):
  pass
class MagneticFieldNotStored(Exception):
  pass
class IncorrectNumberOfPointMasses(Exception):
  pass
class TryingToChangeMassInNonEmptySimulation(Exception):
  pass

class Simulation():
  def __init__(self, infile, tmax, dtmax, analysis=False):
    self.libph = CDLL('./libphantom.so')
    self.infile = infile
    if not analysis:
      self.logfile, self.evfile, self.dumpfile = self.init(infile, tmax, dtmax)
      self.update_time()
    
  def update_time(self):
    try:
      self.dtlast = self.dt
    except AttributeError:
      self.dtlast = 0.
    self.time, self.tmax, self.nsteps, self.nmax, self.dt = self.get_time()
    
  def __iter__(self):
    return self
  
  def next(self):
    self.update_time()
    if self.time >= self.tmax or self.nsteps >= self.nmax:
      raise StopIteration
    else:
      self.init_step()
      self.step()
      self.logfile, self.evfile, self.dumpfile = self.finalize_step(self.infile, self.logfile, self.evfile, self.dumpfile)
    return self.time
    
  def __del__(self):
    self.finalize()
    
  def read_ev_files(self):
    prefix = self.infile.rstrip('.in')
    print "Searching files with prefix '%s'..." % prefix
    evfiles = sorted([f for f in listdir('.') if f.startswith(prefix) and f.endswith('.ev')])
    udist, umass, utime, udens, umagfd = self.get_units()
    wololo = {'time':utime, 'ekin':udist**2*umass/utime**2, 'etherm':udist**2*umass/utime**2, 'emag':udist**2*umass/utime**2, 'epot':udist**2*umass/utime**2, 'etot':udist**2*umass/utime**2, 'totmom':umass*udist/utime, 'angtot':umass*udist**2/utime, 'rhomax':udens, 'rhomean':udens, 'dt':utime, 'rmsmach':1., 'alphamax':1., 'alphau':1., 'alphaBmax':1., 'totmomall':umass*udist/utime, 'angall':umass*udist**2/utime, 'accretedmas':umass, 'totentrop':udist**2*umass/utime**2, 'Macc sink 1':umass, 'Macc sink 2':umass, 'divB max':umagfd, 'divB mean':umagfd, 'hdivB/B max':1./udist, 'hdivB/B av':1./udist}
    def concatdicts(A, B):
      for k in B.keys():
	if k in A.keys():
	  try:
	    A[k] = concatenate([A[k], B[k]])
	  except ValueError: pass
	else:
	  A[k] = B[k]
    output = dict()
    datefirst = path.getmtime(evfiles[0])
    for evfile in evfiles:
      if path.getmtime(evfile) >= datefirst:
	data = self.read_ev_file(evfile)
	concatdicts(output, data)
    for column in output.keys():
      if column in wololo.keys():
	output[column] *= wololo[column] # Unit conversion to cgs
    return output

  def read_ev_file(self, evfile):
    print 'Reading ev file: ', evfile
    with open(evfile, 'r') as f:
      header = f.readline()
      column_names = re.findall(r'\[\d+\s+([^\]]+)\]', header)
    data = loadtxt(evfile,unpack=True)
    output = dict()
    for column_number, column_name in enumerate(column_names):
      output[column_name] = data[column_number]
    return output
  
  def update_from_dump(self, dumpfile):
    self.time, _, _ = self.read_dump(dumpfile)
    
  def inject_or_update_particles(self, ifirst, n, position, velocity, h, u, boundary=False):
    # Inputs
    ifirst_c = c_int(ifirst+1)
    n_c = c_int(n)
    boundary_c = c_int(boundary)
    # Calling subroutine
    self.libph.inject_or_update_particles_(byref(ifirst_c), byref(n_c), c_void_p(position.ctypes.data), c_void_p(velocity.ctypes.data), c_void_p(h.ctypes.data), c_void_p(u.ctypes.data), byref(boundary_c))
    
  def inject_or_update_sphere(self, ifirst, resolution, center, radius, center_velocity, expansion_velocity, h, u, angles, boundary=False):
    # Inputs
    ifirst_c = c_int(ifirst+1)
    resolution_c = c_int(resolution)
    radius_c = c_double(radius)
    expansion_velocity_c = c_double(expansion_velocity)
    h_c = c_double(h)
    u_c = c_double(u)
    boundary_c = c_int(boundary)
    # Calling subroutine
    self.libph.inject_or_update_sphere_(byref(ifirst_c), byref(resolution_c), c_void_p(center.ctypes.data), byref(radius_c), c_void_p(center_velocity.ctypes.data), byref(expansion_velocity_c), byref(h_c), byref(u_c), c_void_p(angles.ctypes.data), byref(boundary_c))
    
  def inject_or_update_sink_particle(self, sink_number, position, velocity, radius, mass):
    # Inputs
    sink_number_c = c_int(sink_number+1)
    radius_c = c_double(radius)
    mass_c = c_double(mass)
    # Calling subroutine
    self.libph.inject_or_update_sink_particle_(byref(sink_number_c), c_void_p(position.ctypes.data), c_void_p(velocity.ctypes.data), byref(radius_c), byref(mass_c))

  def init_step(self):
    # Calling subroutine
    self.libph.init_step_wrapper_()
    
  def step(self):
    # Calling subroutine
    self.libph.step_wrapper_()
    
  def finalize_step(self, infile, logfile, evfile, dumpfile):
    # Inputs
    len_infile_c = c_int(len(infile))
    infile_c = create_string_buffer(infile)
    deflen = 120
    len_logfile_c = c_int(deflen)
    logfile_c = create_string_buffer(logfile.ljust(deflen))
    len_evfile_c = c_int(deflen)
    evfile_c = create_string_buffer(evfile.ljust(deflen))
    len_dumpfile_c = c_int(deflen)
    dumpfile_c = create_string_buffer(dumpfile.ljust(deflen))
    # Calling subroutine
    self.libph.finalize_step_wrapper_(byref(len_infile_c), infile_c, byref(len_logfile_c), logfile_c, byref(len_evfile_c), evfile_c, byref(len_dumpfile_c), dumpfile_c)
    # Returning outputs
    return logfile_c.value.strip(), evfile_c.value.strip(), dumpfile_c.value.strip()

  def inject_particles(self):
    # Calling subroutine
    self.libph.inject_particles_wrapper_()

  def get_time(self):
    # Outputs
    time_c = c_double()
    tmax_c = c_double()
    nsteps_c = c_int()
    nmax_c = c_int()
    dt_c = c_double()
    # Calling subroutine
    self.libph.get_time_(byref(time_c), byref(tmax_c), byref(nsteps_c), byref(nmax_c), byref(dt_c))
    # Returning outputs
    return time_c.value, tmax_c.value, nsteps_c.value, nmax_c.value, dt_c.value

  def set_dt(self, dt):
    # Input
    dt_c = c_double(dt)
    # Calling subroutine
    self.libph.set_dt_(byref(dt_c))

  def init(self, infile, tmax, dtmax):
    # Inputs
    len_infile = c_int(len(infile))
    infile_c = create_string_buffer(infile)
    tmax_c = c_double(tmax)
    dtmax_c = c_double(dtmax)
    # Outputs
    deflen = 120
    logfile_c = create_string_buffer(' '*deflen)
    evfile_c = create_string_buffer(' '*deflen)
    dumpfile_c = create_string_buffer(' '*deflen)
    # Calling subroutine
    self.libph.init_(byref(len_infile), infile_c, logfile_c, evfile_c, dumpfile_c, byref(tmax_c), byref(dtmax_c))
    # Returning outputs
    return logfile_c.value.strip(), evfile_c.value.strip(), dumpfile_c.value.strip()

  def override_tmax_dtmax(self, tmax, dtmax):
    # Inputs
    tmax_c = c_double(tmax)
    dtmax_c = c_double(dtmax)
    # Calling subroutine
    self.libph.override_tmax_dtmax_(byref(tmax_c), byref(dtmax_c))
    
  def finalize(self):
    # Calling subroutine
    self.libph.finalize_()
    
  def get_boundaries(self):
    # Outputs
    xmin_c = c_double()
    xmax_c = c_double()
    ymin_c = c_double()
    ymax_c = c_double()
    zmin_c = c_double()
    zmax_c = c_double()
    # Calling subroutine
    self.libph.get_boundaries_(byref(xmin_c), byref(xmax_c), byref(ymin_c), byref(ymax_c), byref(zmin_c), byref(zmax_c))
    # Returning outputs
    return xmin_c.value, xmax_c.value, ymin_c.value, ymax_c.value, zmin_c.value, zmax_c.value

  def read_dump(self, dumpfile, headeronly=False):
    # Inputs
    len_dumpfile = c_int(len(dumpfile))
    dumpfile_c = create_string_buffer(dumpfile)
    headeronly_c = c_int(headeronly)
    # Outputs
    time_c = c_double()
    hfact_c = c_double()
    massofgas_c = c_double()
    ierr_c = c_int()
    # Calling subroutine
    self.libph.read_dump_wrapper_(byref(len_dumpfile), dumpfile_c, byref(headeronly_c), byref(time_c), byref(hfact_c), byref(massofgas_c), byref(ierr_c))
    ierr = ierr_c.value
    # Handling errors
    if not ierr == 0: raise UnableToReadDumpFile
    # Returning outputs
    return time_c.value, hfact_c.value, massofgas_c.value

  def get_massofgas(self):
    # Outputs
    massofgas_c = c_double()
    # Calling subroutine
    self.libph.get_massofgas_(byref(massofgas_c))
    # Returning outputs
    return massofgas_c.value

  def get_hfact(self):
    # Outputs
    hfact_c = c_double()
    # Calling subroutine
    self.libph.get_hfact_(byref(hfact_c))
    # Returning outputs
    return hfact_c.value
  
  def get_npart(self, nodisabled=False):
    # Inputs
    nodisabled_c = c_int(nodisabled)
    # Outputs
    npart_c = c_int()
    # Calling subroutine
    self.libph.get_npart_(byref(npart_c), byref(nodisabled_c))
    # Returning output
    return npart_c.value

  def get_part_xyzh(self, npart, nodisabled=False):
    # Input
    npart_c = c_int(npart)
    nodisabled_c = c_int(nodisabled)
    # Output
    part_xyzh = zeros((4, npart), dtype=float, order='F')
    ierr_c = c_int()
    # Calling subroutine
    self.libph.get_part_xyzh_(byref(npart_c), c_void_p(part_xyzh.ctypes.data), byref(nodisabled_c), byref(ierr_c))
    ierr = ierr_c.value
    # Handling errors
    if ierr == 1: raise IncorrectNumberOfParticles
    # Returning output
    return part_xyzh

  def get_part_vxyz(self, npart, nodisabled=False):
    # Input
    npart_c = c_int(npart)
    nodisabled_c = c_int(nodisabled)
    # Output
    part_vxyz = zeros((3, npart), dtype=float, order='F')
    ierr_c = c_int()
    # Calling subroutine
    self.libph.get_part_vxyz_(byref(npart_c), c_void_p(part_vxyz.ctypes.data), byref(nodisabled_c), byref(ierr_c))
    ierr = ierr_c.value
    # Handling errors
    if ierr == 1: raise IncorrectNumberOfParticles
    # Returning output
    return part_vxyz
   
  def get_part_bxyz(self, npart, nodisabled=False):
    npart_c      = c_int(npart)
    nodisabled_c = c_int(npart) 
    part_bxyz    = zeros(npart, dtype=float)
    ierr_c       = c_int()
    self.libph.get_part_bxyz_(byref(npart_c), c_void_p(part_bxyz.ctypes.data), byref(nodisabled_c), byref(ierr_c))
    ierr = ierr_c.value
    if ierr == 1: raise IncorrectNumberOfParticles
    if ierr == 2: raise MagneticFieldNotStored
    return part_bxyz
    
  def get_part_u(self, npart, nodisabled=False):
    # Input
    npart_c = c_int(npart)
    nodisabled_c = c_int(npart)
    # Output
    part_u = zeros(npart, dtype=float)
    ierr_c = c_int()
    # Calling subroutine
    self.libph.get_part_u_(byref(npart_c), c_void_p(part_u.ctypes.data), byref(nodisabled_c), byref(ierr_c))
    ierr = ierr_c.value
    # Handling errors
    if ierr == 1: raise IncorrectNumberOfParticles
    if ierr == 2: raise SpecificEnergyNotStored
    # Returning output
    return part_u

  def get_part_temp(self, npart, nodisabled=False):
    # Input
    npart_c = c_int(npart)
    nodisabled_c = c_int(npart)
    # Output
    part_temp = zeros(npart, dtype=float)
    ierr_c = c_int()
    # Calling subroutine
    self.libph.get_part_temp_(byref(npart_c), c_void_p(part_temp.ctypes.data), byref(nodisabled_c), byref(ierr_c))
    ierr = ierr_c.value
    # Handling errors
    if ierr == 1: raise IncorrectNumberOfParticles
    if ierr == 2: raise TemperatureNotStored
    # Returning output
    return part_temp

  def get_nptmass(self):
    # Output
    nptmass_c = c_int()
    # Calling subroutine
    self.libph.get_nptmass_(byref(nptmass_c))
    # Returning output
    return nptmass_c.value

  def get_ptmass_xyzmh(self,nptmass):
    # Input
    nptmass_c = c_int(nptmass)
    # Output
    ptmass_xyzmh = zeros((5,nptmass), dtype=float, order='F')
    ierr_c = c_int()
    # Calling subroutine
    self.libph.get_ptmass_xyzmh_(byref(nptmass_c), c_void_p(ptmass_xyzmh.ctypes.data), byref(ierr_c))
    ierr = ierr_c.value
    # Handling error
    if ierr == 1: raise IncorrectNumberOfPointMasses
    # Returning output
    return ptmass_xyzmh

  def get_ptmass_vxyz(self,nptmass):
    # Input
    nptmass_c = c_int(nptmass)
    # Output
    ptmass_vxyz = zeros((3,nptmass), dtype=float, order='F')
    ierr_c = c_int()
    # Calling subroutine
    self.libph.get_ptmass_vxyz_(byref(nptmass_c), c_void_p(ptmass_vxyz.ctypes.data), byref(ierr_c))
    ierr = ierr_c.value
    # Handling error
    if ierr == 1: raise IncorrectNumberOfPointMasses
    # Returning output
    return ptmass_vxyz
  
  def get_ptmass_spinxyz(self,nptmass):
    # Input
    nptmass_c = c_int(nptmass)
    # Output
    ptmass_spinxyz = zeros((3,nptmass), dtype=float, order='F')
    ierr_c = c_int()
    # Calling subroutine
    self.libph.get_ptmass_spinxyz_(byref(nptmass_c), c_void_p(ptmass_spinxyz.ctypes.data), byref(ierr_c))
    ierr = ierr_c.value
    # Handling error
    if ierr == 1: raise IncorrectNumberOfPointMasses
    # Returning output
    return ptmass_spinxyz

  def set_part_mass(self,newmass):
    # Input
    newmass_c = c_double(newmass)
    # Output
    ierr_c = c_int()
    # Calling subroutine
    self.libph.set_part_mass_(byref(newmass_c), byref(ierr_c))
    ierr = ierr_c.value
    # Handling error
    if ierr == 1: raise TryingToChangeMassInNonEmptySimulation

  def plot_rho_xysec(self, z, xmin, ymin, npx, npy, dx, dy, massofgas, hfact, part_xyzh):
    # Inputs
    z_c = c_double(z)
    xmin_c = c_double(xmin)
    ymin_c = c_double(ymin)
    npx_c = c_int(npx)
    npy_c = c_int(npy)
    dx_c = c_double(dx)
    dy_c = c_double(dy)
    npart_c = c_int(shape(part_xyzh)[1])
    massofgas_c = c_double(massofgas)
    hfact_c = c_double(hfact)
    # Outputs
    pixmap = zeros((npx,npy), dtype=float, order='F')
    # Calling subroutine
    self.libph.plot_rho_xysec_(byref(z_c), byref(xmin_c), byref(ymin_c), byref(npx_c), byref(npy_c), byref(dx_c), byref(dy_c), byref(npart_c), byref(massofgas_c), byref(hfact_c), c_void_p(part_xyzh.ctypes.data), c_void_p(pixmap.ctypes.data))
    # Returning output
    return pixmap

  def plot_log_rho_xysec(self, z, xmin, ymin, npx, npy, dx, dy, massofgas, hfact, part_xyzh):
    # Inputs
    z_c = c_double(z)
    xmin_c = c_double(xmin)
    ymin_c = c_double(ymin)
    npx_c = c_int(npx)
    npy_c = c_int(npy)
    dx_c = c_double(dx)
    dy_c = c_double(dy)
    npart_c = c_int(shape(part_xyzh)[1])
    massofgas_c = c_double(massofgas)
    hfact_c = c_double(hfact)
    # Outputs
    pixmap = zeros((npx,npy), dtype=float, order='F')
    # Calling subroutine
    self.libph.plot_log_rho_xysec_(byref(z_c), byref(xmin_c), byref(ymin_c), byref(npx_c), byref(npy_c), byref(dx_c), byref(dy_c), byref(npart_c), byref(massofgas_c), byref(hfact_c), c_void_p(part_xyzh.ctypes.data), c_void_p(pixmap.ctypes.data))
    # Returning output
    return pixmap
  
  def get_units(self):
    # Outputs
    udist = c_double()
    umass = c_double()
    utime = c_double()
    udens = c_double()
    umagfd = c_double()
    # Calling subroutine
    self.libph.get_units_(byref(udist), byref(umass), byref(utime), byref(udens), byref(umagfd))
    # Returning output
    return udist.value, umass.value, utime.value, udens.value, umagfd.value
  
  def delete_particles_outside_box(self, xmin, xmax, ymin, ymax, zmin, zmax):
    # Inputs
    xmin_c = c_double(xmin)
    xmax_c = c_double(xmax)
    ymin_c = c_double(ymin)
    ymax_c = c_double(ymax)
    zmin_c = c_double(zmin)
    zmax_c = c_double(zmax)
    # Calling subroutine
    self.libph.delete_particles_outside_box_wrapper_(byref(xmin_c), byref(xmax_c), byref(ymin_c), byref(ymax_c), byref(zmin_c), byref(zmax_c))
