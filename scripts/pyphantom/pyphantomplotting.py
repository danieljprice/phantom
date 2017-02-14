#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyphantom
from constants_cgs import R, G, au, year
from pylab import *

from multiprocessing import Process, JoinableQueue, connection, Pipe
import matplotlib.animation as animation

import re
from os import listdir

class RealTimePlotSimulation:
  def __init__(self, SimulationClass, remoteplotter):
    self.remoteplotter = remoteplotter
    self.SimulationClass = SimulationClass
    self.plotter_conn, worker_conn = Pipe()
    self.queue = JoinableQueue(maxsize=1)
    self.axesextent = remoteplotter.axesExtent()
    self.plotter_conn.send(('axes',self.axesextent))
    phantomprocess = Process(target=self.work, args=(self.queue, worker_conn, SimulationClass, ))
    phantomprocess.start()
    anim = animation.FuncAnimation(remoteplotter.fig, self.update_plot, interval=100)
    show()
    # The phantom process is properly stopped
    print 'Killing the job...'
    self.plotter_conn.send(('gokillyourself',))
    phantomprocess.join()
    print 'Done.'

  def update_plot(self, num):
    extent = self.remoteplotter.axesExtent()
    if extent != self.axesextent:
      self.plotter_conn.send(('axes',extent))
      self.axesextent = extent
    if self.queue.empty(): return
    data = self.queue.get()
    self.remoteplotter.plotData(data)
    self.queue.task_done()

  def work(self, queue, worker_conn, SimulationClass):
    s = SimulationClass()
    # Main loop
    for time in s.simulation:
      s.inmainloop()
      # Check for orders
      if worker_conn.poll():
	message = worker_conn.recv()
	if message[0] == 'gokillyourself':
	  print 'Suicide order received. Exiting main loop.'
	  break
	if message[0] == 'axes':
	  extent = message[1]
      if queue.empty():
	# Retreive data
	queue.put(self.remoteplotter.pickableData(s.simulation, extent))
    queue.cancel_join_thread()
    
class PostProcVideo:
  def __init__(self, infile, remoteplotter):
    self.remoteplotter = remoteplotter
    prefix = infile.rstrip('.in')
    dumpfiles = []
    pattern = re.compile(r'%s_(\d+)$' % prefix)
    for filename in listdir('.'):
      match = re.search(pattern, filename)
      if match:
	number = int(match.group(1))
	dumpfiles.append(number)
    dumpfiles = sorted(dumpfiles)
    images = []
    for number in dumpfiles:
      if number != 0:
	img = self.plot(prefix, number)
	if img:
	  images.append(img)
    anim = animation.ArtistAnimation(remoteplotter.fig, images, interval=50, blit=True)
    anim.save('%s.mp4' % prefix, extra_args=['-vcodec', 'libx264'])
      
  def plot(self, prefix, number):
    dumpfile = '%s_%05d' % (prefix, number)
    a = pyphantom.Simulation(prefix+'.in', 0., 0., analysis=True)
    print 'Reading ', dumpfile, '...'
    try:
      a.update_from_dump(dumpfile)
    except pyphantom.UnableToReadDumpFile:
      print 'Skipping small dump.'
      return None
    extent = self.remoteplotter.axesExtent()
    data = self.remoteplotter.pickableData(a, extent)
    return self.remoteplotter.plotData(data)

class RemotePlotter():
  def pickableData(self, simulation, extent):
    return None
  def plotData(self, data):
    pass
  def axesExtent(self):
    return None
  
class RemotePlotter_rho_xysec_ptmass(RemotePlotter):
  npx = npy = 512
  def __init__(self, fig, ptmasses, log=False):
    self.fig = fig
    self.axes = fig.add_subplot(1, 1, 1, axisbg='black')
    self.imgzoomed = self.axes.imshow(zeros((self.npx,self.npy)), extent=(-1., 1., -1., 1.), interpolation='bicubic', cmap='gist_heat')
    self.img = self.axes.imshow(zeros((self.npx,self.npy)), extent=(-1., 1., -1., 1.), interpolation='bicubic', cmap='gist_heat')
    self.axes.set_xlabel('x (au)')
    self.axes.set_ylabel('y (au)')
    self.fig.colorbar(self.imgzoomed)
    self.ptmasses = ptmasses
    self.c = []
    for ptmass in ptmasses:
      c = matplotlib.patches.Circle((0.,0.), 1., fill=False, ec='red', linewidth=2)
      self.axes.add_artist(c)
      self.c.append(c)
    self.log = log

  def pickableData(self, simulation, extent):
    npart = simulation.get_npart()
    part_xyzh = simulation.get_part_xyzh(npart)
    limits = simulation.get_boundaries()
    massofgas = simulation.get_massofgas()
    hfact = simulation.get_hfact()
    nptmass = simulation.get_nptmass()
    ptmass_xyzmh = simulation.get_ptmass_xyzmh(nptmass)
    xmin, xmax, ymin, ymax, _, _ = limits
    dx = (xmax-xmin)/float(self.npx)
    dy = (ymax-ymin)/float(self.npy)
    if self.log:
      pixmap = rot90(simulation.plot_log_rho_xysec(0., xmin, ymin, self.npx, self.npy, dx, dy, massofgas, hfact, part_xyzh))
    else:
      pixmap = rot90(simulation.plot_rho_xysec(0., xmin, ymin, self.npx, self.npy, dx, dy, massofgas, hfact, part_xyzh))
    if extent != None:
      xmin, xmax, ymin, ymax = extent
    dx = (xmax-xmin)/float(self.npx)
    dy = (ymax-ymin)/float(self.npy)
    if self.log:
      pixmap_zoomed = rot90(simulation.plot_log_rho_xysec(0., xmin, ymin, self.npx, self.npy, dx, dy, massofgas, hfact, part_xyzh))
    else:
      pixmap_zoomed = rot90(simulation.plot_rho_xysec(0., xmin, ymin, self.npx, self.npy, dx, dy, massofgas, hfact, part_xyzh))
    time = simulation.time
    extent = xmin, xmax, ymin, ymax
    return (limits, extent, pixmap, pixmap_zoomed, ptmass_xyzmh, time)
  def axesExtent(self):
    xmin, xmax = self.axes.get_xlim()
    ymin, ymax = self.axes.get_ylim()
    return (xmin*au, xmax*au, ymin*au, ymax*au)    
  def plotData(self, data):
    limits, extent, pixmap, pixmap_zoomed, ptmass_xyzmh, time = data
    xmin, xmax, ymin, ymax = extent
    self.imgzoomed.set_array(pixmap_zoomed)
    self.imgzoomed.set_extent((xmin/au, xmax/au, ymin/au, ymax/au))
    xmin, xmax, ymin, ymax, _, _ = limits
    self.img.set_array(pixmap)
    self.img.set_extent((xmin/au, xmax/au, ymin/au, ymax/au))
    cmin = amin(pixmap_zoomed)
    cmax = amax(pixmap_zoomed)
    self.img.set_clim(cmin, cmax)
    self.imgzoomed.set_clim(cmin, cmax)
    for ic, ptmass in enumerate(self.ptmasses):
      x, y, _, _, h = ptmass_xyzmh[:,ptmass]
      self.c[ic].set_radius(h/au)
      self.c[ic].center = x/au, y/au
    if self.log:
      self.axes.set_title(u'log density (g/cm³)   t=%.1f yrs' % (time/year))
    else:
      self.axes.set_title(u'density (g/cm³)   t=%.1f yrs' % (time/year))
    return set([self.imgzoomed, self.img] + self.c)
    
class RemotePlotter_sphericalWindProfile(RemotePlotter):
  def __init__(self, fig):
    self.fig = fig
    self.axes = fig.add_subplot(1, 1, 1)
    self.bg = self.axes.imshow(zeros((1,1)), extent=(-1., 1., -1., 1.), cmap='binary')
    self.vplot, = self.axes.plot([], [], '. ', markersize=.1, rasterized=True, label='SPH')
    self.thplot, = self.axes.plot([], [], 'r', label='Theory')
    self.axes.legend(loc=1)
    self.axes.set_xlabel('r (au)')
    self.axes.set_ylabel('v (km/s)')
  def pickableData(self, simulation, extent):
    npart = simulation.get_npart()
    part_xyzh = simulation.get_part_xyzh(npart)
    part_vxyz = simulation.get_part_vxyz(npart)
    time = simulation.time
    r = sqrt(part_xyzh[0,:]**2+part_xyzh[1,:]**2+part_xyzh[2,:]**2)
    vr = sqrt(part_vxyz[0,:]**2+part_vxyz[1,:]**2+part_vxyz[2,:]**2)
    rmax = amax(r)
    r_th, v_th, _ = simulation.sphericalwind.theoretical_solution(rmax)
    return (npart, r, vr, r_th, v_th, time)
  def axesExtent(self):
    xmin, xmax = self.axes.get_xlim()
    ymin, ymax = self.axes.get_ylim()
    return (xmin*au, xmax*au, ymin*au, ymax*au)
  def plotData(self, data):
    npart, r, vr, r_th, v_th, time = data
    self.vplot.set_data(r/au, vr/1.e5)
    self.thplot.set_data(r_th/au, v_th/1.e5)
    self.axes.set_title(u'Spherical wind velocity profile    t=%.1f yrs' % (time/year))
    self.bg.set_extent((amin(r)/au, amax(r)/au, amin(vr)/1.e5, amax(vr)/1.e5))
    #self.axes.set_xlim(amin(r_th)/au, amax(r_th)/au)
    #self.axes.set_ylim(amin(v_th)/1.e5, amax(v_th)/1.e5)
    return self.bg, self.vplot, self.thplot
    
class RemotePlotter_particles_xy(RemotePlotter):
  def __init__(self, fig):
    self.fig = fig
    self.axes = fig.add_subplot(1, 1, 1)
    self.bg = self.axes.imshow(zeros((1,1)), extent=(-1., 1., -1., 1.), cmap='binary')
    self.xyplot, = self.axes.plot([], [], '. ', markersize=.1, rasterized=True)
    self.axes.set_xlabel('r (au)')
    self.axes.set_ylabel('v (km/s)')
  def pickableData(self, simulation, extent):
    npart = simulation.get_npart()
    part_xyzh = simulation.get_part_xyzh(npart)
    time = simulation.time
    x = part_xyzh[0,:]
    y = part_xyzh[1,:]
    return (npart, x, y, time)
  def axesExtent(self):
    xmin, xmax = self.axes.get_xlim()
    ymin, ymax = self.axes.get_ylim()
    return (xmin*au, xmax*au, ymin*au, ymax*au)
  def plotData(self, data):
    npart, x, y, time = data
    self.xyplot.set_data(x/au, y/au)
    self.axes.set_title(u'SPH particles    t=%.1f yrs' % (time/year))
    self.bg.set_extent((amin(x)/au, amax(x)/au, amin(y)/au, amax(y)/au))
    #self.axes.set_xlim(amin(r_th)/au, amax(r_th)/au)
    #self.axes.set_ylim(amin(v_th)/1.e5, amax(v_th)/1.e5)
    return self.bg, self.xyplot
    
def plot_cross_rho(fig, simulation, limits, part_xyzh, ptmass_xyzmh, massofgas, hfact, time):
  xmin, xmax, ymin, ymax = limits
  npx = npy = 512
  dx = (xmax-xmin)/float(npx)
  dy = (ymax-ymin)/float(npy)
  pixmap = rot90(simulation.plot_rho_xysec(0., xmin, ymin, npx, npy, dx, dy, massofgas, hfact, part_xyzh))
  clf()
  axes = fig.add_subplot(1, 1, 1, axisbg='black')
  img = axes.imshow(pixmap, extent=(xmin/au, xmax/au, ymax/au, ymin/au), interpolation='bicubic', cmap='gist_heat')
  fig.colorbar(img)
  nptmass = shape(ptmass_xyzmh)[1]
  axes.set_title(u'log density (g/cm³)   t=%.1f yrs' % (time/year))
  axes.set_xlabel('x (au)')
  axes.set_ylabel('y (au)')
  for i in range(nptmass):
    x, y, z, m, h = ptmass_xyzmh[:,i]
    c = matplotlib.patches.Circle((x/au,y/au), h/au, fill=False, ec='red', linewidth=2)
    axes.add_artist(c)

def plot_part_xy(fig, part_xyzh, ptmass_xyzmh, time):
  clf()
  x = part_xyzh[0,:]/au
  y = part_xyzh[1,:]/au
  axes = fig.add_subplot(1, 1, 1)
  axes.plot(x, y, '.', markersize=.5)
  axes.set_title(u'SPH particles   t=%.1f yrs' % (time/year))
  axes.set_xlabel('x (au)')
  axes.set_ylabel('y (au)')
  nptmass = shape(ptmass_xyzmh)[1]
  for i in range(nptmass):
    x, y, z, m, h = ptmass_xyzmh[:,i]
    c = matplotlib.patches.Circle((x/au,y/au), h/au, fill=False, ec='red', linewidth=2)
    axes.add_artist(c)

def plot_velocity_profile(fig, part_r, part_vxyz, time):
  vx, vy, vz = [part_vxyz[i,:] for i in range(3)]
  vr = sqrt(vx**2+vy**2+vz**2)
  clf()
  plot(part_r/au, vr/(1e5),'b.', markersize=.5, rasterized=True)
  xlabel('r (au)')
  ylabel('v (km/s)')
  
def plot_temperature_profile(fig, part_r, part_T, time):
  clf()
  plot(part_r/au, part_T, 'b.', markersize=.5, rasterized=True)
  xlabel('r (au)')
  ylabel('T (K)')
  
def plot_density_profile(fig, part_r, part_rho, time):
  clf()
  plot(part_r/au, part_rho, 'b.', markersize=.5, rasterized=True)
  xlabel('r (au)')
  ylabel(u'density (g/cm³)')