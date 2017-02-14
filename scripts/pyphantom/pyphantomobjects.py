#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyphantom
from constants_cgs import R, G, au, year
from pylab import norm, cross, asfortranarray, empty, ones, array, newaxis, zeros, sqrt, ceil
from scipy import integrate

class SinkDoesNotExist(Exception):
  pass

class Star:
  def __init__(self, simulation, location, velocity, mass, radius, n):
    nptmass = simulation.get_nptmass()
    self.simulation = simulation
    self.attached_sink = n
    self.location = location
    self.velocity = velocity
    self.mass = mass
    self.radius = radius
    if n < nptmass:
      self.update_properties
    else:
      simulation.inject_or_update_sink_particle(n, location, velocity, mass, radius)
  def update_properties(self):
    nptmass = self.simulation.get_nptmass()
    if self.attached_sink < 0 or self.attached_sink >= nptmass:
      raise SinkDoesNotExist
    ptmass_xyzmh = self.simulation.get_ptmass_xyzmh(nptmass)
    ptmass_vxyz = self.simulation.get_ptmass_vxyz(nptmass)
    xyzmh = ptmass_xyzmh[:,self.attached_sink]
    self.location = xyzmh[:3]
    self.mass = xyzmh[3]
    self.radius = xyzmh[4]
    self.velocity = ptmass_vxyz[:,self.attached_sink]
    
class SquareWind:
  def __init__(self, simulation, center, normal, width, density, velocity, temperature):
    self.simulation, self.center, self.normal, self.width, self.density, self.velocity, self.temperature = simulation, center, normal, width, density, velocity, temperature
    self.resolution = 8
    self.dist_factor = 1.
    self.shift = 0.
    self.handled_walls = 2
    self.gamma = 5./3.
    self.mu = 1.26
    self.update_common_quantities()
  def update_common_quantities(self):
    self.particles_per_wall = self.resolution**2
    self.mass_of_particles = self.simulation.get_massofgas()
    self.mass_of_walls = self.particles_per_wall*self.mass_of_particles
    self.mass_loss_rate = self.density*self.velocity*self.width**2
    self.number_of_particles = self.resolution**2
    self.dx = self.width/self.resolution
    self.time_between_walls = self.mass_of_walls/self.mass_loss_rate
    self.normal /= norm(self.normal)
    self.center = array(self.center)
    self.specific_energy_to_temperature_ratio = R/(self.mu*(self.gamma-1.))
    self.hfact = self.simulation.get_hfact()
    for v in [[1.,0.,0.], [0.,1.,0.], [0.,0.,1.]]:
      u = cross(v,self.normal)
      if norm(u)>0.:
	break
    self.u = u/norm(u)
    v = cross(u,self.normal)
    self.v = v/norm(v)
  def ideal_particle_mass(self):
    idealmass = self.dist_factor*self.density*(self.width/self.resolution)**3
    return idealmass
  def ideal_resolution(self):
    self.update_common_quantities()
    idealresolution = self.width*(self.dist_factor*self.density/self.mass_of_particles)**(1./3.)
    return int(round(idealresolution))
  def inject(self, simulation, time, dtlast):
    #outer_wall = int((time-dtlast)/self.time_between_walls)+1
    #inner_wall = int(time/self.time_between_walls)
    outer_wall = int(ceil((time-dtlast)/self.time_between_walls))
    inner_wall = int(ceil(time/self.time_between_walls)-1)
    inner_handled_wall = inner_wall+self.handled_walls
    s = ''
    for i in range(inner_handled_wall, inner_wall, -1):
      self.inject_wall(simulation, i, time, self.particles_per_wall*(inner_handled_wall-i), boundary=True)
    for i in range(inner_wall, outer_wall-1, -1):
      npart = simulation.get_npart()
      self.inject_wall(simulation, i, time, npart, boundary=False)
  def inject_wall(self, simulation, nwall, time, first_part, boundary):
    local_time = time-nwall*self.time_between_walls
    N = self.resolution
    vertices = empty((N**2,3))
    v = empty((N**2,3))
    for i in range(N):
      x = (float(i)-float(N-1)/2.)*self.dx
      for j in range(N):
	y = (float(j)-float(N-1)/2.)*self.dx
	vertices[i*N+j,:] = self.center + x*self.u + y*self.v + local_time*self.velocity*self.normal
	v[i*N+j,:] = self.normal*self.velocity
    h = self.dx*self.hfact
    specific_energy = self.temperature * self.specific_energy_to_temperature_ratio
    simulation.inject_or_update_particles(first_part, N**2, asfortranarray(vertices.T), v, ones(N**2)*h, ones(N**2)*specific_energy, boundary)
    
class SphericalWind:
  def __init__(self, simulation, star, velocity, temperature, mass_rate):
    self.simulation = simulation
    self.star = star
    self.velocity = velocity
    self.temperature = temperature
    self.mass_rate = mass_rate
    # Default values
    self.resolution = 4
    self.dist_factor = 1.
    self.shift = 0.
    self.handled_spheres = 3
    self.injection_radius = star.radius
    self.gamma = 5./3.
    self.mu = 1.26
    self.update_common_quantities()
  def update_common_quantities(self):
    resolution = self.resolution
    self.particles_per_sphere = 40*resolution*(resolution-1)+12
    self.mass_of_particles = self.simulation.get_massofgas()
    self.mass_of_spheres = self.mass_of_particles * self.particles_per_sphere
    self.time_between_spheres = self.mass_of_spheres/self.mass_rate
    self.specific_energy_to_temperature_ratio = R/(self.mu*(self.gamma-1.))
    self.hfact = self.simulation.get_hfact()
    phi = (sqrt(5.)+1.)/2. # Golden ratio
    self.neighbour_distance = 2./(float(2*resolution-1)*sqrt(sqrt(5.)*phi))
  def ideal_particle_mass(self):
    self.update_common_quantities()
    idealmass = self.dist_factor*self.neighbour_distance*self.star.radius*self.mass_rate/(self.particles_per_sphere*self.velocity)
    return idealmass
  def ideal_resolution(self):
    self.update_common_quantities()
    idealnumberofparticles = self.dist_factor*self.neighbour_distance*self.star.radius*self.mass_rate/(self.mass_of_particles*self.velocity)
    return int(round((10.+sqrt(10.*(idealnumberofparticles-2)))/20.))
  def inject(self, simulation, time, dtlast):
    self.star.update_properties()
    outer_sphere = int((time-dtlast)/self.time_between_spheres)+1
    inner_sphere = int(time/self.time_between_spheres)
    inner_handled_sphere = inner_sphere+self.handled_spheres
    for i in range(inner_handled_sphere, inner_sphere, -1):
      self.inject_sphere(simulation, i, time, self.particles_per_sphere*(inner_handled_sphere-i), boundary=True)
    for i in range(inner_sphere, outer_sphere-1, -1):
      npart = simulation.get_npart()
      self.inject_sphere(simulation, i, time, npart)
  def inject_sphere(self, simulation, i, time, first_part, boundary=False):
    local_time = time - (float(i)-self.shift)*self.time_between_spheres
    radius, velocity, temperature = self.extrapolate_sphere(local_time)
    specific_energy = temperature*self.specific_energy_to_temperature_ratio
    angles = array([3.05231445647236, 0.397072776282339, 2.27500616856518])*float(i)
    h = radius * self.neighbour_distance * self.hfact
    simulation.inject_or_update_sphere(first_part, self.resolution, self.star.location, radius, self.star.velocity, velocity, h, specific_energy, angles, boundary)
  def extrapolate_sphere(self, local_time):
    equation = integrate.ode(self.theoretical_dv_dt)
    equation.set_integrator('vode', atol=1.e-10, rtol=1.e-10)
    equation.set_initial_value([self.injection_radius, self.velocity], 0.)
    dt = local_time/100.
    while equation.successful() and abs(equation.t) < abs(local_time):
      equation.integrate(equation.t+dt)
    r, v = equation.y
    T = self.temperature * (self.injection_radius**2 * self.velocity / (r**2 * v))**(self.gamma-1.)
    return r, v, T
  def theoretical_dv_dt(self, t, rv):
    r, v = rv
    dv_dr = self.theoretical_dv_dr(r, v)
    return [v, v*dv_dr]
  def theoretical_dv_dr(self, r, v):
    T = self.temperature * (self.injection_radius**2 * self.velocity / (r**2 * v))**(self.gamma-1.)
    vs2 = self.gamma * R * T / self.mu
    return (-G*self.star.mass/r**2+2.*vs2/r)/(v-vs2/v)
  def theoretical_solution(self, rmax):
    N = 1000
    dr = (rmax-self.injection_radius)/float(N)
    equation = integrate.ode(self.theoretical_dv_dr)
    equation.set_integrator('vode', atol=1.e-6, rtol=1.e-6)
    equation.set_initial_value(self.velocity, self.injection_radius)
    r = zeros(N)
    v = zeros(N)
    i = 0
    while equation.successful() and i<N:
      equation.integrate(equation.t+dr)
      r[i],v[i] = equation.t, equation.y
      i += 1
    T = self.temperature * (self.injection_radius**2 * self.velocity / (r**2 * v))**(self.gamma-1.)
    return r, v, T
