#!/usr/bin/env python

class PlasmaConfig:
    def __init__(self):
        self.u_probe = -0.5

    def __str__(self):
        string =  """
# sample config file for plasma2d
# comments start with #

# geometry
coord = CYLINDRICAL

# probe voltage
"""
        string += "u_probe = " + str(self.u_probe) + """		#[V]

#RF trap ?
rf = 0

# probe dimensions
probe_radius = 1e-4		#[m]
probe_length = 1e-2		#[m]

# x component of external field
extern_field = 500.0		#[V/m]

# detect particles?
has_probe = 1

#magnetic field
magnetic_field_file = magfield2.dat
magnetic_field_const = 0
Bt = 0.03

# total number of particles in simulation
n_particles_total = 400000
macroparticle_factor = 2000

# mesh dimension
r_sampl = 400
z_sampl = 200

# working area dimension
r_max = 0.05	#[m]
z_max = 0.4	#[m]

# number of timesteps
niter = 100000

# the factor by which the source area is smaller
# than working area, positive integer
src_fact = 4

# one for selfconsistent simulation, 0 for simulation
# in constant homogenous extern_field
selfconsistent = 0

# simulation state is printed each t_print timesteps
t_print = 200

# after t_equilib timesteps the results are averaged
t_equilib = 1000000

# obvious
do_plot = 0

# if particle_reload==1 simulation state is reloaded
# from directory particle reload_dir
particle_reload = 0
particle_reload_dir = 1e-10s

# method of advancing particles, possible
# values are 
# ADVANCE2, ADVANCE2_EULER, ADVANCE_OLD, ADVANCE_PROB
advancer = ADVANCE2_EULER

# specification of plasma parameters
# energy is defined as 3.0/2.0*k_B*temperature !!!
# or at least should be...

# you can specify temperature or energy, it is the
# same but with different units
# similarly pressure or density specifies the same
[ARGON]
pressure = ${* 133.0 1e-3}		#[Pa]
temperature = 300.0		#[K]
dt = 0				#[s]

[ARGON_POS]
density = 1e15			#[m^-3]
temperature = 300.0		#[K]
dt = 1e-8			#[s]

[O2_POS]
density = 0e15			#[m^-3]
temperature = 300.0		#[K]
dt = 1e-8			#[s]

[O2_NEG]
density = 0e15			#[m^-3]
temperature = 300.0		#[K]
dt = 1e-8			#[s]

[ELECTRON]
density = 1e15			#[m^-3]
energy = 1.5			#[eV]
dt = 5e-11			#[s]

[O2]
pressure = 0.0			#[Pa]
temperature = 300.0		#[K]
dt = 1				#[s]

[O_NEG]
density = 0e15			#[m^-3]
temperature = 300.0		#[K]
dt = 1e-8			#[s]
"""
        return string


voltages = [-1.5, -1.2, -1.1, 1.05, -1.0, -0.95, -0.9, -0.8, -0.7, -0.5, -0.3]

#import os
i=-1
config = PlasmaConfig()
for u in voltages :
    i += 1
    config.u_probe = u
    configfile = "config_%d.txt" % i
    f = open(configfile,"w")
    print >>f, config
    f.close()
    command = "echo 'y'|./test config=config_%d.txt output_dir=%gV" % (i,u)
    #if (i+1) % 2 != 0 :
	#command += " &"	
    print command
    #os.system(command)
