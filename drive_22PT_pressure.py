#!/usr/bin/env python

class PlasmaConfig:
    def __init__(self):
        self.u_probe = -0.5
        self.He_pressure = 1e-3

    def __str__(self):
        string =  """
# sample config file for plasma2d
# comments start with #

# geometry
coord = CARTESIAN

# probe voltage
# or rf trap amplitude
u_probe =  20.0		#[V]
rf_amplitude = 20.0

#RF trap ?
rf = 1

#RF trap frequency
rf_omega = 11.3e7

# probe dimensions
probe_radius = 1e-4		#[m]
probe_length = 1e-2		#[m]

# detect particles?
has_probe = 1

# x component of external field
extern_field = 500.0		#[V/m]

#magnetic field
#magnetic_field_file = MagPole.txt
magnetic_field_const = 1
Bt = 0.03

# total number of particles in simulation
n_particles_total = 40000
macroparticle_factor = 1

# mesh dimension
r_sampl = 100
z_sampl = 100

# working area dimension
r_max = 2e-2	#[m]
z_max = 2e-2	#[m]

# number of timesteps
niter = 200000

# the factor by which the source area is smaller
# than working area, positive integer
src_fact = 4

# one for selfconsistent simulation, 0 for simulation
# in constant homogenous extern_field
selfconsistent = 0

# simulation state is printed each t_print timesteps
t_print = 10000

# after t_equilib timesteps the results are averaged
t_equilib = 100000

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
pressure = ${* 133.0 0.0}		#[Pa]
temperature = 300.0		#[K]
dt = 0				#[s]

[HELIUM]
"""
        string += "pressure = ${* 133.0 " +  str(self.He_pressure) + """}		#[Pa]
temperature = 300.0		#[K]
dt = 0				#[s]

[H_NEG]
density = 1e15			#[m^-3]
temperature = 300.0		#[K]
dt = 2e-9			#[s]

[HYDROGEN]
pressure = ${* 133.0 0.0}		#[Pa]
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
energy = 3.0			#[eV]
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


pressures = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6]
#import os
i=-1
config = PlasmaConfig()
for p in pressures :
    i += 1
    config.He_pressure = p
    configfile = "config_%d.txt" % i
    f = open(configfile,"w")
    print >>f, config
    f.close()
    command = "echo 'y'|./test config=config_%d.txt output_dir=%.0eTorr" % (i,p)
    #if (i+1) % 2 != 0 :
	#command += " &"	
    print command
    #os.system(command)
