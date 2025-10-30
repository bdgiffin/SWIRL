# Append the location of the locally installed SWIRL package to sys.path
import sys
sys.path.append("../install/package/")

# Load the SWIRL package
import SWIRL

# Python package for reading/writing data in the Exodus mesh database format
# NOTE: PYEXODUS V0.1.5 NEEDS TO BE MODIFIED TO WORK CORRECTLY WITH PYTHON 3.12
# https://pypi.org/project/pyexodus/
import pyexodus

# Other needed Python packages
from math import *
import numpy as np
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
import os

def run_SWIRL():
    # ---------------------------------------------------------------------------- #

    # Call C/C++ library API functions from Python:

    # Define randomized spherical particle parameters
    num_particles = 1000
    particle_density         =  0.5 # [kg/m^3] (roughly the density of wood)
    particle_min_diameter    = 0.01 # [m]
    particle_diameter_range  =  1.0 # [m]
    particle_cylinder_radius = 40.0 # [m]
    particle_cylinder_height = 40.0 # [m]
    particle_cylinder_center = [0.0,0.0,0.0] # [m,m,m]
    random_seed = 1
    SWIRL.create_random_particles(num_particles,particle_density,particle_min_diameter,particle_diameter_range,particle_cylinder_radius,particle_cylinder_height,particle_cylinder_center,random_seed)
    
    # Create the parameterized wind field model (Baker Sterling Vortex)
    wind_field_params = np.zeros(12)
    wind_field_params[0]  = 100.0 # [m/s]      Um: reference radial velocity
    wind_field_params[1]  = 1.0   # [m]        rm: reference radius
    wind_field_params[2]  = 4.0   # [m]        zm: reference height
    wind_field_params[3]  = 2.0   #             S: swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
    wind_field_params[4]  = 2.0   #         gamma: 
    wind_field_params[5]  = 1.293 # [kg/m^3] rho0: reference density of air at STP
    wind_field_params[6]  = 10.0  # [m]       xc0: x-position of the vortex center
    wind_field_params[7]  = 0.0   # [m]       yc0: y-position of the vortex center
    wind_field_params[8]  = 0.0   # [m]       zc0: z-position of the vortex center
    wind_field_params[9]  = 0.0   # [m/s]     vxc: x-velocity of the vortex center
    wind_field_params[10] = 0.0   # [m/s]     vyc: y-velocity of the vortex center
    wind_field_params[11] = 0.0   # [m/s]     vzc: z-velocity of the vortex center
    SWIRL.API.define_wind_field(b"BakerSterlingVortex",wind_field_params)

    # RUN analysis -------------------------------------------------------------

    # perform the analysis
    time = 0.0 # [s] starting time
    dt   = 0.001 # [s] time increment
    SWIRL.API.update_state(time) # required for initialization
    SWIRL.output_state(time)
    for step_id in range(1,100):
        time = time + dt
        SWIRL.API.update_state(time)
        SWIRL.output_state(time)

    # finalize the SWIRL module (close the Exodus files)
    SWIRL.finalize()

    # finalize the SWIRL module (restart SWIRL)
    SWIRL.wipe()

# Test running multiple instances of SWIRL
print("Running first analysis")
run_SWIRL()
os.rename("particles.exo", "particles1.exo")

print("Running second analysis")
run_SWIRL()
os.rename("particles.exo", "particles2.exo")
