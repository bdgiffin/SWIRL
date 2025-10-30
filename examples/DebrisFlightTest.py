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

# Post-process results -----------------------------------------------------

# Pre-allocate particle data arrays for use during output
x = np.zeros(num_particles)
y = np.zeros(num_particles)
z = np.zeros(num_particles)

# retrieve the simulation state info at the current time
SWIRL.API.get_particle_positions(x,y,z)

# transform particle positions into cylindrical coordinates
r = np.zeros(num_particles)
for i in range(0,num_particles):
    r[i] = sqrt(x[i]*x[i]+y[i]*y[i])

# compute logarithm of cylindrical particle coordinates
logX = np.zeros((num_particles,2))
for i in range(0,num_particles):
    logX[i,0] = log(max(r[i],1.0e-6))
    logX[i,1] = log(max(z[i],1.0e-6))

# create a bivariate (log) normal distribution
mean_logX = np.mean(logX,axis=0)
print('mean:')
print(mean_logX)
cov_logX  = np.cov(logX,rowvar=0)
print('covariance:')
print(cov_logX)

# randomly sample from the fitted distribution
num_samples = 1000
samples = np.exp(multivariate_normal.rvs(mean=mean_logX, cov=cov_logX, size=num_samples, random_state=None))

# wipe the SWIRL module (reset all internal state data structures)
SWIRL.wipe()

# plot results
xgrid, ygrid = np.mgrid[-1:max(logX[:,0]):.01, -1:max(logX[:,1]):.01]
pos = np.dstack((xgrid, ygrid))
fig = plt.figure()
ax = fig.add_subplot(111)
ax.contourf(np.exp(xgrid), np.exp(ygrid), multivariate_normal.pdf(pos, mean=mean_logX, cov=cov_logX))
plt.plot(r, z, 'o')
plt.plot(samples[:,0], samples[:,1], 'o')
ax.axis('equal')
ax.set_xlim([0, particle_cylinder_radius])
ax.set_ylim([0, particle_cylinder_height])
plt.show()
