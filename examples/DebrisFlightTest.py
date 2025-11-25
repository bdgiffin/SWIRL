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
import time as perftime

# ---------------------------------------------------------------------------- #

# Define generic problem dimensional parameters:

# Define the x,y coordinate center of the vortex
x_vortex = 5.0 # [m]
y_vortex =  0.0 # [m]

# Define the x,y coordinate center of the tower
x_tower = 0.0 # [m]
y_tower = 0.0 # [m]

# Define the base width and height of the tower
base_width_tower = 3.5  # [m]
height_tower     = 20.0 # [m]

# Pre-compute relevant derived dimensional quantities:

# Determine the (radial) distace between the tower and the vortex
distance_to_tower_from_vortex = sqrt(pow(x_tower-x_vortex,2) + pow(y_tower-y_vortex,2))

# Determine the directional position of the center of the tower (as an angle measured from the x-axis) relative to the center of the vortex
direction_to_tower_from_vortex = atan2(y_tower-y_vortex, x_tower-x_vortex)

# Specify the position and dimensions of the cylinder containing randomly generated particles, such that:
#  1) The particle cylinder center is at the vortex center
particle_cylinder_center = [x_vortex,y_vortex,0.0] # [m,m,m]
#  2) The particle cylinder radius is defined to encompass the tower
particle_cylinder_radius = distance_to_tower_from_vortex + base_width_tower # [m]
#  3) The particle cylinder height is greater than the height of the tower
particle_cylinder_height = 1.0*height_tower # [m] 1.5

# Determine the radius of the inner cylinder for particle generation (particles will NOT be initialized within the inner cylinder)
particle_cylinder_inner_radius = max(distance_to_tower_from_vortex - 3.0*base_width_tower, 1.0e-16)

# Define the "exclusion cylinder" inside of which no particles should be initialized (this should closely enclose the tower)
exclusion_cylinder_center = [x_tower,y_tower,0.0] # [m,m,m]
exclusion_cylinder_radius = 0.5*base_width_tower # [m]
exclusion_cylinder_height = height_tower # [m]

# ---------------------------------------------------------------------------- #

# Call C/C++ library API functions from Python:

# Define randomized spherical particle parameters
num_particles = 1000
particle_density        =  0.5 # [kg/m^3] (roughly the density of wood)
particle_min_diameter   = 0.01 # [m]
particle_diameter_range =  1.0 # [m]
random_seed = 1
average_radial_coordinate   = 2.0 # [m] the mean radial coordinate at which particles are likely to be positioned (determined according to equation 36 in Baker & Sterling paper: https://www.sciencedirect.com/science/article/pii/S0167610517301174)
average_vertical_coordinate = 1.0 # [m] the mean vertical coordinate at which particles are likely to be positioned (determined according to equation 37 in Baker & Sterling paper: https://www.sciencedirect.com/science/article/pii/S0167610517301174)
radial_coordinate_stddev    = 0.25*average_radial_coordinate  # [m] The standard deviation characterizing the lognormal distribution of particles' initial radial coordinates   (this example assumes a COV of 0.25)
vertical_coordinate_stddev  = 0.5*average_vertical_coordinate # [m] The standard deviation characterizing the lognormal distribution of particles' initial vertical coordinates (this example assumes a COV of 0.5)
#SWIRL.create_random_particles(num_particles,particle_density,particle_min_diameter,particle_diameter_range,particle_cylinder_radius,particle_cylinder_height,particle_cylinder_center,random_seed)
SWIRL.create_constrained_random_particles(num_particles,particle_density,particle_min_diameter,particle_diameter_range,particle_cylinder_radius,particle_cylinder_height,particle_cylinder_center,random_seed,particle_cylinder_inner_radius,average_radial_coordinate,radial_coordinate_stddev,average_vertical_coordinate,vertical_coordinate_stddev,exclusion_cylinder_center,exclusion_cylinder_radius,exclusion_cylinder_height)

# Create the parameterized wind field model (Baker Sterling Vortex)
wind_field_params = np.zeros(12)
wind_field_params[0]  = 100.0    # [m/s]      Um: reference radial velocity
wind_field_params[1]  = 1.0      # [m]        rm: reference radius
wind_field_params[2]  = 4.0      # [m]        zm: reference height
wind_field_params[3]  = 2.0      #             S: swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
wind_field_params[4]  = 2.0      #         gamma: 
wind_field_params[5]  = 1.293    # [kg/m^3] rho0: reference density of air at STP
wind_field_params[6]  = x_vortex # [m]       xc0: x-position of the vortex center
wind_field_params[7]  = y_vortex # [m]       yc0: y-position of the vortex center
wind_field_params[8]  = 0.0      # [m]       zc0: z-position of the vortex center
wind_field_params[9]  = 0.0      # [m/s]     vxc: x-velocity of the vortex center
wind_field_params[10] = 0.0      # [m/s]     vyc: y-velocity of the vortex center
wind_field_params[11] = 0.0      # [m/s]     vzc: z-velocity of the vortex center
SWIRL.API.define_wind_field(b"BakerSterlingVortex",wind_field_params)

# Modify flag to activate particles
# (setting "particles_active" to "False" will freeze all particles and prevent contact forces from being applied to the structure - this should be done prior to the start of the static analysis initialization)
# (setting "particles_active" to "True" will unfreeze all particles and apply contact forces to the structures - this should be done prior to the start of the dynamic analysis)
SWIRL.API.set_particles_active(False)

# OPTIONAL: define bounding wedge-shaped control volume with periodic inflow/outflow BCs to contain all particles
use_particle_cv = True
if (use_particle_cv):
    # Specify the x,y coordinate center of the wedge-shaped particle control volume (should be the same as the center of the vortex):
    cv_x0 = x_vortex # [m] (x-center of the vortex)
    cv_y0 = y_vortex # [m] (y-center of the vortex)

    # Specify the angular thickness of the wedge-shaped control volume:
    # (this should be determined based on the tower dimensions, as well as the relative position of the tower compared to the vortex)
    # (in this example, the width of the wedged-shaped CV is 5.0 times the base width of the tower, and is capped at a max of 2pi radians)
    cv_dtheta = min(5.0*base_width_tower/(distance_to_tower_from_vortex+1.0e-16), 2.0*pi) # [radians]

    # Specify the angular coordinate of the inflow boundary plane for the wedge-shaped control volume:
    # (this should be determined based on the relative position of the tower compared to the vortex)
    cv_theta_in = direction_to_tower_from_vortex - 0.5*cv_dtheta; # [radians]
    # NOTE: the outflow boundary plane is positioned at an angular coordinate of: cv_theta_out = cv_theta_in + cv_dtheta

    # Define the wedge-shaped control volume in SWIRL (pass relevant parameters to corresponding SWIRL API function):
    SWIRL.API.define_particle_control_volume(cv_x0,cv_y0,cv_theta_in,cv_dtheta)

# RUN analysis -------------------------------------------------------------

# Wall time for performance measurement
start_time = perftime.perf_counter()

# perform the analysis
time = 0.0 # [s] starting time
dt   = 0.001 # [s] time increment
SWIRL.API.update_state(time) # required for initialization
SWIRL.output_state(time)
for step_id in range(1,100):
    time = time + dt
    SWIRL.API.update_state(time)
    SWIRL.output_state(time)

# Modify flag to activate particles
# (setting "particles_active" to "False" will freeze all particles and prevent contact forces from being applied to the structure - this should be done prior to the start of the static analysis initialization)
# (setting "particles_active" to "True" will unfreeze all particles and apply contact forces to the structures - this should be done prior to the start of the dynamic analysis)
SWIRL.API.set_particles_active(True)
for step_id in range(1,100):
    time = time + dt
    SWIRL.API.update_state(time)
    SWIRL.output_state(time)

# Report wall time for performance measurement
end_time = perftime.perf_counter()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time:.6f} seconds")
    
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
