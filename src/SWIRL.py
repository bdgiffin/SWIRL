# Module for calling the C/C++ SWIRL API functions from Python

# Python package for calling C/C++ functions from Python
from ctypes import CDLL, POINTER, byref
from ctypes import c_size_t, c_double, c_int, c_char_p, c_bool

# Package for reading/writing mesh files
import pyexodus

# Other needed Python packages
import sys
import os
import math
import time as timer
from datetime import timedelta
import numpy as np
from argparse import ArgumentParser
from scipy.stats import truncnorm

# ---------------------------------------------------------------------------- #

# Module initialization:

# Load the pre-compiled external C/C++ "shared object" libraries
library_name = os.path.join(os.path.dirname(__file__), 'libSWIRL.so')
API = CDLL(library_name)

# Define types to convert Numpy arrays into C arrays:

# C-type corresponding to 1D numpy array
ND_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.float64, 
                                      ndim=1,
                                      flags="C")
NI_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.int32, 
                                      ndim=1,
                                      flags="C")
NI_POINTER_2 = np.ctypeslib.ndpointer(dtype=np.int32, 
                                      ndim=2,
                                      flags="C")

# Define all C/C++ library API function signatures
API.define_parameter.argtypes = [c_char_p, c_double]
API.define_parameter.restype  = None
API.define_wind_field.argtypes = [c_char_p, ND_POINTER_1]
API.define_wind_field.restype  = None
API.define_particles.argtypes = [c_size_t, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.define_particles.restype  = None
API.set_particles_active.argtypes = [c_bool]
API.set_particles_active.restype  = None
API.define_particle_control_volume.argtypes = [c_double, c_double, c_double, c_double]
API.define_particle_control_volume.restype  = None
API.define_members.argtypes = [c_size_t, c_size_t, NI_POINTER_2, c_size_t, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.define_members.restype  = None
API.update_state.argtypes = [c_double]
API.update_state.restype  = None
API.get_particle_field_data.argtypes = [ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.get_particle_field_data.restype  = None
API.get_particle_positions.argtypes = [ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.get_particle_positions.restype  = None
API.get_wind_field_data.argtypes = [c_size_t, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.get_wind_field_data.restype  = None
API.get_impact_event_metrics.argtypes = [POINTER(c_int), POINTER(c_double), POINTER(c_double)]
API.get_impact_event_metrics.restype  = None
API.get_impact_events.argtypes = [NI_POINTER_1, NI_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.get_impact_events.restype  = None
API.finalize.argtypes = None
API.finalize.restype  = None

# ---------------------------------------------------------------------------- #

# Generate a log-normally distributed collection of particles with variable diameters and positions within a constrained range,
# and initialize the SWIRL object with these particles prior to initialization
def create_constrained_random_particles(n_particles,density,min_diameter,diameter_range,cylinder_radius,cylinder_height,cylinder_center,random_seed,inner_cylinder_radius,average_radial_coordinate,radial_coordinate_stddev,average_vertical_coordinate,vertical_coordinate_stddev,exclusion_cylinder_center,exclusion_cylinder_radius,exclusion_cylinder_height):
    # n_particles     [int]    The total number of spherical debris particles to be defined
    # density         [kg/m^3] The constant mass density assigned to all particles
    # min_diameter    [m]      The minimum particle diameter
    # diameter_range  [m]      The range of random particle diameters
    # cylinder_radius [m]      The radius of the cylinder in which particles will be randomly distributed
    # cylinder_height [m]      The height of the cylinder in which particles will be randomly distributed
    # cylinder_center [m,m,m]  The x,y,z coordinate center of cylinder at the time of initialization
    # random_seed     [int]    The random number generator seed value, to ensure reproducibility
    # Extra variables:
    # inner_cylinder_radius       [m] The radius of the inner cylinder (particles will NOT be initialized within the inner cylinder)
    # average_radial_coordinate   [m] The mean radial coordinate at which particles are likely to be positioned
    # radial_coordinate_stddev    [m] The standard deviation characterizing the lognormal distribution of particles' initial radial coordinates
    # average_vertical_coordinate [m] The mean vertical coordinate at which particles are likely to be positioned
    # vertical_coordinate_stddev  [m] The standard deviation characterizing the lognormal distribution of particles' initial vertical coordinates
    # exclusion_cylinder_center   [m,m,m] The The x,y,z coordinate center of cylinder inside of which no particles should be initialized
    # exclusion_cylinder_radius   [m] The radius of the cylinder in which no particles should be initialized
    # exclusion_cylinder_height   [m] The height of the cylinder in which no particles should be initialized

    global num_particles
    num_particles = n_particles

    # Instantiate a random number generator (rng) initialized with the provided random_seed value
    rng = np.random.default_rng(random_seed)

    # Generate a random collection of spherical particles with variable diameters but constant density
    diameters = min_diameter + diameter_range*rng.random(n_particles) # [m]
    masses = np.zeros(n_particles)
    for i in range(0,n_particles):
        iradius = 0.5*diameters[i]
        masses[i] = density*(4.0/3.0)*math.pi*iradius*iradius*iradius # [kg] (assuming roughly spherical shape)

    # Randomize the initial positions of all particles
    rng_state = np.random.RandomState(random_seed)
    position_x = np.zeros(n_particles)
    position_y = np.zeros(n_particles)
    position_z = np.zeros(n_particles)
    radial_position   = sample_constrained_lognormal(average_radial_coordinate, radial_coordinate_stddev, inner_cylinder_radius, cylinder_radius, n_particles, rng_state)
    circum_position   = 2.0*math.pi*rng.random(n_particles)
    vertical_position = sample_constrained_lognormal(average_vertical_coordinate, vertical_coordinate_stddev, 1.0e-16, cylinder_height, n_particles, rng_state)
    for i in range(0,n_particles):
        position_x[i] = cylinder_center[0] + radial_position[i]*math.cos(circum_position[i])
        position_y[i] = cylinder_center[1] + radial_position[i]*math.sin(circum_position[i])
        position_z[i] = cylinder_center[2] + vertical_position[i]

        # Repeatedly check if particle lies inside of the exclusion cylinder (CAUTION: this could result in an infinite loop if the exclusion cylinder is too big)
        while ((math.sqrt(pow(position_x[i]-exclusion_cylinder_center[0],2)+pow(position_y[i]-exclusion_cylinder_center[1],2)) < exclusion_cylinder_radius) and ((position_z[i] - exclusion_cylinder_center[2]) < exclusion_cylinder_height)):
            # Attempt to re-initialize the particle's position
            iradial_position   = sample_constrained_lognormal(average_radial_coordinate, radial_coordinate_stddev, inner_cylinder_radius, cylinder_radius, 1, rng_state)
            icircum_position   = 2.0*math.pi*rng.random(1)
            ivertical_position = sample_constrained_lognormal(average_vertical_coordinate, vertical_coordinate_stddev, 1.0e-16, cylinder_height, 1, rng_state)
            position_x[i] = cylinder_center[0] + iradial_position*math.cos(icircum_position)
            position_y[i] = cylinder_center[1] + iradial_position*math.sin(icircum_position)
            position_z[i] = cylinder_center[2] + ivertical_position

    # Call SWIRL initialization API function
    API.define_particles(n_particles,masses,diameters,position_x,position_y,position_z)

    # Create the Exodus file containing particle info:
    if (num_particles > 0):

        # create a new Exodus file
        filename = "particles.exo"
        try:
            os.remove(filename)
        except OSError:
            pass
        global exo
        exo = pyexodus.exodus(file=filename, mode='w', array_type='numpy', title='Debris particle trajectory time-history file - produced by SWIRL module', numDims=3, numNodes=num_particles, numElems=num_particles, numBlocks=1, numNodeSets=0, numSideSets=0, io_size=0, compression=None)

        # put node coordinates
        exo.put_coords(xCoords=position_x,yCoords=position_y,zCoords=position_z)

        # put element block info for all particles
        exo.put_elem_blk_info(id=1, elemType='SPHERE', numElems=num_particles, numNodesPerElem=1, numAttrsPerElem=0)
        exo.put_elem_connectivity(id=1, connectivity=np.arange(num_particles), shift_indices=1, chunk_size_in_mb=128)

        # set the number of output node (particle) variables and their names
        num_node_variables = 3 + 3 + 3
        exo.set_node_variable_number(num_node_variables)
        exo.put_node_variable_name("displacement_x", 1)
        exo.put_node_variable_name("displacement_y", 2)
        exo.put_node_variable_name("displacement_z", 3)
        exo.put_node_variable_name("velocity_x",     4)
        exo.put_node_variable_name("velocity_y",     5)
        exo.put_node_variable_name("velocity_z",     6)
        exo.put_node_variable_name("force_x",        7)
        exo.put_node_variable_name("force_y",        8)
        exo.put_node_variable_name("force_z",        9)

        # initialize the total number of time states
        global step_id
        step_id = 0

# ---------------------------------------------------------------------------- #

# Generate a randomized collection of particles with variable diameters and positions,
# and initialize the SWIRL object with these particles prior to initialization
def create_random_particles(n_particles,density,min_diameter,diameter_range,cylinder_radius,cylinder_height,cylinder_center,random_seed):
    # n_particles     [int]    The total number of spherical debris particles to be defined
    # density         [kg/m^3] The constant mass density assigned to all particles
    # min_diameter    [m]      The minimum particle diameter
    # diameter_range  [m]      The range of random particle diameters
    # cylinder_radius [m]      The radius of the cylinder in which particles will be randomly distributed
    # cylinder_height [m]      The height of the cylinder in which particles will be randomly distributed
    # cylinder_center [m,m,m]  The x,y,z coordinate center of cylinder at the time of initialization
    # random_seed     [int]    The random number generator seed value, to ensure reproducibility

    global num_particles
    num_particles = n_particles

    # Instantiate a random number generator (rng) initialized with the provided random_seed value
    rng = np.random.default_rng(random_seed)

    # Generate a random collection of spherical particles with variable diameters but constant density
    diameters = min_diameter + diameter_range*rng.random(n_particles) # [m]
    masses = np.zeros(n_particles)
    for i in range(0,n_particles):
        iradius = 0.5*diameters[i]
        masses[i] = density*(4.0/3.0)*math.pi*iradius*iradius*iradius # [kg] (assuming roughly spherical shape)

    # Randomize the initial positions of all particles
    position_x = np.zeros(n_particles)
    position_y = np.zeros(n_particles)
    for i in range(0,n_particles):
        iradial_position = cylinder_radius*rng.random(1)
        icircum_position = 2.0*math.pi*rng.random(1)
        position_x[i] = cylinder_center[0] + iradial_position[0]*math.cos(icircum_position[0])
        position_y[i] = cylinder_center[1] + iradial_position[0]*math.sin(icircum_position[0])
    position_z = cylinder_center[2] + cylinder_height*rng.random(n_particles)

    # Call SWIRL initialization API function
    API.define_particles(n_particles,masses,diameters,position_x,position_y,position_z)

    # Create the Exodus file containing particle info:
    if (num_particles > 0):

        # create a new Exodus file
        filename = "particles.exo"
        try:
            os.remove(filename)
        except OSError:
            pass
        global exo
        exo = pyexodus.exodus(file=filename, mode='w', array_type='numpy', title='Debris particle trajectory time-history file - produced by SWIRL module', numDims=3, numNodes=num_particles, numElems=num_particles, numBlocks=1, numNodeSets=0, numSideSets=0, io_size=0, compression=None)

        # put node coordinates
        exo.put_coords(xCoords=position_x,yCoords=position_y,zCoords=position_z)

        # put element block info for all particles
        exo.put_elem_blk_info(id=1, elemType='SPHERE', numElems=num_particles, numNodesPerElem=1, numAttrsPerElem=0)
        exo.put_elem_connectivity(id=1, connectivity=np.arange(num_particles), shift_indices=1, chunk_size_in_mb=128)

        # set the number of output node (particle) variables and their names
        num_node_variables = 3 + 3 + 3
        exo.set_node_variable_number(num_node_variables)
        exo.put_node_variable_name("displacement_x", 1)
        exo.put_node_variable_name("displacement_y", 2)
        exo.put_node_variable_name("displacement_z", 3)
        exo.put_node_variable_name("velocity_x",     4)
        exo.put_node_variable_name("velocity_y",     5)
        exo.put_node_variable_name("velocity_z",     6)
        exo.put_node_variable_name("force_x",        7)
        exo.put_node_variable_name("force_y",        8)
        exo.put_node_variable_name("force_z",        9)

        # initialize the total number of time states
        global step_id
        step_id = 0

# ---------------------------------------------------------------------------- #

# Get impact event info: (Nimpacts, max_force, max_impulse)
def impact_event_metrics():
    Nimpacts    = c_int(0)
    max_force   = c_double(0.0)
    max_impulse = c_double(0.0)
    API.get_impact_event_metrics(byref(Nimpacts),byref(max_force),byref(max_impulse))
    return Nimpacts.value, max_force.value, max_impulse.value

# ---------------------------------------------------------------------------- #

# Output information for all discrete impact events
def output_impact_events():
    # Get the total number of impacts
    Nimpacts, max_force, max_impulse = impact_event_metrics()

    # Pre-allocate data arrays for impact events
    particle_ID = np.zeros(Nimpacts,dtype=np.int32)
    segment_ID  = np.zeros(Nimpacts,dtype=np.int32)
    x           = np.zeros(Nimpacts)
    y           = np.zeros(Nimpacts)
    z           = np.zeros(Nimpacts)
    start_time  = np.zeros(Nimpacts)
    end_time    = np.zeros(Nimpacts)
    impulse     = np.zeros(Nimpacts)
    max_force   = np.zeros(Nimpacts)

    # Get data for all impact events
    API.get_impact_events(particle_ID,segment_ID,x,y,z,start_time,end_time,impulse,max_force)

    # Create the Exodus file containing impact info:
    if (Nimpacts > 0):

        # create a new Exodus file
        filename = "impact_events.exo"
        try:
            os.remove(filename)
        except OSError:
            pass
        impact_exo = pyexodus.exodus(file=filename, mode='w', array_type='numpy', title='Impact event file - produced by SWIRL module', numDims=3, numNodes=Nimpacts, numElems=Nimpacts, numBlocks=1, numNodeSets=0, numSideSets=0, io_size=0, compression=None)

        # put node coordinates (locations of impacts)
        impact_exo.put_coords(xCoords=x,yCoords=y,zCoords=z)

        # put element block info for all impacts
        impact_exo.put_elem_blk_info(id=1, elemType='SPHERE', numElems=Nimpacts, numNodesPerElem=1, numAttrsPerElem=0)
        impact_exo.put_elem_connectivity(id=1, connectivity=np.arange(Nimpacts), shift_indices=1, chunk_size_in_mb=128)

        # set the number of output impact event variables and their names
        impact_exo.set_node_variable_number(6)
        impact_exo.put_node_variable_name("particle_ID",1)
        impact_exo.put_node_variable_name("segment_ID", 2)
        impact_exo.put_node_variable_name("start_time", 3)
        impact_exo.put_node_variable_name("end_time",   4)
        impact_exo.put_node_variable_name("impulse",    5)
        impact_exo.put_node_variable_name("max_force",  6)
    
        # create a new (single) output time state
        impact_exo.put_time(1, 0.0)
    
        # write nodal variable values at the current time state
        impact_exo.put_node_variable_values("particle_ID",1, particle_ID.astype(np.float64))
        impact_exo.put_node_variable_values("segment_ID", 1, segment_ID.astype(np.float64))
        impact_exo.put_node_variable_values("start_time", 1, start_time)
        impact_exo.put_node_variable_values("end_time",   1, end_time)
        impact_exo.put_node_variable_values("impulse",    1, impulse)
        impact_exo.put_node_variable_values("max_force",  1, max_force)

        # Close the Exodus file
        impact_exo.close()

# ---------------------------------------------------------------------------- #

# Write data to the Exodus file containing particle info for the current time state
def output_state(time):

    # Write data to the Exodus file for all particles
    if (num_particles > 0):

        # increment the total number of time states
        global step_id
        step_id = step_id + 1

        # Pre-allocate particle data arrays for use during output
        ux = np.zeros(num_particles)
        uy = np.zeros(num_particles)
        uz = np.zeros(num_particles)
        vx = np.zeros(num_particles)
        vy = np.zeros(num_particles)
        vz = np.zeros(num_particles)
        fx = np.zeros(num_particles)
        fy = np.zeros(num_particles)
        fz = np.zeros(num_particles)

        # retrieve the simulation state info at the current time
        API.get_particle_field_data(ux,uy,uz,vx,vy,vz,fx,fy,fz)
    
        # create a new output time state
        exo.put_time(step_id, time)
    
        # write nodal variable values at the current time state
        exo.put_node_variable_values("displacement_x", step_id, ux)
        exo.put_node_variable_values("displacement_y", step_id, uy)
        exo.put_node_variable_values("displacement_z", step_id, uz)
        exo.put_node_variable_values("velocity_x",     step_id, vx)
        exo.put_node_variable_values("velocity_y",     step_id, vy)
        exo.put_node_variable_values("velocity_z",     step_id, vz)
        exo.put_node_variable_values("force_x",        step_id, fx)
        exo.put_node_variable_values("force_y",        step_id, fy)
        exo.put_node_variable_values("force_z",        step_id, fz)

# ---------------------------------------------------------------------------- #

# Close the connection to the Exodus files
def finalize():

    # Close the Exodus file
    if (num_particles > 0):
        exo.close()

# ---------------------------------------------------------------------------- #

# Close the connection to the SWIRL library
def wipe():

    # Call SWIRL finalization API function
    API.finalize()

# ---------------------------------------------------------------------------- #

def sample_constrained_lognormal(mu, sigma, lower_bound, upper_bound, Nsamples, rng_state):
    """
    Samples from a log-normal distribution in a constrained range using scipy.stats.truncnorm.
    """
    # The bounds of the log-normal distribution must be transformed to 
    # the bounds of the underlying normal distribution using the natural logarithm.
    log_lower_bound = np.log(lower_bound)
    log_upper_bound = np.log(upper_bound)

    # Normalize the log-bounds to the standard normal scale (z-scores)
    a = (log_lower_bound - mu) / sigma
    b = (log_upper_bound - mu) / sigma

    # Sample from the truncated normal distribution
    truncated_normal_samples = truncnorm.rvs(a, b, loc=mu, scale=sigma, size=Nsamples, random_state=rng_state)

    # Exponentiate the samples to get the log-normal samples
    truncated_lognormal_samples = np.exp(truncated_normal_samples)
    
    return truncated_lognormal_samples

# ---------------------------------------------------------------------------- #
