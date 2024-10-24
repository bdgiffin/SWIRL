# Module for calling the C/C++ SWIRL API functions from Python

# Python package for calling C/C++ functions from Python
from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double, c_int, c_char_p

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

# ---------------------------------------------------------------------------- #

# Module initialization:

# Load the pre-compiled external C/C++ "shared object" libraries
library_name = "./libSWIRL.so"
API = CDLL(library_name)

# Define types to convert Numpy arrays into C arrays:

# C-type corresponding to 1D numpy array
ND_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.float64, 
                                      ndim=1,
                                      flags="C")
NI_POINTER_2 = np.ctypeslib.ndpointer(dtype=np.int32, 
                                      ndim=2,
                                      flags="C")

# Define all C/C++ library API function signatures
API.define_wind_field.argtypes = [c_char_p, ND_POINTER_1]
API.define_wind_field.restype  = None
API.define_particles.argtypes = [c_size_t, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.define_particles.restype  = None
API.define_members.argtypes = [c_size_t, c_size_t, NI_POINTER_2, c_size_t, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.define_members.restype  = None
API.get_particle_field_data.argtypes = [ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.get_particle_field_data.restype  = None
API.get_wind_field_data.argtypes = [c_size_t, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.get_wind_field_data.restype  = None

# ---------------------------------------------------------------------------- #

# Generate a randomized collection of particles with variable diameters and positions,
# and initialize the SWIRL object with these particles prior to initialization
def create_random_particles(n_particles,density,min_diameter,diameter_range,cylinder_radius,cylinder_height,cylinder_center,random_seed):
    # n_particles     [int]    The total number of spherical debris particles to be defined
    # density         [kg/m^3] The constant mass density assigned to all particles
    # min_diameter    [m]      The minimum particle diameter
    # diameter_range  [m]      The range of random particle diameters
    # cylinder_radius [m]      The diameter of the cylinder in which particles will be randomly distributed
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

# ---------------------------------------------------------------------------- #

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

# Close the connection to the SWIRL library and any Exodus files
def finalize():

    # Close the Exodus file
    if (num_particles > 0):
        exo.close()

# ---------------------------------------------------------------------------- #

