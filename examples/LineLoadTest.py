# Module for OpenSees
# NOTE: THE LOCALLY MODIFIED VERSION OF OPENSEES WITH THE LINELOAD ELEMENT MUST BE USED
from openseespy.opensees import *

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
import os
import math
import time as timer
from datetime import timedelta
import numpy as np
from argparse import ArgumentParser

# ---------------------------------------------------------------------------- #

# read any input arguments the user may have provided
parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="filename",
                    help="input Exodus model file for the frame structure", metavar="FILE")
parser.add_argument("-q", "--quiet",
                    action="store_false", dest="verbose", default=True,
                    help="don't print status messages to stdout")
args = parser.parse_args()

# ---------------------------------------------------------------------------- #

# check to make sure that the user has specified an Exodus file to define the geometry of the structure
if args.filename is None:
    print("ERROR: a valid Exodus model file must be specified to define the frame structure's geometry.")
    quit()

# read data from the Exodus model file for the frame structure
exoin = pyexodus.exodus(file=args.filename, mode='r', array_type='numpy', title=None, numDims=None, numNodes=None, numElems=None, numBlocks=None, numNodeSets=None, numSideSets=None, io_size=0, compression=None)
x_in,y_in,z_in = exoin.get_coords() # [m] (assumed units/dimensions of structure expressed in meters)
n_joints = len(x_in)
connect_in,n_members,n_nodes_per_member = exoin.get_elem_connectivity(id=1)
exoin.close()

# do some model pre-processing: identify the joints with supports
supports = []
for i in range(0,n_joints):
    if (abs(z_in[i]) < 1.0e-6):
        supports.append(i)

# define the cross-sectional area, effective member radius, and mass density assigned to all members
modulus_of_elasticity = 200.0e+9 # [kg*m/s^2] modulus of elasticity of steel (200 GPa)
poissons_ratio        = 0.3      # dimensionless
steel_mass_density    = 7850.0   # [kg/m^3] (mass density of steel)
cross_sectional_area  = 0.00065  # (cross-sectional area 1in^2 = 0.00065m^2)
radius_of_gyration    = 0.05     # (effective radius 2in = 0.05m)

shear_modulus        = 0.5*modulus_of_elasticity/(1.0+poissons_ratio)
mass_per_unit_length = steel_mass_density*cross_sectional_area
polar_moment_of_area = cross_sectional_area*radius_of_gyration*radius_of_gyration
moment_of_area_x     = 0.5*polar_moment_of_area
moment_of_area_y     = 0.5*polar_moment_of_area
        
# lump the nodal masses to the joints of the structure
#lumped_mass = np.zeros(n_joints)
#for i in range(0,n_members):
#    i1 = connect_in[i][0]-1
#    i2 = connect_in[i][1]-1
#    dx = x_in[i2] - x_in[i1]
#    dy = y_in[i2] - y_in[i1]
#    dz = z_in[i2] - z_in[i1]
#    member_length = math.sqrt(dx*dx + dy*dy + dz*dz)
#    member_mass   = steel_mass_density*cross_sectional_area*member_length
#    lumped_mass[i1] = lumped_mass[i1] + 0.5*member_mass
#    lumped_mass[i2] = lumped_mass[i2] + 0.5*member_mass

# ---------------------------------------------------------------------------- #

# Call C/C++ library API functions from Python:

# Define randomized spherical particle parameters
n_particles = 200
particle_density         =  0.5 # [kg/m^3] (roughly the density of wood)
particle_min_diameter    = 0.01 # [m]
particle_diameter_range  =  1.0 # [m]
particle_cylinder_radius = 20.0 # [m]
particle_cylinder_height = 40.0 # [m]
particle_cylinder_center = [10.0,0.0,0.0] # [m,m,m]
random_seed = 1
SWIRL.create_random_particles(n_particles,particle_density,particle_min_diameter,particle_diameter_range,particle_cylinder_radius,particle_cylinder_height,particle_cylinder_center,random_seed)

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

# Create a Rankine vortex
#wind_field_params = np.zeros(12)
#wind_field_params[0]  = 100.0 # [m/s]      Um: reference tangential velocity
#wind_field_params[1]  = 10.0  # [m]        rm: reference outer radius
#wind_field_params[2]  = 1.0   # [m]        rc: reference core radius
#wind_field_params[3]  = 1.0   #             E: decay index
#wind_field_params[4]  = 0.0   #       (unused) 
#wind_field_params[5]  = 1.293 # [kg/m^3] rho0: reference density of air at STP
#wind_field_params[6]  = 10.0  # [m]       xc0: x-position of the vortex center
#wind_field_params[7]  = 0.0   # [m]       yc0: y-position of the vortex center
#wind_field_params[8]  = 0.0   # [m]       zc0: z-position of the vortex center
#wind_field_params[9]  = 0.0   # [m/s]     vxc: x-velocity of the vortex center
#wind_field_params[10] = 0.0   # [m/s]     vyc: y-velocity of the vortex center
#wind_field_params[11] = 0.0   # [m/s]     vzc: z-velocity of the vortex center
#SWIRL.API.define_wind_field(b"RankineVortex",wind_field_params)

# ----------------------------------
# Start of OpenSees model generation
# ----------------------------------

# --------------------------------------------------------------------------------------------------
# Example 3D linear dynamic truss structure subjected to time-varying vortex wind and debris loading
# all units are in Newtons, meters, seconds
#

# SET UP ----------------------------------------------------------------------------

wipe()				               # clear opensees model
model('basic', '-ndm', 3, '-ndf', 6)	       # 3 dimensions, 6 dof per node (3 displacements + 3 rotations)
# file mkdir data 			       # create data directory

# define GEOMETRY -------------------------------------------------------------

# nodal coordinates:
for i in range(0,n_joints):
    node(i+1, x_in[i], y_in[i], z_in[i]) # [m] (node, X, Y, Z)

# Single point constraints -- Boundary Conditions
for isupport in supports:
    fix(isupport+1, 1, 1, 1, 0, 0, 0) # node DX DY DZ RX RY RZ

# define MATERIAL -------------------------------------------------------------

# nodal masses: (only needed if masses are not already computed by the element)
#for i in range(0,n_joints):
#    mass(i+1, lumped_mass[i], lumped_mass[i], lumped_mass[i]) # [kg] node#, Mx My Mz, Mass=Weight/g.

# define materials
#matTag = 1
#uniaxialMaterial("Elastic", matTag, modulus_of_elasticity) # [kg*m/s^2] modulus of elasticity of steel (200 GPa)

# Define SECTION -------------------------------------------------------------

# define section
secTag = 1
section('Elastic', secTag, modulus_of_elasticity, cross_sectional_area, moment_of_area_x, moment_of_area_y, shear_modulus, polar_moment_of_area)

# define geometric transformation (linear, for now, but will eventually need to make this non-linear to capture buckling instabilities)
transfTag = 1
geomTransf('Linear', transfTag, 1.0, 1.0, 1.0)

# Define ELEMENTS -------------------------------------------------------------

# NOTE: The structure should be represented in terms of nonlinear beam-column elements, rather than simple truss elements

# define truss element connectivity
for i in range(0,n_members):
    element('elasticBeamColumn', i+1, int(connect_in[i][0]), int(connect_in[i][1]), secTag, transfTag, '-mass', mass_per_unit_length, '-cMass')
    #element("Truss", i+1, int(connect_in[i][0]), int(connect_in[i][1]), cross_sectional_area, matTag) # [m^2] (Truss, TrussID, node1, node2, area, material)

# define LineLoad elements
for i in range(0,n_members):
    element("LineLoad", i+1+n_members, int(connect_in[i][0]), int(connect_in[i][1]), radius_of_gyration, SWIRL.library_name) # [m] (LineLoad, LineLoadID, node1, node2, radius, library)

# RECORDER -------------------------------------------------------------

# output position-velocity-displacement (PVD) data
output_directory = 'LineLoadTest_PVD'
recorder('PVD', output_directory, 'disp', 'reaction' ,'unbalancedLoad')

# create a matching directory in which to dump the output data
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# DYNAMIC analysis -------------------------------------------------------------

# create TimeSeries
timeSeries("Linear", 1)

# create a plain load pattern
pattern("Plain", 1, 1)

# set damping based on first eigen mode
#freq = eigen('-fullGenLapack', 1)[0]**0.5
#dampRatio = 0.02
#rayleigh(0., 0., 0., 2*dampRatio/freq)

# create the analysis
#wipeAnalysis()			 # clear previously-define analysis parameters
constraints('Plain')    	 # how it handles boundary conditions
numberer("RCM")                  # renumber dof's to minimize band-width (optimization), if you want to
system('BandGeneral')            # how to store and solve the system of equations in the analysis
algorithm('Linear')	         # use Linear algorithm for linear analysis
integrator('Newmark', 0.5, 0.25) # determine the next time step for an analysis
analysis('Transient')            # define type of analysis: time-dependent

# RUN analysis -------------------------------------------------------------

# perform the analysis
time = 0.0 # [s] starting time
dt   = 0.01 # [s] time increment
SWIRL.output_state(time)
for step_id in range(1,100):
    time = time + dt
    analyze(1,dt) # apply 1 time step of size dt in the opensees analysis
    SWIRL.output_state(time)

# Request and report impact event metrics that occurred during the analysis
Nimpacts, max_force, max_impulse = SWIRL.impact_event_metrics()
print("  # discrete impact events  = " + str(Nimpacts))
print("max discrete impact force   = " + str(max_force))
print("max discrete impact impulse = " + str(max_impulse))

# Output all recorded impact events to a separate file for visualization (impact_events.exo)
SWIRL.output_impact_events()

# finalize the SWIRL module (close the Exodus files)
SWIRL.finalize()
