# =============================================================================
# Module for OpenSees
# =============================================================================
# NOTE: THE LOCALLY MODIFIED VERSION OF OPENSEES WITH THE LINELOAD ELEMENT MUST BE USED
import os
import sys
import math
import numpy as np
import time as timer
from datetime import timedelta
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import pyexodus  # https://pypi.org/project/pyexodus/ PYEXODUS V0.1.5 WITH PYTHON 3.12


# NOTE:  NEEDS TO BE MODIFIED TO WORK CORRECTLY 
import openseespy.opensees as op
import vfo.vfo as vfo
import opsvis as opsv


from Modules.SiUnits import * #Importing SI units
# Append the location of the locally installed SWIRL package to sys.path
sys.path.append("../../install/package/")
import SWIRL

# ============================================================================
# Command-Line for providing exo file
# =============================================================================
parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="filename",
                    help="input Exodus model file for the frame structure", metavar="FILE")
parser.add_argument("-q", "--quiet",
                    action="store_false", dest="verbose", default=True,
                    help="don't print status messages to stdout")
args = parser.parse_args()

# # ---------------------------------------------------------------------------- #

# check to make sure that the user has specified an Exodus file to define the geometry of the structure
if args.filename.lower().endswith('.exo'):
    file_format = "exo"
    if args.filename is None:
        print("ERROR: a valid Exodus model file must be specified to define the frame structure's geometry.")
        quit()

    # read data from the Exodus model file for the frame structure
    exoin = pyexodus.exodus(file=args.filename, mode='r', array_type='numpy', title=None, numDims=None, numNodes=None, numElems=None, numBlocks=None, numNodeSets=None, numSideSets=None, io_size=0, compression=None)
    x_in,y_in,z_in = exoin.get_coords() # [m] (assumed units/dimensions of structure expressed in meters)
    n_joints = len(x_in)
    connect_in,n_members,n_nodes_per_member = exoin.get_elem_connectivity(id=1)
    exoin.close()
    #print(connect_in[:5])
    
    # Initialize layer_in to hold the array of connect_in
    layer_in = np.empty((1,), dtype=object)  # Create an array to hold lists or arrays
    
    # Assign connect_in to the first element of layer_in
    layer_in[0] = connect_in
    
    # Print the lengths
    print(len(layer_in), len(layer_in[0]))

# =============================================================================
# Command-Line for reading CSV file for layer use #if using exo file then comment this
# =============================================================================
 # read any input arguments the user may have provided
elif args.filename.lower().endswith('.csv'):
    file_format = "csv"
    def read_points_and_connectivity_from_txt(filename):
        points = []
        connectivity = []
        layer = []  # Will store connectivity for each layer as a list of lists
        current_layer = None  # Track the current layer
        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
                reading_points = False
                reading_connectivity = False

                for line in lines:
                    line = line.strip()
                    if line == "Points:":
                        reading_points = True
                        reading_connectivity = False
                    elif line == "Connectivity:":
                        reading_points = False
                        reading_connectivity = True
                    elif reading_points:
                        coords = line.split(',')
                        if len(coords) == 3:
                            x, y, z = map(float, coords)
                            points.append((x, y, z))
                    elif reading_connectivity:
                        if line.startswith("Layer name:"):
                            layer_name = line.split(":")[1].strip()
                            layer.append([])
                            current_layer = len(layer) - 1  # Index of the current layer
                        else:
                            element_nodes = list(map(int, line.split()))
                            element_nodes = [node_id for node_id in element_nodes]  #Change this for TT because it sarts with 1 
                            connectivity.append(element_nodes)
                            # Append element nodes to the current layer in connectivity2
                            layer[current_layer].append(element_nodes)
        except Exception as e:
            print(f"Error reading from {filename}: {e}")
        return points, connectivity, layer

    txt_filename = args.filename  #For Transmission tower change line 172 element_nodes = [node_id for node_id in element_nodes]]
    points, connect_in, layer_in = read_points_and_connectivity_from_txt(txt_filename)



    n_members = len(connect_in)
    n_nodes_per_member = 2
    x_in = [round(row[0]*m,3) for row in points]
    y_in = [round(row[1]*m ,3 )for row in points]
    z_in = [round(row[2]*m, 3 ) for row in points]
    n_joints = len(x_in)

else:
    print("Provide a vaid file either exo or csv")
    exit()

supports = []
for i in range(0,n_joints):
    if (abs(z_in[i]) < 1.0e-6):
        supports.append(i)


# =============================================================================
# Command-Line for Defining properties of fluid and particles
# =============================================================================
# Call C/C++ library API functions from Python:
# Define randomized spherical particle parameters
n_particles = 0
particle_density         =  0.5*kg/pow(m,3) # [kg/m^3] (roughly the density of wood)
particle_min_diameter    = 0.01*m # [m]
particle_diameter_range  =  1.0*m # [m]
particle_cylinder_radius = 50.0*m # [m]
particle_cylinder_height = 25.25*m # [m]
particle_cylinder_center = [0.0*m,0.0*m,0.0*m] # [m,m,m]
random_seed = 1
SWIRL.create_random_particles(n_particles,particle_density,particle_min_diameter,particle_diameter_range,particle_cylinder_radius,particle_cylinder_height,particle_cylinder_center,random_seed)


#------------|----------------|
# EF Rating  | | 3S Gust(mph) | 
#------------|----------------|
# 0          | 65-85          | 
# 1          | 86-110         |
# 2          | 111-135        |
# 3          | 136-165        | 
# 4          | 166-200        |  
# 5          | Over 200       |    
#------------|----------------|

#------------|--------------|------|-------------------
# Parameters | Distribution | Mean | Standard Deviation
#------------|--------------|------|-------------------
# Vm         | Normal       | 40   | 10
# Um         | Normal       | 20   | 2
# Qm         | Normal       | 10   | 2
# rm         | Normal       | 100  | 20
# Y0         | Uniform      |      | 0 to 500
#------------|--------------|------|-------------------
#                  Debris     
#------------|--------------|------|-------------------
# Parameters | Distribution | Mean | Standard Deviation
#------------|--------------|------|-------------------
# Mass(kg)   | Normal       |0.2   |0.05
# Area(m^2)  | Normal       |0.1   |0.02
#------------|--------------|------|-------------------


# Create the parameterized wind field model (Baker Sterling Vortex)
wind_field_params = np.zeros(12)
wind_field_params[0]  = 60*m/sec # [m/s]      Um: reference radial velocity
Vm = 10*m/sec                                 #Vm : maximum Cirumfrential velocity 
wind_field_params[1]  = 100.0*m   # [m]        rm: reference radius
wind_field_params[2]  = 30.0*m  # [m]        zm: reference height
wind_field_params[3]  =  2  #             S: swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
wind_field_params[4]  = 2.0   #         gamma: 
wind_field_params[5]  = 1.293*kg/pow(m,3) # [kg/m^3] rho0: reference density of air at STP
wind_field_params[6]  = particle_cylinder_center[0]  # [m]       xc0: x-position of the vortex center
wind_field_params[7]  = particle_cylinder_center[1]  # [m]       yc0: y-position of the vortex center
wind_field_params[8]  = particle_cylinder_center[2] # [m]       zc0: z-position of the vortex center
wind_field_params[9]  = 0.0*m/sec   # [m/s]     vxc: x-velocity of the vortex center
wind_field_params[10] = 0.0*m/sec   # [m/s]     vyc: y-velocity of the vortex center
wind_field_params[11] = 0.0*m/sec   # [m/s]     vzc: z-velocity of the vortex center
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
#ParticleDynamics.API.define_wind_field(b"RankineVortex",wind_field_params)

# =============================================================================
# Start of OpenSees model generation ------------------------------------------
# =============================================================================
op.wipe()	
op.wipeAnalysis()
# clear opensees model
op.model('basic', '-ndm', 3, '-ndf', 6)	       # 3 dimensions, 3 dof per node
print(n_joints)
for i in range(0,n_joints):
    op.node(i+1, x_in[i], y_in[i], z_in[i]) # [m] (node, X, Y, Z)

for isupport in supports:
    op.fix(isupport+1, 1, 1, 1, 1, 1, 1) # node DX DY DZ RX RY RZ
    
# =============================================================================
# define MATERIAL -------------------------------------------------------------
# =============================================================================
Fy = 345 * Mpa
Es = 200* Gpa  # Steel Young's Modulus
nu = 0.3
Gs = Es / (2 * (1 + nu))  # Torsional stiffness Modulus
J = 10  # Large torsional stiffness
Bs = 0.01
R0 = 20
cR1 = 0.925
cR2 = 0.15
a1 = 0.39
a2 = 1.0
a3 = 0.029
a4 = 1.0
matIDhard = 1
matType = 'Steel02'
# Function to define uniaxial material in Python
op.uniaxialMaterial(matType, matIDhard, Fy, Es, Bs, R0, cR1, cR2)

# =============================================================================
#Transformation absed on the orientation
# =============================================================================
Trans_Type = "Corotational" 
Transf = [1,2,3,4,5,6,7]

op.geomTransf(Trans_Type, Transf[0], 0, 0, 1)     #main vertical members
op.geomTransf(Trans_Type, Transf[1], 0, 0, 1)  
op.geomTransf(Trans_Type, Transf[2], 0, 0, 1)  
op.geomTransf(Trans_Type, Transf[3], 0, 1, 1)  
op.geomTransf(Trans_Type, Transf[4], 0, 0, 1)  
op.geomTransf(Trans_Type, Transf[5], 0, 0, 1)  #Cross Members
op.geomTransf(Trans_Type, Transf[6], 0, 0, 1)  #Cross Members
# =============================================================================
# Elements size assign
# =============================================================================
secTag = [1,2,3,4,5,6,7]
BreID = 2
Lfiber = 20
Sfiber = 3
Thick = np.array([12*mm,12*mm,10*mm,6*mm,5*mm,5*mm,5*mm]) #should be based on number of layers
Length = np.array([130*mm,110*mm,100*mm,60*mm,45*mm,50*mm,50*mm]) #should be based on number of layers
Ly1= -Thick/2
Hy1= -Thick/2
Ly2= Length-Thick/2
Hy2= Thick/2    

# =============================================================================
# Ql Section ------------------------------------------------------------------
# =============================================================================
QLsection = 310*N/pow(m,3)
steel_mass_density    = 7850.0*kg/(m**3)  # [kg/m^3] (mass density of steel)
shear_modulus        = 0.5*Es /(1.0+nu)
radius_of_gyration    = 0.05*m     # (effective radius 2in = 0.05m)
cross_sectional_area = np.zeros(len(layer_in))
mass_per_unit_length =np.zeros(len(layer_in))
polar_moment_of_area = np.zeros(len(layer_in))
moment_of_area_x = np.zeros(len(layer_in))
moment_of_area_y =np.zeros(len(layer_in))
lumped_mass = np.zeros(n_joints)

for i in range(len(layer_in)):
   cross_sectional_area[i]  = Length[i]*Thick[i]+(Length[i]-Thick[i])*Thick[i]  # (cross-sectional area 1in^2 = 0.00065m^2)
   mass_per_unit_length[i] = steel_mass_density*cross_sectional_area[i]
   polar_moment_of_area[i] = cross_sectional_area[i]*radius_of_gyration*radius_of_gyration
   moment_of_area_x[i]     = 0.5*polar_moment_of_area[i]
   moment_of_area_y[i]     = 0.5*polar_moment_of_area[i]


# Elastic Section Use only when not using Fiber Section
# ElasTag = 1
# for i in range(len(layer_in)):
#     section('Elastic', ElasTag, Es, cross_sectional_area[i], moment_of_area_x[i], moment_of_area_y[i], shear_modulus, polar_moment_of_area[i])


# =============================================================================
#lump the nodal masses to the joints of the structure -------------------------
# =============================================================================
for i in range(len(layer_in)):
        for j in range(len(layer_in[i])):
        # for j in range(0,n_members):
            i1 = layer_in[i][j][0]-1
            i2 = layer_in[i][j][1]-1
            # i1 = connect_in[i][0]-1
            # i2 = connect_in[i][1]-1 
            dx = x_in[i2] - x_in[i1]
            dy = y_in[i2] - y_in[i1]
            dz = z_in[i2] - z_in[i1]
            member_length = math.sqrt(dx*dx + dy*dy + dz*dz)
            member_mass   = steel_mass_density*cross_sectional_area[i]*member_length
            lumped_mass[i1] = lumped_mass[i1] + 0.5*member_mass
            lumped_mass[i2] = lumped_mass[i2] + 0.5*member_mass


# nodal masses: (only needed if masses are not already computed by the element)
op.timeSeries("Linear", 1)
op.pattern("Plain", 1, 1)
for i in range(0,n_joints):
     op.mass(i+1, lumped_mass[i], lumped_mass[i], 0.0) # [kg] node#, Mx My Mz, Mass=Weight/g.
     op.load(i+1, 0.0, 0.0, -lumped_mass[i]*g, 0.0, 0.0 ,0.0)
     
# =============================================================================
# Defining Fiber Section
# =============================================================================
from Modules.Fiber import *  #The fiber module is present in Module file 
for i in range(len(layer_in)):
    FiberCreation(secTag[i],matIDhard,Sfiber,Lfiber,Ly1[i],Hy1[i],Ly2[i],Hy2[i])

# =============================================================================
# Defining Element Section
# =============================================================================
Radius = 0.5 # This is supposed radius for the member to calculate wind and debris forces  
#Mass_Den = 
k = 0
for i in range(len(layer_in)):
        for j in range(len(layer_in[i])):
        #  print(k,layer_in[i][j][0], layer_in[i][j][1], secTag[i], Transf[i])
        # for j in range(0,n_members):
         #print(f"Index {i+1}: {connectivity[i]}")  
        #  element('elasticBeamColumn', j+1, int(connect_in[j][0]), int(connect_in[j][1]), ElasTag, Transf[i], '-mass', mass_per_unit_length[i], '-cMass')
        #  element("LineLoad", j+1+n_members, int(connect_in[j][0]), int(connect_in[j][1]), radius_of_gyration, ParticleDynamics.library_name) # [m] (LineLoad, LineLoadID, node1, node2, radius, library)
        #  element('elasticBeamColumn', k+1, int(connect_in[j][0]), int(connect_in[j][1]), secTag[i], Transf[i])
        #  element("LineLoad", k+1+n_members, int(connect_in[j][0]), int(connect_in[j][1]), radius_of_gyration, ParticleDynamics.library_name) # [m] (LineLoad, LineLoadID, node1, node2, radius, library)
         if file_format =="exo":
            #  op.element('elasticBeamColumn', k+1, int(layer_in[i][j][0]), int(layer_in[i][j][1]), secTag[i], Transf[i])
             op.element('nonlinearBeamColumn', k+1, int(layer_in[i][j][0]), int(layer_in[i][j][1]), 10, secTag[i], Transf[i])
             op.element("LineLoad", k+1+n_members, int(layer_in[i][j][0]), int(layer_in[i][j][1]), radius_of_gyration, SWIRL.library_name) # [m] (LineLoad, LineLoadID, node1, node2, radius, library)
         
         elif file_format =="csv":
            #  op.element('elasticBeamColumn', k+1, layer_in[i][j][0], layer_in[i][j][1], secTag[i], Transf[i])
             op.element('nonlinearBeamColumn', k+1, layer_in[i][j][0], layer_in[i][j][1], 10, secTag[i], Transf[i])
             op.element("LineLoad", k+1+n_members, layer_in[i][j][0], layer_in[i][j][1], radius_of_gyration, SWIRL.library_name) # [m] (LineLoad, LineLoadID, node1, node2, radius, library)
         k+=1
        

# =============================================================================
# RECORDER -------------------------------------------------------------
# =============================================================================
# Ensure dataDir is a valid directory path
dataDir = "Main_Dynamic"
if not os.path.exists(dataDir):
    os.makedirs(dataDir)
    print(f"Directory '{dataDir}' created.")

# array to store node ID and z_in value
push = np.array([[i + 1, z_in[i], x_in[i], y_in[i]] for i in range(n_joints)])
# Sort the array by (z_in values)
condu = push[push[:, 3].argsort()]
for i in range(-3, 3):
    op.mass(75*kg, 75*kg, 75*kg)
    op.load(int(condu[i][0]), 0.0, 0.0, -75*kg*g, 0.0, 0.0 ,0.0)

push = push[push[:, 1].argsort()]
Height = push[-1][1]-push[0][1]
print("Height of the tower is ",Height)
# =============================================================================
# For 3-d Visualization
# =============================================================================

# x = [0] * (len(layer_in) + 1)
# x[0] = 1
# for i in range(len(layer_in)):
#     x[i+1] = x[i] + len(layer_in[i])


# element_ranges = [list(range(x[i], x[i+1])) for i in range(len(layer_in))]
# colors = ["red", "blue", "green", "yellow", "cyan", "magenta", "orange"]


# vfo.plot_model(
#     elementgroups=[element_ranges, colors[:len(element_ranges)]],
#     show_nodes='yes',
#     show_nodetags='yes',
#     show_eletags='no',
#     font_size=15,
#     setview='3D',
#     line_width=3
# )
# exit()

# =============================================================================
# # # Static analysis (To initilize wind load) -------------------------------
# =============================================================================
Static_step = 1
D_Gravity = 1/Static_step

op.constraints('Plain')    
op.numberer('RCM')       
op.system('BandGeneral')      
op.test('NormDispIncr', 1.0e-6, 6)
op.algorithm('Newton')
op.integrator('LoadControl',D_Gravity)
op.analysis('Static')
op.analyze(Static_step)
op.loadConst('-time', 0.0)
print("Static Analysis Complete it initilizes wind load ")


# =============================================================================
# Recorder for max displacement at the top and base reaction
# =============================================================================
free_file = os.path.join(dataDir, "DFree.out")
fixed_file = os.path.join(dataDir, "DFixed.out")
react = os.path.join(dataDir, "RXN.out")

ISDR = [11,12,31,32,33,43,44,67,72,70,75,77,131]
for i in range(len(ISDR)):
    file_name = f"node_{ISDR[i]}_drift.out"  # Or customize the file name as needed
    file_path = os.path.join(dataDir, file_name)
    op.recorder("Node", '-file', file_path, 'time', '-node', int(ISDR[i]), '-precision', 3, '-time', '-dof', 1, 2, 'disp')

op.recorder("Node", '-file', free_file, 'time', '-node', int(push[-1][0]),'-precision',3, '-time' ,'-dof', 1, 'disp')
op.recorder("Node", '-file', fixed_file, 'time', '-node', int(push[0][0]),'-precision',3, '-time' ,'-dof', 1, 'disp')
op.recorder("Node", '-file', react, 'time', '-node', int(push[0][0]), '-precision',3,'-time', '-dof', 1, 'reaction')


output_directory = 'LineLoadTest_PVD'
# output position-velocity-displacement (PVD) data
op.recorder('PVD', output_directory, 'disp', 'reaction' ,'unbalancedLoad')

# create a matching directory in which to dump the output data
if not os.path.exists(output_directory):
    os.makedirs(output_directory)


# from Modules.Dynamic_Analysis import run_analysis
# # Pass the variables as arguments to the function in Dynamic_Analysis
# run_analysis(n_joints)

# =============================================================================
# # # DYNAMIC analysis -------------------------------------------------------------
# =============================================================================
op.wipeAnalysis()	 # clear previously-define analysis parameters
tCurrent = op.getTime()
time = tCurrent # [s] starting time
dt   = 0.1 # [s] time increment
nPts = 30
tFinal = nPts*dt
ok = 0

SWIRL.output_state(time)
op.timeSeries("Linear", 2)
op.pattern("Plain", 2, 2)

#New analysis 
Constrant_Type = "Transformation" 
numberer_Type = "RCM"
system_type = "BandGeneral"
Integrator_type = "Newmark"
N_Gamma = 0.5
N_Beta = 0.25
analysis_type = "Transient"
Tol = math.exp(-5)
maxNumIter = 10
test = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 3:'EnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 3:'ModifiedNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

		      
op.constraints(Constrant_Type)    # how it handles boundary conditions
op.numberer(numberer_Type)        # renumber dof's to minimize band-width (optimization), if you want to
op.system(system_type)            # how to store and solve the system of equations in the analysis
op.test(test[1], Tol, 1000, 0)
op.algorithm(algorithm[3])
op.integrator(Integrator_type,0.5,0.25)
op.analysis(analysis_type)

# while tCurrent < tFinal:
# #    ok = op.analyze(1, .01)     
#     for i in test:
#         for j in algorithm: 
#             if j < 4:
#                 op.algorithm(algorithm[j], '-initial')
                
#             else:
#                 op.algorithm(algorithm[j])
#             while ok == 0 and tCurrent < tFinal:
                    
#                 op.test(test[i], Tol, maxNumIter)        
#                 NewmarkGamma = 0.5
#                 NewmarkBeta = 0.25
#                 op.integrator(Integrator_type, NewmarkGamma, NewmarkBeta)
#                 op.analysis(analysis_type)
#                 ok = op.analyze(1, .01)
#                 SWIRL.output_state(time)
                
#                 if ok == 0 :
#                     tCurrent = op.getTime()                
#                     time.append(tCurrent)
#                     print(test[i], algorithm[j], 'tCurrent=', tCurrent)


for step_id in range(nPts):
    time = time + dt
    print(time)
    op.analyze(1,dt) # apply 1 time step of size dt in the opensees analysis
    SWIRL.API.update_state(time)
    SWIRL.output_state(time)
 
# finalize the ParticleDynamics module (close the Exodus files)
SWIRL.finalize()

# =============================================================================
# # # Plotting Displacement at the top point ----------------------------------
# =============================================================================
time_series1 = np.loadtxt(free_file)
Drift = (time_series1[:, 1]) # 3.0 is the height between node 1 and node 2
Time = (time_series1[:, 0])
print(Drift)
exit()
# plt.figure(figsize=(8,8))
# plt.plot(Time, Drift, color='blue', linewidth=2, label='TimeSeries')
# plt.xlabel('Time (s)')
# plt.ylabel("Displacement at the Top (m)", fontsize=14)
# plt.title("Time-History", fontsize=16)
# plt.show()
# =============================================================================
# # # Plotting ISDR ratio-----------------------------------------------------
# =============================================================================
# Formula for inter-Segmental Drift Ratio
# ISDR(i) = max[{u_i(t)- u_i-1(t)}/{h_i(t)- h_i-1(t)}-Thita_i-1(t)]


# output_file = os.path.join(dataDir, "drift_data_output.out")
# max_file =  os.path.join(dataDir, "drift_max.out")
# max_abs_drift_x = 0
# max_abs_drift_y = 0
# vx_values = []  # Placeholder for storing Vx values

# # Open the output file for writing drift data
# with open(output_file, 'w') as f_out:
#     # f_out.write("Time\tDrift_x\tDrift_y\n")  # Write headers

#     for i in range(len(ISDR) - 1):
#         # File for node i
#         file_name1 = f"node_{ISDR[i]}_drift.out" 
#         file_path1 = os.path.join(dataDir, file_name1)
#         coords1 = op.nodeCoord(ISDR[i])
#         z_value1 = coords1[2]  # Z-coordinate for node i
        
#         # File for node i+1
#         file_name2 = f"node_{ISDR[i+1]}_drift.out" 
#         file_path2 = os.path.join(dataDir, file_name2)
#         coords2 = op.nodeCoord(ISDR[i+1]) 
#         z_value2 = coords2[2]  # Z-coordinate for node i+1
        
#         # Load displacement data for both nodes
#         Disp1 = np.loadtxt(file_path1)
#         Disp2 = np.loadtxt(file_path2)
        
#         # Corrected drift calculation
#         Drift_x = (Disp1[:, 1] - Disp2[:, 1]) / (z_value1 - z_value2)  # X-direction drift
#         Drift_y = (Disp1[:, 2] - Disp2[:, 2]) / (z_value1 - z_value2)  # Y-direction drift
        
#         max_abs_drift_x = max(max_abs_drift_x, np.max(np.abs(Drift_x)))
#         max_abs_drift_y = max(max_abs_drift_y, np.max(np.abs(Drift_y)))
        
        
#         # Write time, Drift_x, and Drift_y for all time steps to file
#         for j in range(len(Disp1)):
#             time = Disp1[j, 0]
#             drift_x = Drift_x[j]
#             drift_y = Drift_y[j]
#             f_out.write(f"{time}\t{drift_x}\t{drift_y}\n")
        

#         print(f"Drift between node {ISDR[i]} and node {ISDR[i+1]} recorded.")

# print(f"Drift data has been saved to {output_file}")

# time_series1 = np.loadtxt(output_file)
# Driftx = (time_series1[:, 1]) # 3.0 is the height between node 1 and node 2
# Driftx = (time_series1[:, 2])
# Time = (time_series1[:, 0])


# # plt.figure(figsize=(8,8))
# # plt.plot(Time, Driftx, color='blue', linewidth=2, label='TimeSeries')
# # plt.xlabel('Time (s)')
# # plt.ylabel("Inter-Segmental Drift Ratio", fontsize=14)
# # plt.title("ISDR time history", fontsize=16)
# # plt.show()
# with open(max_file, 'a') as summary_out:
#     summary_out.write(f"{wind_field_params[0]}\t{max_abs_drift_x}\t{max_abs_drift_y} \n")
    
# print(wind_field_params[0])
# # import os
# # import numpy as np
# # import matplotlib.pyplot as plt
# # file = "Main_Dynamic/Result"
# # max_file = os.path.join(file, "near1000")


# data = np.loadtxt(max_file)
# Vm = data[:, 0]  # First column: Vm
# max_abs_drift_x = data[:, 1]  # Second column: max_abs_drift_x
# max_abs_drift_y = data[:, 2]  # Third column: max_abs_drift_y


# plt.figure(figsize=(8, 6))
# plt.plot(Vm, max_abs_drift_x, label='Max Abs Drift X', marker='o', linestyle='-', color='blue')
# plt.plot(Vm, max_abs_drift_y, label='Max Abs Drift Y', marker='s', linestyle='--', color='green')

# # Add labels, legend, and title
# plt.xlabel("Vm", fontsize=12)
# plt.ylabel("Max Absolute Drift", fontsize=12)
# plt.title("Vm vs Max Absolute Drift", fontsize=14)
# plt.legend()
# plt.grid(True)

# # Show the plot
# plt.show()