#Analysis
# CONSTRAINTS handler(http://opensees.berkeley.edu/OpenSees/manuals/usermanual/617.htm)
# DOF NUMBERER  (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/366.htm)
# SYSTEM (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/371.htm)
# Convergence TEST (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/360.htm)
# Solution ALGORITHM (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/682.htm)
# Static INTEGRATOR (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/689.htm)
#  or use Transient INTEGRATOR: 
# ANALYSIS (http://opensees.berkeley.edu/OpenSees/manuals/usermanual/324.htm)

import openseespy.opensees as op
import math

def run_analysis(n_joints): 
    op.timeSeries("Linear", 1)
    op.pattern("Plain", 1, 1)

    # set damping based on first eigen mode
    #freq = eigen('-fullGenLapack', 1)[0]**0.5
    #dampRatio = 0.02
    #rayleigh(0., 0., 0., 2*dampRatio/freq)

    #New analysis 
    Constrant_Type = "Transformation" 
    numberer_Type = "RCM"
    system_type = "BandGeneral"
    algorith_type = "ModifiedNewton"
    Integrator_type = "Newmark"
    N_Gamma = 0.5
    N_Beta = 0.25
    analysis_type = "Transient"
    Tol = n_joints *1*math.exp(-5)
    maxNumIter = 10
    testTypeDynamic = "NormDispIncr"

    op.wipeAnalysis()			       # clear previously-define analysis parameters
    op.constraints(Constrant_Type)    # how it handles boundary conditions
    op.numberer(numberer_Type)        # renumber dof's to minimize band-width (optimization), if you want to
    op.system(system_type)            # how to store and solve the system of equations in the analysis
    op.test(testTypeDynamic, Tol, maxNumIter) #determine if convergence has been achieved at the end of an iteration step
    op.algorithm(algorith_type)	   # use Linear algorithm for linear analysis
    op.integrator(Integrator_type,N_Gamma,N_Beta)   # determine the next time step for an analysis
    op.analysis(analysis_type)        # define type of analysis: time-dependent

  