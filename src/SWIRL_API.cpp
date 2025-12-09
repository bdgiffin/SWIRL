// ======================================================================== //
// C API functions for interacting with the SWIRL library                   //
// ======================================================================== //

#include "SWIRL.h"
#include "Logger.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <unordered_map>

// global parameter list
std::unordered_map<std::string,double> params;

// global instance of the particle dynamics object
SWIRL SWIRL_instance;

// ======================================================================== //

// Define all C API functions within the following block:
extern "C" {
  
  // ------------------------------------------------------------------------ //

  // Define a named parameter and its corresponding value
  //   name: the string-valued name of the prameter being defined
  //  value: the real-valued parameter constant being defined
  void define_parameter(const char* name, double value) {
    params[name] = value;
  } // define_parameter()
  
  // ------------------------------------------------------------------------ //

  // Define the parameterized wind field within the simulation
  //       type: the string-valued name of the chosen wind field model "type" to instantiate
  // parameters: the condensed string containing the list of all wind field model parameters
  void define_wind_field(const char* type, double* parameters) {
    SWIRL_instance.wind_model = new_WindField(type,parameters);
  } // define_wind_field()
  
  // ------------------------------------------------------------------------ //

  // Define all particles within the simulation
  // n_particles: the total number of compact (spherical) particles to define
  // m: array of particle masses
  // d: array of particle diameters
  // {x,y,z}: arrays of particle initial positions in 3D space
  void define_particles(size_t n_particles, double *m, double *d, double *x, double *y, double *z) {
    SWIRL_instance.debris.define_particles(n_particles,m,d,x,y,z);
  } // define_particles()
  
  // ------------------------------------------------------------------------ //

  // Define all particles within the simulation
  // set_particles_active: boolean flag 
  void set_particles_active(bool particles_active) {
    SWIRL_instance.debris.particles_active = particles_active;
  } // set_particles_active()
  
  // ------------------------------------------------------------------------ //

  // Define a wedge-shaped control volume with periodic inflow/outflow conditions
  // {x0,y0}:   in-plane coordinates of the vertex of the wedge (should coincide with the center of the vortex)
  // theta_in:  the angular coordinate (measured in radians relative to the x-axis) at which the inflow  surface is defined
  // dtheta:    the positive angular dimension (measured in radians relative to theta_in) at which the outflow surface is defined
  void define_particle_control_volume(double x0, double y0, double theta_in, double dtheta) {
    SWIRL_instance.debris.define_control_volume(x0,y0,theta_in,dtheta);
  } // define_particle_control_volume()
  
  // ------------------------------------------------------------------------ //

  // (Legacy) method to define all structural members within the simulation
  // n_particles: the total number of compact (spherical) particles to define
  // m: array of particle masses
  // d: array of particle diameters
  // {x,y,z}: arrays of particle initial positions in 3D space
  void define_members(size_t n_members, size_t n_nodes_per_member, int *connectivity,
		      size_t n_joints, double *x, double *y, double *z) {
    SWIRL_instance.members.define_members(n_members,n_nodes_per_member,connectivity,n_joints,x,y,z);
  } // define_members()
  
  // ------------------------------------------------------------------------ //

  // (Legacy) API method to update the simulation state to the indicated analysis time
  // time: indicated analysis time to which the analysis should be updated
  void update_state(double time) {
    SWIRL_instance.update_state(time,params);
  } // update_state()
  
  // ------------------------------------------------------------------------ //

  // Retrieve current particle field data at the current analysis time
  // {ux,uy,uz}: arrays of vector components of particle displacements at the current time
  // {vx,vy,vz}: arrays of vector components of particle velocities at the current time
  // {fx,fy,fz}: arrays of vector components of particle forces at the current time
  void get_particle_field_data(double *ux, double *uy, double *uz,
	 	               double *vx, double *vy, double *vz,
		               double *fx, double *fy, double *fz) {
    // conditionally initialize the simulation state
    if (!SWIRL_instance.initialized()) SWIRL_instance.initialize(params);
    SWIRL_instance.debris.get_field_data(ux,uy,uz,vx,vy,vz,fx,fy,fz);
  } // get_particle_field_data()
  
  // ------------------------------------------------------------------------ //

  // Retrieve current particle positions at the current analysis time
  // {x,y,z}: arrays of vector components of particle positions at the current time
  void get_particle_positions(double *x, double *y, double *z) {
    // conditionally initialize the simulation state
    if (!SWIRL_instance.initialized()) SWIRL_instance.initialize(params);
    SWIRL_instance.debris.get_positions(x,y,z);
  } // get_particle_positions()
  
  // ------------------------------------------------------------------------ //

  // Retrieve current wind field data at the current analysis time
  //    Npoints: number of sampling points at which to measure wind velocity and density
  // { x, y, z}: arrays of vector components of sampling point positions
  // {vx,vy,vz}: arrays of vector components of velocities at sampling points
  //       rhof: array of fluid densities at sampling points
  void get_wind_field_data(size_t Npoints, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *rhof) {
    // conditionally initialize the simulation state
    if (!SWIRL_instance.initialized()) SWIRL_instance.initialize(params);
    SWIRL_instance.wind_model->get_fluid_velocity_and_density(Npoints,SWIRL_instance.time,x,y,z,vx,vy,vz,rhof);
  } // get_particle_field_data()
  
  // ------------------------------------------------------------------------ //

  // Retrieve metrics concerning discrete impact events
  //    Nimpacts: total number of discrete impact events
  //   max_force: maximum impact force (among all discrete impact events)
  // max_impulse: maximum impact impulse (among all discrete impact events)
  void get_impact_event_metrics(int* Nimpacts, double* max_force, double* max_impulse) {
    SWIRL_instance.members.impacts.get_metrics(Nimpacts,max_force,max_impulse);
  } // get_impact_event_metrics()
  
  // ------------------------------------------------------------------------ //

  // Finalize SWIRL instance to prepare for possible re-initialization
  void finalize(void) {
    SWIRL_instance = SWIRL();
    params = std::unordered_map<std::string,double>();
  } // finalize()
  
// ======================================================================== //

  // The main API function called by OpenSees to initialize the external module
  void OPS_InitializeLineLoad(void) {

    // check to see if this is the first time this function is being called
    if (!SWIRL_instance.initialized()) {
      // initialize the particle dynamics object and define randomized particle positions
      SWIRL_instance.initialize(params);
    }
    
  } // OPS_InitializeLineLoad
  
  // ------------------------------------------------------------------------ //

  // The API function called by OpenSees to initialize a new LineLoad element
  void OPS_DefineLineLoadSegment(int element_tag, double radius, const double* coordinates) {
    // (input)  element_tag:      the unique "tag" identifier for the current element
    // (input)  radius:           the effective radius of the current element
    // (input)  coordinates[2*3]: the nodal coordinates of the current element

    // define a new member in the Structure (if it doesn't already exist)
    DEBUG(std::cout << "Defining new line load: element " << element_tag << ", radius " << radius << ", coords " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " " << coordinates[3] << " " << coordinates[4] << " " << coordinates[5] << std::endl;)
    SWIRL_instance.members.define_member(coordinates, element_tag, radius);
    
  } // OPS_DefineLineLoadSegment
  
  // ------------------------------------------------------------------------ //

  // The main API function called by OpenSees to apply loads to the current LineLoad element at a requested analysis time
  void OPS_ApplyLineLoad(double time, int element_tag, const double* coordinates, double* forces) {
    // (input)  time:             the current analysis time
    // (input)  element_tag:      the unique "tag" identifier for the current element
    // (input)  coordinates[2*3]: the updated nodal coordinates of the current element
    // (output) forces[2*3]:      the forces applied to the nodes of the current element

    // conditionally update the simulation state to the indicated analysis time
    SWIRL_instance.update_state(time,params);

    // get loads applied to the requested element whose tag is specified
    SWIRL_instance.members.get_applied_forces(element_tag,coordinates,forces);
    DEBUG(std::cout << "Applying line load: element " << element_tag << ", time " << time << ", coords " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " " << coordinates[3] << " " << coordinates[4] << " " << coordinates[5] << ", forces " << forces[0] << " " << forces[1] << " " << forces[2] << " " << forces[3] << " " << forces[4] << " " << forces[5] << std::endl;)
    
  } // OPS_ApplyLineLoad
  
} // extern "C"

// ======================================================================== //
