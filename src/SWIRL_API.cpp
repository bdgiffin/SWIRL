// ======================================================================== //
// C API functions for interacting with the SWIRL library        //
// ======================================================================== //

#include "SWIRL.h"
#include "Logger.h"
#include <iostream>
#include <stdio.h>

// global instance of the particle dynamics object
SWIRL SWIRL_instance;

// ======================================================================== //

// Define all C API functions within the following block:
extern "C" {

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

  // Retrieve current particle field data at the current analysis time
  // {ux,uy,uz}: arrays of vector components of particle displacements at the current time
  // {vx,vy,vz}: arrays of vector components of particle velocities at the current time
  // {fx,fy,fz}: arrays of vector components of particle forces at the current time
  void get_particle_field_data(double *ux, double *uy, double *uz,
	 	               double *vx, double *vy, double *vz,
		               double *fx, double *fy, double *fz) {
    // conditionally initialize the simulation state
    if (!SWIRL_instance.initialized()) SWIRL_instance.initialize();
    SWIRL_instance.debris.get_field_data(ux,uy,uz,vx,vy,vz,fx,fy,fz);
  } // get_particle_field_data()
  
  // ------------------------------------------------------------------------ //

  // Retrieve current wind field data at the current analysis time
  //    Npoints: number of sampling points at which to measure wind velocity and density
  // { x, y, z}: arrays of vector components of sampling point positions
  // {vx,vy,vz}: arrays of vector components of velocities at sampling points
  //       rhof: array of fluid densities at sampling points
  void get_wind_field_data(size_t Npoints, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *rhof) {
    // conditionally initialize the simulation state
    if (!SWIRL_instance.initialized()) SWIRL_instance.initialize();
    SWIRL_instance.wind_model->get_fluid_velocity_and_density(Npoints,SWIRL_instance.time,x,y,z,vx,vy,vz,rhof);
  } // get_particle_field_data()
  
// ======================================================================== //

  // The main API function called by OpenSees to initialize the external module
  void OPS_InitializeLineLoad(void) {

    // check to see if this is the first time this function is being called
    if (!SWIRL_instance.initialized()) {
      // initialize the particle dynamics object and define randomized particle positions
      SWIRL_instance.initialize();
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
    SWIRL_instance.update_state(time);

    // get loads applied to the requested element whose tag is specified
    SWIRL_instance.members.get_applied_forces(element_tag,coordinates,forces);
    DEBUG(std::cout << "Applying line load: element " << element_tag << ", time " << time << ", coords " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " " << coordinates[3] << " " << coordinates[4] << " " << coordinates[5] << ", forces " << forces[0] << " " << forces[1] << " " << forces[2] << " " << forces[3] << " " << forces[4] << " " << forces[5] << std::endl;)
    
  } // OPS_ApplyLineLoad
  
} // extern "C"

// ======================================================================== //
