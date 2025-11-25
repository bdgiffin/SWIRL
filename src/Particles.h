#ifndef PARTICLES_H
#define PARTICLES_H

#include "WindField.h"
#include "Logger.h"
#include <vector>
#include <cmath>
#include <memory>

// ======================================================================== //


// A collection of particles
struct Particles {

  // ------------------- Declare public member functions ------------------ //

  // Default empty constructor
  Particles(void) { }
  
  // ---------------------------------------------------------------------- //

  // Method to define all compact particles at once
  void define_particles(size_t n_particles, double *m_in, double *d_in, double *x_in, double *y_in, double *z_in) {
    
    // Store incoming data for all newly defined particles
    num_particles = n_particles;
    mass.resize(num_particles);
    radius.resize(num_particles);
    x0.resize(num_particles);
    y0.resize(num_particles);
    z0.resize(num_particles);

    for (int i=0; i<num_particles; i++) {
      mass[i]   = m_in[i];
      radius[i] = 0.5*d_in[i];
      x0[i] = x_in[i];
      y0[i] = y_in[i];
      z0[i] = z_in[i];
    } // for i=1,...,num_particles
    
  } // define_particles()

  // ---------------------------------------------------------------------- //

  void define_control_volume(double x0_new, double y0_new, double theta_in_new, double dtheta_new) {
    cv_x0 = x0_new;
    cv_y0 = y0_new;
    cv_theta_in = theta_in_new;
    cv_dtheta   = dtheta_new;
    if (cv_dtheta < 0.0) std::cerr << "WARNING: CV dtheta should be positive!" << std::endl;
  } // define_control_volume()
  
  // ---------------------------------------------------------------------- //

  // Method to retrieve field data for all compact particles at the current time
  void get_field_data(double *ux_out, double *uy_out, double *uz_out,
	 	      double *vx_out, double *vy_out, double *vz_out,
		      double *fx_out, double *fy_out, double *fz_out) {

    for (int i=0; i<num_particles; i++) {
      ux_out[i] =  x[i] - x0[i];
      uy_out[i] =  y[i] - y0[i];
      uz_out[i] =  z[i] - z0[i];
      vx_out[i] = vx[i];
      vy_out[i] = vy[i];
      vz_out[i] = vz[i];
      fx_out[i] = fx[i];
      fy_out[i] = fy[i];
      fz_out[i] = fz[i];
    } // for i=1,...,num_particles
    
  } // get_field_data()

  // ---------------------------------------------------------------------- //

  // Method to retrieve the current positions for all compact particles at the current time
  void get_positions(double *x_out, double *y_out, double *z_out) {

    for (int i=0; i<num_particles; i++) {
      x_out[i] = x[i];
      y_out[i] = y[i];
      z_out[i] = z[i];
    } // for i=1,...,num_particles
    
  } // get_positions()

  // ---------------------------------------------------------------------- //

  // Initialize the Particles object (assuming all particles have been defined)
  void initialize(std::unique_ptr<WindField>& wind_model, double initial_time = 0.0, double new_drag_coeff = 0.47, double new_contact_stiff = 1.0e+4, double new_contact_damp_ratio = 0.5) {
    
    DEBUG(std::cout << "Initializing Particles object with " << num_particles << " particles defined" << std::endl;)

    // Set constants common to all particles
    drag_coeff         = new_drag_coeff;         // (default = 0.47 for an assumed spherical particle)
    contact_stiff      = new_contact_stiff;      // [N/m] = [kg/s^2]
    contact_damp_ratio = new_contact_damp_ratio; // (default = 0.5 for 50% of critical damping relative to the contact stiffness)

    // Initialize stored data for all particles
    ghost.resize(num_particles);
    contact_damp.resize(num_particles);
    x.resize(num_particles);
    y.resize(num_particles);
    z.resize(num_particles);
    vx.resize(num_particles);
    vy.resize(num_particles);
    vz.resize(num_particles);
    fx.resize(num_particles);
    fy.resize(num_particles);
    fz.resize(num_particles);

    // determine the fluid velocity and density at the initial positions of all particles
    std::vector<double> rhof(num_particles);
    wind_model->get_fluid_velocity_and_density(num_particles,initial_time,x0.data(),y0.data(),z0.data(),
					       vx.data(),vy.data(),vz.data(),rhof.data());
    
    // initialize the particle's position, and set its initial velocity as a specified fraction of the surrounding fluid velocity
    double initial_velocity_fraction = 0.8;
    for (int i=0; i<num_particles; i++) {
      ghost[i] = true; // make all particles ghost particles upon initialization
      double critical_damp = 2.0*std::sqrt(contact_stiff*mass[i]);
      contact_damp[i] = contact_damp_ratio*critical_damp;
      x[i]   = x0[i];
      y[i]   = y0[i];
      z[i]   = z0[i];
      vx[i] *= initial_velocity_fraction;
      vy[i] *= initial_velocity_fraction;
      vz[i] *= initial_velocity_fraction;
    } // for i=1,...,num_particles
    
    // (Optionally) constrain the particle positions to fall within the CV
    if (cv_dtheta > 0.0) {
      enforce_periodic_control_volume_BCs();
    }
    
  } // initialize()
  
  // ---------------------------------------------------------------------- //

  // Estimate and return the limiting stable time step size for all particles
  double stable_time_step(void) {
    
    // restrict the stable time step based on the smallest stable dt of any particle
    double stable_dt = std::numeric_limits<double>::max();
    for (int i=0; i<num_particles; i++) {
      double frequency = std::sqrt(contact_stiff/mass[i]); // [1/s]
      stable_dt = std::min(stable_dt,0.5/frequency);       // [s]
    } // for i=1,...,num_particles

    return stable_dt;
    
  } // stable_time_step()
  
  // ---------------------------------------------------------------------- //

  // Zero forces acting on all particles
  void zero_forces(void) {

    // Loop over all particles and zero the applied forces
    for (int i=0; i<num_particles; i++) {
      fx[i] = 0.0;
      fy[i] = 0.0;
      fz[i] = 0.0;
    } // for i=1,...,num_particles
    
  } // zero_forces()
  
  // ---------------------------------------------------------------------- //

  // Apply body (gravitational acceleration) forces to all particles
  void apply_gravitational_forces(double gx, double gy, double gz) {

    // Loop over all particles and apply gravitational forces
    for (int i=0; i<num_particles; i++) {
      fx[i] += mass[i]*gx;
      fy[i] += mass[i]*gy;
      fz[i] += mass[i]*gz;
    } // for i=1,...,num_particles
    
  } // apply_gravitational_forces()
  
  // ---------------------------------------------------------------------- //

  // Apply drag forces to all particles
  void apply_drag_forces(std::unique_ptr<WindField>& wind_model, double time) {

    // declare persistent static data arrays
    static std::vector<double> vxf(num_particles);
    static std::vector<double> vyf(num_particles);
    static std::vector<double> vzf(num_particles);
    static std::vector<double> rhof(num_particles);
    
    // determine the fluid velocity and density at the current positions of all particles
    wind_model->get_fluid_velocity_and_density(num_particles,time,x.data(),y.data(),z.data(),
					       vxf.data(),vyf.data(),vzf.data(),rhof.data());

    // Loop over all particles
    for (int i=0; i<num_particles; i++) {

      // determine the velocity of the particle relative to the fluid
      double vx_rel = vx[i] - vxf[i];
      double vy_rel = vy[i] - vyf[i];
      double vz_rel = vz[i] - vzf[i];
      double vmag_rel = std::sqrt(vx_rel*vx_rel + vy_rel*vy_rel + vz_rel*vz_rel);
    
      // compute and apply the drag force
      const double pi = 2.0*std::acos(0.0);
      double area = pi*radius[i]*radius[i];
      double fdrag = -0.5*drag_coeff*area*rhof[i]*vmag_rel;
      fx[i] += fdrag*vx_rel;
      fy[i] += fdrag*vy_rel;
      fz[i] += fdrag*vz_rel;

    } // for i=1,...,num_particles
    
  } // apply_drag_forces()
  
  // ---------------------------------------------------------------------- //

  // Apply forces to particles due to contact with a rigid boundary plane
  void apply_boundary_forces(double xp = 0.0, double yp = 0.0, double zp = 0.0,
			     double nx = 0.0, double ny = 0.0, double nz = 1.0) {
    // {xp,yp,zp}: reference point defined on the boundary plane
    // {nx,ny,nz}: unit surface normal to the boundary plane

    // Loop over all particles
    for (int i=0; i<num_particles; i++) {

      // compute the normal gap functio for the current particle
      double gap = (x[i] - xp)*nx + (y[i] - yp)*ny + (z[i] - zp)*nz;

      // compute the normal gap rate (assuming the boundary plane is not moving) for the current particle
      double dgap_dt = vx[i]*nx + vy[i]*ny + vz[i]*nz;

      // if interpenetration detected, apply contact force to the current particle
      // the contact force includes combined stiffness and damping effects
      if (gap < 0.0) {
	double contact_force = - contact_stiff*gap - contact_damp[i]*dgap_dt;
	fx[i] += contact_force*nx;
	fy[i] += contact_force*ny;
	fz[i] += contact_force*nz;
      }

    } // for i=1,...,num_particles
    
  } // apply_boundary_forces()
  
  // ---------------------------------------------------------------------- //

  // Update the velocities and positions of all particles by integrating the
  // equations of motion over a finite time step of size dt
  void integrate_equations_of_motion(double dt) {

    // loop over all particles
    for (int i=0; i < num_particles; i++) {
      
      // update the current particle's velocity and position
      double dt_inv_mass = dt/mass[i];
      vx[i] += dt_inv_mass*fx[i];
      vy[i] += dt_inv_mass*fy[i];
      vz[i] += dt_inv_mass*fz[i];
      x[i]  += dt*vx[i];
      y[i]  += dt*vy[i];
      z[i]  += dt*vz[i];
      
    } // for(i=0...num_particles)

    // (Optionally) constrain the particle positions to fall within the CV
    if (cv_dtheta > 0.0) {
      enforce_periodic_control_volume_BCs();
    }
    
  } // integrate_equations_of_motion()
  
  // ---------------------------------------------------------------------- //

  // Constraint the positions of all particles to lie within the
  // optionally defined periodic control volume
  void enforce_periodic_control_volume_BCs(void) {
    const double pi = 4.0 * std::atan(1.0);

    // loop over all particles
    for (int i=0; i < num_particles; i++) {

      // determine the coordinates of the current particle measured relative
      // to the vertex of the periodic control volume (the center of the vortex)
      double dx = x[i] - cv_x0;
      double dy = y[i] - cv_y0;

      // determine the relative radial and angular coordinates of the particle
      double r = std::sqrt(dx*dx+dy*dy);
      double theta = std::atan2(dy,dx);

      // determine the relative (positive) angle between the particle's angular coordinate and cv_theta_in
      double dtheta = theta - cv_theta_in;
      if (dtheta < 0.0) dtheta += 2.0*pi; // Add 2*pi if the angle is negative

      // determine whether the particle's angular coordinate lies outside the CV
      double normalized_dtheta = dtheta / cv_dtheta;
      if (normalized_dtheta > 1.0) { // particle lies outside the CV: enforce periodic BC

	// determine the particle's current radial and circumferential velocity components
	// { vR } = [ +cosT +sinT ] * { vx }
	// { vT }   [ -sinT +cosT ]   { vy }
	double cosT = std::cos(theta);
	double sinT = std::sin(theta);
	double vR = +vx[i]*cosT+vy[i]*sinT;
	double vT = -vx[i]*sinT+vy[i]*cosT;

	// confine the particle's relative angular coordinate to lie within the CV
        normalized_dtheta -= int(normalized_dtheta);
	dtheta = normalized_dtheta*cv_dtheta;
	theta = cv_theta_in + dtheta;

	// transform the particle's x,y coordinates
        cosT = std::cos(theta);
	sinT = std::sin(theta);
	dx = r*cosT;
	dy = r*sinT;
	x[i] = cv_x0 + dx;
	y[i] = cv_y0 + dy;

	// transform the particles' x,y velocity
	// [ +cosT -sinT ] * { vR } = { vx }
	// [ +sinT +cosT ]   { vT }   { vy }
	vx[i] = +vR*cosT-vT*sinT;
	vy[i] = +vR*sinT+vT*cosT;
      }
      
    } // for(i=0...num_particles)
    
  } // enforce_periodic_control_volume_BCs()

  // --------------------- Declare public data members -------------------- //

  // Common constants defined for all particles
  int num_particles;            // The total number of particles
  double drag_coeff;            // Drag coefficient for all particles
  double contact_stiff;         // Contact spring stiffness for all particles [N/m] = [kg/s^2]
  double contact_damp_ratio;    // Contact damping ratio for all particles
  bool particles_active = true; // Logical flag indicating if all particles are active/unfrozen (true - default) or inactive/frozen (false)

  // Optional parameters specifying the wedge-shaped control volume
  double cv_x0 = 0.0; // in-plane x-coordinate of the vertex of the wedge (should coincide with the center of the vortex)
  double cv_y0 = 0.0; // in-plane y-coordinate of the vertex of the wedge (should coincide with the center of the vortex)
  double cv_theta_in = 0.0; // [radians] the angular coordinate (measured relative to the x-axis) at which the inflow CV surface is defined
  double cv_dtheta = 0.0; // [radians] the positive angular dimension (measured in radians relative to theta_in) at which the outflow surface is defined

  // Data defined separately for each particle
  std::vector<bool>   ghost;        // The particle's "ghost" status (true if the particle is a ghost - no contact forces applied to structure; false otherwise)
  std::vector<double> mass;         // The masses defined for all particles
  std::vector<double> radius;       // The radii of all particles
  std::vector<double> contact_damp; // The contact damping for all particles
  std::vector<double> x0, y0, z0;   // The initial spatial coordinates of all particles
  std::vector<double>  x,  y,  z;   // The current spatial coordinates of all particles
  std::vector<double> vx, vy, vz;   // The current velocity of all particles
  std::vector<double> fx, fy, fz;   // The current forces applied to all particles
  
  // ---------------------------------------------------------------------- //

}; // Particles


// ======================================================================== //


#endif /* PARTICLES_H */
