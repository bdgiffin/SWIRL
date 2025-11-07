#ifndef WIND_FIELD_H
#define WIND_FIELD_H

#include "Parameters.h"
#include "Logger.h"
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <memory>


// ======================================================================== //


// Abstract base class for a generic wind field model
class WindField {
public:

  // ------------------- Declare public member functions ------------------ //

  WindField(void) {} // Empty default constructor method

  // Default (blank) initialization routine
  virtual void initialize(void) { }

  virtual ~WindField() { } // Declare destructor as virtual

  // Pure virtual method to compute the fluid velocity and density at a specified time,
  // and at multiple evaluation points simultaneously
  virtual void get_fluid_velocity_and_density(int num_points, double time,
				              const double* x, const double* y, const double* z,
				              double* vx, double* vy, double* vz, double* rhof) = 0;
  
  // ---------------------------------------------------------------------- //

}; // WindField


// ======================================================================== //


// Derived class for a vortex wind field (Baker & Sterling, 2017)
// For details on the originally proposed model:
// - https://www.sciencedirect.com/science/article/pii/S0167610517301174
// For details on how translational effects should be incorporated into the model:
// - https://www.sciencedirect.com/science/article/pii/S0167610517301186
class BakerSterlingVortex : public WindField {
public:

  // ------------------- Declare public member functions ------------------ //
  
  ~BakerSterlingVortex() { }

  // Parameterized constructor method
  BakerSterlingVortex(double* parameters) : WindField() {
    DEBUG(std::cout << "Creating new BakerSterlingVortex WindField model" << std::endl;)
    Um    = parameters[0];  // [m/s] reference radial velocity
    rm    = parameters[1];  // [m]   reference radius
    zm    = parameters[2];  // [m]   reference height
    S     = parameters[3];  // swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
    K     = S*(2.0/std::log(2.0));
    gamma = parameters[4];
    delta = zm/rm;
    rho0  = parameters[5];  // [kg/m^3] reference density of air at STP
    xc0   = parameters[6];  // [m]
    yc0   = parameters[7];  // [m]
    zc0   = parameters[8];  // [m]
    vxc   = parameters[9];  // [m/s]
    vyc   = parameters[10]; // [m/s]
    vzc   = parameters[11]; // [m/s]
  } // BakerSterlingVortex()

  // Parameterized constructor method
  BakerSterlingVortex(Parameters& parameters) : WindField() {
    if (parameters.count("initial_center") > 0) {
      std::vector<double> xyzc = parameters["initial_center"];
      xc0 = xyzc[0];
      yc0 = xyzc[1];
      zc0 = xyzc[2];
    }
  } // BakerSterlingVortex()

  // Virtual method implementation to compute the fluid velocity and density at a specified time,
  // and at multiple evaluation points simultaneously
  virtual void get_fluid_velocity_and_density(int num_points, double time,
				              const double* x, const double* y, const double* z,
				              double* vx, double* vy, double* vz, double* rhof) {

    // Define the shifted center of the vortex at the current evaluation time
    double xct = xc0 + vxc*time; // [m]
    double yct = yc0 + vyc*time; // [m]
    double zct = zc0 + vzc*time; // [m]

    // Loop over all evaluation points
    for (int i=0; i<num_points; i++) {

      // Compute normalized radial and height coordinates
      double rp = std::sqrt((x[i]-xct)*(x[i]-xct) + (y[i]-yct)*(y[i]-yct));
      double inv_rp = 1.0/(rp+std::numeric_limits<double>::min());
      double cosp = (x[i]-xct)*inv_rp;
      double sinp = (y[i]-yct)*inv_rp;
      double rbar = rp/rm;
      double zbar = (z[i]-zct)/zm;

      // Ignore points that are sufficiently far away
      if ((rbar > rbar_max) || (zbar > zbar_max)) {
	vx[i]   = 0.0;
	vy[i]   = 0.0;
	vz[i]   = 0.0;
	rhof[i] = 0.0;
	continue;
      }

      // Compute normalized radial, tangential, and vertical velocity of the vortex
      double one_rbar2 = 1.0 + rbar*rbar;
      double one_zbar2 = 1.0 + zbar*zbar;
      double log_one_zbar2 = std::log(one_zbar2);
      double Ubar = -4.0*rbar*zbar/(one_rbar2*one_zbar2);
      double Vbar = K*std::pow(rbar,gamma-1.0)*std::pow(log_one_zbar2,0.5*gamma)/std::pow(one_rbar2,0.5*gamma);
      double Wbar = 4.0*delta*log_one_zbar2/(one_rbar2*one_rbar2);

      // Compute the x,y,z components of the fluid velocity
      vx[i] = vxc + Um*(Ubar*cosp-Vbar*sinp);
      vy[i] = vyc + Um*(Ubar*sinp+Vbar*cosp);
      vz[i] = vzc + Um*Wbar;

      // Assume the density is constant
      rhof[i] = rho0;

      // Check for nan/inf values
      DEBUG(
      if ((inv_rp != inv_rp) || (log_one_zbar2 != log_one_zbar2) || (Ubar != Ubar) || (Vbar != Vbar) || (Wbar != Wbar)) {
	std::cout << "Bad BakerSterlingVortex velocity computed at evaluation point " << i << std::endl;
	std::cout << " time " <<  time << std::endl;
	std::cout << "    x " <<  x[i] << std::endl;
	std::cout << "    y " <<  y[i] << std::endl;
	std::cout << "    z " <<  z[i] << std::endl;
	std::cout << " inv_rp " << inv_rp << std::endl;
	std::cout << " zbar " << zbar << std::endl;
	std::cout << " one_zbar2 " << one_zbar2 << std::endl;
	std::cout << " log_one_zbar2 " << log_one_zbar2 << std::endl;
	std::cout << " Ubar " <<  Ubar << std::endl;
	std::cout << " Vbar " <<  Vbar << std::endl;
	std::cout << " Wbar " <<  Wbar << std::endl;
	std::exit(1);
      }
	    )

    } // for i=1,...,num_points

  } // get_fluid_velocity_and_density()

  // ---------------------------------------------------------------------- //
  
private:
  
  // -------------------- Declare private data members -------------------- //

  // Define reference values and constants for use in dimensionless evaluations
  double Um; // [m/s] reference radial velocity
  double rm; // [m]   reference radius
  double zm; // [m]   reference height
  double S;  // swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
  double K;
  double gamma;
  double delta;
  double rho0; // [kg/m^3] reference density of air at STP

  // Define the center of the vortex, and its translational velocity
  double xc0; // [m]
  double yc0; // [m]
  double zc0; // [m]
  double vxc; // [m/s]
  double vyc; // [m/s]
  double vzc; // [m/s]

  // Define limits
  double rbar_max = std::numeric_limits<double>::max();
  double zbar_max = std::numeric_limits<double>::max();

}; // BakerSterlingVortex


// ======================================================================== //


// Derived class for a Rankine vortex wind field
class RankineVortex : public WindField {
public:

  ~RankineVortex() { }

  // ------------------- Declare public member functions ------------------ //

  // Parameterized constructor method
  RankineVortex(double* parameters) : WindField() {
    DEBUG(std::cout << "Creating new RankineVortex WindField model" << std::endl;)
    Um    = parameters[0];  // [m/s] reference tangential velocity
    rm    = parameters[1];  // [m]   reference radius
    rc    = parameters[2];  // [m]   reference height
    E     = parameters[3];  // swirl ratio (ratio of max circumferential velocity to radial velocity at reference height)
    //      parameters[4];  // (unused)
    rho0  = parameters[5];  // [kg/m^3] reference density of air at STP
    xc0   = parameters[6];  // [m]
    yc0   = parameters[7];  // [m]
    zc0   = parameters[8];  // [m]
    vxc   = parameters[9];  // [m/s]
    vyc   = parameters[10]; // [m/s]
    vzc   = parameters[11]; // [m/s]
  } // RankineVortex()

  // Virtual method implementation to compute the fluid velocity and density at a specified time,
  // and at multiple evaluation points simultaneously
  virtual void get_fluid_velocity_and_density(int num_points, double time,
				              const double* x, const double* y, const double* z,
				              double* vx, double* vy, double* vz, double* rhof) {

    // Define the shifted center of the vortex at the current evaluation time
    double xct = xc0 + vxc*time; // [m]
    double yct = yc0 + vyc*time; // [m]
    //double zct = zc0 + vzc*time; // [m]

    // Loop over all evaluation points
    for (int i=0; i<num_points; i++) {

      // Compute normalized radial and height coordinates
      double rp = std::sqrt((x[i]-xct)*(x[i]-xct) + (y[i]-yct)*(y[i]-yct));
      double inv_rp = 1.0/(rp+std::numeric_limits<double>::min());
      double cosp = (x[i]-xct)*inv_rp;
      double sinp = (y[i]-yct)*inv_rp;
      double rbar = rp/rc;

      // Compute normalized tangential velocity of the vortex
      double vtf;
      if (rp<rm) {
	if (rp<rc) {
	  vtf = Um*rbar;
	} else {
	  vtf = Um*std::pow(1.0/(rbar+std::numeric_limits<double>::min()),E);
	}
      } else {
	vtf = 0.0;
      }

      // Compute the x,y,z components of the fluid velocity
      vx[i] = vxc - vtf*sinp;
      vy[i] = vyc + vtf*cosp;
      vz[i] = vzc;

      // Assume the density is constant
      rhof[i] = rho0;

    } // for i=1,...,num_points

  } // get_fluid_velocity_and_density()

  // ---------------------------------------------------------------------- //
  
private:
  
  // -------------------- Declare private data members -------------------- //

  // Define reference values and constants for use in dimensionless evaluations
  double Um;   // [m/s] reference radial velocity
  double rm;   // [m]   reference outer radius
  double rc;   // [m]   reference core radius
  double E;    // decay index
  double rho0; // [kg/m^3] reference density of air at STP

  // Define the center of the vortex, and its translational velocity
  double xc0; // [m]
  double yc0; // [m]
  double zc0; // [m]
  double vxc; // [m/s]
  double vyc; // [m/s]
  double vzc; // [m/s]

}; // RankineVortex


// ======================================================================== //


// Factory method to create a new WindField model from (generic) parameterized inputs
std::unique_ptr<WindField> new_WindField(const char* type_cstr, double* parameters) {

  // Attempt to create a new wind field model given the passed input parameters
  std::string type(type_cstr);
  if        (type == "BakerSterlingVortex") {
    return std::make_unique<BakerSterlingVortex>(parameters);
  } else if (type == "RankineVortex") {
    return std::make_unique<RankineVortex>(parameters);
  } else {
    std::cerr << "ERROR in `" <<  __func__ << "`; unrecognized WindField type: " << type << std::endl;
    return nullptr;
  }
  
} // new_WindField()


// ======================================================================== //


#endif /* WIND_FIELD_H */
