import math
m = 1  # meters
kg = 1  # kilograms
N = 1  # Newton
KN = 1000 * N  # kilonewtons
sec = 1  # seconds
  

# Dependent units
inch = 0.0254*m  # meters
kip = 4448.22*N  # Newtons (1 kip = 4448.22 N)
kips2in = 175.126836*kg/m  # kg/m (since 1 kip/in = 175.126836 kg/m)
sq_in = inch * inch  # square meters (1 in^2 = 0.00064516 m^2)
ksi = kip / sq_in  # Pascals (force per unit area) Pascals (since 1 ksi = 6.89476 MPa = 6.89476 * 10^6 Pa)
ft = inch/12  # meters (since 1 ft = 0.3048 meters)
mm = 0.001 * m  # meters (since 1 mm = 0.001 meters)
Mpa = N/mm**2
Gpa = Mpa * 10**3


# Constants
g = 9.81 * m / (sec * sec)  # acceleration due to gravity in m/s^2
pi = math.acos(-1)  # value of pi
