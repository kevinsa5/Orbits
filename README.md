Orbits Notes:
=============

Development is done with Matlab R2013b on x64 Linux. 

Numeric Integration Methods
---------------------------

+ Forward Euler 
  * local error O(h^2), global error O(h)

+ Runge-Kutta 4 
  * local error O(h^5), global error O(h^4)

+ Velocity Verlet 
  * local error in position O(h^4), velocity O(h^2)
  * global error in position O(h^2), velocity O(h^2)


Plans
-----

+ Include a graph of total kinetic energy and angular momentum over time to show which methods conserve which quantitites

+ Figure out a better way to do initial conditions for the spacecraft

+ Shore up the hard-coded values for the planets

+ Improve the Euler method so it uses force-symmetry rather than recalculating

+ Collision/absorption checking

+ Include moons
