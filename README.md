Orbits Notes:
=============

Development is done with Matlab R2013b on x64 Linux. 

Numeric Integration Methods
---------------------------
+ Forward Euler 
  * local error O(h<sup>2</sup>), global error O(h)
+ Runge-Kutta 4 
  * local error O(h<sup>5</sup>), global error O(h<sup>4</sup>)
+ Velocity Verlet 
  * local error in position O(h<sup>4</sup>), velocity O(h<sup>2</sup>)
  * global error in position O(h<sup>2</sup>), velocity O(h<sup>2</sup>)

Cool Features
-------------

+ The sun's initial momentum is initialized such that the system's initial net linear momentum is zero and its center of mass is stationary

+ Collisions are checked and absorptions are done based on mass with conservation of linear momentum

+ Starting positions for planets+moons available for a wide selection of dates

Improvements To-Do
-----
+ Rewrite the tight loops in C or Fortran, since RK4 is just abominably slow.

+ Write a shooting method to solve the BVP (I'm here, I want to be there) and find the initial velocity for an unpowered (and powered?) object.

+ Include a graph of total kinetic energy and angular momentum over time to show which methods conserve which quantities

+ Improve the Euler method so it uses force-symmetry rather than recalculating

Screenshots
-----------
![screenshot](pictures/gui.png)

![screenshot](pictures/orbits.png)
