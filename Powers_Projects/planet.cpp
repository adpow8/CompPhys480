#include "planet.h"

planet::planet(double mass, double x,double y, double z, double vx, double vy, double vz){

    // ratio of mass to the sun.
    planetMass = mass;

    // in AU
    position[0] = x;
    position[1] = y;
    position[2] = z;

    // in AU/yr
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    };
