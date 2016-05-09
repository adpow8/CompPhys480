#include "body.h"
#include<cmath>
#include<fstream>

using namespace std;

body::body(){
    mass = 0.0;
    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;
    position[3] = 0.0;
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    velocity[2] = 0.0;
    velocity[3] = 0.0;
}

body::body(double Mass, double x, double y, double z, double vx, double vy, double vz){
    mass = Mass;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    position[3] = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    velocity[3] = sqrt(vx*vx + vy*vy + vz*vz);
}

