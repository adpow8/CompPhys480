#ifndef PLANET_H
#define PLANET_H


class planet
{
public:

    double position[3];
    double velocity[3];
    double planetMass;

    planet();

    planet(double mass, double x,double y, double z, double vx, double vy, double vz);
};

#endif // PLANET_H
