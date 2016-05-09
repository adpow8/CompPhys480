#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "body.h"

class solarsystem {
public:
    int body_count;
    body* system;
    double mass_total;
    double composition[4];
    double comvelocity[4];

    solarsystem();

    ~solarsystem();

    void add(body sagan);
   /* void update(double**&, int);*/
    void relpositions(double**&);
    void forces(double*&, const double**&);
    void verlet(double**&, int, double, int);
    void RK4(double**&, int, double, int);
    void print();
};

#endif // SOLARSYSTEM_H
