#ifndef GALACTICSYSTEM_H
#define GALACTICSYSTEM_H

#include "body.h"

class galacticsystem {
public:
    int body_count;
    body* system;
    double mass_total;
    double composition[4];
    double comvelocity[4];

    galacticsystem();

    ~galacticsystem();

    void add(body sagan);
   /* void update(double**&, int);*/
    void relpositions(double**&);
    void forces(double*&, const double**&);
    void verlet(double**&, int, double);
    void RK4(double**&, int, double);
    void print();
};

#endif // galacticsystem_H
