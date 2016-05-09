# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
#include "planet.h"
# include "solarsystem.h"
# include "armadillo"

using namespace arma;
using namespace std;

int main(int argc, char *argv[])
{


    // initialize planets (data from 2016_Mar_25 00:00:00
    planet Sun(1, 0, 0, 0, 0, 0, 0);
    planet Mercury(1.645e-7, 3.3933e-1, 5.4839e-2, -2.6650e-2, -3.6230, 1.0600e1, 1.1983);
    planet Venus(2.43e-6, 5.4442e-1, -4.8208e-1, -3.8026e-2, 4.8472, 5.4991, -2.0432e-1);
    planet Earth(2e-6, -9.9411e-1, -7.9150e-2, -6.1253e-6, 3.9524e-1, -6.2810, 2.7898e-5);
    planet Mars(3.2e-7, -1.3450, -8.5374e-1, 1.5118e-2, 2.9283, -3.8756, -1.5308e-1);
    planet Jupiter(9.5e-4, -5.3356, 1.0045, 1.1522e-1, -5.4288e-1, -2.5791, 2.2862e-2);
    planet Saturn(2.84e-4, -3.3012, -9.4571, 2.9579e-1, 1.8114, -6.7946e-1, -6.0240e-2);
    planet Uranus(4.34e-5, 1.8755e1, 6.8487, -2.1742e-1, -5.0295e-1, 1.2790, 1.1297e-2);
    planet Neptune(5.1e-5, 2.8043e1, -1.0530e1, -4.2939e-1, 3.9522e-1, 1.0770, -3.1407e-2);
    planet Pluto(6.55e-9, 8.7912, -3.1864e1, 8.6557e-1, 1.1309, 6.6798e-2, -3.3768e-1);

    solarsystem Sol;

    Sol.add(Sun);
    Sol.add(Mercury);
    Sol.add(Venus);
    Sol.add(Earth);
    Sol.add(Mars);
    Sol.add(Jupiter);
    Sol.add(Saturn);
    Sol.add(Uranus);
    Sol.add(Neptune);
    Sol.add(Pluto);

    int elements = Sol.number_planets;
    cout << "number of element" << elements<< endl;
    Sol.RK4solver(Sol.all_planets, 0.001, 100 );

}
