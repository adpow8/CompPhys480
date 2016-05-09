# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include "armadillo"

using namespace arma;
using namespace std;

int main(/*int argc, char *argv[]*/)
{
    // initialize
    double pi = M_PI;

    double ti = 0., tf = 10.;   // time in years
    double n = (tf-ti)*12000.;    // number of iterations
    double h = (tf-ti)/n;       // step length
    // cout << h << endl;
    double phi = pi/2;      // initial direction angle of earth's orbit
    double v0 = 6.39179;    // AU/yr

    double xInit = -9.903414152297834e-01;         // in AU
    double yInit = -7.716244651654738e-02;               // in AU
    double vxInit = 1.082807915018358e-03/365;    // in AU/yr
    //cout << vxInit << endl;
    double vyInit = -1.720133424980736e-02/365;    // in AU/yr
    //cout << vyInit << endl;

    double rInit = sqrt(xInit*xInit + yInit*yInit); // initial dist of earth from sun

    // initialize x/y positions/velocities
    vec x(n+1), y(n+1), vx(n+1), vy(n+1), r(n+1);

    x(0) = xInit;
    y(0) = yInit;
    vx(0) = vxInit;
    vy(0) = vyInit;
    r(0) = rInit;

    // loop to find positions/velocities using Forces

    for (int i = 0; i < n; i++)
    {
        x(i+1) = x(i) + h*vx(i);
        y(i+1) = y(i) + h*vy(i);
        r(i+1) = sqrt( x(i+1)*x(i+1) + y(i+1)*y(i+1) );
        vx(i+1) = vx(i) - (4*pi*pi*x(i)*h) / (r(i)*r(i)*r(i));
        vy(i+1) = vy(i) - (4*pi*pi*y(i)*h) / (r(i)*r(i)*r(i));
        // cout << i << endl;
        //cout << h*vx(i) << " " << h*vy(i) << endl;
    }

    cout << "RESULTS:" << endl;
    cout << "         r:              x:             y:             vx:            vy:  " << endl;

    for (int i=0;i<=n;i+=100)
    {
    cout << setw(15) << setprecision(8) << r(i);
    cout << setw(15) << setprecision(8) << x(i);
    cout << setw(15) << setprecision(8) << y(i);
    cout << setw(15) << setprecision(8) << vx(i);
    cout << setw(15) << setprecision(8) << vy(i) << endl;
    }
}
