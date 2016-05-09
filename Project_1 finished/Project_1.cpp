#include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
#include "armadillo"
#include "time.h"

using namespace arma;
using namespace std;
ofstream ofile;

double Real(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

double f(double x) {return 100*exp(-10*x);}

int main(int argc, char *argv[])
{
    /*
     Program to be used to solve the linear algebra equation
     A[n,n] * x[n] = w[n]
     when A is a tridiagonal matrix.
    */

    // Declaration of initial variables:

    char *outfilename;
    int n;

    // Read in output file and n,
    // abort if there are too few command-line arguments:

    if( argc <= 2 ){
      cout << "Bad Usage: " << argv[0] <<
          " read output file and n (int) on same line" << endl;
      exit(1);
    }
    else{
      outfilename = argv[1]; // first command line argument.
      n = atoi(argv[2])-1;  // second command line argument.
    }

    //Vectors making up the tridiagonal matrix A
    //Zeroth element included to make indexing easy

    vec a(n+1); vec b(n+1); vec c(n+1);

    //Temporary element needed for Gaussian elimination

    vec dTemp(n+1);

    //Constants within the problem

    double h = 1.0/(n+1.0);
    vec x(n+2);
    vec bTwi(n+1);

    bTwi(0) = 0;

    //initialize real solution and approximated solution:

    vec u(n+2); //Analytical solution
    vec v(n+2); //Numerical solution

    //Make indexing easy
    u(0) = 0;
    v(0) = 0;

    //Fill up x-vector

    for (int i = 0; i <= n+1; i++)
    {
        x(i) = i*h;
        // Could print results to check:
        //cout << "x = " << x[i] << " and " << "h^2*f(x) = " << h*h*f(x[i]) << endl;
    }

    //initialize vector components of matrix to be solved

    for (int i = 1; i<= n; i++)
    {
        bTwi(i) = h * h * f( x(i) );
        // Could print here to check:
        //cout << "b_twidd = " << b_twidd[i] << "for x = " << x[i] << endl;

        u(i) = Real(x(i));
        //cout << "u = " << u[i] << " for x = " << x[i] <<  endl;

        a(i) = 1;
        b(i) = 2;
        c(i) = 1;
    }

    a(1) = 0;
    c(n) = 0;

    //Algorithm for finding v:
    // a(i)*v(i-1) + b(i)*v(i) + c(i)*v(i+1) = b_twidd(i)

    //Start timer for homebrew algorithm

    clock_t start, finish, startLU, finishLU;

    start = clock();

    // Row reduction; forward substitution:

    double bTemp = b(1);
    v(1) = bTwi(1) / bTemp;
    for (int i = 2; i <= n; i++)
    {
        //Temporary value needed in next loop
        dTemp(i) = c(i-1) / bTemp;
        //Temporary diagonal element
        bTemp = b(i) - a(i)*dTemp(i);
        // Update right hand side of equation
        v(i) = (bTwi(i) - v(i-1) * a(i)) / bTemp;
    }

    // Row reduction backward substitution:

    for (int i = n-1; i >= 1; i--)
    {
        v(i) -= dTemp(i+1) * v(i+1);
    }

    finish = clock();

    cout << (finish-start)/(double) CLOCKS_PER_SEC << endl;



    return 0;


}
