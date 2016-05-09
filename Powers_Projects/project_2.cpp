#include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
#include "armadillo"

using namespace arma;
using namespace std;
ofstream ofile;

double offdiag (mat A, int *p, int *q, int n);
double oneEpot(const double r);
double twoEpot(const double r, const double omega);
void Jacobi_rotate (mat &A, mat &R, int k, int l, int n);

// main program
int main(int argc, char *argv[])
{
    // determine size of matrix and declare
    int n = 375;
    mat A = zeros<mat>(n,n);
    char *outfilename = argv[1];
    double omega = 0.01827;

    /*if(argc < 2){
        cout << "Bad usage: " << argv[0] <<
                "read output file and n (int) on same line" << endl;
        exit(1);
    }
    else{
        outfilename = argv[1];
        n = atoi(argv[2]);
    }*/
    //set potential
    vec Vone(n);
    vec Vtwo(n);
    for (int i = 0; i < n; i++)
    {
        Vone(i) = oneEpot(i);
    }
    for (int i = 0; i < n; i++)
    {
        Vtwo(i) = twoEpot(i, omega);
    }

    // set up eigenvector matrix
    mat E = zeros<mat>(n,n);

    // set constants
    double h = 1;    // choosing units to make hbar = 1
    double k = 2/(h*h);
    double e = -1/(h*h);

    // fill  Potential, symmetric, tridiagonal, matrix/make symmetric
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                A(i,i) = k + Vtwo(i);
                E(i,j) = 1;
            }
            else if (i == j+1 || j == i+1)
            {
                A(i,j) = e;
            }
        }
    }

    // cout << A << endl;

    double tolerance = 1.0E-10;
    int iterations = 0;
    int maxIterations = 1000000;
    double maxOffDiag = 1;

    //Start timer

    clock_t start, finish;

    start = clock();

    while (maxOffDiag > tolerance && iterations <= maxIterations)
    {
        int p, q;
        // determine initial p and q values
        maxOffDiag = offdiag(A, &p, &q, n);


        // Run Jacobi transformation
        Jacobi_rotate(A, E, p, q, n);
        iterations++;
    }
    finish = clock();

    //check to see if program is working
    /*cout << "RESULTS:" << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout <<"E = " << setw(15) << setprecision(8) << E << endl;
    cout <<"A = " << setw(15) << setprecision(8) << A << endl;*/


    // Open file and write results to file:

    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);

                ofile << setw(15) << setprecision(8) << A << endl;

                ofile << setw(15) << setprecision(8) << E << endl;


    ofile.close();

    cout << "iterations: " << iterations << endl;
    cout << "time (s): " << (finish-start)/(double) CLOCKS_PER_SEC << endl;

}



// find the max point within the matrix not on diagonal
double offdiag (mat A, int *p, int *q, int n)
{
    double max = 0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = i+1; j < n; ++j)
        {
            double aij = fabs(A(i,j));

            if (aij > max)
            {
                max = aij; *p = i; *q = j;
            }
        }
    }
    return max;
} // end of offdiag function

double oneEpot(const double r)
{
    return r*r;
}

double twoEpot(const double r, const double omega)
{
    return omega*omega*r*r + 1/r;
}


// rotate matrix to create tridiagonal
void Jacobi_rotate (mat &A, mat &R, int k, int l, int n)
{
    double s, c;
    if (A(k,l) != 0.0)
    {
        double t, tau;
        tau = (A(l,l) - A(k,k)) / (2*(A(k,l)));

        if (tau >= 0)
        {
            t = 1.0 / (tau + sqrt(1+ tau*tau));
        }
        else
        {
            t = -1.0 / (-tau + sqrt(1+ tau*tau));
        }

        c = 1 / sqrt(1+t*t);
        s = c*t;
    }
    else
    {
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;    // hard coding non-diagonal elements by hand
    A(l,k) = 0.0;    // same thing here

    for (int i = 0; i < n; i++)
    {
        if (i != k && i != l)
        {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }

        // The new eigenvectors:
        r_ik = R(i,k);
        r_il = R(i,l);

        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }

    return;

} // end of Jacobi_rotate function


