#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H
#include <armadillo>
#include "planet.h"
#include <vector>

using std::vector;
using namespace arma;

class solarsystem
{
public:

    int num_planets = 0;

    solarsystem();
    vector<planet> all_planets;
    void add(planet n);
    void print_position(std::ofstream &output, std::ofstream &output2, vector<planet> vec);
    void print_position(std::ofstream &output, std::ofstream &output2, vector<planet> vec, int n);
    void synctroniz(vector<planet> vec, arma::mat &ma);
    void insert_data(vector<planet> vec, arma::mat &ma);
    void RK4solver(vector<planet> vec, double h, double tmax);
    void Verletsolver(vector<planet> vec, double h, double tmax);
    void sum_matrix(mat &result, double coeff_one, mat &first, double coeff_two, mat &second, int n);
    void printmat(mat &ma, int n);
    double force(double x, double y, double z, double Mothers);
    void derivate(mat &dat, mat &de, int n);

};

#endif // SOLARSYSTEM_H
