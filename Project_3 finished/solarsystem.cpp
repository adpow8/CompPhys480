#include "solarsystem.h"
#include "body.h"
#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>

using namespace std;

solarsystem::solarsystem()
{
    body_count = 0;
    system = nullptr;
    mass_total = 0;
    for(int i=0;i<4;i++){
        composition[i]=0.0;
        comvelocity[i]=0.0;
    }
}

solarsystem::~solarsystem()
{
    delete[] system;
}


void solarsystem::add(body sagan)
{
    double mass_old = mass_total;
    body* temp = new body[body_count+1];
    for(int i=0;i<body_count;i++)
    {
        temp[i].mass=system[i].mass;
        temp[i].position[0]=system[i].position[0];
        temp[i].position[1]=system[i].position[1];
        temp[i].position[2]=system[i].position[2];
        temp[i].position[3]=system[i].position[3];
        temp[i].velocity[0]=system[i].velocity[0];
        temp[i].velocity[1]=system[i].velocity[1];
        temp[i].velocity[2]=system[i].velocity[2];
        temp[i].velocity[3]=system[i].velocity[3];
    }

    temp[body_count].mass=sagan.mass;
    temp[body_count].position[0]=sagan.position[0];
    temp[body_count].position[1]=sagan.position[1];
    temp[body_count].position[2]=sagan.position[2];
    temp[body_count].position[3]=sagan.position[3];
    temp[body_count].velocity[0]=sagan.velocity[0];
    temp[body_count].velocity[1]=sagan.velocity[1];
    temp[body_count].velocity[2]=sagan.velocity[2];
    temp[body_count].velocity[3]=sagan.velocity[3];

    body_count++;

    delete[] system;
    mass_total += sagan.mass;
    for(int i=0;i<3;i++)
    {
        composition[i] = (composition[i]*mass_old + sagan.position[i]*sagan.mass)/mass_total;
        comvelocity[i] = (comvelocity[i]*mass_old + sagan.velocity[i]*sagan.mass)/mass_total;
    }
    composition[3] = sqrt(composition[0]*composition[0] + composition[1]*composition[1] + composition[2]*composition[2]);
    comvelocity[3] = sqrt(comvelocity[0]*comvelocity[0] + comvelocity[1]*comvelocity[1] + comvelocity[2]*comvelocity[2]);
    system = temp;
}

/*void solarsystem::update(double**& output, int i)
    {
        int j;
        mass_total = 0.0;

        for(j=0;j<4;j++)
        {
            composition[j] = 0.0;
            comvelocity[j] = 0.0;
        }

        for(j=0;j<body_count;j++)
        {
            system[j].position[0] = output[i][8*j+1];
            system[j].position[1] = output[i][8*j+2];
            system[j].position[2] = output[i][8*j+3];
            system[j].position[3] = output[i][8*j+4];
            system[j].velocity[0] = output[i][8*j+5];
            system[j].velocity[1] = output[i][8*j+6];
            system[j].velocity[2] = output[i][8*j+7];
            system[j].velocity[3] = output[i][8*j+8];
            mass_total += system[j].mass;
        }

        for(j=0;j<body_count;j++)
        {
            composition[0] += system[j].position[0]*system[j].mass;
            composition[1] += system[j].position[1]*system[j].mass;
            composition[2] += system[j].position[2]*system[j].mass;
            comvelocity[0] += system[j].velocity[0]*system[j].mass;
            comvelocity[1] += system[j].velocity[1]*system[j].mass;
            comvelocity[2] += system[j].velocity[2]*system[j].mass;
        }
        comvelocity[0] /= mass_total;
        comvelocity[1] /= mass_total;
        comvelocity[2] /= mass_total;

        composition[3] = sqrt(composition[0]*composition[0] + composition[1]*composition[1] + composition[2]*composition[2]);
        comvelocity[3] = sqrt(comvelocity[0]*comvelocity[0] + comvelocity[1]*comvelocity[1] + comvelocity[2]*comvelocity[2]);
    }*/

void solarsystem::print()
{
    cout << "mass=" << setw(15) << mass_total << endl;
    cout << "body_count=" << setw(15) << body_count << endl;

    for(int i=0;i<4;i++)
    {
        cout << "pos[" << i<< "]=" << setw(15) << composition[i] << endl;
    }

    for(int i=0;i<4;i++)
    {
        cout << "vel[" << i<< "]=" << setw(15) << comvelocity[i] << endl;
    }
}

void solarsystem::relpositions(double**& relcoord)
{
    int i, j;
    int bodies = body_count;

    for(i=0;i<bodies;i++)
    {
        for(j=0;j<bodies;j++)
        {
            if(i!=j)
            {
                relcoord[i][j*4] = system[i].position[0] - system[j].position[0];
                relcoord[i][j*4+1] = system[i].position[1] - system[j].position[1];
                relcoord[i][j*4+2] = system[i].position[2] - system[j].position[2];
                relcoord[i][j*4+3] = sqrt(relcoord[i][j*4]*relcoord[i][j*4]
                        + relcoord[i][j*4+1]*relcoord[i][j*4+1]+relcoord[i][j*4+2]*relcoord[i][j*4+2]);
            }
        }
    }
}

void solarsystem::forces(double*& force, const double**& relcoord){
    int j, k;
    int bodies = body_count;
    double rcubed, fourpisq;
    fourpisq = 4.0*M_PI*M_PI;

    for(j=0;j<bodies;j++){
        force[j*3] = 0.0;
        force[j*3+1] = 0.0;
        force[j*3+2] = 0.0;
    }

    for(j=0;j<bodies;j++)
    {
        for(k=0;k<bodies;k++)
        {
            if(j!=k)
            {
                rcubed = pow(relcoord[j][k*4+3],3.0);
                force[j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
                force[j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
                force[j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
            }
        }

        force[j*3] *= fourpisq;
        force[j*3+1] *= fourpisq;
        force[j*3+2] *= fourpisq;
    }
}

/*void solarsystem::initialize(double**& output, double*& mass, int steps)
            * {
        int i, j;
        int bodies = body_count;
        output = new double*[steps+1];
        mass = new double[bodies];
        for(int i=0;i<steps+1;i++)
        {
            output[i] = new double[bodies*8+1];
        }
        for(i=1;i<steps+1;i++)
        {
            for(j=0;j<bodies;j++)
            {
                output[i][j*8+1] = 0;
                output[i][j*8+2] = 0;
                output[i][j*8+3] = 0;
                output[i][j*8+4] = 0;
                output[i][j*8+5] = 0;
                output[i][j*8+6] = 0;
                output[i][j*8+7] = 0;
                output[i][j*8+8] = 0;
            }
        }

        for(j=0;j<bodies;j++)
        {
            mass[j] = system[j].mass;
            output[0][j*8+1] = system[j].position[0];
            output[0][j*8+2] = system[j].position[1];
            output[0][j*8+3] = system[j].position[2];
            output[0][j*8+4] = system[j].position[3];
            output[0][j*8+5] = system[j].velocity[0];
            output[0][j*8+6] = system[j].velocity[1];
            output[0][j*8+7] = system[j].velocity[2];
            output[0][j*8+8] = system[j].velocity[3];
        }
    }*/

/*************************
      Class function solvers
    *************************/

void solarsystem::RK4(double**& output, int steps, double tmax, int sun_choice){
    int i, j, k, m, n, bodies;
    bodies = body_count;
    double h, halfh, halfhsq, rcubed, fourpisq;
    double* k2, *k3;
    h = (tmax-0.0)/((double) (steps));
    halfh = 0.5*h;
    halfhsq = halfh*h;
    fourpisq = 4.0*M_PI*M_PI;

    //placeholders arrays
    k2 = new double[bodies*6];
    k3 = new double[bodies*6];
    for(i=0;i<bodies*6;i++){
        k2[i] = 0;
        k3[i] = 0;
    }
    //relcoord to provide relative coordinates, x, y, z, r, of body i with respect to body j
    double** relcoord = new double*[bodies];
    for(i=0;i<bodies;i++){
        relcoord[i] = new double[bodies*4];
    }
    //force provides fx,fy,fz for each body at step i and i+1
    double** force = new double*[2];
    for(i=0;i<2;i++){
        force[i] = new double[bodies*3];
    }
    //initialize all matrices to zero
    for(i=0;i<steps+1;i++){
        for(j=0;j<bodies*8+1;j++){
            output[i][j]=0.0;
        }
    }
    for(i=0;i<bodies;i++){
        for(j=0;j<bodies*4;j++){
            relcoord[i][j]=0.0;
        }
    }
    for(i=0;i<2;i++){
        for(j=0;j<bodies*3;j++){
            force[i][j]=0.0;
        }
    }
    //reinitialize relcoord to initial conditions
    //technically double calculating--consider revising
    for(i=0;i<bodies;i++){
        for(j=0;j<bodies;j++){
            if(i!=j){
                relcoord[i][j*4] = system[i].position[0] - system[j].position[0];
                relcoord[i][j*4+1] = system[i].position[1] - system[j].position[1];
                relcoord[i][j*4+2] = system[i].position[2] - system[j].position[2];
                relcoord[i][j*4+3] = sqrt(relcoord[i][j*4]*relcoord[i][j*4] + relcoord[i][j*4+1]*relcoord[i][j*4+1]+relcoord[i][j*4+2]*relcoord[i][j*4+2]);
            }
        }
    }
    //initialize forces (k1v)
    for(j=0;j<bodies;j++){
        for(k=0;k<bodies;k++){
            if(j!=k){
                rcubed = pow(relcoord[j][k*4+3],3.0);
                force[0][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
                force[0][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
                force[0][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
            }
        }
        force[0][j*3] *= fourpisq;
        force[0][j*3+1] *= fourpisq;
        force[0][j*3+2] *= fourpisq;
    }
    //initialize initial conditions for output
    for(j=0;j<bodies;j++){
        output[0][j*8+1] = system[j].position[0];
        output[0][j*8+2] = system[j].position[1];
        output[0][j*8+3] = system[j].position[2];
        output[0][j*8+4] = system[j].position[3];
        output[0][j*8+5] = system[j].velocity[0];
        output[0][j*8+6] = system[j].velocity[1];
        output[0][j*8+7] = system[j].velocity[2];
        output[0][j*8+8] = system[j].velocity[3];
    }

    //Begin RK4 loop
    for(i=0;i<steps;i++){
        output[i+1][0] = output[i][0] + h;
        //calculate first positions [1+1/2]
        // (using k1x)
        for(j=sun_choice;j<bodies;j++){
            output[i+1][j*8+1] = output[i][j*8+1] + halfh*output[i][j*8+5];
            output[i+1][j*8+2] = output[i][j*8+2] + halfh*output[i][j*8+6];
            output[i+1][j*8+3] = output[i][j*8+3] + halfh*output[i][j*8+7];
        }
        //new relative positions at [i+1/2]
        for(m=0;m<bodies;m++){
            for(n=0;n<bodies;n++){
                if(m!=n){
                    relcoord[m][n*4] = output[i+1][m*8+1] - output[i+1][n*8+1];
                    relcoord[m][n*4+1] = output[i+1][m*8+2] - output[i+1][n*8+2];
                    relcoord[m][n*4+2] = output[i+1][m*8+3] - output[i+1][n*8+3];
                    relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
                }
            }
        }
        //new forces with relative positions
        for(j=0;j<bodies*3;j++){
            force[1][j] = 0.0;
        }
        for(j=0;j<bodies;j++){
            for(k=0;k<bodies;k++){
                if(j!=k){
                    rcubed = pow(relcoord[j][k*4+3],3.0);
                    force[1][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
                    force[1][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
                    force[1][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
                }
            }
            force[1][j*3] *= fourpisq;
            force[1][j*3+1] *= fourpisq;
            force[1][j*3+2] *= fourpisq;
            //k2v
            k2[j*6+3] = force[1][j*3];
            k2[j*6+4] = force[1][j*3+1];
            k2[j*6+5] = force[1][j*3+2];
        }
        //k2x using i+1/2 positions
        for(j=0;j<bodies;j++){
            k2[j*6] = output[i][j*8+5] + halfh*force[1][j*3];
            k2[j*6+1] = output[i][j*8+6] + halfh*force[1][j*3+1];
            k2[j*6+2] = output[i][j*8+7] + halfh*force[1][j*3+2];
        }
        //new temp positions at i+1/2 with k2x
        for(j=sun_choice;j<bodies;j++){
            output[i+1][j*8+1] = output[i][j*8+1] + halfh*k2[j*6];
            output[i+1][j*8+2] = output[i][j*8+2] + halfh*k2[j*6+1];
            output[i+1][j*8+3] = output[i][j*8+3] + halfh*k2[j*6+2];
        }
        //new relative positions at [i+1/2] with k2x
        for(m=0;m<bodies;m++){
            for(n=0;n<bodies;n++){
                if(m!=n){
                    relcoord[m][n*4] = output[i+1][m*8+1] - output[i+1][n*8+1];
                    relcoord[m][n*4+1] = output[i+1][m*8+2] - output[i+1][n*8+2];
                    relcoord[m][n*4+2] = output[i+1][m*8+3] - output[i+1][n*8+3];
                    relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
                }
            }
        }
        //new forces with relative positions
        for(j=0;j<bodies*3;j++){
            force[1][j] = 0.0;
        }
        for(j=0;j<bodies;j++){
            for(k=0;k<bodies;k++){
                if(j!=k){
                    rcubed = pow(relcoord[j][k*4+3],3.0);
                    force[1][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
                    force[1][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
                    force[1][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
                }
            }
            force[1][j*3] *= fourpisq;
            force[1][j*3+1] *= fourpisq;
            force[1][j*3+2] *= fourpisq;
            //k3v
            k3[j*6+3] = force[1][j*3];
            k3[j*6+4] = force[1][j*3+1];
            k3[j*6+5] = force[1][j*3+2];
        }
        //k3x using i+1/2 positions
        for(j=0;j<bodies;j++){
            k3[j*6] = output[i][j*8+5] + halfh*force[1][j*3];
            k3[j*6+1] = output[i][j*8+6] + halfh*force[1][j*3+1];
            k3[j*6+2] = output[i][j*8+7] + halfh*force[1][j*3+2];
        }
        //new temp positions at i+1 using k3x
        for(j=sun_choice;j<bodies;j++){
            output[i+1][j*8+1] = output[i][j*8+1] + h*k3[j*6];
            output[i+1][j*8+2] = output[i][j*8+2] + h*k3[j*6+1];
            output[i+1][j*8+3] = output[i][j*8+3] + h*k3[j*6+2];
        }
        //new relative positions at [i+1] with k3x
        for(m=0;m<bodies;m++){
            for(n=0;n<bodies;n++){
                if(m!=n){
                    relcoord[m][n*4] = output[i+1][m*8+1] - output[i+1][n*8+1];
                    relcoord[m][n*4+1] = output[i+1][m*8+2] - output[i+1][n*8+2];
                    relcoord[m][n*4+2] = output[i+1][m*8+3] - output[i+1][n*8+3];
                    relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
                }
            }
        }
        //new forces with relative positions
        for(j=0;j<bodies*3;j++){
            force[1][j] = 0.0;
        }
        for(j=0;j<bodies;j++){
            for(k=0;k<bodies;k++){
                if(j!=k){
                    rcubed = pow(relcoord[j][k*4+3],3.0);
                    force[1][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
                    force[1][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
                    force[1][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
                }
            }
            force[1][j*3] *= fourpisq;
            force[1][j*3+1] *= fourpisq;
            force[1][j*3+2] *= fourpisq;
        }
        //k4x using i+1/2 positions
        for(j=sun_choice;j<bodies;j++){
            output[i+1][j*8+5] = output[i][j*8+5] + h*force[1][j*3];
            output[i+1][j*8+6] = output[i][j*8+6] + h*force[1][j*3+1];
            output[i+1][j*8+7] = output[i][j*8+7] + h*force[1][j*3+2];
        }
        //final positions at i+1
        for(j=sun_choice;j<bodies;j++){
            output[i+1][j*8+1] = output[i][j*8+1] + (h/6.0)*(output[i][j*8+5] + 2.0*k2[j*6] + 2.0*k3[j*6] + output[i+1][j*8+5]);
            output[i+1][j*8+2] = output[i][j*8+2] + (h/6.0)*(output[i][j*8+6] + 2.0*k2[j*6+1] + 2.0*k3[j*6+1] + output[i+1][j*8+6]);
            output[i+1][j*8+3] = output[i][j*8+3] + (h/6.0)*(output[i][j*8+7] + 2.0*k2[j*6+2] + 2.0*k3[j*6+2] + output[i+1][j*8+7]);
            output[i+1][j*8+4] = sqrt(output[i+1][j*8+1]*output[i+1][j*8+1] + output[i+1][j*8+2]*output[i+1][j*8+2] + output[i+1][j*8+3]*output[i+1][j*8+3]);
        }
        //final velocities at i+1
        for(j=sun_choice;j<bodies;j++){
            output[i+1][j*8+5] = output[i][j*8+5] + (h/6.0)*(force[0][j*3] + 2.0*k2[j*6+3] + 2.0*k3[j*6+3] + force[1][j*3]);
            output[i+1][j*8+6] = output[i][j*8+6] + (h/6.0)*(force[0][j*3+1] + 2.0*k2[j*6+4] + 2.0*k3[j*6+4] + force[1][j*3+1]);
            output[i+1][j*8+7] = output[i][j*8+7] + (h/6.0)*(force[0][j*3+2] + 2.0*k2[j*6+5] + 2.0*k3[j*6+5] + force[1][j*3+2]);
            output[i+1][j*8+8] = sqrt(output[i+1][j*8+5]*output[i+1][j*8+5] + output[i+1][j*8+6]*output[i+1][j*8+6] + output[i+1][j*8+7]*output[i+1][j*8+7]);
        }
        //final relative positions at i+1
        for(m=0;m<bodies;m++){
            for(n=0;n<bodies;n++){
                if(m!=n){
                    relcoord[m][n*4] = output[i+1][m*8+1] - output[i+1][n*8+1];
                    relcoord[m][n*4+1] = output[i+1][m*8+2] - output[i+1][n*8+2];
                    relcoord[m][n*4+2] = output[i+1][m*8+3] - output[i+1][n*8+3];
                    relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
                }
            }
        }
        //final forces with relative positions
        for(j=0;j<bodies*3;j++){
            force[1][j] = 0.0;
        }
        for(j=0;j<bodies;j++){
            for(k=0;k<bodies;k++){
                if(j!=k){
                    rcubed = pow(relcoord[j][k*4+3],3.0);
                    force[1][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
                    force[1][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
                    force[1][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
                }
            }
            force[1][j*3] *= fourpisq;
            force[1][j*3+1] *= fourpisq;
            force[1][j*3+2] *= fourpisq;
        }
        for(j=0;j<bodies;j++){
            force[0][j*3] = force[1][j*3];
            force[0][j*3+1] = force[1][j*3+1];
            force[0][j*3+2] = force[1][j*3+2];
        }

    }

    //delete dross--keep output
    for(i=0;i<bodies;i++){
        delete[] relcoord[i];
    }
    delete[] relcoord;
    for(i=0;i<2;i++){
        delete[] force[i];
    }
    delete[] force;
    delete[] k2;
    delete[] k3;
}

void solarsystem::verlet(double**& output, int steps, double tmax, int sun_choice){
    int i, j, k, m, n, bodies;
    bodies = body_count;
    double h, halfh, halfhsq, rcubed, fourpisq;
    h = (tmax-0.0)/((double) (steps));
    halfh = 0.5*h;
    halfhsq = halfh*h;
    fourpisq = 4.0*M_PI*M_PI;
    //output to provide t and x, y, z, vx, vy, vz for each body at each step i
    output = new double*[steps+1];
    for(int i=0;i<steps+1;i++){
        output[i] = new double[bodies*8+1];
    }
    //relcoord to provide relative coordinates, x, y, z, r, of body i with respect to body j
    //technically double counting--consider revising
    double** relcoord = new double*[bodies];
    for(i=0;i<bodies;i++){
        relcoord[i] = new double[bodies*4];
    }
    //force provides fx,fy,fz for each body at step i and i+1
    double** force = new double*[2];
    for(i=0;i<2;i++){
        force[i] = new double[bodies*3];
    }
    //initialize all matrices to zero
    for(i=0;i<steps+1;i++){
        for(j=0;j<bodies*8+1;j++){
            output[i][j]=0.0;
        }
    }
    for(i=0;i<bodies;i++){
        for(j=0;j<bodies*4;j++){
            relcoord[i][j]=0.0;
        }
    }
    for(i=0;i<2;i++){
        for(j=0;j<bodies*3;j++){
            force[i][j]=0.0;
        }
    }

    //reinitialize relcoord to initial conditions
    for(i=0;i<bodies;i++){
        for(j=0;j<bodies;j++){
            if(i!=j){
                relcoord[i][j*4] = system[i].position[0] - system[j].position[0];
                relcoord[i][j*4+1] = system[i].position[1] - system[j].position[1];
                relcoord[i][j*4+2] = system[i].position[2] - system[j].position[2];
                relcoord[i][j*4+3] = sqrt(relcoord[i][j*4]*relcoord[i][j*4] + relcoord[i][j*4+1]*relcoord[i][j*4+1]+relcoord[i][j*4+2]*relcoord[i][j*4+2]);
            }
        }
    }
    //initialize forces
    for(j=0;j<bodies;j++){
        for(k=0;k<bodies;k++){
            if(j!=k){
                rcubed = pow(relcoord[j][k*4+3],3.0);
                force[0][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
                force[0][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
                force[0][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
            }
        }
        force[0][j*3] *= fourpisq;
        force[0][j*3+1] *= fourpisq;
        force[0][j*3+2] *= fourpisq;
    }
    //initialize initial conditions for output
    for(j=0;j<bodies;j++){
        output[0][j*8+1] = system[j].position[0];
        output[0][j*8+2] = system[j].position[1];
        output[0][j*8+3] = system[j].position[2];
        output[0][j*8+4] = system[j].position[3];
        output[0][j*8+5] = system[j].velocity[0];
        output[0][j*8+6] = system[j].velocity[1];
        output[0][j*8+7] = system[j].velocity[2];
        output[0][j*8+8] = system[j].velocity[3];
    }

    //begin solution for-loop
    for(i=0;i<steps;i++){
        //Set next time value
        output[i+1][0] = output[i][0] + h;
        //calculate x,y,z
        for(j=sun_choice;j<bodies;j++){
            output[i+1][j*8+1] = output[i][j*8+1] + h*output[i][j*8+5] + halfhsq*force[0][j*3];
            output[i+1][j*8+2] = output[i][j*8+2] + h*output[i][j*8+6] + halfhsq*force[0][j*3+1];
            output[i+1][j*8+3] = output[i][j*8+3] + h*output[i][j*8+7] + halfhsq*force[0][j*3+2];
            output[i+1][j*8+4] = sqrt(output[i][j*8+1]*output[i][j*8+1] + output[i][j*8+2]*output[i][j*8+2] + output[i][j*8+3]*output[i][j*8+3]);
        }
        //new relcoord
        for(m=0;m<bodies;m++){
            for(n=0;n<bodies;n++){
                if(m!=n){
                    relcoord[m][n*4] = output[i+1][m*8+1] - output[i+1][n*8+1];
                    relcoord[m][n*4+1] = output[i+1][m*8+2] - output[i+1][n*8+2];
                    relcoord[m][n*4+2] = output[i+1][m*8+3] - output[i+1][n*8+3];
                    relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
                }
            }
        }
        //calculate force_i+1
        for(j=0;j<bodies*3;j++){
            force[1][j] = 0.0;
        }
        for(j=0;j<bodies;j++){
            for(k=0;k<bodies;k++){
                if(j!=k){
                    rcubed = pow(relcoord[j][k*4+3],3.0);
                    force[1][j*3] += -system[k].mass*relcoord[j][k*4]/rcubed;
                    force[1][j*3+1] += -system[k].mass*relcoord[j][k*4+1]/rcubed;
                    force[1][j*3+2] += -system[k].mass*relcoord[j][k*4+2]/rcubed;
                }
            }
            force[1][j*3] *= fourpisq;
            force[1][j*3+1] *= fourpisq;
            force[1][j*3+2] *= fourpisq;
        }
        //calculate vx, vy, vz
        for(j=sun_choice;j<bodies;j++){
            output[i+1][j*8+5] = output[i][j*8+5] + halfh*(force[1][j*3] + force[0][j*3]);
            output[i+1][j*8+6] = output[i][j*8+6] + halfh*(force[1][j*3+1] + force[0][j*3+1]);
            output[i+1][j*8+7] = output[i][j*8+7] + halfh*(force[1][j*3+2] + force[0][j*3+2]);
            output[i+1][j*8+8] = sqrt(output[i+1][j*8+5]*output[i+1][j*8+5] + output[i+1][j*8+6]*output[i+1][j*8+6] + output[i+1][j*8+7]*output[i+1][j*8+7]);
        }
        //force_i for next loop is force_i+1 for this loop. No need for double computation.
        for(j=0;j<bodies;j++){
            force[0][j*3] = force[1][j*3];
            force[0][j*3+1] = force[1][j*3+1];
            force[0][j*3+2] = force[1][j*3+2];
        }
    }
    //delete dross--keep output
    for(i=0;i<bodies;i++){
        delete[] relcoord[i];
    }
    delete[] relcoord;
    for(i=0;i<2;i++){
        delete[] force[i];
    }
    delete[] force;

}


