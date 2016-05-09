#ifndef PROJECT3_LIB_H
#define PROJECT3_LIB_H

#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>

using namespace std;

inline void matrix_alloc(double**&, int, int);
inline void matrix_delete(double**&, int);
inline void matrix_resize(double**&, int, int, int, int);
void verlet(double**&, double*&, int, int, double, int);
void verlet_relcor(double**&, double*&, int, int, double, int);
void helionstates(double**&, double**&, double**&, int, int);
void helionstates_dynamic(double**& aphelion, int& aphguess, double**& perihelion, int& periguess, double**& output, int steps, int body, int range, int dec_or_inc);


/***************
Matrix functions
***************/
inline void matrix_alloc(double**& a, int rows, int columns){
    int i, j;
    a = new double*[rows];
    for(i=0;i<rows;i++){
        a[i] = new double[columns];
    }
    for(i=0;i<rows;i++){
        for(j=0;j<columns;j++){
            a[i][j] = 0.0;
        }
    }
}
inline void matrix_delete(double**& a, int rows){
    for(int i=0;i<rows;i++){
        delete a[i];
    }
    delete[] a;
}
inline void matrix_resize(double**& matrix, int oldrows, int oldcol, int newrows, int newcol){
    int i, j;
    double** temp = new double*[newrows];
    for(i=0;i<newrows;i++){
        temp[i] = new double[newcol];
    }
    for(i=0;i<oldrows;i++){
        for(j=0;j<oldcol;j++){
            temp[i][j]=matrix[i][j];
        }
    }

    for(i=0;i<oldrows;i++){
        delete[] matrix[i];
    }
    delete[] matrix;
    matrix = temp;
}

/***************************
helionstates functions to
retrieve aphelion & peri-
helion data
***************************/
void helionstates(double**& output, double**& aphelion, double**& perihelion, int bodies, int steps){
    int i, j;
    //Initialize with beginning values
    for(i=0;i<bodies;i++){
        aphelion[i][0] = output[0][i*8+4];
        perihelion[i][0] = aphelion[i][0];
        aphelion[i][1] = output[0][0];
        perihelion[i][1] = aphelion[i][1];
    }
    //sort through all radii and find the largest and smallest for each body
    for(i=1;i<steps+1;i++){
        for(j=0;j<bodies;j++){
            if(output[i][j*8+4]>aphelion[j][0]){
                aphelion[j][0] = output[i][j*8+4];
                aphelion[j][1] = output[i][0];
            }
            if(output[i][j*8+4]<perihelion[j][0]){
                perihelion[j][0] = output[i][j*8+4];
                perihelion[j][1] = output[i][0];
            }
        }
    }
}

void helionstates_dynamic(double**& aphelion, int& aphguess, double**& perihelion, int& periguess, double**& output, int steps, int body, int range, int dec_or_inc){
    int i, j, stepcount, aphcount, pericount, rposition, holder0, holder1;
    double maxtemp[2], mintemp[2], holder[2];
    rposition = body*8+4;

    for(i=1;i<5;i++){
        aphelion[0][i] = output[0][body*8+i];
        perihelion[0][i] = aphelion[0][i];
    }

    aphelion[0][0] = output[0][0];
    perihelion[0][0] = aphelion[0][0];
    holder[0] = aphelion[0][4];
    holder[1] = 0;

    stepcount = 0;
    aphcount=1;
    pericount=1;

    while(stepcount<steps+1){
        if(aphcount >= aphguess){
            matrix_resize(aphelion, aphguess, 5, aphcount+1, 5);
            aphelion[aphcount][0] = -1;
            for(i=0;i<5;i++){
                aphelion[aphcount][i] = 0;
            }
            aphguess = aphcount+1;
        }if(pericount >= periguess){
            matrix_resize(perihelion, periguess, 5, pericount+1, 5);
            perihelion[pericount][0] = -1;
            for(i=1;i<5;i++){
                perihelion[pericount][i] = 0;
            }
            periguess = pericount+1;
        }

        if(dec_or_inc==0){
            mintemp[0] = holder[0];
            mintemp[1] = holder[1];
            for(i=stepcount;i<stepcount+range;i++){
                if(output[i][rposition]<=mintemp[0]){
                    mintemp[0] = output[i][rposition];
                    mintemp[1] = i;
                }
            }
            if(mintemp[0]<holder[0]){
                holder[0] = mintemp[0];
                holder[1] = mintemp[1];
            }
            else{
                dec_or_inc = 1;
                holder1=holder[1];
                perihelion[pericount][0] = output[holder1][0];
                for(i=1;i<5;i++){
                    perihelion[pericount][i] = output[holder1][body*8+i];
                }
                pericount++;
            }
        }
        else if(dec_or_inc==1){
            maxtemp[0] = holder[0];
            mintemp[1] = holder[1];
            for(i=stepcount;i<stepcount+range;i++){
                if(output[i][rposition]>=maxtemp[0]){
                    maxtemp[0] = output[i][rposition];
                    maxtemp[1] = i;
                }
            }
            if(maxtemp[0]>holder[0]){
                holder[0] = maxtemp[0];
                holder[1] = maxtemp[1];
            }
            else{
                dec_or_inc = 0;
                holder1=holder[1];
                aphelion[aphcount][0] = output[holder1][0];
                for(i=1;i<5;i++){
                    aphelion[aphcount][i] = output[holder1][body*8+i];
                }
                aphcount++;
            }
        }
        stepcount += range;
        if(stepcount + range > steps+1){
            range = (steps+1) - stepcount;
        }
    }
}

#endif // PROJECT3_LIB_H
