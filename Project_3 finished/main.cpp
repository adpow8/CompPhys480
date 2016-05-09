#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include "time.h"
#include "solarsystem.h"
#include "body.h"
#include "project3_lib.h"

using namespace std;

ofstream ofile, ofile2;

int main(int argc, char* argv[])
{
    double** output, **aphelion, **perihelion;
    char* outfilename;

    int i, j, k, steps, system_choice, sun_choice, method, outputtype, bodies;
    int stepsperyear, timesperyear;
    int dec_or_inc, aphguess, periguess;
    double tmax, time, timesort;

    clock_t start, finish;
    string periname, aphname;
    stringstream number;

    if(argc<8)
    {
        cout << "Bad usage. Enter also 'outfilename system_choice immobilize_sun_choice method_choice output_type steps tmax' on same line." << endl;
        cout << "Current options for system_choice include" << endl;
        cout << "\t1: coplanar earth-sun system" << endl << "\t2: coplanar earth-jupiter-sun system" << endl;
        cout << "\t3: full solar system, with planets and luna" << endl;
        cout << "Current options for immobilze_sun_choice include" << endl;
        cout << "\t0: allow sun to evolve in coordinate space" << endl << "\t1: sun is initialized only and is not affected by other bodies" << endl;
        cout << "Current options for method_choice include" << endl;
        cout << "\t1: RK4" << endl << "\t2: Verlet" << endl;
        cout << "Current options for output_type include" << endl;
        cout << "\t1: pure output, steps by 8*bodies + 1" << endl << "\t2: as one, with print maximum aphelion and minimum perihelion to screen" << endl;
        cout << "\t3: 'timesperyear' ouput per earth year, enter after 'tmax'" << endl << "\t4: as three, with helionstates data files created" << endl;
        exit(1);
    }

    else
    {
        outfilename = argv[1];
        system_choice = atoi(argv[2]);
        sun_choice = atoi(argv[3]);
        method = atoi(argv[4]);
        outputtype = atoi(argv[5]);
        steps = atoi(argv[6]);
        tmax = atof(argv[7]);

        if(outputtype==3||outputtype==4)
        {
            if(argc<9)
            {
                cout << "Bad usage. Enter also 'timesperyear' after 'tmax.'" << endl;
                exit(1);
            }
            tmax = atoi(argv[7]);
            timesperyear = atoi(argv[8]);
            if(timesperyear==0)
            {
                cout << "Bad usage. Need a non-zero 'timesperyear.'" << endl;
                exit(1);
            }
            if(timesperyear*tmax>steps)
            {
                cout << "Bad usage. 'timesperyear'*'tmax' is bigger than total 'steps.' Insufficient data points for output." << endl;
                exit(1);
            }
            stepsperyear = steps/(tmax*timesperyear);
        }
    }

    /********************************************
    Here are the major massive bodies present
    within our solar system. Comment out any
    which are not desired for the calculation.
    Current options are
    1: earth-sun system
    2: earth-jupiter-sun system
    3: full solar system, with planets and luna
    ********************************************/

    solarsystem solSystem;

    if(system_choice==1){
        body sun(1,0,0,0,0,0,0);
        body earth(3.003489e-6,1,0,0,0,2.0*M_PI,0);
        solSystem.add(sun);
        solSystem.add(earth);
    }
    if(system_choice==2){
        if(sun_choice==0){
            body sun(1,-3.003489e-6,9.5479194e-4*5.2,0,-9.5479194e-4*0.439*2.0*M_PI,-3.003489e-6*2.0*M_PI,0);
            body earth(3.003489e-6,1,0,0,0,2.0*M_PI,0);
            body jupiter(9.5479194e-4,0,-5.2,0,0.439*2.0*M_PI,0,0);
            solSystem.add(sun);
            solSystem.add(earth);
            solSystem.add(jupiter);
        }
        if(sun_choice==1){
            body sun(1,0,0,0,0,0,0);
            body earth(3.003489e-6,1,0,0,0,2.0*M_PI,0);
            body jupiter(9.5479194e-4,0,-5.2,0,0.439*2.0*M_PI,0,0);
            solSystem.add(sun);
            solSystem.add(earth);
            solSystem.add(jupiter);
        }
    }
    if(system_choice==3){
        body sun(1,3.771551748320805E-03, 1.938413234187417E-03, -1.625928558000791E-04,365*3.428095941073785E-08, 365*6.978886434512168E-06, 365*-9.372671992938156E-09);
        body mercury(1.6601e-7,3.566221110752382E-01,-1.449153604767920E-01,-4.453344798939488E-02, 365*5.324168987021533E-03,365*2.725352157519689E-02, 365*1.737877336062238E-03);
        body venus(2.4478383e-6,4.456189829066335E-01,-5.759198980424926E-01,-3.358283335789598E-02,365*1.593371764525070E-02, 365*1.222110108236829E-02,365*-7.520174121955455E-04);
        body earth(3.003489e-6,-9.906650404586314E-01,4.353612431574581E-02,-1.569714899841466E-04,365*-9.933032846586867E-04,365*-1.724423592380582E-02,365*2.880574493748607E-07);
        body moon(1.230004e-8, -9.917828800119005E-01, 4.587463373520608E-02, -3.563694480647205E-04, 365*-1.530892843498182E-03, 365*-1.746686609878291E-02, 365*2.803200805250705E-05);
        body mars(3.227151e-7, -1.394909885799664E+00, -7.759974369033253E-01, 1.786251125355763E-02, 365*7.324673346484570E-03, 365*-1.102624283521118E-02, 365*-4.109846883854566E-04);
        body jupiter(9.5479194e-4,-5.321136962878863E+00, 1.055810040982731E+00, 1.146123452783326E-01, 365*-1.556597232697437E-03, 365*-7.046207863842619E-03, 365*6.409351264039102E-05);
        body saturn(2.858860e-4,-3.332098484988519E+00, -9.442038663142483E+00, 2.967846224795282E-01, 365*4.954645896896637E-03, 365*-1.873255191158998E-03, 365*-1.647257228743735E-04);
        body uranus(4.366244e-5, 1.876841090026478E+01, 6.826065082612249E+00, -2.177966843356428E-01, 365*-1.372937790621792E-03, 365*3.512867479374193E-03, 365*3.086802162915850E-05);
        body neptune(5.151389e-5, 2.803900494548452E+01, -1.054870089186826E+01, -4.289565554838171E-01, 365*1.084650993604757E-03, 365*2.957157649376530E-03, 365*-8.562727609311126E-05);
        body pluto(6.5812e-9, 8.773368933896196E+00, -3.186331328356860E+01, 8.718065633574812E-01, 365*3.100891963853092E-03, 365*1.939401372093854E-04, 365*-9.194995916567601E-04);

        solSystem.add(sun);
        solSystem.add(mercury);
        solSystem.add(venus);
        solSystem.add(earth);
        solSystem.add(moon);
        solSystem.add(mars);
        solSystem.add(jupiter);
        solSystem.add(saturn);
        solSystem.add(uranus);
        solSystem.add(neptune);
        solSystem.add(pluto);
    }

    bodies = solSystem.body_count;
    matrix_alloc(output, steps+1, bodies*8+1);


    /**********************************************
    Methods
        1: RK4
        2: Verlet
    **********************************************/
    if(method==1){
        start = clock();
        solSystem.RK4(output,steps,tmax,sun_choice);
        finish = clock();
    }
    else if(method==2){
        start = clock();
        solSystem.verlet(output,steps,tmax,sun_choice);
        finish = clock();
    }
    time = (finish - start)/((double) CLOCKS_PER_SEC);


    /*****************************************************************
    Outputfile choice
        1: Pure data file
        2: As 1, but with the largest aphelion,
            and smallest perihelion for each body printed to screen
        3: Output 'timesperyear' data points per year, for long
            time spans
        4: As two, but with data files for each body containing all
            aphelion and perilion positions
    *****************************************************************/

    ofile.open(outfilename);
    ofile.precision(8);

    if(outputtype==1){
        start = clock();
        for(i=0;i<steps+1;i++){
            for(j=0;j<bodies*8+1;j++){
                ofile << setw(15) << output[i][j];
            }
            ofile << endl;
        }
        finish = clock();
    }

    else if(outputtype==2){
        start = clock();
        matrix_alloc(aphelion,bodies,2);
        matrix_alloc(perihelion,bodies,2);
        for(i=0;i<steps+1;i++){
            for(j=0;j<bodies*8+1;j++){
                ofile << setw(15) << output[i][j];
            }
            ofile << endl;
        }
        helionstates(output, aphelion, perihelion, bodies, steps);
        for(i=0;i<bodies;i++){
            cout << setw(15) << aphelion[i][0] << setw(12) << aphelion[i][1];
            cout << setw(30) << perihelion[i][0] << setw(12) << perihelion[i][1];
            cout << endl;
        }
        matrix_delete(aphelion,bodies);
        matrix_delete(perihelion,bodies);
        finish = clock();
    }

    else if(outputtype==3){
        start = clock();
        for(i=0;i<steps+1;i+=stepsperyear){
            for(j=0;j<bodies*8+1;j++){
                ofile << setw(15) << output[i][j];
            }
            ofile << endl;
        }
        finish = clock();
    }

    else if(outputtype==4){
        start = clock();
        for(i=0;i<steps+1;i+=stepsperyear){
            for(j=0;j<bodies*8+1;j++){
                ofile << setw(15) << output[i][j];
            }
            ofile << endl;
        }
        for(k=0;k<bodies;k++){
            periname = "perihelion_body_";
            aphname = "aphelion_body_";
            aphguess = 1;
            periguess = aphguess;
            number << k;
            periname += number.str() + ".dat";
            aphname += number.str() + ".dat";
            number.str("");
            matrix_alloc(aphelion,aphguess,5);
            matrix_alloc(perihelion,periguess,5);
            dec_or_inc=1;
            helionstates_dynamic(aphelion, aphguess, perihelion, periguess, output,steps,k,20,dec_or_inc);

            ofile2.open(periname);
            ofile2.precision(8);
            for(i=0;i<periguess;i++){
                for(j=0;j<5;j++){
                    ofile2 << setw(15) << perihelion[i][j];
                }
                ofile2 << endl;
            }
            ofile2.close();
            ofile2.open(aphname);
            ofile2.precision(8);
            for(i=0;i<aphguess;i++){
                for(j=0;j<5;j++){
                    ofile2 << setw(15) << aphelion[i][j];
                }
                ofile2 << endl;
            }
            ofile2.close();

            matrix_delete(aphelion,aphguess);
            matrix_delete(perihelion,periguess);
        }
        finish = clock();
    }
    timesort = (finish - start)/((double) CLOCKS_PER_SEC);

    //print time for calculation to screen
    cout << "\tcomputation time " << time << " seconds" << endl;
    cout << "\toutput time " << timesort << " seconds" << endl;
    ofile.close();

    matrix_delete(output,steps+1);
    return 0;
}
