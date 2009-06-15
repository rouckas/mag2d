#include "matrix.cpp"
#include "timer.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int main(int argc, char * argv[])
{
    const int size = 256;
    Array3D<double> array(size, size, size);

    t_timer timer;


    cout << "*** PREHEATING ***\n";
    timer.start();
    for(int i=0; i<size; i++)
        for(int j=0; j<size; j++)
            for(int k=0; k<size; k++)
                array[i][j][k] = 0;
    timer.stop();
    cout << timer.gettime() <<endl;
    timer.reset();

    timer.start();
    for(int i=0; i<size; i++)
        for(int j=0; j<size; j++)
            for(int k=0; k<size; k++)
                array[i][j][k] = 0;
    timer.stop();
    cout << timer.gettime() <<endl;
    timer.reset();


    cout << "\n*** TEST START ***\n";
    cout << "access with []: ";
    timer.start();
    for(int i=0; i<size; i++)
        for(int j=0; j<size; j++)
            for(int k=0; k<size; k++)
                array[i][j][k] = 0;
    timer.stop();
    cout << timer.gettime() <<endl;
    timer.reset();

    cout << "access with (): ";
    timer.start();
    for(int i=0; i<size; i++)
        for(int j=0; j<size; j++)
            for(int k=0; k<size; k++)
                array(i, j, k) = 0;
    timer.stop();
    cout << timer.gettime() <<endl;
    timer.reset();

    double t1, t2;
    cout << "\n*** REVERSED INDICES (inhibiting cache) ***\n";
    cout << "[]\t()\tratio\n";
    double ratios[size];
    int niter = 10;
    for(int fn=0; fn<niter; fn++)
    {
        timer.start();
        for(int i=0; i<size; i++)
            for(int j=0; j<size; j++)
                for(int k=0; k<size; k++)
                    array[k][j][i] = 0;
        timer.stop();
        cout << (t1 = timer.gettime()) <<"\t";
        timer.reset();

        timer.start();
        for(int i=0; i<size; i++)
            for(int j=0; j<size; j++)
                for(int k=0; k<size; k++)
                    array(k, j, i) = 0;
        timer.stop();
        cout << (t2 = timer.gettime()) <<"\t";
        timer.reset();

        cout <<  t1/t2 <<endl;
        ratios[fn] = t1/t2;
    }
    double avg, stddev;
    avg=0;
    stddev=0;
    for(int i=0; i<niter; i++)
    {
        avg += ratios[i]/niter;
    }
    for(int i=0; i<niter; i++)
    {
        stddev += ratios[i]*ratios[i]/niter;
    }
    stddev = sqrt(stddev-avg*avg);
    cout << "operator () is by "<< setprecision(2) << (avg-1.0)*100 << " +/- " << setprecision(1) << stddev*100.0 << " % faster\n";
    return 0;
}
