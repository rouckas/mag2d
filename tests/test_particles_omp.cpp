/*
 * Performance test of particle code parallelization
 * Intel(R) Core(TM) i5-2500K CPU @ 3.30GHz
 * 2 times 4 GB 1333 MHz memory modules
 * 100 iterations with 1e6 particles
 *
 * With Particle class containing 7 doubles a 1 bool sizeof(Particle)==64:
 * t_omp ~ 0.81 s       t_noomp ~ 1.49 s
 * speedup ~ 1.8 times with 4 threads. Data rates needed are
 * 100*1e6*64B/0.81s ~ 8 GB/s - we are not far from the theoretical
 * 21 GB/s bandwidth of core i5 processor
 *
 * Probably unsurprisingly, when we remove the sin calculation,
 * the advantage of parallelization is gone and both versions
 * run at t_omp ~ t_noomp ~ 0.78 s which is a further indication
 * of memory limitations
 *
 * With Particle class containing 3 doubles sizeof(Particle)==24:
 * t_omp ~ 0.45 s       t_noomp ~ 1.36 s
 * speedup ~ 3.0 times with threads.
 *
 *
 * 1000 iterations with 1e5 particles
 * With Particle class containing 7 doubles a 1 bool sizeof(Particle)==64:
 * t_omp ~ 0.53 s       t_noomp ~ 1.46 s
 * speedup ~ 2.7 times
 */
#include <cmath>
#include <vector>
#include <iostream>
#include "../src/timer.hpp"
using namespace std;

class Field
{
        double value;
    public:
        Field(double _value) : value(_value) {};
        double get_value() { return value; };
        inline double calculate(double x) { return value/(x);};
};

class Particle
{
    public:
        double x, y, z;

        //this can be commented out
        double vx,vy,vz;
        double time_to_death;
        bool empty;

        Particle() : x(0), y(0), z(0) {};
};

int main()
{
    double j = 0;
    int nparticles = 1000000;
    int niter = 100;
    t_timer timer;
    timer.start();

    cout << "sizes " << sizeof(double) <<" "<< sizeof(bool) <<" "<< sizeof(Particle) <<endl;

    {
        vector<Particle> particles(nparticles);
        for(int k=0; k<niter; k++)
        {
#pragma omp parallel for reduction(+:j)
            for(vector<Particle>::iterator I = particles.begin(); I < particles.end(); ++I)
            {
                I->z = I->y = sin(I->x+0.5);
                j += I->z + I->y + I->x;
            }
        }
    }
    timer.stop();
    cout << timer.get_real_time() <<" "<< timer.get_cpu_time() << endl;
    timer.reset();
    timer.start();
    

    {
        Particle * particles = new Particle[nparticles];
        for(int k=0; k<niter; k++)
        {
#pragma omp parallel for reduction(+:j)
            for(int i = 0; i < nparticles; i++)
            {
                particles[i].z = particles[i].y = sin(particles[i].x+0.5);
                j += particles[i].z + particles[i].y + particles[i].x;
            }
        }
    }
    
    timer.stop();
    cout << timer.get_real_time() <<" "<< timer.get_cpu_time() << endl;
    cout << j << endl;

    return 0;
}
