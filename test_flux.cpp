#include "random.h"
#include "histogram.hpp"
#include <list>
#include <cmath>
#include <iostream>
using namespace std;

//#include "particles.cpp"
const double x_max = 1e-6;	//[m]
const double dt = 1e-8;
const int n_particles = 10000; 
const int niter = 100000;

const double temperature = 300.0;	//[K]
const double rho = 1e15;	//[m-3]
const double k_B = 1.380662e-23;
const double mass = 6.68173e-26;
const double q_e = 1.602189e-19;


class t_particle
{
    public:
	double x,vx;
	double vy,vz;
};
bool outside(const t_particle & particle){ return particle.x > x_max; };
template <class T> T sqr(T & x){ return x*x;};

int main()
{
    list<t_particle> particles;
    t_histogram energy_dist(100,0.0,temperature*k_B/q_e*10.0);
	

    maxwell_flux_seed();
    maxwell1_seed();
    double V = n_particles/rho;
    double dydz = V/x_max;
    double v_max = sqrt(2.0*k_B*temperature/mass);
    double f2 = rho/(2*sqrt(M_PI))*v_max*dt*dydz;

    for(int i=0; i<niter; i++)
    {
	//SOURCE
	int pocet = (int)f2;
	pocet += (uniform()>f2-pocet) ? 0 : 1;
	for(int j=0; j<pocet; j++)
	{
	    t_particle tmpp;
	    tmpp.vy = normal2()*v_max/(M_SQRT2);
	    tmpp.vz = normal2()*v_max/(M_SQRT2);
	    tmpp.vx = maxwell_flux(v_max*5.0);
	    tmpp.x = -tmpp.vx*dt*uniform();
	    particles.push_back(tmpp);
	}

	for(list<t_particle>::iterator I = particles.begin(); I != particles.end(); ++I)
	    I->x += I->vx*dt;
	particles.remove_if(outside);

	for(list<t_particle>::iterator I = particles.begin(); I != particles.end(); I++)
	    energy_dist.add((sqr(I->vx)+sqr(I->vy)+sqr(I->vz))*mass*0.5/q_e);
	if(i%100==0)
	    cout << i << " " << particles.size() << " " << energy_dist.mean()*q_e/k_B*2.0/3.0 << endl;

    }
}

