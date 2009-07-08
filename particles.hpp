#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <string>
#include "fields.hpp"
#include "param.hpp"
#include "histogram.hpp"
#include "random.cpp"
#include "mymath.cpp"
#include "input.hpp"
#include "tabulate.cpp"
#include "speclist.hpp"
#include "util.cpp"

using namespace std;
//#define NTYPES 8
//typedef enum species_type { NONE, ELECTRON, ARGON, ARGON_POS, ARGON_META, O2, O2_POS, O2_NEG };
//const char * species_names[NTYPES] = {"NONE", "ELECTRON", "ARGON", "ARGON_POS", "ARGON_META", "O2", "O2_POS", "O2_NEG" };

class t_particle
{
    public:
    double x, y, z;
    double vx,vy,vz;
    double time_to_death;
    bool empty;
};



class BaseSpecies
//abstract class representing one species in the model
//part of the class with dimension-independent code
{
    private:
	vector<int> empty;
	double source_len;
	ofstream output;
    protected:
	Fields *field;
	int niter;
	BaseSpecies ** species_list;
	Param *p_param;
	t_random * rnd;
    public:
	vector<t_particle> particles;
	species_type type;
	string name;
	double mass;
	double charge;
	double lifetime;
	double temperature;
	double density;
	double v_max;
	double dt;
	double EeV(double v) { return 0.5*mass*v*v/(p_param->q_e); }
	double veV(double EeV) { return sqrt(EeV*(p_param->q_e)/mass*2.0); }
	void set_pressure(double pa){ density = pa/(p_param->k_B*temperature);}
	Field rho;
	void accumulate();



	// diagnostics
	Histogram energy_dist;
	void energy_dist_compute();
	Histogram source_energy_dist;
	void source_energy_dist_compute();
	Histogram probe_energy_dist;
	Histogram probe_angular_dist;
	Histogram probe_angular_normalized_dist;
	double probe_current;
        double probe_charge;
	double probe_current_avg;
	Field rhoAverage;
	int nsampl;
	double probe_current_sum;


	BaseSpecies(int n1, int n2, Param &param, t_random &_rnd, Fields *_field,
	       double mass, BaseSpecies * _species_list[], species_type _type = NONE );

	// generate randomly new particle velocity according to species' velocity
	// distribution
	void rndv(double & vr, double & vz, double & vt)
	{
	    vr = rnd->rnor()*v_max*M_SQRT1_2;
	    vz = rnd->rnor()*v_max*M_SQRT1_2;
	    vt = rnd->rnor()*v_max*M_SQRT1_2;
	    //alternatively we can take random sample from source particles
	    /*
	    int i = rnd->inuni() % source2_particles.size();
	    vr = source2_particles[i].vx;
	    vz = source2_particles[i].vz;
	    vz = source2_particles[i].vz;
	    */
	}

	void rndv(double & vr)
	{
	    vr = rnd->maxwell1(v_max);
	}


	//void save(const char *filename)
	void save(const string filename);
	void load(const string filename);
	void remove(int n)
	{
	    if(particles[n].empty == true)
		throw runtime_error("removing empty particle");
	    else
	    {
		particles[n].empty = true;
		empty.push_back(n);
	    }
	}
	void remove(vector<t_particle>::iterator I)
	{
	    remove(I-particles.begin());
	}
	int insert()
	{
	    if(empty.size()==0)
		throw runtime_error(name + ": particle array full, size = " + int2string(particles.size()));
	    int tmp = empty.back();
	    if(particles[tmp].empty == false)
		throw runtime_error("error inserting particle");
	    empty.pop_back();
	    particles[tmp].empty = false;
	    return tmp;
	}
	int n_particles(){ return particles.size() - empty.size();}
        void resize(unsigned int newsize);

	virtual void scatter(t_particle &particle) = 0;
	vector<t_particle> source2_particles;
	void source5_save(const string filename);
	void source5_load(const string filename);
	int source5_factor;
	virtual void lifetime_init();
    protected:
	void load_CS(const string & fname, vector<vec_interpolate*> & CS, vector<double> & Loss,
	       	const vector<string> & CSnames, const string & tag="", int ncols=3);
	double svmax_find(const vector<vec_interpolate*> & CS, double _vmax, int samples=1000);
	    // this routine assumes input CS in units 1e-16 cm^2, possible source of errors...
    public:
	void print_status( ostream & out = cout);
	void dist_sample();
	void dist_reset();
};

template <int D>
class Species : public BaseSpecies {};

template <>
class Species<CARTESIAN> : public BaseSpecies
{
    public:
        Species(int _n, int n2, Param &param, t_random &_rnd, Fields *_field,
                double _mass, BaseSpecies * _species_list[], species_type _type) :
            BaseSpecies(_n, n2, param, _rnd, _field, _mass, _species_list, _type) {};

	void advance();
        void accumulate();

        void add_particles_on_disk(int nparticles, double centerx, double centery, double radius);
        void add_particle_beam_on_disk(int nparticles, double centerx, double centery, double radius);

	void source();
	void source_old();
	void source5_refresh(unsigned int factor);
};

template <>
class Species<CYLINDRICAL> : public BaseSpecies
{
    public:
        Species(int _n, int n2, Param &param, t_random &_rnd, Fields *_field,
                double _mass, BaseSpecies * _species_list[], species_type _type) :
            BaseSpecies(_n, n2, param, _rnd, _field, _mass, _species_list, _type) {};

	void advance();
        void accumulate();
	void probe_collect(t_particle *I);

        void add_particle_beam_on_disk_cylindrical(int nparticles, double centerz, double radius);
        void add_monoenergetic_particles_on_cylinder_cylindrical(int nparticles, double energy, double centerz, double radius, double height = 0.0);

};
#endif
