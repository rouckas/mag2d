#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "fields.hpp"
#include "Field2D.hpp"
#include "param.hpp"
#include "histogram.hpp"
#include "random.cpp"
#include "mymath.cpp"
#include "input.hpp"
#include "tabulate.cpp"
#include "util.cpp"
#include "parser.hpp"

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

class BaseSpecies;

class TrackedParticle
{
    public:
        TrackedParticle(BaseSpecies *_species, unsigned int _index, Param &param);
        void print();
    private:
        BaseSpecies * species;
        unsigned int index;
        ofstream fw;
};


class Interaction
{
    public:
        string name;
        CollType type;
        double DE;
        double rate;
        double cutoff;
        BaseSpecies * primary;
        BaseSpecies * secondary;
        Param & param;
        vec_interpolate * cross_section_table;
        double sigma_v(double v)
        {
            double result;
            if(type==COULOMB)
                result = coulomb_sigma(EeV(v)) * v;
            else if(cross_section_table != NULL)
                result = (*cross_section_table)(EeV(v)) * v;
            else result = rate;
            return result;
        };
        double EeV(double v_rel);
        double E(double v_rel);
        double v_rel(double E);
        Interaction(InteractionParams * i_params, Param & _param) : name(i_params->name),
            type(i_params->type), DE(i_params->DE), rate(i_params->rate), cutoff(i_params->cutoff), param(_param), cross_section_table(NULL)
        {
            if(i_params->CS_energy.size() > 0)
            {
                cross_section_table = new vec_interpolate(i_params->CS_energy, i_params->CS_value);
            }
            if(type==LANGEVIN)
                rate *= SQR(cutoff);
        };
        ~Interaction()
        {
            if(cross_section_table != NULL)
                delete cross_section_table;
        }
    private:
        double coulomb_sigma(double E);

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
	unsigned long int niter;
        vector<TrackedParticle*> tracked_particles;
	Param * const p_param;
	t_random * const rnd;
    public:
        vector<BaseSpecies *> speclist;
        vector<Interaction *> interactions;
        vector< vector<Interaction *> > interactions_by_species;
        vector<double> rates_by_species;
	vector<t_particle> particles;
	SpeciesType type;
	string name;
        bool particle;
	double mass;
	double charge;
	double lifetime;
	double temperature;
        double polarizability;
	double density;
	double v_max;
	double dt;
        double t;
	double EeV(double v) { return 0.5*mass*v*v/(p_param->q_e); }
	double veV(double EeV) { return sqrt(EeV*(p_param->q_e)/mass*2.0); }
	void set_pressure(double pa){ density = pa/(p_param->k_B*temperature);}
	Field2D rho;
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
	Field2D rhoAverage;
	int nsampl;
	double probe_current_sum;


        BaseSpecies(SpeciesParams * params, Param &param, t_random &_rnd) :
            empty(0), niter(0),
            p_param(&param), rnd(&_rnd),
            particles(0),
            type(params->type),
            name(params->name),
            particle(true), //XXX
            mass(params->mass), charge(params->charge),
            lifetime(INFINITY),
            temperature(params->temperature),
            polarizability(params->polarizability),
            density(params->density),
            dt(params->dt),
            t(0),
            rho(param.x_sampl, param.z_sampl, param.dx, param.dz),
            energy_dist(100,0.0,temperature*param.k_B/param.q_e*10.0),
            source_energy_dist(100,0.0,temperature*param.k_B/param.q_e*10.0),
            probe_energy_dist(100,0.0,temperature*param.k_B/param.q_e*20.0),
            probe_angular_dist(30,0.0,M_PI*0.5),
            probe_angular_normalized_dist(30,0.0,M_PI*0.5),
            probe_current(0),
            probe_charge(0),
            rhoAverage(param.x_sampl, param.z_sampl, param.dx, param.dz),
            nsampl(0)

        {
            v_max = sqrt(2.0*physconst::k_B*temperature/mass);
            output.open( (param.output_dir + "/" + name + ".dat").c_str() );
            rho.reset();
        };


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

        t_particle * random_particle()
        {
            int i = rnd->iuni() % particles.size();
            while(particles[i].empty == true)
            {
                i = rnd->iuni() % particles.size();
            }

            return &(particles[i]);

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
                resize(int(particles.size()*1.1) + 1);
	    int tmp = empty.back();
	    if(particles[tmp].empty == false)
		throw runtime_error("error inserting particle");
	    empty.pop_back();
	    particles[tmp].empty = false;
	    return tmp;
	}
	int n_particles(){ return particles.size() - empty.size();}
        void resize(unsigned int newsize);

	vector<t_particle> source2_particles;
	void source5_save(const string filename);
	void source5_load(const string filename);
	int source5_factor;
	void lifetime_init(double Emax=5.0);
    protected:
	void load_CS(const string & fname, vector<vec_interpolate*> & CS, vector<double> & Loss,
	       	const vector<string> & CSnames, const string & tag="", int ncols=3);
	double svmax_find(const vector<Interaction *> & interactions, double _vmax, int samples=1000);
	    // this routine assumes input CS in units 1e-16 cm^2, possible source of errors...
	void scatter(t_particle &particle);
    public:
	void print_status( ostream & out = cout);
        void print_trace()
        {
            for(vector<TrackedParticle*>::iterator I = tracked_particles.begin(); I != tracked_particles.end(); ++I)
                (*I)->print();
        }
        void print_distribution();
	void dist_sample();
	void dist_reset();
};

template <int D>
class Species : public BaseSpecies {};

template <>
class Species<CARTESIAN> : public BaseSpecies
{
    public:
        Fields * field;
        Species(SpeciesParams * params, Param &param, t_random &rnd, Fields & _field):
            BaseSpecies(params, param, rnd), field(&_field) {};

	void advance();
	void advance_boris();
	void advance_multicoll();
        void accumulate();

        void add_particles_everywhere(int nparticles);
        void add_particles_bessel(int nparticles, double centerx, double centery, double radius);
        void add_particles_on_disk(int nparticles, double centerx, double centery, double radius);
        void add_particle_beam_on_disk(int nparticles, double centerx, double centery, double radius);
        void add_tracked_particle(double x, double y, double vx, double vy, double vz)
        {
            int ii = insert();
            t_particle *pp = &(particles[ii]);
            pp->x = x;
            pp->z = y;
            pp->y = 0;
            pp->vx = vx;
            pp->vy = vz;
            pp->vz = vy;
            TrackedParticle *ptp = new TrackedParticle(this, ii, *p_param);
            tracked_particles.push_back(ptp);
        }

	void source();
	void source_old();
	void source5_refresh(unsigned int factor);
};

template <>
class Species<CYLINDRICAL> : public BaseSpecies
{
    public:
        Fields * field;
        Species(SpeciesParams * params, Param &param, t_random &rnd, Fields & _field):
            BaseSpecies(params, param, rnd), field(&_field) {};

	void advance();
	void advance_boris();
        void accumulate();
	void probe_collect(t_particle *I);

        void add_particle_beam_on_disk_cylindrical(int nparticles, double centerz, double radius);
        void add_monoenergetic_particles_on_cylinder_cylindrical(int nparticles, double energy, double centerz, double radius, double height = 0.0);

};


inline double Interaction::EeV(double v_rel)
{
    double mu = primary->mass*secondary->mass/(primary->mass + secondary->mass);
    double E = 0.5*mu*v_rel*v_rel/param.q_e;
    return E;
};
inline double Interaction::E(double v_rel)
{
    double mu = primary->mass*secondary->mass/(primary->mass + secondary->mass);
    double E = 0.5*mu*v_rel*v_rel;
    return E;
};

inline double Interaction::v_rel(double E)
{
    double mu = primary->mass*secondary->mass/(primary->mass + secondary->mass);
    double v_rel = sqrt(2*E/mu);
    return v_rel;
};

#endif
