#ifndef PARAM_H
#define PARAM_H
#include <GetPot>
#include <string>
#include <stdexcept>
#include "speclist.hpp"
using namespace std;
enum t_advancer { ADVANCE2, ADVANCE_OLD, ADVANCE_PROB, ADVANCE2_EULER };
enum Coord { CARTESIAN, CYLINDRICAL };
class t_spec_param
{
    public:
	double density;
	double energy;
};

class Param
{
    public:
    enum Boundary { FREE, PERIODIC, MIRROR };
    double r_max, z_max;
    int r_sampl, z_sampl;
    double extern_field;
    std::string magnetic_field_file;
    bool magnetic_field_const;
    double Br, Bz, Bt;
    double probe_radius;
    double probe_length;
    double u_probe;
//    int n_particles;
    double n_particles_total;
    double dr, dz, dth;
    double V, dV;
    double idr, idz;
    const static double eps_0 = 8.854187817e-12;
    const static double k_B = 1.380662e-23;
    const static double q_e = 1.602189e-19;
//    double rho;			//[m-3] electron density
    double pressure;		//[Pa]
    double neutral_temperature;	//[K]
    double macroparticle_factor;
    double density[NTYPES];
    double dt[NTYPES];
    double temperature[NTYPES];
    double dt_elon;
    int niter;
    t_advancer advancer;
    Coord coord;
    Boundary boundary;
    int src_fact;
    bool selfconsistent;
    bool rf;
    double rf_omega;
    bool particle_reload;
    int t_print;
    int t_equilib;
    std::string particle_reload_dir;
    std::string output_dir;
    bool do_plot;

    double neutral_density;
    /*
    t_param(double _x_max, double _y_max, int _x_sampl, int _y_sampl, int _n_particles, double _rho,
	    double _pressure=133.0, double _neutral_temperature=300.0) 
	: x_max(_x_max), y_max(_y_max), x_sampl(_x_sampl), y_sampl(_y_sampl), n_particles(_n_particles), rho(_rho),
	pressure(_pressure), neutral_temperature(_neutral_temperature)
    {
	n_particles_total = n_particles*2.0;
	neutral_density = pressure/(k_B*neutral_temperature);

	dx = x_max/(x_sampl-1);
	dy = y_max/(y_sampl-1);
	idx = 1.0/dx;
	idy = 1.0/dy;
	V = n_particles/(rho);
	dV = V/((x_sampl-1.0)*(y_sampl-1.0));
	dz = dV/(dx*dy);
    };
    */
    Param(GetPot & config)
    {
	r_max = config("r_max",1e-2);
	z_max = config("z_max",1e-2);
	r_sampl = config("r_sampl",100);
	z_sampl = config("z_sampl",100);

//	n_particles = config("n_particles",100000);
	n_particles_total = config("n_particles_total",100000);
//	rho = config("rho",1e15);		// XXX rho must correspond to n_particles
	pressure = config("pressure",133.0);
	neutral_temperature = config("neutral_temperature",300.0);

	probe_radius = config("probe_radius",1e-4);
	probe_length = config("probe_length",1e-2);
	u_probe = config("u_probe",-10.0);
	extern_field = config("extern_field",100.0);

        magnetic_field_file = config("magnetic_field_file", "");
        magnetic_field_const = config("magnetic_field_const", 1);
        Br = config("Br", 0.0);
        Bz = config("Bz", 0.0);
        Bt = config("Bt", 0.0);

	niter = config("niter",100000);
	dt_elon = config("dt_elon",1e-11);

	selfconsistent = config("selfconsistent",1);
        rf = config("rf", 1);
        rf_omega = config("rf_omega", 1e8);
        if(selfconsistent && rf)
        {
            printf("error: selfconsistent rf trap not implemented\n");
            exit(1);
        }
	t_print = config("t_print",0);
	t_equilib = config("t_equilib",niter+1);
	particle_reload = config("particle_reload",0);
	particle_reload_dir = config("particle_reload_dir",".");

	src_fact = config("src_fact",20);

	neutral_density = pressure/(k_B*neutral_temperature);
	do_plot = config("do_plot",1);

        string coord_str = config("coord","CYLINDRICAL");
        if(coord_str == "CYLINDRICAL") { coord = CYLINDRICAL; }
        else { coord = CARTESIAN; }

        string boundary_str = config("boundary","FREE");
        if(boundary_str == "FREE") { boundary = FREE; }
        else if(boundary_str == "MIRROR") { boundary = MIRROR; }
        else if(boundary_str == "PERIODIC") { boundary = PERIODIC; }
        else { throw std::runtime_error("Param: unrecognized boundary value " + boundary_str + "\n"); }
        if(coord == CYLINDRICAL && boundary != FREE)
            throw std::runtime_error("Param: only FREE boundary condition in cylindrical coords is implemented\n");
        if(boundary == MIRROR)
            throw std::runtime_error("Param: MIRROR boundary condition not implemented\n");

	//parse plasma parameters
	double particle_density_total = 0.0;
	for(int i=0; i<NTYPES; i++)
	{
	    string spec_name = species_names[i];
	    config.set_prefix( (spec_name+"/").c_str() );
	    
	    temperature[i] = config("temperature",-1.0);
	    if(temperature[i] == -1.0)
	    {
		double tmp = config("energy",0.0);
		temperature[i] = tmp*q_e/k_B*2.0/3.0;
	    }

	    density[i] = config("density",0.0);
	    if(density[i] == 0.0)
	    {
		double tmp = config("pressure",0.0);
		if(temperature[i]!=0.0)
		    density[i] = tmp/(k_B*temperature[i]);
	    }
	    if(is_particle[i]) particle_density_total += density[i];

	    dt[i] = config("dt",1e-11);

	    cout << spec_name << ' ' << density[i] << ' ' << temperature[i] <<endl;
	}

	dr = r_max/(r_sampl-1);
	dz = z_max/(z_sampl-1);
	idr = 1.0/dr;
	idz = 1.0/dz;
	//V = n_particles/(rho);
        config.set_prefix("");
	macroparticle_factor = config("macroparticle_factor",1e4);
        cout << "macroparticle_factor " << macroparticle_factor <<endl;
	V = n_particles_total/particle_density_total;
        // XXX check this, only in cartesian coords:
        if(selfconsistent)
            dV = 1;
        else
            dV = V/((r_sampl-1)*(z_sampl-1));
	cout << "param: V = "<<V<<endl;
	dth = 2*M_PI/macroparticle_factor;
    }
};
const double Param::eps_0;
const double Param::k_B;
const double Param::q_e;
#endif
