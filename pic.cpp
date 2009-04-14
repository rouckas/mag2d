//#include "random.h"
#define GNUPLOT
#include "fields.cpp"
#include "param.cpp"
#ifdef GNUPLOT
#include "gnuplot_i.h"
#endif
#include "argonO2.cpp"
#include "hydrogen.cpp"
#include "helium.cpp"
//#include "oxygen.cpp"
#include "elon.cpp"
//#include "elon.cpp"
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include "random.cpp"
#include "timer.hpp"
using namespace std;

#define NPARTICL_SAFE 0.5
#define SRC_FACT 20

struct rusage picusage2,picusage1;

Species * make_species(species_type type, int n1, int n2, Param &param, t_random &rnd, Fields &field, Species *species_list[])
{
    switch(type)
    {
	case ARGON:
	    return new t_argon_neutral(n1, n2, 0, 0, param, rnd, &field, species_list, 0);
	case ARGON_POS:
	    return new t_argon(n1, n2, 0, 0, param, rnd, &field, species_list, 0);
	case ELECTRON:
	    return new t_elon(n1, n2, 0, 0, param, rnd, &field, species_list, 0);
	case HYDROGEN:
	    return new t_hydrogen_neutral(n1, n2, 0, 0, param, rnd, &field, species_list, 0);
	case H_NEG:
	    return new t_hydrogen_neg(n1, n2, 0, 0, param, rnd, &field, species_list, 0);
	case HELIUM:
	    return new t_helium_neutral(n1, n2, 0, 0, param, rnd, &field, species_list, 0);
	    /*
	case O2:
	    return new t_O2_neutral(n1, n2, 0, 0, param, rnd, &field, species_list, 0);
	case O2_POS:
	    return new t_O2_pos(n1, n2, 0, 0, param, rnd, &field, species_list, 0);
	case O_NEG:
	    return new t_O_neg(n1, n2, 0, 0, param, rnd, &field, species_list, 0);
	    */
	default:
	    cerr << "make_species(): Warning: unrecognized type\n";
	    return NULL;
    }
}
    

class Pic
{
    private:
	t_timer timer;
    public:
	Param & param;
	Species *species_list[NTYPES];
	Fields field;
	t_random rnd;
	int iter;
	Pic(Param & _param) : param(_param), field(param), iter(0)
	{
	    field.reset();

	    for(int i=0; i<NTYPES; i++)
		species_list[i] = NULL;

	    double total_density = 0;
	    for(int i=0; i<NTYPES; i++)
		if(is_particle[i])
		    total_density += param.density[i];

	    species_list[ARGON] = make_species(ARGON,0 ,0 , param, rnd, field, species_list);
	    species_list[HYDROGEN] = make_species(HYDROGEN,0 ,0 , param, rnd, field, species_list);
	    species_list[HELIUM] = make_species(HELIUM,0 ,0 , param, rnd, field, species_list);
	    //total_density = param.density[ARGON_POS] + param.density[ELECTRON] + param.density[O2_POS];

	    //nparticles_spec = int((param.density[ELECTRON]/total_density)*param.n_particles_total);
	    //species_list[ELECTRON] = make_species(ELECTRON, nparticles_spec*(1+NPARTICL_SAFE),nparticles_spec, param, rnd, field, species_list);
	    species_list[ELECTRON] = make_species(ELECTRON,100000,2, param, rnd, field, species_list);
	    species_list[ELECTRON]->source5_refresh(param.src_fact);

	    species_list[H_NEG] = make_species(H_NEG,2000,2, param, rnd, field, species_list);

	    species_list[ELECTRON]->lifetime_init();
	    species_list[H_NEG]->lifetime_init();

	    if(param.particle_reload)
		for(int ii=0; ii<NTYPES; ii++)
		    if( is_particle[ii] && species_list[ii] != NULL )
		    {
			species_list[ii]->load(param.particle_reload_dir + "/particles_" + species_list[ii]->name + ".dat");
			species_list[ii]->source5_load(param.particle_reload_dir + "/particles_source_" + species_list[ii]->name + ".dat");
		    }

	    dist_reset();
	    if(param.do_plot) plot_init();
	    if(!param.selfconsistent)
	    {
		field.boundary_solve();
		field.reset();
	    }
	};
	~Pic()
	{
	    for(int i=0; i<NTYPES; i++) if(species_list[i] != 0)
		delete species_list[i];
	    if(param.do_plot) plot_destroy();
	}
	void advance()
	{
	    timer.start();
	    if(param.selfconsistent)
	    {
		field.boundary_solve();
		field.reset();
		for(int i=0; i<NTYPES; i++) if(is_particle[i] && species_list[i] != NULL)
		{
		    species_list[i]->rho.reset();
		}
	    }

	    for(int i=0; i<NTYPES; i++) if(species_list[i] != 0)
	    {
		species_list[i]->advance();
		//if(param.selfconsistent)
		  //  species_list[i]->source();
	    }

	    for(int i=0; i<NTYPES; i++)
		if( is_particle[i] && species_list[i] != NULL && param.selfconsistent)
		{
		    field.rho.add(species_list[i]->rho);
		}
	    timer.stop();

	    iter++;
	}
        void emit()
        {
            double emission_current = 1e-4;     //[mA]
            double fnemit = 1.0;
            fnemit = param.dt[ELECTRON]/param.q_e * emission_current/param.macroparticle_factor;
            int nemit = (int)fnemit;
            nemit += (rnd.uni()>fnemit-nemit) ? 0 : 1;
            for(int i=0; i<nemit; i++)
            {
                int ind = species_list[ELECTRON]->insert();
                t_particle * p_p = &(species_list[ELECTRON]->particles[ind]);
                p_p->r = 2e-3*sqrt(rnd.uni());
                p_p->z = 0.5e-2;
                p_p->vr = species_list[ELECTRON]->veV(1);
                p_p->vz = species_list[ELECTRON]->veV(1);
                species_list[ELECTRON]->rndv(p_p->vr, p_p->vz, p_p->vt);
                //p_p->vr = 0;
                //p_p->vz = -500000;
                //p_p->vt = 0;
                p_p->time_to_death = species_list[ELECTRON]->lifetime * rnd.rexp();
            }
        }

	int nsampl;
	void dist_reset()
	{
	    nsampl = 0;
	    timer.reset();
	    for(int i=0; i<NTYPES; i++)
		if( is_particle[i] && species_list[i] != NULL )
		    species_list[i]->dist_reset();
	    field.u_reset();
	}
	void dist_sample()
	{
	    nsampl++;
	    for(int i=0; i<NTYPES; i++)
		if( is_particle[i] && species_list[i] != NULL )
		    species_list[i]->dist_sample();
	    field.u_sample();
	}
	void print_status( ostream & out = cout)
	{
	    double probe_current_total = 0;
	    for(int i=0; i<NTYPES; i++)
		if( is_particle[i] && species_list[i] != NULL )
		    probe_current_total += species_list[i]->probe_current;

	    out << iter << " " << timer.gettime() << " " 
		<< probe_current_total/nsampl * param.probe_length/param.dz << " "
                << field.grid.U_trap << " "
                << field.u[0][(int)(param.z_sampl*4.0/7.5)]
		<< endl;

	    for(int i=0; i<NTYPES; i++)
		if( is_particle[i] && species_list[i] != NULL )
		    species_list[i]->print_status();
	    field.u_print( (param.output_dir + "/potential.dat").c_str() );
	}
	void save()
	{
	    for(int ii=0; ii<NTYPES; ii++)
		if( is_particle[ii] && species_list[ii] != NULL )
		{
		    species_list[ii]->save(param.output_dir + "/particles_" + species_list[ii]->name + ".dat");
		    species_list[ii]->source5_save(param.output_dir + "/particles_source_" + species_list[ii]->name + ".dat");
		}
	};

	gnuplot_ctrl *h1;
	gnuplot_ctrl *h2;
	double *mxw;
	void plot_init()
	{
	    h1 = gnuplot_init();
	    h2 = gnuplot_init();
	    gnuplot_cmd(h1,"set pm3d");
	    gnuplot_cmd(h1,"set view map");
	    //gnuplot_cmd(h1,"set view map");
	    gnuplot_cmd(h1,"unset surface");
	    gnuplot_cmd(h2,"set log y");
	    mxw = new double[species_list[ELECTRON]->energy_dist.N_hist()];
	}
	void plot_destroy()
	{
	    gnuplot_close(h1);
	    gnuplot_close(h2);
	    delete [] mxw;
	}
	void plot()
	{
	    gnuplot_resetplot(h1);
	    gnuplot_resetplot(h2);
	    int spindex;
	    for(spindex=0; spindex<NTYPES; spindex++)
		if(is_particle[spindex] && species_list[spindex] != NULL && species_list[spindex]->particles.size() > 0 )
		    break;
	    for(int ii=0; ii<species_list[spindex]->energy_dist.N_hist(); ii++)
	    {
		double f;
		f = species_list[spindex]->energy_dist.position(ii);
		mxw[ii] = exp(-param.q_e*f/(param.k_B*species_list[spindex]->temperature))*
		    sqrt(param.q_e*f/(M_PI*pow(param.k_B*species_list[spindex]->temperature,3)))*2*
		    param.q_e*species_list[spindex]->energy_dist.norm();
	    }
	    gnuplot_setstyle(h2,"lines");
	    gnuplot_plot_x(h2, mxw, species_list[ELECTRON]->energy_dist.N_hist(),"MAxwell");
	    gnuplot_setstyle(h2,"points");
	    for(int i=0; i<NTYPES; i++)
		if( is_particle[i] && 
			species_list[i] != NULL && 
			species_list[i]->particles.size() > 0)
		    gnuplot_plot_x(h2, &(species_list[i]->energy_dist[0]), species_list[i]->energy_dist.N_hist(),species_list[i]->name.c_str());
	    //gnuplot_plot_x(h2, &(species_list[O_NEG]->source_energy_dist[0]), 
		    //species_list[O_NEG]->source_energy_dist.N_hist(),(species_list[O_NEG]->name + "source").c_str());
		usleep(100);

	    //gnuplot_cmd(h1, "set cbrange [10:-5]");
	    ostringstream cmd_str;
	    cmd_str << "set zrange [" <<param.extern_field*param.r_max*0.5 << ":" 
		<< min(-param.extern_field*param.r_max*0.5, param.u_probe) << "]";
	    string cmd(cmd_str.str());
	    //gnuplot_cmd(h1, "set zrange [10:-5]");
	    
	    //gnuplot_cmd(h1, cmd.c_str());
	    //gnuplot_cmd(h1, "set zrange [-1e-15:5e-13]");
	    //gnuplot_splot_grid(h1, field.u[0], param.r_sampl, param.z_sampl, "potencial");
	    species_list[ELECTRON]->rho.reset();
	    species_list[ELECTRON]->accumulate();
	    gnuplot_splot_grid(h1, species_list[ELECTRON]->rho[0], param.r_sampl, param.z_sampl, "rho");

	    ofstream fw("plot.dat");
	    for(int ii=0; ii<param.z_sampl; ii++)
	    {
		for(int jj=0; jj<param.r_sampl; jj++)
		    fw << jj*param.dr <<" "<< ii*param.dz <<" "<< field.u[ii][jj] <<endl;
		fw << endl;
	    }
	}
};
