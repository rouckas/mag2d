//#include "random.h"
#define GNUPLOT
#include "fields.hpp"
#include "param.hpp"
#ifdef GNUPLOT
#include "gnuplot_i.h"
#endif
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include "random.cpp"
#include "timer.hpp"
#include "particles.hpp"
using namespace std;

#define NPARTICL_SAFE 0.5
#define SRC_FACT 20

struct rusage picusage2,picusage1;

template <int D>
class Speclist
{
    public:
        Speclist(string configfile, Param & param, t_random & rnd, Fields & fields)
        {
            vector<SpeciesParams*> vs;
            vector<InteractionParams*> vi;
            config_parse(configfile, vs, vi);

            /*
             * create Species classes
             */
            for(size_t i=0; i<vs.size(); i++)
                data.push_back(new Species<D>(vs[i], param, rnd, fields));

            // prepare room for lists of interactions by species
            for(size_t i=0; i<vs.size(); i++)
            {
                for(size_t j=0; j<vs.size(); j++)
                    data[i]->speclist.push_back(data[j]);
                data[i]->interactions_by_species.resize(vs.size());
                data[i]->rates_by_species.resize(vs.size());
            }

            /*
             * create Interaction classes
             */
            Interaction * pi;
            for(size_t i=0; i<vi.size(); i++)
            {
                pi = new Interaction(vi[i]);

                // find the primary interacting species
                size_t j=0;
                while(j<data.size() && vi[i]->primary != data[j]->name)
                    j++;
                if(j== data.size())
                    throw runtime_error("Speclist::Speclist: unrecognized primary species \"" + vi[i]->primary + "\" of interaction \"" + vi[i]->name + "\"\n");
                pi->primary = data[j];
                data[j]->interactions.push_back(pi);

                // find the secondary interacting species
                size_t k=0;
                while(k<data.size() && vi[i]->secondary != data[k]->name)
                    k++;
                if(k== data.size())
                    throw runtime_error("Speclist::Speclist: unrecognized primary species \"" + vi[i]->secondary + "\" of interaction \"" + vi[i]->name + "\"\n");
                pi->secondary = data[k];
                data[j]->interactions_by_species[k].push_back(pi);
            }

            for(size_t i=0; i<data.size(); i++)
                data[i]->lifetime_init();

        };
        Species<D> * operator [] (size_t i)
        {
            return data[i];
        };
        Species<D> * operator [] (string name)
        {
            size_t j=0;
            while(j<data.size() && name != data[j]->name)
                j++;
            // A map can be used for this
            return data[j];
        };
        ~Speclist()
        {
            for(size_t i=0; i<data.size(); i++)
                delete data[i];
        }
        size_t size()
        {
            return data.size();
        }


        vector<Species<D>*> data;
};

template <int D>
class Pic
{
    private:
	t_timer timer;
    public:
	Param & param;
        map<string, SpeciesType> string2speciestype;
	Fields field;
	t_random rnd;
        Speclist<D> speclist;
	unsigned long int iter;
	Pic(Param & _param) : param(_param), field(param), speclist(_param.species_conf_file, _param, rnd, field), iter(0)
	{
	    field.reset();

	    for(size_t i=0; i<speclist.size(); i++)
            {
                string2speciestype[speclist[i]->name] = (SpeciesType)i;
            }


	    if(param.particle_reload)
		for(size_t ii=0; ii<speclist.size(); ii++)
		    if( speclist[ii]->particle )
		    {
			speclist[ii]->load(param.particle_reload_dir + "/particles_" + speclist[ii]->name + ".dat");
			speclist[ii]->source5_load(param.particle_reload_dir + "/particles_source_" + speclist[ii]->name + ".dat");
		    }

            if(!param.magnetic_field_const)
                field.load_magnetic_field(param.magnetic_field_file.c_str());

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
	    if(param.do_plot) plot_destroy();
	}
        void run_initscript(string filename)
        {
            ifstream fr(filename.c_str());
            string line;
            while(fr.good())
            {
                getline(fr, line);
                if(line[0] == '#') continue;
                istringstream s_line;
                vector<string> toks(0);
                s_line.str(line);
                string tok;
                while( s_line >> tok )
                    toks.push_back(tok);
                if(toks.size()==0) continue;

                if(toks[0] == "add_particles_bessel")
                {
                    if(toks.size() != 6)
                        throw runtime_error("Pic::run_initscript: wrong number of" +
                               string(" parameters (") + int2string(toks.size()) + ") to " + toks[0] + "\n");
                    if(string2speciestype.find(toks[1])==string2speciestype.end())
                        throw runtime_error("Pic::run_initscript: unrecognized species type \"" + toks[1] + "\"\n");
                    SpeciesType species = string2speciestype[toks[1]];
                    int nparticles = string2<int>(toks[2]);
                    double centerx = string2<double>(toks[3]);
                    double centery = string2<double>(toks[4]);
                    double radius = string2<double>(toks[5]);
                    speclist[species]->add_particles_bessel(nparticles, centerx, centery, radius);
                    cout << speclist[species]->name <<" add bessel " <<centerx<<" "<< centery<<" "<<radius<< "\n";
                }
                if(toks[0] == "add_particles_on_disk")
                {
                    if(toks.size() != 6)
                        throw runtime_error("Pic::run_initscript: wrong number of" +
                               string(" parameters (") + int2string(toks.size()) + ") to " + toks[0] + "\n");
                    if(string2speciestype.find(toks[1])==string2speciestype.end())
                        throw runtime_error("Pic::run_initscript: unrecognized species type \"" + toks[1] + "\"\n");
                    SpeciesType species = string2speciestype[toks[1]];
                    int nparticles = string2<int>(toks[2]);
                    double centerx = string2<double>(toks[3]);
                    double centery = string2<double>(toks[4]);
                    double radius = string2<double>(toks[5]);
                    speclist[species]->add_particles_on_disk(nparticles, centerx, centery, radius);
                    cout << speclist[species]->name <<" add on disk " <<centerx<<" "<< centery<<" "<<radius<< "\n";
                }
                else if(toks[0] == "add_tracked_particle")
                {
                    if(toks.size() != 7)
                        throw runtime_error("Pic::run_initscript: wrong number of" +
                               string(" parameters (") + int2string(toks.size()) + ") to " + toks[0] + "\n");
                    if(string2speciestype.find(toks[1])==string2speciestype.end())
                        throw runtime_error("Pic::run_initscript: unrecognized species type \"" + toks[1] + "\"\n");
                    SpeciesType species = string2speciestype[toks[1]];
                    double x = string2<double>(toks[2]);
                    double y = string2<double>(toks[3]);
                    double vx = string2<double>(toks[4]);
                    double vy = string2<double>(toks[5]);
                    double vz = string2<double>(toks[6]);
                    speclist[species]->add_tracked_particle(x, y, vx, vy, vz);
                    cout << speclist[species]->name <<" add tracked particle " <<x<<" "<< y<<" "<<" "<<vx<<" "<<vy<<" "<<vz<<endl;
                }
            }

        }

	void advance()
	{
	    timer.start();
	    if(param.selfconsistent)
	    {
		field.boundary_solve();
                if(param.u_smooth) field.u_smooth();
		field.reset();
		for(size_t i=0; i<speclist.size(); i++)
                    if(speclist[i]->particle && speclist[i]->n_particles() > 0)
		{
		    speclist[i]->rho.reset();
		}
	    }

	    for(size_t i=0; i<speclist.size(); i++)
	    {
		speclist[i]->advance();
		//if(param.selfconsistent)
		  //  species_list[i]->source();
	    }

	    for(size_t i=0; i<speclist.size(); i++)
                if(speclist[i]->particle && param.selfconsistent)
		{
		    field.rho.add(speclist[i]->rho);
		}
	    timer.stop();

	    iter++;
	}
        void emit()
        {
            throw runtime_error("Pic:emit() not implemented\n");
            /*
            double emission_current = 1e-4;     //[mA]
            double fnemit = 1.0;
            fnemit = param.dt[ELECTRON]/param.q_e * emission_current/param.macroparticle_factor;
            int nemit = (int)fnemit;
            nemit += (rnd.uni()>fnemit-nemit) ? 0 : 1;
            for(int i=0; i<nemit; i++)
            {
                int ind = species_list[ELECTRON]->insert();
                t_particle * p_p = &(species_list[ELECTRON]->particles[ind]);
                p_p->x = 2e-3*sqrt(rnd.uni());
                p_p->z = 0.5e-2;
                p_p->vx = species_list[ELECTRON]->veV(1);
                p_p->vz = species_list[ELECTRON]->veV(1);
                species_list[ELECTRON]->rndv(p_p->vx, p_p->vz, p_p->vy);
                //p_p->vx = 0;
                //p_p->vz = -500000;
                //p_p->vy = 0;
                p_p->time_to_death = species_list[ELECTRON]->lifetime * rnd.rexp();
            }
            */
        }

	int nsampl;
	void dist_reset()
	{
	    nsampl = 0;
	    timer.reset();
	    for(size_t i=0; i<speclist.size(); i++)
		if( speclist[i]->particle )
		    speclist[i]->dist_reset();
	    field.u_reset();
	}
	void dist_sample()
	{
	    nsampl++;
	    for(size_t i=0; i<speclist.size(); i++)
		if( speclist[i]->particle )
		    speclist[i]->dist_sample();
	    field.u_sample();
	}
	void print_status( ostream & out = cout)
	{
	    double probe_current_total = 0;
	    for(size_t i=0; i<speclist.size(); i++)
		if( speclist[i]->particle )
		    probe_current_total += speclist[i]->probe_current;

	    out << iter << " " << timer.gettime() << " " 
		<< probe_current_total/nsampl * param.probe_length/param.dz << " "
                << field.grid.U_trap << " "
                << field.u[0][(int)(param.z_sampl*4.0/7.5)]
		<< endl;

	    for(size_t i=0; i<speclist.size(); i++)
		if( speclist[i]->particle )
		    speclist[i]->print_status();
	}
        void print_trace()
        {
	    for(size_t i=0; i<speclist.size(); i++)
		if( speclist[i]->particle)
		    speclist[i]->print_trace();
        }
        void print_distribution()
        {
	    for(size_t i=0; i<speclist.size(); i++)
		if( speclist[i]->particle)
		    speclist[i]->print_distribution();
	    field.u_print( (param.output_dir + "/potential.dat").c_str() );
        }
	void save()
	{
	    for(size_t ii=0; ii<speclist.size(); ii++)
		if( speclist[ii]->particle )
		{
		    speclist[ii]->save(param.output_dir + "/particles_" + speclist[ii]->name + ".dat");
		    speclist[ii]->source5_save(param.output_dir + "/particles_source_" + speclist[ii]->name + ".dat");
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
	    mxw = new double[speclist[ELECTRON]->energy_dist.N_hist()];
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
	    size_t spindex;
	    for(spindex=0; spindex<speclist.size(); spindex++)
		if(speclist[spindex]->particle && speclist[spindex]->particles.size() > 0 )
		    break;
	    for(int ii=0; ii<speclist[spindex]->energy_dist.N_hist(); ii++)
	    {
		double f;
		f = speclist[spindex]->energy_dist.position(ii);
		mxw[ii] = exp(-param.q_e*f/(param.k_B*speclist[spindex]->temperature))*
		    sqrt(param.q_e*f/(M_PI*pow(param.k_B*speclist[spindex]->temperature,3)))*2*
		    param.q_e*speclist[spindex]->energy_dist.norm();
	    }
	    gnuplot_setstyle(h2,"lines");
	    gnuplot_plot_x(h2, mxw, speclist[ELECTRON]->energy_dist.N_hist(),"MAxwell");
	    gnuplot_setstyle(h2,"points");
	    for(size_t i=0; i<speclist.size(); i++)
		if( speclist[i]->particle &&
			speclist[i]->particles.size() > 0)
		    gnuplot_plot_x(h2, &(speclist[i]->energy_dist[0]), speclist[i]->energy_dist.N_hist(),speclist[i]->name.c_str());
	    //gnuplot_plot_x(h2, &(species_list[O_NEG]->source_energy_dist[0]), 
		    //species_list[O_NEG]->source_energy_dist.N_hist(),(species_list[O_NEG]->name + "source").c_str());
		usleep(100);

	    //gnuplot_cmd(h1, "set cbrange [10:-5]");
	    ostringstream cmd_str;
	    cmd_str << "set zrange [" <<param.extern_field*param.x_max*0.5 << ":" 
		<< min(-param.extern_field*param.x_max*0.5, param.u_probe) << "]";
	    string cmd(cmd_str.str());
	    //gnuplot_cmd(h1, "set zrange [10:-5]");
	    
	    //gnuplot_cmd(h1, cmd.c_str());
	    //gnuplot_cmd(h1, "set zrange [-1e-15:5e-13]");
	    //gnuplot_splot_grid(h1, field.u[0], param.x_sampl, param.z_sampl, "potencial");
	    speclist["ELECTRON"]->rho.reset();
	    speclist["ELECTRON"]->accumulate();
	    gnuplot_splot_grid(h1, speclist["ELECTRON"]->rho[0], param.x_sampl, param.z_sampl, "rho");
	}
};
