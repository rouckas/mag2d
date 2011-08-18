#define GNUPLOT
#include "param.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "pic.cpp"

#include <GetPot>
#include <boost/filesystem.hpp>
#include <string>
#include "output.cpp"
using namespace std;
namespace fs = boost::filesystem;


int t_plot = 10;
int t_print = 20;
int t_source_refresh = 0;
int t_save = 0;

int main(int argc, char *argv[])
{

    GetPot cl(argc,argv);

    // load the configuration from file
    string config_file = cl("config","config.txt");
    string species_conf_file = cl("species_conf","species_conf.txt");
    string initscript = cl("initscript", "initscript.txt");
    GetPot config(config_file.c_str());

    // initialize model parameters
    Param param(config);
    param.species_conf_file = species_conf_file;

    // prepare output directory
    param.output_dir = cl("output_dir","output");
    t_output output(param.output_dir);

    // make backup of configuration
    output.backup(config_file, "config.txt");
    output.backup(species_conf_file, "species_conf.txt");
    output.backup(initscript, "initscript.txt");

    // initialize Poisson solver
    //Fields fields(param);
    //fields.accumulate(1.0, 0.0, param.z_max/2.0);
    //fields.boundary_solve();
    //ofstream fw1((param.output_dir+"/u.dat").c_str());
    //fields.u.print(fw1);

    // initialize the PIC model !
    Pic<CARTESIAN> pic(param);

    pic.run_initscript(initscript);


    ofstream fw((param.output_dir+"/out.dat").c_str());
    for(unsigned long int i=1; i<param.niter+1; ++i)
    {

	pic.advance();
        if(i%1000==0)
        {
            //adaptive lifetime and timestep calculation
            double Emax_factor = 3.0;
            double dt_lifetime = 3.0;
            double Emax = pic.speclist["ELECTRON"]->energy_dist.mean()*Emax_factor;
            pic.speclist["ELECTRON"]->lifetime_init(Emax);
            double lifetime = pic.speclist["ELECTRON"]->lifetime;
            double dt = lifetime/dt_lifetime;
            cout << "ELECTRON lifetime "<< lifetime <<" s;    Emax/5 = "<< Emax/5 << " eV;     dt = "<< dt << endl;
            if(dt > pic.speclist["ELECTRON"]->dt)
            {
                pic.speclist["ELECTRON"]->dt = dt;
                pic.speclist["H3+"]->dt = dt;
            }

        }

        if(i%10==0)
        {
            pic.print_trace();
        }
	if(i%param.t_dist_sample==0) pic.dist_sample();

	if(param.t_print_dist != 0 && i%param.t_print_dist == 0)
	{
	    pic.print_distribution();
	    if(param.do_plot) pic.plot();
	}
	if(param.t_print != 0 && i%param.t_print == 0)
	{
	    pic.print_status(fw);
	    cout << "plot " <<i<<endl;
	    if(i<param.t_equilib) pic.dist_reset();
	}

	if(t_save != 0 && i%t_save==0 )
	{
	    pic.save();
	}

	//XXX deprecated
	if(t_source_refresh != 0 && i%t_source_refresh==0)
	{
	    //pic.species_list[ARGON_POS]->source5_refresh(SRC_FACT);
	    //pic.species_list[ELECTRON]->source5_refresh(SRC_FACT);
	}

    }
    //getchar();

    return 0;
}
