#define GNUPLOT
#include "param.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "pic.cpp"

#include <GetPot>
#include <string>
#include "output.cpp"
using namespace std;


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
    pic.rnd.initialize_seed(1234);


    pic.run_initscript(initscript);

    ofstream fw((param.output_dir+"/out.dat").c_str());
    t_timer timer;
    timer.reset();
    for(unsigned long int i=1; i<param.niter+1; ++i)
    {

        for(size_t i=0; i<pic.speclist.size(); i++)
	    {
                for(size_t j=0; j<pic.speclist[i]->particles.size(); j++)
                {
                    if(!pic.speclist[i]->particles[j].empty)
                        pic.speclist[i]->remove(j);
                }
	    }
        pic.run_initscript(initscript);

        timer.start();
	pic.advance();
        timer.stop();
	if(
                (param.t_print_dist != 0 && i%param.t_print_dist == 0) || i == param.niter
          )
	{
	    pic.print_distribution();
	    pic.print_field();
	    if(param.do_plot) pic.plot();
	}
	if(param.t_print != 0 && i%param.t_print == 0)
	{
	    pic.print_status(fw);
	    cout << "plot " <<i<<" "<< timer.get_cpu_time()/i*1000 << " ms / iteration" << endl;
	    if(i<param.t_equilib) pic.dist_reset();
	}

    }

    return 0;
}
