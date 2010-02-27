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
    GetPot config(config_file.c_str());

    // initialize model parameters
    Param param(config);

    // prepare output directory
    param.output_dir = cl("output_dir","output");
    t_output output(param.output_dir);

    // make backup of configuration
    string config_backup = param.output_dir + "/config.txt";
    if(fs::exists(config_backup) && !fs::is_directory(config_backup))
	fs::remove(config_backup);
    fs::copy_file(config_file,config_backup);

    // initialize Poisson solver
    //Fields fields(param);
    //fields.accumulate(1.0, 0.0, param.z_max/2.0);
    //fields.boundary_solve();
    //ofstream fw1((param.output_dir+"/u.dat").c_str());
    //fields.u.print(fw1);

    // initialize the PIC model !
    Pic<CARTESIAN> pic(param);

    string initscript = cl("initscript", "initscript.txt");
    pic.run_initscript(initscript);


    ofstream fw((param.output_dir+"/out.dat").c_str());
    for(unsigned long int i=1; i<param.niter+1; ++i)
    {

	pic.advance();
        if(i%10==0)
        {
            pic.print_trace();
        }
	if(i%param.t_dist_sample==0) pic.dist_sample();

	if(param.t_print != 0 && i%param.t_print == 0)
	{
	    pic.print_status(fw);
	    cout << "plot " <<i<<endl;
	    if(i<param.t_equilib) pic.dist_reset();
	}
	if(param.t_print_dist != 0 && i%param.t_print_dist == 0)
	{
	    pic.print_distribution();
	    if(param.do_plot) pic.plot();
	}
	continue;

	if(t_save != 0 && i%t_save==0 )
	{
	    pic.save();
	    /*
	    for(int ii=0; ii<NTYPES; ii++)
		if( is_particle[ii] && pic.species_list[ii] != NULL )
		{
		    pic.species_list[ii]->save(param.output_dir + "/particles_" + pic.species_list[ii]->name + ".dat");
		    pic.species_list[ii]->source5_save(param.output_dir + "/particles_source_" + pic.species_list[ii]->name + ".dat");
		}
		*/
	}

	//XXX deprecated
	if(t_source_refresh != 0 && i%t_source_refresh==0)
	{
	    pic.species_list[ARGON_POS]->source5_refresh(SRC_FACT);
	    pic.species_list[ELECTRON]->source5_refresh(SRC_FACT);
	}

    }
    //getchar();

    return 0;
}
