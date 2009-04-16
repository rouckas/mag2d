#define GNUPLOT
#include "param.cpp"
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
    ofstream fw1((param.output_dir+"/u.dat").c_str());
    //fields.u.print(fw1);

    // initialize the PIC model !
    Pic pic(param);

    pic.species_list[H_NEG]->add_particles_on_disk(1000, 1e-2, 1e-2, 3e-3);

    ofstream fw((param.output_dir+"/out.dat").c_str());
    ofstream fwt1((param.output_dir+"/traj1.dat").c_str());
    ofstream fwt2((param.output_dir+"/traj2.dat").c_str());
    ofstream fwv1((param.output_dir+"/vel1.dat").c_str());
    ofstream fwv2((param.output_dir+"/vel2.dat").c_str());

    t_particle *p_p;
    p_p = &(pic.species_list[H_NEG]->particles[0]);
    p_p->vz = -fabs(p_p->vz);
    p_p->r = 1e-2;
    p_p->z = 1.1e-2;
    p_p->vz = 1e4;
    p_p->vr = -pic.species_list[H_NEG]->veV(0.01);
    cout <<"veV "<< -pic.species_list[H_NEG]->veV(1) <<endl;
    p_p->vt = 0;

    p_p[1].z = 1e-2;
    p_p[1].r = 1.1e-2;
    p_p[1].vr = 0;
    p_p[1].vt = 0;
    p_p[1].vz = -pic.species_list[H_NEG]->veV(.001);
    p_p[1].vr = 0.3*p_p[1].vz;
    for(int i=1; i<param.niter+1; ++i)
    {

	pic.advance();
        if(i%10==0)
        {
            fwt1 << p_p->r <<' '<< setprecision(10) << p_p->z <<endl;
            fwt2 << p_p[1].r <<' '<< setprecision(10) << p_p[1].z <<endl;
            fwv1 << p_p->vr <<' '<< p_p->vz <<' '<< p_p->vt <<endl;
            fwv2 << setprecision(10) << p_p[1].vr <<' '<< p_p[1].vz <<' '<< p_p[1].vt <<endl;
        }
	if(i%20==0) pic.dist_sample();

	if(param.t_print != 0 && i%param.t_print == 0)
	{
	    pic.print_status(fw);
	    cout << "plot " <<i<<endl;
	    if(param.do_plot) pic.plot();
	    if(i<param.t_equilib) pic.dist_reset();
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
