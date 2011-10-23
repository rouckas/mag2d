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

    // initialize the PIC model !
    Pic<CYLINDRICAL> pic(param);

    pic.species_list[ELECTRON]->resize(10000);
    pic.species_list[ELECTRON]->add_monoenergetic_particles_on_cylinder_cylindrical(10000, 1.0, 5e-2, 4e-3, 6e-2);

    ofstream fw((param.output_dir+"/out.dat").c_str());
    ofstream fwt1((param.output_dir+"/traj1.dat").c_str());
    ofstream fwt2((param.output_dir+"/traj2.dat").c_str());
    ofstream fwv1((param.output_dir+"/vel1.dat").c_str());
    ofstream fwv2((param.output_dir+"/vel2.dat").c_str());

    t_particle *p_p;
    p_p = &(pic.species_list[ELECTRON]->particles[0]);
    //p_p->vz = -fabs(p_p->vz);
    //p_p->r = 1e-2;
    //p_p->z = 1.1e-2;
    //p_p->vz = 1e4;
    //p_p->vx = -pic.species_list[H_NEG]->veV(0.01);
    //cout <<"veV "<< -pic.species_list[H_NEG]->veV(1) <<endl;
    //p_p->vy = 0;
/*
    p_p[1].z = 1e-2;
    p_p[1].x = 1.1e-2;
    p_p[1].vx = 0;
    p_p[1].vy = 0;
    p_p[1].vz = -pic.species_list[H_NEG]->veV(.001);
    p_p[1].vx = 0.3*p_p[1].vz;
    */
    for(unsigned long int i=1; i<param.niter+1; ++i)
    {

	pic.advance();
        if(i%10==11)
        {
            int p1 = 3;
            int p2 = 5;
            fwt1 << p_p[p1].x <<' '<< setprecision(10) << p_p[p1].z <<endl;
            fwt2 << p_p[p2].x <<' '<< setprecision(10) << p_p[p2].z <<endl;
            fwv1 << p_p[p1].vx <<' '<< p_p[p1].vz <<' '<< p_p[p1].vy <<endl;
            fwv2 << setprecision(10) << p_p[p2].vx <<' '<< p_p[p2].vz <<' '<< p_p[p2].vy <<endl;
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

    }
    //getchar();

    return 0;
}
