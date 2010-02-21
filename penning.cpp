#define GNUPLOT
#include "param.hpp"
#include <iostream>
#include <iomanip>
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
    Pic<CYLINDRICAL> pic(param);

    pic.species_list[ELECTRON]->resize(2);

    ofstream fw((param.output_dir+"/out.dat").c_str());
    ofstream fwt1((param.output_dir+"/traj1.dat").c_str());
    ofstream fwt2((param.output_dir+"/traj2.dat").c_str());
    ofstream fwv1((param.output_dir+"/vel1.dat").c_str());
    ofstream fwv2((param.output_dir+"/vel2.dat").c_str());

    t_particle *p_p = &(pic.species_list[ELECTRON]->particles[0]);
    p_p = &(pic.species_list[ELECTRON]->particles[0]);
    p_p->vz = -fabs(p_p->vz);
    p_p->x = 1e-4;
    p_p->z = 4.1e-2;
    p_p->vx = 1e5;
    p_p->vz = 1e4;
    p_p->vx = -pic.species_list[ELECTRON]->veV(1);
    p_p->vz = pic.species_list[ELECTRON]->veV(1);
    cout <<"veV "<< -pic.species_list[ELECTRON]->veV(1) <<endl;
    p_p->vy = 0;
    p_p->empty  = false;

    p_p[1].z = 1e-2;
    p_p[1].x = 1.1e-2;
    p_p[1].vx = 0;
    p_p[1].vy = 0;
    p_p[1].vz = -pic.species_list[ELECTRON]->veV(.001);
    p_p[1].vx = 0.3*p_p[1].vz;
    p_p[1].empty  = false;

    for(unsigned long int i=1; i<param.niter+1; ++i)
    {

	pic.advance();

        double t0 = 30000;
        double t1 = 430000;
        t1 = 0;
        double Ut = 5.0;
        //double dU = Ut/(t1-t0);
        double dU = 3e-5;
        // check the voltage in trap center
        double U_center = pic.field.u[0][(int)(param.z_sampl*4.0/7.5)];

        if(pic.iter<t1) pic.emit();
        //decrease trap potential
        if(pic.iter>t0 && pic.iter<t1 && U_center<0.1) pic.field.grid.penning_trap_simple(pic.field.grid.U_trap+dU);


	fwt1 << p_p->x <<' '<< setprecision(10) << p_p->z <<endl;
	fwt2 << p_p[1].x <<' '<< setprecision(10) << p_p[1].z <<endl;
	fwv1 << p_p->vx <<' '<< p_p->vz <<' '<< p_p->vy <<endl;
	fwv2 << setprecision(10) << p_p[1].vx <<' '<< p_p[1].vz <<' '<< p_p[1].vy <<endl;
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
    }
    //getchar();

    return 0;
}
