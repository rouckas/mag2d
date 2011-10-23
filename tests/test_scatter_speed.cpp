//#include "random.h"
#define GNUPLOT
#include "fields.cpp"
#include "param.cpp"
#ifdef GNUPLOT
#include "gnuplot_i.h"
#endif
#include "argon.cpp"
#include "elonO2.cpp"
//#include "elon.cpp"
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/resource.h>
#include "random.cpp"
#include "input.cpp"
using namespace std;

#define ROWS 100
#define COLS 100
//#define N_PARTICLES 200000
#define N_PARTICLES 1500000
#define N_SOURCE 1000
#define SRC_FACT 20
//#define ADVANCE_OLD

struct rusage usage2,usage1;
#define STOPUJ(cmd,out) getrusage(RUSAGE_SELF,&usage1);\
        cmd\
    getrusage(RUSAGE_SELF,&usage2);\
    out += (double)(usage2.ru_utime.tv_sec-usage1.ru_utime.tv_sec);\
    out += (double)(usage2.ru_utime.tv_usec-usage1.ru_utime.tv_usec)*1e-6;
#define START getrusage(RUSAGE_SELF,&usage1);
#define STOP(out)  getrusage(RUSAGE_SELF,&usage2);\
    out += (double)(usage2.ru_utime.tv_sec-usage1.ru_utime.tv_sec);\
    out += (double)(usage2.ru_utime.tv_usec-usage1.ru_utime.tv_usec)*1e-6;
int t_plot = 10;
int t_print = 10;
int t_source_refresh = 50;
int t_save = 10000;
double dt_elon = 1e-11;
int niter = 620000*1e-11/dt_elon;

int main(int argc, char *argv[])
{

    if(argc >= 2)
    {
	istringstream stream(argv[1]);
	stream >> dt_elon;
	niter = 620000*1e-11/dt_elon;
    }

    t_param param(1e-2,1e-2, ROWS, COLS, N_PARTICLES, 1e15,133*1e0);
    param.probe_radius = 5e-4;
    param.extern_field = 200.0;
    param.u_probe = -5.0;

    t_field field(param);

    //t_random rnd(123);
    t_random rnd(1234);

    const int nspecies = 6;
    t_species * species_list[nspecies];
    for(int i=0; i<nspecies; i++)
	species_list[i] = NULL;

    t_argon argon(N_PARTICLES+50000,N_PARTICLES,300.0,1.22e-15,param, rnd, &field, species_list, dt_elon*1e3);
    species_list[ARGON_POS] = &argon;
    argon.source5_refresh(SRC_FACT);

    t_elon elon(N_PARTICLES+50000,N_PARTICLES, 23209.8,3.214e-13,param, rnd, &field, species_list, dt_elon);
    species_list[ELECTRON] = &elon;
    elon.source5_refresh(SRC_FACT);

    t_argon_neutral argon_neutral(0,0,300.0,0,param,rnd,&field,species_list,0);
    species_list[ARGON] = &argon_neutral;
    // wrong for multicomponent plasma, TODO
    argon_neutral.density = param.neutral_density;

    elon.lifetime_init();

    int tmp = time(NULL);
    /*
       stringstream elon_file;
       elon_file << "elon_" << tmp << ".dat";
       stringstream argon_file;
       argon_file << "argon_" << tmp << ".dat";
       */
    stringstream elon_file;
    elon_file << "elon.dat";
    stringstream argon_file;
    argon_file << "argon.dat";

    argon.save("argon.dat");
    elon.save("elon.dat");
/*
    elon.energy_dist_compute();
    cout << "E0 " << elon.energy_dist.mean_tot() << endl;
    t_elon elon2(N_PARTICLES+50000,N_PARTICLES, 23209.8,3.214e-13,param, rnd, &field, dt_elon);
    elon2.load("elon.dat");
    elon2.energy_dist_compute();
    cout << "E1 " << elon2.energy_dist.mean_tot() << endl;
	argon.load("argon_1196894732.dat");
	elon.load("elon_1196894732.dat");
    elon.energy_dist_compute();
    cout << "E3 " << elon.energy_dist.mean_tot() << endl;
    field.reset();
    field.boundary_solve();
    elon.advance2();
    cout << "E3 " << elon.energy_dist.mean_tot() << endl;
  */  //argon.save("argon.dat");
    //argon.load("argon.dat");
    //elon.load("elon.dat");


    /*
    {
	vector<t_particle> particles(10);
	std::ofstream ofs("filename");
	boost::archive::text_oarchive oa(ofs);
	oa << particles;
    }
    */

#ifdef GNUPLOT
    gnuplot_ctrl *h1;
    gnuplot_ctrl *h2;
    h1 = gnuplot_init();
    h2 = gnuplot_init();
    gnuplot_cmd(h1,"set pm3d");
    //gnuplot_cmd(h1,"set view map");
    cout << "unset surface\n";
    gnuplot_cmd(h1,"unset surface");
    gnuplot_cmd(h2,"set log y");
#endif
    double *mxw = new double[argon.energy_dist.N_hist()];

    cout << "# dt_el=" << elon.dt << endl;
    cout << "# dt_ar=" << argon.dt << endl;
    cout << "# tau_ar=" << argon.lifetime << endl;
    cout << "# tau_el=" << elon.lifetime << endl;
    cout << "# u_probe=" << param.u_probe << endl;
    cout << "# pressure=" << param.pressure << endl;
#ifdef ADVANCE_OLD
    cout << "# advance_old()" << endl;
#else
    cout << "# advance2()" << endl;
#endif
    cout << endl;


    field.reset();
    double probe_current=0;
    double elon_probe_current=0;
    double argon_n_particles = 0,
	   elon_n_particles = 0,
	   argon_energy_dist_mean = 0,
	   elon_energy_dist_mean = 0;
    for(int i=0; i<niter; ++i)
    {
	field.boundary_solve();
	field.reset();
	double time1=0;
	STOPUJ(
		for(unsigned int j=0; j<elon.particles.size(); j++)
		{
		    elon.scatter(elon.particles[j]);
		}
	    ,time1
	    );

	if(i%t_save==0 && i!=0)
	{
///	    argon.load(argon_file.str());
//	    elon.load(elon_file.str());
	    argon.save(argon_file.str());
	    elon.save(elon_file.str());
	}

	if(i%t_source_refresh==0)
	{
	    argon.source5_refresh(SRC_FACT);
	    elon.source5_refresh(SRC_FACT);
	}

	//accumulate values for averaging
	argon_n_particles += argon.n_particles();
	elon_n_particles += elon.n_particles();
	argon_energy_dist_mean += argon.energy_dist.mean_tot();
	elon_energy_dist_mean += elon.energy_dist.mean_tot();
	probe_current += elon.probe_current;
	probe_current += argon.probe_current;
	elon_probe_current += elon.probe_current;

	if(i%t_print == 0 && t_print != 0)
	{
	    argon.energy_dist.reset();
	    elon.energy_dist.reset();
	    argon.energy_dist_compute();
	    elon.energy_dist_compute();
	    cout << i << " " << time1 << " " 
		<< double(argon_n_particles)/t_print << " " 
		<< double(elon_n_particles)/t_print << " "
	       	<< argon_energy_dist_mean/t_print << " " 
		<< elon_energy_dist_mean/t_print << " " 
		<< probe_current/t_print << " "
	       	<< elon.probe_dist.mean_tot() << " " 
		<< argon.probe_dist.mean_tot() << " " 
		<< elon_probe_current/elon.charge << endl;
	    //cout << i << " " << time1 << " " << argon.n_particles()<< " " << elon.n_particles() << " " <<
	//	argon.energy_dist.mean_tot() << " " << elon.energy_dist.mean_tot() << " " << probe_current/t_print << " "
	  //     	<< elon.probe_dist.mean_tot() << " " << argon.probe_dist.mean_tot() << " " << elon_probe_current/elon.charge << endl;
	    cout << elon.hit1 << ' ' << elon.hit2 <<endl;
	    elon.probe_dist.reset();
	    argon.probe_dist.reset();
	    argon_n_particles = 0;
	    elon_n_particles = 0;
	    argon_energy_dist_mean = 0;
	    elon_energy_dist_mean = 0;
	    probe_current=0;
	    elon_probe_current=0;
	}
	
#ifdef GNUPLOT
	if(t_plot != 0 && i%t_plot == 0)
	{
	    gnuplot_resetplot(h1);
	    gnuplot_resetplot(h2);
	    for(int ii=0; ii<argon.energy_dist.N_hist(); ii++)
	    {
		double f;
		f = argon.energy_dist.position(ii);
		mxw[ii] = exp(-param.q_e*f/(param.k_B*argon.temperature))*
		    sqrt(param.q_e*f/(M_PI*pow(param.k_B*argon.temperature,3)))*2*
		    (param.q_e)*(argon.energy_dist.Max()-argon.energy_dist.Min())/
		    argon.energy_dist.N_hist()*argon.energy_dist.N_val();
	    }
	    gnuplot_setstyle(h2,"lines");
	    gnuplot_plot_x(h2, mxw, argon.energy_dist.N_hist(),"MAxwell");
	    gnuplot_setstyle(h2,"points");
	    gnuplot_plot_x(h2, &(elon.energy_dist[0]), elon.energy_dist.N_hist(),"e-");
	    gnuplot_plot_x(h2, &(argon.energy_dist[0]), argon.energy_dist.N_hist(),"Ar+");
	    usleep(100);

	    gnuplot_cmd(h1, "set cbrange [2:-5]");
	    gnuplot_cmd(h1, "set zrange [2:-5]");
	    //gnuplot_cmd(h1, "set zrange [-1e-15:5e-13]");
	    gnuplot_splot_grid(h1, field.u[0], ROWS, COLS, "potencial");

	    ofstream fw("plot.dat");
	    for(int ii=0; ii<ROWS; ii++)
	    {
		for(int jj=0; jj<COLS; jj++)
		    fw << jj*param.dx <<" "<< ii*param.dy <<" "<< field.u[ii][jj] <<endl;
		fw << endl;
	    }

	}
#endif
	//usleep(100);
    }

    //getchar();

#ifdef GNUPLOT
    gnuplot_close(h1);
    gnuplot_close(h2);
#endif

    return 0;
}
