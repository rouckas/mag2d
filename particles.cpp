#ifndef PARTICLES_H
#define PARTICLES_H

#include <list>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "fields.cpp"
#include "param.cpp"
#include "histogram.hpp"
#include "random.cpp"
#include "mymath.cpp"
#include "input.cpp"
#include "tabulate.cpp"
#include "speclist.hpp"

using namespace std;
//#define NTYPES 8
//typedef enum species_type { NONE, ELECTRON, ARGON, ARGON_POS, ARGON_META, O2, O2_POS, O2_NEG };
//const char * species_names[NTYPES] = {"NONE", "ELECTRON", "ARGON", "ARGON_POS", "ARGON_META", "O2", "O2_POS", "O2_NEG" };

class t_particle
{
    public:
    double r,z;
    double vr,vz,vt;
    double time_to_death;
    bool empty;
};



class Species
//abstract class representing one species in the model
{
    private:
	vector<int> empty;
	Fields *field;
	double source_len;
	ofstream output;
	int niter;
    protected:
	Species ** species_list;
	Param *p_param;
	t_random * rnd;
    public:
	vector<t_particle> particles;
	species_type type;
	string name;
	double mass;
	double charge;
	double lifetime;
	double temperature;
	double density;
	double v_max;
	double dt;
	double EeV(double v) { return 0.5*mass*v*v/(p_param->q_e); }
	double veV(double EeV) { return sqrt(EeV*(p_param->q_e)/mass*2.0); }
	void set_pressure(double pa){ density = pa/(p_param->k_B*temperature);}
	Field rho;
	void accumulate();


	
	// diagnostics
	Histogram energy_dist;
	void energy_dist_compute();
	Histogram source_energy_dist;
	void source_energy_dist_compute();
	Histogram probe_energy_dist;
	Histogram probe_angular_dist;
	Histogram probe_angular_normalized_dist;
	double probe_current;
        double probe_charge;
	double probe_current_avg;
	Field rhoAverage;
	int nsampl;
	double probe_current_sum;


	Species(int n1, int n2, double temperature, Param &param, t_random &rnd, 
		Fields *_field, double vsigma_max, double mass, double _dt, Species * _species_list[], species_type _type = NONE );
	Species(int n1, int n2, Param &param, t_random &_rnd, Fields *_field,
	       double mass, Species * _species_list[], species_type _type = NONE );

        void add_particles_on_disk(int nparticles, double centerx, double centery, double radius);
        void add_particle_beam_on_disk(int nparticles, double centerx, double centery, double radius);
        void add_particle_beam_on_disk_cylindrical(int nparticles, double centerz, double radius);
        void add_monoenergetic_particles_on_cylinder_cylindrical(int nparticles, double energy, double centerz, double radius, double height = 0.0);

	// generate randomly new particle velocity according to species' velocity
	// distribution
	void rndv(double & vr, double & vz, double & vt)
	{
	    vr = rnd->rnor()*v_max*M_SQRT1_2;
	    vz = rnd->rnor()*v_max*M_SQRT1_2;
	    vt = rnd->rnor()*v_max*M_SQRT1_2;
	    //alternatively we can take random sample from source particles
	    /*
	    int i = rnd->inuni() % source2_particles.size();
	    vr = source2_particles[i].vr;
	    vz = source2_particles[i].vz;
	    vz = source2_particles[i].vz;
	    */
	}

	void rndv(double & vr)
	{
	    vr = rnd->maxwell1(v_max);
	}


	//void save(const char *filename)
	void save(const string filename)
	{
	    ofstream fw(filename.c_str(),ios::out | ios::binary);
	    if(fw==false)
		throw runtime_error("Species::save(): failed opening file");

	    int tmp;
	    tmp = particles.size();
	    fw.write((char*)&tmp,sizeof(int));

	    int nparticles = particles.size() - empty.size();
	    fw.write((char*)&nparticles,sizeof(int));
	    tmp=0;
	    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); I++)
	    {
		if(I->empty==false)
		{
		    if(++tmp > nparticles)
		    {
			cerr << "t_particle::save(): data incosistent" << endl;
			break;
		    }
		    fw.write((char*)&(*I),sizeof(t_particle));
		}
	    }
	    //cerr << tmp << " " << nparticles << endl;

	}
	void load(const string filename)
	{
	    ifstream fr(filename.c_str(), ios::in | ios::binary);
	    if(fr==false)
		throw runtime_error("Species::load(): failed opening file");

	    int tmp;
	    fr.read((char*)&tmp,sizeof(int));
	    particles.resize(tmp);

	    int nparticles;
	    fr.read((char*)&nparticles,sizeof(int));
	    for(int i=0; i<nparticles; i++)
	    {
		fr.read((char*)&particles[i], sizeof(t_particle));
		if(!fr.good()) cerr << "Species::load(): read error\n";
		// map particles back to working area if they left it during previous
		// nonselfconsistent simulation
		if(particles[i].r < 0 || particles[i].r > p_param->r_max)
		    particles[i].r = p_param->r_max * rnd->uni();
		if(particles[i].z < 0 || particles[i].z > p_param->z_max)
		    particles[i].z = p_param->z_max * rnd->uni();
	    }

	    empty.resize(0);
	    for(unsigned int i=nparticles; i<particles.size(); i++)
	    {
		empty.push_back(i);
		particles[i].empty=true;
	    }

	}

	void remove(int n)
	{
	    if(particles[n].empty == true)
		throw runtime_error("removing empty particle");
	    else
	    {
		particles[n].empty = true;
		empty.push_back(n);
	    }
	}
	void remove(vector<t_particle>::iterator I)
	{
	    remove(I-particles.begin());
	}
	int insert()
	{
	    if(empty.size()==0)
		throw runtime_error("particle array full");
	    int tmp = empty.back();
	    if(particles[tmp].empty == false)
		throw runtime_error("error inserting particle");
	    empty.pop_back();
	    particles[tmp].empty = false;
	    return tmp;
	}
	int n_particles(){ return particles.size() - empty.size();}
        void resize(unsigned int newsize)
        {
            unsigned int oldsize = particles.size();
            if(newsize > oldsize)
            {
                particles.resize(newsize);
                for(unsigned int i = oldsize; i < newsize; i++)
                {
                    particles[i].empty = true;
                    empty.push_back(i);
                }
            }
        }

	virtual void scatter(t_particle &particle) = 0;
	void advance();
	void source();

	vector<t_particle> source2_particles;
	void source_old();
	void source5_refresh(unsigned int factor);
	void source5_save(const string filename)
	{
	    ofstream fw(filename.c_str(),ios::out | ios::binary);
	    if(fw==false)
		throw runtime_error("Species::source_save(): failed opening file");

	    int nparticles = source2_particles.size() ;
	    fw.write((char*)&nparticles,sizeof(int));

	    for(vector<t_particle>::iterator I = source2_particles.begin(); I != source2_particles.end(); I++)
		fw.write((char*)&(*I),sizeof(t_particle));
	}
	void source5_load(const string filename)
	{
	    ifstream fr(filename.c_str(), ios::in | ios::binary);
	    if(fr==false)
	    {
		cerr << " source5_load(): Warning: cannot open file " << filename <<endl;
		return;
	    }

	    int nparticles;
	    fr.read((char*)&nparticles,sizeof(int));
	    source2_particles.resize(nparticles);

	    for(int i=0; i<nparticles; i++)
	    {
		fr.read((char*)&source2_particles[i], sizeof(t_particle));
		if(!fr.good()) cerr << "Species::load(): read error\n";
	    }
	}
   // private:
	int source5_factor;
	virtual void lifetime_init()
	{
	    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); I++) 
		if( I->empty == false)
		    I->time_to_death = rnd->rexp()*lifetime;

	    for(vector<t_particle>::iterator I = source2_particles.begin(); I != source2_particles.end(); I++) 
		if( I->empty == false)
		    I->time_to_death = rnd->rexp()*lifetime;
	}

    protected:
	void load_CS(const string & fname, vector<vec_interpolate*> & CS, vector<double> & Loss,
	       	const vector<string> & CSnames, const string & tag="", int ncols=3)
	{
	    vector<double> Edata,CSdata;
	    double Eloss;
	    vec_interpolate *p_CS;

	    for(vector<string>::const_iterator I=CSnames.begin(); I != CSnames.end(); I++)
	    {
		::load(fname,Edata,CSdata,*I,&Eloss,tag,ncols);
		p_CS = new vec_interpolate(Edata,CSdata);

		CS.push_back(p_CS);
		Loss.push_back(Eloss);
	    }

	}
	double svmax_find(const vector<vec_interpolate*> & CS, double _vmax, int samples=1000)
	    // this routine assumes input CS in units 1e-16 cm^2, possible source of errors...
	{
	    double dv = _vmax/100.0;
	    double svmax=0.0;
	    //cout << "svmax_find: vmax = "<< _vmax<<endl;
	    for(double v=0; v<_vmax; v+=dv)
	    {
		double sv=0;
		for(unsigned int i=0; i<CS.size(); i++)
		    sv += (*CS[i])(EeV(v)) * v;
		if(isnan(sv))continue;
		if(sv>svmax) svmax=sv;
		//cout <<name <<' '<< v <<' '<< (*CS[0])(EeV(v)) + (*CS[1])(EeV(v)) <<' '<< sv << endl;
	    }
	    svmax *= 1e-20;
	    return svmax;
	}
	void probe_collect(t_particle *I);
    public:
	void print_status( ostream & out = cout)
	{
	    output << niter << " " 
		<< n_particles() << " "
		<< energy_dist.mean_tot() << " "
		//<< probe_current_sum/nsampl*p_param->probe_length/p_param->dz << " "
		<< probe_charge << " "
		<< probe_energy_dist.mean_tot() << " " << probe_current/(nsampl*dt)<< endl;
	    //energy_dist.reset();
	    //energy_dist_compute();
	    energy_dist.print( (p_param->output_dir + "/" + name + "_energy_dist.dat").c_str() );
	    probe_energy_dist.print( (p_param->output_dir + "/" + name + "_probe_energy_dist.dat").c_str() );
	    rhoAverage.print( (p_param->output_dir + "/" + name + "_rho.dat").c_str() , 1.0/nsampl);
	    probe_angular_dist.print( (p_param->output_dir + "/" + name + "_probe_angular_dist.dat").c_str() );
	    probe_angular_normalized_dist.print( (p_param->output_dir + "/" + name + "_probe_angular_normalized_dist.dat").c_str() );
	    //exit(0);
	    // print energy distribution to file DONE
	    // print probe energy distribution to file DONE
	    // TODO print probe angular distribution to file
	    // print spatial density distribution to file DONE
	}
	void dist_sample()
	{
	    energy_dist_compute();
	    source_energy_dist_compute();
	    probe_current_sum += probe_current;
	    rhoAverage.add(rho);
	    nsampl++;
	}
	void dist_reset()
	{
	    energy_dist.reset();
	    source_energy_dist.reset();
	    probe_energy_dist.reset();
	    probe_angular_dist.reset();
	    probe_angular_normalized_dist.reset();
	    rhoAverage.reset();
	    nsampl = 0;
	    probe_current_sum = 0;
            probe_current = 0;
	}
 

};

Species::Species(int _n, int n2, Param &param, t_random &_rnd, Fields *_field,
	double _mass, Species * _species_list[], species_type _type)
    : empty(0),
    field(_field), niter(0), p_param(&param), rnd(&_rnd), 
    particles(_n), type(_type), name(species_names[_type]), mass(_mass),
    temperature(param.temperature[_type]), density(param.density[_type]), dt(param.dt[_type]), 
    rho(param),
    energy_dist(100,0.0,temperature*param.k_B/param.q_e*10.0),
    source_energy_dist(100,0.0,temperature*param.k_B/param.q_e*10.0),
    probe_energy_dist(100,0.0,temperature*param.k_B/param.q_e*20.0),
    probe_angular_dist(30,0.0,M_PI*0.5),
    probe_angular_normalized_dist(30,0.0,M_PI*0.5),
    probe_current(0),
    probe_charge(0),
    rhoAverage(param),
    nsampl(0)
{
    species_list = _species_list;

    lifetime = log(0);
    v_max = sqrt(2.0*Param::k_B*temperature/mass);

    empty.reserve(_n);
    for(int i=n2; i<_n; i++)
    {
	empty.push_back(i);
	particles[i].empty = true;
    }

    for(int i=0; i<n2; i++)
    {
	particles[i].empty = false;
	particles[i].r = param.r_max*(rnd->uni());
	particles[i].z = param.z_max*(rnd->uni());
	//XXX BAD in cylindrical coords (maybe not if vt = d(theta)/dt*r):
	particles[i].vr = rnd->rnor()*v_max/(M_SQRT2);
	particles[i].vz = rnd->rnor()*v_max/(M_SQRT2);
	particles[i].vt = rnd->rnor()*v_max/(M_SQRT2);
	particles[i].time_to_death = rnd->rexp()*lifetime;
    }

    // prepare output file
    output.open( (param.output_dir + "/" + name + ".dat").c_str() );
};
void Species::probe_collect(t_particle *I)
{
    if(!(I->z > 395e-3 && I->r < 45e-3)) return;
    //double center_r = p_param->r_max/2.0;
    //double center_z = p_param->z_max/2.0;
    probe_energy_dist.add((SQR(I->vr)+SQR(I->vz)+SQR(I->vt))*mass*0.5/p_param->q_e);
    //compute incidence angle
    //normal vector: n = (r-x0, z-y0)
    // n.v = |n|*|v|*cos(alpha)
    // (alpha) =-acos(n.v/(|n|*|v|))
    //double nx = center_r-I->r;
    //double ny = center_z-I->z;
    //double absV, alpha;
    //double absN = norm(nx,ny);
    //absV = norm(I->vr, I->vz);
    //alpha = acos( ( nx*I->vr + ny*I->vz )/(absN*absV) );
    //probe_angular_dist.add(abs(alpha));

    //absV = norm(I->vr,I->vz,I->vt);
    //alpha = acos( ( nx*I->vr + ny*I->vz )/(absN*absV) );
    //probe_angular_dist.add(abs(alpha));
    //probe_angular_normalized_dist.add(abs(alpha),1.0/sin(alpha));
    //cerr << (SQR(I->vr)+SQR(I->vz)+SQR(I->vz))*mass*0.5/p_param->q_e << " ";;
    //cerr << (I->vr)<<" "<<(I->vz)<<" "<<(I->vz) << endl;;

    probe_current += charge;
    probe_charge += charge;
}

void Species::add_particles_on_disk(int nparticles, double centerx, double centery, double radius)
{
    double x, y;
    for(int i=0; i<nparticles; i++)
    {
        do{
            x = (rnd->uni() - 0.5);
            y = (rnd->uni() - 0.5);
        }while(SQR(x) + SQR(y) > 0.25 );
        x = x*2*radius + centerx;
        y = y*2*radius + centery;
        if(x < 0 || x > p_param->r_max || y < 0 || y > p_param->z_max) continue;

        int ii = insert();
        t_particle * pp = &(particles[ii]);
        pp->r = x;
        pp->z = y;
        rndv(pp->vr, pp->vz, pp->vt);
    }
};

void Species::add_particle_beam_on_disk(int nparticles, double centerx, double centery, double radius)
{
    double x, y;
    for(int i=0; i<nparticles; i++)
    {
        do{
            x = (rnd->uni() - 0.5);
            y = (rnd->uni() - 0.5);
        }while(SQR(x) + SQR(y) > 0.25 );
        x = x*2*radius + centerx;
        y = y*2*radius + centery;
        if(x < 0 || x > p_param->r_max || y < 0 || y > p_param->z_max) continue;

        int ii = insert();
        t_particle * pp = &(particles[ii]);
        pp->r = x;
        pp->z = y;
        pp->vr = pp->vz = 0.0;
        rndv(pp->vt);
    }
};
void Species::add_particle_beam_on_disk_cylindrical(int nparticles, double centerz, double radius)
{
    double x, y, r;
    if(centerz < 0 || centerz > p_param->z_max) return;

    for(int i=0; i<nparticles; i++)
    {
        //silly implementation coming from variant of this method in cartesian coords
        do{
            x = (rnd->uni() - 0.5);
            y = (rnd->uni() - 0.5);
            r = SQR(x) + SQR(y);
        }while(r > 0.25 );
        r = sqrt(r)*radius*2;
        if(r > p_param->r_max) continue;

        int ii = insert();
        t_particle * pp = &(particles[ii]);
        pp->r = r;
        pp->z = centerz;
        pp->vr = pp->vt = 0.0;
        rndv(pp->vz);
    }
};
void Species::add_monoenergetic_particles_on_cylinder_cylindrical(int nparticles, double energy, double centerz, double radius, double height)
{
    double x, y, r;
    if(centerz < 0 || centerz > p_param->z_max) return;

    for(int i=0; i<nparticles; i++)
    {
        //silly implementation coming from variant of this method in cartesian coords
        do{
            x = (rnd->uni() - 0.5);
            y = (rnd->uni() - 0.5);
            r = SQR(x) + SQR(y);
        }while(r > 0.25 );
        r = sqrt(r)*radius*2;
        if(r > p_param->r_max) continue;

        int ii = insert();
        t_particle * pp = &(particles[ii]);
        pp->r = r;
        pp->z = centerz + height*(rnd->uni()-0.5);

        double v = veV(energy);
        rnd->rot(v, pp->vr, pp->vz, pp->vt);
    }
};
void Species::advance()
{
    double fr, fz;	//force vector
    double qmdt = charge/mass*dt;	//auxilliary constant
    //double center_r = p_param->r_max/2.0;
    //double center_z = p_param->z_max/2.0;
    //double sqpp = SQR(p_param->probe_radius);
    double prob = 1.0-exp(-dt/lifetime);
    //double prob = dt/lifetime;
    //probe_current = 0;

    double Bz = 0.03, Br = 0, Bt = 0.0;
    //double tan_theta = charge*B*dt/(2.0*mass);


    if(p_param->coord == CYLINDRICAL) 
	for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); ++I)
	{
	    if(I->empty==true) continue;

	    //compute field at (I->r, I->z)
	    field->E(I->r, I->z, fr, fz);
            field->B(I->r, I->z, Br, Bz, Bt);
	    //XXX osetrit castice mimo prac oblast !!!
	    //fr=fz=0;


	    // advance velocities as in cartesian coords
	    // use HARHA in magnetic field
	    // half acceleration:
	    I->vr -= fr*qmdt/2.0;
	    I->vz -= fz*qmdt/2.0;

	    // rotation
	    //use Boris' algorithm (Birdsall & Langdon pp. 62) for arbitrary B direction
	    double tmp = charge*dt/(2.0*mass);
	    double tr = Br*tmp;
	    double tt = Bt*tmp;
	    double tz = Bz*tmp;

	    //XXX check orientation
	    // use (r, theta, z)
	    double vprime_r = I->vr + I->vt*tz - I->vz*tt;
	    double vprime_t = I->vt + I->vz*tr - I->vr*tz;
	    double vprime_z = I->vz + I->vr*tt - I->vt*tr;

	    tmp = 2.0/(1+SQR(tr)+SQR(tt)+SQR(tz));
	    double sr = tr*tmp;
	    double st = tt*tmp;
	    double sz = tz*tmp;

	    I->vr = I->vr + vprime_t*sz - vprime_z*st;
	    I->vt = I->vt + vprime_z*sr - vprime_r*sz;
	    I->vz = I->vz + vprime_r*st - vprime_t*sr;

	    // half acceleration:
	    I->vr -= fr*qmdt/2.0;
	    I->vz -= fz*qmdt/2.0;

	    //advance position (Birdsall pp. 338):
	    double x2 = I->r + I->vr*dt;
	    double y2 = I->vt*dt;
	    I->r = sqrt(SQR(x2)+SQR(y2));
	    I->z += I->vz * dt;

	    //rotate the speed vector
	    double sa = y2/I->r;
	    double ca = x2/I->r;
	    if(I->r==0)
	    {
		sa = 0;
		ca = 1;
	    }
	    tmp = I->vr;
	    I->vr = ca*I->vr + sa*I->vt;
	    I->vt = -sa*tmp + ca*I->vt;

	    if( rnd->uni() < prob)
	    {
		scatter(*I);
	    }

	    // OKRAJOVE PODMINKY
	    if(I->r > p_param->r_max || I->r < 0 
		    || I->z > p_param->z_max || I->z < 0 )
	    {
		remove(I);
		continue;
	    }
	    if(!field->grid.is_free(I->r, I->z))
	    {
                if(p_param->has_probe)
                    probe_collect(&*I);

                remove(I);
		continue;
	    }

	    // SUMACE NABOJE
	    rho.accumulate(charge, I->r, I->z);
	}
    else
    {
	for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); ++I)
	{
	    // (r, t, z) ~ (x, y, z)
	    if(I->empty==true) continue;

	    //compute field at (I->r, I->z)
	    field->E(I->r, I->z, fr, fz, niter*dt);
            field->B(I->r, I->z, Br, Bz, Bt);
	    //XXX osetrit castice mimo prac oblast !!!
	    //fr=fz=0;


	    // advance velocities as in cartesian coords
	    // use HARHA in magnetic field
	    // half acceleration:
	    I->vr -= fr*qmdt/2.0;
	    I->vz -= fz*qmdt/2.0;

	    // rotation
	    //use Boris' algorithm (Birdsall & Langdon pp. 62) for arbitrary B direction
	    double tmp = charge*dt/(2.0*mass);
	    double tr = Br*tmp;
	    double tt = Bt*tmp;
	    double tz = Bz*tmp;

	    //XXX check orientation
	    // use (r, theta, z)
	    // use (x, y, z) ~ (r, z, theta) !!!
	    // sign change
	    double vprime_r = I->vr - I->vt*tz + I->vz*tt;
	    double vprime_z = I->vz - I->vr*tt + I->vt*tr;
	    double vprime_t = I->vt - I->vz*tr + I->vr*tz;

	    tmp = 2.0/(1+SQR(tr)+SQR(tt)+SQR(tz));
	    double sr = tr*tmp;
	    double st = tt*tmp;
	    double sz = tz*tmp;

	    I->vr = I->vr - vprime_t*sz + vprime_z*st;
	    I->vt = I->vt - vprime_z*sr + vprime_r*sz;
	    I->vz = I->vz - vprime_r*st + vprime_t*sr;

	    // half acceleration:
	    I->vr -= fr*qmdt/2.0;
	    I->vz -= fz*qmdt/2.0;

	    //advance position (Birdsall pp. 338):
	    I->r += I->vr*dt;
	    I->z += I->vz*dt;


	    if( rnd->uni() < prob)
	    {
		scatter(*I);
	    }
	    // OKRAJOVE PODMINKY

	    //if(p_param->selfconsistent)
	    //{
	    if(I->r > p_param->r_max || I->r < 0 
		    || I->z > p_param->z_max || I->z < 0 )
	    {
                if(p_param->boundary == Param::FREE)
                {
                    remove(I);
                    continue;
                }
                else if(p_param->boundary == Param::PERIODIC)
                {
                    I->r = fmod(I->r, p_param->r_max);
                    if(I->r < 0) I->r += p_param->r_max;
                    I->z = fmod(I->z, p_param->z_max);
                    if(I->z < 0) I->z += p_param->z_max;
                }

	    }
	    if(!field->grid.is_free(I->r, I->z))
	    {
		remove(I);
		continue;
	    }

	    /**
	      if(SQR(I->r-center_r)+SQR(I->z-center_z) < sqpp)
	      {
	      probe_collect(&*I);
	      remove(I);
	      probe_current += charge;
	      continue;
	      }
	      */

	    // SUMACE NABOJE
	    rho.accumulate(charge, I->r, I->z);
	    //}
	}
    }
    niter++;
}
void Species::accumulate()
{
    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); ++I)
    {
	if(I->empty==true) continue;
	    // SUMACE NABOJE
	    rho.accumulate(charge, I->r, I->z);
    }
}
void Species::source5_refresh(unsigned int factor)
{
    source5_factor = factor;
    //double N = p_param->rho*p_param->V;		//XXX wrong for multicomponent plasma
    //double N = p_param->density[type]*p_param->V;		//XXX wrong for multicomponent plasma
    double N = density*p_param->V;
    unsigned int n_particles = N/factor;
    double K = 1.0/factor;

    // the second comparison checks if we simulate this species as particles
    if(source2_particles.size() != n_particles && particles.size() != 0)
	source2_particles.resize(n_particles);

    for(unsigned int i=0; i<source2_particles.size(); i++)
    {
	source2_particles[i].empty = false;
	source2_particles[i].r = K*p_param->r_max*rnd->uni();
	source2_particles[i].z = K*p_param->z_max*rnd->uni();
	source2_particles[i].vr = rnd->rnor()*v_max/(M_SQRT2);
	source2_particles[i].vz = rnd->rnor()*v_max/(M_SQRT2);
	source2_particles[i].vt = rnd->rnor()*v_max/(M_SQRT2);
	source2_particles[i].time_to_death = rnd->rexp()*lifetime;
    }
}


void Species::energy_dist_compute()
{
    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); I++)
	if(I->empty==false)
	    if( energy_dist.add((SQR(I->vr)+SQR(I->vz)+SQR(I->vt))*mass*0.5/p_param->q_e) == -1 )
		;//cout << I-particles.begin() << " " << I->vr << " " << I->vz << " " << I->vz << endl;
}
void Species::source_energy_dist_compute()
{
    for(vector<t_particle>::iterator I = source2_particles.begin(); I != source2_particles.end(); I++)
	if(I->empty==false)
	    source_energy_dist.add((SQR(I->vr)+SQR(I->vz)+SQR(I->vt))*mass*0.5/p_param->q_e);
		;//cout << I-source2_particles.begin() << " " << I->vr << " " << I->vz << " " << I->vz << endl;
}










/*
 * not working currently
 */



void Species::source_old()
{

    //double v_max = sqrt(2.0*t_param::k_B*temperature/mass);
    //double f2 = p_param->rho*M_2_SQRTPI*0.25*v_max*dt*p_param->dz;
    double f2 = density/(2*sqrt(M_PI))*v_max*dt*p_param->dz;
    int pocet;
    int k;
    double f1;
    double u;
    vector<t_particle>::iterator I;

    for(k=0; k<4; k++)
    {
	if(k<2) //r==konst
	    f1 = f2*p_param->z_max;
	else
	    f1 = f2*p_param->r_max;
	//cerr << f1*4 << endl;

	//if(f1 > 10 )
	//{
	    if((pocet = (int)( f1 + sqrt(f1)*rnd->rnor() + .5)) < 0)
		pocet = 0;
	//}else {
	    pocet = (int)f1;
	    pocet += (rnd->uni()>f1-pocet) ? 0 : 1;
	//}
	for(int i=0; i<pocet; i++)
	{
	    I = particles.begin() + insert();

	    if(k<2)
	    {
		I->vr = rnd->rnor()*v_max/(M_SQRT2);
		I->vz = rnd->maxwell_flux(v_max*5.0) * (k==0 ? 1 : -1);
		//I->r = I->vr*dt*uniform() + p_param->r_max * (k==0 ? 0 : 1);
		//I->z = p_param->z_max*uniform();
		u = rnd->uni();
		I->r = -I->vr*dt*u + p_param->r_max * (k==0 ? 0 : 1);
		I->z = p_param->z_max*rnd->uni() - I->vz*dt*u;
	    }else
	    {
		I->vr = rnd->rnor()*v_max/(M_SQRT2);
		I->vz = rnd->maxwell_flux(v_max*5.0) * (k==2 ? 1 : -1);
		//I->z = I->vz*dt*uniform() + p_param->z_max * (k==2 ? 0 : 1);
		//I->r = p_param->r_max*uniform();
		u = rnd->uni();
		I->r = -I->vz*dt*u + p_param->z_max * (k==2 ? 0 : 1);
		I->z = p_param->r_max*rnd->uni() - I->vr*dt*u;
	    }
	    I->vz = rnd->rnor()*v_max/M_SQRT2;
	    I->time_to_death = rnd->rexp()*lifetime;

	    I->r += I->vr * dt;
	    I->z += I->vz * dt;
	    if(I->r < p_param->r_max && I->r > 0 
		    && I->z < p_param->z_max && I->z > 0 )
		rho.accumulate(charge, I->r, I->z);
	    /*
	    for(int jj=0; jj<2; jj++)
	    {	
		I->time_to_death -= dt;
		if(I->time_to_death < 0)
		    I->time_to_death += exponential(lifetime);
	    }
	    */

	}
    }
}

void Species::source()
{
    double K = 1.0/source5_factor;
    vector<t_particle>::iterator J;
    int j;
    //double qm2 = charge/mass*0.5;
    double qm = charge/mass;
    double src_z_max = K*p_param->z_max;
    double src_r_max = K*p_param->r_max;
    double fx=-p_param->extern_field;
    double fy=0;
    double qmdt = qm*dt;

    //cout << "source\n";
    //cout << "diff: " << source2_particles.begin() - source2_particles.end() << endl; 
    for(vector<t_particle>::iterator I=source2_particles.begin(); I != source2_particles.end(); I++)
    {

	I->vr -= fx*qmdt;
	I->vz -= fy*qmdt;
	I->r += I->vr * dt;
	I->z += I->vz * dt;

	if( I->time_to_death < dt)
	{
	    scatter(*I);
	    I->time_to_death += rnd->rexp()*lifetime;
	}
	I->time_to_death -= dt;

	if(I->r > src_r_max)
	    while(I->r > src_r_max)
	    {
		I->r -= src_r_max;
		j = insert();
		particles[j] = *I;
		particles[j].z += rand()%source5_factor*src_z_max;
		if(particles[j].z<p_param->z_max && particles[j].z>0 &&
		       	particles[j].r<p_param->r_max && particles[j].r>0)
		    rho.accumulate(charge, particles[j].r, particles[j].z);
		else
		    remove(j);
	    }
	else if(I->r < 0)
	    while(I->r < 0)
	    {
		j = insert();
		particles[j] = *I;
		particles[j].z += rand()%source5_factor*src_z_max;
		particles[j].r += p_param->r_max;
		I->r += src_r_max;
		if(particles[j].z<p_param->z_max && particles[j].z>0 &&
		       	particles[j].r<p_param->r_max && particles[j].r>0)
		    rho.accumulate(charge, particles[j].r, particles[j].z);
		else
		    remove(j);
	    }
	if(I->z > src_z_max)
	    while(I->z > src_z_max)
	    {
		I->z -= src_z_max;
		j = insert();
		particles[j] = *I;
		particles[j].r += rand()%source5_factor*src_r_max;
		if(particles[j].z<p_param->z_max && particles[j].z>0 &&
		       	particles[j].r<p_param->r_max && particles[j].r>0)
		    rho.accumulate(charge, particles[j].r, particles[j].z);
		else
		    remove(j);
	    }
	else if(I->z < 0)
	    while(I->z < 0)
	    {
		j = insert();
		particles[j] = *I;
		particles[j].r += rand()%source5_factor*src_r_max;
		particles[j].z += p_param->z_max;
		I->z += src_z_max;
		if(particles[j].z<p_param->z_max && particles[j].z>0 &&
		       	particles[j].r<p_param->r_max && particles[j].r>0)
		    rho.accumulate(charge, particles[j].r, particles[j].z);
		else
		    remove(j);
	    }
    }
}
#endif
