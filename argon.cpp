#ifndef ARGON_H
#define ARGON_H

#include "particles.cpp"
#include "random.cpp"
#include "tabulate.cpp"
namespace argon{
static const double E_MAX=1.0;
};

class t_argon_neutral : public t_species
{
    void scatter(t_particle &particle){};
    public:
    t_argon_neutral(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, t_species *_species_list[], double _dt=1e-8)
	: t_species(0, 0, temperature, param, random, _field, vsigma_max,  6.68173e-26, _dt, _species_list, ARGON)
    {
	charge = 0.0;
	double pressure = 0.0;
	density = pressure/(param.k_B*temperature);
    };
};

class t_argon_metastable : public t_species
{
    public:
    t_argon_metastable(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, t_species *_species_list[], double _dt=1e-8)
	: t_species(0, 0, temperature, param, random, _field, vsigma_max,  6.68173e-26, _dt, _species_list, ARGON_META)
    {
	charge = 0.0;
    };
    void scatter(t_particle &particle){};
    void advance()
    {
	// TODO insert high energy electrons into the plasma
	// should we insert new electrons or replace them ???? XXX
    };
};



class t_argon : public t_species
{
    //mass = 5e-26;
    //charge = 1.6e-19;
    void scatter(t_particle &particle);
    public:
    t_argon(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, t_species *_species_list[], double _dt=1e-8)
	: t_species(_n, n2, temperature, param, random, _field, vsigma_max,  6.68173e-26, _dt, _species_list, ARGON_POS)
    {
	charge = 1.602189e-19;
	cerr <<"lifetime "<< lifetime << endl;
	if(dt==0)
	    dt = min(lifetime/10.0,1e-8);
	cerr << "Ar+ dt = "<< dt << " ; " << lifetime/dt << " timesteps per collision" << endl;
    };
};


double sigma_in_el(double E)
{
    double mass = 6.68173e-26;
    double charge = 1.602189e-19;
    double v_rel = sqrt(E*charge*2.0/mass);
    return (2e-19/(sqrt(E)*(1+E)) + 3E-19*E/sqr(1+E/3))*v_rel;
}

double sigma_in_ct(double E)
{
    double mass = 6.68173e-26;
    double charge = 1.602189e-19;
    double v_rel = sqrt(E*charge*2.0/mass);
    return ((1.15E-18*pow(E,-0.1)*pow(1+0.015/E,0.6))*v_rel-sigma_in_el(E))*0.5;
}
tabulate sigma_in_el_tab(0, argon::E_MAX, 1000, sigma_in_el);
tabulate sigma_in_ct_tab(0, argon::E_MAX, 1000, sigma_in_ct);


void t_argon::scatter(t_particle &particle)
{

//    rot(particle.vx, particle.vy, particle.vz);
 //   return;
    double const_E = 0.5*mass/charge;

    //druha interagujici castice
    //double vx2 = rnd->rnor()*v_max*M_SQRT1_2;
    //double vy2 = rnd->rnor()*v_max*M_SQRT1_2;
    //double vz2 = rnd->rnor()*v_max*M_SQRT1_2;
    double vx2, vy2, vz2;
    species_list[ARGON]->rndv(vx2, vy2, vz2);

    double sq_v_rel = SQR(vx2-particle.vx) + SQR(vy2-particle.vy) + SQR(vz2-particle.vz);
    double E = const_E * sq_v_rel; //[eV]

    double sigma[3];

    //pruzny rozptyl
    sigma[1] = (sigma_in_el_tab)(E);

    //charge transfer
    sigma[0] = (sigma_in_ct_tab)(E);

    sigma[2] = 1.22e-15 - sigma[0] - sigma[1];

    int i;
    for(i=0;i<3;i++)
	sigma[i]/=1.22e-15;

    double gamma = rnd->uni();
    double tmp=0;
    for(i=0;i<3;i++)
    {
	tmp += sigma[i];
	if(tmp > gamma)
	    break;
    }

    switch(i)
    {
	//charge transfer
	case 0:
	    particle.vx = vx2;
	    particle.vy = vy2;
	    particle.vz = vz2;
	    break;

	case 1:
	    {
		//prepocet v_1 do tezistove soustavy
		double v1_cm_x = (particle.vx - vx2)*0.5;
		double v1_cm_y = (particle.vy - vy2)*0.5;
		double v1_cm_z = (particle.vz - vz2)*0.5;
		double v1_cm = sqrt(SQR(v1_cm_x) + SQR(v1_cm_y) + SQR(v1_cm_z));

		//provedeni nahodne rotace
		rnd->rot(v1_cm,v1_cm_x,v1_cm_y,v1_cm_z);

		//zpetna transformace
		//particle.vx = v1_cm_x + v_cm_x;
		particle.vx = v1_cm_x + (particle.vx + vx2)*0.5;
		particle.vy = v1_cm_y + (particle.vy + vy2)*0.5;
		particle.vz = v1_cm_z + (particle.vz + vz2)*0.5;
	    }
	    break;

	default:
	    break;
    }

}


#endif
