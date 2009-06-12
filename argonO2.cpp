#ifndef ARGON_H
#define ARGON_H

#include "particles.cpp"
#include "random.cpp"
#include "tabulate.cpp"
namespace argon{
static const double E_MIN=1e-3;
static const double E_MAX=1.0;
};

class t_argon_neutral : public Species
{
    void scatter(t_particle &particle){};
    public:
    t_argon_neutral(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, Species *_species_list[], double _dt=1e-8)
	: Species(0, 0, param, random, _field, 6.68173e-26, _species_list, ARGON)
    {
	charge = 0.0;
//	double pressure = 0.0;
//	density = pressure/(param.k_B*temperature);
    };
};

class t_argon_metastable : public Species
{
    public:
    t_argon_metastable(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, Species *_species_list[], double _dt=1e-8)
	: Species(0, 0, param, random, _field, 6.68173e-26, _species_list, ARGON_META)
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



double sigma2_in_el(double E)
{
    return (2e-19/(sqrt(E)*(1+E)) + 3E-19*E/sqr(1+E/3))*1e20;
}

double sigma2_in_ct(double E)
{
    return ((1.15E-18*pow(E,-0.1)*pow(1+0.015/E,0.6))*1e20-sigma2_in_el(E))*0.5;
}
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

class t_argon : public Species
{
    //mass = 5e-26;
    //charge = 1.6e-19;
    void scatter(t_particle &particle);
    double ArArsvmax;
    double ArArFreq;
    double ArO2svmax;
    double ArO2Freq;
    enum coll_type { ELASTIC, CX };
    vector<vec_interpolate*> ArArCS;
    vector<double> ArArLoss;
    vector<coll_type> ArArType;
    //double O2O2Freq;
    //double O2O2svmax;

    public:
    t_argon(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, Species *_species_list[], double _dt=1e-8)
	: Species(_n, n2, param, random, _field, 6.68173e-26, _species_list, ARGON_POS)
    {
	ArArsvmax = 1.22e-15;
	ArO2svmax = 2e-15;	//XXX just a random no-nonsense number

	charge = 1.602189e-19;

	ArArCS.push_back(new vec_interpolate(sigma2_in_el,argon::E_MIN,argon::E_MAX,200) );
	ArArType.push_back(ELASTIC);

	ArArCS.push_back(new vec_interpolate(sigma2_in_ct,argon::E_MIN,argon::E_MAX,200) );
	ArArType.push_back(CX);

	cout << "t_argon: veV(E_MAX) = "<<veV(argon::E_MAX) <<' '<< argon::E_MAX<<endl;
	ArArsvmax = svmax_find(ArArCS,veV(argon::E_MAX));

	cout << "ArArsvmax = "<< ArArsvmax <<endl;

	if(dt==0)
	    dt = min(lifetime/10.0,1e-8);
	cerr << "Ar+ dt = "<< dt << " ; " << lifetime/dt << " timesteps per collision" << endl;
    };
    void lifetime_init()
    {
	ArO2Freq = ArO2svmax * species_list[O2]->density;
	ArArFreq = ArArsvmax * species_list[ARGON]->density;
	cout << species_list[ARGON]->name << ' '<< species_list[ARGON_POS]->density << ' '<< density<<endl;
	lifetime = 1.0/(ArO2Freq + ArArFreq);
	cout << "argon lifetime " << lifetime/dt <<endl;
	Species::lifetime_init();
    }
};




void t_argon::scatter(t_particle &particle)
{

    double const_E = 0.5*mass/charge;

    //druha interagujici castice
    double vr2, vz2, vt2;
    species_list[ARGON]->rndv(vr2, vz2, vt2);

    double sq_v_rel = SQR(vr2-particle.vx) + SQR(vz2-particle.vz) + SQR(vt2-particle.vy);
    double v = sqrt(sq_v_rel);
    double E = const_E * sq_v_rel; //[eV]


    double gamma = rnd->uni() * ArArsvmax * 1e20;
    double tmp = 0;
    unsigned int i;
    for(i=0; i<ArArCS.size(); i++)
    {
	tmp += (*ArArCS[i])(E)*v;
	if(tmp > gamma)
	    break;
    }


    if(i==ArArCS.size()) return; // NULL collision
    switch(ArArType[i])
    {
	//charge transfer
	case CX:
	    particle.vx = vr2;
	    particle.vz = vz2;
	    particle.vy = vt2;
	    break;

	case ELASTIC:
	    {
		//prepocet v_1 do tezistove soustavy
		double v1_cm_x = (particle.vx - vr2)*0.5;
		double v1_cm_y = (particle.vz - vz2)*0.5;
		double v1_cm_z = (particle.vy - vt2)*0.5;
		double v1_cm = sqrt(SQR(v1_cm_x) + SQR(v1_cm_y) + SQR(v1_cm_z));

		//provedeni nahodne rotace
		rnd->rot(v1_cm,v1_cm_x,v1_cm_y,v1_cm_z);

		//zpetna transformace
		//particle.vx = v1_cm_x + v_cm_x;
		particle.vx = v1_cm_x + (particle.vx + vr2)*0.5;
		particle.vz = v1_cm_y + (particle.vz + vz2)*0.5;
		particle.vy = v1_cm_z + (particle.vy + vt2)*0.5;
	    }
	    break;

	default:
	    break;
    }
}


#endif
