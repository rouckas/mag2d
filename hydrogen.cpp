#ifndef HYDROGEN_H
#define HYDROGEN_H

#include "particles.cpp"
#include "random.cpp"
#include "tabulate.cpp"
namespace hydrogen{
static const double E_MIN=1e-3;
static const double E_MAX=1.0;
};

class t_hydrogen_neutral : public Species
{
    void scatter(t_particle &particle){};
    public:
    t_hydrogen_neutral(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, Species *_species_list[], double _dt=1e-8)
	: Species(0, 0, param, random, _field, 6.68173e-26, _species_list, HYDROGEN)
    {
	charge = 0.0;
//	double pressure = 0.0;
//	density = pressure/(param.k_B*temperature);
    };
};

double h_sigma2_in_el(double E)
{
    return (2e-19/(sqrt(E)*(1+E)) + 3E-19*E/sqr(1+E/3))*1e20;
}

double h_sigma2_in_ct(double E)
{
    return ((1.15E-18*pow(E,-0.1)*pow(1+0.015/E,0.6))*1e20-sigma2_in_el(E))*0.5;
}
double h_sigma_in_el(double E)
{
    double mass = 6.68173e-26;
    double charge = 1.602189e-19;
    double v_rel = sqrt(E*charge*2.0/mass);
    return (2e-19/(sqrt(E)*(1+E)) + 3E-19*E/sqr(1+E/3))*v_rel;
}

double h_sigma_in_ct(double E)
{
    double mass = 6.68173e-26;
    double charge = 1.602189e-19;
    double v_rel = sqrt(E*charge*2.0/mass);
    return ((1.15E-18*pow(E,-0.1)*pow(1+0.015/E,0.6))*v_rel-sigma_in_el(E))*0.5;
}
tabulate h_sigma_in_el_tab(0, argon::E_MAX, 1000, sigma_in_el);
tabulate h_sigma_in_ct_tab(0, argon::E_MAX, 1000, sigma_in_ct);

class t_hydrogen_neg : public Species
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
    t_hydrogen_neg(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, Species *_species_list[], double _dt=1e-8)
	: Species(_n, n2, param, random, _field, 1.66e-27, _species_list, H_NEG)
    {
	ArArsvmax = 1.22e-15;
	ArO2svmax = 2e-15;	//XXX just a random no-nonsense number

	charge = -1.602189e-19;

	ArArCS.push_back(new vec_interpolate(sigma2_in_el,argon::E_MIN,argon::E_MAX,200) );
	ArArType.push_back(ELASTIC);

	ArArCS.push_back(new vec_interpolate(sigma2_in_ct,argon::E_MIN,argon::E_MAX,200) );
	ArArType.push_back(CX);

	cout << "t_hydrogen_neg: veV(E_MAX) = "<<veV(argon::E_MAX) <<' '<< argon::E_MAX<<endl;
	ArArsvmax = svmax_find(ArArCS,veV(argon::E_MAX));

	cout << "ArArsvmax = "<< ArArsvmax <<endl;

	if(dt==0)
	    dt = min(lifetime/10.0,1e-8);
    };
    void lifetime_init()
    {
	ArO2Freq = 0;
	ArArFreq = ArArsvmax * species_list[ARGON]->density;
	cout << species_list[ARGON]->name << ' '<< density<<endl;
	lifetime = 1.0/(ArO2Freq + ArArFreq);
	cout << "hydrogen lifetime " << lifetime/dt <<endl;
	Species::lifetime_init();
    }
};




void t_hydrogen_neg::scatter(t_particle &particle)
{

    double const_E = 0.5*mass/charge;

    //druha interagujici castice
    double vr2, vz2, vt2;
    species_list[ARGON]->rndv(vr2, vz2, vt2);

    double sq_v_rel = SQR(vr2-particle.vr) + SQR(vz2-particle.vz) + SQR(vt2-particle.vt);
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
	    particle.vr = vr2;
	    particle.vz = vz2;
	    particle.vt = vt2;
	    break;

	case ELASTIC:
	    {
		//prepocet v_1 do tezistove soustavy
		double v1_cm_x = (particle.vr - vr2)*0.5;
		double v1_cm_y = (particle.vz - vz2)*0.5;
		double v1_cm_z = (particle.vt - vt2)*0.5;
		double v1_cm = sqrt(SQR(v1_cm_x) + SQR(v1_cm_y) + SQR(v1_cm_z));

		//provedeni nahodne rotace
		rnd->rot(v1_cm,v1_cm_x,v1_cm_y,v1_cm_z);

		//zpetna transformace
		//particle.vr = v1_cm_x + v_cm_x;
		particle.vr = v1_cm_x + (particle.vr + vr2)*0.5;
		particle.vz = v1_cm_y + (particle.vz + vz2)*0.5;
		particle.vt = v1_cm_z + (particle.vt + vt2)*0.5;
	    }
	    break;

	default:
	    break;
    }
}


#endif
