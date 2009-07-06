#ifndef ARGON_H
#define ARGON_H

#include "particles.cpp"
#include "random.cpp"
#include "tabulate.cpp"
namespace argon{
static const double E_MIN=1e-3;
static const double E_MAX=1.0;
};

template <int D>
class t_argon_neutral : public Species<D>
{
    void scatter(t_particle &particle){};
    public:
    t_argon_neutral(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, BaseSpecies *_species_list[], double _dt=1e-8)
	: Species<D>(0, 0, param, random, _field, 6.68173e-26, _species_list, ARGON)
    {
	Species<D>::charge = 0.0;
//	double pressure = 0.0;
//	density = pressure/(param.k_B*temperature);
    };
};

template <int D>
class t_argon_metastable : public Species<D>
{
    public:
    t_argon_metastable(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, BaseSpecies *_species_list[], double _dt=1e-8)
	: Species<D>(0, 0, param, random, _field, 6.68173e-26, _species_list, ARGON_META)
    {
	Species<D>::charge = 0.0;
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

template <int D>
class t_argon : public Species<D>
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
	    Fields *_field, BaseSpecies *_species_list[], double _dt=1e-8)
	: Species<D>(_n, n2, param, random, _field, 6.68173e-26, _species_list, ARGON_POS)
    {
	ArArsvmax = 1.22e-15;
	ArO2svmax = 2e-15;	//XXX just a random no-nonsense number

	Species<D>::charge = 1.602189e-19;

	ArArCS.push_back(new vec_interpolate(sigma2_in_el,argon::E_MIN,argon::E_MAX,200) );
	ArArType.push_back(ELASTIC);

	ArArCS.push_back(new vec_interpolate(sigma2_in_ct,argon::E_MIN,argon::E_MAX,200) );
	ArArType.push_back(CX);

	cout << "t_argon: veV(E_MAX) = "<< Species<D>::veV(argon::E_MAX) <<' '<< argon::E_MAX<<endl;
	ArArsvmax = svmax_find(ArArCS, Species<D>::veV(argon::E_MAX));

	cout << "ArArsvmax = "<< ArArsvmax <<endl;

	if(Species<D>::dt==0)
	    Species<D>::dt = min(Species<D>::lifetime/10.0,1e-8);
	cerr << "Ar+ dt = "<< Species<D>::dt << " ; " << Species<D>::lifetime/Species<D>::dt << " timesteps per collision" << endl;
    };
    void lifetime_init()
    {
	ArO2Freq = ArO2svmax * Species<D>::species_list[O2]->density;
	ArArFreq = ArArsvmax * Species<D>::species_list[ARGON]->density;
	cout << Species<D>::species_list[ARGON]->name << ' '<< Species<D>::species_list[ARGON_POS]->density << ' '<< Species<D>::density<<endl;
	Species<D>::lifetime = 1.0/(ArO2Freq + ArArFreq);
	cout << "argon lifetime " << Species<D>::lifetime/Species<D>::dt <<endl;
	Species<D>::lifetime_init();
    }
};




template <int D>
void t_argon<D>::scatter(t_particle &particle)
{

    double const_E = 0.5*Species<D>::mass/Species<D>::charge;

    //druha interagujici castice
    double vr2, vz2, vt2;
    Species<D>::species_list[ARGON]->rndv(vr2, vz2, vt2);

    double sq_v_rel = SQR(vr2-particle.vx) + SQR(vz2-particle.vz) + SQR(vt2-particle.vy);
    double v = sqrt(sq_v_rel);
    double E = const_E * sq_v_rel; //[eV]


    double gamma = Species<D>::rnd->uni() * ArArsvmax * 1e20;
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
		Species<D>::rnd->rot(v1_cm,v1_cm_x,v1_cm_y,v1_cm_z);

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
