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
	: Species(0, 0, param, random, _field, 1.67262158e-27, _species_list, HYDROGEN)
    {
	charge = 0.0;
    };
};

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
tabulate h_sigma_in_el_tab(0, hydrogen::E_MAX, 1000, sigma_in_el);
tabulate h_sigma_in_ct_tab(0, hydrogen::E_MAX, 1000, sigma_in_ct);

class t_hydrogen_neg : public Species
{
    void scatter(t_particle &particle);
    double HHsvmax;
    double HHFreq;
    double HHesvmax;
    double HHeFreq;
    enum coll_type { ELASTIC, CX };
    vector<vec_interpolate*> HHCS;
    vector<double> HHLoss;
    vector<coll_type> HHType;
    vector<vec_interpolate*> HHeCS;
    vector<coll_type> HHeType;
    vector<double> HHk, HHek;

    public:
    t_hydrogen_neg(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, Species *_species_list[], double _dt=1e-8)
	: Species(_n, n2, param, random, _field, 1.67262158e-27, _species_list, H_NEG)
    {
	HHsvmax = 1.22e-15;
	HHesvmax = 1.22e-15;	//XXX just a random no-nonsense number

	charge = -1.602189e-19;

	HHCS.push_back(new vec_interpolate(sigma2_in_el,hydrogen::E_MIN,hydrogen::E_MAX,200) );
        HHk.push_back(1.18e-15);
	HHType.push_back(ELASTIC);

	HHCS.push_back(new vec_interpolate(sigma2_in_ct,hydrogen::E_MIN,hydrogen::E_MAX,200) );
        HHk.push_back(0.0);
	HHType.push_back(CX);

        //initialize the helium collisions same as hydrogen, we don't know it exactly anyway
	HHeCS.push_back(new vec_interpolate(sigma2_in_el,hydrogen::E_MIN,hydrogen::E_MAX,200) );
        HHek.push_back(1.18e-15);
	HHeType.push_back(ELASTIC);

	HHeCS.push_back(new vec_interpolate(sigma2_in_ct,hydrogen::E_MIN,hydrogen::E_MAX,200) );
        HHek.push_back(0.0);
	HHeType.push_back(CX);


        HHsvmax = 0;
        for(unsigned int i=0; i<HHk.size(); i++) HHsvmax += HHk[i];
        HHesvmax = 0;
        for(unsigned int i=0; i<HHek.size(); i++) HHesvmax += HHek[i];
	//HHsvmax = svmax_find(HHCS,veV(hydrogen::E_MAX));
	//HHesvmax = svmax_find(HHeCS,veV(hydrogen::E_MAX));


	if(dt==0)
	    dt = min(lifetime/10.0,1e-8);
    };
    void lifetime_init()
    {
        HHeFreq = HHesvmax * species_list[HELIUM]->density;
	HHFreq = HHsvmax * species_list[HYDROGEN]->density;
	cout << species_list[HYDROGEN]->name << ' '<< density<<endl;
	lifetime = 1.0/(HHeFreq + HHFreq);
	cout << "hydrogen lifetime " << lifetime/dt <<endl;
	Species::lifetime_init();
    }

    void scatter(t_particle &particle, int species, double svmax, vector<vec_interpolate*> & CCS, vector<coll_type> & CType);
    void scatter_k(t_particle &particle, int species, double svmax, vector<double> & Ck, vector<coll_type> & CType);
};




void t_hydrogen_neg::scatter(t_particle &particle)
{
    double freq = HHFreq + HHeFreq;
    if(rnd->uni() < HHFreq/freq)
    {
        //scatter(particle, HYDROGEN, HHsvmax, HHCS, HHType);
        scatter_k(particle, HYDROGEN, HHsvmax, HHk, HHType);
    }
    else
    {
        scatter_k(particle, HELIUM, HHesvmax, HHk, HHeType);
    }

};

void t_hydrogen_neg::scatter_k(t_particle &particle, int species, double svmax, vector<double> & Ck, vector<coll_type> & CType)
{
    double vr2,vz2,vt2;
    species_list[species]->rndv(vr2,vz2,vt2);


    // choose the collisional process randomly
    double gamma = rnd->uni() * svmax;
    double tmp=0;
    unsigned int i;
    for(i=0; i<Ck.size(); i++)
    {
        tmp += Ck[i];
        if(tmp > gamma)
            break;
    }

    // select apropriate collision type
    switch(CType[i])
    {

        case ELASTIC:
            {
                double m2 = species_list[species]->mass;
                double tmp = 1.0/(mass+m2);
                //prepocet v_1 do tezistove soustavy
                double v1_cm_x = (particle.vx - vr2)*m2*tmp;
                double v1_cm_y = (particle.vz - vz2)*m2*tmp;
                double v1_cm_z = (particle.vy - vt2)*m2*tmp;
                //double v1_cm = sqrt(SQR(v1_cm_x) + SQR(v1_cm_y) + SQR(v1_cm_z));

                //provedeni nahodne rotace
                rnd->rot(v1_cm_x,v1_cm_y,v1_cm_z);

                //zpetna transformace
                //particle.vx = v1_cm_x + v_cm_x;
                particle.vx = v1_cm_x + (particle.vx*mass + vr2*m2)*tmp;
                particle.vz = v1_cm_y + (particle.vz*mass + vz2*m2)*tmp;
                particle.vy = v1_cm_z + (particle.vy*mass + vt2*m2)*tmp;
            }
            break;
        case CX :
            particle.vx = vr2;
            particle.vz = vz2;
            particle.vy = vt2;
            break;
    }
};

void t_hydrogen_neg::scatter(t_particle &particle, int species, double svmax, vector<vec_interpolate*> & CCS, vector<coll_type> & CType)
{
    double vr2,vz2,vt2;
    species_list[species]->rndv(vr2,vz2,vt2);

    double v = norm(particle.vx-vr2, particle.vz-vz2, particle.vy-vt2);
    double const_E = -0.5*mass/charge;
    double E = const_E*v*v;

    // choose the collisional process randomly
    double gamma = rnd->uni() * svmax*1e20;
    double tmp=0;
    unsigned int i;
    for(i=0; i<CCS.size(); i++)
    {
        tmp += (*CCS[i])(E)*v;
        if(tmp > gamma)
            break;
    }

    // select apropriate collision type
    switch(CType[i])
    {

        case ELASTIC:
            {
                double m2 = species_list[species]->mass;
                double tmp = 1.0/(mass+m2);
                //prepocet v_1 do tezistove soustavy
                double v1_cm_x = (particle.vx - vr2)*m2*tmp;
                double v1_cm_y = (particle.vz - vz2)*m2*tmp;
                double v1_cm_z = (particle.vy - vt2)*m2*tmp;
                //double v1_cm = sqrt(SQR(v1_cm_x) + SQR(v1_cm_y) + SQR(v1_cm_z));

                //provedeni nahodne rotace
                rnd->rot(v1_cm_x,v1_cm_y,v1_cm_z);

                //zpetna transformace
                //particle.vx = v1_cm_x + v_cm_x;
                particle.vx = v1_cm_x + (particle.vx*mass + vr2*m2)*tmp;
                particle.vz = v1_cm_y + (particle.vz*mass + vz2*m2)*tmp;
                particle.vy = v1_cm_z + (particle.vy*mass + vt2*m2)*tmp;
            }
            break;
        case CX :
            particle.vx = vr2;
            particle.vz = vz2;
            particle.vy = vt2;
            break;
    }
};

#endif
