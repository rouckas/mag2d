#ifndef ELONO2_H
#define ELONO2_H
#include "particles.hpp"
#include "tabulate.cpp"
#include "random.cpp"
#include "input.hpp"
#include "mymath.cpp"

namespace elon{
    double l_debye;
    double Lamb;
};


template <int D>
class t_elon : public Species<D>
{
    public:
    void scatter(t_particle &particle);
    private:
    
    enum coll_type { ELASTIC, EXCITATION, IONIZATION };

    // electron cross sections in O2 plasma
    vector<vec_interpolate*> eO2CS;
    vector<double> eO2Loss;
    vector<coll_type> eO2Type;
    double eO2Freq;
    double eO2svmax;

    // electron cross sections in Ar plasma
    vector<vec_interpolate*> eArCS;
    vector<double> eArLoss;
    vector<coll_type> eArType;
    double eArFreq;
    double eArsvmax;

    // electron cross sections in He plasma
    vector<vec_interpolate*> eHeCS;
    vector<double> eHeLoss;
    vector<coll_type> eHeType;
    double eHeFreq;
    double eHesvmax;

    public:
    t_elon(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random,
	   Fields *_field, BaseSpecies * _species_list[], double _dt=1e-11)
	: Species<D>(_n, n2, param, random, _field, 9.109534e-31, _species_list, ELECTRON) 
    {
	Species<D>::charge = -1.602189e-19;
	//double vsigma_coul_max = 6.371e-10;

	elon::Lamb = Species<D>::charge*Species<D>::charge/(4.0*M_PI*param.eps_0);
	elon::l_debye = param.eps_0*param.k_B*param.temperature[ELECTRON]/(Species<D>::density*Species<D>::charge*Species<D>::charge);
	cerr << "e-  "<< Species<D>::lifetime/Species<D>::dt << " timesteps per collision" << endl;

	// Initialization of electron inelastic cross sections in O2:
	// load the O2 cross sections
	vector<double> Edata,CSdata;
	vector<string> CSnames;
	CSnames.push_back("O2 MOMENTUM-TRANSFER");
	eO2Type.push_back(ELASTIC);
	CSnames.push_back("O2 SINGL LEVEL ROT");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 V=1 LINDER");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 V=2 LINDER");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 V=3 LINDER");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 V=4 LINDER");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 V=1  9V");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 V=1  9V");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 4.5 LOSS");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 6.0 LOSS");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 8.4 LOSS");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 9.97 LOSS");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 SING DELTA");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 B SINGLET SIGMA");
	eO2Type.push_back(EXCITATION);
	CSnames.push_back("O2 IONIZATION");
	eO2Type.push_back(IONIZATION);
	//CSnames.push_back("dissociative attachment");
	
	Species<D>::load_CS("electron.txt", eO2CS, eO2Loss, CSnames);

	// compute the collision frequency with O2
	eO2svmax = svmax_find(eO2CS,Species<D>::v_max*5.0);
	cout << "eO2svmax = "<< eO2svmax<<endl;

	
	// Initialization of electron inelastic cross sections in Ar:
	// load the Ar cross sections
	CSnames.clear();
	CSnames.push_back("ELASTIC MOMENTUM-TRANSFER CROSS");
	eArType.push_back(ELASTIC);
	CSnames.push_back("TOTAL EXCITATION LOW E");
	eArType.push_back(EXCITATION);
	CSnames.push_back("AR IONIZATION");
	eArType.push_back(IONIZATION);
	
	Species<D>::load_CS("electron.txt", eArCS, eArLoss, CSnames,"ARGON");
	
	// compute the collision frequency with Ar
	eArsvmax = svmax_find(eArCS,Species<D>::v_max*5.0);
	cout << "eArsvmax = "<< eArsvmax<<endl;

	// Initialization of electron inelastic cross sections in He:
	// load the He cross sections
	CSnames.clear();
	CSnames.push_back("He MOMENTUM TRANSFER");
	eHeType.push_back(ELASTIC);
	CSnames.push_back("HE EXCITATION");
	eHeType.push_back(EXCITATION);
	CSnames.push_back("HE TOTAL IONIZATION");
	eHeType.push_back(IONIZATION);
	
	Species<D>::load_CS("electron.txt", eHeCS, eHeLoss, CSnames,"HELIUM");
	
	// compute the collision frequency with Ar
	eHesvmax = svmax_find(eHeCS,Species<D>::v_max*5.0);
	cout << "eHesvmax = "<< eHesvmax<<endl;

	hit1=hit2=0;
    };
    void lifetime_init()
    {
	// finally obtain the frequency
        if(Species<D>::species_list[O2] != NULL)
            eO2Freq = eO2svmax * Species<D>::species_list[O2]->density;
        else
            eO2Freq = eO2svmax * Species<D>::p_param->density[O2];

        if(Species<D>::species_list[ARGON] != NULL)
            eArFreq = eArsvmax * Species<D>::species_list[ARGON]->density;
        else
            eArFreq = eArsvmax * Species<D>::p_param->density[ARGON];

        if(Species<D>::species_list[HELIUM] != NULL)
            eHeFreq = eHesvmax * Species<D>::species_list[HELIUM]->density;
        else
            eHeFreq = eHesvmax * Species<D>::p_param->density[HELIUM];
	Species<D>::lifetime = 1.0/(eO2Freq + eArFreq + eHeFreq);
	cout << "electron lifetime " << Species<D>::lifetime/Species<D>::dt <<endl;
	Species<D>::lifetime_init();
    }
    int hit1, hit2;
    void scatter(t_particle &particle, int species, double svmax, vector<vec_interpolate*> & CCS, vector<coll_type> & CType, vector<double> & CLoss);
};
template <int D>
void t_elon<D>::scatter(t_particle &particle)
{
    // first we have to choose the interacting species
    // the e-e collisions are neglected for now
    double freq = eArFreq + eO2Freq + eHeFreq;
    double gamma = Species<D>::rnd->uni() * freq;
    if( (gamma -= eArFreq) < 0 )
    {
        scatter(particle, ARGON, eArsvmax, eArCS, eArType, eArLoss);
    }
    else if( (gamma -= eHeFreq) < 0 )
    {
        scatter(particle, HELIUM, eHesvmax, eHeCS, eHeType, eHeLoss);
    }
    else
    {
        scatter(particle, O2, eO2svmax, eO2CS, eO2Type, eO2Loss);
    }
    return;

    //rozhodneme, jestli dochazi k e-e nebo e-i
    {
	// e-e INTERAKCE

	//Generujeme partnera
	int i = rand()%Species<D>::particles.size();
	while(Species<D>::particles[i].empty == true)
	{
	    i = rand()%Species<D>::particles.size();
	}

	t_particle *p_c2 = &(Species<D>::particles[i]);

	//vypoceteme relativni rychlost elektronu
	double v_rel = (SQR(p_c2->vx-particle.vx) + SQR(p_c2->vz-particle.vz) + SQR(p_c2->vy-particle.vy));
	double E = 0.5*Species<D>::mass*v_rel;
	v_rel = sqrt(v_rel);

	//spoctemez nyni srazkovou freq

	double Lambda = elon::Lamb/(2*E);
	double vsigma = M_PI*SQR(Lambda)*log(elon::l_debye/Lambda)*v_rel;

	// double vsigma = sigma_ee_el_tab(E)*v_rel;

	//normujeme na maximalni srazkovou freq
	vsigma /= 6.371e-10;

	//rozhodneme, jestli dochazi ke srazce
	if(vsigma > Species<D>::rnd->uni())
	{
	    //prepocet v_1 do tezistove soustavy
	    double v1_cm_x = (particle.vx - p_c2->vx)*0.5;
	    double v1_cm_y = (particle.vz - p_c2->vz)*0.5;
	    double v1_cm_z = (particle.vy - p_c2->vy)*0.5;
	    double v1_cm = sqrt(SQR(v1_cm_x) + SQR(v1_cm_y) + SQR(v1_cm_z));

	    //provedeni nahodne rotace
	    Species<D>::rnd->rot(v1_cm,v1_cm_x,v1_cm_y,v1_cm_z);

	    //zpetna transformace
	    //p_c->vr = v1_cm_x + v_cm_x;
	    particle.vx = v1_cm_x + (particle.vx + p_c2->vx)*0.5;
	    particle.vz = v1_cm_y + (particle.vz + p_c2->vz)*0.5;
	    particle.vy = v1_cm_z + (particle.vy + p_c2->vy)*0.5;

	    //rychlost druhe castice je v tezistove soust. opacna:
	    p_c2->vx = -2*v1_cm_x + particle.vx;
	    p_c2->vz = -2*v1_cm_y + particle.vz;
	    p_c2->vy = -2*v1_cm_z + particle.vy;

	}
    }
}

template <int D>
void t_elon<D>::scatter(t_particle &particle, int species, double svmax, vector<vec_interpolate*> & CCS, vector<coll_type> & CType, vector<double> & CLoss)
{
    // compute the electron energy
    double v = norm(particle.vx, particle.vz, particle.vy);
    double const_E = -0.5*Species<D>::mass/Species<D>::charge;
    double E = const_E*v*v;

    // choose the collisional process randomly
    double gamma = Species<D>::rnd->uni() * svmax*1e20;
    double tmp=0;
    unsigned int i;
    for(i=0; i<CCS.size(); i++)
    {
        tmp += (*CCS[i])(E)*v;
        if(tmp > gamma)
            break;
    }

    // select apropriate collision type
    if(i==CCS.size()) return;// NULL collision
    switch(CType[i])
    {
        case IONIZATION :
            E -= CLoss[i];
            E *= 1-Species<D>::rnd->uni(); //TODO ???
            hit1++;
            break;

        case EXCITATION :
            E -= CLoss[i];
            hit2++;
            break;

        case ELASTIC :
            // TODO extend this to low energy electrons
            E *= 1-Species<D>::rnd->uni()*2*Species<D>::mass/Species<D>::species_list[species]->mass;
            break;
    }
    // do the random rotation if any non-null collision occured
    if(i<CCS.size())
        Species<D>::rnd->rot(Species<D>::veV(E),particle.vx,particle.vz,particle.vy);
	
}
#endif
