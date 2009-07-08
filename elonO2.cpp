#ifndef ELONO2_H
#define ELONO2_H
#include "particles.cpp"
#include "tabulate.cpp"
#include "random.cpp"
#include "input.hpp"
#include "mymath.cpp"

#define M_ARP 6.68173e-26       //kg

namespace elon{
    static const double E_MAX=30;
    double l_debye;
    double Lamb;
    double charge;
};


class t_elon : public Species
{
    public:
    void scatter(t_particle &particle);
    private:
    
    enum coll_type { ELASTIC, EXCITATION, IONIZATION };
    // electron inelastic cross sections in O2 plasma
    // it is not effective to use virtual class here.
    // the test took 
    // 0.9 sec with vector<callable*> where
    // CS were implemented as vec_interpolate : public callable
    // 0.5 sec with vector<vec_interpolate*>
    // 0.27 sec with normal functions using equidistant table
    vector<vec_interpolate*> eO2CS;
    vector<double> eO2Loss;
    vector<coll_type> eO2Type;
    double eO2Freq;
    double eO2svmax;
    // electron inelastic cross sections in Ar plasma
    vector<vec_interpolate*> eArCS;
    vector<double> eArLoss;
    vector<coll_type> eArType;
    double eArFreq;
    double eArsvmax;

    public:
    t_elon(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random,
	   Fields *_field, Species * _species_list[], double _dt=1e-11)
	: Species(_n, n2, param, random, _field, 9.109534e-31, _species_list, ELECTRON) 
    {
	charge = -1.602189e-19;
	//double vsigma_coul_max = 6.371e-10;

	elon::Lamb = charge*charge/(4.0*M_PI*param.eps_0);
	elon::l_debye = param.eps_0*param.k_B*param.temperature[ELECTRON]/(density*charge*charge);
	elon::charge = charge;
	cerr << "e-  "<< lifetime/dt << " timesteps per collision" << endl;

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
	
	load_CS("electron.txt", eO2CS, eO2Loss, CSnames);

	// compute the collision frequency with O2
	eO2svmax = svmax_find(eO2CS,v_max*5.0);
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
	
	load_CS("electron.txt", eArCS, eArLoss, CSnames,"ARGON");
	
	// compute the collision frequency with Ar
	eArsvmax = svmax_find(eArCS,v_max*5.0);
	cout << "eArsvmax = "<< eArsvmax<<endl;

	hit1=hit2=0;
    };
    void lifetime_init()
    {
	// finally obtain the frequency
	eO2Freq = eO2svmax * species_list[O2]->density;
	eArFreq = eArsvmax * species_list[ARGON]->density;
	lifetime = 1.0/(eO2Freq + eArFreq);
	cout << "electron lifetime " << lifetime/dt <<endl;
	Species::lifetime_init();
    }
    int hit1, hit2;
};
void t_elon::scatter(t_particle &particle)
{
    // first we have to choose the interacting species
    double freq = eArFreq + eO2Freq;
    if( rnd->uni() < eArFreq/freq )
    {
	// *** collision with argon ***
	// compute the electron energy
	double v = norm(particle.vr, particle.vz, particle.vt);
	double const_E = -0.5*mass/charge;
	double E = const_E*v*v;
	
	// choose the collisional process randomly
	double gamma = rnd->uni() * eArsvmax*1e20;
	double tmp=0;
	unsigned int i;
	// there is some room for improvement in the following
	// cycle. we should check for the most probable collision
	// (null collision) first
	for(i=0; i<eArCS.size(); i++)
	{
	    tmp += (*eArCS[i])(E)*v;
	    if(tmp > gamma)
		break;
	}

	// select apropriate collision type
	if(i==eArCS.size()) return;// NULL collision
	switch(eArType[i])
	{
	    case IONIZATION :
		E -= eArLoss[i];
		E *= 1-rnd->uni(); //TODO ???
		hit1++;
		break;

	    case EXCITATION :
		E -= eArLoss[i];
		hit2++;
		break;

	    case ELASTIC :
		E *= 1-rnd->uni()*2*mass/M_ARP;
		break;
	}
	// do the random rotation if any non-null collision occured
	if(i<eArCS.size())
	    rnd->rot(veV(E),particle.vr,particle.vz,particle.vt);
	
    }
    else
    {
	// *** collision with oxygen ***
	// compute the electron energy
	double v = norm(particle.vr, particle.vz, particle.vt);
	double const_E = -0.5*mass/charge;
	double E = const_E*v*v;
	
	// choose the collisional process randomly
	double gamma = rnd->uni() * eO2svmax*1e20;
	double tmp=0;
	unsigned int i;
	// there is some room for improvement in the following
	// cycle. we should check for the most probable collision
	// (null collision) first
	for(i=0; i<eO2CS.size(); i++)
	{
	    tmp += (*eO2CS[i])(E)*v;
	    if(tmp > gamma)
		break;
	}

	// select apropriate collision type
	if(i==eO2CS.size()) return;// NULL collision
	switch(eO2Type[i])
	{
	    case IONIZATION :
		E -= eO2Loss[i];
		E *= 1-rnd->uni(); //TODO ???
		hit1++;
		break;

	    case EXCITATION :
		E -= eO2Loss[i];
		hit2++;
		break;

	    case ELASTIC :
		E *= 1-rnd->uni()*2*mass/species_list[O2]->mass;
		break;
	}
	// do the random rotation if any non-null collision occured
	if(i<eO2CS.size())
	    rnd->rot(veV(E),particle.vr,particle.vz,particle.vt);
	
    }
    //return;

    //rozhodneme, jestli dochazi k e-e nebo e-i
    if(  true || rnd->uni() > 0.058 )
    {
    }// konec e-i interakce
    else
    {
	// e-e INTERAKCE

	//Generujeme partnera
	int i = rand()%particles.size();
	while(particles[i].empty == true)
	{
	    i = rand()%particles.size();
	}

	t_particle *p_c2 = &(particles[i]);

	//vypoceteme relativni rychlost elektronu
	double v_rel = (SQR(p_c2->vr-particle.vr) + SQR(p_c2->vz-particle.vz) + SQR(p_c2->vt-particle.vt));
	double E = 0.5*mass*v_rel;
	v_rel = sqrt(v_rel);

	//spoctemez nyni srazkovou freq

	double Lambda = elon::Lamb/(2*E);
	double vsigma = M_PI*SQR(Lambda)*log(elon::l_debye/Lambda)*v_rel;

	// double vsigma = sigma_ee_el_tab(E)*v_rel;

	//normujeme na maximalni srazkovou freq
	vsigma /= 6.371e-10;

	//rozhodneme, jestli dochazi ke srazce
	if(vsigma > rnd->uni())
	{
	    //prepocet v_1 do tezistove soustavy
	    double v1_cm_x = (particle.vr - p_c2->vr)*0.5;
	    double v1_cm_y = (particle.vz - p_c2->vz)*0.5;
	    double v1_cm_z = (particle.vt - p_c2->vt)*0.5;
	    double v1_cm = sqrt(SQR(v1_cm_x) + SQR(v1_cm_y) + SQR(v1_cm_z));

	    //provedeni nahodne rotace
	    rnd->rot(v1_cm,v1_cm_x,v1_cm_y,v1_cm_z);

	    //zpetna transformace
	    //p_c->vr = v1_cm_x + v_cm_x;
	    particle.vr = v1_cm_x + (particle.vr + p_c2->vr)*0.5;
	    particle.vz = v1_cm_y + (particle.vz + p_c2->vz)*0.5;
	    particle.vt = v1_cm_z + (particle.vt + p_c2->vt)*0.5;

	    //rychlost druhe castice je v tezistove soust. opacna:
	    p_c2->vr = -2*v1_cm_x + particle.vr;
	    p_c2->vz = -2*v1_cm_y + particle.vz;
	    p_c2->vt = -2*v1_cm_z + particle.vt;

	}
    }
}

#endif
