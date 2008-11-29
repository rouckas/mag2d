#ifndef ELON_H
#define ELON_H
#include "particles.cpp"
#include "tabulate.cpp"
#include "random.cpp"

#define M_ARP 6.68173e-26       //kg

namespace elon{
    static const double E_MAX=30;
    double l_debye;
    double Lamb;
    double charge;
};


class t_elon : public Species
{
    void scatter(t_particle &particle);
    public:
    double eArsvmax;
    t_elon(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random,
	   Fields *_field, Species * _species_list[], double _dt=1e-11)
	: Species(_n, n2, param, random, _field, 9.109534e-31, _species_list, ELECTRON) 
    {
	charge = -1.602189e-19;
	eArsvmax = 3.214e-13;
	//double vsigma_coul_max = 6.371e-10;
	double vsigma_coul_max = 0;
	double N_coul = density*vsigma_coul_max;
	lifetime = 1.0/(vsigma_max*param.neutral_density + N_coul);
	/*
	if(dt==0)
	    dt = lifetime/10.0;
	    */

	elon::Lamb = charge*charge/(4.0*M_PI*param.eps_0);
	elon::l_debye = param.eps_0*param.k_B*temperature/(density*charge*charge);
	elon::charge = charge;
	cerr << "e-  "<< lifetime/dt << " timesteps per collision" << endl;
    };
    void lifetime_init()
    {
	// finally obtain the frequency
	double eArFreq = eArsvmax * species_list[ARGON]->density;
	lifetime = 1.0/eArFreq;
	cout << "electron lifetime " << lifetime/dt <<endl;
	Species::lifetime_init();
    }

};

double sigma_ei_el(double E)
{
    return (fabs(6/pow((1+(E/0.1)+pow(E/0.6,2.0)),3.3) - 
		1.1*pow(E,1.4)/(1.0+pow(E/15.0,1.2))/
		sqrt(1.0+pow(E/5.5,2.5)
		    +pow(E/60.0,4.1))) + 
	    0.05/SQR(1.0+E/10.0) + 0.01*E*E*E/(1.0+pow(E/12.0,6.0)))*1.0E-20;
}
tabulate sigma_ei_el_tab(0, elon::E_MAX, 1000, sigma_ei_el);

double sigma_ei_ex(double E)
{
    if(E<11.55)
	return 0;
    else 
	return 0.85*(4E-22*pow(E-11.5,1.1)*(1.0+pow(E/15.0,2.8))/(1.0+pow(E/23.0,5.5)) 
		+ 2.7E-22*(E-11.5)/pow(1.0+(E/80),1.9));
}
tabulate sigma_ei_ex_tab(0, elon::E_MAX, 1000, sigma_ei_ex);

double sigma_ei_io(double E)
{
    if(E<15.76)
	return 0;
    else 
	return 9.7E-18*(E-15.8)/SQR(70.0+E) + 6E-22*SQR(E-15.8)*exp(-E/9);
}
tabulate sigma_ei_io_tab(0, elon::E_MAX, 1000, sigma_ei_io);

double sigma_ee_el(double E)
{
    double Lambda = elon::Lamb/(2*E);
    return M_PI*SQR(Lambda)*log(elon::l_debye/Lambda);
}
tabulate sigma_ee_el_tab(0, elon::E_MAX*4*(-elon::charge), 1000, sigma_ee_el);
tabulate_periodic sin_tab(0,2*M_PI,1000,sin);
tabulate_periodic cos_tab(0,2*M_PI,1000,cos);


void t_elon::scatter(t_particle &particle)
{
    //rozhodneme, jestli dochazi k e-e nebo e-i
    if(  true || rnd->uni() > 0.058 )
    {
	double sq_v = SQR(particle.vr) + SQR(particle.vz) + SQR(particle.vt);
	double const_E = -0.5*mass/charge;
	double E = const_E*sq_v;
	double v = sqrt(sq_v);
	double vsigma[3];

	//pruzna srazka
	vsigma[0] = sigma_ei_el_tab(E)*v;

	//excitace
	vsigma[1] = sigma_ei_ex_tab(E)*v;

	//ionizace
	vsigma[2] = sigma_ei_io_tab(E)*v;


	double svmax = 3.32e-13;
	int i;
	for(i=0;i<3;i++)
	    vsigma[i] /= svmax;

	double gamma = rnd->uni();
	double tmp=0;
	for(i=0;i<3;i++)
	{
	    tmp += vsigma[i];
	    if(tmp > gamma)
		break;
	}

	switch(i)
	{

	    //pruzny rozptyl
	    case 0:
		sq_v *= 1-rnd->uni()*2*mass/M_ARP;
	//	cerr << (SQR(particle.vr)+SQR(particle.vz)+SQR(particle.vt)) / sq_v << endl;
		break;

		//excitace Ar
	    case 1:
		if(E>11.55)
		    sq_v = (E-11.55)/const_E;
		break;

		//ionizace Ar
	    case 2:
		if(E>15.76)
		{
		    //t_castice *p_castice2;
		    sq_v = (E-15.76)/const_E;
		    sq_v *= 1-rnd->uni(); //XXX ???

		    /*
		    //vytvoreni elektronu
		    p_castice2 = &castice[pop(&prazdne_castice)];
		    p_castice2->typ = ELEKTRON;

		    //vytvoreni iontu
		    p_castice2 = &castice[pop(&prazdne_castice)];
		    p_castice2->typ = ARGON;
		    */

		}
		break;

	    default:
		break;
	}

	if (i<3)
	{
	    v = sqrt(sq_v);
	    //nahodna zmena smeru
	    rnd->rot(v,particle.vr,particle.vz,particle.vt);
	}
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
