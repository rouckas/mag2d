#ifndef OXYGEN_H
#define OXYGEN_H
#include "particles.cpp"
#include "tabulate.cpp"
#include "random.cpp"
#include "input.cpp"
#include "mymath.cpp"

#define M_ARP 6.68173e-26       //kg


//XXX TODO
class t_O2_neutral : public Species
{
    void scatter(t_particle &particle){};
    public:
    t_O2_neutral(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, Species *_species_list[], double _dt=1e-8)
	: Species(0, 0, param, random, _field, 32*1.66e-27, _species_list, O2)
    {
	charge = 0.0;
    };
};

class t_O2_pos : public Species
{
    public:
    void scatter(t_particle &particle);
    private:
    
    enum coll_type { ELASTIC, CX };
    vector<vec_interpolate*> O2O2CS;
    vector<double> O2O2Loss;
    vector<coll_type> O2O2Type;
    double O2O2Freq;
    double O2O2svmax;
    // oxygen cross sections in Ar plasma
    vector<vec_interpolate*> O2ArCS;
    vector<double> O2ArLoss;
    vector<coll_type> O2ArType;
    double O2ArFreq;
    double O2Arsvmax;

    public:
    t_O2_pos(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random,
	   Fields *_field, Species * _species_list[], double _dt=1e-8)
	: Species(_n, n2, param, random, _field, 32*1.66e-27, _species_list, O2_POS) 
    {
	charge = 1.602189e-19;
	

	// load the O2 cross sections
	vector<string> CSnames;
	CSnames.push_back("O2 RESONANT CHARGE TRANSFER");
	O2O2Type.push_back(CX);
	//muzeme vychazet z udaju na
	// http://cccbdb.nist.gov/
	
	load_CS("oxygen.txt", O2O2CS, O2O2Loss, CSnames);

	// compute the collision frequency with O2
	O2O2svmax = svmax_find(O2O2CS,v_max*5.0);
	
	// load the Ar cross sections
	CSnames.clear();
	CSnames.push_back("ELASTIC MOMENTUM-TRANSFER CROSS");
	O2ArType.push_back(ELASTIC);
	
	load_CS("electron.txt", O2ArCS, O2ArLoss, CSnames);

	// compute the collision frequency with Ar
	O2Arsvmax = svmax_find(O2ArCS,v_max*5.0);

	hit1=hit2=0;
    };
    void lifetime_init()
    {
	// finally obtain the frequency
	O2O2Freq = O2O2svmax * species_list[O2]->density;
	O2ArFreq = O2Arsvmax * species_list[ARGON]->density;
	lifetime = 1.0/(O2O2Freq + O2ArFreq);
	cout << "oxygen2 lifetime " << lifetime/dt <<" "<<lifetime<<endl;
	Species::lifetime_init();
    }
    int hit1, hit2;
};
void t_O2_pos::scatter(t_particle &particle)
{
    /*
    cout << "O2 scatter\n";
    */
    //rnd->rot(particle.vr, particle.vz,particle.vt);
    //return;
    // first we have to choose the interacting species
    double freq = O2ArFreq + O2O2Freq;
    if( rnd->uni() < O2ArFreq/freq )
    {
	// *** collision with argon ***
	// generate some random interaction partner
	double vr2,vz2,vt2;
	species_list[ARGON]->rndv(vr2,vz2,vt2);

	double v = norm(particle.vr-vr2, particle.vz-vz2, particle.vt-vt2);
	double const_E = -0.5*mass/charge;
	double E = const_E*v*v;
	
	// choose the collisional process randomly
	double gamma = rnd->uni() * O2Arsvmax*1e20;
	double tmp=0;
	unsigned int i;
	// there is some room for improvement in the following
	// cycle. we should check for the most probable collision
	// (null collision) first
	for(i=0; i<O2ArCS.size(); i++)
	{
	    tmp += (*O2ArCS[i])(E)*v;
	    if(tmp > gamma)
		break;
	}

	// select apropriate collision type
	switch(O2ArType[i])
	{
	    case ELASTIC :
		E *= 1-rnd->uni()*2*mass/species_list[ARGON]->mass;
		break;
	    case CX :
		break;
	}
	// do the random rotation if any non-null collision occured
	
    }
    else
    {
	// collision with O2
	// *** collision with oxygen ***
	// generate some random interaction partner
	double vr2,vz2,vt2;
	species_list[O2]->rndv(vr2,vz2,vt2);

	double v = norm(particle.vr-vr2, particle.vz-vz2, particle.vt-vt2);
	double const_E = -0.5*mass/charge;
	double E = const_E*v*v;
	
	// choose the collisional process randomly
	double gamma = rnd->uni() * O2O2svmax*1e20;
	double tmp=0;
	unsigned int i;
	for(i=0; i<O2O2CS.size(); i++)
	{
	    tmp += (*O2O2CS[i])(E)*v;
	    if(tmp > gamma)
		break;
	}

	// select apropriate collision type
	switch(O2O2Type[i])
	{
	    case ELASTIC :
		cerr << "Warning: oxygen elastic interaction not implemented\n";
		break;
	    case CX :
//		cout << "CX";
		particle.vr = vr2;
		particle.vz = vz2;
		particle.vt = vt2;
		break;
	}
	// do the random rotation if any non-null collision occured
	
    }
    //return;

};




class t_O_neg : public Species
{
    public:
    void scatter(t_particle &particle);
    private:
    
    enum coll_type { ELASTIC, CX };
    vector<vec_interpolate*> OO2CS;
    vector<double> OO2Loss;
    vector<coll_type> OO2Type;
    double OO2Freq;
    double OO2svmax;
    // oxygen cross sections in Ar plasma
    vector<vec_interpolate*> OArCS;
    vector<double> OArLoss;
    vector<coll_type> OArType;
    double OArFreq;
    double OArsvmax;

    public:
    t_O_neg(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random,
	   Fields *_field, Species * _species_list[], double _dt=1e-8)
	: Species(_n, n2, param, random, _field, 16*1.66e-27, _species_list, O_NEG) 
    {
	charge = -1.602189e-19;
	

	// load the O2 cross sections
	vector<string> CSnames;
	CSnames.push_back("O- O2 ELASTIC");
	OO2Type.push_back(ELASTIC);
	//muzeme vychazet z udaju na
	// http://cccbdb.nist.gov/
	
	load_CS("oxygen_atomic_multiplied.txt", OO2CS, OO2Loss, CSnames,"",2);

	// compute the collision frequency with O2
	OO2svmax = svmax_find(OO2CS,v_max*5.0);
	
	// load the Ar cross sections
	CSnames.clear();
	CSnames.push_back("ELASTIC MOMENTUM-TRANSFER CROSS");
	OArType.push_back(ELASTIC);
	
	load_CS("electron.txt", OArCS, OArLoss, CSnames);

	// compute the collision frequency with Ar
	OArsvmax = svmax_find(OArCS,v_max*5.0);

	hit1=hit2=0;
    };
    void lifetime_init()
    {
	// finally obtain the frequency
	OO2Freq = OO2svmax * species_list[O2]->density;
	OArFreq = OArsvmax * species_list[ARGON]->density;
	lifetime = 1.0/(OO2Freq + OArFreq);
	cout << "O- lifetime " << lifetime/dt <<" "<<lifetime<<endl;
	Species::lifetime_init();
    }
    int hit1, hit2;
};
void t_O_neg::scatter(t_particle &particle)
{
    // first we have to choose the interacting species
    double freq = OArFreq + OO2Freq;
    if( rnd->uni() < OArFreq/freq )
    {
	// *** collision with argon ***
	// generate some random interaction partner
	double vr2,vz2,vt2;
	species_list[ARGON]->rndv(vr2,vz2,vt2);

	double v = norm(particle.vr-vr2, particle.vz-vz2, particle.vt-vt2);
	double const_E = -0.5*mass/charge;
	double E = const_E*v*v;
	
	// choose the collisional process randomly
	double gamma = rnd->uni() * OArsvmax*1e20;
	double tmp=0;
	unsigned int i;
	// there is some room for improvement in the following
	// cycle. we should check for the most probable collision
	// (null collision) first
	for(i=0; i<OArCS.size(); i++)
	{
	    tmp += (*OArCS[i])(E)*v;
	    if(tmp > gamma)
		break;
	}

	// select apropriate collision type
	switch(OArType[i])
	{
	    case ELASTIC :
		E *= 1-rnd->uni()*2*mass/species_list[ARGON]->mass;
		break;
	    case CX :
		break;
	}
	// do the random rotation if any non-null collision occured
	
    }
    else
    {
	// collision with O2
	// *** collision with oxygen ***
	// generate some random interaction partner
	double vr2,vz2,vt2;
	species_list[O2]->rndv(vr2,vz2,vt2);

	double v = norm(particle.vr-vr2, particle.vz-vz2, particle.vt-vt2);
	double const_E = -0.5*mass/charge;
	double E = const_E*v*v;
	
	// choose the collisional process randomly
	double gamma = rnd->uni() * OO2svmax*1e20;
	double tmp=0;
	unsigned int i;
	for(i=0; i<OO2CS.size(); i++)
	{
	    tmp += (*OO2CS[i])(E)*v;
	    if(tmp > gamma)
		break;
	}

	// select apropriate collision type
	switch(OO2Type[i])
	{

	    case ELASTIC:
		{
		    double m2 = species_list[O2]->mass;
		    double tmp = 1.0/(mass+m2);
		    //prepocet v_1 do tezistove soustavy
		    double v1_cm_x = (particle.vr - vr2)*m2*tmp;
		    double v1_cm_y = (particle.vz - vz2)*m2*tmp;
		    double v1_cm_z = (particle.vt - vt2)*m2*tmp;
		    //double v1_cm = sqrt(SQR(v1_cm_x) + SQR(v1_cm_y) + SQR(v1_cm_z));

		    //provedeni nahodne rotace
		    rnd->rot(v1_cm_x,v1_cm_y,v1_cm_z);

		    //zpetna transformace
		    //particle.vr = v1_cm_x + v_cm_x;
		    particle.vr = v1_cm_x + (particle.vr*mass + vr2*m2)*tmp;
		    particle.vz = v1_cm_y + (particle.vz*mass + vz2*m2)*tmp;
		    particle.vt = v1_cm_z + (particle.vt*mass + vt2*m2)*tmp;
		}
		break;
	    case CX :
		particle.vr = vr2;
		particle.vz = vz2;
		particle.vt = vt2;
		break;
	}
    }
}

#endif







