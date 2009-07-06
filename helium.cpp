#ifndef HELIUM_H
#define HELIUM_H

#include "particles.cpp"

template <int D>
class t_helium_neutral : public Species<D>
{
    void scatter(t_particle &particle){};
    public:
    t_helium_neutral(int _n, int n2, double temperature, double vsigma_max, Param &param, t_random &random, 
	    Fields *_field, BaseSpecies *_species_list[], double _dt=1e-8)
	: Species<D>(0, 0, param, random, _field, 6.6905e-27, _species_list, HELIUM)
    {
	Species<D>::charge = 0.0;
    };
};
#endif
