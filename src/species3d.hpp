
#ifndef SPECIES3D_H
#define SPECIES3D_H
#include "fields3d.hpp"
#include "particles.hpp"

template <>
class Species<CARTESIAN3D> : public BaseSpecies
{
    public:
        ElMag3D & fields;
        Species(Param &param, t_random &_rnd, ElMag3D & _fields,
                double _mass, BaseSpecies * _species_list[], species_type _type) :
            BaseSpecies(0, 0, param, _rnd, _mass, _species_list, _type), fields(_fields) {};

        void advance();
};
#endif
