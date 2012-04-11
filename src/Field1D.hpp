#ifndef FIELD1D_H
#define FIELD1D_H

#include "Array.hpp"
#include "param.hpp"
#include "mymath.cpp"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

class Field1D : public Array1D<double>
{
    public:
        Field1D(int x_sampl, double dx) : 
            Array1D<double>(x_sampl), xmin(0),
            idx(1.0/dx) {};

        Field1D(Param &param) :
            Array1D<double>(param.x_sampl), xmin(0),
            idx(param.idx),
            xmax(param.x_max) {};

        void resize(int x_sampl, double dx, double xmin=0);
        inline void accumulate(double charge, double x);
        inline double interpolate(double x);
        inline double grad(double x);
        bool hasnan();
        void print( ostream & out = cout, double factor = 1.0);
        void print( const char * filename, double factor = 1.0);

    private:
        double xmin;
        double idx;
        double xmax;
};

inline void Field1D::accumulate(double charge, double x)
{
    x -= xmin;
    int i = (int)(x * idx);
    double u = x*idx - i;

    if(i<0 || i>imax-1)
	throw std::runtime_error("Field1D::accumulate() outside of range\n");

    data[i] += (1-u)*charge;
    data[i+1] += u*charge;
}

inline double Field1D::interpolate(double x)
{
    x -= xmin;
    int i = (int)(x * idx);
    double u = x*idx - i;

    if(i<0 || i>imax-1)
	throw std::runtime_error("Field1D::interpolate() outside of range\n");

    double res = 0;
    res += (1-u)*data[i];
    res += u*data[i+1];

    return res;
}

inline double Field1D::grad(double x)
{
    double g1, g2;
    int i;
    double u;

    // x component of the gradient
    x = clamp(x*idx, 0.5, xmax*idx - 0.5);
    //this should ensure, than no index is outside the
    //array, but floating point arithmetics is a bitch...

    i = (int)(x + 0.5);
    u = x + 0.5 - i;

    g1 = data[i] - data[i-1];
    g2 = data[i+1] - data[i];

    return (g1*u + g2*(1-u))*idx;
}


#endif
