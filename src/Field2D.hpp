#ifndef FIELD2D_H
#define FIELD2D_H

#include "Array.hpp"
#include "util.cpp"
#include "mymath.cpp"
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>


class Field2D : public Array2D<double>
{
    public:
        Field2D(int x_sampl, int z_sampl, double _dx, double _dy, double _xmin=0, double _ymin=0) :
            Array2D<double>(x_sampl, z_sampl), dx(_dx), dy(_dy),
            idx(1./dx), idy(1./dy), xmin(_xmin), ymin(_ymin) {};
        Field2D() : Array2D<double>(), idx(0), idy(0), xmin(0), ymin(0) {};
        void resize(int x_sampl, int z_sampl, double dx, double dz, double _xmin=0, double _ymin=0);
	inline void accumulate(double charge, double x, double y);
        inline double interpolate(double r, double z) const;
        inline void grad(double x, double y, double &grad_x, double &grad_y) const;
        inline double x(int i) { return i*dx+xmin;};
        inline double y(int j) { return j*dy+ymin;};
        bool hasnan();
	void print( ostream & out = std::cout , double factor = 1.0);
	void print( const char * filename , double factor = 1.0);
	void load( const char * filename);
    private:
        double dx, dy;
	double idx, idy;
        double xmin, ymin;
};


inline void Field2D::accumulate(double charge, double x, double y)
{
    x -= xmin;
    y -= ymin;
    int i = (int)(x * idx);
    int j = (int)(y * idy);

    double u = x*idx - i;
    double v = y*idy - j;

    if(i<0 || i>jmax-1 || j<0 || j>lmax-1)
	throw std::runtime_error("Field2D::accumulate() outside of range\n");
    data[i][j] += (1-u)*(1-v)*charge;
    data[i+1][j] += u*(1-v)*charge;
    data[i][j+1] += (1-u)*v*charge;
    data[i+1][j+1] += u*v*charge;

}

inline double Field2D::interpolate(double x, double y) const
{
    x -= xmin;
    y -= ymin;
    int i = (int)(x * idx);
    int j = (int)(y * idy);

    double u = x*idx - i;
    double v = y*idy - j;

    if(i<0 || i>jmax-1 || j<0 || j>lmax-1)
	throw std::runtime_error("Field2D::interpolate() outside of range\n");

    return (1-u)*(1-v)*data[i][j] + u*(1-v)*data[i+1][j] + (1-u)*v*data[i][j+1] + u*v*data[i+1][j+1];
}

inline void Field2D::grad(double x, double y, double &grad_x, double &grad_y) const
{
    int i,j;
    double g1,g2,g3,g4;
    double fx,fy;

    x -= xmin;
    y -= ymin;
    /*
     * nejprve x-ova slozka
     */
    i = (int)(x*idx+0.5);
    j = (int)(y*idy);
    j = min(j,lmax-2);
    if(i>0 && i<jmax-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (data[i][j]-data[i-1][j])*idx;
	g2 = (data[i][j+1]-data[i-1][j+1])*idx;
	g3 = (data[i+1][j+1]-data[i][j+1])*idx;
	g4 = (data[i+1][j]-data[i][j])*idx;

	//interpolace
	fx = x*idx-i+.5;
	fy = y*idy-j;
	grad_x = g1*(1-fx)*(1-fy) + g2*(1-fx)*fy + g3*fx*fy + g4*fx*(1-fy);
    }else if(i==jmax-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (data[i][j]-data[i-1][j])*idx;
	g2 = (data[i][j+1]-data[i-1][j+1])*idx;

	//interpolace
	fy = y*idy-j;
	grad_x = g1*(1-fy) + g2*fy;
    }else if(i==0)
    {
	//vypocet gradientu v mrizovych bodech
	g3 = (data[i+1][j+1]-data[i][j+1])*idx;
	g4 = (data[i+1][j]-data[i][j])*idx;

	//interpolace
	fy = y*idy-j;
	grad_x = g3*fy + g4*(1-fy);
    }



    /*
     * ted y-ova slozka
     */
    i = (int)(x*idx);
    j = (int)(y*idy+0.5);
    i = min(i,jmax-2);


    if(j>0 && j<lmax-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (data[i][j]-data[i][j-1])*idy;
	g2 = (data[i+1][j]-data[i+1][j-1])*idy;
	g3 = (data[i+1][j+1]-data[i+1][j])*idy;
	g4 = (data[i][j+1]-data[i][j])*idy;

	//interpolace
	fx = x*idx-i;
	fy = y*idy-j+0.5;
	grad_y = g1*(1-fx)*(1-fy) + g2*(1-fy)*fx + g3*fx*fy + g4*fy*(1-fx);
    }else if(j==lmax-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (data[i][j]-data[i][j-1])*idy;
	g2 = (data[i+1][j]-data[i+1][j-1])*idy;

	//interpolace
	fx = x*idx-i;
	grad_y = g1*(1-fx) + g2*fx;
    }else if(j==0)
    {
	//vypocet gradientu v mrizovych bodech
	g3 = (data[i+1][j+1]-data[i+1][j])*idy;
	g4 = (data[i][j+1]-data[i][j])*idy;

	//interpolace
	fx = x*idx-i;
	grad_y = g3*fx + g4*(1-fx);
    }

}
#endif
