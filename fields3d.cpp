#ifndef FIELD3D_H
#define FIELD3D_H

#include "matrix.cpp"
#include <iostream>
#include <fstream>
using namespace std;

class Field3D : public Array3D<double>
{
    public:
        Field3D(int x_sampl, int y_sampl, int z_sampl, double dx, double dy, double dz) : 
            Array3D<double>(x_sampl, y_sampl, z_sampl), xmin(0), ymin(0), zmin(0),
            idx(1.0/dx), idy(1.0/dy), idz(1.0/dz) {};

        inline void accumulate(double charge, double x, double y, double z);
        inline double interpolate(double x, double y, double z);
        void print( ostream & out = cout , double factor = 1.0);
        void print( const char * filename , double factor = 1.0);

    private:
        double idx, idy, idz;
        double xmin, ymin, zmin;
};

inline void Field3D::accumulate(double charge, double x, double y, double z)
{
    x -= xmin;
    y -= ymin;
    z -= zmin;
    int i = (int)(x * idx);
    int j = (int)(y * idy);
    int k = (int)(z * idz);

    double u = x*idx - i;
    double v = y*idy - j;
    double w = z*idz - k;

    if(i<0 || i>imax-1 || j<0 || j>jmax-1 || k<0 || k>kmax-1)
	throw std::runtime_error("Field3D::accumulate() outside of range\n");

    data[i][j][k] += (1-u)*(1-v)*(1-w)*charge;
    data[i+1][j][k] += u*(1-v)*(1-w)*charge;
    data[i][j+1][k] += (1-u)*v*(1-w)*charge;
    data[i+1][j+1][k] += u*v*(1-w)*charge;

    data[i][j][k+1] += (1-u)*(1-v)*w*charge;
    data[i+1][j][k+1] += u*(1-v)*w*charge;
    data[i][j+1][k+1] += (1-u)*v*w*charge;
    data[i+1][j+1][k+1] += u*v*w*charge;
}

inline double Field3D::interpolate(double x, double y, double z)
{
    x -= xmin;
    y -= ymin;
    z -= zmin;
    int i = (int)(x * idx);
    int j = (int)(y * idy);
    int k = (int)(z * idz);

    double u = x*idx - i;
    double v = y*idy - j;
    double w = z*idz - k;

    if(i<0 || i>imax-1 || j<0 || j>jmax-1 || k<0 || k>kmax-1)
	throw std::runtime_error("Field3D::interpolate() outside of range\n");

    double res = 0;

    res += (1-u)*(1-v)*(1-w) * data[i][j][k];
    res += u*(1-v)*(1-w) * data[i+1][j][k];
    res += (1-u)*v*(1-w) * data[i][j+1][k];
    res += u*v*(1-w) * data[i+1][j+1][k];

    res += (1-u)*(1-v)*w * data[i][j][k+1];
    res += u*(1-v)*w * data[i+1][j][k+1];
    res += (1-u)*v*w * data[i][j+1][k+1];
    res += u*v*w * data[i+1][j+1][k+1];

    return res;
}


void Field3D::print( ostream & out, double factor)
{
    double dx=1.0/idx, dy=1.0/idy, dz=1.0/idz;
    for(int i=0; i<imax; i++)
    {
        for(int j=0; j<jmax; j++)
        {
            for(int k=0; k<kmax; k++)
                out << i*dx <<"\t"<< j*dy <<"\t"<< k*dz <<"\t"<< data[i][j][k]*factor << endl;
            out << endl;
        }
        out << endl;
    }
}

void Field3D::print( const char * filename , double factor)
{
    ofstream out(filename);
    print(out, factor);
}

#endif
