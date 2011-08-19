#ifndef FIELD3D_H
#define FIELD3D_H

#include "matrix.cpp"
#include "param.hpp"
#include "mymath.cpp"
#include "fields.hpp"
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
/*
class Electrode
{
    protected:
        signed char id;
    public:
        Electrode(int i, double _voltage = 1.0) : id(i), voltage(_voltage) {};
        double voltage;
        virtual void set_mask(Array3D<signed char> &mask) = 0;
        void set_voltage(Array3D<double> &voltage_array, Array3D<signed char> &mask);
};

class Quadrupole : public Electrode
{
    public:
        Quadrupole(int i) : Electrode(i) {};
        void set_mask(Array3D<signed char> &mask)
        {
            mask[mask.imax/2][mask.jmax/2][mask.kmax/2] = id;
            // to be implemented :)
        };
};


class Geometry
{
    private:
        double idx, idy, idz;
    public:
        double x_sampl, y_sampl, z_sampl;
        Array3D<signed char> mask;
        vector<Electrode*> electrodes;
        vector<Array3D<double>*> potentials;
        Field3D voltage;
        bool is_free(double x, double y, double z)
        {
            int i = (int)(x*idx);
            int j = (int)(y*idy);
            int k = (int)(z*idz);
            if(mask[i][j][k]==FREE || mask[i+1][j][k]==FREE || mask[i][j+1][k]==FREE || mask[i+1][j+1][k]==FREE ||
                mask[i][j][k+1]==FREE || mask[i+1][j][k+1]==FREE || mask[i][j+1][k+1]==FREE || mask[i+1][j+1][k+1]==FREE)
                return true;
            else return false;
        }
        Geometry(Param &param);
};

class Solver
{
    private:
        int *Ap;
        int *Ai;
        double *Ax;
        void *Symbolic, *Numeric ;
        Param &param;
        Geometry &geometry;
        double x_sampl, y_sampl, z_sampl;

    public:
        void matrix_init();
        void umfpack_init();
        void solve(Array3D<double> &u, const Array3D<double> &voltage, Array3D<double> &rho);
        void solve_laplace(Array3D<double> &u, const Array3D<double> &voltage);
            // solver for Laplace equation (zero RHS), it assumes that voltage is set to zero in free space
        void save(string filename);
        void load(string filename);
        Solver(Geometry & _geometry, Param & _param, string filename = "");
};

class ElMag3D
{
    public:
        Field3D u, rho, voltage;
        Field3D Bx, By, Bz;
        Geometry geometry;

        ElMag3D(Param &param);
        void E(double x, double y, double z, double & Ex, double & Ey, double & Ez)
        {
            u.grad(x, y, z, Ex, Ey, Ez);
            Ex *= -1.0;
            Ey *= -1.0;
            Ez *= -1.0;
        }
        void B(double x, double y, double z, double & Bx, double & By, double & Bz)
        {
            Bx = this->Bx.interpolate(x, y, z);
            By = this->By.interpolate(x, y, z);
            Bz = this->Bz.interpolate(x, y, z);
        }
        void load_magnetic_field(const char * fname)
        {
            Bx.load_table(fname, 0);
            By.load_table(fname, 1);
            Bz.load_table(fname, 2);
        }
        void solve();
        void potential_sum();
            //TODO do a weighted sum of solutions for single electrodes,
            //for now it just performs sum
    private:
        vector<Field3D*> potentials;

        bool multielectrode;
    public:
        // solver must be initialized after geometry
        Solver solver;

};
*/

inline void Field1D::accumulate(double charge, double x)
{
    x -= xmin;
    int i = (int)(x * idx);
    double u = x*idx - i;

    if(i<0 || i>imax-1)
	throw std::runtime_error("Field3D::accumulate() outside of range\n");

    data[i] += (1-u)*charge;
    data[i+1] += u*charge;
}

inline double Field1D::interpolate(double x)
{
    x -= xmin;
    int i = (int)(x * idx);
    double u = x*idx - i;

    if(i<0 || i>imax-1)
	throw std::runtime_error("Field3D::interpolate() outside of range\n");

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
