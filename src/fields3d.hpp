#ifndef FIELDS3D_H
#define FIELDS3D_H

#include "Array.hpp"
#include "param.hpp"
#include "mymath.cpp"
#include "fields.hpp"
#include "Field3D.hpp"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;


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

#endif
