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

class Field3D : public Array3D<double>
{
    public:
        Field3D(int x_sampl, int y_sampl, int z_sampl, double dx, double dy, double dz) : 
            Array3D<double>(x_sampl, y_sampl, z_sampl), xmin(0), ymin(0), zmin(0),
            idx(1.0/dx), idy(1.0/dy), idz(1.0/dz) {};

        Field3D(Param &param) :
            Array3D<double>(param.x_sampl, param.y_sampl, param.z_sampl), xmin(0), ymin(0), zmin(0),
            idx(param.idx), idy(param.idy), idz(param.idz),
            xmax(param.x_max), ymax(param.y_max), zmax(param.z_max) {};

        void resize(int x_sampl, int y_sampl, int z_sampl, double dx, double dy, double dz,
                double xmin=0, double ymin=0, double zmin=0);
        inline void accumulate(double charge, double x, double y, double z);
        inline double interpolate(double x, double y, double z);
        inline void grad(double x, double y, double z, double & gx, double & gy, double & gz);
        bool hasnan();
        void print( ostream & out = cout, double factor = 1.0);
        void print_vtk( ostream & out = cout);
        void print( const char * filename, string format = "table", double factor = 1.0);
        void load_table( const char * filename, int col = 0);

    private:
        inline double local_interpolate(double g[], double u, double v, double w);
        inline double grad_component(double x, double y, double z, int dirx, int dirz, int diry);
        inline double grad_component_NGP(double x, double y, double z, int dirx, int dirz, int diry);
        double xmin, ymin, zmin;
        double idx, idy, idz;
        double xmax, ymax, zmax;
};

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

inline double Field3D::local_interpolate(double g[], double u, double v, double w)
{
    return (1-u)*(1-v)*(1-w) * g[0] +
        u*(1-v)*(1-w) * g[1] +
        (1-u)*v*(1-w) * g[2] +
        u*v*(1-w) * g[3] +
        (1-u)*(1-v)*w * g[4] +
        u*(1-v)*w * g[5] +
        (1-u)*v*w * g[6] +
        u*v*w * g[7];
}

inline double Field3D::grad_component(double x, double y, double z, int dirx, int diry, int dirz)
{
    double g[8];
    int i, j, k;
    double u, v, w;

    // x component of the gradient
    x = clamp(x*idx, 0.5*dirx, xmax*idx - 0.5*dirx);
    y = clamp(y*idy, 0.5*diry, ymax*idy - 0.5*diry);
    z = clamp(z*idz, 0.5*dirz, zmax*idz - 0.5*dirz);
    //this should ensure, than no index is outside the
    //array, but floating point arithmetics is a bitch...

    i = (int)(x + 0.5*dirx);
    j = (int)(y + 0.5*diry);
    k = (int)(z + 0.5*dirz);
    u = x + 0.5*dirx - i;
    v = y + 0.5*diry - j;
    w = z + 0.5*dirz - k;

    g[0] = data[i][j][k] - data[i-dirx][j-diry][k-dirz];
    g[1] = data[i+1][j][k] - data[i+1-dirx][j-diry][k-dirz];
    g[2] = data[i][j+1][k] - data[i-dirx][j+1-diry][k-dirz];
    g[3] = data[i+1][j+1][k] - data[i+1-dirx][j+1-diry][k-dirz];

    g[4] = data[i][j][k+1] - data[i-dirx][j-diry][k+1-dirz];
    g[5] = data[i+1][j][k+1] - data[i+1-dirx][j-diry][k+1-dirz];
    g[6] = data[i][j+1][k+1] - data[i-dirx][j+1-diry][k+1-dirz];
    g[7] = data[i+1][j+1][k+1] - data[i+1-dirx][j+1-diry][k+1-dirz];

    return local_interpolate(g, u, v, w)*(idx*dirx + idy*diry + idz*dirz);
}

inline double Field3D::grad_component_NGP(double x, double y, double z, int dirx, int diry, int dirz)
{
    int i, j, k;

    x = clamp(x*idx, 0.5*dirx, xmax*idx - 0.5*dirx);
    y = clamp(y*idy, 0.5*diry, ymax*idy - 0.5*diry);
    z = clamp(z*idz, 0.5*dirz, zmax*idz - 0.5*dirz);

    i = (int)(x + 0.5*(1-dirx));
    j = (int)(y + 0.5*(1-diry));
    k = (int)(z + 0.5*(1-dirz));

    return (data[i+dirx][j+diry][k+dirz] - data[i][j][k])*(idx*dirx + idy*diry + idz*dirz);
}

inline void Field3D::grad(double x, double y, double z, double & gx, double & gy, double & gz)
{
    gx = grad_component(x, y, z, 1, 0, 0);
    gy = grad_component(x, y, z, 0, 1, 0);
    gz = grad_component(x, y, z, 0, 0, 1);
}
#endif
