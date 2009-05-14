#ifndef FIELD_H
#define FIELD_H
//#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)>(y)?(y):(x))


#include <suitesparse/umfpack.h>
#include "param.cpp"
#include "matrix.cpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
using namespace std;

//#define SQR(x) ((x)*(x))
inline double SQR(double x){return x*x;};

enum {FIXED, FREE, BOUNDARY};


typedef struct{
            double l, r, t, b;
}t_umf_metrika;



class t_grid
//definuje geometrii experimentu
{
    public:
	int M;
	int N;
	double dr,dz;
	Matrix<char> mask;
	Matrix<t_umf_metrika> metrika;
	Matrix<double> voltage;
	t_grid(Param &param);
	Param *p_param;
        void penning_trap();
        void penning_trap_simple(double trap_voltage = -1.0);
        double U_trap;
        void rf_trap();
        void rf_22PT();
        void MAC_filter();
        void empty();
	void square_electrode(double rmin, double rmax, double zmin, double zmax, double voltage);
	void circle_electrode(double xcenter, double ycenter, double radius, double voltage);
	bool is_free(double r, double z);
};
class Field : public Matrix<double>
{
    public:
	Field(Param &param) : Matrix<double>(param.r_sampl, param.z_sampl),
            idr(param.idr), idz(param.idz), rmin(0), zmin(0) {};
        Field() : Matrix<double>(), idr(0), idz(0), rmin(0), zmin(0) {};
        void resize(double r_sampl, double z_sampl, double dr, double dz, double _rmin=0, double _zmin=0)
        {
            idr = 1.0/dr, idz = 1.0/dz;
            rmin = _rmin, zmin = _zmin;
            Matrix<double>::resize(r_sampl, z_sampl);
        }
	inline void accumulate(double charge, double x, double y);
        inline double interpolate(double r, double z);
	void print( ostream & out = cout , double factor = 1.0)
	{
	    double dr=1.0/idr, dz=1.0/idz;
	    for(int j=0; j<jmax; j++)
	    {
		for(int l=0; l<lmax; l++)
		    out << j*dr <<"\t"<< l*dz <<"\t"<< data[j][l]*factor <<endl;
		out << endl;
	    }
	}
	void print( const char * filename , double factor = 1.0)
	{
	    ofstream out(filename);
	    print(out, factor);
	}
    private:
	double idr, idz;
        double rmin, zmin;
};

class Fields
//stara se o reseni Poissonovy rce na danem gridu
{
    public:
	t_grid grid;
	Field u;
	Field uAvg;
	Field rho;
	Fields(Param &param);
	void solve();
	void boundary_solve();
	void reset();
	void E(double x, double y, double &grad_x, double &grad_y, double time = 0) ;
	void B(double x, double y, double &Br, double &Bz, double &Bt);
        void load_magnetic_field(const char * fname);
	inline void accumulate(double charge, double x, double y);
	~Fields();
    private:
        Field Br, Bz;
	double idr, idz;
	Param *p_param;
	int *Ap;
	int *Ai;
	double *Ax;
	void *Symbolic, *Numeric ;
    public:
	void u_sample();
	void u_reset();
	void u_print(const char * fname);
    private:
	int nsampl;
};
void Fields::u_sample()
{
    uAvg.add(u);
    nsampl++;
}
void Fields::u_reset()
{
    uAvg.reset();
    nsampl = 0;
}
void Fields::u_print(const char * fname)
{
    uAvg.print(fname,1.0/nsampl);
}



t_grid::t_grid(Param &param) :  M(param.r_sampl), N(param.z_sampl), 
    dr(param.dr), dz(param.dz),
    mask(param.r_sampl,param.z_sampl), metrika(param.r_sampl,param.z_sampl), voltage(param.r_sampl,param.z_sampl)
{
    p_param = &param;
			
    //penning_trap_simple(0);
    //rf_22PT();
    MAC_filter();

}
void t_grid::empty()
{
    int i, j;
    /*
     * Vytvoreni sondy
     */
    for(i=0; i<M; i++)
	for(j=0; j<N; j++)
	{
	    if(i==0 || i==M-1 || j==0 || j==N-1)
	    {
		mask[i][j] = FIXED;
		voltage[i][j] = 0.0;
	    }else
	    {
		mask[i][j] = FREE;
	    }
	}
}
void t_grid::rf_22PT()
{
    int i, j;
    /*
     * Vytvoreni sondy
     */
    for(i=0; i<M; i++)
	for(j=0; j<N; j++)
	{
	    if(i==0 || i==M-1 || j==0 || j==N-1)
	    {
		mask[i][j] = FIXED;
		voltage[i][j] = 0.0;
	    }else
	    {
		mask[i][j] = FREE;
	    }
	}

    double xcenter = 1e-2;
    double ycenter = 1e-2;
    double r_22pt = 0.75e-2;
    double r_rod = 0.05e-2;
    int npoles = 22;
    for(int i=0; i<npoles; i++)
    {
        double x = xcenter + sin(2*M_PI*(i+1.0/32)/npoles)*r_22pt;
        double y = ycenter + cos(2*M_PI*(i+1.0/32)/npoles)*r_22pt;
        int sign = i%2==0 ? -1 : 1;
        circle_electrode(x, y, r_rod, p_param->u_probe*sign);
    }

    for(i=2;i<M-2;i++)
	for(j=2;j<N-2;j++)
	{
	    if( ( mask[i-1][j] == FIXED ||
			mask[i+1][j] == FIXED ||
			mask[i][j-1] == FIXED ||
			mask[i][j+1] == FIXED ) &&
		    mask[i][j] != FIXED )
            {
		mask[i][j] = BOUNDARY;
            }
	}
}
void t_grid::rf_trap()
{
    int i, j;
    /*
     * Vytvoreni sondy
     */
    for(i=0; i<M; i++)
	for(j=0; j<N; j++)
	{
	    if(i==0 || i==M-1 || j==0 || j==N-1)
	    {
		mask[i][j] = FIXED;
		voltage[i][j] = 0.0;
	    }else
	    {
		mask[i][j] = FREE;
	    }
	}

    circle_electrode(5e-3, 1e-2, 2e-3, p_param->u_probe);
    circle_electrode(15e-3, 1e-2, 2e-3, p_param->u_probe);
    circle_electrode(1e-2, 5e-3, 2e-3, -p_param->u_probe);
    circle_electrode(1e-2, 15e-3, 2e-3, -p_param->u_probe);

    for(i=2;i<M-2;i++)
	for(j=2;j<N-2;j++)
	{
	    if( ( mask[i-1][j] == FIXED ||
			mask[i+1][j] == FIXED ||
			mask[i][j-1] == FIXED ||
			mask[i][j+1] == FIXED ) &&
		    mask[i][j] != FIXED )
            {
		mask[i][j] = BOUNDARY;
            }
	}
}
void t_grid::MAC_filter()
{
    int i, j;
    /*
     * Vytvoreni sondy
     */
    for(i=0; i<M; i++)
	for(j=0; j<N; j++)
	{
	    if(i==M-1 || j==0 || j==N-1)
	    {
		mask[i][j] = FIXED;
		voltage[i][j] = 0.0;
	    }else
	    {
		mask[i][j] = FREE;
	    }
	}
    square_electrode(5e-3, 4.5e-2, 2.5e-2, 3e-2, -.05);

    square_electrode(3e-2, 3.3e-2, 11e-2, 14e-2, -.4);
    square_electrode(4.5e-2, 4.8e-2, 15e-2, 25e-2, -.5);
    square_electrode(3e-2, 3.3e-2, 26e-2, 29e-2, -.4);

    square_electrode(5e-3, 4.5e-2, 35e-2, 35.3e-2, .0);
    square_electrode(0.0, 4.5e-2, 39.5e-2, 40e-2, 3e3);

    // collector
    for(i=2;i<M-2;i++)
	for(j=2;j<N-2;j++)
	{
	    if( ( mask[i-1][j] == FIXED ||
			mask[i+1][j] == FIXED ||
			mask[i][j-1] == FIXED ||
			mask[i][j+1] == FIXED ) &&
		    mask[i][j] != FIXED )
		mask[i][j] = BOUNDARY;
	}
}
void t_grid::penning_trap()
{
    int i, j;
    /*
     * Vytvoreni sondy
     */
    for(i=0; i<M; i++)
	for(j=0; j<N; j++)
	{
	    if(i==M-1 || j==0 || j==N-1)
	    {
		mask[i][j] = FIXED;
		voltage[i][j] = -i*p_param->dr*p_param->extern_field + p_param->r_max*p_param->extern_field*0.5;
		voltage[i][j] = 0.0;
	    }else
	    {
		mask[i][j] = FREE;
	    }
	}
    square_electrode(1.57e-2/2, 1.67e-2/2, 1e-3, 25e-3, -5);

    // collector
    square_electrode(0, 1.46e-2/2, 15e-3, 16e-3, 10);
    square_electrode(0, 7e-3/2, 12e-3, 19e-3, 10);
    square_electrode(0, 1.9e-3, 1e-3, 12e-3, 10);


    // lenses
    square_electrode(4e-3, 7e-3, 52e-3, 53e-3, -5);
    square_electrode(2.5e-3, 7e-3, 46e-3, 47e-3, 0);

    for(i=2;i<M-2;i++)
	for(j=2;j<N-2;j++)
	{
	    if( ( mask[i-1][j] == FIXED ||
			mask[i+1][j] == FIXED ||
			mask[i][j-1] == FIXED ||
			mask[i][j+1] == FIXED ) &&
		    mask[i][j] != FIXED )
		mask[i][j] = BOUNDARY;
	}
}
void t_grid::penning_trap_simple(double trap_voltage)
{
    U_trap = trap_voltage;
    int i, j;
    for(i=0; i<M; i++)
	for(j=0; j<N; j++)
	{
	    if(i==M-1 || j==0 || j==N-1)
	    {
		mask[i][j] = FIXED;
		voltage[i][j] = 0.0;
	    }else
	    {
		mask[i][j] = FREE;
	    }
	}
    // this trap consists of series of tubes on different potential
    //
    double inner_radius = 1e-2;
    double outer_radius = 1.1e-2;
    // "injection chamber" at zero potential
    square_electrode(inner_radius, outer_radius, 0, 1e-2, -0.5);

    // first closing electrode
    square_electrode(inner_radius, outer_radius, 1.1e-2, 2e-2, 0);

    // trap tube
    square_electrode(inner_radius, outer_radius, 2.1e-2, 6e-2, trap_voltage);

    // second closing electrode
    square_electrode(inner_radius, outer_radius, 6.1e-2, 7.5e-2, -10);


    for(i=2;i<M-2;i++)
	for(j=2;j<N-2;j++)
	{
	    if( ( mask[i-1][j] == FIXED ||
			mask[i+1][j] == FIXED ||
			mask[i][j-1] == FIXED ||
			mask[i][j+1] == FIXED ) &&
		    mask[i][j] != FIXED )
		mask[i][j] = BOUNDARY;
	}
}
void t_grid::square_electrode(double rmin, double rmax, double zmin, double zmax, double _voltage)
{
    for(int i=0; i<M; i++)
	for(int j=0; j<N; j++)
	{
	    double r, z;
	    r = i*dr;
	    z = j*dz;
	    // suboptimal, but simple
	    if(r>rmin && r<rmax && z>zmin && z<zmax) 
	    {
		mask[i][j] = FIXED;
		voltage[i][j] = _voltage;
	    }
	}
}
void t_grid::circle_electrode(double rcenter, double zcenter, double radius, double _voltage)
{
    double sqradius = SQR(radius);
    for(int i=0; i<M; i++)
	for(int j=0; j<N; j++)
	{
	    double r, z;
	    r = i*dr - rcenter;
	    z = j*dz - zcenter;
	    // suboptimal, but simple
	    if(SQR(r) + SQR(z) <= sqradius) 
	    {
		mask[i][j] = FIXED;
		voltage[i][j] = _voltage;
	    }
	}
}
bool t_grid::is_free(double r, double z)
{
    int i = (int)(r*p_param->idr);
    int j = (int)(z*p_param->idz);
    if(mask[i][j]==FREE || mask[i+1][j]==FREE || mask[i][j+1]==FREE || mask[i+1][j+1]==FREE)
	return true;
    else return false;
}

inline void Field::accumulate(double charge, double r, double z)
{
    int i = (int)(r * idr);
    int j = (int)(z * idz);

    double u = r*idr - i;
    double v = z*idz - j;

    if(i<0 || i>jmax-1 || j<0 || j>lmax-1)
	throw std::runtime_error("Field::accumulate() outside of range\n");
    data[i][j] += (1-u)*(1-v)*charge;
    data[i+1][j] += u*(1-v)*charge;
    data[i][j+1] += (1-u)*v*charge;
    data[i+1][j+1] += u*v*charge;

}
inline double Field::interpolate(double r, double z)
{
    int i = (int)(r * idr);
    int j = (int)(z * idz);

    double u = r*idr - i;
    double v = z*idz - j;

    if(i<0 || i>jmax-1 || j<0 || j>lmax-1)
	throw std::runtime_error("Field::interpolate() outside of range\n");

    return (1-u)*(1-v)*data[i][j] + u*(1-v)*data[i+1][j] + (1-u)*v*data[i][j+1] + u*v*data[i+1][j+1];
}
inline void Fields::accumulate(double charge, double r, double z)
{
    int i = (int)(r * idr);
    int j = (int)(z * idz);

    double u = r*idr - i;
    double v = z*idz - j;

    //if(i<0 || i>p_param->r_sampl-1 || j<0 || j>p_param->z_sampl-1)
//	throw std::runtime_error("Field::accumulate() outside of range\n");
    rho[i][j] += (1-u)*(1-v)*charge;
    rho[i+1][j] += u*(1-v)*charge;
    rho[i][j+1] += (1-u)*v*charge;
    rho[i+1][j+1] += u*v*charge;

}
Fields::Fields(Param &param) : grid(param), u(param), uAvg(param), rho(param), nsampl(0)
{
    //if(param.selfconsistent == false)
//	return;
    p_param = &param;
    int i, j, k, l;
    int n;
    idr=p_param->idr;
    idz=p_param->idz;
    double dr = p_param->dr;
    double dz = p_param->dz;

    n = param.r_sampl * param.z_sampl;

    // r_i = i*dr
    Ap = new int [n+1];
    Ai = new int [n*5]; //horni odhad pro petibodove diff schema
    Ax = new double [n*5];

    if( param.coord == CYLINDRICAL )
    {

        l = 0;
        Ap[0] = 0;
        for(i=0; i<param.r_sampl; i++)	//cislo radku - radialni souradnice
            for(j=0; j<param.z_sampl; j++)	//cislo sloupce
            {
                k = j + param.z_sampl*i;
                if(grid.mask[i][j]==FIXED)
                {
                    Ax[l] = 1;
                    Ai[l] = k;
                    Ap[k+1] = Ap[k]+1;
                    l++;
                }
                else if(i==0 && j>0 && j<param.z_sampl-1)
                {
                    double k2, k3;
                    // psi[j,k-1]
                    k2 = 1.0/(dz*dz);
                    Ax[l] = k2;
                    Ai[l] = k-1;

                    // psi[j,k]
                    k3 = 1.0/(dr*dr*0.25);
                    Ax[l+1] = -2.0*k2 - k3;
                    Ai[l+1] = k;

                    // psi[j,k+1]
                    Ax[l+2] = k2;
                    Ai[l+2] = k+1;

                    // psi[j+1,k]
                    Ax[l+3] = k3;
                    Ai[l+3] = k+param.z_sampl;

                    Ap[k+1] = Ap[k]+4;
                    l += 4;
                }
                else if(grid.mask[i][j]!=FIXED )
                {
                    double k1, k2, k3;
                    // psi[j-1,k]
                    k1 = (i-0.5)/(dr*dr*i);
                    Ax[l] = k1;
                    Ai[l] = k-param.z_sampl;

                    // psi[j,k-1]
                    k2 = 1.0/(dz*dz);
                    Ax[l+1] = k2;
                    Ai[l+1] = k-1;

                    // psi[j,k]
                    k3 = (i+0.5)/(dr*dr*i);
                    Ax[l+2] = -2.0*k2 - k1 - k3;
                    Ai[l+2] = k;

                    // psi[j,k+1]
                    Ax[l+3] = k2;
                    Ai[l+3] = k+1;

                    // psi[j+1,k]
                    Ax[l+4] = k3;
                    Ai[l+4] = k+param.z_sampl;

                    Ap[k+1] = Ap[k]+5;
                    l += 5;
                }
            }
    }
    else
    {
        l = 0;
        Ap[0] = 0;
        for(i=0; i<param.r_sampl; i++)		//cislo radku
            for(j=0; j<param.z_sampl; j++)	//cislo sloupce
            {
                k = j + param.z_sampl*i;
                //if( i>0 && i<M-1 && j>0 && j<N-1)
                if(grid.mask[i][j]==FIXED)
                {
                    Ax[l] = 1;
                    Ai[l] = k;
                    Ap[k+1] = Ap[k]+1;
                    l++;
                }
                else
                {
                    Ax[l] = Ax[l+1] = Ax[l+3] = Ax[l+4] = 1.0;  
                    Ax[l+2] = -4;
                    Ai[l] = k-param.z_sampl;
                    Ai[l+1] = k-1;
                    Ai[l+2] = k;
                    Ai[l+3] = k+1;
                    Ai[l+4] = k+param.z_sampl;
                    Ap[k+1] = Ap[k]+5;
                    l += 5;
                }
                /*
                else if(grid.mask[i][j]==BOUNDARY)
                {
                    Ax[l] = 2/((grid.metrika[i][j].t + grid.metrika[i][j].b) * grid.metrika[i][j].t);	//i-1, j
                    Ax[l+1] = 2/((grid.metrika[i][j].l + grid.metrika[i][j].r) * grid.metrika[i][j].l);	//i, j-1
                    Ax[l+3] = 2/((grid.metrika[i][j].l + grid.metrika[i][j].r) * grid.metrika[i][j].r);	//i, j+1
                    Ax[l+4] = 2/((grid.metrika[i][j].t + grid.metrika[i][j].b) * grid.metrika[i][j].b);	//i+1, j
                    Ax[l+2] = -2*( grid.metrika[i][j].r * grid.metrika[i][j].l + grid.metrika[i][j].t * grid.metrika[i][j].b ) / 
                        ( grid.metrika[i][j].r * grid.metrika[i][j].l * grid.metrika[i][j].t * grid.metrika[i][j].b );
                    //printf("[%d][%d] i-1: %g, j-1: %g, j+1: %g, i+1: %g\n", i, j, Ax[l], Ax[l+1], Ax[l+3], Ax[l+4]);
                    Ai[l] = k-param.y_sampl;
                    Ai[l+1] = k-1;
                    Ai[l+2] = k;
                    Ai[l+3] = k+1;
                    Ai[l+4] = k+param.y_sampl;
                    Ap[k+1] = Ap[k]+5;
                    l += 5;
                }
                */
            }
    }
    /*
    for(i=0;i<param.r_sampl;i++)
	for(j=0;j<param.z_sampl;j++)
	{
	    if( SQR(i*p_param->dr-p_param->r_max/2) + SQR(j*p_param->dz-p_param->z_max/2) <= SQR(p_param->probe_radius))
	    {
		u[i][j] = p_param->u_probe;
	    }else
	    {
		u[i][j] = 0;
	    }
	}
*/
    (void) umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL) ;
    (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL) ;
    umfpack_di_free_symbolic (&Symbolic) ;
}
    

void Fields::boundary_solve()
{
    if( p_param->coord == CYLINDRICAL )
    {
        int k;
        for(int i=0; i<grid.M; i++)		//cislo radku
            for(int j=0; j<grid.N; j++)	//cislo sloupce
            {
                k = j + grid.N*i;
                double dr = p_param->dr;
                double dz = p_param->dz;
                if(grid.mask[i][j] == FIXED)
                    rho[i][j] = grid.voltage[i][j];
                else if(i>0)
                    rho[i][j] *= -1.0/p_param->eps_0/(M_PI*dr*dr*2.0*i*dz)*p_param->macroparticle_factor;
                else
                    rho[i][j] *= -1.0/p_param->eps_0/(M_PI*dr*dr*0.25*dz)*p_param->macroparticle_factor;
            }
    }
    else
    {
        int k;
        for(int i=0; i<grid.M; i++)		//cislo radku
            for(int j=0; j<grid.N; j++)	//cislo sloupce
            {
                k = j + grid.N*i;
                if(grid.mask[i][j] == FIXED)
                    rho[i][j] = grid.voltage[i][j];
                else
                    rho[i][j] *= -SQR(p_param->dr)/p_param->eps_0/p_param->dV;
            }
    }
    (void) umfpack_di_solve (UMFPACK_At, Ap, Ai, Ax, u[0], rho[0], Numeric, NULL, NULL) ;
}

void Fields::solve()
{
    (void) umfpack_di_solve (UMFPACK_At, Ap, Ai, Ax, u[0], rho[0], Numeric, NULL, NULL) ;
}
void Fields::reset()
{
    for(int i=0; i<grid.M; i++)		//cislo radku
	for(int j=0; j<grid.N; j++)	//cislo sloupce
	    rho[i][j] = 0;
}


Fields::~Fields()
{
    umfpack_di_free_numeric (&Numeric) ;
    delete [] Ap;
    delete [] Ai;
    delete [] Ax;
}



/*
 * vypocte gradient 2d potencialu bilinearni interpolaci
 */
void Fields::E(double x, double y, double &grad_x, double &grad_y, double time)
{
    int i,j;
    double g1,g2,g3,g4;
    double fx,fy;
    double idr = p_param->idr;
    double idz = p_param->idz;

    /*
     * nejprve x-ova slozka
     */
    i = (int)(x*idr+0.5);
    j = (int)(y*idz);
    j = MIN(j,p_param->z_sampl-2);
    if(i>0 && i<p_param->r_sampl-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (u[i][j]-u[i-1][j])*idr;
	g2 = (u[i][j+1]-u[i-1][j+1])*idr;
	g3 = (u[i+1][j+1]-u[i][j+1])*idr;
	g4 = (u[i+1][j]-u[i][j])*idr;

	//interpolace
	fx = x*idr-i+.5;
	fy = y*idz-j;
	grad_x = g1*(1-fx)*(1-fy) + g2*(1-fx)*fy + g3*fx*fy + g4*fx*(1-fy);
    }else if(i==p_param->r_sampl-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (u[i][j]-u[i-1][j])*idr;
	g2 = (u[i][j+1]-u[i-1][j+1])*idr;

	//interpolace
	fy = y*idz-j;
	grad_x = g1*(1-fy) + g2*fy;
    }else if(i==0)
    {
	//vypocet gradientu v mrizovych bodech
	g3 = (u[i+1][j+1]-u[i][j+1])*idr;
	g4 = (u[i+1][j]-u[i][j])*idr;

	//interpolace
	fy = y*idz-j;
	grad_x = g3*fy + g4*(1-fy);
    }



    /*
     * ted y-ova slozka
     */
    i = (int)(x*idr);
    j = (int)(y*idz+0.5);
    i = MIN(i,p_param->r_sampl-2);


    /*
       =if(pot_mask[i][j]==OKRAJ || pot_mask[i][j+1]==OKRAJ)
       {
       }else
       */
    if(j>0 && j<p_param->z_sampl-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (u[i][j]-u[i][j-1])*idz;
	g2 = (u[i+1][j]-u[i+1][j-1])*idz;
	g3 = (u[i+1][j+1]-u[i+1][j])*idz;
	g4 = (u[i][j+1]-u[i][j])*idz;

	//interpolace
	fx = x*idr-i;
	fy = y*idz-j+0.5;
	grad_y = g1*(1-fx)*(1-fy) + g2*(1-fy)*fx + g3*fx*fy + g4*fy*(1-fx);
    }else if(j==p_param->z_sampl-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (u[i][j]-u[i][j-1])*idz;
	g2 = (u[i+1][j]-u[i+1][j-1])*idz;

	//interpolace
	fx = x*idr-i;
	grad_y = g1*(1-fx) + g2*fx;
    }else if(j==0)
    {
	//vypocet gradientu v mrizovych bodech
	g3 = (u[i+1][j+1]-u[i+1][j])*idz;
	g4 = (u[i][j+1]-u[i][j])*idz;

	//interpolace
	fx = x*idr-i;
	grad_y = g3*fx + g4*(1-fx);
    }

    if(p_param->rf)
    {
        double phase = sin(p_param->rf_omega*time);
        grad_x *= phase;
        grad_y *= phase;
    }

}
void Fields::B(double x, double y, double &_Br, double &_Bz, double &_Bt)
{
    if(p_param->magnetic_field_const)
    {
        /*
           double K = (0.03-0.003)/(4.1-3.0)*1e2;
           double a = 0.03-K*4.1*1e-2;
           Bz = K*I->z + a;
           Br = -K*0.5*I->r;
           if(I->z<3e-2)
           {
           Bz = 0.003;
           Br = 0;
           }
           */
        _Br = p_param->Br;
        _Bz = p_param->Bz;
        _Bt = p_param->Bt;
    }
    else
    {
        _Br = Br.interpolate(x, y);
        _Bz = Bz.interpolate(x, y);
        _Bt = 0.00;
    }
}
string double2string(double x)
{
    ostringstream o;
    o << x;
    return o.str();
}
int double2int(double x, double eps=1e-2)
{
    int res = (int)(x+0.5);
    if(fabs(res-x) > eps)
	throw std::runtime_error("double2int() " + double2string(x) + " is not integer\n");
    return res;
}
void Fields::load_magnetic_field(const char * fname)
{
    vector<double> rvec, zvec, Brvec, Bzvec;
    double tmp1, tmp2, tmp3, tmp4;

    std::ifstream fr(fname);
    string line;
    istringstream s_line;

    //load the numbers from file
    while(fr.good())
    {
        getline(fr, line);
        s_line.clear();
        s_line.str(line);

        if( s_line >> tmp1 >> tmp2 >> tmp3 >> tmp4 )
        {
            rvec.push_back(tmp1);
            zvec.push_back(tmp2);
            Brvec.push_back(tmp3);
            Bzvec.push_back(tmp4);
        }
    }

    //find dr, rmin, rmax
    double dr=0, rmin=rvec[0], rmax=rvec[rvec.size()-1];
    unsigned int i=1;
    while(i<rvec.size() && rvec[i]-rvec[i-1] == 0.0)
        i++;
    dr = rvec[i]-rvec[i-1];
    if(dr<0)
    {
        dr = -dr;
        rmin = rvec[rvec.size()-1];
        rmax = rvec[0];
    }
    double rsampl_d = (rmax-rmin)/dr+1;
    unsigned int rsampl = double2int(rsampl_d);

    //find dz, zmin, zmax
    double dz=0, zmin=zvec[0], zmax=zvec[rvec.size()-1];
    while(i<zvec.size() && zvec[i]-zvec[i-1] == 0.0)
        i++;
    dz = zvec[i]-zvec[i-1];
    if(dz<0)
    {
        dz = -dz;
        zmin = zvec[zvec.size()-1];
        zmax = zvec[0];
    }
    double zsampl_d = (zmax-zmin)/dz+1;
    unsigned int zsampl = double2int(zsampl_d);
    if(rsampl*zsampl != rvec.size())
        throw std::runtime_error("Fields::load_magnetic_field() wrong size of input vector");

    Br.resize(rsampl, zsampl, dr, dz, zmin, zmax);
    Bz.resize(rsampl, zsampl, dr, dz, zmin, zmax);

    for(i=0; i<rsampl; i++)
        for(unsigned int j=0; j<zsampl; j++)
        {
            Br[i][j] = numeric_limits<double>::quiet_NaN();
            Bz[i][j] = numeric_limits<double>::quiet_NaN();
        }

    int ri, zi;
    for(i=0; i<rvec.size(); i++)
    {
        ri = double2int(rvec[i]/dr);
        zi = double2int(zvec[i]/dz);
        Br[ri][zi] = Brvec[i];
        Bz[ri][zi] = Bzvec[i];
    }

    for(i=0; i<rsampl; i++)
        for(unsigned int j=0; j<zsampl; j++)
            if(isnan(Br[i][j]) || isnan(Bz[i][j]))
                throw std::runtime_error("Fields::load_magnetic_field() garbage loaded");
}
#endif
