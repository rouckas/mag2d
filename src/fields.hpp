#ifndef FIELD_H
#define FIELD_H
//#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)>(y)?(y):(x))


#include "param.hpp"
#include "Field2D.hpp"
#include "util.cpp"
#include "mymath.cpp"
#include <sstream>
#include <string>
#include <limits>
using namespace std;

//#define SQR(x) ((x)*(x))
inline double SQR(double x){return x*x;};

enum {FIXED, FIXED_RF, FREE, BOUNDARY};


typedef struct{
            double l, r, t, b;
}t_umf_metrika;



class t_grid
//definuje geometrii experimentu
{
    public:
	int M;
	int N;
	double dx,dz;
	Array2D<char> mask;
	Array2D<t_umf_metrika> metrika;
	Array2D<double> voltage;
	t_grid(Param &param);
	Param *p_param;
        void penning_trap();
        void penning_trap_simple(double trap_voltage = -1.0);
        double U_trap;
        void tube();
        void rf_trap();
        void rf_22PT();
        void rf_8PT();
        void rf_haitrap();
        void MAC_filter();
        void empty();
	void square_electrode(double rmin, double rmax, double zmin, double zmax, double voltage, char mask_type = FIXED);
	void circle_electrode(double xcenter, double ycenter, double radius, double voltage, char mask_type = FIXED);
	bool is_free(double r, double z) const;
};

class Fields
//stara se o reseni Poissonovy rce na danem gridu
{
    public:
	t_grid grid;
	Field2D u;
	Field2D uRF;
        Field2D uTmp;
	Field2D uAvg;
	Field2D rho;
	Fields(Param &param);
	void solve();
	void boundary_solve();
	void boundary_solve_rf();
	void reset();
	void E(double x, double y, double &grad_x, double &grad_y, double time = 0) const;
	void B(double x, double y, double &Br, double &Bz, double &Bt) const;
        void load_magnetic_field(const char * fname);
	inline void accumulate(double charge, double x, double y);
	~Fields();
    private:
        Field2D Br, Bz;
	double idx, idz;
	Param *p_param;
	int *Ap;
	int *Ai;
	double *Ax;
	void *Symbolic, *Numeric ;
    public:
	void u_sample();
	void u_reset();
	void u_print(const char * fname);
        void u_smooth(bool symmetry = 0, double radius = -1., Field2D * field = NULL);
    private:
	int nsampl;
};

inline bool t_grid::is_free(double r, double z) const
{
    int i = (int)(r*p_param->idx);
    int j = (int)(z*p_param->idz);
    if(mask[i][j]==FREE || mask[i+1][j]==FREE || mask[i][j+1]==FREE || mask[i+1][j+1]==FREE)
	return true;
    else return false;
}


inline void Fields::accumulate(double charge, double r, double z)
{
    int i = (int)(r * idx);
    int j = (int)(z * idz);

    double u = r*idx - i;
    double v = z*idz - j;

    //if(i<0 || i>p_param->x_sampl-1 || j<0 || j>p_param->z_sampl-1)
//	throw std::runtime_error("Field::accumulate() outside of range\n");
    rho[i][j] += (1-u)*(1-v)*charge;
    rho[i+1][j] += u*(1-v)*charge;
    rho[i][j+1] += (1-u)*v*charge;
    rho[i+1][j+1] += u*v*charge;

}

/*
 * vypocte gradient 2d potencialu bilinearni interpolaci
 */
inline void Fields::E(double x, double y, double &grad_x, double &grad_y, double time) const
{
    if(p_param->geometry == Param::EMPTY && p_param->selfconsistent == false)
    {
        grad_x = 0.;
        grad_y = p_param->extern_field;
        return;
    }

    u.grad(x, y, grad_x, grad_y);

    if(p_param->rf)
    {
        double grad_rf_x, grad_rf_y;
        uRF.grad(x, y, grad_rf_x, grad_rf_y);

        double phase = mod(p_param->rf_omega*time, 10000000*M_PI);
        phase = p_param->rf_amplitude*cos(phase) +
            p_param->rf_U0;
        grad_x += grad_rf_x*phase;
        grad_y += grad_rf_y*phase;
    }

}

inline void Fields::B(double x, double y, double &_Br, double &_Bz, double &_Bt) const
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
#endif
