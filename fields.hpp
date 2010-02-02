#ifndef FIELD_H
#define FIELD_H
//#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)>(y)?(y):(x))


#include "param.hpp"
#include "matrix.cpp"
#include "util.cpp"
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
	double dx,dz;
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
        void rf_8PT();
        void MAC_filter();
        void empty();
	void square_electrode(double rmin, double rmax, double zmin, double zmax, double voltage);
	void circle_electrode(double xcenter, double ycenter, double radius, double voltage);
	bool is_free(double r, double z);
};
class Field : public Matrix<double>
{
    public:
	Field(Param &param) : Matrix<double>(param.x_sampl, param.z_sampl),
            idx(param.idx), idz(param.idz), rmin(0), zmin(0) {};
        Field() : Matrix<double>(), idx(0), idz(0), rmin(0), zmin(0) {};
        void resize(int x_sampl, int z_sampl, double dx, double dz, double _rmin=0, double _zmin=0);
	inline void accumulate(double charge, double x, double y);
        inline double interpolate(double r, double z);
        bool hasnan();
	void print( ostream & out = cout , double factor = 1.0);
	void print( const char * filename , double factor = 1.0);
    private:
	double idx, idz;
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
    private:
	int nsampl;
};

inline bool t_grid::is_free(double r, double z)
{
    int i = (int)(r*p_param->idx);
    int j = (int)(z*p_param->idz);
    if(mask[i][j]==FREE || mask[i+1][j]==FREE || mask[i][j+1]==FREE || mask[i+1][j+1]==FREE)
	return true;
    else return false;
}

inline void Field::accumulate(double charge, double r, double z)
{
    r -= rmin;
    z -= zmin;
    int i = (int)(r * idx);
    int j = (int)(z * idz);

    double u = r*idx - i;
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
    r -= rmin;
    z -= zmin;
    int i = (int)(r * idx);
    int j = (int)(z * idz);

    double u = r*idx - i;
    double v = z*idz - j;

    if(i<0 || i>jmax-1 || j<0 || j>lmax-1)
	throw std::runtime_error("Field::interpolate() outside of range\n");

    return (1-u)*(1-v)*data[i][j] + u*(1-v)*data[i+1][j] + (1-u)*v*data[i][j+1] + u*v*data[i+1][j+1];
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
inline void Fields::E(double x, double y, double &grad_x, double &grad_y, double time)
{
    int i,j;
    double g1,g2,g3,g4;
    double fx,fy;
    double idx = p_param->idx;
    double idz = p_param->idz;

    /*
     * nejprve x-ova slozka
     */
    i = (int)(x*idx+0.5);
    j = (int)(y*idz);
    j = MIN(j,p_param->z_sampl-2);
    if(i>0 && i<p_param->x_sampl-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (u[i][j]-u[i-1][j])*idx;
	g2 = (u[i][j+1]-u[i-1][j+1])*idx;
	g3 = (u[i+1][j+1]-u[i][j+1])*idx;
	g4 = (u[i+1][j]-u[i][j])*idx;

	//interpolace
	fx = x*idx-i+.5;
	fy = y*idz-j;
	grad_x = g1*(1-fx)*(1-fy) + g2*(1-fx)*fy + g3*fx*fy + g4*fx*(1-fy);
    }else if(i==p_param->x_sampl-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (u[i][j]-u[i-1][j])*idx;
	g2 = (u[i][j+1]-u[i-1][j+1])*idx;

	//interpolace
	fy = y*idz-j;
	grad_x = g1*(1-fy) + g2*fy;
    }else if(i==0)
    {
	//vypocet gradientu v mrizovych bodech
	g3 = (u[i+1][j+1]-u[i][j+1])*idx;
	g4 = (u[i+1][j]-u[i][j])*idx;

	//interpolace
	fy = y*idz-j;
	grad_x = g3*fy + g4*(1-fy);
    }



    /*
     * ted y-ova slozka
     */
    i = (int)(x*idx);
    j = (int)(y*idz+0.5);
    i = MIN(i,p_param->x_sampl-2);


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
	fx = x*idx-i;
	fy = y*idz-j+0.5;
	grad_y = g1*(1-fx)*(1-fy) + g2*(1-fy)*fx + g3*fx*fy + g4*fy*(1-fx);
    }else if(j==p_param->z_sampl-1)
    {
	//vypocet gradientu v mrizovych bodech
	g1 = (u[i][j]-u[i][j-1])*idz;
	g2 = (u[i+1][j]-u[i+1][j-1])*idz;

	//interpolace
	fx = x*idx-i;
	grad_y = g1*(1-fx) + g2*fx;
    }else if(j==0)
    {
	//vypocet gradientu v mrizovych bodech
	g3 = (u[i+1][j+1]-u[i+1][j])*idz;
	g4 = (u[i][j+1]-u[i][j])*idz;

	//interpolace
	fx = x*idx-i;
	grad_y = g3*fx + g4*(1-fx);
    }

    if(p_param->rf)
    {
        double phase = sin(p_param->rf_omega*time);
        grad_x *= phase;
        grad_y *= phase;
    }

}

inline void Fields::B(double x, double y, double &_Br, double &_Bz, double &_Bt)
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
