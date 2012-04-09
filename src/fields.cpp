#include <cmath>
#include <suitesparse/umfpack.h>
#include "fields.hpp"

/*
 *********************** definitions for Field ********************************
 */
void Field::resize(int x_sampl, int z_sampl, double dx, double dz, double _rmin, double _zmin)
{
    idx = 1.0/dx, idz = 1.0/dz;
    rmin = _rmin, zmin = _zmin;
    Matrix<double>::resize(x_sampl, z_sampl);
}

bool Field::hasnan()
{
    for(int j=0; j<jmax; j++)
        for(int l=0; l<lmax; l++)
            if(isnan(data[j][l]))
                return true;
    return false;
}

void Field::print( ostream & out , double factor)
{
    double dx=1.0/idx, dz=1.0/idz;
    for(int j=0; j<jmax; j++)
    {
        for(int l=0; l<lmax; l++)
            out << j*dx <<"\t"<< l*dz <<"\t"<< data[j][l]*factor <<endl;
        out << endl;
    }
}

void Field::print( const char * filename , double factor)
{
    ofstream out(filename);
    print(out, factor);
}


/*
 *********************** definitions for Fields *******************************
 */

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

void Fields::u_smooth()
{
    int ic = u.jmax-1;
    int jc = u.lmax-1;
    if(ic != jc)
        throw std::runtime_error("Fields::u_smooth: smoothing of non-square matrix not implemented\n");


    // this averages the symmetrical points in rotationally symetric
    // geometry on a square cartesian grid
    for(int i=0; i<=ic/2; i++)
        for(int j=0; j<=i; j++)
        {
            double sum = 0;

            //calculate the average
            sum += u[i][j];
            sum += u[ic-i][j];
            sum += u[i][jc-j];
            sum += u[ic-i][jc-j];
            sum += u[jc-j][ic-i];
            sum += u[j][i];
            sum += u[jc-j][i];
            sum += u[j][ic-i];
            sum /= 8.0;

            //assign to corresponding cells
            u[i][j] = sum;
            u[ic-i][j] = sum;
            u[i][jc-j] = sum;
            u[ic-i][jc-j] = sum;
            u[jc-j][ic-i] = sum;
            u[j][i] = sum;
            u[jc-j][i] = sum;
            u[j][ic-i] = sum;
        }
    // TODO remove this magic constant :)
    double radius = p_param->probe_radius;

    // avoid smoothing of the boundary region with sharp gradients
    radius -= 2.0*p_param->dx;

    // precalculate something ;)
    radius = sqr(radius/p_param->dx);

    uTmp.assign(u);

    double icf = ic/2.0;
    double jcf = jc/2.0;
    for(int i=0; i<=ic; i++)
        for(int j=0; j<=jc; j++)
        {
            double r = sqr(i-icf) + sqr(j-jcf);
            if(r > radius) continue;

            // use the HCIC smoothing scheme of Hockney 1971 ?
            // http://dx.doi.org/10.1016/0021-9991(71)90032-5
            // This scheme was becoming unstable with respect to
            // particular spatial frequencies (dx*2)*(dy*2) at low electron
            // temperatures. Following scheme corresponds to 2dx*2dy
            // square particle.
            double sum =
                uTmp[i][j] +
                uTmp[i-1][j]*0.5 +
                uTmp[i+1][j]*0.5 +
                uTmp[i][j-1]*0.5 +
                uTmp[i][j+1]*0.5 +
                uTmp[i-1][j-1]*0.25 +
                uTmp[i+1][j-1]*0.25 +
                uTmp[i-1][j+1]*0.25 +
                uTmp[i+1][j+1]*0.25;

            u[i][j] = sum/4.0;
        }

}

Fields::Fields(Param &param) : grid(param), u(param), uRF(param), uTmp(param), uAvg(param), rho(param), nsampl(0)
{
    //if(param.selfconsistent == false)
//	return;
    p_param = &param;
    int i, j, k, l;
    int n;
    idx=p_param->idx;
    idz=p_param->idz;
    double dx = p_param->dx;
    double dz = p_param->dz;

    n = param.x_sampl * param.z_sampl;

    // r_i = i*dx
    Ap = new int [n+1];
    Ai = new int [n*5]; //horni odhad pro petibodove diff schema
    Ax = new double [n*5];

    if( param.coord == CYLINDRICAL )
    {

        l = 0;
        Ap[0] = 0;
        for(i=0; i<param.x_sampl; i++)	//cislo radku - radialni souradnice
            for(j=0; j<param.z_sampl; j++)	//cislo sloupce
            {
                k = j + param.z_sampl*i;
                if(grid.mask[i][j]==FIXED || grid.mask[i][j]==FIXED_RF)
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
                    k3 = 1.0/(dx*dx*0.25);
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
                else if(grid.mask[i][j]!=FIXED && grid.mask[i][j]!=FIXED_RF)
                {
                    double k1, k2, k3;
                    // psi[j-1,k]
                    k1 = (i-0.5)/(dx*dx*i);
                    Ax[l] = k1;
                    Ai[l] = k-param.z_sampl;

                    // psi[j,k-1]
                    k2 = 1.0/(dz*dz);
                    Ax[l+1] = k2;
                    Ai[l+1] = k-1;

                    // psi[j,k]
                    k3 = (i+0.5)/(dx*dx*i);
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
        for(i=0; i<param.x_sampl; i++)		//cislo radku
            for(j=0; j<param.z_sampl; j++)	//cislo sloupce
            {
                k = j + param.z_sampl*i;
                //if( i>0 && i<M-1 && j>0 && j<N-1)
                if(grid.mask[i][j]==FIXED || grid.mask[i][j]==FIXED_RF)
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
    for(i=0;i<param.x_sampl;i++)
	for(j=0;j<param.z_sampl;j++)
	{
	    if( SQR(i*p_param->dx-p_param->r_max/2) + SQR(j*p_param->dz-p_param->z_max/2) <= SQR(p_param->probe_radius))
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
                double dx = p_param->dx;
                double dz = p_param->dz;
                if(grid.mask[i][j] == FIXED)
                    rho[i][j] = grid.voltage[i][j];
                else if(grid.mask[i][j] == FIXED_RF)
                    rho[i][j] = 0;      //use zero as first approximation
                else if(i>0)
                    rho[i][j] *= -1.0/p_param->eps_0/(M_PI*dx*dx*2.0*i*dz)*p_param->macroparticle_factor;
                else
                    rho[i][j] *= -1.0/p_param->eps_0/(M_PI*dx*dx*0.25*dz)*p_param->macroparticle_factor;
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
                else if(grid.mask[i][j] == FIXED_RF)
                    rho[i][j] = 0;
                else
                    rho[i][j] *= -SQR(p_param->dx)/p_param->eps_0/p_param->dV;
                    //  the SQR(p_param->dx) comes from the d^2/dx^2 operator
            }
    }
    (void) umfpack_di_solve (UMFPACK_At, Ap, Ai, Ax, u[0], rho[0], Numeric, NULL, NULL) ;
}

void Fields::boundary_solve_rf()
{// XXX this should be later merged to Fields::boundary_solve
    if( p_param->coord == CYLINDRICAL )
    {
        int k;
        for(int i=0; i<grid.M; i++)		//cislo radku
            for(int j=0; j<grid.N; j++)	//cislo sloupce
            {
                k = j + grid.N*i;
                double dx = p_param->dx;
                double dz = p_param->dz;
                if(grid.mask[i][j] == FIXED)
                    rho[i][j] = 0;
                else if(grid.mask[i][j] == FIXED_RF)
                    rho[i][j] = grid.voltage[i][j];      //use zero as first approximation
                else if(i>0)
                    rho[i][j] *= -1.0/p_param->eps_0/(M_PI*dx*dx*2.0*i*dz)*p_param->macroparticle_factor;
                else
                    rho[i][j] *= -1.0/p_param->eps_0/(M_PI*dx*dx*0.25*dz)*p_param->macroparticle_factor;
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
                    rho[i][j] = 0;
                else if(grid.mask[i][j] == FIXED_RF)
                    rho[i][j] = grid.voltage[i][j];
                else
                    rho[i][j] *= -SQR(p_param->dx)/p_param->eps_0/p_param->dV;
                    //  the SQR(p_param->dx) comes from the d^2/dx^2 operator
            }
    }
    (void) umfpack_di_solve (UMFPACK_At, Ap, Ai, Ax, uRF[0], rho[0], Numeric, NULL, NULL) ;
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
 *********************** definitions for t_grid *******************************
 */

t_grid::t_grid(Param &param) :  M(param.x_sampl), N(param.z_sampl), 
    dx(param.dx), dz(param.dz),
    mask(param.x_sampl,param.z_sampl), metrika(param.x_sampl,param.z_sampl), voltage(param.x_sampl,param.z_sampl)
{
    p_param = &param;
			
    switch(param.geometry)
    {
        case Param::PENNING_SIMPLE:
            penning_trap_simple();
            break;
        case Param::PENNING:
            penning_trap();
            break;
        case Param::MAC:
            MAC_filter();
            break;
        case Param::RF_QUAD:
            rf_trap();
            break;
        case Param::RF_HAITRAP:
            rf_haitrap();
            break;
        case Param::RF_8PT:
            rf_8PT();
            break;
        case Param::RF_22PT:
            rf_22PT();
            break;
        case Param::TUBE:
            tube();
            break;
        case Param::EMPTY:
        default:
            empty();
            break;
    }

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
        circle_electrode(x, y, r_rod, sign, FIXED_RF);
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

void t_grid::rf_8PT()
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
		voltage[i][j] = p_param->extern_field*dz*(j-N/2);
	    }else
	    {
		mask[i][j] = FREE;
	    }
	}

    double xcenter = 1e-2;
    double ycenter = 1e-2;
    double r_rod = 0.1e-2;
    double r_8pt = 0.3e-2 + r_rod;
    int npoles = 8;
    for(int i=0; i<npoles; i++)
    {
        double x = xcenter + sin(2*M_PI*(i+1.0/32)/npoles)*r_8pt;
        double y = ycenter + cos(2*M_PI*(i+1.0/32)/npoles)*r_8pt;
        int sign = i%2==0 ? -1 : 1;
        circle_electrode(x, y, r_rod, sign, FIXED_RF);
    }

    for(i=2;i<M-2;i++)
	for(j=2;j<N-2;j++)
	{// should implement BOUNDARY for RF electrodes as well, or throw it away completely
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
void t_grid::rf_haitrap()
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
		voltage[i][j] = p_param->extern_field*dz*(j-N/2);
	    }else
	    {
		mask[i][j] = FREE;
	    }
	}

    double xcenter = 1e-2;
    double ycenter = 1e-2;
    double r_rod = 0.01e-2;
    double r_8pt = 0.3e-2 + r_rod;
    int npoles = 8;
    for(int i=0; i<npoles; i++)
    {
        double x = xcenter + sin(2*M_PI*(i+1.0/32)/npoles)*r_8pt;
        double y = ycenter + cos(2*M_PI*(i+1.0/32)/npoles)*r_8pt;
        int sign = i%2==0 ? -1 : 1;
        circle_electrode(x, y, r_rod, sign, FIXED_RF);
    }

    for(i=2;i<M-2;i++)
	for(j=2;j<N-2;j++)
	{// should implement BOUNDARY for RF electrodes as well, or throw it away completely
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

void t_grid::tube()
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

    double radius = p_param->probe_radius;
    double rcenter, zcenter;
    rcenter = zcenter = p_param->x_max/2.0;
    double sqradius = SQR(radius);
    for(i=0; i<M; i++)
	for(j=0; j<N; j++)
	{
	    double r, z;
	    r = i*dx - rcenter;
	    z = j*dz - zcenter;
	    // suboptimal, but simple
	    if(SQR(r) + SQR(z) >= sqradius)
	    {
		mask[i][j] = FIXED;
		voltage[i][j] = 0.0;
	    }
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

    circle_electrode(5e-3, 1e-2, 2e-3, 1.0, FIXED_RF);
    circle_electrode(15e-3, 1e-2, 2e-3, 1.0, FIXED_RF);
    circle_electrode(1e-2, 5e-3, 2e-3, -1.0, FIXED_RF);
    circle_electrode(1e-2, 15e-3, 2e-3, -1.0, FIXED_RF);

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
    square_electrode(5e-3, 4.5e-2, 1e-2, 1.5e-2, -.00);
    square_electrode(5e-3, 7e-3, 2e-2, 8e-2, 0.0);
    square_electrode(5e-3, 4.5e-2, 8.5e-2, 9e-2, -.00);

    double threshold = p_param->u_probe;
    double ofs = 3e-2;
    square_electrode(3e-2, 3.3e-2, 11e-2+ofs, 14e-2+ofs, 0.8*threshold);
    square_electrode(4.5e-2, 4.8e-2, 15e-2+ofs, 25e-2+ofs, threshold);
    square_electrode(3e-2, 3.3e-2, 26e-2+ofs, 29e-2+ofs, 1.0*threshold);
    square_electrode(2.5e-2, 2.8e-2, 29e-2+ofs, 30.5e-2+ofs, 1.0*threshold);

    square_electrode(15e-3, 4.5e-2, 35e-2, 35.3e-2, .0);
    //probe
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
		voltage[i][j] = -i*p_param->dx*p_param->extern_field + p_param->x_max*p_param->extern_field*0.5;
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

void t_grid::square_electrode(double rmin, double rmax, double zmin, double zmax, double _voltage, char mask_type)
{
    for(int i=0; i<M; i++)
	for(int j=0; j<N; j++)
	{
	    double r, z;
	    r = i*dx;
	    z = j*dz;
	    // suboptimal, but simple
	    if(r>rmin && r<rmax && z>zmin && z<zmax) 
	    {
		mask[i][j] = mask_type;
		voltage[i][j] = _voltage;
	    }
	}
}

void t_grid::circle_electrode(double rcenter, double zcenter, double radius, double _voltage, char mask_type)
{
    double sqradius = SQR(radius);
    for(int i=0; i<M; i++)
	for(int j=0; j<N; j++)
	{
	    double r, z;
	    r = i*dx - rcenter;
	    z = j*dz - zcenter;
	    // suboptimal, but simple
	    if(SQR(r) + SQR(z) <= sqradius) 
	    {
		mask[i][j] = mask_type;
		voltage[i][j] = _voltage;
	    }
	}
}
    
void Fields::load_magnetic_field(const char * fname)
{
    vector<double> rvec, zvec, Brvec, Bzvec;
    double tmp1, tmp2, tmp3, tmp4;

    std::ifstream fr(fname);
    if(fr.fail())
    {
        throw runtime_error("Fields::load_magnetic_field(): failed opening file\n");
    }
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

    //find dx, rmin, rmax
    double dx=0, rmin=rvec[0], rmax=rvec[rvec.size()-1];
    unsigned int i=1;
    while(i<rvec.size() && rvec[i]-rvec[i-1] == 0.0)
        i++;
    dx = rvec[i]-rvec[i-1];
    if(dx<0)
    {
        dx = -dx;
        rmin = rvec[rvec.size()-1];
        rmax = rvec[0];
    }
    double rsampl_d = (rmax-rmin)/dx+1;
    unsigned int rsampl = double2int(rsampl_d);

    //find dz, zmin, zmax
    double dz=0, zmin=zvec[0], zmax=zvec[rvec.size()-1];
    i = 1;
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
#ifdef DEBUG
    cout << "rvec.size " << rvec.size() <<endl;
    cout << "rmin = " << rmin <<"  rmax = "<< rmax <<"  dx = "<<dx <<"  rsampl = "<<rsampl<<endl;
    cout << "zmin = " << zmin <<"  zmax = "<< zmax <<"  dz = "<<dz <<"  zsampl = "<<zsampl <<endl;
#endif
    if(rsampl*zsampl != rvec.size())
        throw std::runtime_error("Fields::load_magnetic_field() wrong size of input vector");

    Br.resize(rsampl, zsampl, dx, dz, rmin, zmin);
    Bz.resize(rsampl, zsampl, dx, dz, rmin, zmin);

    for(i=0; i<rsampl; i++)
        for(unsigned int j=0; j<zsampl; j++)
        {
            Br[i][j] = numeric_limits<double>::quiet_NaN();
            Bz[i][j] = numeric_limits<double>::quiet_NaN();
        }

    int ri, zi;
    for(i=0; i<rvec.size(); i++)
    {
        ri = double2int((rvec[i]-rmin)/dx);
        zi = double2int((zvec[i]-zmin)/dz);
        Br[ri][zi] = Brvec[i];
        Bz[ri][zi] = Bzvec[i];
    }

    for(i=0; i<rsampl; i++)
        for(unsigned int j=0; j<zsampl; j++)
            if(isnan(Br[i][j]) || isnan(Bz[i][j]))
                throw std::runtime_error("Fields::load_magnetic_field() garbage loaded");
}
