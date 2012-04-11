#include <suitesparse/umfpack.h>
#include "fields3d.hpp"
using namespace std;

void Electrode::set_voltage(Array3D<double> &voltage_array, Array3D<signed char> &mask)
{
    for(int i = 0; i < voltage_array.imax*voltage_array.jmax*voltage_array.kmax; i++)
        if(mask[0][0][i] == id)
            voltage_array[0][0][i] = voltage;
}

Geometry::Geometry(Param &param) :
    idx(param.idx), idy(param.idy), idz(param.idz),
    x_sampl(param.x_sampl), y_sampl(param.y_sampl),
    z_sampl(param.z_sampl), mask(param.x_sampl, param.y_sampl, param.z_sampl),
    electrodes(0),
    voltage(param.x_sampl, param.y_sampl, param.z_sampl, param.dx, param.dy, param.dz)
{

    //set zero Dirichlet boundary condition
    for(int i = 0; i < param.x_sampl; i++)
        for(int j = 0; j < param.y_sampl; j++)
            for(int k = 0; k < param.z_sampl; k++)
                if( i==0 || i==param.x_sampl-1 || j==0 || j==param.y_sampl-1 || k==0 || k==param.z_sampl)
                    mask[i][j][k] = FIXED;
                else
                    mask[i][j][k] = FREE;

    electrodes.push_back(new Quadrupole(-1));
    electrodes[0]->set_mask(mask);
    electrodes[0]->set_voltage(voltage, mask);

}

void Solver::matrix_init()
{
    int n = x_sampl * y_sampl * z_sampl;
    Ap = new int [n+1];
    Ai = new int [n*7]; //upper limit for 7-point diff scheme
    Ax = new double [n*7];


    int l = 0, m;
    Ap[0] = 0;
    for(int i=0; i<x_sampl; i++)
        for(int j=0; j<y_sampl; j++)
            for(int k=0; k<z_sampl; k++)
            {
                m = y_sampl*z_sampl*i + z_sampl*j + k;
                if(geometry.mask[i][j][k] == FIXED || geometry.mask[i][j][k] < 0)
                    //XXX is the FIXED still needed or is it replaced by electrode id?
                {
                    Ax[l] = 1;
                    Ai[l] = m;
                    Ap[m+1] = Ap[m]+1;
                    l++;
                }
                else
                {
                    Ax[l] = Ax[l+1] = Ax[l+2] = 1.0;
                    Ax[l+4] = Ax[l+5] = Ax[l+6] = 1.0;
                    Ax[l+3] = -6;
                    Ai[l] = m-y_sampl*z_sampl;
                    Ai[l+1] = m-z_sampl;
                    Ai[l+2] = m-1;
                    Ai[l+3] = m;
                    Ai[l+4] = m+1;
                    Ai[l+5] = m+z_sampl;
                    Ai[l+6] = m+y_sampl*z_sampl;
                    Ap[m+1] = Ap[m]+7;
                    l += 7;
                }
            }
};
void Solver::umfpack_init()
{
    int n = x_sampl * y_sampl * z_sampl;
    (void) umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL) ;
    (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL) ;
    umfpack_di_free_symbolic (&Symbolic) ;
};

void Solver::solve(Array3D<double> &u, const Array3D<double> &voltage, Array3D<double> &rho)
{
    for(int i=0; i<x_sampl; i++)
        for(int j=0; j<y_sampl; j++)
            for(int k=0; k<z_sampl; k++)
            {
                if(geometry.mask[i][j][k] == FIXED || geometry.mask[i][j][k] < 0)
                    rho[i][j][k] = voltage[i][j][k];
                else
                    rho[i][j][k] *= -param.macroparticle_factor/param.eps_0;
            }
    (void) umfpack_di_solve (UMFPACK_At, Ap, Ai, Ax, u[0][0], rho[0][0], Numeric, NULL, NULL) ;
}

void Solver::solve_laplace(Array3D<double> &u, const Array3D<double> &voltage)
    // solver for Laplace equation (zero RHS), it assumes that voltage is set to zero in free space
{
    (void) umfpack_di_solve (UMFPACK_At, Ap, Ai, Ax, u[0][0], voltage[0][0], Numeric, NULL, NULL) ;
}
void Solver::save(string filename)
{
    //ugly trick to override constantness of c_str()
    //when passing it as non-const param of save_numeric
    char *fn = (char*)((void*)filename.c_str());

    UF_long status;
    status =  umfpack_di_save_numeric(Numeric, fn);
    if(status == UMFPACK_ERROR_invalid_Numeric_object)
        throw runtime_error("Solver::save: UMFPACK_ERROR_invalid_Numeric_object\n");
    else if(status == UMFPACK_ERROR_file_IO)
        throw runtime_error("Solver::save: UMFPACK_ERROR_file_IO\n");
}
void Solver::load(string filename)
{
    char *fn = (char*)((void*)filename.c_str());
    UF_long status;
    status = umfpack_di_load_numeric(&Numeric, fn);
    if(status == UMFPACK_ERROR_out_of_memory)
        throw runtime_error("Solver::Solver: UMFPACK_ERROR_out_of_memory\n");
    else if(status == UMFPACK_ERROR_file_IO)
        throw runtime_error("Solver::Solver: UMFPACK_ERROR_file_IO\n");
    else if(status != UMFPACK_OK)
        throw runtime_error("Solver::load: Error loading Numeric object\n");
}
Solver::Solver(Geometry & _geometry, Param & _param, string filename) : param(_param), geometry(_geometry),
    x_sampl(geometry.x_sampl), y_sampl(geometry.y_sampl), z_sampl(geometry.z_sampl)
{
    matrix_init();

    if(filename == "")
    {
        umfpack_init();
    }
    else
    {
        load(filename);
    }
}

ElMag3D::ElMag3D(Param &param) :
    u(param.x_sampl, param.y_sampl, param.z_sampl, param.dx, param.dy, param.dz),
    rho(param.x_sampl, param.y_sampl, param.z_sampl, param.dx, param.dy, param.dz),
    voltage(param.x_sampl, param.y_sampl, param.z_sampl, param.dx, param.dy, param.dz),
    Bx(param.x_sampl, param.y_sampl, param.z_sampl, param.dx, param.dy, param.dz),
    By(param.x_sampl, param.y_sampl, param.z_sampl, param.dx, param.dy, param.dz),
    Bz(param.x_sampl, param.y_sampl, param.z_sampl, param.dx, param.dy, param.dz),
    geometry(param), potentials(0), multielectrode(false),
    solver(geometry, param, "solver_numeric.umf")
{
    //do we have multiple independent sets of electrodes in
    //nonselfconsistent simulation
    if(geometry.electrodes.size() > 1 && !param.selfconsistent)
    {
        for(unsigned int i=0; i<geometry.electrodes.size(); i++)
            potentials.push_back(new
                    Field3D(param.x_sampl, param.y_sampl, param.z_sampl, param.dx, param.dy, param.dz)
                    );
        multielectrode = true;
    }
    rho.reset();
}

void ElMag3D::solve()
{
    if(multielectrode)
    {

        for(unsigned int i=0; i<geometry.electrodes.size(); i++)
        {
            voltage.reset();
            geometry.electrodes[i]->set_voltage(voltage, geometry.mask);
            solver.solve_laplace(*(potentials[i]), voltage);
        }
    }
    else
    {
        voltage.reset();
        for(unsigned int i=0; i<geometry.electrodes.size(); i++)
            geometry.electrodes[i]->set_voltage(voltage, geometry.mask);
        solver.solve(u, voltage, rho);
    }
}

void ElMag3D::potential_sum()
    //TODO do a weighted sum of solutions for single electrodes,
    //for now it just performs sum
{
    u.reset();
    for(unsigned int i=0; i<potentials.size(); i++)
        u.add(*(potentials[i]));
}


void get_vector_sampling( vector<double> & vec, double & dx, double & min, double & max, unsigned int & sampl)
{

    //find dx, rmin, rmax
    dx=0, min=vec[0], max=vec[vec.size()-1];
    unsigned int i=1;
    while(i<vec.size() && vec[i]-vec[i-1] == 0.0)
        i++;
    dx = vec[i]-vec[i-1];
    if(dx<0)
    {
        dx = -dx;
        min = vec[vec.size()-1];
        max = vec[0];
    }
    double sampl_d = (max-min)/dx+1;
    sampl = double2int(sampl_d);

}


