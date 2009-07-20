#include "fields3d.hpp"
#include "species3d.hpp"
//#include "matrix.cpp"
#include "output.cpp"
#include "elon.cpp"
#include <GetPot>

class Elon : public Species<CARTESIAN3D>
{
    public:
    Elon(Param &param, t_random &random, ElMag3D &_field, BaseSpecies * _species_list[])
    : Species<CARTESIAN3D>(param, random, _field, 9.109534e-31, _species_list, ELECTRON) {};

    void scatter(t_particle & particle) {};
};

int main(int argc, char * argv[])
{
    GetPot cl(argc, argv);

    // load the configuration from file
    string config_file = cl("config","config.txt");
    GetPot config(config_file.c_str());

    // initialize model parameters
    Param param(config);

    // prepare output directory
    param.output_dir = cl("output_dir","output");
    t_output output(param.output_dir);


    Array3D<double> array(1,1,1);
    array[0][0][0] = 1.0;
    array.resize(10,10,10);
    array.multiply(0.1);
    for(int i=0; i<10; i++)
        for(int j=0; j<10; j++)
            for(int k=0; k<10; k++)
                array[i][j][k] = 0;
    Field3D field(10,10,10,0.1,0.1,0.1);
    field.reset();
    field.accumulate(0.1,0.5, 0.5, 0.5);
    cout << field.interpolate(0.505, 0.505, 0.505) <<endl;

    //Geometry geometry(param);
    //Solver solver(geometry, param);
    ElMag3D elmag(param);
    elmag.solver.save("solver_numeric.umf");
    elmag.solve();
    elmag.u.print("field3d.dat");
    elmag.u.print("field3d.vtk", "vtk");
    elmag.voltage.print("voltage.dat");

    t_random rnd;

    Elon elon(param, rnd, elmag, NULL);
    elon.charge = -1.6e-19;
    elon.resize(1);
    int ii = elon.insert();
    elon.particles[ii].x = param.x_max*0.5;
    elon.particles[ii].y = param.y_max*0.5;
    elon.particles[ii].z = param.z_max*0.7;
    elon.particles[ii].vx = elon.veV(0.03);
    elon.particles[ii].vy = 0.0;
    elon.particles[ii].vz = 0.0;

    t_particle * pp = &(elon.particles[ii]);

    ofstream fw((param.output_dir + "/traj1.dat").c_str());
    ofstream fwv((param.output_dir + "/vel1.dat").c_str());
    for(int i=0; i<param.niter; i++)
    {
        fw << pp->x <<" "<< pp->y <<" "<< pp->z <<endl;
        fwv << pp->vx <<" "<< pp->vy <<" "<< pp->vz <<endl;
        elon.advance();
    }
    return 0;
}
