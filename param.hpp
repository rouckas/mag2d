#ifndef PARAM_H
#define PARAM_H
#include <string>
#include <GetPot>
#include "speclist.hpp"

using namespace std;
enum t_advancer { ADVANCE2, ADVANCE_OLD, ADVANCE_PROB, ADVANCE2_EULER };
enum Coord { CARTESIAN, CYLINDRICAL, CARTESIAN3D };

namespace physconst
{
    extern const double eps_0;
    extern const double k_B;
    extern const double q_e;
}


class Param
{
    public:
        enum Boundary { FREE, PERIODIC, MIRROR };
        enum Geometry { EMPTY, RF_22PT, RF_8PT, RF_QUAD, MAC, PENNING, PENNING_SIMPLE };
        double x_max, y_max, z_max;
        int x_sampl, y_sampl, z_sampl;
        double extern_field;
        bool has_probe;
        std::string magnetic_field_file;
        bool magnetic_field_const;
        double Br, Bz, Bt;
        double probe_radius;
        double probe_length;
        double u_probe;
        //    int n_particles;
        double n_particles_total;
        double dx, dy, dz;
        double V, dV;
        double idx, idy, idz;
        const double eps_0;
        const double k_B;
        const double q_e;
        //    double rho;                       //[m-3] electron density
        double pressure;            //[Pa]
        double neutral_temperature; //[K]
        double macroparticle_factor;
        double density[NTYPES];
        double dt[NTYPES];
        double temperature[NTYPES];
        double dt_elon;
        unsigned long int niter;
        t_advancer advancer;
        Coord coord;
        Boundary boundary;
        Geometry geometry;
        int src_fact;
        bool selfconsistent;
        bool rf;
        double rf_omega;
        bool particle_reload;
        unsigned long int t_print;
        unsigned long int t_print_dist;
        unsigned long int t_dist_sample;
        unsigned long int t_equilib;
        std::string particle_reload_dir;
        std::string output_dir;
        bool do_plot;
        double neutral_density;

        Param(GetPot & config);
};
#endif
