#ifndef PARAM_H
#define PARAM_H
#include <string>
#include <GetPot>

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
        enum Geometry { EMPTY, RF_22PT, RF_8PT, RF_HAITRAP, RF_QUAD, MAC,
            PENNING, PENNING_SIMPLE, TUBE };
        double x_max, y_max, z_max;
        double x_min, y_min, z_min;
        int x_sampl, y_sampl, z_sampl;

        double extern_field;
        std::string electric_field_static_file;
        std::string electric_field_rf_file;
        bool electric_field_from_file;

        std::string magnetic_field_file;
        bool magnetic_field_const;
        double Br, Bz, Bt;

        bool has_probe;
        double probe_radius;
        double probe_length;
        double u_probe;

        // the cell volume is calculated from these
        double n_particles_total;
        double density_total;

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
        double dt_elon;
        unsigned long int niter;
        t_advancer advancer;
        Coord coord;
        Boundary boundary;
        Geometry geometry;
        int src_fact;
        bool selfconsistent;
        bool u_smooth;
        bool rf;
        double rf_amplitude;
        double rf_U0;
        double rf_omega;
        bool particle_reload;
        unsigned long int t_print;
        unsigned long int t_print_dist;
        unsigned long int t_dist_sample;
        unsigned long int t_equilib;
        std::string particle_reload_dir;
        std::string output_dir;
        std::string species_conf_file;
        bool do_plot;
        double neutral_density;

        Param(GetPot & config);
};
#endif
