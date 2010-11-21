#include <map>
#include <stdexcept>
#include "param.hpp"
#include "util.cpp"

namespace physconst
{
    const double eps_0 = 8.854187817e-12;
    const double k_B = 1.380662e-23;
    const double q_e = 1.602189e-19;
}
Param::Param(GetPot & config) : eps_0(physconst::eps_0), k_B(physconst::k_B), q_e(physconst::q_e)
{
    x_max = config("r_max",1e-2);   //backward compatibility
    x_max = config("x_max",x_max);
    y_max = config("y_max",1e-2);
    z_max = config("z_max",1e-2);
    x_sampl = config("r_sampl",100);
    x_sampl = config("x_sampl",x_sampl);
    y_sampl = config("y_sampl",100);
    z_sampl = config("z_sampl",100);

//	n_particles = config("n_particles",100000);
    n_particles_total = config("n_particles_total",100000);
//	rho = config("rho",1e15);		// XXX rho must correspond to n_particles
    pressure = config("pressure",133.0);
    neutral_temperature = config("neutral_temperature",300.0);

    probe_radius = config("probe_radius",1e-4);
    probe_length = config("probe_length",1e-2);
    u_probe = config("u_probe",-10.0);
    extern_field = config("extern_field",100.0);

    has_probe = config("has_probe", 0);

    magnetic_field_file = config("magnetic_field_file", "");
    magnetic_field_const = config("magnetic_field_const", 1);
    Br = config("Br", 0.0);
    Bz = config("Bz", 0.0);
    Bt = config("Bt", 0.0);

    string niter_str = config("niter","100000");
    niter = string2<unsigned long int>(niter_str);
    dt_elon = config("dt_elon",1e-11);

    selfconsistent = config("selfconsistent",1);
    rf = config("rf", 1);
    rf_amplitude = config("rf_amplitude", 10.0);
    rf_U0 = config("rf_U0", 0.0);
    rf_omega = config("rf_omega", 1e8);
    if(selfconsistent && rf)
    {
        printf("error: selfconsistent rf trap not implemented\n");
        exit(1);
    }
    string t_print_str = config("t_print","0");
    t_print = string2<unsigned long int>(t_print_str);
    string t_print_dist_str = config("t_print_dist","0");
    t_print_dist = string2<unsigned long int>(t_print_dist_str);
    t_dist_sample = t_print>10 ? t_print/10 : 1;
    string t_equilib_str = config("t_equilib","niter+1");
    if(t_equilib_str == "niter+1") t_equilib = niter+1;
    else t_equilib = string2<unsigned long int>(t_equilib_str);
    particle_reload = config("particle_reload",0);
    particle_reload_dir = config("particle_reload_dir",".");

    src_fact = config("src_fact",20);

    neutral_density = pressure/(k_B*neutral_temperature);
    do_plot = config("do_plot",1);

    string coord_str = config("coord","CYLINDRICAL");
    if(coord_str == "CYLINDRICAL") { coord = CYLINDRICAL; }
    else if(coord_str == "CARTESIAN") { coord = CARTESIAN; }
    else if(coord_str == "CARTESIAN3D") { coord = CARTESIAN3D; }
    else { throw std::runtime_error("Param: unrecognized coord value " + coord_str + "\n"); }

    string boundary_str = config("boundary","FREE");
    if(boundary_str == "FREE") { boundary = FREE; }
    else if(boundary_str == "MIRROR") { boundary = MIRROR; }
    else if(boundary_str == "PERIODIC") { boundary = PERIODIC; }
    else { throw std::runtime_error("Param: unrecognized boundary value " + boundary_str + "\n"); }
    if(coord == CYLINDRICAL && boundary != FREE)
        throw std::runtime_error("Param: only FREE boundary condition in cylindrical coords is implemented\n");
    if(boundary == MIRROR)
        throw std::runtime_error("Param: MIRROR boundary condition not implemented\n");

    string geometry_str = config("geometry","EMPTY");
    map<string, Geometry> string2geo;
    string2geo["EMPTY"] = EMPTY;
    string2geo["RF_22PT"] = RF_22PT;
    string2geo["RF_8PT"] = RF_8PT;
    string2geo["RF_QUAD"] = RF_QUAD;
    string2geo["MAC"] = MAC;
    string2geo["PENNING"] = PENNING;
    string2geo["PENNING_SIMPLE"] = PENNING_SIMPLE;
    if(string2geo.find(geometry_str) == string2geo.end())
        throw std::runtime_error("Param: unrecognized geometry value " + geometry_str + "\n");
    else
        geometry = string2geo[geometry_str];

    //parse plasma parameters
    double particle_density_total = 1.0;
    //XXX
    cerr << "particle_density_total is not calculated, this is wrong for selfconsistent sim." <<endl;
    /*
    for(int i=0; i<NTYPES; i++)
    {
        if(is_particle[i]) particle_density_total += density[i];
    }
    */

    dx = x_max/(x_sampl-1);
    dy = y_max/(y_sampl-1);
    dz = z_max/(z_sampl-1);
    idx = 1.0/dx;
    idy = 1.0/dy;
    idz = 1.0/dz;
    //V = n_particles/(rho);
    config.set_prefix("");
    macroparticle_factor = config("macroparticle_factor",1e4);
    cout << "macroparticle_factor " << macroparticle_factor <<endl;
    V = n_particles_total/particle_density_total;
    // XXX check this, only in cartesian coords:
    if(selfconsistent)
        dV = 1;
    else
    {
        if(coord == CYLINDRICAL)
            dV = V/((x_sampl-1)*(z_sampl-1));
        else
            dV = V/((x_sampl-1)*(y_sampl-1));
    }
    cout << "param: V = "<<V<<endl;
    if(coord == CYLINDRICAL)
        dy = 2*M_PI/macroparticle_factor;
}
