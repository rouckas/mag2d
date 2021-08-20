#include "particles.hpp"
#include <tr1/cmath>


TrackedParticle::TrackedParticle(BaseSpecies *_species, unsigned int _index, Param &param) :
    species(_species), index(_index)
{
    string fname = param.output_dir + "/" + species->name
        + "_traj_" + int2string(index) + ".dat";
    fw.open(fname.c_str());
};
void TrackedParticle::print()
{
    // cannot store particle pointer directly, because it chanes during particles.resize()
    t_particle * pp = &(species->particles[index]);
    fw << setprecision(10) << pp->x <<' ' << pp->z <<' '<< pp->vx <<' '<< pp->vz <<' '<< pp->vy <<endl;
};

double Interaction::coulomb_sigma(double E)
{
    E *= param.q_e;
    double lambda_D = sqrt(param.eps_0*param.k_B*primary->temperature/(primary->density*primary->charge*primary->charge));
    double Lambda = primary->charge*secondary->charge/(4*M_PI*param.eps_0*E);
    double sigma = M_PI*Lambda*Lambda*log(lambda_D/Lambda);
    return sigma;
};


/*
 * *********************** definitions for BaseSpecies ************************
 */
void BaseSpecies::save(const string filename)
{
    ofstream fw(filename.c_str(),ios::out | ios::binary);
    if(!fw.is_open())
        throw runtime_error("Species::save(): failed opening file");

    int tmp;
    tmp = particles.size();
    fw.write((char*)&tmp,sizeof(int));

    int nparticles = particles.size() - empty.size();
    fw.write((char*)&nparticles,sizeof(int));
    tmp=0;
    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); I++)
    {
        if(I->empty==false)
        {
            if(++tmp > nparticles)
            {
                cerr << "t_particle::save(): data incosistent" << endl;
                break;
            }
            fw.write((char*)&(*I),sizeof(t_particle));
        }
    }
    //cerr << tmp << " " << nparticles << endl;

}

void BaseSpecies::load(const string filename)
{
    ifstream fr(filename.c_str(), ios::in | ios::binary);
    if(!fr.is_open())
        throw runtime_error("Species::load(): failed opening file");

    int tmp;
    fr.read((char*)&tmp,sizeof(int));
    particles.resize(tmp);

    int nparticles;
    fr.read((char*)&nparticles,sizeof(int));
    for(int i=0; i<nparticles; i++)
    {
        fr.read((char*)&particles[i], sizeof(t_particle));
        if(!fr.good()) cerr << "Species::load(): read error\n";
        // map particles back to working area if they left it during previous
        // nonselfconsistent simulation
        cerr << "warning: BaseSpecies::load() needs testing x_min != 0\n";
        if(particles[i].x < p_param->x_min || particles[i].x > p_param->x_max)
            particles[i].x = (p_param->x_max - p_param->x_min) * rnd->uni() + p_param->x_min;
        if(particles[i].z < p_param->z_min || particles[i].z > p_param->z_max)
            particles[i].z = (p_param->z_max - p_param->z_min) * rnd->uni() + p_param->z_min;
    }

    empty.resize(0);
    for(unsigned int i=nparticles; i<particles.size(); i++)
    {
        empty.push_back(i);
        particles[i].empty=true;
    }

}

void BaseSpecies::resize(unsigned int newsize)
{
    unsigned int oldsize = particles.size();
    if(newsize > oldsize)
    {
        particles.resize(newsize);
        for(unsigned int i = oldsize; i < newsize; i++)
        {
            particles[i].empty = true;
            empty.push_back(i);
        }
    }
}

void BaseSpecies::source5_save(const string filename)
{
    ofstream fw(filename.c_str(),ios::out | ios::binary);
    if(!fw.is_open())
        throw runtime_error("Species::source_save(): failed opening file");

    int nparticles = source2_particles.size() ;
    fw.write((char*)&nparticles,sizeof(int));

    for(vector<t_particle>::iterator I = source2_particles.begin(); I != source2_particles.end(); I++)
        fw.write((char*)&(*I),sizeof(t_particle));
}

void BaseSpecies::source5_load(const string filename)
{
    ifstream fr(filename.c_str(), ios::in | ios::binary);
    if(!fr.is_open())
    {
        cerr << " source5_load(): Warning: cannot open file " << filename <<endl;
        return;
    }

    int nparticles;
    fr.read((char*)&nparticles,sizeof(int));
    source2_particles.resize(nparticles);

    for(int i=0; i<nparticles; i++)
    {
        fr.read((char*)&source2_particles[i], sizeof(t_particle));
        if(!fr.good()) cerr << "Species::load(): read error\n";
    }
}

void BaseSpecies::lifetime_init()
{
    /* the energy cutoff should be rather small -
     * using large E_max will calculate lifetime correctly for high energy particles,
     * however, if dt is comparable to lifetime, short lifetime will cause underestimating
     * of collisions even for low energy particles. Therefore, decreasing E_max will improve
     * collisions of low energy particles. The influence on high energy particles might
     * be positive or negative (this is valid for roughly constant cross section)
     */
    double rate = 0;
    for(size_t i=0; i<rates_by_species.size(); i++)
    {
        rates_by_species[i] = svmax_find(interactions_by_species[i], veV(E_max)) * speclist[i]->density;
        rate += rates_by_species[i];
    }

    if(rate > 0.0)
        lifetime = 1.0/rate;
    else
        lifetime = INFINITY;

    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); I++)
        if( I->empty == false)
            I->time_to_death = rnd->rexp()*lifetime;

    for(vector<t_particle>::iterator I = source2_particles.begin(); I != source2_particles.end(); I++)
        if( I->empty == false)
            I->time_to_death = rnd->rexp()*lifetime;
}

void BaseSpecies::load_CS(const string & fname, vector<vec_interpolate*> & CS, vector<double> & Loss,
        const vector<string> & CSnames, const string & tag, int ncols)
{
    vector<double> Edata,CSdata;
    double Eloss;
    vec_interpolate *p_CS;

    for(vector<string>::const_iterator I=CSnames.begin(); I != CSnames.end(); I++)
    {
        ::load(fname,Edata,CSdata,*I,&Eloss,tag,ncols);
        p_CS = new vec_interpolate(Edata,CSdata);

        CS.push_back(p_CS);
        Loss.push_back(Eloss);
    }

}

double BaseSpecies::svmax_find(const vector<Interaction *> & interactions, double _vmax, int samples)
    // this routine assumes input CS in units 1e-16 cm^2, possible source of errors...
{
    double dv = _vmax/samples;
    double svmax=0.0;
    //cout << "svmax_find: vmax = "<< _vmax<<endl;
    for(double v=0; v<_vmax; v+=dv)
    {
        double sv=0;
        for(unsigned int i=0; i<interactions.size(); i++)
            sv += interactions[i]->sigma_v(v);
        if(isnan(sv))continue;
        if(sv>svmax) svmax=sv;
        //cout <<name <<' '<< v <<' '<< (*CS[0])(EeV(v)) + (*CS[1])(EeV(v)) <<' '<< sv << endl;
    }
    return svmax;
}

void BaseSpecies::scatter(t_particle &particle)
{
    //select interacting species
    double gamma = rnd->uni() / lifetime;
    double tmp = 0.0;
    size_t specid;
    // We do not iterate over the last item in rates_species array,
    // it should be select automatically if no other species was
    // chosen before
    for(specid=0; specid < rates_by_species.size()-1; specid++)
    {
        tmp += rates_by_species[specid];
        if(tmp > gamma)
            break;
    }


    // Generate the interacting particle
    // use a random sample from velocity distribution of "continuum" species,
    // otherwise select a random particle from the ensemble
    double vr2,vz2,vt2;
    t_particle * partner = NULL;
    if(speclist[specid]->particles.size() == 0)
        speclist[specid]->rndv(vr2,vz2,vt2);
    else
    {
        partner = speclist[specid]->random_particle();
        vr2 = partner->vx;
        vz2 = partner->vz;
        vt2 = partner->vy;
    }
    //cout << "    v1 = " << vr2 << " " << vz2 << " " << vt2 <<endl;

    double v_rel = norm(particle.vx-vr2, particle.vz-vz2, particle.vy-vt2);
    double m2 = speclist[specid]->mass;

    //select interaction
    gamma = rnd->uni() * rates_by_species[specid];
    tmp = 0.0;
    size_t intid;
    for(intid=0; intid < interactions_by_species[specid].size(); intid++)
    {
        tmp += interactions_by_species[specid][intid]->sigma_v(v_rel) * speclist[specid]->density;

        if(tmp > gamma)
            break;
    }


    // Null collision selected
    if(intid==interactions_by_species[specid].size()) return;

    Interaction * interaction = interactions_by_species[specid][intid];

    //cout << "interaction " << interaction->name << " of " <<
    //    interaction->primary->name << " with " << interaction->secondary->name <<endl;
    //cout << "interaction type " << interaction->type << " ELASTIC: " << ELASTIC <<endl;
    tmp = 1.0/(mass+m2);
    switch(interaction->type)
    {
        case SUPERELASTIC:
            {
                double E = interaction->E(v_rel) + interaction->DE;
                if(E<0)
                {
                    cout << "Negative E!!: Erel "<< interaction->E(v_rel) << "  DE " << interaction->DE <<
                        "   sigmav " << interaction->sigma_v(v_rel) << endl;
                }
                double v_rel = interaction->v_rel(E);
                // cout << E <<" "<< interaction->E(v_rel) << " " << "CRR"<<endl;

                double v_rel_x, v_rel_z, v_rel_y;
                rnd->rot(v_rel, v_rel_x, v_rel_z, v_rel_y);

                // v1_cm = v_rel * m2/(m1+m2)
                // v1 = v1_cm + v_cm
                particle.vx = (v_rel_x*m2 + particle.vx*mass + vr2*m2)*tmp;
                particle.vz = (v_rel_y*m2 + particle.vz*mass + vz2*m2)*tmp;
                particle.vy = (v_rel_z*m2 + particle.vy*mass + vt2*m2)*tmp;
            }
            break;

        case COULOMB:
        case ELASTIC:
            {
                //transform to center of mass system
                double v_cm_x = (particle.vx*mass + vr2*m2)*tmp;
                double v_cm_z = (particle.vz*mass + vz2*m2)*tmp;
                double v_cm_y = (particle.vy*mass + vt2*m2)*tmp;

                //isotropic random rotation
                double v_rel_x, v_rel_z, v_rel_y;
                rnd->rot(v_rel, v_rel_x, v_rel_y, v_rel_z);
                //cout << "    v1 = " << v1_cm_x << " " << v1_cm_y << " " << v1_cm_z ;

                //reverse transformation
                //v1 = v_rel*m2/(m1+m2) + v_cm;
                particle.vx = v_rel_x*m2*tmp + v_cm_x;
                particle.vz = v_rel_z*m2*tmp + v_cm_z;
                particle.vy = v_rel_y*m2*tmp + v_cm_y;
                //cout << "    v1 = " << particle.vx << " " << particle.vy << " " << particle.vz <<endl;
                if(interaction->type == COULOMB)
                {
                    partner->vx = -v_rel_x*mass*tmp + v_cm_x;
                    partner->vz = -v_rel_z*mass*tmp + v_cm_z;
                    partner->vy = -v_rel_y*mass*tmp + v_cm_y;
                }
            }
            break;

        case CX:
            {
                // in the resonant charge exchange, just swap velocities
                particle.vx = vr2;
                particle.vz = vz2;
                particle.vy = vt2;
            }
            break;

        case LANGEVIN:
            {
                //transform to center of mass system
                double v1_cm_x = (particle.vx - vr2)*m2*tmp;
                double v1_cm_y = (particle.vz - vz2)*m2*tmp;
                double v1_cm_z = (particle.vy - vt2)*m2*tmp;


                //generate random reduced impact parameter
                double beta = rnd->radius()*interaction->cutoff;

                if(beta>1.0)
                {
                    // small deflection angle calculation according to
                    // Nanbu, Kitatani, J. Phys. D., 28 (1995) 324
                    double tmp = sqrt(beta*beta*beta*beta - 1.0);
                    double xi0 = sqrt(beta*beta - tmp);
                    double xi1 = sqrt(beta*beta + tmp);
                    double zeta = xi0/xi1;
                    double theta = std::tr1::comp_ellint_1(zeta)*M_SQRT2*beta/xi1;
                    double chi = M_PI-2*theta;

                    rnd->deflect(chi, v1_cm_x, v1_cm_y, v1_cm_z);

                } else {
                    rnd->rot(v1_cm_x,v1_cm_y,v1_cm_z);
                }
                //reverse transformation
                //particle.vx = v1_cm_x + v_cm_x;
                particle.vx = v1_cm_x + (particle.vx*mass + vr2*m2)*tmp;
                particle.vz = v1_cm_y + (particle.vz*mass + vz2*m2)*tmp;
                particle.vy = v1_cm_z + (particle.vy*mass + vt2*m2)*tmp;
            }
            break;

        default:
            throw runtime_error("BaseSpecies::scatter: interaction " + interaction->name + " not implemented\n");
    }
}

void BaseSpecies::print_status( ostream & out)
{
    output << niter << " "
        << n_particles() << " "
        << energy_dist.mean_tot() << " " << t << endl;
        //<< probe_current_sum/nsampl*p_param->probe_length/p_param->dz << " "
        //<< probe_charge << " "
        //<< probe_energy_dist.mean_tot() << " " << probe_current/(nsampl*dt)<< endl;
}

void BaseSpecies::print_distribution()
{
    energy_dist.print( (p_param->output_dir + "/" + name + "_energy_dist.dat").c_str() );
    probe_energy_dist.print( (p_param->output_dir + "/" + name + "_probe_energy_dist.dat").c_str() );
    rhoAverage.print( (p_param->output_dir + "/" + name + "_rho.dat").c_str() , 1.0/nsampl);
    probe_angular_dist.print( (p_param->output_dir + "/" + name + "_probe_angular_dist.dat").c_str() );
    probe_angular_normalized_dist.print( (p_param->output_dir + "/" + name + "_probe_angular_normalized_dist.dat").c_str() );
}

void BaseSpecies::dist_sample()
{
    energy_dist_compute();
    source_energy_dist_compute();
    probe_current_sum += probe_current;
    rhoAverage.add(rho);
    nsampl++;
}

void BaseSpecies::dist_reset()
{
    energy_dist.reset();
    source_energy_dist.reset();
    probe_energy_dist.reset();
    probe_angular_dist.reset();
    probe_angular_normalized_dist.reset();
    rhoAverage.reset();
    nsampl = 0;
    probe_current_sum = 0;
    probe_current = 0;
}

void BaseSpecies::energy_dist_compute()
{
    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); I++)
        if(I->empty==false)
            if( energy_dist.add((SQR(I->vx)+SQR(I->vz)+SQR(I->vy))*mass*0.5/p_param->q_e) == -1 )
                ;//cout << I-particles.begin() << " " << I->vx << " " << I->vz << " " << I->vz << endl;
}

void BaseSpecies::source_energy_dist_compute()
{
    for(vector<t_particle>::iterator I = source2_particles.begin(); I != source2_particles.end(); I++)
        if(I->empty==false)
            source_energy_dist.add((SQR(I->vx)+SQR(I->vz)+SQR(I->vy))*mass*0.5/p_param->q_e);
                ;//cout << I-source2_particles.begin() << " " << I->vx << " " << I->vz << " " << I->vz << endl;
}


/*
 * ****************** definitions for Species<CYLINDRICAL> *******************
 */
void Species<CYLINDRICAL>::probe_collect(t_particle *I)
{
    if(!(I->z > 395e-3 && I->x < 45e-3)) return;
    //double center_r = p_param->x_max/2.0;
    //double center_z = p_param->z_max/2.0;
    probe_energy_dist.add((SQR(I->vx)+SQR(I->vz)+SQR(I->vy))*mass*0.5/p_param->q_e);
    //compute incidence angle
    //normal vector: n = (r-x0, z-y0)
    // n.v = |n|*|v|*cos(alpha)
    // (alpha) =-acos(n.v/(|n|*|v|))
    //double nx = center_r-I->x;
    //double ny = center_z-I->z;
    //double absV, alpha;
    //double absN = norm(nx,ny);
    //absV = norm(I->vx, I->vz);
    //alpha = acos( ( nx*I->vx + ny*I->vz )/(absN*absV) );
    //probe_angular_dist.add(abs(alpha));

    //absV = norm(I->vx,I->vz,I->vy);
    //alpha = acos( ( nx*I->vx + ny*I->vz )/(absN*absV) );
    //probe_angular_dist.add(abs(alpha));
    //probe_angular_normalized_dist.add(abs(alpha),1.0/sin(alpha));
    //cerr << (SQR(I->vx)+SQR(I->vz)+SQR(I->vz))*mass*0.5/p_param->q_e << " ";;
    //cerr << (I->vx)<<" "<<(I->vz)<<" "<<(I->vz) << endl;;

    probe_current += charge;
    probe_charge += charge;
}

void Species<CYLINDRICAL>::add_particle_beam_on_disk_cylindrical(int nparticles, double centerz, double radius)
{
    double x, y, r;
    if(centerz < 0 || centerz > p_param->z_max) return;

    for(int i=0; i<nparticles; i++)
    {
        //silly implementation coming from variant of this method in cartesian coords
        do{
            x = (rnd->uni() - 0.5);
            y = (rnd->uni() - 0.5);
            r = SQR(x) + SQR(y);
        }while(r > 0.25 );
        r = sqrt(r)*radius*2;
        if(r > p_param->x_max) continue;

        int ii = insert();
        t_particle * pp = &(particles[ii]);
        pp->x = r;
        pp->z = centerz;
        pp->vx = pp->vy = 0.0;
        rndv(pp->vz);
    }
}

void Species<CYLINDRICAL>::add_monoenergetic_particles_on_cylinder_cylindrical(int nparticles, double energy, double centerz, double radius, double height)
{
    double x, y, r;
    if(centerz < 0 || centerz > p_param->z_max) return;

    for(int i=0; i<nparticles; i++)
    {
        //silly implementation coming from variant of this method in cartesian coords
        do{
            x = (rnd->uni() - 0.5);
            y = (rnd->uni() - 0.5);
            r = SQR(x) + SQR(y);
        }while(r > 0.25 );
        r = sqrt(r)*radius*2;
        if(r > p_param->x_max) continue;

        int ii = insert();
        t_particle * pp = &(particles[ii]);
        pp->x = r;
        pp->z = centerz + height*(rnd->uni()-0.5);

        double v = veV(energy);
        rnd->rot(v, pp->vx, pp->vz, pp->vy);
    }
}

void Species<CYLINDRICAL>::advance()
{
    switch(p_param->mover)
    {
        case Param::ADVANCE_BORIS:
            Species<CYLINDRICAL>::advance_boris();
            break;
        default:
            throw runtime_error("Species<CYLINDRICAL>: advance method not implemented\n");
    }

}
void Species<CYLINDRICAL>::advance_boris()
{
    double fr, fz;      //force vector
    double Bx, By, Bz;  //magnetic field vector
    double qmdt = charge/mass*dt;       //auxilliary constant
    double prob = 1.0-exp(-dt/lifetime);

    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); ++I)
    {
        if(I->empty==true) continue;

        //compute field at (I->x, I->z)
        field->E(I->x, I->z, fr, fz);
        field->B(I->x, I->z, Bx, Bz, By);
        //XXX osetrit castice mimo prac oblast !!!
        //fr=fz=0;


        // advance velocities as in cartesian coords
        // use HARHA in magnetic field
        // half acceleration:
        I->vx += fr*qmdt/2.0;
        I->vz += fz*qmdt/2.0;

        // rotation
        //use Boris' algorithm (Birdsall & Langdon pp. 62) for arbitrary B direction
        double tmp = charge*dt/(2.0*mass);
        double tx = Bx*tmp;
        double ty = By*tmp;
        double tz = Bz*tmp;

        //XXX check orientation
        // use (r, theta, z)
        double vprime_r = I->vx + I->vy*tz - I->vz*ty;
        double vprime_t = I->vy + I->vz*tx - I->vx*tz;
        double vprime_z = I->vz + I->vx*ty - I->vy*tx;

        tmp = 2.0/(1+SQR(tx)+SQR(ty)+SQR(tz));
        double sx = tx*tmp;
        double sy = ty*tmp;
        double sz = tz*tmp;

        I->vx = I->vx + vprime_t*sz - vprime_z*sy;
        I->vy = I->vy + vprime_z*sx - vprime_r*sz;
        I->vz = I->vz + vprime_r*sy - vprime_t*sx;

        // half acceleration:
        I->vx += fr*qmdt/2.0;
        I->vz += fz*qmdt/2.0;

        //advance position (Birdsall pp. 338):
        double x2 = I->x + I->vx*dt;
        double y2 = I->vy*dt;
        I->x = sqrt(SQR(x2)+SQR(y2));
        I->z += I->vz * dt;

        //rotate the speed vector
        double sa = y2/I->x;
        double ca = x2/I->x;
        if(I->x==0)
        {
            sa = 0;
            ca = 1;
        }
        tmp = I->vx;
        I->vx = ca*I->vx + sa*I->vy;
        I->vy = -sa*tmp + ca*I->vy;

        if( rnd->uni() < prob)
        {
            scatter(*I);
        }

        // OKRAJOVE PODMINKY
        if(I->x > p_param->x_max || I->x < 0
                || I->z > p_param->z_max || I->z < 0 )
        {
            remove(I);
            continue;
        }
        if(!field->grid.is_free(I->x, I->z))
        {
            if(p_param->has_probe)
                probe_collect(&*I);

            remove(I);
            continue;
        }

        // SUMACE NABOJE
        rho.accumulate(charge, I->x, I->z);
    }
    niter++;
    t += dt;
}

void Species<CYLINDRICAL>::accumulate()
{
    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); ++I)
    {
        if(I->empty==true) continue;
            // SUMACE NABOJE
            rho.accumulate(charge, I->x, I->z);
    }
}


/*
 * ******************** Species<CARTESIAN> definitions *********************
 */
void Species<CARTESIAN>::add_particles_everywhere(int nparticles)
{
    for(int i=0; i<nparticles; i++)
    {
        int ii = insert();
        particles[ii].x = (p_param->x_max - p_param->x_min)*(rnd->uni()) + p_param->x_min;
        particles[ii].y = (p_param->y_max - p_param->y_min)*(rnd->uni()) + p_param->y_min;
        particles[ii].z = (p_param->z_max - p_param->z_min)*(rnd->uni()) + p_param->z_min;
        particles[ii].vx = rnd->rnor()*v_max/(M_SQRT2);
        particles[ii].vy = rnd->rnor()*v_max/(M_SQRT2);
        particles[ii].vz = rnd->rnor()*v_max/(M_SQRT2);
        particles[ii].time_to_death = rnd->rexp()*lifetime;
    }
}

void Species<CARTESIAN>::add_particles_bessel(int nparticles, double centerx, double centery, double radius)
{
    double x, y;
    // std::tr1::cyl_bessel_j(0, x)
    for(int i=0; i<nparticles; i++)
    {
        double r=0;
        const double bessel_root = 2.404825557695773;
        do{
            x = (rnd->uni()*2 - 1.0);
            y = (rnd->uni()*2 - 1.0);
            r = sqrt(SQR(x) + SQR(y));
        }while(r > 1.0 || rnd->uni() > std::tr1::cyl_bessel_j(0, r*bessel_root));
        x = x*radius + centerx;
        y = y*radius + centery;
        if(x < p_param->x_min || x > p_param->x_max ||
                y < p_param->z_min || y > p_param->z_max) continue;

        int ii = insert();
        t_particle * pp = &(particles[ii]);
        pp->x = x;
        pp->z = y;
        rndv(pp->vx, pp->vz, pp->vy);
    }
}

void Species<CARTESIAN>::add_particles_on_disk(int nparticles, double centerx, double centery, double radius)
{
    double x, y;
    for(int i=0; i<nparticles; i++)
    {
        do{
            x = (rnd->uni() - 0.5);
            y = (rnd->uni() - 0.5);
        }while(SQR(x) + SQR(y) > 0.25 );
        x = x*2*radius + centerx;
        y = y*2*radius + centery;
        if(x < p_param->x_min || x > p_param->x_max ||
                y < p_param->z_min || y > p_param->z_max) continue;

        int ii = insert();
        t_particle * pp = &(particles[ii]);
        pp->x = x;
        pp->z = y;
        rndv(pp->vx, pp->vz, pp->vy);
    }
}

void Species<CARTESIAN>::add_particle_beam_on_disk(int nparticles, double centerx, double centery, double radius)
{
    double x, y;
    for(int i=0; i<nparticles; i++)
    {
        do{
            x = (rnd->uni() - 0.5);
            y = (rnd->uni() - 0.5);
        }while(SQR(x) + SQR(y) > 0.25 );
        x = x*2*radius + centerx;
        y = y*2*radius + centery;
        if(x < p_param->x_min || x > p_param->x_max ||
                y < p_param->z_min || y > p_param->z_max) continue;

        int ii = insert();
        t_particle * pp = &(particles[ii]);
        pp->x = x;
        pp->z = y;
        pp->vx = pp->vz = 0.0;
        rndv(pp->vy);
    }
}

void Species<CARTESIAN>::advance()
{
    switch(p_param->mover)
    {
        case Param::ADVANCE_BORIS:
            advance_boris();
            break;
        case Param::ADVANCE_MULTICOLL:
            advance_multicoll();
            break;
        default:
            throw runtime_error("Species<CARTESIAN>: advance method not implemented\n");
    }

}
void Species<CARTESIAN>::advance_multicoll()
{
    double fr, fz;      //force vector
    double Bx, By, Bz;  //magnetic field vector

    field->B(0., 0., Bx, Bz, By);
    if(Bx != 0. || Bz != 0. || By != 0.)
        throw runtime_error("Species<CARTESIAN>::advance_multicoll() implemented for zero B only\n");

    if(p_param->selfconsistent || p_param->geometry != Param::EMPTY)
        throw runtime_error("Species<CARTESIAN>::advance_multicoll() implemented for const extern field only:\n\tset selfconsistent=0 and geometry=EMPTY\n");

    field->E(0., 0., fr, fz, 0.);
    const double qm = charge/mass;       //auxilliary constant

    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); ++I)
    {
        // (r, t, z) ~ (x, y, z)
        if(I->empty==true) continue;

        // advance velocities
        // Currently the method is only suitable for simulation of EEDF
        // in constant homogeneous external field. We don't care about positions
        double local_time = 0.;
        while(local_time + I->time_to_death < dt)
        {
            I->vx += fr*qm*I->time_to_death;
            I->vz += fz*qm*I->time_to_death;
            local_time += I->time_to_death;

            scatter(*I);
            I->time_to_death = lifetime*rnd->rexp();
        }

        I->vx += fr*qm*(dt-local_time);
        I->vz += fz*qm*(dt-local_time);
        I->time_to_death -= dt - local_time;
    }
    niter++;
    t += dt;
}

void Species<CARTESIAN>::advance_boris(vector<t_particle> & what, bool extern_fields)
{
    //force vector
    double fx = 0;
    double fz = p_param->extern_field;
    //
    //magnetic field vector
    double Bx = p_param->Br;
    double Bz = p_param->Bz;
    double By = p_param->Bt;

    const double qmdt = charge/mass*dt;       //auxilliary constant
    const double prob = 1.0-exp(-dt/lifetime);

    for(vector<t_particle>::iterator I = what.begin(); I != what.end(); ++I)
    {
        // (r, t, z) ~ (x, y, z)
        if(I->empty==true) continue;

        //compute field at (I->x, I->z)
        if(!extern_fields)
        {
            field->E(I->x, I->z, fx, fz, niter*dt);
            field->B(I->x, I->z, Bx, Bz, By);
        }


        // advance velocities as in cartesian coords
        // use HARHA in magnetic field
        // half acceleration:
        I->vx += fx*qmdt/2.0;
        I->vz += fz*qmdt/2.0;

        // rotation
        //use Boris' algorithm (Birdsall & Langdon pp. 62) for arbitrary B direction
        double tmp = charge*dt/(2.0*mass);
        double tx = Bx*tmp;
        double ty = By*tmp;
        double tz = Bz*tmp;

        //XXX check orientation
        // use (r, theta, z)
        // use (x, y, z) ~ (r, z, theta) !!!
        // sign change
        double vprime_r = I->vx - I->vy*tz + I->vz*ty;
        double vprime_z = I->vz - I->vx*ty + I->vy*tx;
        double vprime_t = I->vy - I->vz*tx + I->vx*tz;

        tmp = 2.0/(1+SQR(tx)+SQR(ty)+SQR(tz));
        double sx = tx*tmp;
        double sy = ty*tmp;
        double sz = tz*tmp;

        I->vx = I->vx - vprime_t*sz + vprime_z*sy;
        I->vy = I->vy - vprime_z*sx + vprime_r*sz;
        I->vz = I->vz - vprime_r*sy + vprime_t*sx;

        // half acceleration:
        I->vx += fx*qmdt/2.0;
        I->vz += fz*qmdt/2.0;

        //advance position (Birdsall pp. 338):
        I->x += I->vx*dt;
        I->z += I->vz*dt;


        if( rnd->uni() < prob)
        {
            scatter(*I);
        }
    }
}

void Species<CARTESIAN>::advance_boris()
{

    advance_boris(particles, false);

    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); ++I)
    {
        if(I->empty) continue;

        // OKRAJOVE PODMINKY
        if(I->x < p_param->x_min || I->x > p_param->x_max ||
                I->z < p_param->z_min || I->z > p_param->z_max)
        {
            if(p_param->boundary == Param::FREE)
            {
                remove(I);
                continue;
            }
            else if(p_param->boundary == Param::PERIODIC)
            {
                I->x = fmod(I->x, p_param->x_max);
                if(I->x < 0) I->x += p_param->x_max;
                I->z = fmod(I->z, p_param->z_max);
                if(I->z < 0) I->z += p_param->z_max;
            }

        }
        if(!p_param->electric_field_from_file)
        {
            if(!field->grid.is_free(I->x, I->z))
            {
                remove(I);
                continue;
            }
        }

        /**
          if(SQR(I->x-center_r)+SQR(I->z-center_z) < sqpp)
          {
          probe_collect(&*I);
          remove(I);
          probe_current += charge;
          continue;
          }
          */

        // SUMACE NABOJE
        rho.accumulate(charge, I->x, I->z);
        //}
    }
    niter++;
    t += dt;
}

void Species<CARTESIAN>::accumulate()
{
    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); ++I)
    {
        if(I->empty==true) continue;
            // SUMACE NABOJE
            rho.accumulate(charge, I->x, I->z);
    }
}

void Species<CARTESIAN>::source5_refresh(unsigned int factor)
{
    source5_factor = factor;
    //double N = p_param->rho*p_param->V;               //XXX wrong for multicomponent plasma
    //double N = p_param->density[type]*p_param->V;             //XXX wrong for multicomponent plasma
    double N = density*p_param->V;
    unsigned int n_particles = N/factor;
    double K = 1.0/factor;

    // the second comparison checks if we simulate this species as particles
    if(source2_particles.size() != n_particles && particles.size() != 0)
        source2_particles.resize(n_particles);

    for(unsigned int i=0; i<source2_particles.size(); i++)
    {
        source2_particles[i].empty = false;
        source2_particles[i].x = K*p_param->x_max*rnd->uni();
        source2_particles[i].z = K*p_param->z_max*rnd->uni();
        source2_particles[i].vx = rnd->rnor()*v_max/(M_SQRT2);
        source2_particles[i].vz = rnd->rnor()*v_max/(M_SQRT2);
        source2_particles[i].vy = rnd->rnor()*v_max/(M_SQRT2);
        source2_particles[i].time_to_death = rnd->rexp()*lifetime;
    }
}




/*
 * not working currently
 */
void Species<CARTESIAN>::source_old()
{

    //double v_max = sqrt(2.0*t_param::k_B*temperature/mass);
    //double f2 = p_param->rho*M_2_SQRTPI*0.25*v_max*dt*p_param->dz;
    double f2 = density/(2*sqrt(M_PI))*v_max*dt*p_param->dz;
    int pocet;
    int k;
    double f1;
    double u;
    vector<t_particle>::iterator I;

    for(k=0; k<4; k++)
    {
        if(k<2) //r==konst
            f1 = f2*p_param->z_max;
        else
            f1 = f2*p_param->x_max;
        //cerr << f1*4 << endl;

        //if(f1 > 10 )
        //{
            if((pocet = (int)( f1 + sqrt(f1)*rnd->rnor() + .5)) < 0)
                pocet = 0;
        //}else {
            pocet = (int)f1;
            pocet += (rnd->uni()>f1-pocet) ? 0 : 1;
        //}
        for(int i=0; i<pocet; i++)
        {
            I = particles.begin() + insert();

            if(k<2)
            {
                I->vx = rnd->rnor()*v_max/(M_SQRT2);
                I->vz = rnd->maxwell_flux(v_max*5.0) * (k==0 ? 1 : -1);
                //I->x = I->vx*dt*uniform() + p_param->x_max * (k==0 ? 0 : 1);
                //I->z = p_param->z_max*uniform();
                u = rnd->uni();
                I->x = -I->vx*dt*u + p_param->x_max * (k==0 ? 0 : 1);
                I->z = p_param->z_max*rnd->uni() - I->vz*dt*u;
            }else
            {
                I->vx = rnd->rnor()*v_max/(M_SQRT2);
                I->vz = rnd->maxwell_flux(v_max*5.0) * (k==2 ? 1 : -1);
                //I->z = I->vz*dt*uniform() + p_param->z_max * (k==2 ? 0 : 1);
                //I->x = p_param->x_max*uniform();
                u = rnd->uni();
                I->x = -I->vz*dt*u + p_param->z_max * (k==2 ? 0 : 1);
                I->z = p_param->x_max*rnd->uni() - I->vx*dt*u;
            }
            I->vz = rnd->rnor()*v_max/M_SQRT2;
            I->time_to_death = rnd->rexp()*lifetime;

            I->x += I->vx * dt;
            I->z += I->vz * dt;
            if(I->x < p_param->x_max && I->x > 0
                    && I->z < p_param->z_max && I->z > 0 )
                rho.accumulate(charge, I->x, I->z);
            /*
            for(int jj=0; jj<2; jj++)
            {
                I->time_to_death -= dt;
                if(I->time_to_death < 0)
                    I->time_to_death += exponential(lifetime);
            }
            */

        }
    }
}

void Species<CARTESIAN>::source()
{
    double K = 1.0/source5_factor;
    vector<t_particle>::iterator J;
    int j;
    //double qm2 = charge/mass*0.5;
    double qm = charge/mass;
    double src_z_max = K*p_param->z_max;
    double src_x_max = K*p_param->x_max;
    double fx=0;
    double fz=p_param->extern_field;
    double qmdt = qm*dt;

    //cout << "source\n";
    //cout << "diff: " << source2_particles.begin() - source2_particles.end() << endl;
    for(vector<t_particle>::iterator I=source2_particles.begin(); I != source2_particles.end(); I++)
    {

        I->vx += fx*qmdt;
        I->vz += fy*qmdt;
        I->x += I->vx * dt;
        I->z += I->vz * dt;

        if( I->time_to_death < dt)
        {
            scatter(*I);
            I->time_to_death += rnd->rexp()*lifetime;
        }
        I->time_to_death -= dt;

        if(I->x > src_x_max)
            while(I->x > src_x_max)
            {
                I->x -= src_x_max;
                j = insert();
                particles[j] = *I;
                particles[j].z += rand()%source5_factor*src_z_max;
                if(particles[j].z<p_param->z_max && particles[j].z>0 &&
                        particles[j].x<p_param->x_max && particles[j].x>0)
                    rho.accumulate(charge, particles[j].x, particles[j].z);
                else
                    remove(j);
            }
        else if(I->x < 0)
            while(I->x < 0)
            {
                j = insert();
                particles[j] = *I;
                particles[j].z += rand()%source5_factor*src_z_max;
                particles[j].x += p_param->x_max;
                I->x += src_x_max;
                if(particles[j].z<p_param->z_max && particles[j].z>0 &&
                        particles[j].x<p_param->x_max && particles[j].x>0)
                    rho.accumulate(charge, particles[j].x, particles[j].z);
                else
                    remove(j);
            }
        if(I->z > src_z_max)
            while(I->z > src_z_max)
            {
                I->z -= src_z_max;
                j = insert();
                particles[j] = *I;
                particles[j].x += rand()%source5_factor*src_x_max;
                if(particles[j].z<p_param->z_max && particles[j].z>0 &&
                        particles[j].x<p_param->x_max && particles[j].x>0)
                    rho.accumulate(charge, particles[j].x, particles[j].z);
                else
                    remove(j);
            }
        else if(I->z < 0)
            while(I->z < 0)
            {
                j = insert();
                particles[j] = *I;
                particles[j].x += rand()%source5_factor*src_x_max;
                particles[j].z += p_param->z_max;
                I->z += src_z_max;
                if(particles[j].z<p_param->z_max && particles[j].z>0 &&
                        particles[j].x<p_param->x_max && particles[j].x>0)
                    rho.accumulate(charge, particles[j].x, particles[j].z);
                else
                    remove(j);
            }
    }
}
