#include "species3d.hpp"

void Species<CARTESIAN3D>::advance()
{
    double Ex, Ey, Ez;
    double Bx, By, Bz;
    double qmdt = charge/mass*dt;       //auxilliary constant
    double prob = 1.0-exp(-dt/lifetime);

    for(vector<t_particle>::iterator I = particles.begin(); I != particles.end(); ++I)
    {
        // (r, t, z) ~ (x, y, z)
        if(I->empty==true) continue;

        fields.E(I->x, I->y, I->z, Ex, Ey, Ez);
        //cout <<niter <<" "<< I->x <<" "<< I->y <<" "<< I->z <<" ";
        //cout << Ex <<" "<< Ey <<" "<< Ez <<endl;
        //fields.B(I->x, I->y, I->z, Bx, By, Bz);
        Bx = By = Bz = 0;

        // advance velocities as in cartesian coords
        // use HARHA in magnetic field
        // half acceleration:
        I->vx += Ex*qmdt/2.0;
        I->vy += Ey*qmdt/2.0;
        I->vz += Ez*qmdt/2.0;

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
        double vprime_x = I->vx - I->vy*tz + I->vz*ty;
        double vprime_y = I->vy - I->vz*tx + I->vx*tz;
        double vprime_z = I->vz - I->vx*ty + I->vy*tx;

        tmp = 2.0/(1+SQR(tx)+SQR(ty)+SQR(tz));
        double sx = tx*tmp;
        double sy = ty*tmp;
        double sz = tz*tmp;

        I->vx = I->vx - vprime_y*sz + vprime_z*sy;
        I->vy = I->vy - vprime_z*sx + vprime_x*sz;
        I->vz = I->vz - vprime_x*sy + vprime_y*sx;

        // half acceleration:
        I->vx += Ex*qmdt/2.0;
        I->vy += Ey*qmdt/2.0;
        I->vz += Ez*qmdt/2.0;

        //advance position (Birdsall pp. 338):
        I->x += I->vx*dt;
        I->y += I->vy*dt;
        I->z += I->vz*dt;


        if( rnd->uni() < prob)
        {
            scatter(*I);
        }

        // OKRAJOVE PODMINKY
        if(I->x > p_param->x_max || I->x < 0
                || I->y > p_param->y_max || I->y < 0 
                || I->z > p_param->z_max || I->z < 0 )
        {
            if(p_param->boundary == Param::FREE)
            {
                remove(I);
                continue;
            }
            else if(p_param->boundary == Param::PERIODIC)
            {
                I->x = mod(I->x, p_param->x_max);
                I->y = mod(I->y, p_param->y_max);
                I->z = mod(I->z, p_param->z_max);
            }

        }
        if(!fields.geometry.is_free(I->x, I->y, I->z))
        {
            remove(I);
            continue;
        }
        fields.rho.accumulate(charge, I->x, I->y, I->z);
    }
    niter++;
}
