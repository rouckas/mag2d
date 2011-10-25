#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

class Field
{
        double value;
    public:
        Field(double _value) : value(_value) {};
        double get_value() { return value; };
        inline double calculate(double x) { return value/(x);};
};

class Particle
{
    public:
        double x, y, z;
        Particle() : x(0), y(0), z(0) {};
};

int main()
{
    double j = 0;
    int nparticles = 1000000;
    int niter = 100;

    {
        vector<Particle> particles(nparticles);
        for(int k=0; k<niter; k++)
        {
#pragma omp parallel for reduction(+:j)
            for(vector<Particle>::iterator I = particles.begin(); I < particles.end(); ++I)
            {
                I->z = I->y = sin(I->x+0.5);
                j += I->z + I->y + I->x;
            }
        }
    }
    

    {
        Particle * particles = new Particle[nparticles];
        for(int k=0; k<niter; k++)
        {
#pragma omp parallel for reduction(+:j)
            for(int i = 0; i < nparticles; i++)
            {
                particles[i].z = particles[i].y = sin(particles[i].x+0.5);
                j += particles[i].z + particles[i].y + particles[i].x;
            }
        }
    }
    
    cout << j << endl;

    return 0;
}
