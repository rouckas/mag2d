#include <iostream>
#include <tr1/cmath>
using namespace std;

int main()
{
    for(double beta=1.0; beta<10; beta += 0.1)
    {
        // small deflection angle calculation according to
        // Nanbu, Kitatani, J. Phys. D., 28 (1995) 324

        double tmp = sqrt(beta*beta*beta*beta - 1.0);
        double xi0 = sqrt(beta*beta - tmp);
        double xi1 = sqrt(beta*beta + tmp);
        double zeta = xi0/xi1;
        double theta = std::tr1::comp_ellint_1(zeta)*M_SQRT2*beta/xi1;
        double chi = M_PI-2*theta;
        double chi0 = - 3*M_PI/(16.0*beta*beta*beta*beta);
        cout << beta <<" "<< chi <<" "<< chi0 <<endl;
    }

    return 0;
}
