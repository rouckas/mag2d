#include <cmath>
#include "../src/random.cpp"
using namespace std;

int main()
{
    t_random rnd;
    for(int k=0; k<10; k++)
    {
        int imax = 200000000;
        double j=0;
#pragma omp parallel for firstprivate(rnd) reduction(+:j)
        for(int i=0; i<imax; i++)
            //double j=sin(double(i)/imax);
            j += rnd.uni();
        cout << j<<endl;
    }

    return 0;
}
