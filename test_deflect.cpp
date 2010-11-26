#define GNUPLOT
#include "random.cpp"
#include <iostream>
#include <fstream>
using namespace std;


int main(int argc, char *argv[])
{
    t_random rnd;

    double x = 1.0, y = 0.0, z = 0.0;
    for(int i=0; i<1000; i++)
    {
        x = 1.0, y = 0.0, z = 0.0;
        rnd.deflect(0.1, x, y, z);
        cout << x <<" "<< y <<" "<< z <<endl;

        x = 1.0, y = 0.0, z = 0.0;
        rnd.deflect(3.141592/2, x, y, z);
        cout << x <<" "<< y <<" "<< z <<endl;

        x = 1.0, y = 0.0, z = 0.0;
        rnd.deflect(3.141592*0.75, x, y, z);
        cout << x <<" "<< y <<" "<< z <<endl;

        x = 0.707, y = 0.707, z = 0.0;
        rnd.deflect(3.141592*0.75, x, y, z);
        cout << x <<" "<< y <<" "<< z <<endl;

    }

    x = 1.0, y = 0.0, z = 0.0;
    for(int i=0; i<1000; i++)
    {
        rnd.deflect(0.03, x, y, z);
        cout << x <<" "<< y <<" "<< z <<endl;
    }

    return 0;
}
