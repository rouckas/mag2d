#include "random.h"
#include <iostream>
using namespace std;

int main()
{
    double x,y,z;
    x=1.0;
    y=z=0.0;
    for(int i=0; i<10000; i++)
    {
	rot(x,y,z);
	cout << x << " " << y << " " << z << endl;
    }
    return 0;
}
