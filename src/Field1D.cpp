#include "Field1D.hpp"
using namespace std;

void Field1D::resize(int x_sampl, double dx, double _xmin)
{
    idx = 1.0/dx;
    xmin = _xmin;
    Array1D<double>::resize(x_sampl);
}

void Field1D::print( ostream & out, double factor)
{
    double dx=1.0/idx;
    for(int i=0; i<imax; i++)
    {
        out << i*dx <<"\t"<< data[i]*factor << endl;
    }
}

bool Field1D::hasnan()
{
    for(int i=0; i<imax; i++)
        if(isnan(data[i]))
            return true;
    return false;
}

void Field1D::print( const char * filename, double factor)
{
    ofstream out(filename);
    print(out, factor);
}

