#ifndef UTIL_H
#define UTIL_H
#include <sstream>
#include <cmath>
#include <stdexcept>
using namespace std;

inline string double2string(double x)
{
    ostringstream o;
    o << x;
    return o.str();
}

inline string int2string(int x)
{
    ostringstream o;
    o << x;
    return o.str();
}

inline int double2int(double x, double eps=1e-2)
{
    int res = (int)(x+0.5);
    if(fabs(res-x) > eps)
	throw std::runtime_error("double2int() " + double2string(x) + " is not integer\n");
    return res;
}

inline double string2double(string str)
{
    istringstream i(str);
    double x;
    i >> x;
    return x;
}
#endif
