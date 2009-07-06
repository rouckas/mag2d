#ifndef UTIL_H
#define UTIL_H
#include <sstream>
using namespace std;
string double2string(double x)
{
    ostringstream o;
    o << x;
    return o.str();
};

string int2string(int x)
{
    ostringstream o;
    o << x;
    return o.str();
};
#endif
