#ifndef UTIL_H
#define UTIL_H
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
