
#ifndef MYMATH_H
#define MYMATH_H

#include <cmath>
#include <cassert>
using namespace std;

#ifdef AMD64
typedef unsigned int uint32;
typedef int int32;
#else
typedef unsigned long int uint32;
typedef long int int32;
#endif

inline void mymath_test()
{
    assert(sizeof(uint32) == 4);
    assert(sizeof(int32) == 4);
}

template<class T>
inline T sqr(T x){return x*x;};

template<class T>
inline T cube(T x){return x*x*x;};

template<class T>
inline T max(T x, T y){return x>y ? x : y ;};

union intdouble{
    long int i;
    double d;
};

// Test of approximate float equality described in
// http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
// Float is casted to int using union, because it is probably
// the only way compatible with g++ strict aliasing rules
// as described under -fstrict-aliasing in g++ manual
inline bool eq(double x, double y, long int ulps = 16)
{
    assert(sizeof(double) == sizeof(long int));
    if (x == y)
        return true;
    intdouble ux, uy;
    ux.d = x;
    uy.d = y;
    if (abs(ux.i-uy.i) <= ulps)
        return true;
    return false;
}

template<class T>
inline T clamp(T x, T xmin, T xmax)
{
    return x>xmax ? xmax : (x<xmin ? xmin : x);
}

//template<class T>
//inline T min(T x, T y){return x<y ? x : y ;};

template<class T>
inline T mod(T x, T y)
{
    if( x>=0.0 && x<=y ) return x;
    return x - y*int(x/y) + (x<0 ? y : 0);
}
inline double norm(double x, double y, double z){return sqrt(sqr(x)+sqr(y)+sqr(z));}
inline double norm(double x, double y){return sqrt(sqr(x)+sqr(y));}
inline void normalize(double &x, double &y, double &z, double norm = 1.0)
{
    double sqnorm = sqr(x) + sqr(y) + sqr(z);
    if(!eq(sqnorm, 1.0))
    {
        sqnorm = sqrt(sqnorm);
        x /= sqnorm;
        y /= sqnorm;
        z /= sqnorm;
    }
}

#endif
