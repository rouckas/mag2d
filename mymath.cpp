
#ifndef MYMATH_H
#define MYMATH_H

#include <cmath>
using namespace std;

#ifdef AMD64
typedef unsigned int uint32;
typedef int int32;
#else
typedef unsigned long int uint32;
typedef long int uint32;
#endif


template<class T>
inline T sqr(T x){return x*x;};

template<class T>
inline T cube(T x){return x*x*x;};

template<class T>
inline T max(T x, T y){return x>y ? x : y ;};

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

#endif
