#ifndef HISTOGRAM_H
#define HISTOGRAM_H 1
#include <fstream>
#include <iostream>
class Histogram
{
    double *data;
    int n_hist;
    double min,max;
    double n_val, n_val_tot;
    double sum, sum_tot;
    
    public:
    double & operator[](int i);
    double position(int i);
    int N_hist();
    double Min();
    double Max();
    double N_val();
    Histogram(int n, double hist_min, double hist_max);
    ~Histogram();
    int add(double f);
    int add(double f,double weight);
    void reset();
    double mean();
    double mean_tot();
    double norm();
    void print(const char * fname);
    void print(std::ostream & out = std::cout);
};
#endif
