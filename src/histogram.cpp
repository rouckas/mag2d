#include "histogram.hpp"
double & Histogram::operator[](int i){return data[i];};
double Histogram::position(int i){return min + (max-min)*(i+.5)/n_hist;};
double Histogram::Min(){return min;};
double Histogram::Max(){return max;};
int Histogram::N_hist(){return n_hist;};
Histogram::Histogram(int n, double hist_min, double hist_max) : sum(0), sum_tot(0)
{
    data = new double[n];
    n_hist = n;
    min = hist_min;
    max = hist_max;
    n_val_tot = n_val = 0;

    for(int i=0;i<n;i++)
	data[i] = 0;
}
Histogram::~Histogram()
{
    delete [] data;
}
int Histogram::add(double f)
{
    int j = -1;
    if(f < max && f > min)
    {
	j = (int)((f-min)*n_hist/(max-min));
	data[j]++;
	n_val++;
	sum += f;
    }
    n_val_tot++;
    sum_tot += f;
    return j;
}
int Histogram::add(double f,double weight)
{
    int j = -1;
    if(f < max && f > min)
    {
	j = (int)((f-min)*n_hist/(max-min));
	data[j] += weight;
	n_val += weight;
	sum += weight*f; // je to spravne?
    }
    n_val_tot += weight;
    sum_tot += weight*f; // je to spravne?
    return j;
}
double Histogram::N_val()
{
    return n_val;
}
void Histogram::reset()
{
    for(int i=0; i<n_hist; i++)
	data[i] = 0;
    n_val = 0;
    sum = 0;
    n_val_tot = 0;
    sum_tot = 0;
}

double Histogram::mean()
{
    return sum/n_val;
}
double Histogram::mean_tot()
{
    return sum_tot/n_val_tot;
}
double Histogram::norm()
{
    return N_val()*(Max()-Min())/N_hist();
}
void Histogram::print(const char * fname)
{
    std::ofstream fw(fname);
    print(fw);
}
void Histogram::print( std::ostream & out )
{
    for(int i=0; i<n_hist; i++)
	out << position(i) << "\t" << data[i]/norm() << std::endl;
}
