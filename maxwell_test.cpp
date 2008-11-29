#include "random.h"
#include "gnuplot_i.h"
#include "histogram.hpp"
#include <cstdlib>
#include <cmath>

const int nsampl = 10000000;
const int nhist = 200;
const double k_B = 1.380662e-23;
template <class T>
inline T sqr(T x){return x*x;};

int main()
{
    double temperature = 300.0;
    double mass = 6.68173e-26;

    double v_max = sqrt(2.0*k_B*temperature/mass);

    srand(1);
    maxwell_flux_seed();
    t_histogram histogram(nhist, -3*v_max, 3*v_max);

    double x[nhist];
    double mxw[nhist];

    for(int i=0; i<nsampl; i++)
	histogram.add(normal2()*v_max/M_SQRT2);
	//histogram.add(maxwell_flux(v_max*5.0));

    for(int i=0; i<nhist; i++)
    {
	x[i] = histogram.position(i);
	mxw[i] = x[i]*mass/(k_B*temperature)*exp(-sqr(x[i]/v_max))*
	    (histogram.Max()-histogram.Min())/histogram.N_hist()*histogram.N_val();
	mxw[i] = sqrt(mass/(2*M_PI*k_B*temperature))*exp(-sqr(x[i]/v_max))*
	    (histogram.Max()-histogram.Min())/histogram.N_hist()*histogram.N_val();
    }
    gnuplot_ctrl *h1;
    h1 = gnuplot_init();
    gnuplot_cmd(h1, "set log y");
    gnuplot_plot_xy(h1, x, &(histogram[0]), histogram.N_hist(), "");
    gnuplot_setstyle(h1, "lines");
    gnuplot_plot_xy(h1, x, mxw, histogram.N_hist(), "");
    getchar();
    gnuplot_close(h1);
    return 0;
}
