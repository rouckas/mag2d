#ifndef TABUL_H
#define TABUL_H
#include <stdexcept>
class callable
{
    public:
	virtual double operator()(double) = 0;
};

class tabulate 
{
    protected:
	double *data;
	double min,max;
	int samples;
	double dx;
	double idx;
    public:
	tabulate(){}
	tabulate(double _min, double _max, int _samples, double (*func)(double)) : min(_min), max(_max), samples(_samples)
	{
	    dx = (max-min)/(samples-1);
	    idx = 1.0/dx;
	    data = new double[samples];
	    for(int i=0; i<samples; i++)
		data[i] = func(i*dx+min);
	}
	tabulate(double _min, double _max, int _samples) : min(_min), max(_max), samples(_samples)
	{
	    idx = 1.0/((max-min)/(samples-1));
	    data = new double[samples];
	}
	~tabulate()
	{
	    delete [] data;
	}
	inline double operator()(double x)
	{
	    if (x>=max) return data[samples-1];
	    if (x<=min) return data[0];

	    int index = int(idx*(x-min));
	    double w = (x-index*dx) * idx;
	    return data[index]*(1-w) + data[index+1]*w;
	    //return data[x>max ? samples-1 : x<min ? 0 : int(idx*(x-min))];
	    //very quick & very dirty
	    // "stair" interpolation :-)
	}
};
class interpolate : public callable
{
    protected:
	int samples;
	double *xdata;
	double *ydata;
	double xmin,xmax;
    public:
	interpolate(){}
	interpolate(int _samples, double *_xdata, double *_ydata) : samples(_samples), xdata(_xdata), ydata(_ydata)
	{
	    xmin = xdata[0];
	    xmax = xdata[samples-1];
	}
	~interpolate()
	{
	}
	double operator()(double x)
	{
	    if (x>=xmax) return ydata[samples-1];
	    if (x<=xmin) return ydata[0];

	    int j1=0;
	    int j2=samples-1;
	    int j3;
	    while((j2-j1)>1)
	    {
		j3=(j1+j2)/2;
		if(x < xdata[j3]) j2=j3;
		else j1=j3;
	    }
	    double w = (x-xdata[j1]) / (xdata[j2]-xdata[j1]);
	    return ydata[j1]*(1-w) + ydata[j2]*w;
	}
};
//class vec_interpolate : public callable
class vec_interpolate
{
    protected:
	int samples;
	vector<double> xdata;
	vector<double> ydata;
	double xmin,xmax;
    public:
	vec_interpolate(vector<double> & _xdata, vector<double> & _ydata) : xdata(_xdata), ydata(_ydata)
	{
	    samples = xdata.size();
	    if(samples < 1) throw runtime_error("vec_interpolate: error: need at least 1 sample\n");
	    xmin = xdata[0];
	    xmax = xdata[samples-1];
	}
	vec_interpolate(double (*func)(double), double _min, double _max, int _samples) : samples(_samples), xmin(_min), xmax(_max)
	{
	    if(samples < 1) throw runtime_error("vec_interpolate: error: need at least 1 sample\n");
	    double dx = (xmax-xmin)/(samples-1);
	    xdata.reserve(samples);
	    ydata.reserve(samples);
	    for(int i=0; i<samples; i++)
	    {
		ydata.push_back(func(i*dx+xmin));
		xdata.push_back(i*dx+xmin);
	    }
	    for(int i=0; i<samples; i++)
	;//	cout << xdata[i] <<' '<<ydata[i]<<endl;
	}
	~vec_interpolate()
	{
	}
	inline double operator()(double x)
	{
	    if (x>=xmax) return ydata[samples-1];
	    if (x<=xmin) return ydata[0];

	    int j1=0;
	    int j2=samples-1;
	    int j3;
	    while((j2-j1)>1)
	    {
		j3=(j1+j2)/2;
		if(x < xdata[j3]) j2=j3;
		else j1=j3;
	    }
	    double w = (x-xdata[j1]) / (xdata[j2]-xdata[j1]);
	    return ydata[j1]*(1-w) + ydata[j2]*w;
	}
};
class tabulate_periodic : public tabulate
{
    public:
	tabulate_periodic(double _min, double _max, int _samples, double (*func)(double)) : tabulate(_min,_max,_samples,func)
	{};
	double operator()(double x)
	{
	    //return data[ int(idx*(x-min)) % samples + (x>min ? 0 : samples)];
	    return data[ x>min && x<max ? int(idx*(x-min)) : int(idx*(x-min)) % samples + (x>min ? 0 : samples)];
	    //very quick && very dirty
	}
};
#endif
