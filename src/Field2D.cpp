#include <cmath>
#include <vector>
#include "Field2D.hpp"

/*
 *********************** definitions for Field2D ********************************
 */
void Field2D::resize(int x_sampl, int y_sampl, double _dx, double _dy, double _xmin, double _ymin)
{
    dx = _dx, dy = _dy;
    idx = 1.0/dx, idy = 1.0/dy;
    xmin = _xmin, ymin = _ymin;
    Array2D<double>::resize(x_sampl, y_sampl);
}

void Field2D::resize(Field2D & f)
{
    resize(f.jmax, f.lmax, f.GetDx(), f.GetDy(), f.GetXMin(), f.GetYMin());
}

bool Field2D::hasnan()
{
    for(int j=0; j<jmax; j++)
        for(int l=0; l<lmax; l++)
            if(std::isnan(data[j][l]))
                return true;
    return false;
}

void Field2D::print( ostream & out , double factor)
{
    for(int j=0; j<jmax; j++)
    {
        for(int l=0; l<lmax; l++)
            out << j*dx+xmin <<"\t"<< l*dy+ymin <<"\t"<< data[j][l]*factor <<endl;
        out << endl;
    }
}

void Field2D::print( const char * filename , double factor)
{
    ofstream out(filename);
    print(out, factor);
}

void Field2D::load(const char * filename)
{
    vector<double> xvec, yvec, fvec;
    double tmp1, tmp2, tmp3;

    std::ifstream fr(filename);
    if(fr.fail())
    {
        throw runtime_error("Field2D::load(): failed opening file\n");
    }
    string line;
    istringstream s_line;

    //load the numbers from file
    while(fr.good())
    {
        getline(fr, line);
        s_line.clear();
        s_line.str(line);

        if( s_line >> tmp1 >> tmp2 >> tmp3)
        {
            xvec.push_back(tmp1);
            yvec.push_back(tmp2);
            fvec.push_back(tmp3);
        }
    }

    //find dx, xmin, xmax
    double dx=0, xmin=xvec[0], xmax=xvec[xvec.size()-1];
    unsigned int i=1;
    while(i<xvec.size() && xvec[i]-xvec[i-1] == 0.0)
        i++;
    dx = xvec[i]-xvec[i-1];
    if(dx<0)
    {
        dx = -dx;
        xmin = xvec[xvec.size()-1];
        xmax = xvec[0];
    }
    double xsampl_d = (xmax-xmin)/dx+1;
    unsigned int xsampl = double2int(xsampl_d, 1e-1);
    // XXX think about recalculating dx based on double2int conversion

    //find dy, ymin, ymax
    double dy=0, ymin=yvec[0], ymax=yvec[xvec.size()-1];
    i = 1;
    while(i<yvec.size() && yvec[i]-yvec[i-1] == 0.0)
        i++;
    dy = yvec[i]-yvec[i-1];
    if(dy<0)
    {
        dy = -dy;
        ymin = yvec[yvec.size()-1];
        ymax = yvec[0];
    }
    double ysampl_d = (ymax-ymin)/dy+1;
    unsigned int ysampl = double2int(ysampl_d, 1e-1);
#ifdef DEBUG
    cout << "xvec.size " << xvec.size() <<endl;
    cout << "xmin = " << xmin <<"  xmax = "<< xmax <<"  dx = "<<dx <<"  xsampl = "<<xsampl<<endl;
    cout << "ymin = " << ymin <<"  ymax = "<< ymax <<"  dy = "<<dy <<"  ysampl = "<<ysampl <<endl;
#endif
    if(xsampl*ysampl != xvec.size())
        throw std::runtime_error("Fields::load_magnetic_field() wrong size of input vector");

    resize(xsampl, ysampl, dx, dy, xmin, ymin);

    for(i=0; i<xsampl; i++)
        for(unsigned int j=0; j<ysampl; j++)
            data[i][j] = numeric_limits<double>::quiet_NaN();

    int xi, yi;
    for(i=0; i<xvec.size(); i++)
    {
        xi = double2int((xvec[i]-xmin)/dx, 1e-1);
        yi = double2int((yvec[i]-ymin)/dy, 1e-1);
        data[xi][yi] = fvec[i];
    }

    for(i=0; i<xsampl; i++)
        for(unsigned int j=0; j<ysampl; j++)
            if(isnan(data[i][j]))
                throw std::runtime_error("Fields::load_magnetic_field() garbage loaded");
}
