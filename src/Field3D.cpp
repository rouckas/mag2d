#include "Field3D.hpp"
#include "util.cpp"
#include <limits>
using namespace std;

void Field3D::resize(int x_sampl, int y_sampl, int z_sampl, double dx, double dy, double dz,
                double _xmin, double _ymin, double _zmin)
{
    idx = 1.0/dx, idy = 1.0/dy, idz = 1.0/dz;
    xmin = _xmin, ymin = _ymin, zmin = _zmin;
    Array3D<double>::resize(x_sampl, y_sampl, z_sampl);
}

void Field3D::print( ostream & out, double factor)
{
    double dx=1.0/idx, dy=1.0/idy, dz=1.0/idz;
    for(int i=0; i<imax; i++)
    {
        for(int j=0; j<jmax; j++)
        {
            for(int k=0; k<kmax; k++)
                out << i*dx <<"\t"<< j*dy <<"\t"<< k*dz <<"\t"<< data[i][j][k]*factor << endl;
            out << endl;
        }
        out << endl;
    }
}

bool Field3D::hasnan()
{
    for(int i=0; i<imax*jmax*kmax; i++)
        if(isnan(base[i]))
            return true;
    return false;
}

void Field3D::print_vtk( ostream & out)
    //method to export the field in a legacy vtk format
{
    double dx=1.0/idx, dy=1.0/idy, dz=1.0/idz;

    // maybe the version could be even lower, but this is
    // based on http://www.vtk.org/VTK/img/file-formats.pdf
    // for vtk version 4.2
    out << "# vtk DataFile Version 2.0\n";
    out << "3D scalar field saved by plasma2d\n";
    out << "ASCII\n";

    //Dataset format
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << imax <<" "<< jmax <<" "<< kmax <<endl;
    out << "ORIGIN " << xmin <<" "<< ymin <<" "<< zmin <<endl;
    out << "SPACING " << dx <<" "<< dy <<" "<< dz <<endl;

    //Dataset attribute format
    out << "POINT_DATA " << imax*jmax*kmax <<endl;
    out << "SCALARS ScalarField double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int i=0; i<imax; i++)
        for(int j=0; j<jmax; j++)
            for(int k=0; k<kmax; k++)
                out << data[i][j][k] << endl;
}

void Field3D::print( const char * filename, string format, double factor)
{
    ofstream out(filename);
    if(format=="table")
        print(out, factor);
    else if(format == "vtk")
        print_vtk(out);
    else
	throw std::runtime_error("Field3D::print() unknown format " + format + "\n");
}

void get_vector_sampling( vector<double> & vec, double & dx, double & min, double & max, unsigned int & sampl)
{

    //find dx, rmin, rmax
    dx=0, min=vec[0], max=vec[vec.size()-1];
    unsigned int i=1;
    while(i<vec.size() && vec[i]-vec[i-1] == 0.0)
        i++;
    dx = vec[i]-vec[i-1];
    if(dx<0)
    {
        dx = -dx;
        min = vec[vec.size()-1];
        max = vec[0];
    }
    double sampl_d = (max-min)/dx+1;
    sampl = double2int(sampl_d);

}


void Field3D::load_table( const char * filename, int col)
{
    vector<double> xvec, yvec, zvec, Bvec;
    double tmp1, tmp2, tmp3, tmp4;

    std::ifstream fr(filename);
    if(fr.fail())
    {
        throw runtime_error("Fields::load_magnetic_field(): failed opening file\n");
    }
    string line;
    istringstream s_line;

    //load the numbers from file
    while(fr.good())
    {
        getline(fr, line);
        s_line.clear();
        s_line.str(line);

        if(!(s_line >> tmp1 >> tmp2 >> tmp3))
            continue;
        bool bad = false;
        for(int i=0; i<col+1; i++)
        {
            if(!(s_line >> tmp4))
                bad = true;
        }
        if(!bad)
        {
            xvec.push_back(tmp1);
            yvec.push_back(tmp2);
            zvec.push_back(tmp3);
            Bvec.push_back(tmp4);
        }
    }

    //find dx, rmin, rmax
    double dx, xmin, xmax;
    double dy, ymin, ymax;
    double dz, zmin, zmax;
    unsigned int xsampl, ysampl, zsampl;
    get_vector_sampling(xvec, dx, xmin, xmax, xsampl);
    get_vector_sampling(yvec, dy, ymin, ymax, ysampl);
    get_vector_sampling(zvec, dz, zmin, zmax, zsampl);

#ifdef DEBUG
    cout << "xvec.size " << xvec.size() <<endl;
    cout << "xmin = " << xmin <<"  xmax = "<< xmax <<"  dx = "<< dx <<"  xsampl = " << xsampl <<endl;
    cout << "ymin = " << ymin <<"  ymax = "<< ymax <<"  dy = "<< dy <<"  ysampl = " << ysampl <<endl;
    cout << "zmin = " << zmin <<"  zmax = "<< zmax <<"  dz = "<< dz <<"  zsampl = " << zsampl <<endl;
#endif
    if(xsampl*ysampl*zsampl != xvec.size())
        throw std::runtime_error("Field3D::load_table() wrong size of input vector");

    this->resize(xsampl, ysampl, zsampl, dx, dy, dz, xmin, ymin, zmin);

    unsigned int i, j, k;
    for(i=0; i<xsampl; i++)
        for(j=0; j<ysampl; j++)
            for(k=0; k<zsampl; k++)
            {
                (*this)(i, j, k) = numeric_limits<double>::quiet_NaN();
            }

    int xi, yi, zi;
    for(i=0; i<xvec.size(); i++)
    {
        xi = double2int((xvec[i]-xmin)/dx);
        yi = double2int((yvec[i]-ymin)/dy);
        zi = double2int((zvec[i]-zmin)/dz);
        (*this)(xi, yi, zi) = Bvec[i];
    }

    if(this->hasnan()) throw std::runtime_error("Field3D::load_table() garbage loaded");
}
