#ifndef MATRIX_H
#define MATRIX_H

#include <cstdlib>
#include <stdexcept>

template<class T>
class Matrix
{
    public:
	int jmax,lmax;
	T **data;
        Matrix() : jmax(0), lmax(0), data(NULL) {}
	Matrix(int _jmax, int _lmax) : jmax(_jmax), lmax(_lmax)
	{

	    data = new T* [jmax];
    
	    data[0] = new T [jmax*lmax];
    
	    for(int j=1; j<jmax; j++) 
		data[j] = data[j-1] + lmax;
	}
        void resize(int _jmax, int _lmax)
        {
            this->~Matrix();

            jmax = _jmax;
            lmax = _lmax;
	    data = new T* [jmax];
	    data[0] = new T [jmax*lmax];
	    for(int j=1; j<jmax; j++) 
		data[j] = data[j-1] + lmax;
        }
	~Matrix()
	{
            if(data != NULL)
            {
                delete [] data[0];
                delete [] data;
            }
	}
	T * & operator [] (int i) const
	{
	    return data[i];
	}
	void add(Matrix const &matrix)
	{
	    if(jmax != matrix.jmax || lmax != matrix.lmax)
		throw std::runtime_error("Matrix.add(): matrix sizes don't match\n");

	    for(int j=0; j<jmax; j++)
		for(int l=0; l<lmax; l++)
		    data[j][l] += matrix[j][l];
	}
	void multiply(double x)
	{
	    for(int j=0; j<jmax; j++)
		for(int l=0; l<lmax; l++)
		    data[j][l] *= x;
	}
	void reset()
	{
	    for(int j=0; j<jmax; j++)
		for(int l=0; l<lmax; l++)
		    data[j][l] = 0;
	}
};

template<class T>
class Array3D
{
    public:
        int imax, jmax,kmax;
        T ***data;
        Array3D() : imax(0), jmax(0), kmax(0), data(NULL) {}
        Array3D(int _imax, int _jmax, int _kmax) : imax(_imax), jmax(_jmax), kmax(_kmax)
        {

            data = new T** [imax];
            data[0] = new T* [imax*jmax];
            data[0][0] = new T [imax*jmax*kmax];

            for(int i=1; i<imax; i++)
                data[i] = data[i-1] + jmax;

            for(int i=0; i<imax; i++)
            {
                for(int j=0; j<jmax; j++)
                {
                    if(i==0 && j==0) continue;
                    data[i][j] = data[i][j-1] + kmax;
                }
            }
        }
        void resize(int _imax, int _jmax, int _kmax)
        {
            this->~Array3D();

            imax = _imax;
            jmax = _jmax;
            kmax = _kmax;
            data = new T** [imax];
            data[0] = new T* [imax*jmax];
            data[0][0] = new T [imax*jmax*kmax];

            for(int i=1; i<imax; i++)
                data[i] = data[i-1] + jmax;

            for(int i=0; i<imax; i++)
            {
                for(int j=0; j<jmax; j++)
                {
                    if(i==0 && j==0) continue;
                    data[i][j] = data[i][j-1] + kmax;
                }
            }
        }
        ~Array3D()
        {
            if(data != NULL)
            {
                delete [] data[0][0];
                delete [] data[0];
                delete [] data;
            }
        }
        T ** & operator [] (int i) const
        {
            return data[i];
        }
        void add(Array3D const &array)
        {
            if(imax != array.imax || jmax != array.jmax || kmax != array.kmax)
                throw std::runtime_error("Array3D.add(): array sizes don't match\n");

            for(int i=0; i<imax; i++)
                for(int j=0; j<jmax; j++)
                    for(int k=0; k<kmax; k++)
                        data[i][j][k] += array[i][j][k];
        }
        void multiply(double x)
        {
            for(int i=0; i<imax; i++)
                for(int j=0; j<jmax; j++)
                    for(int k=0; k<kmax; k++)
                        data[i][j][k] *= x;
        }
        void reset()
        {
            for(int i=0; i<imax; i++)
                for(int j=0; j<jmax; j++)
                    for(int k=0; k<kmax; k++)
                        data[i][j][k] = 0;
        }
};
#endif
