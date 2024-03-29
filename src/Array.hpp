#ifndef MATRIX_H
#define MATRIX_H

#include <cstdlib>
#include <stdexcept>

template<class T>
class Array2D
{
    public:
	int jmax,lmax;
	T **data;
        Array2D() : jmax(0), lmax(0), data(NULL) {}
	Array2D(int _jmax, int _lmax) : jmax(_jmax), lmax(_lmax)
	{

	    data = new T* [jmax];
    
	    data[0] = new T [jmax*lmax];
    
	    for(int j=1; j<jmax; j++) 
		data[j] = data[j-1] + lmax;
	}
        Array2D(Array2D const &matrix) : Array2D(matrix.jmax, matrix.lmax)
        {
            this->assign(matrix);
        }
        void resize(int _jmax, int _lmax)
        {
            this->~Array2D();

            jmax = _jmax;
            lmax = _lmax;
	    data = new T* [jmax];
	    data[0] = new T [jmax*lmax];
	    for(int j=1; j<jmax; j++) 
		data[j] = data[j-1] + lmax;
        }
	~Array2D()
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
	void add(Array2D const &matrix)
	{
	    if(jmax != matrix.jmax || lmax != matrix.lmax)
		throw std::runtime_error("Array2D.add(): matrix sizes don't match\n");

	    for(int j=0; j<jmax; j++)
		for(int l=0; l<lmax; l++)
		    data[j][l] += matrix[j][l];
	}
        void assign(Array2D const &matrix)
	{
	    if(jmax != matrix.jmax || lmax != matrix.lmax)
		throw std::runtime_error("Array2D.assign(): matrix sizes don't match\n");

	    for(int j=0; j<jmax; j++)
		for(int l=0; l<lmax; l++)
		    data[j][l] = matrix[j][l];
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
        T *base;
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
            base = data[0][0];
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
            base = data[0][0];
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
        T & operator () (int i, int j, int k)
        {
            return base[(i*jmax + j)*kmax + k];
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

template<class T>
class Array1D
{
    public:
        int imax;
        T *data;
        Array1D() : imax(0), data(NULL) {}
        Array1D(int _imax) : imax(_imax)
        {
            data = new T [imax];
        }
        void resize(int _imax)
        {
            this->~Array1D();

            imax = _imax;
            data = new T [imax];
        }
        ~Array1D()
        {
            if(data != NULL)
            {
                delete [] data;
            }
        }
        T & operator [] (int i) const
        {
        }
        T & operator () (int i)
        {
            return data[i];
        }
        void add(Array1D const &array)
        {
            if(imax != array.imax)
                throw std::runtime_error("Array1D.add(): array sizes don't match\n");

            for(int i=0; i<imax; i++)
                data[i] += array[i];
        }
        void multiply(double x)
        {
            for(int i=0; i<imax; i++)
                data[i] *= x;
        }
        void reset()
        {
            for(int i=0; i<imax; i++)
                data[i] = 0;
        }
};
#endif
