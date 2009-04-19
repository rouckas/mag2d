#ifndef MATRIX_H
#define MATRIX_H

#include <exception>

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
            if(data != NULL)
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
	    delete [] data[0];
	    delete [] data;
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
#endif
