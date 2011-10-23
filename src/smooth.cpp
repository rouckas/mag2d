#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <vector>
#include "util.cpp"
using namespace std;

int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout << "program smooth, usage: cat input.txt | smooth winsize ncols \
n2 n3 ... ni > output.txt\n   winsize - size of averaging window\n   ncols - \
total number of columns\n   n2 n3 ... ni - columns that should be averaged\n";
        exit(1);
    }

    int winsize = string2double(argv[1]);
    unsigned int ncols = string2double(argv[2]);
    vector<unsigned int> avgcols;
    for(int i=3; i< argc; i++)
    {
        avgcols.push_back(string2double(argv[i]));
        if(avgcols[i-3] >= ncols)
            throw runtime_error("averaging column number greater than total number of columns\n");
    }

    vector<double> columns(ncols);
    vector<double> sums(avgcols.size());



    istream fr(cin.rdbuf());
    ostream fw(cout.rdbuf());
    if(fr.fail())
    {
        throw runtime_error("smoother main(): failed opening file\n");
    }
    string line;
    istringstream s_line;

    //load the numbers from file
    int lineno = 1;
    while(fr.good())
    {
        getline(fr, line);
        s_line.clear();
        s_line.str(line);

        unsigned int i=0;
        while( i < ncols && s_line >> columns[i] ) i++;
        if(i < ncols) continue;
        for(unsigned int j=0; j<avgcols.size();  j++)
            sums[j] += columns[avgcols[j]];

        if(lineno % winsize == 0)
        {
            for(unsigned int j=0; j<avgcols.size();  j++)
            {
                columns[avgcols[j]] = sums[j]/winsize;
                sums[j] = 0;
            }
            for(unsigned int i=0; i<ncols; i++)
                fw << columns[i] << "\t";
            fw << endl;
        }
        lineno++;


    }

    return 0;
}
