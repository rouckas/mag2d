#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
using namespace std;
#define MAXCOLS 10

int main()
{
    double columns[MAXCOLS];
    int ncols = 2;
    double sums[MAXCOLS];
    int avgcols[MAXCOLS];
    int navgcols=2;
    avgcols[0] = 0;
    avgcols[1] = 1;
    int winsize = 10;

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

        int i=0;
        while( i < ncols && s_line >> columns[i] ) i++;
        if(i < ncols) continue;
        for(int j=0; j<navgcols;  j++)
            sums[j] += columns[avgcols[j]];

        if(lineno % winsize == 0)
        {
            for(int j=0; j<navgcols;  j++)
            {
                columns[avgcols[j]] = sums[j]/winsize;
                sums[j] = 0;
            }
            for(int i=0; i<ncols; i++)
                fw << columns[i] << "\t";
            fw << endl;
        }
        lineno++;


    }

    return 0;
}
