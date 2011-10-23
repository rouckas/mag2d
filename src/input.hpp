#ifndef INPUT_H
#define INPUT_H
#include <string>
#include <vector>
using namespace std;

void load( const string  & fname, vector<double> & xdata, vector<double> & ydata, const string & tag, double * Eloss, const string tag2="", int ncols=3);
/* Load the data from Phelps' cross section database and from
 * tabulated files in general. The data are loaded from file fname into
 * vectors xdata and ydata. The "tag2" denotes species and tag2 denotes
 * collisional process.
 */
#endif
