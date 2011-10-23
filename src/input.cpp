#include <sstream>
#include <fstream>
#include "input.hpp"

void load( const string  & fname, vector<double> & xdata, vector<double> & ydata, const string & tag, double * Eloss, const string tag2, int ncols)
{
    ifstream fr(fname.c_str());
    string line;
    *Eloss = 0;

    // search for the species in file
    if (tag2 != "")
	while(fr.good())
	{
	    getline(fr,line);
	    if (line.find(tag2) != string::npos) break;
	}

    //search for the cross section
    while(fr.good())
    {
	getline(fr,line);
	if (line.find(tag) != string::npos) break;
    }

    //skip some useless data
    string data_tag("  ENERGY");
    while(fr.good())
    {
	getline(fr,line);
	// decode the energy loss optionally
	if (line.find("ENERGY LOSS =") != string::npos)
	{
	    istringstream s_line;
	    s_line.str(line);

	    string tmp;
	    s_line >> tmp >> tmp >> tmp >> *Eloss;
	}
	if (line.find(data_tag) != string::npos) break;
    }

    //load the data
    xdata.clear();
    ydata.clear();
    while(fr.good())
    {
	getline(fr,line);
	if(line=="") continue;

	istringstream s_line;
	s_line.str(line);
	double d1,d2,d3;
	for(int i=0; i<ncols-2; i++) if( !(s_line >> d1)) break;
	if( !( s_line >> d2 >> d3 ) ) break;

	xdata.push_back(d2);
	ydata.push_back(d3);
    }
}
