//#include "ion.cpp"
#include "util.cpp"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <stdexcept>

using namespace std;

enum SpeciesType {Neutral, Electron, Ion};
enum CollType {ELASTIC, CX};
class SpeciesParams
{
    public:
        string name;
        SpeciesType type;
        double charge;
        double mass;
        double dt;
        double density;
};

class InteractionParams
{
    public:
        string name;
        CollType type;
        double rate;
        string primary;
        string secondary;

};

void tokenize(string fname, vector< vector<string> > & tokens)
{
    ifstream fr(fname.c_str());
    string line;
    while(fr.good())
    {
        getline(fr, line);
        if(line[0] == '#') continue;
        istringstream s_line;
        vector<string> token_line(0);
        s_line.str(line);
        string token;
        while( s_line >> token )
        {
            token_line.push_back(token);
        }
        if(token_line.size()==0) continue;
        tokens.push_back(token_line);
    }

};

void config_parse(string fname, vector<SpeciesParams*> & species_params, vector<InteractionParams*> & interactions)
{
    vector< vector<string> > tokens;
    enum DecodingState{ DECODE_INIT, DECODE_DEFAULT, DECODE_SPECIES, DECODE_INTERACTION } state;

    tokenize(fname, tokens);
    for(size_t i=0; i< tokens.size(); i++)
    {
        vector<string> line = tokens[i];

        if(line.size() == 0) continue;

        if(line[0] == "DEFAULT")
        {
            state = DECODE_DEFAULT;
            continue;
        }
        if(line[0] == "SPECIES")
        {
            species_params.push_back(new SpeciesParams);
            //TODO specify defaults
            state = DECODE_SPECIES;
            continue;
        }
        if(line[0] == "INTERACTION")
        {
            interactions.push_back(new InteractionParams);
            state = DECODE_INTERACTION;
            continue;
        }

        if(line.size() > 0)
            if(state == DECODE_INIT)
                throw runtime_error("config_parse: unrecognized first config block\n");

        if(state == DECODE_SPECIES)
        {
            SpeciesParams * last = species_params[species_params.size()-1];
            if(line[0] == "NAME")
                last->name = line[1];
            else if(line[0] == "TYPE")
            {
                if(line[1] == "NEUTRAL")
                    last->type = Neutral;
                else if(line[1] == "ELECTRON")
                    last->type = Electron;
                else if(line[1] == "ION")
                    last->type = Ion;
                else
                    throw runtime_error("config_parse: unrecognized first species type \"" + line[1] + "\"");
            }
            else if(line[0] == "MASS")
                last->mass = string2<double>(line[1]);
            else if(line[0] == "CHARGE")
                last->charge = string2<double>(line[1]);
            else if(line[0] == "DENSITY")
                last->density = string2<double>(line[1]);
            else
                throw runtime_error("config_parse: unrecognized species  parameter\"" + line[0] + "\"");
        }

        if(state == DECODE_INTERACTION)
        {
            InteractionParams * last = interactions[interactions.size()-1];
            if(line[0] == "NAME")
                last->name = line[1];
            else if(line[0] == "TYPE")
            {
                if(line[1] == "ELASTIC")
                    last->type = ELASTIC;
                else if(line[1] == "CX")
                    last->type = CX;
                else
                    throw runtime_error("config_parse: unrecognized interaction type \"" + line[1] + "\"");
            }
            else if(line[0] == "RATE")
                last->rate = string2<double>(line[1]);
            else if(line[0] == "PRIMARY")
                last->primary = line[1];
            else if(line[0] == "SECONDARY")
                last->secondary = line[1];
            else
                throw runtime_error("config_parse: unrecognized species  parameter\"" + line[0] + "\"");
        }



    }
}

