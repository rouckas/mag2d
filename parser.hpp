#ifndef PARSER_H
#define PARSER_H

#include "util.cpp"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <stdexcept>

using namespace std;

enum SpeciesType {NEUTRAL, ELECTRON, ION};
enum CollType {ELASTIC, LANGEVIN, CX};
class SpeciesParams
{
    public:
        string name;
        SpeciesType type;
        double charge;
        double mass;
        double dt;
        double density;
        double temperature;
        double polarizability;
};

class InteractionParams
{
    public:
        string name;
        CollType type;
        double rate;
        double cutoff;
        string primary;
        string secondary;
        vector<double> CS_energy;
        vector<double> CS_value;

};

void config_parse(string fname, vector<SpeciesParams*> & species_params, vector<InteractionParams*> & interactions);
#endif
