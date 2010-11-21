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
        double temperature;
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

void config_parse(string fname, vector<SpeciesParams*> & species_params, vector<InteractionParams*> & interactions);
#endif
