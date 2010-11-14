#include "parser.hpp"
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
using namespace std;

class Speclist;
class Species;

class Interaction
{
    public:
        string name;
        CollType type;
        double rate;
        Species * primary;
        Species * secondary;
        Interaction(InteractionParams * params) : name(params->name),
            type(params->type), rate(params->rate)
        {};

};

class Species
{
    public:
        string name;
        double mass;
        double charge;
        double density;
        SpeciesType type;
        double lifetime;
        Speclist * speclist;
        vector<Interaction *> interactions;

        Species(SpeciesParams * params) : name(params->name),
            mass(params->mass), charge(params->charge),
            density(params->density), type(params->type),
            lifetime(INFINITY)
        {};

        void lifetime_init()
        {
            double rate = 0;
            for(size_t i=0; i<interactions.size(); i++)
                rate += interactions[i]->rate * interactions[i]->secondary->density;
            if(rate > 0.0)
                lifetime = 1.0/rate;
            else
                lifetime = INFINITY;
        };
};

class Speclist
{
    public:
        Speclist(string configfile)
        {
            vector<SpeciesParams*> vs;
            vector<InteractionParams*> vi;
            config_parse(configfile, vs, vi);

            /*
             * create Species classes
             */
            for(size_t i=0; i<vs.size(); i++)
                data.push_back(new Species(vs[i]));

            /*
             * create Interaction classes
             */
            Interaction * pi;
            for(size_t i=0; i<vi.size(); i++)
            {
                pi = new Interaction(vi[i]);

                // find the primary interacting species
                size_t j=0;
                while(j<data.size() && vi[i]->primary != data[j]->name)
                    j++;
                if(j== data.size())
                    throw runtime_error("Speclist::Speclist: unrecognized primary species \"" + vi[i]->primary + "\" of interaction \"" + vi[i]->name + "\"\n");
                pi->primary = data[j];
                data[j]->interactions.push_back(pi);

                // find the secondary interacting species
                j=0;
                while(j<data.size() && vi[i]->secondary != data[j]->name)
                    j++;
                if(j== data.size())
                    throw runtime_error("Speclist::Speclist: unrecognized primary species \"" + vi[i]->secondary + "\" of interaction \"" + vi[i]->name + "\"\n");
                pi->secondary = data[j];
            }

            for(size_t i=0; i<data.size(); i++)
                data[i]->lifetime_init();

        }

        vector<Species*> data;
};


int main()
{
    Speclist speclist("species_conf.txt");
    /*
    vector<SpeciesParams*> vs;
    vector<InteractionParams*> vi;
    config_parse("species_conf.txt", vs, vi);
    cout << vs.size() <<endl;
    for(size_t i=0; i<vs.size(); i++)
    {
        cout << vs[i]->name << endl;
        cout << vs[i]->mass << endl;
        cout << vs[i]->charge << endl << endl;
    }
    cout <<endl;
    cout << vi.size() <<endl;
    for(size_t i=0; i<vi.size(); i++)
    {
        cout << vi[i]->name << endl;
        cout << vi[i]->primary << endl;
        cout << vi[i]->secondary << endl;
        cout << vi[i]->rate << endl << endl;
    }
    */
    return 0;
}
