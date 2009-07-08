#include "speclist.hpp"

extern const char * species_names[NTYPES] = {"NONE", "ELECTRON", "ARGON", "ARGON_POS", "ARGON_META", "O2", "O2_POS", "O2_NEG", "O_NEG", "HYDROGEN", "H_NEG", "HELIUM"};
const bool is_particle[NTYPES] = {false, true, false, true, false, false, true, true , true, false, true, false};
