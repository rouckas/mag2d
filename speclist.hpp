#ifndef SPECLIST_H
#define SPECLIST_H

#define NTYPES 9
enum species_type { NONE, ELECTRON, ARGON, ARGON_POS, ARGON_META, O2, O2_POS, O2_NEG ,O_NEG};
const char * species_names[NTYPES] = {"NONE", "ELECTRON", "ARGON", "ARGON_POS", "ARGON_META", "O2", "O2_POS", "O2_NEG", "O_NEG"};
const bool is_particle[NTYPES] = {false, true, false, true, false, false, true, true , true };

#endif
