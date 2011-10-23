#ifndef SPECLIST_H
#define SPECLIST_H

#define NTYPES 12
enum species_type { NONE, ELECTRON, ARGON, ARGON_POS, ARGON_META, O2, O2_POS, O2_NEG ,O_NEG, HYDROGEN, H_NEG, HELIUM};
extern const char * species_names[NTYPES];
extern const bool is_particle[NTYPES];

#endif
