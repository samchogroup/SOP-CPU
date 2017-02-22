
#ifndef IO_H
#define IO_H

void read_input(const char* const);
void record_traj(char*,char*);
void save_coords(char*,char*);
void load_coords(char*,char*);
void save_unccoords(char*);
void save_vels(char*);
void load_vels(char*);
void print_sim_params();

#endif
