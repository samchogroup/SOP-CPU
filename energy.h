#ifndef ENERGY_H
#define ENERGY_H

void set_potential();
void energy_eval();
void stacking_energy();
void bond_energy();
void fene_energy();
void angular_energy();
void soft_sphere_angular_energy();
void vdw_energy();

void set_forces();
void clear_forces();
void force_eval();
void random_force();
void stacking_forces();
void bond_forces();
void fene_forces();
void angular_forces();
void soft_sphere_angular_forces();
void vdw_forces();

#endif /* ENERGY_H */
