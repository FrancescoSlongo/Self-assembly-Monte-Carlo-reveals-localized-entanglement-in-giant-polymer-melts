///////////////////////////////////////////
// efornasa, February 2026               //
// Header file with functions definition //
// for SAMC sampling.                    //
///////////////////////////////////////////

// functions.h

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdint.h>
#include <stdbool.h>


typedef struct { int dimension, size, n_linear_chains, n_dumps, n_steps_per_dump; } Parameters;

typedef struct {
  int dimension, size, n_sites, n_edges;
  int **site, **edge, **nearest_neighbor, **adjacent_edge, **incident_edge;
} HC_Lattice;

typedef struct {
  int n_spins, n_linear_chains, *end;
	char *spin;
} Configuration;


Parameters *construct_parameters(const char *input_parameters_filename);

HC_Lattice *construct_lattice(int dimension, int size);
HC_Lattice *construct_lattice_from_input_file(const char *input_lattice_filename);
void free_lattice(HC_Lattice *lattice);

Configuration *initialize_configuration(HC_Lattice *lattice, int n_linear_chains);
Configuration *initialize_configuration_from_input_file(const char *input_file, int n_spins, int n_linear_chains);
void free_configuration(Configuration *config);

void metropolis_move(Configuration *config, HC_Lattice *lattice, long *idum);

void EVOLUTION(Configuration *config, HC_Lattice *lattice, int n_dumps, int n_steps_per_dump, long *idum);


// Random generator functions
float ran2(long *idum);
int rand_int(int max_value, long *idum);


#endif
