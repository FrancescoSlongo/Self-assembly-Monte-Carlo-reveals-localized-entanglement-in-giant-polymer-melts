//////////////////////////////
// efornasa, March 2026     //
// SAMC sampling functions. //
//////////////////////////////

// functions.h

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdint.h>
#include <stdbool.h>


typedef struct { size_t dimension, lattice_side, n_linear_chains, n_dumps, n_steps_per_dump; } Parameters;
typedef struct { size_t dimension, lattice_side, n_sites, n_linear_chains, *endpoint, **connected_site; } Configuration;


// SAMC functions
Parameters *generate_parameters(const char *input_parameters_filename);

size_t *get_coordinates_from_site_index(size_t dimension, size_t lattice_side, size_t site_index);
size_t get_site_index_from_coordinates(size_t dimension, size_t lattice_side, size_t *x);
size_t **generate_nearest_neighbors_list(size_t dimension, size_t lattice_side);

void build_snake_recursive(size_t dimension, size_t lattice_side, size_t depth, size_t *coordinate, size_t *path, size_t *idx, int forward);
Configuration *initialize_configuration(Parameters *param, size_t **nearest_neighbor);
Configuration *copy_configuration_from_config_file(const char *config_filename, Parameters *param);
void free_configuration(Configuration *config);

void SAMC_MOVE(Configuration *config, size_t **nearest_neighbor, long *idum);

void samc_trajectory(Configuration *config, size_t **nearest_neighbor, Parameters *param, long *idum, bool output);


// Random integer generator
size_t rand_int(size_t max_value, long *idum);


// Memory allocation functions
size_t **alloc_2d_ull(size_t rows, size_t cols);
void free_2d_ull(size_t **a);


// Chain coordinates extraction functions
void print_unwrapped_configuration(Configuration *config);
static int periodic_step(size_t a, size_t b, size_t lattice_side);


#endif
