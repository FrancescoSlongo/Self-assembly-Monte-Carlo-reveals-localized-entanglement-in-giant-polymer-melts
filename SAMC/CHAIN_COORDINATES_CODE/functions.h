////////////////////////////////////////////////////////////
// efornasa, February 2026                                //
// Header file with functions definition to extract chain //
// coordinates from configurations generated with SAMC.   //
////////////////////////////////////////////////////////////

// function.h

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdbool.h>


typedef struct {
  int dimension, size, n_sites, n_edges;
  int **site, **edge, **nearest_neighbor, **adjacent_edge, **incident_edge;
} HC_Lattice;

typedef struct {
  bool linear;
  int length;
  int *chain_array;
  int **chain_coordinate;
  double *Rcm;
	double Rg, Ree;
} Chain;

typedef struct {
  int n_linear_chains, n_rings;
  Chain **linear_chain, **ring;
} Configuration;


// Function prototypes
HC_Lattice *construct_lattice(int dimension, int size);
HC_Lattice *construct_lattice_from_input_file(const char *lattice_filename);
void free_lattice(HC_Lattice *lattice);

Chain *extract_linear_chain(HC_Lattice *lattice, char *edge_string, int endpoint);
Chain *extract_ring(HC_Lattice *lattice, char *edge_string);
int find_next_site(HC_Lattice *lattice, char *edge_string, int start_site, int *curr_bond);
void append_element_to_int_array(int **array, int *size, int element);
void free_chain(Chain *chain);

Configuration *construct_configuration_chains(HC_Lattice *lattice, char *edge_string, int *end, int n_linear_chains);
void free_configuration(Configuration *config);

#endif
