//////////////////////////////////////
// efornasa, February 2026          //
// Main function for SAMC sampling. //
//////////////////////////////////////

// main.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "functions.h"



int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <input_parameters_filename> <input_lattice_filename>\n", argv[0]);
    return EXIT_FAILURE;
  }
  
  char *input_parameters_filename = argv[1];
  char *input_lattice_filename = argv[2];

  // Construct parameters
  Parameters *param = construct_parameters(input_parameters_filename);
  if (param) printf("\nInput parameters:\nlattice dimension = %d\nlattice size = %d\nn_linear_chains = %d\nn_dumps = %d\nn_steps_per_dump = %d\n\n", param->dimension, param->size, param->n_linear_chains, param->n_dumps, param->n_steps_per_dump);
  else {
    fprintf(stderr, "Failed to construct input parameters.\n\n");
    return EXIT_FAILURE;
  }
  
  // Construct lattice
  HC_Lattice *lattice = construct_lattice_from_input_file(input_lattice_filename);
  if (lattice) printf("Lattice:\ndimension = %d\nsize = %d\nn_sites = %d\nn_edges = %d\n\n", lattice->dimension, lattice->size, lattice->n_sites, lattice->n_edges);
  else {
    fprintf(stderr, "Failed to construct lattice.\n\n");
    return EXIT_FAILURE;
  }

  // Initialize configuration
  Configuration *config = initialize_configuration(lattice, param->n_linear_chains);
  int i;
  if (config) {
    // Print explicitly the configuration in terms of edge occupation: spin[i]=0(1) if edge i is unoccupied(occupied) by a bond
    //printf("Starting configuration:\n");
    //for (i = 0; i < config->n_spins; i++) printf("%c", config->spin[i]);
		//printf("\n\n");
    printf("End-points:");
		for (i = 0; i < 2 * config->n_linear_chains; i++) printf(" %d", config->end[i]);
		printf("\n\n");
  } else {
    fprintf(stderr, "Failed to initialize configuration.\n\n");
    return EXIT_FAILURE;
  }

	// Initialize rand state
  long seed = -1; // set random generator seed
  long *idum = &seed;

  //////////////////////////////////////////////////////////////////////////
  // TRAJECTORY
  EVOLUTION(config, lattice, param->n_dumps, param->n_steps_per_dump, idum);
  //////////////////////////////////////////////////////////////////////////

  // Free memory
  free(param);
  free_configuration(config);
  free_lattice(lattice);
 
  return EXIT_SUCCESS;
}
