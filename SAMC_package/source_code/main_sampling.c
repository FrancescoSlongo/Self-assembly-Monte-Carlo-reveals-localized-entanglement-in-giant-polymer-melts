//////////////////////////////////////
// efornasa, March 2026             //
// Main function for SAMC sampling. //
//////////////////////////////////////

// main_sampling.c

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "functions.h"

const bool output = true; // output flag  



int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <input_parameters_filename>\n", argv[0]);
    return EXIT_FAILURE;
  }

  const char *input_parameters_filename = argv[1];

  // Construct parameters
  Parameters *param = generate_parameters(input_parameters_filename);
  if (param) printf("\nInput parameters:\ndimension = %zu\nlattice_side = %zu\nn_linear_chains = %zu\nn_dumps = %zu\nn_steps_per_dump = %zu\n\n", param->dimension, param->lattice_side, param->n_linear_chains, param->n_dumps, param->n_steps_per_dump);
  else {
    fprintf(stderr, "Failed to generate input parameters.\n\n");
    return EXIT_FAILURE;
  }

	// Generate nearest-neighbors list
	size_t **nearest_neighbor = generate_nearest_neighbors_list(param->dimension, param->lattice_side);
	if (nearest_neighbor) printf("Nearest-neighbors list generated.\n");
  else {
    fprintf(stderr, "Failed to generate nearest-neighbors list.\n\n");
    return EXIT_FAILURE;
  }

  // Initialize configuration
  Configuration *config = initialize_configuration(param, nearest_neighbor);
  size_t i;
  if (config) {
    printf("End-points:");
		for (i = 0; i < 2 * config->n_linear_chains; i++) printf(" %zu", config->endpoint[i]);
		printf("\n\n");
  } else {
    fprintf(stderr, "Failed to initialize configuration.\n\n");
    return EXIT_FAILURE;
  }

	// Initialize rand state
  long seed = -1; // set random generator seed
  long *idum = &seed;

  ///////////////////////////////////////////////////////////////
  // TRAJECTORY
  samc_trajectory(config, nearest_neighbor, param, idum, output);
  ///////////////////////////////////////////////////////////////

  // Free memory
	free_2d_ull(nearest_neighbor);
  free_configuration(config);
  free(param);
 
  return EXIT_SUCCESS;
}
