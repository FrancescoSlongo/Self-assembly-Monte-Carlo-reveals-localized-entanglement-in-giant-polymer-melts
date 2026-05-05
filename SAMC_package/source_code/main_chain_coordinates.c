////////////////////////////////////
// efornasa, March 2026           //
// Extract chain coordinates from //
// SAMC configurations.           //
////////////////////////////////////

// main_chain_coordinates.c

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "functions.h"



int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <parameters_filename> <config_filename>\n", argv[0]);
    return EXIT_FAILURE;
  }
  
  const char *parameters_filename = argv[1];
  const char *config_filename = argv[2];

  // Construct parameters
  Parameters *param = generate_parameters(parameters_filename);
  if (!param) {
    fprintf(stderr, "Failed to generate input parameters.\n");
    return EXIT_FAILURE;
  }

	// Open configuration file
	FILE *file = fopen(config_filename, "r");
	if (!file) {
    fprintf(stderr, "Failed to open %s file.\n", config_filename);
		free(param);
    return EXIT_FAILURE;
  }

  // Copy configuration
  Configuration *config = copy_configuration_from_config_file(config_filename, param);
  size_t i;
  if (!config) {
    fprintf(stderr, "Failed to copy configuration from %s file.\n", config_filename);
    return EXIT_FAILURE;
  }

	// EXTRACT CHAIN COORDINATES
	print_unwrapped_configuration(config);

  // Free memory
  free_configuration(config);
  free(param);
 
  return EXIT_SUCCESS;
}
