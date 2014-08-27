/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'main.c' is part of RASPA-2.0

 *************************************************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include "simulation.h"
#include "run.h"

int main(int argc, char **argv)
{
  int c;
  char *input = NULL, *input_crystal = NULL, *raspa_dir = NULL, *output = NULL;
  bool stream = false;

  // set default for the inputs
  input = strdup("simulation.input");
  input_crystal = strdup("");

  // set default RASPA_DIR
  raspa_dir = getenv("HOME");
  strcat(raspa_dir,"/RASPA/simulations");

  // get the raspa install directory from environement if defined
  if(getenv("RASPA_DIR")&&(strlen(getenv("RASPA_DIR"))>0))
    raspa_dir=getenv("RASPA_DIR");
  // allow multiple versions of raspa to coexist peacefully by using different dirs
  if(getenv("RASPA2_DIR")&&(strlen(getenv("RASPA2_DIR"))>0))
    raspa_dir=getenv("RASPA2_DIR");

  // parse command-line options (":" means the option has an argument)
  while((c=getopt(argc,argv,"a:vhsc:i:d:"))!=-1)
  {
    switch(c)
    {
      case 'a':
        strcpy(FileNameAppend,optarg);
        break;
      case 'h':
        printf("usage: simulate [-hv] [-ifile] [-ddir] [-s [-idata] [-cdata]]\n");
        printf("\t-h help\n");
        printf("\t-v version\n");
        printf("\t-s enables streaming (inputs must be stream, not filepath)\n");
        printf("\t-i the input (either file or stream)\n");
        printf("\t-c if streaming, the crystal structure (as a stream)\n");
        printf("\t-d the raspa directory\n");
        printf("\t-a appends the string to output-files\n");
        return 0;
      case 'i': // set the input-filename
        input = strdup(optarg);
        break;
      case 'd':  // set the raspa-directory
        raspa_dir = strdup(optarg);
        break;
      case 's': // Toggle the streaming boolean
        stream = true;
        break;
      case 'c': // If streaming, pass the input molecule here
        input_crystal = strdup(optarg);
        break;
      case 'v':
        fprintf(stderr, "RASPA 1.6-2 (2012)\n");
        return 0;
      default:
        return 1;
        break;
    }
  }
  output = run(input, input_crystal, raspa_dir, stream);

  // This prints the output, which can be piped into other applications.
  if (stream)
    printf("%s\n", output);
  free(output);

  return 0;
}
