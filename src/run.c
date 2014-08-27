/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'run.c' is part of RASPA-2.0

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
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "utils.h"
#include "molecule.h"
#include "output.h"
#include "mc_moves.h"
#include "statistics.h"
#include "potentials.h"
#include "cbmc.h"
#include "monte_carlo.h"
#include "grids.h"
#include "ewald.h"
#include "input.h"
#include "recrossing.h"
#include "minimization_simulation.h"
#include "molecular_dynamics.h"
#include "spectra.h"
#include "phonon_disperion.h"
#include "spacegroup.h"
#include "inter_energy.h"
#include "numerical.h"
#include "matrix.h"
#include "vector.h"
#include "movies.h"
#include "status.h"
#include "framework.h"

extern bool STREAM = false;
extern char *INPUT_CRYSTAL = NULL;
extern char *PORE_SIZE_DISTRIBUTION_OUTPUT = NULL;
extern size_t PORE_SIZE_DISTRIBUTION_OUTPUT_SIZE = NULL;
extern char **FILE_CONTENTS = NULL;
extern size_t *FILE_SIZES = NULL;

/**
 * The core logic is separated from main to simplify wrapper functionality
 */
char* run(char *inputData, char *inputCrystal, char *raspaDir, bool stream)
{
  int i = 0;
  size_t chars = 0;
  REAL energy, force_factor;
  char *output = NULL, *temp = NULL, *delimiter = NULL;

  // There are a lot of globals kicking around.
  INPUT_CRYSTAL = strdup(inputCrystal);
  RASPA_DIRECTORY = strdup(raspaDir);
  STREAM = stream;

  if (STREAM)
    ReadInput(inputData);
  else
    ReadInputFile(inputData);

  // write the initial positions to files
  if (!STREAM)
  {
    WriteFrameworkDefinitionCSSR("initial");
    WriteFrameworkDefinitionGulp("initial");
    WriteFrameworkDefinitionVASP("initial");
    WriteFrameworkDefinitionPDB("initial");
    WriteFrameworkDefinitionTinker("initial");
    WriteFrameworkDefinitionMOL("initial");
    WriteFrameworkDefinitionCIF("initial");
    WriteSnapshotTinker("initial");

    if(CreateTinkerInput)
    {
      WriteFrameworkDefinitionTinker("tinker");
      WriteTinkerParameterFile();
      WriteTinkerKeyFile();
    }
  }

  // compute powder diffraction pattern
  if(ComputePowderDiffractionPattern)
    PowderDiffraction();

  // run a simulation-type specified in the input-file
  switch(SimulationType)
  {
    case NUMERICAL:
      CheckStatusNumerically();
      break;
    case MONTE_CARLO:
      MonteCarloSimulation();
      break;
    case MOLECULAR_DYNAMICS:
      MolecularDynamicsSimulation();
      break;
    case SPECTRA:
      VibrationalAnalysis();
      break;
    case PHONON_DISPERSION:
      PhononDisperionCurves();
      break;
    case MINIMIZATION:
      MinimalizationSimulation();
      //PhononDisperionCurves();
      break;
    case GLOBAL_MINIMIZATION:
      GlobalMinimumSimulation();
      break;
    case MAKE_GRID:
      OpenOutputFile();
      PrintPreSimulationStatus();
      Framework[0].FrameworkModel=FULL;
      CurrentSystem=0;
      MakeGrid();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case MAKE_ASCI_GRID:
      OpenOutputFile();
      PrintPreSimulationStatus();
      Framework[0].FrameworkModel=FULL;
      CurrentSystem=0;
      MakeASCIGrid();
      PrintPostSimulationStatus();
      CloseOutputFile();
    case VISUALIZATION:
      OpenOutputFile();
      PrintPreSimulationStatus();
      FreeEnergyProfile3D();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case PORE_SIZE_DISTRIBUTION:
      OpenOutputFile();
      PrintPreSimulationStatus();
      ComputePoreSizeDistribution();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case BARRIER_RECROSSING:
      OpenOutputFile();
      PrintPreSimulationStatus();
      BarrierRecrossing();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case TEST_SPACEGROUP:
      TestSpacegroup();
      break;
    case PLOT_POTENTIAL:
      OpenOutputFile();
      PrintPreSimulationStatus();
      CurrentSystem=0;
      for(i=0;i<=1000;i++)
      {
        PotentialGradient(ReturnPseudoAtomNumber("CH4_sp3"),ReturnPseudoAtomNumber("CH4_sp3"),SQR(i*CutOffVDW/1000),&energy,&force_factor,1.0);
        //TestFunction(i*CutOffChargeCharge/1000,&energy,&force_factor);
        fprintf(stderr, "%g %g %g\n",i*CutOffVDW/1000,energy,force_factor);
      }
      CloseOutputFile();
      break;
    case STATUS:
      OpenOutputFile();
      PrintPreSimulationStatus();
      CurrentSystem=0;
      Status();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case GENERATE_FRAMEWORK:
      OpenOutputFile();
      PrintPreSimulationStatus();
      CurrentSystem=0;
      GenerateFramework();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
    case USER:
      OpenOutputFile();
      PrintPreSimulationStatus();
      CurrentSystem=0;
      VibrationalAnalysis();
      PrintPostSimulationStatus();
      CloseOutputFile();
      break;
  }

  // If streaming, merge the simulation contents into one string
  if (STREAM)
  {
    // Just returning relevant data on PSD sim. Note: Should we return full
    // simulation data, or is the overhead not worth it?
    if (SimulationType == PORE_SIZE_DISTRIBUTION)
    {
      output = calloc(PORE_SIZE_DISTRIBUTION_OUTPUT_SIZE, sizeof(char));
      strcat(output, PORE_SIZE_DISTRIBUTION_OUTPUT);
      free(PORE_SIZE_DISTRIBUTION_OUTPUT);
    }
    // Returning full output on standard case.
    else
    {
      delimiter = "\n\nEND OF SIMULATION\n\n";

      for (i = 0; i < NumberOfSystems; i++)
          chars += FILE_SIZES[i];
      chars += strlen(delimiter) * (NumberOfSystems - 1) + 1;
      output = calloc(chars, sizeof(char));

      for (i = 0; i < NumberOfSystems-1; i++)
      {
        strcat(output, FILE_CONTENTS[i]);
        strcat(output, delimiter);
      }
      strcat(output, FILE_CONTENTS[NumberOfSystems - 1]);

      for (i = 0; i < NumberOfSystems-1; i++)
        free(FILE_CONTENTS[i]);
    }

  // Write the final positions to files
  } else {
    WriteFrameworkDefinitionCSSR("final");
    WriteFrameworkDefinitionGulp("final");
    WriteFrameworkDefinitionVASP("final");
    WriteFrameworkDefinitionPDB("final");
    WriteFrameworkDefinitionMOL("final");
    WriteFrameworkDefinitionCIF("final");
    WriteSnapshotTinker("final");
  }

  free(FILE_CONTENTS);
  free(FILE_SIZES);
  free(INPUT_CRYSTAL);

  return output;
}