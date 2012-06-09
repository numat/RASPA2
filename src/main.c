/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'main.c' is part of RASPA.

    RASPA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RASPA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
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

int main(int argc, char **argv)
{
  int c,i;
  REAL energy,force_factor;
  char inputfilename[256],raspa_dir[256];

  // set default name for the input-file
  strcpy(inputfilename,"simulation.input");
  
  // set deault RASPA_DIR
  strcpy(raspa_dir,getenv("HOME"));
  strcat(raspa_dir,"/RASPA/simulations");
  RASPA_DIRECTORY=raspa_dir;
  
  // get the raspa install directory from environement if defined
  if(getenv("RASPA_DIR")&&(strlen(getenv("RASPA_DIR"))>0))
    RASPA_DIRECTORY=getenv("RASPA_DIR");

  // parse command-line options (":" means the option has an argument)
  while((c=getopt(argc,argv,"a:vhi:d:"))!=-1)
  {
    switch(c)
    {
      case 'a':
        strcpy(FileNameAppend,optarg);
        break;
      case 'h':
        printf("usage: simulate [-hv] [-ifile] [-ddir]\n");
        printf("\t-h help\n");
        printf("\t-v version\n");
        printf("\t-i the name of the input-file\n");
        printf("\t-d the raspa directory\n");
        printf("\t-a appends the string to output-files\n");
        return 0;
      case 'i': // set the input-filename
        strcpy(inputfilename,optarg);
        break;
      case 'd':  // set the raspa-directory
        strcpy(raspa_dir,optarg);
        RASPA_DIRECTORY=raspa_dir;
        break;
      default:
        return 1;
        break;
    }
  }

  // read input file 'simulation.res'
  ReadInputFile(inputfilename);

  // write the initial positions to files
  WriteFrameworkDefinitionCSSR("initial");
  WriteFrameworkDefinitionGulp("initial");
  WriteFrameworkDefinitionVASP("initial");
  WriteFrameworkDefinitionPDB("initial");
  WriteFrameworkDefinitionTinker("initial");
  WriteFrameworkDefinitionMOL("initial");
  WriteFrameworkDefinitionCIF("initial");

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
        PotentialGradient(ReturnPseudoAtomNumber("CH4_sp3"),ReturnPseudoAtomNumber("CH4_sp3"),SQR(i*CutOffVDW/1000),&energy,&force_factor);
        //TestFunction(i*CutOffChargeCharge/1000,&energy,&force_factor);
        printf("%g %g %g\n",i*CutOffVDW/1000,energy,force_factor);
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

  // write the final positions to files
  WriteFrameworkDefinitionCSSR("final");
  WriteFrameworkDefinitionGulp("final");
  WriteFrameworkDefinitionVASP("final");
  WriteFrameworkDefinitionPDB("final");
  WriteFrameworkDefinitionTinker("final");
  WriteFrameworkDefinitionMOL("final");
  WriteFrameworkDefinitionCIF("final");

  return 0;
}
