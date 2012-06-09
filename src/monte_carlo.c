/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'monte_carlo.c' is part of RASPA.

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
#include <time.h>
#include "simulation.h"
#include "molecule.h"
#include "framework_energy.h"
#include "framework.h"
#include "utils.h"
#include "molecule.h"
#include "input.h"
#include "output.h"
#include "mc_moves.h"
#include "statistics.h"
#include "potentials.h"
#include "spacegroup.h"
#include "cbmc.h"
#include "grids.h"
#include "movies.h"
#include "sample.h"
#include "recrossing.h"
#include "monte_carlo.h"
#include "integration.h"
#include "thermo_baro_stats.h"
#include "equations_of_state.h"

/*********************************************************************************************************
 * Name       | MonteCarloSimulation                                                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Main routine to do a Monte-Carlo (MC) simulation.                                        *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void MonteCarloSimulation(void)
{
  int i,j;
  int NumberOfSystemMoves;
  int NumberOfParticleMoves;
  int NumberOfSteps;
  int ran_int;
  REAL ran;

  // for a crash-recovery we skip the initialization/equilibration and jump to the
  // position right after the crash-file was written in the production run
  if(ContinueAfterCrash)
  {
    if(SimulationStage==POSITION_INITIALIZATION) goto ContinueAfterCrashLabel1;
    else goto ContinueAfterCrashLabel2;
  }
  
  // allocate memory for sampling routines
  SampleRadialDistributionFunction(ALLOCATE);
  SampleNumberOfMoleculesHistogram(ALLOCATE);
  SamplePositionHistogram(ALLOCATE);
  SampleFreeEnergyProfile(ALLOCATE);
  SamplePoreSizeDistribution(ALLOCATE); 
  SampleEndToEndDistanceHistogram(ALLOCATE);
  SampleEnergyHistogram(ALLOCATE);
  SampleThermoDynamicsFactor(ALLOCATE);
  SampleFrameworkSpacingHistogram(ALLOCATE);
  SampleResidenceTimes(ALLOCATE);
  SampleDistanceHistogram(ALLOCATE);
  SampleBendAngleHistogram(ALLOCATE);
  SampleDihedralAngleHistogram(ALLOCATE);
  SampleAngleBetweenPlanesHistogram(ALLOCATE);
  SampleMoleculePropertyHistogram(ALLOCATE);

  SampleInfraRedSpectra(ALLOCATE);
  SampleMeanSquaredDisplacementOrderN(ALLOCATE);
  SampleVelocityAutoCorrelationFunctionOrderN(ALLOCATE);
  SampleRotationalVelocityAutoCorrelationFunctionOrderN(ALLOCATE);
  SampleMolecularOrientationAutoCorrelationFunctionOrderN(ALLOCATE);
  SampleMeanSquaredDisplacement(ALLOCATE);
  SampleVelocityAutoCorrelationFunction(ALLOCATE);

  SampleDensityProfile3DVTKGrid(ALLOCATE);
  SampleCationAndAdsorptionSites(ALLOCATE);
  SampleDcTSTConfigurationFiles(ALLOCATE);
  SamplePDBMovies(ALLOCATE,-1);

  // loop over all the pressures of the isotherm
  for(CurrentIsothermPressure=0;CurrentIsothermPressure<NumberOfIsothermPressures;CurrentIsothermPressure++)
  {
    // compute gas properties
    // here the pressure is converted to fugacities
    ComputeGasPropertiesForAllSystems();

    // open output-file for systems
    OpenOutputFile();

    // print simulation settings to the output-file
    PrintPreSimulationStatus();

    // initialize the energies and compute the total energies for all systems
    InitializesEnergiesAllSystems();
    InitializesEnergyAveragesAllSystems();

    InitializeSmallMCStatisticsAllSystems();
    InitializeMCMovesStatisticsAllSystems();

    // compute total energy for all systems
    CalculateTotalEnergyAllSystems();

    // initialization to reach equilibration of positions (no averages are computed yet)
    SimulationStage=POSITION_INITIALIZATION;
    for(CurrentCycle=0;CurrentCycle<NumberOfInitializationCycles;CurrentCycle++)
    {
      if((CurrentCycle>0)&&(WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
        WriteBinaryRestartFiles();
      
      // a label to jump to for a restart, everything before here is skipped
      // this is the point where the previous binary restart file was written
      ContinueAfterCrashLabel1: ;

      // Print at 'PrintEvery' intervals the status and a restart-file
      if((CurrentCycle%PrintEvery)==0) 
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          PrintIntervalStatusInit(CurrentCycle,NumberOfInitializationCycles,OutputFilePtr[CurrentSystem]);
          PrintRestartFile();
        }
      }

      // moves
      for(i=0;i<NumberOfSystems;i++)
      { 
        // choose system at random
        CurrentSystem=(int)(RandomNumber()*(REAL)NumberOfSystems);

        NumberOfSystemMoves=12;
        NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem]);
        NumberOfSteps=(NumberOfSystemMoves+NumberOfParticleMoves)*NumberOfComponents;

        for(j=0;j<NumberOfSteps;j++)
        {
          ran_int=(int)(RandomNumber()*NumberOfSteps);
          switch(ran_int)
          {
            case 0:
              if(RandomNumber()<ProbabilityParallelTemperingMove) ParallelTemperingMove();
              break;
            case 1:
              if(RandomNumber()<ProbabilityHyperParallelTemperingMove) HyperParallelTemperingMove();
              break;
            case 2:
              if(RandomNumber()<ProbabilityParallelMolFractionMove) ParallelMolFractionMove();
              break;
            case 3:
              if(RandomNumber()<ProbabilityChiralInversionMove) ChiralInversionMove();
              break;
            case 4:
              if(RandomNumber()<ProbabilityHybridNVEMove) HybridNVEMove();
              break;
            case 5:
              if(RandomNumber()<ProbabilityHybridNPHMove) HybridNPHMove();
              break;
            case 6:
              if(RandomNumber()<ProbabilityHybridNPHPRMove) HybridNPHPRMove();
              break;
            case 7:
              if(RandomNumber()<ProbabilityVolumeChangeMove) VolumeMove();
              break;
            case 8:
              if(RandomNumber()<ProbabilityBoxShapeChangeMove) BoxShapeChangeMove();
              break;
            case 9:
              if(RandomNumber()<ProbabilityGibbsVolumeChangeMove) GibbsVolumeMove();
              break;
            case 10:
              if(RandomNumber()<ProbabilityFrameworkChangeMove) FrameworkChangeMove();
              break;
            case 11:
              if(RandomNumber()<ProbabilityFrameworkShiftMove) FrameworkShiftMove();
              break;
            default:
              // choose component at random
              CurrentComponent=(int)(RandomNumber()*(REAL)NumberOfComponents);

              // choose the Monte Carlo move at random
              ran=RandomNumber();

              if(ran<Components[CurrentComponent].ProbabilityTranslationMove)
                TranslationMove();
              else if(ran<Components[CurrentComponent].ProbabilityRandomTranslationMove)
                RandomTranslationMove();
              else if(ran<Components[CurrentComponent].ProbabilityRotationMove)
                RotationMove();
              else if(ran<Components[CurrentComponent].ProbabilityPartialReinsertionMove)
                PartialReinsertionMove();
              else if(ran<Components[CurrentComponent].ProbabilityReinsertionMove)
                ReinsertionMove();
              else if(ran<Components[CurrentComponent].ProbabilityReinsertionInPlaceMove)
                ReinsertionInPlaceMove();
              else if(ran<Components[CurrentComponent].ProbabilityReinsertionInPlaneMove)
                ReinsertionInPlaneMove();
              else if(ran<Components[CurrentComponent].ProbabilityIdentityChangeMove)
                IdentityChangeMove();
              else if(ran<Components[CurrentComponent].ProbabilitySwapMove)
              {
                if(RandomNumber()<0.5) SwapAddMove();
                else SwapRemoveMove();
              }
              else if(ran<Components[CurrentComponent].ProbabilityWidomMove)
                ;
              else if(ran<Components[CurrentComponent].ProbabilitySurfaceAreaMove)
                ;
              else if(ran<Components[CurrentComponent].ProbabilityGibbsSwapChangeMove)
                GibbsParticleTransferMove();
              else if(ran<Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove)
                GibbsIdentityChangeMove();
              break;
          }
        }
      }

      if(CurrentCycle%PrintEvery==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          OptimizeVolumeChangeAcceptence();
          OptimizeTranslationAcceptence();
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            OptimizeFrameworkChangeAcceptence();
          OptimizeFrameworkShiftAcceptence();
          RescaleMaximumRotationAnglesSmallMC();
        }
      }
    }


    // after equilibration, recompute all the energies
    InitializesEnergiesAllSystems();
    InitializesEnergyAveragesAllSystems();

    InitializeSmallMCStatisticsAllSystems();
    InitializeMCMovesStatisticsAllSystems();

    CalculateTotalEnergyAllSystems();

    // initialize sampling-routines at the start of the production run
    SampleRadialDistributionFunction(INITIALIZE);
    SampleNumberOfMoleculesHistogram(INITIALIZE);
    SamplePositionHistogram(INITIALIZE);
    SampleFreeEnergyProfile(INITIALIZE);
    SamplePoreSizeDistribution(INITIALIZE);
    SampleEndToEndDistanceHistogram(INITIALIZE);
    SampleEnergyHistogram(INITIALIZE);
    SampleThermoDynamicsFactor(INITIALIZE);
    SampleFrameworkSpacingHistogram(INITIALIZE);
    SampleResidenceTimes(INITIALIZE);
    SampleDistanceHistogram(INITIALIZE);
    SampleBendAngleHistogram(INITIALIZE);
    SampleDihedralAngleHistogram(INITIALIZE);
    SampleAngleBetweenPlanesHistogram(INITIALIZE);
    SampleMoleculePropertyHistogram(INITIALIZE);
    SampleDensityProfile3DVTKGrid(INITIALIZE);
    SampleCationAndAdsorptionSites(INITIALIZE);
    SampleDcTSTConfigurationFiles(INITIALIZE);
    SamplePDBMovies(INITIALIZE,-1);

    SimulationStage=PRODUCTION;
    for(CurrentCycle=0;CurrentCycle<NumberOfCycles;CurrentCycle++)
    {

      // detect erroneous chirality changes
      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        CheckChiralityMolecules();

      // sample energy average and system/particle properties
      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      {
        UpdateEnergyAveragesCurrentSystem();

        SampleRadialDistributionFunction(SAMPLE);
        SampleNumberOfMoleculesHistogram(SAMPLE);
        SamplePositionHistogram(SAMPLE);
        SampleFreeEnergyProfile(SAMPLE);
        SamplePoreSizeDistribution(SAMPLE);
        SampleEndToEndDistanceHistogram(SAMPLE);
        SampleEnergyHistogram(SAMPLE);
        SampleThermoDynamicsFactor(SAMPLE);
        SampleFrameworkSpacingHistogram(SAMPLE);
        SampleResidenceTimes(SAMPLE);
        SampleDistanceHistogram(SAMPLE);
        SampleBendAngleHistogram(SAMPLE);
        SampleDihedralAngleHistogram(SAMPLE);
        SampleAngleBetweenPlanesHistogram(SAMPLE);
        SampleMoleculePropertyHistogram(SAMPLE);
        SampleDensityProfile3DVTKGrid(SAMPLE);
        SampleCationAndAdsorptionSites(SAMPLE);
        SampleDcTSTConfigurationFiles(SAMPLE);
        SamplePDBMovies(SAMPLE,-1);
      }

      if(CurrentCycle%PrintPropertiesEvery==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
          PrintPropertyStatus(CurrentCycle,NumberOfCycles,OutputFilePtr[CurrentSystem]);
      }

      // Print at 'PrintEvery' intervals the status and a restart-file
      if((CurrentCycle%PrintEvery)==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          PrintIntervalStatus(CurrentCycle,NumberOfCycles,OutputFilePtr[CurrentSystem]);
          PrintRestartFile();
        }
      }

      // moves
      for(i=0;i<NumberOfSystems;i++)
      {
        // choose system at random
        CurrentSystem=(int)(RandomNumber()*(REAL)NumberOfSystems);

        NumberOfSystemMoves=12;
        NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem]);
        NumberOfSteps=(NumberOfSystemMoves+NumberOfParticleMoves)*NumberOfComponents;

        // loop over the MC 'steps' per MC 'cycle'
        for(j=0;j<NumberOfSteps;j++)
        {
          // choose any of the MC moves randomly
          ran_int=(int)(RandomNumber()*NumberOfSteps);
          switch(ran_int)
          {
            case 0:
              if(RandomNumber()<ProbabilityParallelTemperingMove) ParallelTemperingMove();
              break;
            case 1:
              if(RandomNumber()<ProbabilityHyperParallelTemperingMove) HyperParallelTemperingMove();
              break;
            case 2:
              if(RandomNumber()<ProbabilityParallelMolFractionMove) ParallelMolFractionMove();
              break;
            case 3:
              if(RandomNumber()<ProbabilityChiralInversionMove) ChiralInversionMove();
              break;
            case 4:
              if(RandomNumber()<ProbabilityHybridNVEMove) HybridNVEMove();
              break;
            case 5:
              if(RandomNumber()<ProbabilityHybridNPHMove) HybridNPHMove();
              break;
            case 6:
              if(RandomNumber()<ProbabilityHybridNPHPRMove) HybridNPHPRMove();
              break;
            case 7:
              if(RandomNumber()<ProbabilityVolumeChangeMove) VolumeMove();
              break;
            case 8:
              if(RandomNumber()<ProbabilityBoxShapeChangeMove) BoxShapeChangeMove();
              break;
            case 9:
              if(RandomNumber()<ProbabilityGibbsVolumeChangeMove) GibbsVolumeMove();
              break;
            case 10:
              if(RandomNumber()<ProbabilityFrameworkChangeMove) FrameworkChangeMove();
              break;
            case 11:
              if(RandomNumber()<ProbabilityFrameworkShiftMove) FrameworkShiftMove();
              break;
            default:
              // choose component at random
              CurrentComponent=(int)(RandomNumber()*(REAL)NumberOfComponents);

              // choose the Monte Carlo move at random
              ran=RandomNumber();
              if(ran<Components[CurrentComponent].ProbabilityTranslationMove)
                TranslationMove();
              else if(ran<Components[CurrentComponent].ProbabilityRandomTranslationMove)
                RandomTranslationMove();
              else if(ran<Components[CurrentComponent].ProbabilityRotationMove)
                RotationMove();
              else if(ran<Components[CurrentComponent].ProbabilityPartialReinsertionMove)
                PartialReinsertionMove();          
              else if(ran<Components[CurrentComponent].ProbabilityReinsertionMove)
                ReinsertionMove();          
              else if(ran<Components[CurrentComponent].ProbabilityReinsertionInPlaceMove)
                ReinsertionInPlaceMove();          
              else if(ran<Components[CurrentComponent].ProbabilityReinsertionInPlaneMove)
                ReinsertionInPlaneMove();
              else if(ran<Components[CurrentComponent].ProbabilityIdentityChangeMove)
                IdentityChangeMove();          
              else if(ran<Components[CurrentComponent].ProbabilitySwapMove)
              {
                if(RandomNumber()<0.5) SwapAddMove();
                else SwapRemoveMove();
              }
              else if(ran<Components[CurrentComponent].ProbabilityWidomMove)
                WidomMove();
              else if(ran<Components[CurrentComponent].ProbabilitySurfaceAreaMove)
                SurfaceAreaMove();
              else if(ran<Components[CurrentComponent].ProbabilityGibbsSwapChangeMove)
                GibbsParticleTransferMove();
              else if(ran<Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove)
                GibbsIdentityChangeMove();
              break;
          }
        }
      }

      if(CurrentCycle%PrintEvery==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          OptimizeVolumeChangeAcceptence();
          OptimizeTranslationAcceptence();
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            OptimizeFrameworkChangeAcceptence();
          OptimizeFrameworkShiftAcceptence();
          RescaleMaximumRotationAnglesSmallMC();
        }
      }

      if((CurrentCycle>0)&&(WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
        WriteBinaryRestartFiles();

      // a label to jump to for a restart, everything before here is skipped
      // this is the point where the previous binary restart file was written
      ContinueAfterCrashLabel2: ;

      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      {
        SampleRadialDistributionFunction(PRINT);
        SampleNumberOfMoleculesHistogram(PRINT);
        SamplePositionHistogram(PRINT);
        SampleFreeEnergyProfile(PRINT);
        SamplePoreSizeDistribution(PRINT);
        SampleEndToEndDistanceHistogram(PRINT);
        SampleEnergyHistogram(PRINT);
        SampleThermoDynamicsFactor(PRINT);
        SampleFrameworkSpacingHistogram(PRINT);
        SampleResidenceTimes(PRINT);
        SampleDistanceHistogram(PRINT);
        SampleBendAngleHistogram(PRINT);
        SampleDihedralAngleHistogram(PRINT);
        SampleAngleBetweenPlanesHistogram(PRINT);
        SampleMoleculePropertyHistogram(PRINT);
        SampleDensityProfile3DVTKGrid(PRINT);
        SampleCationAndAdsorptionSites(PRINT);
        SampleDcTSTConfigurationFiles(PRINT);
        if(Movies[CurrentSystem]&&(CurrentCycle%WriteMoviesEvery[CurrentSystem]==0))
          SamplePDBMovies(PRINT,-1);
      }
    }

    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      PrintRestartFile();

    PrintPostSimulationStatus();

    CloseOutputFile();
  }

  // set current prssure to the last one
  CurrentIsothermPressure=NumberOfIsothermPressures-1;

  // write last binary restart-file
  // sampled properties can be remade from the binary restart-file
  if(WriteBinaryRestartFileEvery>0)
     WriteBinaryRestartFiles();

  // finalize output
  SampleRadialDistributionFunction(FINALIZE);
  SampleNumberOfMoleculesHistogram(FINALIZE);
  SamplePositionHistogram(FINALIZE);
  SampleFreeEnergyProfile(FINALIZE);
  SamplePoreSizeDistribution(FINALIZE);
  SampleEndToEndDistanceHistogram(FINALIZE);
  SampleEnergyHistogram(FINALIZE);
  SampleThermoDynamicsFactor(FINALIZE);
  SampleFrameworkSpacingHistogram(FINALIZE);
  SampleResidenceTimes(FINALIZE);
  SampleDistanceHistogram(FINALIZE);
  SampleBendAngleHistogram(FINALIZE);
  SampleDihedralAngleHistogram(FINALIZE);
  SampleAngleBetweenPlanesHistogram(FINALIZE);
  SampleMoleculePropertyHistogram(FINALIZE);
  SampleDensityProfile3DVTKGrid(FINALIZE);
  SampleCationAndAdsorptionSites(FINALIZE);
  SampleDcTSTConfigurationFiles(FINALIZE);
  SamplePDBMovies(FINALIZE,-1);
}


/*********************************************************************************************************
 * Name       | ComputePoreSizeDistribution                                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Main routine to compute the Pore-Size Distribution function (PSD).                       *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void ComputePoreSizeDistribution(void)
{
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    ComputePSDHistogram[CurrentSystem]=TRUE;

  SamplePoreSizeDistribution(ALLOCATE);
  SamplePoreSizeDistribution(INITIALIZE);

  for(CurrentCycle=0;CurrentCycle<NumberOfCycles;CurrentCycle++)
  {
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      if((CurrentCycle%PrintEvery)==0)
      {
        fprintf(OutputFilePtr[CurrentSystem],"Iteration: %lld\n",CurrentCycle);
        fflush(OutputFilePtr[CurrentSystem]);
      }
      SamplePoreSizeDistribution(SAMPLE);
      SamplePoreSizeDistribution(PRINT);
    }
  }

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
     SamplePoreSizeDistribution(PRINT);

  SamplePoreSizeDistribution(FINALIZE);
}
