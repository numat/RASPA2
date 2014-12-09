/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2015 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'monte_carlo.c' is part of RASPA-2.0

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
 *************************************************************************************************************/

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
#include "ewald.h"

void DebugEnergyStatus(void);

/*********************************************************************************************************
 * Name       | MonteCarloSimulation                                                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Main routine to do a Monte-Carlo (MC) simulation.                                        *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void MonteCarloSimulation(void)
{
  int i,j;
  int NumberOfParticleMoves;
  int NumberOfSteps;
  int SelectedSystem;
  double cpu_start,cpu_end;
  double cpu_before,cpu_after;
  REAL ran;


  // for a crash-recovery we skip the initialization/equilibration and jump to the
  // position right after the crash-file was written in the production run
  if(ContinueAfterCrash)
  {
    if(SimulationStage==POSITION_INITIALIZATION) goto ContinueAfterCrashLabel1;
    else if(SimulationStage==CF_WANG_LANDAU_EQUILIBRATION) goto ContinueAfterCrashLabel2;
    else if(SimulationStage==PRODUCTION) goto ContinueAfterCrashLabel3;
    else goto ContinueAfterCrashLabel4;
  }

  // allocate memory for sampling routines
  SampleRadialDistributionFunction(ALLOCATE);
  SampleProjectedLengthsDistributionFunction(ALLOCATE);
  SampleProjectedAnglesDistributionFunction(ALLOCATE);
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
  SampleBondOrientationAutoCorrelationFunctionOrderN(ALLOCATE);
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

    CFWangLandauIteration(INITIALIZE);
    CFRXMCWangLandauIteration(INITIALIZE);

    // initialization to reach equilibration of positions (no averages are computed yet)
    SimulationStage=POSITION_INITIALIZATION;
    for(CurrentCycle=0;CurrentCycle<NumberOfInitializationCycles;CurrentCycle++)
    {
      if((WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
        WriteBinaryRestartFiles();

      // a label to jump to for a restart, everything before here is skipped
      // this is the point where the previous binary restart file was written
      ContinueAfterCrashLabel1: ;

      cpu_start=get_cpu_time();

      // detect erroneous chirality changes
      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        CheckChiralityMolecules();

      // Print at 'PrintEvery' intervals the status and a restart-file
      if((CurrentCycle%PrintEvery)==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          PrintIntervalStatusInit(CurrentCycle,NumberOfInitializationCycles,OutputFilePtr[CurrentSystem]);
          PrintRestartFile();
        }
      }

      for(i=0;i<NumberOfSystems;i++)
      {
        // choose a random system
        SelectedSystem=(int)(RandomNumber()*(REAL)NumberOfSystems);

        NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[SelectedSystem]+NumberOfCationMolecules[SelectedSystem]);
        NumberOfSteps=NumberOfParticleMoves*(NumberOfComponents==0?1:NumberOfComponents);

        for(j=0;j<NumberOfSteps;j++)
        {
          // set the selected system
          CurrentSystem=SelectedSystem;

          // choose a random component
          CurrentComponent=(int)(RandomNumber()*(REAL)NumberOfComponents);

          // choose any of the MC moves randomly with the selected probability
          ran=RandomNumber();

          if(ran<Components[CurrentComponent].ProbabilityTranslationMove)
            TranslationMove();
          else if(ran<Components[CurrentComponent].ProbabilityRandomTranslationMove)
            RandomTranslationMove();
          else if(ran<Components[CurrentComponent].ProbabilityRotationMove)
            RotationMove();
          else if(ran<Components[CurrentComponent].ProbabilityRandomRotationMove)
            RandomRotationMove();
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
          else if(ran<Components[CurrentComponent].ProbabilityCFSwapLambdaMove)
            CFSwapLambaMove();
          else if(ran<Components[CurrentComponent].ProbabilityCBCFSwapLambdaMove)
            CBCFSwapLambaMove();
          else if(ran<Components[CurrentComponent].ProbabilityWidomMove)
            ;
          else if(ran<Components[CurrentComponent].ProbabilitySurfaceAreaMove)
            ;
          else if(ran<Components[CurrentComponent].ProbabilityGibbsChangeMove)
            GibbsParticleTransferMove();
          else if(ran<Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove)
            GibbsIdentityChangeMove();
          else if(ran<Components[CurrentComponent].ProbabilityCFGibbsChangeMove)
            CFGibbsParticleTransferMove();
          else if(ran<Components[CurrentComponent].ProbabilityCBCFGibbsChangeMove)
            CBCFGibbsParticleTransferMove();
          else if(ran<Components[CurrentComponent].ProbabilityParallelTemperingMove)
            ParallelTemperingMove();
          else if(ran<Components[CurrentComponent].ProbabilityHyperParallelTemperingMove)
            HyperParallelTemperingMove();
          else if(ran<Components[CurrentComponent].ProbabilityParallelMolFractionMove)
            ParallelMolFractionMove();
          else if(ran<Components[CurrentComponent].ProbabilityChiralInversionMove)
            ChiralInversionMove();
          else if(ran<Components[CurrentComponent].ProbabilityHybridNVEMove)
            HybridNVEMove();
          else if(ran<Components[CurrentComponent].ProbabilityHybridNPHMove)
            HybridNPHMove();
          else if(ran<Components[CurrentComponent].ProbabilityHybridNPHPRMove)
            HybridNPHPRMove();
          else if(ran<Components[CurrentComponent].ProbabilityVolumeChangeMove)
            VolumeMove();
          else if(ran<Components[CurrentComponent].ProbabilityBoxShapeChangeMove)
            BoxShapeChangeMove();
          else if(ran<Components[CurrentComponent].ProbabilityGibbsVolumeChangeMove)
            GibbsVolumeMove();
          else if(ran<Components[CurrentComponent].ProbabilityFrameworkChangeMove)
            FrameworkChangeMove();
          else if(ran<Components[CurrentComponent].ProbabilityFrameworkShiftMove)
            FrameworkShiftMove();
          else if(ran<Components[CurrentComponent].ProbabilityCFCRXMCLambdaChangeMove)
            CFCRXMCLambdaChangeMove();

          #ifdef DEBUG
            DebugEnergyStatus();
          #endif
        }
      }

      if(CurrentCycle%OptimizeAcceptenceEvery==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          OptimizeVolumeChangeAcceptence();
          OptimizeGibbsVolumeChangeAcceptence();
          OptimizeTranslationAcceptence();
          OptimizeRotationAcceptence();
          OptimizeFrameworkChangeAcceptence();
          OptimizeFrameworkShiftAcceptence();
          OptimizeCFLambdaChangeAcceptence();
          OptimizeCFGibbsLambdaChangeAcceptence();
          OptimizeCBCFLambdaChangeAcceptence();
          OptimizeCBCFGibbsLambdaChangeAcceptence();
          OptimizeRXMCLambdaChangeAcceptence();
          RescaleMaximumRotationAnglesSmallMC();
        }
      }
      cpu_end=get_cpu_time();
      CpuTotal+=(cpu_end-cpu_start);
      CpuTimeInitialization+=(cpu_end-cpu_start);

    }



    // initialize the energies and compute the total energies for all systems
    InitializesEnergiesAllSystems();
    InitializesEnergyAveragesAllSystems();

    InitializeSmallMCStatisticsAllSystems();
    InitializeMCMovesStatisticsAllSystems();

    // compute total energy for all systems
    CalculateTotalEnergyAllSystems();

    if(NumberOfEquilibrationCycles>0)
    {
      CFWangLandauIteration(INITIALIZE);
      CFRXMCWangLandauIteration(INITIALIZE);

      // initialization to reach equilibration of positions (no averages are computed yet)
      SimulationStage=CF_WANG_LANDAU_EQUILIBRATION;
      for(CurrentCycle=0;CurrentCycle<NumberOfEquilibrationCycles;CurrentCycle++)
      {
        if((WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
          WriteBinaryRestartFiles();

        // a label to jump to for a restart, everything before here is skipped
        // this is the point where the previous binary restart file was written
        ContinueAfterCrashLabel2: ;

        cpu_start=get_cpu_time();

        // detect erroneous chirality changes
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
          CheckChiralityMolecules();

        // Print at 'PrintEvery' intervals the status and a restart-file
        if((CurrentCycle%PrintEvery)==0)
        {
          for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
          {
            PrintIntervalStatusEquilibration(CurrentCycle,NumberOfEquilibrationCycles,OutputFilePtr[CurrentSystem]);
            PrintRestartFile();
          }
        }

        // moves
        for(i=0;i<NumberOfSystems;i++)
        {
          // choose a random system
          SelectedSystem=(int)(RandomNumber()*(REAL)NumberOfSystems);

          NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[SelectedSystem]+NumberOfCationMolecules[SelectedSystem]);
          NumberOfSteps=NumberOfParticleMoves*(NumberOfComponents==0?1:NumberOfComponents);

          for(j=0;j<NumberOfSteps;j++)
          {
            // set the selected system
            CurrentSystem=SelectedSystem;

            // choose a random component
            CurrentComponent=(int)(RandomNumber()*(REAL)NumberOfComponents);

            // choose any of the MC moves randomly with the selected probability
            ran=RandomNumber();

            if(ran<Components[CurrentComponent].ProbabilityTranslationMove)
              TranslationMove();
            else if(ran<Components[CurrentComponent].ProbabilityRandomTranslationMove)
              RandomTranslationMove();
            else if(ran<Components[CurrentComponent].ProbabilityRotationMove)
              RotationMove();
            else if(ran<Components[CurrentComponent].ProbabilityRandomRotationMove)
              RandomRotationMove();
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
            else if(ran<Components[CurrentComponent].ProbabilityCFSwapLambdaMove)
            {
              CFSwapLambaMove();
              CFWangLandauIteration(SAMPLE);
            }
            else if(ran<Components[CurrentComponent].ProbabilityCBCFSwapLambdaMove)
            {
              CFWangLandauIteration(SAMPLE);
              CBCFSwapLambaMove();
            }
            else if(ran<Components[CurrentComponent].ProbabilityWidomMove)
              ;
            else if(ran<Components[CurrentComponent].ProbabilitySurfaceAreaMove)
              ;
            else if(ran<Components[CurrentComponent].ProbabilityGibbsChangeMove)
              GibbsParticleTransferMove();
            else if(ran<Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove)
              GibbsIdentityChangeMove();
            else if(ran<Components[CurrentComponent].ProbabilityCFGibbsChangeMove)
            {
              CFWangLandauIteration(SAMPLE);
              CFGibbsParticleTransferMove();
            }
            else if(ran<Components[CurrentComponent].ProbabilityCBCFGibbsChangeMove)
            {
              CFWangLandauIteration(SAMPLE);
              CBCFGibbsParticleTransferMove();
            }
            else if(ran<Components[CurrentComponent].ProbabilityParallelTemperingMove)
              ParallelTemperingMove();
            else if(ran<Components[CurrentComponent].ProbabilityHyperParallelTemperingMove)
              HyperParallelTemperingMove();
            else if(ran<Components[CurrentComponent].ProbabilityParallelMolFractionMove)
              ParallelMolFractionMove();
            else if(ran<Components[CurrentComponent].ProbabilityChiralInversionMove)
              ChiralInversionMove();
            else if(ran<Components[CurrentComponent].ProbabilityHybridNVEMove)
              HybridNVEMove();
            else if(ran<Components[CurrentComponent].ProbabilityHybridNPHMove)
              HybridNPHMove();
            else if(ran<Components[CurrentComponent].ProbabilityHybridNPHPRMove)
              HybridNPHPRMove();
            else if(ran<Components[CurrentComponent].ProbabilityVolumeChangeMove)
              VolumeMove();
            else if(ran<Components[CurrentComponent].ProbabilityBoxShapeChangeMove)
              BoxShapeChangeMove();
            else if(ran<Components[CurrentComponent].ProbabilityGibbsVolumeChangeMove)
              GibbsVolumeMove();
            else if(ran<Components[CurrentComponent].ProbabilityFrameworkChangeMove)
              FrameworkChangeMove();
            else if(ran<Components[CurrentComponent].ProbabilityFrameworkShiftMove)
              FrameworkShiftMove();
            else if(ran<Components[CurrentComponent].ProbabilityCFCRXMCLambdaChangeMove)
            {
              CFCRXMCLambdaChangeMove();
              CFRXMCWangLandauIteration(SAMPLE);
            }

            #ifdef DEBUG
              DebugEnergyStatus();
            #endif
          }
        }

        if(CurrentCycle%OptimizeAcceptenceEvery==0)
        {
          for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
          {
            OptimizeVolumeChangeAcceptence();
            OptimizeGibbsVolumeChangeAcceptence();
            OptimizeTranslationAcceptence();
            OptimizeRotationAcceptence();
            OptimizeFrameworkChangeAcceptence();
            OptimizeFrameworkShiftAcceptence();
            OptimizeCFLambdaChangeAcceptence();
            OptimizeCFGibbsLambdaChangeAcceptence();
            OptimizeCBCFLambdaChangeAcceptence();
            OptimizeCBCFGibbsLambdaChangeAcceptence();
            OptimizeRXMCLambdaChangeAcceptence();
            RescaleMaximumRotationAnglesSmallMC();
          }
        }

        if(CurrentCycle%CFWangLandauEvery==0)
        {
          CFWangLandauIteration(PRINT);
          CFRXMCWangLandauIteration(PRINT);
        }

        cpu_end=get_cpu_time();
        CpuTotal+=(cpu_end-cpu_start);
        CpuTimeEquilibration+=(cpu_end-cpu_start);
      }

      CFWangLandauIteration(FINALIZE);
      CFRXMCWangLandauIteration(FINALIZE);

      // after equilibration, recompute all the energies
      InitializesEnergiesAllSystems();
      InitializesEnergyAveragesAllSystems();

      InitializeSmallMCStatisticsAllSystems();
      InitializeMCMovesStatisticsAllSystems();

      CalculateTotalEnergyAllSystems();
    }

    // initialize sampling-routines at the start of the production run
    SampleRadialDistributionFunction(INITIALIZE);
    SampleProjectedLengthsDistributionFunction(INITIALIZE);
    SampleProjectedAnglesDistributionFunction(INITIALIZE);
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
      cpu_start=get_cpu_time();

      // detect erroneous chirality changes
      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        CheckChiralityMolecules();

      // sample energy average and system/particle properties
      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      {
        UpdateEnergyAveragesCurrentSystem();

        SampleRadialDistributionFunction(SAMPLE);
        SampleProjectedLengthsDistributionFunction(SAMPLE);
        SampleProjectedAnglesDistributionFunction(SAMPLE);
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
        // choose a random system
        SelectedSystem=(int)(RandomNumber()*(REAL)NumberOfSystems);

        NumberOfParticleMoves=MAX2(MinimumInnerCycles,NumberOfAdsorbateMolecules[SelectedSystem]+NumberOfCationMolecules[SelectedSystem]);
        NumberOfSteps=NumberOfParticleMoves*(NumberOfComponents==0?1:NumberOfComponents);

        // loop over the MC 'steps' per MC 'cycle'
        for(j=0;j<NumberOfSteps;j++)
        {
          // set the selected system
          CurrentSystem=SelectedSystem;

          // choose a random component
          CurrentComponent=(int)(RandomNumber()*(REAL)NumberOfComponents);

          // choose any of the MC moves randomly with the selected probability
          ran=RandomNumber();

          if(ran<Components[CurrentComponent].ProbabilityTranslationMove)
          {
            cpu_before=get_cpu_time();
            TranslationMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeTranslationMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityRandomTranslationMove)
          {
            cpu_before=get_cpu_time();
            RandomTranslationMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeRandomTranslationMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityRotationMove)
          {
            cpu_before=get_cpu_time();
            RotationMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeRotationMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityRandomRotationMove)
          {
            cpu_before=get_cpu_time();
            RandomRotationMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeRandomRotationMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityPartialReinsertionMove)
          {
            cpu_before=get_cpu_time();
            PartialReinsertionMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimePartialReinsertionMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityReinsertionMove)
          {
            cpu_before=get_cpu_time();
            ReinsertionMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeReinsertionMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityReinsertionInPlaceMove)
          {
            cpu_before=get_cpu_time();
            ReinsertionInPlaceMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeReinsertionInPlaceMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityReinsertionInPlaneMove)
          {
            cpu_before=get_cpu_time();
            ReinsertionInPlaneMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeReinsertionInPlaneMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityIdentityChangeMove)
          {
            cpu_before=get_cpu_time();
            IdentityChangeMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeIdentityChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilitySwapMove)
          {
            if(RandomNumber()<0.5) 
            {
              cpu_before=get_cpu_time();
              SwapAddMove();
              cpu_after=get_cpu_time();
              Components[CurrentComponent].CpuTimeSwapMoveInsertion[CurrentSystem]+=(cpu_after-cpu_before);
            }
            else 
            {
              cpu_before=get_cpu_time();
              SwapRemoveMove();
              cpu_after=get_cpu_time();
              Components[CurrentComponent].CpuTimeSwapMoveDeletion[CurrentSystem]+=(cpu_after-cpu_before);
            }
          }
          else if(ran<Components[CurrentComponent].ProbabilityCFSwapLambdaMove)
          {
            cpu_before=get_cpu_time();
            CFSwapLambaMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCFSwapLambdaMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCBCFSwapLambdaMove)
          {
            cpu_before=get_cpu_time();
            CBCFSwapLambaMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCBCFSwapLambdaMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityWidomMove)
          {
            cpu_before=get_cpu_time();
            WidomMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeWidomMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilitySurfaceAreaMove)
          {
            cpu_before=get_cpu_time();
            SurfaceAreaMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeSurfaceAreaMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityGibbsChangeMove)
          {
            cpu_before=get_cpu_time();
            GibbsParticleTransferMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeGibbsChangeMove[0]+=0.5*(cpu_after-cpu_before);
            Components[CurrentComponent].CpuTimeGibbsChangeMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove)
          {
            cpu_before=get_cpu_time();
            GibbsIdentityChangeMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeGibbsIdentityChangeMove[0]+=0.5*(cpu_after-cpu_before);
            Components[CurrentComponent].CpuTimeGibbsIdentityChangeMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCFGibbsChangeMove)
          {
            cpu_before=get_cpu_time();
            CFGibbsParticleTransferMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCFGibbsChangeMove[0]+=0.5*(cpu_after-cpu_before);
            Components[CurrentComponent].CpuTimeCFGibbsChangeMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCBCFGibbsChangeMove)
          {
            cpu_before=get_cpu_time();
            CBCFGibbsParticleTransferMove();
            cpu_after=get_cpu_time();
            Components[CurrentComponent].CpuTimeCBCFGibbsChangeMove[0]+=0.5*(cpu_after-cpu_before);
            Components[CurrentComponent].CpuTimeCBCFGibbsChangeMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityParallelTemperingMove) 
          {
            // note: cpu-timings done in the move itself
            ParallelTemperingMove();
          }
          else if(ran<Components[CurrentComponent].ProbabilityHyperParallelTemperingMove) 
          {
            // note: cpu-timings done in the move itself
            HyperParallelTemperingMove();
          }
          else if(ran<Components[CurrentComponent].ProbabilityParallelMolFractionMove) 
          {
            // note: cpu-timings done in the move itself
            ParallelMolFractionMove();
          }
          else if(ran<Components[CurrentComponent].ProbabilityChiralInversionMove) 
          {
            cpu_before=get_cpu_time();
            ChiralInversionMove();
            cpu_after=get_cpu_time();
            CpuTimeChiralInversionMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityHybridNVEMove) 
          {
            cpu_before=get_cpu_time();
            HybridNVEMove();
            cpu_after=get_cpu_time();
            CpuTimeHybridNVEMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityHybridNPHMove) 
          {
            cpu_before=get_cpu_time();
            HybridNPHMove();
            cpu_after=get_cpu_time();
            CpuTimeHybridNPHMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityHybridNPHPRMove) 
          {
            cpu_before=get_cpu_time();
            HybridNPHPRMove();
            cpu_after=get_cpu_time();
            CpuTimeHybridNPHPRMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityVolumeChangeMove) 
          {
            cpu_before=get_cpu_time();
            VolumeMove();
            cpu_after=get_cpu_time();
            CpuTimeVolumeChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityBoxShapeChangeMove) 
          {
            cpu_before=get_cpu_time();
            BoxShapeChangeMove();
            cpu_after=get_cpu_time();
            CpuTimeBoxShapeChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityGibbsVolumeChangeMove) 
          {
            cpu_before=get_cpu_time();
            GibbsVolumeMove();
            cpu_after=get_cpu_time();
            CpuTimeGibbsVolumeChangeMove[0]+=0.5*(cpu_after-cpu_before);
            CpuTimeGibbsVolumeChangeMove[1]+=0.5*(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityFrameworkChangeMove) 
          {
            cpu_before=get_cpu_time();
            FrameworkChangeMove();
            cpu_after=get_cpu_time();
            CpuTimeFrameworkChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityFrameworkShiftMove) 
          {
            cpu_before=get_cpu_time();
            FrameworkShiftMove();
            cpu_after=get_cpu_time();
            CpuTimeFrameworkShiftMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          else if(ran<Components[CurrentComponent].ProbabilityCFCRXMCLambdaChangeMove)
          {
            cpu_before=get_cpu_time();
            CFCRXMCLambdaChangeMove();
            cpu_after=get_cpu_time();
            CpuCFCRXMCLambdaChangeMove[CurrentSystem]+=(cpu_after-cpu_before);
          }
          
          #ifdef DEBUG
            DebugEnergyStatus();
          #endif
        }
      }


      if(CurrentCycle%OptimizeAcceptenceEvery==0)
      {
        for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
        {
          OptimizeVolumeChangeAcceptence();
          OptimizeGibbsVolumeChangeAcceptence();
          OptimizeTranslationAcceptence();
          OptimizeRotationAcceptence();
          OptimizeFrameworkChangeAcceptence();
          OptimizeFrameworkShiftAcceptence();
          OptimizeCFLambdaChangeAcceptence();
          OptimizeCFGibbsLambdaChangeAcceptence();
          OptimizeCBCFLambdaChangeAcceptence();
          OptimizeCBCFGibbsLambdaChangeAcceptence();
          OptimizeRXMCLambdaChangeAcceptence();
          RescaleMaximumRotationAnglesSmallMC();
        }
      }

      cpu_end=get_cpu_time();
      CpuTotal+=(cpu_end-cpu_start);
      CpuTimeProductionRun+=(cpu_end-cpu_start);

      if((WriteBinaryRestartFileEvery>0)&&(CurrentCycle%WriteBinaryRestartFileEvery==0))
        WriteBinaryRestartFiles();

      // a label to jump to for a restart, everything before here is skipped
      // this is the point where the previous binary restart file was written
      ContinueAfterCrashLabel3: ;

      for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      {
        SampleRadialDistributionFunction(PRINT);
        SampleProjectedLengthsDistributionFunction(PRINT);
        SampleProjectedAnglesDistributionFunction(PRINT);
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

    SimulationStage=FINISHED;

    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      OptimizeVolumeChangeAcceptence();
      OptimizeGibbsVolumeChangeAcceptence();
      OptimizeTranslationAcceptence();
      OptimizeRotationAcceptence();
      OptimizeFrameworkChangeAcceptence();
      OptimizeFrameworkShiftAcceptence();
      OptimizeCFLambdaChangeAcceptence();
      OptimizeCFGibbsLambdaChangeAcceptence();
      OptimizeCBCFLambdaChangeAcceptence();
      OptimizeCBCFGibbsLambdaChangeAcceptence();
      OptimizeRXMCLambdaChangeAcceptence();
      RescaleMaximumRotationAnglesSmallMC();
    }


    // write last binary restart-file
    // sampled properties can be remade from the binary restart-file
    if(WriteBinaryRestartFileEvery>0)
       WriteBinaryRestartFiles();

    ContinueAfterCrashLabel4: ;

    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      SampleRadialDistributionFunction(PRINT);
      SampleProjectedLengthsDistributionFunction(PRINT);
      SampleProjectedAnglesDistributionFunction(PRINT);
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
    }

    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      PrintRestartFile();

    PrintPostSimulationStatus();

    CloseOutputFile();
  }

  // set current prssure to the last one
  CurrentIsothermPressure=NumberOfIsothermPressures-1;


  // finalize output
  SampleRadialDistributionFunction(FINALIZE);
  SampleProjectedLengthsDistributionFunction(FINALIZE);
  SampleProjectedAnglesDistributionFunction(FINALIZE);
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
 * Name       | CheckEnergyOfSystem                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Routine to check the energy state of the system                                          *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void DebugEnergyStatus(void)
{
  int i;

  REAL StoredUHostBond,StoredUHostUreyBradley,StoredUHostBend,StoredUHostInversionBend;
  REAL StoredUHostTorsion,StoredUHostImproperTorsion,StoredUHostBondBond;
  REAL StoredUHostBendBend,StoredUHostBondBend,StoredUHostBondTorsion,StoredUHostBendTorsion;

  REAL StoredUAdsorbateBond,StoredUAdsorbateUreyBradley,StoredUAdsorbateBend,StoredUAdsorbateInversionBend;
  REAL StoredUAdsorbateTorsion,StoredUAdsorbateImproperTorsion,StoredUAdsorbateBondBond;
  REAL StoredUAdsorbateBendBend,StoredUAdsorbateBondBend,StoredUAdsorbateBondTorsion,StoredUAdsorbateBendTorsion;
  REAL StoredUAdsorbateIntraVDW,StoredUAdsorbateIntraChargeCharge;
  REAL StoredUAdsorbateIntraChargeBondDipole,StoredUAdsorbateIntraBondDipoleBondDipole;

  REAL StoredUCationBond,StoredUCationUreyBradley,StoredUCationBend,StoredUCationInversionBend;
  REAL StoredUCationTorsion,StoredUCationImproperTorsion,StoredUCationBondBond;
  REAL StoredUCationBendBend,StoredUCationBondBend,StoredUCationBondTorsion,StoredUCationBendTorsion;
  REAL StoredUCationIntraVDW,StoredUCationIntraChargeCharge;
  REAL StoredUCationIntraChargeBondDipole,StoredUCationIntraBondDipoleBondDipole;

  REAL StoredUHostHost,StoredUHostHostVDW,StoredUHostHostChargeChargeReal;
  REAL StoredUHostHostChargeBondDipoleReal,StoredUHostHostBondDipoleBondDipoleReal;
  REAL StoredUHostHostChargeChargeFourier,StoredUHostHostCoulomb;
  REAL StoredUHostHostChargeBondDipoleFourier,StoredUHostHostBondDipoleBondDipoleFourier;
  REAL StoredUHostAdsorbate,StoredUHostAdsorbateVDW,StoredUHostAdsorbateChargeChargeReal;
  REAL StoredUHostAdsorbateChargeBondDipoleReal,StoredUHostAdsorbateBondDipoleBondDipoleReal;
  REAL StoredUHostAdsorbateChargeChargeFourier,StoredUHostAdsorbateCoulomb;
  REAL StoredUHostAdsorbateChargeBondDipoleFourier,StoredUHostAdsorbateBondDipoleBondDipoleFourier;
  REAL StoredUHostCation,StoredUHostCationVDW,StoredUHostCationChargeChargeReal;
  REAL StoredUHostCationChargeBondDipoleReal,StoredUHostCationBondDipoleBondDipoleReal;
  REAL StoredUHostCationChargeChargeFourier,StoredUHostCationCoulomb;
  REAL StoredUHostCationChargeBondDipoleFourier,StoredUHostCationBondDipoleBondDipoleFourier;

  REAL StoredUAdsorbateAdsorbate,StoredUAdsorbateAdsorbateVDW,StoredUAdsorbateAdsorbateChargeChargeReal;
  REAL StoredUAdsorbateAdsorbateChargeBondDipoleReal,StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal;
  REAL StoredUAdsorbateAdsorbateChargeChargeFourier,StoredUAdsorbateAdsorbateCoulomb;
  REAL StoredUAdsorbateAdsorbateChargeBondDipoleFourier,StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier;
  REAL StoredUAdsorbateCation,StoredUAdsorbateCationVDW,StoredUAdsorbateCationChargeChargeReal;
  REAL StoredUAdsorbateCationChargeBondDipoleReal,StoredUAdsorbateCationBondDipoleBondDipoleReal;
  REAL StoredUAdsorbateCationChargeChargeFourier,StoredUAdsorbateCationCoulomb;
  REAL StoredUAdsorbateCationChargeBondDipoleFourier,StoredUAdsorbateCationBondDipoleBondDipoleFourier;
  REAL StoredUCationCation,StoredUCationCationVDW,StoredUCationCationChargeChargeReal;
  REAL StoredUCationCationChargeBondDipoleReal,StoredUCationCationBondDipoleBondDipoleReal;
  REAL StoredUCationCationChargeChargeFourier,StoredUCationCationCoulomb;
  REAL StoredUCationCationChargeBondDipoleFourier,StoredUCationCationBondDipoleBondDipoleFourier;
  REAL StoredUTotal,StoredUTailCorrection;

  REAL UHostPolarizationStored,UAdsorbatePolarizationStored,UCationPolarizationStored;
  REAL UHostBackPolarizationStored,UAdsorbateBackPolarizationStored,UCationBackPolarizationStored;


  // store all energies
  StoredUTotal=UTotal[CurrentSystem];
  StoredUTailCorrection=UTailCorrection[CurrentSystem];

  StoredUHostBond=UHostBond[CurrentSystem];
  StoredUHostUreyBradley=UHostUreyBradley[CurrentSystem];
  StoredUHostBend=UHostBend[CurrentSystem];
  StoredUHostInversionBend=UHostInversionBend[CurrentSystem];
  StoredUHostTorsion=UHostTorsion[CurrentSystem];
  StoredUHostImproperTorsion=UHostImproperTorsion[CurrentSystem];
  StoredUHostBondBond=UHostBondBond[CurrentSystem];
  StoredUHostBendBend=UHostBendBend[CurrentSystem];
  StoredUHostBondBend=UHostBondBend[CurrentSystem];
  StoredUHostBondTorsion=UHostBondTorsion[CurrentSystem];
  StoredUHostBendTorsion=UHostBendTorsion[CurrentSystem];

  StoredUAdsorbateBond=UAdsorbateBond[CurrentSystem];
  StoredUAdsorbateUreyBradley=UAdsorbateUreyBradley[CurrentSystem];
  StoredUAdsorbateBend=UAdsorbateBend[CurrentSystem];
  StoredUAdsorbateInversionBend=UAdsorbateInversionBend[CurrentSystem];
  StoredUAdsorbateTorsion=UAdsorbateTorsion[CurrentSystem];
  StoredUAdsorbateImproperTorsion=UAdsorbateImproperTorsion[CurrentSystem];
  StoredUAdsorbateBondBond=UAdsorbateBondBond[CurrentSystem];
  StoredUAdsorbateBendBend=UAdsorbateBendBend[CurrentSystem];
  StoredUAdsorbateBondBend=UAdsorbateBondBend[CurrentSystem];
  StoredUAdsorbateBondTorsion=UAdsorbateBondTorsion[CurrentSystem];
  StoredUAdsorbateBendTorsion=UAdsorbateBendTorsion[CurrentSystem];
  StoredUAdsorbateIntraVDW=UAdsorbateIntraVDW[CurrentSystem];
  StoredUAdsorbateIntraChargeCharge=UAdsorbateIntraChargeCharge[CurrentSystem];
  StoredUAdsorbateIntraChargeBondDipole=UAdsorbateIntraChargeBondDipole[CurrentSystem];
  StoredUAdsorbateIntraBondDipoleBondDipole=UAdsorbateIntraBondDipoleBondDipole[CurrentSystem];

  StoredUCationBond=UCationBond[CurrentSystem];
  StoredUCationUreyBradley=UCationUreyBradley[CurrentSystem];
  StoredUCationBend=UCationBend[CurrentSystem];
  StoredUCationInversionBend=UCationInversionBend[CurrentSystem];
  StoredUCationTorsion=UCationTorsion[CurrentSystem];
  StoredUCationImproperTorsion=UCationImproperTorsion[CurrentSystem];
  StoredUCationBondBond=UCationBondBond[CurrentSystem];
  StoredUCationBendBend=UCationBendBend[CurrentSystem];
  StoredUCationBondBend=UCationBondBend[CurrentSystem];
  StoredUCationBondTorsion=UCationBondTorsion[CurrentSystem];
  StoredUCationBendTorsion=UCationBendTorsion[CurrentSystem];
  StoredUCationIntraVDW=UCationIntraVDW[CurrentSystem];
  StoredUCationIntraChargeCharge=UCationIntraChargeCharge[CurrentSystem];
  StoredUCationIntraChargeBondDipole=UCationIntraChargeBondDipole[CurrentSystem];
  StoredUCationIntraBondDipoleBondDipole=UCationIntraBondDipoleBondDipole[CurrentSystem];

  StoredUHostHost=UHostHost[CurrentSystem];
  StoredUHostHostVDW=UHostHostVDW[CurrentSystem];
  StoredUHostHostChargeChargeReal=UHostHostChargeChargeReal[CurrentSystem];
  StoredUHostHostChargeChargeFourier=UHostHostChargeChargeFourier[CurrentSystem];
  StoredUHostHostChargeBondDipoleReal=UHostHostChargeBondDipoleReal[CurrentSystem];
  StoredUHostHostChargeBondDipoleFourier=UHostHostChargeBondDipoleFourier[CurrentSystem];
  StoredUHostHostBondDipoleBondDipoleReal=UHostHostBondDipoleBondDipoleReal[CurrentSystem];
  StoredUHostHostBondDipoleBondDipoleFourier=UHostHostBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUHostHostCoulomb=UHostHostCoulomb[CurrentSystem];

  StoredUHostAdsorbate=UHostAdsorbate[CurrentSystem];
  StoredUHostAdsorbateVDW=UHostAdsorbateVDW[CurrentSystem];
  StoredUHostAdsorbateChargeChargeReal=UHostAdsorbateChargeChargeReal[CurrentSystem];
  StoredUHostAdsorbateChargeBondDipoleReal=UHostAdsorbateChargeBondDipoleReal[CurrentSystem];
  StoredUHostAdsorbateBondDipoleBondDipoleReal=UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
  StoredUHostAdsorbateChargeChargeFourier=UHostAdsorbateChargeChargeFourier[CurrentSystem];
  StoredUHostAdsorbateChargeBondDipoleFourier=UHostAdsorbateChargeBondDipoleFourier[CurrentSystem];
  StoredUHostAdsorbateBondDipoleBondDipoleFourier=UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUHostAdsorbateCoulomb=UHostAdsorbateCoulomb[CurrentSystem];

  StoredUHostCation=UHostCation[CurrentSystem];
  StoredUHostCationVDW=UHostCationVDW[CurrentSystem];
  StoredUHostCationChargeChargeReal=UHostCationChargeChargeReal[CurrentSystem];
  StoredUHostCationChargeBondDipoleReal=UHostCationChargeBondDipoleReal[CurrentSystem];
  StoredUHostCationBondDipoleBondDipoleReal=UHostCationBondDipoleBondDipoleReal[CurrentSystem];
  StoredUHostCationChargeChargeFourier=UHostCationChargeChargeFourier[CurrentSystem];
  StoredUHostCationChargeBondDipoleFourier=UHostCationChargeBondDipoleFourier[CurrentSystem];
  StoredUHostCationBondDipoleBondDipoleFourier=UHostCationBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUHostCationCoulomb=UHostCationCoulomb[CurrentSystem];

  StoredUAdsorbateAdsorbate=UAdsorbateAdsorbate[CurrentSystem];
  StoredUAdsorbateAdsorbateVDW=UAdsorbateAdsorbateVDW[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeChargeReal=UAdsorbateAdsorbateChargeChargeReal[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeBondDipoleReal=UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem];
  StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal=UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeChargeFourier=UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeBondDipoleFourier=UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier=UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateCoulomb=UAdsorbateAdsorbateCoulomb[CurrentSystem];

  StoredUAdsorbateCation=UAdsorbateCation[CurrentSystem];
  StoredUAdsorbateCationVDW=UAdsorbateCationVDW[CurrentSystem];
  StoredUAdsorbateCationChargeChargeReal=UAdsorbateCationChargeChargeReal[CurrentSystem];
  StoredUAdsorbateCationChargeBondDipoleReal=UAdsorbateCationChargeBondDipoleReal[CurrentSystem];
  StoredUAdsorbateCationBondDipoleBondDipoleReal=UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem];
  StoredUAdsorbateCationChargeChargeFourier=UAdsorbateCationChargeChargeFourier[CurrentSystem];
  StoredUAdsorbateCationChargeBondDipoleFourier=UAdsorbateCationChargeBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateCationBondDipoleBondDipoleFourier=UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateCationCoulomb=UAdsorbateCationCoulomb[CurrentSystem];

  StoredUCationCation=UCationCation[CurrentSystem];
  StoredUCationCationVDW=UCationCationVDW[CurrentSystem];
  StoredUCationCationChargeChargeReal=UCationCationChargeChargeReal[CurrentSystem];
  StoredUCationCationChargeBondDipoleReal=UCationCationChargeBondDipoleReal[CurrentSystem];
  StoredUCationCationBondDipoleBondDipoleReal=UCationCationBondDipoleBondDipoleReal[CurrentSystem];
  StoredUCationCationChargeChargeFourier=UCationCationChargeChargeFourier[CurrentSystem];
  StoredUCationCationChargeBondDipoleFourier=UCationCationChargeBondDipoleFourier[CurrentSystem];
  StoredUCationCationBondDipoleBondDipoleFourier=UCationCationBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUCationCationCoulomb=UCationCationCoulomb[CurrentSystem];

  UHostPolarizationStored=UHostPolarization[CurrentSystem];
  UAdsorbatePolarizationStored=UAdsorbatePolarization[CurrentSystem];
  UCationPolarizationStored=UCationPolarization[CurrentSystem];

  UHostBackPolarizationStored=UHostBackPolarization[CurrentSystem];
  UAdsorbateBackPolarizationStored=UAdsorbateBackPolarization[CurrentSystem];
  UCationBackPolarizationStored=UCationBackPolarization[CurrentSystem];

  // store the structure-factors for the Ewald-summations
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    SaveCurrentEwaldStructureFactors(0,CurrentSystem);

 // NOTE: CalculateEnergy does not recompute the positions from the center-of-mass and orientation (otherwise the type needs to be switched too).
  CalculateEnergy();

  if(fabs(StoredUTotal-UTotal[CurrentSystem])<1e-4)
  {
    fprintf(stderr, "Energy status system [%d]: okay, %18.10f vs. %18.10f\n",CurrentSystem,StoredUTotal,UTotal[CurrentSystem]);
  }
  else
  {
    for(i=0;i<NumberOfSystems;i++)
    {
      fprintf(stderr, "\n\n");
      fprintf(stderr, "ERROR [energy status]: current energy (%18.10f) and true energy (%18.10f) are different in system %d\n",StoredUTotal,UTotal[CurrentSystem],CurrentSystem);
      if(fabs(StoredUTotal)>1e-8) fprintf(stderr, "UTotal: %18.10f %18.10f\n",StoredUTotal,UTotal[CurrentSystem]);
      if(fabs(StoredUTailCorrection)>1e-8) fprintf(stderr, "UTailCorrection: %18.10f %18.10f\n",StoredUTailCorrection,UTailCorrection[CurrentSystem]);

      if(fabs(StoredUHostBond)>1e-8) fprintf(stderr, "UHostBond: %18.10f %18.10f\n",StoredUHostBond,UHostBond[CurrentSystem]);
      if(fabs(StoredUHostUreyBradley)>1e-8) fprintf(stderr, "UHostUreyBradley: %18.10f %18.10f\n",StoredUHostUreyBradley,UHostUreyBradley[CurrentSystem]);
      if(fabs(StoredUHostBend)>1e-8) fprintf(stderr, "UHostBend: %18.10f %18.10f\n",StoredUHostBend,UHostBend[CurrentSystem]);
      if(fabs(StoredUHostInversionBend)>1e-8) fprintf(stderr, "UHostInversionBend: %18.10f %18.10f\n",StoredUHostInversionBend,UHostInversionBend[CurrentSystem]);
      if(fabs(StoredUHostTorsion)>1e-8) fprintf(stderr, "UHostTorsion: %18.10f %18.10f\n",StoredUHostTorsion,UHostTorsion[CurrentSystem]);
      if(fabs(StoredUHostImproperTorsion)>1e-8) fprintf(stderr, "UHostImproperTorsion: %18.10f %18.10f\n",StoredUHostImproperTorsion,UHostImproperTorsion[CurrentSystem]);
      if(fabs(StoredUHostBondBond)>1e-8) fprintf(stderr, "UHostBondBond: %18.10f %18.10f\n",StoredUHostBondBond,UHostBondBond[CurrentSystem]);
      if(fabs(StoredUHostBendBend)>1e-8) fprintf(stderr, "UHostBendBend: %18.10f %18.10f\n",StoredUHostBendBend,UHostBendBend[CurrentSystem]);
      if(fabs(StoredUHostBondBend)>1e-8) fprintf(stderr, "UHostBondBend: %18.10f %18.10f\n",StoredUHostBondBend,UHostBondBend[CurrentSystem]);
      if(fabs(StoredUHostBondTorsion)>1e-8) fprintf(stderr, "UHostBondTorsion: %18.10f %18.10f\n",StoredUHostBondTorsion,UHostBondTorsion[CurrentSystem]);
      if(fabs(StoredUHostBendTorsion)>1e-8) fprintf(stderr, "UHostBendTorsion: %18.10f %18.10f\n",StoredUHostBendTorsion,UHostBendTorsion[CurrentSystem]);

      if(fabs(StoredUAdsorbateBond)>1e-8) fprintf(stderr, "UAdsorbateBond: %18.10f %18.10f\n",StoredUAdsorbateBond,UAdsorbateBond[CurrentSystem]);
      if(fabs(StoredUAdsorbateUreyBradley)>1e-8) fprintf(stderr, "UAdsorbateUreyBradley: %18.10f %18.10f\n",StoredUAdsorbateUreyBradley,UAdsorbateUreyBradley[CurrentSystem]);
      if(fabs(StoredUAdsorbateBend)>1e-8) fprintf(stderr, "UAdsorbateBend: %18.10f %18.10f\n",StoredUAdsorbateBend,UAdsorbateBend[CurrentSystem]);
      if(fabs(StoredUAdsorbateInversionBend)>1e-8) fprintf(stderr, "UAdsorbateInversionBend: %18.10f %18.10f\n",StoredUAdsorbateInversionBend,UAdsorbateInversionBend[CurrentSystem]);
      if(fabs(StoredUAdsorbateTorsion)>1e-8) fprintf(stderr, "UAdsorbateTorsion: %18.10f %18.10f\n",StoredUAdsorbateTorsion,UAdsorbateTorsion[CurrentSystem]);
      if(fabs(StoredUAdsorbateImproperTorsion)>1e-8) fprintf(stderr, "UAdsorbateImproperTorsion: %18.10f %18.10f\n",StoredUAdsorbateImproperTorsion,UAdsorbateImproperTorsion[CurrentSystem]);
      if(fabs(StoredUAdsorbateBondBond)>1e-8) fprintf(stderr, "UAdsorbateBondBond: %18.10f %18.10f\n",StoredUAdsorbateBondBond,UAdsorbateBondBond[CurrentSystem]);
      if(fabs(StoredUAdsorbateBendBend)>1e-8) fprintf(stderr, "UAdsorbateBendBend: %18.10f %18.10f\n",StoredUAdsorbateBendBend,UAdsorbateBendBend[CurrentSystem]);
      if(fabs(StoredUAdsorbateBondBend)>1e-8) fprintf(stderr, "UAdsorbateBondBend: %18.10f %18.10f\n",StoredUAdsorbateBondBend,UAdsorbateBondBend[CurrentSystem]);
      if(fabs(StoredUAdsorbateBondTorsion)>1e-8) fprintf(stderr, "UAdsorbateBondTorsion: %18.10f %18.10f\n",StoredUAdsorbateBondTorsion,UAdsorbateBondTorsion[CurrentSystem]);
      if(fabs(StoredUAdsorbateBendTorsion)>1e-8) fprintf(stderr, "UAdsorbateBendTorsion: %18.10f %18.10f\n",StoredUAdsorbateBendTorsion,UAdsorbateBendTorsion[CurrentSystem]);
      if(fabs(StoredUAdsorbateIntraVDW)>1e-8) fprintf(stderr, "UAdsorbateIntraVDW: %18.10f %18.10f\n",StoredUAdsorbateIntraVDW,UAdsorbateIntraVDW[CurrentSystem]);
      if(fabs(StoredUAdsorbateIntraChargeCharge)>1e-8) fprintf(stderr, "UAdsorbateIntraChargeCharge: %18.10f %18.10f\n",StoredUAdsorbateIntraChargeCharge,UAdsorbateIntraChargeCharge[CurrentSystem]);
      if(fabs(StoredUAdsorbateIntraChargeBondDipole)>1e-8) fprintf(stderr, "UAdsorbateIntraChargeBondDipole: %18.10f %18.10f\n",StoredUAdsorbateIntraChargeBondDipole,UAdsorbateIntraChargeBondDipole[CurrentSystem]);
      if(fabs(StoredUAdsorbateIntraBondDipoleBondDipole)>1e-8) fprintf(stderr, "UAdsorbateIntraBondDipoleBondDipole: %18.10f %18.10f\n",StoredUAdsorbateIntraBondDipoleBondDipole,UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]);

      if(fabs(StoredUCationBond)>1e-8) fprintf(stderr, "UCationBond: %18.10f %18.10f\n",StoredUCationBond,UCationBond[CurrentSystem]);
      if(fabs(StoredUCationUreyBradley)>1e-8) fprintf(stderr, "UCationUreyBradley: %18.10f %18.10f\n",StoredUCationUreyBradley,UCationUreyBradley[CurrentSystem]);
      if(fabs(StoredUCationBend)>1e-8) fprintf(stderr, "UCationBend: %18.10f %18.10f\n",StoredUCationBend,UCationBend[CurrentSystem]);
      if(fabs(StoredUCationInversionBend)>1e-8) fprintf(stderr, "UCationInversionBend: %18.10f %18.10f\n",StoredUCationInversionBend,UCationInversionBend[CurrentSystem]);
      if(fabs(StoredUCationTorsion)>1e-8) fprintf(stderr, "UCationTorsion: %18.10f %18.10f\n",StoredUCationTorsion,UCationTorsion[CurrentSystem]);
      if(fabs(StoredUCationImproperTorsion)>1e-8) fprintf(stderr, "UCationImproperTorsion: %18.10f %18.10f\n",StoredUCationImproperTorsion,UCationImproperTorsion[CurrentSystem]);
      if(fabs(StoredUCationBondBond)>1e-8) fprintf(stderr, "UCationBondBond: %18.10f %18.10f\n",StoredUCationBondBond,UCationBondBond[CurrentSystem]);
      if(fabs(StoredUCationBendBend)>1e-8) fprintf(stderr, "UCationBendBend: %18.10f %18.10f\n",StoredUCationBendBend,UCationBendBend[CurrentSystem]);
      if(fabs(StoredUCationBondBend)>1e-8) fprintf(stderr, "UCationBondBend: %18.10f %18.10f\n",StoredUCationBondBend,UCationBondBend[CurrentSystem]);
      if(fabs(StoredUCationBondTorsion)>1e-8) fprintf(stderr, "UCationBondTorsion: %18.10f %18.10f\n",StoredUCationBondTorsion,UCationBondTorsion[CurrentSystem]);
      if(fabs(StoredUCationBendTorsion)>1e-8) fprintf(stderr, "UCationBendTorsion: %18.10f %18.10f\n",StoredUCationBendTorsion,UCationBendTorsion[CurrentSystem]);
      if(fabs(StoredUCationIntraVDW)>1e-8) fprintf(stderr, "UCationIntraVDW: %18.10f %18.10f\n",StoredUCationIntraVDW,UCationIntraVDW[CurrentSystem]);
      if(fabs(StoredUCationIntraChargeCharge)>1e-8) fprintf(stderr, "UCationIntraChargeCharge: %18.10f %18.10f\n",StoredUCationIntraChargeCharge,UCationIntraChargeCharge[CurrentSystem]);
      if(fabs(StoredUCationIntraChargeBondDipole)>1e-8) fprintf(stderr, "UCationIntraChargeBondDipole: %18.10f %18.10f\n",StoredUCationIntraChargeBondDipole,UCationIntraChargeBondDipole[CurrentSystem]);
      if(fabs(StoredUCationIntraBondDipoleBondDipole)>1e-8) fprintf(stderr, "UCationIntraBondDipoleBondDipole: %18.10f %18.10f\n",StoredUCationIntraBondDipoleBondDipole,UCationIntraBondDipoleBondDipole[CurrentSystem]);

      if(fabs(StoredUHostHost)>1e-8) fprintf(stderr, "UHostHost: %18.10f %18.10f\n",StoredUHostHost,UHostHost[CurrentSystem]);
      if(fabs(StoredUHostHostVDW)>1e-8) fprintf(stderr, "UHostHostVDW: %18.10f %18.10f\n",StoredUHostHostVDW,UHostHostVDW[CurrentSystem]);
      if(fabs(StoredUHostHostChargeChargeReal)>1e-8) fprintf(stderr, "UHostHostChargeChargeReal: %18.10f %18.10f\n",StoredUHostHostChargeChargeReal,UHostHostChargeChargeReal[CurrentSystem]);
      if(fabs(StoredUHostHostChargeChargeFourier)>1e-8) fprintf(stderr, "UHostHostChargeChargeFourier: %18.10f %18.10f\n",StoredUHostHostChargeChargeFourier,UHostHostChargeChargeFourier[CurrentSystem]);
      if(fabs(StoredUHostHostChargeBondDipoleReal)>1e-8) fprintf(stderr, "UHostHostChargeBondDipoleReal: %18.10f %18.10f\n",StoredUHostHostChargeBondDipoleReal,UHostHostChargeBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUHostHostChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UHostHostChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUHostHostChargeBondDipoleFourier,UHostHostChargeBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUHostHostBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UHostHostBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUHostHostBondDipoleBondDipoleReal,UHostHostBondDipoleBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUHostHostBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UHostHostBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUHostHostBondDipoleBondDipoleFourier,UHostHostBondDipoleBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUHostHostCoulomb)>1e-8) fprintf(stderr, "UHostHostCoulomb: %18.10f %18.10f\n",StoredUHostHostCoulomb,UHostHostCoulomb[CurrentSystem]);

      if(fabs(StoredUHostAdsorbate)>1e-8) fprintf(stderr, "UHostAdsorbate: %18.10f %18.10f\n",StoredUHostAdsorbate,UHostAdsorbate[CurrentSystem]);
      if(fabs(StoredUHostAdsorbateVDW)>1e-8) fprintf(stderr, "UHostAdsorbateVDW: %18.10f %18.10f\n",StoredUHostAdsorbateVDW,UHostAdsorbateVDW[CurrentSystem]);
      if(fabs(StoredUHostAdsorbateChargeChargeReal)>1e-8) fprintf(stderr, "UHostAdsorbateChargeChargeReal: %18.10f %18.10f\n",StoredUHostAdsorbateChargeChargeReal,UHostAdsorbateChargeChargeReal[CurrentSystem]);
      if(fabs(StoredUHostAdsorbateChargeBondDipoleReal)>1e-8) fprintf(stderr, "UHostAdsorbateChargeBondDipoleReal: %18.10f %18.10f\n",StoredUHostAdsorbateChargeBondDipoleReal,UHostAdsorbateChargeBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUHostAdsorbateBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UHostAdsorbateBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUHostAdsorbateBondDipoleBondDipoleReal,UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUHostAdsorbateChargeChargeFourier)>1e-8) fprintf(stderr, "UHostAdsorbateChargeChargeFourier: %18.10f %18.10f\n",StoredUHostAdsorbateChargeChargeFourier,UHostAdsorbateChargeChargeFourier[CurrentSystem]);
      if(fabs(StoredUHostAdsorbateChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UHostAdsorbateChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUHostAdsorbateChargeBondDipoleFourier,UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUHostAdsorbateBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UHostAdsorbateBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUHostAdsorbateBondDipoleBondDipoleFourier,UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUHostAdsorbateCoulomb)>1e-8) fprintf(stderr, "UHostAdsorbateCoulomb: %18.10f %18.10f\n",StoredUHostAdsorbateCoulomb,UHostAdsorbateCoulomb[CurrentSystem]);

      if(fabs(StoredUHostCation)>1e-8) fprintf(stderr, "UHostCation: %18.10f %18.10f\n",StoredUHostCation,UHostCation[CurrentSystem]);
      if(fabs(StoredUHostCationVDW)>1e-8) fprintf(stderr, "UHostCationVDW: %18.10f %18.10f\n",StoredUHostCationVDW,UHostCationVDW[CurrentSystem]);
      if(fabs(StoredUHostCationChargeChargeReal)>1e-8) fprintf(stderr, "UHostCationChargeChargeReal: %18.10f %18.10f\n",StoredUHostCationChargeChargeReal,UHostCationChargeChargeReal[CurrentSystem]);
      if(fabs(StoredUHostCationChargeBondDipoleReal)>1e-8) fprintf(stderr, "UHostCationChargeBondDipoleReal: %18.10f %18.10f\n",StoredUHostCationChargeBondDipoleReal,UHostCationChargeBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUHostCationBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UHostCationBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUHostCationBondDipoleBondDipoleReal,UHostCationBondDipoleBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUHostCationChargeChargeFourier)>1e-8) fprintf(stderr, "UHostCationChargeChargeFourier: %18.10f %18.10f\n",StoredUHostCationChargeChargeFourier,UHostCationChargeChargeFourier[CurrentSystem]);
      if(fabs(StoredUHostCationChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UHostCationChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUHostCationChargeBondDipoleFourier,UHostCationChargeBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUHostCationBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UHostCationBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUHostCationBondDipoleBondDipoleFourier,UHostCationBondDipoleBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUHostCationCoulomb)>1e-8) fprintf(stderr, "UHostCationCoulomb: %18.10f %18.10f\n",StoredUHostCationCoulomb,UHostCationCoulomb[CurrentSystem]);

      if(fabs(StoredUAdsorbateAdsorbate)>1e-8) fprintf(stderr, "UAdsorbateAdsorbate: %18.10f %18.10f\n",StoredUAdsorbateAdsorbate,UAdsorbateAdsorbate[CurrentSystem]);
      if(fabs(StoredUAdsorbateAdsorbateVDW)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateVDW: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateVDW,UAdsorbateAdsorbateVDW[CurrentSystem]);
      if(fabs(StoredUAdsorbateAdsorbateChargeChargeReal)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateChargeChargeReal: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateChargeChargeReal,UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]);
      if(fabs(StoredUAdsorbateAdsorbateChargeBondDipoleReal)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateChargeBondDipoleReal: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateChargeBondDipoleReal,UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal,UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUAdsorbateAdsorbateChargeChargeFourier)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateChargeChargeFourier: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateChargeChargeFourier,UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]);
      if(fabs(StoredUAdsorbateAdsorbateChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateChargeBondDipoleFourier,UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier,UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUAdsorbateAdsorbateCoulomb)>1e-8) fprintf(stderr, "UAdsorbateAdsorbateCoulomb: %18.10f %18.10f\n",StoredUAdsorbateAdsorbateCoulomb,UAdsorbateAdsorbateCoulomb[CurrentSystem]);

      if(fabs(StoredUAdsorbateCation)>1e-8) fprintf(stderr, "UAdsorbateCation: %18.10f %18.10f\n",StoredUAdsorbateCation,UAdsorbateCation[CurrentSystem]);
      if(fabs(StoredUAdsorbateCationVDW)>1e-8) fprintf(stderr, "UAdsorbateCationVDW: %18.10f %18.10f\n",StoredUAdsorbateCationVDW,UAdsorbateCationVDW[CurrentSystem]);
      if(fabs(StoredUAdsorbateCationChargeChargeReal)>1e-8) fprintf(stderr, "UAdsorbateCationChargeChargeReal: %18.10f %18.10f\n",StoredUAdsorbateCationChargeChargeReal,UAdsorbateCationChargeChargeReal[CurrentSystem]);
      if(fabs(StoredUAdsorbateCationChargeBondDipoleReal)>1e-8) fprintf(stderr, "UAdsorbateCationChargeBondDipoleReal: %18.10f %18.10f\n",StoredUAdsorbateCationChargeBondDipoleReal,UAdsorbateCationChargeBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUAdsorbateCationBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UAdsorbateCationBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUAdsorbateCationBondDipoleBondDipoleReal,UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUAdsorbateCationChargeChargeFourier)>1e-8) fprintf(stderr, "UAdsorbateCationChargeChargeFourier: %18.10f %18.10f\n",StoredUAdsorbateCationChargeChargeFourier,UAdsorbateCationChargeChargeFourier[CurrentSystem]);
      if(fabs(StoredUAdsorbateCationChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UAdsorbateCationChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUAdsorbateCationChargeBondDipoleFourier,UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUAdsorbateCationBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UAdsorbateCationBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUAdsorbateCationBondDipoleBondDipoleFourier,UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUAdsorbateCationCoulomb)>1e-8) fprintf(stderr, "UAdsorbateCationCoulomb: %18.10f %18.10f\n",StoredUAdsorbateCationCoulomb,UAdsorbateCationCoulomb[CurrentSystem]);

      if(fabs(StoredUCationCation)>1e-8) fprintf(stderr, "UCationCation: %18.10f %18.10f\n",StoredUCationCation,UCationCation[CurrentSystem]);
      if(fabs(StoredUCationCationVDW)>1e-8) fprintf(stderr, "UCationCationVDW: %18.10f %18.10f\n",StoredUCationCationVDW,UCationCationVDW[CurrentSystem]);
      if(fabs(StoredUCationCationChargeChargeReal)>1e-8) fprintf(stderr, "UCationCationChargeChargeReal: %18.10f %18.10f\n",StoredUCationCationChargeChargeReal,UCationCationChargeChargeReal[CurrentSystem]);
      if(fabs(StoredUCationCationChargeBondDipoleReal)>1e-8) fprintf(stderr, "UCationCationChargeBondDipoleReal: %18.10f %18.10f\n",StoredUCationCationChargeBondDipoleReal,UCationCationChargeBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUCationCationBondDipoleBondDipoleReal)>1e-8) fprintf(stderr, "UCationCationBondDipoleBondDipoleReal: %18.10f %18.10f\n",StoredUCationCationBondDipoleBondDipoleReal,UCationCationBondDipoleBondDipoleReal[CurrentSystem]);
      if(fabs(StoredUCationCationChargeChargeFourier)>1e-8) fprintf(stderr, "UCationCationChargeChargeFourier: %18.10f %18.10f\n",StoredUCationCationChargeChargeFourier,UCationCationChargeChargeFourier[CurrentSystem]);
      if(fabs(StoredUCationCationChargeBondDipoleFourier)>1e-8) fprintf(stderr, "UCationCationChargeBondDipoleFourier: %18.10f %18.10f\n",StoredUCationCationChargeBondDipoleFourier,UCationCationChargeBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUCationCationBondDipoleBondDipoleFourier)>1e-8) fprintf(stderr, "UCationCationBondDipoleBondDipoleFourier: %18.10f %18.10f\n",StoredUCationCationBondDipoleBondDipoleFourier,UCationCationBondDipoleBondDipoleFourier[CurrentSystem]);
      if(fabs(StoredUCationCationCoulomb)>1e-8) fprintf(stderr, "UCationCationCoulomb: %18.10f %18.10f\n",StoredUCationCationCoulomb,UCationCationCoulomb[CurrentSystem]);

      if(fabs(UHostPolarizationStored)>1e-8) fprintf(stderr, "UHostPolarization: %18.10f %18.10f\n",UHostPolarizationStored,UHostPolarization[CurrentSystem]);
      if(fabs(UAdsorbatePolarizationStored)>1e-8) fprintf(stderr, "UAdsorbatePolarization: %18.10f %18.10f\n",UAdsorbatePolarizationStored,UAdsorbatePolarization[CurrentSystem]);
      if(fabs(UCationPolarizationStored)>1e-8) fprintf(stderr, "UCationPolarization: %18.10f %18.10f\n",UCationPolarizationStored,UCationPolarization[CurrentSystem]);

      if(fabs(UHostBackPolarizationStored)>1e-8) fprintf(stderr, "UHostBackPolarization: %18.10f %18.10f\n",UHostBackPolarizationStored,UHostBackPolarization[CurrentSystem]);
      if(fabs(UAdsorbateBackPolarizationStored)>1e-8) fprintf(stderr, "UAdsorbateBackPolarization: %18.10f %18.10f\n",UAdsorbateBackPolarizationStored,UAdsorbateBackPolarization[CurrentSystem]);
      if(fabs(UCationBackPolarizationStored)>1e-8) fprintf(stderr, "UCationBackPolarization: %18.10f %18.10f\n",UCationBackPolarizationStored,UCationBackPolarization[CurrentSystem]);
      fclose(OutputFilePtr[i]);
    }
    exit(0);
  }

  UHostBond[CurrentSystem]=StoredUHostBond;
  UHostUreyBradley[CurrentSystem]=StoredUHostUreyBradley;
  UHostBend[CurrentSystem]=StoredUHostBend;
  UHostInversionBend[CurrentSystem]=StoredUHostInversionBend;
  UHostTorsion[CurrentSystem]=StoredUHostTorsion;
  UHostImproperTorsion[CurrentSystem]=StoredUHostImproperTorsion;
  UHostBondBond[CurrentSystem]=StoredUHostBondBond;
  UHostBendBend[CurrentSystem]=StoredUHostBendBend;
  UHostBondBend[CurrentSystem]=StoredUHostBondBend;
  UHostBondTorsion[CurrentSystem]=StoredUHostBondTorsion;
  UHostBendTorsion[CurrentSystem]=StoredUHostBendTorsion;

  UAdsorbateBond[CurrentSystem]=StoredUAdsorbateBond;
  UAdsorbateUreyBradley[CurrentSystem]=StoredUAdsorbateUreyBradley;
  UAdsorbateBend[CurrentSystem]=StoredUAdsorbateBend;
  UAdsorbateInversionBend[CurrentSystem]=StoredUAdsorbateInversionBend;
  UAdsorbateTorsion[CurrentSystem]=StoredUAdsorbateTorsion;
  UAdsorbateImproperTorsion[CurrentSystem]=StoredUAdsorbateImproperTorsion;
  UAdsorbateBondBond[CurrentSystem]=StoredUAdsorbateBondBond;
  UAdsorbateBendBend[CurrentSystem]=StoredUAdsorbateBendBend;
  UAdsorbateBondTorsion[CurrentSystem]=StoredUAdsorbateBondTorsion;
  UAdsorbateBondBend[CurrentSystem]=StoredUAdsorbateBondBend;
  UAdsorbateBendTorsion[CurrentSystem]=StoredUAdsorbateBendTorsion;
  UAdsorbateIntraVDW[CurrentSystem]=StoredUAdsorbateIntraVDW;
  UAdsorbateIntraChargeCharge[CurrentSystem]=StoredUAdsorbateIntraChargeCharge;
  UAdsorbateIntraChargeBondDipole[CurrentSystem]=StoredUAdsorbateIntraChargeBondDipole;
  UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]=StoredUAdsorbateIntraBondDipoleBondDipole;

  UCationBond[CurrentSystem]=StoredUCationBond;
  UCationUreyBradley[CurrentSystem]=StoredUCationUreyBradley;
  UCationBend[CurrentSystem]=StoredUCationBend;
  UCationInversionBend[CurrentSystem]=StoredUCationInversionBend;
  UCationTorsion[CurrentSystem]=StoredUCationTorsion;
  UCationImproperTorsion[CurrentSystem]=StoredUCationImproperTorsion;
  UCationBondBond[CurrentSystem]=StoredUCationBondBond;
  UCationBendBend[CurrentSystem]=StoredUCationBendBend;
  UCationBondBend[CurrentSystem]=StoredUCationBondBend;
  UCationBondTorsion[CurrentSystem]=StoredUCationBondTorsion;
  UCationBendTorsion[CurrentSystem]=StoredUCationBendTorsion;
  UCationIntraVDW[CurrentSystem]=StoredUCationIntraVDW;
  UCationIntraChargeCharge[CurrentSystem]=StoredUCationIntraChargeCharge;
  UCationIntraChargeBondDipole[CurrentSystem]=StoredUCationIntraChargeBondDipole;
  UCationIntraBondDipoleBondDipole[CurrentSystem]=StoredUCationIntraBondDipoleBondDipole;

  UHostHost[CurrentSystem]=StoredUHostHost;
  UHostHostVDW[CurrentSystem]=StoredUHostHostVDW;
  UHostHostChargeChargeReal[CurrentSystem]=StoredUHostHostChargeChargeReal;
  UHostHostChargeBondDipoleReal[CurrentSystem]=StoredUHostHostChargeBondDipoleReal;
  UHostHostBondDipoleBondDipoleReal[CurrentSystem]=StoredUHostHostBondDipoleBondDipoleReal;
  UHostHostChargeChargeFourier[CurrentSystem]=StoredUHostHostChargeChargeFourier;
  UHostHostChargeBondDipoleFourier[CurrentSystem]=StoredUHostHostChargeBondDipoleFourier;
  UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=StoredUHostHostBondDipoleBondDipoleFourier;
  UHostHostCoulomb[CurrentSystem]=StoredUHostHostCoulomb;

  UHostAdsorbate[CurrentSystem]=StoredUHostAdsorbate;
  UHostAdsorbateVDW[CurrentSystem]=StoredUHostAdsorbateVDW;
  UHostAdsorbateChargeChargeReal[CurrentSystem]=StoredUHostAdsorbateChargeChargeReal;
  UHostAdsorbateChargeBondDipoleReal[CurrentSystem]=StoredUHostAdsorbateChargeBondDipoleReal;
  UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=StoredUHostAdsorbateBondDipoleBondDipoleReal;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=StoredUHostAdsorbateChargeChargeFourier;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=StoredUHostAdsorbateChargeBondDipoleFourier;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=StoredUHostAdsorbateBondDipoleBondDipoleFourier;
  UHostAdsorbateCoulomb[CurrentSystem]=StoredUHostAdsorbateCoulomb;

  UHostCation[CurrentSystem]=StoredUHostCation;
  UHostCationVDW[CurrentSystem]=StoredUHostCationVDW;
  UHostCationChargeChargeReal[CurrentSystem]=StoredUHostCationChargeChargeReal;
  UHostCationChargeBondDipoleReal[CurrentSystem]=StoredUHostCationChargeBondDipoleReal;
  UHostCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUHostCationBondDipoleBondDipoleReal;
  UHostCationChargeChargeFourier[CurrentSystem]=StoredUHostCationChargeChargeFourier;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=StoredUHostCationChargeBondDipoleFourier;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUHostCationBondDipoleBondDipoleFourier;
  UHostCationCoulomb[CurrentSystem]=StoredUHostCationCoulomb;

  UAdsorbateAdsorbate[CurrentSystem]=StoredUAdsorbateAdsorbate;
  UAdsorbateAdsorbateVDW[CurrentSystem]=StoredUAdsorbateAdsorbateVDW;
  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=StoredUAdsorbateAdsorbateChargeChargeReal;
  UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=StoredUAdsorbateAdsorbateChargeBondDipoleReal;
  UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal;
  UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=StoredUAdsorbateAdsorbateChargeChargeFourier;
  UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=StoredUAdsorbateAdsorbateChargeBondDipoleFourier;
  UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier;
  UAdsorbateAdsorbateCoulomb[CurrentSystem]=StoredUAdsorbateAdsorbateCoulomb;

  UAdsorbateCation[CurrentSystem]=StoredUAdsorbateCation;
  UAdsorbateCationVDW[CurrentSystem]=StoredUAdsorbateCationVDW;
  UAdsorbateCationChargeChargeReal[CurrentSystem]=StoredUAdsorbateCationChargeChargeReal;
  UAdsorbateCationChargeBondDipoleReal[CurrentSystem]=StoredUAdsorbateCationChargeBondDipoleReal;
  UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUAdsorbateCationBondDipoleBondDipoleReal;
  UAdsorbateCationChargeChargeFourier[CurrentSystem]=StoredUAdsorbateCationChargeChargeFourier;
  UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=StoredUAdsorbateCationChargeBondDipoleFourier;
  UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUAdsorbateCationBondDipoleBondDipoleFourier;
  UAdsorbateCationCoulomb[CurrentSystem]=StoredUAdsorbateCationCoulomb;

  UCationCation[CurrentSystem]=StoredUCationCation;
  UCationCationVDW[CurrentSystem]=StoredUCationCationVDW;
  UCationCationChargeChargeReal[CurrentSystem]=StoredUCationCationChargeChargeReal;
  UCationCationChargeBondDipoleReal[CurrentSystem]=StoredUCationCationChargeBondDipoleReal;
  UCationCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUCationCationBondDipoleBondDipoleReal;
  UCationCationChargeChargeFourier[CurrentSystem]=StoredUCationCationChargeChargeFourier;
  UCationCationChargeBondDipoleFourier[CurrentSystem]=StoredUCationCationChargeBondDipoleFourier;
  UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUCationCationBondDipoleBondDipoleFourier;
  UCationCationCoulomb[CurrentSystem]=StoredUCationCationCoulomb;

  UHostPolarization[CurrentSystem]=UHostPolarizationStored;
  UAdsorbatePolarization[CurrentSystem]=UHostPolarizationStored;
  UCationPolarization[CurrentSystem]=UCationPolarizationStored;

  UHostBackPolarization[CurrentSystem]=UHostBackPolarizationStored;
  UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationStored;
  UCationBackPolarization[CurrentSystem]=UCationBackPolarizationStored;

  UTotal[CurrentSystem]=StoredUTotal;
  UTailCorrection[CurrentSystem]=StoredUTailCorrection;

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
