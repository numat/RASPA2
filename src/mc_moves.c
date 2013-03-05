/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'mc_moves.c' is part of RASPA.

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "potentials.h"
#include "molecule.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "framework.h"
#include "simulation.h"
#include "cbmc.h"
#include "ewald.h"
#include "utils.h"
#include "mc_moves.h"
#include "inter_energy.h"
#include "inter_force.h"
#include "integration.h"
#include "statistics.h"
#include "movies.h"
#include "output.h"
#include "recrossing.h"
#include "grids.h"
#include "thermo_baro_stats.h"
#include "spacegroup.h"
#include "rigid.h"
#include "numerical.h"
#include "internal_energy.h"
#include "sample.h"

// this module contains the Monte Carlo moves
//  1) molecule translation-change
//  2) molecule random translation
//  3) molecule random rotation
//  4) molecule insertion
//  5) molecule deletion
//  6) molcule full regrow
//  7) molecule partial regrow
//  8) molecule regrow in plane
//  9) molecule identity change
// 10) system: parallell tempering
// 11) system: volume change
// 12) system: Gibbs volume change
// 13) molecule: Gibbs particle swap
// 14) framework: atom-translation

// It is assumed that mc-moves operate on at most two systems, e.g. parallel replica and Gibbs-volume move.
// This implies that some temperary storage is created for at most two systems in e.g. the Ewald part.

// Notes
// =============================================================================================
// *) Any mc-move that changes the box-shape/volume and uses the routine 'EwaldFourierForce' [e.g. CalculateForce()]
//    needs to call 'PrecomputeTotalEwaldContributions()' on acceptance to recompute the stored wavevectors.
//    This applies e.g. to MC/MD hybrid moves.
// *) Any mc-move that changes the box-shape/volume should call 'SetupKVectors()' after the new box has been computed.
//    The current KVectors are saved using 'SaveCurrentKVectors()' and restored using 'RetrieveStoredKVectors()' 
//    if the move is not accepted.
// *) MC-moves store the current structure factors using 'SaveCurrentEwaldStructureFactors()'. If the move is not accepted the
//    structure factors are restored using 'RetrieveStoredEwaldStructureFactors'.


// the trial positions for translation and rotation
VECTOR **TrialPosition;
VECTOR **TrialAnisotropicPosition;
VECTOR **ElectricFieldAtTrialPosition;
VECTOR **ReferenceElectricFieldAtTrialPosition;
VECTOR **InducedDipoleAtTrialPosition;
static VECTOR *cord;

// delta energies
static REAL UDeltaPolarization;
static REAL *UHostPolarizationNew;
static REAL *UAdsorbatePolarizationNew;
static REAL *UCationPolarizationNew;
static REAL *UPolarizationNew;

static REAL *UHostBackPolarizationNew;
static REAL *UAdsorbateBackPolarizationNew;
static REAL *UCationBackPolarizationNew;
static REAL *UBackPolarizationNew;

REAL *UHostVDWDelta;
REAL *UHostChargeChargeRealDelta;
REAL *UHostChargeBondDipoleRealDelta;
REAL *UHostBondDipoleBondDipoleRealDelta;
REAL *UHostFourierDelta;

REAL *UCationVDWDelta;
REAL *UAdsorbateVDWDelta;

REAL *UCationChargeChargeRealDelta;
REAL *UAdsorbateChargeChargeRealDelta;

REAL *UCationChargeBondDipoleRealDelta;
REAL *UAdsorbateChargeBondDipoleRealDelta;

REAL *UCationBondDipoleBondDipoleRealDelta;
REAL *UAdsorbateBondDipoleBondDipoleRealDelta;

REAL *UCationFourierDelta;
REAL *UAdsorbateFourierDelta;


// the acceptance-ratio for translation
REAL TargetAccRatioTranslation;
REAL TargetAccRatioLambdaChange;

static VECTOR **TranslationAttempts;
static VECTOR **TranslationAccepted;
static VECTOR **TotalTranslationAttempts;
static VECTOR **TotalTranslationAccepted;
VECTOR **MaximumTranslation;

static REAL **RandomTranslationAttempts;
static REAL **RandomTranslationAccepted;

static REAL **RotationAttempts;
static REAL **RotationAccepted;
REAL **MaximumRotation;

static REAL **SwapAddAttempts;
static REAL (**SwapAddAccepted)[2];

static REAL **SwapRemoveAttempts;
static REAL (**SwapRemoveAccepted)[2];

static REAL (**CFSwapLambdaAttempts)[3];
static REAL (**CFSwapLambdaAccepted)[3];
static REAL **TotalCFSwapLambdaAttempts;
static REAL **TotalCFSwapLambdaAccepted;
REAL **MaximumCFLambdaChange;

static REAL (**CFCBSwapLambdaAttempts)[3];
static REAL (**CFCBSwapLambdaAccepted)[3];
static REAL **TotalCFCBSwapLambdaAttempts;
static REAL **TotalCFCBSwapLambdaAccepted;
REAL **MaximumCFCBLambdaChange;

static REAL **ReinsertionAttempts;
static REAL (**ReinsertionAccepted)[2];

static REAL **PartialReinsertionAttempts;
static REAL (**PartialReinsertionAccepted)[2];

static VECTOR **TranslationInPlaneAttempts;
static VECTOR **TranslationInPlaneAccepted;
static VECTOR **TotalTranslationInPlaneAttempts;
static VECTOR **TotalTranslationInPlaneAccepted;
VECTOR **MaximumTranslationInPlane;
static REAL **ReinsertionInPlaneAttempts;
static REAL (**ReinsertionInPlaneAccepted)[2];

static VECTOR **TranslationInPlaceAttempts;
static VECTOR **TranslationInPlaceAccepted;
static VECTOR **TotalTranslationInPlaceAttempts;
static VECTOR **TotalTranslationInPlaceAccepted;
static REAL **ReinsertionInPlaceAttempts;
static REAL (**ReinsertionInPlaceAccepted)[2];

static REAL ***IdentityChangeAttempts;
static REAL (***IdentityChangeAccepted)[2];

static REAL **ParallelTemperingAttempts;
static REAL **ParallelTemperingAccepted;

static REAL **HyperParallelTemperingAttempts;
static REAL **HyperParallelTemperingAccepted;

static REAL **ParallelMolFractionAttempts;
static REAL **ParallelMolFractionAccepted;
int ParallelMolFractionComponentA;
int ParallelMolFractionComponentB;

static REAL *ChiralInversionAttempts;
static REAL *ChiralInversionAccepted;

static REAL *VolumeChangeAttempts;
static REAL *VolumeChangeAccepted;
static REAL *TotalVolumeChangeAttempts;
static REAL *TotalVolumeChangeAccepted;
REAL *MaximumVolumeChange;
REAL TargetAccRatioVolumeChange;

static REAL_MATRIX3x3 *BoxShapeChangeAttempts;
static REAL_MATRIX3x3 *BoxShapeChangeAccepted;
REAL_MATRIX3x3 *MaximumBoxShapeChange;
REAL TargetAccRatioBoxShapeChange;

static REAL *GibbsVolumeChangeAttempts;
static REAL *GibbsVolumeChangeAccepted;
REAL *MaximumGibbsVolumeChange;
REAL TargetAccRatioGibbsVolumeChange;

static REAL **GibbsSwapAttempts;
static REAL **GibbsSwapAccepted;

static REAL **FrameworkChangeAttempts;
static REAL **FrameworkChangeAccepted;
static REAL **TotalFrameworkChangeAttempts;
static REAL **TotalFrameworkChangeAccepted;
REAL **FrameworkMaximumTranslation;

static VECTOR **FrameworkShiftAttempts;
static VECTOR **FrameworkShiftAccepted;
static VECTOR **TotalFrameworkShiftAttempts;
static VECTOR **TotalFrameworkShiftAccepted;
VECTOR **FrameworkMaximumShiftTranslation;


static REAL ***GibbsIdentityChangeAttempts;
static REAL (***GibbsIdentityChangeAccepted)[2];

int NumberOfHybridNVESteps;
static REAL *HybridNVEDrift;
static REAL *HybridNVEDriftCount;

static REAL *HybridNVEAttempts;
static REAL *HybridNVEAccepted;

static REAL *HybridNVEStartTemperature;
static REAL *HybridNVEStartTranslationalTemperature;
static REAL *HybridNVEStartRotationalTemperature;
static REAL *HybridNVEStartTemperatureFramework;
static REAL *HybridNVEStartTemperatureAdsorbate;
static REAL *HybridNVEStartTemperatureCation;
static REAL *HybridNVEStartTemperatureCount;
static REAL *HybridNVEStartTemperatureTranslationCount;
static REAL *HybridNVEStartTemperatureRotationCount;
static REAL *HybridNVEStartTemperatureFrameworkCount;
static REAL *HybridNVEStartTemperatureAdsorbateCount;
static REAL *HybridNVEStartTemperatureCationCount;

static REAL *HybridNVEEndTemperature;
static REAL *HybridNVEEndTranslationalTemperature;
static REAL *HybridNVEEndRotationalTemperature;
static REAL *HybridNVEEndTemperatureFramework;
static REAL *HybridNVEEndTemperatureAdsorbate;
static REAL *HybridNVEEndTemperatureCation;
static REAL *HybridNVEEndTemperatureCount;
static REAL *HybridNVEEndTemperatureTranslationCount;
static REAL *HybridNVEEndTemperatureRotationCount;
static REAL *HybridNVEEndTemperatureFrameworkCount;
static REAL *HybridNVEEndTemperatureAdsorbateCount;
static REAL *HybridNVEEndTemperatureCationCount;

int NumberOfHybridNPHSteps;
static REAL *HybridNPHDrift;
static REAL *HybridNPHDriftCount;

static REAL *HybridNPHAttempts;
static REAL *HybridNPHAccepted;

static REAL *HybridNPHStartTemperature;
static REAL *HybridNPHStartCellTemperature;
static REAL *HybridNPHStartTranslationalTemperature;
static REAL *HybridNPHStartRotationalTemperature;
static REAL *HybridNPHStartTemperatureFramework;
static REAL *HybridNPHStartTemperatureAdsorbate;
static REAL *HybridNPHStartTemperatureCation;
static REAL *HybridNPHStartTemperatureCount;
static REAL *HybridNPHStartTemperatureTranslationCount;
static REAL *HybridNPHStartTemperatureRotationCount;
static REAL *HybridNPHStartTemperatureFrameworkCount;
static REAL *HybridNPHStartTemperatureAdsorbateCount;
static REAL *HybridNPHStartTemperatureCationCount;

static REAL *HybridNPHEndTemperature;
static REAL *HybridNPHEndCellTemperature;
static REAL *HybridNPHEndTranslationalTemperature;
static REAL *HybridNPHEndRotationalTemperature;
static REAL *HybridNPHEndTemperatureFramework;
static REAL *HybridNPHEndTemperatureAdsorbate;
static REAL *HybridNPHEndTemperatureCation;
static REAL *HybridNPHEndTemperatureCount;
static REAL *HybridNPHEndTemperatureTranslationCount;
static REAL *HybridNPHEndTemperatureRotationCount;
static REAL *HybridNPHEndTemperatureFrameworkCount;
static REAL *HybridNPHEndTemperatureAdsorbateCount;
static REAL *HybridNPHEndTemperatureCationCount;

int NumberOfHybridNPHPRSteps;
static REAL *HybridNPHPRDrift;
static REAL *HybridNPHPRDriftCount;

static REAL *HybridNPHPRAttempts;
static REAL *HybridNPHPRAccepted;

static REAL *HybridNPHPRStartTemperature;
static REAL *HybridNPHPRStartCellTemperature;
static REAL *HybridNPHPRStartTranslationalTemperature;
static REAL *HybridNPHPRStartRotationalTemperature;
static REAL *HybridNPHPRStartTemperatureFramework;
static REAL *HybridNPHPRStartTemperatureAdsorbate;
static REAL *HybridNPHPRStartTemperatureCation;
static REAL *HybridNPHPRStartTemperatureCount;
static REAL *HybridNPHPRStartTemperatureTranslationCount;
static REAL *HybridNPHPRStartTemperatureRotationCount;
static REAL *HybridNPHPRStartTemperatureFrameworkCount;
static REAL *HybridNPHPRStartTemperatureAdsorbateCount;
static REAL *HybridNPHPRStartTemperatureCationCount;

static REAL *HybridNPHPREndTemperature;
static REAL *HybridNPHPREndCellTemperature;
static REAL *HybridNPHPREndTranslationalTemperature;
static REAL *HybridNPHPREndRotationalTemperature;
static REAL *HybridNPHPREndTemperatureFramework;
static REAL *HybridNPHPREndTemperatureAdsorbate;
static REAL *HybridNPHPREndTemperatureCation;
static REAL *HybridNPHPREndTemperatureCount;
static REAL *HybridNPHPREndTemperatureTranslationCount;
static REAL *HybridNPHPREndTemperatureRotationCount;
static REAL *HybridNPHPREndTemperatureFrameworkCount;
static REAL *HybridNPHPREndTemperatureAdsorbateCount;
static REAL *HybridNPHPREndTemperatureCationCount;

//----------------------------------------------------------------------------------------
// CFC-RXMC Parameters
//----------------------------------------------------------------------------------------

REAL *MaximumCFCRXMCLambdaChange;
REAL TargetAccRatioCFCRXMCLambdaChange;

//----------------------------------------------------------------------------------------

int CFWangLandauEvery;
static REAL ***CFLambdaHistogram;

void CheckEnergy(void);


void InitializeMCMovesStatisticsAllSystems(void)
{
  int i,j,k;

  for(j=0;j<NumberOfSystems;j++)
  {
    for(i=0;i<NumberOfComponents;i++)
    {
      TranslationAttempts[j][i].x=0.0;
      TranslationAttempts[j][i].y=0.0;
      TranslationAttempts[j][i].z=0.0;
      TranslationAccepted[j][i].x=0.0;
      TranslationAccepted[j][i].y=0.0;
      TranslationAccepted[j][i].z=0.0;
      TotalTranslationAttempts[j][i].x=0.0;
      TotalTranslationAttempts[j][i].y=0.0;
      TotalTranslationAttempts[j][i].z=0.0;
      TotalTranslationAccepted[j][i].x=0.0;
      TotalTranslationAccepted[j][i].y=0.0;
      TotalTranslationAccepted[j][i].z=0.0;

      RandomTranslationAttempts[j][i]=0.0;
      RandomTranslationAccepted[j][i]=0.0;

      RotationAttempts[j][i]=0.0;
      RotationAccepted[j][i]=0.0;;

      ReinsertionAttempts[j][i]=0.0;
      ReinsertionAccepted[j][i][0]=0.0;
      ReinsertionAccepted[j][i][1]=0.0;

      SwapAddAttempts[j][i]=0.0;
      SwapAddAccepted[j][i][0]=0.0;
      SwapAddAccepted[j][i][1]=0.0;

      SwapRemoveAttempts[j][i]=0.0;
      SwapRemoveAccepted[j][i][0]=0.0;
      SwapRemoveAccepted[j][i][1]=0.0;

      CFSwapLambdaAttempts[j][i][0]=0.0;
      CFSwapLambdaAttempts[j][i][1]=0.0;
      CFSwapLambdaAttempts[j][i][2]=0.0;
      CFSwapLambdaAccepted[j][i][0]=0.0;
      CFSwapLambdaAccepted[j][i][1]=0.0;
      CFSwapLambdaAccepted[j][i][2]=0.0;
      TotalCFSwapLambdaAttempts[j][i]=0.0;
      TotalCFSwapLambdaAccepted[j][i]=0.0;

      CFCBSwapLambdaAttempts[j][i][0]=0.0;
      CFCBSwapLambdaAttempts[j][i][1]=0.0;
      CFCBSwapLambdaAttempts[j][i][2]=0.0;
      CFCBSwapLambdaAccepted[j][i][0]=0.0;
      CFCBSwapLambdaAccepted[j][i][1]=0.0;
      CFCBSwapLambdaAccepted[j][i][2]=0.0;
      TotalCFCBSwapLambdaAttempts[j][i]=0.0;
      TotalCFCBSwapLambdaAccepted[j][i]=0.0;

      PartialReinsertionAttempts[j][i]=0.0;
      PartialReinsertionAccepted[j][i][0]=0.0;
      PartialReinsertionAccepted[j][i][1]=0.0;

      TranslationInPlaneAttempts[j][i].x=0.0;
      TranslationInPlaneAttempts[j][i].y=0.0;
      TranslationInPlaneAttempts[j][i].z=0.0;
      TranslationInPlaneAccepted[j][i].x=0.0;
      TranslationInPlaneAccepted[j][i].y=0.0;
      TranslationInPlaneAccepted[j][i].z=0.0;
      TotalTranslationInPlaneAttempts[j][i].x=0.0;
      TotalTranslationInPlaneAttempts[j][i].y=0.0;
      TotalTranslationInPlaneAttempts[j][i].z=0.0;
      TotalTranslationInPlaneAccepted[j][i].x=0.0;
      TotalTranslationInPlaneAccepted[j][i].y=0.0;
      TotalTranslationInPlaneAccepted[j][i].z=0.0;
      MaximumTranslationInPlane[j][i].x=0.0;
      MaximumTranslationInPlane[j][i].y=0.0;
      MaximumTranslationInPlane[j][i].z=0.0;
      ReinsertionInPlaneAttempts[j][i]=0.0;
      ReinsertionInPlaneAccepted[j][i][0]=0.0;
      ReinsertionInPlaneAccepted[j][i][1]=0.0;

      TranslationInPlaceAttempts[j][i].x=0.0;
      TranslationInPlaceAttempts[j][i].y=0.0;
      TranslationInPlaceAttempts[j][i].z=0.0;
      TranslationInPlaceAccepted[j][i].x=0.0;
      TranslationInPlaceAccepted[j][i].y=0.0;
      TranslationInPlaceAccepted[j][i].z=0.0;
      TotalTranslationInPlaceAttempts[j][i].x=0.0;
      TotalTranslationInPlaceAttempts[j][i].y=0.0;
      TotalTranslationInPlaceAttempts[j][i].z=0.0;
      TotalTranslationInPlaceAccepted[j][i].x=0.0;
      TotalTranslationInPlaceAccepted[j][i].y=0.0;
      TotalTranslationInPlaceAccepted[j][i].z=0.0;
      ReinsertionInPlaceAttempts[j][i]=0.0;
      ReinsertionInPlaceAccepted[j][i][0]=0.0;
      ReinsertionInPlaceAccepted[j][i][1]=0.0;

      for(k=0;k<NumberOfComponents;k++)
      {
        IdentityChangeAttempts[j][i][k]=0.0;
        IdentityChangeAccepted[j][i][k][0]=0.0;
        IdentityChangeAccepted[j][i][k][1]=0.0;

        GibbsIdentityChangeAttempts[j][i][k]=0.0;
        GibbsIdentityChangeAccepted[j][i][k][0]=0.0;
        GibbsIdentityChangeAccepted[j][i][k][1]=0.0;

        GibbsSwapAttempts[k][j]=0.0;
        GibbsSwapAccepted[k][j]=0.0;
      }
    }

    for(i=0;i<NumberOfSystems;i++)
    {
      ParallelTemperingAttempts[i][j]=0.0;
      ParallelTemperingAccepted[i][j]=0.0;

      HyperParallelTemperingAttempts[i][j]=0.0;
      HyperParallelTemperingAccepted[i][j]=0.0;

      ParallelMolFractionAttempts[i][j]=0.0;
      ParallelMolFractionAccepted[i][j]=0.0;
    }

    ChiralInversionAttempts[j]=0.0;
    ChiralInversionAccepted[j]=0.0;

    VolumeChangeAttempts[j]=0.0;
    VolumeChangeAccepted[j]=0.0;
    TotalVolumeChangeAttempts[j]=0.0;
    TotalVolumeChangeAccepted[j]=0.0;

    BoxShapeChangeAttempts[j].ax=0.0; BoxShapeChangeAttempts[j].bx=0.0; BoxShapeChangeAttempts[j].cx=0.0;
    BoxShapeChangeAttempts[j].ay=0.0; BoxShapeChangeAttempts[j].by=0.0; BoxShapeChangeAttempts[j].cy=0.0;
    BoxShapeChangeAttempts[j].az=0.0; BoxShapeChangeAttempts[j].bz=0.0; BoxShapeChangeAttempts[j].cz=0.0;

    BoxShapeChangeAccepted[j].ax=0.0; BoxShapeChangeAccepted[j].bx=0.0; BoxShapeChangeAccepted[j].cx=0.0;
    BoxShapeChangeAccepted[j].ay=0.0; BoxShapeChangeAccepted[j].by=0.0; BoxShapeChangeAccepted[j].cy=0.0;
    BoxShapeChangeAccepted[j].az=0.0; BoxShapeChangeAccepted[j].bz=0.0; BoxShapeChangeAccepted[j].cz=0.0;

    GibbsVolumeChangeAttempts[j]=0.0;
    GibbsVolumeChangeAccepted[j]=0.0;

    HybridNVEDrift[j]=0.0;
    HybridNVEDriftCount[j]=0.0;

    HybridNVEAttempts[j]=0.0;
    HybridNVEAccepted[j]=0.0;

    HybridNVEStartTemperature[j]=0.0;
    HybridNVEStartTranslationalTemperature[j]=0.0;
    HybridNVEStartRotationalTemperature[j]=0.0;
    HybridNVEStartTemperatureFramework[j]=0.0;
    HybridNVEStartTemperatureAdsorbate[j]=0.0;
    HybridNVEStartTemperatureCation[j]=0.0;
    HybridNVEStartTemperatureCount[j]=0.0;
    HybridNVEStartTemperatureTranslationCount[j]=0.0;
    HybridNVEStartTemperatureRotationCount[j]=0.0;
    HybridNVEStartTemperatureFrameworkCount[j]=0.0;
    HybridNVEStartTemperatureAdsorbateCount[j]=0.0;
    HybridNVEStartTemperatureCationCount[j]=0.0;

    HybridNVEEndTemperature[j]=0.0;
    HybridNVEEndTranslationalTemperature[j]=0.0;
    HybridNVEEndRotationalTemperature[j]=0.0;
    HybridNVEEndTemperatureFramework[j]=0.0;
    HybridNVEEndTemperatureAdsorbate[j]=0.0;
    HybridNVEEndTemperatureCation[j]=0.0;
    HybridNVEEndTemperatureCount[j]=0.0;
    HybridNVEEndTemperatureTranslationCount[j]=0.0;
    HybridNVEEndTemperatureRotationCount[j]=0.0;
    HybridNVEEndTemperatureFrameworkCount[j]=0.0;
    HybridNVEEndTemperatureAdsorbateCount[j]=0.0;
    HybridNVEEndTemperatureCationCount[j]=0.0;

    HybridNPHDrift[j]=0.0;
    HybridNPHDriftCount[j]=0.0;

    HybridNPHAttempts[j]=0.0;
    HybridNPHAccepted[j]=0.0;

    HybridNPHStartTemperature[j]=0.0;
    HybridNPHStartCellTemperature[j]=0.0;
    HybridNPHStartTranslationalTemperature[j]=0.0;
    HybridNPHStartRotationalTemperature[j]=0.0;
    HybridNPHStartTemperatureFramework[j]=0.0;
    HybridNPHStartTemperatureAdsorbate[j]=0.0;
    HybridNPHStartTemperatureCation[j]=0.0;
    HybridNPHStartTemperatureCount[j]=0.0;
    HybridNPHStartTemperatureTranslationCount[j]=0.0;
    HybridNPHStartTemperatureRotationCount[j]=0.0;
    HybridNPHStartTemperatureFrameworkCount[j]=0.0;
    HybridNPHStartTemperatureAdsorbateCount[j]=0.0;
    HybridNPHStartTemperatureCationCount[j]=0.0;

    HybridNPHEndTemperature[j]=0.0;
    HybridNPHEndCellTemperature[j]=0.0;
    HybridNPHEndTranslationalTemperature[j]=0.0;
    HybridNPHEndRotationalTemperature[j]=0.0;
    HybridNPHEndTemperatureFramework[j]=0.0;
    HybridNPHEndTemperatureAdsorbate[j]=0.0;
    HybridNPHEndTemperatureCation[j]=0.0;
    HybridNPHEndTemperatureCount[j]=0.0;
    HybridNPHEndTemperatureTranslationCount[j]=0.0;
    HybridNPHEndTemperatureRotationCount[j]=0.0;
    HybridNPHEndTemperatureFrameworkCount[j]=0.0;
    HybridNPHEndTemperatureAdsorbateCount[j]=0.0;
    HybridNPHEndTemperatureCationCount[j]=0.0;

    HybridNPHPRDrift[j]=0.0;
    HybridNPHPRDriftCount[j]=0.0;

    HybridNPHPRAttempts[j]=0.0;
    HybridNPHPRAccepted[j]=0.0;

    HybridNPHPRStartTemperature[j]=0.0;
    HybridNPHPRStartCellTemperature[j]=0.0;
    HybridNPHPRStartTranslationalTemperature[j]=0.0;
    HybridNPHPRStartRotationalTemperature[j]=0.0;
    HybridNPHPRStartTemperatureFramework[j]=0.0;
    HybridNPHPRStartTemperatureAdsorbate[j]=0.0;
    HybridNPHPRStartTemperatureCation[j]=0.0;
    HybridNPHPRStartTemperatureCount[j]=0.0;
    HybridNPHPRStartTemperatureTranslationCount[j]=0.0;
    HybridNPHPRStartTemperatureRotationCount[j]=0.0;
    HybridNPHPRStartTemperatureFrameworkCount[j]=0.0;
    HybridNPHPRStartTemperatureAdsorbateCount[j]=0.0;
    HybridNPHPRStartTemperatureCationCount[j]=0.0;

    HybridNPHPREndTemperature[j]=0.0;
    HybridNPHPREndCellTemperature[j]=0.0;
    HybridNPHPREndTranslationalTemperature[j]=0.0;
    HybridNPHPREndRotationalTemperature[j]=0.0;
    HybridNPHPREndTemperatureFramework[j]=0.0;
    HybridNPHPREndTemperatureAdsorbate[j]=0.0;
    HybridNPHPREndTemperatureCation[j]=0.0;
    HybridNPHPREndTemperatureCount[j]=0.0;
    HybridNPHPREndTemperatureTranslationCount[j]=0.0;
    HybridNPHPREndTemperatureRotationCount[j]=0.0;
    HybridNPHPREndTemperatureFrameworkCount[j]=0.0;
    HybridNPHPREndTemperatureAdsorbateCount[j]=0.0;
    HybridNPHPREndTemperatureCationCount[j]=0.0;
  }
}

int ComputeNewPolarizationEnergy(int New,int excl_ads,int excl_cation)
{
  int i,j,k,f1;
  int Type;
  REAL_MATRIX3x3 pol_factor;
  VECTOR electric_field,induced_dipole;
  VECTOR electric_field1,electric_field2;

  UHostBackPolarizationNew[CurrentSystem]=0.0;
  UAdsorbateBackPolarizationNew[CurrentSystem]=0.0;
  UCationBackPolarizationNew[CurrentSystem]=0.0;
  UBackPolarizationNew[CurrentSystem]=0.0;

  OVERLAP=FALSE;

  if((BackPolarization)||(Framework[CurrentSystem].FrameworkModel==FLEXIBLE))
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField.x=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField.y=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField.z=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.x=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.y=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.z=0.0;
      }
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField.z=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField.z=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }

  if(New)
  {
   for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      ReferenceElectricFieldAtTrialPosition[CurrentSystem][i].x=0.0;
      ReferenceElectricFieldAtTrialPosition[CurrentSystem][i].y=0.0;
      ReferenceElectricFieldAtTrialPosition[CurrentSystem][i].z=0.0;
      ElectricFieldAtTrialPosition[CurrentSystem][i].x=0.0;
      ElectricFieldAtTrialPosition[CurrentSystem][i].y=0.0;
      ElectricFieldAtTrialPosition[CurrentSystem][i].z=0.0;
    }
  }


  // calculate New 'trial' electric field
  CalculateTotalInterChargeChargeCoulombElectricFieldMC(New,excl_ads,excl_cation);
  CalculateFrameworkChargeChargeElectricFieldMC(New,excl_ads,excl_cation);
  ComputeStaticElectricFieldEwald(New,excl_ads,excl_cation);

  if((BackPolarization)||(Framework[CurrentSystem].FrameworkModel==FLEXIBLE))
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Type=Framework[CurrentSystem].Atoms[f1][i].Type;
        pol_factor=PseudoAtoms[Type].Polarization;
        electric_field=Framework[CurrentSystem].Atoms[f1][i].ElectricField;
        Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField=electric_field;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.x=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.y=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.z=0.0;
        Framework[CurrentSystem].Atoms[f1][i].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
        Framework[CurrentSystem].Atoms[f1][i].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
        Framework[CurrentSystem].Atoms[f1][i].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;
      }
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      pol_factor=PseudoAtoms[Type].Polarization;
      electric_field=Adsorbates[CurrentSystem][i].Atoms[j].ElectricField;
      Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
      Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
      Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;

      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField=electric_field;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][i].Atoms[j].Type;
      pol_factor=PseudoAtoms[Type].Polarization;
      electric_field=Cations[CurrentSystem][i].Atoms[j].ElectricField;
      Cations[CurrentSystem][i].Atoms[j].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
      Cations[CurrentSystem][i].Atoms[j].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
      Cations[CurrentSystem][i].Atoms[j].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;

      Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField=electric_field;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
    }
  }

  if(New)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      Type=Components[CurrentComponent].Type[i];
      pol_factor=PseudoAtoms[Type].Polarization;
      electric_field=ElectricFieldAtTrialPosition[CurrentSystem][i];
      InducedDipoleAtTrialPosition[CurrentSystem][i].x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
      InducedDipoleAtTrialPosition[CurrentSystem][i].y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
      InducedDipoleAtTrialPosition[CurrentSystem][i].z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;

      ReferenceElectricFieldAtTrialPosition[CurrentSystem][i]=electric_field;
      ElectricFieldAtTrialPosition[CurrentSystem][i].x=0.0;
      ElectricFieldAtTrialPosition[CurrentSystem][i].y=0.0;
      ElectricFieldAtTrialPosition[CurrentSystem][i].z=0.0;
    } 
  }

  // sum the energy
  UHostPolarizationNew[CurrentSystem]=0.0;
  UAdsorbatePolarizationNew[CurrentSystem]=0.0;
  UCationPolarizationNew[CurrentSystem]=0.0;
  UPolarizationNew[CurrentSystem]=0.0;

  if((BackPolarization)||(Framework[CurrentSystem].FrameworkModel==FLEXIBLE))
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          Type=Framework[CurrentSystem].Atoms[f1][i].Type;
          electric_field=Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField;
          induced_dipole=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;
          UHostPolarizationNew[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
        }
      }
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      induced_dipole=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;
      electric_field=Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField;
      UAdsorbatePolarizationNew[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][i].Atoms[j].Type;
      induced_dipole=Cations[CurrentSystem][i].Atoms[j].InducedDipole;
      electric_field=Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField;
      UCationPolarizationNew[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }

  if(New)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      Type=Components[CurrentComponent].Type[i];
      induced_dipole=InducedDipoleAtTrialPosition[CurrentSystem][i];
      electric_field=ReferenceElectricFieldAtTrialPosition[CurrentSystem][i];
      UPolarizationNew[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }

  if(!BackPolarization) return 0;

  // iterate to convergence
  //
  //

  for(k=0;k<NumberOfBackPolarizationSteps;k++)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.x=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.y=0.0;
        Framework[CurrentSystem].Atoms[f1][i].ElectricField.z=0.0;
      }
    }

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
        Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
        Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
      }
    }

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Cations[CurrentSystem][i].Atoms[j].ElectricField.x=0.0;
        Cations[CurrentSystem][i].Atoms[j].ElectricField.y=0.0;
        Cations[CurrentSystem][i].Atoms[j].ElectricField.z=0.0;
      }
    }

    if(New)
    {
      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
      {
        ElectricFieldAtTrialPosition[CurrentSystem][i].x=0.0;
        ElectricFieldAtTrialPosition[CurrentSystem][i].y=0.0;
        ElectricFieldAtTrialPosition[CurrentSystem][i].z=0.0;
      }
    }

    // compute electric field generate by the induced dipoles
    CalculateInterElectricFieldFromInducedDipoleMC(New,excl_ads,excl_cation);
    CalculateFrameworkMoleculeElectricFieldFromInducedDipoleMC(New,excl_ads,excl_cation);
    ComputeElectricFieldFromInducedDipolesEwaldMC(New,excl_ads,excl_cation);


    // computer New induced dipoles from electric field
    if((BackPolarization)||(Framework[CurrentSystem].FrameworkModel==FLEXIBLE))
    {
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          Type=Framework[CurrentSystem].Atoms[f1][i].Type;
          pol_factor=PseudoAtoms[Type].Polarization;
          electric_field1=Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField;
          electric_field2=Framework[CurrentSystem].Atoms[f1][i].ElectricField;
          electric_field.x=electric_field1.x+electric_field2.x;
          electric_field.y=electric_field1.y+electric_field2.y;
          electric_field.z=electric_field1.z+electric_field2.z;
          Framework[CurrentSystem].Atoms[f1][i].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
          Framework[CurrentSystem].Atoms[f1][i].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
          Framework[CurrentSystem].Atoms[f1][i].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;
        }
      }
    }

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        pol_factor=PseudoAtoms[Type].Polarization;
        electric_field1=Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField;
        electric_field2=Adsorbates[CurrentSystem][i].Atoms[j].ElectricField;
        electric_field.x=electric_field1.x+electric_field2.x;
        electric_field.y=electric_field1.y+electric_field2.y;
        electric_field.z=electric_field1.z+electric_field2.z;
        Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
        Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
        Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;
      }
    }

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Type=Cations[CurrentSystem][i].Atoms[j].Type;
        pol_factor=PseudoAtoms[Type].Polarization;
        electric_field1=Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField;
        electric_field2=Cations[CurrentSystem][i].Atoms[j].ElectricField;
        electric_field.x=electric_field1.x+electric_field2.x;
        electric_field.y=electric_field1.y+electric_field2.y;
        electric_field.z=electric_field1.z+electric_field2.z;
        Cations[CurrentSystem][i].Atoms[j].InducedDipole.x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
        Cations[CurrentSystem][i].Atoms[j].InducedDipole.y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
        Cations[CurrentSystem][i].Atoms[j].InducedDipole.z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;
      }
    }
    if(New)
    {
      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
      {
        Type=Components[CurrentComponent].Type[i];
        pol_factor=PseudoAtoms[Type].Polarization;
        electric_field1=ReferenceElectricFieldAtTrialPosition[CurrentSystem][i];
        electric_field2=ElectricFieldAtTrialPosition[CurrentSystem][i];
        electric_field.x=electric_field1.x+electric_field2.x;
        electric_field.y=electric_field1.y+electric_field2.y;
        electric_field.z=electric_field1.z+electric_field2.z;
        InducedDipoleAtTrialPosition[CurrentSystem][i].x=pol_factor.ax*electric_field.x+pol_factor.ay*electric_field.y+pol_factor.az*electric_field.z;
        InducedDipoleAtTrialPosition[CurrentSystem][i].y=pol_factor.bx*electric_field.x+pol_factor.by*electric_field.y+pol_factor.bz*electric_field.z;
        InducedDipoleAtTrialPosition[CurrentSystem][i].z=pol_factor.cx*electric_field.x+pol_factor.cy*electric_field.y+pol_factor.cz*electric_field.z;
      }
    }
  }

  // sum the energy
  UHostBackPolarizationNew[CurrentSystem]=0.0;
  UAdsorbateBackPolarizationNew[CurrentSystem]=0.0;
  UCationBackPolarizationNew[CurrentSystem]=0.0;
  UBackPolarizationNew[CurrentSystem]=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Type=Framework[CurrentSystem].Atoms[f1][i].Type;
      pol_factor=PseudoAtoms[Type].Polarization;
      electric_field=Framework[CurrentSystem].Atoms[f1][i].ReferenceElectricField;
      induced_dipole=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;
      UHostBackPolarizationNew[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      induced_dipole=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;
      electric_field=Adsorbates[CurrentSystem][i].Atoms[j].ReferenceElectricField;
      UAdsorbateBackPolarizationNew[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][i].Atoms[j].Type;
      induced_dipole=Cations[CurrentSystem][i].Atoms[j].InducedDipole;
      electric_field=Cations[CurrentSystem][i].Atoms[j].ReferenceElectricField;
      UCationBackPolarizationNew[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }  
  if(New)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      Type=Components[CurrentComponent].Type[i];
      induced_dipole=InducedDipoleAtTrialPosition[CurrentSystem][i];
      electric_field=ReferenceElectricFieldAtTrialPosition[CurrentSystem][i];
      UBackPolarizationNew[CurrentSystem]-=0.5*(induced_dipole.x*electric_field.x+induced_dipole.y*electric_field.y+induced_dipole.z*electric_field.z);
    }
  }

  // subtract the "static" polarization to get the iterative part only
  UHostBackPolarizationNew[CurrentSystem]-=UHostPolarizationNew[CurrentSystem];
  UAdsorbateBackPolarizationNew[CurrentSystem]-=UAdsorbatePolarizationNew[CurrentSystem];
  UCationBackPolarizationNew[CurrentSystem]-=UCationPolarizationNew[CurrentSystem];
  UBackPolarizationNew[CurrentSystem]-=UPolarizationNew[CurrentSystem];

  return 0;
}

/*********************************************************************************************************
 * Name       | GetDisplacementVector                                                                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Get the displacement vector for a translation move.                                      *
 * Parameters | -                                                                                        *
 * Note       |                                                                                          *
 * Used       | in the translation moves to sample e.g. on a plane for dcTST                             *
 *********************************************************************************************************/

VECTOR GetDisplacementVector(REAL vNew)
{
  VECTOR displacement,s;
  REAL choice;

  // initialize displacement to zero
  displacement.x=0.0;
  displacement.y=0.0;
  displacement.z=0.0;

  switch(Components[CurrentComponent].TranslationDirection)
  {
    case XYZ_DIR:
      choice=Dimension*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
      else if(choice<2.0)
        displacement.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
      else
        displacement.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].z;
      break;
    case XY_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
      else if(choice<2.0)
        displacement.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
      break;
    case XZ_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
      else if(choice<2.0)
        displacement.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].z;
      break;
    case YZ_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
      else if(choice<2.0)
        displacement.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].z;
      break;
    case X_DIR:
      displacement.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
      break;
    case Y_DIR:
      displacement.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
      break;
    case Z_DIR:
      displacement.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].z;
      break;
    case ABC_DIR:
      choice=Dimension*RandomNumber();
      if(choice<1.0)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=0.0;
        s.z=0.0;
      }
      else if(choice<2.0)
      {
        s.x=0.0;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=0.0;
      }
      else
      {
        s.x=0.0;
        s.y=0.0;
        s.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].z;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case AB_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=0.0;
        s.z=0.0;
      }
      else
      {
        s.x=0.0;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=0.0;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case AC_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=0.0;
        s.z=0.0;
      }
      else
      {
        s.x=0.0;
        s.y=0.0;
        s.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].z;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case BC_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=0.0;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=0.0;
      }
      else
      {
        s.x=0.0;
        s.y=0.0;
        s.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].z;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case ORTHOGONAL_TO_AB_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.z=0.0;
      }
      else
      {
        s.x=0.0;
        s.y=0.0;
        s.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].z;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case ORTHOGONAL_TO_BC_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=0.0;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
      }
      else
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=0.0;
        s.z=0.0;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case ORTHOGONAL_TO_AC_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=0.0;
        s.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
      }
      else
      {
        s.x=0.0;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=0.0;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case ORTHOGONAL_TO_O_AB_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=-vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.z=0.0;
      }
      else
      {
        s.x=0.0;
        s.y=0.0;
        s.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].z;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case ORTHOGONAL_TO_O_AC_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=0.0;
        s.z=-vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
      }
      else
      {
        s.x=0.0;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=0.0;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case ORTHOGONAL_TO_O_BC_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=0.0;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=-vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
      }
      else
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=0.0;
        s.z=0.0;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case ORTHOGONAL_TO_A_BC_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.z=0.0;
      }
      else
      {
        s.x=-vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=-2.0*vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case ORTHOGONAL_TO_B_AC_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.z=0.0;
      }
      else
      {
        s.x=-vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=2.0*vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case ORTHOGONAL_TO_C_AB_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=0.0;
        s.z=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
      }
      else
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.y=-2.0*vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=-vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
    case ORTHOGONAL_TO_O_ABC_DIR:
      choice=RandomNumber();
      if(choice<0.5)
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.y=-vNew*MaximumTranslation[CurrentSystem][CurrentComponent].x;
        s.z=0.0;
      }
      else
      {
        s.x=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.y=vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
        s.z=-2.0*vNew*MaximumTranslation[CurrentSystem][CurrentComponent].y;
      }
      displacement.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      displacement.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      displacement.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
  }

  return displacement;
}

/*********************************************************************************************************
 * Name       | TranslationMoveAdsorbate                                                                 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Translation Monte Carlo move of an adsorbate molecule                                    *
 * Parameters | -                                                                                        *
 * Note       | The translation moves is the most often and basic used move to change the system.        *
 *            | Before calling this function, a component has been randomly chosen. A molecule of this   *
 *            | type is chosen randomly and an attempt is made to displace the molecule. The internal    *
 *            | configuration of the molecule remains unchanged.                                         *
 *            | The displacement can be restricted to specified directions, e.g. a plane.                *
 *********************************************************************************************************/

int TranslationMoveAdsorbate(void)
{
  int i,nr_atoms;
  REAL vNew;
  VECTOR displacement;
  POINT pos;
  int StartingBead;
  REAL BiasingWeight=1.0;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfAdsorbateMolecules[CurrentSystem]==0) return 0;
  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentCationMolecule=-1;
  CurrentAdsorbateMolecule=SelectRandomMoleculeOfType(CurrentComponent);

  // calculate a random displacement
  vNew=2.0*RandomNumber()-1.0;

  // choose a possible displacement
  displacement=GetDisplacementVector(vNew);

  // update statistics and register a translation move
  if(fabs(displacement.x)>0.0)
    TranslationAttempts[CurrentSystem][CurrentComponent].x+=1.0; 
  if(fabs(displacement.y)>0.0)
    TranslationAttempts[CurrentSystem][CurrentComponent].y+=1.0; 
  if(fabs(displacement.z)>0.0)
    TranslationAttempts[CurrentSystem][CurrentComponent].z+=1.0; 

  // calculate the energy of the current configuration translated
  nr_atoms=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].NumberOfAtoms;
  for(i=0;i<nr_atoms;i++)
  {
    TrialPosition[CurrentSystem][i].x=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position.x+displacement.x;
    TrialPosition[CurrentSystem][i].y=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position.y+displacement.y;
    TrialPosition[CurrentSystem][i].z=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position.z+displacement.z;
  }

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<nr_atoms;i++)
  {
    CFVDWScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFChargeScalingParameter;
  }

  // calculate anisotropic sites
  CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);

  for(i=0;i<nr_atoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  // check if still a valid position
  if(Components[CurrentComponent].RestrictMoves)
  {
    StartingBead=Components[CurrentComponent].StartingBead;
    pos=TrialPosition[CurrentSystem][StartingBead];

    if(!ValidCartesianPoint(CurrentComponent,pos)) return 0;
  }

  // compute inter-molecular energy differences
  CalculateInterVDWEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeChargeEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute energy differences framework-adsorbate
  CalculateFrameworkAdsorbateVDWEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkAdsorbateChargeChargeEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute the energy differenes in Fourier-space
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,0);

  if(ComputePolarization)
  { 
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
  }

  StartingBead=Components[CurrentComponent].StartingBead;

  switch(Components[CurrentComponent].Biased)
  {
    case UMBRELLA:
      BiasingWeight=BiasingPotentialUmbrella(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialUmbrella(CurrentComponent,Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position);
      break;
    case RUIZ_MONTERO:
      BiasingWeight=BiasingPotentialRuizMontero(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialRuizMontero(CurrentComponent,Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position);
      break;
    case NO_BIASING:
      BiasingWeight=1.0;
      break;
  }

  DeltaU=UHostVDWDelta[CurrentSystem]+UAdsorbateVDWDelta[CurrentSystem]+UCationVDWDelta[CurrentSystem]+
         UHostChargeChargeRealDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
         UHostChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
         UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+
         UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
         UDeltaPolarization;


  if(RandomNumber()<BiasingWeight*exp(-Beta[CurrentSystem]*DeltaU))
  {
    #ifdef DEBUG
      printf("Translation adsorbate accepted\n");
    #endif
    if(fabs(displacement.x)>0.0)
      TranslationAccepted[CurrentSystem][CurrentComponent].x+=1.0;
    if(fabs(displacement.y)>0.0)
      TranslationAccepted[CurrentSystem][CurrentComponent].y+=1.0;
    if(fabs(displacement.z)>0.0)
      TranslationAccepted[CurrentSystem][CurrentComponent].z+=1.0;

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UCationVDWDelta[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWDelta[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleRealDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                 UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                     UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem];

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }

    nr_atoms=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassAdsorbate(CurrentAdsorbateMolecule);

    UTotal[CurrentSystem]+=DeltaU;
  }
  else
  {
    #ifdef DEBUG
      printf("Translation adsorbate rejected\n");
    #endif
  }

  return 0;
}

/*********************************************************************************************************
 * Name       | TranslationMoveCation                                                                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Translation Monte Carlo move of a cation molecule                                        *
 * Parameters | -                                                                                        *
 * Note       | The translation moves is the most often and basic used move to change the system.        *
 *            | Before calling this function, a component has been randomly chosen. A molecule of this   *
 *            | type is chosen randomly and an attempt is made to displace the molecule. The internal    *
 *            | configuration of the molecule remains unchanged.                                         *
 *            | The displacement can be restricted to specified directions, e.g. a plane.                *
 *********************************************************************************************************/

int TranslationMoveCation(void)
{
  int i,nr_atoms;
  REAL vNew;
  VECTOR displacement;
  POINT pos;
  int StartingBead;
  REAL BiasingWeight=1.0;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfCationMolecules[CurrentSystem]==0) return 0;
  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentAdsorbateMolecule=-1;
  CurrentCationMolecule=SelectRandomMoleculeOfType(CurrentComponent);

  // calculate a random displacement
  vNew=2.0*RandomNumber()-1.0;

  // choose a possible displacement
  displacement=GetDisplacementVector(vNew);

  // update statistics and register a translation move
  if(fabs(displacement.x)>0.0)
    TranslationAttempts[CurrentSystem][CurrentComponent].x+=1.0; 
  if(fabs(displacement.y)>0.0)
    TranslationAttempts[CurrentSystem][CurrentComponent].y+=1.0; 
  if(fabs(displacement.z)>0.0)
    TranslationAttempts[CurrentSystem][CurrentComponent].z+=1.0; 

  // calculate the energy of the current configuration translated
  nr_atoms=Cations[CurrentSystem][CurrentCationMolecule].NumberOfAtoms;
  for(i=0;i<nr_atoms;i++)
  {
    TrialPosition[CurrentSystem][i].x=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position.x+displacement.x;
    TrialPosition[CurrentSystem][i].y=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position.y+displacement.y;
    TrialPosition[CurrentSystem][i].z=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position.z+displacement.z;
  }

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<nr_atoms;i++)
  {
    CFVDWScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFChargeScalingParameter;
  }

  // calculate anisotropic sites
  CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);

  for(i=0;i<nr_atoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  // check if still a valid position
  if(Components[CurrentComponent].RestrictMoves)
  {
    StartingBead=Components[CurrentComponent].StartingBead;
    pos=TrialPosition[CurrentSystem][StartingBead];

    if(!ValidCartesianPoint(CurrentComponent,pos)) return 0;
  }

  // compute inter-molecular energy differences
  CalculateInterVDWEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeChargeEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeBondDipoleEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterBondDipoleBondDipoleEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute energy differences framework-cation
  CalculateFrameworkCationVDWEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkCationChargeChargeEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkCationChargeBondDipoleEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute the energy differenes in Fourier-space
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
     CalculateEwaldFourierCation(TRUE,TRUE,CurrentCationMolecule,0);

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,-1,CurrentCationMolecule);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       UAdsorbatePolarizationNew[CurrentSystem]-UAdsorbatePolarization[CurrentSystem]+
                       (UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       UAdsorbateBackPolarizationNew[CurrentSystem]-UAdsorbateBackPolarization[CurrentSystem]+
                       (UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UCationBackPolarization[CurrentSystem];
  }

  StartingBead=Components[CurrentComponent].StartingBead;

  switch(Components[CurrentComponent].Biased)
  {
    case UMBRELLA:
      BiasingWeight=BiasingPotentialUmbrella(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialUmbrella(CurrentComponent,Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBead].Position);
      break;
    case RUIZ_MONTERO:
      BiasingWeight=BiasingPotentialRuizMontero(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialRuizMontero(CurrentComponent,Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBead].Position);
      break;
    case NO_BIASING:
      BiasingWeight=1.0;
      break;
  }

  DeltaU=UHostVDWDelta[CurrentSystem]+UAdsorbateVDWDelta[CurrentSystem]+UCationVDWDelta[CurrentSystem]+
         UHostChargeChargeRealDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
         UHostChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
         UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+
         UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
         UDeltaPolarization;

  if(RandomNumber()<BiasingWeight*exp(-Beta[CurrentSystem]*DeltaU))
  {
    if(fabs(displacement.x)>0.0)
      TranslationAccepted[CurrentSystem][CurrentComponent].x+=1.0;
    if(fabs(displacement.y)>0.0)
      TranslationAccepted[CurrentSystem][CurrentComponent].y+=1.0;
    if(fabs(displacement.z)>0.0)
      TranslationAccepted[CurrentSystem][CurrentComponent].z+=1.0;

    UCationCation[CurrentSystem]+=UCationVDWDelta[CurrentSystem];
    UCationCationVDW[CurrentSystem]+=UCationVDWDelta[CurrentSystem];
    UHostCation[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UHostCationVDW[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UCationCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleRealDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                           UCationChargeBondDipoleRealDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                    UCationChargeBondDipoleRealDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(0);
    }

    nr_atoms=Cations[CurrentSystem][CurrentCationMolecule].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassCation(CurrentCationMolecule);

    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

void TranslationMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    TranslationMoveCation();
  else
    TranslationMoveAdsorbate();
}

/*********************************************************************************************************
 * Name       | OptimizeTranslationAcceptence                                                            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Adjusts the maximum displacement on-the-fly to obtain an acceptance-rate of 50%.         *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void OptimizeTranslationAcceptence(void)
{
  int i;
  VECTOR ratio,vandr;

  for(i=0;i<NumberOfComponents;i++)
  {
    if(TranslationAttempts[CurrentSystem][i].x>0.0)
      ratio.x=TranslationAccepted[CurrentSystem][i].x/TranslationAttempts[CurrentSystem][i].x;
    else
      ratio.x=0.0;
    if(TranslationAttempts[CurrentSystem][i].y>0.0)
      ratio.y=TranslationAccepted[CurrentSystem][i].y/TranslationAttempts[CurrentSystem][i].y;
    else
      ratio.y=0.0;
    if(TranslationAttempts[CurrentSystem][i].z>0.0)
      ratio.z=TranslationAccepted[CurrentSystem][i].z/TranslationAttempts[CurrentSystem][i].z;
    else
      ratio.z=0.0;

    vandr.x=ratio.x/TargetAccRatioTranslation;
    if(vandr.x>1.5) vandr.x=1.5;
    else if(vandr.x<0.5) vandr.x=0.5;
    MaximumTranslation[CurrentSystem][i].x*=vandr.x;
    if(MaximumTranslation[CurrentSystem][i].x<0.01)
       MaximumTranslation[CurrentSystem][i].x=0.01;
    if(MaximumTranslation[CurrentSystem][i].x>1.0)
       MaximumTranslation[CurrentSystem][i].x=1.0;

    vandr.y=ratio.y/TargetAccRatioTranslation;
    if(vandr.y>1.5) vandr.y=1.5;
    else if(vandr.y<0.5) vandr.y=0.5;
    MaximumTranslation[CurrentSystem][i].y*=vandr.y;
    if(MaximumTranslation[CurrentSystem][i].y<0.01)
       MaximumTranslation[CurrentSystem][i].y=0.01;
    if(MaximumTranslation[CurrentSystem][i].y>1.0)
       MaximumTranslation[CurrentSystem][i].y=1.0;

    vandr.z=ratio.z/TargetAccRatioTranslation;
     if(vandr.z>1.5) vandr.z=1.5;
    else if(vandr.z<0.5) vandr.z=0.5;
    MaximumTranslation[CurrentSystem][i].z*=vandr.z;
    if(MaximumTranslation[CurrentSystem][i].z<0.01)
       MaximumTranslation[CurrentSystem][i].z=0.01;
    if(MaximumTranslation[CurrentSystem][i].z>1.0)
       MaximumTranslation[CurrentSystem][i].z=1.0;

    TotalTranslationAttempts[CurrentSystem][i].x+=TranslationAttempts[CurrentSystem][i].x;
    TotalTranslationAccepted[CurrentSystem][i].x+=TranslationAccepted[CurrentSystem][i].x;
    TranslationAttempts[CurrentSystem][i].x=TranslationAccepted[CurrentSystem][i].x=0.0;

    TotalTranslationAttempts[CurrentSystem][i].y+=TranslationAttempts[CurrentSystem][i].y;
    TotalTranslationAccepted[CurrentSystem][i].y+=TranslationAccepted[CurrentSystem][i].y;
    TranslationAttempts[CurrentSystem][i].y=TranslationAccepted[CurrentSystem][i].y=0.0;

    TotalTranslationAttempts[CurrentSystem][i].z+=TranslationAttempts[CurrentSystem][i].z;
    TotalTranslationAccepted[CurrentSystem][i].z+=TranslationAccepted[CurrentSystem][i].z;
    TranslationAttempts[CurrentSystem][i].z=TranslationAccepted[CurrentSystem][i].z=0.0;
  }
}

void PrintTranslationStatistics(FILE *FilePtr)
{
  int i,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfTranslationMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the translation move:\n");
    fprintf(FilePtr,"======================================\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      fprintf(FilePtr,"Component %d [%s]\n",i,Components[i].Name);
      fprintf(FilePtr,"\ttotal        %f %f %f\n",
        (double)TotalTranslationAttempts[CurrentSystem][i].x,
        (double)TotalTranslationAttempts[CurrentSystem][i].y,
        (double)TotalTranslationAttempts[CurrentSystem][i].z);
      fprintf(FilePtr,"\tsuccesfull   %f %f %f\n",
        (double)TotalTranslationAccepted[CurrentSystem][i].x,
        (double)TotalTranslationAccepted[CurrentSystem][i].y,
        (double)TotalTranslationAccepted[CurrentSystem][i].z);

      fprintf(FilePtr,"\taccepted   %f %f %f\n",
        (double)(TotalTranslationAttempts[CurrentSystem][i].x>(REAL)0.0?
           TotalTranslationAccepted[CurrentSystem][i].x/TotalTranslationAttempts[CurrentSystem][i].x:(REAL)0.0),
        (double)(TotalTranslationAttempts[CurrentSystem][i].y>(REAL)0.0?
            TotalTranslationAccepted[CurrentSystem][i].y/TotalTranslationAttempts[CurrentSystem][i].y:(REAL)0.0),
        (double)(TotalTranslationAttempts[CurrentSystem][i].z>(REAL)0.0?
            TotalTranslationAccepted[CurrentSystem][i].z/TotalTranslationAttempts[CurrentSystem][i].z:(REAL)0.0));
      fprintf(FilePtr,"\tdisplacement %f %f %f\n",
        (double)MaximumTranslation[CurrentSystem][i].x,
        (double)MaximumTranslation[CurrentSystem][i].y,
        (double)MaximumTranslation[CurrentSystem][i].z);
      fprintf(FilePtr,"\n");
    }
    fprintf(FilePtr,"\n");
  }
  else
   fprintf(FilePtr,"Translation move was OFF for all components\n\n");
}


/*********************************************************************************************************
 * Name       | RandomTranslationMoveAdsorbate                                                           *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Random translation Monte Carlo move of an adsorbate molecule                             *
 * Parameters | -                                                                                        *
 * Note       | The random translation move translates the molecule to a random position in the box.     *
 *            | Before calling this function, a component has been randomly chosen. A molecule of this   *
 *            | type is chosen randomly and an attempt is made to displace the molecule. The internal    *
 *            | configuration of the molecule remains unchanged.                                         *
 *********************************************************************************************************/

int RandomTranslationMoveAdsorbate(void)
{
  int i,nr_atoms;
  VECTOR displacement;
  POINT pos;
  int StartingBead;
  REAL BiasingWeight=1.0;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfAdsorbateMolecules[CurrentSystem]==0) return 0;
  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentAdsorbateMolecule=SelectRandomMoleculeOfType(CurrentComponent);
  CurrentCationMolecule=-1;

  RandomTranslationAttempts[CurrentSystem][CurrentComponent]+=1.0;

  // calculate a random displacement
  displacement.x=(RandomNumber()-0.5)*Box[CurrentSystem].ax;
  displacement.y=(RandomNumber()-0.5)*Box[CurrentSystem].by;
  displacement.z=(RandomNumber()-0.5)*Box[CurrentSystem].cz;

  // calculate the energy of the current configuration translated
  nr_atoms=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].NumberOfAtoms;
  for(i=0;i<nr_atoms;i++)
  {
    TrialPosition[CurrentSystem][i].x=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position.x+displacement.x;
    TrialPosition[CurrentSystem][i].y=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position.y+displacement.y;
    TrialPosition[CurrentSystem][i].z=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position.z+displacement.z;
  }

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<nr_atoms;i++)
  {
    CFVDWScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFChargeScalingParameter;
  }

  // calculate anisotropic sites
  CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);

  for(i=0;i<nr_atoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }


  // check if still a valid position
  if(Components[CurrentComponent].RestrictMoves)
  {
    StartingBead=Components[CurrentComponent].StartingBead;
    pos=TrialPosition[CurrentSystem][StartingBead];

    if(!ValidCartesianPoint(CurrentComponent,pos)) return 0;
  }

  // compute inter-molecular energy differences
  CalculateInterVDWEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeChargeEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute energy differences framework-adsorbate
  CalculateFrameworkAdsorbateVDWEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkAdsorbateChargeChargeEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute the energy differenes in Fourier-space
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
     CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,0);

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
  }

  StartingBead=Components[CurrentComponent].StartingBead;

  switch(Components[CurrentComponent].Biased)
  {
    case UMBRELLA:
      BiasingWeight=BiasingPotentialUmbrella(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialUmbrella(CurrentComponent,Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position);
      break;
    case RUIZ_MONTERO:
      BiasingWeight=BiasingPotentialRuizMontero(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialRuizMontero(CurrentComponent,Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position);
      break;
    case NO_BIASING:
      BiasingWeight=1.0;
      break;
  }

  DeltaU=UHostVDWDelta[CurrentSystem]+UAdsorbateVDWDelta[CurrentSystem]+UCationVDWDelta[CurrentSystem]+
         UHostChargeChargeRealDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
         UHostChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
         UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+
         UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
         UDeltaPolarization;

  if(RandomNumber()<BiasingWeight*exp(-Beta[CurrentSystem]*DeltaU))
  {
    RandomTranslationAccepted[CurrentSystem][CurrentComponent]+=1.0;

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UCationVDWDelta[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWDelta[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleRealDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                 UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                     UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem];

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }

    nr_atoms=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassAdsorbate(CurrentAdsorbateMolecule);

    UTotal[CurrentSystem]+=DeltaU;
  }

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassCation(i);

  return 0;
}

/*********************************************************************************************************
 * Name       | RandomTranslationMoveCation                                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Random translation Monte Carlo move of a cation molecule                                 *
 * Parameters | -                                                                                        *
 * Note       | The random translation move translates the molecule to a random position in the box.     *
 *            | Before calling this function, a component has been randomly chosen. A molecule of this   *
 *            | type is chosen randomly and an attempt is made to displace the molecule. The internal    *
 *            | configuration of the molecule remains unchanged.                                         *
 *********************************************************************************************************/

int RandomTranslationMoveCation(void)
{
  int i,nr_atoms;
  VECTOR displacement;
  POINT pos;
  int StartingBead;
  REAL BiasingWeight=1.0;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfCationMolecules[CurrentSystem]==0) return 0;
  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentAdsorbateMolecule=-1;
  CurrentCationMolecule=SelectRandomMoleculeOfType(CurrentComponent);

  RandomTranslationAttempts[CurrentSystem][CurrentComponent]+=1.0;

  // calculate a random displacement
  displacement.x=(RandomNumber()-0.5)*Box[CurrentSystem].ax;
  displacement.y=(RandomNumber()-0.5)*Box[CurrentSystem].by;
  displacement.z=(RandomNumber()-0.5)*Box[CurrentSystem].cz;

  // calculate the energy of the current configuration translated
  nr_atoms=Cations[CurrentSystem][CurrentCationMolecule].NumberOfAtoms;
  for(i=0;i<nr_atoms;i++)
  {
    TrialPosition[CurrentSystem][i].x=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position.x+displacement.x;
    TrialPosition[CurrentSystem][i].y=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position.y+displacement.y;
    TrialPosition[CurrentSystem][i].z=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position.z+displacement.z;
  }

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<nr_atoms;i++)
  {
    CFVDWScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFChargeScalingParameter;
  }

  // calculate anisotropic sites
  CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);

  for(i=0;i<nr_atoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  // check if still a valid position
  if(Components[CurrentComponent].RestrictMoves)
  {
    StartingBead=Components[CurrentComponent].StartingBead;
    pos=TrialPosition[CurrentSystem][StartingBead];

    if(!ValidCartesianPoint(CurrentComponent,pos)) return 0;
  }

  // compute inter-molecular energy differences
  CalculateInterVDWEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeChargeEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeBondDipoleEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterBondDipoleBondDipoleEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute energy differences framework-cation
  CalculateFrameworkCationVDWEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkCationChargeChargeEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkCationChargeBondDipoleEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute the energy differenes in Fourier-space
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
     CalculateEwaldFourierCation(TRUE,TRUE,CurrentCationMolecule,0);

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,-1,CurrentCationMolecule);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       UAdsorbatePolarizationNew[CurrentSystem]-UAdsorbatePolarization[CurrentSystem]+
                       (UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       UAdsorbateBackPolarizationNew[CurrentSystem]-UAdsorbateBackPolarization[CurrentSystem]+
                       (UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UCationBackPolarization[CurrentSystem];
  }

  StartingBead=Components[CurrentComponent].StartingBead;

  BiasingWeight=1.0;
  switch(Components[CurrentComponent].Biased)
  {
    case UMBRELLA:
      BiasingWeight=BiasingPotentialUmbrella(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialUmbrella(CurrentComponent,Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBead].Position);
      break;
    case RUIZ_MONTERO:
      BiasingWeight=BiasingPotentialRuizMontero(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialRuizMontero(CurrentComponent,Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBead].Position);
      break;
    case NO_BIASING:
      BiasingWeight=1.0;
      break;
  }

  DeltaU=UHostVDWDelta[CurrentSystem]+UAdsorbateVDWDelta[CurrentSystem]+UCationVDWDelta[CurrentSystem]+
         UHostChargeChargeRealDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
         UHostChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
         UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+
         UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
         UDeltaPolarization;

  if(RandomNumber()<BiasingWeight*exp(-Beta[CurrentSystem]*DeltaU))
  {
    RandomTranslationAccepted[CurrentSystem][CurrentComponent]+=1.0;

    UAdsorbateCation[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UHostCation[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UHostCationVDW[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UCationCation[CurrentSystem]+=UCationVDWDelta[CurrentSystem];
    UCationCationVDW[CurrentSystem]+=UCationVDWDelta[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UCationCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleRealDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                           UCationChargeBondDipoleRealDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                    UCationChargeBondDipoleRealDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(0);
    }

    nr_atoms=Cations[CurrentSystem][CurrentCationMolecule].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassCation(CurrentCationMolecule);

    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

void RandomTranslationMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    RandomTranslationMoveCation();
  else
    RandomTranslationMoveAdsorbate();
}

void PrintRandomTranslationStatistics(FILE *FilePtr)
{
  int i,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfRandomTranslationMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the random translation move:\n");
    fprintf(FilePtr,"============================================\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      fprintf(FilePtr,"Component %d [%s]\n",i,Components[i].Name);
      fprintf(FilePtr,"\ttotal        %lf\n",
        (double)RandomTranslationAttempts[CurrentSystem][i]);
      fprintf(FilePtr,"\tsuccesfull   %lf\n",
        (double)RandomTranslationAccepted[CurrentSystem][i]);
      fprintf(FilePtr,"\taccepted   %lf\n",
        (double)(RandomTranslationAttempts[CurrentSystem][i]>(REAL)0.0?
          RandomTranslationAccepted[CurrentSystem][i]/RandomTranslationAttempts[CurrentSystem][i]:(REAL)0.0));
      fprintf(FilePtr,"\n");
    }
    fprintf(FilePtr,"\n");
  }
  else
   fprintf(FilePtr,"Random translation move was OFF for all components\n\n");
}


/*********************************************************************************************************
 * Name       | RotationMoveAdsorbate                                                                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Rotation Monte Carlo move of an adsorbate molecule                                       *
 * Parameters | -                                                                                        *
 * Note       | The rotation move randomly rotates the molecule around a randomly chosen vector.         *
 *            | Before calling this function, a component has been randomly chosen. A molecule of this   *
 *            | type is chosen randomly and an attempt is made to rotate the molecule. The internal      *
 *            | configuration of the molecule remains unchanged.                                         *
 *********************************************************************************************************/


int RotationMoveAdsorbate(void)
{
  int i,nr_atoms,start;
  int StartingBead;
  REAL DeltaU;
  VECTOR pos,posA,posB,posC,posD,posE;

  // return if the component currently has zero molecules
  if(NumberOfAdsorbateMolecules[CurrentSystem]==0) return 0;
  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentAdsorbateMolecule=SelectRandomMoleculeOfType(CurrentComponent);
  CurrentCationMolecule=-1;

  RotationAttempts[CurrentSystem][CurrentComponent]+=1.0; 

  start=Components[CurrentComponent].StartingBead;
  nr_atoms=Components[CurrentComponent].NumberOfAtoms;

  for(i=0;i<nr_atoms;i++)
  {
    cord[i].x=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position.x-
              Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[start].Position.x;
    cord[i].y=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position.y-
              Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[start].Position.y;
    cord[i].z=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position.z-
              Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[start].Position.z;
  }

  RandomArrayRotationMatrix(cord,nr_atoms);

  for(i=0;i<nr_atoms;i++)
  {
    TrialPosition[CurrentSystem][i].x=cord[i].x+Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[start].Position.x;
    TrialPosition[CurrentSystem][i].y=cord[i].y+Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[start].Position.y;
    TrialPosition[CurrentSystem][i].z=cord[i].z+Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[start].Position.z;
  }

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<nr_atoms;i++)
  {
    CFVDWScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFChargeScalingParameter;
  }

  // calculate anisotropic sites
  CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);

  for(i=0;i<nr_atoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  if(Components[CurrentComponent].RestrictEnantionface)
  {
    posA=Components[CurrentComponent].EnantiofaceAtoms[0]->Position;
    if((Components[CurrentComponent].EnantiofaceAtomDefinitions[0][0]==ADSORBATE)&&
       (Components[CurrentComponent].EnantiofaceAtomDefinitions[0][1]==CurrentAdsorbateMolecule))
       posA=TrialPosition[CurrentSystem][Components[CurrentComponent].EnantiofaceAtomDefinitions[0][2]];

    posB=Components[CurrentComponent].EnantiofaceAtoms[1]->Position;
    if((Components[CurrentComponent].EnantiofaceAtomDefinitions[1][0]==ADSORBATE)&&
       (Components[CurrentComponent].EnantiofaceAtomDefinitions[1][1]==CurrentAdsorbateMolecule))
       posB=TrialPosition[CurrentSystem][Components[CurrentComponent].EnantiofaceAtomDefinitions[1][2]];

    posC=Components[CurrentComponent].EnantiofaceAtoms[2]->Position;
    if((Components[CurrentComponent].EnantiofaceAtomDefinitions[2][0]==ADSORBATE)&&
       (Components[CurrentComponent].EnantiofaceAtomDefinitions[2][1]==CurrentAdsorbateMolecule))
       posC=TrialPosition[CurrentSystem][Components[CurrentComponent].EnantiofaceAtomDefinitions[2][2]];

    posD=Components[CurrentComponent].EnantiofaceAtoms[3]->Position;
    if((Components[CurrentComponent].EnantiofaceAtomDefinitions[3][0]==ADSORBATE)&&
       (Components[CurrentComponent].EnantiofaceAtomDefinitions[3][1]==CurrentAdsorbateMolecule))
       posD=TrialPosition[CurrentSystem][Components[CurrentComponent].EnantiofaceAtomDefinitions[3][2]];

    posE=Components[CurrentComponent].EnantiofaceAtoms[4]->Position;
    if((Components[CurrentComponent].EnantiofaceAtomDefinitions[4][0]==ADSORBATE)&&
       (Components[CurrentComponent].EnantiofaceAtomDefinitions[4][1]==CurrentAdsorbateMolecule))
       posE=TrialPosition[CurrentSystem][Components[CurrentComponent].EnantiofaceAtomDefinitions[4][2]];

    if(Components[CurrentComponent].Enantioface!=CheckEnantioFace(posA,posB,posC,posD,posE))
      return 0;
  }

  // check if still a valid position
  if(Components[CurrentComponent].RestrictMoves)
  {
    StartingBead=Components[CurrentComponent].StartingBead;
    pos=TrialPosition[CurrentSystem][StartingBead];

    if(!ValidCartesianPoint(CurrentComponent,pos)) return 0;
  }

  // compute inter-molecular energy differences
  CalculateInterVDWEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeChargeEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute energy differences framework-adsorbate
  CalculateFrameworkAdsorbateVDWEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkAdsorbateChargeChargeEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute the energy differenes in Fourier-space
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
     CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,0);

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
  }

  DeltaU=UHostVDWDelta[CurrentSystem]+UAdsorbateVDWDelta[CurrentSystem]+UCationVDWDelta[CurrentSystem]+
         UHostChargeChargeRealDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
         UHostChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
         UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+
         UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
         UDeltaPolarization;

  if(RandomNumber()<exp(-Beta[CurrentSystem]*DeltaU))
  {
    RotationAccepted[CurrentSystem][CurrentComponent]+=1.0;

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UCationVDWDelta[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWDelta[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleRealDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                 UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                     UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem];

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }

    nr_atoms=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassAdsorbate(CurrentAdsorbateMolecule);
    ComputeQuaternionAdsorbate(CurrentAdsorbateMolecule);

    UTotal[CurrentSystem]+=DeltaU;
  }

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassCation(i);


  return 0;
}

/*********************************************************************************************************
 * Name       | RotationMoveCation                                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Rotation Monte Carlo move of a cation molecule                                           *
 * Parameters | -                                                                                        *
 * Note       | The rotation move randomly rotates the molecule around a randomly chosen vector.         *
 *            | Before calling this function, a component has been randomly chosen. A molecule of this   *
 *            | type is chosen randomly and an attempt is made to rotate the molecule. The internal      *
 *            | configuration of the molecule remains unchanged.                                         *
 *********************************************************************************************************/

int RotationMoveCation(void)
{
  int i,nr_atoms,start;
  int StartingBead;
  REAL DeltaU;
  VECTOR pos,posA,posB,posC,posD,posE;

  // return if the component currently has zero molecules
  if(NumberOfCationMolecules[CurrentSystem]==0) return 0;
  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentAdsorbateMolecule=-1;
  CurrentCationMolecule=SelectRandomMoleculeOfType(CurrentComponent);

  RotationAttempts[CurrentSystem][CurrentComponent]+=1.0; 

  start=Components[CurrentComponent].StartingBead;
  nr_atoms=Components[CurrentComponent].NumberOfAtoms;

  for(i=0;i<nr_atoms;i++)
  {
    cord[i].x=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position.x-
              Cations[CurrentSystem][CurrentCationMolecule].Atoms[start].Position.x;
    cord[i].y=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position.y-
              Cations[CurrentSystem][CurrentCationMolecule].Atoms[start].Position.y;
    cord[i].z=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position.z-
              Cations[CurrentSystem][CurrentCationMolecule].Atoms[start].Position.z;
  }

  RandomArrayRotationMatrix(cord,nr_atoms);

  for(i=0;i<nr_atoms;i++)
  {
    TrialPosition[CurrentSystem][i].x=cord[i].x+Cations[CurrentSystem][CurrentCationMolecule].Atoms[start].Position.x;
    TrialPosition[CurrentSystem][i].y=cord[i].y+Cations[CurrentSystem][CurrentCationMolecule].Atoms[start].Position.y;
    TrialPosition[CurrentSystem][i].z=cord[i].z+Cations[CurrentSystem][CurrentCationMolecule].Atoms[start].Position.z;
  }

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<nr_atoms;i++)
  {
    CFVDWScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFChargeScalingParameter;
  }

  // calculate anisotropic sites
  CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);

  for(i=0;i<nr_atoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  if(Components[CurrentComponent].RestrictEnantionface)
  {
    posA=Components[CurrentComponent].EnantiofaceAtoms[0]->Position;
    if((Components[CurrentComponent].EnantiofaceAtomDefinitions[0][0]==ADSORBATE)&&
       (Components[CurrentComponent].EnantiofaceAtomDefinitions[0][1]==CurrentAdsorbateMolecule))
       posA=TrialPosition[CurrentSystem][Components[CurrentComponent].EnantiofaceAtomDefinitions[0][2]];

    posB=Components[CurrentComponent].EnantiofaceAtoms[1]->Position;
    if((Components[CurrentComponent].EnantiofaceAtomDefinitions[1][0]==ADSORBATE)&&
       (Components[CurrentComponent].EnantiofaceAtomDefinitions[1][1]==CurrentAdsorbateMolecule))
       posB=TrialPosition[CurrentSystem][Components[CurrentComponent].EnantiofaceAtomDefinitions[1][2]];

    posC=Components[CurrentComponent].EnantiofaceAtoms[2]->Position;
    if((Components[CurrentComponent].EnantiofaceAtomDefinitions[2][0]==ADSORBATE)&&
       (Components[CurrentComponent].EnantiofaceAtomDefinitions[2][1]==CurrentAdsorbateMolecule))
       posC=TrialPosition[CurrentSystem][Components[CurrentComponent].EnantiofaceAtomDefinitions[2][2]];

    posD=Components[CurrentComponent].EnantiofaceAtoms[3]->Position;
    if((Components[CurrentComponent].EnantiofaceAtomDefinitions[3][0]==ADSORBATE)&&
       (Components[CurrentComponent].EnantiofaceAtomDefinitions[3][1]==CurrentAdsorbateMolecule))
       posD=TrialPosition[CurrentSystem][Components[CurrentComponent].EnantiofaceAtomDefinitions[3][2]];

    posE=Components[CurrentComponent].EnantiofaceAtoms[4]->Position;
    if((Components[CurrentComponent].EnantiofaceAtomDefinitions[4][0]==ADSORBATE)&&
       (Components[CurrentComponent].EnantiofaceAtomDefinitions[4][1]==CurrentAdsorbateMolecule))
       posE=TrialPosition[CurrentSystem][Components[CurrentComponent].EnantiofaceAtomDefinitions[4][2]];

    if(Components[CurrentComponent].Enantioface!=CheckEnantioFace(posA,posB,posC,posD,posE))
      return 0;
  }

  // check if still a valid position
  if(Components[CurrentComponent].RestrictMoves)
  {
    StartingBead=Components[CurrentComponent].StartingBead;
    pos=TrialPosition[CurrentSystem][StartingBead];

    if(!ValidCartesianPoint(CurrentComponent,pos)) return 0;
  }

  // compute inter-molecular energy differences
  CalculateInterVDWEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeChargeEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterChargeBondDipoleEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateInterBondDipoleBondDipoleEnergyDifferenceCation(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute energy differences framework-cation
  CalculateFrameworkCationVDWEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkCationChargeChargeEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkCationChargeBondDipoleEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference(CurrentCationMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) return 0;

  // compute the energy differenes in Fourier-space
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
     CalculateEwaldFourierCation(TRUE,TRUE,CurrentCationMolecule,0);

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,-1,CurrentCationMolecule);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       UAdsorbatePolarizationNew[CurrentSystem]-UAdsorbatePolarization[CurrentSystem]+
                       (UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       UAdsorbateBackPolarizationNew[CurrentSystem]-UAdsorbateBackPolarization[CurrentSystem]+
                       (UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UCationBackPolarization[CurrentSystem];
  }

  DeltaU=UHostVDWDelta[CurrentSystem]+UAdsorbateVDWDelta[CurrentSystem]+UCationVDWDelta[CurrentSystem]+
         UHostChargeChargeRealDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
         UHostChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
         UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+
         UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
         UDeltaPolarization;

  if(RandomNumber()<exp(-Beta[CurrentSystem]*DeltaU))
  {
    RotationAccepted[CurrentSystem][CurrentComponent]+=1.0;

    UAdsorbateCation[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UHostCation[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UHostCationVDW[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UCationCation[CurrentSystem]+=UCationVDWDelta[CurrentSystem];
    UCationCationVDW[CurrentSystem]+=UCationVDWDelta[CurrentSystem];


    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UCationCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleRealDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                           UCationChargeBondDipoleRealDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                    UCationChargeBondDipoleRealDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(0);
    }

    nr_atoms=Cations[CurrentSystem][CurrentCationMolecule].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassCation(CurrentCationMolecule);

    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

void RotationMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    RotationMoveCation();
  else
    RotationMoveAdsorbate();
}

void PrintRotationStatistics(FILE *FilePtr)
{
  int i,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfRotationMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the rotation move:\n");
    fprintf(FilePtr,"======================================\n");
    for(i=0;i<NumberOfComponents;i++)
      fprintf(FilePtr,"Component [%s] total tried: %lf accepted: %lf (%lf [%%])\n",
        Components[i].Name,
        (double)RotationAttempts[CurrentSystem][i],
        (double)RotationAccepted[CurrentSystem][i], 
        (double)(RotationAccepted[CurrentSystem][i]>(REAL)0.0?
          100.0*RotationAccepted[CurrentSystem][i]/RotationAttempts[CurrentSystem][i]:(REAL)0.0));
    fprintf(FilePtr,"\n");
  }
  else
   fprintf(FilePtr,"Rotation move was OFF for all components\n\n");
}


/*********************************************************************************************************
 * Name       | PartialReinsertionAdsorbateMove                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Partial reinsertion Monte Carlo move of an adsorbate molecule                            *
 * Parameters | -                                                                                        *
 * Note       | The partial regrow move keeps a selected part of a molecule fixed while regenerating the *
 *            | remainder. This selection is described in the definition of the molecule file. Before    *
 *            | calling this function, a component has been randomly chosen. A molecule of this type is  *
 *            | chosen randomly and an attempt is made to regenerate the selected atoms. The internal    *
 *            | configuration of the newly generated part is modified while the configuration of the     *
 *            | fixed part remains unchanged.                                                            *
 *********************************************************************************************************/

int PartialReinsertionAdsorbateMove(void)
{
  int i,d,nr_atoms;
  REAL RosenbluthOld,RosenbluthNew;
  REAL PreFactor;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfAdsorbateMolecules[CurrentSystem]==0) return 0;

  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  if(Components[CurrentComponent].NumberOfConfigMoves==0) return 0;

  // choose a random molecule of this component
  CurrentCationMolecule=-1;
  CurrentAdsorbateMolecule=SelectRandomMoleculeOfType(CurrentComponent);

  PartialReinsertionAttempts[CurrentSystem][CurrentComponent]+=1.0;

  // select which config move
  d=(int)(RandomNumber()*(REAL)Components[CurrentComponent].NumberOfConfigMoves);

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFChargeScalingParameter;
  }

  // calculate the New Rosenbluth factor
  NumberOfBeadsAlreadyPlaced=Components[CurrentComponent].NumberOfUnchangedAtomsConfig[d];
  for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
    BeadsAlreadyPlaced[i]=Components[CurrentComponent].UnchangedAtomsConfig[d][i];
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    NewPosition[CurrentSystem][i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;

  RosenbluthNew=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  // calculate the Old Rosenbluth factor
  NumberOfBeadsAlreadyPlaced=Components[CurrentComponent].NumberOfUnchangedAtomsConfig[d];
  for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
    BeadsAlreadyPlaced[i]=Components[CurrentComponent].UnchangedAtomsConfig[d][i];
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  PartialReinsertionAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  PreFactor=1.0;
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,0);
    PreFactor*=exp(-Beta[CurrentSystem]*(
         UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    PreFactor*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }

  if(RandomNumber()<PreFactor*RosenbluthNew/RosenbluthOld)
  {
    PartialReinsertionAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UAdsorbateBond[CurrentSystem]+=UBondNew[CurrentSystem]-UBondOld[CurrentSystem];
    UAdsorbateUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem]-UUreyBradleyOld[CurrentSystem];
    UAdsorbateBend[CurrentSystem]+=UBendNew[CurrentSystem]-UBendOld[CurrentSystem];
    UAdsorbateBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem]-UBendBendOld[CurrentSystem];
    UAdsorbateInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem]-UInversionBendOld[CurrentSystem];
    UAdsorbateTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem]-UTorsionOld[CurrentSystem];
    UAdsorbateImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem]-UImproperTorsionOld[CurrentSystem];
    UAdsorbateBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem]-UBondBondOld[CurrentSystem];
    UAdsorbateBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem]-UBondBendOld[CurrentSystem];
    UAdsorbateBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem]-UBondTorsionOld[CurrentSystem];
    UAdsorbateBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem]-UBendTorsionOld[CurrentSystem];
    UAdsorbateIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem]-UIntraVDWOld[CurrentSystem];

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UAdsorbateIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem];
      UAdsorbateIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem];
      UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                                 UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                                 UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                          UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                            UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                            UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                            UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                     UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                     UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                     UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                     UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                              UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                              UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                       UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                       UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];


      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }

    nr_atoms=Components[CurrentComponent].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassAdsorbate(CurrentAdsorbateMolecule);

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
          -UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
          -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
          -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
          -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
          -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
          -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
          -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
          -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
           UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization;


    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

/*********************************************************************************************************
 * Name       | PartialReinsertionCationMove                                                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Partial reinsertion Monte Carlo move of a cation molecule                                *
 * Parameters | -                                                                                        *
 * Note       | The partial regrow move keeps a selected part of a molecule fixed while regenerating the *
 *            | remainder. This selection is described in the definition of the molecule file. Before    *
 *            | calling this function, a component has been randomly chosen. A molecule of this type is  *
 *            | chosen randomly and an attempt is made to regenerate the selected atoms. The internal    *
 *            | configuration of the newly generated part is modified while the configuration of the     *
 *            | fixed part remains unchanged.                                                            *
 *********************************************************************************************************/

int PartialReinsertionCationMove(void)
{
  int i,d,nr_atoms;
  REAL RosenbluthOld,RosenbluthNew;
  REAL PreFactor;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfCationMolecules[CurrentSystem]==0) return 0;

  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  if(Components[CurrentComponent].NumberOfConfigMoves==0) return 0;

  // choose a random molecule of this component
  CurrentCationMolecule=SelectRandomMoleculeOfType(CurrentComponent);
  CurrentAdsorbateMolecule=-1;

  PartialReinsertionAttempts[CurrentSystem][CurrentComponent]+=1.0;

  // select which config move
  d=(int)(RandomNumber()*(REAL)Components[CurrentComponent].NumberOfConfigMoves);

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFChargeScalingParameter;
  }

  // calculate the New Rosenbluth factor
  NumberOfBeadsAlreadyPlaced=Components[CurrentComponent].NumberOfUnchangedAtomsConfig[d];
  for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
    BeadsAlreadyPlaced[i]=Components[CurrentComponent].UnchangedAtomsConfig[d][i];
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    NewPosition[CurrentSystem][i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position;
  RosenbluthNew=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  // calculate the Old Rosenbluth factor
  NumberOfBeadsAlreadyPlaced=Components[CurrentComponent].NumberOfUnchangedAtomsConfig[d];
  for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
    BeadsAlreadyPlaced[i]=Components[CurrentComponent].UnchangedAtomsConfig[d][i];
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  PartialReinsertionAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  PreFactor=1.0;
  if(ChargeMethod!=NONE)
  {
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      CalculateEwaldFourierCation(TRUE,TRUE,CurrentCationMolecule,0);
      PreFactor*=exp(-Beta[CurrentSystem]*(
           UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
    }
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       UAdsorbatePolarizationNew[CurrentSystem]-UAdsorbatePolarization[CurrentSystem]+
                       (UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       UAdsorbateBackPolarizationNew[CurrentSystem]-UAdsorbateBackPolarization[CurrentSystem]+
                       (UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UCationBackPolarization[CurrentSystem];
    PreFactor*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }


  if(RandomNumber()<PreFactor*RosenbluthNew/RosenbluthOld)
  {
    PartialReinsertionAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UCationBond[CurrentSystem]+=UBondNew[CurrentSystem]-UBondOld[CurrentSystem];
    UCationUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem]-UUreyBradleyOld[CurrentSystem];
    UCationBend[CurrentSystem]+=UBendNew[CurrentSystem]-UBendOld[CurrentSystem];
    UCationBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem]-UBendBendOld[CurrentSystem];
    UCationInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem]-UInversionBendOld[CurrentSystem];
    UCationTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem]-UTorsionOld[CurrentSystem];
    UCationImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem]-UImproperTorsionOld[CurrentSystem];
    UCationBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem]-UBondBondOld[CurrentSystem];
    UCationBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem]-UBondBendOld[CurrentSystem];
    UCationBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem]-UBondTorsionOld[CurrentSystem];
    UCationBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem]-UBendTorsionOld[CurrentSystem];
    UCationIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem]-UIntraVDWOld[CurrentSystem];

    UCationCation[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UCationCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UHostCation[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];
    UHostCationVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UCationIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem];
      UCationIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem];
      UCationIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UCationCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                           UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                           UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                           UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                           UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                    UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                    UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                    UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                    UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                         UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                         UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                         UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                  UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                  UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                  UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                              UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                              UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                       UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                       UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(0);
    }

    nr_atoms=Components[CurrentComponent].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassCation(CurrentCationMolecule);

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
          -UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
          -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
          -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
          -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
          -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
          -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
          -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
          -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
           UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization;

    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

void PartialReinsertionMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    PartialReinsertionCationMove();
  else
   PartialReinsertionAdsorbateMove();
}

void PrintPartialReinsertionStatistics(FILE *FilePtr)
{
  int i,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfPartialReinsertionMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the partial reinsertion move:\n");
    fprintf(FilePtr,"============================================\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      if(PartialReinsertionAttempts[CurrentSystem][i]>0.0)
      {
        fprintf(FilePtr,"Component [%s] total tried: %lf succesfull growth: %lf (%lf [%%]) accepted: %lf (%lf [%%])\n",
          Components[i].Name,
          (double)PartialReinsertionAttempts[CurrentSystem][i],
          (double)PartialReinsertionAccepted[CurrentSystem][i][0],
          (double)(PartialReinsertionAttempts[CurrentSystem][i]>(REAL)0.0?
            100.0*PartialReinsertionAccepted[CurrentSystem][i][0]/PartialReinsertionAttempts[CurrentSystem][i]:(REAL)0.0),
          (double)PartialReinsertionAccepted[CurrentSystem][i][1],
          (double)(PartialReinsertionAttempts[CurrentSystem][i]>(REAL)0.0?
            100.0*PartialReinsertionAccepted[CurrentSystem][i][1]/PartialReinsertionAttempts[CurrentSystem][i]:(REAL)0.0));
      }
    }
    fprintf(FilePtr,"\n");
  }
  else
   fprintf(FilePtr,"Partial reinsertion move was OFF for all components\n\n");
}



/*********************************************************************************************************
 * Name       | ReinsertionAdsorbateMove                                                                 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Reinsertion Monte Carlo move of an adsorbate molecule                                    *
 * Parameters | -                                                                                        *
 * Note       | The reinsertion move regenerates a randomly chosen molecule at a random position in the  *
 *            | simulation cell. Before calling this function, a component has been randomly chosen.     *
 *********************************************************************************************************/

int ReinsertionAdsorbateMove(void)
{
  int i,nr_atoms,StartingBead;
  REAL RosenbluthOld,RosenbluthNew;
  REAL PreFactor,BiasingWeight=1.0;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfAdsorbateMolecules[CurrentSystem]==0) return 0;

  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentCationMolecule=-1;
  CurrentAdsorbateMolecule=SelectRandomMoleculeOfType(CurrentComponent);

  ReinsertionAttempts[CurrentSystem][CurrentComponent]+=1.0;

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFChargeScalingParameter;
  }

  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNew=GrowMolecule(CBMC_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  ReinsertionAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  NumberOfBeadsAlreadyPlaced=0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_RETRACE_REINSERTION);

  PreFactor=1.0;
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,0);
    PreFactor*=exp(-Beta[CurrentSystem]*(
         UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    PreFactor*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }

  StartingBead=Components[CurrentComponent].StartingBead;
  BiasingWeight=1.0;
  switch(Components[CurrentComponent].Biased)
  {
    case UMBRELLA:
      BiasingWeight=BiasingPotentialUmbrella(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialUmbrella(CurrentComponent,Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position);
      break;
    case RUIZ_MONTERO:
      BiasingWeight=BiasingPotentialRuizMontero(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialRuizMontero(CurrentComponent,Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position);
      break;
    case NO_BIASING:
      BiasingWeight=1.0;
      break;
  }

  if(RandomNumber()<BiasingWeight*PreFactor*RosenbluthNew/RosenbluthOld)
  {
    ReinsertionAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UAdsorbateBond[CurrentSystem]+=UBondNew[CurrentSystem]-UBondOld[CurrentSystem];
    UAdsorbateUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem]-UUreyBradleyOld[CurrentSystem];
    UAdsorbateBend[CurrentSystem]+=UBendNew[CurrentSystem]-UBendOld[CurrentSystem];
    UAdsorbateBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem]-UBendBendOld[CurrentSystem];
    UAdsorbateInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem]-UInversionBendOld[CurrentSystem];
    UAdsorbateTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem]-UTorsionOld[CurrentSystem];
    UAdsorbateImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem]-UImproperTorsionOld[CurrentSystem];
    UAdsorbateBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem]-UBondBondOld[CurrentSystem];
    UAdsorbateBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem]-UBondBendOld[CurrentSystem];
    UAdsorbateBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem]-UBondTorsionOld[CurrentSystem];
    UAdsorbateBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem]-UBendTorsionOld[CurrentSystem];
    UAdsorbateIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem]-UIntraVDWOld[CurrentSystem];

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UAdsorbateIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem];
      UAdsorbateIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem];
      UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                                 UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                                 UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                          UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                            UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                            UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                            UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                     UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                     UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                     UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                     UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                              UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                              UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                       UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                       UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }

    nr_atoms=Components[CurrentComponent].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassAdsorbate(CurrentAdsorbateMolecule);


    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
          -UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
          -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
          -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
          -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
          -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
          -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
          -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
          -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
           UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization;

    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

/********************************************************************************************************* 
 * Name       | ReinsertionAdsorbateMove                                                                 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Reinsertion Monte Carlo move of an adsorbate molecule                                    *
 * Parameters | -                                                                                        *
 * Note       | The reinsertion move regenerates a randomly chosen molecule at a random position in the  *
 *            | simulation cell. Before calling this function, a component has been randomly chosen      *
 *********************************************************************************************************/

int ReinsertionCationMove(void)
{
  int i,nr_atoms;
  REAL RosenbluthOld,RosenbluthNew;
  REAL PreFactor;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfCationMolecules[CurrentSystem]==0) return 0;

  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentCationMolecule=SelectRandomMoleculeOfType(CurrentComponent);
  CurrentAdsorbateMolecule=-1;

  ReinsertionAttempts[CurrentSystem][CurrentComponent]+=1.0;

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFChargeScalingParameter;
  }

  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNew=GrowMolecule(CBMC_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_RETRACE_REINSERTION);

  PreFactor=1.0;
  if(ChargeMethod!=NONE)
  {
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      CalculateEwaldFourierCation(TRUE,TRUE,CurrentCationMolecule,0);
      PreFactor*=exp(-Beta[CurrentSystem]*(
          UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
          UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
          UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
    }
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,-1,CurrentCationMolecule);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       UAdsorbatePolarizationNew[CurrentSystem]-UAdsorbatePolarization[CurrentSystem]+
                       (UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       UAdsorbateBackPolarizationNew[CurrentSystem]-UAdsorbateBackPolarization[CurrentSystem]+
                       (UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UCationBackPolarization[CurrentSystem];
    PreFactor*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }

  ReinsertionAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  if(RandomNumber()<PreFactor*RosenbluthNew/RosenbluthOld)
  {
    ReinsertionAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UCationBond[CurrentSystem]+=UBondNew[CurrentSystem]-UBondOld[CurrentSystem];
    UCationUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem]-UUreyBradleyOld[CurrentSystem];
    UCationBend[CurrentSystem]+=UBendNew[CurrentSystem]-UBendOld[CurrentSystem];
    UCationBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem]-UBendBendOld[CurrentSystem];
    UCationInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem]-UInversionBendOld[CurrentSystem];
    UCationTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem]-UTorsionOld[CurrentSystem];
    UCationImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem]-UImproperTorsionOld[CurrentSystem];
    UCationBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem]-UBondBondOld[CurrentSystem];
    UCationBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem]-UBondBendOld[CurrentSystem];
    UCationBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem]-UBondTorsionOld[CurrentSystem];
    UCationBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem]-UBendTorsionOld[CurrentSystem];
    UCationIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem]-UIntraVDWOld[CurrentSystem];

    UCationCation[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UCationCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UHostCation[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];
    UHostCationVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UCationIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem];
      UCationIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem];
      UCationIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UCationCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                           UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                           UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                           UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                           UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                    UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                    UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                    UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                    UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                         UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                         UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                         UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                  UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                  UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                  UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                              UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                              UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                       UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                       UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(0);
    }

    nr_atoms=Components[CurrentComponent].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassCation(CurrentCationMolecule);

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
          -UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
          -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
          -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
          -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
          -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
          -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
          -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
          -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
           UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization;



    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

void ReinsertionMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    ReinsertionCationMove();
  else
    ReinsertionAdsorbateMove();
}

void PrintReinsertionStatistics(FILE *FilePtr)
{
  int i,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfReinsertionMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the Reinsertion move:\n");
    fprintf(FilePtr,"====================================\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      if(ReinsertionAttempts[CurrentSystem][i]>0.0)
      {
        fprintf(FilePtr,"Component [%s] total tried: %lf succesfull growth: %lf (%lf [%%]) accepted: %lf (%lf [%%])\n",
          Components[i].Name,
          (double)ReinsertionAttempts[CurrentSystem][i],
          (double)ReinsertionAccepted[CurrentSystem][i][0],
          (double)(ReinsertionAttempts[CurrentSystem][i]>(REAL)0.0?
            100.0*ReinsertionAccepted[CurrentSystem][i][0]/ReinsertionAttempts[CurrentSystem][i]:(REAL)0.0),
          (double)ReinsertionAccepted[CurrentSystem][i][1],
          (double)(ReinsertionAttempts[CurrentSystem][i]>(REAL)0.0?
            100.0*ReinsertionAccepted[CurrentSystem][i][1]/ReinsertionAttempts[CurrentSystem][i]:(REAL)0.0));
      }
    }
    fprintf(FilePtr,"\n");
  }
  else
   fprintf(FilePtr,"Reinsertion move was OFF for all components\n\n");
}


/*********************************************************************************************************
 * Name       | ReinsertionInPlaceAdsorbateMove                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Reinsertion Monte Carlo move of an adsorbate molecule at the starting-bead position      *
 * Parameters | -                                                                                        *
 * Note       | The reinsertion-in-place move fully regenerates a randomly chosen molecule at the        *
 *            | position of its starting bead.                                                           *
 *            | Before calling this function, a component has been randomly chosen.                      *
 *********************************************************************************************************/

int ReinsertionInPlaceAdsorbateMove(void)
{
  int i,nr_atoms,StartingBead;
  REAL RosenbluthOld,RosenbluthNew;
  REAL PreFactor,BiasingWeight=1.0;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfAdsorbateMolecules[CurrentSystem]==0) return 0;

  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentCationMolecule=-1;
  CurrentAdsorbateMolecule=SelectRandomMoleculeOfType(CurrentComponent);

  ReinsertionInPlaceAttempts[CurrentSystem][CurrentComponent]+=1.0;

  NumberOfBeadsAlreadyPlaced=0;

  StartingBead=Components[CurrentComponent].StartingBead;
  NewPosition[CurrentSystem][StartingBead]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position;
  RosenbluthNew=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  ReinsertionInPlaceAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFChargeScalingParameter;
  }

  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_PARTIAL_INSERTION);

  PreFactor=1.0;
  if(ChargeMethod!=NONE)
  {
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,0);
      PreFactor*=exp(-Beta[CurrentSystem]*(
          UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
          UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
          UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
    }
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       UAdsorbatePolarizationNew[CurrentSystem]-UAdsorbatePolarization[CurrentSystem]+
                       (UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       UAdsorbateBackPolarizationNew[CurrentSystem]-UAdsorbateBackPolarization[CurrentSystem]+
                       (UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UCationBackPolarization[CurrentSystem];
    PreFactor*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }


  StartingBead=Components[CurrentComponent].StartingBead;
  BiasingWeight=1.0;
  switch(Components[CurrentComponent].Biased)
  {
    case UMBRELLA:
      BiasingWeight=BiasingPotentialUmbrella(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialUmbrella(CurrentComponent,Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position);
      break;
    case RUIZ_MONTERO:
      BiasingWeight=BiasingPotentialRuizMontero(CurrentComponent,TrialPosition[CurrentSystem][StartingBead])/
                    BiasingPotentialRuizMontero(CurrentComponent,Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position);
      break;
    case NO_BIASING:
      BiasingWeight=1.0;
      break;
  }

  if(RandomNumber()<BiasingWeight*PreFactor*RosenbluthNew/RosenbluthOld)
  {
    ReinsertionInPlaceAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UAdsorbateBond[CurrentSystem]+=UBondNew[CurrentSystem]-UBondOld[CurrentSystem];
    UAdsorbateUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem]-UUreyBradleyOld[CurrentSystem];
    UAdsorbateBend[CurrentSystem]+=UBendNew[CurrentSystem]-UBendOld[CurrentSystem];
    UAdsorbateInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem]-UInversionBendOld[CurrentSystem];
    UAdsorbateTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem]-UTorsionOld[CurrentSystem];
    UAdsorbateImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem]-UImproperTorsionOld[CurrentSystem];
    UAdsorbateBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem]-UBondBondOld[CurrentSystem];
    UAdsorbateBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem]-UBondBendOld[CurrentSystem];
    UAdsorbateBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem]-UBondTorsionOld[CurrentSystem];
    UAdsorbateBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem]-UBendTorsionOld[CurrentSystem];
    UAdsorbateIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem]-UIntraVDWOld[CurrentSystem];

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UAdsorbateIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem];
      UAdsorbateIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem];
      UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                                 UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                                 UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                          UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                            UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                            UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                            UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                     UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                     UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                     UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                     UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                              UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                              UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                       UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                       UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }

    nr_atoms=Components[CurrentComponent].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassAdsorbate(CurrentAdsorbateMolecule);

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
          -UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
          -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
          -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
          -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
          -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
          -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
          -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
          -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
           UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization;

    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

/*********************************************************************************************************
 * Name       | ReinsertionInPlaceCationMove                                                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Reinsertion Monte Carlo move of an adsorbate molecule at the starting-bead position      *
 * Parameters | -                                                                                        *
 * Note       | The reinsertion-in-place move fully regenerates a randomly chosen molecule at the        *
 *            | position of its starting bead.                                                           *
 *            | Before calling this function, a component has been randomly chosen.                      *
 *********************************************************************************************************/

int ReinsertionInPlaceCationMove(void)
{
  int i,nr_atoms,StartingBead;
  REAL RosenbluthOld,RosenbluthNew;
  REAL PreFactor;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfCationMolecules[CurrentSystem]==0) return 0;

  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentCationMolecule=SelectRandomMoleculeOfType(CurrentComponent);
  CurrentAdsorbateMolecule=-1;

  ReinsertionInPlaceAttempts[CurrentSystem][CurrentComponent]+=1.0;

  NumberOfBeadsAlreadyPlaced=0;
  StartingBead=Components[CurrentComponent].StartingBead;
  NewPosition[CurrentSystem][StartingBead]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBead].Position;
  RosenbluthNew=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFChargeScalingParameter;
  }

  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_PARTIAL_INSERTION);

  PreFactor=1.0;
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CalculateEwaldFourierCation(TRUE,TRUE,CurrentCationMolecule,0);
    PreFactor*=exp(-Beta[CurrentSystem]*(
        UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
        UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
        UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  ReinsertionInPlaceAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  if(RandomNumber()<PreFactor*RosenbluthNew/RosenbluthOld)
  {
    ReinsertionInPlaceAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UCationBond[CurrentSystem]+=UBondNew[CurrentSystem]-UBondOld[CurrentSystem];
    UCationUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem]-UUreyBradleyOld[CurrentSystem];
    UCationBend[CurrentSystem]+=UBendNew[CurrentSystem]-UBendOld[CurrentSystem];
    UCationBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem]-UBendBendOld[CurrentSystem];
    UCationInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem]-UInversionBendOld[CurrentSystem];
    UCationTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem]-UTorsionOld[CurrentSystem];
    UCationImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem]-UImproperTorsionOld[CurrentSystem];
    UCationBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem]-UBondBondOld[CurrentSystem];
    UCationBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem]-UBondBendOld[CurrentSystem];
    UCationBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem]-UBondTorsionOld[CurrentSystem];
    UCationBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem]-UBendTorsionOld[CurrentSystem];
    UCationIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem]-UIntraVDWOld[CurrentSystem];

    UCationCation[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UCationCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UHostCation[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];
    UHostCationVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UCationIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem];
      UCationIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem];
      UCationIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UCationCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                           UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                           UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                           UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                           UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                    UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                    UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                    UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                    UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                         UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                         UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                         UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                  UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                  UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                  UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];


      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                              UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                              UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                       UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                       UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];


      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(0);
    }

    nr_atoms=Components[CurrentComponent].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassCation(CurrentCationMolecule);

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
          -UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
          -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
          -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
          -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
          -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
          -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
          -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
          -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
           UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization;


    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

void ReinsertionInPlaceMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    ReinsertionInPlaceCationMove();
  else
    ReinsertionInPlaceAdsorbateMove();
}


void PrintReinsertionInPlaceStatistics(FILE *FilePtr)
{
  int i,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfReinsertionInPlaceMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the reinsertion-in-place move:\n");
    fprintf(FilePtr,"=============================================\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      if(ReinsertionInPlaceAttempts[CurrentSystem][i]>0.0)
      {
        fprintf(FilePtr,"Component [%s] total tried: %lf succesfull growth: %lf (%lf [%%]) accepted: %lf (%lf [%%])\n",
          Components[i].Name,
          (double)ReinsertionInPlaceAttempts[CurrentSystem][i],
          (double)ReinsertionInPlaceAccepted[CurrentSystem][i][0],
          (double)(ReinsertionInPlaceAttempts[CurrentSystem][i]>(REAL)0.0?
            100.0*ReinsertionInPlaceAccepted[CurrentSystem][i][0]/ReinsertionInPlaceAttempts[CurrentSystem][i]:(REAL)0.0),
          (double)ReinsertionInPlaceAccepted[CurrentSystem][i][1],
          (double)(ReinsertionInPlaceAttempts[CurrentSystem][i]>(REAL)0.0?
            100.0*ReinsertionInPlaceAccepted[CurrentSystem][i][1]/ReinsertionInPlaceAttempts[CurrentSystem][i]:(REAL)0.0));
      }
      fprintf(FilePtr,"\ttotal        %lf %lf %lf\n",
        (double)TotalTranslationInPlaceAttempts[CurrentSystem][i].x,
        (double)TotalTranslationInPlaceAttempts[CurrentSystem][i].y,
        (double)TotalTranslationInPlaceAttempts[CurrentSystem][i].z);
      fprintf(FilePtr,"\tsuccesfull   %lf %lf %lf\n",
        (double)TotalTranslationInPlaneAccepted[CurrentSystem][i].x,
        (double)TotalTranslationInPlaneAccepted[CurrentSystem][i].y,
        (double)TotalTranslationInPlaneAccepted[CurrentSystem][i].z);
      fprintf(FilePtr,"\taccepted   %lf %lf %lf\n",
        (double)(TotalTranslationInPlaceAttempts[CurrentSystem][i].x>(REAL)0.0?
          TotalTranslationInPlaceAccepted[CurrentSystem][i].x/TotalTranslationInPlaceAttempts[CurrentSystem][i].x:(REAL)0.0),
        (double)(TotalTranslationInPlaceAttempts[CurrentSystem][i].y>(REAL)0.0?
          TotalTranslationInPlaceAccepted[CurrentSystem][i].y/TotalTranslationInPlaceAttempts[CurrentSystem][i].y:(REAL)0.0),
        (double)(TotalTranslationInPlaceAttempts[CurrentSystem][i].z>(REAL)0.0?
          TotalTranslationInPlaceAccepted[CurrentSystem][i].z/TotalTranslationInPlaceAttempts[CurrentSystem][i].z:(REAL)0.0));
      fprintf(FilePtr,"\n");
    }
    fprintf(FilePtr,"\n");
  }
  else
   fprintf(FilePtr,"Reinsertion-in-place move was OFF for all components\n\n");
}


/*********************************************************************************************************
 * Name       | ReinsertionInPlaneAdsorbateMove                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Reinsertion Monte Carlo move of an adsorbate molecule at a plane                         *
 * Parameters | -                                                                                        *
 * Note       | The reinsertion-in-place move fully regenerates a randomly chosen molecule at the a      *
 *            | plane (position of its starting bead).                                                   *
 *            | Before calling this function, a component has been randomly chosen.                      *
 *********************************************************************************************************/

int ReinsertionInPlaneAdsorbateMove(void)
{
  int i,nr_atoms,StartingBead;
  REAL RosenbluthOld,RosenbluthNew;
  REAL PreFactor;
  VECTOR displacement,rotated_displacement;
  REAL vNew,choice;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfAdsorbateMolecules[CurrentSystem]==0) return 0;

  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentCationMolecule=-1;
  CurrentAdsorbateMolecule=SelectRandomMoleculeOfType(CurrentComponent);

  ReinsertionInPlaneAttempts[CurrentSystem][CurrentComponent]+=1.0;


  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].CFChargeScalingParameter;
  }

  NumberOfBeadsAlreadyPlaced=0;

  // calculate a random displacement
  vNew=2.0*RandomNumber()-1.0;

  displacement.x=0.0;
  displacement.y=0.0;
  displacement.z=0.0;

  switch(Components[CurrentComponent].TranslationDirection)
  {
    case XYZ_DIR:
      choice=3.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].x;
      else if(choice<2.0)
        displacement.y=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].y;
      else
        displacement.z=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].z;
      break;
    case XY_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].x;
      else if(choice<2.0)
        displacement.y=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].y;
      break;
    case XZ_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].x;
      else if(choice<2.0)
        displacement.z=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].z;
      break;
    case YZ_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.y=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].y;
      else if(choice<2.0)
        displacement.z=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].z;
      break;
    case X_DIR:
      displacement.x=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].x;
      break;
    case Y_DIR:
      displacement.y=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].y;
      break;
    case Z_DIR:
      displacement.z=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].z;
      break;
  }

  if(fabs(displacement.x)>0.0)
    TranslationInPlaneAttempts[CurrentSystem][CurrentComponent].x+=1.0;
  if(fabs(displacement.y)>0.0)
    TranslationInPlaneAttempts[CurrentSystem][CurrentComponent].y+=1.0;
  if(fabs(displacement.z)>0.0)
    TranslationInPlaneAttempts[CurrentSystem][CurrentComponent].z+=1.0;

  rotated_displacement=RotateVectorAboutX(displacement,BarrierAngle[CurrentSystem].x);
  rotated_displacement=RotateVectorAboutY(rotated_displacement,BarrierAngle[CurrentSystem].y);
  rotated_displacement=RotateVectorAboutZ(rotated_displacement,BarrierAngle[CurrentSystem].z);

  StartingBead=Components[CurrentComponent].StartingBead;
  NewPosition[CurrentSystem][StartingBead].x=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position.x+rotated_displacement.x;
  NewPosition[CurrentSystem][StartingBead].y=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position.y+rotated_displacement.y;
  NewPosition[CurrentSystem][StartingBead].z=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBead].Position.z+rotated_displacement.z;
  RosenbluthNew=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  ReinsertionInPlaneAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_PARTIAL_INSERTION);

  PreFactor=1.0;
  if(ChargeMethod!=NONE)
  {
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,0);
      PreFactor*=exp(-Beta[CurrentSystem]*(
           UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
    }
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    PreFactor*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }

  if(RandomNumber()<PreFactor*RosenbluthNew/RosenbluthOld)
  {
    if(fabs(displacement.x)>0.0)
      TranslationInPlaneAccepted[CurrentSystem][CurrentComponent].x+=1.0;
    if(fabs(displacement.y)>0.0)
      TranslationInPlaneAccepted[CurrentSystem][CurrentComponent].y+=1.0;
    if(fabs(displacement.z)>0.0)
      TranslationInPlaneAccepted[CurrentSystem][CurrentComponent].z+=1.0;

    ReinsertionInPlaneAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UAdsorbateBond[CurrentSystem]+=UBondNew[CurrentSystem]-UBondOld[CurrentSystem];
    UAdsorbateUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem]-UUreyBradleyOld[CurrentSystem];
    UAdsorbateBend[CurrentSystem]+=UBendNew[CurrentSystem]-UBendOld[CurrentSystem];
    UAdsorbateBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem]-UBendBendOld[CurrentSystem];
    UAdsorbateInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem]-UInversionBendOld[CurrentSystem];
    UAdsorbateTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem]-UTorsionOld[CurrentSystem];
    UAdsorbateImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem]-UImproperTorsionOld[CurrentSystem];
    UAdsorbateBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem]-UBondBondOld[CurrentSystem];
    UAdsorbateBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem]-UBondBendOld[CurrentSystem];
    UAdsorbateBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem]-UBondTorsionOld[CurrentSystem];
    UAdsorbateBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem]-UBendTorsionOld[CurrentSystem];
    UAdsorbateIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem]-UIntraVDWOld[CurrentSystem];

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UAdsorbateIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem];
      UAdsorbateIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem];
      UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                                 UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                                 UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                          UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                            UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                            UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                            UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                     UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                     UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                     UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                     UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                              UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                              UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                       UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                       UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];


      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }

    nr_atoms=Components[CurrentComponent].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassAdsorbate(CurrentAdsorbateMolecule);

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
          -UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
          -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
          -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
          -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
          -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
          -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
          -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
          -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
           UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization;

    UTotal[CurrentSystem]+=DeltaU;

  }
  return 0;
}


/*********************************************************************************************************
 * Name       | ReinsertionInPlaneCationMove                                                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Reinsertion Monte Carlo move of an cation molecule at a plane                            *
 * Parameters | -                                                                                        *
 * Note       | The reinsertion-in-place move fully regenerates a randomly chosen molecule at the a      *
 *            | plane (position of its starting bead).                                                   *
 *            | Before calling this function, a component has been randomly chosen.                      *
 *********************************************************************************************************/
int ReinsertionInPlaneCationMove(void)
{
  int i,nr_atoms,StartingBead;
  REAL RosenbluthOld,RosenbluthNew;
  REAL PreFactor;
  VECTOR displacement,rotated_displacement;
  REAL vNew,choice;
  REAL DeltaU;

  // return if the component currently has zero molecules
  if(NumberOfCationMolecules[CurrentSystem]==0) return 0;

  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]==0) return 0;

  // choose a random molecule of this component
  CurrentCationMolecule=SelectRandomMoleculeOfType(CurrentComponent);
  CurrentAdsorbateMolecule=-1;

  ReinsertionInPlaneAttempts[CurrentSystem][CurrentComponent]+=1.0;

  // set Continuous Fraction (CF) atomic scaling-factors, translation is at fixed CF-'lambda'
  // if no CF is used, then these scaling factors are unity
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScaling[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].CFChargeScalingParameter;
  }

  NumberOfBeadsAlreadyPlaced=0;

  // calculate a random displacement
  vNew=2.0*RandomNumber()-1.0;

  displacement.x=0.0;
  displacement.y=0.0;
  displacement.z=0.0;

  switch(Components[CurrentComponent].TranslationDirection)
  {
    case XYZ_DIR:
      choice=3.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].x;
      else if(choice<2.0)
        displacement.y=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].y;
      else
        displacement.z=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].z;
      break;
    case XY_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].x;
      else if(choice<2.0)
        displacement.y=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].y;
      break;
    case XZ_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].x;
      else if(choice<2.0)
        displacement.z=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].z;
      break;
    case YZ_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.y=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].y;
      else if(choice<2.0)
        displacement.z=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].z;
      break;
    case X_DIR:
      displacement.x=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].x;
      break;
    case Y_DIR:
      displacement.y=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].y;
      break;
    case Z_DIR:
      displacement.z=vNew*MaximumTranslationInPlane[CurrentSystem][CurrentComponent].z;
      break;
  }

  if(fabs(displacement.x)>0.0)
    TranslationInPlaneAttempts[CurrentSystem][CurrentComponent].x+=1.0;
  if(fabs(displacement.y)>0.0)
    TranslationInPlaneAttempts[CurrentSystem][CurrentComponent].y+=1.0;
  if(fabs(displacement.z)>0.0)
    TranslationInPlaneAttempts[CurrentSystem][CurrentComponent].z+=1.0;

  rotated_displacement=RotateVectorAboutX(displacement,BarrierAngle[CurrentSystem].x);
  rotated_displacement=RotateVectorAboutY(rotated_displacement,BarrierAngle[CurrentSystem].y);
  rotated_displacement=RotateVectorAboutZ(rotated_displacement,BarrierAngle[CurrentSystem].z);

  StartingBead=Components[CurrentComponent].StartingBead;
  NewPosition[CurrentSystem][StartingBead].x=Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBead].Position.x+rotated_displacement.x;
  NewPosition[CurrentSystem][StartingBead].y=Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBead].Position.y+rotated_displacement.y;
  NewPosition[CurrentSystem][StartingBead].z=Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBead].Position.z+rotated_displacement.z;

  RosenbluthNew=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  ReinsertionInPlaneAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_PARTIAL_INSERTION);

  PreFactor=1.0;
  if(ChargeMethod!=NONE)
  {
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentCationMolecule,0);
      PreFactor*=exp(-Beta[CurrentSystem]*(
           UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
    }
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentCationMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UCationPolarization[CurrentSystem]+
                       UAdsorbatePolarizationNew[CurrentSystem]-UAdsorbatePolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UCationBackPolarization[CurrentSystem]+
                       UAdsorbateBackPolarizationNew[CurrentSystem]-UAdsorbateBackPolarization[CurrentSystem];
    PreFactor*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }

  if(RandomNumber()<PreFactor*RosenbluthNew/RosenbluthOld)
  {
    if(fabs(displacement.x)>0.0)
      TranslationInPlaneAccepted[CurrentSystem][CurrentComponent].x+=1.0;
    if(fabs(displacement.y)>0.0)
      TranslationInPlaneAccepted[CurrentSystem][CurrentComponent].y+=1.0;
    if(fabs(displacement.z)>0.0)
      TranslationInPlaneAccepted[CurrentSystem][CurrentComponent].z+=1.0;

    ReinsertionInPlaneAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UCationBond[CurrentSystem]+=UBondNew[CurrentSystem]-UBondOld[CurrentSystem];
    UCationUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem]-UUreyBradleyOld[CurrentSystem];
    UCationBend[CurrentSystem]+=UBendNew[CurrentSystem]-UBendOld[CurrentSystem];
    UCationBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem]-UBendBendOld[CurrentSystem];
    UCationInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem]-UInversionBendOld[CurrentSystem];
    UCationTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem]-UTorsionOld[CurrentSystem];
    UCationImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem]-UImproperTorsionOld[CurrentSystem];
    UCationBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem]-UBondBondOld[CurrentSystem];
    UCationBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem]-UBondBendOld[CurrentSystem];
    UCationBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem]-UBondTorsionOld[CurrentSystem];
    UCationBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem]-UBendTorsionOld[CurrentSystem];
    UCationIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem]-UIntraVDWOld[CurrentSystem];

    UCationCation[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UCationCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UHostCation[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];
    UHostCationVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UCationIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem];
      UCationIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem];
      UCationIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UCationCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                           UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                           UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                           UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                           UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                    UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                    UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                    UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                    UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                         UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                         UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                         UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                  UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                  UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                  UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                              UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                              UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                       UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                       UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];


      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(0);
    }

    nr_atoms=Components[CurrentComponent].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    }

    UpdateGroupCenterOfMassCation(CurrentCationMolecule);

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UCationVDWNew[CurrentSystem]+UAdsorbateVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UCationChargeChargeNew[CurrentSystem]+UAdsorbateChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UCationChargeBondDipoleNew[CurrentSystem]+UAdsorbateChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UCationBondDipoleBondDipoleNew[CurrentSystem]+UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
          -UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
          -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
          -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
          -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
          -UCationVDWOld[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
          -UCationChargeChargeOld[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
          -UCationChargeBondDipoleOld[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
          -UCationBondDipoleBondDipoleOld[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
           UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization;

    UTotal[CurrentSystem]+=DeltaU;

  }
  return 0;
}

void ReinsertionInPlaneMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    ReinsertionInPlaneCationMove();
  else
    ReinsertionInPlaneAdsorbateMove();
}


void OptimizeTranslationInPlaneAcceptence(void)
{
  int i;
  VECTOR ratio,vandr;

  for(i=0;i<NumberOfComponents;i++)
  {
    if(TranslationInPlaneAttempts[CurrentSystem][i].x>0.0)
      ratio.x=TranslationInPlaneAccepted[CurrentSystem][i].x/TranslationInPlaneAttempts[CurrentSystem][i].x;
    else
      ratio.x=0.0;
    if(TranslationInPlaneAttempts[CurrentSystem][i].y>0.0)
      ratio.y=TranslationInPlaneAccepted[CurrentSystem][i].y/TranslationInPlaneAttempts[CurrentSystem][i].y;
    else
      ratio.y=0.0;
    if(TranslationInPlaneAttempts[CurrentSystem][i].z>0.0)
      ratio.z=TranslationInPlaneAccepted[CurrentSystem][i].z/TranslationInPlaneAttempts[CurrentSystem][i].z;
    else
      ratio.z=0.0;

    vandr.x=ratio.x/TargetAccRatioTranslation;
    if(vandr.x>1.5) vandr.x=1.5;
    else if(vandr.x<0.5) vandr.x=0.5;
    MaximumTranslationInPlane[CurrentSystem][i].x*=vandr.x;
    if(MaximumTranslationInPlane[CurrentSystem][i].x<0.01)
       MaximumTranslationInPlane[CurrentSystem][i].x=0.01;
    if(MaximumTranslationInPlane[CurrentSystem][i].x>CutOffVDW)
       MaximumTranslationInPlane[CurrentSystem][i].x=CutOffVDW;

    vandr.y=ratio.y/TargetAccRatioTranslation;
    if(vandr.y>1.5) vandr.y=1.5;
    else if(vandr.y<0.5) vandr.y=0.5;
    MaximumTranslationInPlane[CurrentSystem][i].y*=vandr.y;
    if(MaximumTranslationInPlane[CurrentSystem][i].y<0.01)
       MaximumTranslationInPlane[CurrentSystem][i].y=0.01;
    if(MaximumTranslationInPlane[CurrentSystem][i].y>CutOffVDW)
       MaximumTranslationInPlane[CurrentSystem][i].y=CutOffVDW;

    vandr.z=ratio.z/TargetAccRatioTranslation;
     if(vandr.z>1.5) vandr.z=1.5;
    else if(vandr.z<0.5) vandr.z=0.5;
    MaximumTranslationInPlane[CurrentSystem][i].z*=vandr.z;
    if(MaximumTranslationInPlane[CurrentSystem][i].z<0.01)
       MaximumTranslationInPlane[CurrentSystem][i].z=0.01;
    if(MaximumTranslationInPlane[CurrentSystem][i].z>CutOffVDW)
       MaximumTranslationInPlane[CurrentSystem][i].z=CutOffVDW;

    TotalTranslationInPlaneAttempts[CurrentSystem][i].x+=TranslationInPlaneAttempts[CurrentSystem][i].x;
    TotalTranslationInPlaneAccepted[CurrentSystem][i].x+=TranslationInPlaneAccepted[CurrentSystem][i].x;
    TranslationInPlaneAttempts[CurrentSystem][i].x=TranslationInPlaneAccepted[CurrentSystem][i].x=0.0;

    TotalTranslationInPlaneAttempts[CurrentSystem][i].y+=TranslationInPlaneAttempts[CurrentSystem][i].y;
    TotalTranslationInPlaneAccepted[CurrentSystem][i].y+=TranslationInPlaneAccepted[CurrentSystem][i].y;
    TranslationInPlaneAttempts[CurrentSystem][i].y=TranslationInPlaneAccepted[CurrentSystem][i].y=0.0;

    TotalTranslationInPlaneAttempts[CurrentSystem][i].z+=TranslationInPlaneAttempts[CurrentSystem][i].z;
    TotalTranslationInPlaneAccepted[CurrentSystem][i].z+=TranslationInPlaneAccepted[CurrentSystem][i].z;
    TranslationInPlaneAttempts[CurrentSystem][i].z=TranslationInPlaneAccepted[CurrentSystem][i].z=0.0;
  }
}

void PrintReinsertionInPlaneStatistics(FILE *FilePtr)
{
  int i,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfReinsertionInPlaneMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the regrow in plane move:\n");
    fprintf(FilePtr,"========================================\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      if(ReinsertionInPlaneAttempts[CurrentSystem][i]>0.0)
      {
        fprintf(FilePtr,"Component [%s] total tried: %lf succesfull growth: %lf (%lf [%%]) accepted: %lf (%lf [%%])\n",
          Components[i].Name,
          (double)ReinsertionInPlaneAttempts[CurrentSystem][i],
          (double)ReinsertionInPlaneAccepted[CurrentSystem][i][0],
          (double)(ReinsertionInPlaneAttempts[CurrentSystem][i]>(REAL)0.0?
            100.0*ReinsertionInPlaneAccepted[CurrentSystem][i][0]/ReinsertionInPlaneAttempts[CurrentSystem][i]:(REAL)0.0),
          (double)ReinsertionInPlaneAccepted[CurrentSystem][i][1],
          (double)(ReinsertionInPlaneAttempts[CurrentSystem][i]>(REAL)0.0?
            100.0*ReinsertionInPlaneAccepted[CurrentSystem][i][1]/ReinsertionInPlaneAttempts[CurrentSystem][i]:(REAL)0.0));
      }
      fprintf(FilePtr,"\ttotal        %lf %lf %lf\n",
        (double)TotalTranslationInPlaneAttempts[CurrentSystem][i].x,
        (double)TotalTranslationInPlaneAttempts[CurrentSystem][i].y,
        (double)TotalTranslationInPlaneAttempts[CurrentSystem][i].z);
      fprintf(FilePtr,"\tsuccesfull   %lf %lf %lf\n",
        (double)TotalTranslationInPlaneAccepted[CurrentSystem][i].x,
        (double)TotalTranslationInPlaneAccepted[CurrentSystem][i].y,
        (double)TotalTranslationInPlaneAccepted[CurrentSystem][i].z);
      fprintf(FilePtr,"\taccepted   %lf %lf %lf\n",
        (double)(TotalTranslationInPlaneAttempts[CurrentSystem][i].x>(REAL)0.0?
          TotalTranslationInPlaneAccepted[CurrentSystem][i].x/TotalTranslationInPlaneAttempts[CurrentSystem][i].x:(REAL)0.0),
        (double)(TotalTranslationInPlaneAttempts[CurrentSystem][i].y>(REAL)0.0?
          TotalTranslationInPlaneAccepted[CurrentSystem][i].y/TotalTranslationInPlaneAttempts[CurrentSystem][i].y:(REAL)0.0),
        (double)(TotalTranslationInPlaneAttempts[CurrentSystem][i].z>(REAL)0.0?
          TotalTranslationInPlaneAccepted[CurrentSystem][i].z/TotalTranslationInPlaneAttempts[CurrentSystem][i].z:(REAL)0.0));
      fprintf(FilePtr,"\tdisplacement %lf %lf %lf\n",
        (double)MaximumTranslationInPlane[CurrentSystem][i].x,
        (double)MaximumTranslationInPlane[CurrentSystem][i].y,
        (double)MaximumTranslationInPlane[CurrentSystem][i].z);
      fprintf(FilePtr,"\n");
    }
    fprintf(FilePtr,"\n");
  }
  else
   fprintf(FilePtr,"Reinsertion-in-plane move was OFF for all components\n\n");
}


// Identity-switch Monte-Carlo move in the grand-canonical ensemble (Adsorbate-version)
// A molecule of type 'A' is randomly selected and an attempt is made to change its identity to 'B'
// The move is called semi-grand ensemble, but is also a special case of the Gibbs ensemble.
// The move is called 'swotch' in original paper of M. G. Martin and J. I. Siepmann
// JACS, 1997, 8921-8924

int IdentityChangeAdsorbateMove(void)
{
  int i,d,nr_atoms,type;
  int OldComponent,NewComponent;
  int StartingBeadOld,StartingBeadNew;
  REAL RosenbluthOld,RosenbluthNew,TailCorrectionDiff;
  REAL RosenbluthIdealOld,RosenbluthIdealNew;
  REAL PreFactor,PartialFugacityNew,PartialFugacityOld;
  REAL DeltaU;

  // Choose the 'Old' and 'New'-component
  OldComponent=CurrentComponent;
  d=(int)(RandomNumber()*(REAL)Components[OldComponent].NumberOfIdentityChanges);
  NewComponent=Components[OldComponent].IdentityChanges[d];

  if(Components[NewComponent].ExtraFrameworkMolecule) return 0;

  // garanty detailed balance
  if(RandomNumber()<0.5)
  {
    d=OldComponent;
    OldComponent=NewComponent;
    NewComponent=d;
  }

  // return if no molecules of the 'Old'-component are currently present
  if(Components[OldComponent].NumberOfMolecules[CurrentSystem]<=Components[OldComponent].CFMoleculePresent[CurrentSystem]?1:0) return 0;

  // choose a random molecule of the 'Old'-component
  CurrentCationMolecule=-1;
  CurrentAdsorbateMolecule=SelectRandomMoleculeOfTypeExcludingFractionalMolecule(OldComponent);

  // register an attempt to change the 'Old'-component to the 'New'-component
  IdentityChangeAttempts[CurrentSystem][OldComponent][NewComponent]+=1.0;

  // set Continuous Fraction (CF) atomic scaling-factors to unity
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  StartingBeadOld=Components[OldComponent].StartingBead;
  StartingBeadNew=Components[NewComponent].StartingBead;

  // the starting position of the 'Old'-component is taken as the starting position for the 'New'-component
  NewPosition[CurrentSystem][StartingBeadNew]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBeadOld].Position;

  // grow the 'New'-component
  NumberOfBeadsAlreadyPlaced=0;
  CurrentComponent=NewComponent;

  // select which config move
/*
  d=(int)(RandomNumber()*(REAL)Components[CurrentComponent].NumberOfIdentityConfigMoves);
  NumberOfBeadsAlreadyPlaced=Components[CurrentComponent].NumberOfUnchangedAtomsIdentityConfig[d];
  for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
    BeadsAlreadyPlaced[i]=Components[CurrentComponent].UnchangedAtomsIdentityConfig[d][i];
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    NewPosition[CurrentSystem][i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
    //TrialAnisotropicPosition[CurrentSystem][i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].AnisotropicPosition;
  }
*/


  RosenbluthNew=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[NewComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(NewPosition[CurrentSystem][i]))
      return 0;
  }

  // retrace the 'Old'-component
  NumberOfBeadsAlreadyPlaced=0;
  CurrentComponent=OldComponent;

/*
  NumberOfBeadsAlreadyPlaced=Components[CurrentComponent].NumberOfUnchangedAtomsIdentityConfig[d];
  for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
    BeadsAlreadyPlaced[i]=Components[CurrentComponent].UnchangedAtomsIdentityConfig[d][i];
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    NewPosition[CurrentSystem][i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
    //TrialAnisotropicPosition[CurrentSystem][i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].AnisotropicPosition;
  }
*/
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_PARTIAL_INSERTION);


  // compute the tail-correction difference
  TailCorrectionDiff=TailMolecularEnergyDifferenceAddRemove(NewComponent,OldComponent);
  RosenbluthNew*=exp(-Beta[CurrentSystem]*TailCorrectionDiff);

  // compute the Ewald-correction difference
  PreFactor=1.0;
  if(ChargeMethod!=NONE)
  {
    CurrentComponent=NewComponent;
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,0);
      PreFactor*=exp(-Beta[CurrentSystem]*(
           UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
    }
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    PreFactor*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }

  // register an succesfull growth/retrace before acceptance
  IdentityChangeAccepted[CurrentSystem][OldComponent][NewComponent][0]+=1.0;

  RosenbluthIdealNew=Components[NewComponent].IdealGasRosenbluthWeight[CurrentSystem];
  RosenbluthIdealOld=Components[OldComponent].IdealGasRosenbluthWeight[CurrentSystem];

  PartialFugacityNew=Components[NewComponent].FugacityCoefficient[CurrentSystem]*Components[NewComponent].PartialPressure[CurrentSystem];
  PartialFugacityOld=Components[OldComponent].FugacityCoefficient[CurrentSystem]*Components[OldComponent].PartialPressure[CurrentSystem];

  // acceptance-rule
  if(RandomNumber()<PreFactor*((RosenbluthNew/RosenbluthIdealNew)*PartialFugacityNew*
                    (REAL)Components[OldComponent].NumberOfMolecules[CurrentSystem]/
                   ((RosenbluthOld/RosenbluthIdealOld)*PartialFugacityOld*
                    ((REAL)Components[NewComponent].NumberOfMolecules[CurrentSystem]+(REAL)1.0))))
  {
    // register an succesfull growth/retrace after acceptance
    IdentityChangeAccepted[CurrentSystem][OldComponent][NewComponent][1]+=1.0;

    // register the changes in the energies
    UAdsorbateBond[CurrentSystem]+=UBondNew[CurrentSystem]-UBondOld[CurrentSystem];
    UAdsorbateUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem]-UUreyBradleyOld[CurrentSystem];
    UAdsorbateBend[CurrentSystem]+=UBendNew[CurrentSystem]-UBendOld[CurrentSystem];
    UAdsorbateBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem]-UBendBendOld[CurrentSystem];
    UAdsorbateInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem]-UInversionBendOld[CurrentSystem];
    UAdsorbateTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem]-UTorsionOld[CurrentSystem];
    UAdsorbateImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem]-UImproperTorsionOld[CurrentSystem];
    UAdsorbateBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem]-UBondBondOld[CurrentSystem];
    UAdsorbateBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem]-UBondBendOld[CurrentSystem];
    UAdsorbateBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem]-UBondTorsionOld[CurrentSystem];
    UAdsorbateBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem]-UBendTorsionOld[CurrentSystem];
    UAdsorbateIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem]-UIntraVDWOld[CurrentSystem];

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    UTailCorrection[CurrentSystem]+=TailCorrectionDiff;

    if(ChargeMethod!=NONE)
    {
      UAdsorbateIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem];
      UAdsorbateIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem];
      UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                                 UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                                 UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                          UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                            UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                            UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                            UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                     UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                     UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                     UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                     UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                              UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                              UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                       UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                       UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];

      NetChargeAdsorbates[CurrentSystem]+=NetChargeAdsorbateDelta;
      NetChargeSystem[CurrentSystem]+=NetChargeAdsorbateDelta;

      // register the changes in the stored structure factors
      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }

    // register the changes in the types for the 'Old'-component atoms
    nr_atoms=Components[OldComponent].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      type=Components[OldComponent].Type[i];
      NumberOfPseudoAtomsType[CurrentSystem][type]--;
    }

    // overwrite the 'Old'-component with the 'New'-component
    nr_atoms=Components[NewComponent].NumberOfAtoms;
    Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].NumberOfAtoms=nr_atoms;
    Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Type=NewComponent;

    // realloc memory (the size of the molecule could be different)
    Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms=(ATOM*)
          realloc(Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms,nr_atoms*sizeof(ATOM));
    if(Components[CurrentComponent].NumberOfGroups>0)
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Groups=(GROUP*)
          realloc(Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Groups,
                  Components[CurrentComponent].NumberOfGroups*sizeof(GROUP));

    for(i=0;i<nr_atoms;i++)
    {
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
      type=Components[NewComponent].Type[i];
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Type=type;
      Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Charge=Components[NewComponent].Charge[i];
      NumberOfPseudoAtomsType[CurrentSystem][type]++;
    }

    // register the changes in the amount of molecules of each type
    Components[NewComponent].NumberOfMolecules[CurrentSystem]++;
    Components[OldComponent].NumberOfMolecules[CurrentSystem]--;

    NumberOfAtomsPerSystem[CurrentSystem]-=Components[OldComponent].NumberOfAtoms;
    NumberOfChargesPerSystem[CurrentSystem]-=Components[OldComponent].NumberOfCharges;
    NumberOfBondDipolesPerSystem[CurrentSystem]-=Components[OldComponent].NumberOfBondDipoles;

    NumberOfAtomsPerSystem[CurrentSystem]+=Components[NewComponent].NumberOfAtoms;
    NumberOfChargesPerSystem[CurrentSystem]+=Components[NewComponent].NumberOfCharges;
    NumberOfBondDipolesPerSystem[CurrentSystem]+=Components[NewComponent].NumberOfBondDipoles;

    // update the center of mass of the molecule
    UpdateGroupCenterOfMassAdsorbate(CurrentAdsorbateMolecule);

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
          -UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
          -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
          -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
          -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
          -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
          -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
          -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
          -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
           UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization+TailCorrectionDiff;

    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

// Identity-switch Monte-Carlo move in the grand-canonical ensemble (Cation-version)
// A molecule of type 'A' is randomly selected and an attempt is made to change its identity to 'B'
// The move is called semi-grand ensemble, but is also a special case of the Gibbs ensemble.
// The move is called 'swotch' in original paper of M. G. Martin and J. I. Siepmann
// JACS, 1997, 8921-8924

int IdentityChangeCationMove(void)
{
  int i,d,nr_atoms,type;
  int OldComponent,NewComponent;
  int StartingBeadOld,StartingBeadNew;
  REAL RosenbluthOld,RosenbluthNew,UTailOld,UTailNew;
  REAL RosenbluthIdealOld,RosenbluthIdealNew;
  REAL PreFactor,PartialFugacityNew,PartialFugacityOld;
  REAL DeltaU;

  OldComponent=CurrentComponent;
  d=(int)(RandomNumber()*(REAL)Components[OldComponent].NumberOfIdentityChanges);
  NewComponent=Components[OldComponent].IdentityChanges[d];

  if(!Components[NewComponent].ExtraFrameworkMolecule) return 0;

  // garanty detailed balance
  if(RandomNumber()<0.5)
  {
    d=OldComponent;
    OldComponent=NewComponent;
    NewComponent=d;
  }

  // return if no molecules of the 'Old'-component are currently present
  if(Components[OldComponent].NumberOfMolecules[CurrentSystem]<=Components[OldComponent].CFMoleculePresent[CurrentSystem]?1:0) return 0;

  // choose a random molecule of this component
  CurrentAdsorbateMolecule=-1;
  CurrentCationMolecule=SelectRandomMoleculeOfTypeExcludingFractionalMolecule(OldComponent);

  // register an attempt to change the 'Old'-component to the 'New'-component
  IdentityChangeAttempts[CurrentSystem][OldComponent][NewComponent]+=1.0;

  // set Continuous Fraction (CF) atomic scaling-factors to unity
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  StartingBeadOld=Components[OldComponent].StartingBead;
  StartingBeadNew=Components[NewComponent].StartingBead;

  // the starting position of the 'Old'-component is taken as the starting position
  NewPosition[CurrentSystem][StartingBeadNew]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBeadOld].Position;

  // grow the 'New'-component
  NumberOfBeadsAlreadyPlaced=0;
  CurrentComponent=NewComponent;
  RosenbluthNew=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[NewComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(NewPosition[CurrentSystem][i]))
      return 0;
  }

  // compute the tail-correction difference
  UTailNew=TailMolecularEnergyDifferenceAdd();
  RosenbluthNew*=exp(-Beta[CurrentSystem]*UTailNew);

  // retrace the 'Old'-component
  NumberOfBeadsAlreadyPlaced=0;
  CurrentComponent=OldComponent;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_PARTIAL_INSERTION);

  // compute the tail-correction difference
  UTailOld=TailMolecularEnergyDifferenceRemove();
  RosenbluthOld*=exp(-Beta[CurrentSystem]*UTailOld);

  PreFactor=1.0;
  if(ChargeMethod!=NONE)
  {
    CurrentComponent=NewComponent;
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      CalculateEwaldFourierCation(TRUE,TRUE,CurrentCationMolecule,0);
      PreFactor*=exp(-Beta[CurrentSystem]*(
         UHostCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+
         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
    }
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,-1,CurrentCationMolecule);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       UAdsorbatePolarizationNew[CurrentSystem]-UAdsorbatePolarization[CurrentSystem]+
                       (UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       UAdsorbateBackPolarizationNew[CurrentSystem]-UAdsorbateBackPolarization[CurrentSystem]+
                       (UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UCationBackPolarization[CurrentSystem];

    PreFactor*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }


  // register an succesfull growth/retrace before acceptance
  IdentityChangeAccepted[CurrentSystem][OldComponent][NewComponent][0]+=1.0;

  RosenbluthIdealNew=Components[NewComponent].IdealGasRosenbluthWeight[CurrentSystem];
  RosenbluthIdealOld=Components[OldComponent].IdealGasRosenbluthWeight[CurrentSystem];

  PartialFugacityNew=Components[NewComponent].FugacityCoefficient[CurrentSystem]*Components[NewComponent].PartialPressure[CurrentSystem];
  PartialFugacityOld=Components[OldComponent].FugacityCoefficient[CurrentSystem]*Components[OldComponent].PartialPressure[CurrentSystem];

  // acceptance-rule
  if(RandomNumber()<PreFactor*((RosenbluthNew/RosenbluthIdealNew)*PartialFugacityNew*
                               (REAL)Components[OldComponent].NumberOfMolecules[CurrentSystem]/
                              ((RosenbluthOld/RosenbluthIdealOld)*PartialFugacityOld*
                               ((REAL)Components[NewComponent].NumberOfMolecules[CurrentSystem]+(REAL)1.0))))
  {
    // register an succesfull growth/retrace after acceptance
    IdentityChangeAccepted[CurrentSystem][OldComponent][NewComponent][1]+=1.0;

    // register the changes in the energies
    UCationBond[CurrentSystem]+=UBondNew[CurrentSystem]-UBondOld[CurrentSystem];
    UCationUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem]-UUreyBradleyOld[CurrentSystem];
    UCationBend[CurrentSystem]+=UBendNew[CurrentSystem]-UBendOld[CurrentSystem];
    UCationBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem]-UBendBendOld[CurrentSystem];
    UCationInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem]-UInversionBendOld[CurrentSystem];
    UCationTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem]-UTorsionOld[CurrentSystem];
    UCationImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem]-UImproperTorsionOld[CurrentSystem];
    UCationBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem]-UBondBondOld[CurrentSystem];
    UCationBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem]-UBondBendOld[CurrentSystem];
    UCationBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem]-UBondTorsionOld[CurrentSystem];
    UCationBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem]-UBendTorsionOld[CurrentSystem];
    UCationIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem]-UIntraVDWOld[CurrentSystem];

    UCationCation[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UCationCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem]-UCationVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem]-UAdsorbateVDWOld[CurrentSystem];
    UHostCation[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];
    UHostCationVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem]-UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];

    UTailCorrection[CurrentSystem]+=UTailNew-UTailOld;

    if(ChargeMethod!=NONE)
    {
      UCationIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem];
      UCationIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem];
      UCationIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UCationCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                           UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                           UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                           UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                           UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                    UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                    UCationChargeChargeNew[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]+
                                    UCationChargeBondDipoleNew[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]+
                                    UCationBondDipoleBondDipoleNew[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                         UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                         UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                         UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                  UHostChargeChargeNew[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]+
                                  UHostChargeBondDipoleNew[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]+
                                  UHostBondDipoleBondDipoleNew[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                              UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                              UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateChargeChargeNew[CurrentSystem]-UAdsorbateChargeChargeOld[CurrentSystem]+
                                       UAdsorbateChargeBondDipoleNew[CurrentSystem]-UAdsorbateChargeBondDipoleOld[CurrentSystem]+
                                       UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      NetChargeCations[CurrentSystem]+=NetChargeCationDelta;
      NetChargeSystem[CurrentSystem]+=NetChargeCationDelta;

      // register the changes in the stored structure factors
      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(0);
    }

    // register the changes in the types for the 'Old'-component atoms
    nr_atoms=Components[OldComponent].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      type=Components[OldComponent].Type[i];
      NumberOfPseudoAtomsType[CurrentSystem][type]--;
    }

    // overwrite the 'Old'-component with the 'New'-component
    nr_atoms=Components[NewComponent].NumberOfAtoms;
    Cations[CurrentSystem][CurrentCationMolecule].NumberOfAtoms=nr_atoms;
    Cations[CurrentSystem][CurrentCationMolecule].Type=NewComponent;

    // realloc memory (the size of the molecule could be different)
    Cations[CurrentSystem][CurrentCationMolecule].Atoms=(ATOM*)
          realloc(Cations[CurrentSystem][CurrentCationMolecule].Atoms,nr_atoms*sizeof(ATOM));
    if(Components[CurrentComponent].NumberOfGroups>0)
      Cations[CurrentSystem][CurrentCationMolecule].Groups=(GROUP*)
          realloc(Cations[CurrentSystem][CurrentCationMolecule].Groups,
                  Components[CurrentComponent].NumberOfGroups*sizeof(GROUP));

    for(i=0;i<nr_atoms;i++)
    {
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position=TrialPosition[CurrentSystem][i];
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];

      type=Components[NewComponent].Type[i];
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Type=type;
      Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Charge=Components[NewComponent].Charge[i];
      NumberOfPseudoAtomsType[CurrentSystem][type]++;
    }

    // register the changes in the amount of molecules of each type
    Components[NewComponent].NumberOfMolecules[CurrentSystem]++;
    Components[OldComponent].NumberOfMolecules[CurrentSystem]--;

    NumberOfAtomsPerSystem[CurrentSystem]-=Components[OldComponent].NumberOfAtoms;
    NumberOfChargesPerSystem[CurrentSystem]-=Components[OldComponent].NumberOfCharges;
    NumberOfBondDipolesPerSystem[CurrentSystem]-=Components[OldComponent].NumberOfBondDipoles;

    NumberOfAtomsPerSystem[CurrentSystem]+=Components[NewComponent].NumberOfAtoms;
    NumberOfChargesPerSystem[CurrentSystem]+=Components[NewComponent].NumberOfCharges;
    NumberOfBondDipolesPerSystem[CurrentSystem]+=Components[NewComponent].NumberOfBondDipoles;

    // update the center of mass of the molecule
    UpdateGroupCenterOfMassCation(CurrentCationMolecule);

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
          -UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
          -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
          -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
          -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
          -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
          -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
          -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
          -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
           UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization+UTailNew-UTailOld;

    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

void IdentityChangeMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    IdentityChangeCationMove();
  else
    IdentityChangeAdsorbateMove();
}

void PrintIdentityChangeStatistics(FILE *FilePtr)
{
  int i,j,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfIdentityChangeMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the identity change move:\n");
    fprintf(FilePtr,"======================================\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      for(j=0;j<NumberOfComponents;j++)
      {
        if(IdentityChangeAttempts[CurrentSystem][i][j]>0.0)
        {
          fprintf(FilePtr,"Component [%s]->[%s] total tried: %lf succesfull growth: %lf (%lf [%%]) accepted: %lf (%lf [%%])\n",
            Components[i].Name,
            Components[j].Name,
            (double)IdentityChangeAttempts[CurrentSystem][i][j],
            (double)IdentityChangeAccepted[CurrentSystem][i][j][0],
            (double)(IdentityChangeAttempts[CurrentSystem][i][j]>(REAL)0.0?
              100.0*IdentityChangeAccepted[CurrentSystem][i][j][0]/IdentityChangeAttempts[CurrentSystem][i][j]:(REAL)0.0),
            (double)IdentityChangeAccepted[CurrentSystem][i][j][1],
            (double)(IdentityChangeAttempts[CurrentSystem][i][j]>(REAL)0.0?
              100.0*IdentityChangeAccepted[CurrentSystem][i][j][1]/IdentityChangeAttempts[CurrentSystem][i][j]:(REAL)0.0));
        }
      }
    }
    fprintf(FilePtr,"\n");
  }
  else
   fprintf(FilePtr,"Identity change move was OFF for all components\n\n");

}


// Swap Monte-Carlo move
//
//
//
//

int SwapAddAdsorbateMove(void)
{
  int i;
  REAL RosenbluthNew,PartialFugacity,UTailNew;
  REAL RosenbluthIdealNew;
  REAL DeltaU;

  SwapAddAttempts[CurrentSystem][CurrentComponent]+=1.0;
  CurrentAdsorbateMolecule=NumberOfAdsorbateMolecules[CurrentSystem];
  CurrentCationMolecule=NumberOfCationMolecules[CurrentSystem];

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be added here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNew=GrowMolecule(CBMC_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }


  SwapAddAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  UTailNew=TailMolecularEnergyDifferenceAdd();
  RosenbluthNew*=exp(-Beta[CurrentSystem]*UTailNew);

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    RosenbluthNew*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }


  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CalculateEwaldFourierAdsorbate(TRUE,FALSE,NumberOfAdsorbateMolecules[CurrentSystem],0);

    RosenbluthNew*=exp(-Beta[CurrentSystem]*(
        UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
        UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
        UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  // get partial pressure for this component
  PartialFugacity=Components[CurrentComponent].FugacityCoefficient[CurrentSystem]*
                  Components[CurrentComponent].PartialPressure[CurrentSystem];

  RosenbluthIdealNew=Components[CurrentComponent].IdealGasRosenbluthWeight[CurrentSystem];

  // acceptence rule
  if(RandomNumber()<((RosenbluthNew/RosenbluthIdealNew)*Beta[CurrentSystem]*PartialFugacity*Volume[CurrentSystem]/
                    (1.0+Components[CurrentComponent].NumberOfMolecules[CurrentSystem])))
  {
    SwapAddAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UAdsorbateBond[CurrentSystem]+=UBondNew[CurrentSystem];
    UAdsorbateUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem];
    UAdsorbateBend[CurrentSystem]+=UBendNew[CurrentSystem];
    UAdsorbateBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem];
    UAdsorbateInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem];
    UAdsorbateTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem];
    UAdsorbateImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem];
    UAdsorbateBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem];
    UAdsorbateBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem];
    UAdsorbateBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem];
    UAdsorbateBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem];
    UAdsorbateIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem];

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UCationVDWNew[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UHostVDWNew[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem];

    UTailCorrection[CurrentSystem]+=UTailNew;

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UAdsorbateIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem];
      UAdsorbateIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem];
      UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem];

      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]+
                                                 UAdsorbateChargeBondDipoleNew[CurrentSystem]+
                                                 UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+
                                                 UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]+
                                          UAdsorbateChargeBondDipoleNew[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+
                                          UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]+
                                            UHostChargeBondDipoleNew[CurrentSystem]+
                                            UHostBondDipoleBondDipoleNew[CurrentSystem]+
                                            UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]+
                                     UHostChargeBondDipoleNew[CurrentSystem]+
                                     UHostBondDipoleBondDipoleNew[CurrentSystem]+
                                     UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                     UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]+
                                              UCationChargeBondDipoleNew[CurrentSystem]+
                                              UCationBondDipoleBondDipoleNew[CurrentSystem]+
                                              UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]+
                                       UCationChargeBondDipoleNew[CurrentSystem]+
                                       UCationBondDipoleBondDipoleNew[CurrentSystem]+
                                       UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

      NetChargeAdsorbates[CurrentSystem]+=NetChargeAdsorbateDelta;
      NetChargeSystem[CurrentSystem]+=NetChargeAdsorbateDelta;

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }

    InsertAdsorbateMolecule();

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
           UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization+UTailNew;

    UTotal[CurrentSystem]+=DeltaU;

  }
  return 0;
}

void SwapAddMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    SwapAddCationMove();
  else
    SwapAddAdsorbateMove();
}

int SwapRemoveAdsorbateMove(void)
{
  int i;
  REAL RosenbluthOld,PartialFugacity,UTailOld;
  REAL RosenbluthIdealOld;
  REAL DeltaU;

  SwapRemoveAttempts[CurrentSystem][CurrentComponent]+=1.0;

  // return if the component currently has zero molecules
  if(NumberOfAdsorbateMolecules[CurrentSystem]==0) return 0;
  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]<=Components[CurrentComponent].CFMoleculePresent[CurrentSystem]?1:0) return 0;

  CurrentAdsorbateMolecule=SelectRandomMoleculeOfTypeExcludingFractionalMolecule(CurrentComponent);
  CurrentCationMolecule=-1;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be removed here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // calculate the Old Rosenbluth factor
  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_DELETION);
  if (OVERLAP) return 0;

  SwapRemoveAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  UTailOld=TailMolecularEnergyDifferenceRemove();
  RosenbluthOld*=exp(-Beta[CurrentSystem]*UTailOld);

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CalculateEwaldFourierAdsorbate(FALSE,TRUE,CurrentAdsorbateMolecule,0);

    RosenbluthOld*=exp(Beta[CurrentSystem]*(
         UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(FALSE,CurrentAdsorbateMolecule,-1);
    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    RosenbluthOld*=exp(Beta[CurrentSystem]*UDeltaPolarization);
  }

  PartialFugacity=Components[CurrentComponent].FugacityCoefficient[CurrentSystem]*
                  Components[CurrentComponent].PartialPressure[CurrentSystem];

  RosenbluthIdealOld=Components[CurrentComponent].IdealGasRosenbluthWeight[CurrentSystem];

  // acceptence rule
  if(RandomNumber()<((REAL)Components[CurrentComponent].NumberOfMolecules[CurrentSystem]/
                    ((RosenbluthOld/RosenbluthIdealOld)*Beta[CurrentSystem]*PartialFugacity*Volume[CurrentSystem])))
  {
    SwapRemoveAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UAdsorbateBond[CurrentSystem]-=UBondOld[CurrentSystem];
    UAdsorbateUreyBradley[CurrentSystem]-=UUreyBradleyOld[CurrentSystem];
    UAdsorbateBend[CurrentSystem]-=UBendOld[CurrentSystem];
    UAdsorbateBendBend[CurrentSystem]-=UBendBendOld[CurrentSystem];
    UAdsorbateInversionBend[CurrentSystem]-=UInversionBendOld[CurrentSystem];
    UAdsorbateTorsion[CurrentSystem]-=UTorsionOld[CurrentSystem];
    UAdsorbateImproperTorsion[CurrentSystem]-=UImproperTorsionOld[CurrentSystem];
    UAdsorbateBondBond[CurrentSystem]-=UBondBondOld[CurrentSystem];
    UAdsorbateBondBend[CurrentSystem]-=UBondBendOld[CurrentSystem];
    UAdsorbateBondTorsion[CurrentSystem]-=UBondTorsionOld[CurrentSystem];
    UAdsorbateBendTorsion[CurrentSystem]-=UBendTorsionOld[CurrentSystem];
    UAdsorbateIntraVDW[CurrentSystem]-=UIntraVDWOld[CurrentSystem];

    UAdsorbateAdsorbate[CurrentSystem]-=UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]-=UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]-=UCationVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]-=UCationVDWOld[CurrentSystem];
    UHostAdsorbate[CurrentSystem]-=UHostVDWOld[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]-=UHostVDWOld[CurrentSystem];

    UTailCorrection[CurrentSystem]-=UTailOld;

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UAdsorbateIntraChargeCharge[CurrentSystem]-=UIntraChargeChargeOld[CurrentSystem];
      UAdsorbateIntraChargeBondDipole[CurrentSystem]-=UIntraChargeBondDipoleOld[CurrentSystem];
      UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]-=UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]-=UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]-=UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]-=UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]
                                                +UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
                                                +UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                                -UAdsorbateChargeChargeOld[CurrentSystem]
                                                -UAdsorbateChargeBondDipoleOld[CurrentSystem]
                                                -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]
                                         +UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
                                         +UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                         -UAdsorbateChargeChargeOld[CurrentSystem]
                                         -UAdsorbateChargeBondDipoleOld[CurrentSystem]
                                         -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]-=UHostChargeChargeOld[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]-=UHostChargeBondDipoleOld[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]-=UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]
                                           +UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
                                           +UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                           -UHostChargeChargeOld[CurrentSystem]
                                           -UHostChargeBondDipoleOld[CurrentSystem]
                                           -UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]
                                    +UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
                                    +UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                    -UHostChargeChargeOld[CurrentSystem]
                                    -UHostChargeBondDipoleOld[CurrentSystem]
                                    -UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]-=UCationChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=UCationChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]-=UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]
                                             +UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]
                                             +UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                             -UCationChargeChargeOld[CurrentSystem]
                                             -UCationChargeBondDipoleOld[CurrentSystem]
                                             -UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]
                                      +UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]
                                      +UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                      -UCationChargeChargeOld[CurrentSystem]
                                      -UCationChargeBondDipoleOld[CurrentSystem]
                                      -UCationBondDipoleBondDipoleOld[CurrentSystem];

      NetChargeAdsorbates[CurrentSystem]+=NetChargeAdsorbateDelta;
      NetChargeSystem[CurrentSystem]+=NetChargeAdsorbateDelta;

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }

    RemoveAdsorbateMolecule();

    DeltaU=-UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
           -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
           -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
           -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
           -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
           -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
           -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
            UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
            UDeltaPolarization-UTailOld;

    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

void SwapRemoveMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    SwapRemoveCationMove();
  else
    SwapRemoveAdsorbateMove();
}


int SwapAddCationMove(void)
{
  int i;
  REAL RosenbluthNew,PartialFugacity,UTailNew;
  REAL RosenbluthIdealNew;
  REAL DeltaU;

  SwapAddAttempts[CurrentSystem][CurrentComponent]+=1.0;
  CurrentAdsorbateMolecule=NumberOfAdsorbateMolecules[CurrentSystem];
  CurrentCationMolecule=NumberOfCationMolecules[CurrentSystem];

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be added here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNew=GrowMolecule(CBMC_INSERTION);
  if (OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  SwapAddAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  UTailNew=TailMolecularEnergyDifferenceAdd();
  RosenbluthNew*=exp(-Beta[CurrentSystem]*UTailNew);

  if(ChargeMethod!=NONE)
  {
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      CalculateEwaldFourierCation(TRUE,FALSE,NumberOfCationMolecules[CurrentSystem],0);
      RosenbluthNew*=exp(-Beta[CurrentSystem]*(
         UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
    }
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,-1,CurrentCationMolecule);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       UAdsorbatePolarizationNew[CurrentSystem]-UAdsorbatePolarization[CurrentSystem]+
                       (UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       UAdsorbateBackPolarizationNew[CurrentSystem]-UAdsorbateBackPolarization[CurrentSystem]+
                       (UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UCationBackPolarization[CurrentSystem];
    RosenbluthNew*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }

  // get partial pressure for this component
  PartialFugacity=Components[CurrentComponent].FugacityCoefficient[CurrentSystem]*
                  Components[CurrentComponent].PartialPressure[CurrentSystem];

  RosenbluthIdealNew=Components[CurrentComponent].IdealGasRosenbluthWeight[CurrentSystem];

  // acceptence rule
  if(RandomNumber()<((RosenbluthNew/RosenbluthIdealNew)*Beta[CurrentSystem]*PartialFugacity*Volume[CurrentSystem]/
                    ((REAL)1.0+(REAL)Components[CurrentComponent].NumberOfMolecules[CurrentSystem])))
  {
    SwapAddAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UCationBond[CurrentSystem]+=UBondNew[CurrentSystem];
    UCationUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem];
    UCationBend[CurrentSystem]+=UBendNew[CurrentSystem];
    UCationBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem];
    UCationInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem];
    UCationTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem];
    UCationImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem];
    UCationBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem];
    UCationBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem];
    UCationBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem];
    UCationBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem];
    UCationIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem];

    UCationCation[CurrentSystem]+=UCationVDWNew[CurrentSystem];
    UCationCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem];
    UHostCation[CurrentSystem]+=UHostVDWNew[CurrentSystem];
    UHostCationVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];

    UTailCorrection[CurrentSystem]+=UTailNew;

    if(ChargeMethod!=NONE)
    {
      UCationIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem];
      UCationIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem];
      UCationIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem];

      UCationCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]+
                                           UCationChargeBondDipoleNew[CurrentSystem]+
                                           UCationBondDipoleBondDipoleNew[CurrentSystem]+
                                           UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                           UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]+
                                    UCationChargeBondDipoleNew[CurrentSystem]+
                                    UCationBondDipoleBondDipoleNew[CurrentSystem]+
                                    UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                    UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                         UHostChargeChargeNew[CurrentSystem]+
                                         UHostChargeBondDipoleNew[CurrentSystem]+
                                         UHostBondDipoleBondDipoleNew[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                  UHostChargeChargeNew[CurrentSystem]+
                                  UHostChargeBondDipoleNew[CurrentSystem]+
                                  UHostBondDipoleBondDipoleNew[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateChargeChargeNew[CurrentSystem]+
                                              UAdsorbateChargeBondDipoleNew[CurrentSystem]+
                                              UAdsorbateBondDipoleBondDipoleNew[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateChargeChargeNew[CurrentSystem]+
                                       UAdsorbateChargeBondDipoleNew[CurrentSystem]+
                                       UAdsorbateBondDipoleBondDipoleNew[CurrentSystem];

      NetChargeCations[CurrentSystem]+=NetChargeCationDelta;
      NetChargeSystem[CurrentSystem]+=NetChargeCationDelta;

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(0);
    }

    InsertCationMolecule();

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
           UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization+UTailNew;

    UTotal[CurrentSystem]+=DeltaU;
  }
  return 0;
}

int SwapRemoveCationMove(void)
{
  int i;
  REAL RosenbluthOld,PartialFugacity,UTailOld;
  REAL RosenbluthIdealOld;
  REAL DeltaU;

  SwapRemoveAttempts[CurrentSystem][CurrentComponent]+=1.0;

  // return if the component currently has zero molecules
  if(NumberOfCationMolecules[CurrentSystem]==0) return 0;
  if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]<=Components[CurrentComponent].CFMoleculePresent[CurrentSystem]?1:0) return 0;

  // choose a random molecule of this component
  CurrentCationMolecule=SelectRandomMoleculeOfTypeExcludingFractionalMolecule(CurrentComponent);
  CurrentAdsorbateMolecule=-1;

  if(Cations[CurrentSystem][CurrentCationMolecule].Type!=CurrentComponent) printf("ERROR !\n"),exit(0);

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be removed here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // calculate the Old Rosenbluth factor
  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_DELETION);
  if (OVERLAP) return 0;
  SwapRemoveAccepted[CurrentSystem][CurrentComponent][0]+=1.0;

  UTailOld=TailMolecularEnergyDifferenceRemove();
  RosenbluthOld*=exp(-Beta[CurrentSystem]*UTailOld);

  if(ChargeMethod!=NONE)
  {
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      CalculateEwaldFourierCation(FALSE,TRUE,CurrentCationMolecule,0);
      RosenbluthOld*=exp(Beta[CurrentSystem]*(
         UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
    }
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(FALSE,-1,CurrentCationMolecule);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       UAdsorbatePolarizationNew[CurrentSystem]-UAdsorbatePolarization[CurrentSystem]+
                       (UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       UAdsorbateBackPolarizationNew[CurrentSystem]-UAdsorbateBackPolarization[CurrentSystem]+
                       (UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UCationBackPolarization[CurrentSystem];
    RosenbluthOld*=exp(Beta[CurrentSystem]*UDeltaPolarization);
  }


  PartialFugacity=Components[CurrentComponent].FugacityCoefficient[CurrentSystem]*
                  Components[CurrentComponent].PartialPressure[CurrentSystem];

  RosenbluthIdealOld=Components[CurrentComponent].IdealGasRosenbluthWeight[CurrentSystem];

  // acceptence rule
  if(RandomNumber()<((REAL)Components[CurrentComponent].NumberOfMolecules[CurrentSystem]/
                    ((RosenbluthOld/RosenbluthIdealOld)*Beta[CurrentSystem]*PartialFugacity*Volume[CurrentSystem])))
  {
    SwapRemoveAccepted[CurrentSystem][CurrentComponent][1]+=1.0;

    UCationBond[CurrentSystem]-=UBondOld[CurrentSystem];
    UCationUreyBradley[CurrentSystem]-=UUreyBradleyOld[CurrentSystem];
    UCationBend[CurrentSystem]-=UBendOld[CurrentSystem];
    UCationBendBend[CurrentSystem]-=UBendBendOld[CurrentSystem];
    UCationInversionBend[CurrentSystem]-=UInversionBendOld[CurrentSystem];
    UCationTorsion[CurrentSystem]-=UTorsionOld[CurrentSystem];
    UCationImproperTorsion[CurrentSystem]-=UImproperTorsionOld[CurrentSystem];
    UCationBondBond[CurrentSystem]-=UBondBondOld[CurrentSystem];
    UCationBondBend[CurrentSystem]-=UBondBendOld[CurrentSystem];
    UCationBondTorsion[CurrentSystem]-=UBondTorsionOld[CurrentSystem];
    UCationBendTorsion[CurrentSystem]-=UBendTorsionOld[CurrentSystem];
    UCationIntraVDW[CurrentSystem]-=UIntraVDWOld[CurrentSystem];

    UCationCation[CurrentSystem]-=UCationVDWOld[CurrentSystem];
    UCationCationVDW[CurrentSystem]-=UCationVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]-=UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]-=UAdsorbateVDWOld[CurrentSystem];
    UHostCation[CurrentSystem]-=UHostVDWOld[CurrentSystem];
    UHostCationVDW[CurrentSystem]-=UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];

    UTailCorrection[CurrentSystem]-=UTailOld;

    if(ChargeMethod!=NONE)
    {
      UCationIntraChargeCharge[CurrentSystem]-=UIntraChargeChargeOld[CurrentSystem];
      UCationIntraChargeBondDipole[CurrentSystem]-=UIntraChargeBondDipoleOld[CurrentSystem];
      UCationIntraBondDipoleBondDipole[CurrentSystem]-=UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UCationCationChargeChargeReal[CurrentSystem]-=UCationChargeChargeOld[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]-=UCationChargeBondDipoleOld[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]-=UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]
                                          +UCationCationChargeBondDipoleFourierDelta[CurrentSystem]
                                          +UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                          -UCationChargeChargeOld[CurrentSystem]
                                          -UCationChargeBondDipoleOld[CurrentSystem]
                                          -UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]
                                   +UCationCationChargeBondDipoleFourierDelta[CurrentSystem]
                                   +UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                   -UCationChargeChargeOld[CurrentSystem]
                                   -UCationChargeBondDipoleOld[CurrentSystem]
                                   -UCationBondDipoleBondDipoleOld[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]-=UHostChargeChargeOld[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]-=UHostChargeBondDipoleOld[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]-=UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]
                                        +UHostCationChargeBondDipoleFourierDelta[CurrentSystem]
                                        +UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                        -UHostChargeChargeOld[CurrentSystem]
                                        -UHostChargeBondDipoleOld[CurrentSystem]
                                        -UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]
                                 +UHostCationChargeBondDipoleFourierDelta[CurrentSystem]
                                 +UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                 -UHostChargeChargeOld[CurrentSystem]
                                 -UHostChargeBondDipoleOld[CurrentSystem]
                                 -UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]-=UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]-=UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]
                                             +UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]
                                             +UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                             -UAdsorbateChargeChargeOld[CurrentSystem]
                                             -UAdsorbateChargeBondDipoleOld[CurrentSystem]
                                             -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]
                                      +UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]
                                      +UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                      -UAdsorbateChargeChargeOld[CurrentSystem]
                                      -UAdsorbateChargeBondDipoleOld[CurrentSystem]
                                      -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      NetChargeCations[CurrentSystem]+=NetChargeCationDelta;
      NetChargeSystem[CurrentSystem]+=NetChargeCationDelta;

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(0);
    }

    RemoveCationMolecule();

    DeltaU=-UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
           -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
           -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
           -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
           -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
           -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
           -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
           -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
            UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
            UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
            UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
            UDeltaPolarization-UTailOld;

    UTotal[CurrentSystem]+=DeltaU;

  }
  return 0;
}

void PrintSwapAddStatistics(FILE *FilePtr)
{
  int i,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfSwapMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the swap addition move:\n");
    fprintf(FilePtr,"======================================\n");
    for(i=0;i<NumberOfComponents;i++)
      fprintf(FilePtr,"Component [%s] total tried: %lf succesfull growth: %lf (%lf [%%]) accepted: %lf (%lf [%%])\n",
        Components[i].Name,
        (double)SwapAddAttempts[CurrentSystem][i],
        (double)SwapAddAccepted[CurrentSystem][i][0],
        (double)(SwapAddAttempts[CurrentSystem][i]>(REAL)0.0?
          100.0*SwapAddAccepted[CurrentSystem][i][0]/SwapAddAttempts[CurrentSystem][i]:(REAL)0.0),
        (double)SwapAddAccepted[CurrentSystem][i][1],
        (double)(SwapAddAttempts[CurrentSystem][i]>(REAL)0.0?
          100.0*SwapAddAccepted[CurrentSystem][i][1]/SwapAddAttempts[CurrentSystem][i]:(REAL)0.0));
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"Swap addition move was OFF for all components\n\n");
}

void PrintSwapRemoveStatistics(FILE *FilePtr)
{
  int i,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfSwapMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the swap deletion move:\n");
    fprintf(FilePtr,"======================================\n");
    for(i=0;i<NumberOfComponents;i++)
      fprintf(FilePtr,"Component [%s] total tried: %lf succesfull growth: %lf (%lf [%%]) accepted: %lf (%lf [%%])\n",
        Components[i].Name,
        (double)SwapRemoveAttempts[CurrentSystem][i],
        (double)SwapRemoveAccepted[CurrentSystem][i][0],
        (double)(SwapRemoveAttempts[CurrentSystem][i]>(REAL)0.0?
          100.0*SwapRemoveAccepted[CurrentSystem][i][0]/SwapRemoveAttempts[CurrentSystem][i]:(REAL)0.0),
        (double)SwapRemoveAccepted[CurrentSystem][i][1],
        (double)(SwapRemoveAttempts[CurrentSystem][i]>(REAL)0.0?
          100.0*SwapRemoveAccepted[CurrentSystem][i][1]/SwapRemoveAttempts[CurrentSystem][i]:(REAL)0.0));
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"Swap deletion move was OFF for all components\n\n");
}


/*********************************************************************************************************
 * Name       | WidomAdsorbateMove                                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to measure the Rosenbluth weight                                        *
 * Parameters | -                                                                                        *
 * Note       |                                                                                          *
 *********************************************************************************************************/

REAL WidomAdsorbateMove(void)
{
  int i;
  REAL RosenbluthNew,UTailNew,DeltaU;

  CurrentAdsorbateMolecule=NumberOfAdsorbateMolecules[CurrentSystem];
  CurrentCationMolecule=NumberOfCationMolecules[CurrentSystem];

  WidomRosenbluthFactorCount[CurrentSystem][CurrentComponent][Block]+=1.0;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (probed as an integer molecule)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNew=GrowMolecule(CBMC_INSERTION);

  if(OVERLAP) 
    return 0.0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }

  UTailNew=TailMolecularEnergyDifferenceAdd();
  RosenbluthNew*=exp(-Beta[CurrentSystem]*UTailNew);

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CalculateEwaldFourierAdsorbate(TRUE,FALSE,NumberOfAdsorbateMolecules[CurrentSystem],0);
    RosenbluthNew*=exp(-Beta[CurrentSystem]*(
        UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
        UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
        UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    RosenbluthNew*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }


  DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
         UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
         UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
         UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
         UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
         UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
         UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
         UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
         UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
         UDeltaPolarization+UTailNew;

  WidomEnergyDifferenceAverage[CurrentSystem][CurrentComponent][Block]+=(UTotal[CurrentSystem]+DeltaU)*RosenbluthNew;
  WidomRosenbluthFactorAverage[CurrentSystem][CurrentComponent][Block]+=RosenbluthNew;

  WidomEnergyFrameworkAverage[CurrentSystem][CurrentComponent][Block]+=UTotal[CurrentSystem];
  WidomEnergyFrameworkCount[CurrentSystem][CurrentComponent][Block]+=1.0;

  return RosenbluthNew;
}


/*********************************************************************************************************
 * Name       | WidomCationMove                                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to measure the Rosenbluth weight                                        *
 * Parameters | -                                                                                        *
 * Note       |                                                                                          *
 *********************************************************************************************************/

REAL WidomCationMove(void)
{
  int i;
  REAL RosenbluthNew,UTailNew,DeltaU;
    
  CurrentCationMolecule=NumberOfCationMolecules[CurrentSystem];
  CurrentAdsorbateMolecule=NumberOfAdsorbateMolecules[CurrentSystem];
    
  WidomRosenbluthFactorCount[CurrentSystem][CurrentComponent][Block]+=1.0;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (probed as an integer molecule)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }
    
  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNew=GrowMolecule(CBMC_INSERTION);
  if(OVERLAP) return 0;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(BlockedPocket(TrialPosition[CurrentSystem][i]))
      return 0;
  }
    
  UTailNew=TailMolecularEnergyDifferenceAdd();
  RosenbluthNew*=exp(-Beta[CurrentSystem]*UTailNew);

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CalculateEwaldFourierCation(TRUE,FALSE,NumberOfCationMolecules[CurrentSystem],0);
    RosenbluthNew*=exp(-Beta[CurrentSystem]*(
          UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
          UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
          UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,-1,CurrentCationMolecule);
    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    RosenbluthNew*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }

  DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
         UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
         UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
         UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
         UHostVDWNew[CurrentSystem]+UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+
         UHostChargeChargeNew[CurrentSystem]+UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+
         UHostChargeBondDipoleNew[CurrentSystem]+UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+
         UHostBondDipoleBondDipoleNew[CurrentSystem]+UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+
         UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
         UDeltaPolarization+UTailNew;

  WidomEnergyDifferenceAverage[CurrentSystem][CurrentComponent][Block]+=(UTotal[CurrentSystem]+DeltaU)*RosenbluthNew;
  WidomRosenbluthFactorAverage[CurrentSystem][CurrentComponent][Block]+=RosenbluthNew;

  WidomEnergyFrameworkAverage[CurrentSystem][CurrentComponent][Block]+=UTotal[CurrentSystem];
  WidomEnergyFrameworkCount[CurrentSystem][CurrentComponent][Block]+=1.0;

  return RosenbluthNew;
}

void WidomMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    WidomCationMove();
  else
    WidomAdsorbateMove();
}


/*********************************************************************************************************
 * Name       | VolumeMove                                                                               *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to isotropically change the volume of the simulation box.               *
 * Parameters | -                                                                                        *
 * Note       | CF scaling factors are used when set.                                                    *
 *********************************************************************************************************/

int VolumeMove(void)
{
  int i,j,f1;
  int NumberOfMolecules;
  REAL vol,det,scale;
  REAL_MATRIX3x3 StoredBox;
  REAL Pressure,StoredVolume;
  VECTOR com,d;

  CurrentSystem=(int)(RandomNumber()*NumberOfSystems);

  NumberOfMolecules=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    NumberOfMolecules+=Framework[CurrentSystem].TotalNumberOfAtoms;
  if(NumberOfMolecules==0) return 0;

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
  REAL StoredUTotal,StoredUIon,StoredUTailCorrection;

  REAL UHostPolarizationStored,UAdsorbatePolarizationStored,UCationPolarizationStored;
  REAL UHostBackPolarizationStored,UAdsorbateBackPolarizationStored,UCationBackPolarizationStored;

  Pressure=therm_baro_stats.ExternalPressure[CurrentSystem][0];

  VolumeChangeAttempts[CurrentSystem]+=1.0;

  StoredBox=Box[CurrentSystem];
  StoredVolume=Volume[CurrentSystem];

  StoredUTotal=UTotal[CurrentSystem];
  StoredUIon=UIon[CurrentSystem];
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


  // store the positions of the framework
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
       for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
         Framework[CurrentSystem].Atoms[f1][i].ReferencePosition=Framework[CurrentSystem].Atoms[f1][i].Position;

  // store the positions of the adsorbates
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++) 
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition=Adsorbates[CurrentSystem][i].Atoms[j].Position;

  // store the positions of the cations
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++) 
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      Cations[CurrentSystem][i].Atoms[j].ReferencePosition=Cations[CurrentSystem][i].Atoms[j].Position;

  vol=exp(log(Volume[CurrentSystem])+MaximumVolumeChange[CurrentSystem]*(2.0*RandomNumber()-1.0));

  scale=pow(vol/Volume[CurrentSystem],(REAL)1.0/3.0);
  Box[CurrentSystem].ax*=scale;  Box[CurrentSystem].bx*=scale;  Box[CurrentSystem].cx*=scale;
  Box[CurrentSystem].ay*=scale;  Box[CurrentSystem].by*=scale;  Box[CurrentSystem].cy*=scale;
  Box[CurrentSystem].az*=scale;  Box[CurrentSystem].bz*=scale;  Box[CurrentSystem].cz*=scale;
  Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

  if(MIN3(BoxProperties[CurrentSystem].cx,BoxProperties[CurrentSystem].cy,
          BoxProperties[CurrentSystem].cz)<2.0*CutOffVDW)
  {
    for(i=0;i<NumberOfSystems;i++)
    {
      fprintf(OutputFilePtr[i],"\n");
      fprintf(OutputFilePtr[i],"ERROR: (System (%d) Cutoff smaller than half of one of the perpendicular boxlengths !!!\n",CurrentSystem);
      fprintf(OutputFilePtr[i],"       Cutoff: %lf perpendicular boxlengths: %lf %lf %lf\n",(double)CutOffVDW,
              (double)BoxProperties[CurrentSystem].cx,(double)BoxProperties[CurrentSystem].cy,(double)BoxProperties[CurrentSystem].cz);
      fflush(OutputFilePtr[i]);
    }
    exit(0);
  }

  // store the structure-factors for the Ewald-summations
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    SaveCurrentEwaldStructureFactors(0,CurrentSystem);
    SaveCurrentKVectors(0,CurrentSystem);
    SetupKVectors();
  }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Framework[CurrentSystem].Atoms[f1][i].Position.x*=scale;
      Framework[CurrentSystem].Atoms[f1][i].Position.y*=scale;
      Framework[CurrentSystem].Atoms[f1][i].Position.z*=scale;
    }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++) 
  {
    com=GetAdsorbateCenterOfMass(i);

    d.x=com.x*(scale-1.0);
    d.y=com.y*(scale-1.0);
    d.z=com.z*(scale-1.0);

    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Position.x+=d.x;
      Adsorbates[CurrentSystem][i].Atoms[j].Position.y+=d.y;
      Adsorbates[CurrentSystem][i].Atoms[j].Position.z+=d.z;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    com=GetCationCenterOfMass(i);

    d.x=com.x*(scale-1.0);
    d.y=com.y*(scale-1.0);
    d.z=com.z*(scale-1.0);

    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Position.x+=d.x;
      Cations[CurrentSystem][i].Atoms[j].Position.y+=d.y;
      Cations[CurrentSystem][i].Atoms[j].Position.z+=d.z;
    }
  }


  CalculateEnergy();

  if(RandomNumber()<exp(((REAL)NumberOfMolecules+(REAL)1.0)*log(Volume[CurrentSystem]/StoredVolume)
        -(Pressure*(Volume[CurrentSystem]-StoredVolume)+
          (UTotal[CurrentSystem]-StoredUTotal))*Beta[CurrentSystem]))
  {
    VolumeChangeAccepted[CurrentSystem]+=1.0;

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++) 
      UpdateGroupCenterOfMassAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++) 
      UpdateGroupCenterOfMassCation(i);

    EwaldEnergyIon();
  }
  else
  {
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
    UTailCorrection[CurrentSystem]=StoredUTailCorrection;

    UHostPolarization[CurrentSystem]=UHostPolarizationStored;
    UAdsorbatePolarization[CurrentSystem]=UHostPolarizationStored;
    UCationPolarization[CurrentSystem]=UCationPolarizationStored;

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationStored;
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationStored;
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationStored;

    UTotal[CurrentSystem]=StoredUTotal;
    UIon[CurrentSystem]=StoredUIon;
    Box[CurrentSystem]=StoredBox;
    Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
    CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
          Framework[CurrentSystem].Atoms[f1][i].Position=Framework[CurrentSystem].Atoms[f1][i].ReferencePosition;

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        Adsorbates[CurrentSystem][i].Atoms[j].Position=Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition;

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        Cations[CurrentSystem][i].Atoms[j].Position=Cations[CurrentSystem][i].Atoms[j].ReferencePosition;

    CalculateAnisotropicSites();

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      RetrieveStoredEwaldStructureFactors(0,CurrentSystem);
      RetrieveStoredKVectors(0,CurrentSystem);
    }
  }
  return 0;
}

void OptimizeVolumeChangeAcceptence(void)
{
  REAL ratio,vandr;

  if(TotalVolumeChangeAttempts[CurrentSystem]>0.0)
    ratio=TotalVolumeChangeAccepted[CurrentSystem]/TotalVolumeChangeAttempts[CurrentSystem];
  else
    ratio=0.0;

  vandr=ratio/TargetAccRatioVolumeChange;
  if(vandr>1.5) vandr=1.5;
  else if(vandr<0.5) vandr=0.5;
  MaximumVolumeChange[CurrentSystem]*=vandr;
  if(MaximumVolumeChange[CurrentSystem]<0.0005)
     MaximumVolumeChange[CurrentSystem]=0.0005;
  if(MaximumVolumeChange[CurrentSystem]>0.5)
     MaximumVolumeChange[CurrentSystem]=0.5;

  TotalVolumeChangeAttempts[CurrentSystem]+=VolumeChangeAttempts[CurrentSystem];
  TotalVolumeChangeAccepted[CurrentSystem]+=VolumeChangeAccepted[CurrentSystem];
  VolumeChangeAttempts[CurrentSystem]=VolumeChangeAccepted[CurrentSystem]=0.0;
}


void PrintVolumeChangeStatistics(FILE *FilePtr)
{
  if(ProbabilityVolumeChangeMove>0.0)
  {
    fprintf(FilePtr,"Performance of the volume change move:\n");
    fprintf(FilePtr,"======================================\n");
    if(TotalVolumeChangeAttempts[CurrentSystem]>0.0)
    {
       fprintf(FilePtr,"total tried: %lf accepted: %lf (%lf [%%])\n",
          (double)TotalVolumeChangeAttempts[CurrentSystem],
          (double)TotalVolumeChangeAccepted[CurrentSystem],
          (double)(TotalVolumeChangeAttempts[CurrentSystem]>(REAL)0.0?
            100.0*TotalVolumeChangeAccepted[CurrentSystem]/TotalVolumeChangeAttempts[CurrentSystem]:(REAL)0.0));
       fprintf(FilePtr,"\tmaximum volume change %lf\n",(double)MaximumVolumeChange[CurrentSystem]);
    }
    fprintf(FilePtr,"\n");
  }
  else
   fprintf(FilePtr,"Volume move was OFF\n\n");
}


/*********************************************************************************************************
 * Name       | BoxShapeChangeMove                                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to change the shape and size of the simulation box.                     *
 * Parameters | -                                                                                        *
 * Note       | CF scaling factors are used when set.                                                    *
 *********************************************************************************************************/

int BoxShapeChangeMove(void)
{
  int i,j,f1;
  int NumberOfMolecules,ShapeChange;
  REAL det;
  REAL_MATRIX3x3 StoredBox;
  REAL Pressure,StoredVolume;
  VECTOR com,d,pos,s;

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
  REAL StoredUTotal,StoredUIon,StoredUTailCorrection;

  REAL UHostPolarizationStored,UAdsorbatePolarizationStored,UCationPolarizationStored;
  REAL UHostBackPolarizationStored,UAdsorbateBackPolarizationStored,UCationBackPolarizationStored;

  CurrentSystem=(int)(RandomNumber()*NumberOfSystems);

  Pressure=therm_baro_stats.ExternalPressure[CurrentSystem][0];

  StoredBox=Box[CurrentSystem];
  StoredVolume=Volume[CurrentSystem];

  StoredUTotal=UTotal[CurrentSystem];
  StoredUIon=UIon[CurrentSystem];
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

  // store the positions of the framework
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Framework[CurrentSystem].Atoms[f1][i].ReferencePosition=Framework[CurrentSystem].Atoms[f1][i].Position;
      Framework[CurrentSystem].Atoms[f1][i].ReferenceAnisotropicPosition=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;
    }

  // store the positions of the adsorbates
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++) 
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition=Adsorbates[CurrentSystem][i].Atoms[j].Position;
      Adsorbates[CurrentSystem][i].Atoms[j].ReferenceAnisotropicPosition=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
    }

  // store the positions of the cations
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++) 
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].ReferencePosition=Cations[CurrentSystem][i].Atoms[j].Position;
      Cations[CurrentSystem][i].Atoms[j].ReferenceAnisotropicPosition=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
    }

  // store the structure-factors for the Ewald-summations
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    SaveCurrentEwaldStructureFactors(0,CurrentSystem);
    SaveCurrentKVectors(0,CurrentSystem);
    SetupKVectors();
  }

  // choose which index to alter
  ShapeChange=(int)((REAL)6.0*RandomNumber());

  switch(ShapeChange)
  {
    case 0: 
      BoxShapeChangeAttempts[CurrentSystem].ax+=1.0; 
      Box[CurrentSystem].ax+=MaximumBoxShapeChange[CurrentSystem].ax*(2.0*RandomNumber()-1.0);
      break;
    case 1: 
      BoxShapeChangeAttempts[CurrentSystem].bx+=1.0; 
      Box[CurrentSystem].bx+=MaximumBoxShapeChange[CurrentSystem].bx*(2.0*RandomNumber()-1.0);
      break;
    case 2: 
      BoxShapeChangeAttempts[CurrentSystem].by+=1.0; 
      Box[CurrentSystem].by+=MaximumBoxShapeChange[CurrentSystem].by*(2.0*RandomNumber()-1.0);
      break;
    case 3: 
      BoxShapeChangeAttempts[CurrentSystem].cx+=1.0; 
      Box[CurrentSystem].cx+=MaximumBoxShapeChange[CurrentSystem].cx*(2.0*RandomNumber()-1.0);
      break;
    case 4: 
      BoxShapeChangeAttempts[CurrentSystem].cy+=1.0; 
      Box[CurrentSystem].cy+=MaximumBoxShapeChange[CurrentSystem].cy*(2.0*RandomNumber()-1.0);
      break;
    case 5: 
      BoxShapeChangeAttempts[CurrentSystem].cz+=1.0; 
      Box[CurrentSystem].cz+=MaximumBoxShapeChange[CurrentSystem].cz*(2.0*RandomNumber()-1.0);
      break;
  }

  // compute New box properties and volume
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

  if(MIN3(BoxProperties[CurrentSystem].cx,BoxProperties[CurrentSystem].cy,
          BoxProperties[CurrentSystem].cz)<2.0*CutOffVDW)
  {
    for(i=0;i<NumberOfSystems;i++)
    {
      fprintf(OutputFilePtr[i],"\n");
      fprintf(OutputFilePtr[i],"ERROR: (System (%d) Cutoff smaller than half of one of the perpendicular boxlengths !!!\n",CurrentSystem);
      fprintf(OutputFilePtr[i],"       Cutoff: %lf perpendicular boxlengths: %lf %lf %lf\n",(double)CutOffVDW,
              (double)BoxProperties[CurrentSystem].cx,(double)BoxProperties[CurrentSystem].cy,(double)BoxProperties[CurrentSystem].cz);
      fflush(OutputFilePtr[i]);
    }
    exit(0);
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          pos=Framework[CurrentSystem].Atoms[f1][i].Position;

          // convert from xyz to abc
          s.x=InverseBox[CurrentSystem].ax*pos.x+InverseBox[CurrentSystem].bx*pos.y+InverseBox[CurrentSystem].cx*pos.z;
          s.y=InverseBox[CurrentSystem].ay*pos.x+InverseBox[CurrentSystem].by*pos.y+InverseBox[CurrentSystem].cy*pos.z;
          s.z=InverseBox[CurrentSystem].az*pos.x+InverseBox[CurrentSystem].bz*pos.y+InverseBox[CurrentSystem].cz*pos.z;

          // convert from abc to xyz
          pos.x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
          pos.y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
          pos.z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;

          Framework[CurrentSystem].Atoms[f1][i].Position=pos;
        }
      }
    }
  }
  else
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      com=GetIndividualFrameworkCenterOfMass(f1);

      // convert from xyz to abc
      s.x=InverseBox[CurrentSystem].ax*com.x+InverseBox[CurrentSystem].bx*com.y+InverseBox[CurrentSystem].cx*com.z;
      s.y=InverseBox[CurrentSystem].ay*com.x+InverseBox[CurrentSystem].by*com.y+InverseBox[CurrentSystem].cy*com.z;
      s.z=InverseBox[CurrentSystem].az*com.x+InverseBox[CurrentSystem].bz*com.y+InverseBox[CurrentSystem].cz*com.z;

      // convert from abc to xyz
      pos.x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
      pos.y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
      pos.z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;

      d.x=com.x-pos.x;
      d.y=com.y-pos.y;
      d.z=com.z-pos.z;

      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Framework[CurrentSystem].Atoms[f1][i].Position.x+=d.x;
        Framework[CurrentSystem].Atoms[f1][i].Position.y+=d.y;
        Framework[CurrentSystem].Atoms[f1][i].Position.z+=d.z;
      }
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++) 
  {
    com=GetAdsorbateCenterOfMass(i);

    // convert from xyz to abc
    s.x=InverseBox[CurrentSystem].ax*com.x+InverseBox[CurrentSystem].bx*com.y+InverseBox[CurrentSystem].cx*com.z;
    s.y=InverseBox[CurrentSystem].ay*com.x+InverseBox[CurrentSystem].by*com.y+InverseBox[CurrentSystem].cy*com.z;
    s.z=InverseBox[CurrentSystem].az*com.x+InverseBox[CurrentSystem].bz*com.y+InverseBox[CurrentSystem].cz*com.z;

    // convert from abc to xyz
    pos.x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
    pos.y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
    pos.z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;

    d.x=com.x-pos.x;
    d.y=com.y-pos.y;
    d.z=com.z-pos.z;

    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Position.x+=d.x;
      Adsorbates[CurrentSystem][i].Atoms[j].Position.y+=d.y;
      Adsorbates[CurrentSystem][i].Atoms[j].Position.z+=d.z;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    com=GetCationCenterOfMass(i);

    // convert from xyz to abc
    s.x=InverseBox[CurrentSystem].ax*com.x+InverseBox[CurrentSystem].bx*com.y+InverseBox[CurrentSystem].cx*com.z;
    s.y=InverseBox[CurrentSystem].ay*com.x+InverseBox[CurrentSystem].by*com.y+InverseBox[CurrentSystem].cy*com.z;
    s.z=InverseBox[CurrentSystem].az*com.x+InverseBox[CurrentSystem].bz*com.y+InverseBox[CurrentSystem].cz*com.z;

    // convert from abc to xyz
    pos.x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
    pos.y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
    pos.z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;

    d.x=com.x-pos.x;
    d.y=com.y-pos.y;
    d.z=com.z-pos.z;

    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Position.x+=d.x;
      Cations[CurrentSystem][i].Atoms[j].Position.y+=d.y;
      Cations[CurrentSystem][i].Atoms[j].Position.z+=d.z;
    }
  }

  Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);

  CalculateEnergy();

  NumberOfMolecules=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    NumberOfMolecules+=Framework[CurrentSystem].TotalNumberOfAtoms;
  else
    NumberOfMolecules+=Framework[CurrentSystem].NumberOfFrameworks;

  if(RandomNumber()<exp(((REAL)NumberOfMolecules+(REAL)1.0)*log(Volume[CurrentSystem]/StoredVolume)
        -(Pressure*(Volume[CurrentSystem]-StoredVolume)+
          (UTotal[CurrentSystem]-StoredUTotal))*Beta[CurrentSystem]))

  {
    switch(ShapeChange)
    {
      case 0: BoxShapeChangeAccepted[CurrentSystem].ax+=1.0; break;
      case 1: BoxShapeChangeAccepted[CurrentSystem].bx+=1.0; break;
      case 2: BoxShapeChangeAccepted[CurrentSystem].by+=1.0; break;
      case 3: BoxShapeChangeAccepted[CurrentSystem].cx+=1.0; break;
      case 4: BoxShapeChangeAccepted[CurrentSystem].cy+=1.0; break;
      case 5: BoxShapeChangeAccepted[CurrentSystem].cz+=1.0; break;
    }

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++) 
      UpdateGroupCenterOfMassAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++) 
      UpdateGroupCenterOfMassCation(i);

    Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);

    AlphaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bx);
    BetaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].by);
    GammaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bz);

    EwaldEnergyIon();
  }
  else
  {
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
    UTailCorrection[CurrentSystem]=StoredUTailCorrection;

    UHostPolarization[CurrentSystem]=UHostPolarizationStored;
    UAdsorbatePolarization[CurrentSystem]=UHostPolarizationStored;
    UCationPolarization[CurrentSystem]=UCationPolarizationStored;

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationStored;
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationStored;
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationStored;

    UTotal[CurrentSystem]=StoredUTotal;
    UIon[CurrentSystem]=StoredUIon;
    Box[CurrentSystem]=StoredBox;
    Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
    CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Framework[CurrentSystem].Atoms[f1][i].Position=Framework[CurrentSystem].Atoms[f1][i].ReferencePosition;
        Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition=Framework[CurrentSystem].Atoms[f1][i].ReferenceAnisotropicPosition;
      }

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Adsorbates[CurrentSystem][i].Atoms[j].Position=Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition;
        Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition=Adsorbates[CurrentSystem][i].Atoms[j].ReferenceAnisotropicPosition;
      }

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Cations[CurrentSystem][i].Atoms[j].Position=Cations[CurrentSystem][i].Atoms[j].ReferencePosition;
        Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition=Cations[CurrentSystem][i].Atoms[j].ReferenceAnisotropicPosition;
      }

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      RetrieveStoredEwaldStructureFactors(0,CurrentSystem);
      RetrieveStoredKVectors(0,CurrentSystem);
    }
  }
  CalculateAnisotropicSites();
  return 0;
}

void PrintBoxShapeChangeStatistics(FILE *FilePtr)
{
  if(ProbabilityBoxShapeChangeMove>0.0)
  {
    fprintf(FilePtr,"Performance of the box shape change move:\n");
    fprintf(FilePtr,"=========================================\n");
    if(BoxShapeChangeAttempts[CurrentSystem].ax>0.0)
    {
       fprintf(FilePtr,"Box[ax] total tried: %lf accepted: %lf (%lf [%%])\n",
          (double)BoxShapeChangeAttempts[CurrentSystem].ax,
          (double)BoxShapeChangeAccepted[CurrentSystem].ax,
          (double)(BoxShapeChangeAttempts[CurrentSystem].ax>(REAL)0.0?
            100.0*BoxShapeChangeAccepted[CurrentSystem].ax/BoxShapeChangeAttempts[CurrentSystem].ax:(REAL)0.0));
       fprintf(FilePtr,"\tmaximum volume change %lf\n",(double)MaximumBoxShapeChange[CurrentSystem].ax);

       fprintf(FilePtr,"Box[bx] total tried: %lf accepted: %lf (%lf [%%])\n",
          (double)BoxShapeChangeAttempts[CurrentSystem].bx,
          (double)BoxShapeChangeAccepted[CurrentSystem].bx,
          (double)(BoxShapeChangeAttempts[CurrentSystem].bx>(REAL)0.0?
            100.0*BoxShapeChangeAccepted[CurrentSystem].bx/BoxShapeChangeAttempts[CurrentSystem].bx:(REAL)0.0));
       fprintf(FilePtr,"\tmaximum volume change %lf\n",(double)MaximumBoxShapeChange[CurrentSystem].bx);
       fprintf(FilePtr,"Box[by] total tried: %lf accepted: %lf (%lf [%%])\n",
          (double)BoxShapeChangeAttempts[CurrentSystem].by,
          (double)BoxShapeChangeAccepted[CurrentSystem].by,
          (double)(BoxShapeChangeAttempts[CurrentSystem].by>(REAL)0.0?
            100.0*BoxShapeChangeAccepted[CurrentSystem].by/BoxShapeChangeAttempts[CurrentSystem].by:(REAL)0.0));
       fprintf(FilePtr,"\tmaximum volume change %lf\n",(double)MaximumBoxShapeChange[CurrentSystem].by);

       fprintf(FilePtr,"Box[cx] total tried: %lf accepted: %lf (%lf [%%])\n",
          (double)BoxShapeChangeAttempts[CurrentSystem].cx,
          (double)BoxShapeChangeAccepted[CurrentSystem].cx,
          (double)(BoxShapeChangeAttempts[CurrentSystem].cx>(REAL)0.0?
            100.0*BoxShapeChangeAccepted[CurrentSystem].cx/BoxShapeChangeAttempts[CurrentSystem].cx:(REAL)0.0));
       fprintf(FilePtr,"\tmaximum volume change %lf\n",(double)MaximumBoxShapeChange[CurrentSystem].cx);
       fprintf(FilePtr,"Box[cy] total tried: %lf accepted: %lf (%lf [%%])\n",
          (double)BoxShapeChangeAttempts[CurrentSystem].cy,
          (double)BoxShapeChangeAccepted[CurrentSystem].cy,
          (double)(BoxShapeChangeAttempts[CurrentSystem].cy>(REAL)0.0?
            100.0*BoxShapeChangeAccepted[CurrentSystem].cy/BoxShapeChangeAttempts[CurrentSystem].cy:(REAL)0.0));
       fprintf(FilePtr,"\tmaximum volume change %lf\n",(double)MaximumBoxShapeChange[CurrentSystem].cy);
       fprintf(FilePtr,"Box[cz] total tried: %lf accepted: %lf (%lf [%%])\n",
          (double)BoxShapeChangeAttempts[CurrentSystem].cz,
          (double)BoxShapeChangeAccepted[CurrentSystem].cz,
          (double)(BoxShapeChangeAttempts[CurrentSystem].cz>(REAL)0.0?
            100.0*BoxShapeChangeAccepted[CurrentSystem].cz/BoxShapeChangeAttempts[CurrentSystem].cz:(REAL)0.0));
       fprintf(FilePtr,"\tmaximum volume change %lf\n",(double)MaximumBoxShapeChange[CurrentSystem].cz);
    }
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"Box shape change move was OFF\n\n");
}


/*********************************************************************************************************
 * Name       | ParallelTemperingMove                                                                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to swap two systems that differ in temperature.                         *
 * Parameters | -                                                                                        *
 * Note       | CF scaling factors are used when set.                                                    *
 *********************************************************************************************************/

int ParallelTemperingMove(void)
{
  int i,SystemA,SystemB,f1;
  REAL PartialFugacityA,PartialFugacityB;
  int temp_int;
  int *temp_int_pointer;
  REAL temp_real;
  ADSORBATE_MOLECULE *adsorbate_pointer;
  CATION_MOLECULE *cation_pointer;
  ATOM *temp_framework_pointer;

  SystemA=(int)(((REAL)NumberOfSystems-(REAL)1.0)*RandomNumber());
  SystemB=SystemA+1;

  ParallelTemperingAttempts[SystemA][SystemB]+=1.0;
  ParallelTemperingAttempts[SystemB][SystemA]+=1.0;

  PartialFugacityA=Components[CurrentComponent].FugacityCoefficient[SystemA]*
                   Components[CurrentComponent].PartialPressure[SystemA];
  PartialFugacityB=Components[CurrentComponent].FugacityCoefficient[SystemB]*
                   Components[CurrentComponent].PartialPressure[SystemB];

  if(RandomNumber()<exp((Beta[SystemB]-Beta[SystemA])*(UTotal[SystemB]-UTotal[SystemA])))
  {
    ParallelTemperingAccepted[SystemA][SystemB]+=1.0;
    ParallelTemperingAccepted[SystemB][SystemA]+=1.0;

    // swap adsorbates-data for system A and B
    SWAP(NumberOfAdsorbateMolecules[SystemA],NumberOfAdsorbateMolecules[SystemB],temp_int);
    SWAP(NumberOfCationMolecules[SystemA],NumberOfCationMolecules[SystemB],temp_int);
    SWAP(MaxNumberOfAdsorbateMolecules[SystemA],MaxNumberOfAdsorbateMolecules[SystemB],temp_int);
    SWAP(MaxNumberOfCationMolecules[SystemA],MaxNumberOfCationMolecules[SystemB],temp_int);
    SWAP(NumberOfAtomsPerSystem[SystemA],NumberOfAtomsPerSystem[SystemB],temp_int);
    SWAP(NumberOfBondDipolesPerSystem[SystemA],NumberOfBondDipolesPerSystem[SystemB],temp_int);

    SWAP(NumberOfPseudoAtomsCount[SystemA],NumberOfPseudoAtomsCount[SystemB],temp_int_pointer);
    SWAP(NumberOfPseudoAtomsType[SystemA],NumberOfPseudoAtomsType[SystemB],temp_int_pointer);

    // swap the adsorbates and cations of the systems
    SWAP(Adsorbates[SystemA],Adsorbates[SystemB],adsorbate_pointer);
    SWAP(Cations[SystemA],Cations[SystemB],cation_pointer);

    for(i=0;i<NumberOfComponents;i++)
      SWAP(Components[i].NumberOfMolecules[SystemA],Components[i].NumberOfMolecules[SystemB],temp_int);

    for(f1=0;f1<Framework[SystemA].NumberOfFrameworks;f1++)
      if((Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE))
        SWAP(Framework[SystemA].Atoms[f1],Framework[SystemB].Atoms[f1],temp_framework_pointer);

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      SwapEwaldSystem(SystemA,SystemB);

    // swap energies
    SWAP(UTotal[SystemA],UTotal[SystemB],temp_real);
    SWAP(UIon[SystemA],UIon[SystemB],temp_real);
    SWAP(UTailCorrection[SystemA],UTailCorrection[SystemB],temp_real);

    SWAP(UHostPolarization[SystemA],UHostPolarization[SystemB],temp_real);
    SWAP(UAdsorbatePolarization[SystemA],UAdsorbatePolarization[SystemB],temp_real);
    SWAP(UCationPolarization[SystemA],UCationPolarization[SystemB],temp_real);
    SWAP(UHostBackPolarization[SystemA],UHostBackPolarization[SystemB],temp_real);
    SWAP(UAdsorbateBackPolarization[SystemA],UAdsorbateBackPolarization[SystemB],temp_real);
    SWAP(UCationBackPolarization[SystemA],UCationBackPolarization[SystemB],temp_real);

    SWAP(UHostBond[SystemA],UHostBond[SystemB],temp_real);
    SWAP(UHostUreyBradley[SystemA],UHostUreyBradley[SystemB],temp_real);
    SWAP(UHostBend[SystemA],UHostBend[SystemB],temp_real);
    SWAP(UHostInversionBend[SystemA],UHostInversionBend[SystemB],temp_real);
    SWAP(UHostTorsion[SystemA],UHostTorsion[SystemB],temp_real);
    SWAP(UHostImproperTorsion[SystemA],UHostImproperTorsion[SystemB],temp_real);
    SWAP(UHostBondBond[SystemA],UHostBondBond[SystemB],temp_real);
    SWAP(UHostBendBend[SystemA],UHostBendBend[SystemB],temp_real);
    SWAP(UHostBondBend[SystemA],UHostBondBend[SystemB],temp_real);
    SWAP(UHostBondTorsion[SystemA],UHostBondTorsion[SystemB],temp_real);
    SWAP(UHostBendTorsion[SystemA],UHostBendTorsion[SystemB],temp_real);

    SWAP(UAdsorbateBond[SystemA],UAdsorbateBond[SystemB],temp_real);
    SWAP(UAdsorbateUreyBradley[SystemA],UAdsorbateUreyBradley[SystemB],temp_real);
    SWAP(UAdsorbateBend[SystemA],UAdsorbateBend[SystemB],temp_real);
    SWAP(UAdsorbateInversionBend[SystemA],UAdsorbateInversionBend[SystemB],temp_real);
    SWAP(UAdsorbateTorsion[SystemA],UAdsorbateTorsion[SystemB],temp_real);
    SWAP(UAdsorbateImproperTorsion[SystemA],UAdsorbateImproperTorsion[SystemB],temp_real);
    SWAP(UAdsorbateBondBond[SystemA],UAdsorbateBondBond[SystemB],temp_real);
    SWAP(UAdsorbateBendBend[SystemA],UAdsorbateBendBend[SystemB],temp_real);
    SWAP(UAdsorbateBondBend[SystemA],UAdsorbateBondBend[SystemB],temp_real);
    SWAP(UAdsorbateBondTorsion[SystemA],UAdsorbateBondTorsion[SystemB],temp_real);
    SWAP(UAdsorbateBendTorsion[SystemA],UAdsorbateBendTorsion[SystemB],temp_real);
    SWAP(UAdsorbateIntraVDW[SystemA],UAdsorbateIntraVDW[SystemB],temp_real);
    SWAP(UAdsorbateIntraChargeCharge[SystemA],UAdsorbateIntraChargeCharge[SystemB],temp_real);
    SWAP(UAdsorbateIntraChargeBondDipole[SystemA],UAdsorbateIntraChargeBondDipole[SystemB],temp_real);
    SWAP(UAdsorbateIntraBondDipoleBondDipole[SystemA],UAdsorbateIntraBondDipoleBondDipole[SystemB],temp_real);

    SWAP(UCationBond[SystemA],UCationBond[SystemB],temp_real);
    SWAP(UCationUreyBradley[SystemA],UCationUreyBradley[SystemB],temp_real);
    SWAP(UCationBend[SystemA],UCationBend[SystemB],temp_real);
    SWAP(UCationInversionBend[SystemA],UCationInversionBend[SystemB],temp_real);
    SWAP(UCationTorsion[SystemA],UCationTorsion[SystemB],temp_real);
    SWAP(UCationImproperTorsion[SystemA],UCationImproperTorsion[SystemB],temp_real);
    SWAP(UCationBondBond[SystemA],UCationBondBond[SystemB],temp_real);
    SWAP(UCationBendBend[SystemA],UCationBendBend[SystemB],temp_real);
    SWAP(UCationBondBend[SystemA],UCationBondBend[SystemB],temp_real);
    SWAP(UCationBondTorsion[SystemA],UCationBondTorsion[SystemB],temp_real);
    SWAP(UCationBendTorsion[SystemA],UCationBendTorsion[SystemB],temp_real);
    SWAP(UCationIntraVDW[SystemA],UCationIntraVDW[SystemB],temp_real);
    SWAP(UCationIntraChargeCharge[SystemA],UCationIntraChargeCharge[SystemB],temp_real);
    SWAP(UCationIntraChargeBondDipole[SystemA],UCationIntraChargeBondDipole[SystemB],temp_real);
    SWAP(UCationIntraBondDipoleBondDipole[SystemA],UCationIntraBondDipoleBondDipole[SystemB],temp_real);

    SWAP(UHostHost[SystemA],UHostHost[SystemB],temp_real);
    SWAP(UHostHostVDW[SystemA],UHostHostVDW[SystemB],temp_real);
    SWAP(UHostHostChargeChargeReal[SystemA],UHostHostChargeChargeReal[SystemB],temp_real);
    SWAP(UHostHostChargeBondDipoleReal[SystemA],UHostHostChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UHostHostBondDipoleBondDipoleReal[SystemA],UHostHostBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UHostHostChargeChargeFourier[SystemA],UHostHostChargeChargeFourier[SystemB],temp_real);
    SWAP(UHostHostChargeBondDipoleFourier[SystemA],UHostHostChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostHostBondDipoleBondDipoleFourier[SystemA],UHostHostBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostHostCoulomb[SystemA],UHostHostCoulomb[SystemB],temp_real);

    SWAP(UHostAdsorbate[SystemA],UHostAdsorbate[SystemB],temp_real);
    SWAP(UHostAdsorbateVDW[SystemA],UHostAdsorbateVDW[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeChargeReal[SystemA],UHostAdsorbateChargeChargeReal[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeBondDipoleReal[SystemA],UHostAdsorbateChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UHostAdsorbateBondDipoleBondDipoleReal[SystemA],UHostAdsorbateBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeChargeFourier[SystemA],UHostAdsorbateChargeChargeFourier[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeBondDipoleFourier[SystemA],UHostAdsorbateChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostAdsorbateBondDipoleBondDipoleFourier[SystemA],UHostAdsorbateBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostAdsorbateCoulomb[SystemA],UHostAdsorbateCoulomb[SystemB],temp_real);

    SWAP(UHostCation[SystemA],UHostCation[SystemB],temp_real);
    SWAP(UHostCationVDW[SystemA],UHostCationVDW[SystemB],temp_real);
    SWAP(UHostCationChargeChargeReal[SystemA],UHostCationChargeChargeReal[SystemB],temp_real);
    SWAP(UHostCationChargeBondDipoleReal[SystemA],UHostCationChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UHostCationBondDipoleBondDipoleReal[SystemA],UHostCationBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UHostCationChargeChargeFourier[SystemA],UHostCationChargeChargeFourier[SystemB],temp_real);
    SWAP(UHostCationChargeBondDipoleFourier[SystemA],UHostCationChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostCationBondDipoleBondDipoleFourier[SystemA],UHostCationBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostCationCoulomb[SystemA],UHostCationCoulomb[SystemB],temp_real);

    SWAP(UAdsorbateAdsorbate[SystemA],UAdsorbateAdsorbate[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateVDW[SystemA],UAdsorbateAdsorbateVDW[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeChargeReal[SystemA],UAdsorbateAdsorbateChargeChargeReal[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeBondDipoleReal[SystemA],UAdsorbateAdsorbateChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateBondDipoleBondDipoleReal[SystemA],UAdsorbateAdsorbateBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeChargeFourier[SystemA],UAdsorbateAdsorbateChargeChargeFourier[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeBondDipoleFourier[SystemA],UAdsorbateAdsorbateChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateBondDipoleBondDipoleFourier[SystemA],UAdsorbateAdsorbateBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateCoulomb[SystemA],UAdsorbateAdsorbateCoulomb[SystemB],temp_real);

    SWAP(UAdsorbateCation[SystemA],UAdsorbateCation[SystemB],temp_real);
    SWAP(UAdsorbateCationVDW[SystemA],UAdsorbateCationVDW[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeChargeReal[SystemA],UAdsorbateCationChargeChargeReal[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeBondDipoleReal[SystemA],UAdsorbateCationChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateCationBondDipoleBondDipoleReal[SystemA],UAdsorbateCationBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeChargeFourier[SystemA],UAdsorbateCationChargeChargeFourier[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeBondDipoleFourier[SystemA],UAdsorbateCationChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateCationBondDipoleBondDipoleFourier[SystemA],UAdsorbateCationBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateCationCoulomb[SystemA],UAdsorbateCationCoulomb[SystemB],temp_real);

    SWAP(UCationCation[SystemA],UCationCation[SystemB],temp_real);
    SWAP(UCationCationVDW[SystemA],UCationCationVDW[SystemB],temp_real);
    SWAP(UCationCationChargeChargeReal[SystemA],UCationCationChargeChargeReal[SystemB],temp_real);
    SWAP(UCationCationChargeBondDipoleReal[SystemA],UCationCationChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UCationCationBondDipoleBondDipoleReal[SystemA],UCationCationBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UCationCationChargeChargeFourier[SystemA],UCationCationChargeChargeFourier[SystemB],temp_real);
    SWAP(UCationCationChargeBondDipoleFourier[SystemA],UCationCationChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UCationCationBondDipoleBondDipoleFourier[SystemA],UCationCationBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UCationCationCoulomb[SystemA],UCationCationCoulomb[SystemB],temp_real);
  }
  return 0;
}

void PrintParallelTemperingStatistics(FILE *FilePtr)
{
  int i;

  if(ProbabilityParallelTemperingMove>0.0)
  {
    fprintf(FilePtr,"Performance of the parallel tempering move:\n");
    fprintf(FilePtr,"===========================================\n");
    for(i=0;i<NumberOfSystems;i++)
    {
      if(ParallelTemperingAttempts[CurrentSystem][i]>0.0)
      {
         fprintf(FilePtr,"System [%d]<->[%d] total tried: %lf accepted: %lf (%lf [%%])\n",
            CurrentSystem,
            i,
            (double)ParallelTemperingAttempts[CurrentSystem][i],
            (double)ParallelTemperingAccepted[CurrentSystem][i],
            (double)(ParallelTemperingAttempts[CurrentSystem][i]>(REAL)0.0?
              100.0*ParallelTemperingAccepted[CurrentSystem][i]/ParallelTemperingAttempts[CurrentSystem][i]:(REAL)0.0));
    }
    fprintf(FilePtr,"\n");
    }
  }
  else
    fprintf(FilePtr,"Parallel tempering move was OFF\n\n");
}

/*********************************************************************************************************
 * Name       | HyperParallelTemperingMove                                                               *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to swap two systems that differ in temperature and chemical potential.  *
 * Parameters | -                                                                                        *
 * Note       | CF scaling factors are used when set.                                                    *
 * Ref        | "Hyper-parallel tempering Monte Carlo: Application to the Lennard-Jones fluid and the    *
 *            |  restricted primitive model",  G. Yan and J.J. de Pablo, JCP, 111(21): 9509-9516, 1999   *
 *********************************************************************************************************/

int HyperParallelTemperingMove(void)
{
  int i,SystemA,SystemB,f1;
  REAL PartialFugacityA,PartialFugacityB;
  int temp_int;
  int *temp_int_pointer;
  REAL temp_real;
  ADSORBATE_MOLECULE *adsorbate_pointer;
  CATION_MOLECULE *cation_pointer;
  ATOM *temp_framework_pointer;

  SystemA=(int)(((REAL)NumberOfSystems-(REAL)1.0)*RandomNumber());
  SystemB=SystemA+1;

  HyperParallelTemperingAttempts[SystemA][SystemB]+=1.0;
  HyperParallelTemperingAttempts[SystemB][SystemA]+=1.0;

  PartialFugacityA=Components[CurrentComponent].FugacityCoefficient[SystemA]*
                   Components[CurrentComponent].PartialPressure[SystemA];
  PartialFugacityB=Components[CurrentComponent].FugacityCoefficient[SystemB]*
                   Components[CurrentComponent].PartialPressure[SystemB];

  if(RandomNumber()<exp((Beta[SystemB]-Beta[SystemA])*(UTotal[SystemB]-UTotal[SystemA])-
     log((Beta[SystemB]*PartialFugacityB)/(Beta[SystemA]*PartialFugacityA))*
     (Components[CurrentComponent].NumberOfMolecules[SystemB]-Components[CurrentComponent].NumberOfMolecules[SystemA])))

  {
    HyperParallelTemperingAccepted[SystemA][SystemB]+=1.0;
    HyperParallelTemperingAccepted[SystemB][SystemA]+=1.0;

    // swap adsorbates-data for system A and B
    SWAP(NumberOfAdsorbateMolecules[SystemA],NumberOfAdsorbateMolecules[SystemB],temp_int);
    SWAP(NumberOfCationMolecules[SystemA],NumberOfCationMolecules[SystemB],temp_int);
    SWAP(MaxNumberOfAdsorbateMolecules[SystemA],MaxNumberOfAdsorbateMolecules[SystemB],temp_int);
    SWAP(MaxNumberOfCationMolecules[SystemA],MaxNumberOfCationMolecules[SystemB],temp_int);
    SWAP(NumberOfAtomsPerSystem[SystemA],NumberOfAtomsPerSystem[SystemB],temp_int);
    SWAP(NumberOfBondDipolesPerSystem[SystemA],NumberOfBondDipolesPerSystem[SystemB],temp_int);

    SWAP(NumberOfPseudoAtomsCount[SystemA],NumberOfPseudoAtomsCount[SystemB],temp_int_pointer);
    SWAP(NumberOfPseudoAtomsType[SystemA],NumberOfPseudoAtomsType[SystemB],temp_int_pointer);

    // swap the adsorbates and cations of the systems
    SWAP(Adsorbates[SystemA],Adsorbates[SystemB],adsorbate_pointer);
    SWAP(Cations[SystemA],Cations[SystemB],cation_pointer);

    for(i=0;i<NumberOfComponents;i++)
      SWAP(Components[i].NumberOfMolecules[SystemA],Components[i].NumberOfMolecules[SystemB],temp_int);

    for(f1=0;f1<Framework[SystemA].NumberOfFrameworks;f1++)
      if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
        SWAP(Framework[SystemA].Atoms[f1],Framework[SystemB].Atoms[f1],temp_framework_pointer);

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      SwapEwaldSystem(SystemA,SystemB);

    // swap energies
    SWAP(UTotal[SystemA],UTotal[SystemB],temp_real);
    SWAP(UIon[SystemA],UIon[SystemB],temp_real);
    SWAP(UTailCorrection[SystemA],UTailCorrection[SystemB],temp_real);

    SWAP(UHostPolarization[SystemA],UHostPolarization[SystemB],temp_real);
    SWAP(UAdsorbatePolarization[SystemA],UAdsorbatePolarization[SystemB],temp_real);
    SWAP(UCationPolarization[SystemA],UCationPolarization[SystemB],temp_real);
    SWAP(UHostBackPolarization[SystemA],UHostBackPolarization[SystemB],temp_real);
    SWAP(UAdsorbateBackPolarization[SystemA],UAdsorbateBackPolarization[SystemB],temp_real);
    SWAP(UCationBackPolarization[SystemA],UCationBackPolarization[SystemB],temp_real);

    SWAP(UHostBond[SystemA],UHostBond[SystemB],temp_real);
    SWAP(UHostUreyBradley[SystemA],UHostUreyBradley[SystemB],temp_real);
    SWAP(UHostBend[SystemA],UHostBend[SystemB],temp_real);
    SWAP(UHostInversionBend[SystemA],UHostInversionBend[SystemB],temp_real);
    SWAP(UHostTorsion[SystemA],UHostTorsion[SystemB],temp_real);
    SWAP(UHostImproperTorsion[SystemA],UHostImproperTorsion[SystemB],temp_real);
    SWAP(UHostBondBond[SystemA],UHostBondBond[SystemB],temp_real);
    SWAP(UHostBendBend[SystemA],UHostBendBend[SystemB],temp_real);
    SWAP(UHostBondBend[SystemA],UHostBondBend[SystemB],temp_real);
    SWAP(UHostBondTorsion[SystemA],UHostBondTorsion[SystemB],temp_real);
    SWAP(UHostBendTorsion[SystemA],UHostBendTorsion[SystemB],temp_real);

    SWAP(UAdsorbateBond[SystemA],UAdsorbateBond[SystemB],temp_real);
    SWAP(UAdsorbateUreyBradley[SystemA],UAdsorbateUreyBradley[SystemB],temp_real);
    SWAP(UAdsorbateBend[SystemA],UAdsorbateBend[SystemB],temp_real);
    SWAP(UAdsorbateInversionBend[SystemA],UAdsorbateInversionBend[SystemB],temp_real);
    SWAP(UAdsorbateTorsion[SystemA],UAdsorbateTorsion[SystemB],temp_real);
    SWAP(UAdsorbateImproperTorsion[SystemA],UAdsorbateImproperTorsion[SystemB],temp_real);
    SWAP(UAdsorbateBondBond[SystemA],UAdsorbateBondBond[SystemB],temp_real);
    SWAP(UAdsorbateBendBend[SystemA],UAdsorbateBendBend[SystemB],temp_real);
    SWAP(UAdsorbateBondBend[SystemA],UAdsorbateBondBend[SystemB],temp_real);
    SWAP(UAdsorbateBondTorsion[SystemA],UAdsorbateBondTorsion[SystemB],temp_real);
    SWAP(UAdsorbateBendTorsion[SystemA],UAdsorbateBendTorsion[SystemB],temp_real);
    SWAP(UAdsorbateIntraVDW[SystemA],UAdsorbateIntraVDW[SystemB],temp_real);
    SWAP(UAdsorbateIntraChargeCharge[SystemA],UAdsorbateIntraChargeCharge[SystemB],temp_real);
    SWAP(UAdsorbateIntraChargeBondDipole[SystemA],UAdsorbateIntraChargeBondDipole[SystemB],temp_real);
    SWAP(UAdsorbateIntraBondDipoleBondDipole[SystemA],UAdsorbateIntraBondDipoleBondDipole[SystemB],temp_real);

    SWAP(UCationBond[SystemA],UCationBond[SystemB],temp_real);
    SWAP(UCationUreyBradley[SystemA],UCationUreyBradley[SystemB],temp_real);
    SWAP(UCationBend[SystemA],UCationBend[SystemB],temp_real);
    SWAP(UCationInversionBend[SystemA],UCationInversionBend[SystemB],temp_real);
    SWAP(UCationTorsion[SystemA],UCationTorsion[SystemB],temp_real);
    SWAP(UCationImproperTorsion[SystemA],UCationImproperTorsion[SystemB],temp_real);
    SWAP(UCationBondBond[SystemA],UCationBondBond[SystemB],temp_real);
    SWAP(UCationBendBend[SystemA],UCationBendBend[SystemB],temp_real);
    SWAP(UCationBondBend[SystemA],UCationBondBend[SystemB],temp_real);
    SWAP(UCationBondTorsion[SystemA],UCationBondTorsion[SystemB],temp_real);
    SWAP(UCationBendTorsion[SystemA],UCationBendTorsion[SystemB],temp_real);
    SWAP(UCationIntraVDW[SystemA],UCationIntraVDW[SystemB],temp_real);
    SWAP(UCationIntraChargeCharge[SystemA],UCationIntraChargeCharge[SystemB],temp_real);
    SWAP(UCationIntraChargeBondDipole[SystemA],UCationIntraChargeBondDipole[SystemB],temp_real);
    SWAP(UCationIntraBondDipoleBondDipole[SystemA],UCationIntraBondDipoleBondDipole[SystemB],temp_real);

    SWAP(UHostHost[SystemA],UHostHost[SystemB],temp_real);
    SWAP(UHostHostVDW[SystemA],UHostHostVDW[SystemB],temp_real);
    SWAP(UHostHostChargeChargeReal[SystemA],UHostHostChargeChargeReal[SystemB],temp_real);
    SWAP(UHostHostChargeBondDipoleReal[SystemA],UHostHostChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UHostHostBondDipoleBondDipoleReal[SystemA],UHostHostBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UHostHostChargeChargeFourier[SystemA],UHostHostChargeChargeFourier[SystemB],temp_real);
    SWAP(UHostHostChargeBondDipoleFourier[SystemA],UHostHostChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostHostBondDipoleBondDipoleFourier[SystemA],UHostHostBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostHostCoulomb[SystemA],UHostHostCoulomb[SystemB],temp_real);

    SWAP(UHostAdsorbate[SystemA],UHostAdsorbate[SystemB],temp_real);
    SWAP(UHostAdsorbateVDW[SystemA],UHostAdsorbateVDW[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeChargeReal[SystemA],UHostAdsorbateChargeChargeReal[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeBondDipoleReal[SystemA],UHostAdsorbateChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UHostAdsorbateBondDipoleBondDipoleReal[SystemA],UHostAdsorbateBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeChargeFourier[SystemA],UHostAdsorbateChargeChargeFourier[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeBondDipoleFourier[SystemA],UHostAdsorbateChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostAdsorbateBondDipoleBondDipoleFourier[SystemA],UHostAdsorbateBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostAdsorbateCoulomb[SystemA],UHostAdsorbateCoulomb[SystemB],temp_real);

    SWAP(UHostCation[SystemA],UHostCation[SystemB],temp_real);
    SWAP(UHostCationVDW[SystemA],UHostCationVDW[SystemB],temp_real);
    SWAP(UHostCationChargeChargeReal[SystemA],UHostCationChargeChargeReal[SystemB],temp_real);
    SWAP(UHostCationChargeBondDipoleReal[SystemA],UHostCationChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UHostCationBondDipoleBondDipoleReal[SystemA],UHostCationBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UHostCationChargeChargeFourier[SystemA],UHostCationChargeChargeFourier[SystemB],temp_real);
    SWAP(UHostCationChargeBondDipoleFourier[SystemA],UHostCationChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostCationBondDipoleBondDipoleFourier[SystemA],UHostCationBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostCationCoulomb[SystemA],UHostCationCoulomb[SystemB],temp_real);

    SWAP(UAdsorbateAdsorbate[SystemA],UAdsorbateAdsorbate[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateVDW[SystemA],UAdsorbateAdsorbateVDW[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeChargeReal[SystemA],UAdsorbateAdsorbateChargeChargeReal[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeBondDipoleReal[SystemA],UAdsorbateAdsorbateChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateBondDipoleBondDipoleReal[SystemA],UAdsorbateAdsorbateBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeChargeFourier[SystemA],UAdsorbateAdsorbateChargeChargeFourier[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeBondDipoleFourier[SystemA],UAdsorbateAdsorbateChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateBondDipoleBondDipoleFourier[SystemA],UAdsorbateAdsorbateBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateCoulomb[SystemA],UAdsorbateAdsorbateCoulomb[SystemB],temp_real);

    SWAP(UAdsorbateCation[SystemA],UAdsorbateCation[SystemB],temp_real);
    SWAP(UAdsorbateCationVDW[SystemA],UAdsorbateCationVDW[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeChargeReal[SystemA],UAdsorbateCationChargeChargeReal[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeBondDipoleReal[SystemA],UAdsorbateCationChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateCationBondDipoleBondDipoleReal[SystemA],UAdsorbateCationBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeChargeFourier[SystemA],UAdsorbateCationChargeChargeFourier[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeBondDipoleFourier[SystemA],UAdsorbateCationChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateCationBondDipoleBondDipoleFourier[SystemA],UAdsorbateCationBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateCationCoulomb[SystemA],UAdsorbateCationCoulomb[SystemB],temp_real);

    SWAP(UCationCation[SystemA],UCationCation[SystemB],temp_real);
    SWAP(UCationCationVDW[SystemA],UCationCationVDW[SystemB],temp_real);
    SWAP(UCationCationChargeChargeReal[SystemA],UCationCationChargeChargeReal[SystemB],temp_real);
    SWAP(UCationCationChargeBondDipoleReal[SystemA],UCationCationChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UCationCationBondDipoleBondDipoleReal[SystemA],UCationCationBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UCationCationChargeChargeFourier[SystemA],UCationCationChargeChargeFourier[SystemB],temp_real);
    SWAP(UCationCationChargeBondDipoleFourier[SystemA],UCationCationChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UCationCationBondDipoleBondDipoleFourier[SystemA],UCationCationBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UCationCationCoulomb[SystemA],UCationCationCoulomb[SystemB],temp_real);
  }
  return 0;
}

void PrintHyperParallelTemperingStatistics(FILE *FilePtr)
{
  int i;

  if(ProbabilityHyperParallelTemperingMove>0.0)
  {
    fprintf(FilePtr,"Performance of the hyper parallel tempering move:\n");
    fprintf(FilePtr,"=================================================\n");
    for(i=0;i<NumberOfSystems;i++)
    {
      if(HyperParallelTemperingAttempts[CurrentSystem][i]>0.0)
      {
         fprintf(FilePtr,"System [%d]<->[%d] total tried: %lf accepted: %lf (%lf [%%])\n",
            CurrentSystem,
            i,
            (double)HyperParallelTemperingAttempts[CurrentSystem][i],
            (double)HyperParallelTemperingAccepted[CurrentSystem][i],
            (double)(HyperParallelTemperingAttempts[CurrentSystem][i]>(REAL)0.0?
              100.0*HyperParallelTemperingAccepted[CurrentSystem][i]/HyperParallelTemperingAttempts[CurrentSystem][i]:(REAL)0.0));
    }
    fprintf(FilePtr,"\n");
    }
  }
  else
    fprintf(FilePtr,"Hyper parallel tempering move was OFF\n\n");
}


/*********************************************************************************************************
 * Name       | ParallelMolFractionMove                                                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to swap two systems that differ in mol-fractions                        *
 * Parameters | -                                                                                        *
 * Note       |                                                                                          *
 *********************************************************************************************************/

int ParallelMolFractionMove(void)
{
  int i,SystemA,SystemB,f1;
  int temp_int;
  int *temp_int_pointer;
  REAL temp_real;
  CATION_MOLECULE *cation_pointer;
  ADSORBATE_MOLECULE *adsorbate_pointer;
  ATOM *temp_framework_pointer;
  REAL QSA,QSB;
  int NRA,NRB,NSA,NSB;

  SystemA=(int)(((REAL)NumberOfSystems-(REAL)1.0)*RandomNumber());
  SystemB=SystemA+1;

  ParallelMolFractionAttempts[SystemA][SystemB]+=1.0;
  ParallelMolFractionAttempts[SystemB][SystemA]+=1.0;

  QSA=Components[ParallelMolFractionComponentA].MolFraction[SystemA];
  QSB=Components[ParallelMolFractionComponentA].MolFraction[SystemB];

  NSA=Components[ParallelMolFractionComponentA].NumberOfMolecules[SystemA];
  NSB=Components[ParallelMolFractionComponentA].NumberOfMolecules[SystemB];

  NRA=Components[ParallelMolFractionComponentB].NumberOfMolecules[SystemA];
  NRB=Components[ParallelMolFractionComponentB].NumberOfMolecules[SystemB];

  if(RandomNumber()<pow(QSA,NSB-NSA)*pow(1.0-QSA,NRB-NRA)*pow(QSB,NSA-NSB)*pow(1.0-QSB,NRA-NRB))
  {
    ParallelMolFractionAccepted[SystemA][SystemB]+=1.0;
    ParallelMolFractionAccepted[SystemB][SystemA]+=1.0;

    // swap adsorbates-data for system A and B
    SWAP(NumberOfAdsorbateMolecules[SystemA],NumberOfAdsorbateMolecules[SystemB],temp_int);
    SWAP(NumberOfCationMolecules[SystemA],NumberOfCationMolecules[SystemB],temp_int);
    SWAP(MaxNumberOfAdsorbateMolecules[SystemA],MaxNumberOfAdsorbateMolecules[SystemB],temp_int);
    SWAP(MaxNumberOfCationMolecules[SystemA],MaxNumberOfCationMolecules[SystemB],temp_int);
    SWAP(NumberOfAtomsPerSystem[SystemA],NumberOfAtomsPerSystem[SystemB],temp_int);
    SWAP(NumberOfBondDipolesPerSystem[SystemA],NumberOfBondDipolesPerSystem[SystemB],temp_int);

    SWAP(NumberOfPseudoAtomsCount[SystemA],NumberOfPseudoAtomsCount[SystemB],temp_int_pointer);
    SWAP(NumberOfPseudoAtomsType[SystemA],NumberOfPseudoAtomsType[SystemB],temp_int_pointer);

    // swap the adsorbates and cations of the systems
    SWAP(Adsorbates[SystemA],Adsorbates[SystemB],adsorbate_pointer);
    SWAP(Cations[SystemA],Cations[SystemB],cation_pointer);

    for(i=0;i<NumberOfComponents;i++)
      SWAP(Components[i].NumberOfMolecules[SystemA],Components[i].NumberOfMolecules[SystemB],temp_int);

    for(f1=0;f1<Framework[SystemA].NumberOfFrameworks;f1++)
      if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
        SWAP(Framework[SystemA].Atoms[f1],Framework[SystemB].Atoms[f1],temp_framework_pointer);

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      SwapEwaldSystem(SystemA,SystemB);

    // swap energies
    SWAP(UTotal[SystemA],UTotal[SystemB],temp_real);
    SWAP(UIon[SystemA],UIon[SystemB],temp_real);
    SWAP(UTailCorrection[SystemA],UTailCorrection[SystemB],temp_real);

    SWAP(UHostPolarization[SystemA],UHostPolarization[SystemB],temp_real);
    SWAP(UAdsorbatePolarization[SystemA],UAdsorbatePolarization[SystemB],temp_real);
    SWAP(UCationPolarization[SystemA],UCationPolarization[SystemB],temp_real);
    SWAP(UHostBackPolarization[SystemA],UHostBackPolarization[SystemB],temp_real);
    SWAP(UAdsorbateBackPolarization[SystemA],UAdsorbateBackPolarization[SystemB],temp_real);
    SWAP(UCationBackPolarization[SystemA],UCationBackPolarization[SystemB],temp_real);

    SWAP(UHostBond[SystemA],UHostBond[SystemB],temp_real);
    SWAP(UHostUreyBradley[SystemA],UHostUreyBradley[SystemB],temp_real);
    SWAP(UHostBend[SystemA],UHostBend[SystemB],temp_real);
    SWAP(UHostInversionBend[SystemA],UHostInversionBend[SystemB],temp_real);
    SWAP(UHostTorsion[SystemA],UHostTorsion[SystemB],temp_real);
    SWAP(UHostImproperTorsion[SystemA],UHostImproperTorsion[SystemB],temp_real);
    SWAP(UHostBondBond[SystemA],UHostBondBond[SystemB],temp_real);
    SWAP(UHostBendBend[SystemA],UHostBendBend[SystemB],temp_real);
    SWAP(UHostBondBend[SystemA],UHostBondBend[SystemB],temp_real);
    SWAP(UHostBondTorsion[SystemA],UHostBondTorsion[SystemB],temp_real);
    SWAP(UHostBendTorsion[SystemA],UHostBendTorsion[SystemB],temp_real);

    SWAP(UAdsorbateBond[SystemA],UAdsorbateBond[SystemB],temp_real);
    SWAP(UAdsorbateUreyBradley[SystemA],UAdsorbateUreyBradley[SystemB],temp_real);
    SWAP(UAdsorbateBend[SystemA],UAdsorbateBend[SystemB],temp_real);
    SWAP(UAdsorbateInversionBend[SystemA],UAdsorbateInversionBend[SystemB],temp_real);
    SWAP(UAdsorbateTorsion[SystemA],UAdsorbateTorsion[SystemB],temp_real);
    SWAP(UAdsorbateImproperTorsion[SystemA],UAdsorbateImproperTorsion[SystemB],temp_real);
    SWAP(UAdsorbateBondBond[SystemA],UAdsorbateBondBond[SystemB],temp_real);
    SWAP(UAdsorbateBendBend[SystemA],UAdsorbateBendBend[SystemB],temp_real);
    SWAP(UAdsorbateBondBend[SystemA],UAdsorbateBondBend[SystemB],temp_real);
    SWAP(UAdsorbateBondTorsion[SystemA],UAdsorbateBondTorsion[SystemB],temp_real);
    SWAP(UAdsorbateBendTorsion[SystemA],UAdsorbateBendTorsion[SystemB],temp_real);
    SWAP(UAdsorbateIntraVDW[SystemA],UAdsorbateIntraVDW[SystemB],temp_real);
    SWAP(UAdsorbateIntraChargeCharge[SystemA],UAdsorbateIntraChargeCharge[SystemB],temp_real);
    SWAP(UAdsorbateIntraChargeBondDipole[SystemA],UAdsorbateIntraChargeBondDipole[SystemB],temp_real);
    SWAP(UAdsorbateIntraBondDipoleBondDipole[SystemA],UAdsorbateIntraBondDipoleBondDipole[SystemB],temp_real);

    SWAP(UCationBond[SystemA],UCationBond[SystemB],temp_real);
    SWAP(UCationUreyBradley[SystemA],UCationUreyBradley[SystemB],temp_real);
    SWAP(UCationBend[SystemA],UCationBend[SystemB],temp_real);
    SWAP(UCationInversionBend[SystemA],UCationInversionBend[SystemB],temp_real);
    SWAP(UCationTorsion[SystemA],UCationTorsion[SystemB],temp_real);
    SWAP(UCationImproperTorsion[SystemA],UCationImproperTorsion[SystemB],temp_real);
    SWAP(UCationBondBond[SystemA],UCationBondBond[SystemB],temp_real);
    SWAP(UCationBendBend[SystemA],UCationBendBend[SystemB],temp_real);
    SWAP(UCationBondBend[SystemA],UCationBondBend[SystemB],temp_real);
    SWAP(UCationBondTorsion[SystemA],UCationBondTorsion[SystemB],temp_real);
    SWAP(UCationBendTorsion[SystemA],UCationBendTorsion[SystemB],temp_real);
    SWAP(UCationIntraVDW[SystemA],UCationIntraVDW[SystemB],temp_real);
    SWAP(UCationIntraChargeCharge[SystemA],UCationIntraChargeCharge[SystemB],temp_real);
    SWAP(UCationIntraChargeBondDipole[SystemA],UCationIntraChargeBondDipole[SystemB],temp_real);
    SWAP(UCationIntraBondDipoleBondDipole[SystemA],UCationIntraBondDipoleBondDipole[SystemB],temp_real);

    SWAP(UHostHost[SystemA],UHostHost[SystemB],temp_real);
    SWAP(UHostHostVDW[SystemA],UHostHostVDW[SystemB],temp_real);
    SWAP(UHostHostChargeChargeReal[SystemA],UHostHostChargeChargeReal[SystemB],temp_real);
    SWAP(UHostHostChargeBondDipoleReal[SystemA],UHostHostChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UHostHostBondDipoleBondDipoleReal[SystemA],UHostHostBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UHostHostChargeChargeFourier[SystemA],UHostHostChargeChargeFourier[SystemB],temp_real);
    SWAP(UHostHostChargeBondDipoleFourier[SystemA],UHostHostChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostHostBondDipoleBondDipoleFourier[SystemA],UHostHostBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostHostCoulomb[SystemA],UHostHostCoulomb[SystemB],temp_real);

    SWAP(UHostAdsorbate[SystemA],UHostAdsorbate[SystemB],temp_real);
    SWAP(UHostAdsorbateVDW[SystemA],UHostAdsorbateVDW[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeChargeReal[SystemA],UHostAdsorbateChargeChargeReal[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeBondDipoleReal[SystemA],UHostAdsorbateChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UHostAdsorbateBondDipoleBondDipoleReal[SystemA],UHostAdsorbateBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeChargeFourier[SystemA],UHostAdsorbateChargeChargeFourier[SystemB],temp_real);
    SWAP(UHostAdsorbateChargeBondDipoleFourier[SystemA],UHostAdsorbateChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostAdsorbateBondDipoleBondDipoleFourier[SystemA],UHostAdsorbateBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostAdsorbateCoulomb[SystemA],UHostAdsorbateCoulomb[SystemB],temp_real);

    SWAP(UHostCation[SystemA],UHostCation[SystemB],temp_real);
    SWAP(UHostCationVDW[SystemA],UHostCationVDW[SystemB],temp_real);
    SWAP(UHostCationChargeChargeReal[SystemA],UHostCationChargeChargeReal[SystemB],temp_real);
    SWAP(UHostCationChargeBondDipoleReal[SystemA],UHostCationChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UHostCationBondDipoleBondDipoleReal[SystemA],UHostCationBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UHostCationChargeChargeFourier[SystemA],UHostCationChargeChargeFourier[SystemB],temp_real);
    SWAP(UHostCationChargeBondDipoleFourier[SystemA],UHostCationChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostCationBondDipoleBondDipoleFourier[SystemA],UHostCationBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UHostCationCoulomb[SystemA],UHostCationCoulomb[SystemB],temp_real);

    SWAP(UAdsorbateAdsorbate[SystemA],UAdsorbateAdsorbate[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateVDW[SystemA],UAdsorbateAdsorbateVDW[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeChargeReal[SystemA],UAdsorbateAdsorbateChargeChargeReal[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeBondDipoleReal[SystemA],UAdsorbateAdsorbateChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateBondDipoleBondDipoleReal[SystemA],UAdsorbateAdsorbateBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeChargeFourier[SystemA],UAdsorbateAdsorbateChargeChargeFourier[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateChargeBondDipoleFourier[SystemA],UAdsorbateAdsorbateChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateBondDipoleBondDipoleFourier[SystemA],UAdsorbateAdsorbateBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateAdsorbateCoulomb[SystemA],UAdsorbateAdsorbateCoulomb[SystemB],temp_real);

    SWAP(UAdsorbateCation[SystemA],UAdsorbateCation[SystemB],temp_real);
    SWAP(UAdsorbateCationVDW[SystemA],UAdsorbateCationVDW[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeChargeReal[SystemA],UAdsorbateCationChargeChargeReal[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeBondDipoleReal[SystemA],UAdsorbateCationChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateCationBondDipoleBondDipoleReal[SystemA],UAdsorbateCationBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeChargeFourier[SystemA],UAdsorbateCationChargeChargeFourier[SystemB],temp_real);
    SWAP(UAdsorbateCationChargeBondDipoleFourier[SystemA],UAdsorbateCationChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateCationBondDipoleBondDipoleFourier[SystemA],UAdsorbateCationBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UAdsorbateCationCoulomb[SystemA],UAdsorbateCationCoulomb[SystemB],temp_real);

    SWAP(UCationCation[SystemA],UCationCation[SystemB],temp_real);
    SWAP(UCationCationVDW[SystemA],UCationCationVDW[SystemB],temp_real);
    SWAP(UCationCationChargeChargeReal[SystemA],UCationCationChargeChargeReal[SystemB],temp_real);
    SWAP(UCationCationChargeBondDipoleReal[SystemA],UCationCationChargeBondDipoleReal[SystemB],temp_real);
    SWAP(UCationCationBondDipoleBondDipoleReal[SystemA],UCationCationBondDipoleBondDipoleReal[SystemB],temp_real);
    SWAP(UCationCationChargeChargeFourier[SystemA],UCationCationChargeChargeFourier[SystemB],temp_real);
    SWAP(UCationCationChargeBondDipoleFourier[SystemA],UCationCationChargeBondDipoleFourier[SystemB],temp_real);
    SWAP(UCationCationBondDipoleBondDipoleFourier[SystemA],UCationCationBondDipoleBondDipoleFourier[SystemB],temp_real);
    SWAP(UCationCationCoulomb[SystemA],UCationCationCoulomb[SystemB],temp_real);
  }
  return 0;
}

void PrintParallelMolFractionStatistics(FILE *FilePtr)
{
  int i;

  if(ProbabilityParallelMolFractionMove>0.0)
  {
    fprintf(FilePtr,"Performance of the parallel mol-fraction move:\n");
    fprintf(FilePtr,"==============================================\n");
    for(i=0;i<NumberOfSystems;i++)
    {
      if(ParallelMolFractionAttempts[CurrentSystem][i]>0.0)
      {
         fprintf(FilePtr,"System [%d]<->[%d] total tried: %lf accepted: %lf (%lf [%%])\n",
            CurrentSystem,
            i,
            (double)ParallelMolFractionAttempts[CurrentSystem][i],
            (double)ParallelMolFractionAccepted[CurrentSystem][i],
            (double)(ParallelMolFractionAttempts[CurrentSystem][i]>(REAL)0.0?
              100.0*ParallelMolFractionAccepted[CurrentSystem][i]/ParallelMolFractionAttempts[CurrentSystem][i]:(REAL)0.0));
      }
    }
  }
  else
    fprintf(FilePtr,"Parallel mol-fraction move was OFF\n");
  fprintf(FilePtr,"\n");
}


/*********************************************************************************************************
 * Name       | ChiralInversionMove                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to change all chiral R into S, and S into R                             *
 * Parameters | -                                                                                        *
 * Note       | This move is only allowed for systems with an inversion center at (1/2,1/2,1/2), e.g.    *
 *            | MFI and LTA. It modifies all adsorbate and cations posiitons as s'=1-s.                  *
 *********************************************************************************************************/

int ChiralInversionMove(void)
{
  int i,j,m,Type;
  int NumberOfMolecules;
  VECTOR pos;
  int NS,NR;
  REAL QS,QR;

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

  CurrentSystem=(int)(RandomNumber()*NumberOfSystems);

  NumberOfMolecules=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    NumberOfMolecules+=Framework[CurrentSystem].TotalNumberOfAtoms;
  if(NumberOfMolecules==0) return 0;

  ChiralInversionAttempts[CurrentSystem]+=1.0;

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

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++) 
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfGroups;i++)
    {
      Adsorbates[CurrentSystem][m].Groups[i].CenterOfMassReferencePosition=Adsorbates[CurrentSystem][m].Groups[i].CenterOfMassPosition;
      pos=ChiralInversion(Framework[CurrentSystem].SpaceGroupIdentifier[0],Adsorbates[CurrentSystem][m].Groups[i].CenterOfMassPosition);
      Adsorbates[CurrentSystem][m].Groups[i].CenterOfMassPosition=pos;
    }
    for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][m].Atoms[j].ReferencePosition=Adsorbates[CurrentSystem][m].Atoms[j].Position;
      pos=ChiralInversion(Framework[CurrentSystem].SpaceGroupIdentifier[0],Adsorbates[CurrentSystem][m].Atoms[j].Position);
      Adsorbates[CurrentSystem][m].Atoms[j].Position=pos;
    }
  }
 
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(i=0;i<Components[Type].NumberOfGroups;i++)
    {
      Cations[CurrentSystem][m].Groups[i].CenterOfMassReferencePosition=Cations[CurrentSystem][m].Groups[i].CenterOfMassPosition;
      pos=ChiralInversion(Framework[CurrentSystem].SpaceGroupIdentifier[0],Cations[CurrentSystem][m].Groups[i].CenterOfMassPosition);
      Cations[CurrentSystem][m].Groups[i].CenterOfMassPosition=pos;
    }
    for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][m].Atoms[j].ReferencePosition=Cations[CurrentSystem][m].Atoms[j].Position;
      pos=ChiralInversion(Framework[CurrentSystem].SpaceGroupIdentifier[0],Cations[CurrentSystem][m].Atoms[j].Position);
      Cations[CurrentSystem][m].Atoms[j].Position=pos;
    }
  }

  // NOTE: CalculateEnergy does not recompute the positions from the center-of-mass and orientation (otherwise the type needs to be switched too).
  CalculateEnergy();

  if(fabs(StoredUTotal-UTotal[CurrentSystem])>1e-3)
  {
    for(i=0;i<NumberOfSystems;i++)
    {
      fprintf(OutputFilePtr[i],"\n\n");
      fprintf(OutputFilePtr[i],"ERROR [Chiral inversion move]: energy before (%g) and after (%g) are different in system %d\n",StoredUTotal,UTotal[CurrentSystem],CurrentSystem);
      fclose(OutputFilePtr[i]);
    }
    exit(0);
  }

  QS=Components[ParallelMolFractionComponentA].MolFraction[CurrentSystem];
  QR=Components[ParallelMolFractionComponentB].MolFraction[CurrentSystem];
  NS=Components[ParallelMolFractionComponentA].NumberOfMolecules[CurrentSystem];
  NR=Components[ParallelMolFractionComponentB].NumberOfMolecules[CurrentSystem];

  if(RandomNumber()<pow(QS,NR)*pow(QR,NS)/(pow(QS,NS)*pow(QR,NR)))
  {
    ChiralInversionAccepted[CurrentSystem]+=1.0;

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++) 
    {
      if(Adsorbates[CurrentSystem][i].Type==ParallelMolFractionComponentA)
      {
        Adsorbates[CurrentSystem][i].Type=ParallelMolFractionComponentB;
        Components[ParallelMolFractionComponentA].NumberOfMolecules[CurrentSystem]--;
        Components[ParallelMolFractionComponentB].NumberOfMolecules[CurrentSystem]++;
      }
      else if(Adsorbates[CurrentSystem][i].Type==ParallelMolFractionComponentB)
      {
        Adsorbates[CurrentSystem][i].Type=ParallelMolFractionComponentA;
        Components[ParallelMolFractionComponentA].NumberOfMolecules[CurrentSystem]++;
        Components[ParallelMolFractionComponentB].NumberOfMolecules[CurrentSystem]--;
      }
    }

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++) 
    {
      if(Cations[CurrentSystem][i].Type==ParallelMolFractionComponentA)
      {
        Cations[CurrentSystem][i].Type=ParallelMolFractionComponentB;
        Components[ParallelMolFractionComponentA].NumberOfMolecules[CurrentSystem]--;
        Components[ParallelMolFractionComponentB].NumberOfMolecules[CurrentSystem]++;
      }
      else if(Cations[CurrentSystem][i].Type==ParallelMolFractionComponentB)
      {
        Cations[CurrentSystem][i].Type=ParallelMolFractionComponentA;
        Components[ParallelMolFractionComponentA].NumberOfMolecules[CurrentSystem]++;
        Components[ParallelMolFractionComponentB].NumberOfMolecules[CurrentSystem]--;
      }
    }
  }
  else
  {
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

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(i=0;i<Components[Type].NumberOfGroups;i++)
        Adsorbates[CurrentSystem][m].Groups[i].CenterOfMassPosition=Adsorbates[CurrentSystem][m].Groups[i].CenterOfMassReferencePosition;
      for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
        Adsorbates[CurrentSystem][m].Atoms[j].Position=Adsorbates[CurrentSystem][m].Atoms[j].ReferencePosition;
    }

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      Type=Cations[CurrentSystem][m].Type;
      for(i=0;i<Components[Type].NumberOfGroups;i++)
        Cations[CurrentSystem][m].Groups[i].CenterOfMassPosition=Cations[CurrentSystem][m].Groups[i].CenterOfMassReferencePosition;
      for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
        Cations[CurrentSystem][m].Atoms[j].Position=Cations[CurrentSystem][m].Atoms[j].ReferencePosition;
    }

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      RetrieveStoredEwaldStructureFactors(0,CurrentSystem);

    CalculateAnisotropicSites();
  }
  return 0;
}

void PrintChiralInversionStatistics(FILE *FilePtr)
{
  if(ProbabilityChiralInversionMove>0.0)
  {
    fprintf(FilePtr,"Performance of the chiral-inversion move:\n");
    fprintf(FilePtr,"=========================================\n");
    if(ChiralInversionAttempts[CurrentSystem]>0.0)
    {
       fprintf(FilePtr,"total tried: %lf accepted: %lf (%lf [%%])\n",
          (double)ChiralInversionAttempts[CurrentSystem],
          (double)ChiralInversionAccepted[CurrentSystem],
          (double)(ChiralInversionAttempts[CurrentSystem]>(REAL)0.0?
            100.0*ChiralInversionAccepted[CurrentSystem]/ChiralInversionAttempts[CurrentSystem]:(REAL)0.0));
    }
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"Chiral inversion move was OFF\n\n");
}



/*********************************************************************************************************
 * Name       | FrameworkChangeMove                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to change the positions of frameworkatoms                               *
 * Parameters | -                                                                                        *
 * Note       | The number of framework atoms that are randomly chosen is the minimum of 500 and the     *
 *            | number of framework atoms.                                                               *
 *********************************************************************************************************/

// FIX: CF compatibility

int FrameworkChangeMove(void)
{
  int i,index,k,Type,f1,A,B;
  REAL StoredUHostBond,StoredUHostBend,StoredUHostUreyBradley,StoredUHostTorsion,StoredUHostImproperTorsion,StoredUHostBendTorsion;

  REAL StoredUHostHost,StoredUHostHostVDW,StoredUHostHostCoulomb;
  REAL StoredUHostHostChargeChargeReal,StoredUHostHostChargeBondDipoleReal,StoredUHostHostBondDipoleBondDipoleReal;
  REAL StoredUHostHostChargeChargeFourier,StoredUHostHostChargeBondDipoleFourier,StoredUHostHostBondDipoleBondDipoleFourier;

  REAL StoredUHostAdsorbate,StoredUHostAdsorbateVDW,StoredUHostAdsorbateCoulomb;
  REAL StoredUHostAdsorbateChargeChargeReal,StoredUHostAdsorbateChargeBondDipoleReal,StoredUHostAdsorbateBondDipoleBondDipoleReal;
  REAL StoredUHostAdsorbateChargeChargeFourier,StoredUHostAdsorbateChargeBondDipoleFourier,StoredUHostAdsorbateBondDipoleBondDipoleFourier;

  REAL StoredUHostCation,StoredUHostCationVDW,StoredUHostCationCoulomb;
  REAL StoredUHostCationChargeChargeReal,StoredUHostCationChargeBondDipoleReal,StoredUHostCationBondDipoleBondDipoleReal;
  REAL StoredUHostCationChargeChargeFourier,StoredUHostCationChargeBondDipoleFourier,StoredUHostCationBondDipoleBondDipoleFourier;

  REAL UHostBondNew,UHostUreyBradleyNew,UHostBendNew,UHostInversionBendNew,UHostTorsionNew,UHostImproperTorsionNew;
  REAL UHostBondBondNew,UHostBondBendNew,UHostBendBendNew,UHostBondTorsionNew,UHostBendTorsionNew;
  REAL UHostBondOld,UHostUreyBradleyOld,UHostBendOld,UHostInversionBendOld,UHostTorsionOld,UHostImproperTorsionOld;
  REAL UHostBondBondOld,UHostBondBendOld,UHostBendBendOld,UHostBondTorsionOld,UHostBendTorsionOld;
  REAL UAdsorbateVDWNew,UAdsorbateVDWOld,UAdsorbateChargeChargeNew,UAdsorbateChargeChargeOld;
  REAL UCationVDWNew,UCationVDWOld,UCationChargeChargeNew,UCationChargeChargeOld;
  REAL UAdsorbateChargeBondDipoleNew,UAdsorbateChargeBondDipoleOld;
  REAL UAdsorbateBondDipoleChargeNew,UAdsorbateBondDipoleChargeOld;
  REAL UAdsorbateBondDipoleBondDipoleNew,UAdsorbateBondDipoleBondDipoleOld;
  REAL UCationChargeBondDipoleNew,UCationChargeBondDipoleOld;
  REAL UCationBondDipoleChargeNew,UCationBondDipoleChargeOld;
  REAL UCationBondDipoleBondDipoleNew,UCationBondDipoleBondDipoleOld;
  REAL UHostHostVDWNew,UHostHostVDWOld;
  REAL Magnitude,temp,temp2,DeltaU;
  VECTOR posA,posB;
  int AnisotropicNeighboringAtoms;
  int Atoms[20];
  int NumberOfAtoms,TypeA;


  // choose a system at random
  CurrentSystem=(int)(RandomNumber()*(REAL)NumberOfSystems);

  if(Framework[CurrentSystem].FrameworkModel!=FLEXIBLE) return 0;

  for(k=0;k<MIN2(500,Framework[CurrentSystem].TotalNumberOfAtoms);k++)
  {
    // choose a framework at random
    CurrentFramework=(int)(RandomNumber()*Framework[CurrentSystem].NumberOfFrameworks);

    if(Framework[CurrentSystem].FrameworkModels[CurrentFramework]!=FLEXIBLE) continue;

    // choose an atom within that framework at random
    index=(int)(RandomNumber()*(REAL)Framework[CurrentSystem].NumberOfAtoms[CurrentFramework]);

    UHostHostVDWNew=UHostHostVDWOld=0.0;
    UAdsorbateVDWNew=UAdsorbateVDWOld=0.0;
    UAdsorbateChargeChargeNew=UAdsorbateChargeChargeOld=0.0;
    UAdsorbateChargeBondDipoleNew=UAdsorbateChargeBondDipoleOld=0.0;
    UAdsorbateBondDipoleChargeNew=UAdsorbateBondDipoleChargeOld=0.0;
    UAdsorbateBondDipoleBondDipoleNew=UAdsorbateBondDipoleBondDipoleOld=0.0;
    UCationVDWNew=UCationVDWOld=0.0;
    UCationChargeChargeNew=UCationChargeChargeOld=0.0;
    UCationChargeBondDipoleNew=UCationChargeBondDipoleOld=0.0;
    UCationBondDipoleChargeNew=UCationBondDipoleChargeOld=0.0;
    UCationBondDipoleBondDipoleNew=UCationBondDipoleBondDipoleOld=0.0;

    UHostVDWDelta[CurrentSystem]=UHostChargeChargeRealDelta[CurrentSystem]=UHostChargeBondDipoleRealDelta[CurrentSystem]=UHostBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;
    UHostHostChargeChargeFourierDelta[CurrentSystem]=UHostHostChargeBondDipoleFourierDelta[CurrentSystem]=UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;


    Type=Framework[CurrentSystem].Atoms[CurrentFramework][index].Type;
    FrameworkChangeAttempts[CurrentSystem][Type]+=1.0;

    // store the current energies
    StoredUHostBond=UHostBond[CurrentSystem];
    StoredUHostBend=UHostBend[CurrentSystem];
    StoredUHostUreyBradley=UHostUreyBradley[CurrentSystem];
    StoredUHostTorsion=UHostTorsion[CurrentSystem];
    StoredUHostImproperTorsion=UHostImproperTorsion[CurrentSystem];
    StoredUHostBendTorsion=UHostBendTorsion[CurrentSystem];

    StoredUHostHost=UHostHost[CurrentSystem];
    StoredUHostHostVDW=UHostHostVDW[CurrentSystem];
    StoredUHostHostCoulomb=UHostHostCoulomb[CurrentSystem];
    StoredUHostHostChargeChargeReal=UHostHostChargeChargeReal[CurrentSystem];
    StoredUHostHostChargeBondDipoleReal=UHostHostChargeBondDipoleReal[CurrentSystem];
    StoredUHostHostBondDipoleBondDipoleReal=UHostHostBondDipoleBondDipoleReal[CurrentSystem];
    StoredUHostHostChargeChargeFourier=UHostHostChargeChargeFourier[CurrentSystem];
    StoredUHostHostChargeBondDipoleFourier=UHostHostChargeBondDipoleFourier[CurrentSystem];
    StoredUHostHostBondDipoleBondDipoleFourier=UHostHostBondDipoleBondDipoleFourier[CurrentSystem];

    StoredUHostAdsorbate=UHostAdsorbate[CurrentSystem];
    StoredUHostAdsorbateVDW=UHostAdsorbateVDW[CurrentSystem];
    StoredUHostAdsorbateCoulomb=UHostAdsorbateCoulomb[CurrentSystem];
    StoredUHostAdsorbateChargeChargeReal=UHostAdsorbateChargeChargeReal[CurrentSystem];
    StoredUHostAdsorbateChargeBondDipoleReal=UHostAdsorbateChargeBondDipoleReal[CurrentSystem];
    StoredUHostAdsorbateBondDipoleBondDipoleReal=UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
    StoredUHostAdsorbateChargeChargeFourier=UHostAdsorbateChargeChargeFourier[CurrentSystem];
    StoredUHostAdsorbateChargeBondDipoleFourier=UHostAdsorbateChargeBondDipoleFourier[CurrentSystem];
    StoredUHostAdsorbateBondDipoleBondDipoleFourier=UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];

    StoredUHostCation=UHostCation[CurrentSystem];
    StoredUHostCationVDW=UHostCationVDW[CurrentSystem];
    StoredUHostCationCoulomb=UHostCationCoulomb[CurrentSystem];
    StoredUHostCationChargeChargeReal=UHostCationChargeChargeReal[CurrentSystem];
    StoredUHostCationChargeBondDipoleReal=UHostCationChargeBondDipoleReal[CurrentSystem];
    StoredUHostCationBondDipoleBondDipoleReal=UHostCationBondDipoleBondDipoleReal[CurrentSystem];
    StoredUHostCationChargeChargeFourier=UHostCationChargeChargeFourier[CurrentSystem];
    StoredUHostCationChargeBondDipoleFourier=UHostCationChargeBondDipoleFourier[CurrentSystem];
    StoredUHostCationBondDipoleBondDipoleFourier=UHostCationBondDipoleBondDipoleFourier[CurrentSystem];


    // save the selected atom
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Framework[CurrentSystem].Atoms[f1][i].ReferencePosition=Framework[CurrentSystem].Atoms[f1][i].Position;
        Framework[CurrentSystem].Atoms[f1][i].ReferenceAnisotropicPosition=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;
      }

    UHostBondOld=CalculateFrameworkBondEnergy(FALSE,CurrentFramework,index);
    UHostUreyBradleyOld=CalculateFrameworkUreyBradleyEnergy(FALSE,CurrentFramework,index);
    UHostBendOld=CalculateFrameworkBendEnergy(FALSE,CurrentFramework,index);
    UHostInversionBendOld=CalculateFrameworkInversionBendEnergy(FALSE,CurrentFramework,index);
    UHostTorsionOld=CalculateFrameworkTorsionEnergy(FALSE,CurrentFramework,index);
    UHostImproperTorsionOld=CalculateFrameworkImproperTorsionEnergy(FALSE,CurrentFramework,index);
    UHostBondBondOld=CalculateFrameworkBondBondEnergy(FALSE,CurrentFramework,index);
    UHostBondBendOld=CalculateFrameworkBondBendEnergy(FALSE,CurrentFramework,index);
    UHostBendBendOld=CalculateFrameworkBendBendEnergy(FALSE,CurrentFramework,index);
    UHostBondTorsionOld=CalculateFrameworkBondTorsionEnergy(FALSE,CurrentFramework,index);
    UHostBendTorsionOld=CalculateFrameworkBendTorsionEnergy(FALSE,CurrentFramework,index);

    // displace the selected atom
    Framework[CurrentSystem].Atoms[CurrentFramework][index].Position.x+=0.1*(2.0*RandomNumber()-1.0);
    Framework[CurrentSystem].Atoms[CurrentFramework][index].Position.y+=0.1*(2.0*RandomNumber()-1.0);
    Framework[CurrentSystem].Atoms[CurrentFramework][index].Position.z+=0.1*(2.0*RandomNumber()-1.0);

    CalculateAnisotropicSites();

    UHostBondNew=CalculateFrameworkBondEnergy(FALSE,CurrentFramework,index);
    UHostUreyBradleyNew=CalculateFrameworkUreyBradleyEnergy(FALSE,CurrentFramework,index);
    UHostBendNew=CalculateFrameworkBendEnergy(FALSE,CurrentFramework,index);
    UHostInversionBendNew=CalculateFrameworkInversionBendEnergy(FALSE,CurrentFramework,index);
    UHostTorsionNew=CalculateFrameworkTorsionEnergy(FALSE,CurrentFramework,index);
    UHostImproperTorsionNew=CalculateFrameworkImproperTorsionEnergy(FALSE,CurrentFramework,index);
    UHostBondBondNew=CalculateFrameworkBondBondEnergy(FALSE,CurrentFramework,index);
    UHostBondBendNew=CalculateFrameworkBondBendEnergy(FALSE,CurrentFramework,index);
    UHostBendBendNew=CalculateFrameworkBendBendEnergy(FALSE,CurrentFramework,index);
    UHostBondTorsionNew=CalculateFrameworkBondTorsionEnergy(FALSE,CurrentFramework,index);
    UHostBendTorsionNew=CalculateFrameworkBendTorsionEnergy(FALSE,CurrentFramework,index);

    NumberOfAtoms=1;
    Atoms[0]=index;

    Type=Framework[CurrentSystem].Atoms[CurrentFramework][index].Type;
    AnisotropicNeighboringAtoms=PseudoAtoms[Type].AnisotropicCorrection;
    for(i=0;i<Framework[CurrentSystem].Connectivity[CurrentFramework][index];i++)
    {
      TypeA=Framework[CurrentSystem].Atoms[CurrentFramework][Framework[CurrentSystem].Neighbours[CurrentFramework][index][i]].Type;
      if(PseudoAtoms[TypeA].AnisotropicCorrection) AnisotropicNeighboringAtoms=TRUE;
    }

    // add neighbors for anisotropic sites
    if(AnisotropicNeighboringAtoms)
    {
      for(i=0;i<Framework[CurrentSystem].Connectivity[CurrentFramework][index];i++)
        Atoms[NumberOfAtoms++]=Framework[CurrentSystem].Neighbours[CurrentFramework][index][i];
    }

    for(i=0;i<NumberOfAtoms;i++)
    {
      index=Atoms[i];
      TypeA=Framework[CurrentSystem].Atoms[CurrentFramework][index].Type;

      UAdsorbateVDWNew+=CalculateInterVDWEnergyAdsorbateAtPosition(Framework[CurrentSystem].Atoms[CurrentFramework][index].AnisotropicPosition,TypeA,-1,1.0);
      UAdsorbateVDWOld+=CalculateInterVDWEnergyAdsorbateAtPosition(Framework[CurrentSystem].Atoms[CurrentFramework][index].ReferenceAnisotropicPosition,TypeA,-1,1.0);

      UCationVDWNew+=CalculateInterVDWEnergyCationAtPosition(Framework[CurrentSystem].Atoms[CurrentFramework][index].AnisotropicPosition,TypeA,-1,1.0);
      UCationVDWOld+=CalculateInterVDWEnergyCationAtPosition(Framework[CurrentSystem].Atoms[CurrentFramework][index].ReferenceAnisotropicPosition,TypeA,-1,1.0);

      UHostHostVDWNew+=CalculateEnergyDifferenceFrameworkMoveVDW(index,Framework[CurrentSystem].Atoms[CurrentFramework][index].AnisotropicPosition,TypeA);
      UHostHostVDWOld+=CalculateEnergyDifferenceFrameworkMoveVDW(index,Framework[CurrentSystem].Atoms[CurrentFramework][index].ReferenceAnisotropicPosition,TypeA);
    }

    index=Atoms[0];

    CalculateInterChargeEnergyAdsorbateAtPosition(Framework[CurrentSystem].Atoms[CurrentFramework][index].Position,Type,&UAdsorbateChargeChargeNew,&UAdsorbateChargeBondDipoleNew,-1,1.0);
    CalculateInterChargeEnergyAdsorbateAtPosition(Framework[CurrentSystem].Atoms[CurrentFramework][index].ReferencePosition,Type,&UAdsorbateChargeChargeOld,&UAdsorbateChargeBondDipoleOld,-1,1.0);

    CalculateInterChargeEnergyCationAtPosition(Framework[CurrentSystem].Atoms[CurrentFramework][index].Position,Type,&UCationChargeChargeNew,&UCationChargeBondDipoleNew,-1,1.0);
    CalculateInterChargeEnergyCationAtPosition(Framework[CurrentSystem].Atoms[CurrentFramework][index].ReferencePosition,Type,&UCationChargeChargeOld,&UCationChargeBondDipoleOld,-1,1.0);

    CalculateEnergyDifferenceFrameworkMoveCharge(index);

    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
    {
      Magnitude=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
      A=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
      B=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;
      if((A==index)||(B==index))
      {
        posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
        posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;
        CalculateInterBondDipoleEnergyAdsorbateAtPosition(posA,posB,Magnitude,&temp,&temp2,-1);
        UAdsorbateBondDipoleChargeNew+=temp;
        UAdsorbateBondDipoleBondDipoleNew+=temp2;
        CalculateInterBondDipoleEnergyCationAtPosition(posA,posB,Magnitude,&temp,&temp2,-1);
        UCationBondDipoleChargeNew+=temp;
        UCationBondDipoleBondDipoleNew+=temp2;


        posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].ReferencePosition;
        posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].ReferencePosition;
        CalculateInterBondDipoleEnergyAdsorbateAtPosition(posA,posB,Magnitude,&temp,&temp2,-1);
        UAdsorbateBondDipoleChargeOld+=temp;
        UAdsorbateBondDipoleBondDipoleOld+=temp2;
        CalculateInterBondDipoleEnergyCationAtPosition(posA,posB,Magnitude,&temp,&temp2,-1);
        UCationBondDipoleChargeOld+=temp;
        UCationBondDipoleBondDipoleOld+=temp2;
      }
    }
   

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
       CalculateEwaldFourierFrameworkAtomTranslate(index);


    DeltaU=(UHostBondNew-UHostBondOld)+(UHostUreyBradleyNew-UHostUreyBradleyOld)+(UHostBendNew-UHostBendOld)+
           (UHostInversionBendNew-UHostInversionBendOld)+(UHostTorsionNew-UHostTorsionOld)+(UHostImproperTorsionNew-UHostImproperTorsionOld)+
           (UHostBondBondNew-UHostBondBondOld)+(UHostBondBendNew-UHostBondBendOld)+(UHostBendBendNew-UHostBendBendOld)+
           (UHostBondTorsionNew-UHostBondTorsionOld)+(UHostBendTorsionNew-UHostBendTorsionOld)+
           (UHostHostVDWNew-UHostHostVDWOld)+(UAdsorbateVDWNew-UAdsorbateVDWOld)+(UCationVDWNew-UCationVDWOld)+
           UHostHostChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
           UHostChargeBondDipoleRealDelta[CurrentSystem]+UHostHostChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           (UAdsorbateChargeChargeNew-UAdsorbateChargeChargeOld)+UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
           (UAdsorbateChargeBondDipoleNew-UAdsorbateChargeBondDipoleOld)+UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
           (UAdsorbateBondDipoleChargeNew-UAdsorbateBondDipoleChargeOld)+
           (UAdsorbateBondDipoleBondDipoleNew-UAdsorbateBondDipoleBondDipoleOld)+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           (UCationChargeChargeNew-UCationChargeChargeOld)+UHostCationChargeChargeFourierDelta[CurrentSystem]+
           (UCationChargeBondDipoleNew-UCationChargeBondDipoleOld)+UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
           (UCationBondDipoleChargeNew-UCationBondDipoleChargeOld)+
           (UCationBondDipoleBondDipoleNew-UCationBondDipoleBondDipoleOld)+UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];


    if(RandomNumber()<exp(-Beta[CurrentSystem]*DeltaU))
    {
      FrameworkChangeAccepted[CurrentSystem][Type]+=1.0;

      UTotal[CurrentSystem]+=DeltaU;

      UHostHostVDW[CurrentSystem]+=UHostHostVDWNew-UHostHostVDWOld;
      UHostHost[CurrentSystem]+=UHostHostVDWNew-UHostHostVDWOld;

      UHostAdsorbateVDW[CurrentSystem]+=(UAdsorbateVDWNew-UAdsorbateVDWOld);
      UHostAdsorbate[CurrentSystem]+=(UAdsorbateVDWNew-UAdsorbateVDWOld);

      UHostCationVDW[CurrentSystem]+=(UCationVDWNew-UCationVDWOld);
      UHostCation[CurrentSystem]+=(UCationVDWNew-UCationVDWOld);

      UHostBond[CurrentSystem]+=(UHostBondNew-UHostBondOld);
      UHostUreyBradley[CurrentSystem]+=(UHostUreyBradleyNew-UHostUreyBradleyOld);
      UHostBend[CurrentSystem]+=(UHostBendNew-UHostBendOld);
      UHostInversionBend[CurrentSystem]+=(UHostInversionBendNew-UHostInversionBendOld);
      UHostTorsion[CurrentSystem]+=(UHostTorsionNew-UHostTorsionOld);
      UHostImproperTorsion[CurrentSystem]+=(UHostImproperTorsionNew-UHostImproperTorsionOld);
      UHostBondBond[CurrentSystem]+=(UHostBondBondNew-UHostBondBondOld);
      UHostBondBend[CurrentSystem]+=(UHostBondBendNew-UHostBondBendOld);
      UHostBendBend[CurrentSystem]+=(UHostBendBendNew-UHostBendBendOld);
      UHostBondTorsion[CurrentSystem]+=(UHostBondTorsionNew-UHostBondTorsionOld);
      UHostBendTorsion[CurrentSystem]+=(UHostBendTorsionNew-UHostBendTorsionOld);

      if(ChargeMethod!=NONE)
      {
        UHostAdsorbateChargeChargeReal[CurrentSystem]+=(UAdsorbateChargeChargeNew-UAdsorbateChargeChargeOld);
        UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=(UAdsorbateChargeBondDipoleNew-UAdsorbateChargeBondDipoleOld)+(UAdsorbateBondDipoleChargeNew-UAdsorbateBondDipoleChargeOld);
        UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=(UAdsorbateBondDipoleBondDipoleNew-UAdsorbateBondDipoleBondDipoleOld);
        UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
        UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
        UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
        UHostAdsorbateCoulomb[CurrentSystem]+=(UAdsorbateChargeChargeNew-UAdsorbateChargeChargeOld)+UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                              (UAdsorbateChargeBondDipoleNew-UAdsorbateChargeBondDipoleOld)+UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                              (UAdsorbateBondDipoleChargeNew-UAdsorbateBondDipoleChargeOld)+
                                              (UAdsorbateBondDipoleBondDipoleNew-UAdsorbateBondDipoleBondDipoleOld)+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
        UHostAdsorbate[CurrentSystem]+=(UAdsorbateChargeChargeNew-UAdsorbateChargeChargeOld)+UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                       (UAdsorbateChargeBondDipoleNew-UAdsorbateChargeBondDipoleOld)+UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                       (UAdsorbateBondDipoleChargeNew-UAdsorbateBondDipoleChargeOld)+
                                       (UAdsorbateBondDipoleBondDipoleNew-UAdsorbateBondDipoleBondDipoleOld)+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

        UHostCationChargeChargeReal[CurrentSystem]+=(UCationChargeChargeNew-UCationChargeChargeOld);
        UHostCationChargeBondDipoleReal[CurrentSystem]+=(UCationChargeBondDipoleNew-UCationChargeBondDipoleOld)+(UCationBondDipoleChargeNew-UCationBondDipoleChargeOld);
        UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=(UCationBondDipoleBondDipoleNew-UCationBondDipoleBondDipoleOld);
        UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
        UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
        UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
        UHostCationCoulomb[CurrentSystem]+=(UCationChargeChargeNew-UCationChargeChargeOld)+UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                           (UCationChargeBondDipoleNew-UCationChargeBondDipoleOld)+UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           (UCationBondDipoleChargeNew-UCationBondDipoleChargeOld)+
                                           (UCationBondDipoleBondDipoleNew-UCationBondDipoleBondDipoleOld)+UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
        UHostCation[CurrentSystem]+=(UCationChargeChargeNew-UCationChargeChargeOld)+UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                    (UCationChargeBondDipoleNew-UCationChargeBondDipoleOld)+UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    (UCationBondDipoleChargeNew-UCationBondDipoleChargeOld)+
                                    (UCationBondDipoleBondDipoleNew-UCationBondDipoleBondDipoleOld)+UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];


        UHostHostChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem];
        UHostHostChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDelta[CurrentSystem];
        UHostHostBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
        UHostHostChargeChargeFourier[CurrentSystem]+=UHostHostChargeChargeFourierDelta[CurrentSystem];
        UHostHostChargeBondDipoleFourier[CurrentSystem]+=UHostHostChargeBondDipoleFourierDelta[CurrentSystem];
        UHostHostBondDipoleBondDipoleFourier[CurrentSystem]+=UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem];
        UHostHostCoulomb[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem]+UHostHostChargeChargeFourierDelta[CurrentSystem]+
                                         UHostChargeBondDipoleRealDelta[CurrentSystem]+UHostHostChargeBondDipoleFourierDelta[CurrentSystem]+
                                         UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem];
        UHostHost[CurrentSystem]+=UHostHostChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                  UHostChargeBondDipoleRealDelta[CurrentSystem]+UHostHostChargeBondDipoleFourierDelta[CurrentSystem]+
                                  UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem];

        if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
          AcceptEwaldFrameworkMove(0);
      }
    }
    else
    {
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          Framework[CurrentSystem].Atoms[f1][i].Position=Framework[CurrentSystem].Atoms[f1][i].ReferencePosition;
          Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition=Framework[CurrentSystem].Atoms[f1][i].ReferenceAnisotropicPosition;
        }

      UHostBond[CurrentSystem]=StoredUHostBond;
      UHostBend[CurrentSystem]=StoredUHostBend;
      UHostUreyBradley[CurrentSystem]=StoredUHostUreyBradley;
      UHostTorsion[CurrentSystem]=StoredUHostTorsion;
      UHostImproperTorsion[CurrentSystem]=StoredUHostImproperTorsion;
      UHostBendTorsion[CurrentSystem]=StoredUHostBendTorsion;
  
      UHostHost[CurrentSystem]=StoredUHostHost;
      UHostHostVDW[CurrentSystem]=StoredUHostHostVDW;
      UHostHostCoulomb[CurrentSystem]=StoredUHostHostCoulomb;
      UHostHostChargeChargeReal[CurrentSystem]=StoredUHostHostChargeChargeReal;
      UHostHostChargeBondDipoleReal[CurrentSystem]=StoredUHostHostChargeBondDipoleReal;
      UHostHostBondDipoleBondDipoleReal[CurrentSystem]=StoredUHostHostBondDipoleBondDipoleReal;
      UHostHostChargeChargeFourier[CurrentSystem]=StoredUHostHostChargeChargeFourier;
      UHostHostChargeBondDipoleFourier[CurrentSystem]=StoredUHostHostChargeBondDipoleFourier;
      UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=StoredUHostHostBondDipoleBondDipoleFourier;

      UHostAdsorbate[CurrentSystem]=StoredUHostAdsorbate;
      UHostAdsorbateVDW[CurrentSystem]=StoredUHostAdsorbateVDW;
      UHostAdsorbateCoulomb[CurrentSystem]=StoredUHostAdsorbateCoulomb;
      UHostAdsorbateChargeChargeReal[CurrentSystem]=StoredUHostAdsorbateChargeChargeReal;
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]=StoredUHostAdsorbateChargeBondDipoleReal;
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=StoredUHostAdsorbateBondDipoleBondDipoleReal;
      UHostAdsorbateChargeChargeFourier[CurrentSystem]=StoredUHostAdsorbateChargeChargeFourier;
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=StoredUHostAdsorbateChargeBondDipoleFourier;
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=StoredUHostAdsorbateBondDipoleBondDipoleFourier;

      UHostCation[CurrentSystem]=StoredUHostCation;
      UHostCationVDW[CurrentSystem]=StoredUHostCationVDW;
      UHostCationCoulomb[CurrentSystem]=StoredUHostCationCoulomb;
      UHostCationChargeChargeReal[CurrentSystem]=StoredUHostCationChargeChargeReal;
      UHostCationChargeBondDipoleReal[CurrentSystem]=StoredUHostCationChargeBondDipoleReal;
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUHostCationBondDipoleBondDipoleReal;
      UHostCationChargeChargeFourier[CurrentSystem]=StoredUHostCationChargeChargeFourier;
      UHostCationChargeBondDipoleFourier[CurrentSystem]=StoredUHostCationChargeBondDipoleFourier;
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUHostCationBondDipoleBondDipoleFourier;
    }
    CalculateAnisotropicSites();
  }
  return 0;
}

void OptimizeFrameworkChangeAcceptence(void)
{
  int i;
  REAL ratio,vandr;

  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    if(FrameworkChangeAttempts[CurrentSystem][i]>0.0)
      ratio=FrameworkChangeAccepted[CurrentSystem][i]/FrameworkChangeAttempts[CurrentSystem][i];
    else
      ratio=0.0;

    vandr=ratio/TargetAccRatioTranslation;
    if(vandr>1.5) vandr=1.5;
    else if(vandr<0.5) vandr=0.5;
    FrameworkMaximumTranslation[CurrentSystem][i]*=vandr;
    if(FrameworkMaximumTranslation[CurrentSystem][i]<0.0005)
       FrameworkMaximumTranslation[CurrentSystem][i]=0.0005;
    if(FrameworkMaximumTranslation[CurrentSystem][i]>0.1)
       FrameworkMaximumTranslation[CurrentSystem][i]=0.1;

    TotalFrameworkChangeAttempts[CurrentSystem][i]+=FrameworkChangeAttempts[CurrentSystem][i];
    TotalFrameworkChangeAccepted[CurrentSystem][i]+=FrameworkChangeAccepted[CurrentSystem][i];
    FrameworkChangeAttempts[CurrentSystem][i]=FrameworkChangeAccepted[CurrentSystem][i]=0.0;
  }
}

void PrintFrameworkStatistics(FILE *FilePtr)
{
  int i;

  if(ProbabilityFrameworkChangeMove>0.0)
  {
    fprintf(FilePtr,"Performance of the framework change move:\n");
    fprintf(FilePtr,"=========================================\n");
    for(i=0;i<NumberOfPseudoAtoms;i++)
    {
      if(TotalFrameworkChangeAttempts[CurrentSystem][i]>0.0)
      {
        fprintf(FilePtr,"\tatom type %d [%s]\n",i,PseudoAtoms[i].Name);
        fprintf(FilePtr,"\t-----------------------------------------\n");
        fprintf(FilePtr,"\t\ttotal        %lf\n",(double)TotalFrameworkChangeAttempts[CurrentSystem][i]);
        fprintf(FilePtr,"\t\tsuccesfull   %lf\n",(double)TotalFrameworkChangeAccepted[CurrentSystem][i]);
        fprintf(FilePtr,"\t\taccepted   %lf (%lf [%%])\n",
          (double)TotalFrameworkChangeAccepted[CurrentSystem][i],
          (double)(TotalFrameworkChangeAttempts[CurrentSystem][i]>(REAL)0.0?
            100.0*TotalFrameworkChangeAccepted[CurrentSystem][i]/TotalFrameworkChangeAttempts[CurrentSystem][i]:(REAL)0.0));
        fprintf(FilePtr,"\t\tdisplacement %lf\n",(double)FrameworkMaximumTranslation[CurrentSystem][i]);
        fprintf(FilePtr,"\n");
      }
    }
  }
  else
    fprintf(FilePtr,"Framework change move was OFF\n\n");
}


/*********************************************************************************************************
 * Name       | FrameworkShiftMove                                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to change the relative positions of a framework.                        *
 * Parameters | -                                                                                        *
 * Note       | A framework (1,..,n) is randomly chosen and randomly displaced in either x, y, or z.     *
 *            | Framework 0 remains fixed.                                                               *
 *********************************************************************************************************/

// Framework[].Fixed: true) eveything is fixed false) some parts are flexible
// Framework[].Framework[].Fixed: true) framework is fixed, false) framework is flexible

// FIX: CF compatibility

int FrameworkShiftMove(void)
{
  int i;
  int nr_atoms;
  REAL choice,vNew;
  VECTOR displacement;
  REAL DeltaU;

  // no or a single framework do not need to be shifted
  if(Framework[CurrentSystem].NumberOfFrameworks<=1) return 0;

  // choose the framework randomly using the first framework as the fixed reference framework
  CurrentFramework=1+(int)((Framework[CurrentSystem].NumberOfFrameworks-1)*RandomNumber());

  // calculate a random displacement
  vNew=2.0*RandomNumber()-1.0;

  // initialize displacement to zero
  displacement.x=0.0;
  displacement.y=0.0;
  displacement.z=0.0;

  // choose a possible displacement
  switch(Framework[CurrentSystem].TranslationDirection)
  {
    default:
    case XYZ_DIR:
      choice=3.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentFramework].x;
      else if(choice<2.0)
        displacement.y=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentFramework].y;
      else
        displacement.z=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentFramework].z;
      break;
    case XY_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentFramework].x;
      else if(choice<2.0)
        displacement.y=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentFramework].y;
      break;
    case XZ_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.x=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentComponent].x;
      else if(choice<2.0)
        displacement.z=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentComponent].z;
      break;
    case YZ_DIR:
      choice=2.0*RandomNumber();
      if(choice<1.0)
        displacement.y=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentFramework].y;
      else if(choice<2.0)
        displacement.z=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentFramework].z;
      break;
    case X_DIR:
      displacement.x=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentFramework].x;
      break;
    case Y_DIR:
      displacement.y=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentFramework].y;
      break;
    case Z_DIR:
      displacement.z=vNew*FrameworkMaximumShiftTranslation[CurrentSystem][CurrentFramework].z;
      break;
  }

  // update statistics and register a translation move
  if(fabs(displacement.x)>0.0)
    FrameworkShiftAttempts[CurrentSystem][CurrentFramework].x+=1.0;
  if(fabs(displacement.y)>0.0)
    FrameworkShiftAttempts[CurrentSystem][CurrentFramework].y+=1.0;
  if(fabs(displacement.z)>0.0)
    FrameworkShiftAttempts[CurrentSystem][CurrentFramework].z+=1.0;

  nr_atoms=Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];
  for(i=0;i<nr_atoms;i++)
  {
    Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferencePosition=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
    Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferenceAnisotropicPosition=Framework[CurrentSystem].Atoms[CurrentFramework][i].AnisotropicPosition;
    Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.x+=displacement.x;
    Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.y+=displacement.y;
    Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.z+=displacement.z;
  }
  CalculateAnisotropicSites();

  CalculateFrameworkShiftEnergyDifferenceAdsorbateVDW();
  CalculateFrameworkShiftEnergyDifferenceAdsorbateCharge();

  CalculateFrameworkShiftEnergyDifferenceCationVDW();
  CalculateFrameworkShiftEnergyDifferenceCationCharge();

  CalculateFrameworkEnergyDifferenceShiftedFramework();

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
     CalculateEwaldFourierFrameworkDisplacement();

  DeltaU=UHostVDWDelta[CurrentSystem]+UAdsorbateVDWDelta[CurrentSystem]+UCationVDWDelta[CurrentSystem]+
         UHostChargeChargeRealDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+
         UAdsorbateChargeChargeRealDelta[CurrentSystem]+UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+
         UCationChargeChargeRealDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+
         UHostHostChargeChargeFourierDelta[CurrentSystem]+UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostCationChargeChargeFourierDelta[CurrentSystem]+
         UHostHostChargeBondDipoleFourierDelta[CurrentSystem]+UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
         UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

  if(RandomNumber()<exp(-Beta[CurrentSystem]*DeltaU))
  {
    if(fabs(displacement.x)>0.0)
      FrameworkShiftAccepted[CurrentSystem][CurrentFramework].x+=1.0;
    if(fabs(displacement.y)>0.0)
      FrameworkShiftAccepted[CurrentSystem][CurrentFramework].y+=1.0;
    if(fabs(displacement.z)>0.0)
      FrameworkShiftAccepted[CurrentSystem][CurrentFramework].z+=1.0;

    UHostHostVDW[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
    UHostHost[CurrentSystem]+=UHostVDWDelta[CurrentSystem];

    UHostAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];

    UHostCationVDW[CurrentSystem]+=UCationVDWDelta[CurrentSystem];
    UHostCation[CurrentSystem]+=UCationVDWDelta[CurrentSystem];

    if(ChargeMethod!=NONE)
    {
      UHostHostChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem];
      UHostHostChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDelta[CurrentSystem];
      UHostHostBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostHostChargeChargeFourier[CurrentSystem]+=UHostHostChargeChargeFourierDelta[CurrentSystem];
      UHostHostChargeBondDipoleFourier[CurrentSystem]+=UHostHostChargeBondDipoleFourierDelta[CurrentSystem];
      UHostHostBondDipoleBondDipoleFourier[CurrentSystem]+=UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostHostCoulomb[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem]+UHostHostChargeChargeFourierDelta[CurrentSystem]+
                                       UHostChargeBondDipoleRealDelta[CurrentSystem]+UHostHostChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostHost[CurrentSystem]+=UHostHostChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                UHostChargeBondDipoleRealDelta[CurrentSystem]+UHostHostChargeBondDipoleFourierDelta[CurrentSystem]+
                                UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem]+UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                            UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                            UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                     UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                     UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleRealDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostCationCoulomb[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem]+UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                         UCationChargeBondDipoleRealDelta[CurrentSystem]+UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                         UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                  UCationChargeBondDipoleRealDelta[CurrentSystem]+UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
    }

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      AcceptEwaldFrameworkDiplacementMove();

      if(Framework[CurrentSystem].FrameworkModels[CurrentFramework]!=FLEXIBLE)
        AcceptEwaldFrameworkDiplacementMoveRigid();
    }

    UTotal[CurrentSystem]+=DeltaU;
  }
  else
  {
    for(i=0;i<nr_atoms;i++)
    {
      Framework[CurrentSystem].Atoms[CurrentFramework][i].Position=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferencePosition;
      Framework[CurrentSystem].Atoms[CurrentFramework][i].AnisotropicPosition=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferenceAnisotropicPosition;
    }
  }

  return 0;
}

void OptimizeFrameworkShiftAcceptence(void)
{
  int i;
  VECTOR ratio,vandr;

  for(i=0;i<Framework[CurrentSystem].NumberOfFrameworks;i++)
  {
    if(FrameworkShiftAttempts[CurrentSystem][i].x>0.0)
      ratio.x=FrameworkShiftAccepted[CurrentSystem][i].x/FrameworkShiftAttempts[CurrentSystem][i].x;
    else
      ratio.x=0.0;
    if(FrameworkShiftAttempts[CurrentSystem][i].y>0.0)
      ratio.y=FrameworkShiftAccepted[CurrentSystem][i].y/FrameworkShiftAttempts[CurrentSystem][i].y;
    else
      ratio.y=0.0;
    if(FrameworkShiftAttempts[CurrentSystem][i].z>0.0)
      ratio.z=FrameworkShiftAccepted[CurrentSystem][i].z/FrameworkShiftAttempts[CurrentSystem][i].z;
    else
      ratio.z=0.0;

    vandr.x=ratio.x/TargetAccRatioTranslation;
    if(vandr.x>1.5) vandr.x=1.5;
    else if(vandr.x<0.5) vandr.x=0.5;
    FrameworkMaximumShiftTranslation[CurrentSystem][i].x*=vandr.x;
    if(FrameworkMaximumShiftTranslation[CurrentSystem][i].x<0.01)
       FrameworkMaximumShiftTranslation[CurrentSystem][i].x=0.01;
    if(FrameworkMaximumShiftTranslation[CurrentSystem][i].x>1.0)
       FrameworkMaximumShiftTranslation[CurrentSystem][i].x=1.0;

    vandr.y=ratio.y/TargetAccRatioTranslation;
    if(vandr.y>1.5) vandr.y=1.5;
    else if(vandr.y<0.5) vandr.y=0.5;
    FrameworkMaximumShiftTranslation[CurrentSystem][i].y*=vandr.y;
    if(FrameworkMaximumShiftTranslation[CurrentSystem][i].y<0.01)
       FrameworkMaximumShiftTranslation[CurrentSystem][i].y=0.01;
    if(FrameworkMaximumShiftTranslation[CurrentSystem][i].y>1.0)
       FrameworkMaximumShiftTranslation[CurrentSystem][i].y=1.0;

    vandr.z=ratio.z/TargetAccRatioTranslation;
     if(vandr.z>1.5) vandr.z=1.5;
    else if(vandr.z<0.5) vandr.z=0.5;
    FrameworkMaximumShiftTranslation[CurrentSystem][i].z*=vandr.z;
    if(FrameworkMaximumShiftTranslation[CurrentSystem][i].z<0.01)
       FrameworkMaximumShiftTranslation[CurrentSystem][i].z=0.01;
    if(FrameworkMaximumShiftTranslation[CurrentSystem][i].z>1.0)
       FrameworkMaximumShiftTranslation[CurrentSystem][i].z=1.0;

    TotalFrameworkShiftAttempts[CurrentSystem][i].x+=FrameworkShiftAttempts[CurrentSystem][i].x;
    TotalFrameworkShiftAccepted[CurrentSystem][i].x+=FrameworkShiftAccepted[CurrentSystem][i].x;
    FrameworkShiftAttempts[CurrentSystem][i].x=FrameworkShiftAccepted[CurrentSystem][i].x=0.0;

    TotalFrameworkShiftAttempts[CurrentSystem][i].y+=FrameworkShiftAttempts[CurrentSystem][i].y;
    TotalFrameworkShiftAccepted[CurrentSystem][i].y+=FrameworkShiftAccepted[CurrentSystem][i].y;
    FrameworkShiftAttempts[CurrentSystem][i].y=FrameworkShiftAccepted[CurrentSystem][i].y=0.0;

    TotalFrameworkShiftAttempts[CurrentSystem][i].z+=FrameworkShiftAttempts[CurrentSystem][i].z;
    TotalFrameworkShiftAccepted[CurrentSystem][i].z+=FrameworkShiftAccepted[CurrentSystem][i].z;
    FrameworkShiftAttempts[CurrentSystem][i].z=FrameworkShiftAccepted[CurrentSystem][i].z=0.0;
  }
}


void PrintFrameworkShiftStatistics(FILE *FilePtr)
{
  int i;

  if(ProbabilityFrameworkShiftMove>0.0)
  {
    fprintf(FilePtr,"Performance of the framework shift move:\n");
    fprintf(FilePtr,"=========================================\n");
    for(i=0;i<Framework[CurrentSystem].NumberOfFrameworks;i++)
    {
      fprintf(FilePtr,"Framework %d [%s]\n",i,Framework[CurrentSystem].Name[i]);
      fprintf(FilePtr,"\ttotal        %f %f %f\n",
        (double)TotalFrameworkShiftAttempts[CurrentSystem][i].x,
        (double)TotalFrameworkShiftAttempts[CurrentSystem][i].y,
        (double)TotalFrameworkShiftAttempts[CurrentSystem][i].z);
      fprintf(FilePtr,"\tsuccesfull   %f %f %f\n",
        (double)TotalFrameworkShiftAccepted[CurrentSystem][i].x,
        (double)TotalFrameworkShiftAccepted[CurrentSystem][i].y,
        (double)TotalFrameworkShiftAccepted[CurrentSystem][i].z);

      fprintf(FilePtr,"\taccepted   %f %f %f\n",
        (double)(TotalFrameworkShiftAttempts[CurrentSystem][i].x>(REAL)0.0?
           TotalFrameworkShiftAccepted[CurrentSystem][i].x/TotalFrameworkShiftAttempts[CurrentSystem][i].x:(REAL)0.0),
        (double)(TotalFrameworkShiftAttempts[CurrentSystem][i].y>(REAL)0.0?
            TotalFrameworkShiftAccepted[CurrentSystem][i].y/TotalFrameworkShiftAttempts[CurrentSystem][i].y:(REAL)0.0),
        (double)(TotalFrameworkShiftAttempts[CurrentSystem][i].z>(REAL)0.0?
            TotalFrameworkShiftAccepted[CurrentSystem][i].z/TotalFrameworkShiftAttempts[CurrentSystem][i].z:(REAL)0.0));
      fprintf(FilePtr,"\tdisplacement %f %f %f\n",
        (double)FrameworkMaximumShiftTranslation[CurrentSystem][i].x,
        (double)FrameworkMaximumShiftTranslation[CurrentSystem][i].y,
        (double)FrameworkMaximumShiftTranslation[CurrentSystem][i].z);

        fprintf(FilePtr,"\n");
    }
  }
  else
    fprintf(FilePtr,"Framework shift move was OFF\n\n");
}

/*********************************************************************************************************
 * Name       | GibbsParticleTransferAdsorbateMove                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to do a Gibbs particle transfer from one box to the other.              *
 * Parameters | -                                                                                        *
 * Note       | The ideal gas Rosenbluth weights cancel.                                                 *
 *********************************************************************************************************/

int GibbsParticleTransferAdsorbateMove(void)
{
  int i,A,B;
  REAL RosenbluthNew,UTailNew;
  REAL RosenbluthOld,UTailOld;
  REAL NetChargeDeltaNew,NetChargeDeltaOld;
  int CurrentSystemStored;
  REAL DeltaU;

  CurrentSystemStored=CurrentSystem;

  // choose the box random
  if(RandomNumber()<0.5)
  {
    A=0;
    B=1;
  }
  else
  {
    A=1;
    B=0;
  }

  GibbsSwapAttempts[CurrentComponent][A]+=1.0;

  // return if there are no particles of the chosen component in box B
  if(Components[CurrentComponent].NumberOfMolecules[B]<=Components[CurrentComponent].CFMoleculePresent[B]?1:0) return 0;

  // first grow a particle of the chosen type in box A
  CurrentSystem=A;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be added here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  CurrentAdsorbateMolecule=NumberOfAdsorbateMolecules[CurrentSystem];
  CurrentCationMolecule=NumberOfCationMolecules[CurrentSystem];
  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNew=GrowMolecule(CBMC_INSERTION);
  if (OVERLAP) return 0;

  UTailNew=TailMolecularEnergyDifferenceAdd();
  RosenbluthNew*=exp(-Beta[CurrentSystem]*UTailNew);

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CalculateEwaldFourierAdsorbate(TRUE,FALSE,NumberOfAdsorbateMolecules[CurrentSystem],A);
    NetChargeDeltaNew=NetChargeAdsorbateDelta;
    RosenbluthNew*=exp(-Beta[CurrentSystem]*(
        UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
        UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
        UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    RosenbluthNew*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }

  // second retrace a particle of the chosen type in box B
  CurrentSystem=B;

  // choose a random molecule of this component
  CurrentAdsorbateMolecule=SelectRandomMoleculeOfTypeExcludingFractionalMolecule(CurrentComponent);
  CurrentCationMolecule=-1;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be removed here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // calculate the Old Rosenbluth factor
  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_DELETION);

  UTailOld=TailMolecularEnergyDifferenceRemove();
  RosenbluthOld*=exp(-Beta[CurrentSystem]*UTailOld);

  if(ChargeMethod!=NONE)
  {
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      CalculateEwaldFourierAdsorbate(FALSE,TRUE,CurrentAdsorbateMolecule,B);
      NetChargeDeltaOld=NetChargeAdsorbateDelta;
      RosenbluthOld*=exp(Beta[CurrentSystem]*(
          UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
          UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
          UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
    }
  }


  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(FALSE,CurrentAdsorbateMolecule,-1);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    RosenbluthOld*=exp(Beta[CurrentSystem]*UDeltaPolarization);
  }


  // acceptence rule
  if(RandomNumber()<((RosenbluthNew*(REAL)(Components[CurrentComponent].NumberOfMolecules[B])*Volume[A])/
                    (RosenbluthOld*((REAL)(1.0+Components[CurrentComponent].NumberOfMolecules[A]))*Volume[B])))
  {
    GibbsSwapAccepted[CurrentComponent][A]+=1.0;

    CurrentSystem=A;

    UAdsorbateBond[CurrentSystem]+=UBondNew[CurrentSystem];
    UAdsorbateUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem];
    UAdsorbateBend[CurrentSystem]+=UBendNew[CurrentSystem];
    UAdsorbateBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem];
    UAdsorbateInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem];
    UAdsorbateTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem];
    UAdsorbateImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem];
    UAdsorbateBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem];
    UAdsorbateBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem];
    UAdsorbateBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem];
    UAdsorbateBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem];
    UAdsorbateIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem];
    UAdsorbateIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem];
    UAdsorbateIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem];
    UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem];

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UCationVDWNew[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem];
    UHostAdsorbate[CurrentSystem]+=UHostVDWNew[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    UTailCorrection[CurrentSystem]+=UTailNew;

    if(ChargeMethod!=NONE)
    {
      UAdsorbateIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem];
      UAdsorbateIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem];
      UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem];

      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]+
                                                 UAdsorbateChargeBondDipoleNew[CurrentSystem]+
                                                 UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+
                                                 UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                 UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]+
                                          UAdsorbateChargeBondDipoleNew[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+
                                          UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]+
                                            UHostChargeBondDipoleNew[CurrentSystem]+
                                            UHostBondDipoleBondDipoleNew[CurrentSystem]+
                                            UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]+
                                     UHostChargeBondDipoleNew[CurrentSystem]+
                                     UHostBondDipoleBondDipoleNew[CurrentSystem]+
                                     UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                     UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]+
                                              UCationChargeBondDipoleNew[CurrentSystem]+
                                              UCationBondDipoleBondDipoleNew[CurrentSystem]+
                                              UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]+
                                       UCationChargeBondDipoleNew[CurrentSystem]+
                                       UCationBondDipoleBondDipoleNew[CurrentSystem]+
                                       UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

      NetChargeAdsorbates[CurrentSystem]+=NetChargeDeltaNew;
      NetChargeSystem[CurrentSystem]+=NetChargeDeltaNew;

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(A);
    }

    InsertAdsorbateMolecule();

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
           UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization+UTailNew;

    UTotal[CurrentSystem]+=DeltaU;

    CurrentSystem=B;

    UAdsorbateBond[CurrentSystem]-=UBondOld[CurrentSystem];
    UAdsorbateUreyBradley[CurrentSystem]-=UUreyBradleyOld[CurrentSystem];
    UAdsorbateBend[CurrentSystem]-=UBendOld[CurrentSystem];
    UAdsorbateBendBend[CurrentSystem]-=UBendBendOld[CurrentSystem];
    UAdsorbateInversionBend[CurrentSystem]-=UInversionBendOld[CurrentSystem];
    UAdsorbateTorsion[CurrentSystem]-=UTorsionOld[CurrentSystem];
    UAdsorbateImproperTorsion[CurrentSystem]-=UImproperTorsionOld[CurrentSystem];
    UAdsorbateBondBond[CurrentSystem]-=UBondBondOld[CurrentSystem];
    UAdsorbateBondBend[CurrentSystem]-=UBondBendOld[CurrentSystem];
    UAdsorbateBondTorsion[CurrentSystem]-=UBondTorsionOld[CurrentSystem];
    UAdsorbateBendTorsion[CurrentSystem]-=UBendTorsionOld[CurrentSystem];
    UAdsorbateIntraVDW[CurrentSystem]-=UIntraVDWOld[CurrentSystem];
    UAdsorbateIntraChargeCharge[CurrentSystem]-=UIntraChargeChargeOld[CurrentSystem];
    UAdsorbateIntraChargeBondDipole[CurrentSystem]-=UIntraChargeBondDipoleOld[CurrentSystem];
    UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]-=UIntraBondDipoleBondDipoleOld[CurrentSystem];

    UAdsorbateAdsorbate[CurrentSystem]-=UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]-=UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]-=UCationVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]-=UCationVDWOld[CurrentSystem];
    UHostAdsorbate[CurrentSystem]-=UHostVDWOld[CurrentSystem];
    UHostAdsorbateVDW[CurrentSystem]-=UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    UTailCorrection[CurrentSystem]-=UTailOld;

    if(ChargeMethod!=NONE)
    {
      UAdsorbateIntraChargeCharge[CurrentSystem]-=UIntraChargeChargeOld[CurrentSystem];
      UAdsorbateIntraChargeBondDipole[CurrentSystem]-=UIntraChargeBondDipoleOld[CurrentSystem];
      UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]-=UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]-=UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]-=UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]-=UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]
                                                +UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
                                                +UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                                -UAdsorbateChargeChargeOld[CurrentSystem]
                                                -UAdsorbateChargeBondDipoleOld[CurrentSystem]
                                                -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]
                                         +UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
                                         +UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                         -UAdsorbateChargeChargeOld[CurrentSystem]
                                         -UAdsorbateChargeBondDipoleOld[CurrentSystem]
                                         -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      UHostAdsorbateChargeChargeReal[CurrentSystem]-=UHostChargeChargeOld[CurrentSystem];
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]-=UHostChargeBondDipoleOld[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]-=UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]
                                           +UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
                                           +UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                           -UHostChargeChargeOld[CurrentSystem]
                                           -UHostChargeBondDipoleOld[CurrentSystem]
                                           -UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]
                                    +UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
                                    +UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                    -UHostChargeChargeOld[CurrentSystem]
                                    -UHostChargeBondDipoleOld[CurrentSystem]
                                    -UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]-=UCationChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=UCationChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]-=UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]
                                             +UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]
                                             +UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                             -UCationChargeChargeOld[CurrentSystem]
                                             -UCationChargeBondDipoleOld[CurrentSystem]
                                             -UCationBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]
                                      +UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]
                                      +UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                      -UCationChargeChargeOld[CurrentSystem]
                                      -UCationChargeBondDipoleOld[CurrentSystem]
                                      -UCationBondDipoleBondDipoleOld[CurrentSystem];

      NetChargeAdsorbates[CurrentSystem]+=NetChargeDeltaOld;
      NetChargeSystem[CurrentSystem]+=NetChargeDeltaOld;

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(B);
    }

    RemoveAdsorbateMolecule();

    DeltaU=-UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
           -UImproperTorsionOld[CurrentSystem]-UBondBondNew[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
           -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
           -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
           -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
           -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
           -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
           -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
            UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
            UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
            UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
            UDeltaPolarization-UTailOld;

    UTotal[CurrentSystem]+=DeltaU;
  }
  CurrentSystem=CurrentSystemStored;
  return 0;
}

/*********************************************************************************************************
 * Name       | GibbsParticleTransferCationMove                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to do a Gibbs particle transfer from one box to the other.              *
 * Parameters | -                                                                                        *
 * Note       | The ideal gas Rosenbluth weights cancel.                                                 *
 *********************************************************************************************************/

int GibbsParticleTransferCationMove(void)
{
  int i,A,B;
  REAL RosenbluthNew,UTailNew;
  REAL RosenbluthOld,UTailOld;
  REAL NetChargeDeltaNew,NetChargeDeltaOld;
  int CurrentSystemStored;
  REAL DeltaU;

  CurrentSystemStored=CurrentSystem;

  // choose the box random
  if(RandomNumber()<0.5)
  {
    A=0;
    B=1;
  }
  else
  {
    A=1;
    B=0;
  }

  GibbsSwapAttempts[CurrentComponent][A]+=1.0;

  // return if there are no particles of the chosen component in box B
  if(Components[CurrentComponent].NumberOfMolecules[B]<=Components[CurrentComponent].CFMoleculePresent[B]?1:0) return 0;

  // first grow a particle of the chosen type in box A
  CurrentSystem=A;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be added here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  CurrentAdsorbateMolecule=NumberOfAdsorbateMolecules[CurrentSystem];
  CurrentCationMolecule=NumberOfCationMolecules[CurrentSystem];
  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNew=GrowMolecule(CBMC_INSERTION);
  if (OVERLAP) return 0;

  UTailNew=TailMolecularEnergyDifferenceAdd();
  RosenbluthNew*=exp(-Beta[CurrentSystem]*UTailNew);

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CalculateEwaldFourierCation(TRUE,FALSE,NumberOfCationMolecules[CurrentSystem],A);
    NetChargeDeltaNew=NetChargeCationDelta;
    RosenbluthNew*=exp(-Beta[CurrentSystem]*(
        UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
        UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
        UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(TRUE,-1,CurrentCationMolecule);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    RosenbluthNew*=exp(-Beta[CurrentSystem]*UDeltaPolarization);
  }

  // second retrace a particle of the chosen type in box B
  CurrentSystem=B;

  // choose a random molecule of this component
  CurrentAdsorbateMolecule=-1;
  CurrentCationMolecule=SelectRandomMoleculeOfTypeExcludingFractionalMolecule(CurrentComponent);

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be removed here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // calculate the Old Rosenbluth factor
  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position;
  RosenbluthOld=RetraceMolecule(CBMC_DELETION);

  UTailOld=TailMolecularEnergyDifferenceRemove();
  RosenbluthOld*=exp(-Beta[CurrentSystem]*UTailOld);

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CalculateEwaldFourierCation(FALSE,TRUE,CurrentCationMolecule,B);
    NetChargeDeltaOld=NetChargeCationDelta;
    RosenbluthOld*=exp(Beta[CurrentSystem]*(
       UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
       UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
       UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  if(ComputePolarization)
  {
    ComputeNewPolarizationEnergy(FALSE,-1,CurrentCationMolecule);
    if(OVERLAP) return 0;

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
                       (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
                       UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
                       UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
                       (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
                       UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
    RosenbluthOld*=exp(Beta[CurrentSystem]*UDeltaPolarization);
  }

  // acceptence rule
  if(RandomNumber()<((RosenbluthNew*(REAL)(Components[CurrentComponent].NumberOfMolecules[B])*Volume[A])/
                    (RosenbluthOld*((REAL)(1.0+Components[CurrentComponent].NumberOfMolecules[A]))*Volume[B])))
  {
    GibbsSwapAccepted[CurrentComponent][A]+=1.0;

    CurrentSystem=A;

    UCationBond[CurrentSystem]+=UBondNew[CurrentSystem];
    UCationUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem];
    UCationBend[CurrentSystem]+=UBendNew[CurrentSystem];
    UCationBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem];
    UCationInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem];
    UCationTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem];
    UCationImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem];
    UCationBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem];
    UCationBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem];
    UCationBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem];
    UCationBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem];
    UCationIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem];
    UCationIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem];
    UCationIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem];
    UCationIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem];

    UCationCation[CurrentSystem]+=UCationVDWNew[CurrentSystem];
    UCationCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem];
    UAdsorbateCation[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem];
    UHostCation[CurrentSystem]+=UHostVDWNew[CurrentSystem];
    UHostCationVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    UTailCorrection[CurrentSystem]+=UTailNew;

    if(ChargeMethod!=NONE)
    {
      UCationIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem];
      UCationIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem];
      UCationIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem];

      UCationCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]+
                                           UCationChargeBondDipoleNew[CurrentSystem]+
                                           UCationBondDipoleBondDipoleNew[CurrentSystem]+
                                           UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                           UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                           UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]+
                                    UCationChargeBondDipoleNew[CurrentSystem]+
                                    UCationBondDipoleBondDipoleNew[CurrentSystem]+
                                    UCationCationChargeChargeFourierDelta[CurrentSystem]+
                                    UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                    UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                         UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                         UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                         UHostChargeChargeNew[CurrentSystem]+
                                         UHostChargeBondDipoleNew[CurrentSystem]+
                                         UHostBondDipoleBondDipoleNew[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]+
                                  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                  UHostChargeChargeNew[CurrentSystem]+
                                  UHostChargeBondDipoleNew[CurrentSystem]+
                                  UHostBondDipoleBondDipoleNew[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                              UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateChargeChargeNew[CurrentSystem]+
                                              UAdsorbateChargeBondDipoleNew[CurrentSystem]+
                                              UAdsorbateBondDipoleBondDipoleNew[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
                                       UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
                                       UAdsorbateChargeChargeNew[CurrentSystem]+
                                       UAdsorbateChargeBondDipoleNew[CurrentSystem]+
                                       UAdsorbateBondDipoleBondDipoleNew[CurrentSystem];

      NetChargeCations[CurrentSystem]+=NetChargeDeltaNew;
      NetChargeSystem[CurrentSystem]+=NetChargeDeltaNew;

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(A);
    }

    InsertCationMolecule();

    DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
           UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
           UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
           UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
           UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
           UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
           UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
           UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
           UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
           UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
           UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
           UDeltaPolarization+UTailNew;

    UTotal[CurrentSystem]+=DeltaU;

    CurrentSystem=B;

    UCationBond[CurrentSystem]-=UBondOld[CurrentSystem];
    UCationUreyBradley[CurrentSystem]-=UUreyBradleyOld[CurrentSystem];
    UCationBend[CurrentSystem]-=UBendOld[CurrentSystem];
    UCationBendBend[CurrentSystem]-=UBendBendOld[CurrentSystem];
    UCationInversionBend[CurrentSystem]-=UInversionBendOld[CurrentSystem];
    UCationTorsion[CurrentSystem]-=UTorsionOld[CurrentSystem];
    UCationImproperTorsion[CurrentSystem]-=UImproperTorsionOld[CurrentSystem];
    UCationBondBond[CurrentSystem]-=UBondBondOld[CurrentSystem];
    UCationBondBend[CurrentSystem]-=UBondBendOld[CurrentSystem];
    UCationBondTorsion[CurrentSystem]-=UBondTorsionOld[CurrentSystem];
    UCationBendTorsion[CurrentSystem]-=UBendTorsionOld[CurrentSystem];
    UCationIntraVDW[CurrentSystem]-=UIntraVDWOld[CurrentSystem];
    UCationIntraChargeCharge[CurrentSystem]-=UIntraChargeChargeOld[CurrentSystem];
    UCationIntraChargeBondDipole[CurrentSystem]-=UIntraChargeBondDipoleOld[CurrentSystem];
    UCationIntraBondDipoleBondDipole[CurrentSystem]-=UIntraBondDipoleBondDipoleOld[CurrentSystem];

    UCationCation[CurrentSystem]-=UCationVDWOld[CurrentSystem];
    UCationCationVDW[CurrentSystem]-=UCationVDWOld[CurrentSystem];
    UAdsorbateCation[CurrentSystem]-=UAdsorbateVDWOld[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]-=UAdsorbateVDWOld[CurrentSystem];
    UHostCation[CurrentSystem]-=UHostVDWOld[CurrentSystem];
    UHostCationVDW[CurrentSystem]-=UHostVDWOld[CurrentSystem];

    UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
    UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

    UTailCorrection[CurrentSystem]-=UTailOld;

    if(ChargeMethod!=NONE)
    {
      UCationIntraChargeCharge[CurrentSystem]-=UIntraChargeChargeOld[CurrentSystem];
      UCationIntraChargeBondDipole[CurrentSystem]-=UIntraChargeBondDipoleOld[CurrentSystem];
      UCationIntraBondDipoleBondDipole[CurrentSystem]-=UIntraBondDipoleBondDipoleOld[CurrentSystem];

      UCationCationChargeChargeReal[CurrentSystem]-=UCationChargeChargeOld[CurrentSystem];
      UCationCationChargeBondDipoleReal[CurrentSystem]-=UCationChargeBondDipoleOld[CurrentSystem];
      UCationCationBondDipoleBondDipoleReal[CurrentSystem]-=UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCationChargeChargeFourier[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem];
      UCationCationChargeBondDipoleFourier[CurrentSystem]+=UCationCationChargeBondDipoleFourierDelta[CurrentSystem];
      UCationCationBondDipoleBondDipoleFourier[CurrentSystem]+=UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UCationCationCoulomb[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]
                                          +UCationCationChargeBondDipoleFourierDelta[CurrentSystem]
                                          +UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                          -UCationChargeChargeOld[CurrentSystem]
                                          -UCationChargeBondDipoleOld[CurrentSystem]
                                          -UCationBondDipoleBondDipoleOld[CurrentSystem];
      UCationCation[CurrentSystem]+=UCationCationChargeChargeFourierDelta[CurrentSystem]
                                   +UCationCationChargeBondDipoleFourierDelta[CurrentSystem]
                                   +UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                   -UCationChargeChargeOld[CurrentSystem]
                                   -UCationChargeBondDipoleOld[CurrentSystem]
                                   -UCationBondDipoleBondDipoleOld[CurrentSystem];

      UHostCationChargeChargeReal[CurrentSystem]-=UHostChargeChargeOld[CurrentSystem];
      UHostCationChargeBondDipoleReal[CurrentSystem]-=UHostChargeBondDipoleOld[CurrentSystem];
      UHostCationBondDipoleBondDipoleReal[CurrentSystem]-=UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCationChargeChargeFourier[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem];
      UHostCationChargeBondDipoleFourier[CurrentSystem]+=UHostCationChargeBondDipoleFourierDelta[CurrentSystem];
      UHostCationBondDipoleBondDipoleFourier[CurrentSystem]+=UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UHostCationCoulomb[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]
                                        +UHostCationChargeBondDipoleFourierDelta[CurrentSystem]
                                        +UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                        -UHostChargeChargeOld[CurrentSystem]
                                        -UHostChargeBondDipoleOld[CurrentSystem]
                                        -UHostBondDipoleBondDipoleOld[CurrentSystem];
      UHostCation[CurrentSystem]+=UHostCationChargeChargeFourierDelta[CurrentSystem]
                                 +UHostCationChargeBondDipoleFourierDelta[CurrentSystem]
                                 +UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                 -UHostChargeChargeOld[CurrentSystem]
                                 -UHostChargeBondDipoleOld[CurrentSystem]
                                 -UHostBondDipoleBondDipoleOld[CurrentSystem];

      UAdsorbateCationChargeChargeReal[CurrentSystem]-=UAdsorbateChargeChargeOld[CurrentSystem];
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=UAdsorbateChargeBondDipoleOld[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]-=UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]
                                             +UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]
                                             +UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                             -UAdsorbateChargeChargeOld[CurrentSystem]
                                             -UAdsorbateChargeBondDipoleOld[CurrentSystem]
                                             -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]
                                      +UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]
                                      +UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
                                      -UAdsorbateChargeChargeOld[CurrentSystem]
                                      -UAdsorbateChargeBondDipoleOld[CurrentSystem]
                                      -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

      NetChargeCations[CurrentSystem]+=NetChargeDeltaOld;
      NetChargeSystem[CurrentSystem]+=NetChargeDeltaOld;

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldCationMove(B);
    }

    RemoveCationMolecule();

    DeltaU=-UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
           -UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
           -UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]
           -UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
           -UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
           -UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
           -UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
           -UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
            UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
            UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
            UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
            UDeltaPolarization-UTailOld;

    UTotal[CurrentSystem]+=DeltaU;
  }
  CurrentSystem=CurrentSystemStored;
  return 0;
}

void GibbsParticleTransferMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    GibbsParticleTransferCationMove();
  else
    GibbsParticleTransferAdsorbateMove();
}

void PrintGibbsSwapStatistics(FILE *FilePtr)
{
  int i,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfGibbsSwapChangeMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the Gibbs swap move:\n");
    fprintf(FilePtr,"============================================\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      if(GibbsSwapAttempts[i][CurrentSystem]>0.0)
      {
        fprintf(FilePtr,"Component [%s] transfered to system [%d] total tried: %lf accepted: %lf (%lf [%%])\n",
          Components[i].Name,
          CurrentSystem,
          (double)GibbsSwapAttempts[i][CurrentSystem],
          (double)GibbsSwapAccepted[i][CurrentSystem],
          (double)(GibbsSwapAttempts[i][CurrentSystem]>(REAL)0.0?
            100.0*GibbsSwapAccepted[i][CurrentSystem]/GibbsSwapAttempts[i][CurrentSystem]:(REAL)0.0));
      }
    }
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"Gibbs swap move was OFF for all components\n\n");
}


/*********************************************************************************************************
 * Name       | GibbsVolumeMove                                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to change the volumes of the boxes in the Gibbs-ensemble.               *
 * Parameters | -                                                                                        *
 * Note       | The Gibbs ensemble Monte Carlo technique allows direct simulation of phase equilibria in *
 *            | fluids. Gibbs ensemnle simulations are performed in two separate microscopic regions,    *
 *            | each with periodic boundary conditions.                                                  *
 *            | An n-component system at constant temperature T, total volume V, and total number        *
 *            | of particles N is divided into two regions, with volumes V_I and V_II=V-V_I              *
 *            | The "GibbsVolume" move changes the volume of both boxes, while keeping the total volume  *
 *            | constant. The volumes are changed by a random change in ln(V_I/V_II).                    *
 *********************************************************************************************************/

int GibbsVolumeMove(void)
{
  int i,j;
  int A,B;
  int NumberOfMoleculesA,NumberOfMoleculesB;
  REAL det,scale,expdv;
  REAL_MATRIX3x3 StoredBox[2];
  REAL StoredVolume[2];
  REAL TotalVolume;
  VECTOR com,d;
  REAL StoredUAdsorbateAdsorbate[2],StoredUAdsorbateAdsorbateVDW[2],StoredUAdsorbateAdsorbateChargeChargeReal[2];
  REAL StoredUAdsorbateAdsorbateChargeBondDipoleReal[2],StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal[2];
  REAL StoredUAdsorbateAdsorbateChargeChargeFourier[2],StoredUAdsorbateAdsorbateCoulomb[2];
  REAL StoredUAdsorbateAdsorbateChargeBondDipoleFourier[2],StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier[2];
  REAL StoredUAdsorbateCation[2],StoredUAdsorbateCationVDW[2],StoredUAdsorbateCationChargeChargeReal[2];
  REAL StoredUAdsorbateCationChargeBondDipoleReal[2],StoredUAdsorbateCationBondDipoleBondDipoleReal[2];
  REAL StoredUAdsorbateCationChargeChargeFourier[2],StoredUAdsorbateCationCoulomb[2];
  REAL StoredUAdsorbateCationChargeBondDipoleFourier[2],StoredUAdsorbateCationBondDipoleBondDipoleFourier[2];
  REAL StoredUCationCation[2],StoredUCationCationVDW[2],StoredUCationCationChargeChargeReal[2];
  REAL StoredUCationCationChargeBondDipoleReal[2],StoredUCationCationBondDipoleBondDipoleReal[2];
  REAL StoredUCationCationChargeChargeFourier[2],StoredUCationCationCoulomb[2];
  REAL StoredUCationCationChargeBondDipoleFourier[2],StoredUCationCationBondDipoleBondDipoleFourier[2];
  REAL StoredUTotal[2],StoredUIon[2],StoredUTail[2];
  int CurrentSystemStored;

  CurrentSystemStored=CurrentSystem;

  A=0;
  B=1;

  // store state of system A
  CurrentSystem=A;
  GibbsVolumeChangeAttempts[CurrentSystem]+=1.0;

  StoredBox[CurrentSystem]=Box[CurrentSystem];
  StoredVolume[CurrentSystem]=Volume[CurrentSystem];

  StoredUTotal[CurrentSystem]=UTotal[CurrentSystem];
  StoredUIon[CurrentSystem]=UIon[CurrentSystem];
  StoredUTail[CurrentSystem]=UTailCorrection[CurrentSystem];

  StoredUAdsorbateAdsorbate[CurrentSystem]=UAdsorbateAdsorbate[CurrentSystem];
  StoredUAdsorbateAdsorbateVDW[CurrentSystem]=UAdsorbateAdsorbateVDW[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=UAdsorbateAdsorbateChargeChargeReal[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem];
  StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateCoulomb[CurrentSystem]=UAdsorbateAdsorbateCoulomb[CurrentSystem];

  StoredUAdsorbateCation[CurrentSystem]=UAdsorbateCation[CurrentSystem];
  StoredUAdsorbateCationVDW[CurrentSystem]=UAdsorbateCationVDW[CurrentSystem];
  StoredUAdsorbateCationChargeChargeReal[CurrentSystem]=UAdsorbateCationChargeChargeReal[CurrentSystem];
  StoredUAdsorbateCationChargeBondDipoleReal[CurrentSystem]=UAdsorbateCationChargeBondDipoleReal[CurrentSystem];
  StoredUAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem];
  StoredUAdsorbateCationChargeChargeFourier[CurrentSystem]=UAdsorbateCationChargeChargeFourier[CurrentSystem];
  StoredUAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=UAdsorbateCationChargeBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateCationCoulomb[CurrentSystem]=UAdsorbateCationCoulomb[CurrentSystem];

  StoredUCationCation[CurrentSystem]=UCationCation[CurrentSystem];
  StoredUCationCationVDW[CurrentSystem]=UCationCationVDW[CurrentSystem];
  StoredUCationCationChargeChargeReal[CurrentSystem]=UCationCationChargeChargeReal[CurrentSystem];
  StoredUCationCationChargeBondDipoleReal[CurrentSystem]=UCationCationChargeBondDipoleReal[CurrentSystem];
  StoredUCationCationBondDipoleBondDipoleReal[CurrentSystem]=UCationCationBondDipoleBondDipoleReal[CurrentSystem];
  StoredUCationCationChargeChargeFourier[CurrentSystem]=UCationCationChargeChargeFourier[CurrentSystem];
  StoredUCationCationChargeBondDipoleFourier[CurrentSystem]=UCationCationChargeBondDipoleFourier[CurrentSystem];
  StoredUCationCationBondDipoleBondDipoleFourier[CurrentSystem]=UCationCationBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUCationCationCoulomb[CurrentSystem]=UCationCationCoulomb[CurrentSystem];

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++) 
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition=Adsorbates[CurrentSystem][i].Atoms[j].Position;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++) 
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      Cations[CurrentSystem][i].Atoms[j].ReferencePosition=Cations[CurrentSystem][i].Atoms[j].Position;

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    SaveCurrentEwaldStructureFactors(CurrentSystem,CurrentSystem);
    SaveCurrentKVectors(CurrentSystem,CurrentSystem);
    SetupKVectors();
  }

  // store state of system B
  CurrentSystem=B;
  GibbsVolumeChangeAttempts[CurrentSystem]+=1.0;

  StoredBox[CurrentSystem]=Box[CurrentSystem];
  StoredUIon[CurrentSystem]=UIon[CurrentSystem];
  StoredVolume[CurrentSystem]=Volume[CurrentSystem];

  StoredUTotal[CurrentSystem]=UTotal[CurrentSystem];
  StoredUTail[CurrentSystem]=UTailCorrection[CurrentSystem];

  StoredUAdsorbateAdsorbate[CurrentSystem]=UAdsorbateAdsorbate[CurrentSystem];
  StoredUAdsorbateAdsorbateVDW[CurrentSystem]=UAdsorbateAdsorbateVDW[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=UAdsorbateAdsorbateChargeChargeReal[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem];
  StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateAdsorbateCoulomb[CurrentSystem]=UAdsorbateAdsorbateCoulomb[CurrentSystem];

  StoredUAdsorbateCation[CurrentSystem]=UAdsorbateCation[CurrentSystem];
  StoredUAdsorbateCationVDW[CurrentSystem]=UAdsorbateCationVDW[CurrentSystem];
  StoredUAdsorbateCationChargeChargeReal[CurrentSystem]=UAdsorbateCationChargeChargeReal[CurrentSystem];
  StoredUAdsorbateCationChargeBondDipoleReal[CurrentSystem]=UAdsorbateCationChargeBondDipoleReal[CurrentSystem];
  StoredUAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem];
  StoredUAdsorbateCationChargeChargeFourier[CurrentSystem]=UAdsorbateCationChargeChargeFourier[CurrentSystem];
  StoredUAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=UAdsorbateCationChargeBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUAdsorbateCationCoulomb[CurrentSystem]=UAdsorbateCationCoulomb[CurrentSystem];

  StoredUCationCation[CurrentSystem]=UCationCation[CurrentSystem];
  StoredUCationCationVDW[CurrentSystem]=UCationCationVDW[CurrentSystem];
  StoredUCationCationChargeChargeReal[CurrentSystem]=UCationCationChargeChargeReal[CurrentSystem];
  StoredUCationCationChargeBondDipoleReal[CurrentSystem]=UCationCationChargeBondDipoleReal[CurrentSystem];
  StoredUCationCationBondDipoleBondDipoleReal[CurrentSystem]=UCationCationBondDipoleBondDipoleReal[CurrentSystem];
  StoredUCationCationChargeChargeFourier[CurrentSystem]=UCationCationChargeChargeFourier[CurrentSystem];
  StoredUCationCationChargeBondDipoleFourier[CurrentSystem]=UCationCationChargeBondDipoleFourier[CurrentSystem];
  StoredUCationCationBondDipoleBondDipoleFourier[CurrentSystem]=UCationCationBondDipoleBondDipoleFourier[CurrentSystem];
  StoredUCationCationCoulomb[CurrentSystem]=UCationCationCoulomb[CurrentSystem];

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition=Adsorbates[CurrentSystem][i].Atoms[j].Position;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      Cations[CurrentSystem][i].Atoms[j].ReferencePosition=Cations[CurrentSystem][i].Atoms[j].Position;

  // determine New box-volumes leaving the total volume constant
  TotalVolume=Volume[A]+Volume[B];
  expdv=exp(log(Volume[A]/Volume[B])+MaximumGibbsVolumeChange[A]*(2.0*RandomNumber()-1.0));
  Volume[A]=expdv*TotalVolume/(1.0+expdv);
  Volume[B]=TotalVolume-Volume[A];

  // compute total energy for state A
  CurrentSystem=A;
  scale=pow(Volume[CurrentSystem]/StoredVolume[CurrentSystem],(REAL)1.0/3.0);
  Box[CurrentSystem].ax*=scale;  Box[CurrentSystem].bx*=scale;  Box[CurrentSystem].cx*=scale;
  Box[CurrentSystem].ay*=scale;  Box[CurrentSystem].by*=scale;  Box[CurrentSystem].cy*=scale;
  Box[CurrentSystem].az*=scale;  Box[CurrentSystem].bz*=scale;  Box[CurrentSystem].cz*=scale;
  Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

  if(MIN3(BoxProperties[CurrentSystem].cx,BoxProperties[CurrentSystem].cy,
          BoxProperties[CurrentSystem].cz)<2.0*CutOffVDW)
  {
    for(i=0;i<NumberOfSystems;i++)
    {
      fprintf(OutputFilePtr[i],"\n");
      fprintf(OutputFilePtr[i],"Gibbs ERROR 1: (System %d) Cutoff smaller than half of one of the perpendicular boxlengths !!!\n",CurrentSystem);
      fprintf(OutputFilePtr[i],"       Cutoff: %lf perpendicular boxlengths: %lf %lf %lf\n",(double)CutOffVDW,
               (double)BoxProperties[CurrentSystem].cx,(double)BoxProperties[CurrentSystem].cy,(double)BoxProperties[CurrentSystem].cz);
      fflush(OutputFilePtr[i]);
    }
    exit(0);
  }

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    SaveCurrentEwaldStructureFactors(CurrentSystem,CurrentSystem);
    SaveCurrentKVectors(CurrentSystem,CurrentSystem);
    SetupKVectors();
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    UpdateGroupCenterOfMassAdsorbate(i);
    com=GetAdsorbateCenterOfMass(i);

    d.x=com.x*(scale-1.0);
    d.y=com.y*(scale-1.0);
    d.z=com.z*(scale-1.0);

    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Position.x+=d.x;
      Adsorbates[CurrentSystem][i].Atoms[j].Position.y+=d.y;
      Adsorbates[CurrentSystem][i].Atoms[j].Position.z+=d.z;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    UpdateGroupCenterOfMassCation(i);
    com=GetCationCenterOfMass(i);

    d.x=com.x*(scale-1.0);
    d.y=com.y*(scale-1.0);
    d.z=com.z*(scale-1.0);

    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Position.x+=d.x;
      Cations[CurrentSystem][i].Atoms[j].Position.y+=d.y;
      Cations[CurrentSystem][i].Atoms[j].Position.z+=d.z;
    }
  }


  CalculateEnergy();

  // compute total energy for state B
  CurrentSystem=B;
  scale=pow(Volume[CurrentSystem]/StoredVolume[CurrentSystem],(REAL)1.0/3.0);
  Box[CurrentSystem].ax*=scale;  Box[CurrentSystem].bx*=scale;  Box[CurrentSystem].cx*=scale;
  Box[CurrentSystem].ay*=scale;  Box[CurrentSystem].by*=scale;  Box[CurrentSystem].cy*=scale;
  Box[CurrentSystem].az*=scale;  Box[CurrentSystem].bz*=scale;  Box[CurrentSystem].cz*=scale;
  Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

  if(MIN3(BoxProperties[CurrentSystem].cx,BoxProperties[CurrentSystem].cy,
          BoxProperties[CurrentSystem].cz)<2.0*CutOffVDW)
  {
    for(i=0;i<NumberOfSystems;i++)
    {
      fprintf(OutputFilePtr[i],"\n");
      fprintf(OutputFilePtr[i],"Gibbs ERROR 2: (System %d) Cutoff smaller than half of one of the perpendicular boxlengths !!!\n",CurrentSystem);
      fprintf(OutputFilePtr[i],"       Cutoff: %lf perpendicular boxlengths: %lf %lf %lf\n",(double)CutOffVDW,
              (double)BoxProperties[CurrentSystem].cx,(double)BoxProperties[CurrentSystem].cy,(double)BoxProperties[CurrentSystem].cz);
      fflush(OutputFilePtr[i]);
    }
    exit(0);
  }

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    SaveCurrentEwaldStructureFactors(CurrentSystem,CurrentSystem);
    SaveCurrentKVectors(CurrentSystem,CurrentSystem);
    SetupKVectors();
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    UpdateGroupCenterOfMassAdsorbate(i);
    com=GetAdsorbateCenterOfMass(i);

    d.x=com.x*(scale-1.0);
    d.y=com.y*(scale-1.0);
    d.z=com.z*(scale-1.0);

    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Position.x+=d.x;
      Adsorbates[CurrentSystem][i].Atoms[j].Position.y+=d.y;
      Adsorbates[CurrentSystem][i].Atoms[j].Position.z+=d.z;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    UpdateGroupCenterOfMassCation(i);
    com=GetCationCenterOfMass(i);

    d.x=com.x*(scale-1.0);
    d.y=com.y*(scale-1.0);
    d.z=com.z*(scale-1.0);

    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Position.x+=d.x;
      Cations[CurrentSystem][i].Atoms[j].Position.y+=d.y;
      Cations[CurrentSystem][i].Atoms[j].Position.z+=d.z;
    }
  }

  CalculateEnergy();

  NumberOfMoleculesA=NumberOfAdsorbateMolecules[A]+NumberOfCationMolecules[A];
  NumberOfMoleculesB=NumberOfAdsorbateMolecules[B]+NumberOfCationMolecules[B];

  if(RandomNumber()<exp(-Beta[A]*((UTotal[A]-StoredUTotal[A])+(UTotal[B]-StoredUTotal[B]))+
        ((NumberOfMoleculesA+1.0)*log(Volume[A]/StoredVolume[A]))+
        ((NumberOfMoleculesB+1.0)*log(Volume[B]/StoredVolume[B]))))
  {
    GibbsVolumeChangeAccepted[A]+=1.0;
    GibbsVolumeChangeAccepted[B]+=1.0;

    CurrentSystem=A;
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassCation(i);

    CurrentSystem=B;
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassCation(i);
  }
  else
  {
    // restore the state of system A
    CurrentSystem=A;
    UAdsorbateAdsorbate[CurrentSystem]=StoredUAdsorbateAdsorbate[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]=StoredUAdsorbateAdsorbateVDW[CurrentSystem];
    UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=StoredUAdsorbateAdsorbateChargeChargeReal[CurrentSystem];
    UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=StoredUAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem];
    UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
    UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=StoredUAdsorbateAdsorbateChargeChargeFourier[CurrentSystem];
    UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=StoredUAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem];
    UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
    UAdsorbateAdsorbateCoulomb[CurrentSystem]=StoredUAdsorbateAdsorbateCoulomb[CurrentSystem];

    UAdsorbateCation[CurrentSystem]=StoredUAdsorbateCation[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]=StoredUAdsorbateCationVDW[CurrentSystem];
    UAdsorbateCationChargeChargeReal[CurrentSystem]=StoredUAdsorbateCationChargeChargeReal[CurrentSystem];
    UAdsorbateCationChargeBondDipoleReal[CurrentSystem]=StoredUAdsorbateCationChargeBondDipoleReal[CurrentSystem];
    UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem];
    UAdsorbateCationChargeChargeFourier[CurrentSystem]=StoredUAdsorbateCationChargeChargeFourier[CurrentSystem];
    UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=StoredUAdsorbateCationChargeBondDipoleFourier[CurrentSystem];
    UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem];
    UAdsorbateCationCoulomb[CurrentSystem]=StoredUAdsorbateCationCoulomb[CurrentSystem];

    UCationCation[CurrentSystem]=StoredUCationCation[CurrentSystem];
    UCationCationVDW[CurrentSystem]=StoredUCationCationVDW[CurrentSystem];
    UCationCationChargeChargeReal[CurrentSystem]=StoredUCationCationChargeChargeReal[CurrentSystem];
    UCationCationChargeBondDipoleReal[CurrentSystem]=StoredUCationCationChargeBondDipoleReal[CurrentSystem];
    UCationCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUCationCationBondDipoleBondDipoleReal[CurrentSystem];
    UCationCationChargeChargeFourier[CurrentSystem]=StoredUCationCationChargeChargeFourier[CurrentSystem];
    UCationCationChargeBondDipoleFourier[CurrentSystem]=StoredUCationCationChargeBondDipoleFourier[CurrentSystem];
    UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUCationCationBondDipoleBondDipoleFourier[CurrentSystem];
    UCationCationCoulomb[CurrentSystem]=StoredUCationCationCoulomb[CurrentSystem];

    UTotal[CurrentSystem]=StoredUTotal[CurrentSystem];
    UIon[CurrentSystem]=StoredUIon[CurrentSystem];
    UTailCorrection[CurrentSystem]=StoredUTail[CurrentSystem];

    Box[CurrentSystem]=StoredBox[CurrentSystem];
    Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
    CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        Adsorbates[CurrentSystem][i].Atoms[j].Position=Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition;

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        Cations[CurrentSystem][i].Atoms[j].Position=Cations[CurrentSystem][i].Atoms[j].ReferencePosition;

    CalculateAnisotropicSites();

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      RetrieveStoredEwaldStructureFactors(CurrentSystem,CurrentSystem);
      RetrieveStoredKVectors(CurrentSystem,CurrentSystem);
    }

    // restore the state of system B
    CurrentSystem=B;
    UAdsorbateAdsorbate[CurrentSystem]=StoredUAdsorbateAdsorbate[CurrentSystem];
    UAdsorbateAdsorbateVDW[CurrentSystem]=StoredUAdsorbateAdsorbateVDW[CurrentSystem];
    UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=StoredUAdsorbateAdsorbateChargeChargeReal[CurrentSystem];
    UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=StoredUAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem];
    UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=StoredUAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem];
    UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=StoredUAdsorbateAdsorbateChargeChargeFourier[CurrentSystem];
    UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=StoredUAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem];
    UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=StoredUAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem];
    UAdsorbateAdsorbateCoulomb[CurrentSystem]=StoredUAdsorbateAdsorbateCoulomb[CurrentSystem];

    UAdsorbateCation[CurrentSystem]=StoredUAdsorbateCation[CurrentSystem];
    UAdsorbateCationVDW[CurrentSystem]=StoredUAdsorbateCationVDW[CurrentSystem];
    UAdsorbateCationChargeChargeReal[CurrentSystem]=StoredUAdsorbateCationChargeChargeReal[CurrentSystem];
    UAdsorbateCationChargeBondDipoleReal[CurrentSystem]=StoredUAdsorbateCationChargeBondDipoleReal[CurrentSystem];
    UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem];
    UAdsorbateCationChargeChargeFourier[CurrentSystem]=StoredUAdsorbateCationChargeChargeFourier[CurrentSystem];
    UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=StoredUAdsorbateCationChargeBondDipoleFourier[CurrentSystem];
    UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem];
    UAdsorbateCationCoulomb[CurrentSystem]=StoredUAdsorbateCationCoulomb[CurrentSystem];

    UCationCation[CurrentSystem]=StoredUCationCation[CurrentSystem];
    UCationCationVDW[CurrentSystem]=StoredUCationCationVDW[CurrentSystem];
    UCationCationChargeChargeReal[CurrentSystem]=StoredUCationCationChargeChargeReal[CurrentSystem];
    UCationCationChargeBondDipoleReal[CurrentSystem]=StoredUCationCationChargeBondDipoleReal[CurrentSystem];
    UCationCationBondDipoleBondDipoleReal[CurrentSystem]=StoredUCationCationBondDipoleBondDipoleReal[CurrentSystem];
    UCationCationChargeChargeFourier[CurrentSystem]=StoredUCationCationChargeChargeFourier[CurrentSystem];
    UCationCationChargeBondDipoleFourier[CurrentSystem]=StoredUCationCationChargeBondDipoleFourier[CurrentSystem];
    UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=StoredUCationCationBondDipoleBondDipoleFourier[CurrentSystem];
    UCationCationCoulomb[CurrentSystem]=StoredUCationCationCoulomb[CurrentSystem];

    UTotal[CurrentSystem]=StoredUTotal[CurrentSystem];
    UIon[CurrentSystem]=StoredUIon[CurrentSystem];
    UTailCorrection[CurrentSystem]=StoredUTail[CurrentSystem];

    Box[CurrentSystem]=StoredBox[CurrentSystem];
    Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
    CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        Adsorbates[CurrentSystem][i].Atoms[j].Position=Adsorbates[CurrentSystem][i].Atoms[j].ReferencePosition;

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        Cations[CurrentSystem][i].Atoms[j].Position=Cations[CurrentSystem][i].Atoms[j].ReferencePosition;

    CalculateAnisotropicSites();

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      RetrieveStoredEwaldStructureFactors(CurrentSystem,CurrentSystem);
      RetrieveStoredKVectors(CurrentSystem,CurrentSystem);
    }
  }
  CurrentSystem=CurrentSystemStored;
  return 0;
}

void PrintGibbsVolumeChangeStatistics(FILE *FilePtr)
{
  if(ProbabilityGibbsVolumeChangeMove)
  {
    fprintf(FilePtr,"Performance of the Gibbs volume change move:\n");
    fprintf(FilePtr,"============================================\n");
    if(GibbsVolumeChangeAttempts[CurrentSystem]>0.0)
    {
       fprintf(FilePtr,"total tried: %lf accepted: %lf (%lf [%%])\n",
          (double)GibbsVolumeChangeAttempts[CurrentSystem],
          (double)GibbsVolumeChangeAccepted[CurrentSystem],
          (double)(GibbsVolumeChangeAttempts[CurrentSystem]>(REAL)0.0?
            100.0*GibbsVolumeChangeAccepted[CurrentSystem]/GibbsVolumeChangeAttempts[CurrentSystem]:(REAL)0.0));
       fprintf(FilePtr,"\tmaximum volume change %lf\n",(double)MaximumGibbsVolumeChange[CurrentSystem]);
    }
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"Gibbs volume change move was OFF\n\n");
}

/*********************************************************************************************************
 * Name       | GibbsIdentityChangeAdsorbateMove                                                         *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Monte Carlo move to reinsert a particle with another identity in the Gibbs-ensemble.     *
 * Parameters | -                                                                                        *
 * Note       |                                                                                          *
 *********************************************************************************************************/

// Gibbs identity-change Monte-Carlo move
// ==========================================================
// A very efficient trial move for the simulation of mixtures is the use of identity
// switch. In this move, a particle is regrown in the same box at the same
// position but with a different identity. Suppose that we have a mixture of n
// components and that we use the following algorithm to select particles
// for an identity switch in the Gibbs ensemble
//   1. Out of the n components, two components i!=j are selected at random
//   2. At random, it is selected to switch the identity of component i
//      in box A or in box B, and the identity of the component j in the other box.
//   3. In each box, a particle is selected at random which matches the desired
//      identity.
//  The acceptance rule was originally derived by Panagiotopoulos et al. Martin et
//  al. have used this trial moves using the CBMC algorithm to compute
//  vapor-liquid equilibria of binary mixtures of linear alkanes.

// A.Z. Panagiotopoulos,
// "Exact calculations of fluid-phase equilibria by Monte Carlo simulations in a new ensemble",
// International Journal of Thermophysiics, volume 10, issue 2, pages 447-457, 1989

// M.G. Martin and J.I. Siepmann
// "Predicting Mulitcomponent Phase Equilibria and Free Energies of Transfer for Alkens by Molecular Simulation",
// Journal of American Chemical Society, volume 119, issue 38, pages 8921-8924, 1997

// A.Z. Panagiotopoulos, N. Quirke, M. Stapleton, and D.J. Tildesley
// "Phase-Equilibria by Simulation in the Gibbs Ensemble - Alternative Derivation, Generalization and
//  Application to Mixture and Membrane Equilibria",
// Molecular Physics, volume 63, issue 4, pages 527-545, 1988

// A.Z. Panagiotopoulos
// "Direct Determination of Fluid-Phase Equilbriua by Simulation in the Gibbs Ensemble - A Review",
// Molecular Simulation, volume 9, issue 1, pages 1-23, 1992

int GibbsIdentityChangeAdsorbateMove(void)
{
  int i,BoxI,BoxII,d,count,nr_atoms,type;
  int ComponentA,ComponentB;
  int StartingBeadOld,StartingBeadNew;
  int AdsorbateMoleculeA,AdsorbateMoleculeB;
  REAL RosenbluthNewA,RosenbluthNewB,RosenbluthOldA,RosenbluthOldB;
  REAL UTailCorrectionDifferenceA,UTailCorrectionDifferenceB,PreFactorA,PreFactorB;
  REAL RosenbluthIdealNew,RosenbluthIdealOld;
  int CurrentSystemStored;
  GROUP *group_temp;
  ATOM *atom_temp;

  CurrentSystemStored=CurrentSystem;

  // 'A'-component was already randomly chosen (as CurrentComponent), 
  // 'B'-component randomly chosen from the list defined by 'A'
  ComponentA=CurrentComponent;
  d=(int)(RandomNumber()*(REAL)Components[ComponentA].NumberOfGibbsIdentityChanges);
  ComponentB=Components[ComponentA].GibbsIdentityChanges[d];

  if(Components[ComponentB].ExtraFrameworkMolecule) return 0;

  // garantue detail balance, i.e. choose the box random
  if(RandomNumber()<0.5)
  {
    BoxI=0;
    BoxII=1;
  }
  else
  {
    BoxI=1;
    BoxII=0;
  }

  GibbsIdentityChangeAttempts[BoxI][ComponentA][ComponentB]+=1.0;
  GibbsIdentityChangeAttempts[BoxII][ComponentB][ComponentA]+=1.0;

  // return when particles of the chosen conditions are not present
  if((Components[ComponentA].NumberOfMolecules[BoxI]==0)||
     (Components[ComponentB].NumberOfMolecules[BoxII]==0)) return 0;

  // choose a random molecule of the 'Old'-component in box I
  d=(int)(RandomNumber()*(REAL)Components[ComponentA].NumberOfMolecules[BoxI]);
  count=-1;
  CurrentCationMolecule=-1;
  AdsorbateMoleculeA=-1;
  do   // search for d-th molecule of the right type
    if(Adsorbates[BoxI][++AdsorbateMoleculeA].Type==ComponentA) count++;
  while(d!=count);

  // choose a random molecule of the 'New'-component in box II
  d=(int)(RandomNumber()*(REAL)Components[ComponentB].NumberOfMolecules[BoxII]);
  count=-1;
  CurrentCationMolecule=-1;
  AdsorbateMoleculeB=-1;
  do   // search for d-th molecule of the right type
    if(Adsorbates[BoxII][++AdsorbateMoleculeB].Type==ComponentB) count++;
  while(d!=count);

  // in box A an attempt is made to change a molecule of type A into B
  CurrentSystem=BoxI;
  CurrentAdsorbateMolecule=AdsorbateMoleculeA;
  CurrentComponent=ComponentB;

  StartingBeadNew=Components[ComponentB].StartingBead;
  StartingBeadOld=Components[ComponentA].StartingBead;

  // the starting position of the 'Old'-component is taken as the starting position for the 'New'-component
  NewPosition[CurrentSystem][StartingBeadNew]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBeadOld].Position;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be added here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // grow the 'New'-component
  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNewA=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  // in box B an attempt is made to change a molecule of type B into A
  CurrentSystem=BoxII;
  CurrentAdsorbateMolecule=AdsorbateMoleculeB;
  CurrentComponent=ComponentA;

  StartingBeadNew=Components[ComponentA].StartingBead;
  StartingBeadOld=Components[ComponentB].StartingBead;

  // the starting position of the 'Old'-component is taken as the starting position for the 'New'-component
  NewPosition[CurrentSystem][StartingBeadNew]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[StartingBeadOld].Position;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be added here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // grow the 'New'-component
  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNewB=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  // in box A an attempt is made to change a molecule of type A into B
  CurrentSystem=BoxI;
  CurrentAdsorbateMolecule=AdsorbateMoleculeA;
  CurrentComponent=ComponentA;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be removed here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // retrace the 'Old'-component
  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
  RosenbluthOldA=RetraceMolecule(CBMC_PARTIAL_INSERTION);

  // compute the tail-correction and Ewald-correction difference in Box I
  UTailCorrectionDifferenceA=TailMolecularEnergyDifferenceAddRemove(ComponentB,ComponentA);
  PreFactorA=exp(-Beta[CurrentSystem]*UTailCorrectionDifferenceA);

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CurrentComponent=ComponentB;
    CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,BoxI);
    PreFactorA*=exp(-Beta[CurrentSystem]*(
       UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
       UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
       UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  // in box B an attempt is made to change a molecule of type B into A
  CurrentSystem=BoxII;
  CurrentAdsorbateMolecule=AdsorbateMoleculeB;
  CurrentComponent=ComponentB;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be removed here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // retrace the 'Old'-component
  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;
  RosenbluthOldB=RetraceMolecule(CBMC_PARTIAL_INSERTION);

  // compute the tail-correction and Ewald-correction difference in Box II
  UTailCorrectionDifferenceB=TailMolecularEnergyDifferenceAddRemove(ComponentA,ComponentB);
  PreFactorB=exp(-Beta[CurrentSystem]*UTailCorrectionDifferenceB);

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CurrentComponent=ComponentA;
    CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,BoxII);
    PreFactorB*=exp(-Beta[CurrentSystem]*(
       UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
       UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
       UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  GibbsIdentityChangeAccepted[BoxI][ComponentA][ComponentB][0]+=1.0;
  GibbsIdentityChangeAccepted[BoxII][ComponentB][ComponentA][0]+=1.0;

  RosenbluthIdealNew=Components[ComponentA].IdealGasRosenbluthWeight[0];
  RosenbluthIdealOld=Components[ComponentB].IdealGasRosenbluthWeight[0];

  RosenbluthIdealNew=RosenbluthIdealOld=1.0;


  // acceptance-rule
  if(RandomNumber()<((PreFactorA*PreFactorB*RosenbluthNewA*RosenbluthNewB/(RosenbluthIdealNew*RosenbluthIdealOld*RosenbluthOldA*RosenbluthOldB))*
                    ((REAL)Components[ComponentA].NumberOfMolecules[BoxI]*(REAL)Components[ComponentB].NumberOfMolecules[BoxII])/
                    (((REAL)Components[ComponentA].NumberOfMolecules[BoxII]+1.0)*((REAL)Components[ComponentB].NumberOfMolecules[BoxI]+1.0))))
  {
    // register an succesfull growth/retrace after acceptance
    GibbsIdentityChangeAccepted[BoxI][ComponentA][ComponentB][1]+=1.0;
    GibbsIdentityChangeAccepted[BoxII][ComponentB][ComponentA][1]+=1.0;

    // register the changes in the energies for system A
    UAdsorbateBond[BoxI]+=UBondNew[BoxI]-UBondOld[BoxI];
    UAdsorbateUreyBradley[BoxI]+=UUreyBradleyNew[BoxI]-UUreyBradleyOld[BoxI];
    UAdsorbateBend[BoxI]+=UBendNew[BoxI]-UBendOld[BoxI];
    UAdsorbateBendBend[BoxI]+=UBendBendNew[BoxI]-UBendBendOld[BoxI];
    UAdsorbateInversionBend[BoxI]+=UInversionBendNew[BoxI]-UInversionBendOld[BoxI];
    UAdsorbateTorsion[BoxI]+=UTorsionNew[BoxI]-UTorsionOld[BoxI];
    UAdsorbateImproperTorsion[BoxI]+=UImproperTorsionNew[BoxI]-UImproperTorsionOld[BoxI];
    UAdsorbateBondBond[BoxI]+=UBondBondNew[BoxI]-UBondBondOld[BoxI];
    UAdsorbateBondBend[BoxI]+=UBondBendNew[BoxI]-UBondBendOld[BoxI];
    UAdsorbateBondTorsion[BoxI]+=UBondTorsionNew[BoxI]-UBondTorsionOld[BoxI];
    UAdsorbateBendTorsion[BoxI]+=UBendTorsionNew[BoxI]-UBendTorsionOld[BoxI];
    UAdsorbateIntraVDW[BoxI]+=UIntraVDWNew[BoxI]-UIntraVDWOld[BoxI];
    UAdsorbateIntraChargeCharge[BoxI]+=UIntraChargeChargeNew[BoxI]-UIntraChargeChargeOld[BoxI];
    UAdsorbateIntraChargeBondDipole[BoxI]+=UIntraChargeBondDipoleNew[BoxI]-UIntraChargeBondDipoleOld[BoxI];
    UAdsorbateIntraBondDipoleBondDipole[BoxI]+=UIntraBondDipoleBondDipoleNew[BoxI]-UIntraBondDipoleBondDipoleOld[BoxI];

    UAdsorbateAdsorbate[BoxI]+=UAdsorbateVDWNew[BoxI]-UAdsorbateVDWOld[BoxI];
    UAdsorbateAdsorbateVDW[BoxI]+=UAdsorbateVDWNew[BoxI]-UAdsorbateVDWOld[BoxI];
    UAdsorbateCation[BoxI]+=UCationVDWNew[BoxI]-UCationVDWOld[BoxI];
    UAdsorbateCationVDW[BoxI]+=UCationVDWNew[BoxI]-UCationVDWOld[BoxI];
    UHostAdsorbate[BoxI]+=UHostVDWNew[BoxI]-UHostVDWOld[BoxI];
    UHostAdsorbateVDW[BoxI]+=UHostVDWNew[BoxI]-UHostVDWOld[BoxI];

    UTailCorrection[BoxI]+=UTailCorrectionDifferenceA;

    if(ChargeMethod!=NONE)
    {
      UAdsorbateAdsorbateChargeChargeReal[BoxI]+=UAdsorbateChargeChargeNew[BoxI]-UAdsorbateChargeChargeOld[BoxI];
      UAdsorbateAdsorbateChargeBondDipoleReal[BoxI]+=UAdsorbateChargeBondDipoleNew[BoxI]-UAdsorbateChargeBondDipoleOld[BoxI];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[BoxI]+=UAdsorbateBondDipoleBondDipoleNew[BoxI]-UAdsorbateBondDipoleBondDipoleOld[BoxI];
      UAdsorbateAdsorbateChargeChargeFourier[BoxI]+=UAdsorbateAdsorbateChargeChargeFourierDelta[BoxI];
      UAdsorbateAdsorbateChargeBondDipoleFourier[BoxI]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[BoxI];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[BoxI]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[BoxI];
      UAdsorbateAdsorbateCoulomb[BoxI]+=UAdsorbateAdsorbateChargeChargeFourierDelta[BoxI]+UAdsorbateChargeChargeNew[BoxI]-UAdsorbateChargeChargeOld[BoxI]+
                                        UAdsorbateAdsorbateChargeBondDipoleFourierDelta[BoxI]+UAdsorbateChargeBondDipoleNew[BoxI]-UAdsorbateChargeBondDipoleOld[BoxI]+
                                        UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[BoxI]+UAdsorbateBondDipoleBondDipoleNew[BoxI]-UAdsorbateBondDipoleBondDipoleOld[BoxI];
      UAdsorbateAdsorbate[BoxI]+=UAdsorbateAdsorbateChargeChargeFourierDelta[BoxI]+UAdsorbateChargeChargeNew[BoxI]-UAdsorbateChargeChargeOld[BoxI]+
                                 UAdsorbateAdsorbateChargeBondDipoleFourierDelta[BoxI]+UAdsorbateChargeBondDipoleNew[BoxI]-UAdsorbateChargeBondDipoleOld[BoxI]+
                                 UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[BoxI]+UAdsorbateBondDipoleBondDipoleNew[BoxI]-UAdsorbateBondDipoleBondDipoleOld[BoxI];

      UHostAdsorbateChargeChargeReal[BoxI]+=UHostChargeChargeNew[BoxI]-UHostChargeChargeOld[BoxI];
      UHostAdsorbateChargeBondDipoleReal[BoxI]+=UHostChargeBondDipoleNew[BoxI]-UHostChargeBondDipoleOld[BoxI];
      UHostAdsorbateBondDipoleBondDipoleReal[BoxI]+=UHostBondDipoleBondDipoleNew[BoxI]-UHostBondDipoleBondDipoleOld[BoxI];
      UHostAdsorbateChargeChargeFourier[BoxI]+=UHostAdsorbateChargeChargeFourierDelta[BoxI];
      UHostAdsorbateChargeBondDipoleFourier[BoxI]+=UHostAdsorbateChargeBondDipoleFourierDelta[BoxI];
      UHostAdsorbateBondDipoleBondDipoleFourier[BoxI]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[BoxI];
      UHostAdsorbateCoulomb[BoxI]+=UHostAdsorbateChargeChargeFourierDelta[BoxI]+UHostChargeChargeNew[BoxI]-UHostChargeChargeOld[BoxI]+
                                   UHostAdsorbateChargeBondDipoleFourierDelta[BoxI]+UHostChargeBondDipoleNew[BoxI]-UHostChargeBondDipoleOld[BoxI]+
                                   UHostAdsorbateBondDipoleBondDipoleFourierDelta[BoxI]+UHostBondDipoleBondDipoleNew[BoxI]-UHostBondDipoleBondDipoleOld[BoxI];
      UHostAdsorbate[BoxI]+=UHostAdsorbateChargeChargeFourierDelta[BoxI]+UHostChargeChargeNew[BoxI]-UHostChargeChargeOld[BoxI]+
                            UHostAdsorbateChargeBondDipoleFourierDelta[BoxI]+UHostChargeBondDipoleNew[BoxI]-UHostChargeBondDipoleOld[BoxI]+
                            UHostAdsorbateBondDipoleBondDipoleFourierDelta[BoxI]+UHostBondDipoleBondDipoleNew[BoxI]-UHostBondDipoleBondDipoleOld[BoxI];

      UAdsorbateCationChargeChargeReal[BoxI]+=UCationChargeChargeNew[BoxI]-UCationChargeChargeOld[BoxI];
      UAdsorbateCationChargeBondDipoleReal[BoxI]+=UCationChargeBondDipoleNew[BoxI]-UCationChargeBondDipoleOld[BoxI];
      UAdsorbateCationBondDipoleBondDipoleReal[BoxI]+=UCationBondDipoleBondDipoleNew[BoxI]-UCationBondDipoleBondDipoleOld[BoxI];
      UAdsorbateCationChargeChargeFourier[BoxI]+=UAdsorbateCationChargeChargeFourierDelta[BoxI];
      UAdsorbateCationChargeBondDipoleFourier[BoxI]+=UAdsorbateCationChargeBondDipoleFourierDelta[BoxI];
      UAdsorbateCationBondDipoleBondDipoleFourier[BoxI]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxI];
      UAdsorbateCationCoulomb[BoxI]+=UAdsorbateCationChargeChargeFourierDelta[BoxI]+UCationChargeChargeNew[BoxI]-UCationChargeChargeOld[BoxI]+
                                     UAdsorbateCationChargeBondDipoleFourierDelta[BoxI]+UCationChargeBondDipoleNew[BoxI]-UCationChargeBondDipoleOld[BoxI]+
                                     UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxI]+UCationBondDipoleBondDipoleNew[BoxI]-UCationBondDipoleBondDipoleOld[BoxI];
      UAdsorbateCation[BoxI]+=UAdsorbateCationChargeChargeFourierDelta[BoxI]+UCationChargeChargeNew[BoxI]-UCationChargeChargeOld[BoxI]+
                              UAdsorbateCationChargeBondDipoleFourierDelta[BoxI]+UCationChargeBondDipoleNew[BoxI]-UCationChargeBondDipoleOld[BoxI]+
                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxI]+UCationBondDipoleBondDipoleNew[BoxI]-UCationBondDipoleBondDipoleOld[BoxI];

      NetChargeAdsorbates[CurrentSystem]+=NetChargeAdsorbateDelta;

      // register the changes in the stored structure factors
      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      {
        CurrentSystem=BoxI;
        AcceptEwaldAdsorbateMove(BoxI);
      }
    }

    // register the changes in the energies for system B
    UAdsorbateBond[BoxII]+=UBondNew[BoxII]-UBondOld[BoxII];
    UAdsorbateUreyBradley[BoxII]+=UUreyBradleyNew[BoxII]-UUreyBradleyOld[BoxII];
    UAdsorbateBend[BoxII]+=UBendNew[BoxII]-UBendOld[BoxII];
    UAdsorbateBendBend[BoxII]+=UBendBendNew[BoxII]-UBendBendOld[BoxII];
    UAdsorbateInversionBend[BoxII]+=UInversionBendNew[BoxII]-UInversionBendOld[BoxII];
    UAdsorbateTorsion[BoxII]+=UTorsionNew[BoxII]-UTorsionOld[BoxII];
    UAdsorbateImproperTorsion[BoxII]+=UImproperTorsionNew[BoxII]-UImproperTorsionOld[BoxII];
    UAdsorbateBondBond[BoxII]+=UBondBondNew[BoxII]-UBondBondOld[BoxII];
    UAdsorbateBondBend[BoxII]+=UBondBendNew[BoxII]-UBondBendOld[BoxII];
    UAdsorbateBondTorsion[BoxII]+=UBondTorsionNew[BoxII]-UBondTorsionOld[BoxII];
    UAdsorbateBendTorsion[BoxII]+=UBendTorsionNew[BoxII]-UBendTorsionOld[BoxII];
    UAdsorbateIntraVDW[BoxII]+=UIntraVDWNew[BoxII]-UIntraVDWOld[BoxII];
    UAdsorbateIntraChargeCharge[BoxII]+=UIntraChargeChargeNew[BoxII]-UIntraChargeChargeOld[BoxII];
    UAdsorbateIntraChargeBondDipole[BoxII]+=UIntraChargeBondDipoleNew[BoxII]-UIntraChargeBondDipoleOld[BoxII];
    UAdsorbateIntraBondDipoleBondDipole[BoxII]+=UIntraBondDipoleBondDipoleNew[BoxII]-UIntraBondDipoleBondDipoleOld[BoxII];

    UAdsorbateAdsorbate[BoxII]+=UAdsorbateVDWNew[BoxII]-UAdsorbateVDWOld[BoxII];
    UAdsorbateAdsorbateVDW[BoxII]+=UAdsorbateVDWNew[BoxII]-UAdsorbateVDWOld[BoxII];
    UAdsorbateCation[BoxII]+=UCationVDWNew[BoxII]-UCationVDWOld[BoxII];
    UAdsorbateCationVDW[BoxII]+=UCationVDWNew[BoxII]-UCationVDWOld[BoxII];
    UHostAdsorbate[BoxII]+=UHostVDWNew[BoxII]-UHostVDWOld[BoxII];
    UHostAdsorbateVDW[BoxII]+=UHostVDWNew[BoxII]-UHostVDWOld[BoxII];

    UTailCorrection[BoxII]+=UTailCorrectionDifferenceB;

   if(ChargeMethod!=NONE)
    {
      UAdsorbateAdsorbateChargeChargeReal[BoxII]+=UAdsorbateChargeChargeNew[BoxII]-UAdsorbateChargeChargeOld[BoxII];
      UAdsorbateAdsorbateChargeBondDipoleReal[BoxII]+=UAdsorbateChargeBondDipoleNew[BoxII]-UAdsorbateChargeBondDipoleOld[BoxII];
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[BoxII]+=UAdsorbateBondDipoleBondDipoleNew[BoxII]-UAdsorbateBondDipoleBondDipoleOld[BoxII];
      UAdsorbateAdsorbateChargeChargeFourier[BoxII]+=UAdsorbateAdsorbateChargeChargeFourierDelta[BoxII];
      UAdsorbateAdsorbateChargeBondDipoleFourier[BoxII]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[BoxII];
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[BoxII]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[BoxII];
      UAdsorbateAdsorbateCoulomb[BoxII]+=UAdsorbateAdsorbateChargeChargeFourierDelta[BoxII]+UAdsorbateChargeChargeNew[BoxII]-UAdsorbateChargeChargeOld[BoxII]+
                                         UAdsorbateAdsorbateChargeBondDipoleFourierDelta[BoxII]+UAdsorbateChargeBondDipoleNew[BoxII]-UAdsorbateChargeBondDipoleOld[BoxII]+
                                         UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[BoxII]+UAdsorbateBondDipoleBondDipoleNew[BoxII]-UAdsorbateBondDipoleBondDipoleOld[BoxII];
      UAdsorbateAdsorbate[BoxII]+=UAdsorbateAdsorbateChargeChargeFourierDelta[BoxII]+UAdsorbateChargeChargeNew[BoxII]-UAdsorbateChargeChargeOld[BoxII]+
                                  UAdsorbateAdsorbateChargeBondDipoleFourierDelta[BoxII]+UAdsorbateChargeBondDipoleNew[BoxII]-UAdsorbateChargeBondDipoleOld[BoxII]+
                                  UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[BoxII]+UAdsorbateBondDipoleBondDipoleNew[BoxII]-UAdsorbateBondDipoleBondDipoleOld[BoxII];

      UHostAdsorbateChargeChargeReal[BoxII]+=UHostChargeChargeNew[BoxII]-UHostChargeChargeOld[BoxII];
      UHostAdsorbateChargeBondDipoleReal[BoxII]+=UHostChargeBondDipoleNew[BoxII]-UHostChargeBondDipoleOld[BoxII];
      UHostAdsorbateBondDipoleBondDipoleReal[BoxII]+=UHostBondDipoleBondDipoleNew[BoxII]-UHostBondDipoleBondDipoleOld[BoxII];
      UHostAdsorbateChargeChargeFourier[BoxII]+=UHostAdsorbateChargeChargeFourierDelta[BoxII];
      UHostAdsorbateChargeBondDipoleFourier[BoxII]+=UHostAdsorbateChargeBondDipoleFourierDelta[BoxII];
      UHostAdsorbateBondDipoleBondDipoleFourier[BoxII]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[BoxII];
      UHostAdsorbateCoulomb[BoxII]+=UHostAdsorbateChargeChargeFourierDelta[BoxII]+UHostChargeChargeNew[BoxII]-UHostChargeChargeOld[BoxII]+
                                    UHostAdsorbateChargeBondDipoleFourierDelta[BoxII]+UHostChargeBondDipoleNew[BoxII]-UHostChargeBondDipoleOld[BoxII]+
                                    UHostAdsorbateBondDipoleBondDipoleFourierDelta[BoxII]+UHostBondDipoleBondDipoleNew[BoxII]-UHostBondDipoleBondDipoleOld[BoxII];
      UHostAdsorbate[BoxII]+=UHostAdsorbateChargeChargeFourierDelta[BoxII]+UHostChargeChargeNew[BoxII]-UHostChargeChargeOld[BoxII]+
                             UHostAdsorbateChargeBondDipoleFourierDelta[BoxII]+UHostChargeBondDipoleNew[BoxII]-UHostChargeBondDipoleOld[BoxII]+
                             UHostAdsorbateBondDipoleBondDipoleFourierDelta[BoxII]+UHostBondDipoleBondDipoleNew[BoxII]-UHostBondDipoleBondDipoleOld[BoxII];

      UAdsorbateCationChargeChargeReal[BoxII]+=UCationChargeChargeNew[BoxII]-UCationChargeChargeOld[BoxII];
      UAdsorbateCationChargeBondDipoleReal[BoxII]+=UCationChargeBondDipoleNew[BoxII]-UCationChargeBondDipoleOld[BoxII];
      UAdsorbateCationBondDipoleBondDipoleReal[BoxII]+=UCationBondDipoleBondDipoleNew[BoxII]-UCationBondDipoleBondDipoleOld[BoxII];
      UAdsorbateCationChargeChargeFourier[BoxII]+=UAdsorbateCationChargeChargeFourierDelta[BoxII];
      UAdsorbateCationChargeBondDipoleFourier[BoxII]+=UAdsorbateCationChargeBondDipoleFourierDelta[BoxII];
      UAdsorbateCationBondDipoleBondDipoleFourier[BoxII]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxII];
      UAdsorbateCationCoulomb[BoxII]+=UAdsorbateCationChargeChargeFourierDelta[BoxII]+UCationChargeChargeNew[BoxII]-UCationChargeChargeOld[BoxII]+
                                      UAdsorbateCationChargeBondDipoleFourierDelta[BoxII]+UCationChargeBondDipoleNew[BoxII]-UCationChargeBondDipoleOld[BoxII]+
                                      UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxII]+UCationBondDipoleBondDipoleNew[BoxII]-UCationBondDipoleBondDipoleOld[BoxII];
      UAdsorbateCation[BoxII]+=UAdsorbateCationChargeChargeFourierDelta[BoxII]+UCationChargeChargeNew[BoxII]-UCationChargeChargeOld[BoxII]+
                                      UAdsorbateCationChargeBondDipoleFourierDelta[BoxII]+UCationChargeBondDipoleNew[BoxII]-UCationChargeBondDipoleOld[BoxII]+
                                      UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxII]+UCationBondDipoleBondDipoleNew[BoxII]-UCationBondDipoleBondDipoleOld[BoxII];

      // register the changes in the stored structure factors
      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      {
        CurrentSystem=BoxII;
        AcceptEwaldAdsorbateMove(BoxII);
      }
    }

    // swap the alocate memory for the atoms and the groups
    SWAP(Adsorbates[BoxI][AdsorbateMoleculeA].Atoms,Adsorbates[BoxII][AdsorbateMoleculeB].Atoms,atom_temp);
    SWAP(Adsorbates[BoxI][AdsorbateMoleculeA].Groups,Adsorbates[BoxII][AdsorbateMoleculeB].Groups,group_temp);

    // register the changes in the types for the 'Old'-component atoms for system A
    nr_atoms=Components[ComponentA].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      type=Components[ComponentA].Type[i];
      NumberOfPseudoAtomsType[BoxI][type]--;
    }

    // overwrite the 'Old'-component with the 'New'-component for system A
    nr_atoms=Components[ComponentB].NumberOfAtoms;
    Adsorbates[BoxI][AdsorbateMoleculeA].NumberOfAtoms=nr_atoms;
    Adsorbates[BoxI][AdsorbateMoleculeA].Type=ComponentB;
    for(i=0;i<nr_atoms;i++)
    {
      Adsorbates[BoxI][AdsorbateMoleculeA].Atoms[i].Position=TrialPosition[BoxI][i];
      type=Components[ComponentB].Type[i];
      Adsorbates[BoxI][AdsorbateMoleculeA].Atoms[i].Type=type;
      Adsorbates[BoxI][AdsorbateMoleculeA].Atoms[i].Charge=Components[ComponentB].Charge[i];
      NumberOfPseudoAtomsType[BoxI][type]++;
    }

    // register the changes in the types for the 'Old'-component atoms for system B
    nr_atoms=Components[ComponentB].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      type=Components[ComponentB].Type[i];
      NumberOfPseudoAtomsType[BoxII][type]--;
    }

    // overwrite the 'New'-component with the 'Old'-component for system A
    nr_atoms=Components[ComponentA].NumberOfAtoms;
    Adsorbates[BoxII][AdsorbateMoleculeB].NumberOfAtoms=nr_atoms;
    Adsorbates[BoxII][AdsorbateMoleculeB].Type=ComponentA;
    for(i=0;i<nr_atoms;i++)
    {
      Adsorbates[BoxII][AdsorbateMoleculeB].Atoms[i].Position=TrialPosition[BoxII][i];
      type=Components[ComponentA].Type[i];
      Adsorbates[BoxII][AdsorbateMoleculeB].Atoms[i].Type=type;
      Adsorbates[BoxII][AdsorbateMoleculeB].Atoms[i].Charge=Components[ComponentA].Charge[i];
      NumberOfPseudoAtomsType[BoxII][type]++;
    }


    // register the changes in the amount of molecules of each type
    Components[ComponentA].NumberOfMolecules[BoxI]--;
    Components[ComponentB].NumberOfMolecules[BoxI]++;

    Components[ComponentA].NumberOfMolecules[BoxII]++;
    Components[ComponentB].NumberOfMolecules[BoxII]--;

    // update the center of mass of the molecule
    CurrentSystem=BoxI;
    UpdateGroupCenterOfMassAdsorbate(AdsorbateMoleculeA);

    CurrentSystem=BoxII;
    UpdateGroupCenterOfMassAdsorbate(AdsorbateMoleculeB);

    // recompute the total energy
    UTotal[BoxI]=
         // inter
         UHostHost[CurrentSystem]+UCationCation[CurrentSystem]+UAdsorbateAdsorbate[CurrentSystem]+
         UHostCation[CurrentSystem]+UHostAdsorbate[CurrentSystem]+UAdsorbateCation[CurrentSystem]+
         // intra framework
         UHostBond[CurrentSystem]+UHostUreyBradley[CurrentSystem]+UHostBend[CurrentSystem]+
         UHostInversionBend[CurrentSystem]+UHostTorsion[CurrentSystem]+UHostImproperTorsion[CurrentSystem]+
         UHostBondBond[CurrentSystem]+UHostBondBend[CurrentSystem]+UHostBendBend[CurrentSystem]+
         UHostBondTorsion[CurrentSystem]+UHostBendTorsion[CurrentSystem]+
         // intra adsorbates
         UAdsorbateBond[CurrentSystem]+UAdsorbateUreyBradley[CurrentSystem]+UAdsorbateBend[CurrentSystem]+
         UAdsorbateInversionBend[CurrentSystem]+UAdsorbateTorsion[CurrentSystem]+UAdsorbateImproperTorsion[CurrentSystem]+
         UAdsorbateBondBond[CurrentSystem]+UAdsorbateBondBend[CurrentSystem]+UAdsorbateBendBend[CurrentSystem]+
         UAdsorbateBondTorsion[CurrentSystem]+UAdsorbateBendTorsion[CurrentSystem]+
         UAdsorbateIntraVDW[CurrentSystem]+UAdsorbateIntraChargeCharge[CurrentSystem]+
         UAdsorbateIntraChargeBondDipole[CurrentSystem]+UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+
         // intra cations
         UCationBond[CurrentSystem]+UCationUreyBradley[CurrentSystem]+UCationBend[CurrentSystem]+
         UCationInversionBend[CurrentSystem]+UCationTorsion[CurrentSystem]+UCationImproperTorsion[CurrentSystem]+
         UCationBondBond[CurrentSystem]+UCationBondBend[CurrentSystem]+UCationBendBend[CurrentSystem]+
         UCationBondTorsion[CurrentSystem]+UCationBendTorsion[CurrentSystem]+
         UCationIntraVDW[CurrentSystem]+UCationIntraChargeCharge[CurrentSystem]+
         UCationIntraChargeBondDipole[CurrentSystem]+UCationIntraBondDipoleBondDipole[CurrentSystem]+
         // polarization
         UHostPolarization[CurrentSystem]+UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem]+
         // tailcorrection
         UTailCorrection[CurrentSystem];

    UTotal[BoxII]=
         // inter
         UHostHost[CurrentSystem]+UCationCation[CurrentSystem]+UAdsorbateAdsorbate[CurrentSystem]+
         UHostCation[CurrentSystem]+UHostAdsorbate[CurrentSystem]+UAdsorbateCation[CurrentSystem]+
         // intra framework
         UHostBond[CurrentSystem]+UHostUreyBradley[CurrentSystem]+UHostBend[CurrentSystem]+
         UHostInversionBend[CurrentSystem]+UHostTorsion[CurrentSystem]+UHostImproperTorsion[CurrentSystem]+
         UHostBondBond[CurrentSystem]+UHostBondBend[CurrentSystem]+UHostBendBend[CurrentSystem]+
         UHostBondTorsion[CurrentSystem]+UHostBendTorsion[CurrentSystem]+
         // intra adsorbates
         UAdsorbateBond[CurrentSystem]+UAdsorbateUreyBradley[CurrentSystem]+UAdsorbateBend[CurrentSystem]+
         UAdsorbateInversionBend[CurrentSystem]+UAdsorbateTorsion[CurrentSystem]+UAdsorbateImproperTorsion[CurrentSystem]+
         UAdsorbateBondBond[CurrentSystem]+UAdsorbateBondBend[CurrentSystem]+UAdsorbateBendBend[CurrentSystem]+
         UAdsorbateBondTorsion[CurrentSystem]+UAdsorbateBendTorsion[CurrentSystem]+
         UAdsorbateIntraVDW[CurrentSystem]+UAdsorbateIntraChargeCharge[CurrentSystem]+
         UAdsorbateIntraChargeBondDipole[CurrentSystem]+UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+
         // intra cations
         UCationBond[CurrentSystem]+UCationUreyBradley[CurrentSystem]+UCationBend[CurrentSystem]+
         UCationInversionBend[CurrentSystem]+UCationTorsion[CurrentSystem]+UCationImproperTorsion[CurrentSystem]+
         UCationBondBond[CurrentSystem]+UCationBondBend[CurrentSystem]+UCationBendBend[CurrentSystem]+
         UCationBondTorsion[CurrentSystem]+UCationBendTorsion[CurrentSystem]+
         UCationIntraVDW[CurrentSystem]+UCationIntraChargeCharge[CurrentSystem]+
         UCationIntraChargeBondDipole[CurrentSystem]+UCationIntraBondDipoleBondDipole[CurrentSystem]+
         // polarization
         UHostPolarization[CurrentSystem]+UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem]+
         // tailcorrection
         UTailCorrection[CurrentSystem];
  }
  CurrentSystem=CurrentSystemStored;

  return 0;
}

int GibbsIdentityChangeCationMove(void)
{
  int i,BoxI,BoxII,d,count,nr_atoms,type;
  int ComponentA,ComponentB;
  int StartingBeadOld,StartingBeadNew;
  int CationMoleculeA,CationMoleculeB;
  REAL RosenbluthNewA,RosenbluthNewB,RosenbluthOldA,RosenbluthOldB;
  REAL UTailCorrectionDifferenceA,UTailCorrectionDifferenceB,PreFactorA,PreFactorB;
  REAL RosenbluthIdealNew,RosenbluthIdealOld;
  int CurrentSystemStored;
  GROUP *group_temp;
  ATOM *atom_temp;

  CurrentSystemStored=CurrentSystem;

  // 'A'-component was already randomly chosen (as CurrentComponent), 
  // 'B'-component randomly chosen from the list defined by 'A'
  ComponentA=CurrentComponent;
  d=(int)(RandomNumber()*(REAL)Components[ComponentA].NumberOfGibbsIdentityChanges);
  ComponentB=Components[ComponentA].GibbsIdentityChanges[d];

  if(!Components[ComponentB].ExtraFrameworkMolecule) return 0;

  // garantue detail balance, i.e. choose the box random
  if(RandomNumber()<0.5)
  {
    BoxI=0;
    BoxII=1;
  }
  else
  {
    BoxI=1;
    BoxII=0;
  }

  GibbsIdentityChangeAttempts[BoxI][ComponentA][ComponentB]+=1.0;
  GibbsIdentityChangeAttempts[BoxII][ComponentB][ComponentA]+=1.0;

  // return when particles of the chosen conditions are not present
  if((Components[ComponentA].NumberOfMolecules[BoxI]==0)||
     (Components[ComponentB].NumberOfMolecules[BoxII]==0)) return 0;

  // choose a random molecule of the 'Old'-component in box I
  d=(int)(RandomNumber()*(REAL)Components[ComponentA].NumberOfMolecules[BoxI]);
  count=-1;
  CurrentAdsorbateMolecule=-1;
  CationMoleculeA=-1;
  do   // search for d-th molecule of the right type
    if(Cations[BoxI][++CationMoleculeA].Type==ComponentA) count++;
  while(d!=count);

  // choose a random molecule of the 'New'-component in box II
  d=(int)(RandomNumber()*(REAL)Components[ComponentB].NumberOfMolecules[BoxII]);
  count=-1;
  CurrentAdsorbateMolecule=-1;
  CationMoleculeB=-1;
  do   // search for d-th molecule of the right type
    if(Cations[BoxII][++CationMoleculeB].Type==ComponentB) count++;
  while(d!=count);

  // in box A an attempt is made to change a molecule of type A into B
  CurrentSystem=BoxI;
  CurrentCationMolecule=CationMoleculeA;
  CurrentComponent=ComponentB;

  StartingBeadNew=Components[ComponentB].StartingBead;
  StartingBeadOld=Components[ComponentA].StartingBead;

  // the starting position of the 'Old'-component is taken as the starting position for the 'New'-component
  NewPosition[CurrentSystem][StartingBeadNew]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBeadOld].Position;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be added here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // grow the 'New'-component
  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNewA=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  // in box B an attempt is made to change a molecule of type B into A
  CurrentSystem=BoxII;
  CurrentCationMolecule=CationMoleculeB;
  CurrentComponent=ComponentA;

  StartingBeadNew=Components[ComponentA].StartingBead;
  StartingBeadOld=Components[ComponentB].StartingBead;

  // the starting position of the 'Old'-component is taken as the starting position for the 'New'-component
  NewPosition[CurrentSystem][StartingBeadNew]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[StartingBeadOld].Position;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be added here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // grow the 'New'-component
  NumberOfBeadsAlreadyPlaced=0;
  RosenbluthNewB=GrowMolecule(CBMC_PARTIAL_INSERTION);
  if (OVERLAP) return 0;

  // in box A an attempt is made to change a molecule of type A into B
  CurrentSystem=BoxI;
  CurrentCationMolecule=CationMoleculeA;
  CurrentComponent=ComponentA;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be added here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // retrace the 'Old'-component
  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position;
  RosenbluthOldA=RetraceMolecule(CBMC_PARTIAL_INSERTION);

  // compute the tail-correction and Ewald-correction difference in Box I
  UTailCorrectionDifferenceA=TailMolecularEnergyDifferenceAddRemove(ComponentB,ComponentA);
  PreFactorA=exp(-Beta[CurrentSystem]*UTailCorrectionDifferenceA);
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CurrentComponent=ComponentB;
    CalculateEwaldFourierCation(TRUE,TRUE,CurrentCationMolecule,BoxI);
    PreFactorA*=exp(-Beta[CurrentSystem]*(
        UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
        UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
        UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  // in box B an attempt is made to change a molecule of type B into A
  CurrentSystem=BoxII;
  CurrentCationMolecule=CationMoleculeB;
  CurrentComponent=ComponentB;

  // set Continuous Fraction (CF) atomic scaling-factors to unity (only integer molecules can be added here)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScaling[i]=1.0;
    CFChargeScaling[i]=1.0;
  }

  // retrace the 'Old'-component
  NumberOfBeadsAlreadyPlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    OldPosition[i]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[i].Position;
  RosenbluthOldB=RetraceMolecule(CBMC_PARTIAL_INSERTION);

  // compute the tail-correction and Ewald-correction difference in Box II
  UTailCorrectionDifferenceB=TailMolecularEnergyDifferenceAddRemove(ComponentA,ComponentB);
  PreFactorB=exp(-Beta[CurrentSystem]*UTailCorrectionDifferenceB);
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    CurrentComponent=ComponentA;
    CalculateEwaldFourierCation(TRUE,TRUE,CurrentCationMolecule,BoxII);
    PreFactorB*=exp(-Beta[CurrentSystem]*(
        UHostCationChargeChargeFourierDelta[CurrentSystem]+UCationCationChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
        UHostCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationCationChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
        UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]));
  }

  GibbsIdentityChangeAccepted[BoxI][ComponentA][ComponentB][0]+=1.0;
  GibbsIdentityChangeAccepted[BoxII][ComponentB][ComponentA][0]+=1.0;

  RosenbluthIdealNew=Components[ComponentA].IdealGasRosenbluthWeight[0];
  RosenbluthIdealOld=Components[ComponentB].IdealGasRosenbluthWeight[0];

  RosenbluthIdealNew=RosenbluthIdealOld=1.0;


  // acceptance-rule
  if(RandomNumber()<((PreFactorA*PreFactorB*RosenbluthNewA*RosenbluthNewB/(RosenbluthIdealNew*RosenbluthIdealOld*RosenbluthOldA*RosenbluthOldB))*
                    ((REAL)Components[BoxI].NumberOfMolecules[BoxI]*(REAL)Components[BoxII].NumberOfMolecules[BoxII])/
                    (((REAL)Components[BoxI].NumberOfMolecules[BoxII]+1.0)*((REAL)Components[BoxII].NumberOfMolecules[BoxI]+1.0))))
  {
    // register an succesfull growth/retrace after acceptance
    GibbsIdentityChangeAccepted[BoxI][ComponentA][ComponentB][1]+=1.0;
    GibbsIdentityChangeAccepted[BoxII][ComponentB][ComponentA][1]+=1.0;

    // register the changes in the energies for system A
    UCationBond[BoxI]+=UBondNew[BoxI]-UBondOld[BoxI];
    UCationUreyBradley[BoxI]+=UUreyBradleyNew[BoxI]-UUreyBradleyOld[BoxI];
    UCationBend[BoxI]+=UBendNew[BoxI]-UBendOld[BoxI];
    UCationBendBend[BoxI]+=UBendBendNew[BoxI]-UBendBendOld[BoxI];
    UCationInversionBend[BoxI]+=UInversionBendNew[BoxI]-UInversionBendOld[BoxI];
    UCationTorsion[BoxI]+=UTorsionNew[BoxI]-UTorsionOld[BoxI];
    UCationImproperTorsion[BoxI]+=UImproperTorsionNew[BoxI]-UImproperTorsionOld[BoxI];
    UCationBondBond[BoxI]+=UBondBondNew[BoxI]-UBondBondOld[BoxI];
    UCationBondBend[BoxI]+=UBondBendNew[BoxI]-UBondBendOld[BoxI];
    UCationBondTorsion[BoxI]+=UBondTorsionNew[BoxI]-UBondTorsionOld[BoxI];
    UCationBendTorsion[BoxI]+=UBendTorsionNew[BoxI]-UBendTorsionOld[BoxI];
    UCationIntraVDW[BoxI]+=UIntraVDWNew[BoxI]-UIntraVDWOld[BoxI];
    UCationIntraChargeCharge[BoxI]+=UIntraChargeChargeNew[BoxI]-UIntraChargeChargeOld[BoxI];
    UCationIntraChargeBondDipole[BoxI]+=UIntraChargeBondDipoleNew[BoxI]-UIntraChargeBondDipoleOld[BoxI];
    UCationIntraBondDipoleBondDipole[BoxI]+=UIntraBondDipoleBondDipoleNew[BoxI]-UIntraBondDipoleBondDipoleOld[BoxI];

    UCationCation[BoxI]+=UCationVDWNew[BoxI]-UCationVDWOld[BoxI];
    UCationCationVDW[BoxI]+=UCationVDWNew[BoxI]-UCationVDWOld[BoxI];
    UAdsorbateCation[BoxI]+=UAdsorbateVDWNew[BoxI]-UAdsorbateVDWOld[BoxI];
    UAdsorbateCationVDW[BoxI]+=UAdsorbateVDWNew[BoxI]-UAdsorbateVDWOld[BoxI];
    UHostCation[BoxI]+=UHostVDWNew[BoxI]-UHostVDWOld[BoxI];
    UHostCationVDW[BoxI]+=UHostVDWNew[BoxI]-UHostVDWOld[BoxI];

    UTailCorrection[BoxI]+=UTailCorrectionDifferenceA;

    if(ChargeMethod!=NONE)
    {
      UCationCationChargeChargeReal[BoxI]+=UCationChargeChargeNew[BoxI]-UCationChargeChargeOld[BoxI];
      UCationCationChargeBondDipoleReal[BoxI]+=UCationChargeBondDipoleNew[BoxI]-UCationChargeBondDipoleOld[BoxI];
      UCationCationBondDipoleBondDipoleReal[BoxI]+=UCationBondDipoleBondDipoleNew[BoxI]-UCationBondDipoleBondDipoleOld[BoxI];
      UCationCationChargeChargeFourier[BoxI]+=UCationCationChargeChargeFourierDelta[BoxI];
      UCationCationChargeBondDipoleFourier[BoxI]+=UCationCationChargeBondDipoleFourierDelta[BoxI];
      UCationCationBondDipoleBondDipoleFourier[BoxI]+=UCationCationBondDipoleBondDipoleFourierDelta[BoxI];
      UCationCationCoulomb[BoxI]+=UCationCationChargeChargeFourierDelta[BoxI]+UCationChargeChargeNew[BoxI]-UCationChargeChargeOld[BoxI]+
                                  UCationCationChargeBondDipoleFourierDelta[BoxI]+UCationChargeBondDipoleNew[BoxI]-UCationChargeBondDipoleOld[BoxI]+
                                  UCationCationBondDipoleBondDipoleFourierDelta[BoxI]+UCationBondDipoleBondDipoleNew[BoxI]-UCationBondDipoleBondDipoleOld[BoxI];
      UCationCation[BoxI]+=UCationCationChargeChargeFourierDelta[BoxI]+UCationChargeChargeNew[BoxI]-UCationChargeChargeOld[BoxI]+
                           UCationCationChargeBondDipoleFourierDelta[BoxI]+UCationChargeBondDipoleNew[BoxI]-UCationChargeBondDipoleOld[BoxI]+
                           UCationCationBondDipoleBondDipoleFourierDelta[BoxI]+UCationBondDipoleBondDipoleNew[BoxI]-UCationBondDipoleBondDipoleOld[BoxI];

      UHostCationChargeChargeReal[BoxI]+=UHostChargeChargeNew[BoxI]-UHostChargeChargeOld[BoxI];
      UHostCationChargeBondDipoleReal[BoxI]+=UHostChargeBondDipoleNew[BoxI]-UHostChargeBondDipoleOld[BoxI];
      UHostCationBondDipoleBondDipoleReal[BoxI]+=UHostBondDipoleBondDipoleNew[BoxI]-UHostBondDipoleBondDipoleOld[BoxI];
      UHostCationChargeChargeFourier[BoxI]+=UHostCationChargeChargeFourierDelta[BoxI];
      UHostCationChargeBondDipoleFourier[BoxI]+=UHostCationChargeBondDipoleFourierDelta[BoxI];
      UHostCationBondDipoleBondDipoleFourier[BoxI]+=UHostCationBondDipoleBondDipoleFourierDelta[BoxI];
      UHostCationCoulomb[BoxI]+=UHostCationChargeChargeFourierDelta[BoxI]+UHostChargeChargeNew[BoxI]-UHostChargeChargeOld[BoxI]+
                                UHostCationChargeBondDipoleFourierDelta[BoxI]+UHostChargeBondDipoleNew[BoxI]-UHostChargeBondDipoleOld[BoxI]+
                                UHostCationBondDipoleBondDipoleFourierDelta[BoxI]+UHostBondDipoleBondDipoleNew[BoxI]-UHostBondDipoleBondDipoleOld[BoxI];
      UHostCation[BoxI]+=UHostCationChargeChargeFourierDelta[BoxI]+UHostChargeChargeNew[BoxI]-UHostChargeChargeOld[BoxI]+
                         UHostCationChargeBondDipoleFourierDelta[BoxI]+UHostChargeBondDipoleNew[BoxI]-UHostChargeBondDipoleOld[BoxI]+
                         UHostCationBondDipoleBondDipoleFourierDelta[BoxI]+UHostBondDipoleBondDipoleNew[BoxI]-UHostBondDipoleBondDipoleOld[BoxI];

      UAdsorbateCationChargeChargeReal[BoxI]+=UAdsorbateChargeChargeNew[BoxI]-UAdsorbateChargeChargeOld[BoxI];
      UAdsorbateCationChargeBondDipoleReal[BoxI]+=UAdsorbateChargeBondDipoleNew[BoxI]-UAdsorbateChargeBondDipoleOld[BoxI];
      UAdsorbateCationBondDipoleBondDipoleReal[BoxI]+=UAdsorbateBondDipoleBondDipoleNew[BoxI]-UAdsorbateBondDipoleBondDipoleOld[BoxI];
      UAdsorbateCationChargeChargeFourier[BoxI]+=UAdsorbateCationChargeChargeFourierDelta[BoxI];
      UAdsorbateCationChargeBondDipoleFourier[BoxI]+=UAdsorbateCationChargeBondDipoleFourierDelta[BoxI];
      UAdsorbateCationBondDipoleBondDipoleFourier[BoxI]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxI];
      UAdsorbateCationCoulomb[BoxI]+=UAdsorbateCationChargeChargeFourierDelta[BoxI]+UAdsorbateChargeChargeNew[BoxI]-UAdsorbateChargeChargeOld[BoxI]+
                                     UAdsorbateCationChargeBondDipoleFourierDelta[BoxI]+UAdsorbateChargeBondDipoleNew[BoxI]-UAdsorbateChargeBondDipoleOld[BoxI]+
                                     UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxI]+UAdsorbateBondDipoleBondDipoleNew[BoxI]-UAdsorbateBondDipoleBondDipoleOld[BoxI];
      UAdsorbateCation[BoxI]+=UAdsorbateCationChargeChargeFourierDelta[BoxI]+UAdsorbateChargeChargeNew[BoxI]-UAdsorbateChargeChargeOld[BoxI]+
                              UAdsorbateCationChargeBondDipoleFourierDelta[BoxI]+UAdsorbateChargeBondDipoleNew[BoxI]-UAdsorbateChargeBondDipoleOld[BoxI]+
                              UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxI]+UAdsorbateBondDipoleBondDipoleNew[BoxI]-UAdsorbateBondDipoleBondDipoleOld[BoxI];

      NetChargeCations[CurrentSystem]+=NetChargeCationDelta;

      // register the changes in the stored structure factors
      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      {
        CurrentSystem=BoxI;
        AcceptEwaldCationMove(BoxI);
      }
    }

    // register the changes in the energies for system B
    UCationBond[BoxII]+=UBondNew[BoxII]-UBondOld[BoxII];
    UCationUreyBradley[BoxII]+=UUreyBradleyNew[BoxII]-UUreyBradleyOld[BoxII];
    UCationBend[BoxII]+=UBendNew[BoxII]-UBendOld[BoxII];
    UCationBendBend[BoxII]+=UBendBendNew[BoxII]-UBendBendOld[BoxII];
    UCationInversionBend[BoxII]+=UInversionBendNew[BoxII]-UInversionBendOld[BoxII];
    UCationTorsion[BoxII]+=UTorsionNew[BoxII]-UTorsionOld[BoxII];
    UCationImproperTorsion[BoxII]+=UImproperTorsionNew[BoxII]-UImproperTorsionOld[BoxII];
    UCationBondBond[BoxII]+=UBondBondNew[BoxII]-UBondBondOld[BoxII];
    UCationBondBend[BoxII]+=UBondBendNew[BoxII]-UBondBendOld[BoxII];
    UCationBondTorsion[BoxII]+=UBondTorsionNew[BoxII]-UBondTorsionOld[BoxII];
    UCationBendTorsion[BoxII]+=UBendTorsionNew[BoxII]-UBendTorsionOld[BoxII];
    UCationIntraVDW[BoxII]+=UIntraVDWNew[BoxII]-UIntraVDWOld[BoxII];
    UCationIntraChargeCharge[BoxII]+=UIntraChargeChargeNew[BoxII]-UIntraChargeChargeOld[BoxII];
    UCationIntraChargeBondDipole[BoxII]+=UIntraChargeBondDipoleNew[BoxII]-UIntraChargeBondDipoleOld[BoxII];
    UCationIntraBondDipoleBondDipole[BoxII]+=UIntraBondDipoleBondDipoleNew[BoxII]-UIntraBondDipoleBondDipoleOld[BoxII];

    UCationCation[BoxII]+=UCationVDWNew[BoxII]-UCationVDWOld[BoxII];
    UCationCationVDW[BoxII]+=UCationVDWNew[BoxII]-UCationVDWOld[BoxII];
    UAdsorbateCation[BoxII]+=UAdsorbateVDWNew[BoxII]-UAdsorbateVDWOld[BoxII];
    UAdsorbateCationVDW[BoxII]+=UAdsorbateVDWNew[BoxII]-UAdsorbateVDWOld[BoxII];
    UHostCation[BoxII]+=UHostVDWNew[BoxII]-UHostVDWOld[BoxII];
    UHostCationVDW[BoxII]+=UHostVDWNew[BoxII]-UHostVDWOld[BoxII];

    UTailCorrection[BoxII]+=UTailCorrectionDifferenceB;

   if(ChargeMethod!=NONE)
    {
      UCationCationChargeChargeReal[BoxII]+=UCationChargeChargeNew[BoxII]-UCationChargeChargeOld[BoxII];
      UCationCationChargeBondDipoleReal[BoxII]+=UCationChargeBondDipoleNew[BoxII]-UCationChargeBondDipoleOld[BoxII];
      UCationCationBondDipoleBondDipoleReal[BoxII]+=UCationBondDipoleBondDipoleNew[BoxII]-UCationBondDipoleBondDipoleOld[BoxII];
      UCationCationChargeChargeFourier[BoxII]+=UCationCationChargeChargeFourierDelta[BoxII];
      UCationCationChargeBondDipoleFourier[BoxII]+=UCationCationChargeBondDipoleFourierDelta[BoxII];
      UCationCationBondDipoleBondDipoleFourier[BoxII]+=UCationCationBondDipoleBondDipoleFourierDelta[BoxII];
      UCationCationCoulomb[BoxII]+=UCationCationChargeChargeFourierDelta[BoxII]+UCationChargeChargeNew[BoxII]-UCationChargeChargeOld[BoxII]+
                                   UCationCationChargeBondDipoleFourierDelta[BoxII]+UCationChargeBondDipoleNew[BoxII]-UCationChargeBondDipoleOld[BoxII]+
                                   UCationCationBondDipoleBondDipoleFourierDelta[BoxII]+UCationBondDipoleBondDipoleNew[BoxII]-UCationBondDipoleBondDipoleOld[BoxII];
      UCationCation[BoxII]+=UCationCationChargeChargeFourierDelta[BoxII]+UCationChargeChargeNew[BoxII]-UCationChargeChargeOld[BoxII]+
                            UCationCationChargeBondDipoleFourierDelta[BoxII]+UCationChargeBondDipoleNew[BoxII]-UCationChargeBondDipoleOld[BoxII]+
                            UCationCationBondDipoleBondDipoleFourierDelta[BoxII]+UCationBondDipoleBondDipoleNew[BoxII]-UCationBondDipoleBondDipoleOld[BoxII];

      UHostCationChargeChargeReal[BoxII]+=UHostChargeChargeNew[BoxII]-UHostChargeChargeOld[BoxII];
      UHostCationChargeBondDipoleReal[BoxII]+=UHostChargeBondDipoleNew[BoxII]-UHostChargeBondDipoleOld[BoxII];
      UHostCationBondDipoleBondDipoleReal[BoxII]+=UHostBondDipoleBondDipoleNew[BoxII]-UHostBondDipoleBondDipoleOld[BoxII];
      UHostCationChargeChargeFourier[BoxII]+=UHostCationChargeChargeFourierDelta[BoxII];
      UHostCationChargeBondDipoleFourier[BoxII]+=UHostCationChargeBondDipoleFourierDelta[BoxII];
      UHostCationBondDipoleBondDipoleFourier[BoxII]+=UHostCationBondDipoleBondDipoleFourierDelta[BoxII];
      UHostCationCoulomb[BoxII]+=UHostCationChargeChargeFourierDelta[BoxII]+UHostChargeChargeNew[BoxII]-UHostChargeChargeOld[BoxII]+
                                 UHostCationChargeBondDipoleFourierDelta[BoxII]+UHostChargeBondDipoleNew[BoxII]-UHostChargeBondDipoleOld[BoxII]+
                                 UHostCationBondDipoleBondDipoleFourierDelta[BoxII]+UHostBondDipoleBondDipoleNew[BoxII]-UHostBondDipoleBondDipoleOld[BoxII];
      UHostCation[BoxII]+=UHostCationChargeChargeFourierDelta[BoxII]+UHostChargeChargeNew[BoxII]-UHostChargeChargeOld[BoxII]+
                          UHostCationChargeBondDipoleFourierDelta[BoxII]+UHostChargeBondDipoleNew[BoxII]-UHostChargeBondDipoleOld[BoxII]+
                          UHostCationBondDipoleBondDipoleFourierDelta[BoxII]+UHostBondDipoleBondDipoleNew[BoxII]-UHostBondDipoleBondDipoleOld[BoxII];

      UAdsorbateCationChargeChargeReal[BoxII]+=UAdsorbateChargeChargeNew[BoxII]-UAdsorbateChargeChargeOld[BoxII];
      UAdsorbateCationChargeBondDipoleReal[BoxII]+=UAdsorbateChargeBondDipoleNew[BoxII]-UAdsorbateChargeBondDipoleOld[BoxII];
      UAdsorbateCationBondDipoleBondDipoleReal[BoxII]+=UAdsorbateBondDipoleBondDipoleNew[BoxII]-UAdsorbateBondDipoleBondDipoleOld[BoxII];
      UAdsorbateCationChargeChargeFourier[BoxII]+=UAdsorbateCationChargeChargeFourierDelta[BoxII];
      UAdsorbateCationChargeBondDipoleFourier[BoxII]+=UAdsorbateCationChargeBondDipoleFourierDelta[BoxII];
      UAdsorbateCationBondDipoleBondDipoleFourier[BoxII]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxII];
      UAdsorbateCationCoulomb[BoxII]+=UAdsorbateCationChargeChargeFourierDelta[BoxII]+UAdsorbateChargeChargeNew[BoxII]-UAdsorbateChargeChargeOld[BoxII]+
                                      UAdsorbateCationChargeBondDipoleFourierDelta[BoxII]+UAdsorbateChargeBondDipoleNew[BoxII]-UAdsorbateChargeBondDipoleOld[BoxII]+
                                      UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxII]+UAdsorbateBondDipoleBondDipoleNew[BoxII]-UAdsorbateBondDipoleBondDipoleOld[BoxII];
      UAdsorbateCation[BoxII]+=UAdsorbateCationChargeChargeFourierDelta[BoxII]+UAdsorbateChargeChargeNew[BoxII]-UAdsorbateChargeChargeOld[BoxII]+
                               UAdsorbateCationChargeBondDipoleFourierDelta[BoxII]+UAdsorbateChargeBondDipoleNew[BoxII]-UAdsorbateChargeBondDipoleOld[BoxII]+
                               UAdsorbateCationBondDipoleBondDipoleFourierDelta[BoxII]+UAdsorbateBondDipoleBondDipoleNew[BoxII]-UAdsorbateBondDipoleBondDipoleOld[BoxII];

      // register the changes in the stored structure factors
      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      {
        CurrentSystem=BoxII;
        AcceptEwaldCationMove(BoxII);
      }
    }

    // swap the alocate memory for the atoms and the groups
    SWAP(Cations[BoxI][CationMoleculeA].Atoms,Cations[BoxII][CationMoleculeB].Atoms,atom_temp);
    SWAP(Cations[BoxI][CationMoleculeA].Groups,Cations[BoxII][CationMoleculeB].Groups,group_temp);

    // register the changes in the types for the 'Old'-component atoms for system A
    nr_atoms=Components[ComponentA].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      type=Components[ComponentA].Type[i];
      NumberOfPseudoAtomsType[BoxI][type]--;
    }

    // overwrite the 'Old'-component with the 'New'-component for system A
    nr_atoms=Components[ComponentB].NumberOfAtoms;
    Cations[BoxI][CationMoleculeA].NumberOfAtoms=nr_atoms;
    Cations[BoxI][CationMoleculeA].Type=ComponentB;
    for(i=0;i<nr_atoms;i++)
    {
      Cations[BoxI][CationMoleculeA].Atoms[i].Position=TrialPosition[BoxI][i];
      type=Components[ComponentB].Type[i];
      Cations[BoxI][CationMoleculeA].Atoms[i].Type=type;
      Cations[BoxI][CationMoleculeA].Atoms[i].Charge=Components[ComponentB].Charge[i];
      NumberOfPseudoAtomsType[BoxI][type]++;
    }

    // register the changes in the types for the 'Old'-component atoms for system B
    nr_atoms=Components[ComponentB].NumberOfAtoms;
    for(i=0;i<nr_atoms;i++)
    {
      type=Components[ComponentB].Type[i];
      NumberOfPseudoAtomsType[BoxII][type]--;
    }

    // overwrite the 'New'-component with the 'Old'-component for system A
    nr_atoms=Components[ComponentA].NumberOfAtoms;
    Cations[BoxII][CationMoleculeB].NumberOfAtoms=nr_atoms;
    Cations[BoxII][CationMoleculeB].Type=ComponentA;
    for(i=0;i<nr_atoms;i++)
    {
      Cations[BoxII][CationMoleculeB].Atoms[i].Position=TrialPosition[BoxII][i];
      type=Components[ComponentA].Type[i];
      Cations[BoxII][CationMoleculeB].Atoms[i].Type=type;
      Cations[BoxII][CationMoleculeB].Atoms[i].Charge=Components[ComponentA].Charge[i];
      NumberOfPseudoAtomsType[BoxII][type]++;
    }


    // register the changes in the amount of molecules of each type
    Components[ComponentA].NumberOfMolecules[BoxI]--;
    Components[ComponentB].NumberOfMolecules[BoxI]++;

    Components[ComponentA].NumberOfMolecules[BoxII]++;
    Components[ComponentB].NumberOfMolecules[BoxII]--;

    // update the center of mass of the molecule
    CurrentSystem=BoxI;
    UpdateGroupCenterOfMassCation(CationMoleculeA);

    CurrentSystem=BoxII;
    UpdateGroupCenterOfMassCation(CationMoleculeB);

    // recompute the total energy
    UTotal[BoxI]=
         // inter
         UHostHost[CurrentSystem]+UCationCation[CurrentSystem]+UAdsorbateAdsorbate[CurrentSystem]+
         UHostCation[CurrentSystem]+UHostAdsorbate[CurrentSystem]+UAdsorbateCation[CurrentSystem]+
         // intra framework
         UHostBond[CurrentSystem]+UHostUreyBradley[CurrentSystem]+UHostBend[CurrentSystem]+
         UHostInversionBend[CurrentSystem]+UHostTorsion[CurrentSystem]+UHostImproperTorsion[CurrentSystem]+
         UHostBondBond[CurrentSystem]+UHostBondBend[CurrentSystem]+UHostBendBend[CurrentSystem]+
         UHostBondTorsion[CurrentSystem]+UHostBendTorsion[CurrentSystem]+
         // intra adsorbates
         UAdsorbateBond[CurrentSystem]+UAdsorbateUreyBradley[CurrentSystem]+UAdsorbateBend[CurrentSystem]+
         UAdsorbateInversionBend[CurrentSystem]+UAdsorbateTorsion[CurrentSystem]+UAdsorbateImproperTorsion[CurrentSystem]+
         UAdsorbateBondBond[CurrentSystem]+UAdsorbateBondBend[CurrentSystem]+UAdsorbateBendBend[CurrentSystem]+
         UAdsorbateBondTorsion[CurrentSystem]+UAdsorbateBendTorsion[CurrentSystem]+
         UAdsorbateIntraVDW[CurrentSystem]+UAdsorbateIntraChargeCharge[CurrentSystem]+
         UAdsorbateIntraChargeBondDipole[CurrentSystem]+UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+
         // intra cations
         UCationBond[CurrentSystem]+UCationUreyBradley[CurrentSystem]+UCationBend[CurrentSystem]+
         UCationInversionBend[CurrentSystem]+UCationTorsion[CurrentSystem]+UCationImproperTorsion[CurrentSystem]+
         UCationBondBond[CurrentSystem]+UCationBondBend[CurrentSystem]+UCationBendBend[CurrentSystem]+
         UCationBondTorsion[CurrentSystem]+UCationBendTorsion[CurrentSystem]+
         UCationIntraVDW[CurrentSystem]+UCationIntraChargeCharge[CurrentSystem]+
         UCationIntraChargeBondDipole[CurrentSystem]+UCationIntraBondDipoleBondDipole[CurrentSystem]+
         // polarization
         UHostPolarization[CurrentSystem]+UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem]+
         // tailcorrection
         UTailCorrection[CurrentSystem];

    UTotal[BoxII]=
         // inter
         UHostHost[CurrentSystem]+UCationCation[CurrentSystem]+UAdsorbateAdsorbate[CurrentSystem]+
         UHostCation[CurrentSystem]+UHostAdsorbate[CurrentSystem]+UAdsorbateCation[CurrentSystem]+
         // intra framework
         UHostBond[CurrentSystem]+UHostUreyBradley[CurrentSystem]+UHostBend[CurrentSystem]+
         UHostInversionBend[CurrentSystem]+UHostTorsion[CurrentSystem]+UHostImproperTorsion[CurrentSystem]+
         UHostBondBond[CurrentSystem]+UHostBondBend[CurrentSystem]+UHostBendBend[CurrentSystem]+
         UHostBondTorsion[CurrentSystem]+UHostBendTorsion[CurrentSystem]+
         // intra adsorbates
         UAdsorbateBond[CurrentSystem]+UAdsorbateUreyBradley[CurrentSystem]+UAdsorbateBend[CurrentSystem]+
         UAdsorbateInversionBend[CurrentSystem]+UAdsorbateTorsion[CurrentSystem]+UAdsorbateImproperTorsion[CurrentSystem]+
         UAdsorbateBondBond[CurrentSystem]+UAdsorbateBondBend[CurrentSystem]+UAdsorbateBendBend[CurrentSystem]+
         UAdsorbateBondTorsion[CurrentSystem]+UAdsorbateBendTorsion[CurrentSystem]+
         UAdsorbateIntraVDW[CurrentSystem]+UAdsorbateIntraChargeCharge[CurrentSystem]+
         UAdsorbateIntraChargeBondDipole[CurrentSystem]+UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+
         // intra cations
         UCationBond[CurrentSystem]+UCationUreyBradley[CurrentSystem]+UCationBend[CurrentSystem]+
         UCationInversionBend[CurrentSystem]+UCationTorsion[CurrentSystem]+UCationImproperTorsion[CurrentSystem]+
         UCationBondBond[CurrentSystem]+UCationBondBend[CurrentSystem]+UCationBendBend[CurrentSystem]+
         UCationBondTorsion[CurrentSystem]+UCationBendTorsion[CurrentSystem]+
         UCationIntraVDW[CurrentSystem]+UCationIntraChargeCharge[CurrentSystem]+
         UCationIntraChargeBondDipole[CurrentSystem]+UCationIntraBondDipoleBondDipole[CurrentSystem]+
         // polarization
         UHostPolarization[CurrentSystem]+UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem]+
         // tailcorrection
         UTailCorrection[CurrentSystem];
  }
  CurrentSystem=CurrentSystemStored;

  return 0;
}


void GibbsIdentityChangeMove(void)
{
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
    GibbsIdentityChangeCationMove();
  else
    GibbsIdentityChangeAdsorbateMove();
}


void PrintGibbsIdentityChangeStatistics(FILE *FilePtr)
{
  int i,j,k,MoveUsed;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfGibbsIdentityChangeMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the Gibbs identity change move:\n");
    fprintf(FilePtr,"==============================================================\n");

    for(j=0;j<NumberOfComponents;j++)
    {
      for(k=0;k<NumberOfComponents;k++)
      {
        fprintf(FilePtr,"System [%d] ([%s]->[%s]) total tried: %lf succesfull growth: %lf (%lf [%%]) accepted: %lf (%lf [%%])\n",
             CurrentSystem,
             Components[j].Name,
             Components[k].Name,
             (double)GibbsIdentityChangeAttempts[CurrentSystem][j][k],
             (double)GibbsIdentityChangeAccepted[CurrentSystem][j][k][0],
             (double)(GibbsIdentityChangeAttempts[CurrentSystem][j][k]>(REAL)0.0?
               100.0*GibbsIdentityChangeAccepted[CurrentSystem][j][k][0]/GibbsIdentityChangeAttempts[CurrentSystem][j][k]:(REAL)0.0),
             (double)GibbsIdentityChangeAccepted[CurrentSystem][j][k][1],
             (double)(GibbsIdentityChangeAttempts[CurrentSystem][j][k]>(REAL)0.0?
               100.0*GibbsIdentityChangeAccepted[CurrentSystem][j][k][1]/GibbsIdentityChangeAttempts[CurrentSystem][j][k]:(REAL)0.0));
      }
    }
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"Gibbs identity change move was OFF for all components\n\n");
}

/*********************************************************************************************************
 * Name       | HybrideNVEMove                                                                           *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | A short molecular dynamics trajetory is used as a Monte Carlo move.                      *
 * Parameters | -                                                                                        *
 * Acc. rule  | RandomNumber()<exp(-Beta dU)                                                             *
 * Note       | All atoms are given a random velocity from the Maxwell-Boltzmann distribution. The       *
 *            | system is integrated usin a time-step DeltaT for N steps in the NVE ensemble. The dU in  *
 *            | the acceptance rule is the change in total energy (kinetic+potential energy) which       *
 *            | corresponds to the integrayion error. Accurate integration is not the main priority.     *
 *            | The algorithm is theoretical exact, irrespective of Dt and N. The hybrid-move is by far  *
 *            | the most efficient move for flexible frameworks. Not only the framework is moved, but    *
 *            | also all molecules in the system.                                                        *
 * Refs.      | S. Duane, A.D. Kennedy, B.J. Pendleton, and D. Roeth,                                    *
 *            | Physcial Letters B, Vol. 195, pages 217-, 1987.                                          *
 *            | B. Mehlig, D. Heermanm, and B. Forrest                                                   *
 *            | Physical Review B, Vol. 45, pages 679-, 1992.                                            *
 *            | S. chempath, L.A. Clark, and R.Q. Snurr                                                  *
 *            | Journal of Chemical Physics, Vol. 118(16), pages 7365-7643, 2003.                        *
 *********************************************************************************************************/

void HybridNVEMove(void)
{
  int i;
  REAL Drift;
  REAL ReferenceEnergy;

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

  REAL StoredUKinetic,StoredUHostKinetic,StoredUAdsorbateTranslationalKinetic;
  REAL StoredUCationTranslationalKinetic,StoredUAdsorbateRotationalKinetic;
  REAL StoredUCationRotationalKinetic,StoredUAdsorbateKinetic,StoredUCationKinetic;

  CurrentSystem=(int)(RandomNumber()*NumberOfSystems);

  HybridNVEAttempts[CurrentSystem]+=1.0;

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


  // store the positions of the framework
  SaveFrameworkPositionsToReferenceValues();
  SaveAdsorbateAtomPositionsToReferenceValues();
  SaveCationAtomPositionsToReferenceValues();

  // store the structure-factors for the Ewald-summations
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    SaveCurrentEwaldStructureFactors(0,CurrentSystem);

  Ensemble[CurrentSystem]=NVE;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    InitializeVelocityAdsorbate(i);

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    InitializeVelocityCation(i);

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    InitializeFrameworkVelocities();

  StoredUKinetic=UKinetic[CurrentSystem];
  StoredUHostKinetic=UHostKinetic[CurrentSystem];
  StoredUAdsorbateTranslationalKinetic=UAdsorbateTranslationalKinetic[CurrentSystem];
  StoredUCationTranslationalKinetic=UCationTranslationalKinetic[CurrentSystem];
  StoredUAdsorbateRotationalKinetic=UAdsorbateRotationalKinetic[CurrentSystem];
  StoredUCationRotationalKinetic=UCationRotationalKinetic[CurrentSystem];
  StoredUAdsorbateKinetic=UAdsorbateKinetic[CurrentSystem];
  StoredUCationKinetic=UCationKinetic[CurrentSystem];

  InitializeForces();

  ReferenceEnergy=ConservedEnergy[CurrentSystem];
  Drift=0.0;

  // integrated the system 'NumberOfHybridNVESteps' steps
  for(i=0;i<NumberOfHybridNVESteps;i++)
  {
    // evolve the system a full time-step
    Integration();

    // update the drift in the energy
    Drift+=fabs((ConservedEnergy[CurrentSystem]-ReferenceEnergy)/ReferenceEnergy);
  }

  if((RandomNumber()<exp(-Beta[CurrentSystem]*(ConservedEnergy[CurrentSystem]-ReferenceEnergy)))&&(finite(Drift)))
  {
    HybridNVEAccepted[CurrentSystem]+=1.0;

    // register the starting temperatures
    if(DegreesOfFreedom[CurrentSystem]>0)
    {
      HybridNVEStartTemperature[CurrentSystem]+=2.0*StoredUKinetic/(K_B*DegreesOfFreedom[CurrentSystem]);
      HybridNVEStartTemperatureCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomTranslation[CurrentSystem]>0)
    {
      HybridNVEStartTranslationalTemperature[CurrentSystem]+=2.0*(StoredUHostKinetic+
         StoredUAdsorbateTranslationalKinetic+StoredUCationTranslationalKinetic)/
                                                (K_B*DegreesOfFreedomTranslation[CurrentSystem]);
      HybridNVEStartTemperatureTranslationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomRotation[CurrentSystem]>0)
    {
      HybridNVEStartRotationalTemperature[CurrentSystem]+=2.0*(StoredUAdsorbateRotationalKinetic+
            StoredUCationRotationalKinetic)/(K_B*DegreesOfFreedomRotation[CurrentSystem]);
      HybridNVEStartTemperatureRotationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomFramework[CurrentSystem]>0)
    {
      HybridNVEStartTemperatureFramework[CurrentSystem]+=2.0*StoredUHostKinetic/(K_B*DegreesOfFreedomFramework[CurrentSystem]);
      HybridNVEStartTemperatureFrameworkCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomAdsorbates[CurrentSystem]>0)
    {
      HybridNVEStartTemperatureAdsorbate[CurrentSystem]+=2.0*StoredUAdsorbateKinetic/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem]);
      HybridNVEStartTemperatureAdsorbateCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomCations[CurrentSystem]>0)
    {
      HybridNVEStartTemperatureCation[CurrentSystem]+=2.0*StoredUCationKinetic/(K_B*DegreesOfFreedomCations[CurrentSystem]);
      HybridNVEStartTemperatureCationCount[CurrentSystem]+=1.0;
    }

    // register the end temperatures
    HybridNVEDrift[CurrentSystem]+=Drift;
    HybridNVEDriftCount[CurrentSystem]+=1.0;

    if(DegreesOfFreedom[CurrentSystem]>0)
    {
      HybridNVEEndTemperature[CurrentSystem]+=2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem]);
      HybridNVEEndTemperatureCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomTranslation[CurrentSystem]>0)
    {
      HybridNVEEndTranslationalTemperature[CurrentSystem]+=2.0*(UHostKinetic[CurrentSystem]+
          UAdsorbateTranslationalKinetic[CurrentSystem]+UCationTranslationalKinetic[CurrentSystem])/
                                              (K_B*DegreesOfFreedomTranslation[CurrentSystem]);
      HybridNVEEndTemperatureTranslationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomRotation[CurrentSystem]>0)
    {
      HybridNVEEndRotationalTemperature[CurrentSystem]+=2.0*(UAdsorbateRotationalKinetic[CurrentSystem]+
             UCationRotationalKinetic[CurrentSystem])/(K_B*DegreesOfFreedomRotation[CurrentSystem]);
      HybridNVEEndTemperatureRotationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomFramework[CurrentSystem]>0)
    {
      HybridNVEEndTemperatureFramework[CurrentSystem]+=2.0*UHostKinetic[CurrentSystem]/(K_B*DegreesOfFreedomFramework[CurrentSystem]);
      HybridNVEEndTemperatureFrameworkCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomAdsorbates[CurrentSystem]>0)
    {
      HybridNVEEndTemperatureAdsorbate[CurrentSystem]+=2.0*UAdsorbateKinetic[CurrentSystem]/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem]);
      HybridNVEEndTemperatureAdsorbateCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomCations[CurrentSystem]>0)
    {
      HybridNVEEndTemperatureCation[CurrentSystem]+=2.0*UCationKinetic[CurrentSystem]/(K_B*DegreesOfFreedomCations[CurrentSystem]);
      HybridNVEEndTemperatureCationCount[CurrentSystem]+=1.0;
    }

    // register the changes in the stored structure factors
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      PrecomputeTotalEwaldContributions();
  }
  else
  {
    // restore all the energy to the Old state
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
    UTailCorrection[CurrentSystem]=StoredUTailCorrection;

    UHostPolarization[CurrentSystem]=UHostPolarizationStored;
    UAdsorbatePolarization[CurrentSystem]=UHostPolarizationStored;
    UCationPolarization[CurrentSystem]=UCationPolarizationStored;

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationStored;
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationStored;
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationStored;

    UTotal[CurrentSystem]=StoredUTotal;
    UTailCorrection[CurrentSystem]=StoredUTailCorrection;

    // restore all the positions to the Old state
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      RetrieveStoredEwaldStructureFactors(0,CurrentSystem);

    //ADDED
    CalculateAnisotropicSites();
  }
}

void PrintHybridNVEStatistics(FILE *FilePtr)
{
  if(ProbabilityHybridNVEMove)
  {
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Performance of the hybrid MCMD in the NVE-ensemble:\n");
    fprintf(FilePtr,"==============================================================\n");

    fprintf(FilePtr,"total tried: %lf accepted: %lf (%lf [%%])\n\n",
        (double)HybridNVEAttempts[CurrentSystem],
        (double)HybridNVEAccepted[CurrentSystem],
        (double)(HybridNVEAttempts[CurrentSystem]>(REAL)0.0?
          100.0*HybridNVEAccepted[CurrentSystem]/HybridNVEAttempts[CurrentSystem]:(REAL)0.0));

    fprintf(FilePtr,"total amount of MD-time simulated: %18.10lf [ps]\n\n",
        (double)((REAL)NumberOfHybridNVESteps*DeltaT*HybridNVEAccepted[CurrentSystem]));

    if(HybridNVEAccepted[CurrentSystem]>0.0)
    {
      fprintf(FilePtr,"\tAverage drift in the energy:               % 18.10lf\n\n",
                  (double)(HybridNVEDrift[CurrentSystem]/HybridNVEDriftCount[CurrentSystem]));

      if(HybridNVEStartTemperatureCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature:               % 18.10lf\n",
           (double)(HybridNVEStartTemperature[CurrentSystem]/HybridNVEStartTemperatureCount[CurrentSystem]));
      if(HybridNVEStartTemperatureTranslationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature (translation): % 18.10lf\n",
           (double)(HybridNVEStartTranslationalTemperature[CurrentSystem]/HybridNVEStartTemperatureTranslationCount[CurrentSystem]));
      if(HybridNVEStartTemperatureRotationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature (rotation): % 18.10lf\n",
           (double)(HybridNVEStartRotationalTemperature[CurrentSystem]/HybridNVEStartTemperatureRotationCount[CurrentSystem]));
      if(HybridNVEStartTemperatureFrameworkCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature framework : % 18.10lf\n",
           (double)(HybridNVEStartTemperatureFramework[CurrentSystem]/HybridNVEStartTemperatureFrameworkCount[CurrentSystem]));
      if(HybridNVEStartTemperatureAdsorbateCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature adsorbates: % 18.10lf\n",
           (double)(HybridNVEStartTemperatureAdsorbate[CurrentSystem]/HybridNVEStartTemperatureAdsorbateCount[CurrentSystem]));
      if(HybridNVEStartTemperatureCationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature cations: % 18.10lf\n\n",
           (double)(HybridNVEStartTemperatureCation[CurrentSystem]/HybridNVEStartTemperatureCationCount[CurrentSystem]));

      if(HybridNVEEndTemperatureCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature: % 18.10lf\n",
          (double)(HybridNVEEndTemperature[CurrentSystem]/HybridNVEEndTemperatureCount[CurrentSystem]));
      if(HybridNVEEndTemperatureTranslationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature (translation): % 18.10lf\n",
           (double)(HybridNVEEndTranslationalTemperature[CurrentSystem]/HybridNVEEndTemperatureTranslationCount[CurrentSystem]));
      if(HybridNVEEndTemperatureRotationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature (rotation): % 18.10lf\n",
           (double)(HybridNVEEndRotationalTemperature[CurrentSystem]/HybridNVEEndTemperatureRotationCount[CurrentSystem]));
      if(HybridNVEEndTemperatureFrameworkCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature framework : % 18.10lf\n",
           (double)(HybridNVEEndTemperatureFramework[CurrentSystem]/HybridNVEEndTemperatureFrameworkCount[CurrentSystem]));
      if(HybridNVEEndTemperatureAdsorbateCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature adsorbates: % 18.10lf\n",
           (double)(HybridNVEEndTemperatureAdsorbate[CurrentSystem]/HybridNVEEndTemperatureAdsorbateCount[CurrentSystem]));
      if(HybridNVEEndTemperatureCationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature cations: % 18.10lf\n\n",
           (double)(HybridNVEEndTemperatureCation[CurrentSystem]/HybridNVEEndTemperatureCationCount[CurrentSystem]));
    }
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"Hybrid MC/MD move in the NVE-ensemble was OFF\n\n");
}



void HybridNPHMove(void)
{
  int i;
  REAL Drift,det;
  REAL ReferenceEnergy;

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
  REAL StoredUTotal,StoredUIon,StoredUTailCorrection;

  REAL UHostPolarizationStored,UAdsorbatePolarizationStored,UCationPolarizationStored;
  REAL UHostBackPolarizationStored,UAdsorbateBackPolarizationStored,UCationBackPolarizationStored;

  REAL StoredUKinetic,StoredUHostKinetic,StoredUAdsorbateTranslationalKinetic;
  REAL StoredUCationTranslationalKinetic,StoredUAdsorbateRotationalKinetic;
  REAL StoredUCationRotationalKinetic,StoredUAdsorbateKinetic,StoredUCationKinetic;
  REAL StoredCellTemperature;

  REAL StoredVolume;
  REAL_MATRIX3x3 StoredBox;

  CurrentSystem=(int)(RandomNumber()*NumberOfSystems);

  HybridNPHAttempts[CurrentSystem]+=1.0;

  StoredBox=Box[CurrentSystem];
  StoredVolume=Volume[CurrentSystem];

  StoredUTotal=UTotal[CurrentSystem];
  StoredUIon=UIon[CurrentSystem];
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

  // store the positions of the framework
  SaveFrameworkPositionsToReferenceValues();
  SaveAdsorbateAtomPositionsToReferenceValues();
  SaveCationAtomPositionsToReferenceValues();

  // store the structure-factors for the Ewald-summations
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    SaveCurrentEwaldStructureFactors(0,CurrentSystem);
    SaveCurrentKVectors(0,CurrentSystem);
  }

  Ensemble[CurrentSystem]=NPH;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    InitializeVelocityAdsorbate(i);

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    InitializeVelocityCation(i);

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    InitializeFrameworkVelocities();

  InitializeForces();
  InitializeNoseHooverCurrentSystem();
  InitializeBoxVelocities();

  StoredCellTemperature=GetCellTemperature();
  StoredUKinetic=UKinetic[CurrentSystem];
  StoredUHostKinetic=UHostKinetic[CurrentSystem];
  StoredUAdsorbateTranslationalKinetic=UAdsorbateTranslationalKinetic[CurrentSystem];
  StoredUCationTranslationalKinetic=UCationTranslationalKinetic[CurrentSystem];
  StoredUAdsorbateRotationalKinetic=UAdsorbateRotationalKinetic[CurrentSystem];
  StoredUCationRotationalKinetic=UCationRotationalKinetic[CurrentSystem];
  StoredUAdsorbateKinetic=UAdsorbateKinetic[CurrentSystem];
  StoredUCationKinetic=UCationKinetic[CurrentSystem];

  ReferenceEnergy=ConservedEnergy[CurrentSystem];
  Drift=0.0;

  // integrated the system 'NumberOfHybridNPHSteps' steps
  for(i=0;i<NumberOfHybridNPHSteps;i++)
  {
    // evolve the system a full time-step
    Integration();

    // update the drift in the energy
    Drift+=fabs((ConservedEnergy[CurrentSystem]-ReferenceEnergy)/ReferenceEnergy);
  }

  //if((RandomNumber()<exp(-Beta[CurrentSystem]*(ConservedEnergy[CurrentSystem]-ReferenceEnergy))*SQR(Volume[CurrentSystem]/StoredVolume))&&(finite(Drift)))
  if((RandomNumber()<exp(-Beta[CurrentSystem]*(ConservedEnergy[CurrentSystem]-ReferenceEnergy)))&&(finite(Drift)))
  {
    HybridNPHAccepted[CurrentSystem]+=1.0;

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassCation(i);

    // register the starting temperatures
    if(DegreesOfFreedom[CurrentSystem]>0)
    {
      HybridNPHStartCellTemperature[CurrentSystem]+=StoredCellTemperature;
      HybridNPHStartTemperature[CurrentSystem]+=2.0*StoredUKinetic/(K_B*DegreesOfFreedom[CurrentSystem]);
      HybridNPHStartTemperatureCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomTranslation[CurrentSystem]>0)
    {
      HybridNPHStartTranslationalTemperature[CurrentSystem]+=2.0*(StoredUHostKinetic+
         StoredUAdsorbateTranslationalKinetic+StoredUCationTranslationalKinetic)/
                                                (K_B*DegreesOfFreedomTranslation[CurrentSystem]);
      HybridNPHStartTemperatureTranslationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomRotation[CurrentSystem]>0)
    {
      HybridNPHStartRotationalTemperature[CurrentSystem]+=2.0*(StoredUAdsorbateRotationalKinetic+
            StoredUCationRotationalKinetic)/(K_B*DegreesOfFreedomRotation[CurrentSystem]);
      HybridNPHStartTemperatureRotationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomFramework[CurrentSystem]>0)
    {
      HybridNPHStartTemperatureFramework[CurrentSystem]+=2.0*StoredUHostKinetic/(K_B*DegreesOfFreedomFramework[CurrentSystem]);
      HybridNPHStartTemperatureFrameworkCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomAdsorbates[CurrentSystem]>0)
    {
      HybridNPHStartTemperatureAdsorbate[CurrentSystem]+=2.0*StoredUAdsorbateKinetic/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem]);
      HybridNPHStartTemperatureAdsorbateCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomCations[CurrentSystem]>0)
    {
      HybridNPHStartTemperatureCation[CurrentSystem]+=2.0*StoredUCationKinetic/(K_B*DegreesOfFreedomCations[CurrentSystem]);
      HybridNPHStartTemperatureCationCount[CurrentSystem]+=1.0;
    }

    // register the end temperatures
    HybridNPHDrift[CurrentSystem]+=Drift;
    HybridNPHDriftCount[CurrentSystem]+=1.0;

    if(DegreesOfFreedom[CurrentSystem]>0)
    {
      HybridNPHEndCellTemperature[CurrentSystem]+=GetCellTemperature();
      HybridNPHEndTemperature[CurrentSystem]+=2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem]);
      HybridNPHEndTemperatureCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomTranslation[CurrentSystem]>0)
    {
      HybridNPHEndTranslationalTemperature[CurrentSystem]+=2.0*(UHostKinetic[CurrentSystem]+
          UAdsorbateTranslationalKinetic[CurrentSystem]+UCationTranslationalKinetic[CurrentSystem])/
                                              (K_B*DegreesOfFreedomTranslation[CurrentSystem]);
      HybridNPHEndTemperatureTranslationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomRotation[CurrentSystem]>0)
    {
      HybridNPHEndRotationalTemperature[CurrentSystem]+=2.0*(UAdsorbateRotationalKinetic[CurrentSystem]+
             UCationRotationalKinetic[CurrentSystem])/(K_B*DegreesOfFreedomRotation[CurrentSystem]);
      HybridNPHEndTemperatureRotationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomFramework[CurrentSystem]>0)
    {
      HybridNPHEndTemperatureFramework[CurrentSystem]+=2.0*UHostKinetic[CurrentSystem]/(K_B*DegreesOfFreedomFramework[CurrentSystem]);
      HybridNPHEndTemperatureFrameworkCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomAdsorbates[CurrentSystem]>0)
    {
      HybridNPHEndTemperatureAdsorbate[CurrentSystem]+=2.0*UAdsorbateKinetic[CurrentSystem]/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem]);
      HybridNPHEndTemperatureAdsorbateCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomCations[CurrentSystem]>0)
    {
      HybridNPHEndTemperatureCation[CurrentSystem]+=2.0*UCationKinetic[CurrentSystem]/(K_B*DegreesOfFreedomCations[CurrentSystem]);
      HybridNPHEndTemperatureCationCount[CurrentSystem]+=1.0;
    }
  }
  else
  {
    // restore all the energy to the Old state
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
    UTailCorrection[CurrentSystem]=StoredUTailCorrection;

    UHostPolarization[CurrentSystem]=UHostPolarizationStored;
    UAdsorbatePolarization[CurrentSystem]=UHostPolarizationStored;
    UCationPolarization[CurrentSystem]=UCationPolarizationStored;

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationStored;
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationStored;
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationStored;

    UTotal[CurrentSystem]=StoredUTotal;
    UIon[CurrentSystem]=StoredUIon;
    UTailCorrection[CurrentSystem]=StoredUTailCorrection;


    Box[CurrentSystem]=StoredBox;
    Volume[CurrentSystem]=StoredVolume;
    Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
    CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&det);

    // restore all the positions to the Old state
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassCation(i);

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      RetrieveStoredEwaldStructureFactors(0,CurrentSystem);
      RetrieveStoredKVectors(0,CurrentSystem);
    }
  }
}

void PrintHybridNPHStatistics(FILE *FilePtr)
{
  if(ProbabilityHybridNPHMove)
  {
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Performance of the hybrid MCMD move in the NPH-ensemble:\n");
    fprintf(FilePtr,"==============================================================\n");

    fprintf(FilePtr,"total tried: %lf accepted: %lf (%lf [%%])\n\n",
        (double)HybridNPHAttempts[CurrentSystem],
        (double)HybridNPHAccepted[CurrentSystem],
        (double)(HybridNPHAttempts[CurrentSystem]>(REAL)0.0?
          100.0*HybridNPHAccepted[CurrentSystem]/HybridNPHAttempts[CurrentSystem]:(REAL)0.0));

    fprintf(FilePtr,"total amount of MD-time simulated: %18.10lf [ps]\n\n",
        (double)((REAL)NumberOfHybridNPHSteps*DeltaT*HybridNPHAccepted[CurrentSystem]));

    if(HybridNPHAccepted[CurrentSystem]>0.0)
    {
      fprintf(FilePtr,"\tAverage drift in the energy:               % 18.10lf\n\n",
                  (double)(HybridNPHDrift[CurrentSystem]/HybridNPHDriftCount[CurrentSystem]));

      if(HybridNPHStartTemperatureCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature:               % 18.10lf\n",
           (double)(HybridNPHStartTemperature[CurrentSystem]/HybridNPHStartTemperatureCount[CurrentSystem]));
      if(HybridNPHStartTemperatureCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature cell:               % 18.10lf\n",
           (double)(HybridNPHStartCellTemperature[CurrentSystem]/HybridNPHStartTemperatureCount[CurrentSystem]));
      if(HybridNPHStartTemperatureTranslationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature (translation): % 18.10lf\n",
           (double)(HybridNPHStartTranslationalTemperature[CurrentSystem]/HybridNPHStartTemperatureTranslationCount[CurrentSystem]));
      if(HybridNPHStartTemperatureRotationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature (rotation): % 18.10lf\n",
           (double)(HybridNPHStartRotationalTemperature[CurrentSystem]/HybridNPHStartTemperatureRotationCount[CurrentSystem]));
      if(HybridNPHStartTemperatureFrameworkCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature framework : % 18.10lf\n",
           (double)(HybridNPHStartTemperatureFramework[CurrentSystem]/HybridNPHStartTemperatureFrameworkCount[CurrentSystem]));
      if(HybridNPHStartTemperatureAdsorbateCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature adsorbates: % 18.10lf\n",
           (double)(HybridNPHStartTemperatureAdsorbate[CurrentSystem]/HybridNPHStartTemperatureAdsorbateCount[CurrentSystem]));
      if(HybridNPHStartTemperatureCationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature cations: % 18.10lf\n\n",
           (double)(HybridNPHStartTemperatureCation[CurrentSystem]/HybridNPHStartTemperatureCationCount[CurrentSystem]));

      if(HybridNPHEndTemperatureCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature: % 18.10lf\n",
          (double)(HybridNPHEndTemperature[CurrentSystem]/HybridNPHEndTemperatureCount[CurrentSystem]));
      if(HybridNPHEndTemperatureCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature (cell): % 18.10lf\n",
          (double)(HybridNPHEndCellTemperature[CurrentSystem]/HybridNPHEndTemperatureCount[CurrentSystem]));
      if(HybridNPHEndTemperatureTranslationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature (translation): % 18.10lf\n",
           (double)(HybridNPHEndTranslationalTemperature[CurrentSystem]/HybridNPHEndTemperatureTranslationCount[CurrentSystem]));
      if(HybridNPHEndTemperatureRotationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature (rotation): % 18.10lf\n",
           (double)(HybridNPHEndRotationalTemperature[CurrentSystem]/HybridNPHEndTemperatureRotationCount[CurrentSystem]));
      if(HybridNPHEndTemperatureFrameworkCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature framework : % 18.10lf\n",
           (double)(HybridNPHEndTemperatureFramework[CurrentSystem]/HybridNPHEndTemperatureFrameworkCount[CurrentSystem]));
      if(HybridNPHEndTemperatureAdsorbateCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature adsorbates: % 18.10lf\n",
           (double)(HybridNPHEndTemperatureAdsorbate[CurrentSystem]/HybridNPHEndTemperatureAdsorbateCount[CurrentSystem]));
      if(HybridNPHEndTemperatureCationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature cations: % 18.10lf\n\n",
           (double)(HybridNPHEndTemperatureCation[CurrentSystem]/HybridNPHEndTemperatureCationCount[CurrentSystem]));
    }
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"Hybrid MC/MD in the NPH-ensemble move was OFF\n\n");
}


void HybridNPHPRMove(void)
{
  int i;
  REAL Drift,det;
  REAL ReferenceEnergy;

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
  REAL StoredUTotal,StoredUIon,StoredUTailCorrection;

  REAL UHostPolarizationStored,UAdsorbatePolarizationStored,UCationPolarizationStored;
  REAL UHostBackPolarizationStored,UAdsorbateBackPolarizationStored,UCationBackPolarizationStored;

  REAL StoredUKinetic,StoredUHostKinetic,StoredUAdsorbateTranslationalKinetic;
  REAL StoredUCationTranslationalKinetic,StoredUAdsorbateRotationalKinetic;
  REAL StoredUCationRotationalKinetic,StoredUAdsorbateKinetic,StoredUCationKinetic;
  REAL StoredCellTemperature;

  REAL StoredVolume;
  REAL_MATRIX3x3 StoredBox;

  CurrentSystem=(int)(RandomNumber()*NumberOfSystems);

  HybridNPHPRAttempts[CurrentSystem]+=1.0;

  StoredBox=Box[CurrentSystem];
  StoredVolume=Volume[CurrentSystem];

  StoredUTotal=UTotal[CurrentSystem];
  StoredUIon=UIon[CurrentSystem];
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


  // store the positions of the framework
  SaveFrameworkPositionsToReferenceValues();
  SaveAdsorbateAtomPositionsToReferenceValues();
  SaveCationAtomPositionsToReferenceValues();

  // store the structure-factors for the Ewald-summations
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    SaveCurrentEwaldStructureFactors(0,CurrentSystem);
    SaveCurrentKVectors(0,CurrentSystem);
  }

  Ensemble[CurrentSystem]=NPHPR;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    InitializeVelocityAdsorbate(i);

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    InitializeVelocityCation(i);

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    InitializeFrameworkVelocities();

  InitializeNoseHooverCurrentSystem();
  InitializeForces();

  StoredUKinetic=UKinetic[CurrentSystem];
  StoredUHostKinetic=UHostKinetic[CurrentSystem];
  StoredUAdsorbateTranslationalKinetic=UAdsorbateTranslationalKinetic[CurrentSystem];
  StoredUCationTranslationalKinetic=UCationTranslationalKinetic[CurrentSystem];
  StoredUAdsorbateRotationalKinetic=UAdsorbateRotationalKinetic[CurrentSystem];
  StoredUCationRotationalKinetic=UCationRotationalKinetic[CurrentSystem];
  StoredUAdsorbateKinetic=UAdsorbateKinetic[CurrentSystem];
  StoredUCationKinetic=UCationKinetic[CurrentSystem];
  StoredCellTemperature=GetCellTemperature();

  ReferenceEnergy=ConservedEnergy[CurrentSystem];
  Drift=0.0;

  // integrated the system 'NumberOfHybridNPHPRSteps' steps
  for(i=0;i<NumberOfHybridNPHPRSteps;i++)
  {
    // evolve the system a full time-step
    Integration();

    // update the drift in the energy
    Drift+=fabs((ConservedEnergy[CurrentSystem]-ReferenceEnergy)/ReferenceEnergy);
  }

  if((RandomNumber()<exp(-Beta[CurrentSystem]*(ConservedEnergy[CurrentSystem]-ReferenceEnergy))*SQR(Volume[CurrentSystem]/StoredVolume))&&(finite(Drift)))
  {
    HybridNPHPRAccepted[CurrentSystem]+=1.0;

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassCation(i);

    // register the starting temperatures
    if(DegreesOfFreedom[CurrentSystem]>0)
    {
      HybridNPHPRStartTemperature[CurrentSystem]+=2.0*StoredUKinetic/(K_B*DegreesOfFreedom[CurrentSystem]);
      HybridNPHPRStartCellTemperature[CurrentSystem]+=StoredCellTemperature;
      HybridNPHPRStartTemperatureCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomTranslation[CurrentSystem]>0)
    {
      HybridNPHPRStartTranslationalTemperature[CurrentSystem]+=2.0*(StoredUHostKinetic+
         StoredUAdsorbateTranslationalKinetic+StoredUCationTranslationalKinetic)/
                                                (K_B*DegreesOfFreedomTranslation[CurrentSystem]);
      HybridNPHPRStartTemperatureTranslationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomRotation[CurrentSystem]>0)
    {
      HybridNPHPRStartRotationalTemperature[CurrentSystem]+=2.0*(StoredUAdsorbateRotationalKinetic+
            StoredUCationRotationalKinetic)/(K_B*DegreesOfFreedomRotation[CurrentSystem]);
      HybridNPHPRStartTemperatureRotationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomFramework[CurrentSystem]>0)
    {
      HybridNPHPRStartTemperatureFramework[CurrentSystem]+=2.0*StoredUHostKinetic/(K_B*DegreesOfFreedomFramework[CurrentSystem]);
      HybridNPHPRStartTemperatureFrameworkCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomAdsorbates[CurrentSystem]>0)
    {
      HybridNPHPRStartTemperatureAdsorbate[CurrentSystem]+=2.0*StoredUAdsorbateKinetic/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem]);
      HybridNPHPRStartTemperatureAdsorbateCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomCations[CurrentSystem]>0)
    {
      HybridNPHPRStartTemperatureCation[CurrentSystem]+=2.0*StoredUCationKinetic/(K_B*DegreesOfFreedomCations[CurrentSystem]);
      HybridNPHPRStartTemperatureCationCount[CurrentSystem]+=1.0;
    }

    // register the end temperatures
    HybridNPHPRDrift[CurrentSystem]+=Drift;
    HybridNPHPRDriftCount[CurrentSystem]+=1.0;

    if(DegreesOfFreedom[CurrentSystem]>0)
    {
      HybridNPHPREndTemperature[CurrentSystem]+=2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem]);
      HybridNPHPREndCellTemperature[CurrentSystem]+=GetCellTemperature();
      HybridNPHPREndTemperatureCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomTranslation[CurrentSystem]>0)
    {
      HybridNPHPREndTranslationalTemperature[CurrentSystem]+=2.0*(UHostKinetic[CurrentSystem]+
          UAdsorbateTranslationalKinetic[CurrentSystem]+UCationTranslationalKinetic[CurrentSystem])/
                                              (K_B*DegreesOfFreedomTranslation[CurrentSystem]);
      HybridNPHPREndTemperatureTranslationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomRotation[CurrentSystem]>0)
    {
      HybridNPHPREndRotationalTemperature[CurrentSystem]+=2.0*(UAdsorbateRotationalKinetic[CurrentSystem]+
             UCationRotationalKinetic[CurrentSystem])/(K_B*DegreesOfFreedomRotation[CurrentSystem]);
      HybridNPHPREndTemperatureRotationCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomFramework[CurrentSystem]>0)
    {
      HybridNPHPREndTemperatureFramework[CurrentSystem]+=2.0*UHostKinetic[CurrentSystem]/(K_B*DegreesOfFreedomFramework[CurrentSystem]);
      HybridNPHPREndTemperatureFrameworkCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomAdsorbates[CurrentSystem]>0)
    {
      HybridNPHPREndTemperatureAdsorbate[CurrentSystem]+=2.0*UAdsorbateKinetic[CurrentSystem]/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem]);
      HybridNPHPREndTemperatureAdsorbateCount[CurrentSystem]+=1.0;
    }

    if(DegreesOfFreedomCations[CurrentSystem]>0)
    {
      HybridNPHPREndTemperatureCation[CurrentSystem]+=2.0*UCationKinetic[CurrentSystem]/(K_B*DegreesOfFreedomCations[CurrentSystem]);
      HybridNPHPREndTemperatureCationCount[CurrentSystem]+=1.0;
    }
  }
  else
  {
    // restore all the energy to the Old state
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
    UTailCorrection[CurrentSystem]=StoredUTailCorrection;

    UHostPolarization[CurrentSystem]=UHostPolarizationStored;
    UAdsorbatePolarization[CurrentSystem]=UHostPolarizationStored;
    UCationPolarization[CurrentSystem]=UCationPolarizationStored;

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationStored;
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationStored;
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationStored;

    UTotal[CurrentSystem]=StoredUTotal;
    UIon[CurrentSystem]=StoredUIon;
    UTailCorrection[CurrentSystem]=StoredUTailCorrection;

    Box[CurrentSystem]=StoredBox;
    Volume[CurrentSystem]=StoredVolume;
    Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
    CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&det);

    // restore all the positions to the Old state
    PlaceFrameworkInBoxFromReferenceValues();
    PlaceAdsorbateAtomsInBoxFromReferenceValues();
    PlaceCationAtomsInBoxFromReferenceValues();

    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassAdsorbate(i);

    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassCation(i);

    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      RetrieveStoredEwaldStructureFactors(0,CurrentSystem);
      RetrieveStoredKVectors(0,CurrentSystem);
    }
  }
}

void PrintHybridNPHPRStatistics(FILE *FilePtr)
{
  if(ProbabilityHybridNPHPRMove)
  {
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Performance of the hybrid MCMD move in the NPH-ensemble (Parrinello-Rahman):\n");
    fprintf(FilePtr,"============================================================================\n");

    fprintf(FilePtr,"total tried: %lf accepted: %lf (%lf [%%])\n\n",
        (double)HybridNPHPRAttempts[CurrentSystem],
        (double)HybridNPHPRAccepted[CurrentSystem],
        (double)(HybridNPHPRAttempts[CurrentSystem]>(REAL)0.0?
          100.0*HybridNPHPRAccepted[CurrentSystem]/HybridNPHPRAttempts[CurrentSystem]:(REAL)0.0));

    fprintf(FilePtr,"total amount of MD-time simulated: %18.10lf [ps]\n\n",
        (double)((REAL)NumberOfHybridNPHPRSteps*DeltaT*HybridNPHPRAccepted[CurrentSystem]));

    if(HybridNPHPRAccepted[CurrentSystem]>0.0)
    {
      fprintf(FilePtr,"\tAverage drift in the energy:               % 18.10lf\n\n",
                  (double)(HybridNPHPRDrift[CurrentSystem]/HybridNPHPRDriftCount[CurrentSystem]));

      if(HybridNPHPRStartTemperatureCount[CurrentSystem]>0)
      {
        fprintf(FilePtr,"\tAverage begin temperature:               % 18.10lf\n",
           (double)(HybridNPHPRStartTemperature[CurrentSystem]/HybridNPHPRStartTemperatureCount[CurrentSystem]));
        fprintf(FilePtr,"\tAverage begin temperature (cell):        % 18.10lf\n",
           (double)(HybridNPHPRStartCellTemperature[CurrentSystem]/HybridNPHPRStartTemperatureCount[CurrentSystem]));
      }
      if(HybridNPHPRStartTemperatureTranslationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature (translation): % 18.10lf\n",
           (double)(HybridNPHPRStartTranslationalTemperature[CurrentSystem]/HybridNPHPRStartTemperatureTranslationCount[CurrentSystem]));
      if(HybridNPHPRStartTemperatureRotationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature (rotation): % 18.10lf\n",
           (double)(HybridNPHPRStartRotationalTemperature[CurrentSystem]/HybridNPHPRStartTemperatureRotationCount[CurrentSystem]));
      if(HybridNPHPRStartTemperatureFrameworkCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature framework : % 18.10lf\n",
           (double)(HybridNPHPRStartTemperatureFramework[CurrentSystem]/HybridNPHPRStartTemperatureFrameworkCount[CurrentSystem]));
      if(HybridNPHPRStartTemperatureAdsorbateCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature adsorbates: % 18.10lf\n",
           (double)(HybridNPHPRStartTemperatureAdsorbate[CurrentSystem]/HybridNPHPRStartTemperatureAdsorbateCount[CurrentSystem]));
      if(HybridNPHPRStartTemperatureCationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage begin temperature cations: % 18.10lf\n\n",
           (double)(HybridNPHPRStartTemperatureCation[CurrentSystem]/HybridNPHPRStartTemperatureCationCount[CurrentSystem]));

      if(HybridNPHPREndTemperatureCount[CurrentSystem]>0)
      {
        fprintf(FilePtr,"\tAverage end temperature: % 18.10lf\n",
          (double)(HybridNPHPREndTemperature[CurrentSystem]/HybridNPHPREndTemperatureCount[CurrentSystem]));
        fprintf(FilePtr,"\tAverage end temperature (cell): % 18.10lf\n",
          (double)(HybridNPHPREndCellTemperature[CurrentSystem]/HybridNPHPREndTemperatureCount[CurrentSystem]));
      }
      if(HybridNPHPREndTemperatureTranslationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature (translation): % 18.10lf\n",
           (double)(HybridNPHPREndTranslationalTemperature[CurrentSystem]/HybridNPHPREndTemperatureTranslationCount[CurrentSystem]));
      if(HybridNPHPREndTemperatureRotationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature (rotation): % 18.10lf\n",
           (double)(HybridNPHPREndRotationalTemperature[CurrentSystem]/HybridNPHPREndTemperatureRotationCount[CurrentSystem]));
      if(HybridNPHPREndTemperatureFrameworkCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature framework : % 18.10lf\n",
           (double)(HybridNPHPREndTemperatureFramework[CurrentSystem]/HybridNPHPREndTemperatureFrameworkCount[CurrentSystem]));
      if(HybridNPHPREndTemperatureAdsorbateCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature adsorbates: % 18.10lf\n",
           (double)(HybridNPHPREndTemperatureAdsorbate[CurrentSystem]/HybridNPHPREndTemperatureAdsorbateCount[CurrentSystem]));
      if(HybridNPHPREndTemperatureCationCount[CurrentSystem]>0)
        fprintf(FilePtr,"\tAverage end temperature cations: % 18.10lf\n\n",
           (double)(HybridNPHPREndTemperatureCation[CurrentSystem]/HybridNPHPREndTemperatureCationCount[CurrentSystem]));
    }
    fprintf(FilePtr,"\n");
  }
  else
    fprintf(FilePtr,"Hybrid MC/MD in the NPH-ensemble (Parrinello-Rahman) move was OFF\n\n");
}



// Surface area move
//
// The surface area of the framework and cations is probed

void SurfaceAreaMove(void)
{
  int i,j,k;
  int typeA,typeB,Overlap,StartingBead;
  int NumberOfTrials;
  REAL equilibrium_distance,well_depth_factor,temp;
  VECTOR vec,posA;
  REAL total,counted;

  total=0.0;
  counted=0.0;
  equilibrium_distance=0.0;

  StartingBead=Components[CurrentComponent].StartingBead;
  typeA=Components[CurrentComponent].Type[StartingBead];
  well_depth_factor=Framework[CurrentSystem].SurfaceAreaProbeDistance;

  for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
    {
      typeB=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
      if(PseudoAtoms[typeB].Interaction)
      {
        equilibrium_distance=well_depth_factor*PotentialParms[typeA][typeB][1];

        total=0.0;
        counted=0.0;
        NumberOfTrials=Framework[CurrentSystem].SurfaceAreaSamplingPointsPerShere;
        for(k=0;k<NumberOfTrials;k++)
        {
          total+=1.0;

          vec=RandomNumberOnUnitSphere();
          posA.x=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.x+vec.x*equilibrium_distance;
          posA.y=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.y+vec.y*equilibrium_distance;
          posA.z=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.z+vec.z*equilibrium_distance;

          if(!BlockedPocket(posA))
          {
            Overlap=CheckSurfaceAreaOverlap(typeA,posA,CurrentFramework,i,-1,-1);
            if(!Overlap) counted+=1.0;
          }
        }
        temp=(counted/total)*4.0*M_PI*SQR(equilibrium_distance);
        SurfaceAreaFrameworkAverage[CurrentSystem][Block]+=temp;
        SurfaceAreaFrameworksAverage[CurrentSystem][CurrentFramework][Block]+=temp;
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeB=Cations[CurrentSystem][i].Atoms[j].Type;
      if(PseudoAtoms[typeB].Interaction)
      {
        equilibrium_distance=well_depth_factor*PotentialParms[typeA][typeB][1];

        total=0.0;
        counted=0.0;
        NumberOfTrials=Framework[CurrentSystem].SurfaceAreaSamplingPointsPerShere;
        for(k=0;k<NumberOfTrials;k++)
        {
          total+=1.0;

          vec=RandomNumberOnUnitSphere();
          posA.x=Cations[CurrentSystem][i].Atoms[j].Position.x+vec.x*equilibrium_distance;
          posA.y=Cations[CurrentSystem][i].Atoms[j].Position.y+vec.y*equilibrium_distance;
          posA.z=Cations[CurrentSystem][i].Atoms[j].Position.z+vec.z*equilibrium_distance;

          if((!BlockedPocket(posA))&&(ValidCartesianPoint(CurrentComponent,posA)))
          {
            Overlap=CheckSurfaceAreaOverlap(typeA,posA,-1,-1,i,j);
            if(!Overlap) counted+=1.0;
          }
        }
      }
    }
    temp=(counted/total)*4.0*M_PI*SQR(equilibrium_distance);
    SurfaceAreaFrameworkAverage[CurrentSystem][Block]+=temp;
    SurfaceAreaCationsAverage[CurrentSystem][Block]+=temp;
  }
  SurfaceAreaCount[CurrentSystem][Block]+=1.0;
}


enum {CF_INSERT_MOVE,CF_DELETE_MOVE,CF_MOVE};

void CFWangLandauIteration(int Switch)
{
  int i,j,k;
  int FractionalMolecule,index;
  int condition,MostProbableIndex;
  REAL Lambda,MostProbableValue,shift;

  switch(Switch)
  {
    case INITIALIZE:
      for(i=0;i<NumberOfComponents;i++)
        for(j=0;j<Components[i].CFLambdaHistogramSize;j++)
          CFLambdaHistogram[CurrentSystem][i][j]=0.0;
      break;
    case SAMPLE:
      FractionalMolecule=Components[CurrentComponent].FractionalMolecule[CurrentSystem];
      Lambda=Adsorbates[CurrentSystem][FractionalMolecule].Atoms[0].CFVDWScalingParameter;
      index=(int)(Components[CurrentComponent].CFLambdaHistogramSize*Lambda);
      if(index==Components[CurrentComponent].CFLambdaHistogramSize) index--;
      Components[CurrentComponent].CFBiasingFactors[CurrentSystem][index]-=Components[CurrentComponent].CFWangLandauScalingFactor[CurrentSystem];
      break;
    case PRINT:
      for(i=0;i<NumberOfSystems;i++)
      {
        for(j=0;j<NumberOfComponents;j++)
        {
          // find highest probability
          MostProbableIndex=0;
          MostProbableValue=CFLambdaHistogram[i][j][0];
          for(k=1;k<Components[j].CFLambdaHistogramSize;k++)
          {
            if(Components[j].CFBiasingFactors[i][k]>MostProbableValue)
            {
              MostProbableIndex=k;
              MostProbableValue=CFLambdaHistogram[i][j][k];
            }
          }

          // check if all bins are at least 70% of this value
          condition=TRUE;
          for(k=0;k<Components[j].CFLambdaHistogramSize;k++)
          {
            if(CFLambdaHistogram[i][j][k]<0.8*MostProbableValue)
            {
              condition=FALSE;
              break;
            }
          }

          if(condition)
          {
            Components[j].CFWangLandauScalingFactor[CurrentSystem]*=0.5;
            for(k=0;k<Components[j].CFLambdaHistogramSize;k++)
             CFLambdaHistogram[i][j][k]=0.0;
          }
        }
      }
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
        for(j=0;j<NumberOfComponents;j++)
        {
          shift=Components[j].CFBiasingFactors[i][0];
          for(k=0;k<Components[j].CFLambdaHistogramSize;k++)
          {
            CFLambdaHistogram[i][j][k]=0.0;
            Components[j].CFBiasingFactors[i][k]-=shift;
          }
        }
      break;
  }
}


/*********************************************************************************************************
 * Name       | CFSwapLambaMove                                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Grand-canonical Continuous Fraction Monte Carlo move of an adsorbate molecule            *
 * Parameters | -                                                                                        *
 * Note       | The Continous Fraction move changes the lambda-scaling parameter of a fractional         *
 *            | molecule.                                                                                *
 * Refs.      | W. Shi and E.J. Maginn, "Continuous Fractional Component Monte Carlo: An Adaptive        *
 *            | Biasing Method for Open System Atomistic Simulations", Journal of Chemical Theory and    *
 *            | Computation, 3, 1451-1463 (2007)                                                         *
 *********************************************************************************************************/

int CFSwapLambaMove(void)
{
  int i;
  REAL vNew;
  int MoveType;
  REAL RosenbluthNew,RosenbluthOld;
  int FractionalMolecule,index_old,index_new;
  REAL DeltaU,DeltaUFirstStep;
  REAL UAdsorbateVDWDeltaFirstStep,UHostVDWDeltaFirstStep,UCationVDWDeltaFirstStep;
  REAL UHostPolarizationNewFirstStep,UAdsorbatePolarizationNewFirstStep;
  REAL UPolarizationNewFirstStep,UCationPolarizationNewFirstStep;
  REAL UHostBackPolarizationNewFirstStep,UAdsorbateBackPolarizationNewFirstStep;
  REAL UBackPolarizationNewFirstStep,UCationBackPolarizationNewFirstStep;
  REAL UAdsorbateChargeChargeRealDeltaFirstStep,UAdsorbateBondDipoleBondDipoleRealDeltaFirstStep;
  REAL UAdsorbateAdsorbateChargeChargeFourierDeltaFirstStep,UAdsorbateChargeBondDipoleRealDeltaFirstStep;
  REAL UAdsorbateAdsorbateChargeBondDipoleFourierDeltaFirstStep,UAdsorbateAdsorbateBondDipoleBondDipoleFourierDeltaFirstStep;
  REAL UHostChargeChargeRealDeltaFirstStep,UHostChargeBondDipoleRealDeltaFirstStep;
  REAL UHostBondDipoleBondDipoleRealDeltaFirstStep,UHostAdsorbateChargeChargeFourierDeltaFirstStep;
  REAL UHostAdsorbateChargeBondDipoleFourierDeltaFirstStep,UHostAdsorbateBondDipoleBondDipoleFourierDeltaFirstStep;
  REAL UCationChargeChargeRealDeltaFirstStep,UCationChargeBondDipoleRealDeltaFirstStep;
  REAL UCationBondDipoleBondDipoleRealDeltaFirstStep,UAdsorbateCationChargeChargeFourierDeltaFirstStep;
  REAL UAdsorbateCationChargeBondDipoleFourierDeltaFirstStep,UAdsorbateCationBondDipoleBondDipoleFourierDeltaFirstStep;
  REAL UAdsorbateBondFirstStep,UAdsorbateUreyBradleyFirstStep,UAdsorbateBendFirstStep,UAdsorbateBendBendFirstStep;
  REAL UAdsorbateInversionBendFirstStep,UAdsorbateTorsionFirstStep,UAdsorbateImproperTorsionFirstStep;
  REAL UAdsorbateBondBondFirstStep,UAdsorbateBondBendFirstStep,UAdsorbateBondTorsionFirstStep;
  REAL UAdsorbateBendTorsionFirstStep,UAdsorbateIntraVDWFirstStep;
  REAL UAdsorbateIntraChargeChargeFirstStep,UAdsorbateIntraChargeBondDipoleFirstStep,UAdsorbateIntraBondDipoleBondDipoleFirstStep;
  int SelectedRetraceMolecule;
  REAL PartialFugacity,LambdaNew,LambdaOld,BiasNew,BiasOld,UTailDelta;

  FractionalMolecule=Components[CurrentComponent].FractionalMolecule[CurrentSystem];
  CurrentAdsorbateMolecule=FractionalMolecule;

  LambdaOld=Adsorbates[CurrentSystem][FractionalMolecule].Atoms[0].CFVDWScalingParameter;
  index_old=(int)(Components[CurrentComponent].CFLambdaHistogramSize*LambdaOld);
  if(index_old==Components[CurrentComponent].CFLambdaHistogramSize) index_old--;
  BiasOld=Components[CurrentComponent].CFBiasingFactors[CurrentSystem][index_old];


  CFLambdaHistogram[CurrentSystem][CurrentComponent][index_old]+=1.0;

  // determine a random change in Lambda
  vNew=(2.0*RandomNumber()-1.0)*MaximumCFLambdaChange[CurrentSystem][CurrentComponent];

  LambdaNew=LambdaOld+vNew;
  if(LambdaNew>1.0) 
  {
    MoveType=CF_INSERT_MOVE;
    LambdaNew-=1.0;
  }
  else if(LambdaNew<0.0) 
  {
    MoveType=CF_DELETE_MOVE;
    LambdaNew+=1.0;
  }
  else
    MoveType=CF_MOVE;
  index_new=(int)(Components[CurrentComponent].CFLambdaHistogramSize*LambdaNew);
  if(index_new==Components[CurrentComponent].CFLambdaHistogramSize) index_new--;
  BiasNew=Components[CurrentComponent].CFBiasingFactors[CurrentSystem][index_new];

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    CFVDWScalingStored[i]=Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFVDWScalingParameter;
    CFChargeScalingStored[i]=Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFChargeScalingParameter;
  }

  switch(MoveType)
  {
    case CF_INSERT_MOVE:
      CFSwapLambdaAttempts[CurrentSystem][CurrentComponent][1]++;
      // compute the energy difference for making the fractional molecule an integer molecule
      // (in the second step a new molecule will be grown containing the remainder of lambda)
      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
      {
        CFVDWScaling[i]=1.0;
        CFChargeScaling[i]=1.0;
      }
      break;
    case CF_DELETE_MOVE:
      CFSwapLambdaAttempts[CurrentSystem][CurrentComponent][2]++;
      if(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]<=1) return 0;

      // compute the energy difference for removing the fractional molecule
      // (in the second step a new fractional molecule will be retraced containing the remainder of lambda)
      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
      {
        CFVDWScaling[i]=0.0;
        CFChargeScaling[i]=0.0;
      }
      break;
    default:
      CFSwapLambdaAttempts[CurrentSystem][CurrentComponent][0]++;
      // compute the energy difference for the change in lambda of the fractional molecule
      // (there is no second step for this case)
      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
      {
        CFVDWScaling[i]=LambdaNew;
        CFChargeScaling[i]=pow(LambdaNew,5);
      }
      break;
  }

  // calculate the energy of the current configuration with a change in lambda (but with the same positions)
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    TrialPosition[CurrentSystem][i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;

  // calculate anisotropic sites
  CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);

  // compute inter-molecular energy differences
  CalculateInterVDWEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) goto l1;

  CalculateInterChargeChargeEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) goto l1;

  CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) goto l1;

  CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) goto l1;

  // compute energy differences framework-adsorbate
  CalculateFrameworkAdsorbateVDWEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) goto l1;

  CalculateFrameworkAdsorbateChargeChargeEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) goto l1;

  CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) goto l1;

  CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);
  if(OVERLAP) goto l1;

  // compute the energy differenes in Fourier-space
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    CalculateEwaldFourierAdsorbate(TRUE,TRUE,CurrentAdsorbateMolecule,0);

  if(ComputePolarization)
  { 
    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);

    UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
          (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
          UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
          UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
          (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
          UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
  }

  DeltaU=UHostVDWDelta[CurrentSystem]+UAdsorbateVDWDelta[CurrentSystem]+UCationVDWDelta[CurrentSystem]+
         UHostChargeChargeRealDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
         UCationChargeChargeRealDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
         UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
         UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+
         UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
         UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
         UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
         UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
         UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
         UDeltaPolarization;

  // save energy-changes of the first step

  UAdsorbateVDWDeltaFirstStep=UAdsorbateVDWDelta[CurrentSystem];
  UHostVDWDeltaFirstStep=UHostVDWDelta[CurrentSystem];
  UCationVDWDeltaFirstStep=UCationVDWDelta[CurrentSystem];

  UHostPolarizationNewFirstStep=UHostPolarizationNew[CurrentSystem];
  UAdsorbatePolarizationNewFirstStep=UAdsorbatePolarizationNew[CurrentSystem];
  UPolarizationNewFirstStep=UPolarizationNew[CurrentSystem];
  UCationPolarizationNewFirstStep=UCationPolarizationNew[CurrentSystem];
  UHostBackPolarizationNewFirstStep=UHostBackPolarizationNew[CurrentSystem];
  UAdsorbateBackPolarizationNewFirstStep=UAdsorbateBackPolarizationNew[CurrentSystem];
  UBackPolarizationNewFirstStep=UBackPolarizationNew[CurrentSystem];
  UCationBackPolarizationNewFirstStep=UCationBackPolarizationNew[CurrentSystem];

  UAdsorbateChargeChargeRealDeltaFirstStep=UAdsorbateChargeChargeRealDelta[CurrentSystem];
  UAdsorbateBondDipoleBondDipoleRealDeltaFirstStep=UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
  UAdsorbateAdsorbateChargeChargeFourierDeltaFirstStep=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
  UAdsorbateChargeBondDipoleRealDeltaFirstStep=UAdsorbateChargeBondDipoleRealDelta[CurrentSystem];
  UAdsorbateAdsorbateChargeBondDipoleFourierDeltaFirstStep=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
  UAdsorbateAdsorbateBondDipoleBondDipoleFourierDeltaFirstStep=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

  UHostChargeChargeRealDeltaFirstStep=UHostChargeChargeRealDelta[CurrentSystem];
  UHostChargeBondDipoleRealDeltaFirstStep=UHostChargeBondDipoleRealDelta[CurrentSystem];
  UHostBondDipoleBondDipoleRealDeltaFirstStep=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
  UHostAdsorbateChargeChargeFourierDeltaFirstStep=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
  UHostAdsorbateChargeBondDipoleFourierDeltaFirstStep=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
  UHostAdsorbateBondDipoleBondDipoleFourierDeltaFirstStep=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

  UCationChargeChargeRealDeltaFirstStep=UCationChargeChargeRealDelta[CurrentSystem];
  UCationChargeBondDipoleRealDeltaFirstStep=UCationChargeBondDipoleRealDelta[CurrentSystem];
  UCationBondDipoleBondDipoleRealDeltaFirstStep=UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
  UAdsorbateCationChargeChargeFourierDeltaFirstStep=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
  UAdsorbateCationChargeBondDipoleFourierDeltaFirstStep=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
  UAdsorbateCationBondDipoleBondDipoleFourierDeltaFirstStep=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

  DeltaUFirstStep=DeltaU;


  // get partial pressure for this component
  PartialFugacity=Components[CurrentComponent].FugacityCoefficient[CurrentSystem]*
                  Components[CurrentComponent].PartialPressure[CurrentSystem];

  RosenbluthNew=1.0;
  RosenbluthOld=1.0;
  switch(MoveType)
  {
    case CF_INSERT_MOVE:
      // a new molecule will be grown containing the remainder of lambda
      // the current fractional molecule is fully present

      GrowReservoirMolecule();

      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
      {
        CFVDWScaling[i]=LambdaNew;
        CFChargeScaling[i]=pow(LambdaNew,5);

        Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFVDWScalingParameter=1.0;
        Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFChargeScalingParameter=1.0;
      }

      // move blocked pocket to cbmc.c!
      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
      {
        if(BlockedPocket(NewPosition[CurrentSystem][i]))
          return 0;
      }

      // calculate the energy of the current configuration with a change in lambda (but with the same positions)
      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
        TrialPosition[CurrentSystem][i]=NewPosition[CurrentSystem][i];

      // calculate anisotropic sites
      CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);

      // compute inter-molecular energy differences
      CalculateInterVDWEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,FALSE);
      if(OVERLAP) goto l1;

      CalculateInterChargeChargeEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,FALSE);
      if(OVERLAP) goto l1;

      CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,FALSE);
      if(OVERLAP) goto l1;

      CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,FALSE);
      if(OVERLAP) goto l1;

      // compute energy differences framework-adsorbate
      CalculateFrameworkAdsorbateVDWEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,FALSE);
      if(OVERLAP) goto l1;

      CalculateFrameworkAdsorbateChargeChargeEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,FALSE);
      if(OVERLAP) goto l1;

      CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,FALSE);
      if(OVERLAP) goto l1;

      CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,FALSE);
      if(OVERLAP) goto l1;

      // compute the energy differenes in Fourier-space
      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        CalculateEwaldFourierAdsorbate2(TRUE,FALSE,-1,0); // HERE

      UTailDelta=TailMolecularEnergyDifference(CurrentComponent,CurrentComponent,TRUE,FALSE);

      if(ComputePolarization)
      { 
        ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);

        UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
              (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
              UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
              UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
              (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
              UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
      }

      DeltaU=UHostVDWDelta[CurrentSystem]+UAdsorbateVDWDelta[CurrentSystem]+UCationVDWDelta[CurrentSystem]+
             UHostChargeChargeRealDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
             UCationChargeChargeRealDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
             UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
             UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+
             UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
             UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
             UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
             UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
             UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
             UDeltaPolarization+UTailDelta;

      RosenbluthNew=exp(-Beta[CurrentSystem]*DeltaU)*Beta[CurrentSystem]*PartialFugacity*Volume[CurrentSystem]/
                    (Components[CurrentComponent].NumberOfMolecules[CurrentSystem]);
      break;
    case CF_DELETE_MOVE:
      // a new fractional molecule will be retraced containing the remainder of lambda
      // the current fractional molecule is fully removed

      // choose a random molecule of this component (but not the current fractional molecule)
      CurrentAdsorbateMolecule=SelectRandomMoleculeOfTypeExcludingFractionalMolecule(CurrentComponent);

      SelectedRetraceMolecule=CurrentAdsorbateMolecule;

      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
      {
        // current fractional particle is removed
        Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFVDWScalingParameter=0.0;
        Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFChargeScalingParameter=0.0;

        CFVDWScalingStored2[i]=Adsorbates[CurrentSystem][SelectedRetraceMolecule].Atoms[i].CFVDWScalingParameter;
        CFChargeScalingStored2[i]=Adsorbates[CurrentSystem][SelectedRetraceMolecule].Atoms[i].CFChargeScalingParameter;
      }

      // a newly selected fractional particle will get these lambdas
      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
      {
        Adsorbates[CurrentSystem][SelectedRetraceMolecule].Atoms[i].CFVDWScalingParameter=1.0;
        Adsorbates[CurrentSystem][SelectedRetraceMolecule].Atoms[i].CFChargeScalingParameter=1.0;

        CFVDWScaling[i]=LambdaNew;
        CFChargeScaling[i]=pow(LambdaNew,5);
      }

      // calculate the energy of the current configuration with a change in lambda (but with the same positions)
      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
        TrialPosition[CurrentSystem][i]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[i].Position;

      // calculate anisotropic sites
      CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);

      // compute inter-molecular energy differences
      CalculateInterVDWEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);

      CalculateInterChargeChargeEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);

      CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);

      CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);

      // compute energy differences framework-adsorbate
      CalculateFrameworkAdsorbateVDWEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);

      CalculateFrameworkAdsorbateChargeChargeEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);

      CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);

      CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference(CurrentAdsorbateMolecule,CurrentComponent,TRUE,TRUE);

      // compute the energy differenes in Fourier-space
      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        CalculateEwaldFourierAdsorbate2(TRUE,TRUE,CurrentAdsorbateMolecule,0);

      UTailDelta=TailMolecularEnergyDifference(CurrentComponent,CurrentComponent,FALSE,TRUE);

      if(ComputePolarization)
      { 
        ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);

        UDeltaPolarization=UHostPolarizationNew[CurrentSystem]-UHostPolarization[CurrentSystem]+
              (UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem])-UAdsorbatePolarization[CurrentSystem]+
              UCationPolarizationNew[CurrentSystem]-UCationPolarization[CurrentSystem]+
              UHostBackPolarizationNew[CurrentSystem]-UHostBackPolarization[CurrentSystem]+
              (UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem])-UAdsorbateBackPolarization[CurrentSystem]+
              UCationBackPolarizationNew[CurrentSystem]-UCationBackPolarization[CurrentSystem];
      }

      DeltaU=UHostVDWDelta[CurrentSystem]+UAdsorbateVDWDelta[CurrentSystem]+UCationVDWDelta[CurrentSystem]+
             UHostChargeChargeRealDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
             UCationChargeChargeRealDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
             UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
             UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+
             UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
             UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
             UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
             UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+
             UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
             UDeltaPolarization+UTailDelta;

      for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
      {
        Adsorbates[CurrentSystem][SelectedRetraceMolecule].Atoms[i].CFVDWScalingParameter=CFVDWScaling[i];
        Adsorbates[CurrentSystem][SelectedRetraceMolecule].Atoms[i].CFChargeScalingParameter=CFChargeScaling[i];
      }
      RosenbluthOld=(REAL)(Components[CurrentComponent].NumberOfMolecules[CurrentSystem]-1)/
                    (exp(Beta[CurrentSystem]*DeltaU)*Beta[CurrentSystem]*PartialFugacity*Volume[CurrentSystem]);
      break;
    default:
      break;
  }

  if(RandomNumber()<exp(-Beta[CurrentSystem]*DeltaUFirstStep)*exp(BiasNew-BiasOld)*RosenbluthNew*RosenbluthOld)
  {
    #ifdef DEBUG
      printf("Lambda-move adsorbate accepted\n");
    #endif

    UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWDeltaFirstStep;
    UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWDeltaFirstStep;
    UHostAdsorbate[CurrentSystem]+=UHostVDWDeltaFirstStep;
    UHostAdsorbateVDW[CurrentSystem]+=UHostVDWDeltaFirstStep;
    UAdsorbateCation[CurrentSystem]+=UCationVDWDeltaFirstStep;
    UAdsorbateCationVDW[CurrentSystem]+=UCationVDWDeltaFirstStep;

    UHostPolarization[CurrentSystem]=UHostPolarizationNewFirstStep;
    UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNewFirstStep+UPolarizationNewFirstStep;
    UCationPolarization[CurrentSystem]=UCationPolarizationNewFirstStep;

    UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNewFirstStep;
    UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNewFirstStep+UBackPolarizationNewFirstStep;
    UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNewFirstStep;

    if(ChargeMethod!=NONE)
    {
      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeRealDeltaFirstStep;
      UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleRealDeltaFirstStep;
      UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDeltaFirstStep;
      UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleRealDeltaFirstStep;
      UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDeltaFirstStep;
      UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDeltaFirstStep;
      UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateChargeChargeRealDeltaFirstStep+UAdsorbateAdsorbateChargeChargeFourierDeltaFirstStep+
                                                 UAdsorbateChargeBondDipoleRealDeltaFirstStep+UAdsorbateAdsorbateChargeBondDipoleFourierDeltaFirstStep+
                                                 UAdsorbateBondDipoleBondDipoleRealDeltaFirstStep+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDeltaFirstStep;
      UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDeltaFirstStep+UAdsorbateChargeChargeRealDeltaFirstStep+
                                          UAdsorbateChargeBondDipoleRealDeltaFirstStep+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                          UAdsorbateBondDipoleBondDipoleRealDeltaFirstStep+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDeltaFirstStep;



      UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDeltaFirstStep;
      UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDeltaFirstStep;
      UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDeltaFirstStep;
      UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDeltaFirstStep;
      UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDeltaFirstStep;
      UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDeltaFirstStep;
      UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDeltaFirstStep+UHostChargeChargeRealDeltaFirstStep+
                                            UHostAdsorbateChargeBondDipoleFourierDeltaFirstStep+UHostChargeBondDipoleRealDeltaFirstStep+
                                            UHostAdsorbateBondDipoleBondDipoleFourierDeltaFirstStep+UHostBondDipoleBondDipoleRealDeltaFirstStep;
      UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDeltaFirstStep+UHostChargeChargeRealDeltaFirstStep+
                                     UHostAdsorbateChargeBondDipoleFourierDeltaFirstStep+UHostChargeBondDipoleRealDeltaFirstStep+
                                     UHostAdsorbateBondDipoleBondDipoleFourierDeltaFirstStep+UHostBondDipoleBondDipoleRealDeltaFirstStep;



      UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeRealDeltaFirstStep;
      UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleRealDeltaFirstStep;
      UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleRealDeltaFirstStep;
      UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDeltaFirstStep;
      UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDeltaFirstStep;
      UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDeltaFirstStep;
      UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDeltaFirstStep+UCationChargeChargeRealDeltaFirstStep+
                                              UAdsorbateCationChargeBondDipoleFourierDeltaFirstStep+UCationChargeBondDipoleRealDeltaFirstStep+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDeltaFirstStep+UCationBondDipoleBondDipoleRealDeltaFirstStep;
      UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDeltaFirstStep+UCationChargeChargeRealDeltaFirstStep+
                                              UAdsorbateCationChargeBondDipoleFourierDeltaFirstStep+UCationChargeBondDipoleRealDeltaFirstStep+
                                              UAdsorbateCationBondDipoleBondDipoleFourierDeltaFirstStep+UCationBondDipoleBondDipoleRealDeltaFirstStep;

      if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
        AcceptEwaldAdsorbateMove(0);
    }


    UTotal[CurrentSystem]+=DeltaUFirstStep;

    switch(MoveType)
    {
      case CF_INSERT_MOVE:
        CFSwapLambdaAccepted[CurrentSystem][CurrentComponent][1]++;

        InsertAdsorbateMolecule();

        FractionalMolecule=NumberOfAdsorbateMolecules[CurrentSystem]-1;
        Components[CurrentComponent].FractionalMolecule[CurrentSystem]=FractionalMolecule;

        UAdsorbateBondFirstStep=CalculateBondEnergyAdsorbate(FractionalMolecule);
        UAdsorbateUreyBradleyFirstStep=CalculateUreyBradleyEnergyAdsorbate(FractionalMolecule);
        UAdsorbateBendFirstStep=CalculateBendEnergyAdsorbate(FractionalMolecule);
        UAdsorbateBendBendFirstStep=CalculateBendBendEnergyAdsorbate(FractionalMolecule);
        UAdsorbateInversionBendFirstStep=CalculateInversionBendEnergyAdsorbate(FractionalMolecule);
        UAdsorbateTorsionFirstStep=CalculateTorsionEnergyAdsorbate(FractionalMolecule);
        UAdsorbateImproperTorsionFirstStep=CalculateImproperTorsionEnergyAdsorbate(FractionalMolecule);
        UAdsorbateBondBondFirstStep=CalculateBondBondEnergyAdsorbate(FractionalMolecule);
        UAdsorbateBondBendFirstStep=CalculateBondBendEnergyAdsorbate(FractionalMolecule);
        UAdsorbateBondTorsionFirstStep=CalculateBondTorsionEnergyAdsorbate(FractionalMolecule);
        UAdsorbateBendTorsionFirstStep=CalculateBendTorsionEnergyAdsorbate(FractionalMolecule);
        UAdsorbateIntraVDWFirstStep=CalculateIntraVDWEnergyAdsorbate(FractionalMolecule);

        if(ChargeMethod!=NONE)
        {
          UAdsorbateIntraChargeChargeFirstStep=CalculateIntraChargeChargeEnergyAdsorbate(FractionalMolecule);
          UAdsorbateIntraChargeBondDipoleFirstStep=CalculateIntraChargeBondDipoleEnergyAdsorbate(FractionalMolecule);
          UAdsorbateIntraBondDipoleBondDipoleFirstStep=CalculateIntraBondDipoleBondDipoleEnergyAdsorbate(FractionalMolecule);
        }

        UAdsorbateBond[CurrentSystem]+=UAdsorbateBondFirstStep;
        UAdsorbateUreyBradley[CurrentSystem]+=UAdsorbateUreyBradleyFirstStep;
        UAdsorbateBend[CurrentSystem]+=UAdsorbateBendFirstStep;
        UAdsorbateBendBend[CurrentSystem]+=UAdsorbateBendBendFirstStep;
        UAdsorbateInversionBend[CurrentSystem]+=UAdsorbateInversionBendFirstStep;
        UAdsorbateTorsion[CurrentSystem]+=UAdsorbateTorsionFirstStep;
        UAdsorbateImproperTorsion[CurrentSystem]+=UAdsorbateImproperTorsionFirstStep;
        UAdsorbateBondBond[CurrentSystem]+=UAdsorbateBondBondFirstStep;
        UAdsorbateBondBend[CurrentSystem]+=UAdsorbateBondBendFirstStep;
        UAdsorbateBondTorsion[CurrentSystem]+=UAdsorbateBondTorsionFirstStep;
        UAdsorbateBendTorsion[CurrentSystem]+=UAdsorbateBendTorsionFirstStep;
        UAdsorbateIntraVDW[CurrentSystem]+=UAdsorbateIntraVDWFirstStep;
        UTotal[CurrentSystem]+=UAdsorbateBondFirstStep+UAdsorbateUreyBradleyFirstStep+UAdsorbateBendFirstStep+UAdsorbateBendBendFirstStep+
              UAdsorbateInversionBendFirstStep+UAdsorbateTorsionFirstStep+UAdsorbateImproperTorsionFirstStep+UAdsorbateBondBondFirstStep+
              UAdsorbateBondBendFirstStep+UAdsorbateBondTorsionFirstStep+UAdsorbateBendTorsionFirstStep+UAdsorbateIntraVDWFirstStep;

        if(ChargeMethod!=NONE)
        {
          UAdsorbateIntraChargeCharge[CurrentSystem]+=UAdsorbateIntraChargeChargeFirstStep;
          UAdsorbateIntraChargeBondDipole[CurrentSystem]+=UAdsorbateIntraChargeBondDipoleFirstStep;
          UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=UAdsorbateIntraBondDipoleBondDipoleFirstStep;
          UTotal[CurrentSystem]+=UAdsorbateIntraChargeChargeFirstStep+UAdsorbateIntraChargeBondDipoleFirstStep+UAdsorbateIntraBondDipoleBondDipoleFirstStep;
        }



        UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
        UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
        UHostAdsorbate[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
        UHostAdsorbateVDW[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
        UAdsorbateCation[CurrentSystem]+=UCationVDWDelta[CurrentSystem];
        UAdsorbateCationVDW[CurrentSystem]+=UCationVDWDelta[CurrentSystem];

        UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
        UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
        UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

        UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
        UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
        UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

        if(ChargeMethod!=NONE)
        {
          UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem];
          UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
          UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
          UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleRealDelta[CurrentSystem];
          UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
          UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
          UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                     UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                     UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
          UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];



          UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem];
          UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDelta[CurrentSystem];
          UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
          UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
          UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
          UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
          UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                                UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                                UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
          UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                         UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                         UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];



          UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem];
          UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleRealDelta[CurrentSystem];
          UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
          UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
          UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
          UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
          UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                                  UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
                                                  UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
          UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                                  UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
                                                  UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem];

        }

        UTailCorrection[CurrentSystem]+=UTailDelta;

        UTotal[CurrentSystem]+=DeltaU;
        break;
      case CF_DELETE_MOVE:
        CFSwapLambdaAccepted[CurrentSystem][CurrentComponent][2]++;

        // remove old fractional molecule
        CurrentAdsorbateMolecule=FractionalMolecule;

        // compute the internal energy of the molecule that will be deleted
        UAdsorbateBondFirstStep=CalculateBondEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateUreyBradleyFirstStep=CalculateUreyBradleyEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateBendFirstStep=CalculateBendEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateBendBendFirstStep=CalculateBendBendEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateInversionBendFirstStep=CalculateInversionBendEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateTorsionFirstStep=CalculateTorsionEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateImproperTorsionFirstStep=CalculateImproperTorsionEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateBondBondFirstStep=CalculateBondBondEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateBondBendFirstStep=CalculateBondBendEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateBondTorsionFirstStep=CalculateBondTorsionEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateBendTorsionFirstStep=CalculateBendTorsionEnergyAdsorbate(CurrentAdsorbateMolecule);
        UAdsorbateIntraVDWFirstStep=CalculateIntraVDWEnergyAdsorbate(CurrentAdsorbateMolecule);

        if(ChargeMethod!=NONE)
        {
          UAdsorbateIntraChargeChargeFirstStep=CalculateIntraChargeChargeEnergyAdsorbate(CurrentAdsorbateMolecule);
          UAdsorbateIntraChargeBondDipoleFirstStep=CalculateIntraChargeBondDipoleEnergyAdsorbate(CurrentAdsorbateMolecule);
          UAdsorbateIntraBondDipoleBondDipoleFirstStep=CalculateIntraBondDipoleBondDipoleEnergyAdsorbate(CurrentAdsorbateMolecule);
        }


        UAdsorbateBond[CurrentSystem]-=UAdsorbateBondFirstStep;
        UAdsorbateUreyBradley[CurrentSystem]-=UAdsorbateUreyBradleyFirstStep;
        UAdsorbateBend[CurrentSystem]-=UAdsorbateBendFirstStep;
        UAdsorbateBendBend[CurrentSystem]-=UAdsorbateBendBendFirstStep;
        UAdsorbateInversionBend[CurrentSystem]-=UAdsorbateInversionBendFirstStep;
        UAdsorbateTorsion[CurrentSystem]-=UAdsorbateTorsionFirstStep;
        UAdsorbateImproperTorsion[CurrentSystem]-=UAdsorbateImproperTorsionFirstStep;
        UAdsorbateBondBond[CurrentSystem]-=UAdsorbateBondBondFirstStep;
        UAdsorbateBondBend[CurrentSystem]-=UAdsorbateBondBendFirstStep;
        UAdsorbateBondTorsion[CurrentSystem]-=UAdsorbateBondTorsionFirstStep;
        UAdsorbateBendTorsion[CurrentSystem]-=UAdsorbateBendTorsionFirstStep;
        UAdsorbateIntraVDW[CurrentSystem]-=UAdsorbateIntraVDWFirstStep;
        UTotal[CurrentSystem]-=UAdsorbateBondFirstStep+UAdsorbateUreyBradleyFirstStep+UAdsorbateBendFirstStep+UAdsorbateBendBendFirstStep+
              UAdsorbateInversionBendFirstStep+UAdsorbateTorsionFirstStep+UAdsorbateImproperTorsionFirstStep+UAdsorbateBondBondFirstStep+
              UAdsorbateBondBendFirstStep+UAdsorbateBondTorsionFirstStep+UAdsorbateBendTorsionFirstStep+UAdsorbateIntraVDWFirstStep;

        if(ChargeMethod!=NONE)
        {
          UAdsorbateIntraChargeCharge[CurrentSystem]-=UAdsorbateIntraChargeChargeFirstStep;
          UAdsorbateIntraChargeBondDipole[CurrentSystem]-=UAdsorbateIntraChargeBondDipoleFirstStep;
          UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]-=UAdsorbateIntraBondDipoleBondDipoleFirstStep;
          UTotal[CurrentSystem]-=UAdsorbateIntraChargeChargeFirstStep+UAdsorbateIntraChargeBondDipoleFirstStep+UAdsorbateIntraBondDipoleBondDipoleFirstStep;
        }

        Components[CurrentComponent].FractionalMolecule[CurrentSystem]=SelectedRetraceMolecule;

        // delete the current fractional molecule
        RemoveAdsorbateMolecule();

        UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
        UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWDelta[CurrentSystem];
        UHostAdsorbate[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
        UHostAdsorbateVDW[CurrentSystem]+=UHostVDWDelta[CurrentSystem];
        UAdsorbateCation[CurrentSystem]+=UCationVDWDelta[CurrentSystem];
        UAdsorbateCationVDW[CurrentSystem]+=UCationVDWDelta[CurrentSystem];

        UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
        UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
        UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

        UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
        UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
        UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

        if(ChargeMethod!=NONE)
        {
          UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem];
          UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem];
          UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
          UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleRealDelta[CurrentSystem];
          UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
          UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
          UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateChargeChargeRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
                                                     UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                                     UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
          UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateChargeChargeRealDelta[CurrentSystem]+
                                              UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
                                              UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];



          UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeRealDelta[CurrentSystem];
          UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleRealDelta[CurrentSystem];
          UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
          UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
          UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
          UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
          UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                                UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                                UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];
          UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UHostChargeChargeRealDelta[CurrentSystem]+
                                         UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UHostChargeBondDipoleRealDelta[CurrentSystem]+
                                         UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UHostBondDipoleBondDipoleRealDelta[CurrentSystem];



          UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeRealDelta[CurrentSystem];
          UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleRealDelta[CurrentSystem];
          UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
          UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
          UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
          UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
          UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                                  UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
                                                  UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem];
          UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+UCationChargeChargeRealDelta[CurrentSystem]+
                                                  UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+UCationChargeBondDipoleRealDelta[CurrentSystem]+
                                                  UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+UCationBondDipoleBondDipoleRealDelta[CurrentSystem];

        }

        UTailCorrection[CurrentSystem]+=UTailDelta;

        UTotal[CurrentSystem]+=DeltaU;
        break;
      default:
        CFSwapLambdaAccepted[CurrentSystem][CurrentComponent][0]++;
        for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
        {
          Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFVDWScalingParameter=CFVDWScaling[i];
          Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFChargeScalingParameter=CFChargeScaling[i];
        }
        break;
    }

  }
  else
  {
l1:;
    #ifdef DEBUG
      printf("Lambda-move adsorbate rejected\n");
    #endif
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFVDWScalingParameter=CFVDWScalingStored[i];
      Adsorbates[CurrentSystem][FractionalMolecule].Atoms[i].CFChargeScalingParameter=CFChargeScalingStored[i];
      if(MoveType==CF_DELETE_MOVE)
      {
        Adsorbates[CurrentSystem][SelectedRetraceMolecule].Atoms[i].CFVDWScalingParameter=CFVDWScalingStored2[i];
        Adsorbates[CurrentSystem][SelectedRetraceMolecule].Atoms[i].CFChargeScalingParameter=CFChargeScalingStored2[i];
      }
    }
  }

  return 0;

}

/*********************************************************************************************************
 * Name       | OptimizeTranslationAcceptence                                                            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Adjusts the maximum displacement on-the-fly to obtain an acceptance-rate of 50%.         *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void OptimizeCFLambdaChangeAcceptence(void)
{
  int i;
  REAL ratio,vandr;

  for(i=0;i<NumberOfComponents;i++)
  {
    if(CFSwapLambdaAttempts[CurrentSystem][i][0]>0.0)
      ratio=CFSwapLambdaAccepted[CurrentSystem][i][0]/CFSwapLambdaAttempts[CurrentSystem][i][0];
    else
      ratio=0.0;

    vandr=ratio/TargetAccRatioTranslation;
    if(vandr>1.0) vandr=1.0;
    else if(vandr<0.5) vandr=0.5;
    MaximumCFLambdaChange[CurrentSystem][i]*=vandr;
    if(MaximumCFLambdaChange[CurrentSystem][i]<0.01)
       MaximumCFLambdaChange[CurrentSystem][i]=0.01;
    if(MaximumCFLambdaChange[CurrentSystem][i]>1.0)
       MaximumCFLambdaChange[CurrentSystem][i]=1.0;

    TotalCFSwapLambdaAttempts[CurrentSystem][i]+=CFSwapLambdaAttempts[CurrentSystem][i][0];
    TotalCFSwapLambdaAccepted[CurrentSystem][i]+=CFSwapLambdaAccepted[CurrentSystem][i][0];
    CFSwapLambdaAttempts[CurrentSystem][i][0]=CFSwapLambdaAccepted[CurrentSystem][i][0]=0.0;
  }
}


void PrintCFSwapLambdaStatistics(FILE *FilePtr)
{
  int i,k,MoveUsed;
  REAL total;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfCFSwapLambdaMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the CF swap lambda move:\n");
    fprintf(FilePtr,"=======================================\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      if(Components[i].FractionOfCFSwapLambdaMove>0.0)
      {
        fprintf(FilePtr,"Component [%s] total tried: %lf constant-lambda accepted: %lf (%lf [%%])\n",
          Components[i].Name,
          (double)CFSwapLambdaAttempts[CurrentSystem][i][0],
          (double)CFSwapLambdaAccepted[CurrentSystem][i][0],
          (double)(CFSwapLambdaAttempts[CurrentSystem][i][0]>(REAL)0.0?
            100.0*CFSwapLambdaAccepted[CurrentSystem][i][0]/CFSwapLambdaAttempts[CurrentSystem][i][0]:(REAL)0.0));
        fprintf(FilePtr,"               total tried: %lf insert-lambda accepted: %lf (%lf [%%])\n",
          (double)CFSwapLambdaAttempts[CurrentSystem][i][1],
          (double)CFSwapLambdaAccepted[CurrentSystem][i][1],
          (double)(CFSwapLambdaAttempts[CurrentSystem][i][1]>(REAL)0.0?
            100.0*CFSwapLambdaAccepted[CurrentSystem][i][1]/CFSwapLambdaAttempts[CurrentSystem][i][1]:(REAL)0.0));
        fprintf(FilePtr,"               total tried: %lf remove-lambda accepted: %lf (%lf [%%])\n",
          (double)CFSwapLambdaAttempts[CurrentSystem][i][2],
          (double)CFSwapLambdaAccepted[CurrentSystem][i][2],
          (double)(CFSwapLambdaAttempts[CurrentSystem][i][2]>(REAL)0.0?
            100.0*CFSwapLambdaAccepted[CurrentSystem][i][2]/CFSwapLambdaAttempts[CurrentSystem][i][2]:(REAL)0.0));

        total=0.0;
        for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
          total+=CFLambdaHistogram[CurrentSystem][i][k];

        fprintf(FilePtr,"\n\tLambda probabilities:\n");
        fprintf(FilePtr,"\t---------------------\n");
        for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
          fprintf(FilePtr,"\tLambda [ %4f - %4f ]: %18.10f (biasing factor: %18.10f)\n",
            (REAL)k/Components[i].CFLambdaHistogramSize,(REAL)(k+1)/Components[i].CFLambdaHistogramSize,
            100.0*CFLambdaHistogram[CurrentSystem][i][k]/total,Components[i].CFBiasingFactors[CurrentSystem][k]);
      }
    }
    fprintf(FilePtr,"\n\n");
  }
  else
    fprintf(FilePtr,"CF swap lambda move was OFF for all components\n\n");
}

int CFCBSwapLambaMove(void)
{
  return 0;
}


void PrintCFCBSwapLambdaStatistics(FILE *FilePtr)
{
  int i,k,MoveUsed;
  REAL total;

  MoveUsed=FALSE;
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].FractionOfCFCBSwapLambdaMove>0.0)
    {
      MoveUsed=TRUE;
      break;
    }
  }

  if(MoveUsed)
  {
    fprintf(FilePtr,"Performance of the CFCB swap lambda move:\n");
    fprintf(FilePtr,"=========================================\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      if(Components[i].FractionOfCFCBSwapLambdaMove>0.0)
        fprintf(FilePtr,"Component [%s] total tried: %lf constant-lambda accepted: %lf (%lf [%%])\n\t\tinsert accepted %lf (%lf [%%])\n\t\tremoved accepted %lf (%lf [%%])\n",
          Components[i].Name,
          (double)CFCBSwapLambdaAttempts[CurrentSystem][i][0],
          (double)CFCBSwapLambdaAccepted[CurrentSystem][i][0],
          (double)(CFCBSwapLambdaAttempts[CurrentSystem][i][0]>(REAL)0.0?
            100.0*CFCBSwapLambdaAccepted[CurrentSystem][i][0]/CFCBSwapLambdaAttempts[CurrentSystem][i][0]:(REAL)0.0),
          (double)CFCBSwapLambdaAccepted[CurrentSystem][i][1],
          (double)(CFCBSwapLambdaAttempts[CurrentSystem][i][1]>(REAL)0.0?
            100.0*CFCBSwapLambdaAccepted[CurrentSystem][i][1]/CFCBSwapLambdaAttempts[CurrentSystem][i][1]:(REAL)0.0),
          (double)CFCBSwapLambdaAccepted[CurrentSystem][i][2],
          (double)(CFCBSwapLambdaAttempts[CurrentSystem][i][2]>(REAL)0.0?
            100.0*CFCBSwapLambdaAccepted[CurrentSystem][i][2]/CFCBSwapLambdaAttempts[CurrentSystem][i][2]:(REAL)0.0));

      total=0.0;
      for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
        total+=CFLambdaHistogram[CurrentSystem][i][k];

      fprintf(FilePtr,"\tLambda probabilities:\n"); 
      fprintf(FilePtr,"\t---------------------\n"); 
      for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
        fprintf(FilePtr,"\tLambda [%g-%g]: %18.10f\n",(REAL)k/Components[i].CFLambdaHistogramSize,
        (REAL)(k+1)/Components[i].CFLambdaHistogramSize,CFLambdaHistogram[CurrentSystem][i][k]);
        
    }
    fprintf(FilePtr,"\n\n");
  }
  else
    fprintf(FilePtr,"CFCB swap lambda move was OFF for all components\n\n");
}

int CFGibbsParticleTransferMove(void)
{
  return 0;
}

int CBCFGibbsParticleTransferMove(void)
{
  return 0;
}


void WriteRestartMcMoves(FILE *FilePtr)
{
  int i,j;
  REAL Check;

  for(i=0;i<NumberOfSystems;i++)
  {
    fwrite(TranslationAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(TranslationAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(TotalTranslationAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(TotalTranslationAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(MaximumTranslation[i],sizeof(VECTOR),NumberOfComponents,FilePtr);

    fwrite(RandomTranslationAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(RandomTranslationAccepted[i],sizeof(REAL),NumberOfComponents,FilePtr);

    fwrite(RotationAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(RotationAccepted[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(MaximumRotation[i],sizeof(REAL),NumberOfComponents,FilePtr);

    fwrite(SwapAddAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(SwapAddAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fwrite(SwapRemoveAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(SwapRemoveAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fwrite(CFSwapLambdaAttempts[i],sizeof(REAL[3]),NumberOfComponents,FilePtr);
    fwrite(CFSwapLambdaAccepted[i],sizeof(REAL[3]),NumberOfComponents,FilePtr);
    fwrite(TotalCFSwapLambdaAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(TotalCFSwapLambdaAccepted[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(MaximumCFLambdaChange[i],sizeof(REAL),NumberOfComponents,FilePtr);

    fwrite(CFCBSwapLambdaAttempts[i],sizeof(REAL[3]),NumberOfComponents,FilePtr);
    fwrite(CFCBSwapLambdaAccepted[i],sizeof(REAL[3]),NumberOfComponents,FilePtr);
    fwrite(TotalCFCBSwapLambdaAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(TotalCFCBSwapLambdaAccepted[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(MaximumCFCBLambdaChange[i],sizeof(REAL),NumberOfComponents,FilePtr);

    fwrite(ReinsertionAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(ReinsertionAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fwrite(PartialReinsertionAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(PartialReinsertionAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fwrite(ReinsertionInPlaceAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(ReinsertionInPlaceAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fwrite(TranslationInPlaceAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(TranslationInPlaceAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(TotalTranslationInPlaceAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(TotalTranslationInPlaceAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);

    fwrite(ReinsertionInPlaneAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fwrite(ReinsertionInPlaneAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fwrite(TranslationInPlaneAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(TranslationInPlaneAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(TotalTranslationInPlaneAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(TotalTranslationInPlaneAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fwrite(MaximumTranslationInPlane[i],sizeof(VECTOR),NumberOfComponents,FilePtr);

    for(j=0;j<NumberOfComponents;j++)
    {
      fwrite(IdentityChangeAttempts[i][j],sizeof(REAL),NumberOfComponents,FilePtr);
      fwrite(IdentityChangeAccepted[i][j],sizeof(REAL[2]),NumberOfComponents,FilePtr);

      fwrite(GibbsIdentityChangeAttempts[i][j],sizeof(REAL),NumberOfComponents,FilePtr);
      fwrite(GibbsIdentityChangeAccepted[i][j],sizeof(REAL[2]),NumberOfComponents,FilePtr);

      fwrite(CFLambdaHistogram[i][j],sizeof(REAL),Components[j].CFLambdaHistogramSize,FilePtr);
    }

    fwrite(ParallelTemperingAttempts[i],sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(ParallelTemperingAccepted[i],sizeof(REAL),NumberOfSystems,FilePtr);

    fwrite(HyperParallelTemperingAttempts[i],sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(HyperParallelTemperingAccepted[i],sizeof(REAL),NumberOfSystems,FilePtr);

    fwrite(ParallelMolFractionAttempts[i],sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(ParallelMolFractionAccepted[i],sizeof(REAL),NumberOfSystems,FilePtr);
  }
  fwrite(&ParallelMolFractionComponentA,sizeof(int),1,FilePtr);
  fwrite(&ParallelMolFractionComponentB,sizeof(int),1,FilePtr);

  fwrite(ChiralInversionAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(ChiralInversionAccepted,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(VolumeChangeAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(VolumeChangeAccepted,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(TotalVolumeChangeAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(TotalVolumeChangeAccepted,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(MaximumVolumeChange,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(BoxShapeChangeAttempts,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(BoxShapeChangeAccepted,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(MaximumBoxShapeChange,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);

  fwrite(GibbsVolumeChangeAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(GibbsVolumeChangeAccepted,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(MaximumGibbsVolumeChange,sizeof(REAL),NumberOfSystems,FilePtr);

  for(i=0;i<NumberOfComponents;i++)
  {
    fwrite(GibbsSwapAttempts[i],sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(GibbsSwapAccepted[i],sizeof(REAL),NumberOfSystems,FilePtr);
  }

  for(i=0;i<NumberOfSystems;i++)
  {
    fwrite(FrameworkShiftAttempts[i],sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
    fwrite(FrameworkShiftAccepted[i],sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
    fwrite(TotalFrameworkShiftAttempts[i],sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
    fwrite(TotalFrameworkShiftAccepted[i],sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
    fwrite(FrameworkMaximumShiftTranslation[i],sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
  }

  fwrite(&NumberOfHybridNVESteps,sizeof(int),1,FilePtr);
  fwrite(HybridNVEDrift,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEDriftCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEAccepted,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(HybridNVEStartTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEStartTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(HybridNVEEndTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNVEEndTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);


  fwrite(&NumberOfHybridNPHSteps,sizeof(int),1,FilePtr);
  fwrite(HybridNPHDrift,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHDriftCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHAccepted,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(HybridNPHStartTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartCellTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHStartTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(HybridNPHEndTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndCellTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHEndTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(&NumberOfHybridNPHPRSteps,sizeof(int),1,FilePtr);
  fwrite(HybridNPHPRDrift,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRDriftCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRAccepted,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(HybridNPHPRStartTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartCellTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPRStartTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(HybridNPHPREndTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndCellTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(HybridNPHPREndTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(&CFWangLandauEvery,sizeof(int),1,FilePtr);
  fwrite(&TargetAccRatioLambdaChange,sizeof(REAL),1,FilePtr);

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void AllocateMCMovesMemory(void)
{
  int i,j;

  TrialPosition=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TrialAnisotropicPosition=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  ElectricFieldAtTrialPosition=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  ReferenceElectricFieldAtTrialPosition=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  InducedDipoleAtTrialPosition=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  for(i=0;i<NumberOfSystems;i++)
  {
    TrialPosition[i]=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));
    TrialAnisotropicPosition[i]=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));
    ElectricFieldAtTrialPosition[i]=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));
    ReferenceElectricFieldAtTrialPosition[i]=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));
    InducedDipoleAtTrialPosition[i]=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));
  }

  cord=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));

  // delta energies
  UHostVDWDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostChargeChargeRealDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostChargeBondDipoleRealDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBondDipoleBondDipoleRealDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UCationVDWDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateVDWDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UCationChargeChargeRealDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateChargeChargeRealDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UCationChargeBondDipoleRealDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateChargeBondDipoleRealDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UCationBondDipoleBondDipoleRealDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBondDipoleBondDipoleRealDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UCationFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostPolarizationNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbatePolarizationNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationPolarizationNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UPolarizationNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostBackPolarizationNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBackPolarizationNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBackPolarizationNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBackPolarizationNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  TranslationAttempts=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TranslationAccepted=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalTranslationAttempts=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalTranslationAccepted=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  MaximumTranslation=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));

  RandomTranslationAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  RandomTranslationAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  RotationAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  RotationAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  MaximumRotation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  SwapAddAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  SwapAddAccepted=(REAL(**)[2])calloc(NumberOfSystems,sizeof(REAL(*)[2]));

  SwapRemoveAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  SwapRemoveAccepted=(REAL(**)[2])calloc(NumberOfSystems,sizeof(REAL(*)[2]));

  CFSwapLambdaAttempts=(REAL(**)[3])calloc(NumberOfSystems,sizeof(REAL(*)[3]));
  CFSwapLambdaAccepted=(REAL(**)[3])calloc(NumberOfSystems,sizeof(REAL(*)[3]));
  TotalCFSwapLambdaAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TotalCFSwapLambdaAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  MaximumCFLambdaChange=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  CFCBSwapLambdaAttempts=(REAL(**)[3])calloc(NumberOfSystems,sizeof(REAL(*)[3]));
  CFCBSwapLambdaAccepted=(REAL(**)[3])calloc(NumberOfSystems,sizeof(REAL(*)[3]));
  TotalCFCBSwapLambdaAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TotalCFCBSwapLambdaAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  MaximumCFCBLambdaChange=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  ReinsertionAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ReinsertionAccepted=(REAL(**)[2])calloc(NumberOfSystems,sizeof(REAL(*)[2]));

  PartialReinsertionAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  PartialReinsertionAccepted=(REAL(**)[2])calloc(NumberOfSystems,sizeof(REAL(*)[2]));

  ReinsertionInPlaceAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ReinsertionInPlaceAccepted=(REAL(**)[2])calloc(NumberOfSystems,sizeof(REAL(*)[2]));

  TranslationInPlaceAttempts=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TranslationInPlaceAccepted=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalTranslationInPlaceAttempts=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalTranslationInPlaceAccepted=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));

  ReinsertionInPlaneAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ReinsertionInPlaneAccepted=(REAL(**)[2])calloc(NumberOfSystems,sizeof(REAL(*)[2]));

  TranslationInPlaneAttempts=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TranslationInPlaneAccepted=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalTranslationInPlaneAttempts=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalTranslationInPlaneAccepted=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  MaximumTranslationInPlane=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));

  IdentityChangeAttempts=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  IdentityChangeAccepted=(REAL(***)[2])calloc(NumberOfSystems,sizeof(REAL(**)[2]));

  GibbsIdentityChangeAttempts=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  GibbsIdentityChangeAccepted=(REAL(***)[2])calloc(NumberOfSystems,sizeof(REAL(**)[2]));

  CFLambdaHistogram=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

  ParallelTemperingAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ParallelTemperingAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  HyperParallelTemperingAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  HyperParallelTemperingAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  ParallelMolFractionAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  ParallelMolFractionAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  for(i=0;i<NumberOfSystems;i++)
  {
    TranslationAttempts[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    TranslationAccepted[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    TotalTranslationAttempts[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    TotalTranslationAccepted[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    MaximumTranslation[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));

    RandomTranslationAttempts[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    RandomTranslationAccepted[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));

    RotationAttempts[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    RotationAccepted[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    MaximumRotation[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));

    SwapAddAttempts[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    SwapAddAccepted[i]=(REAL(*)[2])calloc(NumberOfComponents,sizeof(REAL[2]));

    SwapRemoveAttempts[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    SwapRemoveAccepted[i]=(REAL(*)[2])calloc(NumberOfComponents,sizeof(REAL[2]));

    CFSwapLambdaAttempts[i]=(REAL(*)[3])calloc(NumberOfComponents,sizeof(REAL[3]));
    CFSwapLambdaAccepted[i]=(REAL(*)[3])calloc(NumberOfComponents,sizeof(REAL[3]));
    TotalCFSwapLambdaAttempts[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    TotalCFSwapLambdaAccepted[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    MaximumCFLambdaChange[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));

    CFCBSwapLambdaAttempts[i]=(REAL(*)[3])calloc(NumberOfComponents,sizeof(REAL[3]));
    CFCBSwapLambdaAccepted[i]=(REAL(*)[3])calloc(NumberOfComponents,sizeof(REAL[3]));
    TotalCFCBSwapLambdaAttempts[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    TotalCFCBSwapLambdaAccepted[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    MaximumCFCBLambdaChange[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));

    ReinsertionAttempts[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    ReinsertionAccepted[i]=(REAL(*)[2])calloc(NumberOfComponents,sizeof(REAL[2]));

    PartialReinsertionAttempts[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    PartialReinsertionAccepted[i]=(REAL(*)[2])calloc(NumberOfComponents,sizeof(REAL[2]));

    ReinsertionInPlaceAttempts[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    ReinsertionInPlaceAccepted[i]=(REAL(*)[2])calloc(NumberOfComponents,sizeof(REAL[2]));

    TranslationInPlaceAttempts[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    TranslationInPlaceAccepted[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    TotalTranslationInPlaceAttempts[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    TotalTranslationInPlaceAccepted[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));

    ReinsertionInPlaneAttempts[i]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
    ReinsertionInPlaneAccepted[i]=(REAL(*)[2])calloc(NumberOfComponents,sizeof(REAL[2]));

    TranslationInPlaneAttempts[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    TranslationInPlaneAccepted[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    TotalTranslationInPlaneAttempts[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    TotalTranslationInPlaneAccepted[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));
    MaximumTranslationInPlane[i]=(VECTOR*)calloc(NumberOfComponents,sizeof(VECTOR));

    IdentityChangeAttempts[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    IdentityChangeAccepted[i]=(REAL(**)[2])calloc(NumberOfComponents,sizeof(REAL(*)[2]));

    GibbsIdentityChangeAttempts[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    GibbsIdentityChangeAccepted[i]=(REAL(**)[2])calloc(NumberOfComponents,sizeof(REAL(*)[2]));

    CFLambdaHistogram[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

    for(j=0;j<NumberOfComponents;j++)
    {
      IdentityChangeAttempts[i][j]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
      IdentityChangeAccepted[i][j]=(REAL(*)[2])calloc(NumberOfComponents,sizeof(REAL[2]));

      GibbsIdentityChangeAttempts[i][j]=(REAL*)calloc(NumberOfComponents,sizeof(REAL));
      GibbsIdentityChangeAccepted[i][j]=(REAL(*)[2])calloc(NumberOfComponents,sizeof(REAL[2]));

      CFLambdaHistogram[i][j]=(REAL*)calloc(Components[j].CFLambdaHistogramSize,sizeof(REAL));
    }

    ParallelTemperingAttempts[i]=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    ParallelTemperingAccepted[i]=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

    HyperParallelTemperingAttempts[i]=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    HyperParallelTemperingAccepted[i]=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

    ParallelMolFractionAttempts[i]=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    ParallelMolFractionAccepted[i]=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  }

  ChiralInversionAttempts=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  ChiralInversionAccepted=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  VolumeChangeAttempts=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  VolumeChangeAccepted=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  TotalVolumeChangeAttempts=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  TotalVolumeChangeAccepted=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  MaximumVolumeChange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  BoxShapeChangeAttempts=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  BoxShapeChangeAccepted=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  MaximumBoxShapeChange=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));

  GibbsVolumeChangeAttempts=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  GibbsVolumeChangeAccepted=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  MaximumGibbsVolumeChange=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  GibbsSwapAttempts=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
  GibbsSwapAccepted=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
  for(i=0;i<NumberOfComponents;i++)
  {
    GibbsSwapAttempts[i]=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    GibbsSwapAccepted[i]=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  }

  FrameworkChangeAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  FrameworkChangeAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TotalFrameworkChangeAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TotalFrameworkChangeAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  FrameworkMaximumTranslation=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  for(i=0;i<NumberOfSystems;i++)
  {
    FrameworkChangeAttempts[i]=(REAL*)calloc(NumberOfPseudoAtoms,sizeof(REAL));
    FrameworkChangeAccepted[i]=(REAL*)calloc(NumberOfPseudoAtoms,sizeof(REAL));
    TotalFrameworkChangeAttempts[i]=(REAL*)calloc(NumberOfPseudoAtoms,sizeof(REAL));
    TotalFrameworkChangeAccepted[i]=(REAL*)calloc(NumberOfPseudoAtoms,sizeof(REAL));
    FrameworkMaximumTranslation[i]=(REAL*)calloc(NumberOfPseudoAtoms,sizeof(REAL));
  }

  FrameworkShiftAttempts=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  FrameworkShiftAccepted=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalFrameworkShiftAttempts=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalFrameworkShiftAccepted=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  FrameworkMaximumShiftTranslation=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));

  for(i=0;i<NumberOfSystems;i++)
  {
    FrameworkShiftAttempts[i]=(VECTOR*)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR));
    FrameworkShiftAccepted[i]=(VECTOR*)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR));
    TotalFrameworkShiftAttempts[i]=(VECTOR*)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR));
    TotalFrameworkShiftAccepted[i]=(VECTOR*)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR));
    FrameworkMaximumShiftTranslation[i]=(VECTOR*)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR));
  }

  HybridNVEDrift=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEDriftCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNVEAttempts=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAccepted=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNVEStartTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartTranslationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartRotationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartTemperatureFramework=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartTemperatureAdsorbate=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartTemperatureCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartTemperatureCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartTemperatureTranslationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartTemperatureRotationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartTemperatureFrameworkCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartTemperatureAdsorbateCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEStartTemperatureCationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNVEEndTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndTranslationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndRotationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndTemperatureFramework=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndTemperatureAdsorbate=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndTemperatureCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndTemperatureCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndTemperatureTranslationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndTemperatureRotationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndTemperatureFrameworkCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndTemperatureAdsorbateCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEEndTemperatureCationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));


  HybridNPHDrift=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHDriftCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNPHAttempts=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHAccepted=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNPHStartTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartCellTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartTranslationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartRotationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartTemperatureFramework=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartTemperatureAdsorbate=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartTemperatureCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartTemperatureCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartTemperatureTranslationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartTemperatureRotationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartTemperatureFrameworkCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartTemperatureAdsorbateCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHStartTemperatureCationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNPHEndTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndCellTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndTranslationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndRotationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndTemperatureFramework=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndTemperatureAdsorbate=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndTemperatureCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndTemperatureCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndTemperatureTranslationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndTemperatureRotationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndTemperatureFrameworkCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndTemperatureAdsorbateCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHEndTemperatureCationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));


  HybridNPHPRDrift=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRDriftCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNPHPRAttempts=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRAccepted=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNPHPRStartTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartCellTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartTranslationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartRotationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartTemperatureFramework=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartTemperatureAdsorbate=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartTemperatureCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartTemperatureCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartTemperatureTranslationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartTemperatureRotationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartTemperatureFrameworkCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartTemperatureAdsorbateCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPRStartTemperatureCationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNPHPREndTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndCellTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndTranslationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndRotationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndTemperatureFramework=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndTemperatureAdsorbate=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndTemperatureCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndTemperatureCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndTemperatureTranslationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndTemperatureRotationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndTemperatureFrameworkCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndTemperatureAdsorbateCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNPHPREndTemperatureCationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UDeltaPolarization=0.0;
}

void ReadRestartMcMoves(FILE *FilePtr)
{
  int i,j;
  REAL Check;

  AllocateMCMovesMemory();

  for(i=0;i<NumberOfSystems;i++)
  {
    fread(TranslationAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(TranslationAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(TotalTranslationAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(TotalTranslationAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(MaximumTranslation[i],sizeof(VECTOR),NumberOfComponents,FilePtr);

    fread(RandomTranslationAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(RandomTranslationAccepted[i],sizeof(REAL),NumberOfComponents,FilePtr);

    fread(RotationAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(RotationAccepted[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(MaximumRotation[i],sizeof(REAL),NumberOfComponents,FilePtr);

    fread(SwapAddAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(SwapAddAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fread(SwapRemoveAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(SwapRemoveAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fread(CFSwapLambdaAttempts[i],sizeof(REAL[3]),NumberOfComponents,FilePtr);
    fread(CFSwapLambdaAccepted[i],sizeof(REAL[3]),NumberOfComponents,FilePtr);
    fread(TotalCFSwapLambdaAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(TotalCFSwapLambdaAccepted[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(MaximumCFLambdaChange[i],sizeof(REAL),NumberOfComponents,FilePtr);

    fread(CFCBSwapLambdaAttempts[i],sizeof(REAL[3]),NumberOfComponents,FilePtr);
    fread(CFCBSwapLambdaAccepted[i],sizeof(REAL[3]),NumberOfComponents,FilePtr);
    fread(TotalCFCBSwapLambdaAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(TotalCFCBSwapLambdaAccepted[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(MaximumCFCBLambdaChange[i],sizeof(REAL),NumberOfComponents,FilePtr);

    fread(ReinsertionAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(ReinsertionAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fread(PartialReinsertionAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(PartialReinsertionAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fread(ReinsertionInPlaceAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(ReinsertionInPlaceAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fread(TranslationInPlaceAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(TranslationInPlaceAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(TotalTranslationInPlaceAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(TotalTranslationInPlaceAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);

    fread(ReinsertionInPlaneAttempts[i],sizeof(REAL),NumberOfComponents,FilePtr);
    fread(ReinsertionInPlaneAccepted[i],sizeof(REAL[2]),NumberOfComponents,FilePtr);

    fread(TranslationInPlaneAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(TranslationInPlaneAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(TotalTranslationInPlaneAttempts[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(TotalTranslationInPlaneAccepted[i],sizeof(VECTOR),NumberOfComponents,FilePtr);
    fread(MaximumTranslationInPlane[i],sizeof(VECTOR),NumberOfComponents,FilePtr);

    for(j=0;j<NumberOfComponents;j++)
    {
      fread(IdentityChangeAttempts[i][j],sizeof(REAL),NumberOfComponents,FilePtr);
      fread(IdentityChangeAccepted[i][j],sizeof(REAL[2]),NumberOfComponents,FilePtr);

      fread(GibbsIdentityChangeAttempts[i][j],sizeof(REAL),NumberOfComponents,FilePtr);
      fread(GibbsIdentityChangeAccepted[i][j],sizeof(REAL[2]),NumberOfComponents,FilePtr);

      fread(CFLambdaHistogram[i][j],sizeof(REAL),Components[j].CFLambdaHistogramSize,FilePtr);
    }

    fread(ParallelTemperingAttempts[i],sizeof(REAL),NumberOfSystems,FilePtr);
    fread(ParallelTemperingAccepted[i],sizeof(REAL),NumberOfSystems,FilePtr);

    fread(HyperParallelTemperingAttempts[i],sizeof(REAL),NumberOfSystems,FilePtr);
    fread(HyperParallelTemperingAccepted[i],sizeof(REAL),NumberOfSystems,FilePtr);

    fread(ParallelMolFractionAttempts[i],sizeof(REAL),NumberOfSystems,FilePtr);
    fread(ParallelMolFractionAccepted[i],sizeof(REAL),NumberOfSystems,FilePtr);
  }
  fread(&ParallelMolFractionComponentA,sizeof(int),1,FilePtr);
  fread(&ParallelMolFractionComponentB,sizeof(int),1,FilePtr);

  fread(ChiralInversionAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(ChiralInversionAccepted,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(VolumeChangeAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(VolumeChangeAccepted,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(TotalVolumeChangeAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(TotalVolumeChangeAccepted,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(MaximumVolumeChange,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(BoxShapeChangeAttempts,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(BoxShapeChangeAccepted,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(MaximumBoxShapeChange,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);

  fread(GibbsVolumeChangeAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(GibbsVolumeChangeAccepted,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(MaximumGibbsVolumeChange,sizeof(REAL),NumberOfSystems,FilePtr);
  
  for(i=0;i<NumberOfComponents;i++)
  {
    fread(GibbsSwapAttempts[i],sizeof(REAL),NumberOfSystems,FilePtr);
    fread(GibbsSwapAccepted[i],sizeof(REAL),NumberOfSystems,FilePtr);
  }

  for(i=0;i<NumberOfSystems;i++)
  {
    fread(FrameworkShiftAttempts[i],sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
    fread(FrameworkShiftAccepted[i],sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
    fread(TotalFrameworkShiftAttempts[i],sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
    fread(TotalFrameworkShiftAccepted[i],sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
    fread(FrameworkMaximumShiftTranslation[i],sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
  }

  fread(&NumberOfHybridNVESteps,sizeof(int),1,FilePtr);
  fread(HybridNVEDrift,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEDriftCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEAccepted,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(HybridNVEStartTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEStartTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(HybridNVEEndTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNVEEndTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);


  fread(&NumberOfHybridNPHSteps,sizeof(int),1,FilePtr);
  fread(HybridNPHDrift,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHDriftCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHAccepted,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(HybridNPHStartTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartCellTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHStartTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(HybridNPHEndTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndCellTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHEndTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);


  fread(&NumberOfHybridNPHPRSteps,sizeof(int),1,FilePtr);
  fread(HybridNPHPRDrift,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRDriftCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRAttempts,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRAccepted,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(HybridNPHPRStartTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartCellTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPRStartTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(HybridNPHPREndTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndCellTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndTranslationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndRotationalTemperature,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndTemperatureFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndTemperatureAdsorbate,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndTemperatureCation,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndTemperatureCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndTemperatureTranslationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndTemperatureRotationCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndTemperatureFrameworkCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndTemperatureAdsorbateCount,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(HybridNPHPREndTemperatureCationCount,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(&CFWangLandauEvery,sizeof(int),1,FilePtr);
  fread(&TargetAccRatioLambdaChange,sizeof(REAL),1,FilePtr);

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    printf("Error in binary restart-file (ReadRestartMcMoves)\n");
    exit(0);
  }
}

