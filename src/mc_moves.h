/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'mc_moves.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef MC_MOVES_H
#define MC_MOVES_H

#include "simulation.h"
#include "molecule.h"
#include "framework.h"

// delta energies
extern REAL *UHostVDWDelta;
extern REAL *UHostChargeChargeRealDelta;
extern REAL *UHostChargeBondDipoleRealDelta;
extern REAL *UHostBondDipoleBondDipoleRealDelta;
extern REAL *UHostFourierDelta;

extern REAL *UCationVDWDelta;
extern REAL *UAdsorbateVDWDelta;

extern REAL *UCationChargeChargeRealDelta;
extern REAL *UAdsorbateChargeChargeRealDelta;

extern REAL *UCationChargeBondDipoleRealDelta;
extern REAL *UAdsorbateChargeBondDipoleRealDelta;

extern REAL *UCationBondDipoleBondDipoleRealDelta;
extern REAL *UAdsorbateBondDipoleBondDipoleRealDelta;

extern REAL *UCationFourierDelta;
extern REAL *UAdsorbateFourierDelta;


extern VECTOR **MaximumTranslation;
extern VECTOR **MaximumTranslationInPlane;
extern REAL **MaximumRotation;

extern VECTOR **TrialPosition;
extern VECTOR **TrialAnisotropicPosition;
extern VECTOR **ElectricFieldAtTrialPosition;
extern VECTOR **InducedDipoleAtTrialPosition;
extern REAL TargetAccRatioTranslation;

extern REAL *MaximumVolumeChange;
extern REAL TargetAccRatioVolumeChange;

extern REAL_MATRIX3x3 *MaximumBoxShapeChange;
extern REAL TargetAccRatioBoxShapeChange;

extern REAL *MaximumGibbsVolumeChange;
extern REAL TargetAccRatioGibbsVolumeChange;

extern REAL **FrameworkMaximumTranslation;
extern VECTOR **FrameworkMaximumShiftTranslation;

extern int ParallelMolFractionComponentA;
extern int ParallelMolFractionComponentB;

extern int NumberOfHybridNVESteps;
extern int NumberOfHybridNPHSteps;
extern int NumberOfHybridNPHPRSteps;

//----------------------------------------------------------------------------------------
// CFC-RXMC Parameters
//----------------------------------------------------------------------------------------

extern REAL *MaximumCFCRXMCLambdaChange;
extern REAL TargetAccRatioCFCRXMCLambdaChange;

extern int CFLambdaHistogramSize;
extern int CFWangLandauEvery;
extern REAL **MaximumCFLambdaChange;

//----------------------------------------------------------------------------------------

void InitializeMCMovesStatisticsAllSystems(void);

void TranslationMove(void);
void RandomTranslationMove(void);
void RotationMove(void);
void PartialReinsertionMove(void);
void ReinsertionMove(void);
void ReinsertionInPlaceMove(void);
void ReinsertionInPlaneMove(void);
void SwapAddMove(void);
void SwapRemoveMove(void);
void IdentityChangeMove(void);
void WidomMove(void);
void GibbsParticleTransferMove(void);
void GibbsIdentityChangeMove(void);

int TranslationMoveAdsorbate(void);
int RandomTranslationMoveAdsorbate(void);
int RotationMoveAdsorbate(void);
int SwapAddAdsorbateMove(void);
int SwapRemoveAdsorbateMove(void);
int ReinsertionAdsorbateMove(void);
int ReinsertionInPlaceAdsorbateMove(void);
int ReinsertionInPlaneAdsorbateMove(void);
int IdentityChangeAdsorbateMove(void);
int PartialReinsertionAdsorbateMove(void);

void OptimizeTranslationAcceptence(void);
void OptimizeTranslationInPlaneAcceptence(void);
void OptimizeVolumeChangeAcceptence(void);
void OptimizeFrameworkChangeAcceptence(void);
void OptimizeFrameworkShiftAcceptence(void);

int TranslationMoveCation(void);
int RandomTranslationMoveCation(void);
int RotationMoveCation(void);
int SwapAddCationMove(void);
int SwapRemoveCationMove(void);
int ReinsertionCationMove(void);
int ReinsertionInPlaceCationMove(void);
int ReinsertionInPlaneCationMove(void);
int PartialReinsertionCationMove(void);
int IdentityChangeCationMove(void);

int ParallelTemperingMove(void);
int HyperParallelTemperingMove(void);
int ParallelMolFractionMove(void);
int ChiralInversionMove(void);
int VolumeMove(void);
int BoxShapeChangeMove(void);
int GibbsParticleTransferAdsorbateMove(void);
int GibbsParticleTransferCationMove(void);
int GibbsVolumeMove(void);
int GibbsIdentityChangeAdsorbateMove(void);
int FrameworkChangeMove(void);
int FrameworkShiftMove(void);

void PrintGibbsSwapStatistics(FILE *FilePtr);
void PrintTranslationStatistics(FILE *FilePtr);
void PrintRandomTranslationStatistics(FILE *FilePtr);
void PrintRotationStatistics(FILE *FilePtr);
void PrintSwapAddStatistics(FILE *FilePtr);
void PrintSwapRemoveStatistics(FILE *FilePtr);
void PrintReinsertionStatistics(FILE *FilePtr);
void PrintReinsertionInPlaneStatistics(FILE *FilePtr);
void PrintReinsertionInPlaceStatistics(FILE *FilePtr);
void PrintPartialReinsertionStatistics(FILE *FilePtr);
void PrintIdentityChangeStatistics(FILE *FilePtr);
void PrintParallelTemperingStatistics(FILE *FilePtr);
void PrintHyperParallelTemperingStatistics(FILE *FilePtr);
void PrintParallelMolFractionStatistics(FILE *FilePtr);
void PrintChiralInversionStatistics(FILE *FilePtr);
void PrintVolumeChangeStatistics(FILE *FilePtr);
void PrintBoxShapeChangeStatistics(FILE *FilePtr);
void PrintGibbsVolumeChangeStatistics(FILE *FilePtr);
void PrintFrameworkStatistics(FILE *FilePtr);
void PrintFrameworkShiftStatistics(FILE *FilePtr);
void PrintGibbsIdentityChangeStatistics(FILE *FilePtr);
void PrintCFSwapLambdaStatistics(FILE *FilePtr);
void PrintCFCBSwapLambdaStatistics(FILE *FilePtr);

REAL WidomAdsorbateMove(void);
REAL WidomCationMove(void);

void HybridNVEMove(void);
void PrintHybridNVEStatistics(FILE *FilePtr);

void HybridNPHMove(void);
void PrintHybridNPHStatistics(FILE *FilePtr);

void HybridNPHPRMove(void);
void PrintHybridNPHPRStatistics(FILE *FilePtr);

void SurfaceAreaMove(void);

void CFWangLandauIteration(int Switch);
int CFSwapLambaMove(void);
void OptimizeCFLambdaChangeAcceptence(void);
int CFCBSwapLambaMove(void);

void WriteRestartMcMoves(FILE *FilePtr);
void AllocateMCMovesMemory(void);
void ReadRestartMcMoves(FILE *FilePtr);

#endif
