/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'ewald.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef EWALD_H
#define EWALD_H

#include "simulation.h"
#include "vector.h"
#include "utils.h"

extern int *NumberOfKVectors;
extern int OmitEwaldFourier;

extern REAL *Alpha;
extern INT_VECTOR3 *kvec;
extern REAL *ReciprocalCutOffSquared;
extern REAL EwaldPrecision;
extern REAL DielectricConstantOfTheMedium;

extern REAL *UAdsorbateAdsorbateChargeChargeFourierDelta;
extern REAL *UCationCationChargeChargeFourierDelta;
extern REAL *UAdsorbateCationChargeChargeFourierDelta;
extern REAL *UHostCationChargeChargeFourierDelta;
extern REAL *UHostAdsorbateChargeChargeFourierDelta;
extern REAL *UHostHostChargeChargeFourierDelta;

extern REAL *UAdsorbateAdsorbateChargeBondDipoleFourierDelta;
extern REAL *UCationCationChargeBondDipoleFourierDelta;
extern REAL *UAdsorbateCationChargeBondDipoleFourierDelta;
extern REAL *UHostCationChargeBondDipoleFourierDelta;
extern REAL *UHostAdsorbateChargeBondDipoleFourierDelta;
extern REAL *UHostHostChargeBondDipoleFourierDelta;

extern REAL *UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta;
extern REAL *UCationCationBondDipoleBondDipoleFourierDelta;
extern REAL *UAdsorbateCationBondDipoleBondDipoleFourierDelta;
extern REAL *UHostCationBondDipoleBondDipoleFourierDelta;
extern REAL *UHostAdsorbateBondDipoleBondDipoleFourierDelta;
extern REAL *UHostHostBondDipoleBondDipoleFourierDelta;


// net charges
extern REAL *NetChargeSystem;
extern REAL *NetChargeFramework;
extern REAL *NetChargeCations;
extern REAL *NetChargeAdsorbates;

extern REAL NetChargeCationDelta;
extern REAL NetChargeAdsorbateDelta;

int AllocateNewEwald(void);
void SetupKVectors(void);
int PrecomputeFixedEwaldContributions(void);
int PrecomputeTotalEwaldContributions(void);
void InitializeEwald(REAL precision,int Automatic);
int EwaldEnergyIon(void);
int EwaldFourierEnergy(void);
int EwaldFourierForce(void);
int EwaldFourierBornTerm(void);
int CalculateEwaldFourierAdsorbate(int New,int old,int mol,int store);
int CalculateEwaldFourierAdsorbateCF(int New,int old,int mol,int store);
int CalculateEwaldFourierCation(int New,int old,int mol,int store);
int CalculateEwaldFourierCationCF(int New,int old,int mol,int store);
void AcceptEwaldAdsorbateMove(int A);
void AcceptEwaldCationMove(int A);
void SaveCurrentKVectors(int A,int B);
void RetrieveStoredKVectors(int A,int B);
void SaveCurrentEwaldStructureFactors(int A,int B);
void RetrieveStoredEwaldStructureFactors(int A,int B);
void SwapEwaldSystem(int A,int B);

int AllocateEwaldMemory(void);
int ReallocateEwaldChargeMemory(void);
int ReallocateEwaldBondDipoleMemory(void);

void AcceptEwaldFrameworkDiplacementMove(void);
void AcceptEwaldFrameworkDiplacementMoveRigid(void);

int EwaldFourierStaticElectricField(void);

int CalculateEwaldFourierFrameworkDisplacement(void);
int CalculateEwaldFourierFrameworkAtomTranslate(int index);
void AcceptEwaldFrameworkMove(int A);

// RXMC
int CalculateEwaldFourierAdsorbateRXMC(int reaction,REAL Lambda1,REAL Lambda2,REAL LambdaNew,int **LambdaRetraceMolecules,int direction,int store);
int CalculateEwaldFourierAdsorbateRXMC2(int reaction,REAL LambdaNew,int direction,int store);

int ComputeStaticElectricFieldEwaldAdsorbateDifference(int NewMolecule,int OldMolecule,int mol,int store);

void CalculateEwaldFourierBornTerm(REAL *Energy,REAL* Gradient,REAL_MATRIX3x3 *StrainDerivativeTensor);

void CalculateEwaldFourierCrossTerms(REAL *Energy,REAL* Gradient,REAL_MATRIX3x3 *StrainDerivativeTensor,REAL_MATRIX CrossTerms);

REAL WolfEnergyFull(void);

int ComputeStaticElectricFieldEwald(int New,int excl_ads,int excl_cation);

void WriteRestartEwald(FILE *FilePtr);
void ReadRestartEwald(FILE *FilePtr);

void ComputeElectricFieldFromInducedDipolesEwald(void);
void ComputeElectricFieldFromInducedDipolesEwaldMC(int New,int excl_ads,int excl_cation);

void ComputeInducedDipolesForcesEwald(void);

int CalculateEwaldFourierDerivatives(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainFirstDerivative,int ComputeGradient,int ComputeHessian);

REAL EwaldFourierElectrostaticPotentialAdsorbate(int m,int l);
REAL EwaldFourierElectrostaticPotentialCation(int m,int l);

int CalculateEwaldFourierPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainFirstDerivative,int ComputeGradient,int ComputeHessian);

#endif
