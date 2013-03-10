/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'cbmc.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef CBMC_H
#define CBMC_H

#include "molecule.h"

enum {LJ_BIASING,LJ_AND_REAL_BIASING};

enum {CBMC_INSERTION,CBMC_PARTIAL_INSERTION,CBMC_DELETION,CBMC_RETRACE_REINSERTION};

extern int BiasingMethod;

#define MAX_NUMBER_OF_TRIAL_POSITIONS 2000

extern int NumberOfTrialPositions;
extern int NumberOfTrialPositionsForTheFirstBead;
extern int NumberOfTrialMovesPerOpenBead;
extern int NumberOfTrialPositionsTorsion;
extern REAL TargetAccRatioSmallMCScheme;
extern REAL EnergyOverlapCriteria;
extern REAL MinimumRosenbluthFactor;

extern int MaxNumberOfBeads;
extern int MaxNumberOfBonds;
extern int MaxNumberOfBondDipoles;
extern int MaxNumberOfBends;
extern int MaxNumberOfBendBends;
extern int MaxNumberOfInversionBends;
extern int MaxNumberOfUreyBradleys;
extern int MaxNumberOfTorsions;
extern int MaxNumberOfImproperTorsions;
extern int MaxNumberOfOutOfPlanes;
extern int MaxNumberOfBondBonds;
extern int MaxNumberOfBondBends;
extern int MaxNumberOfBondTorsions;
extern int MaxNumberOfBendTorsions;
extern int MaxNumberOfIntraVDW;
extern int MaxNumberOfIntraChargeCharge;
extern int MaxNumberOfIntraChargeBondDipole;
extern int MaxNumberOfIntraBondDipoleBondDipole;

extern REAL *UBondNew;
extern REAL *UUreyBradleyNew;
extern REAL *UBendNew;
extern REAL *UBendBendNew;
extern REAL *UInversionBendNew;
extern REAL *UTorsionNew;
extern REAL *UImproperTorsionNew;
extern REAL *UBondBondNew;
extern REAL *UBondBendNew;
extern REAL *UBondTorsionNew;
extern REAL *UBendTorsionNew;
extern REAL *UIntraVDWNew;
extern REAL *UIntraChargeChargeNew;
extern REAL *UIntraChargeBondDipoleNew;
extern REAL *UIntraBondDipoleBondDipoleNew;

extern REAL *UHostVDWNew;
extern REAL *UAdsorbateVDWNew;
extern REAL *UCationVDWNew;
extern REAL *UHostChargeChargeNew;
extern REAL *UAdsorbateChargeChargeNew;
extern REAL *UCationChargeChargeNew;
extern REAL *UHostChargeBondDipoleNew;
extern REAL *UAdsorbateChargeBondDipoleNew;
extern REAL *UCationChargeBondDipoleNew;
extern REAL *UHostBondDipoleBondDipoleNew;
extern REAL *UAdsorbateBondDipoleBondDipoleNew;
extern REAL *UCationBondDipoleBondDipoleNew;

extern REAL *UBondOld;
extern REAL *UUreyBradleyOld;
extern REAL *UBendOld;
extern REAL *UBendBendOld;
extern REAL *UInversionBendOld;
extern REAL *UTorsionOld;
extern REAL *UImproperTorsionOld;
extern REAL *UBondBondOld;
extern REAL *UBondBendOld;
extern REAL *UBondTorsionOld;
extern REAL *UBendTorsionOld;
extern REAL *UIntraVDWOld;
extern REAL *UIntraChargeChargeOld;
extern REAL *UIntraChargeBondDipoleOld;
extern REAL *UIntraBondDipoleBondDipoleOld;

extern REAL *UHostVDWOld;
extern REAL *UAdsorbateVDWOld;
extern REAL *UCationVDWOld;
extern REAL *UHostChargeChargeOld;
extern REAL *UAdsorbateChargeChargeOld;
extern REAL *UCationChargeChargeOld;
extern REAL *UHostChargeBondDipoleOld;
extern REAL *UAdsorbateChargeBondDipoleOld;
extern REAL *UCationChargeBondDipoleOld;
extern REAL *UHostBondDipoleBondDipoleOld;
extern REAL *UAdsorbateBondDipoleBondDipoleOld;
extern REAL *UCationBondDipoleBondDipoleOld;

extern int NumberOfBeadsAlreadyPlaced;                         // number of atoms that are already placed
extern int *BeadsAlreadyPlaced;
extern VECTOR **NewPosition;
extern VECTOR *OldPosition;
extern VECTOR **NewVelocity;
extern VECTOR **NewForce;
extern VECTOR FirstBeadPosition;
extern VECTOR **TrialPositions;
extern REAL *CFVDWScaling;
extern REAL *CFChargeScaling;
extern int OVERLAP;           

void CalculateAnisotropicTrialPositions(int TypeMolA,VECTOR *TrialPosition,VECTOR *TrialAnisotropicPosition);

REAL GrowMolecule(int Iicode);
REAL RetraceMolecule(int Iicode);

void GrowReservoirMolecule(void);

void MakeInitialAdsorbate(void);
void MakeInitialAdsorbates(int n,int type);
void MakeInitialCation(void);
void MakeInitialCations(int n,int type);

void RescaleMaximumRotationAnglesSmallMC(void);
void InitializeSmallMCStatisticsAllSystems(void);
void PrintSmallMCAddStatistics(FILE *FilePtr);

void CheckConfigMoves(void);

void WriteRestartCBMC(FILE *FilePtr);
void AllocateCBMCMemory(void);
void ReadRestartCBMC(FILE *FilePtr);

#endif

