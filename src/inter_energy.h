/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'inter_energy.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef INTER_ENERGY_H
#define INTER_ENERGY_H

#include "utils.h"
#include "simulation.h"

REAL CalculateInterVDWSelfEnergyCorrectionNew(void);
REAL CalculateInterVDWSelfEnergyCorrectionAdsorbateOld(int mol);
REAL CalculateInterVDWSelfEnergyCorrectionCationOld(int mol);
REAL CalculateInterChargeChargeSelfEnergyCorrectionNew(void);
REAL CalculateInterChargeChargeSelfEnergyCorrectionAdsorbateOld(int mol);
REAL CalculateInterChargeChargeSelfEnergyCorrectionCationOld(int mol);

REAL CalculateInterVDWEnergyAdsorbateAtPosition(POINT posA,int typeA,int exclude,REAL scaling);
int CalculateInterChargeEnergyAdsorbateAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,int exclude,REAL scaling);
int CalculateInterBondDipoleEnergyAdsorbateAtPosition(POINT posA1,POINT posA2,REAL DipoleMagnitudeA,REAL *UChargeBondDipole,REAL *UBondDipoleBondDipole,int exclude);

REAL CalculateInterVDWEnergyCationAtPosition(POINT posA,int typeA,int exclude,REAL scaling);
int CalculateInterChargeEnergyCationAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,int exclude,REAL scaling);
int CalculateInterBondDipoleEnergyCationAtPosition(POINT posA1,POINT posA2,REAL DipoleMagnitudeA,REAL *UChargeBondDipole,REAL *UBondDipoleBondDipole,int exclude);

int CalculateTotalInterVDWEnergy(void);
int CalculateTotalInterChargeChargeCoulombEnergy(void);
int CalculateTotalInterReplicaChargeChargeCoulombEnergy(void);
int CalculateTotalInterChargeBondDipoleCoulombEnergy(void);
int CalculateTotalInterBondDipoleBondDipoleCoulombEnergy(void);

int CalculateInterVDWEnergyDifferenceAdsorbate(int m,int comp,int New,int Old);
int CalculateInterVDWEnergyDifferenceCation(int m,int comp,int New,int Old);

int CalculateInterChargeChargeEnergyDifferenceAdsorbate(int m,int comp,int New,int Old);
int CalculateInterChargeChargeEnergyDifferenceCation(int m,int comp,int New,int Old);

int CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(int m,int comp,int New,int Old);
int CalculateInterChargeBondDipoleEnergyDifferenceCation(int m,int comp,int New,int Old);

int CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(int m,int comp,int New,int Old);
int CalculateInterBondDipoleBondDipoleEnergyDifferenceCation(int m,int comp,int New,int Old);

REAL CalculateInterChargeElectrostaticPotential(POINT posA);

REAL CalculateInterVDWEnergyCorrectionAdsorbate(VECTOR* Positions,VECTOR *AnisotropicPositions,int exclude);
REAL CalculateInterVDWEnergyCorrectionCation(VECTOR* Positions,VECTOR *AnisotropicPositions,int exclude);

#endif
