/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_energy.h' is part of RASPA.

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

#ifndef FRAMEWORK_ENERGY_H
#define FRAMEWORK_ENERGY_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <molecule.h>

REAL CalculateFrameworkVDWEnergyAtPosition(POINT posA,int typeA);
REAL CalculateFrameworkChargeChargeEnergyAtPosition(POINT pos,int typeA);

void CalculateFrameworkChargeEnergyAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole);
void CalculateFrameworkBondDipoleEnergyAtPosition(VECTOR posA1,VECTOR posA2,REAL DipoleMagnitudeA,REAL *UChargeBondDipole,REAL *UBondDipoleBondDipole);

REAL CalculateFrameworkBondEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkUreyBradleyEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBendEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkInversionBendEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkTorsionEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkImproperTorsionEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBondBondEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBondBendEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBendBendEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBondTorsionEnergy(int flag,int f1,int atom_id);
REAL CalculateFrameworkBendTorsionEnergy(int flag,int f1,int atom_id);

int CalculateFrameworkIntraVDWEnergy(void);
int CalculateFrameworkIntraReplicaVDWEnergy(void);
int CalculateFrameworkIntraChargeChargeEnergy(void);
int CalculateFrameworkIntraReplicaChargeChargeEnergy(void);
int CalculateFrameworkIntraChargeBondDipoleEnergy(void);
int CalculateFrameworkIntraBondDipoleBondDipoleEnergy(void);

void CalculateFrameworkAdsorbateVDWEnergy(void);
void CalculateFrameworkAdsorbateReplicaVDWEnergy(void);
void CalculateFrameworkCationVDWEnergy(void);
void CalculateFrameworkCationReplicaVDWEnergy(void);

int CalculateFrameworkAdsorbateChargeChargeEnergy(void);
int CalculateFrameworkAdsorbateReplicaChargeChargeEnergy(void);
int CalculateFrameworkCationChargeChargeEnergy(void);
int CalculateFrameworkCationReplicaChargeChargeEnergy(void);

int CalculateFrameworkAdsorbateChargeBondDipoleEnergy(void);
int CalculateFrameworkCationChargeBondDipoleEnergy(void);

int CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergy(void);
int CalculateFrameworkCationBondDipoleBondDipoleEnergy(void);

void CalculateFrameworkShiftEnergyDifferenceAdsorbateVDW(void);
void CalculateFrameworkShiftEnergyDifferenceAdsorbateCharge(void);
void CalculateFrameworkShiftEnergyDifferenceCationVDW(void);
void CalculateFrameworkShiftEnergyDifferenceCationCharge(void);

REAL CalculateEnergyDifferenceFrameworkMoveVDW(int atom_id,VECTOR posA,int typeA);
void CalculateEnergyDifferenceFrameworkMoveCharge(int index);
void CalculateFrameworkEnergyDifferenceShiftedFramework(void);

int CalculateFrameworkAdsorbateVDWEnergyDifference(int m);
int CalculateFrameworkCationVDWEnergyDifference(int m);

void CalculateFrameworkAdsorbateChargeChargeEnergyDifference(int m);
void CalculateFrameworkCationChargeChargeEnergyDifference(int m);

int CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference(int m);
int CalculateFrameworkCationChargeBondDipoleEnergyDifference(int m);

int CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference(int m);
int CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference(int m);

REAL CalculateFrameworkElectrostaticPotential(POINT posA);

REAL CalculateFrameworkVDWEnergyCorrection(VECTOR* Positions,VECTOR *AnisotropicPositions);

#endif
