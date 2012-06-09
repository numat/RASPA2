/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'inter_energy.h' is part of RASPA.

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

#ifndef INTER_ENERGY_H
#define INTER_ENERGY_H

#include "utils.h"
#include "simulation.h"

REAL CalculateInterVDWEnergyAdsorbateAtPosition(POINT posA,int typeA,int exclude);
int CalculateInterChargeEnergyAdsorbateAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,int exclude);
int CalculateInterBondDipoleEnergyAdsorbateAtPosition(POINT posA1,POINT posA2,REAL DipoleMagnitudeA,REAL *UChargeBondDipole,REAL *UBondDipoleBondDipole,int exclude);

REAL CalculateInterVDWEnergyCationAtPosition(POINT posA,int typeA,int exclude);
int CalculateInterChargeEnergyCationAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,int exclude);
int CalculateInterBondDipoleEnergyCationAtPosition(POINT posA1,POINT posA2,REAL DipoleMagnitudeA,REAL *UChargeBondDipole,REAL *UBondDipoleBondDipole,int exclude);

int CalculateTotalInterVDWEnergy(void);
int CalculateTotalInterReplicaVDWEnergy(void);
int CalculateTotalInterChargeChargeCoulombEnergy(void);
int CalculateTotalInterReplicaChargeChargeCoulombEnergy(void);
int CalculateTotalInterChargeBondDipoleCoulombEnergy(void);
int CalculateTotalInterBondDipoleBondDipoleCoulombEnergy(void);

int CalculateInterVDWEnergyDifferenceAdsorbate(int m);
int CalculateInterVDWEnergyDifferenceCation(int m);

int CalculateInterChargeChargeEnergyDifferenceAdsorbate(int m);
int CalculateInterChargeChargeEnergyDifferenceCation(int m);

int CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(int m);
int CalculateInterChargeBondDipoleEnergyDifferenceCation(int m);

int CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(int m);
int CalculateInterBondDipoleBondDipoleEnergyDifferenceCation(int m);

REAL CalculateInterChargeElectrostaticPotential(POINT posA);

REAL CalculateInterVDWEnergyCorrectionAdsorbate(VECTOR* Positions,VECTOR *AnisotropicPositions,int exclude);
REAL CalculateInterVDWEnergyCorrectionCation(VECTOR* Positions,VECTOR *AnisotropicPositions,int exclude);

#endif
