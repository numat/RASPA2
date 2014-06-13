/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'integration.c' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef INTEGRATION_H
#define INTEGRATION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

void ComputeKineticEnergySystem(void);

void CalculateEnergy(void);
void CalculateTotalEnergyAllSystems(void);

void InitializeForces(void);
void InitializeForcesAllSystems(void);
void CalculateForce(void);

void CalculateMolecularExcessPressure(void);

void Integration(void);

int CalculateElectricField(void);

void IntegrationLeapFrog(void);

void IntegrationLeapFrogAtomicPositions(void);
void IntegrationLeapFrogAtomicVelocities(void);

void IntegrationLeapFrogAtomicPositionsAndQuaternions(void);
void IntegrationLeapFrogAtomicVelocitiesAndQuaternions(void);

void TestElectrostaticPotential(void);

#endif
