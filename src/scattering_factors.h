/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'scattering_factors.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef SCATTERING_FACTORS_H
#define SCATTERING_FACTORS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

REAL GetAtomicMass(char *Name);
REAL GetCovalentRadius(char *Name);
REAL GetCovalentRadiusExtended(int type,char *Name);
REAL GetAtomPolarization(char *Name);

void SetDiffractionWaveLength(REAL Wavelength);
REAL ScatteringFactor(int index,REAL stol2);
COMPLEX GetAnomalousScatteringFactor(int index);
int GetScatteringNumber(char *Name);
int GetAnomalousScatteringNumber(char *Name);

#endif
