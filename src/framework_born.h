/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_born.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef FRAMEWORK_BORN_H
#define FRAMEWORK_BORN_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <molecule.h>

void CalculateFrameworkBondBornTerm(void);
void CalculateFrameworkUreyBradleyBornTerm(void);
void CalculateFrameworkBendBornTerm(void);
void CalculateFrameworkTorsionBornTerm(void);
void CalculateFrameworkImproperTorsionBornTerm(void);

int CalculateFrameworkIntraVDWBornTerm(void);
int CalculateFrameworkIntraReplicaVDWBornTerm(void);
int CalculateFrameworkIntraChargeChargeBornTerm(void);
int CalculateFrameworkIntraReplicaChargeChargeBornTerm(void);

void FrameworkAdsorbateVDWBornTerm(void);
void FrameworkAdsorbateChargeChargeBornTerm(void);
void FrameworkCationVDWBornTerm(void);
void FrameworkCationChargeChargeBornTerm(void);

void FrameworkAdsorbateReplicaVDWBornTerm(void);
void FrameworkAdsorbateReplicaChargeChargeBornTerm(void);
void FrameworkCationReplicaVDWBornTerm(void);
void FrameworkCationReplicaChargeChargeBornTerm(void);
#endif
