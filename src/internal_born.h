/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'internal_born.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef INTERNAL_BORN_H
#define INTERNAL_BORN_H

void CalculateAdsorbateBondBornTerm(void);
void CalculateAdsorbateUreyBradleyBornTerm(void);
void CalculateAdsorbateBendBornTerm(void);
void CalculateAdsorbateTorsionBornTerm(void);
void CalculateAdsorbateImproperTorsionBornTerm(void);
void CalculateAdsorbateIntraVDWBornTerm(void);
void CalculateAdsorbateIntraCoulombBornTerm(void);

void CalculateCationBondBornTerm(void);
void CalculateCationUreyBradleyBornTerm(void);
void CalculateCationBendBornTerm(void);
void CalculateCationTorsionBornTerm(void);
void CalculateCationImproperTorsionBornTerm(void);
void CalculateCationIntraVDWBornTerm(void);
void CalculateCationIntraCoulombBornTerm(void);

#endif
