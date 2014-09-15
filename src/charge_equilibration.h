/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'charge_equilibration.c' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef CHARGE_EQUILIBRATION_H
#define CHARGE_EQUILIBRATION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

enum{DING_YAZAYDIN_DUBBELDAM,RAPPE_GODDARD,WILMER_SNURR};

extern int ChargeEquilibrationPeriodic;
extern int ChargeEquilibrationEwald;
void ChargeEquilibration(void);

#endif
