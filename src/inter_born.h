/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'inter_born.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef INTER_BORN_H
#define INTER_BORN_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

void ComputeInterVDWBornTerm(void);
int ComputeInterChargeChargeBornTerm(void);

#endif

