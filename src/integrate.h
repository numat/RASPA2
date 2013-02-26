/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'integrate.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "cubic_spline_1d.h"

//double QuadratureSimpson(CUBIC_SPLINE spline,double f(REAL),REAL a,REAL b);

REAL IntegrateGeneralizedSimpsonOrderN(REAL *data,REAL *count,REAL h,int n);

REAL IntegrateGeneralizedSimpsonConventional(REAL *data,REAL h,int n);
REAL IntegrateGeneralizedSimpsonConventionalVectorX(VECTOR *data,REAL delta,int n);
REAL IntegrateGeneralizedSimpsonConventionalVectorY(VECTOR *data,REAL delta,int n);
REAL IntegrateGeneralizedSimpsonConventionalVectorZ(VECTOR *data,REAL delta,int n);

#endif
