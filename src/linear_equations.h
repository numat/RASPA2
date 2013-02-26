/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'linear_equations.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef LINEAR_EQUATIONS_H
#define LINEAR_EQUATIONS_H

#include "matrix.h"

int Inverse(REAL_MATRIX c);
int SolveLinearSystem(REAL_FORTRAN_MATRIX *a,REAL_FORTRAN_MATRIX *b);

int trdiag(int n,REAL *lower,REAL *diag,REAL *upper,REAL *b,int rep);
int tzdiag(int n,REAL *lower,REAL *diag,REAL *upper,
           REAL *lowrow,REAL *ricol,REAL *b,int rep);

#endif
