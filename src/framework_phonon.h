/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_phonon.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef FRAMEWORK_PHONON_H
#define FRAMEWORK_PHONON_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <molecule.h>

void ComputeFrameworkAdsorbateVDWPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkAdsorbateChargeChargePhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkCationVDWPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkCationChargeChargePhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);

int ComputeFrameworkIntraVDWPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
int ComputeFrameworkIntraChargeChargePhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkBondPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkBendPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkTorsionPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkImproperTorsionPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkInversionBendPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivativeTensor,int ComputeGradient,int ComputeHessian);

#endif
