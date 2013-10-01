/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_hessian.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef FRAMEWORK_HESSIAN_H
#define FRAMEWORK_HESSIAN_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <molecule.h>

void CalculateDerivativesAtPositionVDW(VECTOR,int,REAL *,VECTOR *,REAL_MATRIX3x3 *,REAL *);
REAL CalculateDerivativesAtPositionReal(VECTOR,int,REAL *,VECTOR *,REAL_MATRIX3x3 *,REAL *);

void ComputeFrameworkAdsorbateVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkAdsorbateChargeChargeHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkCationVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkCationChargeChargeHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);

int ComputeFrameworkIntraVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
int ComputeFrameworkIntraChargeChargeHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkBondHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkUreyBradleyHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkImproperTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeFrameworkInversionBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivativeTensor,int ComputeGradient,int ComputeHessian);

void ComputeFrameworkBendTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivativeTensor,int ComputeGradient,int ComputeHessian);

#endif
