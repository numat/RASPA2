/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_phonon.h' is part of RASPA.

    RASPA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RASPA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************************************/

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
