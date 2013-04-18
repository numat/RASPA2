/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'numerical.c' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef NUMERICAL_H
#define NUMERICAL_H

void SaveFrameworkPositionsToReferenceValues(void);
void PlaceFrameworkInBoxFromReferenceValues(void);
void SaveAdsorbateAtomPositionsToReferenceValues(void);
void PlaceAdsorbateAtomsInBoxFromReferenceValues(void);
void SaveCationAtomPositionsToReferenceValues(void);
void PlaceCationAtomsInBoxFromReferenceValues(void);

REAL_MATRIX3x3 ComputeStrainDerivativeNumerically(void);
REAL_MATRIX9x9 ComputeStrainSecondDerivativeNumerically(void);
void ComputeHessianMatrixNumerically(REAL_MATRIX HessianMatrix);
void ConstructCrossTermNumerically(void);
void ComputeCrossTermNumerically(REAL_MATRIX CrossTerm);

void ComputeGradientsNumerically(REAL *Gradients);

void CheckStatusNumerically(void);
REAL_MATRIX9x9 ComputeRelaxationTerm(int NumberOfPositionVariables,int NumberOfBoxVariables);
void ConstrucGeneralizedHessianMatrix(REAL_MATRIX GeneralizedHessianMatrix);
void ComputeNormalModeDerivativeNumerically(REAL_MATRIX GeneralizedHessianMatrix,REAL *Eigenvalues,REAL *Positions,
               CUBIC_SPLINE *splines,REAL_MATRIX3x3 StoredBox,REAL_MATRIX3x3 StoredInverseBox);

void ComputeCrossTermNumericallyMinimalSet(REAL_MATRIX CrossTerm);

void ComputeThirdOrderElasticConstantsNumerically(int NumberOfPositionVariables,int NumberOfBoxVariables,REAL_MATRIX6x6x6 *VoigtMatrixThirdOrder);

void ComputeEnergyGradientHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX3x3 *Strain_Derivative_Tensor,
                                  REAL_MATRIX HessianMatrix,REAL_MATRIX CrossTerm);


void AddRemainderOfCrossTermNumerically(REAL_MATRIX HessianMatrix);
void AddRemainderOfBornTermNumerically(REAL_MATRIX HessianMatrix);

 #endif
