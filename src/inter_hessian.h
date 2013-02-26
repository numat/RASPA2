/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'inter_hessian.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef INTER_HESSIAN_H
#define INTER_HESSIAN_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

void PreComputeRotationDerivatives(void);
void HessianAtomicPositionPosition(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,
                                                 REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor);
void HessianCenterOfMassOrientation(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_i2,INT_VECTOR3 index_j,INT_VECTOR3 index_j2,
              int index1,int index2,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor);
void HessianOrientationOrientation(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_i2,INT_VECTOR3 index_j,INT_VECTOR3 index_j2,
             int index1,int index2,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor);

void HessianCenterOfMassStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,
      REAL f1,REAL f2,VECTOR dr,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB,int RigidA,int RigidB);
void HessianOrientationStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i2,INT_VECTOR3 index_j2,int index1,int index2,
                  REAL f1,REAL f2,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB,VECTOR dr);

void HessianAtomicStrainStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,
                       REAL f1,REAL f2,VECTOR dr,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB);


void GradientStrainI(REAL *Gradient,REAL f1,VECTOR dr,VECTOR posA,VECTOR comA);
void GradientStrainJ(REAL *Gradient,REAL f1,VECTOR dr,VECTOR posB,VECTOR comB);


void ComputeInterVDWMolecularHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeInterChargeChargeMolecularHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);

void CalculateBondConstraintExclusionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);

#endif

