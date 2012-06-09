/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_hessian.c' is part of RASPA.

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "simulation.h"
#include "ewald.h"
#include "potentials.h"
#include "utils.h"
#include "output.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "spacegroup.h"
#include "spectra.h"
#include "rigid.h"
#include "minimization.h"


void CalculateDerivativesAtPositionVDW(VECTOR pos,int typeA,REAL *value,VECTOR *first_derivative,
                                       REAL_MATRIX3x3 *second_derivative,REAL *third_derivative)
{
  int i,f;
  int typeB;
  VECTOR dr;
  REAL rr,F,DF,DDF,DDDF;

  *value=0;
  first_derivative->x=0.0;
  first_derivative->y=0.0;
  first_derivative->z=0.0;
  second_derivative->ax=second_derivative->bx=second_derivative->cx=0.0;
  second_derivative->ay=second_derivative->by=second_derivative->cy=0.0;
  second_derivative->az=second_derivative->bz=second_derivative->cz=0.0;
  *third_derivative=0.0;

  for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f];i++)
    {
      typeB=Framework[CurrentSystem].Atoms[f][i].Type;
      dr.x=pos.x-Framework[CurrentSystem].Atoms[f][i].AnisotropicPosition.x;
      dr.y=pos.y-Framework[CurrentSystem].Atoms[f][i].AnisotropicPosition.y;
      dr.z=pos.z-Framework[CurrentSystem].Atoms[f][i].AnisotropicPosition.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      if(rr<CutOffVDWSquared)
      {
        PotentialThirdDerivative(typeA,typeB,rr,&F,&DF,&DDF,&DDDF);

        *value+=F;

        first_derivative->x+=dr.x*DF;
        first_derivative->y+=dr.y*DF;
        first_derivative->z+=dr.z*DF;

        // add contribution to the second derivatives (Hessian matrix)
        second_derivative->ax+=DDF*dr.x*dr.x+DF; second_derivative->bx+=DDF*dr.y*dr.x;    second_derivative->cx+=DDF*dr.z*dr.x;
        second_derivative->ay+=DDF*dr.x*dr.y;    second_derivative->by+=DDF*dr.y*dr.y+DF; second_derivative->cy+=DDF*dr.z*dr.y;
        second_derivative->az+=DDF*dr.x*dr.z;    second_derivative->bz+=DDF*dr.y*dr.z;    second_derivative->cz+=DDF*dr.z*dr.z+DF;

        *third_derivative+=DDDF*dr.x*dr.y*dr.z;
      }
    }
  }
}

void CalculateDerivativesAtPositionReal(VECTOR pos,int typeA,REAL *value,VECTOR *first_derivative,
                                       REAL_MATRIX3x3 *second_derivative,REAL *third_derivative)
{
  int i,f;
  int typeB;
  VECTOR dr;
  REAL ChargeA,ChargeB;
  REAL r,rr,F,DF,DDF,DDDF;

  *value=0;
  first_derivative->x=0.0;
  first_derivative->y=0.0;
  first_derivative->z=0.0;
  second_derivative->ax=second_derivative->bx=second_derivative->cx=0.0;
  second_derivative->ay=second_derivative->by=second_derivative->cy=0.0;
  second_derivative->az=second_derivative->bz=second_derivative->cz=0.0;
  *third_derivative=0.0;

  ChargeA=1.0;
  for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f];i++)
    {
      typeB=Framework[CurrentSystem].Atoms[f][i].Type;
      ChargeB=Framework[CurrentSystem].Atoms[f][i].Charge;
      dr.x=pos.x-Framework[CurrentSystem].Atoms[f][i].Position.x;
      dr.y=pos.y-Framework[CurrentSystem].Atoms[f][i].Position.y;
      dr.z=pos.z-Framework[CurrentSystem].Atoms[f][i].Position.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      if(rr<CutOffChargeChargeSquared)
      {
        r=sqrt(rr);

        F=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(erfc(Alpha[CurrentSystem]*r)/r);

        DF=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
            (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)+erfc(Alpha[CurrentSystem]*r))/
            (r*rr);

        DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
             (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*r*(3.0+2.0*SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)+
              3.0*erfc(Alpha[CurrentSystem]*r))/(rr*rr*r);

        DDDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
             (-2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*r*(15.0+10.0*SQR(Alpha[CurrentSystem])*rr+4.0*SQR(SQR(Alpha[CurrentSystem])*rr))/sqrt(M_PI)
              -15.0*erfc(Alpha[CurrentSystem]*r))/(rr*rr*rr*r);

        *value+=F;

        first_derivative->x+=dr.x*DF;
        first_derivative->y+=dr.y*DF;
        first_derivative->z+=dr.z*DF;

        // add contribution to the second derivatives (Hessian matrix)
        second_derivative->ax+=DDF*dr.x*dr.x+DF; second_derivative->bx+=DDF*dr.y*dr.x;    second_derivative->cx+=DDF*dr.z*dr.x;
        second_derivative->ay+=DDF*dr.x*dr.y;    second_derivative->by+=DDF*dr.y*dr.y+DF; second_derivative->cy+=DDF*dr.z*dr.y;
        second_derivative->az+=DDF*dr.x*dr.z;    second_derivative->bz+=DDF*dr.y*dr.z;    second_derivative->cz+=DDF*dr.z*dr.z+DF;

        *third_derivative+=DDDF*dr.x*dr.y*dr.z;
      }
    }
  }
}


// Hessian: Center of mass - Center of mass
// ========================================
static inline void HessianAtomicPositionPosition(REAL_MATRIX HessianMatrix,int index_i,int index_j,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor)
{
  REAL_MATRIX3x3 Hessian;

  Hessian.ax=f2*dr.x*dr.x+f1; Hessian.bx=f2*dr.y*dr.x;    Hessian.cx=f2*dr.z*dr.x;
  Hessian.ay=f2*dr.x*dr.y;    Hessian.by=f2*dr.y*dr.y+f1; Hessian.cy=f2*dr.z*dr.y;
  Hessian.az=f2*dr.x*dr.z;    Hessian.bz=f2*dr.y*dr.z;    Hessian.cz=f2*dr.z*dr.z+f1;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  { 
    if(index_i>=0)
    {
      HessianMatrix.element[index_i][index_i]+=ReplicaFactor*Hessian.ax;
      HessianMatrix.element[index_i][index_i+1]+=ReplicaFactor*Hessian.ay;
      HessianMatrix.element[index_i][index_i+2]+=ReplicaFactor*Hessian.az;
      HessianMatrix.element[index_i+1][index_i+1]+=ReplicaFactor*Hessian.by;
      HessianMatrix.element[index_i+1][index_i+2]+=ReplicaFactor*Hessian.bz;
      HessianMatrix.element[index_i+2][index_i+2]+=ReplicaFactor*Hessian.cz;
    }
  }

  if(index_j>=0)
  { 
    HessianMatrix.element[index_j][index_j]+=ReplicaFactor*Hessian.ax;
    HessianMatrix.element[index_j][index_j+1]+=ReplicaFactor*Hessian.ay;
    HessianMatrix.element[index_j][index_j+2]+=ReplicaFactor*Hessian.az;
    HessianMatrix.element[index_j+1][index_j+1]+=ReplicaFactor*Hessian.by;
    HessianMatrix.element[index_j+1][index_j+2]+=ReplicaFactor*Hessian.bz;
    HessianMatrix.element[index_j+2][index_j+2]+=ReplicaFactor*Hessian.cz;
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    if((index_i>=0)&&(index_j>=0))
    {
      HessianMatrix.element[index_i][index_j]-=Hessian.ax;
      HessianMatrix.element[index_i][index_j+1]-=Hessian.ay;
      HessianMatrix.element[index_i][index_j+2]-=Hessian.az;
      HessianMatrix.element[index_i+1][index_j]-=Hessian.ay;
      HessianMatrix.element[index_i+1][index_j+1]-=Hessian.by;
      HessianMatrix.element[index_i+1][index_j+2]-=Hessian.bz;
      HessianMatrix.element[index_i+2][index_j]-=Hessian.az;
      HessianMatrix.element[index_i+2][index_j+1]-=Hessian.bz;
      HessianMatrix.element[index_i+2][index_j+2]-=Hessian.cz;
    }
  }
}


// Hessian: Center of mass - Orientation
// =====================================
static inline void HessianCenterOfMassOrientationJ(REAL_MATRIX HessianMatrix,int index_i,int index_j,int index_j2,int index2,REAL f1,REAL f2,VECTOR dr)
{
  REAL_MATRIX3x3 Hessian;
  VECTOR vecj1,vecj2,vecj3;

  vecj1=DVecX[index2];
  vecj2=DVecY[index2];
  vecj3=DVecZ[index2];

  Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
  Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
  Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

  Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
  Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
  Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

  Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
  Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
  Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

  if((index_j>=0)&&(index_j2>=0))
  {
    HessianMatrix.element[index_j][index_j2]+=Hessian.ax;
    HessianMatrix.element[index_j][index_j2+1]+=Hessian.ay;
    HessianMatrix.element[index_j][index_j2+2]+=Hessian.az;

    HessianMatrix.element[index_j+1][index_j2]+=Hessian.bx;
    HessianMatrix.element[index_j+1][index_j2+1]+=Hessian.by;
    HessianMatrix.element[index_j+1][index_j2+2]+=Hessian.bz;

    HessianMatrix.element[index_j+2][index_j2]+=Hessian.cx;
    HessianMatrix.element[index_j+2][index_j2+1]+=Hessian.cy;
    HessianMatrix.element[index_j+2][index_j2+2]+=Hessian.cz;
  }

  Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
  Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
  Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

  Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
  Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
  Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

  Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
  Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
  Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    if((index_i>=0)&&(index_j2>=0))
    {
      HessianMatrix.element[index_i][index_j2]-=Hessian.ax;
      HessianMatrix.element[index_i][index_j2+1]-=Hessian.ay;
      HessianMatrix.element[index_i][index_j2+2]-=Hessian.az;

      HessianMatrix.element[index_i+1][index_j2]-=Hessian.bx;
      HessianMatrix.element[index_i+1][index_j2+1]-=Hessian.by;
      HessianMatrix.element[index_i+1][index_j2+2]-=Hessian.bz;

      HessianMatrix.element[index_i+2][index_j2]-=Hessian.cx;
      HessianMatrix.element[index_i+2][index_j2+1]-=Hessian.cy;
      HessianMatrix.element[index_i+2][index_j2+2]-=Hessian.cz;
    }
  }
}

// Hessian: Orientation - Orientation
// ==================================
static inline void HessianOrientationOrientationJ(REAL_MATRIX HessianMatrix,int index_j2,int index2,REAL f1,REAL f2,VECTOR dr)
{
  REAL_MATRIX3x3 Hessian;
  VECTOR vecj1,vecj2,vecj3;
  VECTOR DDvecJAX,DDvecJBY,DDvecJCZ,DDvecJAY,DDvecJAZ,DDvecJBZ;

  vecj1=DVecX[index2];
  vecj2=DVecY[index2];
  vecj3=DVecZ[index2];

  DDvecJAX=DDVecAX[index2];
  DDvecJBY=DDVecBY[index2];
  DDvecJCZ=DDVecCZ[index2];
  DDvecJAY=DDVecAY[index2];
  DDvecJAZ=DDVecAZ[index2];
  DDvecJBZ=DDVecBZ[index2];

  if(index_j2>=0)
  {
    Hessian.ax=-f1*(dr.x*DDvecJAX.x+dr.y*DDvecJAX.y+dr.z*DDvecJAX.z)
               +f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)
               +f1*(vecj1.x*vecj1.x+vecj1.y*vecj1.y+vecj1.z*vecj1.z);

    Hessian.by=-f1*(dr.x*DDvecJBY.x+dr.y*DDvecJBY.y+dr.z*DDvecJBY.z)+
               f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)+
               f1*(vecj2.x*vecj2.x+vecj2.y*vecj2.y+vecj2.z*vecj2.z);

    Hessian.cz=-f1*(dr.x*DDvecJCZ.x+dr.y*DDvecJCZ.y+dr.z*DDvecJCZ.z)+
               f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
               f1*(vecj3.x*vecj3.x+vecj3.y*vecj3.y+vecj3.z*vecj3.z);

    Hessian.ay=-f1*(dr.x*DDvecJAY.x+dr.y*DDvecJAY.y+dr.z*DDvecJAY.z)+
               f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)+
               f1*(vecj1.x*vecj2.x+vecj1.y*vecj2.y+vecj1.z*vecj2.z);

    Hessian.az=-f1*(dr.x*DDvecJAZ.x+dr.y*DDvecJAZ.y+dr.z*DDvecJAZ.z)+
               f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
               f1*(vecj1.x*vecj3.x+vecj1.y*vecj3.y+vecj1.z*vecj3.z);

    Hessian.bz=-f1*(dr.x*DDvecJBZ.x+dr.y*DDvecJBZ.y+dr.z*DDvecJBZ.z)+
               f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
               f1*(vecj2.x*vecj3.x+vecj2.y*vecj3.y+vecj2.z*vecj3.z);

    HessianMatrix.element[index_j2][index_j2]+=Hessian.ax;
    HessianMatrix.element[index_j2+1][index_j2+1]+=Hessian.by;
    HessianMatrix.element[index_j2+2][index_j2+2]+=Hessian.cz;
    HessianMatrix.element[index_j2][index_j2+1]+=Hessian.ay;
    HessianMatrix.element[index_j2][index_j2+2]+=Hessian.az;
    HessianMatrix.element[index_j2+1][index_j2+2]+=Hessian.bz;
  }
}

// Hessian: Center of mass - Strain (part I)
// =========================================
static inline void HessianAtomicPositionStrain(REAL_MATRIX HessianMatrix,int index_i,int index_j,REAL f1,REAL f2,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        if(index_i>=0)
        {
          HessianMatrix.element[index_i][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;     // xx x + yy x + zz x
          HessianMatrix.element[index_i+1][n]+=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;   // xx y + yy y + zz y
          HessianMatrix.element[index_i+2][n]+=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;   // xx z + yy z + zz z
        }
      }

      if(index_j>=0)
      {
        HessianMatrix.element[index_j][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;
        HessianMatrix.element[index_j+1][n]-=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;
        HessianMatrix.element[index_j+2][n]-=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            if(index_i>=0)
            {
              HessianMatrix.element[index_i][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              HessianMatrix.element[index_i][n+1]+=f2*dr.y*dr.y*dr.x;
              HessianMatrix.element[index_i][n+2]+=f2*dr.z*dr.z*dr.x;

              HessianMatrix.element[index_i+1][n]+=f2*dr.x*dr.x*dr.y;
              HessianMatrix.element[index_i+1][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
              HessianMatrix.element[index_i+1][n+2]+=f2*dr.z*dr.z*dr.y;

              HessianMatrix.element[index_i+2][n]+=f2*dr.x*dr.x*dr.z;
              HessianMatrix.element[index_i+2][n+1]+=f2*dr.y*dr.y*dr.z;
              HessianMatrix.element[index_i+2][n+2]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
            }
          }

          if(index_j>=0)
          {
            HessianMatrix.element[index_j][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
            HessianMatrix.element[index_j][n+1]-=f2*dr.y*dr.y*dr.x;
            HessianMatrix.element[index_j][n+2]-=f2*dr.z*dr.z*dr.x;

            HessianMatrix.element[index_j+1][n]-=f2*dr.x*dr.x*dr.y;
            HessianMatrix.element[index_j+1][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
            HessianMatrix.element[index_j+1][n+2]-=f2*dr.z*dr.z*dr.y;

            HessianMatrix.element[index_j+2][n]-=f2*dr.x*dr.x*dr.z;
            HessianMatrix.element[index_j+2][n+1]-=f2*dr.y*dr.y*dr.z;
            HessianMatrix.element[index_j+2][n+2]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
          }
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            if(index_i>=0)
            {
              HessianMatrix.element[index_i][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
              HessianMatrix.element[index_i][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;
              HessianMatrix.element[index_i][n+2]+=f2*dr.x*dr.z*dr.x+f1*dr.z;
              HessianMatrix.element[index_i][n+3]+=f2*dr.y*dr.y*dr.x;
              HessianMatrix.element[index_i][n+4]+=f2*dr.y*dr.z*dr.x;
              HessianMatrix.element[index_i][n+5]+=f2*dr.z*dr.z*dr.x;

              HessianMatrix.element[index_i+1][n]+=f2*dr.x*dr.x*dr.y;
              HessianMatrix.element[index_i+1][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;
              HessianMatrix.element[index_i+1][n+2]+=f2*dr.x*dr.z*dr.y;
              HessianMatrix.element[index_i+1][n+3]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
              HessianMatrix.element[index_i+1][n+4]+=f2*dr.y*dr.z*dr.y+f1*dr.z;
              HessianMatrix.element[index_i+1][n+5]+=f2*dr.z*dr.z*dr.y;

              HessianMatrix.element[index_i+2][n]+=f2*dr.x*dr.x*dr.z;
              HessianMatrix.element[index_i+2][n+1]+=f2*dr.x*dr.y*dr.z;
              HessianMatrix.element[index_i+2][n+2]+=f2*dr.x*dr.z*dr.z+f1*dr.x;
              HessianMatrix.element[index_i+2][n+3]+=f2*dr.y*dr.y*dr.z;
              HessianMatrix.element[index_i+2][n+4]+=f2*dr.y*dr.z*dr.z+f1*dr.y;
              HessianMatrix.element[index_i+2][n+5]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
            }
          }

          if(index_j>=0)
          {
            HessianMatrix.element[index_j][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
            HessianMatrix.element[index_j][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;
            HessianMatrix.element[index_j][n+2]-=f2*dr.x*dr.z*dr.x+f1*dr.z;
            HessianMatrix.element[index_j][n+3]-=f2*dr.y*dr.y*dr.x;
            HessianMatrix.element[index_j][n+4]-=f2*dr.y*dr.z*dr.x;
            HessianMatrix.element[index_j][n+5]-=f2*dr.z*dr.z*dr.x;

            HessianMatrix.element[index_j+1][n]-=f2*dr.x*dr.x*dr.y;
            HessianMatrix.element[index_j+1][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;
            HessianMatrix.element[index_j+1][n+2]-=f2*dr.x*dr.z*dr.y;
            HessianMatrix.element[index_j+1][n+3]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
            HessianMatrix.element[index_j+1][n+4]-=f2*dr.y*dr.z*dr.y+f1*dr.z;
            HessianMatrix.element[index_j+1][n+5]-=f2*dr.z*dr.z*dr.y;

            HessianMatrix.element[index_j+2][n]-=f2*dr.x*dr.x*dr.z;
            HessianMatrix.element[index_j+2][n+1]-=f2*dr.x*dr.y*dr.z;
            HessianMatrix.element[index_j+2][n+2]-=f2*dr.x*dr.z*dr.z+f1*dr.x;
            HessianMatrix.element[index_j+2][n+3]-=f2*dr.y*dr.y*dr.z;
            HessianMatrix.element[index_j+2][n+4]-=f2*dr.y*dr.z*dr.z+f1*dr.y;
            HessianMatrix.element[index_j+2][n+5]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
          }
          break;                                   
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i>=0)
                {
                  HessianMatrix.element[index_i][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
                  HessianMatrix.element[index_i][n+1]+=f2*dr.y*dr.y*dr.x;
                  HessianMatrix.element[index_i][n+2]+=f2*dr.y*dr.z*dr.x;
                  HessianMatrix.element[index_i][n+3]+=f2*dr.z*dr.z*dr.x;

                  HessianMatrix.element[index_i+1][n]+=f2*dr.x*dr.x*dr.y;
                  HessianMatrix.element[index_i+1][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
                  HessianMatrix.element[index_i+1][n+2]+=f2*dr.y*dr.z*dr.y+f1*dr.z;
                  HessianMatrix.element[index_i+1][n+3]+=f2*dr.z*dr.z*dr.y;

                  HessianMatrix.element[index_i+2][n]+=f2*dr.x*dr.x*dr.z;
                  HessianMatrix.element[index_i+2][n+1]+=f2*dr.y*dr.y*dr.z;
                  HessianMatrix.element[index_i+2][n+2]+=f2*dr.y*dr.z*dr.z+f1*dr.y;
                  HessianMatrix.element[index_i+2][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
                }
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
                HessianMatrix.element[index_j][n+1]-=f2*dr.y*dr.y*dr.x;
                HessianMatrix.element[index_j][n+2]-=f2*dr.y*dr.z*dr.x;
                HessianMatrix.element[index_j][n+3]-=f2*dr.z*dr.z*dr.x;

                HessianMatrix.element[index_j+1][n]-=f2*dr.x*dr.x*dr.y;
                HessianMatrix.element[index_j+1][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
                HessianMatrix.element[index_j+1][n+2]-=f2*dr.y*dr.z*dr.y+f1*dr.z;
                HessianMatrix.element[index_j+1][n+3]-=f2*dr.z*dr.z*dr.y;

                HessianMatrix.element[index_j+2][n]-=f2*dr.x*dr.x*dr.z;
                HessianMatrix.element[index_j+2][n+1]-=f2*dr.y*dr.y*dr.z;
                HessianMatrix.element[index_j+2][n+2]-=f2*dr.y*dr.z*dr.z+f1*dr.y;
                HessianMatrix.element[index_j+2][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i>=0)
                {
                  HessianMatrix.element[index_i][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
                  HessianMatrix.element[index_i][n+1]+=f2*dr.x*dr.z*dr.x+f1*dr.z;
                  HessianMatrix.element[index_i][n+2]+=f2*dr.y*dr.y*dr.x;
                  HessianMatrix.element[index_i][n+3]+=f2*dr.z*dr.z*dr.x;

                  HessianMatrix.element[index_i+1][n]+=f2*dr.x*dr.x*dr.y;
                  HessianMatrix.element[index_i+1][n+1]+=f2*dr.x*dr.z*dr.y;
                  HessianMatrix.element[index_i+1][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
                  HessianMatrix.element[index_i+1][n+3]+=f2*dr.z*dr.z*dr.y;

                  HessianMatrix.element[index_i+2][n]+=f2*dr.x*dr.x*dr.z;
                  HessianMatrix.element[index_i+2][n+1]+=f2*dr.x*dr.z*dr.z+f1*dr.x;
                  HessianMatrix.element[index_i+2][n+2]+=f2*dr.y*dr.y*dr.z;
                  HessianMatrix.element[index_i+2][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
                }
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
                HessianMatrix.element[index_j][n+1]-=f2*dr.x*dr.z*dr.x+f1*dr.z;
                HessianMatrix.element[index_j][n+2]-=f2*dr.y*dr.y*dr.x;
                HessianMatrix.element[index_j][n+3]-=f2*dr.z*dr.z*dr.x;

                HessianMatrix.element[index_j+1][n]-=f2*dr.x*dr.x*dr.y;
                HessianMatrix.element[index_j+1][n+1]-=f2*dr.x*dr.z*dr.y;
                HessianMatrix.element[index_j+1][n+2]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
                HessianMatrix.element[index_j+1][n+3]-=f2*dr.z*dr.z*dr.y;

                HessianMatrix.element[index_j+2][n]-=f2*dr.x*dr.x*dr.z;
                HessianMatrix.element[index_j+2][n+1]-=f2*dr.x*dr.z*dr.z+f1*dr.x;
                HessianMatrix.element[index_j+2][n+2]-=f2*dr.y*dr.y*dr.z;
                HessianMatrix.element[index_j+2][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i>=0)
                {
                  HessianMatrix.element[index_i][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
                  HessianMatrix.element[index_i][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;
                  HessianMatrix.element[index_i][n+2]+=f2*dr.y*dr.y*dr.x;
                  HessianMatrix.element[index_i][n+3]+=f2*dr.z*dr.z*dr.x;

                  HessianMatrix.element[index_i+1][n]+=f2*dr.x*dr.x*dr.y;
                  HessianMatrix.element[index_i+1][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;
                  HessianMatrix.element[index_i+1][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
                  HessianMatrix.element[index_i+1][n+3]+=f2*dr.z*dr.z*dr.y;

                  HessianMatrix.element[index_i+2][n]+=f2*dr.x*dr.x*dr.z;
                  HessianMatrix.element[index_i+2][n+1]+=f2*dr.x*dr.y*dr.z;
                  HessianMatrix.element[index_i+2][n+2]+=f2*dr.y*dr.y*dr.z;
                  HessianMatrix.element[index_i+2][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
                }
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
                HessianMatrix.element[index_j][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;
                HessianMatrix.element[index_j][n+2]-=f2*dr.y*dr.y*dr.x;
                HessianMatrix.element[index_j][n+3]-=f2*dr.z*dr.z*dr.x;

                HessianMatrix.element[index_j+1][n]-=f2*dr.x*dr.x*dr.y;
                HessianMatrix.element[index_j+1][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;
                HessianMatrix.element[index_j+1][n+2]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
                HessianMatrix.element[index_j+1][n+3]-=f2*dr.z*dr.z*dr.y;

                HessianMatrix.element[index_j+2][n]-=f2*dr.x*dr.x*dr.z;
                HessianMatrix.element[index_j+2][n+1]-=f2*dr.x*dr.y*dr.z;
                HessianMatrix.element[index_j+2][n+2]-=f2*dr.y*dr.y*dr.z;
                HessianMatrix.element[index_j+2][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
              }
              break;
          }
          break;  
        default:  
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}

// Hessian: Center of mass - Strain (part II)
// ==========================================
static inline void HessianCenterOfMassStrainJ(REAL_MATRIX HessianMatrix,int index_i,int index_j,REAL f1,REAL f2,VECTOR dr,VECTOR posB,VECTOR comB)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        if(index_i>=0)
        {
          HessianMatrix.element[index_i][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1)+(posB.y-comB.y)*(f2*dr.x*dr.y)+(posB.z-comB.z)*(f2*dr.x*dr.z);
          HessianMatrix.element[index_i+1][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x)+(posB.y-comB.y)*(f2*dr.y*dr.y+f1)+(posB.z-comB.z)*(f2*dr.y*dr.z);
          HessianMatrix.element[index_i+2][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x)+(posB.y-comB.y)*(f2*dr.z*dr.y)+(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
        }
      }

      if(index_j>=0)
      {
        HessianMatrix.element[index_j][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1)+(posB.y-comB.y)*(f2*dr.x*dr.y)+(posB.z-comB.z)*(f2*dr.x*dr.z);
        HessianMatrix.element[index_j+1][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x)+(posB.y-comB.y)*(f2*dr.y*dr.y+f1)+(posB.z-comB.z)*(f2*dr.y*dr.z);
        HessianMatrix.element[index_j+2][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x)+(posB.y-comB.y)*(f2*dr.z*dr.y)+(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            if(index_i>=0)
            {
              HessianMatrix.element[index_i][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
              HessianMatrix.element[index_i][n+1]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
              HessianMatrix.element[index_i][n+2]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

              HessianMatrix.element[index_i+1][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
              HessianMatrix.element[index_i+1][n+1]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
              HessianMatrix.element[index_i+1][n+2]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

              HessianMatrix.element[index_i+2][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
              HessianMatrix.element[index_i+2][n+1]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
              HessianMatrix.element[index_i+2][n+2]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
            }
          }

          if(index_j>=0)
          {
            HessianMatrix.element[index_j][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
            HessianMatrix.element[index_j][n+1]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
            HessianMatrix.element[index_j][n+2]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

            HessianMatrix.element[index_j+1][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
            HessianMatrix.element[index_j+1][n+1]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
            HessianMatrix.element[index_j+1][n+2]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

            HessianMatrix.element[index_j+2][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
            HessianMatrix.element[index_j+2][n+1]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
            HessianMatrix.element[index_j+2][n+2]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
          }
          break;
        case REGULAR:
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            if(index_i>=0)
            {
              HessianMatrix.element[index_i][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
              HessianMatrix.element[index_i][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.y));
              HessianMatrix.element[index_i][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.z));
              HessianMatrix.element[index_i][n+3]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
              HessianMatrix.element[index_i][n+4]+=0.5*((posB.z-comB.z)*(f2*dr.x*dr.y)+(posB.y-comB.y)*(f2*dr.x*dr.z));
              HessianMatrix.element[index_i][n+5]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

              HessianMatrix.element[index_i+1][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
              HessianMatrix.element[index_i+1][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.y+f1));
              HessianMatrix.element[index_i+1][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.z));
              HessianMatrix.element[index_i+1][n+3]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
              HessianMatrix.element[index_i+1][n+4]+=0.5*((posB.z-comB.z)*(f2*dr.y*dr.y+f1)+(posB.y-comB.y)*(f2*dr.y*dr.z));
              HessianMatrix.element[index_i+1][n+5]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

              HessianMatrix.element[index_i+2][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
              HessianMatrix.element[index_i+2][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.y));
              HessianMatrix.element[index_i+2][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.z+f1));
              HessianMatrix.element[index_i+2][n+3]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
              HessianMatrix.element[index_i+2][n+4]+=0.5*((posB.z-comB.z)*(f2*dr.z*dr.y)+(posB.y-comB.y)*(f2*dr.z*dr.z+f1));
              HessianMatrix.element[index_i+2][n+5]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
            }
          }

          if(index_j>=0)
          {
            HessianMatrix.element[index_j][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
            HessianMatrix.element[index_j][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.y));
            HessianMatrix.element[index_j][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.z));
            HessianMatrix.element[index_j][n+3]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
            HessianMatrix.element[index_j][n+4]-=0.5*((posB.z-comB.z)*(f2*dr.x*dr.y)+(posB.y-comB.y)*(f2*dr.x*dr.z));
            HessianMatrix.element[index_j][n+5]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

            HessianMatrix.element[index_j+1][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
            HessianMatrix.element[index_j+1][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.y+f1));
            HessianMatrix.element[index_j+1][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.z));
            HessianMatrix.element[index_j+1][n+3]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
            HessianMatrix.element[index_j+1][n+4]-=0.5*((posB.z-comB.z)*(f2*dr.y*dr.y+f1)+(posB.y-comB.y)*(f2*dr.y*dr.z));
            HessianMatrix.element[index_j+1][n+5]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

            HessianMatrix.element[index_j+2][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
            HessianMatrix.element[index_j+2][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.y));
            HessianMatrix.element[index_j+2][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.z+f1));
            HessianMatrix.element[index_j+2][n+3]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
            HessianMatrix.element[index_j+2][n+4]-=0.5*((posB.z-comB.z)*(f2*dr.z*dr.y)+(posB.y-comB.y)*(f2*dr.z*dr.z+f1));
            HessianMatrix.element[index_j+2][n+5]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
          }
          break;                                   
        case REGULAR_UPPER_TRIANGLE:
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            if(index_i>=0)
            {
              HessianMatrix.element[index_i  ][n  ]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
              HessianMatrix.element[index_i  ][n+1]+=(posB.y-comB.y)*(f2*dr.x*dr.x+f1);
              HessianMatrix.element[index_i  ][n+2]+=(posB.z-comB.z)*(f2*dr.x*dr.x+f1);
              HessianMatrix.element[index_i  ][n+3]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
              HessianMatrix.element[index_i  ][n+4]+=(posB.z-comB.z)*(f2*dr.x*dr.y);
              HessianMatrix.element[index_i  ][n+5]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

              HessianMatrix.element[index_i+1][n  ]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
              HessianMatrix.element[index_i+1][n+1]+=(posB.y-comB.y)*(f2*dr.y*dr.x);
              HessianMatrix.element[index_i+1][n+2]+=(posB.z-comB.z)*(f2*dr.y*dr.x);
              HessianMatrix.element[index_i+1][n+3]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
              HessianMatrix.element[index_i+1][n+4]+=(posB.z-comB.z)*(f2*dr.y*dr.y+f1);
              HessianMatrix.element[index_i+1][n+5]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

              HessianMatrix.element[index_i+2][n  ]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
              HessianMatrix.element[index_i+2][n+1]+=(posB.y-comB.y)*(f2*dr.z*dr.x);
              HessianMatrix.element[index_i+2][n+2]+=(posB.z-comB.z)*(f2*dr.z*dr.x);
              HessianMatrix.element[index_i+2][n+3]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
              HessianMatrix.element[index_i+2][n+4]+=(posB.z-comB.z)*(f2*dr.z*dr.y);
              HessianMatrix.element[index_i+2][n+5]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
            }
          }

          if(index_j>=0)
          {
            HessianMatrix.element[index_j  ][n  ]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
            HessianMatrix.element[index_j  ][n+1]-=(posB.y-comB.y)*(f2*dr.x*dr.x+f1);
            HessianMatrix.element[index_j  ][n+2]-=(posB.z-comB.z)*(f2*dr.x*dr.x+f1);
            HessianMatrix.element[index_j  ][n+3]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
            HessianMatrix.element[index_j  ][n+4]-=(posB.z-comB.z)*(f2*dr.x*dr.y);
            HessianMatrix.element[index_j  ][n+5]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

            HessianMatrix.element[index_j+1][n  ]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
            HessianMatrix.element[index_j+1][n+1]-=(posB.y-comB.y)*(f2*dr.y*dr.x);
            HessianMatrix.element[index_j+1][n+2]-=(posB.z-comB.z)*(f2*dr.y*dr.x);
            HessianMatrix.element[index_j+1][n+3]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
            HessianMatrix.element[index_j+1][n+4]-=(posB.z-comB.z)*(f2*dr.y*dr.y+f1);
            HessianMatrix.element[index_j+1][n+5]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

            HessianMatrix.element[index_j+2][n  ]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
            HessianMatrix.element[index_j+2][n+1]-=(posB.y-comB.y)*(f2*dr.z*dr.x);
            HessianMatrix.element[index_j+2][n+2]-=(posB.z-comB.z)*(f2*dr.z*dr.x);
            HessianMatrix.element[index_j+2][n+3]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
            HessianMatrix.element[index_j+2][n+4]-=(posB.z-comB.z)*(f2*dr.z*dr.y);
            HessianMatrix.element[index_j+2][n+5]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
          }
          break;                                   
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i>=0)
                {
                  HessianMatrix.element[index_i][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                  HessianMatrix.element[index_i][n+1]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                  HessianMatrix.element[index_i][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.x*dr.y)+(posB.y-comB.y)*(f2*dr.x*dr.z));
                  HessianMatrix.element[index_i][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

                  HessianMatrix.element[index_i+1][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                  HessianMatrix.element[index_i+1][n+1]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                  HessianMatrix.element[index_i+1][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.y*dr.y+f1)+(posB.y-comB.y)*(f2*dr.y*dr.z));
                  HessianMatrix.element[index_i+1][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

                  HessianMatrix.element[index_i+2][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                  HessianMatrix.element[index_i+2][n+1]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                  HessianMatrix.element[index_i+2][n+2]+=0.5*((posB.z-comB.z)*(f2*dr.z*dr.y)+(posB.y-comB.y)*(f2*dr.z*dr.z+f1));
                  HessianMatrix.element[index_i+2][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
                }
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                HessianMatrix.element[index_j][n+1]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
                HessianMatrix.element[index_j][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.x*dr.y)+(posB.y-comB.y)*(f2*dr.x*dr.z));
                HessianMatrix.element[index_j][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

                HessianMatrix.element[index_j+1][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
                HessianMatrix.element[index_j+1][n+1]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                HessianMatrix.element[index_j+1][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.y*dr.y+f1)+(posB.y-comB.y)*(f2*dr.y*dr.z));
                HessianMatrix.element[index_j+1][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

                HessianMatrix.element[index_j+2][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
                HessianMatrix.element[index_j+2][n+1]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
                HessianMatrix.element[index_j+2][n+2]-=0.5*((posB.z-comB.z)*(f2*dr.z*dr.y)+(posB.y-comB.y)*(f2*dr.z*dr.z+f1));
                HessianMatrix.element[index_j+2][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i>=0)
                {
                  HessianMatrix.element[index_i][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                  HessianMatrix.element[index_i][n+1]+=0.5*((posB.z-comB.z)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.z));
                  HessianMatrix.element[index_i][n+2]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                  HessianMatrix.element[index_i][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

                  HessianMatrix.element[index_i+1][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                  HessianMatrix.element[index_i+1][n+1]+=0.5*((posB.z-comB.z)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.z));
                  HessianMatrix.element[index_i+1][n+2]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                  HessianMatrix.element[index_i+1][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

                  HessianMatrix.element[index_i+2][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                  HessianMatrix.element[index_i+2][n+1]+=0.5*((posB.z-comB.z)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.z+f1));
                  HessianMatrix.element[index_i+2][n+2]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                  HessianMatrix.element[index_i+2][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
                }
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                HessianMatrix.element[index_j][n+1]-=0.5*((posB.z-comB.z)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.z));
                HessianMatrix.element[index_j][n+2]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
                HessianMatrix.element[index_j][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

                HessianMatrix.element[index_j+1][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
                HessianMatrix.element[index_j+1][n+1]-=0.5*((posB.z-comB.z)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.z));
                HessianMatrix.element[index_j+1][n+2]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                HessianMatrix.element[index_j+1][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

                HessianMatrix.element[index_j+2][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
                HessianMatrix.element[index_j+2][n+1]-=0.5*((posB.z-comB.z)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.z+f1));
                HessianMatrix.element[index_j+2][n+2]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
                HessianMatrix.element[index_j+2][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i>=0)
                {
                  HessianMatrix.element[index_i][n]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                  HessianMatrix.element[index_i][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.y));
                  HessianMatrix.element[index_i][n+2]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                  HessianMatrix.element[index_i][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);

                  HessianMatrix.element[index_i+1][n]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                  HessianMatrix.element[index_i+1][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.y+f1));
                  HessianMatrix.element[index_i+1][n+2]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                  HessianMatrix.element[index_i+1][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);

                  HessianMatrix.element[index_i+2][n]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                  HessianMatrix.element[index_i+2][n+1]+=0.5*((posB.y-comB.y)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.y));
                  HessianMatrix.element[index_i+2][n+2]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                  HessianMatrix.element[index_i+2][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
                }
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j][n]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                HessianMatrix.element[index_j][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.x*dr.x+f1)+(posB.x-comB.x)*(f2*dr.x*dr.y));
                HessianMatrix.element[index_j][n+2]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
                HessianMatrix.element[index_j][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);

                HessianMatrix.element[index_j+1][n]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
                HessianMatrix.element[index_j+1][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.y*dr.x)+(posB.x-comB.x)*(f2*dr.y*dr.y+f1));
                HessianMatrix.element[index_j+1][n+2]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                HessianMatrix.element[index_j+1][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);

                HessianMatrix.element[index_j+2][n]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
                HessianMatrix.element[index_j+2][n+1]-=0.5*((posB.y-comB.y)*(f2*dr.z*dr.x)+(posB.x-comB.x)*(f2*dr.z*dr.y));
                HessianMatrix.element[index_j+2][n+2]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
                HessianMatrix.element[index_j+2][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i>=0)
                {
                  HessianMatrix.element[index_i  ][n  ]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                  HessianMatrix.element[index_i  ][n+1]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                  HessianMatrix.element[index_i  ][n+2]+=(posB.z-comB.z)*(f2*dr.x*dr.y);
                  HessianMatrix.element[index_i  ][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);
    
                  HessianMatrix.element[index_i+1][n  ]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                  HessianMatrix.element[index_i+1][n+1]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                  HessianMatrix.element[index_i+1][n+2]+=(posB.z-comB.z)*(f2*dr.y*dr.y+f1);
                  HessianMatrix.element[index_i+1][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);
    
                  HessianMatrix.element[index_i+2][n  ]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                  HessianMatrix.element[index_i+2][n+1]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                  HessianMatrix.element[index_i+2][n+2]+=(posB.z-comB.z)*(f2*dr.z*dr.y);
                  HessianMatrix.element[index_i+2][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
                }
              }
    
              if(index_j>=0)
              {
                HessianMatrix.element[index_j  ][n  ]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                HessianMatrix.element[index_j  ][n+1]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
                HessianMatrix.element[index_j  ][n+2]-=(posB.z-comB.z)*(f2*dr.x*dr.y);
                HessianMatrix.element[index_j  ][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);
    
                HessianMatrix.element[index_j+1][n  ]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
                HessianMatrix.element[index_j+1][n+1]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                HessianMatrix.element[index_j+1][n+2]-=(posB.z-comB.z)*(f2*dr.y*dr.y+f1);
                HessianMatrix.element[index_j+1][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);
    
                HessianMatrix.element[index_j+2][n  ]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
                HessianMatrix.element[index_j+2][n+1]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
                HessianMatrix.element[index_j+2][n+2]-=(posB.z-comB.z)*(f2*dr.z*dr.y);
                HessianMatrix.element[index_j+2][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i>=0)
                {
                  HessianMatrix.element[index_i  ][n  ]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                  HessianMatrix.element[index_i  ][n+1]+=(posB.z-comB.z)*(f2*dr.x*dr.x+f1);
                  HessianMatrix.element[index_i  ][n+2]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                  HessianMatrix.element[index_i  ][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);
    
                  HessianMatrix.element[index_i+1][n  ]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                  HessianMatrix.element[index_i+1][n+1]+=(posB.z-comB.z)*(f2*dr.y*dr.x);
                  HessianMatrix.element[index_i+1][n+2]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                  HessianMatrix.element[index_i+1][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);
    
                  HessianMatrix.element[index_i+2][n  ]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                  HessianMatrix.element[index_i+2][n+1]+=(posB.z-comB.z)*(f2*dr.z*dr.x);
                  HessianMatrix.element[index_i+2][n+2]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                  HessianMatrix.element[index_i+2][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
                }
              }
    
              if(index_j>=0)
              {
                HessianMatrix.element[index_j  ][n  ]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                HessianMatrix.element[index_j  ][n+1]-=(posB.z-comB.z)*(f2*dr.x*dr.x+f1);
                HessianMatrix.element[index_j  ][n+2]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
                HessianMatrix.element[index_j  ][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);
    
                HessianMatrix.element[index_j+1][n  ]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
                HessianMatrix.element[index_j+1][n+1]-=(posB.z-comB.z)*(f2*dr.y*dr.x);
                HessianMatrix.element[index_j+1][n+2]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                HessianMatrix.element[index_j+1][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);
    
                HessianMatrix.element[index_j+2][n  ]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
                HessianMatrix.element[index_j+2][n+1]-=(posB.z-comB.z)*(f2*dr.z*dr.x);
                HessianMatrix.element[index_j+2][n+2]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
                HessianMatrix.element[index_j+2][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
              {
                if(index_i>=0)
                {
                  HessianMatrix.element[index_i  ][n  ]+=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                  HessianMatrix.element[index_i  ][n+1]+=(posB.y-comB.y)*(f2*dr.x*dr.x+f1);
                  HessianMatrix.element[index_i  ][n+2]+=(posB.y-comB.y)*(f2*dr.x*dr.y);
                  HessianMatrix.element[index_i  ][n+3]+=(posB.z-comB.z)*(f2*dr.x*dr.z);
    
                  HessianMatrix.element[index_i+1][n  ]+=(posB.x-comB.x)*(f2*dr.y*dr.x);
                  HessianMatrix.element[index_i+1][n+1]+=(posB.y-comB.y)*(f2*dr.y*dr.x);
                  HessianMatrix.element[index_i+1][n+2]+=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                  HessianMatrix.element[index_i+1][n+3]+=(posB.z-comB.z)*(f2*dr.y*dr.z);
    
                  HessianMatrix.element[index_i+2][n  ]+=(posB.x-comB.x)*(f2*dr.z*dr.x);
                  HessianMatrix.element[index_i+2][n+1]+=(posB.y-comB.y)*(f2*dr.z*dr.x);
                  HessianMatrix.element[index_i+2][n+2]+=(posB.y-comB.y)*(f2*dr.z*dr.y);
                  HessianMatrix.element[index_i+2][n+3]+=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
                }
              }
    
              if(index_j>=0)
              {
                HessianMatrix.element[index_j  ][n  ]-=(posB.x-comB.x)*(f2*dr.x*dr.x+f1);
                HessianMatrix.element[index_j  ][n+1]-=(posB.y-comB.y)*(f2*dr.x*dr.x+f1);
                HessianMatrix.element[index_j  ][n+2]-=(posB.y-comB.y)*(f2*dr.x*dr.y);
                HessianMatrix.element[index_j  ][n+3]-=(posB.z-comB.z)*(f2*dr.x*dr.z);
    
                HessianMatrix.element[index_j+1][n  ]-=(posB.x-comB.x)*(f2*dr.y*dr.x);
                HessianMatrix.element[index_j+1][n+1]-=(posB.y-comB.y)*(f2*dr.y*dr.x);
                HessianMatrix.element[index_j+1][n+2]-=(posB.y-comB.y)*(f2*dr.y*dr.y+f1);
                HessianMatrix.element[index_j+1][n+3]-=(posB.z-comB.z)*(f2*dr.y*dr.z);
    
                HessianMatrix.element[index_j+2][n  ]-=(posB.x-comB.x)*(f2*dr.z*dr.x);
                HessianMatrix.element[index_j+2][n+1]-=(posB.y-comB.y)*(f2*dr.z*dr.x);
                HessianMatrix.element[index_j+2][n+2]-=(posB.y-comB.y)*(f2*dr.z*dr.y);
                HessianMatrix.element[index_j+2][n+3]-=(posB.z-comB.z)*(f2*dr.z*dr.z+f1);
              }
              break;
          }
          break;  
        default:  
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }


}

// Hessian: Orientation - Strain (part I)
// ======================================
static inline void HessianOrientationStrainJ(REAL_MATRIX HessianMatrix,int index_j2,int index2,REAL f1,REAL f2,VECTOR dr)
{
  int n;
  REAL_MATRIX3x3 Hessian;
  VECTOR vecj1,vecj2,vecj3;

  n=NumberOfCoordinatesMinimizationVariables;

  vecj1=DVecX[index2];
  vecj2=DVecY[index2];
  vecj3=DVecZ[index2];

  Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
  Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
  Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

  Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
  Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
  Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

  Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
  Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
  Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_j2>=0)
      {
        HessianMatrix.element[index_j2][n]-=Hessian.ax*dr.x+Hessian.bx*dr.y+Hessian.cx*dr.z;
        HessianMatrix.element[index_j2+1][n]-=Hessian.ay*dr.x+Hessian.by*dr.y+Hessian.cy*dr.z;
        HessianMatrix.element[index_j2+2][n]-=Hessian.az*dr.x+Hessian.bz*dr.y+Hessian.cz*dr.z;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_j2>=0)
          { 
            HessianMatrix.element[index_j2][n]-=Hessian.ax*dr.x;
            HessianMatrix.element[index_j2][n+1]-=Hessian.bx*dr.y;
            HessianMatrix.element[index_j2][n+2]-=Hessian.cx*dr.z;
            
            HessianMatrix.element[index_j2+1][n]-=Hessian.ay*dr.x;
            HessianMatrix.element[index_j2+1][n+1]-=Hessian.by*dr.y;
            HessianMatrix.element[index_j2+1][n+2]-=Hessian.cy*dr.z;
          
            HessianMatrix.element[index_j2+2][n]-=Hessian.az*dr.x;
            HessianMatrix.element[index_j2+2][n+1]-=Hessian.bz*dr.y;
            HessianMatrix.element[index_j2+2][n+2]-=Hessian.cz*dr.z;
          } 
          break;
        case REGULAR:
          if(index_j2>=0)
          {
            HessianMatrix.element[index_j2][n]-=Hessian.ax*dr.x;
            HessianMatrix.element[index_j2][n+1]-=0.5*(Hessian.ax*dr.y+Hessian.bx*dr.x);
            HessianMatrix.element[index_j2][n+2]-=0.5*(Hessian.ax*dr.z+Hessian.cx*dr.x);
            HessianMatrix.element[index_j2][n+3]-=Hessian.bx*dr.y;
            HessianMatrix.element[index_j2][n+4]-=0.5*(Hessian.bx*dr.z+Hessian.cx*dr.y);
            HessianMatrix.element[index_j2][n+5]-=Hessian.cx*dr.z;

            HessianMatrix.element[index_j2+1][n]-=Hessian.ay*dr.x;
            HessianMatrix.element[index_j2+1][n+1]-=0.5*(Hessian.ay*dr.y+Hessian.by*dr.x);
            HessianMatrix.element[index_j2+1][n+2]-=0.5*(Hessian.ay*dr.z+Hessian.cy*dr.x);
            HessianMatrix.element[index_j2+1][n+3]-=Hessian.by*dr.y;
            HessianMatrix.element[index_j2+1][n+4]-=0.5*(Hessian.by*dr.z+Hessian.cy*dr.y);
            HessianMatrix.element[index_j2+1][n+5]-=Hessian.cy*dr.z;

            HessianMatrix.element[index_j2+2][n]-=Hessian.az*dr.x;
            HessianMatrix.element[index_j2+2][n+1]-=0.5*(Hessian.az*dr.y+Hessian.bz*dr.x);
            HessianMatrix.element[index_j2+2][n+2]-=0.5*(Hessian.az*dr.z+Hessian.cz*dr.x);
            HessianMatrix.element[index_j2+2][n+3]-=Hessian.bz*dr.y;
            HessianMatrix.element[index_j2+2][n+4]-=0.5*(Hessian.bz*dr.z+Hessian.cz*dr.y);
            HessianMatrix.element[index_j2+2][n+5]-=Hessian.cz*dr.z;
          }
          break;                                   
        case REGULAR_UPPER_TRIANGLE:
          if(index_j2>=0)
          {
            HessianMatrix.element[index_j2  ][n  ]-=Hessian.ax*dr.x;
            HessianMatrix.element[index_j2  ][n+1]-=Hessian.ax*dr.y;
            HessianMatrix.element[index_j2  ][n+2]-=Hessian.ax*dr.z;
            HessianMatrix.element[index_j2  ][n+3]-=Hessian.bx*dr.y;
            HessianMatrix.element[index_j2  ][n+4]-=Hessian.bx*dr.z;
            HessianMatrix.element[index_j2  ][n+5]-=Hessian.cx*dr.z;

            HessianMatrix.element[index_j2+1][n  ]-=Hessian.ay*dr.x;
            HessianMatrix.element[index_j2+1][n+1]-=Hessian.ay*dr.y;
            HessianMatrix.element[index_j2+1][n+2]-=Hessian.ay*dr.z;
            HessianMatrix.element[index_j2+1][n+3]-=Hessian.by*dr.y;
            HessianMatrix.element[index_j2+1][n+4]-=Hessian.by*dr.z;
            HessianMatrix.element[index_j2+1][n+5]-=Hessian.cy*dr.z;

            HessianMatrix.element[index_j2+2][n  ]-=Hessian.az*dr.x;
            HessianMatrix.element[index_j2+2][n+1]-=Hessian.az*dr.y;
            HessianMatrix.element[index_j2+2][n+2]-=Hessian.az*dr.z;
            HessianMatrix.element[index_j2+2][n+3]-=Hessian.bz*dr.y;
            HessianMatrix.element[index_j2+2][n+4]-=Hessian.bz*dr.z;
            HessianMatrix.element[index_j2+2][n+5]-=Hessian.cz*dr.z;
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2][n]-=Hessian.ax*dr.x;
                HessianMatrix.element[index_j2][n+1]-=Hessian.bx*dr.y;
                HessianMatrix.element[index_j2][n+2]-=0.5*(Hessian.bx*dr.z+Hessian.cx*dr.y);
                HessianMatrix.element[index_j2][n+3]-=Hessian.cx*dr.z;

                HessianMatrix.element[index_j2+1][n]-=Hessian.ay*dr.x;
                HessianMatrix.element[index_j2+1][n+1]-=Hessian.by*dr.y;
                HessianMatrix.element[index_j2+1][n+2]-=0.5*(Hessian.by*dr.z+Hessian.cy*dr.y);
                HessianMatrix.element[index_j2+1][n+3]-=Hessian.cy*dr.z;

                HessianMatrix.element[index_j2+2][n]-=Hessian.az*dr.x;
                HessianMatrix.element[index_j2+2][n+1]-=Hessian.bz*dr.y;
                HessianMatrix.element[index_j2+2][n+2]-=0.5*(Hessian.bz*dr.z+Hessian.cz*dr.y);
                HessianMatrix.element[index_j2+2][n+3]-=Hessian.cz*dr.z;
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2][n]-=Hessian.ax*dr.x;
                HessianMatrix.element[index_j2][n+1]-=0.5*(Hessian.ax*dr.z+Hessian.cx*dr.x);
                HessianMatrix.element[index_j2][n+2]-=Hessian.bx*dr.y;
                HessianMatrix.element[index_j2][n+3]-=Hessian.cx*dr.z;

                HessianMatrix.element[index_j2+1][n]-=Hessian.ay*dr.x;
                HessianMatrix.element[index_j2+1][n+1]-=0.5*(Hessian.ay*dr.z+Hessian.cy*dr.x);
                HessianMatrix.element[index_j2+1][n+2]-=Hessian.by*dr.y;
                HessianMatrix.element[index_j2+1][n+3]-=Hessian.cy*dr.z;

                HessianMatrix.element[index_j2+2][n]-=Hessian.az*dr.x;
                HessianMatrix.element[index_j2+2][n+1]-=0.5*(Hessian.az*dr.z+Hessian.cz*dr.x);
                HessianMatrix.element[index_j2+2][n+2]-=Hessian.bz*dr.y;
                HessianMatrix.element[index_j2+2][n+3]-=Hessian.cz*dr.z;
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2][n]-=Hessian.ax*dr.x;
                HessianMatrix.element[index_j2][n+1]-=0.5*(Hessian.ax*dr.y+Hessian.bx*dr.x);
                HessianMatrix.element[index_j2][n+2]-=Hessian.bx*dr.y;
                HessianMatrix.element[index_j2][n+3]-=Hessian.cx*dr.z;

                HessianMatrix.element[index_j2+1][n]-=Hessian.ay*dr.x;
                HessianMatrix.element[index_j2+1][n+1]-=0.5*(Hessian.ay*dr.y+Hessian.by*dr.x);
                HessianMatrix.element[index_j2+1][n+2]-=Hessian.by*dr.y;
                HessianMatrix.element[index_j2+1][n+3]-=Hessian.cy*dr.z;

                HessianMatrix.element[index_j2+2][n]-=Hessian.az*dr.x;
                HessianMatrix.element[index_j2+2][n+1]-=0.5*(Hessian.az*dr.y+Hessian.bz*dr.x);
                HessianMatrix.element[index_j2+2][n+2]-=Hessian.bz*dr.y;
                HessianMatrix.element[index_j2+2][n+3]-=Hessian.cz*dr.z;
              }
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2  ][n  ]-=Hessian.ax*dr.x;
                HessianMatrix.element[index_j2  ][n+1]-=Hessian.bx*dr.y;
                HessianMatrix.element[index_j2  ][n+2]-=Hessian.bx*dr.z;
                HessianMatrix.element[index_j2  ][n+3]-=Hessian.cx*dr.z;

                HessianMatrix.element[index_j2+1][n  ]-=Hessian.ay*dr.x;
                HessianMatrix.element[index_j2+1][n+1]-=Hessian.by*dr.y;
                HessianMatrix.element[index_j2+1][n+2]-=Hessian.by*dr.z;
                HessianMatrix.element[index_j2+1][n+3]-=Hessian.cy*dr.z;

                HessianMatrix.element[index_j2+2][n  ]-=Hessian.az*dr.x;
                HessianMatrix.element[index_j2+2][n+1]-=Hessian.bz*dr.y;
                HessianMatrix.element[index_j2+2][n+2]-=Hessian.bz*dr.z;
                HessianMatrix.element[index_j2+2][n+3]-=Hessian.cz*dr.z;
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2  ][n  ]-=Hessian.ax*dr.x;
                HessianMatrix.element[index_j2  ][n+1]-=Hessian.ax*dr.z;
                HessianMatrix.element[index_j2  ][n+2]-=Hessian.bx*dr.y;
                HessianMatrix.element[index_j2  ][n+3]-=Hessian.cx*dr.z;

                HessianMatrix.element[index_j2+1][n  ]-=Hessian.ay*dr.x;
                HessianMatrix.element[index_j2+1][n+1]-=Hessian.ay*dr.z;
                HessianMatrix.element[index_j2+1][n+2]-=Hessian.by*dr.y;
                HessianMatrix.element[index_j2+1][n+3]-=Hessian.cy*dr.z;

                HessianMatrix.element[index_j2+2][n  ]-=Hessian.az*dr.x;
                HessianMatrix.element[index_j2+2][n+1]-=Hessian.az*dr.z;
                HessianMatrix.element[index_j2+2][n+2]-=Hessian.bz*dr.y;
                HessianMatrix.element[index_j2+2][n+3]-=Hessian.cz*dr.z;
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2  ][n  ]-=Hessian.ax*dr.x;
                HessianMatrix.element[index_j2  ][n+1]-=Hessian.ax*dr.y;
                HessianMatrix.element[index_j2  ][n+2]-=Hessian.bx*dr.y;
                HessianMatrix.element[index_j2  ][n+3]-=Hessian.cx*dr.z;

                HessianMatrix.element[index_j2+1][n  ]-=Hessian.ay*dr.x;
                HessianMatrix.element[index_j2+1][n+1]-=Hessian.ay*dr.y;
                HessianMatrix.element[index_j2+1][n+2]-=Hessian.by*dr.y;
                HessianMatrix.element[index_j2+1][n+3]-=Hessian.cy*dr.z;

                HessianMatrix.element[index_j2+2][n  ]-=Hessian.az*dr.x;
                HessianMatrix.element[index_j2+2][n+1]-=Hessian.az*dr.y;
                HessianMatrix.element[index_j2+2][n+2]-=Hessian.bz*dr.y;
                HessianMatrix.element[index_j2+2][n+3]-=Hessian.cz*dr.z;
              }
              break;
          }
          break;  
        default:  
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}

// Hessian: Orientation - Strain (part II)
// =======================================
static inline void HessianOrientationStrainJ_FR(REAL_MATRIX HessianMatrix,int index_j2,int index2,REAL f1,REAL f2,
                                                VECTOR posB,VECTOR comB,VECTOR dr)
{
  int n;
  REAL_MATRIX3x3 Hessian;
  VECTOR vecj1,vecj2,vecj3;

  n=NumberOfCoordinatesMinimizationVariables;

  vecj1=DVecX[index2];
  vecj2=DVecY[index2];
  vecj3=DVecZ[index2];

  Hessian.ax=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.x+f1*vecj1.x;
  Hessian.ay=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.x+f1*vecj2.x;
  Hessian.az=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.x+f1*vecj3.x;

  Hessian.bx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.y+f1*vecj1.y;
  Hessian.by=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.y+f1*vecj2.y;
  Hessian.bz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.y+f1*vecj3.y;

  Hessian.cx=f2*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)*dr.z+f1*vecj1.z;
  Hessian.cy=f2*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)*dr.z+f1*vecj2.z;
  Hessian.cz=f2*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)*dr.z+f1*vecj3.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_j2>=0)
      {
        HessianMatrix.element[index_j2][n]-=(posB.x-comB.x)*Hessian.ax+(posB.y-comB.y)*Hessian.bx+(posB.z-comB.z)*Hessian.cx;
        HessianMatrix.element[index_j2+1][n]-=(posB.x-comB.x)*Hessian.ay+(posB.y-comB.y)*Hessian.by+(posB.z-comB.z)*Hessian.cy;
        HessianMatrix.element[index_j2+2][n]-=(posB.x-comB.x)*Hessian.az+(posB.y-comB.y)*Hessian.bz+(posB.z-comB.z)*Hessian.cz;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_j2>=0)
          {
            HessianMatrix.element[index_j2][n]-=(posB.x-comB.x)*Hessian.ax;
            HessianMatrix.element[index_j2][n+1]-=(posB.y-comB.y)*Hessian.bx;
            HessianMatrix.element[index_j2][n+2]-=(posB.z-comB.z)*Hessian.cx;

            HessianMatrix.element[index_j2+1][n]-=(posB.x-comB.x)*Hessian.ay;
            HessianMatrix.element[index_j2+1][n+1]-=(posB.y-comB.y)*Hessian.by;
            HessianMatrix.element[index_j2+1][n+2]-=(posB.z-comB.z)*Hessian.cy;

            HessianMatrix.element[index_j2+2][n]-=(posB.x-comB.x)*Hessian.az;
            HessianMatrix.element[index_j2+2][n+1]-=(posB.y-comB.y)*Hessian.bz;
            HessianMatrix.element[index_j2+2][n+2]-=(posB.z-comB.z)*Hessian.cz;
          }
          break;
        case REGULAR:
          if(index_j2>=0)
          {
            HessianMatrix.element[index_j2][n]-=(posB.x-comB.x)*Hessian.ax;
            HessianMatrix.element[index_j2][n+1]-=0.5*((posB.y-comB.y)*Hessian.ax+(posB.x-comB.x)*Hessian.bx);
            HessianMatrix.element[index_j2][n+2]-=0.5*((posB.z-comB.z)*Hessian.ax+(posB.x-comB.x)*Hessian.cx);
            HessianMatrix.element[index_j2][n+3]-=(posB.y-comB.y)*Hessian.bx;
            HessianMatrix.element[index_j2][n+4]-=0.5*((posB.z-comB.z)*Hessian.bx+(posB.y-comB.y)*Hessian.cx);
            HessianMatrix.element[index_j2][n+5]-=(posB.z-comB.z)*Hessian.cx;

            HessianMatrix.element[index_j2+1][n]-=(posB.x-comB.x)*Hessian.ay;
            HessianMatrix.element[index_j2+1][n+1]-=0.5*((posB.y-comB.y)*Hessian.ay+(posB.x-comB.x)*Hessian.by);
            HessianMatrix.element[index_j2+1][n+2]-=0.5*((posB.z-comB.z)*Hessian.ay+(posB.x-comB.x)*Hessian.cy);
            HessianMatrix.element[index_j2+1][n+3]-=(posB.y-comB.y)*Hessian.by;
            HessianMatrix.element[index_j2+1][n+4]-=0.5*((posB.z-comB.z)*Hessian.by+(posB.y-comB.y)*Hessian.cy);
            HessianMatrix.element[index_j2+1][n+5]-=(posB.z-comB.z)*Hessian.cy;

            HessianMatrix.element[index_j2+2][n]-=(posB.x-comB.x)*Hessian.az;
            HessianMatrix.element[index_j2+2][n+1]-=0.5*((posB.y-comB.y)*Hessian.az+(posB.x-comB.x)*Hessian.bz);
            HessianMatrix.element[index_j2+2][n+2]-=0.5*((posB.z-comB.z)*Hessian.az+(posB.x-comB.x)*Hessian.cz);
            HessianMatrix.element[index_j2+2][n+3]-=(posB.y-comB.y)*Hessian.bz;
            HessianMatrix.element[index_j2+2][n+4]-=0.5*((posB.z-comB.z)*Hessian.bz+(posB.y-comB.y)*Hessian.cz);
            HessianMatrix.element[index_j2+2][n+5]-=(posB.z-comB.z)*Hessian.cz;
          }
          break;                                   
        case REGULAR_UPPER_TRIANGLE:
          if(index_j2>=0)
          {
            HessianMatrix.element[index_j2  ][n  ]-=(posB.x-comB.x)*Hessian.ax;
            HessianMatrix.element[index_j2  ][n+1]-=(posB.y-comB.y)*Hessian.ax;
            HessianMatrix.element[index_j2  ][n+2]-=(posB.z-comB.z)*Hessian.ax;
            HessianMatrix.element[index_j2  ][n+3]-=(posB.y-comB.y)*Hessian.bx;
            HessianMatrix.element[index_j2  ][n+4]-=(posB.z-comB.z)*Hessian.bx;
            HessianMatrix.element[index_j2  ][n+5]-=(posB.z-comB.z)*Hessian.cx;

            HessianMatrix.element[index_j2+1][n  ]-=(posB.x-comB.x)*Hessian.ay;
            HessianMatrix.element[index_j2+1][n+1]-=(posB.y-comB.y)*Hessian.ay;
            HessianMatrix.element[index_j2+1][n+2]-=(posB.z-comB.z)*Hessian.ay;
            HessianMatrix.element[index_j2+1][n+3]-=(posB.y-comB.y)*Hessian.by;
            HessianMatrix.element[index_j2+1][n+4]-=(posB.z-comB.z)*Hessian.by;
            HessianMatrix.element[index_j2+1][n+5]-=(posB.z-comB.z)*Hessian.cy;

            HessianMatrix.element[index_j2+2][n  ]-=(posB.x-comB.x)*Hessian.az;
            HessianMatrix.element[index_j2+2][n+1]-=(posB.y-comB.y)*Hessian.az;
            HessianMatrix.element[index_j2+2][n+2]-=(posB.z-comB.z)*Hessian.az;
            HessianMatrix.element[index_j2+2][n+3]-=(posB.y-comB.y)*Hessian.bz;
            HessianMatrix.element[index_j2+2][n+4]-=(posB.z-comB.z)*Hessian.bz;
            HessianMatrix.element[index_j2+2][n+5]-=(posB.z-comB.z)*Hessian.cz;
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2][n]-=(posB.x-comB.x)*Hessian.ax;
                HessianMatrix.element[index_j2][n+1]-=(posB.y-comB.y)*Hessian.bx;
                HessianMatrix.element[index_j2][n+2]-=0.5*((posB.z-comB.z)*Hessian.bx+(posB.y-comB.y)*Hessian.cx);
                HessianMatrix.element[index_j2][n+3]-=(posB.z-comB.z)*Hessian.cx;

                HessianMatrix.element[index_j2+1][n]-=(posB.x-comB.x)*Hessian.ay;
                HessianMatrix.element[index_j2+1][n+1]-=(posB.y-comB.y)*Hessian.by;
                HessianMatrix.element[index_j2+1][n+2]-=0.5*((posB.z-comB.z)*Hessian.by+(posB.y-comB.y)*Hessian.cy);
                HessianMatrix.element[index_j2+1][n+3]-=(posB.z-comB.z)*Hessian.cy;

                HessianMatrix.element[index_j2+2][n]-=(posB.x-comB.x)*Hessian.az;
                HessianMatrix.element[index_j2+2][n+1]-=(posB.y-comB.y)*Hessian.bz;
                HessianMatrix.element[index_j2+2][n+2]-=0.5*((posB.z-comB.z)*Hessian.bz+(posB.y-comB.y)*Hessian.cz);
                HessianMatrix.element[index_j2+2][n+3]-=(posB.z-comB.z)*Hessian.cz;
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2][n]-=(posB.x-comB.x)*Hessian.ax;
                HessianMatrix.element[index_j2][n+1]-=0.5*((posB.z-comB.z)*Hessian.ax+(posB.x-comB.x)*Hessian.cx);
                HessianMatrix.element[index_j2][n+2]-=(posB.y-comB.y)*Hessian.bx;
                HessianMatrix.element[index_j2][n+3]-=(posB.z-comB.z)*Hessian.cx;

                HessianMatrix.element[index_j2+1][n]-=(posB.x-comB.x)*Hessian.ay;
                HessianMatrix.element[index_j2+1][n+1]-=0.5*((posB.z-comB.z)*Hessian.ay+(posB.x-comB.x)*Hessian.cy);
                HessianMatrix.element[index_j2+1][n+2]-=(posB.y-comB.y)*Hessian.by;
                HessianMatrix.element[index_j2+1][n+3]-=(posB.z-comB.z)*Hessian.cy;

                HessianMatrix.element[index_j2+2][n]-=(posB.x-comB.x)*Hessian.az;
                HessianMatrix.element[index_j2+2][n+1]-=0.5*((posB.z-comB.z)*Hessian.az+(posB.x-comB.x)*Hessian.cz);
                HessianMatrix.element[index_j2+2][n+2]-=(posB.y-comB.y)*Hessian.bz;
                HessianMatrix.element[index_j2+2][n+3]-=(posB.z-comB.z)*Hessian.cz;
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2][n]-=(posB.x-comB.x)*Hessian.ax;
                HessianMatrix.element[index_j2][n+1]-=0.5*((posB.y-comB.y)*Hessian.ax+(posB.x-comB.x)*Hessian.bx);
                HessianMatrix.element[index_j2][n+2]-=(posB.y-comB.y)*Hessian.bx;
                HessianMatrix.element[index_j2][n+3]-=(posB.z-comB.z)*Hessian.cx;

                HessianMatrix.element[index_j2+1][n]-=(posB.x-comB.x)*Hessian.ay;
                HessianMatrix.element[index_j2+1][n+1]-=0.5*((posB.y-comB.y)*Hessian.ay+(posB.x-comB.x)*Hessian.by);
                HessianMatrix.element[index_j2+1][n+2]-=(posB.y-comB.y)*Hessian.by;
                HessianMatrix.element[index_j2+1][n+3]-=(posB.z-comB.z)*Hessian.cy;

                HessianMatrix.element[index_j2+2][n]-=(posB.x-comB.x)*Hessian.az;
                HessianMatrix.element[index_j2+2][n+1]-=0.5*((posB.y-comB.y)*Hessian.az+(posB.x-comB.x)*Hessian.bz);
                HessianMatrix.element[index_j2+2][n+2]-=(posB.y-comB.y)*Hessian.bz;
                HessianMatrix.element[index_j2+2][n+3]-=(posB.z-comB.z)*Hessian.cz;
              }
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2  ][n  ]-=(posB.x-comB.x)*Hessian.ax;
                HessianMatrix.element[index_j2  ][n+1]-=(posB.y-comB.y)*Hessian.bx;
                HessianMatrix.element[index_j2  ][n+2]-=(posB.z-comB.z)*Hessian.bx;
                HessianMatrix.element[index_j2  ][n+3]-=(posB.z-comB.z)*Hessian.cx;

                HessianMatrix.element[index_j2+1][n  ]-=(posB.x-comB.x)*Hessian.ay;
                HessianMatrix.element[index_j2+1][n+1]-=(posB.y-comB.y)*Hessian.by;
                HessianMatrix.element[index_j2+1][n+2]-=(posB.z-comB.z)*Hessian.by;
                HessianMatrix.element[index_j2+1][n+3]-=(posB.z-comB.z)*Hessian.cy;

                HessianMatrix.element[index_j2+2][n  ]-=(posB.x-comB.x)*Hessian.az;
                HessianMatrix.element[index_j2+2][n+1]-=(posB.y-comB.y)*Hessian.bz;
                HessianMatrix.element[index_j2+2][n+2]-=(posB.z-comB.z)*Hessian.bz;
                HessianMatrix.element[index_j2+2][n+3]-=(posB.z-comB.z)*Hessian.cz;
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2  ][n  ]-=(posB.x-comB.x)*Hessian.ax;
                HessianMatrix.element[index_j2  ][n+1]-=(posB.z-comB.z)*Hessian.ax;
                HessianMatrix.element[index_j2  ][n+2]-=(posB.y-comB.y)*Hessian.bx;
                HessianMatrix.element[index_j2  ][n+3]-=(posB.z-comB.z)*Hessian.cx;

                HessianMatrix.element[index_j2+1][n  ]-=(posB.x-comB.x)*Hessian.ay;
                HessianMatrix.element[index_j2+1][n+1]-=(posB.z-comB.z)*Hessian.ay;
                HessianMatrix.element[index_j2+1][n+2]-=(posB.y-comB.y)*Hessian.by;
                HessianMatrix.element[index_j2+1][n+3]-=(posB.z-comB.z)*Hessian.cy;

                HessianMatrix.element[index_j2+2][n  ]-=(posB.x-comB.x)*Hessian.az;
                HessianMatrix.element[index_j2+2][n+1]-=(posB.z-comB.z)*Hessian.az;
                HessianMatrix.element[index_j2+2][n+2]-=(posB.y-comB.y)*Hessian.bz;
                HessianMatrix.element[index_j2+2][n+3]-=(posB.z-comB.z)*Hessian.cz;
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_j2>=0)
              {
                HessianMatrix.element[index_j2  ][n  ]-=(posB.x-comB.x)*Hessian.ax;
                HessianMatrix.element[index_j2  ][n+1]-=(posB.y-comB.y)*Hessian.ax;
                HessianMatrix.element[index_j2  ][n+2]-=(posB.y-comB.y)*Hessian.bx;
                HessianMatrix.element[index_j2  ][n+3]-=(posB.z-comB.z)*Hessian.cx;

                HessianMatrix.element[index_j2+1][n  ]-=(posB.x-comB.x)*Hessian.ay;
                HessianMatrix.element[index_j2+1][n+1]-=(posB.y-comB.y)*Hessian.ay;
                HessianMatrix.element[index_j2+1][n+2]-=(posB.y-comB.y)*Hessian.by;
                HessianMatrix.element[index_j2+1][n+3]-=(posB.z-comB.z)*Hessian.cy;

                HessianMatrix.element[index_j2+2][n  ]-=(posB.x-comB.x)*Hessian.az;
                HessianMatrix.element[index_j2+2][n+1]-=(posB.y-comB.y)*Hessian.az;
                HessianMatrix.element[index_j2+2][n+2]-=(posB.y-comB.y)*Hessian.bz;
                HessianMatrix.element[index_j2+2][n+3]-=(posB.z-comB.z)*Hessian.cz;
              }
              break;
          }
          break;  
        default:  
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}


// Hessian: Strain - Strain 
// ========================
static inline void HessianAtomicStrainStrain(REAL_MATRIX HessianMatrix,int index_i,int index_j,REAL f1,REAL f2,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x+2.0*f2*dr.x*dr.x*dr.y*dr.y+2.0*f2*dr.x*dr.x*dr.z*dr.z+
                            f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y+2.0*f2*dr.y*dr.y*dr.z*dr.z+f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x; // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                  // xxyy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.z*dr.z;                  // xxzz
          HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
          HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz
          HessianMatrix.element[n+2][n+2]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
          break;
        case REGULAR:
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;              // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;                  // xxxy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;                  // xxxz
          HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.y*dr.y;                               // xxyy
          HessianMatrix.element[n  ][n+4]+=f2*dr.x*dr.x*dr.y*dr.z;                               // xxyz
          HessianMatrix.element[n  ][n+5]+=f2*dr.x*dr.x*dr.z*dr.z;                               // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+0.5*f1*(dr.x*dr.x+dr.y*dr.y);  // xyxy
          HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.x*dr.z+0.5*f1*dr.y*dr.z;              // xyxz
          HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.y*dr.y+f1*dr.x*dr.y;                  // xyyy
          HessianMatrix.element[n+1][n+4]+=f2*dr.x*dr.y*dr.y*dr.z+0.5*f1*(dr.x*dr.z);            // xyyz
          HessianMatrix.element[n+1][n+5]+=f2*dr.x*dr.y*dr.z*dr.z;                               // xyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.x*dr.z*dr.x*dr.z+0.5*f1*(dr.x*dr.x+dr.z*dr.z);  // xzxz
          HessianMatrix.element[n+2][n+3]+=f2*dr.x*dr.z*dr.y*dr.y;                               // xzyy
          HessianMatrix.element[n+2][n+4]+=f2*dr.x*dr.z*dr.y*dr.z+0.5*f1*dr.x*dr.y;              // xzyz
          HessianMatrix.element[n+2][n+5]+=f2*dr.x*dr.z*dr.z*dr.z+f1*dr.x*dr.z;                  // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;                  // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

          HessianMatrix.element[n+4][n+4]+=f2*dr.y*dr.z*dr.y*dr.z+0.5*f1*(dr.y*dr.y+dr.z*dr.z);  // yzyz
          HessianMatrix.element[n+4][n+5]+=f2*dr.y*dr.z*dr.z*dr.z+f1*dr.y*dr.z;                  // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
          break;                                   
        case REGULAR_UPPER_TRIANGLE:
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x; // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;     // xxxy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;     // xxxz
          HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.y*dr.y;                  // xxyy
          HessianMatrix.element[n  ][n+4]+=f2*dr.x*dr.x*dr.y*dr.z;                  // xxyz
          HessianMatrix.element[n  ][n+5]+=f2*dr.x*dr.x*dr.z*dr.z;                  // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+f1*dr.y*dr.y;     // xyxy
          HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.x*dr.z+f1*dr.y*dr.z;     // xyxz
          HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.y*dr.y;                  // xyyy
          HessianMatrix.element[n+1][n+4]+=f2*dr.x*dr.y*dr.y*dr.z;                  // xyyz
          HessianMatrix.element[n+1][n+5]+=f2*dr.x*dr.y*dr.z*dr.z;                  // xyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.x*dr.z*dr.x*dr.z+f1*dr.z*dr.z;     // xzxz
          HessianMatrix.element[n+2][n+3]+=f2*dr.x*dr.z*dr.y*dr.y;                  // xzyy
          HessianMatrix.element[n+2][n+4]+=f2*dr.x*dr.z*dr.y*dr.z;                  // xzyz
          HessianMatrix.element[n+2][n+5]+=f2*dr.x*dr.z*dr.z*dr.z;                  // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;     // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

          HessianMatrix.element[n+4][n+4]+=f2*dr.y*dr.z*dr.y*dr.z+f1*dr.z*dr.z;     // yzyz
          HessianMatrix.element[n+4][n+5]+=f2*dr.y*dr.z*dr.z*dr.z;                  // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
          break;                          
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;              // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                               // xxyy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.z;                               // xxyz
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;                  // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.z*dr.y*dr.z+0.5*f1*(dr.y*dr.y+dr.z*dr.z);  // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.z*dr.z*dr.z+f1*dr.y*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;              // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;                  // xxxz
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                               // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.z*dr.x*dr.z+0.5*f1*(dr.x*dr.x+dr.z*dr.z);  // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.z*dr.y*dr.y;                               // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.z*dr.z*dr.z+f1*dr.x*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;              // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;                  // xxxy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                               // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+0.5*f1*(dr.x*dr.x+dr.y*dr.y);  // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.y*dr.y+f1*dr.x*dr.y;                  // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.z*dr.z;                               // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
          }
          break;  
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.z;                  // xxyz
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;     // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.z*dr.y*dr.z+f1*dr.z*dr.z;     // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.z*dr.z*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;     // xxxz
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.z*dr.x*dr.z+f1*dr.z*dr.z;     // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.z*dr.y*dr.y;                  // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.z*dr.z*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;     // xxxy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+f1*dr.y*dr.y;     // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.y*dr.y;                  // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.z*dr.z;                  // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
          }
          break;  
        default:  
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}

// Hessian: Strain - Strain (part II)
// ==================================
static inline void HessianCorrectionStrainStrain(REAL_MATRIX HessianMatrix,REAL f1,REAL f2,VECTOR dr,VECTOR posB,VECTOR comB)
{
  int n;
  VECTOR dB,drJ;

  n=NumberOfCoordinatesMinimizationVariables;
  dB.x=posB.x-comB.x;
  dB.y=posB.y-comB.y;
  dB.z=posB.z-comB.z;

  drJ.x=dr.x+dB.x;
  drJ.y=dr.y+dB.y;
  drJ.z=dr.z+dB.z;


  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // atomic correction
      HessianMatrix.element[n][n]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxy
      HessianMatrix.element[n][n]+=2.0*f2*dr.x*dB.x*dr.y*dr.y;              // xxyy
      HessianMatrix.element[n][n]+=2.0*f2*dr.x*dB.x*dr.z*dr.z;              // xxzz

      HessianMatrix.element[n][n]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
      HessianMatrix.element[n][n]+=2.0*f2*dr.y*dB.y*dr.z*dr.z;              // yyzz

      HessianMatrix.element[n][n]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

      // rigid correction
      HessianMatrix.element[n][n]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x;   // xxxy
      HessianMatrix.element[n][n]+=2.0*f2*dr.x*drJ.x*dr.y*dB.y;             // xxyy
      HessianMatrix.element[n][n]+=2.0*f2*dr.x*drJ.x*dr.z*dB.z;             // xxzz

      HessianMatrix.element[n][n]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y;   // yyyy
      HessianMatrix.element[n][n]+=2.0*f2*dr.y*drJ.y*dr.z*dB.z;             // yyzz

      HessianMatrix.element[n][n]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z;   // zzzz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          // atomic correction
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxy
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dB.x*dr.y*dr.y;                  // xxyy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dB.x*dr.z*dr.z;                  // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
          HessianMatrix.element[n+1][n+2]+=f2*dr.y*dB.y*dr.z*dr.z;                  // yyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

          // rigid correction
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x;   // xxxy
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*drJ.x*dr.y*dB.y;                 // xxyy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*drJ.x*dr.z*dB.z;                 // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y;   // yyyy
          HessianMatrix.element[n+1][n+2]+=f2*dr.y*drJ.y*dr.z*dB.z;                 // yyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z;   // zzzz
          break;
        case REGULAR:
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dB.x*dr.x+2.0*f1*dB.x*dr.x; // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dB.x*dr.y+f1*dB.x*dr.y;     // xxxy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.x*dB.x*dr.z+f1*dB.x*dr.z;     // xxxz
          HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.y*dB.x*dr.y;                  // xxyy
          HessianMatrix.element[n  ][n+4]+=f2*dr.x*dr.y*dB.x*dr.z;                  // xxyz
          HessianMatrix.element[n  ][n+5]+=f2*dr.x*dr.z*dB.x*dr.z;                  // xxzz

          HessianMatrix.element[n+1][n+1]+=0.5*(f2*(dr.x*dB.y*dr.x*dr.y+dr.y*dB.x*dr.x*dr.y)+f1*(dB.x*dr.x+dB.y*dr.y));   // xyxy
          HessianMatrix.element[n+1][n+2]+=0.5*(f2*(dr.x*dB.y*dr.x*dr.z+dr.y*dB.x*dr.x*dr.z)+f1*(dB.y*dr.z));             // xyxz
          HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.x*dB.y*dr.y*dr.y+dr.y*dB.x*dr.y*dr.y)+f1*(dB.x*dr.y);               // xyyy
          HessianMatrix.element[n+1][n+4]+=0.5*f2*(dr.x*dB.y*dr.y*dr.z+dr.y*dB.x*dr.y*dr.z)+0.5*f1*(dB.x*dr.z);           // xyyz
          HessianMatrix.element[n+1][n+5]+=0.5*f2*(dr.x*dB.y*dr.z*dr.z+dr.y*dB.x*dr.z*dr.z);                              // xyzz

          HessianMatrix.element[n+2][n+2]+=0.5*(f2*(dr.x*dB.z*dr.x*dr.z+dr.z*dB.x*dr.x*dr.z)+f1*(dB.x*dr.x+dB.z*dr.z));   // xzxz
          HessianMatrix.element[n+2][n+3]+=0.5*f2*(dr.x*dB.z*dr.y*dr.y+dB.x*dr.z*dr.y*dr.y);                              // xzyy
          HessianMatrix.element[n+2][n+4]+=0.5*(f2*(dr.x*dB.z*dr.y*dr.z+dB.x*dr.z*dr.y*dr.z)+f1*dB.x*dr.y);               // xzyz
          HessianMatrix.element[n+2][n+5]+=0.5*f2*(dr.x*dB.z*dr.z*dr.z+dB.x*dr.z*dr.z*dr.z)+f1*dB.x*dr.z;                 // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dr.y*dB.y*dr.y+2.0*f1*dB.y*dr.y;   // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dB.y*dr.y*dr.z+f1*dB.y*dr.z;       // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dB.y*dr.z*dr.z;                    // yyzz

          HessianMatrix.element[n+4][n+4]+=0.5*(f2*(dr.y*dB.z*dr.y*dr.z+dB.y*dr.z*dr.y*dr.z)+f1*(dB.y*dr.y+dB.z*dr.z));   // yzyz
          HessianMatrix.element[n+4][n+5]+=0.5*f2*(dr.y*dB.z*dr.z*dr.z+dB.y*dr.z*dr.z*dr.z)+f1*dB.y*dr.z;                 // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z;     

          HessianMatrix.element[n][n]+=f2*dr.x*dB.x*drJ.x*dr.x+f1*dB.x*drJ.x;
          HessianMatrix.element[n][n+1]+=0.5*f2*(dr.x*dB.y+dr.y*dB.x)*drJ.x*dr.x+0.5*f1*dB.y*(drJ.x);
          HessianMatrix.element[n][n+2]+=0.5*f2*(dr.x*dB.z+dr.z*dB.x)*drJ.x*dr.x+0.5*f1*dB.z*(drJ.x);
          HessianMatrix.element[n][n+3]+=0.5*(dr.y*dB.y+dr.y*dB.y)*f2*drJ.x*dr.x;
          HessianMatrix.element[n][n+4]+=0.5*(dr.y*dB.z+dr.z*dB.y)*f2*drJ.x*dr.x;
          HessianMatrix.element[n][n+5]+=0.5*(dr.z*dB.z+dr.z*dB.z)*f2*drJ.x*dr.x;

          HessianMatrix.element[n+1][n+1]+=0.25*f2*(dr.x*dB.y+dr.y*dB.x)*(drJ.x*dr.y+drJ.y*dr.x)+
                                           0.25*f1*(dB.y*drJ.y+dB.x*drJ.x);
          HessianMatrix.element[n+1][n+2]+=0.25*f2*(dr.x*dB.z+dr.z*dB.x)*(drJ.x*dr.y+drJ.y*dr.x)+0.25*f1*dB.z*drJ.y;
          HessianMatrix.element[n+1][n+3]+=0.25*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.x*dr.y+drJ.y*dr.x)
                                           +0.25*f1*(dB.y*drJ.x+dB.y*drJ.x);
          HessianMatrix.element[n+1][n+4]+=0.25*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.x*dr.y+drJ.y*dr.x)+0.25*f1*dB.z*drJ.x;
          HessianMatrix.element[n+1][n+5]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.x*dr.y+drJ.y*dr.x);

          HessianMatrix.element[n+2][n+2]+=0.25*f2*(dr.x*dB.z+dr.z*dB.x)*(drJ.x*dr.z+drJ.z*dr.x)+
                                           0.25*f1*(dB.z*drJ.z+dB.x*drJ.x);
          HessianMatrix.element[n+2][n+3]+=0.25*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.x*dr.z+drJ.z*dr.x);
          HessianMatrix.element[n+2][n+4]+=0.25*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.x*dr.z+drJ.z*dr.x)+0.25*f1*dB.y*drJ.x;
          HessianMatrix.element[n+2][n+5]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.x*dr.z+drJ.z*dr.x)+
                                           0.25*f1*(dB.z*drJ.x+dB.z*drJ.x);

          HessianMatrix.element[n+3][n+3]+=0.5*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.y*dr.y)+0.5*f1*(dB.y*drJ.y+dB.y*drJ.y);
          HessianMatrix.element[n+3][n+4]+=0.5*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.y*dr.y)+0.5*f1*dB.z*drJ.y;
          HessianMatrix.element[n+3][n+5]+=0.5*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.y);

          HessianMatrix.element[n+4][n+4]+=0.25*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.y*dr.z+drJ.z*dr.y)+
                                           0.25*f1*(dB.z*drJ.z+dB.y*drJ.y);
          HessianMatrix.element[n+4][n+5]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.z+drJ.z*dr.y)+
                                           0.25*f1*(dB.z*drJ.y+dB.z*drJ.y);

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dB.z*drJ.z*dr.z+f1*dB.z*drJ.z;
          break;                                   
        case REGULAR_UPPER_TRIANGLE:
          // atomic correction
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*dB.x*dr.x*dr.y+f1*dB.x*dr.y;     // xxxy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*dB.x*dr.x*dr.z+f1*dB.x*dr.z;     // xxxz
          HessianMatrix.element[n  ][n+3]+=f2*dr.x*dB.x*dr.y*dr.y;                  // xxyy
          HessianMatrix.element[n  ][n+4]+=f2*dr.x*dB.x*dr.y*dr.z;                  // xxyz
          HessianMatrix.element[n  ][n+5]+=f2*dr.x*dB.x*dr.z*dr.z;                  // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.x*dB.y*dr.x*dr.y+f1*dB.y*dr.y;     // xyxy
          HessianMatrix.element[n+1][n+2]+=f2*dr.x*dB.y*dr.x*dr.z+f1*dB.y*dr.z;     // xyxz
          HessianMatrix.element[n+1][n+3]+=f2*dr.x*dB.y*dr.y*dr.y;                  // xyyy
          HessianMatrix.element[n+1][n+4]+=f2*dr.x*dB.y*dr.y*dr.z;                  // xyyz
          HessianMatrix.element[n+1][n+5]+=f2*dr.x*dB.y*dr.z*dr.z;                  // xyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.x*dB.z*dr.x*dr.z+f1*dB.z*dr.z;     // xzxz
          HessianMatrix.element[n+2][n+3]+=f2*dr.x*dB.z*dr.y*dr.y;                  // xzyy
          HessianMatrix.element[n+2][n+4]+=f2*dr.x*dB.z*dr.y*dr.z;                  // xzyz
          HessianMatrix.element[n+2][n+5]+=f2*dr.x*dB.z*dr.z*dr.z;                  // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dB.y*dr.y*dr.z+f1*dB.y*dr.z;     // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dB.y*dr.z*dr.z;                  // yyzz

          HessianMatrix.element[n+4][n+4]+=f2*dr.y*dB.z*dr.y*dr.z+f1*dB.z*dr.z;     // yzyz
          HessianMatrix.element[n+4][n+5]+=f2*dr.y*dB.z*dr.z*dr.z;                  // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

          // rigid correction
          HessianMatrix.element[n  ][n  ]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x; // xxxx
          HessianMatrix.element[n  ][n+1]+=f2*dr.x*drJ.x*dr.x*dB.y+f1*drJ.x*dB.y; // xxxy
          HessianMatrix.element[n  ][n+2]+=f2*dr.x*drJ.x*dr.x*dB.z+f1*drJ.x*dB.z; // xxxz
          HessianMatrix.element[n  ][n+3]+=f2*dr.x*drJ.x*dr.y*dB.y;               // xxyy
          HessianMatrix.element[n  ][n+4]+=f2*dr.x*drJ.x*dr.y*dB.z;               // xxyz
          HessianMatrix.element[n  ][n+5]+=f2*dr.x*drJ.x*dr.z*dB.z;               // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.x*drJ.y*dr.x*dB.y+f1*drJ.y*dB.y; // xyxy
          HessianMatrix.element[n+1][n+2]+=f2*dr.x*drJ.y*dr.x*dB.z+f1*drJ.y*dB.z; // xyxz
          HessianMatrix.element[n+1][n+3]+=f2*dr.x*drJ.y*dr.y*dB.y;               // xyyy
          HessianMatrix.element[n+1][n+4]+=f2*dr.x*drJ.y*dr.y*dB.z;               // xyyz
          HessianMatrix.element[n+1][n+5]+=f2*dr.x*drJ.y*dr.z*dB.z;               // xyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.x*drJ.z*dr.x*dB.z+f1*drJ.z*dB.z; // xzxz
          HessianMatrix.element[n+2][n+3]+=f2*dr.x*drJ.z*dr.y*dB.y;               // xzyy
          HessianMatrix.element[n+2][n+4]+=f2*dr.x*drJ.z*dr.y*dB.z;               // xzyz
          HessianMatrix.element[n+2][n+5]+=f2*dr.x*drJ.z*dr.z*dB.z;               // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y; // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*drJ.y*dr.y*dB.z+f1*drJ.y*dB.z; // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*drJ.y*dr.z*dB.z;               // yyzz

          HessianMatrix.element[n+4][n+4]+=f2*dr.y*drJ.z*dr.y*dB.z+f1*drJ.z*dB.z; // yzyz
          HessianMatrix.element[n+4][n+5]+=f2*dr.y*drJ.z*dr.z*dB.z;               // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z; // zzzz
          break;                                   
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dB.x*dr.x+2.0*f1*dB.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.y*dB.x*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.y*dB.x*dr.z;                  // xxyz
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.z*dB.x*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dB.y*dr.y+2.0*f1*dB.y*dr.y;   // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dB.y*dr.y*dr.z+f1*dB.y*dr.z;       // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+2][n+2]+=0.5*(f2*(dr.y*dB.z*dr.y*dr.z+dB.y*dr.z*dr.y*dr.z)+f1*(dB.y*dr.y+dB.z*dr.z));   // yzyz
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(dr.y*dB.z*dr.z*dr.z+dB.y*dr.z*dr.z*dr.z)+f1*dB.y*dr.z;                 // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z;     

              HessianMatrix.element[n][n]+=f2*dr.x*dB.x*drJ.x*dr.x+f1*dB.x*drJ.x;
              HessianMatrix.element[n][n+1]+=0.5*(dr.y*dB.y+dr.y*dB.y)*f2*drJ.x*dr.x;
              HessianMatrix.element[n][n+2]+=0.5*(dr.y*dB.z+dr.z*dB.y)*f2*drJ.x*dr.x;
              HessianMatrix.element[n][n+3]+=0.5*(dr.z*dB.z+dr.z*dB.z)*f2*drJ.x*dr.x;

              HessianMatrix.element[n+1][n+1]+=0.5*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.y*dr.y)+0.5*f1*(dB.y*drJ.y+dB.y*drJ.y);
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.y*dr.y)+0.5*f1*dB.z*drJ.y;
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.y);

              HessianMatrix.element[n+2][n+2]+=0.25*f2*(dr.y*dB.z+dr.z*dB.y)*(drJ.y*dr.z+drJ.z*dr.y)+
                                               0.25*f1*(dB.z*drJ.z+dB.y*drJ.y);
              HessianMatrix.element[n+2][n+3]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.z+drJ.z*dr.y)+
                                               0.25*f1*(dB.z*drJ.y+dB.z*drJ.y);

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*drJ.z*dr.z+f1*dB.z*drJ.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dB.x*dr.x+2.0*f1*dB.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dB.x*dr.z+f1*dB.x*dr.z;     // xxxz
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.y*dB.x*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.z*dB.x*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=0.5*(f2*(dr.x*dB.z*dr.x*dr.z+dr.z*dB.x*dr.x*dr.z)+f1*(dB.x*dr.x+dB.z*dr.z));   // xzxz
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(dr.x*dB.z*dr.y*dr.y+dB.x*dr.z*dr.y*dr.y);                              // xzyy
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.x*dB.z*dr.z*dr.z+dB.x*dr.z*dr.z*dr.z)+f1*dB.x*dr.z;                 // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dB.y*dr.y+2.0*f1*dB.y*dr.y;   // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z;     

              HessianMatrix.element[n][n]+=f2*dr.x*dB.x*drJ.x*dr.x+f1*dB.x*drJ.x;
              HessianMatrix.element[n][n+1]+=0.5*f2*(dr.x*dB.z+dr.z*dB.x)*drJ.x*dr.x+0.5*f1*dB.z*(drJ.x);
              HessianMatrix.element[n][n+2]+=0.5*(dr.y*dB.y+dr.y*dB.y)*f2*drJ.x*dr.x;
              HessianMatrix.element[n][n+3]+=0.5*(dr.z*dB.z+dr.z*dB.z)*f2*drJ.x*dr.x;

              HessianMatrix.element[n+1][n+1]+=0.25*f2*(dr.x*dB.z+dr.z*dB.x)*(drJ.x*dr.z+drJ.z*dr.x)+
                                               0.25*f1*(dB.z*drJ.z+dB.x*drJ.x);
              HessianMatrix.element[n+1][n+2]+=0.25*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.x*dr.z+drJ.z*dr.x);
              HessianMatrix.element[n+1][n+3]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.x*dr.z+drJ.z*dr.x)+
                                               0.25*f1*(dB.z*drJ.x+dB.z*drJ.x);

              HessianMatrix.element[n+2][n+2]+=0.5*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.y*dr.y)+0.5*f1*(dB.y*drJ.y+dB.y*drJ.y);
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.y);

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*drJ.z*dr.z+f1*dB.z*drJ.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dr.x*dB.x*dr.x+2.0*f1*dB.x*dr.x; // xxxx
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dr.x*dB.x*dr.y+f1*dB.x*dr.y;     // xxxy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dr.y*dB.x*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dr.z*dB.x*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=0.5*(f2*(dr.x*dB.y*dr.x*dr.y+dr.y*dB.x*dr.x*dr.y)+f1*(dB.x*dr.x+dB.y*dr.y));   // xyxy
              HessianMatrix.element[n+1][n+2]+=0.5*f2*(dr.x*dB.y*dr.y*dr.y+dr.y*dB.x*dr.y*dr.y)+f1*(dB.x*dr.y);               // xyyy
              HessianMatrix.element[n+1][n+3]+=0.5*f2*(dr.x*dB.y*dr.z*dr.z+dr.y*dB.x*dr.z*dr.z);                              // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dB.y*dr.y+2.0*f1*dB.y*dr.y;   // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                    // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z;     

              HessianMatrix.element[n][n]+=f2*dr.x*dB.x*drJ.x*dr.x+f1*dB.x*drJ.x;
              HessianMatrix.element[n][n+1]+=0.5*f2*(dr.x*dB.y+dr.y*dB.x)*drJ.x*dr.x+0.5*f1*dB.y*(drJ.x);
              HessianMatrix.element[n][n+2]+=0.5*(dr.y*dB.y+dr.y*dB.y)*f2*drJ.x*dr.x;
              HessianMatrix.element[n][n+3]+=0.5*(dr.z*dB.z+dr.z*dB.z)*f2*drJ.x*dr.x;

              HessianMatrix.element[n+1][n+1]+=0.25*f2*(dr.x*dB.y+dr.y*dB.x)*(drJ.x*dr.y+drJ.y*dr.x)+
                                               0.25*f1*(dB.y*drJ.y+dB.x*drJ.x);
              HessianMatrix.element[n+1][n+2]+=0.25*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.x*dr.y+drJ.y*dr.x)
                                               +0.25*f1*(dB.y*drJ.x+dB.y*drJ.x);
              HessianMatrix.element[n+1][n+3]+=0.25*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.x*dr.y+drJ.y*dr.x);

              HessianMatrix.element[n+2][n+2]+=0.5*f2*(dr.y*dB.y+dr.y*dB.y)*(drJ.y*dr.y)+0.5*f1*(dB.y*drJ.y+dB.y*drJ.y);
              HessianMatrix.element[n+2][n+3]+=0.5*f2*(dr.z*dB.z+dr.z*dB.z)*(drJ.y*dr.y);

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*drJ.z*dr.z+f1*dB.z*drJ.z;
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // atomic correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dB.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dB.x*dr.y*dr.z;                  // xxyz
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dB.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dB.y*dr.y*dr.z+f1*dB.y*dr.z;     // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dB.z*dr.y*dr.z+f1*dB.z*dr.z;     // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dB.z*dr.z*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

              // rigid correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*drJ.x*dr.y*dB.y;               // xxyy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*drJ.x*dr.y*dB.z;               // xxyz
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*drJ.x*dr.z*dB.z;               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y; // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*drJ.y*dr.y*dB.z+f1*drJ.y*dB.z; // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*drJ.y*dr.z*dB.z;               // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*drJ.z*dr.y*dB.z+f1*drJ.z*dB.z; // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*drJ.z*dr.z*dB.z;               // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z; // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // atomic correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dB.x*dr.x*dr.z+f1*dB.x*dr.z;     // xxxz
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dB.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dB.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dB.z*dr.x*dr.z+f1*dB.z*dr.z;     // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dB.z*dr.y*dr.y;                  // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dB.z*dr.z*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

              // rigid correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*drJ.x*dr.x*dB.z+f1*drJ.x*dB.z; // xxxz
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*drJ.x*dr.y*dB.y;               // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*drJ.x*dr.z*dB.z;               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*drJ.z*dr.x*dB.z+f1*drJ.z*dB.z; // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*drJ.z*dr.y*dB.y;               // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*drJ.z*dr.z*dB.z;               // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*drJ.y*dr.z*dB.z;               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z; // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // atomic correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*dB.x*dr.x*dr.x+2.0*f1*dB.x*dr.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*dB.x*dr.x*dr.y+f1*dB.x*dr.y;     // xxxy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*dB.x*dr.y*dr.y;                  // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*dB.x*dr.z*dr.z;                  // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dB.y*dr.x*dr.y+f1*dB.y*dr.y;     // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dB.y*dr.y*dr.y;                  // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dB.y*dr.z*dr.z;                  // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dB.y*dr.y*dr.y+2.0*f1*dB.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dB.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dB.z*dr.z*dr.z+2.0*f1*dB.z*dr.z; // zzzz

              // rigid correction
              HessianMatrix.element[n  ][n  ]+=f2*dr.x*drJ.x*dr.x*dB.x+f1*drJ.x*dB.x; // xxxy
              HessianMatrix.element[n  ][n+1]+=f2*dr.x*drJ.x*dr.x*dB.y+f1*drJ.x*dB.y; // xxxy
              HessianMatrix.element[n  ][n+2]+=f2*dr.x*drJ.x*dr.y*dB.y;               // xxyy
              HessianMatrix.element[n  ][n+3]+=f2*dr.x*drJ.x*dr.z*dB.z;               // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*drJ.y*dr.x*dB.y+f1*drJ.y*dB.y; // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*drJ.y*dr.y*dB.y;               // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*drJ.y*dr.z*dB.z;               // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*drJ.y*dr.y*dB.y+f1*drJ.y*dB.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*drJ.y*dr.z*dB.z;               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*drJ.z*dr.z*dB.z+f1*drJ.z*dB.z; // zzzz
              break;
          }
          break;  
        default:  
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

static inline void GradientStrain(REAL *Gradient,REAL f1,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]+=f1*(dr.x*dr.x+dr.y*dr.y+dr.z*dr.z); // xx + yy + zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n  ]+=f1*dr.x*dr.x; // xx
          Gradient[n+1]+=f1*dr.y*dr.y; // yy
          Gradient[n+2]+=f1*dr.z*dr.z; // zz
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n  ]+=f1*dr.x*dr.x;  // xx
          Gradient[n+1]+=f1*dr.x*dr.y;  // xy
          Gradient[n+2]+=f1*dr.x*dr.z;  // xz
          Gradient[n+3]+=f1*dr.y*dr.y;  // yy
          Gradient[n+4]+=f1*dr.y*dr.z;  // yz
          Gradient[n+5]+=f1*dr.z*dr.z;  // zz
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]+=f1*dr.x*dr.x;  // xx
              Gradient[n+1]+=f1*dr.y*dr.y;  // yy
              Gradient[n+2]+=f1*dr.y*dr.z;  // yz
              Gradient[n+3]+=f1*dr.z*dr.z;  // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]+=f1*dr.x*dr.x;  // xx
              Gradient[n+1]+=f1*dr.x*dr.z;  // xz
              Gradient[n+2]+=f1*dr.y*dr.y;  // yy
              Gradient[n+3]+=f1*dr.z*dr.z;  // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]+=f1*dr.x*dr.x;  // xx
              Gradient[n+1]+=f1*dr.x*dr.y;  // xy
              Gradient[n+2]+=f1*dr.y*dr.y;  // yy
              Gradient[n+3]+=f1*dr.z*dr.z;  // zz
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}

static inline void GradientStrainJ(REAL *Gradient,REAL f1,VECTOR dr,VECTOR posB,VECTOR comB)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n  ]+=f1*(dr.x*(posB.x-comB.x)+dr.y*(posB.y-comB.y)+dr.z*(posB.z-comB.z)); // xx + yy + zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n  ]+=f1*dr.x*(posB.x-comB.x); // xx
          Gradient[n+1]+=f1*dr.y*(posB.y-comB.y); // yy
          Gradient[n+2]+=f1*dr.z*(posB.z-comB.z); // zz
          break;
        case REGULAR:
          Gradient[n  ]+=f1*dr.x*(posB.x-comB.x);                            // xx
          Gradient[n+1]+=f1*0.5*(dr.x*(posB.y-comB.y)+dr.y*(posB.x-comB.x)); // xy
          Gradient[n+2]+=f1*0.5*(dr.x*(posB.z-comB.z)+dr.z*(posB.x-comB.x)); // xz
          Gradient[n+3]+=f1*(posB.y-comB.y)*dr.y;                            // yy
          Gradient[n+4]+=f1*0.5*(dr.y*(posB.z-comB.z)+dr.z*(posB.y-comB.y)); // yz
          Gradient[n+5]+=f1*dr.z*(posB.z-comB.z);                            // zz
          break;
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n  ]+=f1*dr.x*(posB.x-comB.x); // xx
          Gradient[n+1]+=f1*dr.x*(posB.y-comB.y); // xy
          Gradient[n+2]+=f1*dr.x*(posB.z-comB.z); // xz
          Gradient[n+3]+=f1*dr.y*(posB.y-comB.y); // yy
          Gradient[n+4]+=f1*dr.y*(posB.z-comB.z); // yz
          Gradient[n+5]+=f1*dr.z*(posB.z-comB.z); // zz
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x);                            // xx
              Gradient[n+1]+=f1*(posB.y-comB.y)*dr.y;                            // yy
              Gradient[n+2]+=f1*0.5*(dr.y*(posB.z-comB.z)+dr.z*(posB.y-comB.y)); // yz
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z);                            // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x);                            // xx
              Gradient[n+1]+=f1*0.5*(dr.x*(posB.z-comB.z)+dr.z*(posB.x-comB.x)); // xz
              Gradient[n+2]+=f1*(posB.y-comB.y)*dr.y;                            // yy
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z);                            // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x);                            // xx
              Gradient[n+1]+=f1*0.5*(dr.x*(posB.y-comB.y)+dr.y*(posB.x-comB.x)); // xy
              Gradient[n+2]+=f1*(posB.y-comB.y)*dr.y;                            // yy
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z);                            // zz
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x); // xx
              Gradient[n+1]+=f1*dr.y*(posB.y-comB.y); // yy
              Gradient[n+2]+=f1*dr.y*(posB.z-comB.z); // yz
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z); // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x); // xx
              Gradient[n+1]+=f1*dr.x*(posB.z-comB.z); // xz
              Gradient[n+2]+=f1*dr.y*(posB.y-comB.y); // yy
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z); // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]+=f1*dr.x*(posB.x-comB.x); // xx
              Gradient[n+1]+=f1*dr.x*(posB.y-comB.y); // xy
              Gradient[n+2]+=f1*dr.y*(posB.y-comB.y); // yy
              Gradient[n+3]+=f1*dr.z*(posB.z-comB.z); // zz
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

void ComputeFrameworkAdsorbateVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                         REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr;
  REAL energy,f1,f2;
  VECTOR posA,posB,dr;
  int index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR pos,comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;

  f1=f2=0.0;
  index2=0;
  // first loop over adsorbate molecules

  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
    {
      index_i=Framework[CurrentSystem].Atoms[fr][i].HessianIndex;
      posA=Framework[CurrentSystem].Atoms[fr][i].Position;
      typeA=Framework[CurrentSystem].Atoms[fr][i].Type;
      ChargeA=Framework[CurrentSystem].Atoms[fr][i].Charge;

        // second loop over adsorbates
      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
            {
              TypeMolB=Adsorbates[CurrentSystem][J].Type;
              for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
              {
                for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                {
                  j=Components[TypeMolB].Groups[jg].Atoms[ja];
 
                  if(Components[TypeMolB].Groups[jg].Rigid)
                  {
                    index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                    index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                    index2=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                  }
                  else
                  {
                    index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                    index_j2=-1;
                  }

                  typeB=Adsorbates[CurrentSystem][J].Atoms[j].Type;
                  posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;
                  ChargeB=Adsorbates[CurrentSystem][J].Atoms[j].Charge;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffVDWSquared)
                  {
                    PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2);

                    // add contribution to the energy
                    *Energy+=energy;

                    if((index_i<0)&&(index_j<0)&&(index_j2<0)) continue;

                    StrainDerivative->ax+=f1*dr.x*dr.x;
                    StrainDerivative->bx+=f1*dr.x*dr.y;
                    StrainDerivative->cx+=f1*dr.x*dr.z;

                    StrainDerivative->ay+=f1*dr.y*dr.x;
                    StrainDerivative->by+=f1*dr.y*dr.y;
                    StrainDerivative->cy+=f1*dr.y*dr.z;

                    StrainDerivative->az+=f1*dr.z*dr.x;
                    StrainDerivative->bz+=f1*dr.z*dr.y;
                    StrainDerivative->cz+=f1*dr.z*dr.z;

                    if(Components[TypeMolB].Groups[jg].Rigid)
                    {
                      comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                      comB.x+=ReplicaShift[ncell].x;
                      comB.y+=ReplicaShift[ncell].y;
                      comB.z+=ReplicaShift[ncell].z;

                      pos=Components[TypeMolB].Positions[j];

                      temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                      temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                      temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                      StrainDerivative->ax+=(posB.x-comB.x)*f1*dr.x;
                      StrainDerivative->ay+=temp1;
                      StrainDerivative->az+=temp2;
                      StrainDerivative->bx+=temp1;
                      StrainDerivative->by+=(posB.y-comB.y)*f1*dr.y;
                      StrainDerivative->bz+=temp3;
                      StrainDerivative->cx+=temp2;
                      StrainDerivative->cy+=temp3;
                      StrainDerivative->cz+=(posB.z-comB.z)*f1*dr.z;
                    }

                    // add contribution to the first derivatives
                    if(ComputeGradient)
                    {
                      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      {
                        if(index_i>=0)
                        {
                          Gradient[index_i]+=f1*dr.x;
                          Gradient[index_i+1]+=f1*dr.y;
                          Gradient[index_i+2]+=f1*dr.z;
                        }
                      }
  
                      if(index_j>=0)
                      {
                        Gradient[index_j]-=f1*dr.x;
                        Gradient[index_j+1]-=f1*dr.y;
                        Gradient[index_j+2]-=f1*dr.z;
                      }

                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                        // add contribution to the first derivatives
                        if(index_j2>=0)
                        {
                          Gradient[index_j2]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                          Gradient[index_j2+1]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                          Gradient[index_j2+2]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                        }
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,1.0);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,f1,f2,dr);
                      HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,f1,f2,dr);


                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianCenterOfMassStrainJ(HessianMatrix,index_i,index_j,f1,f2,dr,posB,comB);
                        HessianOrientationStrainJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianOrientationStrainJ_FR(HessianMatrix,index_j2,index2,f1,f2,posB,comB,dr);

                        HessianCorrectionStrainStrain(HessianMatrix,f1,f2,dr,posB,comB);
                      }
                    } 
                  }
                }
              }
            }
            ncell++;
          }
    }
  }
}

void ComputeFrameworkAdsorbateChargeChargeHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                                  REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr;
  REAL f1,f2;
  VECTOR posA,posB,dr;
  int index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;

  f1=f2=0.0;
  index2=0;
  // first loop over adsorbate molecules
  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
    {
      index_i=Framework[CurrentSystem].Atoms[fr][i].HessianIndex;
      posA=Framework[CurrentSystem].Atoms[fr][i].Position;
      typeA=Framework[CurrentSystem].Atoms[fr][i].Type;
      ChargeA=Framework[CurrentSystem].Atoms[fr][i].Charge;

      // second loop over adsorbates
      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
            {
              TypeMolB=Adsorbates[CurrentSystem][J].Type;
              for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
              {
                for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                {
                  j=Components[TypeMolB].Groups[jg].Atoms[ja];

                  if(Components[TypeMolB].Groups[jg].Rigid)
                  {
                    index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                    index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                    index2=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                  }
                  else
                  {
                    index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                    index_j2=-1;
                  }

                  typeB=Adsorbates[CurrentSystem][J].Atoms[j].Type;
                  posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;
                  ChargeB=Adsorbates[CurrentSystem][J].Atoms[j].Charge;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffChargeChargeSquared)
                  {
                    switch(ChargeMethod)
                    {
                      case NONE:
                        f1=f2=0.0;
                        break;
                      case SHIFTED_COULOMB:
                        *Energy+=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r-InverseCutOffChargeCharge);
                        f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                        f2=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                        break;
                      case TRUNCATED_COULOMB:
                        *Energy+=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r);
                        f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                        f2=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                        break;
                      case EWALD:
                      default:
                        *Energy+=COULOMBIC_CONVERSION_FACTOR*erfc(Alpha[CurrentSystem]*r)*
                                                ChargeA*ChargeB/r;
                        f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                            (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                            (r*rr);

                        f2=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                           (((3.0*erfc(Alpha[CurrentSystem]*r)/(r*rr))+
                           (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                           (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));
                        break;
                    }

                    if((index_i<0)&&(index_j<0)) continue;

                    StrainDerivative->ax+=f1*dr.x*dr.x;
                    StrainDerivative->bx+=f1*dr.x*dr.y;
                    StrainDerivative->cx+=f1*dr.x*dr.z;

                    StrainDerivative->ay+=f1*dr.y*dr.x;
                    StrainDerivative->by+=f1*dr.y*dr.y;
                    StrainDerivative->cy+=f1*dr.y*dr.z;

                    StrainDerivative->az+=f1*dr.z*dr.x;
                    StrainDerivative->bz+=f1*dr.z*dr.y;
                    StrainDerivative->cz+=f1*dr.z*dr.z;

                    if(Components[TypeMolB].Groups[jg].Rigid)
                    {
                      comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                      comB.x+=ReplicaShift[ncell].x;
                      comB.y+=ReplicaShift[ncell].y;
                      comB.z+=ReplicaShift[ncell].z;

                      temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                      temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                      temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                      StrainDerivative->ax+=(posB.x-comB.x)*f1*dr.x;
                      StrainDerivative->ay+=temp1;
                      StrainDerivative->az+=temp2;
                      StrainDerivative->bx+=temp1;
                      StrainDerivative->by+=(posB.y-comB.y)*f1*dr.y;
                      StrainDerivative->bz+=temp3;
                      StrainDerivative->cx+=temp2;
                      StrainDerivative->cy+=temp3;
                      StrainDerivative->cz+=(posB.z-comB.z)*f1*dr.z;
                    }

                    // add contribution to the first derivatives
                    if(ComputeGradient)
                    {
                      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      {
                        if(index_i>=0)
                        {
                          Gradient[index_i]+=f1*dr.x;
                          Gradient[index_i+1]+=f1*dr.y;
                          Gradient[index_i+2]+=f1*dr.z;
                        }
                      }

                      if(index_j>=0)
                      {
                        Gradient[index_j]-=f1*dr.x;
                        Gradient[index_j+1]-=f1*dr.y;
                        Gradient[index_j+2]-=f1*dr.z;
                      }
  
                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                         // add contribution to the first derivatives
                        if(index_j2>=0)
                        {
                          Gradient[index_j2]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                          Gradient[index_j2+1]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                          Gradient[index_j2+2]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                        }
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,1.0);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,f1,f2,dr);
                      HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,f1,f2,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianCenterOfMassStrainJ(HessianMatrix,index_i,index_j,f1,f2,dr,posB,comB);
                        HessianOrientationStrainJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianOrientationStrainJ_FR(HessianMatrix,index_j2,index2,f1,f2,posB,comB,dr);

                        HessianCorrectionStrainStrain(HessianMatrix,f1,f2,dr,posB,comB);
                      }
                    }
                  }
                }
              }
            }
            ncell++;
          }
    }
  }
}

void ComputeFrameworkCationVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                      REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr;
  REAL energy,f1,f2;
  VECTOR posA,posB,dr;
  int index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;

  f1=f2=0.0;
  index2=0;
  // first loop over adsorbate molecules

  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
    {
      index_i=Framework[CurrentSystem].Atoms[fr][i].HessianIndex;
      posA=Framework[CurrentSystem].Atoms[fr][i].Position;
      typeA=Framework[CurrentSystem].Atoms[fr][i].Type;
      ChargeA=Framework[CurrentSystem].Atoms[fr][i].Charge;

      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
            {
              TypeMolB=Cations[CurrentSystem][J].Type;
              for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
              {
                for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                {
                  j=Components[TypeMolB].Groups[jg].Atoms[ja];

                  if(Components[TypeMolB].Groups[jg].Rigid)
                  {
                    index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                    index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                    index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                  }
                  else
                  {
                    index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                    index_j2=-1;
                  }

                  typeB=Cations[CurrentSystem][J].Atoms[j].Type;
                  posB=Cations[CurrentSystem][J].Atoms[j].Position;
                  ChargeB=Cations[CurrentSystem][J].Atoms[j].Charge;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffVDWSquared)
                  {
                    PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2);

                    // add contribution to the energy
                    *Energy+=energy;

                    StrainDerivative->ax+=f1*dr.x*dr.x;
                    StrainDerivative->bx+=f1*dr.x*dr.y;
                    StrainDerivative->cx+=f1*dr.x*dr.z;

                    StrainDerivative->ay+=f1*dr.y*dr.x;
                    StrainDerivative->by+=f1*dr.y*dr.y;
                    StrainDerivative->cy+=f1*dr.y*dr.z;

                    StrainDerivative->az+=f1*dr.z*dr.x;
                    StrainDerivative->bz+=f1*dr.z*dr.y;
                    StrainDerivative->cz+=f1*dr.z*dr.z;

                    if(Components[TypeMolB].Groups[jg].Rigid)
                    {
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                      comB.x+=ReplicaShift[ncell].x;
                      comB.y+=ReplicaShift[ncell].y;
                      comB.z+=ReplicaShift[ncell].z;

                      temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                      temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                      temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                      StrainDerivative->ax+=(posB.x-comB.x)*f1*dr.x;
                      StrainDerivative->ay+=temp1;
                      StrainDerivative->az+=temp2;
                      StrainDerivative->bx+=temp1;
                      StrainDerivative->by+=(posB.y-comB.y)*f1*dr.y;
                      StrainDerivative->bz+=temp3;
                      StrainDerivative->cx+=temp2;
                      StrainDerivative->cy+=temp3;
                      StrainDerivative->cz+=(posB.z-comB.z)*f1*dr.z;
                    }

                    // add contribution to the first derivatives
                    if(ComputeGradient)
                    {
                      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      {
                        if(index_i>=0)
                        {
                          Gradient[index_i]+=f1*dr.x;
                          Gradient[index_i+1]+=f1*dr.y;
                          Gradient[index_i+2]+=f1*dr.z;
                        }
                      }

                      if(index_j>=0)
                      {
                        Gradient[index_j]-=f1*dr.x;
                        Gradient[index_j+1]-=f1*dr.y;
                        Gradient[index_j+2]-=f1*dr.z;
                      }

                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                        // add contribution to the first derivatives
                        if(index_j2>=0)
                        {
                          Gradient[index_j2]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                          Gradient[index_j2+1]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                          Gradient[index_j2+2]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                        }
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,1.0);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,f1,f2,dr);
                      HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,f1,f2,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianCenterOfMassStrainJ(HessianMatrix,index_i,index_j,f1,f2,dr,posB,comB);
                        HessianOrientationStrainJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianOrientationStrainJ_FR(HessianMatrix,index_j2,index2,f1,f2,posB,comB,dr);

                        HessianCorrectionStrainStrain(HessianMatrix,f1,f2,dr,posB,comB);
                      }
                    }
                  }
                }
              }
            }
            ncell++;
          }
    }
  }
}

void ComputeFrameworkCationChargeChargeHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                               REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr;
  REAL f1,f2;
  VECTOR posA,posB,dr;
  int index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;

  f1=f2=0.0;
  index2=0;
  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
    {
      index_i=Framework[CurrentSystem].Atoms[fr][i].HessianIndex;
      posA=Framework[CurrentSystem].Atoms[fr][i].Position;
      typeA=Framework[CurrentSystem].Atoms[fr][i].Type;
      ChargeA=Framework[CurrentSystem].Atoms[fr][i].Charge;

      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
            {
              TypeMolB=Cations[CurrentSystem][J].Type;
              for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
              {
                for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                {
                  j=Components[TypeMolB].Groups[jg].Atoms[ja];

                  if(Components[TypeMolB].Groups[jg].Rigid)
                  {
                    index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                    index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                    index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                  }
                  else
                  {
                    index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                    index_j2=-1;
                  }

                  typeB=Cations[CurrentSystem][J].Atoms[j].Type;
                  posB=Cations[CurrentSystem][J].Atoms[j].Position;
                  ChargeB=Cations[CurrentSystem][J].Atoms[j].Charge;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffChargeChargeSquared)
                  {
                    switch(ChargeMethod)
                    {
                      case NONE:
                        f1=f2=0.0;
                        break;
                      case SHIFTED_COULOMB:
                        *Energy+=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r-InverseCutOffChargeCharge);
                        f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                        f2=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                        break;
                      case TRUNCATED_COULOMB:
                        *Energy+=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r);
                        f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                        f2=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                        break;
                      case EWALD:
                      default:
                        *Energy+=COULOMBIC_CONVERSION_FACTOR*erfc(Alpha[CurrentSystem]*r)*
                                                ChargeA*ChargeB/r;
                        f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                            (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                            (r*rr);

                        f2=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                           (((3.0*erfc(Alpha[CurrentSystem]*r)/(r*rr))+
                           (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                           (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));
                        break;
                    }

                    if((index_j<0)&&(index_j2<0)) continue;

                    StrainDerivative->ax+=f1*dr.x*dr.x;
                    StrainDerivative->bx+=f1*dr.x*dr.y;
                    StrainDerivative->cx+=f1*dr.x*dr.z;

                    StrainDerivative->ay+=f1*dr.y*dr.x;
                    StrainDerivative->by+=f1*dr.y*dr.y;
                    StrainDerivative->cy+=f1*dr.y*dr.z;

                    StrainDerivative->az+=f1*dr.z*dr.x;
                    StrainDerivative->bz+=f1*dr.z*dr.y;
                    StrainDerivative->cz+=f1*dr.z*dr.z;

                    if(Components[TypeMolB].Groups[jg].Rigid)
                    {
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                      comB.x+=ReplicaShift[ncell].x;
                      comB.y+=ReplicaShift[ncell].y;
                      comB.z+=ReplicaShift[ncell].z;

                      temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                      temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                      temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                      StrainDerivative->ax+=(posB.x-comB.x)*f1*dr.x;
                      StrainDerivative->ay+=temp1;
                      StrainDerivative->az+=temp2;
                      StrainDerivative->bx+=temp1;
                      StrainDerivative->by+=(posB.y-comB.y)*f1*dr.y;
                      StrainDerivative->bz+=temp3;
                      StrainDerivative->cx+=temp2;
                      StrainDerivative->cy+=temp3;
                      StrainDerivative->cz+=(posB.z-comB.z)*f1*dr.z;
                    }

                    // add contribution to the first derivatives
                    if(ComputeGradient)
                    {
                      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      {
                        if(index_i>=0)
                        {
                          Gradient[index_i]+=f1*dr.x;
                          Gradient[index_i+1]+=f1*dr.y;
                          Gradient[index_i+2]+=f1*dr.z;
                        }
                      }

                      if(index_j>=0)
                      {
                        Gradient[index_j]-=f1*dr.x;
                        Gradient[index_j+1]-=f1*dr.y;
                        Gradient[index_j+2]-=f1*dr.z;
                      }

                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                        // add contribution to the first derivatives
                        if(index_j2>=0)
                        {
                          Gradient[index_j2]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                          Gradient[index_j2+1]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                          Gradient[index_j2+2]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                        }
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,f1,f2,dr,1.0);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,f1,f2,dr);
                      HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,f1,f2,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianCenterOfMassStrainJ(HessianMatrix,index_i,index_j,f1,f2,dr,posB,comB);
                        HessianOrientationStrainJ(HessianMatrix,index_j2,index2,f1,f2,dr);
                        HessianOrientationStrainJ_FR(HessianMatrix,index_j2,index2,f1,f2,posB,comB,dr);

                        HessianCorrectionStrainStrain(HessianMatrix,f1,f2,dr,posB,comB);
                      }
                    }
                  }
                }
              }
            }
            ncell++;
          }
    }
  }
}


void ComputeFrameworkIntraVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                     REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,j,typeA,typeB,f1,f2,start;
  REAL ChargeA,ChargeB,energy,DF,DDF;
  REAL rr,r;
  VECTOR posA,posB,dr;
  int index_i,index_j;
  int ncell,k1,k2,k3,indj;
  REAL ReplicaFactor;

  if(!InternalFrameworkLennardJonesInteractions) return;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        ChargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;
        index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if((f1==f2)&&(ncell==0)) start=i+1;
              else start=0;
              for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
              {
                indj=j+ncell*Framework[CurrentSystem].NumberOfAtoms[f2];
                if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][indj],0))
                {
                  typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
                  posB=Framework[CurrentSystem].Atoms[f2][j].Position;
                  ChargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;
                  index_j=Framework[CurrentSystem].Atoms[f2][j].HessianIndex;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffVDWSquared)
                  {
                    if(ncell==0) ReplicaFactor=1.0;
                    else ReplicaFactor=0.5;

                    PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF);

                    // add contribution to the energy
                    *Energy+=ReplicaFactor*energy;

                    if((index_i<0)&&(index_j)<0) continue;

                    if(ComputeGradient)
                    {
                      // add contribution to the first derivatives
                      if(index_i>=0)
                      {
                        Gradient[index_i]+=DF*dr.x;
                        Gradient[index_i+1]+=DF*dr.y;
                        Gradient[index_i+2]+=DF*dr.z;
                      }

                      if(ncell==0)
                      {
                        if(index_j>=0)
                        {
                          Gradient[index_j]-=DF*dr.x;
                          Gradient[index_j+1]-=DF*dr.y;
                          Gradient[index_j+2]-=DF*dr.z;
                        }
                      }

                      GradientStrain(Gradient,ReplicaFactor*DF,dr);
                    }

                    // add contribution to the strain derivative tensor
                    StrainDerivative->ax+=ReplicaFactor*DF*dr.x*dr.x;
                    StrainDerivative->bx+=ReplicaFactor*DF*dr.y*dr.x;
                    StrainDerivative->cx+=ReplicaFactor*DF*dr.z*dr.x;

                    StrainDerivative->ay+=ReplicaFactor*DF*dr.x*dr.y;
                    StrainDerivative->by+=ReplicaFactor*DF*dr.y*dr.y;
                    StrainDerivative->cy+=ReplicaFactor*DF*dr.z*dr.y;

                    StrainDerivative->az+=ReplicaFactor*DF*dr.x*dr.z;
                    StrainDerivative->bz+=ReplicaFactor*DF*dr.y*dr.z;
                    StrainDerivative->cz+=ReplicaFactor*DF*dr.z*dr.z;

                    if(ComputeHessian)
                    {
                      // add contribution to the second derivatives (Hessian matrix)
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,ReplicaFactor);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,ReplicaFactor*DF,ReplicaFactor*DDF,dr);
                      HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,ReplicaFactor*DF,ReplicaFactor*DDF,dr);
                    }
                  }
                }
              }
              ncell++;
            }
      }
    }
  }
}

void ComputeFrameworkIntraChargeChargeHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,
                                              REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,j,typeA,typeB,f1,f2,start;
  REAL ChargeA,ChargeB,DF,DDF;
  REAL rr,r;
  VECTOR posA,posB,dr;
  int index_i,index_j;
  int ncell,k1,k2,k3,indj;
  REAL ReplicaFactor;

  if(ChargeMethod==NONE) return;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        ChargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;
        index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if((f1==f2)&&(ncell==0)) start=i+1;
              else start=0;
              for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
              {
                indj=j+ncell*Framework[CurrentSystem].NumberOfAtoms[f2];
                if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][indj],0))
                {
                  typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
                  posB=Framework[CurrentSystem].Atoms[f2][j].Position;
                  ChargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;
                  index_j=Framework[CurrentSystem].Atoms[f2][j].HessianIndex;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffChargeChargeSquared)
                  {
                    if(ncell==0) ReplicaFactor=1.0;
                    else ReplicaFactor=0.5;

                    switch(ChargeMethod)
                    {
                      case NONE:
                        DF=DDF=0.0;
                        break;
                      case SHIFTED_COULOMB:
                        *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r-InverseCutOffChargeCharge);
                        DF=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                        DDF=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                        break;
                      case TRUNCATED_COULOMB:
                        *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r);
                        DF=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                        DDF=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                        break;
                      case EWALD:
                      default:
                        *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*erfc(Alpha[CurrentSystem]*r)*
                                                ChargeA*ChargeB/r;
                        DF=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                            (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                            (r*rr);

                        DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                            (((3.0*erfc(Alpha[CurrentSystem]*r)/(r*rr))+
                             (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                             (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));
                        break;
                    }

                    if((index_i<0)&&(index_j)<0) continue;

                    if(ComputeGradient)
                    {
                      // add contribution to the first derivatives
                      if(index_i>=0)
                      {
                        Gradient[index_i]+=DF*dr.x;
                        Gradient[index_i+1]+=DF*dr.y;
                        Gradient[index_i+2]+=DF*dr.z;
                      }

                      if(ncell==0)
                      {
                        if(index_j>=0)
                        {
                          Gradient[index_j]-=DF*dr.x;
                          Gradient[index_j+1]-=DF*dr.y;
                          Gradient[index_j+2]-=DF*dr.z;
                        }
                      }

                      GradientStrain(Gradient,ReplicaFactor*DF,dr);
                    }

                    // add contribution to the strain derivative tensor
                    StrainDerivative->ax+=ReplicaFactor*DF*dr.x*dr.x;
                    StrainDerivative->bx+=ReplicaFactor*DF*dr.y*dr.x;
                    StrainDerivative->cx+=ReplicaFactor*DF*dr.z*dr.x;

                    StrainDerivative->ay+=ReplicaFactor*DF*dr.x*dr.y;
                    StrainDerivative->by+=ReplicaFactor*DF*dr.y*dr.y;
                    StrainDerivative->cy+=ReplicaFactor*DF*dr.z*dr.y;

                    StrainDerivative->az+=ReplicaFactor*DF*dr.x*dr.z;
                    StrainDerivative->bz+=ReplicaFactor*DF*dr.y*dr.z;
                    StrainDerivative->cz+=ReplicaFactor*DF*dr.z*dr.z;

                    if(ComputeHessian)
                    {
                      // add contribution to the second derivatives (Hessian matrix)
                      HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,ReplicaFactor);
                      HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,ReplicaFactor*DF,ReplicaFactor*DDF,dr);
                      HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,ReplicaFactor*DF,ReplicaFactor*DDF,dr);
                    }
                  }
                }
              }
              ncell++;
            }
      }
    }
  }
}



void ComputeFrameworkBondHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i;   // loop variable
  int f1;  // loop over all frameworks
  REAL r;  // distance
  REAL rr; // distance squared
  REAL temp,temp2; // temporary
  REAL exp_term;   // temporary
  REAL U;  // energy of a specific interaction
  REAL DF; // first derivative
  REAL DDF;  // second derivative
  VECTOR dr; // atoms separation vector
  int A,B;   // atom-indices
  int index_i,index_j; // indices of the Hessian
  REAL *parms;  // pointer to potential parameter

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBonds[f1];i++)
    {
      A=Framework[CurrentSystem].Bonds[f1][i].A;
      B=Framework[CurrentSystem].Bonds[f1][i].B;

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;

      dr.x=Framework[CurrentSystem].Atoms[f1][A].Position.x-
           Framework[CurrentSystem].Atoms[f1][B].Position.x;
      dr.y=Framework[CurrentSystem].Atoms[f1][A].Position.y-
           Framework[CurrentSystem].Atoms[f1][B].Position.y;
      dr.z=Framework[CurrentSystem].Atoms[f1][A].Position.z-
           Framework[CurrentSystem].Atoms[f1][B].Position.z;

      // apply boundary condition
      dr=ApplyBoundaryCondition(dr);

      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      parms=(REAL*)&Framework[CurrentSystem].BondArguments[f1][i];

      switch(Framework[CurrentSystem].BondType[f1][i])
      {
        case HARMONIC_BOND:
          // 0.5*p0*SQR(r-p1);
          // ===============================================
          // p_0/k_B [K/A^2]   force constant
          // p_1     [A]       reference bond distance
          U=0.5*parms[0]*SQR(r-parms[1]);
          DF=parms[0]*(r-parms[1])/r;
          DDF=(parms[0]*parms[1])/(r*rr);
          //DDF=(parms[0]/r-DF)/rr;
          break;
        case CORE_SHELL_SPRING:
          U=0.5*parms[0]*SQR(r);
          DF=parms[0];
          DDF=0.0;
          break;
        case MORSE_BOND:
          // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
          // ===============================================
          // p_0/k_B [K]       force constant
          // p_1     [A^-1]    parameter
          // p_2     [A]       reference bond distance
          temp=exp(parms[1]*(parms[2]-r));
          U=parms[0]*(SQR(1.0-temp)-1.0);
          DF=2.0*parms[0]*parms[1]*(1.0-temp)*temp/r;
          DDF=2.0*parms[0]*parms[1]*temp*((1.0+2.0*parms[1]*r)*temp-parms[1]*r-1.0)/(r*rr);
          break;
        case LJ_12_6_BOND:
          // A/r_ij^12-B/r_ij^6
          // ===============================================
          // p_0/k_B [K A^12]
          // p_1/k_B [K A^6]
          temp=CUBE(1.0/rr);
          U=parms[0]*SQR(temp)-parms[1]*temp;
          DF=6.0*(parms[1]*temp-2.0*parms[0]*SQR(temp))/rr;
          DDF=24.0*(7.0*parms[0]*SQR(temp)-2.0*parms[1]*temp)/SQR(rr);
          break;
        case LENNARD_JONES_BOND:
          // 4*p_0*((p_1/r)^12-(p_1/r)^6)
          // ===============================================
          // p_0/k_B [K]
          // p_1     [A]
          temp=CUBE(parms[1]/rr);
          U=4.0*parms[0]*(temp*(temp-1.0));
          DF=24.0*parms[0]*(temp*(1.0-2.0*temp))/rr;
          DDF=96.0*parms[0]*(temp*(7.0*temp-2.0))/SQR(rr);
          break;
        case BUCKINGHAM_BOND:
          // p_0*exp(-p_1 r)-p_2/r^6
          // ===============================================
          // p_0/k_B [K]
          // p_1     [A^-1]
          // p_2/k_B [K A^6]
          temp=parms[2]*CUBE(1.0/rr);
          exp_term=parms[0]*exp(-parms[1]*r);
          U=-temp+exp_term;
          DF=(6.0/rr)*temp-parms[1]*exp_term/r;
          DDF=(-48.0*temp/rr+parms[1]*(1.0+parms[1]*r)*exp_term/r)/rr;
          break;
        case RESTRAINED_HARMONIC_BOND:
          // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
          // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
          // ===============================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2     [A]
          temp=r-parms[1];
          U=0.5*parms[0]*SQR(MIN2(fabs(temp),parms[2]))
                +parms[0]*parms[2]*MAX2(fabs(temp)-parms[2],(REAL)0.0);
          DF=-parms[0]*(SIGN(MIN2(fabs(temp),parms[2]),temp))/r;
          DDF=fabs(temp)<parms[2]?-parms[0]*parms[1]/(r*rr):parms[0]*SIGN(parms[2],temp)/(r*rr);
          break;
        case QUARTIC_BOND:
          // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
          // ===========================================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2/k_B [K/A^3]
          // p_3/k_B [K/A^4]
          temp=r-parms[1];
          temp2=SQR(r-parms[1]);
          U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
          DF=temp*(parms[0]+parms[2]*temp+parms[3]*temp2)/r;
          DDF=2.0*parms[3]+(parms[2]-3.0*parms[1]*parms[3])/r+(parms[1]*(parms[0]+parms[1]*(parms[1]*parms[3]-parms[2])))/(r*rr);
          break;
        case CFF_QUARTIC_BOND:
          // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
          // ===============================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2/k_B [K/A^3]
          // p_3/k_B [K/A^4]
          temp=r-parms[1];
          temp2=SQR(r-parms[1]);
          U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
          DF=temp*(2.0*parms[0]+3.0*parms[2]*temp+4.0*parms[3]*temp2)/r;
          DDF=8.0*parms[3]+(3.0*parms[2]-3.0*parms[1]*4.0*parms[3])/r+(parms[1]*(2.0*parms[0]+parms[1]*(parms[1]*4.0*parms[3]-3.0*parms[2])))/(r*rr);
          break;
        case MM3_BOND:
          // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
          // =================================================================
          // p_0     [mdyne/A molecule]
          // p_1     [A]
          temp=r-parms[1];
          temp2=SQR(temp);
          U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
          DF=parms[0]*(2.0+2.55*(4.0*2.55*(7.0/12.0)*temp-3.0)*temp)*temp/r;
          DDF=(parms[0]*(SQR(2.55)*4.0*7.0*temp2*(parms[1]+2.0*r)+12.0*(2.0*parms[1]+2.55*3.0*(SQR(parms[1])-SQR(r)))))/(12.0*SQR(r)*r);
          break;
        case RIGID_BOND:
          U=DF=DDF=0.0;
          break;
        case FIXED_BOND:
          U=DF=DDF=0.0;
          break;
        case MEASURE_BOND:
          U=DF=DDF=0.0;
          break;
        default:
          printf("Undefined Bond potential in routine 'CalculateFrameworkBondHessian' ('framework_hessian.c')\n");
          exit(0);
          break;
      }

      // add contribution to the energy
      *Energy+=U;

      if((index_i<0)&&(index_j<0)) continue;

      if(ComputeGradient)
      {
        // add contribution to the first derivatives
        if(index_i>=0)
        {
          Gradient[index_i]+=DF*dr.x;
          Gradient[index_i+1]+=DF*dr.y;
          Gradient[index_i+2]+=DF*dr.z;
        }

        if(index_j>=0)
        {
          Gradient[index_j]-=DF*dr.x;
          Gradient[index_j+1]-=DF*dr.y;
          Gradient[index_j+2]-=DF*dr.z;
        }

        GradientStrain(Gradient,DF,dr);
      }

      // add contribution to the strain derivative tensor
      StrainDerivative->ax+=dr.x*DF*dr.x;
      StrainDerivative->bx+=dr.y*DF*dr.x;
      StrainDerivative->cx+=dr.z*DF*dr.x;

      StrainDerivative->ay+=dr.x*DF*dr.y;
      StrainDerivative->by+=dr.y*DF*dr.y;
      StrainDerivative->cy+=dr.z*DF*dr.y;

      StrainDerivative->az+=dr.x*DF*dr.z;
      StrainDerivative->bz+=dr.y*DF*dr.z;
      StrainDerivative->cz+=dr.z*DF*dr.z;

      if(ComputeHessian)
      {
        // add contribution to the second derivatives (Hessian matrix)
        HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
        HessianAtomicPositionStrain(HessianMatrix,index_i,index_j,DF,DDF,dr);
        HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,DF,DDF,dr);
      }
    }
  }
}



static inline void HessianBendStrainPosition(int index_i,int index_j,int index_k,REAL_MATRIX HessianMatrix,VECTOR vec_u,VECTOR vec_v,REAL u,REAL v,
           REAL rab,REAL rbc,VECTOR Rab,VECTOR Rbc,VECTOR dtA,VECTOR dtC,REAL DF,REAL DDF,REAL_MATRIX3x3 S,REAL CosTheta)
{
  int n;
  VECTOR veci,veck;
  VECTOR vec_ddtax_i,vec_ddtax_k;
  VECTOR vec_ddtay_i,vec_ddtay_k;
  VECTOR vec_ddtaz_i,vec_ddtaz_k;
  VECTOR vec_ddtcx_i,vec_ddtcx_k;
  VECTOR vec_ddtcy_i,vec_ddtcy_k;
  VECTOR vec_ddtcz_i,vec_ddtcz_k;

  n=NumberOfCoordinatesMinimizationVariables;

  veci.x=(1.0/CUBE(u))*(vec_u.x*vec_u.x*Rbc.x+vec_u.y*vec_u.x*Rbc.y+vec_u.z*vec_u.x*Rbc.z)-Rbc.x/u;
  veci.y=(1.0/CUBE(u))*(vec_u.x*vec_u.y*Rbc.x+vec_u.y*vec_u.y*Rbc.y+vec_u.z*vec_u.y*Rbc.z)-Rbc.y/u;
  veci.z=(1.0/CUBE(u))*(vec_u.x*vec_u.z*Rbc.x+vec_u.y*vec_u.z*Rbc.y+vec_u.z*vec_u.z*Rbc.z)-Rbc.z/u;

  veck.x=(1.0/CUBE(v))*(vec_v.x*vec_v.x*Rab.x+vec_v.y*vec_v.x*Rab.y+vec_v.z*vec_v.x*Rab.z)-Rab.x/v;
  veck.y=(1.0/CUBE(v))*(vec_v.x*vec_v.y*Rab.x+vec_v.y*vec_v.y*Rab.y+vec_v.z*vec_v.y*Rab.z)-Rab.y/v;
  veck.z=(1.0/CUBE(v))*(vec_v.x*vec_v.z*Rab.x+vec_v.y*vec_v.z*Rab.y+vec_v.z*vec_v.z*Rab.z)-Rab.z/v;


  // FINAL PART 1: derivative dtA.x=(Rbc.x-CosTheta*Rab.x)/rab
  vec_ddtax_i.x=(veci.x*Rab.x+CosTheta*((1.0/CUBE(u))*(vec_u.x*vec_u.x)-1/u))/rab+CosTheta*Rab.x*(1.0/CUBE(u))*vec_u.x-(1.0/CUBE(u))*vec_u.x*Rbc.x;
  vec_ddtax_i.y=(veci.y*Rab.x+CosTheta*(1.0/CUBE(u))*(vec_u.x*vec_u.y))/rab+CosTheta*Rab.x*(1.0/CUBE(u))*vec_u.y-(1.0/CUBE(u))*vec_u.y*Rbc.x;
  vec_ddtax_i.z=(veci.z*Rab.x+CosTheta*(1.0/CUBE(u))*(vec_u.x*vec_u.z))/rab+CosTheta*Rab.x*(1.0/CUBE(u))*vec_u.z-(1.0/CUBE(u))*vec_u.z*Rbc.x;

  vec_ddtax_k.x=veck.x*Rab.x/rab-((1.0/CUBE(v))*(vec_v.x*vec_v.x)-1/v)/rab;
  vec_ddtax_k.y=veck.y*Rab.x/rab-(1.0/CUBE(v))*(vec_v.x*vec_v.y)/rab;
  vec_ddtax_k.z=veck.z*Rab.x/rab-(1.0/CUBE(v))*(vec_v.x*vec_v.z)/rab;

  // FINAL PART 1: derivative dtA.y=(Rbc.y-CosTheta*Rab.y)/rab
  vec_ddtay_i.x=(veci.x*Rab.y+CosTheta*(1.0/CUBE(u))*(vec_u.y*vec_u.x))/rab+CosTheta*Rab.y*(1.0/CUBE(u))*vec_u.x-(1.0/CUBE(u))*vec_u.x*Rbc.y;
  vec_ddtay_i.y=(veci.y*Rab.y+CosTheta*((1.0/CUBE(u))*(vec_u.y*vec_u.y)-1.0/u))/rab+CosTheta*Rab.y*(1.0/CUBE(u))*vec_u.y-(1.0/CUBE(u))*vec_u.y*Rbc.y;
  vec_ddtay_i.z=(veci.z*Rab.y+CosTheta*(1.0/CUBE(u))*(vec_u.y*vec_u.z))/rab+CosTheta*Rab.y*(1.0/CUBE(u))*vec_u.z-(1.0/CUBE(u))*vec_u.z*Rbc.y;

  vec_ddtay_k.x=veck.x*Rab.y/rab-(1.0/CUBE(v))*(vec_v.y*vec_v.x)/rab;
  vec_ddtay_k.y=veck.y*Rab.y/rab-((1.0/CUBE(v))*(vec_v.y*vec_v.y)-1.0/v)/rab;
  vec_ddtay_k.z=veck.z*Rab.y/rab-(1.0/CUBE(v))*(vec_v.y*vec_v.z)/rab;

  // FINAL PART 1: derivative dtA.z=(Rbc.z-CosTheta*Rab.z)/rab
  vec_ddtaz_i.x=(veci.x*Rab.z+CosTheta*(1.0/CUBE(u))*(vec_u.z*vec_u.x))/rab+CosTheta*Rab.z*(1.0/CUBE(u))*vec_u.x-(1.0/CUBE(u))*vec_u.x*Rbc.z;
  vec_ddtaz_i.y=(veci.y*Rab.z+CosTheta*(1.0/CUBE(u))*(vec_u.z*vec_u.y))/rab+CosTheta*Rab.z*(1.0/CUBE(u))*vec_u.y-(1.0/CUBE(u))*vec_u.y*Rbc.z;
  vec_ddtaz_i.z=(veci.z*Rab.z+CosTheta*((1.0/CUBE(u))*(vec_u.z*vec_u.z)-1.0/u))/rab+CosTheta*Rab.z*(1.0/CUBE(u))*vec_u.z-(1.0/CUBE(u))*vec_u.z*Rbc.z;

  vec_ddtaz_k.x=veck.x*Rab.z/rab-(1.0/CUBE(v))*(vec_v.z*vec_v.x)/rab;
  vec_ddtaz_k.y=veck.y*Rab.z/rab-(1.0/CUBE(v))*(vec_v.z*vec_v.y)/rab;
  vec_ddtaz_k.z=veck.z*Rab.z/rab-((1.0/CUBE(v))*(vec_v.z*vec_v.z)-1.0/v)/rab;

  // FINAL PART 2: derivative dtC.x=(Rab.x-CosTheta*Rbc.x)/rbc
  vec_ddtcx_i.x=veci.x*Rbc.x/rbc-((1.0/CUBE(u))*(vec_u.x*vec_u.x)-1/u)/rbc;
  vec_ddtcx_i.y=veci.y*Rbc.x/rbc-(1.0/CUBE(u))*(vec_u.x*vec_u.y)/rbc;
  vec_ddtcx_i.z=veci.z*Rbc.x/rbc-(1.0/CUBE(u))*(vec_u.x*vec_u.z)/rbc;

  vec_ddtcx_k.x=(veck.x*Rbc.x+CosTheta*((1.0/CUBE(v))*(vec_v.x*vec_v.x)-1/v))/rbc+CosTheta*Rbc.x*(1.0/CUBE(v))*vec_v.x-(1.0/CUBE(v))*vec_v.x*Rab.x;
  vec_ddtcx_k.y=(veck.y*Rbc.x+CosTheta*(1.0/CUBE(v))*(vec_v.x*vec_v.y))/rbc+CosTheta*Rbc.x*(1.0/CUBE(v))*vec_v.y-(1.0/CUBE(v))*vec_v.y*Rab.x;
  vec_ddtcx_k.z=(veck.z*Rbc.x+CosTheta*(1.0/CUBE(v))*(vec_v.x*vec_v.z))/rbc+CosTheta*Rbc.x*(1.0/CUBE(v))*vec_v.z-(1.0/CUBE(v))*vec_v.z*Rab.x;

  // FINAL PART 2: derivative dtC.y=(Rab.y-CosTheta*Rbc.y)/rbc
  vec_ddtcy_i.x=veci.x*Rbc.y/rbc-(1.0/CUBE(u))*(vec_u.y*vec_u.x)/rbc;
  vec_ddtcy_i.y=veci.y*Rbc.y/rbc-((1.0/CUBE(u))*(vec_u.y*vec_u.y)-1.0/u)/rbc;
  vec_ddtcy_i.z=veci.z*Rbc.y/rbc-(1.0/CUBE(u))*(vec_u.y*vec_u.z)/rbc;

  vec_ddtcy_k.x=(veck.x*Rbc.y+CosTheta*(1.0/CUBE(v))*(vec_v.y*vec_v.x))/rbc+CosTheta*Rbc.y*(1.0/CUBE(v))*vec_v.x-(1.0/CUBE(v))*vec_v.x*Rab.y;
  vec_ddtcy_k.y=(veck.y*Rbc.y+CosTheta*((1.0/CUBE(v))*(vec_v.y*vec_v.y)-1.0/v))/rbc+CosTheta*Rbc.y*(1.0/CUBE(v))*vec_v.y-(1.0/CUBE(v))*vec_v.y*Rab.y;
  vec_ddtcy_k.z=(veck.z*Rbc.y+CosTheta*(1.0/CUBE(v))*(vec_v.y*vec_v.z))/rbc+CosTheta*Rbc.y*(1.0/CUBE(v))*vec_v.z-(1.0/CUBE(v))*vec_v.z*Rab.y;

  // FINAL PART 2: derivative dtC.z=(Rab.z-CosTheta*Rbc.z)/rbc
  vec_ddtcz_i.x=veci.x*Rbc.z/rbc-(1.0/CUBE(u))*(vec_u.z*vec_u.x)/rbc;
  vec_ddtcz_i.y=veci.y*Rbc.z/rbc-(1.0/CUBE(u))*(vec_u.z*vec_u.y)/rbc;
  vec_ddtcz_i.z=veci.z*Rbc.z/rbc-((1.0/CUBE(u))*(vec_u.z*vec_u.z)-1.0/u)/rbc;

  vec_ddtcz_k.x=(veck.x*Rbc.z+CosTheta*(1.0/CUBE(v))*(vec_v.z*vec_v.x))/rbc+CosTheta*Rbc.z*(1.0/CUBE(v))*vec_v.x-(1.0/CUBE(v))*vec_v.x*Rab.z;
  vec_ddtcz_k.y=(veck.y*Rbc.z+CosTheta*(1.0/CUBE(v))*(vec_v.z*vec_v.y))/rbc+CosTheta*Rbc.z*(1.0/CUBE(v))*vec_v.y-(1.0/CUBE(v))*vec_v.y*Rab.z;
  vec_ddtcz_k.z=(veck.z*Rbc.z+CosTheta*((1.0/CUBE(v))*(vec_v.z*vec_v.z)-1.0/v))/rbc+CosTheta*Rbc.z*(1.0/CUBE(v))*vec_v.z-(1.0/CUBE(v))*vec_v.z*Rab.z;


  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // contribution is zero (i.e. cancels out)
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          // S.ax=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [0]
          veci.x=vec_ddtax_i.x*DF*vec_u.x+dtA.x*(DDF*dtA.x*vec_u.x+DF) +vec_ddtcx_i.x*DF*vec_v.x+dtC.x*DDF*dtA.x*vec_v.x;
          veci.y=vec_ddtax_i.y*DF*vec_u.x+dtA.x*DDF*dtA.y*vec_u.x      +vec_ddtcx_i.y*DF*vec_v.x+dtC.x*DDF*dtA.y*vec_v.x;
          veci.z=vec_ddtax_i.z*DF*vec_u.x+dtA.x*DDF*dtA.z*vec_u.x      +vec_ddtcx_i.z*DF*vec_v.x+dtC.x*DDF*dtA.z*vec_v.x;

          veck.x=vec_ddtax_k.x*DF*vec_u.x+dtA.x*DDF*dtC.x*vec_u.x      +vec_ddtcx_k.x*DF*vec_v.x+dtC.x*(DDF*dtC.x*vec_v.x+DF);
          veck.y=vec_ddtax_k.y*DF*vec_u.x+dtA.x*DDF*dtC.y*vec_u.x      +vec_ddtcx_k.y*DF*vec_v.x+dtC.x*DDF*dtC.y*vec_v.x;
          veck.z=vec_ddtax_k.z*DF*vec_u.x+dtA.x*DDF*dtC.z*vec_u.x      +vec_ddtcx_k.z*DF*vec_v.x+dtC.x*DDF*dtC.z*vec_v.x;

          HessianMatrix.element[index_i][n]+=veci.x;
          HessianMatrix.element[index_i+1][n]+=veci.y;
          HessianMatrix.element[index_i+2][n]+=veci.z;

          HessianMatrix.element[index_j][n]-=(veci.x+veck.x);
          HessianMatrix.element[index_j+1][n]-=(veci.y+veck.y);
          HessianMatrix.element[index_j+2][n]-=(veci.z+veck.z);

          HessianMatrix.element[index_k][n]+=veck.x;
          HessianMatrix.element[index_k+1][n]+=veck.y;
          HessianMatrix.element[index_k+2][n]+=veck.z;

          veci.x=vec_ddtay_i.x*DF*vec_u.y+dtA.y*DDF*dtA.x*vec_u.y      +vec_ddtcy_i.x*DF*vec_v.y+dtC.y*DDF*dtA.x*vec_v.y;
          veci.y=vec_ddtay_i.y*DF*vec_u.y+dtA.y*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcy_i.y*DF*vec_v.y+dtC.y*DDF*dtA.y*vec_v.y;
          veci.z=vec_ddtay_i.z*DF*vec_u.y+dtA.y*DDF*dtA.z*vec_u.y      +vec_ddtcy_i.z*DF*vec_v.y+dtC.y*DDF*dtA.z*vec_v.y;

          veck.x=vec_ddtay_k.x*DF*vec_u.y+dtA.y*DDF*dtC.x*vec_u.y      +vec_ddtcy_k.x*DF*vec_v.y+dtC.y*DDF*dtC.x*vec_v.y;
          veck.y=vec_ddtay_k.y*DF*vec_u.y+dtA.y*DDF*dtC.y*vec_u.y      +vec_ddtcy_k.y*DF*vec_v.y+dtC.y*(DDF*dtC.y*vec_v.y+DF);
          veck.z=vec_ddtay_k.z*DF*vec_u.y+dtA.y*DDF*dtC.z*vec_u.y      +vec_ddtcy_k.z*DF*vec_v.y+dtC.y*DDF*dtC.z*vec_v.y;

          HessianMatrix.element[index_i][n+1]+=veci.x;
          HessianMatrix.element[index_i+1][n+1]+=veci.y;
          HessianMatrix.element[index_i+2][n+1]+=veci.z;

          HessianMatrix.element[index_j][n+1]-=(veci.x+veck.x);
          HessianMatrix.element[index_j+1][n+1]-=(veci.y+veck.y);
          HessianMatrix.element[index_j+2][n+1]-=(veci.z+veck.z);

          HessianMatrix.element[index_k][n+1]+=veck.x;
          HessianMatrix.element[index_k+1][n+1]+=veck.y;
          HessianMatrix.element[index_k+2][n+1]+=veck.z;

          veci.x=vec_ddtaz_i.x*DF*vec_u.z+dtA.z*DDF*dtA.x*vec_u.z      +vec_ddtcz_i.x*DF*vec_v.z+dtC.z*DDF*dtA.x*vec_v.z;
          veci.y=vec_ddtaz_i.y*DF*vec_u.z+dtA.z*DDF*dtA.y*vec_u.z      +vec_ddtcz_i.y*DF*vec_v.z+dtC.z*DDF*dtA.y*vec_v.z;
          veci.z=vec_ddtaz_i.z*DF*vec_u.z+dtA.z*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcz_i.z*DF*vec_v.z+dtC.z*DDF*dtA.z*vec_v.z;

          veck.x=vec_ddtaz_k.x*DF*vec_u.z+dtA.z*DDF*dtC.x*vec_u.z      +vec_ddtcz_k.x*DF*vec_v.z+dtC.z*DDF*dtC.x*vec_v.z;
          veck.y=vec_ddtaz_k.y*DF*vec_u.z+dtA.z*DDF*dtC.y*vec_u.z      +vec_ddtcz_k.y*DF*vec_v.z+dtC.z*DDF*dtC.y*vec_v.z;
          veck.z=vec_ddtaz_k.z*DF*vec_u.z+dtA.z*DDF*dtC.z*vec_u.z      +vec_ddtcz_k.z*DF*vec_v.z+dtC.z*(DDF*dtC.z*vec_v.z+DF);

          HessianMatrix.element[index_i][n+2]+=veci.x;
          HessianMatrix.element[index_i+1][n+2]+=veci.y;
          HessianMatrix.element[index_i+2][n+2]+=veci.z;

          HessianMatrix.element[index_j][n+2]-=(veci.x+veck.x);
          HessianMatrix.element[index_j+1][n+2]-=(veci.y+veck.y);
          HessianMatrix.element[index_j+2][n+2]-=(veci.z+veck.z);

          HessianMatrix.element[index_k][n+2]+=veck.x;
          HessianMatrix.element[index_k+1][n+2]+=veck.y;
          HessianMatrix.element[index_k+2][n+2]+=veck.z;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          // S.ax=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [0]
          veci.x=vec_ddtax_i.x*DF*vec_u.x+dtA.x*(DDF*dtA.x*vec_u.x+DF) +vec_ddtcx_i.x*DF*vec_v.x+dtC.x*DDF*dtA.x*vec_v.x;
          veci.y=vec_ddtax_i.y*DF*vec_u.x+dtA.x*DDF*dtA.y*vec_u.x      +vec_ddtcx_i.y*DF*vec_v.x+dtC.x*DDF*dtA.y*vec_v.x;
          veci.z=vec_ddtax_i.z*DF*vec_u.x+dtA.x*DDF*dtA.z*vec_u.x      +vec_ddtcx_i.z*DF*vec_v.x+dtC.x*DDF*dtA.z*vec_v.x;

          veck.x=vec_ddtax_k.x*DF*vec_u.x+dtA.x*DDF*dtC.x*vec_u.x      +vec_ddtcx_k.x*DF*vec_v.x+dtC.x*(DDF*dtC.x*vec_v.x+DF);
          veck.y=vec_ddtax_k.y*DF*vec_u.x+dtA.x*DDF*dtC.y*vec_u.x      +vec_ddtcx_k.y*DF*vec_v.x+dtC.x*DDF*dtC.y*vec_v.x;
          veck.z=vec_ddtax_k.z*DF*vec_u.x+dtA.x*DDF*dtC.z*vec_u.x      +vec_ddtcx_k.z*DF*vec_v.x+dtC.x*DDF*dtC.z*vec_v.x;

          HessianMatrix.element[index_i][n]+=veci.x;
          HessianMatrix.element[index_i+1][n]+=veci.y;
          HessianMatrix.element[index_i+2][n]+=veci.z;

          HessianMatrix.element[index_j][n]-=(veci.x+veck.x);
          HessianMatrix.element[index_j+1][n]-=(veci.y+veck.y);
          HessianMatrix.element[index_j+2][n]-=(veci.z+veck.z);

          HessianMatrix.element[index_k][n]+=veck.x;
          HessianMatrix.element[index_k+1][n]+=veck.y;
          HessianMatrix.element[index_k+2][n]+=veck.z;

          // S.ay=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [1] and [3]
          veci.x=vec_ddtax_i.x*DF*vec_u.y+dtA.x*DDF*dtA.x*vec_u.y      +vec_ddtcx_i.x*DF*vec_v.y+dtC.x*DDF*dtA.x*vec_v.y;
          veci.y=vec_ddtax_i.y*DF*vec_u.y+dtA.x*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcx_i.y*DF*vec_v.y+dtC.x*DDF*dtA.y*vec_v.y;
          veci.z=vec_ddtax_i.z*DF*vec_u.y+dtA.x*DDF*dtA.z*vec_u.y      +vec_ddtcx_i.z*DF*vec_v.y+dtC.x*DDF*dtA.z*vec_v.y;

          veck.x=vec_ddtax_k.x*DF*vec_u.y+dtA.x*DDF*dtC.x*vec_u.y      +vec_ddtcx_k.x*DF*vec_v.y+dtC.x*DDF*dtC.x*vec_v.y;
          veck.y=vec_ddtax_k.y*DF*vec_u.y+dtA.x*DDF*dtC.y*vec_u.y      +vec_ddtcx_k.y*DF*vec_v.y+dtC.x*(DDF*dtC.y*vec_v.y+DF);
          veck.z=vec_ddtax_k.z*DF*vec_u.y+dtA.x*DDF*dtC.z*vec_u.y      +vec_ddtcx_k.z*DF*vec_v.y+dtC.x*DDF*dtC.z*vec_v.y;

          HessianMatrix.element[index_i][n+1]+=veci.x;
          HessianMatrix.element[index_i+1][n+1]+=veci.y;
          HessianMatrix.element[index_i+2][n+1]+=veci.z;

          HessianMatrix.element[index_j][n+1]-=(veci.x+veck.x);
          HessianMatrix.element[index_j+1][n+1]-=(veci.y+veck.y);
          HessianMatrix.element[index_j+2][n+1]-=(veci.z+veck.z);

          HessianMatrix.element[index_k][n+1]+=veck.x;
          HessianMatrix.element[index_k+1][n+1]+=veck.y;
          HessianMatrix.element[index_k+2][n+1]+=veck.z;

          // S.cx=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;   [2] and [6]
          veci.x=vec_ddtax_i.x*DF*vec_u.z+dtA.x*DDF*dtA.x*vec_u.z      +vec_ddtcx_i.x*DF*vec_v.z+dtC.x*DDF*dtA.x*vec_v.z;
          veci.y=vec_ddtax_i.y*DF*vec_u.z+dtA.x*DDF*dtA.y*vec_u.z      +vec_ddtcx_i.y*DF*vec_v.z+dtC.x*DDF*dtA.y*vec_v.z;
          veci.z=vec_ddtax_i.z*DF*vec_u.z+dtA.x*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcx_i.z*DF*vec_v.z+dtC.x*DDF*dtA.z*vec_v.z;

          veck.x=vec_ddtax_k.x*DF*vec_u.z+dtA.x*DDF*dtC.x*vec_u.z      +vec_ddtcx_k.x*DF*vec_v.z+dtC.x*DDF*dtC.x*vec_v.z;
          veck.y=vec_ddtax_k.y*DF*vec_u.z+dtA.x*DDF*dtC.y*vec_u.z      +vec_ddtcx_k.y*DF*vec_v.z+dtC.x*DDF*dtC.y*vec_v.z;
          veck.z=vec_ddtax_k.z*DF*vec_u.z+dtA.x*DDF*dtC.z*vec_u.z      +vec_ddtcx_k.z*DF*vec_v.z+dtC.x*(DDF*dtC.z*vec_v.z+DF);

          HessianMatrix.element[index_i][n+2]+=veci.x;
          HessianMatrix.element[index_i+1][n+2]+=veci.y;
          HessianMatrix.element[index_i+2][n+2]+=veci.z;

          HessianMatrix.element[index_j][n+2]-=(veci.x+veck.x);
          HessianMatrix.element[index_j+1][n+2]-=(veci.y+veck.y);
          HessianMatrix.element[index_j+2][n+2]-=(veci.z+veck.z);

          HessianMatrix.element[index_k][n+2]+=veck.x;
          HessianMatrix.element[index_k+1][n+2]+=veck.y;
          HessianMatrix.element[index_k+2][n+2]+=veck.z;

          // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [4]
          veci.x=vec_ddtay_i.x*DF*vec_u.y+dtA.y*DDF*dtA.x*vec_u.y      +vec_ddtcy_i.x*DF*vec_v.y+dtC.y*DDF*dtA.x*vec_v.y;
          veci.y=vec_ddtay_i.y*DF*vec_u.y+dtA.y*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcy_i.y*DF*vec_v.y+dtC.y*DDF*dtA.y*vec_v.y;
          veci.z=vec_ddtay_i.z*DF*vec_u.y+dtA.y*DDF*dtA.z*vec_u.y      +vec_ddtcy_i.z*DF*vec_v.y+dtC.y*DDF*dtA.z*vec_v.y;

          veck.x=vec_ddtay_k.x*DF*vec_u.y+dtA.y*DDF*dtC.x*vec_u.y      +vec_ddtcy_k.x*DF*vec_v.y+dtC.y*DDF*dtC.x*vec_v.y;
          veck.y=vec_ddtay_k.y*DF*vec_u.y+dtA.y*DDF*dtC.y*vec_u.y      +vec_ddtcy_k.y*DF*vec_v.y+dtC.y*(DDF*dtC.y*vec_v.y+DF);
          veck.z=vec_ddtay_k.z*DF*vec_u.y+dtA.y*DDF*dtC.z*vec_u.y      +vec_ddtcy_k.z*DF*vec_v.y+dtC.y*DDF*dtC.z*vec_v.y;

          HessianMatrix.element[index_i][n+3]+=veci.x;
          HessianMatrix.element[index_i+1][n+3]+=veci.y;
          HessianMatrix.element[index_i+2][n+3]+=veci.z;

          HessianMatrix.element[index_j][n+3]-=(veci.x+veck.x);
          HessianMatrix.element[index_j+1][n+3]-=(veci.y+veck.y);
          HessianMatrix.element[index_j+2][n+3]-=(veci.z+veck.z);

          HessianMatrix.element[index_k][n+3]+=veck.x;
          HessianMatrix.element[index_k+1][n+3]+=veck.y;
          HessianMatrix.element[index_k+2][n+3]+=veck.z;

          // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
          veci.x=vec_ddtay_i.x*DF*vec_u.z+dtA.y*DDF*dtA.x*vec_u.z      +vec_ddtcy_i.x*DF*vec_v.z+dtC.y*DDF*dtA.x*vec_v.z;
          veci.y=vec_ddtay_i.y*DF*vec_u.z+dtA.y*DDF*dtA.y*vec_u.z      +vec_ddtcy_i.y*DF*vec_v.z+dtC.y*DDF*dtA.y*vec_v.z;
          veci.z=vec_ddtay_i.z*DF*vec_u.z+dtA.y*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcy_i.z*DF*vec_v.z+dtC.y*DDF*dtA.z*vec_v.z;

          veck.x=vec_ddtay_k.x*DF*vec_u.z+dtA.y*DDF*dtC.x*vec_u.z      +vec_ddtcy_k.x*DF*vec_v.z+dtC.y*DDF*dtC.x*vec_v.z;
          veck.y=vec_ddtay_k.y*DF*vec_u.z+dtA.y*DDF*dtC.y*vec_u.z      +vec_ddtcy_k.y*DF*vec_v.z+dtC.y*DDF*dtC.y*vec_v.z;
          veck.z=vec_ddtay_k.z*DF*vec_u.z+dtA.y*DDF*dtC.z*vec_u.z      +vec_ddtcy_k.z*DF*vec_v.z+dtC.y*(DDF*dtC.z*vec_v.z+DF);

          HessianMatrix.element[index_i][n+4]+=veci.x;
          HessianMatrix.element[index_i+1][n+4]+=veci.y;
          HessianMatrix.element[index_i+2][n+4]+=veci.z;

          HessianMatrix.element[index_j][n+4]-=(veci.x+veck.x);
          HessianMatrix.element[index_j+1][n+4]-=(veci.y+veck.y);
          HessianMatrix.element[index_j+2][n+4]-=(veci.z+veck.z);

          HessianMatrix.element[index_k][n+4]+=veck.x;
          HessianMatrix.element[index_k+1][n+4]+=veck.y;
          HessianMatrix.element[index_k+2][n+4]+=veck.z;

          // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
          veci.x=vec_ddtaz_i.x*DF*vec_u.z+dtA.z*DDF*dtA.x*vec_u.z      +vec_ddtcz_i.x*DF*vec_v.z+dtC.z*DDF*dtA.x*vec_v.z;
          veci.y=vec_ddtaz_i.y*DF*vec_u.z+dtA.z*DDF*dtA.y*vec_u.z      +vec_ddtcz_i.y*DF*vec_v.z+dtC.z*DDF*dtA.y*vec_v.z;
          veci.z=vec_ddtaz_i.z*DF*vec_u.z+dtA.z*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcz_i.z*DF*vec_v.z+dtC.z*DDF*dtA.z*vec_v.z;

          veck.x=vec_ddtaz_k.x*DF*vec_u.z+dtA.z*DDF*dtC.x*vec_u.z      +vec_ddtcz_k.x*DF*vec_v.z+dtC.z*DDF*dtC.x*vec_v.z;
          veck.y=vec_ddtaz_k.y*DF*vec_u.z+dtA.z*DDF*dtC.y*vec_u.z      +vec_ddtcz_k.y*DF*vec_v.z+dtC.z*DDF*dtC.y*vec_v.z;
          veck.z=vec_ddtaz_k.z*DF*vec_u.z+dtA.z*DDF*dtC.z*vec_u.z      +vec_ddtcz_k.z*DF*vec_v.z+dtC.z*(DDF*dtC.z*vec_v.z+DF);

          HessianMatrix.element[index_i][n+5]+=veci.x;
          HessianMatrix.element[index_i+1][n+5]+=veci.y;
          HessianMatrix.element[index_i+2][n+5]+=veci.z;

          HessianMatrix.element[index_j][n+5]-=(veci.x+veck.x);
          HessianMatrix.element[index_j+1][n+5]-=(veci.y+veck.y);
          HessianMatrix.element[index_j+2][n+5]-=(veci.z+veck.z);

          HessianMatrix.element[index_k][n+5]+=veck.x;
          HessianMatrix.element[index_k+1][n+5]+=veck.y;
          HessianMatrix.element[index_k+2][n+5]+=veck.z;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // S.ax=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [0]
              veci.x=vec_ddtax_i.x*DF*vec_u.x+dtA.x*(DDF*dtA.x*vec_u.x+DF) +vec_ddtcx_i.x*DF*vec_v.x+dtC.x*DDF*dtA.x*vec_v.x;
              veci.y=vec_ddtax_i.y*DF*vec_u.x+dtA.x*DDF*dtA.y*vec_u.x      +vec_ddtcx_i.y*DF*vec_v.x+dtC.x*DDF*dtA.y*vec_v.x;
              veci.z=vec_ddtax_i.z*DF*vec_u.x+dtA.x*DDF*dtA.z*vec_u.x      +vec_ddtcx_i.z*DF*vec_v.x+dtC.x*DDF*dtA.z*vec_v.x;

              veck.x=vec_ddtax_k.x*DF*vec_u.x+dtA.x*DDF*dtC.x*vec_u.x      +vec_ddtcx_k.x*DF*vec_v.x+dtC.x*(DDF*dtC.x*vec_v.x+DF);
              veck.y=vec_ddtax_k.y*DF*vec_u.x+dtA.x*DDF*dtC.y*vec_u.x      +vec_ddtcx_k.y*DF*vec_v.x+dtC.x*DDF*dtC.y*vec_v.x;
              veck.z=vec_ddtax_k.z*DF*vec_u.x+dtA.x*DDF*dtC.z*vec_u.x      +vec_ddtcx_k.z*DF*vec_v.x+dtC.x*DDF*dtC.z*vec_v.x;

              HessianMatrix.element[index_i][n]+=veci.x;
              HessianMatrix.element[index_i+1][n]+=veci.y;
              HessianMatrix.element[index_i+2][n]+=veci.z;

              HessianMatrix.element[index_j][n]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n]+=veck.x;
              HessianMatrix.element[index_k+1][n]+=veck.y;
              HessianMatrix.element[index_k+2][n]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [4]
              veci.x=vec_ddtay_i.x*DF*vec_u.y+dtA.y*DDF*dtA.x*vec_u.y      +vec_ddtcy_i.x*DF*vec_v.y+dtC.y*DDF*dtA.x*vec_v.y;
              veci.y=vec_ddtay_i.y*DF*vec_u.y+dtA.y*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcy_i.y*DF*vec_v.y+dtC.y*DDF*dtA.y*vec_v.y;
              veci.z=vec_ddtay_i.z*DF*vec_u.y+dtA.y*DDF*dtA.z*vec_u.y      +vec_ddtcy_i.z*DF*vec_v.y+dtC.y*DDF*dtA.z*vec_v.y;

              veck.x=vec_ddtay_k.x*DF*vec_u.y+dtA.y*DDF*dtC.x*vec_u.y      +vec_ddtcy_k.x*DF*vec_v.y+dtC.y*DDF*dtC.x*vec_v.y;
              veck.y=vec_ddtay_k.y*DF*vec_u.y+dtA.y*DDF*dtC.y*vec_u.y      +vec_ddtcy_k.y*DF*vec_v.y+dtC.y*(DDF*dtC.y*vec_v.y+DF);
              veck.z=vec_ddtay_k.z*DF*vec_u.y+dtA.y*DDF*dtC.z*vec_u.y      +vec_ddtcy_k.z*DF*vec_v.y+dtC.y*DDF*dtC.z*vec_v.y;

              HessianMatrix.element[index_i][n+1]+=veci.x;
              HessianMatrix.element[index_i+1][n+1]+=veci.y;
              HessianMatrix.element[index_i+2][n+1]+=veci.z;

              HessianMatrix.element[index_j][n+1]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n+1]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n+1]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n+1]+=veck.x;
              HessianMatrix.element[index_k+1][n+1]+=veck.y;
              HessianMatrix.element[index_k+2][n+1]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
              veci.x=vec_ddtay_i.x*DF*vec_u.z+dtA.y*DDF*dtA.x*vec_u.z      +vec_ddtcy_i.x*DF*vec_v.z+dtC.y*DDF*dtA.x*vec_v.z;
              veci.y=vec_ddtay_i.y*DF*vec_u.z+dtA.y*DDF*dtA.y*vec_u.z      +vec_ddtcy_i.y*DF*vec_v.z+dtC.y*DDF*dtA.y*vec_v.z;
              veci.z=vec_ddtay_i.z*DF*vec_u.z+dtA.y*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcy_i.z*DF*vec_v.z+dtC.y*DDF*dtA.z*vec_v.z;

              veck.x=vec_ddtay_k.x*DF*vec_u.z+dtA.y*DDF*dtC.x*vec_u.z      +vec_ddtcy_k.x*DF*vec_v.z+dtC.y*DDF*dtC.x*vec_v.z;
              veck.y=vec_ddtay_k.y*DF*vec_u.z+dtA.y*DDF*dtC.y*vec_u.z      +vec_ddtcy_k.y*DF*vec_v.z+dtC.y*DDF*dtC.y*vec_v.z;
              veck.z=vec_ddtay_k.z*DF*vec_u.z+dtA.y*DDF*dtC.z*vec_u.z      +vec_ddtcy_k.z*DF*vec_v.z+dtC.y*(DDF*dtC.z*vec_v.z+DF);

              HessianMatrix.element[index_i][n+2]+=veci.x;
              HessianMatrix.element[index_i+1][n+2]+=veci.y;
              HessianMatrix.element[index_i+2][n+2]+=veci.z;

              HessianMatrix.element[index_j][n+2]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n+2]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n+2]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n+2]+=veck.x;
              HessianMatrix.element[index_k+1][n+2]+=veck.y;
              HessianMatrix.element[index_k+2][n+2]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
              veci.x=vec_ddtaz_i.x*DF*vec_u.z+dtA.z*DDF*dtA.x*vec_u.z      +vec_ddtcz_i.x*DF*vec_v.z+dtC.z*DDF*dtA.x*vec_v.z;
              veci.y=vec_ddtaz_i.y*DF*vec_u.z+dtA.z*DDF*dtA.y*vec_u.z      +vec_ddtcz_i.y*DF*vec_v.z+dtC.z*DDF*dtA.y*vec_v.z;
              veci.z=vec_ddtaz_i.z*DF*vec_u.z+dtA.z*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcz_i.z*DF*vec_v.z+dtC.z*DDF*dtA.z*vec_v.z;

              veck.x=vec_ddtaz_k.x*DF*vec_u.z+dtA.z*DDF*dtC.x*vec_u.z      +vec_ddtcz_k.x*DF*vec_v.z+dtC.z*DDF*dtC.x*vec_v.z;
              veck.y=vec_ddtaz_k.y*DF*vec_u.z+dtA.z*DDF*dtC.y*vec_u.z      +vec_ddtcz_k.y*DF*vec_v.z+dtC.z*DDF*dtC.y*vec_v.z;
              veck.z=vec_ddtaz_k.z*DF*vec_u.z+dtA.z*DDF*dtC.z*vec_u.z      +vec_ddtcz_k.z*DF*vec_v.z+dtC.z*(DDF*dtC.z*vec_v.z+DF);

              HessianMatrix.element[index_i][n+3]+=veci.x;
              HessianMatrix.element[index_i+1][n+3]+=veci.y;
              HessianMatrix.element[index_i+2][n+3]+=veci.z;

              HessianMatrix.element[index_j][n+3]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n+3]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n+3]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n+3]+=veck.x;
              HessianMatrix.element[index_k+1][n+3]+=veck.y;
              HessianMatrix.element[index_k+2][n+3]+=veck.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              // S.ax=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [0]
              veci.x=vec_ddtax_i.x*DF*vec_u.x+dtA.x*(DDF*dtA.x*vec_u.x+DF) +vec_ddtcx_i.x*DF*vec_v.x+dtC.x*DDF*dtA.x*vec_v.x;
              veci.y=vec_ddtax_i.y*DF*vec_u.x+dtA.x*DDF*dtA.y*vec_u.x      +vec_ddtcx_i.y*DF*vec_v.x+dtC.x*DDF*dtA.y*vec_v.x;
              veci.z=vec_ddtax_i.z*DF*vec_u.x+dtA.x*DDF*dtA.z*vec_u.x      +vec_ddtcx_i.z*DF*vec_v.x+dtC.x*DDF*dtA.z*vec_v.x;

              veck.x=vec_ddtax_k.x*DF*vec_u.x+dtA.x*DDF*dtC.x*vec_u.x      +vec_ddtcx_k.x*DF*vec_v.x+dtC.x*(DDF*dtC.x*vec_v.x+DF);
              veck.y=vec_ddtax_k.y*DF*vec_u.x+dtA.x*DDF*dtC.y*vec_u.x      +vec_ddtcx_k.y*DF*vec_v.x+dtC.x*DDF*dtC.y*vec_v.x;
              veck.z=vec_ddtax_k.z*DF*vec_u.x+dtA.x*DDF*dtC.z*vec_u.x      +vec_ddtcx_k.z*DF*vec_v.x+dtC.x*DDF*dtC.z*vec_v.x;

              HessianMatrix.element[index_i][n]+=veci.x;
              HessianMatrix.element[index_i+1][n]+=veci.y;
              HessianMatrix.element[index_i+2][n]+=veci.z;

              HessianMatrix.element[index_j][n]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n]+=veck.x;
              HessianMatrix.element[index_k+1][n]+=veck.y;
              HessianMatrix.element[index_k+2][n]+=veck.z;

              // S.cx=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;   [2] and [6]
              veci.x=vec_ddtax_i.x*DF*vec_u.z+dtA.x*DDF*dtA.x*vec_u.z      +vec_ddtcx_i.x*DF*vec_v.z+dtC.x*DDF*dtA.x*vec_v.z;
              veci.y=vec_ddtax_i.y*DF*vec_u.z+dtA.x*DDF*dtA.y*vec_u.z      +vec_ddtcx_i.y*DF*vec_v.z+dtC.x*DDF*dtA.y*vec_v.z;
              veci.z=vec_ddtax_i.z*DF*vec_u.z+dtA.x*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcx_i.z*DF*vec_v.z+dtC.x*DDF*dtA.z*vec_v.z;

              veck.x=vec_ddtax_k.x*DF*vec_u.z+dtA.x*DDF*dtC.x*vec_u.z      +vec_ddtcx_k.x*DF*vec_v.z+dtC.x*DDF*dtC.x*vec_v.z;
              veck.y=vec_ddtax_k.y*DF*vec_u.z+dtA.x*DDF*dtC.y*vec_u.z      +vec_ddtcx_k.y*DF*vec_v.z+dtC.x*DDF*dtC.y*vec_v.z;
              veck.z=vec_ddtax_k.z*DF*vec_u.z+dtA.x*DDF*dtC.z*vec_u.z      +vec_ddtcx_k.z*DF*vec_v.z+dtC.x*(DDF*dtC.z*vec_v.z+DF);

              HessianMatrix.element[index_i][n+1]+=veci.x;
              HessianMatrix.element[index_i+1][n+1]+=veci.y;
              HessianMatrix.element[index_i+2][n+1]+=veci.z;

              HessianMatrix.element[index_j][n+1]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n+1]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n+1]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n+1]+=veck.x;
              HessianMatrix.element[index_k+1][n+1]+=veck.y;
              HessianMatrix.element[index_k+2][n+1]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [4]
              veci.x=vec_ddtay_i.x*DF*vec_u.y+dtA.y*DDF*dtA.x*vec_u.y      +vec_ddtcy_i.x*DF*vec_v.y+dtC.y*DDF*dtA.x*vec_v.y;
              veci.y=vec_ddtay_i.y*DF*vec_u.y+dtA.y*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcy_i.y*DF*vec_v.y+dtC.y*DDF*dtA.y*vec_v.y;
              veci.z=vec_ddtay_i.z*DF*vec_u.y+dtA.y*DDF*dtA.z*vec_u.y      +vec_ddtcy_i.z*DF*vec_v.y+dtC.y*DDF*dtA.z*vec_v.y;

              veck.x=vec_ddtay_k.x*DF*vec_u.y+dtA.y*DDF*dtC.x*vec_u.y      +vec_ddtcy_k.x*DF*vec_v.y+dtC.y*DDF*dtC.x*vec_v.y;
              veck.y=vec_ddtay_k.y*DF*vec_u.y+dtA.y*DDF*dtC.y*vec_u.y      +vec_ddtcy_k.y*DF*vec_v.y+dtC.y*(DDF*dtC.y*vec_v.y+DF);
              veck.z=vec_ddtay_k.z*DF*vec_u.y+dtA.y*DDF*dtC.z*vec_u.y      +vec_ddtcy_k.z*DF*vec_v.y+dtC.y*DDF*dtC.z*vec_v.y;

              HessianMatrix.element[index_i][n+2]+=veci.x;
              HessianMatrix.element[index_i+1][n+2]+=veci.y;
              HessianMatrix.element[index_i+2][n+2]+=veci.z;

              HessianMatrix.element[index_j][n+2]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n+2]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n+2]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n+2]+=veck.x;
              HessianMatrix.element[index_k+1][n+2]+=veck.y;
              HessianMatrix.element[index_k+2][n+2]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
              veci.x=vec_ddtaz_i.x*DF*vec_u.z+dtA.z*DDF*dtA.x*vec_u.z      +vec_ddtcz_i.x*DF*vec_v.z+dtC.z*DDF*dtA.x*vec_v.z;
              veci.y=vec_ddtaz_i.y*DF*vec_u.z+dtA.z*DDF*dtA.y*vec_u.z      +vec_ddtcz_i.y*DF*vec_v.z+dtC.z*DDF*dtA.y*vec_v.z;
              veci.z=vec_ddtaz_i.z*DF*vec_u.z+dtA.z*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcz_i.z*DF*vec_v.z+dtC.z*DDF*dtA.z*vec_v.z;

              veck.x=vec_ddtaz_k.x*DF*vec_u.z+dtA.z*DDF*dtC.x*vec_u.z      +vec_ddtcz_k.x*DF*vec_v.z+dtC.z*DDF*dtC.x*vec_v.z;
              veck.y=vec_ddtaz_k.y*DF*vec_u.z+dtA.z*DDF*dtC.y*vec_u.z      +vec_ddtcz_k.y*DF*vec_v.z+dtC.z*DDF*dtC.y*vec_v.z;
              veck.z=vec_ddtaz_k.z*DF*vec_u.z+dtA.z*DDF*dtC.z*vec_u.z      +vec_ddtcz_k.z*DF*vec_v.z+dtC.z*(DDF*dtC.z*vec_v.z+DF);

              HessianMatrix.element[index_i][n+3]+=veci.x;
              HessianMatrix.element[index_i+1][n+3]+=veci.y;
              HessianMatrix.element[index_i+2][n+3]+=veci.z;

              HessianMatrix.element[index_j][n+3]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n+3]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n+3]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n+3]+=veck.x;
              HessianMatrix.element[index_k+1][n+3]+=veck.y;
              HessianMatrix.element[index_k+2][n+3]+=veck.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // S.ax=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [0]
              veci.x=vec_ddtax_i.x*DF*vec_u.x+dtA.x*(DDF*dtA.x*vec_u.x+DF) +vec_ddtcx_i.x*DF*vec_v.x+dtC.x*DDF*dtA.x*vec_v.x;
              veci.y=vec_ddtax_i.y*DF*vec_u.x+dtA.x*DDF*dtA.y*vec_u.x      +vec_ddtcx_i.y*DF*vec_v.x+dtC.x*DDF*dtA.y*vec_v.x;
              veci.z=vec_ddtax_i.z*DF*vec_u.x+dtA.x*DDF*dtA.z*vec_u.x      +vec_ddtcx_i.z*DF*vec_v.x+dtC.x*DDF*dtA.z*vec_v.x;

              veck.x=vec_ddtax_k.x*DF*vec_u.x+dtA.x*DDF*dtC.x*vec_u.x      +vec_ddtcx_k.x*DF*vec_v.x+dtC.x*(DDF*dtC.x*vec_v.x+DF);
              veck.y=vec_ddtax_k.y*DF*vec_u.x+dtA.x*DDF*dtC.y*vec_u.x      +vec_ddtcx_k.y*DF*vec_v.x+dtC.x*DDF*dtC.y*vec_v.x;
              veck.z=vec_ddtax_k.z*DF*vec_u.x+dtA.x*DDF*dtC.z*vec_u.x      +vec_ddtcx_k.z*DF*vec_v.x+dtC.x*DDF*dtC.z*vec_v.x;

              HessianMatrix.element[index_i][n]+=veci.x;
              HessianMatrix.element[index_i+1][n]+=veci.y;
              HessianMatrix.element[index_i+2][n]+=veci.z;

              HessianMatrix.element[index_j][n]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n]+=veck.x;
              HessianMatrix.element[index_k+1][n]+=veck.y;
              HessianMatrix.element[index_k+2][n]+=veck.z;

              // S.ay=vec_u.x*DF*dtA.x+vec_v.x*DF*dtC.x;    [1] and [3]
              veci.x=vec_ddtax_i.x*DF*vec_u.y+dtA.x*DDF*dtA.x*vec_u.y      +vec_ddtcx_i.x*DF*vec_v.y+dtC.x*DDF*dtA.x*vec_v.y;
              veci.y=vec_ddtax_i.y*DF*vec_u.y+dtA.x*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcx_i.y*DF*vec_v.y+dtC.x*DDF*dtA.y*vec_v.y;
              veci.z=vec_ddtax_i.z*DF*vec_u.y+dtA.x*DDF*dtA.z*vec_u.y      +vec_ddtcx_i.z*DF*vec_v.y+dtC.x*DDF*dtA.z*vec_v.y;

              veck.x=vec_ddtax_k.x*DF*vec_u.y+dtA.x*DDF*dtC.x*vec_u.y      +vec_ddtcx_k.x*DF*vec_v.y+dtC.x*DDF*dtC.x*vec_v.y;
              veck.y=vec_ddtax_k.y*DF*vec_u.y+dtA.x*DDF*dtC.y*vec_u.y      +vec_ddtcx_k.y*DF*vec_v.y+dtC.x*(DDF*dtC.y*vec_v.y+DF);
              veck.z=vec_ddtax_k.z*DF*vec_u.y+dtA.x*DDF*dtC.z*vec_u.y      +vec_ddtcx_k.z*DF*vec_v.y+dtC.x*DDF*dtC.z*vec_v.y;

              HessianMatrix.element[index_i][n+1]+=veci.x;
              HessianMatrix.element[index_i+1][n+1]+=veci.y;
              HessianMatrix.element[index_i+2][n+1]+=veci.z;

              HessianMatrix.element[index_j][n+1]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n+1]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n+1]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n+1]+=veck.x;
              HessianMatrix.element[index_k+1][n+1]+=veck.y;
              HessianMatrix.element[index_k+2][n+1]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [4]
              veci.x=vec_ddtay_i.x*DF*vec_u.y+dtA.y*DDF*dtA.x*vec_u.y      +vec_ddtcy_i.x*DF*vec_v.y+dtC.y*DDF*dtA.x*vec_v.y;
              veci.y=vec_ddtay_i.y*DF*vec_u.y+dtA.y*(DDF*dtA.y*vec_u.y+DF) +vec_ddtcy_i.y*DF*vec_v.y+dtC.y*DDF*dtA.y*vec_v.y;
              veci.z=vec_ddtay_i.z*DF*vec_u.y+dtA.y*DDF*dtA.z*vec_u.y      +vec_ddtcy_i.z*DF*vec_v.y+dtC.y*DDF*dtA.z*vec_v.y;

              veck.x=vec_ddtay_k.x*DF*vec_u.y+dtA.y*DDF*dtC.x*vec_u.y      +vec_ddtcy_k.x*DF*vec_v.y+dtC.y*DDF*dtC.x*vec_v.y;
              veck.y=vec_ddtay_k.y*DF*vec_u.y+dtA.y*DDF*dtC.y*vec_u.y      +vec_ddtcy_k.y*DF*vec_v.y+dtC.y*(DDF*dtC.y*vec_v.y+DF);
              veck.z=vec_ddtay_k.z*DF*vec_u.y+dtA.y*DDF*dtC.z*vec_u.y      +vec_ddtcy_k.z*DF*vec_v.y+dtC.y*DDF*dtC.z*vec_v.y;

              HessianMatrix.element[index_i][n+2]+=veci.x;
              HessianMatrix.element[index_i+1][n+2]+=veci.y;
              HessianMatrix.element[index_i+2][n+2]+=veci.z;

              HessianMatrix.element[index_j][n+2]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n+2]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n+2]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n+2]+=veck.x;
              HessianMatrix.element[index_k+1][n+2]+=veck.y;
              HessianMatrix.element[index_k+2][n+2]+=veck.z;

              // vec_u.y*DF*dtA.y+rbc*vec_c.y*DF*dtC.y  [5] [7]
              veci.x=vec_ddtaz_i.x*DF*vec_u.z+dtA.z*DDF*dtA.x*vec_u.z      +vec_ddtcz_i.x*DF*vec_v.z+dtC.z*DDF*dtA.x*vec_v.z;
              veci.y=vec_ddtaz_i.y*DF*vec_u.z+dtA.z*DDF*dtA.y*vec_u.z      +vec_ddtcz_i.y*DF*vec_v.z+dtC.z*DDF*dtA.y*vec_v.z;
              veci.z=vec_ddtaz_i.z*DF*vec_u.z+dtA.z*(DDF*dtA.z*vec_u.z+DF) +vec_ddtcz_i.z*DF*vec_v.z+dtC.z*DDF*dtA.z*vec_v.z;

              veck.x=vec_ddtaz_k.x*DF*vec_u.z+dtA.z*DDF*dtC.x*vec_u.z      +vec_ddtcz_k.x*DF*vec_v.z+dtC.z*DDF*dtC.x*vec_v.z;
              veck.y=vec_ddtaz_k.y*DF*vec_u.z+dtA.z*DDF*dtC.y*vec_u.z      +vec_ddtcz_k.y*DF*vec_v.z+dtC.z*DDF*dtC.y*vec_v.z;
              veck.z=vec_ddtaz_k.z*DF*vec_u.z+dtA.z*DDF*dtC.z*vec_u.z      +vec_ddtcz_k.z*DF*vec_v.z+dtC.z*(DDF*dtC.z*vec_v.z+DF);

              HessianMatrix.element[index_i][n+3]+=veci.x;
              HessianMatrix.element[index_i+1][n+3]+=veci.y;
              HessianMatrix.element[index_i+2][n+3]+=veci.z;

              HessianMatrix.element[index_j][n+3]-=(veci.x+veck.x);
              HessianMatrix.element[index_j+1][n+3]-=(veci.y+veck.y);
              HessianMatrix.element[index_j+2][n+3]-=(veci.z+veck.z);

              HessianMatrix.element[index_k][n+3]+=veck.x;
              HessianMatrix.element[index_k+1][n+3]+=veck.y;
              HessianMatrix.element[index_k+2][n+3]+=veck.z;
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

// strain-strain

static inline void HessianBendStrainStrain(REAL_MATRIX HessianMatrix,VECTOR vec_u,VECTOR vec_v,REAL u,REAL v,REAL DF,REAL DDF,REAL_MATRIX3x3 S,REAL CosTheta)
{
  int n;
  REAL_MATRIX3x3 T;

  n=NumberOfCoordinatesMinimizationVariables;

  T.ax=(vec_u.x*vec_v.x+vec_v.x*vec_u.x)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.x*vec_u.x)/SQR(u)+(vec_v.x*vec_v.x)/SQR(v));
  T.ay=(vec_u.x*vec_v.y+vec_v.x*vec_u.y)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.x*vec_u.y)/SQR(u)+(vec_v.x*vec_v.y)/SQR(v));
  T.az=(vec_u.x*vec_v.z+vec_v.x*vec_u.z)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.x*vec_u.z)/SQR(u)+(vec_v.x*vec_v.z)/SQR(v));

  T.bx=(vec_u.y*vec_v.x+vec_v.y*vec_u.x)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.y*vec_u.x)/SQR(u)+(vec_v.y*vec_v.x)/SQR(v));
  T.by=(vec_u.y*vec_v.y+vec_v.y*vec_u.y)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.y*vec_u.y)/SQR(u)+(vec_v.y*vec_v.y)/SQR(v));
  T.bz=(vec_u.y*vec_v.z+vec_v.y*vec_u.z)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.y*vec_u.z)/SQR(u)+(vec_v.y*vec_v.z)/SQR(v));

  T.cx=(vec_u.z*vec_v.x+vec_v.z*vec_u.x)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.z*vec_u.x)/SQR(u)+(vec_v.z*vec_v.x)/SQR(v));
  T.cy=(vec_u.z*vec_v.y+vec_v.z*vec_u.y)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.z*vec_u.y)/SQR(u)+(vec_v.z*vec_v.y)/SQR(v));
  T.cz=(vec_u.z*vec_v.z+vec_v.z*vec_u.z)/(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)-((vec_u.z*vec_u.z)/SQR(u)+(vec_v.z*vec_v.z)/SQR(v));

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                    +T.ax*T.ax)+2.0*S.ax;
      HessianMatrix.element[n][n]+=2.0*(DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                    +T.ax*T.by));
      HessianMatrix.element[n][n]+=2.0*(DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                    -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                    +T.ax*T.cz));
      HessianMatrix.element[n][n]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                    -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                    +T.by*T.by)+2.0*S.by;
      HessianMatrix.element[n][n]+=2.0*(DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                    -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                    +T.by*T.cz));
      HessianMatrix.element[n][n]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                    -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                    2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                    +T.cz*T.cz)+2.0*S.cz;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                        +T.ax*T.ax)+2.0*S.ax;
          HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.ax*T.by);
          HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.cz);

          HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.by*T.by)+2.0*S.by;
          HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.by*T.cz);

          HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.cz*T.cz)+2.0*S.cz;
          break;
        case REGULAR:
         HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                        +T.ax*T.ax)+2.0*S.ax;
          HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.ay+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.y)/SQR(SQR(v)))
                        +T.ax*T.ay)+S.ay;
          HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.az)+S.az;
          HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.ax*T.by);
          HessianMatrix.element[n][n+4]+=DDF*CosTheta*T.ax*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.bz);
          HessianMatrix.element[n][n+5]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.cz);

          HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.ay*CosTheta*T.ay+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.y)/SQR(SQR(v)))
                        +T.ay*T.ay)+0.5*(S.ax+S.by);
          HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.ay*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.az)+0.5*(S.bz);
          HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.ay*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.ay*T.by)+S.ay;
          HessianMatrix.element[n+1][n+4]+=DDF*CosTheta*T.ay*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.bz)+0.5*S.az;
          HessianMatrix.element[n+1][n+5]+=DDF*CosTheta*T.ay*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.cz);

          HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.az*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.az)+0.5*(S.ax+S.cz);
          HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.az*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.az*T.by);
          HessianMatrix.element[n+2][n+4]+=DDF*CosTheta*T.az*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.bz)+0.5*S.ay;
          HessianMatrix.element[n+2][n+5]+=DDF*CosTheta*T.az*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.cz)+S.az;

          HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.by*T.by)+2.0*S.by;
          HessianMatrix.element[n+3][n+4]+=DDF*CosTheta*T.by*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.by*T.bz)+S.bz;
          HessianMatrix.element[n+3][n+5]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.by*T.cz);

          HessianMatrix.element[n+4][n+4]+=DDF*CosTheta*T.bz*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.bz*T.bz)+0.5*(S.by+S.cz);
          HessianMatrix.element[n+4][n+5]+=DDF*CosTheta*T.bz*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.bz*T.cz)+S.bz;

          HessianMatrix.element[n+5][n+5]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.cz*T.cz)+2.0*S.cz;
          break;
        case REGULAR_UPPER_TRIANGLE:
          HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                        +T.ax*T.ax)+2.0*S.ax;
          HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.ay+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.y)/SQR(SQR(v)))
                        +T.ax*T.ay)+S.ay;
          HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.az)+S.az;
          HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.ax*T.by);
          HessianMatrix.element[n][n+4]+=DDF*CosTheta*T.ax*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.bz);
          HessianMatrix.element[n][n+5]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.ax*T.cz);

          HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.ay*CosTheta*T.ay+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.y)/SQR(SQR(v)))
                        +T.ay*T.ay)+S.by;
          HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.ay*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.az)+S.bz;
          HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.ay*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.ay*T.by);
          HessianMatrix.element[n+1][n+4]+=DDF*CosTheta*T.ay*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.bz);
          HessianMatrix.element[n+1][n+5]+=DDF*CosTheta*T.ay*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.ay*T.cz);

          HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.az*CosTheta*T.az+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.x*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.az)+S.cz;
          HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.az*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.az*T.by);
          HessianMatrix.element[n+2][n+4]+=DDF*CosTheta*T.az*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.bz);
          HessianMatrix.element[n+2][n+5]+=DDF*CosTheta*T.az*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.x*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.az*T.cz);

          HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                        +T.by*T.by)+2.0*S.by;
          HessianMatrix.element[n+3][n+4]+=DDF*CosTheta*T.by*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.by*T.bz)+S.bz;
          HessianMatrix.element[n+3][n+5]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.by*T.cz);

          HessianMatrix.element[n+4][n+4]+=DDF*CosTheta*T.bz*CosTheta*T.bz+DF*CosTheta*(
                        -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                        +T.bz*T.bz)+S.cz;
          HessianMatrix.element[n+4][n+5]+=DDF*CosTheta*T.bz*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.y*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.bz*T.cz);

          HessianMatrix.element[n+5][n+5]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                        -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                        2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                        +T.cz*T.cz)+2.0*S.cz;
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.bz);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.by*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.bz)+S.bz;
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.bz*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.bz*T.bz)+0.5*(S.by+S.cz);
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.bz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.bz*T.cz)+S.bz;

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.az+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.az)+S.az;
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.az*CosTheta*T.az+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.x*vec_v.z)/SQR(SQR(v)))
                            +T.az*T.az)+0.5*(S.ax+S.cz);
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.az*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.az*T.by);
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.az*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.az*T.cz)+S.az;

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.ay+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.ay)+S.ay;
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.ay*CosTheta*T.ay+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.y)/SQR(SQR(v)))
                            +T.ay*T.ay)+0.5*(S.ax+S.by);
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.ay*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ay*T.by)+S.ay;
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.ay*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ay*T.cz);

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.bz);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.by*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.bz)+S.bz;
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.bz*CosTheta*T.bz+DF*CosTheta*(
                            -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.z*vec_u.y*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.y*vec_v.z)/SQR(SQR(v)))
                            +T.bz*T.bz)+S.cz;
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.bz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.bz*T.cz);

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.az+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.az)+S.az;
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.az*CosTheta*T.az+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.x*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.x*vec_v.z)/SQR(SQR(v)))
                            +T.az*T.az)+S.cz;
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.az*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.az*T.by);
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.az*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.az*T.cz);

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=DDF*CosTheta*T.ax*CosTheta*T.ax+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.x)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.x)/SQR(SQR(v)))
                            +T.ax*T.ax)+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=DDF*CosTheta*T.ax*CosTheta*T.ay+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.x*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.ay)+S.ay;
              HessianMatrix.element[n][n+2]+=DDF*CosTheta*T.ax*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ax*T.by);
              HessianMatrix.element[n][n+3]+=DDF*CosTheta*T.ax*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.x*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.x*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ax*T.cz);

              HessianMatrix.element[n+1][n+1]+=DDF*CosTheta*T.ay*CosTheta*T.ay+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.x*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.x*vec_v.y)/SQR(SQR(v)))
                            +T.ay*T.ay)+S.by;
              HessianMatrix.element[n+1][n+2]+=DDF*CosTheta*T.ay*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.ay*T.by);
              HessianMatrix.element[n+1][n+3]+=DDF*CosTheta*T.ay*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.x*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.x*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.ay*T.cz);

              HessianMatrix.element[n+2][n+2]+=DDF*CosTheta*T.by*CosTheta*T.by+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.y*vec_u.y)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.y*vec_v.y)/SQR(SQR(v)))
                            +T.by*T.by)+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=DDF*CosTheta*T.by*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.y*vec_u.y*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.y*vec_v.y*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.by*T.cz);

              HessianMatrix.element[n+3][n+3]+=DDF*CosTheta*T.cz*CosTheta*T.cz+DF*CosTheta*(
                            -(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)/SQR(vec_u.x*vec_v.x+vec_u.y*vec_v.y+vec_u.z*vec_v.z)+
                            2.0*((vec_u.z*vec_u.z*vec_u.z*vec_u.z)/SQR(SQR(u))+(vec_v.z*vec_v.z*vec_v.z*vec_v.z)/SQR(SQR(v)))
                            +T.cz*T.cz)+2.0*S.cz;
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}


static inline void GradientStrainBend(REAL *Gradient,REAL_MATRIX3x3 S)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]+=S.ax+S.by+S.cz;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n]+=S.ax;
          Gradient[n+1]+=S.by;
          Gradient[n+2]+=S.cz;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n]+=S.ax;
          Gradient[n+1]+=S.bx;
          Gradient[n+2]+=S.cx;
          Gradient[n+3]+=S.by;
          Gradient[n+4]+=S.cy;
          Gradient[n+5]+=S.cz;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.by;
              Gradient[n+2]+=S.cy;
              Gradient[n+3]+=S.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.cx;
              Gradient[n+2]+=S.by;
              Gradient[n+3]+=S.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.bx;
              Gradient[n+2]+=S.by;
              Gradient[n+3]+=S.cz;
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}



void ComputeFrameworkBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,f1,A,B,C;
  REAL *parms,U;
  REAL CosTheta,Theta,SinTheta,temp,temp2;
  REAL rab,rbc,rac;
  POINT posA,posB,posC;
  VECTOR Rab,Rbc,Rac;
  int index_i,index_j,index_k;
  REAL DTDX,DF,DDF;
  REAL_MATRIX3x3 D2I,D2K,D2IK;
  VECTOR dtA,dtB,dtC;
  REAL_MATRIX3x3 S;
  VECTOR vec_u,vec_v;
  REAL u,v;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBends[f1];i++)
    {
      A=Framework[CurrentSystem].Bends[f1][i].A;
      B=Framework[CurrentSystem].Bends[f1][i].B;
      C=Framework[CurrentSystem].Bends[f1][i].C;
      parms=Framework[CurrentSystem].BendArguments[f1][i];

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;
      index_k=Framework[CurrentSystem].Atoms[f1][C].HessianIndex;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      Rab=ApplyBoundaryCondition(Rab);
      vec_u=Rab;
      rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
      u=rab;
      Rab.x/=rab;
      Rab.y/=rab;
      Rab.z/=rab;

      Rbc.x=posC.x-posB.x;
      Rbc.y=posC.y-posB.y;
      Rbc.z=posC.z-posB.z;
      Rbc=ApplyBoundaryCondition(Rbc);
      vec_v=Rbc;
      rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
      v=rbc;
      Rbc.x/=rbc;
      Rbc.y/=rbc;
      Rbc.z/=rbc;

      Rac.x=posC.x-posA.x;
      Rac.y=posC.y-posA.y;
      Rac.z=posC.z-posA.z;
      Rac=ApplyBoundaryCondition(Rac);
      rac=sqrt(SQR(Rac.x)+SQR(Rac.y)+SQR(Rac.z));
      Rac.x/=rac;
      Rac.y/=rac;
      Rac.z/=rac;

      CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
      Theta=acos(CosTheta);
      SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));
      DTDX=-1.0/SinTheta;


      switch(Framework[CurrentSystem].BendType[f1][i])
      {
        case HARMONIC_BEND:
          // (1/2)p_0*(theta-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          temp=Theta-parms[1];
          U=0.5*parms[0]*SQR(temp);
          DF=parms[0]*temp*DTDX;
          DDF=parms[0]*SQR(DTDX)+parms[0]*temp*CosTheta*CUBE(DTDX);
          break;
        case CORE_SHELL_BEND:
          // (1/2)p_0*(theta-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          temp=Theta-parms[1];
          U=0.5*parms[0]*SQR(temp);
          DF=parms[0]*temp*DTDX;
          DDF=parms[0]*SQR(DTDX)+parms[0]*temp*CosTheta*CUBE(DTDX);
          break;
        case QUARTIC_BEND:
          // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
          // ======================================================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // p_2/k_B [K/rad^3]
          // p_3/k_B [K/rad^4]
          temp=(Theta-parms[1]);
          temp2=SQR(temp);
          U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
          DF=(parms[0]*temp+parms[2]*temp2+parms[3]*temp*temp2)*DTDX;
          DDF=parms[0]*SQR(DTDX)+parms[0]*temp*CosTheta*CUBE(DTDX)
              +2.0*parms[2]*temp*SQR(DTDX)+parms[2]*temp2*CosTheta*CUBE(DTDX)
              +3.0*parms[3]*temp2*SQR(DTDX)+parms[3]*temp*temp2*CosTheta*CUBE(DTDX);
          break;
        case CFF_QUARTIC_BEND:
          // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
          // =====================================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // p_2/k_B [K/rad^3]
          // p_3/k_B [K/rad^4]
          temp=(Theta-parms[1]);
          temp2=SQR(temp);
          U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
          DF=(2.0*parms[0]*temp+3.0*parms[2]*temp2+4.0*parms[3]*temp*temp2)*DTDX;
          DDF=2.0*parms[0]*SQR(DTDX)+2.0*parms[0]*temp*CosTheta*CUBE(DTDX)
              +6.0*parms[2]*temp*SQR(DTDX)+3.0*parms[2]*temp2*CosTheta*CUBE(DTDX)
              +12.0*parms[3]*temp2*SQR(DTDX)+4.0*parms[3]*temp*temp2*CosTheta*CUBE(DTDX);
          break;
        case HARMONIC_COSINE_BEND:
          // (1/2)*p_0*(cos(theta)-cos(p_1))^2
          // ===============================================
          // p_0/k_B [K]
          // p_1     [degrees]
          U=0.5*parms[0]*SQR(CosTheta-parms[1]);
          DF=parms[0]*(CosTheta-parms[1]);
          DDF=parms[0];
          break;
        case COSINE_BEND:
          // p_0*(1+cos(p_1*theta-p_2))
          // ===============================================
          // p_0/k_B [K]
          // p_1     [-]
          // p_2     [degrees]
          temp=parms[1]*Theta-parms[2];
          U=parms[0]*(1.0+cos(temp));
          DF=-(parms[0]*parms[1]*sin(temp))*DTDX;
          DDF=-parms[0]*parms[1]*(parms[1]*cos(temp)+sin(temp)*CosTheta*DTDX)*SQR(DTDX);
          break;
        case MM3_BEND:
        case MM3_IN_PLANE_BEND:
          // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
          // =================================================================================================
          // p_0/k_B [mdyne A/rad^2]
          // p_1     [degrees]
          temp=RAD2DEG*(Theta-parms[1]);
          temp2=SQR(temp);
          U=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
          DF=parms[0]*RAD2DEG*(2.0-(3.0*0.014-(4.0*5.6e-5-(5.0*7.0e-7-6.0*2.2e-8*temp)*temp)*temp)*temp)*temp*DTDX;
          printf("TO BE DONE!\n");
          exit(0);
          break;
        case FIXED_BEND:
          U=DF=DDF=0.0;
          break;
        case MEASURE_BEND:
          U=DF=DDF=0.0;
          break;
        default:
          U=DF=DDF=0.0;
          printf("Undefined Bend potential\n");
          exit(0);
          break;
      }

      // add contribution to the energy
      *Energy+=U;

      if((index_i<0)&&(index_j<0)&&(index_k<0)) continue;

      // Calculate the components of the derivatives.
      dtA.x=(Rbc.x-CosTheta*Rab.x)/rab;
      dtA.y=(Rbc.y-CosTheta*Rab.y)/rab;
      dtA.z=(Rbc.z-CosTheta*Rab.z)/rab;

      dtC.x=(Rab.x-CosTheta*Rbc.x)/rbc;
      dtC.y=(Rab.y-CosTheta*Rbc.y)/rbc;
      dtC.z=(Rab.z-CosTheta*Rbc.z)/rbc;

      dtB.x=-(dtA.x+dtC.x);
      dtB.y=-(dtA.y+dtC.y);
      dtB.z=-(dtA.z+dtC.z);

      S.ax=rab*Rab.x*DF*dtA.x+rbc*Rbc.x*DF*dtC.x;
      S.bx=rab*Rab.y*DF*dtA.x+rbc*Rbc.y*DF*dtC.x;
      S.cx=rab*Rab.z*DF*dtA.x+rbc*Rbc.z*DF*dtC.x;

      S.ay=rab*Rab.x*DF*dtA.y+rbc*Rbc.x*DF*dtC.y;
      S.by=rab*Rab.y*DF*dtA.y+rbc*Rbc.y*DF*dtC.y;
      S.cy=rab*Rab.z*DF*dtA.y+rbc*Rbc.z*DF*dtC.y;

      S.az=rab*Rab.x*DF*dtA.z+rbc*Rbc.x*DF*dtC.z;
      S.bz=rab*Rab.y*DF*dtA.z+rbc*Rbc.y*DF*dtC.z;
      S.cz=rab*Rab.z*DF*dtA.z+rbc*Rbc.z*DF*dtC.z;

      StrainDerivative->ax+=S.ax;
      StrainDerivative->bx+=S.bx;
      StrainDerivative->cx+=S.cx;

      StrainDerivative->ay+=S.ay;
      StrainDerivative->by+=S.by;
      StrainDerivative->cy+=S.cy;

      StrainDerivative->az+=S.az;
      StrainDerivative->bz+=S.bz;
      StrainDerivative->cz+=S.cz;

      // add contribution to the first derivatives
      if(ComputeGradient)
      {
        if(index_i>=0)
        {
          Gradient[index_i]+=DF*dtA.x;
          Gradient[index_i+1]+=DF*dtA.y;
          Gradient[index_i+2]+=DF*dtA.z;
        }

        if(index_j>=0)
        {
          Gradient[index_j]+=DF*dtB.x;
          Gradient[index_j+1]+=DF*dtB.y;
          Gradient[index_j+2]+=DF*dtB.z;
        }

        if(index_k>=0)
        {
          Gradient[index_k]+=DF*dtC.x;
          Gradient[index_k+1]+=DF*dtC.y;
          Gradient[index_k+2]+=DF*dtC.z;
        }

        GradientStrainBend(Gradient,S);
      }

      if(ComputeHessian)
      {
        // Calculate some diagonal Hessian terms for A
        D2I.ax=-DF*(2.0*dtA.x*Rab.x+CosTheta*(1.0-Rab.x*Rab.x)/rab)/rab;
        D2I.by=-DF*(2.0*dtA.y*Rab.y+CosTheta*(1.0-Rab.y*Rab.y)/rab)/rab;
        D2I.cz=-DF*(2.0*dtA.z*Rab.z+CosTheta*(1.0-Rab.z*Rab.z)/rab)/rab;
        D2I.ay=DF*(CosTheta*Rab.x*Rab.y/rab-dtA.x*Rab.y-dtA.y*Rab.x)/rab;
        D2I.az=DF*(CosTheta*Rab.x*Rab.z/rab-dtA.x*Rab.z-dtA.z*Rab.x)/rab;
        D2I.bz=DF*(CosTheta*Rab.y*Rab.z/rab-dtA.y*Rab.z-dtA.z*Rab.y)/rab;

        // Calculate some diagonal Hessian terms for C
        D2K.ax=-DF*(2.0*dtC.x*Rbc.x+CosTheta*(1.0-Rbc.x*Rbc.x)/rbc)/rbc;
        D2K.by=-DF*(2.0*dtC.y*Rbc.y+CosTheta*(1.0-Rbc.y*Rbc.y)/rbc)/rbc;
        D2K.cz=-DF*(2.0*dtC.z*Rbc.z+CosTheta*(1.0-Rbc.z*Rbc.z)/rbc)/rbc;
        D2K.ay=DF*(CosTheta*Rbc.x*Rbc.y/rbc-dtC.x*Rbc.y-dtC.y*Rbc.x)/rbc;
        D2K.az=DF*(CosTheta*Rbc.x*Rbc.z/rbc-dtC.x*Rbc.z-dtC.z*Rbc.x)/rbc;
        D2K.bz=DF*(CosTheta*Rbc.y*Rbc.z/rbc-dtC.y*Rbc.z-dtC.z*Rbc.y)/rbc;

        // Calculate the AC off-diagonal terms.
        D2IK.ax=DF*(CosTheta*Rab.x*Rbc.x-Rab.x*Rab.x-Rbc.x*Rbc.x+1.0)/(rab*rbc);
        D2IK.ay=DF*(CosTheta*Rab.x*Rbc.y-Rab.x*Rab.y-Rbc.x*Rbc.y)/(rab*rbc);
        D2IK.az=DF*(CosTheta*Rab.x*Rbc.z-Rab.x*Rab.z-Rbc.x*Rbc.z)/(rab*rbc);
        D2IK.bx=DF*(CosTheta*Rab.y*Rbc.x-Rab.y*Rab.x-Rbc.y*Rbc.x)/(rab*rbc);
        D2IK.by=DF*(CosTheta*Rab.y*Rbc.y-Rab.y*Rab.y-Rbc.y*Rbc.y+1.0)/(rab*rbc);
        D2IK.bz=DF*(CosTheta*Rab.y*Rbc.z-Rab.y*Rab.z-Rbc.y*Rbc.z)/(rab*rbc);
        D2IK.cx=DF*(CosTheta*Rab.z*Rbc.x-Rab.z*Rab.x-Rbc.z*Rbc.x)/(rab*rbc);
        D2IK.cy=DF*(CosTheta*Rab.z*Rbc.y-Rab.z*Rab.y-Rbc.z*Rbc.y)/(rab*rbc);
        D2IK.cz=DF*(CosTheta*Rab.z*Rbc.z-Rab.z*Rab.z-Rbc.z*Rbc.z+1.0)/(rab*rbc);

        // Calculate AA-block of the Hessian
        if(index_i>=0)
        {
          HessianMatrix.element[index_i][index_i]+=DDF*dtA.x*dtA.x+D2I.ax;
          HessianMatrix.element[index_i][index_i+1]+=DDF*dtA.x*dtA.y+D2I.ay;
          HessianMatrix.element[index_i+1][index_i+1]+=DDF*dtA.y*dtA.y+D2I.by;
          HessianMatrix.element[index_i][index_i+2]+=DDF*dtA.x*dtA.z+D2I.az;
          HessianMatrix.element[index_i+1][index_i+2]+=DDF*dtA.y*dtA.z+D2I.bz;
          HessianMatrix.element[index_i+2][index_i+2]+=DDF*dtA.z*dtA.z+D2I.cz;
        }

        // Calculate BB-block of the Hessian
        if(index_j>=0)
        {
          HessianMatrix.element[index_j][index_j]+=DDF*dtB.x*dtB.x+(D2I.ax+D2K.ax+2.0*D2IK.ax);
          HessianMatrix.element[index_j][index_j+1]+=DDF*dtB.x*dtB.y+(D2I.ay+D2K.ay+D2IK.ay+D2IK.bx);
          HessianMatrix.element[index_j+1][index_j+1]+=DDF*dtB.y*dtB.y+(D2I.by+D2K.by+2.0*D2IK.by);
          HessianMatrix.element[index_j][index_j+2]+=DDF*dtB.x*dtB.z+(D2I.az+D2K.az+D2IK.az+D2IK.cx);
          HessianMatrix.element[index_j+1][index_j+2]+=DDF*dtB.y*dtB.z+(D2I.bz+D2K.bz+D2IK.bz+D2IK.cy);
          HessianMatrix.element[index_j+2][index_j+2]+=DDF*dtB.z*dtB.z+(D2I.cz+D2K.cz+2.0*D2IK.cz);
        }

        // Calculate the CC-block of the Hessian
        if(index_k>=0)
        {
          HessianMatrix.element[index_k][index_k]+=DDF*dtC.x*dtC.x+D2K.ax;
          HessianMatrix.element[index_k][index_k+1]+=DDF*dtC.x*dtC.y+D2K.ay;
          HessianMatrix.element[index_k+1][index_k+1]+=DDF*dtC.y*dtC.y+D2K.by;
          HessianMatrix.element[index_k][index_k+2]+=DDF*dtC.x*dtC.z+D2K.az;
          HessianMatrix.element[index_k+1][index_k+2]+=DDF*dtC.y*dtC.z+D2K.bz;
          HessianMatrix.element[index_k+2][index_k+2]+=DDF*dtC.z*dtC.z+D2K.cz;
        }

        // Calculate the AB-block of the Hessian
        if((index_i>=0)&&(index_j>=0))
        {
          if(index_i<index_j)
          {
            HessianMatrix.element[index_i][index_j]+=DDF*dtA.x*dtB.x-D2I.ax-D2IK.ax;
            HessianMatrix.element[index_i][index_j+1]+=DDF*dtA.x*dtB.y-D2I.ay-D2IK.ay;
            HessianMatrix.element[index_i][index_j+2]+=DDF*dtA.x*dtB.z-D2I.az-D2IK.az;
            HessianMatrix.element[index_i+1][index_j]+=DDF*dtA.y*dtB.x-D2I.ay-D2IK.bx;
            HessianMatrix.element[index_i+1][index_j+1]+=DDF*dtA.y*dtB.y-D2I.by-D2IK.by;
            HessianMatrix.element[index_i+1][index_j+2]+=DDF*dtA.y*dtB.z-D2I.bz-D2IK.bz;
            HessianMatrix.element[index_i+2][index_j]+=DDF*dtA.z*dtB.x-D2I.az-D2IK.cx;
            HessianMatrix.element[index_i+2][index_j+1]+=DDF*dtA.z*dtB.y-D2I.bz-D2IK.cy;
            HessianMatrix.element[index_i+2][index_j+2]+=DDF*dtA.z*dtB.z-D2I.cz-D2IK.cz;
          }
          else
          {
            HessianMatrix.element[index_j][index_i]+=DDF*dtA.x*dtB.x-D2I.ax-D2IK.ax;
            HessianMatrix.element[index_j][index_i+1]+=DDF*dtA.y*dtB.x-D2I.ay-D2IK.bx;
            HessianMatrix.element[index_j][index_i+2]+=DDF*dtA.z*dtB.x-D2I.az-D2IK.cx;
            HessianMatrix.element[index_j+1][index_i]+=DDF*dtA.x*dtB.y-D2I.ay-D2IK.ay;
            HessianMatrix.element[index_j+1][index_i+1]+=DDF*dtA.y*dtB.y-D2I.by-D2IK.by;
            HessianMatrix.element[index_j+1][index_i+2]+=DDF*dtA.z*dtB.y-D2I.bz-D2IK.cy;
            HessianMatrix.element[index_j+2][index_i]+=DDF*dtA.x*dtB.z-D2I.az-D2IK.az;
            HessianMatrix.element[index_j+2][index_i+1]+=DDF*dtA.y*dtB.z-D2I.bz-D2IK.bz;
            HessianMatrix.element[index_j+2][index_i+2]+=DDF*dtA.z*dtB.z-D2I.cz-D2IK.cz;
          }
        }

        // Calculate the AC-block of the Hessian
        if((index_i>=0)&&(index_k>=0))
        {
          if(index_i<index_k)
          {
            HessianMatrix.element[index_i][index_k]+=DDF*dtA.x*dtC.x+D2IK.ax;
            HessianMatrix.element[index_i][index_k+1]+=DDF*dtA.x*dtC.y+D2IK.ay;
            HessianMatrix.element[index_i][index_k+2]+=DDF*dtA.x*dtC.z+D2IK.az;
            HessianMatrix.element[index_i+1][index_k]+=DDF*dtA.y*dtC.x+D2IK.bx;
            HessianMatrix.element[index_i+1][index_k+1]+=DDF*dtA.y*dtC.y+D2IK.by;
            HessianMatrix.element[index_i+1][index_k+2]+=DDF*dtA.y*dtC.z+D2IK.bz;
            HessianMatrix.element[index_i+2][index_k]+=DDF*dtA.z*dtC.x+D2IK.cx;
            HessianMatrix.element[index_i+2][index_k+1]+=DDF*dtA.z*dtC.y+D2IK.cy;
            HessianMatrix.element[index_i+2][index_k+2]+=DDF*dtA.z*dtC.z+D2IK.cz;
          }
          else
          {
            HessianMatrix.element[index_k][index_i]+=DDF*dtA.x*dtC.x+D2IK.ax;
            HessianMatrix.element[index_k][index_i+1]+=DDF*dtA.y*dtC.x+D2IK.bx;
            HessianMatrix.element[index_k][index_i+2]+=DDF*dtA.z*dtC.x+D2IK.cx;
            HessianMatrix.element[index_k+1][index_i]+=DDF*dtA.x*dtC.y+D2IK.ay;
            HessianMatrix.element[index_k+1][index_i+1]+=DDF*dtA.y*dtC.y+D2IK.by;
            HessianMatrix.element[index_k+1][index_i+2]+=DDF*dtA.z*dtC.y+D2IK.cy;
            HessianMatrix.element[index_k+2][index_i]+=DDF*dtA.x*dtC.z+D2IK.az;
            HessianMatrix.element[index_k+2][index_i+1]+=DDF*dtA.y*dtC.z+D2IK.bz;
            HessianMatrix.element[index_k+2][index_i+2]+=DDF*dtA.z*dtC.z+D2IK.cz;
          }
        }

        // Calculate the BC-block of the Hessian
        if((index_j>=0)&&(index_k>=0))
        {
          if(index_k<index_j)
          {
            HessianMatrix.element[index_k][index_j]+=DDF*dtC.x*dtB.x-D2K.ax-D2IK.ax;
            HessianMatrix.element[index_k][index_j+1]+=DDF*dtC.x*dtB.y-D2K.ay-D2IK.bx;
            HessianMatrix.element[index_k][index_j+2]+=DDF*dtC.x*dtB.z-D2K.az-D2IK.cx;
            HessianMatrix.element[index_k+1][index_j]+=DDF*dtC.y*dtB.x-D2K.ay-D2IK.ay;
            HessianMatrix.element[index_k+1][index_j+1]+=DDF*dtC.y*dtB.y-D2K.by-D2IK.by;
            HessianMatrix.element[index_k+1][index_j+2]+=DDF*dtC.y*dtB.z-D2K.bz-D2IK.cy;
            HessianMatrix.element[index_k+2][index_j]+=DDF*dtC.z*dtB.x-D2K.az-D2IK.az;
            HessianMatrix.element[index_k+2][index_j+1]+=DDF*dtC.z*dtB.y-D2K.bz-D2IK.bz;
            HessianMatrix.element[index_k+2][index_j+2]+=DDF*dtC.z*dtB.z-D2K.cz-D2IK.cz;
          }
          else
          {
            HessianMatrix.element[index_j][index_k]+=DDF*dtC.x*dtB.x-D2K.ax-D2IK.ax;
            HessianMatrix.element[index_j][index_k+1]+=DDF*dtC.y*dtB.x-D2K.ay-D2IK.ay;
            HessianMatrix.element[index_j][index_k+2]+=DDF*dtC.z*dtB.x-D2K.az-D2IK.az;
            HessianMatrix.element[index_j+1][index_k]+=DDF*dtC.x*dtB.y-D2K.ay-D2IK.bx;
            HessianMatrix.element[index_j+1][index_k+1]+=DDF*dtC.y*dtB.y-D2K.by-D2IK.by;
            HessianMatrix.element[index_j+1][index_k+2]+=DDF*dtC.z*dtB.y-D2K.bz-D2IK.bz;
            HessianMatrix.element[index_j+2][index_k]+=DDF*dtC.x*dtB.z-D2K.az-D2IK.cx;
            HessianMatrix.element[index_j+2][index_k+1]+=DDF*dtC.y*dtB.z-D2K.bz-D2IK.cy;
            HessianMatrix.element[index_j+2][index_k+2]+=DDF*dtC.z*dtB.z-D2K.cz-D2IK.cz;
          }
        }

        HessianBendStrainPosition(index_i,index_j,index_k,HessianMatrix,vec_u,vec_v,u,v,
                                  rab,rbc,Rab,Rbc,dtA,dtC,DF,DDF,S,CosTheta);

        HessianBendStrainStrain(HessianMatrix,vec_u,vec_v,u,v,DF,DDF,S,CosTheta);
      }
    }
  }
}

static inline void HessianTorsionStrainPosition(REAL_MATRIX HessianMatrix,int index_i,int index_j,int index_k,int index_l,
             VECTOR vec_u,VECTOR vec_v,VECTOR vec_w,REAL u,REAL v,REAL w,VECTOR fa,VECTOR fc,VECTOR fd,
             VECTOR Dab,VECTOR Dcb,VECTOR Ddc,REAL rbc,REAL_MATRIX3x3 D2I,REAL_MATRIX3x3 D2K,REAL_MATRIX3x3 D2L,
             REAL_MATRIX3x3 D2IJ,REAL_MATRIX3x3 D2IK,REAL_MATRIX3x3 D2IL,REAL_MATRIX3x3 D2JK,REAL_MATRIX3x3 D2JL,REAL_MATRIX3x3 D2KL,
             VECTOR dtA,VECTOR dtB,VECTOR dtC,VECTOR dtD,REAL DDF,REAL_MATRIX3x3 S,REAL CosPhi)
{
  int n;
  VECTOR veci,vecj,veck,vecl;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // contribution is zero (i.e. cancels out)
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          veci.x=Dab.x*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.x*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.x*(DDF*dtD.x*dtA.x+D2IL.ax)+fa.x;
          veci.y=Dab.x*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.x*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.x*(DDF*dtD.x*dtA.y+D2IL.bx);
          veci.z=Dab.x*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.x*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.x*(DDF*dtD.x*dtA.z+D2IL.cx);

          vecj.x=Dab.x*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.x*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.x*(DDF*dtD.x*dtB.x+D2JL.ax)-fa.x-fc.x-fd.x;
          vecj.y=Dab.x*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.x*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.x*(DDF*dtD.x*dtB.y+D2JL.bx);
          vecj.z=Dab.x*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.x*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.x*(DDF*dtD.x*dtB.z+D2JL.cx);

          veck.x=Dab.x*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.x*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.x*(DDF*dtD.x*dtC.x+D2KL.ax)+fc.x+fd.x-fd.x;
          veck.y=Dab.x*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.x*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.x*(DDF*dtD.x*dtC.y+D2KL.bx);
          veck.z=Dab.x*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.x*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.x*(DDF*dtD.x*dtC.z+D2KL.cx);

          vecl.x=Dab.x*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.x*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.x*(DDF*dtD.x*dtD.x+D2L.ax)+fd.x;
          vecl.y=Dab.x*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.x*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.x*(DDF*dtD.x*dtD.y+D2L.ay);
          vecl.z=Dab.x*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.x*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.x*(DDF*dtD.x*dtD.z+D2L.az);

          HessianMatrix.element[index_i][n]+=veci.x;
          HessianMatrix.element[index_i+1][n]+=veci.y;
          HessianMatrix.element[index_i+2][n]+=veci.z;
          HessianMatrix.element[index_j][n]+=vecj.x;
          HessianMatrix.element[index_j+1][n]+=vecj.y;
          HessianMatrix.element[index_j+2][n]+=vecj.z;
          HessianMatrix.element[index_k][n]+=veck.x;
          HessianMatrix.element[index_k+1][n]+=veck.y;
          HessianMatrix.element[index_k+2][n]+=veck.z;
          HessianMatrix.element[index_l][n]+=vecl.x;
          HessianMatrix.element[index_l+1][n]+=vecl.y;
          HessianMatrix.element[index_l+2][n]+=vecl.z;

          veci.x=Dab.y*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.y*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.y*(DDF*dtD.y*dtA.x+D2IL.ay);
          veci.y=Dab.y*(DDF*dtA.y*dtA.y+D2I.by)+fa.y+Dcb.y*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.y*(DDF*dtD.y*dtA.y+D2IL.by);
          veci.z=Dab.y*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.y*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.y*(DDF*dtD.y*dtA.z+D2IL.cy);

          vecj.x=Dab.y*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.y*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.y*(DDF*dtD.y*dtB.x+D2JL.ay);
          vecj.y=Dab.y*(DDF*dtA.y*dtB.y+D2IJ.by)-fa.y+Dcb.y*rbc*(DDF*dtC.y*dtB.y+D2JK.by)-fc.y+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.y+D2JL.by)-fd.y+Ddc.y*(DDF*dtD.y*dtB.y+D2JL.by);
          vecj.z=Dab.y*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.y*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.y*(DDF*dtD.y*dtB.z+D2JL.cy);

          veck.x=Dab.y*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.y*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.y*(DDF*dtD.y*dtC.x+D2KL.ay);
          veck.y=Dab.y*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.y*rbc*(DDF*dtC.y*dtC.y+D2K.by)+fc.y+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+fd.y+Ddc.y*(DDF*dtD.y*dtC.y+D2KL.by)-fd.y;
          veck.z=Dab.y*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.y*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.y*(DDF*dtD.y*dtC.z+D2KL.cy);

          vecl.x=Dab.y*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.y*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.y*(DDF*dtD.y*dtD.x+D2L.ay);
          vecl.y=Dab.y*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.y*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.y*(DDF*dtD.y*dtD.y+D2L.by)+fd.y;
          vecl.z=Dab.y*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.y*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.y*(DDF*dtD.y*dtD.z+D2L.bz);

          HessianMatrix.element[index_i][n+1]+=veci.x;
          HessianMatrix.element[index_i+1][n+1]+=veci.y;
          HessianMatrix.element[index_i+2][n+1]+=veci.z;
          HessianMatrix.element[index_j][n+1]+=vecj.x;
          HessianMatrix.element[index_j+1][n+1]+=vecj.y;
          HessianMatrix.element[index_j+2][n+1]+=vecj.z;
          HessianMatrix.element[index_k][n+1]+=veck.x;
          HessianMatrix.element[index_k+1][n+1]+=veck.y;
          HessianMatrix.element[index_k+2][n+1]+=veck.z;
          HessianMatrix.element[index_l][n+1]+=vecl.x;
          HessianMatrix.element[index_l+1][n+1]+=vecl.y;
          HessianMatrix.element[index_l+2][n+1]+=vecl.z;

          veci.x=Dab.z*(DDF*dtA.z*dtA.x+D2I.az)+Dcb.z*rbc*(DDF*dtC.z*dtA.x+D2IK.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.x+D2IL.az)+Ddc.z*(DDF*dtD.z*dtA.x+D2IL.az);
          veci.y=Dab.z*(DDF*dtA.z*dtA.y+D2I.bz)+Dcb.z*rbc*(DDF*dtC.z*dtA.y+D2IK.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.y+D2IL.bz)+Ddc.z*(DDF*dtD.z*dtA.y+D2IL.bz);
          veci.z=Dab.z*(DDF*dtA.z*dtA.z+D2I.cz)+fa.z+Dcb.z*rbc*(DDF*dtC.z*dtA.z+D2IK.cz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.z+D2IL.cz)+Ddc.z*(DDF*dtD.z*dtA.z+D2IL.cz);

          vecj.x=Dab.z*(DDF*dtA.z*dtB.x+D2IJ.cx)+Dcb.z*rbc*(DDF*dtC.z*dtB.x+D2JK.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.x+D2JL.az)+Ddc.z*(DDF*dtD.z*dtB.x+D2JL.az);
          vecj.y=Dab.z*(DDF*dtA.z*dtB.y+D2IJ.cy)+Dcb.z*rbc*(DDF*dtC.z*dtB.y+D2JK.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.y+D2JL.bz)+Ddc.z*(DDF*dtD.z*dtB.y+D2JL.bz);
          vecj.z=Dab.z*(DDF*dtA.z*dtB.z+D2IJ.cz)-fa.z+Dcb.z*rbc*(DDF*dtC.z*dtB.z+D2JK.cz)-fc.z+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.z+D2JL.cz)-fd.z+Ddc.z*(DDF*dtD.z*dtB.z+D2JL.cz);

          veck.x=Dab.z*(DDF*dtA.z*dtC.x+D2IK.cx)+Dcb.z*rbc*(DDF*dtC.z*dtC.x+D2K.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.x+D2KL.az)+Ddc.z*(DDF*dtD.z*dtC.x+D2KL.az);
          veck.y=Dab.z*(DDF*dtA.z*dtC.y+D2IK.cy)+Dcb.z*rbc*(DDF*dtC.z*dtC.y+D2K.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.y+D2KL.bz)+Ddc.z*(DDF*dtD.z*dtC.y+D2KL.bz);
          veck.z=Dab.z*(DDF*dtA.z*dtC.z+D2IK.cz)+Dcb.z*rbc*(DDF*dtC.z*dtC.z+D2K.cz)+fc.z+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.z+D2KL.cz)+fd.z+Ddc.z*(DDF*dtD.z*dtC.z+D2KL.cz)-fd.z;

          vecl.x=Dab.z*(DDF*dtA.z*dtD.x+D2IL.cx)+Dcb.z*rbc*(DDF*dtC.z*dtD.x+D2KL.cx)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.x+D2L.az)+Ddc.z*(DDF*dtD.z*dtD.x+D2L.az);
          vecl.y=Dab.z*(DDF*dtA.z*dtD.y+D2IL.cy)+Dcb.z*rbc*(DDF*dtC.z*dtD.y+D2KL.cy)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.y+D2L.bz)+Ddc.z*(DDF*dtD.z*dtD.y+D2L.bz);
          vecl.z=Dab.z*(DDF*dtA.z*dtD.z+D2IL.cz)+Dcb.z*rbc*(DDF*dtC.z*dtD.z+D2KL.cz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.z+D2L.cz)+Ddc.z*(DDF*dtD.z*dtD.z+D2L.cz)+fd.z;

          HessianMatrix.element[index_i][n+2]+=veci.x;
          HessianMatrix.element[index_i+1][n+2]+=veci.y;
          HessianMatrix.element[index_i+2][n+2]+=veci.z;
          HessianMatrix.element[index_j][n+2]+=vecj.x;
          HessianMatrix.element[index_j+1][n+2]+=vecj.y;
          HessianMatrix.element[index_j+2][n+2]+=vecj.z;
          HessianMatrix.element[index_k][n+2]+=veck.x;
          HessianMatrix.element[index_k+1][n+2]+=veck.y;
          HessianMatrix.element[index_k+2][n+2]+=veck.z;
          HessianMatrix.element[index_l][n+2]+=vecl.x;
          HessianMatrix.element[index_l+1][n+2]+=vecl.y;
          HessianMatrix.element[index_l+2][n+2]+=vecl.z;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          //S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
          veci.x=Dab.x*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.x*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.x*(DDF*dtD.x*dtA.x+D2IL.ax)+fa.x;
          veci.y=Dab.x*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.x*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.x*(DDF*dtD.x*dtA.y+D2IL.bx);
          veci.z=Dab.x*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.x*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.x*(DDF*dtD.x*dtA.z+D2IL.cx);

          vecj.x=Dab.x*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.x*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.x*(DDF*dtD.x*dtB.x+D2JL.ax)-fa.x-fc.x-fd.x;
          vecj.y=Dab.x*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.x*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.x*(DDF*dtD.x*dtB.y+D2JL.bx);
          vecj.z=Dab.x*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.x*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                 Dcb.x*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.x*(DDF*dtD.x*dtB.z+D2JL.cx);

          veck.x=Dab.x*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.x*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.x*(DDF*dtD.x*dtC.x+D2KL.ax)+fc.x+fd.x-fd.x;
          veck.y=Dab.x*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.x*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.x*(DDF*dtD.x*dtC.y+D2KL.bx);
          veck.z=Dab.x*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.x*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                 Dcb.x*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.x*(DDF*dtD.x*dtC.z+D2KL.cx);

          vecl.x=Dab.x*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.x*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.x*(DDF*dtD.x*dtD.x+D2L.ax)+fd.x;
          vecl.y=Dab.x*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.x*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.x*(DDF*dtD.x*dtD.y+D2L.ay);
          vecl.z=Dab.x*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.x*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                 Dcb.x*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.x*(DDF*dtD.x*dtD.z+D2L.az);

          HessianMatrix.element[index_i][n]+=veci.x;
          HessianMatrix.element[index_i+1][n]+=veci.y;
          HessianMatrix.element[index_i+2][n]+=veci.z;
          HessianMatrix.element[index_j][n]+=vecj.x;
          HessianMatrix.element[index_j+1][n]+=vecj.y;
          HessianMatrix.element[index_j+2][n]+=vecj.z;
          HessianMatrix.element[index_k][n]+=veck.x;
          HessianMatrix.element[index_k+1][n]+=veck.y;
          HessianMatrix.element[index_k+2][n]+=veck.z;
          HessianMatrix.element[index_l][n]+=vecl.x;
          HessianMatrix.element[index_l+1][n]+=vecl.y;
          HessianMatrix.element[index_l+2][n]+=vecl.z;

          //S.bx=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
          veci.x=Dab.y*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.y*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                 Dcb.y*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.y*(DDF*dtD.x*dtA.x+D2IL.ax);
          veci.y=Dab.y*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.y*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                 Dcb.y*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.y*(DDF*dtD.x*dtA.y+D2IL.bx)+fa.x;
          veci.z=Dab.y*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.y*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                 Dcb.y*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.y*(DDF*dtD.x*dtA.z+D2IL.cx);

          vecj.x=Dab.y*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.y*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                 Dcb.y*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.y*(DDF*dtD.x*dtB.x+D2JL.ax);
          vecj.y=Dab.y*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.y*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                 Dcb.y*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.y*(DDF*dtD.x*dtB.y+D2JL.bx)-fa.x-fc.x-fd.x;
          vecj.z=Dab.y*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.y*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                 Dcb.y*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.y*(DDF*dtD.x*dtB.z+D2JL.cx);

          veck.x=Dab.y*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.y*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                 Dcb.y*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.y*(DDF*dtD.x*dtC.x+D2KL.ax);
          veck.y=Dab.y*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.y*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                 Dcb.y*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.y*(DDF*dtD.x*dtC.y+D2KL.bx)+fc.x+fd.x-fd.x;
          veck.z=Dab.y*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.y*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                 Dcb.y*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.y*(DDF*dtD.x*dtC.z+D2KL.cx);

          vecl.x=Dab.y*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.y*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                 Dcb.y*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.y*(DDF*dtD.x*dtD.x+D2L.ax);
          vecl.y=Dab.y*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.y*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                 Dcb.y*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.y*(DDF*dtD.x*dtD.y+D2L.ay)+fd.x;
          vecl.z=Dab.y*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.y*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                 Dcb.y*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.y*(DDF*dtD.x*dtD.z+D2L.az);

          HessianMatrix.element[index_i][n+1]+=veci.x;
          HessianMatrix.element[index_i+1][n+1]+=veci.y;
          HessianMatrix.element[index_i+2][n+1]+=veci.z;
          HessianMatrix.element[index_j][n+1]+=vecj.x;
          HessianMatrix.element[index_j+1][n+1]+=vecj.y;
          HessianMatrix.element[index_j+2][n+1]+=vecj.z;
          HessianMatrix.element[index_k][n+1]+=veck.x;
          HessianMatrix.element[index_k+1][n+1]+=veck.y;
          HessianMatrix.element[index_k+2][n+1]+=veck.z;
          HessianMatrix.element[index_l][n+1]+=vecl.x;
          HessianMatrix.element[index_l+1][n+1]+=vecl.y;
          HessianMatrix.element[index_l+2][n+1]+=vecl.z;


          veci.x=Dab.z*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.z*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                Dcb.z*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.z*(DDF*dtD.x*dtA.x+D2IL.ax);
          veci.y=Dab.z*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.z*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                Dcb.z*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.z*(DDF*dtD.x*dtA.y+D2IL.bx);
          veci.z=Dab.z*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.z*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                 Dcb.z*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.z*(DDF*dtD.x*dtA.z+D2IL.cx)+fa.x;

          vecj.x=Dab.z*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.z*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                 Dcb.z*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.z*(DDF*dtD.x*dtB.x+D2JL.ax);
          vecj.y=Dab.z*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.z*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                 Dcb.z*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.z*(DDF*dtD.x*dtB.y+D2JL.bx);
          vecj.z=Dab.z*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.z*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                 Dcb.z*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.z*(DDF*dtD.x*dtB.z+D2JL.cx)-fa.x-fc.x-fd.x;

          veck.x=Dab.z*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.z*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                 Dcb.z*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.z*(DDF*dtD.x*dtC.x+D2KL.ax);
          veck.y=Dab.z*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.z*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                 Dcb.z*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.z*(DDF*dtD.x*dtC.y+D2KL.bx);
          veck.z=Dab.z*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.z*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                 Dcb.z*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.z*(DDF*dtD.x*dtC.z+D2KL.cx)+fc.x+fd.x-fd.x;

          vecl.x=Dab.z*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.z*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                 Dcb.z*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.z*(DDF*dtD.x*dtD.x+D2L.ax);
          vecl.y=Dab.z*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.z*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                 Dcb.z*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.z*(DDF*dtD.x*dtD.y+D2L.ay);
          vecl.z=Dab.z*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.z*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                 Dcb.z*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.z*(DDF*dtD.x*dtD.z+D2L.az)+fd.x;

          HessianMatrix.element[index_i][n+2]+=veci.x;
          HessianMatrix.element[index_i+1][n+2]+=veci.y;
          HessianMatrix.element[index_i+2][n+2]+=veci.z;
          HessianMatrix.element[index_j][n+2]+=vecj.x;
          HessianMatrix.element[index_j+1][n+2]+=vecj.y;
          HessianMatrix.element[index_j+2][n+2]+=vecj.z;
          HessianMatrix.element[index_k][n+2]+=veck.x;
          HessianMatrix.element[index_k+1][n+2]+=veck.y;
          HessianMatrix.element[index_k+2][n+2]+=veck.z;
          HessianMatrix.element[index_l][n+2]+=vecl.x;
          HessianMatrix.element[index_l+1][n+2]+=vecl.y;
          HessianMatrix.element[index_l+2][n+2]+=vecl.z;


          veci.x=Dab.y*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.y*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.y*(DDF*dtD.y*dtA.x+D2IL.ay);
          veci.y=Dab.y*(DDF*dtA.y*dtA.y+D2I.by)+fa.y+Dcb.y*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.y*(DDF*dtD.y*dtA.y+D2IL.by);
          veci.z=Dab.y*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.y*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                 Dcb.y*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.y*(DDF*dtD.y*dtA.z+D2IL.cy);

          vecj.x=Dab.y*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.y*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.y*(DDF*dtD.y*dtB.x+D2JL.ay);
          vecj.y=Dab.y*(DDF*dtA.y*dtB.y+D2IJ.by)-fa.y+Dcb.y*rbc*(DDF*dtC.y*dtB.y+D2JK.by)-fc.y+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.y+D2JL.by)-fd.y+Ddc.y*(DDF*dtD.y*dtB.y+D2JL.by);
          vecj.z=Dab.y*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.y*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                 Dcb.y*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.y*(DDF*dtD.y*dtB.z+D2JL.cy);

          veck.x=Dab.y*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.y*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.y*(DDF*dtD.y*dtC.x+D2KL.ay);
          veck.y=Dab.y*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.y*rbc*(DDF*dtC.y*dtC.y+D2K.by)+fc.y+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+fd.y+Ddc.y*(DDF*dtD.y*dtC.y+D2KL.by)-fd.y;
          veck.z=Dab.y*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.y*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                 Dcb.y*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.y*(DDF*dtD.y*dtC.z+D2KL.cy);

          vecl.x=Dab.y*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.y*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.y*(DDF*dtD.y*dtD.x+D2L.ay);
          vecl.y=Dab.y*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.y*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.y*(DDF*dtD.y*dtD.y+D2L.by)+fd.y;
          vecl.z=Dab.y*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.y*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                 Dcb.y*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.y*(DDF*dtD.y*dtD.z+D2L.bz);

          HessianMatrix.element[index_i][n+3]+=veci.x;
          HessianMatrix.element[index_i+1][n+3]+=veci.y;
          HessianMatrix.element[index_i+2][n+3]+=veci.z;
          HessianMatrix.element[index_j][n+3]+=vecj.x;
          HessianMatrix.element[index_j+1][n+3]+=vecj.y;
          HessianMatrix.element[index_j+2][n+3]+=vecj.z;
          HessianMatrix.element[index_k][n+3]+=veck.x;
          HessianMatrix.element[index_k+1][n+3]+=veck.y;
          HessianMatrix.element[index_k+2][n+3]+=veck.z;
          HessianMatrix.element[index_l][n+3]+=vecl.x;
          HessianMatrix.element[index_l+1][n+3]+=vecl.y;
          HessianMatrix.element[index_l+2][n+3]+=vecl.z;


          veci.x=Dab.z*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.z*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                 Dcb.z*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.z*(DDF*dtD.y*dtA.x+D2IL.ay);
          veci.y=Dab.z*(DDF*dtA.y*dtA.y+D2I.by)+Dcb.z*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                 Dcb.z*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.z*(DDF*dtD.y*dtA.y+D2IL.by);
          veci.z=Dab.z*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.z*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                 Dcb.z*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.z*(DDF*dtD.y*dtA.z+D2IL.cy)+fa.y;

          vecj.x=Dab.z*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.z*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                 Dcb.z*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.z*(DDF*dtD.y*dtB.x+D2JL.ay);
          vecj.y=Dab.z*(DDF*dtA.y*dtB.y+D2IJ.by)+Dcb.z*rbc*(DDF*dtC.y*dtB.y+D2JK.by)+
                 Dcb.z*rbc*(DDF*dtD.y*dtB.y+D2JL.by)+Ddc.z*(DDF*dtD.y*dtB.y+D2JL.by);
          vecj.z=Dab.z*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.z*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                 Dcb.z*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.z*(DDF*dtD.y*dtB.z+D2JL.cy)-fa.y-fc.y-fd.y;

          veck.x=Dab.z*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.z*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                 Dcb.z*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.z*(DDF*dtD.y*dtC.x+D2KL.ay);
          veck.y=Dab.z*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.z*rbc*(DDF*dtC.y*dtC.y+D2K.by)+
                 Dcb.z*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+Ddc.z*(DDF*dtD.y*dtC.y+D2KL.by);
          veck.z=Dab.z*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.z*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                 Dcb.z*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.z*(DDF*dtD.y*dtC.z+D2KL.cy)+fc.y+fd.y-fd.y;

          vecl.x=Dab.z*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.z*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                 Dcb.z*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.z*(DDF*dtD.y*dtD.x+D2L.ay);
          vecl.y=Dab.z*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.z*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                 Dcb.z*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.z*(DDF*dtD.y*dtD.y+D2L.by);
          vecl.z=Dab.z*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.z*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                 Dcb.z*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.z*(DDF*dtD.y*dtD.z+D2L.bz)+fd.y;

          HessianMatrix.element[index_i][n+4]+=veci.x;
          HessianMatrix.element[index_i+1][n+4]+=veci.y;
          HessianMatrix.element[index_i+2][n+4]+=veci.z;
          HessianMatrix.element[index_j][n+4]+=vecj.x;
          HessianMatrix.element[index_j+1][n+4]+=vecj.y;
          HessianMatrix.element[index_j+2][n+4]+=vecj.z;
          HessianMatrix.element[index_k][n+4]+=veck.x;
          HessianMatrix.element[index_k+1][n+4]+=veck.y;
          HessianMatrix.element[index_k+2][n+4]+=veck.z;
          HessianMatrix.element[index_l][n+4]+=vecl.x;
          HessianMatrix.element[index_l+1][n+4]+=vecl.y;
          HessianMatrix.element[index_l+2][n+4]+=vecl.z;

          veci.x=Dab.z*(DDF*dtA.z*dtA.x+D2I.az)+Dcb.z*rbc*(DDF*dtC.z*dtA.x+D2IK.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.x+D2IL.az)+Ddc.z*(DDF*dtD.z*dtA.x+D2IL.az);
          veci.y=Dab.z*(DDF*dtA.z*dtA.y+D2I.bz)+Dcb.z*rbc*(DDF*dtC.z*dtA.y+D2IK.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.y+D2IL.bz)+Ddc.z*(DDF*dtD.z*dtA.y+D2IL.bz);
          veci.z=Dab.z*(DDF*dtA.z*dtA.z+D2I.cz)+fa.z+Dcb.z*rbc*(DDF*dtC.z*dtA.z+D2IK.cz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtA.z+D2IL.cz)+Ddc.z*(DDF*dtD.z*dtA.z+D2IL.cz);

          vecj.x=Dab.z*(DDF*dtA.z*dtB.x+D2IJ.cx)+Dcb.z*rbc*(DDF*dtC.z*dtB.x+D2JK.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.x+D2JL.az)+Ddc.z*(DDF*dtD.z*dtB.x+D2JL.az);
          vecj.y=Dab.z*(DDF*dtA.z*dtB.y+D2IJ.cy)+Dcb.z*rbc*(DDF*dtC.z*dtB.y+D2JK.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.y+D2JL.bz)+Ddc.z*(DDF*dtD.z*dtB.y+D2JL.bz);
          vecj.z=Dab.z*(DDF*dtA.z*dtB.z+D2IJ.cz)-fa.z+Dcb.z*rbc*(DDF*dtC.z*dtB.z+D2JK.cz)-fc.z+
                 Dcb.z*rbc*(DDF*dtD.z*dtB.z+D2JL.cz)-fd.z+Ddc.z*(DDF*dtD.z*dtB.z+D2JL.cz);

          veck.x=Dab.z*(DDF*dtA.z*dtC.x+D2IK.cx)+Dcb.z*rbc*(DDF*dtC.z*dtC.x+D2K.az)+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.x+D2KL.az)+Ddc.z*(DDF*dtD.z*dtC.x+D2KL.az);
          veck.y=Dab.z*(DDF*dtA.z*dtC.y+D2IK.cy)+Dcb.z*rbc*(DDF*dtC.z*dtC.y+D2K.bz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.y+D2KL.bz)+Ddc.z*(DDF*dtD.z*dtC.y+D2KL.bz);
          veck.z=Dab.z*(DDF*dtA.z*dtC.z+D2IK.cz)+Dcb.z*rbc*(DDF*dtC.z*dtC.z+D2K.cz)+fc.z+
                 Dcb.z*rbc*(DDF*dtD.z*dtC.z+D2KL.cz)+fd.z+Ddc.z*(DDF*dtD.z*dtC.z+D2KL.cz)-fd.z;

          vecl.x=Dab.z*(DDF*dtA.z*dtD.x+D2IL.cx)+Dcb.z*rbc*(DDF*dtC.z*dtD.x+D2KL.cx)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.x+D2L.az)+Ddc.z*(DDF*dtD.z*dtD.x+D2L.az);
          vecl.y=Dab.z*(DDF*dtA.z*dtD.y+D2IL.cy)+Dcb.z*rbc*(DDF*dtC.z*dtD.y+D2KL.cy)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.y+D2L.bz)+Ddc.z*(DDF*dtD.z*dtD.y+D2L.bz);
          vecl.z=Dab.z*(DDF*dtA.z*dtD.z+D2IL.cz)+Dcb.z*rbc*(DDF*dtC.z*dtD.z+D2KL.cz)+
                 Dcb.z*rbc*(DDF*dtD.z*dtD.z+D2L.cz)+Ddc.z*(DDF*dtD.z*dtD.z+D2L.cz)+fd.z;

          HessianMatrix.element[index_i][n+5]+=veci.x;
          HessianMatrix.element[index_i+1][n+5]+=veci.y;
          HessianMatrix.element[index_i+2][n+5]+=veci.z;
          HessianMatrix.element[index_j][n+5]+=vecj.x;
          HessianMatrix.element[index_j+1][n+5]+=vecj.y;
          HessianMatrix.element[index_j+2][n+5]+=vecj.z;
          HessianMatrix.element[index_k][n+5]+=veck.x;
          HessianMatrix.element[index_k+1][n+5]+=veck.y;
          HessianMatrix.element[index_k+2][n+5]+=veck.z;
          HessianMatrix.element[index_l][n+5]+=vecl.x;
          HessianMatrix.element[index_l+1][n+5]+=vecl.y;
          HessianMatrix.element[index_l+2][n+5]+=vecl.z;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              //S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
              veci.x=Dab.x*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.x*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.x*(DDF*dtD.x*dtA.x+D2IL.ax)+fa.x;
              veci.y=Dab.x*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.x*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.x*(DDF*dtD.x*dtA.y+D2IL.bx);
              veci.z=Dab.x*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.x*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.x*(DDF*dtD.x*dtA.z+D2IL.cx);

              vecj.x=Dab.x*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.x*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.x*(DDF*dtD.x*dtB.x+D2JL.ax)-fa.x-fc.x-fd.x;
              vecj.y=Dab.x*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.x*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.x*(DDF*dtD.x*dtB.y+D2JL.bx);
              vecj.z=Dab.x*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.x*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.x*(DDF*dtD.x*dtB.z+D2JL.cx);

              veck.x=Dab.x*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.x*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.x*(DDF*dtD.x*dtC.x+D2KL.ax)+fc.x+fd.x-fd.x;
              veck.y=Dab.x*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.x*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.x*(DDF*dtD.x*dtC.y+D2KL.bx);
              veck.z=Dab.x*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.x*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.x*(DDF*dtD.x*dtC.z+D2KL.cx);

              vecl.x=Dab.x*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.x*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.x*(DDF*dtD.x*dtD.x+D2L.ax)+fd.x;
              vecl.y=Dab.x*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.x*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.x*(DDF*dtD.x*dtD.y+D2L.ay);
              vecl.z=Dab.x*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.x*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.x*(DDF*dtD.x*dtD.z+D2L.az);

              HessianMatrix.element[index_i][n]+=veci.x;
              HessianMatrix.element[index_i+1][n]+=veci.y;
              HessianMatrix.element[index_i+2][n]+=veci.z;
              HessianMatrix.element[index_j][n]+=vecj.x;
              HessianMatrix.element[index_j+1][n]+=vecj.y;
              HessianMatrix.element[index_j+2][n]+=vecj.z;
              HessianMatrix.element[index_k][n]+=veck.x;
              HessianMatrix.element[index_k+1][n]+=veck.y;
              HessianMatrix.element[index_k+2][n]+=veck.z;
              HessianMatrix.element[index_l][n]+=vecl.x;
              HessianMatrix.element[index_l+1][n]+=vecl.y;
              HessianMatrix.element[index_l+2][n]+=vecl.z;

              veci.x=Dab.y*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.y*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.y*(DDF*dtD.y*dtA.x+D2IL.ay);
              veci.y=Dab.y*(DDF*dtA.y*dtA.y+D2I.by)+fa.y+Dcb.y*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.y*(DDF*dtD.y*dtA.y+D2IL.by);
              veci.z=Dab.y*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.y*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.y*(DDF*dtD.y*dtA.z+D2IL.cy);

              vecj.x=Dab.y*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.y*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.y*(DDF*dtD.y*dtB.x+D2JL.ay);
              vecj.y=Dab.y*(DDF*dtA.y*dtB.y+D2IJ.by)-fa.y+Dcb.y*rbc*(DDF*dtC.y*dtB.y+D2JK.by)-fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.y+D2JL.by)-fd.y+Ddc.y*(DDF*dtD.y*dtB.y+D2JL.by);
              vecj.z=Dab.y*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.y*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.y*(DDF*dtD.y*dtB.z+D2JL.cy);

              veck.x=Dab.y*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.y*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.y*(DDF*dtD.y*dtC.x+D2KL.ay);
              veck.y=Dab.y*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.y*rbc*(DDF*dtC.y*dtC.y+D2K.by)+fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+fd.y+Ddc.y*(DDF*dtD.y*dtC.y+D2KL.by)-fd.y;
              veck.z=Dab.y*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.y*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.y*(DDF*dtD.y*dtC.z+D2KL.cy);

              vecl.x=Dab.y*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.y*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.y*(DDF*dtD.y*dtD.x+D2L.ay);
              vecl.y=Dab.y*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.y*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.y*(DDF*dtD.y*dtD.y+D2L.by)+fd.y;
              vecl.z=Dab.y*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.y*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.y*(DDF*dtD.y*dtD.z+D2L.bz);

              HessianMatrix.element[index_i][n+1]+=veci.x;
              HessianMatrix.element[index_i+1][n+1]+=veci.y;
              HessianMatrix.element[index_i+2][n+1]+=veci.z;
              HessianMatrix.element[index_j][n+1]+=vecj.x;
              HessianMatrix.element[index_j+1][n+1]+=vecj.y;
              HessianMatrix.element[index_j+2][n+1]+=vecj.z;
              HessianMatrix.element[index_k][n+1]+=veck.x;
              HessianMatrix.element[index_k+1][n+1]+=veck.y;
              HessianMatrix.element[index_k+2][n+1]+=veck.z;
              HessianMatrix.element[index_l][n+1]+=vecl.x;
              HessianMatrix.element[index_l+1][n+1]+=vecl.y;
              HessianMatrix.element[index_l+2][n+1]+=vecl.z;


              veci.x=Dab.z*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.z*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                     Dcb.z*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.z*(DDF*dtD.y*dtA.x+D2IL.ay);
              veci.y=Dab.z*(DDF*dtA.y*dtA.y+D2I.by)+Dcb.z*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                     Dcb.z*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.z*(DDF*dtD.y*dtA.y+D2IL.by);
              veci.z=Dab.z*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.z*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                     Dcb.z*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.z*(DDF*dtD.y*dtA.z+D2IL.cy)+fa.y;

              vecj.x=Dab.z*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.z*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                     Dcb.z*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.z*(DDF*dtD.y*dtB.x+D2JL.ay);
              vecj.y=Dab.z*(DDF*dtA.y*dtB.y+D2IJ.by)+Dcb.z*rbc*(DDF*dtC.y*dtB.y+D2JK.by)+
                     Dcb.z*rbc*(DDF*dtD.y*dtB.y+D2JL.by)+Ddc.z*(DDF*dtD.y*dtB.y+D2JL.by);
              vecj.z=Dab.z*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.z*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                     Dcb.z*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.z*(DDF*dtD.y*dtB.z+D2JL.cy)-fa.y-fc.y-fd.y;

              veck.x=Dab.z*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.z*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                     Dcb.z*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.z*(DDF*dtD.y*dtC.x+D2KL.ay);
              veck.y=Dab.z*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.z*rbc*(DDF*dtC.y*dtC.y+D2K.by)+
                     Dcb.z*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+Ddc.z*(DDF*dtD.y*dtC.y+D2KL.by);
              veck.z=Dab.z*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.z*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                     Dcb.z*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.z*(DDF*dtD.y*dtC.z+D2KL.cy)+fc.y+fd.y-fd.y;

              vecl.x=Dab.z*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.z*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                     Dcb.z*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.z*(DDF*dtD.y*dtD.x+D2L.ay);
              vecl.y=Dab.z*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.z*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                     Dcb.z*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.z*(DDF*dtD.y*dtD.y+D2L.by);
              vecl.z=Dab.z*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.z*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                     Dcb.z*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.z*(DDF*dtD.y*dtD.z+D2L.bz)+fd.y;

              HessianMatrix.element[index_i][n+2]+=veci.x;
              HessianMatrix.element[index_i+1][n+2]+=veci.y;
              HessianMatrix.element[index_i+2][n+2]+=veci.z;
              HessianMatrix.element[index_j][n+2]+=vecj.x;
              HessianMatrix.element[index_j+1][n+2]+=vecj.y;
              HessianMatrix.element[index_j+2][n+2]+=vecj.z;
              HessianMatrix.element[index_k][n+2]+=veck.x;
              HessianMatrix.element[index_k+1][n+2]+=veck.y;
              HessianMatrix.element[index_k+2][n+2]+=veck.z;
              HessianMatrix.element[index_l][n+2]+=vecl.x;
              HessianMatrix.element[index_l+1][n+2]+=vecl.y;
              HessianMatrix.element[index_l+2][n+2]+=vecl.z;

              veci.x=Dab.z*(DDF*dtA.z*dtA.x+D2I.az)+Dcb.z*rbc*(DDF*dtC.z*dtA.x+D2IK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.x+D2IL.az)+Ddc.z*(DDF*dtD.z*dtA.x+D2IL.az);
              veci.y=Dab.z*(DDF*dtA.z*dtA.y+D2I.bz)+Dcb.z*rbc*(DDF*dtC.z*dtA.y+D2IK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.y+D2IL.bz)+Ddc.z*(DDF*dtD.z*dtA.y+D2IL.bz);
              veci.z=Dab.z*(DDF*dtA.z*dtA.z+D2I.cz)+fa.z+Dcb.z*rbc*(DDF*dtC.z*dtA.z+D2IK.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.z+D2IL.cz)+Ddc.z*(DDF*dtD.z*dtA.z+D2IL.cz);

              vecj.x=Dab.z*(DDF*dtA.z*dtB.x+D2IJ.cx)+Dcb.z*rbc*(DDF*dtC.z*dtB.x+D2JK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.x+D2JL.az)+Ddc.z*(DDF*dtD.z*dtB.x+D2JL.az);
              vecj.y=Dab.z*(DDF*dtA.z*dtB.y+D2IJ.cy)+Dcb.z*rbc*(DDF*dtC.z*dtB.y+D2JK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.y+D2JL.bz)+Ddc.z*(DDF*dtD.z*dtB.y+D2JL.bz);
              vecj.z=Dab.z*(DDF*dtA.z*dtB.z+D2IJ.cz)-fa.z+Dcb.z*rbc*(DDF*dtC.z*dtB.z+D2JK.cz)-fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.z+D2JL.cz)-fd.z+Ddc.z*(DDF*dtD.z*dtB.z+D2JL.cz);

              veck.x=Dab.z*(DDF*dtA.z*dtC.x+D2IK.cx)+Dcb.z*rbc*(DDF*dtC.z*dtC.x+D2K.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.x+D2KL.az)+Ddc.z*(DDF*dtD.z*dtC.x+D2KL.az);
              veck.y=Dab.z*(DDF*dtA.z*dtC.y+D2IK.cy)+Dcb.z*rbc*(DDF*dtC.z*dtC.y+D2K.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.y+D2KL.bz)+Ddc.z*(DDF*dtD.z*dtC.y+D2KL.bz);
              veck.z=Dab.z*(DDF*dtA.z*dtC.z+D2IK.cz)+Dcb.z*rbc*(DDF*dtC.z*dtC.z+D2K.cz)+fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.z+D2KL.cz)+fd.z+Ddc.z*(DDF*dtD.z*dtC.z+D2KL.cz)-fd.z;

              vecl.x=Dab.z*(DDF*dtA.z*dtD.x+D2IL.cx)+Dcb.z*rbc*(DDF*dtC.z*dtD.x+D2KL.cx)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.x+D2L.az)+Ddc.z*(DDF*dtD.z*dtD.x+D2L.az);
              vecl.y=Dab.z*(DDF*dtA.z*dtD.y+D2IL.cy)+Dcb.z*rbc*(DDF*dtC.z*dtD.y+D2KL.cy)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.y+D2L.bz)+Ddc.z*(DDF*dtD.z*dtD.y+D2L.bz);
              vecl.z=Dab.z*(DDF*dtA.z*dtD.z+D2IL.cz)+Dcb.z*rbc*(DDF*dtC.z*dtD.z+D2KL.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.z+D2L.cz)+Ddc.z*(DDF*dtD.z*dtD.z+D2L.cz)+fd.z;

              HessianMatrix.element[index_i][n+3]+=veci.x;
              HessianMatrix.element[index_i+1][n+3]+=veci.y;
              HessianMatrix.element[index_i+2][n+3]+=veci.z;
              HessianMatrix.element[index_j][n+3]+=vecj.x;
              HessianMatrix.element[index_j+1][n+3]+=vecj.y;
              HessianMatrix.element[index_j+2][n+3]+=vecj.z;
              HessianMatrix.element[index_k][n+3]+=veck.x;
              HessianMatrix.element[index_k+1][n+3]+=veck.y;
              HessianMatrix.element[index_k+2][n+3]+=veck.z;
              HessianMatrix.element[index_l][n+3]+=vecl.x;
              HessianMatrix.element[index_l+1][n+3]+=vecl.y;
              HessianMatrix.element[index_l+2][n+3]+=vecl.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              //S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
              veci.x=Dab.x*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.x*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.x*(DDF*dtD.x*dtA.x+D2IL.ax)+fa.x;
              veci.y=Dab.x*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.x*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.x*(DDF*dtD.x*dtA.y+D2IL.bx);
              veci.z=Dab.x*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.x*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.x*(DDF*dtD.x*dtA.z+D2IL.cx);

              vecj.x=Dab.x*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.x*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.x*(DDF*dtD.x*dtB.x+D2JL.ax)-fa.x-fc.x-fd.x;
              vecj.y=Dab.x*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.x*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.x*(DDF*dtD.x*dtB.y+D2JL.bx);
              vecj.z=Dab.x*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.x*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.x*(DDF*dtD.x*dtB.z+D2JL.cx);

              veck.x=Dab.x*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.x*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.x*(DDF*dtD.x*dtC.x+D2KL.ax)+fc.x+fd.x-fd.x;
              veck.y=Dab.x*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.x*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.x*(DDF*dtD.x*dtC.y+D2KL.bx);
              veck.z=Dab.x*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.x*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.x*(DDF*dtD.x*dtC.z+D2KL.cx);

              vecl.x=Dab.x*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.x*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.x*(DDF*dtD.x*dtD.x+D2L.ax)+fd.x;
              vecl.y=Dab.x*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.x*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.x*(DDF*dtD.x*dtD.y+D2L.ay);
              vecl.z=Dab.x*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.x*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.x*(DDF*dtD.x*dtD.z+D2L.az);

              HessianMatrix.element[index_i][n]+=veci.x;
              HessianMatrix.element[index_i+1][n]+=veci.y;
              HessianMatrix.element[index_i+2][n]+=veci.z;
              HessianMatrix.element[index_j][n]+=vecj.x;
              HessianMatrix.element[index_j+1][n]+=vecj.y;
              HessianMatrix.element[index_j+2][n]+=vecj.z;
              HessianMatrix.element[index_k][n]+=veck.x;
              HessianMatrix.element[index_k+1][n]+=veck.y;
              HessianMatrix.element[index_k+2][n]+=veck.z;
              HessianMatrix.element[index_l][n]+=vecl.x;
              HessianMatrix.element[index_l+1][n]+=vecl.y;
              HessianMatrix.element[index_l+2][n]+=vecl.z;

              veci.x=Dab.z*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.z*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                    Dcb.z*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.z*(DDF*dtD.x*dtA.x+D2IL.ax);
              veci.y=Dab.z*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.z*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                    Dcb.z*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.z*(DDF*dtD.x*dtA.y+D2IL.bx);
              veci.z=Dab.z*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.z*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                     Dcb.z*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.z*(DDF*dtD.x*dtA.z+D2IL.cx)+fa.x;

              vecj.x=Dab.z*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.z*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                     Dcb.z*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.z*(DDF*dtD.x*dtB.x+D2JL.ax);
              vecj.y=Dab.z*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.z*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                     Dcb.z*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.z*(DDF*dtD.x*dtB.y+D2JL.bx);
              vecj.z=Dab.z*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.z*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                     Dcb.z*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.z*(DDF*dtD.x*dtB.z+D2JL.cx)-fa.x-fc.x-fd.x;

              veck.x=Dab.z*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.z*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                     Dcb.z*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.z*(DDF*dtD.x*dtC.x+D2KL.ax);
              veck.y=Dab.z*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.z*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                     Dcb.z*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.z*(DDF*dtD.x*dtC.y+D2KL.bx);
              veck.z=Dab.z*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.z*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                     Dcb.z*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.z*(DDF*dtD.x*dtC.z+D2KL.cx)+fc.x+fd.x-fd.x;

              vecl.x=Dab.z*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.z*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                     Dcb.z*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.z*(DDF*dtD.x*dtD.x+D2L.ax);
              vecl.y=Dab.z*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.z*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                     Dcb.z*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.z*(DDF*dtD.x*dtD.y+D2L.ay);
              vecl.z=Dab.z*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.z*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                     Dcb.z*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.z*(DDF*dtD.x*dtD.z+D2L.az)+fd.x;

              HessianMatrix.element[index_i][n+1]+=veci.x;
              HessianMatrix.element[index_i+1][n+1]+=veci.y;
              HessianMatrix.element[index_i+2][n+1]+=veci.z;
              HessianMatrix.element[index_j][n+1]+=vecj.x;
              HessianMatrix.element[index_j+1][n+1]+=vecj.y;
              HessianMatrix.element[index_j+2][n+1]+=vecj.z;
              HessianMatrix.element[index_k][n+1]+=veck.x;
              HessianMatrix.element[index_k+1][n+1]+=veck.y;
              HessianMatrix.element[index_k+2][n+1]+=veck.z;
              HessianMatrix.element[index_l][n+1]+=vecl.x;
              HessianMatrix.element[index_l+1][n+1]+=vecl.y;
              HessianMatrix.element[index_l+2][n+1]+=vecl.z;


              veci.x=Dab.y*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.y*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.y*(DDF*dtD.y*dtA.x+D2IL.ay);
              veci.y=Dab.y*(DDF*dtA.y*dtA.y+D2I.by)+fa.y+Dcb.y*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.y*(DDF*dtD.y*dtA.y+D2IL.by);
              veci.z=Dab.y*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.y*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.y*(DDF*dtD.y*dtA.z+D2IL.cy);

              vecj.x=Dab.y*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.y*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.y*(DDF*dtD.y*dtB.x+D2JL.ay);
              vecj.y=Dab.y*(DDF*dtA.y*dtB.y+D2IJ.by)-fa.y+Dcb.y*rbc*(DDF*dtC.y*dtB.y+D2JK.by)-fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.y+D2JL.by)-fd.y+Ddc.y*(DDF*dtD.y*dtB.y+D2JL.by);
              vecj.z=Dab.y*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.y*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.y*(DDF*dtD.y*dtB.z+D2JL.cy);

              veck.x=Dab.y*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.y*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.y*(DDF*dtD.y*dtC.x+D2KL.ay);
              veck.y=Dab.y*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.y*rbc*(DDF*dtC.y*dtC.y+D2K.by)+fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+fd.y+Ddc.y*(DDF*dtD.y*dtC.y+D2KL.by)-fd.y;
              veck.z=Dab.y*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.y*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.y*(DDF*dtD.y*dtC.z+D2KL.cy);

              vecl.x=Dab.y*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.y*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.y*(DDF*dtD.y*dtD.x+D2L.ay);
              vecl.y=Dab.y*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.y*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.y*(DDF*dtD.y*dtD.y+D2L.by)+fd.y;
              vecl.z=Dab.y*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.y*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.y*(DDF*dtD.y*dtD.z+D2L.bz);

              HessianMatrix.element[index_i][n+2]+=veci.x;
              HessianMatrix.element[index_i+1][n+2]+=veci.y;
              HessianMatrix.element[index_i+2][n+2]+=veci.z;
              HessianMatrix.element[index_j][n+2]+=vecj.x;
              HessianMatrix.element[index_j+1][n+2]+=vecj.y;
              HessianMatrix.element[index_j+2][n+2]+=vecj.z;
              HessianMatrix.element[index_k][n+2]+=veck.x;
              HessianMatrix.element[index_k+1][n+2]+=veck.y;
              HessianMatrix.element[index_k+2][n+2]+=veck.z;
              HessianMatrix.element[index_l][n+2]+=vecl.x;
              HessianMatrix.element[index_l+1][n+2]+=vecl.y;
              HessianMatrix.element[index_l+2][n+2]+=vecl.z;


              veci.x=Dab.z*(DDF*dtA.z*dtA.x+D2I.az)+Dcb.z*rbc*(DDF*dtC.z*dtA.x+D2IK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.x+D2IL.az)+Ddc.z*(DDF*dtD.z*dtA.x+D2IL.az);
              veci.y=Dab.z*(DDF*dtA.z*dtA.y+D2I.bz)+Dcb.z*rbc*(DDF*dtC.z*dtA.y+D2IK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.y+D2IL.bz)+Ddc.z*(DDF*dtD.z*dtA.y+D2IL.bz);
              veci.z=Dab.z*(DDF*dtA.z*dtA.z+D2I.cz)+fa.z+Dcb.z*rbc*(DDF*dtC.z*dtA.z+D2IK.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.z+D2IL.cz)+Ddc.z*(DDF*dtD.z*dtA.z+D2IL.cz);

              vecj.x=Dab.z*(DDF*dtA.z*dtB.x+D2IJ.cx)+Dcb.z*rbc*(DDF*dtC.z*dtB.x+D2JK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.x+D2JL.az)+Ddc.z*(DDF*dtD.z*dtB.x+D2JL.az);
              vecj.y=Dab.z*(DDF*dtA.z*dtB.y+D2IJ.cy)+Dcb.z*rbc*(DDF*dtC.z*dtB.y+D2JK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.y+D2JL.bz)+Ddc.z*(DDF*dtD.z*dtB.y+D2JL.bz);
              vecj.z=Dab.z*(DDF*dtA.z*dtB.z+D2IJ.cz)-fa.z+Dcb.z*rbc*(DDF*dtC.z*dtB.z+D2JK.cz)-fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.z+D2JL.cz)-fd.z+Ddc.z*(DDF*dtD.z*dtB.z+D2JL.cz);

              veck.x=Dab.z*(DDF*dtA.z*dtC.x+D2IK.cx)+Dcb.z*rbc*(DDF*dtC.z*dtC.x+D2K.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.x+D2KL.az)+Ddc.z*(DDF*dtD.z*dtC.x+D2KL.az);
              veck.y=Dab.z*(DDF*dtA.z*dtC.y+D2IK.cy)+Dcb.z*rbc*(DDF*dtC.z*dtC.y+D2K.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.y+D2KL.bz)+Ddc.z*(DDF*dtD.z*dtC.y+D2KL.bz);
              veck.z=Dab.z*(DDF*dtA.z*dtC.z+D2IK.cz)+Dcb.z*rbc*(DDF*dtC.z*dtC.z+D2K.cz)+fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.z+D2KL.cz)+fd.z+Ddc.z*(DDF*dtD.z*dtC.z+D2KL.cz)-fd.z;

              vecl.x=Dab.z*(DDF*dtA.z*dtD.x+D2IL.cx)+Dcb.z*rbc*(DDF*dtC.z*dtD.x+D2KL.cx)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.x+D2L.az)+Ddc.z*(DDF*dtD.z*dtD.x+D2L.az);
              vecl.y=Dab.z*(DDF*dtA.z*dtD.y+D2IL.cy)+Dcb.z*rbc*(DDF*dtC.z*dtD.y+D2KL.cy)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.y+D2L.bz)+Ddc.z*(DDF*dtD.z*dtD.y+D2L.bz);
              vecl.z=Dab.z*(DDF*dtA.z*dtD.z+D2IL.cz)+Dcb.z*rbc*(DDF*dtC.z*dtD.z+D2KL.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.z+D2L.cz)+Ddc.z*(DDF*dtD.z*dtD.z+D2L.cz)+fd.z;

              HessianMatrix.element[index_i][n+3]+=veci.x;
              HessianMatrix.element[index_i+1][n+3]+=veci.y;
              HessianMatrix.element[index_i+2][n+3]+=veci.z;
              HessianMatrix.element[index_j][n+3]+=vecj.x;
              HessianMatrix.element[index_j+1][n+3]+=vecj.y;
              HessianMatrix.element[index_j+2][n+3]+=vecj.z;
              HessianMatrix.element[index_k][n+3]+=veck.x;
              HessianMatrix.element[index_k+1][n+3]+=veck.y;
              HessianMatrix.element[index_k+2][n+3]+=veck.z;
              HessianMatrix.element[index_l][n+3]+=vecl.x;
              HessianMatrix.element[index_l+1][n+3]+=vecl.y;
              HessianMatrix.element[index_l+2][n+3]+=vecl.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              //S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
              veci.x=Dab.x*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.x*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.x*(DDF*dtD.x*dtA.x+D2IL.ax)+fa.x;
              veci.y=Dab.x*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.x*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.x*(DDF*dtD.x*dtA.y+D2IL.bx);
              veci.z=Dab.x*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.x*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.x*(DDF*dtD.x*dtA.z+D2IL.cx);

              vecj.x=Dab.x*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.x*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.x*(DDF*dtD.x*dtB.x+D2JL.ax)-fa.x-fc.x-fd.x;
              vecj.y=Dab.x*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.x*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.x*(DDF*dtD.x*dtB.y+D2JL.bx);
              vecj.z=Dab.x*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.x*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                     Dcb.x*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.x*(DDF*dtD.x*dtB.z+D2JL.cx);

              veck.x=Dab.x*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.x*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.x*(DDF*dtD.x*dtC.x+D2KL.ax)+fc.x+fd.x-fd.x;
              veck.y=Dab.x*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.x*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.x*(DDF*dtD.x*dtC.y+D2KL.bx);
              veck.z=Dab.x*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.x*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.x*(DDF*dtD.x*dtC.z+D2KL.cx);

              vecl.x=Dab.x*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.x*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.x*(DDF*dtD.x*dtD.x+D2L.ax)+fd.x;
              vecl.y=Dab.x*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.x*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.x*(DDF*dtD.x*dtD.y+D2L.ay);
              vecl.z=Dab.x*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.x*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                     Dcb.x*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.x*(DDF*dtD.x*dtD.z+D2L.az);

              HessianMatrix.element[index_i][n]+=veci.x;
              HessianMatrix.element[index_i+1][n]+=veci.y;
              HessianMatrix.element[index_i+2][n]+=veci.z;
              HessianMatrix.element[index_j][n]+=vecj.x;
              HessianMatrix.element[index_j+1][n]+=vecj.y;
              HessianMatrix.element[index_j+2][n]+=vecj.z;
              HessianMatrix.element[index_k][n]+=veck.x;
              HessianMatrix.element[index_k+1][n]+=veck.y;
              HessianMatrix.element[index_k+2][n]+=veck.z;
              HessianMatrix.element[index_l][n]+=vecl.x;
              HessianMatrix.element[index_l+1][n]+=vecl.y;
              HessianMatrix.element[index_l+2][n]+=vecl.z;

              //S.bx=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
              veci.x=Dab.y*(DDF*dtA.x*dtA.x+D2I.ax)+Dcb.y*rbc*(DDF*dtC.x*dtA.x+D2IK.ax)+
                     Dcb.y*rbc*(DDF*dtD.x*dtA.x+D2IL.ax)+Ddc.y*(DDF*dtD.x*dtA.x+D2IL.ax);
              veci.y=Dab.y*(DDF*dtA.x*dtA.y+D2I.ay)+Dcb.y*rbc*(DDF*dtC.x*dtA.y+D2IK.bx)+
                     Dcb.y*rbc*(DDF*dtD.x*dtA.y+D2IL.bx)+Ddc.y*(DDF*dtD.x*dtA.y+D2IL.bx)+fa.x;
              veci.z=Dab.y*(DDF*dtA.x*dtA.z+D2I.az)+Dcb.y*rbc*(DDF*dtC.x*dtA.z+D2IK.cx)+
                     Dcb.y*rbc*(DDF*dtD.x*dtA.z+D2IL.cx)+Ddc.y*(DDF*dtD.x*dtA.z+D2IL.cx);

              vecj.x=Dab.y*(DDF*dtA.x*dtB.x+D2IJ.ax)+Dcb.y*rbc*(DDF*dtC.x*dtB.x+D2JK.ax)+
                     Dcb.y*rbc*(DDF*dtD.x*dtB.x+D2JL.ax)+Ddc.y*(DDF*dtD.x*dtB.x+D2JL.ax);
              vecj.y=Dab.y*(DDF*dtA.x*dtB.y+D2IJ.ay)+Dcb.y*rbc*(DDF*dtC.x*dtB.y+D2JK.bx)+
                     Dcb.y*rbc*(DDF*dtD.x*dtB.y+D2JL.bx)+Ddc.y*(DDF*dtD.x*dtB.y+D2JL.bx)-fa.x-fc.x-fd.x;
              vecj.z=Dab.y*(DDF*dtA.x*dtB.z+D2IJ.az)+Dcb.y*rbc*(DDF*dtC.x*dtB.z+D2JK.cx)+
                     Dcb.y*rbc*(DDF*dtD.x*dtB.z+D2JL.cx)+Ddc.y*(DDF*dtD.x*dtB.z+D2JL.cx);

              veck.x=Dab.y*(DDF*dtA.x*dtC.x+D2IK.ax)+Dcb.y*rbc*(DDF*dtC.x*dtC.x+D2K.ax)+
                     Dcb.y*rbc*(DDF*dtD.x*dtC.x+D2KL.ax)+Ddc.y*(DDF*dtD.x*dtC.x+D2KL.ax);
              veck.y=Dab.y*(DDF*dtA.x*dtC.y+D2IK.ay)+Dcb.y*rbc*(DDF*dtC.x*dtC.y+D2K.ay)+
                     Dcb.y*rbc*(DDF*dtD.x*dtC.y+D2KL.bx)+Ddc.y*(DDF*dtD.x*dtC.y+D2KL.bx)+fc.x+fd.x-fd.x;
              veck.z=Dab.y*(DDF*dtA.x*dtC.z+D2IK.az)+Dcb.y*rbc*(DDF*dtC.x*dtC.z+D2K.az)+
                     Dcb.y*rbc*(DDF*dtD.x*dtC.z+D2KL.cx)+Ddc.y*(DDF*dtD.x*dtC.z+D2KL.cx);

              vecl.x=Dab.y*(DDF*dtA.x*dtD.x+D2IL.ax)+Dcb.y*rbc*(DDF*dtC.x*dtD.x+D2KL.ax)+
                     Dcb.y*rbc*(DDF*dtD.x*dtD.x+D2L.ax)+Ddc.y*(DDF*dtD.x*dtD.x+D2L.ax);
              vecl.y=Dab.y*(DDF*dtA.x*dtD.y+D2IL.ay)+Dcb.y*rbc*(DDF*dtC.x*dtD.y+D2KL.ay)+
                     Dcb.y*rbc*(DDF*dtD.x*dtD.y+D2L.ay)+Ddc.y*(DDF*dtD.x*dtD.y+D2L.ay)+fd.x;
              vecl.z=Dab.y*(DDF*dtA.x*dtD.z+D2IL.az)+Dcb.y*rbc*(DDF*dtC.x*dtD.z+D2KL.az)+
                     Dcb.y*rbc*(DDF*dtD.x*dtD.z+D2L.az)+Ddc.y*(DDF*dtD.x*dtD.z+D2L.az);

              HessianMatrix.element[index_i][n+1]+=veci.x;
              HessianMatrix.element[index_i+1][n+1]+=veci.y;
              HessianMatrix.element[index_i+2][n+1]+=veci.z;
              HessianMatrix.element[index_j][n+1]+=vecj.x;
              HessianMatrix.element[index_j+1][n+1]+=vecj.y;
              HessianMatrix.element[index_j+2][n+1]+=vecj.z;
              HessianMatrix.element[index_k][n+1]+=veck.x;
              HessianMatrix.element[index_k+1][n+1]+=veck.y;
              HessianMatrix.element[index_k+2][n+1]+=veck.z;
              HessianMatrix.element[index_l][n+1]+=vecl.x;
              HessianMatrix.element[index_l+1][n+1]+=vecl.y;
              HessianMatrix.element[index_l+2][n+1]+=vecl.z;

              veci.x=Dab.y*(DDF*dtA.y*dtA.x+D2I.ay)+Dcb.y*rbc*(DDF*dtC.y*dtA.x+D2IK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.x+D2IL.ay)+Ddc.y*(DDF*dtD.y*dtA.x+D2IL.ay);
              veci.y=Dab.y*(DDF*dtA.y*dtA.y+D2I.by)+fa.y+Dcb.y*rbc*(DDF*dtC.y*dtA.y+D2IK.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.y+D2IL.by)+Ddc.y*(DDF*dtD.y*dtA.y+D2IL.by);
              veci.z=Dab.y*(DDF*dtA.y*dtA.z+D2I.bz)+Dcb.y*rbc*(DDF*dtC.y*dtA.z+D2IK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtA.z+D2IL.cy)+Ddc.y*(DDF*dtD.y*dtA.z+D2IL.cy);

              vecj.x=Dab.y*(DDF*dtA.y*dtB.x+D2IJ.bx)+Dcb.y*rbc*(DDF*dtC.y*dtB.x+D2JK.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.x+D2JL.ay)+Ddc.y*(DDF*dtD.y*dtB.x+D2JL.ay);
              vecj.y=Dab.y*(DDF*dtA.y*dtB.y+D2IJ.by)-fa.y+Dcb.y*rbc*(DDF*dtC.y*dtB.y+D2JK.by)-fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.y+D2JL.by)-fd.y+Ddc.y*(DDF*dtD.y*dtB.y+D2JL.by);
              vecj.z=Dab.y*(DDF*dtA.y*dtB.z+D2IJ.bz)+Dcb.y*rbc*(DDF*dtC.y*dtB.z+D2JK.cy)+
                     Dcb.y*rbc*(DDF*dtD.y*dtB.z+D2JL.cy)+Ddc.y*(DDF*dtD.y*dtB.z+D2JL.cy);

              veck.x=Dab.y*(DDF*dtA.y*dtC.x+D2IK.bx)+Dcb.y*rbc*(DDF*dtC.y*dtC.x+D2K.ay)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.x+D2KL.ay)+Ddc.y*(DDF*dtD.y*dtC.x+D2KL.ay);
              veck.y=Dab.y*(DDF*dtA.y*dtC.y+D2IK.by)+Dcb.y*rbc*(DDF*dtC.y*dtC.y+D2K.by)+fc.y+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.y+D2KL.by)+fd.y+Ddc.y*(DDF*dtD.y*dtC.y+D2KL.by)-fd.y;
              veck.z=Dab.y*(DDF*dtA.y*dtC.z+D2IK.bz)+Dcb.y*rbc*(DDF*dtC.y*dtC.z+D2K.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtC.z+D2KL.cy)+Ddc.y*(DDF*dtD.y*dtC.z+D2KL.cy);

              vecl.x=Dab.y*(DDF*dtA.y*dtD.x+D2IL.bx)+Dcb.y*rbc*(DDF*dtC.y*dtD.x+D2KL.bx)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.x+D2L.ay)+Ddc.y*(DDF*dtD.y*dtD.x+D2L.ay);
              vecl.y=Dab.y*(DDF*dtA.y*dtD.y+D2IL.by)+Dcb.y*rbc*(DDF*dtC.y*dtD.y+D2KL.by)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.y+D2L.by)+Ddc.y*(DDF*dtD.y*dtD.y+D2L.by)+fd.y;
              vecl.z=Dab.y*(DDF*dtA.y*dtD.z+D2IL.bz)+Dcb.y*rbc*(DDF*dtC.y*dtD.z+D2KL.bz)+
                     Dcb.y*rbc*(DDF*dtD.y*dtD.z+D2L.bz)+Ddc.y*(DDF*dtD.y*dtD.z+D2L.bz);

              HessianMatrix.element[index_i][n+2]+=veci.x;
              HessianMatrix.element[index_i+1][n+2]+=veci.y;
              HessianMatrix.element[index_i+2][n+2]+=veci.z;
              HessianMatrix.element[index_j][n+2]+=vecj.x;
              HessianMatrix.element[index_j+1][n+2]+=vecj.y;
              HessianMatrix.element[index_j+2][n+2]+=vecj.z;
              HessianMatrix.element[index_k][n+2]+=veck.x;
              HessianMatrix.element[index_k+1][n+2]+=veck.y;
              HessianMatrix.element[index_k+2][n+2]+=veck.z;
              HessianMatrix.element[index_l][n+2]+=vecl.x;
              HessianMatrix.element[index_l+1][n+2]+=vecl.y;
              HessianMatrix.element[index_l+2][n+2]+=vecl.z;


              veci.x=Dab.z*(DDF*dtA.z*dtA.x+D2I.az)+Dcb.z*rbc*(DDF*dtC.z*dtA.x+D2IK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.x+D2IL.az)+Ddc.z*(DDF*dtD.z*dtA.x+D2IL.az);
              veci.y=Dab.z*(DDF*dtA.z*dtA.y+D2I.bz)+Dcb.z*rbc*(DDF*dtC.z*dtA.y+D2IK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.y+D2IL.bz)+Ddc.z*(DDF*dtD.z*dtA.y+D2IL.bz);
              veci.z=Dab.z*(DDF*dtA.z*dtA.z+D2I.cz)+fa.z+Dcb.z*rbc*(DDF*dtC.z*dtA.z+D2IK.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtA.z+D2IL.cz)+Ddc.z*(DDF*dtD.z*dtA.z+D2IL.cz);

              vecj.x=Dab.z*(DDF*dtA.z*dtB.x+D2IJ.cx)+Dcb.z*rbc*(DDF*dtC.z*dtB.x+D2JK.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.x+D2JL.az)+Ddc.z*(DDF*dtD.z*dtB.x+D2JL.az);
              vecj.y=Dab.z*(DDF*dtA.z*dtB.y+D2IJ.cy)+Dcb.z*rbc*(DDF*dtC.z*dtB.y+D2JK.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.y+D2JL.bz)+Ddc.z*(DDF*dtD.z*dtB.y+D2JL.bz);
              vecj.z=Dab.z*(DDF*dtA.z*dtB.z+D2IJ.cz)-fa.z+Dcb.z*rbc*(DDF*dtC.z*dtB.z+D2JK.cz)-fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtB.z+D2JL.cz)-fd.z+Ddc.z*(DDF*dtD.z*dtB.z+D2JL.cz);

              veck.x=Dab.z*(DDF*dtA.z*dtC.x+D2IK.cx)+Dcb.z*rbc*(DDF*dtC.z*dtC.x+D2K.az)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.x+D2KL.az)+Ddc.z*(DDF*dtD.z*dtC.x+D2KL.az);
              veck.y=Dab.z*(DDF*dtA.z*dtC.y+D2IK.cy)+Dcb.z*rbc*(DDF*dtC.z*dtC.y+D2K.bz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.y+D2KL.bz)+Ddc.z*(DDF*dtD.z*dtC.y+D2KL.bz);
              veck.z=Dab.z*(DDF*dtA.z*dtC.z+D2IK.cz)+Dcb.z*rbc*(DDF*dtC.z*dtC.z+D2K.cz)+fc.z+
                     Dcb.z*rbc*(DDF*dtD.z*dtC.z+D2KL.cz)+fd.z+Ddc.z*(DDF*dtD.z*dtC.z+D2KL.cz)-fd.z;

              vecl.x=Dab.z*(DDF*dtA.z*dtD.x+D2IL.cx)+Dcb.z*rbc*(DDF*dtC.z*dtD.x+D2KL.cx)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.x+D2L.az)+Ddc.z*(DDF*dtD.z*dtD.x+D2L.az);
              vecl.y=Dab.z*(DDF*dtA.z*dtD.y+D2IL.cy)+Dcb.z*rbc*(DDF*dtC.z*dtD.y+D2KL.cy)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.y+D2L.bz)+Ddc.z*(DDF*dtD.z*dtD.y+D2L.bz);
              vecl.z=Dab.z*(DDF*dtA.z*dtD.z+D2IL.cz)+Dcb.z*rbc*(DDF*dtC.z*dtD.z+D2KL.cz)+
                     Dcb.z*rbc*(DDF*dtD.z*dtD.z+D2L.cz)+Ddc.z*(DDF*dtD.z*dtD.z+D2L.cz)+fd.z;

              HessianMatrix.element[index_i][n+3]+=veci.x;
              HessianMatrix.element[index_i+1][n+3]+=veci.y;
              HessianMatrix.element[index_i+2][n+3]+=veci.z;
              HessianMatrix.element[index_j][n+3]+=vecj.x;
              HessianMatrix.element[index_j+1][n+3]+=vecj.y;
              HessianMatrix.element[index_j+2][n+3]+=vecj.z;
              HessianMatrix.element[index_k][n+3]+=veck.x;
              HessianMatrix.element[index_k+1][n+3]+=veck.y;
              HessianMatrix.element[index_k+2][n+3]+=veck.z;
              HessianMatrix.element[index_l][n+3]+=vecl.x;
              HessianMatrix.element[index_l+1][n+3]+=vecl.y;
              HessianMatrix.element[index_l+2][n+3]+=vecl.z;
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}


static inline void HessianTorsionStrainStrain(REAL_MATRIX HessianMatrix,VECTOR vec_u,VECTOR vec_v,VECTOR vec_w,REAL u,REAL v,REAL w,REAL DF,REAL DDF,REAL_MATRIX3x3 S,REAL CosPhi)
{
  int n;
  REAL_MATRIX3x3 dme,dne,dmn,T;
  REAL m2,n2,dot_mn;
  VECTOR vec_m,vec_n;

  n=NumberOfCoordinatesMinimizationVariables;

  dme.ax=2.0*((SQR(u))*vec_v.x*vec_v.x+(SQR(v))*vec_u.x*vec_u.x-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x));
  dme.ay=2.0*((SQR(u))*vec_v.x*vec_v.y+(SQR(v))*vec_u.x*vec_u.y-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y));
  dme.az=2.0*((SQR(u))*vec_v.x*vec_v.z+(SQR(v))*vec_u.x*vec_u.z-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z));
  dme.bx=2.0*((SQR(u))*vec_v.y*vec_v.x+(SQR(v))*vec_u.y*vec_u.x-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.y*vec_u.x+vec_u.y*vec_v.x));
  dme.by=2.0*((SQR(u))*vec_v.y*vec_v.y+(SQR(v))*vec_u.y*vec_u.y-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y));
  dme.bz=2.0*((SQR(u))*vec_v.y*vec_v.z+(SQR(v))*vec_u.y*vec_u.z-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z));
  dme.cx=2.0*((SQR(u))*vec_v.z*vec_v.x+(SQR(v))*vec_u.z*vec_u.x-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x));
  dme.cy=2.0*((SQR(u))*vec_v.z*vec_v.y+(SQR(v))*vec_u.z*vec_u.y-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.z*vec_u.y+vec_u.z*vec_v.y));
  dme.cz=2.0*((SQR(u))*vec_v.z*vec_v.z+(SQR(v))*vec_u.z*vec_u.z-(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z));

  dne.ax=2.0*((SQR(u))*vec_w.x*vec_w.x+(SQR(w))*vec_u.x*vec_u.x-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x));
  dne.ay=2.0*((SQR(u))*vec_w.x*vec_w.y+(SQR(w))*vec_u.x*vec_u.y-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y));
  dne.az=2.0*((SQR(u))*vec_w.x*vec_w.z+(SQR(w))*vec_u.x*vec_u.z-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z));
  dne.bx=2.0*((SQR(u))*vec_w.y*vec_w.x+(SQR(w))*vec_u.y*vec_u.x-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.y*vec_u.x+vec_u.y*vec_w.x));
  dne.by=2.0*((SQR(u))*vec_w.y*vec_w.y+(SQR(w))*vec_u.y*vec_u.y-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y));
  dne.bz=2.0*((SQR(u))*vec_w.y*vec_w.z+(SQR(w))*vec_u.y*vec_u.z-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z));
  dne.cx=2.0*((SQR(u))*vec_w.z*vec_w.x+(SQR(w))*vec_u.z*vec_u.x-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.z*vec_u.x+vec_u.z*vec_w.x));
  dne.cy=2.0*((SQR(u))*vec_w.z*vec_w.y+(SQR(w))*vec_u.z*vec_u.y-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.z*vec_u.y+vec_u.z*vec_w.y));
  dne.cz=2.0*((SQR(u))*vec_w.z*vec_w.z+(SQR(w))*vec_u.z*vec_u.z-(vec_w.x*vec_u.x+vec_w.y*vec_u.y+vec_w.z*vec_u.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z));

  dmn.ax=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.x*vec_w.x+vec_w.x*vec_u.x)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x);
  dmn.ay=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.x*vec_w.y+vec_w.x*vec_u.y)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y);
  dmn.az=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.x*vec_w.z+vec_w.x*vec_u.z)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z);

  dmn.bx=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.y*vec_w.x+vec_w.y*vec_u.x)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.y*vec_u.x+vec_u.y*vec_v.x)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.y*vec_u.x+vec_u.y*vec_u.x)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.y*vec_w.x+vec_w.y*vec_v.x);
  dmn.by=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.y*vec_w.y+vec_w.y*vec_u.y)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y);
  dmn.bz=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.y*vec_w.z+vec_w.y*vec_u.z)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z);

  dmn.cx=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.z*vec_w.x+vec_w.z*vec_u.x)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.z*vec_u.x+vec_u.z*vec_v.x)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.z*vec_u.x+vec_u.z*vec_u.x)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.z*vec_w.x+vec_w.z*vec_v.x);
  dmn.cy=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.z*vec_w.y+vec_w.z*vec_u.y)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.z*vec_u.y+vec_u.z*vec_v.y)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.z*vec_u.y+vec_u.z*vec_u.y)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.z*vec_w.y+vec_w.z*vec_v.y);
  dmn.cz=(vec_v.x*vec_u.x+vec_v.y*vec_u.y+vec_v.z*vec_u.z)*(vec_u.z*vec_w.z+vec_w.z*vec_u.z)
            +(vec_u.x*vec_w.x+vec_u.y*vec_w.y+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
            -(vec_v.x*vec_w.x+vec_v.y*vec_w.y+vec_v.z*vec_w.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)
            -(vec_u.x*vec_u.x+vec_u.y*vec_u.y+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z);

  vec_m.x=vec_v.y*vec_u.z-vec_v.z*vec_u.y;
  vec_m.y=vec_v.z*vec_u.x-vec_v.x*vec_u.z;
  vec_m.z=vec_v.x*vec_u.y-vec_v.y*vec_u.x;
  m2=SQR(vec_m.x)+SQR(vec_m.y)+SQR(vec_m.z);

  vec_n.x=vec_u.y*vec_w.z-vec_u.z*vec_w.y;
  vec_n.y=vec_u.z*vec_w.x-vec_u.x*vec_w.z;
  vec_n.z=vec_u.x*vec_w.y-vec_u.y*vec_w.x;
  n2=SQR(vec_n.x)+SQR(vec_n.y)+SQR(vec_n.z);

  dot_mn=(vec_m.x*vec_n.x+vec_m.y*vec_n.y+vec_m.z*vec_n.z);

  T.ax=0.5*CosPhi*(2.0*dmn.ax/dot_mn-dme.ax/m2-dne.ax/n2);
  T.ay=0.5*CosPhi*(2.0*dmn.ay/dot_mn-dme.ay/m2-dne.ay/n2);
  T.az=0.5*CosPhi*(2.0*dmn.az/dot_mn-dme.az/m2-dne.az/n2);

  T.bx=0.5*CosPhi*(2.0*dmn.bx/dot_mn-dme.bx/m2-dne.bx/n2);
  T.by=0.5*CosPhi*(2.0*dmn.by/dot_mn-dme.by/m2-dne.by/n2);
  T.bz=0.5*CosPhi*(2.0*dmn.bz/dot_mn-dme.bz/m2-dne.bz/n2);

  T.cx=0.5*CosPhi*(2.0*dmn.cx/dot_mn-dme.cx/m2-dne.cx/n2);
  T.cy=0.5*CosPhi*(2.0*dmn.cy/dot_mn-dme.cy/m2-dne.cy/n2);
  T.cz=0.5*CosPhi*(2.0*dmn.cz/dot_mn-dme.cz/m2-dne.cz/n2);

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
          HessianMatrix.element[n][n]+=
             DDF*T.ax*T.ax+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
             +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
             +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
          HessianMatrix.element[n][n]+=
             2.0*(DDF*T.ax*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ax*T.by);
          HessianMatrix.element[n][n]+=
             2.0*(DDF*T.ax*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.cz);
          HessianMatrix.element[n][n]+=
             DDF*T.by*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
          HessianMatrix.element[n][n]+=
             2.0*(DDF*T.by*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.by*T.cz);
          HessianMatrix.element[n][n]+=
             DDF*T.cz*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
             ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          HessianMatrix.element[n][n]+=
             DDF*T.ax*T.ax+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
             +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
             +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
          HessianMatrix.element[n][n+1]+=
             DDF*T.ax*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ax*T.by;
          HessianMatrix.element[n][n+2]+=
             DDF*T.ax*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.cz;

          HessianMatrix.element[n+1][n+1]+=
             DDF*T.by*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
          HessianMatrix.element[n+1][n+2]+=
             DDF*T.by*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.by*T.cz;
          HessianMatrix.element[n+2][n+2]+=
             DDF*T.cz*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
             ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
          break;
        case REGULAR:
          HessianMatrix.element[n][n]+=
             DDF*T.ax*T.ax+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
             +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
             +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
          HessianMatrix.element[n][n+1]+=
             DDF*T.ax*T.ay+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.ay/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ay)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
             +(dme.ax*dme.ay/SQR(m2)+dne.ax*dne.ay/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ax*T.ay+S.ay;
          HessianMatrix.element[n][n+2]+=
             DDF*T.ax*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.ax*dme.az/SQR(m2)+dne.ax*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.az+S.az;
          HessianMatrix.element[n][n+3]+=
             DDF*T.ax*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ax*T.by;
          HessianMatrix.element[n][n+4]+=
             DDF*T.ax*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.bz)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.ax*dme.bz/SQR(m2)+dne.ax*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.bz;
          HessianMatrix.element[n][n+5]+=
             DDF*T.ax*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.cz;

          HessianMatrix.element[n+1][n+1]+=
             DDF*T.ay*T.ay+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.ay/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.ay)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
             +(dme.ay*dme.ay/SQR(m2)+dne.ay*dne.ay/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ay*T.ay+0.5*(S.ax+S.by);
          HessianMatrix.element[n+1][n+2]+=
             DDF*T.ay*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.ay*dme.az/SQR(m2)+dne.ay*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.z+vec_u.x*vec_u.y*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.z+vec_u.x*vec_u.y*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.az+0.5*S.bz;
          HessianMatrix.element[n+1][n+3]+=
             DDF*T.ay*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.ay*dme.by/SQR(m2)+dne.ay*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ay*T.by+S.ay;
          HessianMatrix.element[n+1][n+4]+=
             DDF*T.ay*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.bz)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.ay*dme.bz/SQR(m2)+dne.ay*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.z+vec_u.x*vec_u.y*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.z+vec_u.x*vec_u.y*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.bz+0.5*S.az;
          HessianMatrix.element[n+1][n+5]+=
             DDF*T.ay*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.ay*dme.cz/SQR(m2)+dne.ay*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.cz;

          HessianMatrix.element[n+2][n+2]+=DDF*T.az*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.az*dme.az/SQR(m2)+dne.az*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.az+0.5*(S.ax+S.cz);
          HessianMatrix.element[n+2][n+3]+=
             DDF*T.az*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.az*dme.by/SQR(m2)+dne.az*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.az*T.by;
          HessianMatrix.element[n+2][n+4]+=
             DDF*T.az*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.bz)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.az*dme.bz/SQR(m2)+dne.az*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.z+vec_u.x*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.z+vec_u.x*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.bz+0.5*S.ay;
          HessianMatrix.element[n+2][n+5]+=
             DDF*T.az*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.az*dme.cz/SQR(m2)+dne.az*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.cz+S.az;

          HessianMatrix.element[n+3][n+3]+=
             DDF*T.by*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
          HessianMatrix.element[n+3][n+4]+=
             DDF*T.by*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.bz)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.by*dme.bz/SQR(m2)+dne.by*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.by*T.bz+S.bz;
          HessianMatrix.element[n+3][n+5]+=
             DDF*T.by*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.by*T.cz;

          HessianMatrix.element[n+4][n+4]+=
             DDF*T.bz*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.bz*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.bz)+2.0*dot_mn*
             ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.bz*dme.bz/SQR(m2)+dne.bz*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.bz*T.bz+0.5*(S.by+S.cz);
          HessianMatrix.element[n+4][n+5]+=
             DDF*T.bz*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.bz*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.cz)+2.0*dot_mn*
             ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.bz*dme.cz/SQR(m2)+dne.bz*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.bz*T.cz+S.bz;

          HessianMatrix.element[n+5][n+5]+=
             DDF*T.cz*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
             ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
          break;
        case REGULAR_UPPER_TRIANGLE:
          HessianMatrix.element[n][n]+=
             DDF*T.ax*T.ax+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
             +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
             +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
          HessianMatrix.element[n][n+1]+=
             DDF*T.ax*T.ay+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.ay/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ay)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
             +(dme.ax*dme.ay/SQR(m2)+dne.ax*dne.ay/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ax*T.ay+S.ay;
          HessianMatrix.element[n][n+2]+=
             DDF*T.ax*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.ax*dme.az/SQR(m2)+dne.ax*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.az+S.az;
          HessianMatrix.element[n][n+3]+=
             DDF*T.ax*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ax*T.by;
          HessianMatrix.element[n][n+4]+=
             DDF*T.ax*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.bz)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.ax*dme.bz/SQR(m2)+dne.ax*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.bz;
          HessianMatrix.element[n][n+5]+=
             DDF*T.ax*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ax*T.cz;

          HessianMatrix.element[n+1][n+1]+=
             DDF*T.ay*T.ay+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.ay/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.ay)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
             +(dme.ay*dme.ay/SQR(m2)+dne.ay*dne.ay/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ay*T.ay+S.by;
          HessianMatrix.element[n+1][n+2]+=
             DDF*T.ay*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.ay*dme.az/SQR(m2)+dne.ay*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.z+vec_u.x*vec_u.y*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.z+vec_u.x*vec_u.y*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.az+S.bz;
          HessianMatrix.element[n+1][n+3]+=
             DDF*T.ay*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.ay*dme.by/SQR(m2)+dne.ay*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.ay*T.by;
          HessianMatrix.element[n+1][n+4]+=
             DDF*T.ay*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.bz)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.ay*dme.bz/SQR(m2)+dne.ay*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.z+vec_u.x*vec_u.y*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.z+vec_u.x*vec_u.y*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.bz;
          HessianMatrix.element[n+1][n+5]+=
             DDF*T.ay*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.ay*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.ay*dme.cz/SQR(m2)+dne.ay*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.ay*T.cz;

          HessianMatrix.element[n+2][n+2]+=DDF*T.az*T.az+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.az/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.az)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
             +(dme.az*dme.az/SQR(m2)+dne.az*dne.az/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.az+S.cz;
          HessianMatrix.element[n+2][n+3]+=
             DDF*T.az*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.by)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.az*dme.by/SQR(m2)+dne.az*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.az*T.by;
          HessianMatrix.element[n+2][n+4]+=
             DDF*T.az*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.bz)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.az*dme.bz/SQR(m2)+dne.az*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.z+vec_u.x*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.z+vec_u.x*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.bz;
          HessianMatrix.element[n+2][n+5]+=
             DDF*T.az*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.az*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.cz)+2.0*dot_mn*
             ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.az*dme.cz/SQR(m2)+dne.az*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.az*T.cz;

          HessianMatrix.element[n+3][n+3]+=
             DDF*T.by*T.by+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.by/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
             +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
             +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
          HessianMatrix.element[n+3][n+4]+=
             DDF*T.by*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.bz)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.by*dme.bz/SQR(m2)+dne.by*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.by*T.bz+S.bz;
          HessianMatrix.element[n+3][n+5]+=
             DDF*T.by*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.by*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
             ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.by*T.cz;

          HessianMatrix.element[n+4][n+4]+=
             DDF*T.bz*T.bz+
             DF*0.5*CosPhi*(
             -4.0*dmn.bz*dmn.bz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.bz)+2.0*dot_mn*
             ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
             -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
             +(dme.bz*dme.bz/SQR(m2)+dne.bz*dne.bz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
             +DF*(1.0/CosPhi)*T.bz*T.bz+S.cz;
          HessianMatrix.element[n+4][n+5]+=
             DDF*T.bz*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.bz*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.cz)+2.0*dot_mn*
             ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.bz*dme.cz/SQR(m2)+dne.bz*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.bz*T.cz;

          HessianMatrix.element[n+5][n+5]+=
             DDF*T.cz*T.cz+
             DF*0.5*CosPhi*(
             -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
             +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
             ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
             -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
             +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
             -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
             -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
             +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n]+=
                 DDF*T.ax*T.ax+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
                 +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
                 +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=
                 DDF*T.ax*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ax*T.by;
              HessianMatrix.element[n][n+2]+=
                 DDF*T.ax*T.bz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.bz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.bz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
                 +(dme.ax*dme.bz/SQR(m2)+dne.ax*dne.bz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.bz;
              HessianMatrix.element[n][n+3]+=
                 DDF*T.ax*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.cz;

              HessianMatrix.element[n+1][n+1]+=
                 DDF*T.by*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
              HessianMatrix.element[n+1][n+2]+=
                 DDF*T.by*T.bz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.bz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.bz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
                 +(dme.by*dme.bz/SQR(m2)+dne.by*dne.bz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.by*T.bz+S.bz;
              HessianMatrix.element[n+1][n+3]+=
                 DDF*T.by*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.by*T.cz;

              HessianMatrix.element[n+2][n+2]+=
                 DDF*T.bz*T.bz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.bz*dmn.bz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.bz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
                 -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
                 +(dme.bz*dme.bz/SQR(m2)+dne.bz*dne.bz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.bz*T.bz+0.5*(S.by+S.cz);
              HessianMatrix.element[n+2][n+3]+=
                 DDF*T.bz*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.bz*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.cz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.bz*dme.cz/SQR(m2)+dne.bz*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.bz*T.cz+S.bz;

              HessianMatrix.element[n+3][n+3]+=
                 DDF*T.cz*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
                 ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=
                 DDF*T.ax*T.ax+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
                 +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
                 +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=
                 DDF*T.ax*T.az+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.az/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.az)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
                 +(dme.ax*dme.az/SQR(m2)+dne.ax*dne.az/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.az+S.az;
              HessianMatrix.element[n][n+2]+=
                 DDF*T.ax*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ax*T.by;
              HessianMatrix.element[n][n+3]+=
                 DDF*T.ax*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.cz;

              HessianMatrix.element[n+1][n+1]+=DDF*T.az*T.az+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.az/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.az)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
                 +(dme.az*dme.az/SQR(m2)+dne.az*dne.az/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.az*T.az+0.5*(S.ax+S.cz);
              HessianMatrix.element[n+1][n+2]+=
                 DDF*T.az*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.az*dme.by/SQR(m2)+dne.az*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.az*T.by;
              HessianMatrix.element[n+1][n+3]+=
                 DDF*T.az*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.az*dme.cz/SQR(m2)+dne.az*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.az*T.cz+S.az;

              HessianMatrix.element[n+2][n+2]+=
                 DDF*T.by*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=
                 DDF*T.by*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.by*T.cz;

              HessianMatrix.element[n+3][n+3]+=
                 DDF*T.cz*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
                 ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=
                 DDF*T.ax*T.ax+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
                 +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
                 +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=
                 DDF*T.ax*T.ay+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.ay/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ay)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
                 +(dme.ax*dme.ay/SQR(m2)+dne.ax*dne.ay/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ax*T.ay+S.ay;
              HessianMatrix.element[n][n+2]+=
                 DDF*T.ax*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ax*T.by;
              HessianMatrix.element[n][n+3]+=
                 DDF*T.ax*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.cz;

              HessianMatrix.element[n+1][n+1]+=
                 DDF*T.ay*T.ay+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.ay/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.ay)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
                 +(dme.ay*dme.ay/SQR(m2)+dne.ay*dne.ay/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ay*T.ay+0.5*(S.ax+S.by);
              HessianMatrix.element[n+1][n+2]+=
                 DDF*T.ay*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.ay*dme.by/SQR(m2)+dne.ay*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ay*T.by+S.ay;
              HessianMatrix.element[n+1][n+3]+=
                 DDF*T.ay*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.ay*dme.cz/SQR(m2)+dne.ay*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ay*T.cz;

              HessianMatrix.element[n+2][n+2]+=
                 DDF*T.by*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=
                 DDF*T.by*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.by*T.cz;

              HessianMatrix.element[n+3][n+3]+=
                 DDF*T.cz*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
                 ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n]+=
                 DDF*T.ax*T.ax+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
                 +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
                 +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=
                 DDF*T.ax*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ax*T.by;
              HessianMatrix.element[n][n+2]+=
                 DDF*T.ax*T.bz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.bz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.bz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
                 +(dme.ax*dme.bz/SQR(m2)+dne.ax*dne.bz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_v.y*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.z+vec_u.x*vec_u.x*vec_w.y*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.bz;
              HessianMatrix.element[n][n+3]+=
                 DDF*T.ax*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.cz;

              HessianMatrix.element[n+1][n+1]+=
                 DDF*T.by*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
              HessianMatrix.element[n+1][n+2]+=
                 DDF*T.by*T.bz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.bz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.bz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
                 +(dme.by*dme.bz/SQR(m2)+dne.by*dne.bz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.z+vec_u.y*vec_u.y*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.by*T.bz+S.bz;
              HessianMatrix.element[n+1][n+3]+=
                 DDF*T.by*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.by*T.cz;

              HessianMatrix.element[n+2][n+2]+=
                 DDF*T.bz*T.bz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.bz*dmn.bz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.bz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)
                 -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.y*vec_u.z+vec_u.y*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.y*vec_w.z+vec_w.y*vec_v.z)))
                 +(dme.bz*dme.bz/SQR(m2)+dne.bz*dne.bz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_v.y*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.y*vec_u.z+vec_u.y*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.y*vec_u.z+vec_u.y*vec_u.z*vec_w.y*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.bz*T.bz+S.cz;
              HessianMatrix.element[n+2][n+3]+=
                 DDF*T.bz*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.bz*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.bz*dmn.cz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.y*vec_w.z+vec_w.y*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.z+vec_u.y*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.bz*dme.cz/SQR(m2)+dne.bz*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.z+vec_u.y*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.z*vec_u.z*vec_u.z+vec_u.y*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.z+vec_u.y*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.bz*T.cz;

              HessianMatrix.element[n+3][n+3]+=
                 DDF*T.cz*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
                 ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=
                 DDF*T.ax*T.ax+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
                 +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
                 +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=
                 DDF*T.ax*T.az+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.az/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.az)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
                 +(dme.ax*dme.az/SQR(m2)+dne.ax*dne.az/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.z+vec_u.x*vec_u.x*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.az+S.az;
              HessianMatrix.element[n][n+2]+=
                 DDF*T.ax*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ax*T.by;
              HessianMatrix.element[n][n+3]+=
                 DDF*T.ax*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.cz;

              HessianMatrix.element[n+1][n+1]+=DDF*T.az*T.az+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.az/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.az)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.x*vec_u.z+vec_u.x*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.x*vec_w.z+vec_w.x*vec_v.z)))
                 +(dme.az*dme.az/SQR(m2)+dne.az*dne.az/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_v.x*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.x*vec_u.z+vec_u.x*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.x*vec_u.z+vec_u.x*vec_u.z*vec_w.x*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.az*T.az+S.cz;
              HessianMatrix.element[n+1][n+2]+=
                 DDF*T.az*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.az*dme.by/SQR(m2)+dne.az*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.y*vec_u.y+vec_u.x*vec_u.z*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.az*T.by;
              HessianMatrix.element[n+1][n+3]+=
                 DDF*T.az*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.az*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.az*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.z+vec_w.x*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.z+vec_u.x*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.az*dme.cz/SQR(m2)+dne.az*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.z+vec_u.x*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.z*vec_u.z*vec_u.z+vec_u.x*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.z+vec_u.x*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.az*T.cz;

              HessianMatrix.element[n+2][n+2]+=
                 DDF*T.by*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=
                 DDF*T.by*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.by*T.cz;

              HessianMatrix.element[n+3][n+3]+=
                 DDF*T.cz*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
                 ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=
                 DDF*T.ax*T.ax+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.ax/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ax)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.x+vec_u.x*vec_u.x)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.x+vec_w.x*vec_v.x)))
                 +(dme.ax*dme.ax/SQR(m2)+dne.ax*dne.ax/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_v.x*vec_v.x)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.x+vec_u.x*vec_v.x))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.x+vec_u.x*vec_u.x*vec_w.x*vec_w.x)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)))
                 +DF*(1.0/CosPhi)*T.ax*T.ax+2.0*S.ax;
              HessianMatrix.element[n][n+1]+=
                 DDF*T.ax*T.ay+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.ay/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.ay)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
                 +(dme.ax*dme.ay/SQR(m2)+dne.ax*dne.ay/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.x*vec_u.y+vec_u.x*vec_u.x*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ax*T.ay+S.ay;
              HessianMatrix.element[n][n+2]+=
                 DDF*T.ax*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.ax*dme.by/SQR(m2)+dne.ax*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.y*vec_u.y+vec_u.x*vec_u.x*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ax*T.by;
              HessianMatrix.element[n][n+3]+=
                 DDF*T.ax*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ax*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ax*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.x+vec_w.x*vec_v.x)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.x+vec_u.x*vec_u.x)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.ax*dme.cz/SQR(m2)+dne.ax*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.x+vec_u.x*vec_v.x)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.x*vec_u.z*vec_u.z+vec_u.x*vec_u.x*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.x+vec_u.x*vec_w.x)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ax*T.cz;

              HessianMatrix.element[n+1][n+1]+=
                 DDF*T.ay*T.ay+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.ay/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.ay)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.x*vec_u.y+vec_u.x*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.x*vec_w.y+vec_w.x*vec_v.y)))
                 +(dme.ay*dme.ay/SQR(m2)+dne.ay*dne.ay/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_v.x*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.x*vec_u.y+vec_u.x*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.x*vec_u.y+vec_u.x*vec_u.y*vec_w.x*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ay*T.ay+S.by;
              HessianMatrix.element[n+1][n+2]+=
                 DDF*T.ay*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.by)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.ay*dme.by/SQR(m2)+dne.ay*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.y*vec_u.y+vec_u.x*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.ay*T.by;
              HessianMatrix.element[n+1][n+3]+=
                 DDF*T.ay*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.ay*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.ay*dmn.cz)+2.0*dot_mn*
                 ((vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.x*vec_w.y+vec_w.x*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.x*vec_u.y+vec_u.x*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.ay*dme.cz/SQR(m2)+dne.ay*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.x*vec_v.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.x*vec_u.y+vec_u.x*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.x*vec_w.y*vec_u.z*vec_u.z+vec_u.x*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.x*vec_u.y+vec_u.x*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.ay*T.cz;

              HessianMatrix.element[n+2][n+2]+=
                 DDF*T.by*T.by+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.by/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.by)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.y*vec_u.y+vec_u.y*vec_u.y)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.y*vec_w.y+vec_w.y*vec_v.y)))
                 +(dme.by*dme.by/SQR(m2)+dne.by*dne.by/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_v.y*vec_v.y)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.y*vec_u.y+vec_u.y*vec_v.y))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.y*vec_u.y+vec_u.y*vec_u.y*vec_w.y*vec_w.y)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)))
                 +DF*(1.0/CosPhi)*T.by*T.by+2.0*S.by;
              HessianMatrix.element[n+2][n+3]+=
                 DDF*T.by*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.by*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.by*dmn.cz)+2.0*dot_mn*
                 ((vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.y*vec_w.y+vec_w.y*vec_v.y)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.y*vec_u.y+vec_u.y*vec_u.y)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.by*dme.cz/SQR(m2)+dne.by*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.y*vec_v.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_v.z*vec_v.z)-2.0*(vec_v.y*vec_u.y+vec_u.y*vec_v.y)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.y*vec_w.y*vec_u.z*vec_u.z+vec_u.y*vec_u.y*vec_w.z*vec_w.z)-2.0*(vec_w.y*vec_u.y+vec_u.y*vec_w.y)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.by*T.cz;

              HessianMatrix.element[n+3][n+3]+=
                 DDF*T.cz*T.cz+
                 DF*0.5*CosPhi*(
                 -4.0*dmn.cz*dmn.cz/SQR(dot_mn)
                 +(1.0/SQR(dot_mn))*(2.0*(dmn.cz*dmn.cz)+2.0*dot_mn*
                 ((vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)+(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)
                 -(vec_v.z*vec_w.z+vec_w.z*vec_v.z)*(vec_u.z*vec_u.z+vec_u.z*vec_u.z)-(vec_u.z*vec_u.z+vec_u.z*vec_u.z)*(vec_v.z*vec_w.z+vec_w.z*vec_v.z)))
                 +(dme.cz*dme.cz/SQR(m2)+dne.cz*dne.cz/SQR(n2))
                 -(1.0/m2)*(4.0*(vec_v.z*vec_v.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_v.z*vec_v.z)-2.0*(vec_v.z*vec_u.z+vec_u.z*vec_v.z)*(vec_v.z*vec_u.z+vec_u.z*vec_v.z))
                 -(1.0/n2)*(4.0*(vec_w.z*vec_w.z*vec_u.z*vec_u.z+vec_u.z*vec_u.z*vec_w.z*vec_w.z)-2.0*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)*(vec_w.z*vec_u.z+vec_u.z*vec_w.z)))
                 +DF*(1.0/CosPhi)*T.cz*T.cz+2.0*S.cz;
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}


static inline void GradientStrainTorsion(REAL *Gradient,REAL_MATRIX3x3 S)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]+=S.ax+S.by+S.cz;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n]+=S.ax;
          Gradient[n+1]+=S.by;
          Gradient[n+2]+=S.cz;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n]+=S.ax;
          Gradient[n+1]+=S.bx;
          Gradient[n+2]+=S.cx;
          Gradient[n+3]+=S.by;
          Gradient[n+4]+=S.cy;
          Gradient[n+5]+=S.cz;
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.by;
              Gradient[n+2]+=S.cy;
              Gradient[n+3]+=S.cz;
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.cx;
              Gradient[n+2]+=S.by;
              Gradient[n+3]+=S.cz;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]+=S.ax;
              Gradient[n+1]+=S.bx;
              Gradient[n+2]+=S.by;
              Gradient[n+3]+=S.cz;
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}



void ComputeFrameworkTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL *parms;
  int index_i,index_j,index_k,index_l;
  REAL_MATRIX3x3 D2I,D2L,D2IL,D2IJ,D2IK,D2JL;
  REAL_MATRIX3x3 D2KL,D2J,D2K,D2JK,S;
  VECTOR DIL,DDJ,DDK,DEJ,DEK;
  REAL DF,DDF;
  VECTOR vec_u,vec_v,vec_w;
  REAL u,w,v;
  VECTOR fa,fb,fc,fd;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[f1];i++)
    {
      A=Framework[CurrentSystem].Torsions[f1][i].A;
      B=Framework[CurrentSystem].Torsions[f1][i].B;
      C=Framework[CurrentSystem].Torsions[f1][i].C;
      D=Framework[CurrentSystem].Torsions[f1][i].D;
      parms=Framework[CurrentSystem].TorsionArguments[f1][i];

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD=Framework[CurrentSystem].Atoms[f1][D].Position;

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;
      index_k=Framework[CurrentSystem].Atoms[f1][C].HessianIndex;
      index_l=Framework[CurrentSystem].Atoms[f1][D].HessianIndex;

      vec_u.x=posB.x-posC.x;
      vec_u.y=posB.y-posC.y;
      vec_u.z=posB.z-posC.z;
      vec_u=ApplyBoundaryCondition(vec_u);
      u=sqrt(SQR(vec_u.x)+SQR(vec_u.y)+SQR(vec_u.z));

      vec_v.x=posD.x-posC.x;
      vec_v.y=posD.y-posC.y;
      vec_v.z=posD.z-posC.z;
      vec_v=ApplyBoundaryCondition(vec_v);
      v=sqrt(SQR(vec_v.x)+SQR(vec_v.y)+SQR(vec_v.z));

      vec_w.x=posA.x-posB.x;
      vec_w.y=posA.y-posB.y;
      vec_w.z=posA.z-posB.z;
      vec_w=ApplyBoundaryCondition(vec_w);
      w=sqrt(SQR(vec_w.x)+SQR(vec_w.y)+SQR(vec_w.z));

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      Dab=ApplyBoundaryCondition(Dab);

      Dcb.x=posC.x-posB.x;
      Dcb.y=posC.y-posB.y;
      Dcb.z=posC.z-posB.z;
      Dcb=ApplyBoundaryCondition(Dcb);
      rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
      Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

      Ddc.x=posD.x-posC.x;
      Ddc.y=posD.y-posC.y;
      Ddc.z=posD.z-posC.z;
      Ddc=ApplyBoundaryCondition(Ddc);

      dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
      dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

      dr.x=Dab.x-dot_ab*Dcb.x;
      dr.y=Dab.y-dot_ab*Dcb.y;
      dr.z=Dab.z-dot_ab*Dcb.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Ddc.x-dot_cd*Dcb.x;
      ds.y=Ddc.y-dot_cd*Dcb.y;
      ds.z=Ddc.z-dot_cd*Dcb.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      // compute Cos(Phi)
      // Phi is defined in protein convention Phi(trans)=Pi
      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      // Ensure CosPhi is between -1 and 1.
      CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
      CosPhi2=SQR(CosPhi);

      switch(Framework[CurrentSystem].TorsionType[f1][i])
      {
        case HARMONIC_DIHEDRAL:
          // (1/2)*p_0*(phi-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // potential defined in terms of 'phi' and therefore contains a singularity
          // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
          // same direction as Rbc, and negative otherwise
          Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
          Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
          Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
          Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
          Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
          Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
          sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
          Phi=SIGN(acos(CosPhi),sign);
          SinPhi=sin(Phi);
          Phi-=parms[1];
          Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
          SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
          U=0.5*parms[0]*SQR(Phi);
          DF=-parms[0]*Phi/SinPhi;
          DDF=-parms[0]*(Phi*CosPhi-SinPhi)/CUBE(SinPhi);
          break;
        case HARMONIC_COSINE_DIHEDRAL:
          // (1/2)*p_0*(cos(phi)-cos(p_1))^2
          // ===============================================
          // p_0/k_B [K]
          // p_1     [degrees]
          U=0.5*parms[0]*SQR(CosPhi-parms[1]);
          DF=parms[0]*(CosPhi-parms[1]);
          DDF=parms[0];
          break;
        case THREE_COSINE_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
          DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[1]-6.0*parms[2]*CosPhi);
          break;
        case MM3_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0     [kcal/mol]
          // p_1     [kcal/mol]
          // p_2     [kcal/mol]
          U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
          DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[1]-6.0*parms[2]*CosPhi);
          break;
        case CFF_DIHEDRAL:
          // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
          DF=-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2);
          DDF=-4.0*(parms[1]+6.0*parms[2]*CosPhi);
          break;
        case CFF_DIHEDRAL2:
          // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
          DF=parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi);
          DDF=4.0*(parms[1]+6.0*parms[2]*CosPhi);
          break;
        case SIX_COSINE_DIHEDRAL:
          // Prod_i=0^5 p_i*cos(phi)^i
          // =========================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
          // between the first and last atoms of the dihedral, and Phi'=Phi-pi is defined accoording to the
          // polymer convention Phi'(trans)=0.
          U=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
                 parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
          DF=-parms[1]+2.0*parms[2]*CosPhi-3.0*parms[3]*CosPhi2+4.0*parms[4]*CosPhi2*CosPhi-5.0*parms[5]*SQR(CosPhi2);
          DDF=2.0*parms[2]-6.0*parms[3]*CosPhi+12.0*parms[4]*CosPhi2-20.0*parms[5]*CosPhi2*CosPhi;
          break;
        case TRAPPE_DIHEDRAL:
          // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
          // =============================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
          DF=parms[1]-4.0*parms[2]*CosPhi+3.0*parms[3]*(4.0*CosPhi2-1.0);
          DDF=-4.0*(parms[2]-6.0*parms[3]*CosPhi);
          break;
        case CVFF_DIHEDRAL:
          // p_0*(1+cos(p_1*phi-p_2))
          // ========================
          // p_0/k_B [K]
          // p_1     [-]
          // p_2     [degrees]
          // potential defined in terms of 'phi' and therefore contains a singularity
          // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
          // same direction as Dbc, and negative otherwise
          Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
          Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
          Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
          Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
          Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
          Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
          sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
          Phi=SIGN(acos(CosPhi),sign);
          SinPhi=sin(Phi);
          U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
          if(fabs(SinPhi)<1.0e-8)
          {
            SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
            DF=parms[0]*parms[1]*cos(parms[2])*SIGN(parms[1],sin(parms[1]*Phi)*Phi)+parms[0]*parms[1]*cos(parms[1]*Phi)*sin(parms[2])/SinPhi;
          }
          else
            DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
          DDF=-parms[0]*parms[1]*(parms[1]*cos(-parms[1]*Phi+parms[2])+sin(-parms[1]*Phi+parms[2])*CosPhi/SinPhi)/SQR(SinPhi);
          break;
        case OPLS_DIHEDRAL:
          // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
          // =================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
          DF=0.5*parms[1]-2.0*parms[2]*CosPhi+1.5*parms[3]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[2]-6.0*parms[3]*CosPhi);
          break;
        case FOURIER_SERIES_DIHEDRAL:
          // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
          // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
          // =======================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          U=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
                 2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi2+4.0*(parms[2]-5.0*parms[4])*CosPhi2*CosPhi-
                 8.0*(parms[3]-6.0*parms[5])*SQR(CosPhi2)+16.0*parms[4]*SQR(CosPhi2)*CosPhi-32.0*parms[5]*CUBE(CosPhi2));
          DF=0.5*(parms[0]-3.0*parms[2]+5.0*parms[4])-2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi+
             6.0*(parms[2]-5.0*parms[4])*CosPhi2-16.0*(parms[3]-6.0*parms[5])*CosPhi2*CosPhi+40.0*parms[4]*SQR(CosPhi2)-96.0*parms[5]*CosPhi2*CUBE(CosPhi);
          DDF=-2.0*parms[1]+8.0*parms[3]-18.0*parms[5]+12.0*(parms[2]-5.0*parms[4])*CosPhi
              -48.0*(parms[3]-6.0*parms[5])*CosPhi2+160.0*parms[4]*CosPhi2*CosPhi-480.0*parms[5]*SQR(CosPhi2);
          break;
        case FOURIER_SERIES_DIHEDRAL2:
          // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
          // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
          // =======================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          U=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
            2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
            2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
          DF=0.5*parms[0]+parms[2]*(6.0*CosPhi2-1.5)+parms[4]*(2.5-30.0*CosPhi2+40.0*SQR(CosPhi2))+
             CosPhi*(-2.0*parms[1]+parms[3]*(16.0*CosPhi2-8.0)+parms[5]*(18.0-96.0*CosPhi2+96.0*SQR(CosPhi2)));
          DDF=-2.0*parms[1]+12.0*parms[2]*CosPhi+parms[3]*(48.0*CosPhi2-8.0)+parms[4]*CosPhi*(160.0*CosPhi2-60.0)
              +parms[5]*(18.0-288.0*CosPhi2+480.0*SQR(CosPhi2));
          break;
        default:
          printf("Undefined Torsion potential in 'framework_hessian.c'\n");
          U=DF=DDF=0.0;
          exit(0);
          break;
      }

      *Energy+=U;

      if((index_i<0)&&(index_j<0)&&(index_k<0)&&(index_l<0)) continue;

     // Calculate the first derivative vectors.
      d=dot_ab/rbc;
      e=dot_cd/rbc;

      dtA.x=(ds.x-CosPhi*dr.x)/r;
      dtA.y=(ds.y-CosPhi*dr.y)/r;
      dtA.z=(ds.z-CosPhi*dr.z)/r;

      dtD.x=(dr.x-CosPhi*ds.x)/s;
      dtD.y=(dr.y-CosPhi*ds.y)/s;
      dtD.z=(dr.z-CosPhi*ds.z)/s;

      dtB.x=dtA.x*(d-1.0)+e*dtD.x;
      dtB.y=dtA.y*(d-1.0)+e*dtD.y;
      dtB.z=dtA.z*(d-1.0)+e*dtD.z;

      dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
      dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
      dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

      // forces are oppositely directed to the gradient
      fa.x=DF*dtA.x;
      fa.y=DF*dtA.y;
      fa.z=DF*dtA.z;

      fb.x=DF*dtB.x;
      fb.y=DF*dtB.y;
      fb.z=DF*dtB.z;

      fc.x=DF*dtC.x;
      fc.y=DF*dtC.y;
      fc.z=DF*dtC.z;

      fd.x=DF*dtD.x;
      fd.y=DF*dtD.y;
      fd.z=DF*dtD.z;

      // add contribution to the strain derivative tensor
      // Note: rbc is here because the vector was normalized before
      S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
      S.bx=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
      S.cx=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

      S.ay=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
      S.by=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
      S.cy=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

      S.az=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
      S.bz=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
      S.cz=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;

      StrainDerivative->ax+=S.ax;
      StrainDerivative->bx+=S.bx;
      StrainDerivative->cx+=S.cx;

      StrainDerivative->ay+=S.ay;
      StrainDerivative->by+=S.by;
      StrainDerivative->cy+=S.cy;

      StrainDerivative->az+=S.az;
      StrainDerivative->bz+=S.bz;
      StrainDerivative->cz+=S.cz;

      if(ComputeGradient)
      {
        if(index_i>=0)
        {
          Gradient[index_i]+=fa.x;
          Gradient[index_i+1]+=fa.y;
          Gradient[index_i+2]+=fa.z;
        }

        if(index_j>=0)
        {
          Gradient[index_j]+=fb.x;
          Gradient[index_j+1]+=fb.y;
          Gradient[index_j+2]+=fb.z;
        }

        if(index_k>=0)
        {
          Gradient[index_k]+=fc.x;
          Gradient[index_k+1]+=fc.y;
          Gradient[index_k+2]+=fc.z;
        }

        if(index_l>=0)
        {
          Gradient[index_l]+=fd.x;
          Gradient[index_l+1]+=fd.y;
          Gradient[index_l+2]+=fd.z;
        }

        GradientStrainTorsion(Gradient,S);
      }

      if(ComputeHessian)
      {
        // Calculate the derivatives of DOTIJ and DOTLK.
        DIL.x=DF*Dcb.x/rbc;
        DIL.y=DF*Dcb.y/rbc;
        DIL.z=DF*Dcb.z/rbc;
        DDJ.x=DF*((2.0*d-1.0)*Dcb.x-Dab.x/rbc)/rbc;
        DDJ.y=DF*((2.0*d-1.0)*Dcb.y-Dab.y/rbc)/rbc;
        DDJ.z=DF*((2.0*d-1.0)*Dcb.z-Dab.z/rbc)/rbc;
        DDK.x=-DF*(2.0*d*Dcb.x-Dab.x/rbc)/rbc;
        DDK.y=-DF*(2.0*d*Dcb.y-Dab.y/rbc)/rbc;
        DDK.z=-DF*(2.0*d*Dcb.z-Dab.z/rbc)/rbc;
        DEJ.x=DF*(2.0*e*Dcb.x-Ddc.x/rbc)/rbc;
        DEJ.y=DF*(2.0*e*Dcb.y-Ddc.y/rbc)/rbc;
        DEJ.z=DF*(2.0*e*Dcb.z-Ddc.z/rbc)/rbc;
        DEK.x=-DF*((2.0*e+1.0)*Dcb.x-Ddc.x/rbc)/rbc;
        DEK.y=-DF*((2.0*e+1.0)*Dcb.y-Ddc.y/rbc)/rbc;
        DEK.z=-DF*((2.0*e+1.0)*Dcb.z-Ddc.z/rbc)/rbc;

        // Calculate some diagonal Hessian terms for I.
        D2I.ax=DF*(CosPhi*(dr.x*dr.x+Dcb.x*Dcb.x-1.0)/r-2.0*dr.x*dtA.x)/r;
        D2I.by=DF*(CosPhi*(dr.y*dr.y+Dcb.y*Dcb.y-1.0)/r-2.0*dr.y*dtA.y)/r;
        D2I.cz=DF*(CosPhi*(dr.z*dr.z+Dcb.z*Dcb.z-1.0)/r-2.0*dr.z*dtA.z)/r;
        D2I.ay=DF*(CosPhi*(dr.x*dr.y+Dcb.x*Dcb.y)/r-dr.x*dtA.y-dr.y*dtA.x)/r;
        D2I.az=DF*(CosPhi*(dr.x*dr.z+Dcb.x*Dcb.z)/r-dr.x*dtA.z-dr.z*dtA.x)/r;
        D2I.bz=DF*(CosPhi*(dr.y*dr.z+Dcb.y*Dcb.z)/r-dr.y*dtA.z-dr.z*dtA.y)/r;

        // Calculate some diagonal Hessian terms for L.
        D2L.ax=DF*(CosPhi*(ds.x*ds.x+Dcb.x*Dcb.x-1.0)/s-2.0*ds.x*dtD.x)/s;
        D2L.by=DF*(CosPhi*(ds.y*ds.y+Dcb.y*Dcb.y-1.0)/s-2.0*ds.y*dtD.y)/s;
        D2L.cz=DF*(CosPhi*(ds.z*ds.z+Dcb.z*Dcb.z-1.0)/s-2.0*ds.z*dtD.z)/s;
        D2L.ay=DF*(CosPhi*(ds.x*ds.y+Dcb.x*Dcb.y)/s-ds.x*dtD.y-ds.y*dtD.x)/s;
        D2L.az=DF*(CosPhi*(ds.x*ds.z+Dcb.x*Dcb.z)/s-ds.x*dtD.z-ds.z*dtD.x)/s;
        D2L.bz=DF*(CosPhi*(ds.y*ds.z+Dcb.y*Dcb.z)/s-ds.y*dtD.z-ds.z*dtD.y)/s;

        // Calculate the IL off-diagonal terms.
        D2IL.ax=DF*(CosPhi*dr.x*ds.x-dr.x*dr.x-ds.x*ds.x-Dcb.x*Dcb.x+1.0)/(r*s);
        D2IL.ay=DF*(CosPhi*dr.x*ds.y-dr.x*dr.y-ds.x*ds.y-Dcb.x*Dcb.y)/(r*s);
        D2IL.az=DF*(CosPhi*dr.x*ds.z-dr.x*dr.z-ds.x*ds.z-Dcb.x*Dcb.z)/(r*s);
        D2IL.bx=DF*(CosPhi*dr.y*ds.x-dr.y*dr.x-ds.y*ds.x-Dcb.y*Dcb.x)/(r*s);
        D2IL.by=DF*(CosPhi*dr.y*ds.y-dr.y*dr.y-ds.y*ds.y-Dcb.y*Dcb.y+1.0)/(r*s);
        D2IL.bz=DF*(CosPhi*dr.y*ds.z-dr.y*dr.z-ds.y*ds.z-Dcb.y*Dcb.z)/(r*s);
        D2IL.cx=DF*(CosPhi*dr.z*ds.x-dr.z*dr.x-ds.z*ds.x-Dcb.z*Dcb.x)/(r*s);
        D2IL.cy=DF*(CosPhi*dr.z*ds.y-dr.z*dr.y-ds.z*ds.y-Dcb.z*Dcb.y)/(r*s);
        D2IL.cz=DF*(CosPhi*dr.z*ds.z-dr.z*dr.z-ds.z*ds.z-Dcb.z*Dcb.z+1.0)/(r*s);

        // Calculate the IJ off-diagonal terms.
        D2IJ.ax=D2I.ax*(d-1.0)+D2IL.ax*e+DIL.x*dtA.x;
        D2IJ.ay=D2I.ay*(d-1.0)+D2IL.ay*e+DIL.x*dtA.y;
        D2IJ.az=D2I.az*(d-1.0)+D2IL.az*e+DIL.x*dtA.z;
        D2IJ.bx=D2I.ay*(d-1.0)+D2IL.bx*e+DIL.y*dtA.x;
        D2IJ.by=D2I.by*(d-1.0)+D2IL.by*e+DIL.y*dtA.y;
        D2IJ.bz=D2I.bz*(d-1.0)+D2IL.bz*e+DIL.y*dtA.z;
        D2IJ.cx=D2I.az*(d-1.0)+D2IL.cx*e+DIL.z*dtA.x;
        D2IJ.cy=D2I.bz*(d-1.0)+D2IL.cy*e+DIL.z*dtA.y;
        D2IJ.cz=D2I.cz*(d-1.0)+D2IL.cz*e+DIL.z*dtA.z;

        // Calculate the IK off-diagonal terms.
        D2IK.ax=-D2I.ax*d-D2IL.ax*(e+1.0)-DIL.x*dtA.x;
        D2IK.ay=-D2I.ay*d-D2IL.ay*(e+1.0)-DIL.x*dtA.y;
        D2IK.az=-D2I.az*d-D2IL.az*(e+1.0)-DIL.x*dtA.z;
        D2IK.bx=-D2I.ay*d-D2IL.bx*(e+1.0)-DIL.y*dtA.x;
        D2IK.by=-D2I.by*d-D2IL.by*(e+1.0)-DIL.y*dtA.y;
        D2IK.bz=-D2I.bz*d-D2IL.bz*(e+1.0)-DIL.y*dtA.z;
        D2IK.cx=-D2I.az*d-D2IL.cx*(e+1.0)-DIL.z*dtA.x;
        D2IK.cy=-D2I.bz*d-D2IL.cy*(e+1.0)-DIL.z*dtA.y;
        D2IK.cz=-D2I.cz*d-D2IL.cz*(e+1.0)-DIL.z*dtA.z;

        // Calculate the JL off-diagonal terms.
        D2JL.ax=D2IL.ax*(d-1.0)+D2L.ax*e+DIL.x*dtD.x;
        D2JL.ay=D2IL.ay*(d-1.0)+D2L.ay*e+DIL.y*dtD.x;
        D2JL.az=D2IL.az*(d-1.0)+D2L.az*e+DIL.z*dtD.x;
        D2JL.bx=D2IL.bx*(d-1.0)+D2L.ay*e+DIL.x*dtD.y;
        D2JL.by=D2IL.by*(d-1.0)+D2L.by*e+DIL.y*dtD.y;
        D2JL.bz=D2IL.bz*(d-1.0)+D2L.bz*e+DIL.z*dtD.y;
        D2JL.cx=D2IL.cx*(d-1.0)+D2L.az*e+DIL.x*dtD.z;
        D2JL.cy=D2IL.cy*(d-1.0)+D2L.bz*e+DIL.y*dtD.z;
        D2JL.cz=D2IL.cz*(d-1.0)+D2L.cz*e+DIL.z*dtD.z;

        // Calculate the KL off-diagonal terms.
        D2KL.ax=-D2IL.ax*d-D2L.ax*(e+1.0)-DIL.x*dtD.x;
        D2KL.ay=-D2IL.ay*d-D2L.ay*(e+1.0)-DIL.y*dtD.x;
        D2KL.az=-D2IL.az*d-D2L.az*(e+1.0)-DIL.z*dtD.x;
        D2KL.bx=-D2IL.bx*d-D2L.ay*(e+1.0)-DIL.x*dtD.y;
        D2KL.by=-D2IL.by*d-D2L.by*(e+1.0)-DIL.y*dtD.y;
        D2KL.bz=-D2IL.bz*d-D2L.bz*(e+1.0)-DIL.z*dtD.y;
        D2KL.cx=-D2IL.cx*d-D2L.az*(e+1.0)-DIL.x*dtD.z;
        D2KL.cy=-D2IL.cy*d-D2L.bz*(e+1.0)-DIL.y*dtD.z;
        D2KL.cz=-D2IL.cz*d-D2L.cz*(e+1.0)-DIL.z*dtD.z;

        // Calculate the JJ diagonal terms.
        D2J.ax=D2IJ.ax*(d-1.0)+D2JL.ax*e+DDJ.x*dtA.x+DEJ.x*dtD.x;
        D2J.by=D2IJ.by*(d-1.0)+D2JL.by*e+DDJ.y*dtA.y+DEJ.y*dtD.y;
        D2J.cz=D2IJ.cz*(d-1.0)+D2JL.cz*e+DDJ.z*dtA.z+DEJ.z*dtD.z;
        D2J.ay=D2IJ.ay*(d-1.0)+D2JL.bx*e+DDJ.y*dtA.x+DEJ.y*dtD.x;
        D2J.az=D2IJ.az*(d-1.0)+D2JL.cx*e+DDJ.z*dtA.x+DEJ.z*dtD.x;
        D2J.bz=D2IJ.bz*(d-1.0)+D2JL.cy*e+DDJ.z*dtA.y+DEJ.z*dtD.y;

        // Calculate the KK diagonal terms.
        D2K.ax=-D2KL.ax*(e+1.0)-D2IK.ax*d-DDK.x*dtA.x-DEK.x*dtD.x;
        D2K.by=-D2KL.by*(e+1.0)-D2IK.by*d-DDK.y*dtA.y-DEK.y*dtD.y;
        D2K.cz=-D2KL.cz*(e+1.0)-D2IK.cz*d-DDK.z*dtA.z-DEK.z*dtD.z;
        D2K.ay=-D2KL.ay*(e+1.0)-D2IK.bx*d-DDK.x*dtA.y-DEK.x*dtD.y;
        D2K.az=-D2KL.az*(e+1.0)-D2IK.cx*d-DDK.x*dtA.z-DEK.x*dtD.z;
        D2K.bz=-D2KL.bz*(e+1.0)-D2IK.cy*d-DDK.y*dtA.z-DEK.y*dtD.z;

        // Calculate the JK off-diagonal terms.
        D2JK.ax=-D2IJ.ax*d-D2JL.ax*(e+1.0)-DDJ.x*dtA.x-DEJ.x*dtD.x;
        D2JK.ay=-D2IJ.bx*d-D2JL.ay*(e+1.0)-DDJ.x*dtA.y-DEJ.x*dtD.y;
        D2JK.az=-D2IJ.cx*d-D2JL.az*(e+1.0)-DDJ.x*dtA.z-DEJ.x*dtD.z;
        D2JK.bx=-D2IJ.ay*d-D2JL.bx*(e+1.0)-DDJ.y*dtA.x-DEJ.y*dtD.x;
        D2JK.by=-D2IJ.by*d-D2JL.by*(e+1.0)-DDJ.y*dtA.y-DEJ.y*dtD.y;
        D2JK.bz=-D2IJ.cy*d-D2JL.bz*(e+1.0)-DDJ.y*dtA.z-DEJ.y*dtD.z;
        D2JK.cx=-D2IJ.az*d-D2JL.cx*(e+1.0)-DDJ.z*dtA.x-DEJ.z*dtD.x;
        D2JK.cy=-D2IJ.bz*d-D2JL.cy*(e+1.0)-DDJ.z*dtA.y-DEJ.z*dtD.y;
        D2JK.cz=-D2IJ.cz*d-D2JL.cz*(e+1.0)-DDJ.z*dtA.z-DEJ.z*dtD.z;

        // Calculate the diagonal blocks of the hessian.
        if(index_i>=0)
        {
          HessianMatrix.element[index_i][index_i]+=DDF*dtA.x*dtA.x+D2I.ax;
          HessianMatrix.element[index_i][index_i+1]+=DDF*dtA.x*dtA.y+D2I.ay;
          HessianMatrix.element[index_i+1][index_i+1]+=DDF*dtA.y*dtA.y+D2I.by;
          HessianMatrix.element[index_i][index_i+2]+=DDF*dtA.x*dtA.z+D2I.az;
          HessianMatrix.element[index_i+1][index_i+2]+=DDF*dtA.y*dtA.z+D2I.bz;
          HessianMatrix.element[index_i+2][index_i+2]+=DDF*dtA.z*dtA.z+D2I.cz;
        }

        if(index_j>=0)
        {
          HessianMatrix.element[index_j][index_j]+=DDF*dtB.x*dtB.x+D2J.ax;
          HessianMatrix.element[index_j][index_j+1]+=DDF*dtB.x*dtB.y+D2J.ay;
          HessianMatrix.element[index_j+1][index_j+1]+=DDF*dtB.y*dtB.y+D2J.by;
          HessianMatrix.element[index_j][index_j+2]+=DDF*dtB.x*dtB.z+D2J.az;
          HessianMatrix.element[index_j+1][index_j+2]+=DDF*dtB.y*dtB.z+D2J.bz;
          HessianMatrix.element[index_j+2][index_j+2]+=DDF*dtB.z*dtB.z+D2J.cz;
        }

        if(index_k>=0)
        {
          HessianMatrix.element[index_k][index_k]+=DDF*dtC.x*dtC.x+D2K.ax;
          HessianMatrix.element[index_k][index_k+1]+=DDF*dtC.x*dtC.y+D2K.ay;
          HessianMatrix.element[index_k+1][index_k+1]+=DDF*dtC.y*dtC.y+D2K.by;
          HessianMatrix.element[index_k][index_k+2]+=DDF*dtC.x*dtC.z+D2K.az;
          HessianMatrix.element[index_k+1][index_k+2]+=DDF*dtC.y*dtC.z+D2K.bz;
          HessianMatrix.element[index_k+2][index_k+2]+=DDF*dtC.z*dtC.z+D2K.cz;
        }

        if(index_l>=0)
        {
          HessianMatrix.element[index_l][index_l]+=DDF*dtD.x*dtD.x+D2L.ax;
          HessianMatrix.element[index_l][index_l+1]+=DDF*dtD.x*dtD.y+D2L.ay;
          HessianMatrix.element[index_l+1][index_l+1]+=DDF*dtD.y*dtD.y+D2L.by;
          HessianMatrix.element[index_l][index_l+2]+=DDF*dtD.x*dtD.z+D2L.az;
          HessianMatrix.element[index_l+1][index_l+2]+=DDF*dtD.y*dtD.z+D2L.bz;
          HessianMatrix.element[index_l+2][index_l+2]+=DDF*dtD.z*dtD.z+D2L.cz;
        }

        if((index_i>=0)&&(index_j>=0))
        {
          if(index_i<index_j)
          {
            HessianMatrix.element[index_i][index_j]+=DDF*dtA.x*dtB.x+D2IJ.ax;
            HessianMatrix.element[index_i][index_j+1]+=DDF*dtA.x*dtB.y+D2IJ.ay;
            HessianMatrix.element[index_i][index_j+2]+=DDF*dtA.x*dtB.z+D2IJ.az;
            HessianMatrix.element[index_i+1][index_j]+=DDF*dtA.y*dtB.x+D2IJ.bx;
            HessianMatrix.element[index_i+1][index_j+1]+=DDF*dtA.y*dtB.y+D2IJ.by;
            HessianMatrix.element[index_i+1][index_j+2]+=DDF*dtA.y*dtB.z+D2IJ.bz;
            HessianMatrix.element[index_i+2][index_j]+=DDF*dtA.z*dtB.x+D2IJ.cx;
            HessianMatrix.element[index_i+2][index_j+1]+=DDF*dtA.z*dtB.y+D2IJ.cy;
            HessianMatrix.element[index_i+2][index_j+2]+=DDF*dtA.z*dtB.z+D2IJ.cz;
          }
          else
          {
            HessianMatrix.element[index_j][index_i]+=DDF*dtA.x*dtB.x+D2IJ.ax;
            HessianMatrix.element[index_j][index_i+1]+=DDF*dtA.y*dtB.x+D2IJ.bx;
            HessianMatrix.element[index_j][index_i+2]+=DDF*dtA.z*dtB.x+D2IJ.cx;
            HessianMatrix.element[index_j+1][index_i]+=DDF*dtA.x*dtB.y+D2IJ.ay;
            HessianMatrix.element[index_j+1][index_i+1]+=DDF*dtA.y*dtB.y+D2IJ.by;
            HessianMatrix.element[index_j+1][index_i+2]+=DDF*dtA.z*dtB.y+D2IJ.cy;
            HessianMatrix.element[index_j+2][index_i]+=DDF*dtA.x*dtB.z+D2IJ.az;
            HessianMatrix.element[index_j+2][index_i+1]+=DDF*dtA.y*dtB.z+D2IJ.bz;
            HessianMatrix.element[index_j+2][index_i+2]+=DDF*dtA.z*dtB.z+D2IJ.cz;
          }
        }

        if((index_i>=0)&&(index_k>=0))
        {
          if(index_i<index_k)
          {
            HessianMatrix.element[index_i][index_k]+=DDF*dtA.x*dtC.x+D2IK.ax;
            HessianMatrix.element[index_i][index_k+1]+=DDF*dtA.x*dtC.y+D2IK.ay;
            HessianMatrix.element[index_i][index_k+2]+=DDF*dtA.x*dtC.z+D2IK.az;
            HessianMatrix.element[index_i+1][index_k]+=DDF*dtA.y*dtC.x+D2IK.bx;
            HessianMatrix.element[index_i+1][index_k+1]+=DDF*dtA.y*dtC.y+D2IK.by;
            HessianMatrix.element[index_i+1][index_k+2]+=DDF*dtA.y*dtC.z+D2IK.bz;
            HessianMatrix.element[index_i+2][index_k]+=DDF*dtA.z*dtC.x+D2IK.cx;
            HessianMatrix.element[index_i+2][index_k+1]+=DDF*dtA.z*dtC.y+D2IK.cy;
            HessianMatrix.element[index_i+2][index_k+2]+=DDF*dtA.z*dtC.z+D2IK.cz;
          }
          else
          {
            HessianMatrix.element[index_k][index_i]+=DDF*dtA.x*dtC.x+D2IK.ax;
            HessianMatrix.element[index_k][index_i+1]+=DDF*dtA.y*dtC.x+D2IK.bx;
            HessianMatrix.element[index_k][index_i+2]+=DDF*dtA.z*dtC.x+D2IK.cx;
            HessianMatrix.element[index_k+1][index_i]+=DDF*dtA.x*dtC.y+D2IK.ay;
            HessianMatrix.element[index_k+1][index_i+1]+=DDF*dtA.y*dtC.y+D2IK.by;
            HessianMatrix.element[index_k+1][index_i+2]+=DDF*dtA.z*dtC.y+D2IK.cy;
            HessianMatrix.element[index_k+2][index_i]+=DDF*dtA.x*dtC.z+D2IK.az;
            HessianMatrix.element[index_k+2][index_i+1]+=DDF*dtA.y*dtC.z+D2IK.bz;
            HessianMatrix.element[index_k+2][index_i+2]+=DDF*dtA.z*dtC.z+D2IK.cz;
          }
        }

        if((index_i>=0)&&(index_l>=0))
        {
          if(index_i<index_l)
          {
            HessianMatrix.element[index_i][index_l]+=DDF*dtA.x*dtD.x+D2IL.ax;
            HessianMatrix.element[index_i][index_l+1]+=DDF*dtA.x*dtD.y+D2IL.ay; 
            HessianMatrix.element[index_i][index_l+2]+=DDF*dtA.x*dtD.z+D2IL.az; 
            HessianMatrix.element[index_i+1][index_l]+=DDF*dtA.y*dtD.x+D2IL.bx; 
            HessianMatrix.element[index_i+1][index_l+1]+=DDF*dtA.y*dtD.y+D2IL.by;
            HessianMatrix.element[index_i+1][index_l+2]+=DDF*dtA.y*dtD.z+D2IL.bz;
            HessianMatrix.element[index_i+2][index_l]+=DDF*dtA.z*dtD.x+D2IL.cx;
            HessianMatrix.element[index_i+2][index_l+1]+=DDF*dtA.z*dtD.y+D2IL.cy;
            HessianMatrix.element[index_i+2][index_l+2]+=DDF*dtA.z*dtD.z+D2IL.cz;
          }
          else
          {
            HessianMatrix.element[index_l][index_i]+=DDF*dtA.x*dtD.x+D2IL.ax;
            HessianMatrix.element[index_l][index_i+1]+=DDF*dtA.y*dtD.x+D2IL.bx;
            HessianMatrix.element[index_l][index_i+2]+=DDF*dtA.z*dtD.x+D2IL.cx;
            HessianMatrix.element[index_l+1][index_i]+=DDF*dtA.x*dtD.y+D2IL.ay;
            HessianMatrix.element[index_l+1][index_i+1]+=DDF*dtA.y*dtD.y+D2IL.by;
            HessianMatrix.element[index_l+1][index_i+2]+=DDF*dtA.z*dtD.y+D2IL.cy;
            HessianMatrix.element[index_l+2][index_i]+=DDF*dtA.x*dtD.z+D2IL.az;
            HessianMatrix.element[index_l+2][index_i+1]+=DDF*dtA.y*dtD.z+D2IL.bz;
            HessianMatrix.element[index_l+2][index_i+2]+=DDF*dtA.z*dtD.z+D2IL.cz;
          }
        }

        if((index_j>=0)&&(index_k>=0))
        {
          if(index_j<index_k)
          {
            HessianMatrix.element[index_j][index_k]+=DDF*dtB.x*dtC.x+D2JK.ax;
            HessianMatrix.element[index_j][index_k+1]+=DDF*dtB.x*dtC.y+D2JK.ay; 
            HessianMatrix.element[index_j][index_k+2]+=DDF*dtB.x*dtC.z+D2JK.az;
            HessianMatrix.element[index_j+1][index_k]+=DDF*dtB.y*dtC.x+D2JK.bx;
            HessianMatrix.element[index_j+1][index_k+1]+=DDF*dtB.y*dtC.y+D2JK.by;
            HessianMatrix.element[index_j+1][index_k+2]+=DDF*dtB.y*dtC.z+D2JK.bz;
            HessianMatrix.element[index_j+2][index_k]+=DDF*dtB.z*dtC.x+D2JK.cx;
            HessianMatrix.element[index_j+2][index_k+1]+=DDF*dtB.z*dtC.y+D2JK.cy;
            HessianMatrix.element[index_j+2][index_k+2]+=DDF*dtB.z*dtC.z+D2JK.cz;
          }
          else
          {
            HessianMatrix.element[index_k][index_j]+=DDF*dtB.x*dtC.x+D2JK.ax;
            HessianMatrix.element[index_k][index_j+1]+=DDF*dtB.y*dtC.x+D2JK.bx;
            HessianMatrix.element[index_k][index_j+2]+=DDF*dtB.z*dtC.x+D2JK.cx;
            HessianMatrix.element[index_k+1][index_j]+=DDF*dtB.x*dtC.y+D2JK.ay;
            HessianMatrix.element[index_k+1][index_j+1]+=DDF*dtB.y*dtC.y+D2JK.by;
            HessianMatrix.element[index_k+1][index_j+2]+=DDF*dtB.z*dtC.y+D2JK.cy;
            HessianMatrix.element[index_k+2][index_j]+=DDF*dtB.x*dtC.z+D2JK.az;
            HessianMatrix.element[index_k+2][index_j+1]+=DDF*dtB.y*dtC.z+D2JK.bz;
            HessianMatrix.element[index_k+2][index_j+2]+=DDF*dtB.z*dtC.z+D2JK.cz;
          }
        }

        if((index_j>=0)&&(index_l>=0))
        {
          if(index_j<index_l)
          {
            HessianMatrix.element[index_j][index_l]+=DDF*dtB.x*dtD.x+D2JL.ax;
            HessianMatrix.element[index_j][index_l+1]+=DDF*dtB.x*dtD.y+D2JL.ay;
            HessianMatrix.element[index_j][index_l+2]+=DDF*dtB.x*dtD.z+D2JL.az;
            HessianMatrix.element[index_j+1][index_l]+=DDF*dtB.y*dtD.x+D2JL.bx;
            HessianMatrix.element[index_j+1][index_l+1]+=DDF*dtB.y*dtD.y+D2JL.by;
            HessianMatrix.element[index_j+1][index_l+2]+=DDF*dtB.y*dtD.z+D2JL.bz;
            HessianMatrix.element[index_j+2][index_l]+=DDF*dtB.z*dtD.x+D2JL.cx;
            HessianMatrix.element[index_j+2][index_l+1]+=DDF*dtB.z*dtD.y+D2JL.cy;
            HessianMatrix.element[index_j+2][index_l+2]+=DDF*dtB.z*dtD.z+D2JL.cz;
          }
          else
          {
            HessianMatrix.element[index_l][index_j]+=DDF*dtB.x*dtD.x+D2JL.ax;
            HessianMatrix.element[index_l][index_j+1]+=DDF*dtB.y*dtD.x+D2JL.bx;
            HessianMatrix.element[index_l][index_j+2]+=DDF*dtB.z*dtD.x+D2JL.cx;
            HessianMatrix.element[index_l+1][index_j]+=DDF*dtB.x*dtD.y+D2JL.ay;
            HessianMatrix.element[index_l+1][index_j+1]+=DDF*dtB.y*dtD.y+D2JL.by;
            HessianMatrix.element[index_l+1][index_j+2]+=DDF*dtB.z*dtD.y+D2JL.cy;
            HessianMatrix.element[index_l+2][index_j]+=DDF*dtB.x*dtD.z+D2JL.az;
            HessianMatrix.element[index_l+2][index_j+1]+=DDF*dtB.y*dtD.z+D2JL.bz;
            HessianMatrix.element[index_l+2][index_j+2]+=DDF*dtB.z*dtD.z+D2JL.cz;
          }
        }

        if((index_k>=0)&&(index_l>=0))
        {
          if(index_k<index_l)
          {
            HessianMatrix.element[index_k][index_l]+=DDF*dtC.x*dtD.x+D2KL.ax;
            HessianMatrix.element[index_k][index_l+1]+=DDF*dtC.x*dtD.y+D2KL.ay;
            HessianMatrix.element[index_k][index_l+2]+=DDF*dtC.x*dtD.z+D2KL.az;
            HessianMatrix.element[index_k+1][index_l]+=DDF*dtC.y*dtD.x+D2KL.bx;
            HessianMatrix.element[index_k+1][index_l+1]+=DDF*dtC.y*dtD.y+D2KL.by;
            HessianMatrix.element[index_k+1][index_l+2]+=DDF*dtC.y*dtD.z+D2KL.bz;
            HessianMatrix.element[index_k+2][index_l]+=DDF*dtC.z*dtD.x+D2KL.cx;
            HessianMatrix.element[index_k+2][index_l+1]+=DDF*dtC.z*dtD.y+D2KL.cy;
            HessianMatrix.element[index_k+2][index_l+2]+=DDF*dtC.z*dtD.z+D2KL.cz;
          }
          else
          {
            HessianMatrix.element[index_l][index_k]+=DDF*dtC.x*dtD.x+D2KL.ax;
            HessianMatrix.element[index_l][index_k+1]+=DDF*dtC.y*dtD.x+D2KL.bx;
            HessianMatrix.element[index_l][index_k+2]+=DDF*dtC.z*dtD.x+D2KL.cx;
            HessianMatrix.element[index_l+1][index_k]+=DDF*dtC.x*dtD.y+D2KL.ay;
            HessianMatrix.element[index_l+1][index_k+1]+=DDF*dtC.y*dtD.y+D2KL.by;
            HessianMatrix.element[index_l+1][index_k+2]+=DDF*dtC.z*dtD.y+D2KL.cy;
            HessianMatrix.element[index_l+2][index_k]+=DDF*dtC.x*dtD.z+D2KL.az;
            HessianMatrix.element[index_l+2][index_k+1]+=DDF*dtC.y*dtD.z+D2KL.bz;
            HessianMatrix.element[index_l+2][index_k+2]+=DDF*dtC.z*dtD.z+D2KL.cz;
          }
        }

        HessianTorsionStrainPosition(HessianMatrix,index_i,index_j,index_k,index_l,
             vec_u,vec_v,vec_w,u,v,w,fa,fc,fd,Dab,Dcb,Ddc,rbc,D2I,D2K,D2L,
             D2IJ,D2IK,D2IL,D2JK,D2JL,D2KL,dtA,dtB,dtC,dtD,DDF,S,CosPhi);

        HessianTorsionStrainStrain(HessianMatrix,vec_u,vec_v,vec_w,u,v,w,DF,DDF,S,CosPhi);
      }
    }
  }
}

void ComputeFrameworkImproperTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL d,e,rbc;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL U,CosPhi,Phi,CosPhi2,SinPhi;
  VECTOR dtA,dtB,dtC,dtD,Pb,Pc;
  REAL *parms;
  int index_i,index_j,index_k,index_l;
  REAL_MATRIX3x3 D2I,D2L,D2IL,D2IJ,D2IK,D2JL;
  REAL_MATRIX3x3 D2KL,D2J,D2K,D2JK,S;
  VECTOR DIL,DDJ,DDK,DEJ,DEK;
  REAL DF,DDF;
  REAL u,w,v;
  VECTOR vec_u,vec_v,vec_w;
  VECTOR fa,fb,fc,fd;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfImproperTorsions[f1];i++)
    {
      A=Framework[CurrentSystem].ImproperTorsions[f1][i].A;
      B=Framework[CurrentSystem].ImproperTorsions[f1][i].B;
      C=Framework[CurrentSystem].ImproperTorsions[f1][i].C;
      D=Framework[CurrentSystem].ImproperTorsions[f1][i].D;
      parms=Framework[CurrentSystem].ImproperTorsionArguments[f1][i];

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD=Framework[CurrentSystem].Atoms[f1][D].Position;

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;
      index_k=Framework[CurrentSystem].Atoms[f1][C].HessianIndex;
      index_l=Framework[CurrentSystem].Atoms[f1][D].HessianIndex;

      vec_u.x=posB.x-posC.x;
      vec_u.y=posB.y-posC.y;
      vec_u.z=posB.z-posC.z;
      vec_u=ApplyBoundaryCondition(vec_u);
      u=sqrt(SQR(vec_u.x)+SQR(vec_u.y)+SQR(vec_u.z));

      vec_v.x=posD.x-posC.x;
      vec_v.y=posD.y-posC.y;
      vec_v.z=posD.z-posC.z;
      vec_v=ApplyBoundaryCondition(vec_v);
      v=sqrt(SQR(vec_v.x)+SQR(vec_v.y)+SQR(vec_v.z));

      vec_w.x=posA.x-posB.x;
      vec_w.y=posA.y-posB.y;
      vec_w.z=posA.z-posB.z;
      vec_w=ApplyBoundaryCondition(vec_w);
      w=sqrt(SQR(vec_w.x)+SQR(vec_w.y)+SQR(vec_w.z));

      Dab.x=posA.x-posB.x;
      Dab.y=posA.y-posB.y;
      Dab.z=posA.z-posB.z;
      Dab=ApplyBoundaryCondition(Dab);

      Dcb.x=posC.x-posB.x;
      Dcb.y=posC.y-posB.y;
      Dcb.z=posC.z-posB.z;
      Dcb=ApplyBoundaryCondition(Dcb);
      rbc=sqrt(SQR(Dcb.x)+SQR(Dcb.y)+SQR(Dcb.z));
      Dcb.x/=rbc; Dcb.y/=rbc; Dcb.z/=rbc;

      Ddc.x=posD.x-posC.x;
      Ddc.y=posD.y-posC.y;
      Ddc.z=posD.z-posC.z;
      Ddc=ApplyBoundaryCondition(Ddc);

      dot_ab=Dab.x*Dcb.x+Dab.y*Dcb.y+Dab.z*Dcb.z;
      dot_cd=Ddc.x*Dcb.x+Ddc.y*Dcb.y+Ddc.z*Dcb.z;

      dr.x=Dab.x-dot_ab*Dcb.x;
      dr.y=Dab.y-dot_ab*Dcb.y;
      dr.z=Dab.z-dot_ab*Dcb.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      dr.x/=r; dr.y/=r; dr.z/=r;

      ds.x=Ddc.x-dot_cd*Dcb.x;
      ds.y=Ddc.y-dot_cd*Dcb.y;
      ds.z=Ddc.z-dot_cd*Dcb.z;
      s=sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z));
      ds.x/=s; ds.y/=s; ds.z/=s;

      // compute Cos(Phi)
      // Phi is defined in protein convention Phi(trans)=Pi
      CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;

      // Ensure CosPhi is between -1 and 1.
      CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
      CosPhi2=SQR(CosPhi);

      switch(Framework[CurrentSystem].ImproperTorsionType[f1][i])
      {
        case HARMONIC_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(phi-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // potential defined in terms of 'phi' and therefore contains a singularity
          // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
          // same direction as Rbc, and negative otherwise
          Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
          Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
          Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
          Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
          Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
          Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
          sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
          Phi=SIGN(acos(CosPhi),sign);
          SinPhi=sin(Phi);
          Phi-=parms[1];
          Phi-=NINT(Phi/(2.0*M_PI))*2.0*M_PI;
          SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi); // remove singularity
          U=0.5*parms[0]*SQR(Phi);
          DF=-parms[0]*Phi/SinPhi;
          DDF=-parms[0]*(Phi*CosPhi-SinPhi)/CUBE(SinPhi);
          break;
        case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(cos(phi)-cos(p_1))^2
          // ===============================================
          // p_0/k_B [K]
          // p_1     [degrees]
          U=0.5*parms[0]*SQR(CosPhi-parms[1]);
          DF=parms[0]*(CosPhi-parms[1]);
          DDF=parms[0];
          break;
        case THREE_COSINE_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
          DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[1]-6.0*parms[2]*CosPhi);
          break;
        case MM3_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0     [kcal/mol]
          // p_1     [kcal/mol]
          // p_2     [kcal/mol]
          U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
          DF=0.5*parms[0]-2.0*parms[1]*CosPhi+1.5*parms[2]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[1]-6.0*parms[2]*CosPhi);
          break;
        case CFF_IMPROPER_DIHEDRAL:
          // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
          DF=-parms[0]-4.0*parms[1]*CosPhi+3.0*parms[2]*(1.0-4.0*CosPhi2);
          DDF=-4.0*(parms[1]+6.0*parms[2]*CosPhi);
          break;
        case CFF_IMPROPER_DIHEDRAL2:
          // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
          DF=parms[0]-3.0*parms[2]+4.0*CosPhi*(parms[1]+3.0*parms[2]*CosPhi);
          DDF=4.0*(parms[1]+6.0*parms[2]*CosPhi);
          break;
        case SIX_COSINE_IMPROPER_DIHEDRAL:
          // Prod_i=0^5 p_i*cos(phi)^i
          // =========================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
          // between the first and last atoms of the dihedral, and Phi'=Phi-pi is defined accoording to the
          // polymer convention Phi'(trans)=0.
          U=parms[0]-parms[1]*CosPhi+parms[2]*CosPhi2-parms[3]*CosPhi*CosPhi2+
                 parms[4]*SQR(CosPhi2)-parms[5]*SQR(CosPhi2)*CosPhi;
          DF=-parms[1]+2.0*parms[2]*CosPhi-3.0*parms[3]*CosPhi2+4.0*parms[4]*CosPhi2*CosPhi-5.0*parms[5]*SQR(CosPhi2);
          DDF=2.0*parms[2]-6.0*parms[3]*CosPhi+12.0*parms[4]*CosPhi2-20.0*parms[5]*CosPhi2*CosPhi;
          break;
        case TRAPPE_IMPROPER_DIHEDRAL:
          U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
          DF=parms[1]-4.0*parms[2]*CosPhi+3.0*parms[3]*(4.0*CosPhi2-1.0);
          DDF=-4.0*(parms[2]-6.0*parms[3]*CosPhi);
          // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
          // =============================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          break;
        case CVFF_IMPROPER_DIHEDRAL:
          // p_0*(1+cos(p_1*phi-p_2))
          // ========================
          // p_0/k_B [K]
          // p_1     [-]
          // p_2     [degrees]
          // potential defined in terms of 'phi' and therefore contains a singularity
          // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
          // same direction as Dbc, and negative otherwise
          Pb.x=Dab.z*Dcb.y-Dab.y*Dcb.z;
          Pb.y=Dab.x*Dcb.z-Dab.z*Dcb.x;
          Pb.z=Dab.y*Dcb.x-Dab.x*Dcb.y;
          Pc.x=Dcb.y*Ddc.z-Dcb.z*Ddc.y;
          Pc.y=Dcb.z*Ddc.x-Dcb.x*Ddc.z;
          Pc.z=Dcb.x*Ddc.y-Dcb.y*Ddc.x;
          sign=(Dcb.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dcb.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                +Dcb.z*(Pc.y*Pb.x-Pc.x*Pb.y));
          Phi=SIGN(acos(CosPhi),sign);
          SinPhi=sin(Phi);
          U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
          if(fabs(SinPhi)<1.0e-8)
          {
            SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
            DF=parms[0]*parms[1]*cos(parms[2])*SIGN(parms[1],sin(parms[1]*Phi)*Phi)+parms[0]*parms[1]*cos(parms[1]*Phi)*sin(parms[2])/SinPhi;
          }
          else
            DF=(parms[0]*parms[1]*sin(parms[1]*Phi-parms[2]))/SinPhi;
          DDF=-parms[0]*parms[1]*(parms[1]*cos(-parms[1]*Phi+parms[2])+sin(-parms[1]*Phi+parms[2])*CosPhi/SinPhi)/SQR(SinPhi);
          break;
        case OPLS_IMPROPER_DIHEDRAL:
          // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
          // =================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
          DF=0.5*parms[1]-2.0*parms[2]*CosPhi+1.5*parms[3]*(4.0*CosPhi2-1.0);
          DDF=-2.0*(parms[2]-6.0*parms[3]*CosPhi);
          break;
        case FOURIER_SERIES_IMPROPER_DIHEDRAL:
          // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
          // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
          // =======================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          U=0.5*(parms[0]+2.0*parms[1]+parms[2]+parms[4]+2.0*parms[5]+(parms[0]-3.0*parms[2]+5.0*parms[4])*CosPhi-
                 2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi2+4.0*(parms[2]-5.0*parms[4])*CosPhi2*CosPhi-
                 8.0*(parms[3]-6.0*parms[5])*SQR(CosPhi2)+16.0*parms[4]*SQR(CosPhi2)*CosPhi-32.0*parms[5]*CUBE(CosPhi2));
          DF=0.5*(parms[0]-3.0*parms[2]+5.0*parms[4])-2.0*(parms[1]-4.0*parms[3]+9.0*parms[5])*CosPhi+
             6.0*(parms[2]-5.0*parms[4])*CosPhi2-16.0*(parms[3]-6.0*parms[5])*CosPhi2*CosPhi+40.0*parms[4]*SQR(CosPhi2)-96.0*parms[5]*CosPhi2*CUBE(CosPhi);
          DDF=-2.0*parms[1]+8.0*parms[3]-18.0*parms[5]+12.0*(parms[2]-5.0*parms[4])*CosPhi
              -48.0*(parms[3]-6.0*parms[5])*CosPhi2+160.0*parms[4]*CosPhi2*CosPhi-480.0*parms[5]*SQR(CosPhi2);
          break;
        case FOURIER_SERIES_IMPROPER_DIHEDRAL2:
          // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
          // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
          // =======================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          U=0.5*(parms[2]+2.0*parms[3]+parms[4]-3.0*parms[2]*CosPhi+5.0*parms[4]*CosPhi+parms[0]*(1.0+CosPhi)+
            2.0*(parms[1]-parms[1]*CosPhi2+CosPhi2*(parms[5]*SQR(3.0-4.0*CosPhi2)+4.0*parms[3]*(CosPhi2-1.0)+
            2.0*CosPhi*(parms[2]+parms[4]*(4.0*CosPhi2-5.0)))));
          DF=0.5*parms[0]+parms[2]*(6.0*CosPhi2-1.5)+parms[4]*(2.5-30.0*CosPhi2+40.0*SQR(CosPhi2))+
             CosPhi*(-2.0*parms[1]+parms[3]*(16.0*CosPhi2-8.0)+parms[5]*(18.0-96.0*CosPhi2+96.0*SQR(CosPhi2)));
          DDF=-2.0*parms[1]+12.0*parms[2]*CosPhi+parms[3]*(48.0*CosPhi2-8.0)+parms[4]*CosPhi*(160.0*CosPhi2-60.0)
              +parms[5]*(18.0-288.0*CosPhi2+480.0*SQR(CosPhi2));
          break;
        case FIXED_IMPROPER_DIHEDRAL:
          U=DF=DDF=0.0;
          break;
        default:
          printf("Undefined Improper Torsion potential in 'framework_hessian.c'\n");
          U=DF=DDF=0.0;
          exit(0);
          break;
      }

      *Energy+=U;

      if((index_i<0)&&(index_j<0)&&(index_k<0)&&(index_l<0)) continue;

     // Calculate the first derivative vectors.
      d=dot_ab/rbc;
      e=dot_cd/rbc;

      dtA.x=(ds.x-CosPhi*dr.x)/r;
      dtA.y=(ds.y-CosPhi*dr.y)/r;
      dtA.z=(ds.z-CosPhi*dr.z)/r;

      dtD.x=(dr.x-CosPhi*ds.x)/s;
      dtD.y=(dr.y-CosPhi*ds.y)/s;
      dtD.z=(dr.z-CosPhi*ds.z)/s;

      dtB.x=dtA.x*(d-1.0)+e*dtD.x;
      dtB.y=dtA.y*(d-1.0)+e*dtD.y;
      dtB.z=dtA.z*(d-1.0)+e*dtD.z;

      dtC.x=-dtD.x*(e+1.0)-d*dtA.x;
      dtC.y=-dtD.y*(e+1.0)-d*dtA.y;
      dtC.z=-dtD.z*(e+1.0)-d*dtA.z;

      // forces are oppositely directed to the gradient
      fa.x=DF*dtA.x;
      fa.y=DF*dtA.y;
      fa.z=DF*dtA.z;

      fb.x=DF*dtB.x;
      fb.y=DF*dtB.y;
      fb.z=DF*dtB.z;

      fc.x=DF*dtC.x;
      fc.y=DF*dtC.y;
      fc.z=DF*dtC.z;

      fd.x=DF*dtD.x;
      fd.y=DF*dtD.y;
      fd.z=DF*dtD.z;

      // add contribution to the strain derivative tensor
      // Note: rbc is here because the vector was normalized before
      S.ax=Dab.x*fa.x+Dcb.x*rbc*(fc.x+fd.x)+Ddc.x*fd.x;
      S.bx=Dab.y*fa.x+Dcb.y*rbc*(fc.x+fd.x)+Ddc.y*fd.x;
      S.cx=Dab.z*fa.x+Dcb.z*rbc*(fc.x+fd.x)+Ddc.z*fd.x;

      S.ay=Dab.x*fa.y+Dcb.x*rbc*(fc.y+fd.y)+Ddc.x*fd.y;
      S.by=Dab.y*fa.y+Dcb.y*rbc*(fc.y+fd.y)+Ddc.y*fd.y;
      S.cy=Dab.z*fa.y+Dcb.z*rbc*(fc.y+fd.y)+Ddc.z*fd.y;

      S.az=Dab.x*fa.z+Dcb.x*rbc*(fc.z+fd.z)+Ddc.x*fd.z;
      S.bz=Dab.y*fa.z+Dcb.y*rbc*(fc.z+fd.z)+Ddc.y*fd.z;
      S.cz=Dab.z*fa.z+Dcb.z*rbc*(fc.z+fd.z)+Ddc.z*fd.z;

      StrainDerivative->ax+=S.ax;
      StrainDerivative->bx+=S.bx;
      StrainDerivative->cx+=S.cx;

      StrainDerivative->ay+=S.ay;
      StrainDerivative->by+=S.by;
      StrainDerivative->cy+=S.cy;

      StrainDerivative->az+=S.az;
      StrainDerivative->bz+=S.bz;
      StrainDerivative->cz+=S.cz;

      if(ComputeGradient)
      {
        if(index_i>=0)
        {
          Gradient[index_i]+=fa.x;
          Gradient[index_i+1]+=fa.y;
          Gradient[index_i+2]+=fa.z;
        }

        if(index_j>=0)
        {
          Gradient[index_j]+=fb.x;
          Gradient[index_j+1]+=fb.y;
          Gradient[index_j+2]+=fb.z;
        }

        if(index_k>=0)
        {
          Gradient[index_k]+=fc.x;
          Gradient[index_k+1]+=fc.y;
          Gradient[index_k+2]+=fc.z;
        }

        if(index_l>=0)
        {
          Gradient[index_l]+=fd.x;
          Gradient[index_l+1]+=fd.y;
          Gradient[index_l+2]+=fd.z;
        }

        GradientStrainTorsion(Gradient,S);
      }

      if(ComputeHessian)
      {
        // Calculate the derivatives of DOTIJ and DOTLK.
        DIL.x=DF*Dcb.x/rbc;
        DIL.y=DF*Dcb.y/rbc;
        DIL.z=DF*Dcb.z/rbc;
        DDJ.x=DF*((2.0*d-1.0)*Dcb.x-Dab.x/rbc)/rbc;
        DDJ.y=DF*((2.0*d-1.0)*Dcb.y-Dab.y/rbc)/rbc;
        DDJ.z=DF*((2.0*d-1.0)*Dcb.z-Dab.z/rbc)/rbc;
        DDK.x=-DF*(2.0*d*Dcb.x-Dab.x/rbc)/rbc;
        DDK.y=-DF*(2.0*d*Dcb.y-Dab.y/rbc)/rbc;
        DDK.z=-DF*(2.0*d*Dcb.z-Dab.z/rbc)/rbc;
        DEJ.x=DF*(2.0*e*Dcb.x-Ddc.x/rbc)/rbc;
        DEJ.y=DF*(2.0*e*Dcb.y-Ddc.y/rbc)/rbc;
        DEJ.z=DF*(2.0*e*Dcb.z-Ddc.z/rbc)/rbc;
        DEK.x=-DF*((2.0*e+1.0)*Dcb.x-Ddc.x/rbc)/rbc;
        DEK.y=-DF*((2.0*e+1.0)*Dcb.y-Ddc.y/rbc)/rbc;
        DEK.z=-DF*((2.0*e+1.0)*Dcb.z-Ddc.z/rbc)/rbc;

        // Calculate some diagonal Hessian terms for I.
        D2I.ax=DF*(CosPhi*(dr.x*dr.x+Dcb.x*Dcb.x-1.0)/r-2.0*dr.x*dtA.x)/r;
        D2I.by=DF*(CosPhi*(dr.y*dr.y+Dcb.y*Dcb.y-1.0)/r-2.0*dr.y*dtA.y)/r;
        D2I.cz=DF*(CosPhi*(dr.z*dr.z+Dcb.z*Dcb.z-1.0)/r-2.0*dr.z*dtA.z)/r;
        D2I.ay=DF*(CosPhi*(dr.x*dr.y+Dcb.x*Dcb.y)/r-dr.x*dtA.y-dr.y*dtA.x)/r;
        D2I.az=DF*(CosPhi*(dr.x*dr.z+Dcb.x*Dcb.z)/r-dr.x*dtA.z-dr.z*dtA.x)/r;
        D2I.bz=DF*(CosPhi*(dr.y*dr.z+Dcb.y*Dcb.z)/r-dr.y*dtA.z-dr.z*dtA.y)/r;

        // Calculate some diagonal Hessian terms for L.
        D2L.ax=DF*(CosPhi*(ds.x*ds.x+Dcb.x*Dcb.x-1.0)/s-2.0*ds.x*dtD.x)/s;
        D2L.by=DF*(CosPhi*(ds.y*ds.y+Dcb.y*Dcb.y-1.0)/s-2.0*ds.y*dtD.y)/s;
        D2L.cz=DF*(CosPhi*(ds.z*ds.z+Dcb.z*Dcb.z-1.0)/s-2.0*ds.z*dtD.z)/s;
        D2L.ay=DF*(CosPhi*(ds.x*ds.y+Dcb.x*Dcb.y)/s-ds.x*dtD.y-ds.y*dtD.x)/s;
        D2L.az=DF*(CosPhi*(ds.x*ds.z+Dcb.x*Dcb.z)/s-ds.x*dtD.z-ds.z*dtD.x)/s;
        D2L.bz=DF*(CosPhi*(ds.y*ds.z+Dcb.y*Dcb.z)/s-ds.y*dtD.z-ds.z*dtD.y)/s;

        // Calculate the IL off-diagonal terms.
        D2IL.ax=DF*(CosPhi*dr.x*ds.x-dr.x*dr.x-ds.x*ds.x-Dcb.x*Dcb.x+1.0)/(r*s);
        D2IL.ay=DF*(CosPhi*dr.x*ds.y-dr.x*dr.y-ds.x*ds.y-Dcb.x*Dcb.y)/(r*s);
        D2IL.az=DF*(CosPhi*dr.x*ds.z-dr.x*dr.z-ds.x*ds.z-Dcb.x*Dcb.z)/(r*s);
        D2IL.bx=DF*(CosPhi*dr.y*ds.x-dr.y*dr.x-ds.y*ds.x-Dcb.y*Dcb.x)/(r*s);
        D2IL.by=DF*(CosPhi*dr.y*ds.y-dr.y*dr.y-ds.y*ds.y-Dcb.y*Dcb.y+1.0)/(r*s);
        D2IL.bz=DF*(CosPhi*dr.y*ds.z-dr.y*dr.z-ds.y*ds.z-Dcb.y*Dcb.z)/(r*s);
        D2IL.cx=DF*(CosPhi*dr.z*ds.x-dr.z*dr.x-ds.z*ds.x-Dcb.z*Dcb.x)/(r*s);
        D2IL.cy=DF*(CosPhi*dr.z*ds.y-dr.z*dr.y-ds.z*ds.y-Dcb.z*Dcb.y)/(r*s);
        D2IL.cz=DF*(CosPhi*dr.z*ds.z-dr.z*dr.z-ds.z*ds.z-Dcb.z*Dcb.z+1.0)/(r*s);

        // Calculate the IJ off-diagonal terms.
        D2IJ.ax=D2I.ax*(d-1.0)+D2IL.ax*e+DIL.x*dtA.x;
        D2IJ.ay=D2I.ay*(d-1.0)+D2IL.ay*e+DIL.x*dtA.y;
        D2IJ.az=D2I.az*(d-1.0)+D2IL.az*e+DIL.x*dtA.z;
        D2IJ.bx=D2I.ay*(d-1.0)+D2IL.bx*e+DIL.y*dtA.x;
        D2IJ.by=D2I.by*(d-1.0)+D2IL.by*e+DIL.y*dtA.y;
        D2IJ.bz=D2I.bz*(d-1.0)+D2IL.bz*e+DIL.y*dtA.z;
        D2IJ.cx=D2I.az*(d-1.0)+D2IL.cx*e+DIL.z*dtA.x;
        D2IJ.cy=D2I.bz*(d-1.0)+D2IL.cy*e+DIL.z*dtA.y;
        D2IJ.cz=D2I.cz*(d-1.0)+D2IL.cz*e+DIL.z*dtA.z;

        // Calculate the IK off-diagonal terms.
        D2IK.ax=-D2I.ax*d-D2IL.ax*(e+1.0)-DIL.x*dtA.x;
        D2IK.ay=-D2I.ay*d-D2IL.ay*(e+1.0)-DIL.x*dtA.y;
        D2IK.az=-D2I.az*d-D2IL.az*(e+1.0)-DIL.x*dtA.z;
        D2IK.bx=-D2I.ay*d-D2IL.bx*(e+1.0)-DIL.y*dtA.x;
        D2IK.by=-D2I.by*d-D2IL.by*(e+1.0)-DIL.y*dtA.y;
        D2IK.bz=-D2I.bz*d-D2IL.bz*(e+1.0)-DIL.y*dtA.z;
        D2IK.cx=-D2I.az*d-D2IL.cx*(e+1.0)-DIL.z*dtA.x;
        D2IK.cy=-D2I.bz*d-D2IL.cy*(e+1.0)-DIL.z*dtA.y;
        D2IK.cz=-D2I.cz*d-D2IL.cz*(e+1.0)-DIL.z*dtA.z;

        // Calculate the JL off-diagonal terms.
        D2JL.ax=D2IL.ax*(d-1.0)+D2L.ax*e+DIL.x*dtD.x;
        D2JL.ay=D2IL.ay*(d-1.0)+D2L.ay*e+DIL.y*dtD.x;
        D2JL.az=D2IL.az*(d-1.0)+D2L.az*e+DIL.z*dtD.x;
        D2JL.bx=D2IL.bx*(d-1.0)+D2L.ay*e+DIL.x*dtD.y;
        D2JL.by=D2IL.by*(d-1.0)+D2L.by*e+DIL.y*dtD.y;
        D2JL.bz=D2IL.bz*(d-1.0)+D2L.bz*e+DIL.z*dtD.y;
        D2JL.cx=D2IL.cx*(d-1.0)+D2L.az*e+DIL.x*dtD.z;
        D2JL.cy=D2IL.cy*(d-1.0)+D2L.bz*e+DIL.y*dtD.z;
        D2JL.cz=D2IL.cz*(d-1.0)+D2L.cz*e+DIL.z*dtD.z;

        // Calculate the KL off-diagonal terms.
        D2KL.ax=-D2IL.ax*d-D2L.ax*(e+1.0)-DIL.x*dtD.x;
        D2KL.ay=-D2IL.ay*d-D2L.ay*(e+1.0)-DIL.y*dtD.x;
        D2KL.az=-D2IL.az*d-D2L.az*(e+1.0)-DIL.z*dtD.x;
        D2KL.bx=-D2IL.bx*d-D2L.ay*(e+1.0)-DIL.x*dtD.y;
        D2KL.by=-D2IL.by*d-D2L.by*(e+1.0)-DIL.y*dtD.y;
        D2KL.bz=-D2IL.bz*d-D2L.bz*(e+1.0)-DIL.z*dtD.y;
        D2KL.cx=-D2IL.cx*d-D2L.az*(e+1.0)-DIL.x*dtD.z;
        D2KL.cy=-D2IL.cy*d-D2L.bz*(e+1.0)-DIL.y*dtD.z;
        D2KL.cz=-D2IL.cz*d-D2L.cz*(e+1.0)-DIL.z*dtD.z;

        // Calculate the JJ diagonal terms.
        D2J.ax=D2IJ.ax*(d-1.0)+D2JL.ax*e+DDJ.x*dtA.x+DEJ.x*dtD.x;
        D2J.by=D2IJ.by*(d-1.0)+D2JL.by*e+DDJ.y*dtA.y+DEJ.y*dtD.y;
        D2J.cz=D2IJ.cz*(d-1.0)+D2JL.cz*e+DDJ.z*dtA.z+DEJ.z*dtD.z;
        D2J.ay=D2IJ.ay*(d-1.0)+D2JL.bx*e+DDJ.y*dtA.x+DEJ.y*dtD.x;
        D2J.az=D2IJ.az*(d-1.0)+D2JL.cx*e+DDJ.z*dtA.x+DEJ.z*dtD.x;
        D2J.bz=D2IJ.bz*(d-1.0)+D2JL.cy*e+DDJ.z*dtA.y+DEJ.z*dtD.y;

        // Calculate the KK diagonal terms.
        D2K.ax=-D2KL.ax*(e+1.0)-D2IK.ax*d-DDK.x*dtA.x-DEK.x*dtD.x;
        D2K.by=-D2KL.by*(e+1.0)-D2IK.by*d-DDK.y*dtA.y-DEK.y*dtD.y;
        D2K.cz=-D2KL.cz*(e+1.0)-D2IK.cz*d-DDK.z*dtA.z-DEK.z*dtD.z;
        D2K.ay=-D2KL.ay*(e+1.0)-D2IK.bx*d-DDK.x*dtA.y-DEK.x*dtD.y;
        D2K.az=-D2KL.az*(e+1.0)-D2IK.cx*d-DDK.x*dtA.z-DEK.x*dtD.z;
        D2K.bz=-D2KL.bz*(e+1.0)-D2IK.cy*d-DDK.y*dtA.z-DEK.y*dtD.z;

        // Calculate the JK off-diagonal terms.
        D2JK.ax=-D2IJ.ax*d-D2JL.ax*(e+1.0)-DDJ.x*dtA.x-DEJ.x*dtD.x;
        D2JK.ay=-D2IJ.bx*d-D2JL.ay*(e+1.0)-DDJ.x*dtA.y-DEJ.x*dtD.y;
        D2JK.az=-D2IJ.cx*d-D2JL.az*(e+1.0)-DDJ.x*dtA.z-DEJ.x*dtD.z;
        D2JK.bx=-D2IJ.ay*d-D2JL.bx*(e+1.0)-DDJ.y*dtA.x-DEJ.y*dtD.x;
        D2JK.by=-D2IJ.by*d-D2JL.by*(e+1.0)-DDJ.y*dtA.y-DEJ.y*dtD.y;
        D2JK.bz=-D2IJ.cy*d-D2JL.bz*(e+1.0)-DDJ.y*dtA.z-DEJ.y*dtD.z;
        D2JK.cx=-D2IJ.az*d-D2JL.cx*(e+1.0)-DDJ.z*dtA.x-DEJ.z*dtD.x;
        D2JK.cy=-D2IJ.bz*d-D2JL.cy*(e+1.0)-DDJ.z*dtA.y-DEJ.z*dtD.y;
        D2JK.cz=-D2IJ.cz*d-D2JL.cz*(e+1.0)-DDJ.z*dtA.z-DEJ.z*dtD.z;

        // Calculate the diagonal blocks of the hessian.
        if(index_i>=0)
        {
          HessianMatrix.element[index_i][index_i]+=DDF*dtA.x*dtA.x+D2I.ax;
          HessianMatrix.element[index_i][index_i+1]+=DDF*dtA.x*dtA.y+D2I.ay;
          HessianMatrix.element[index_i+1][index_i+1]+=DDF*dtA.y*dtA.y+D2I.by;
          HessianMatrix.element[index_i][index_i+2]+=DDF*dtA.x*dtA.z+D2I.az;
          HessianMatrix.element[index_i+1][index_i+2]+=DDF*dtA.y*dtA.z+D2I.bz;
          HessianMatrix.element[index_i+2][index_i+2]+=DDF*dtA.z*dtA.z+D2I.cz;
        }

        if(index_j>=0)
        {
          HessianMatrix.element[index_j][index_j]+=DDF*dtB.x*dtB.x+D2J.ax;
          HessianMatrix.element[index_j][index_j+1]+=DDF*dtB.x*dtB.y+D2J.ay;
          HessianMatrix.element[index_j+1][index_j+1]+=DDF*dtB.y*dtB.y+D2J.by;
          HessianMatrix.element[index_j][index_j+2]+=DDF*dtB.x*dtB.z+D2J.az;
          HessianMatrix.element[index_j+1][index_j+2]+=DDF*dtB.y*dtB.z+D2J.bz;
          HessianMatrix.element[index_j+2][index_j+2]+=DDF*dtB.z*dtB.z+D2J.cz;
        }

        if(index_k>=0)
        {
          HessianMatrix.element[index_k][index_k]+=DDF*dtC.x*dtC.x+D2K.ax;
          HessianMatrix.element[index_k][index_k+1]+=DDF*dtC.x*dtC.y+D2K.ay;
          HessianMatrix.element[index_k+1][index_k+1]+=DDF*dtC.y*dtC.y+D2K.by;
          HessianMatrix.element[index_k][index_k+2]+=DDF*dtC.x*dtC.z+D2K.az;
          HessianMatrix.element[index_k+1][index_k+2]+=DDF*dtC.y*dtC.z+D2K.bz;
          HessianMatrix.element[index_k+2][index_k+2]+=DDF*dtC.z*dtC.z+D2K.cz;
        }

        if(index_l>=0)
        {
          HessianMatrix.element[index_l][index_l]+=DDF*dtD.x*dtD.x+D2L.ax;
          HessianMatrix.element[index_l][index_l+1]+=DDF*dtD.x*dtD.y+D2L.ay;
          HessianMatrix.element[index_l+1][index_l+1]+=DDF*dtD.y*dtD.y+D2L.by;
          HessianMatrix.element[index_l][index_l+2]+=DDF*dtD.x*dtD.z+D2L.az;
          HessianMatrix.element[index_l+1][index_l+2]+=DDF*dtD.y*dtD.z+D2L.bz;
          HessianMatrix.element[index_l+2][index_l+2]+=DDF*dtD.z*dtD.z+D2L.cz;
        }

        if((index_i>=0)&&(index_j>=0))
        {
          if(index_i<index_j)
          {
            HessianMatrix.element[index_i][index_j]+=DDF*dtA.x*dtB.x+D2IJ.ax;
            HessianMatrix.element[index_i][index_j+1]+=DDF*dtA.x*dtB.y+D2IJ.ay;
            HessianMatrix.element[index_i][index_j+2]+=DDF*dtA.x*dtB.z+D2IJ.az;
            HessianMatrix.element[index_i+1][index_j]+=DDF*dtA.y*dtB.x+D2IJ.bx;
            HessianMatrix.element[index_i+1][index_j+1]+=DDF*dtA.y*dtB.y+D2IJ.by;
            HessianMatrix.element[index_i+1][index_j+2]+=DDF*dtA.y*dtB.z+D2IJ.bz;
            HessianMatrix.element[index_i+2][index_j]+=DDF*dtA.z*dtB.x+D2IJ.cx;
            HessianMatrix.element[index_i+2][index_j+1]+=DDF*dtA.z*dtB.y+D2IJ.cy;
            HessianMatrix.element[index_i+2][index_j+2]+=DDF*dtA.z*dtB.z+D2IJ.cz;
          }
          else
          {
            HessianMatrix.element[index_j][index_i]+=DDF*dtA.x*dtB.x+D2IJ.ax;
            HessianMatrix.element[index_j][index_i+1]+=DDF*dtA.y*dtB.x+D2IJ.bx;
            HessianMatrix.element[index_j][index_i+2]+=DDF*dtA.z*dtB.x+D2IJ.cx;
            HessianMatrix.element[index_j+1][index_i]+=DDF*dtA.x*dtB.y+D2IJ.ay;
            HessianMatrix.element[index_j+1][index_i+1]+=DDF*dtA.y*dtB.y+D2IJ.by;
            HessianMatrix.element[index_j+1][index_i+2]+=DDF*dtA.z*dtB.y+D2IJ.cy;
            HessianMatrix.element[index_j+2][index_i]+=DDF*dtA.x*dtB.z+D2IJ.az;
            HessianMatrix.element[index_j+2][index_i+1]+=DDF*dtA.y*dtB.z+D2IJ.bz;
            HessianMatrix.element[index_j+2][index_i+2]+=DDF*dtA.z*dtB.z+D2IJ.cz;
          }
        }

        if((index_i>=0)&&(index_k>=0))
        {
          if(index_i<index_k)
          {
            HessianMatrix.element[index_i][index_k]+=DDF*dtA.x*dtC.x+D2IK.ax;
            HessianMatrix.element[index_i][index_k+1]+=DDF*dtA.x*dtC.y+D2IK.ay;
            HessianMatrix.element[index_i][index_k+2]+=DDF*dtA.x*dtC.z+D2IK.az;
            HessianMatrix.element[index_i+1][index_k]+=DDF*dtA.y*dtC.x+D2IK.bx;
            HessianMatrix.element[index_i+1][index_k+1]+=DDF*dtA.y*dtC.y+D2IK.by;
            HessianMatrix.element[index_i+1][index_k+2]+=DDF*dtA.y*dtC.z+D2IK.bz;
            HessianMatrix.element[index_i+2][index_k]+=DDF*dtA.z*dtC.x+D2IK.cx;
            HessianMatrix.element[index_i+2][index_k+1]+=DDF*dtA.z*dtC.y+D2IK.cy;
            HessianMatrix.element[index_i+2][index_k+2]+=DDF*dtA.z*dtC.z+D2IK.cz;
          }
          else
          {
            HessianMatrix.element[index_k][index_i]+=DDF*dtA.x*dtC.x+D2IK.ax;
            HessianMatrix.element[index_k][index_i+1]+=DDF*dtA.y*dtC.x+D2IK.bx;
            HessianMatrix.element[index_k][index_i+2]+=DDF*dtA.z*dtC.x+D2IK.cx;
            HessianMatrix.element[index_k+1][index_i]+=DDF*dtA.x*dtC.y+D2IK.ay;
            HessianMatrix.element[index_k+1][index_i+1]+=DDF*dtA.y*dtC.y+D2IK.by;
            HessianMatrix.element[index_k+1][index_i+2]+=DDF*dtA.z*dtC.y+D2IK.cy;
            HessianMatrix.element[index_k+2][index_i]+=DDF*dtA.x*dtC.z+D2IK.az;
            HessianMatrix.element[index_k+2][index_i+1]+=DDF*dtA.y*dtC.z+D2IK.bz;
            HessianMatrix.element[index_k+2][index_i+2]+=DDF*dtA.z*dtC.z+D2IK.cz;
          }
        }

        if((index_i>=0)&&(index_l>=0))
        {
          if(index_i<index_l)
          {
            HessianMatrix.element[index_i][index_l]+=DDF*dtA.x*dtD.x+D2IL.ax;
            HessianMatrix.element[index_i][index_l+1]+=DDF*dtA.x*dtD.y+D2IL.ay; 
            HessianMatrix.element[index_i][index_l+2]+=DDF*dtA.x*dtD.z+D2IL.az; 
            HessianMatrix.element[index_i+1][index_l]+=DDF*dtA.y*dtD.x+D2IL.bx; 
            HessianMatrix.element[index_i+1][index_l+1]+=DDF*dtA.y*dtD.y+D2IL.by;
            HessianMatrix.element[index_i+1][index_l+2]+=DDF*dtA.y*dtD.z+D2IL.bz;
            HessianMatrix.element[index_i+2][index_l]+=DDF*dtA.z*dtD.x+D2IL.cx;
            HessianMatrix.element[index_i+2][index_l+1]+=DDF*dtA.z*dtD.y+D2IL.cy;
            HessianMatrix.element[index_i+2][index_l+2]+=DDF*dtA.z*dtD.z+D2IL.cz;
          }
          else
          {
            HessianMatrix.element[index_l][index_i]+=DDF*dtA.x*dtD.x+D2IL.ax;
            HessianMatrix.element[index_l][index_i+1]+=DDF*dtA.y*dtD.x+D2IL.bx;
            HessianMatrix.element[index_l][index_i+2]+=DDF*dtA.z*dtD.x+D2IL.cx;
            HessianMatrix.element[index_l+1][index_i]+=DDF*dtA.x*dtD.y+D2IL.ay;
            HessianMatrix.element[index_l+1][index_i+1]+=DDF*dtA.y*dtD.y+D2IL.by;
            HessianMatrix.element[index_l+1][index_i+2]+=DDF*dtA.z*dtD.y+D2IL.cy;
            HessianMatrix.element[index_l+2][index_i]+=DDF*dtA.x*dtD.z+D2IL.az;
            HessianMatrix.element[index_l+2][index_i+1]+=DDF*dtA.y*dtD.z+D2IL.bz;
            HessianMatrix.element[index_l+2][index_i+2]+=DDF*dtA.z*dtD.z+D2IL.cz;
          }
        }

        if((index_j>=0)&&(index_k>=0))
        {
          if(index_j<index_k)
          {
            HessianMatrix.element[index_j][index_k]+=DDF*dtB.x*dtC.x+D2JK.ax;
            HessianMatrix.element[index_j][index_k+1]+=DDF*dtB.x*dtC.y+D2JK.ay; 
            HessianMatrix.element[index_j][index_k+2]+=DDF*dtB.x*dtC.z+D2JK.az;
            HessianMatrix.element[index_j+1][index_k]+=DDF*dtB.y*dtC.x+D2JK.bx;
            HessianMatrix.element[index_j+1][index_k+1]+=DDF*dtB.y*dtC.y+D2JK.by;
            HessianMatrix.element[index_j+1][index_k+2]+=DDF*dtB.y*dtC.z+D2JK.bz;
            HessianMatrix.element[index_j+2][index_k]+=DDF*dtB.z*dtC.x+D2JK.cx;
            HessianMatrix.element[index_j+2][index_k+1]+=DDF*dtB.z*dtC.y+D2JK.cy;
            HessianMatrix.element[index_j+2][index_k+2]+=DDF*dtB.z*dtC.z+D2JK.cz;
          }
          else
          {
            HessianMatrix.element[index_k][index_j]+=DDF*dtB.x*dtC.x+D2JK.ax;
            HessianMatrix.element[index_k][index_j+1]+=DDF*dtB.y*dtC.x+D2JK.bx;
            HessianMatrix.element[index_k][index_j+2]+=DDF*dtB.z*dtC.x+D2JK.cx;
            HessianMatrix.element[index_k+1][index_j]+=DDF*dtB.x*dtC.y+D2JK.ay;
            HessianMatrix.element[index_k+1][index_j+1]+=DDF*dtB.y*dtC.y+D2JK.by;
            HessianMatrix.element[index_k+1][index_j+2]+=DDF*dtB.z*dtC.y+D2JK.cy;
            HessianMatrix.element[index_k+2][index_j]+=DDF*dtB.x*dtC.z+D2JK.az;
            HessianMatrix.element[index_k+2][index_j+1]+=DDF*dtB.y*dtC.z+D2JK.bz;
            HessianMatrix.element[index_k+2][index_j+2]+=DDF*dtB.z*dtC.z+D2JK.cz;
          }
        }

        if((index_j>=0)&&(index_l>=0))
        {
          if(index_j<index_l)
          {
            HessianMatrix.element[index_j][index_l]+=DDF*dtB.x*dtD.x+D2JL.ax;
            HessianMatrix.element[index_j][index_l+1]+=DDF*dtB.x*dtD.y+D2JL.ay;
            HessianMatrix.element[index_j][index_l+2]+=DDF*dtB.x*dtD.z+D2JL.az;
            HessianMatrix.element[index_j+1][index_l]+=DDF*dtB.y*dtD.x+D2JL.bx;
            HessianMatrix.element[index_j+1][index_l+1]+=DDF*dtB.y*dtD.y+D2JL.by;
            HessianMatrix.element[index_j+1][index_l+2]+=DDF*dtB.y*dtD.z+D2JL.bz;
            HessianMatrix.element[index_j+2][index_l]+=DDF*dtB.z*dtD.x+D2JL.cx;
            HessianMatrix.element[index_j+2][index_l+1]+=DDF*dtB.z*dtD.y+D2JL.cy;
            HessianMatrix.element[index_j+2][index_l+2]+=DDF*dtB.z*dtD.z+D2JL.cz;
          }
          else
          {
            HessianMatrix.element[index_l][index_j]+=DDF*dtB.x*dtD.x+D2JL.ax;
            HessianMatrix.element[index_l][index_j+1]+=DDF*dtB.y*dtD.x+D2JL.bx;
            HessianMatrix.element[index_l][index_j+2]+=DDF*dtB.z*dtD.x+D2JL.cx;
            HessianMatrix.element[index_l+1][index_j]+=DDF*dtB.x*dtD.y+D2JL.ay;
            HessianMatrix.element[index_l+1][index_j+1]+=DDF*dtB.y*dtD.y+D2JL.by;
            HessianMatrix.element[index_l+1][index_j+2]+=DDF*dtB.z*dtD.y+D2JL.cy;
            HessianMatrix.element[index_l+2][index_j]+=DDF*dtB.x*dtD.z+D2JL.az;
            HessianMatrix.element[index_l+2][index_j+1]+=DDF*dtB.y*dtD.z+D2JL.bz;
            HessianMatrix.element[index_l+2][index_j+2]+=DDF*dtB.z*dtD.z+D2JL.cz;
          }
        }

        if((index_k>=0)&&(index_l>=0))
        {
          if(index_k<index_l)
          {
            HessianMatrix.element[index_k][index_l]+=DDF*dtC.x*dtD.x+D2KL.ax;
            HessianMatrix.element[index_k][index_l+1]+=DDF*dtC.x*dtD.y+D2KL.ay;
            HessianMatrix.element[index_k][index_l+2]+=DDF*dtC.x*dtD.z+D2KL.az;
            HessianMatrix.element[index_k+1][index_l]+=DDF*dtC.y*dtD.x+D2KL.bx;
            HessianMatrix.element[index_k+1][index_l+1]+=DDF*dtC.y*dtD.y+D2KL.by;
            HessianMatrix.element[index_k+1][index_l+2]+=DDF*dtC.y*dtD.z+D2KL.bz;
            HessianMatrix.element[index_k+2][index_l]+=DDF*dtC.z*dtD.x+D2KL.cx;
            HessianMatrix.element[index_k+2][index_l+1]+=DDF*dtC.z*dtD.y+D2KL.cy;
            HessianMatrix.element[index_k+2][index_l+2]+=DDF*dtC.z*dtD.z+D2KL.cz;
          }
          else
          {
            HessianMatrix.element[index_l][index_k]+=DDF*dtC.x*dtD.x+D2KL.ax;
            HessianMatrix.element[index_l][index_k+1]+=DDF*dtC.y*dtD.x+D2KL.bx;
            HessianMatrix.element[index_l][index_k+2]+=DDF*dtC.z*dtD.x+D2KL.cx;
            HessianMatrix.element[index_l+1][index_k]+=DDF*dtC.x*dtD.y+D2KL.ay;
            HessianMatrix.element[index_l+1][index_k+1]+=DDF*dtC.y*dtD.y+D2KL.by;
            HessianMatrix.element[index_l+1][index_k+2]+=DDF*dtC.z*dtD.y+D2KL.cy;
            HessianMatrix.element[index_l+2][index_k]+=DDF*dtC.x*dtD.z+D2KL.az;
            HessianMatrix.element[index_l+2][index_k+1]+=DDF*dtC.y*dtD.z+D2KL.bz;
            HessianMatrix.element[index_l+2][index_k+2]+=DDF*dtC.z*dtD.z+D2KL.cz;
          }
        }

        HessianTorsionStrainPosition(HessianMatrix,index_i,index_j,index_k,index_l,
               vec_u,vec_v,vec_w,u,v,w,fa,fc,fd,Dab,Dcb,Ddc,rbc,D2I,D2K,D2L,
               D2IJ,D2IK,D2IL,D2JK,D2JL,D2KL,dtA,dtB,dtC,dtD,DDF,S,CosPhi);

        HessianTorsionStrainStrain(HessianMatrix,vec_u,vec_v,vec_w,u,v,w,DF,DDF,S,CosPhi);
      }
    }
  }
}


// Numerical routines
// TODO: convert them to analytical expressions

void CalculateFrameworkInversionBendForces(int f1,int i,VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,REAL *Energy,VECTOR *fa,VECTOR *fb,VECTOR *fc,VECTOR *fd,REAL_MATRIX3x3 *strain_derivative)
{
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2,rrbc,rbc2,rbd2,rad2,rac2,dot;
  REAL CosChi,Chi,energy;
  REAL temp,temp2;
  VECTOR Rab,Rbc,Rbd,Rad,Rac;
  REAL term,dedcos;
  VECTOR dccdia,dccdic,dccdid;
  VECTOR deedia,deedic,deedid;

  parms=Framework[CurrentSystem].InversionBendArguments[f1][i];

        Rab.x=posA.x-posB.x;
        Rab.y=posA.y-posB.y;
        Rab.z=posA.z-posB.z;
        Rab=ApplyBoundaryCondition(Rab);
        rab2=Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z;
        rrab=sqrt(rab2);
  
        Rbc.x=posC.x-posB.x;
        Rbc.y=posC.y-posB.y;
        Rbc.z=posC.z-posB.z;
        Rbc=ApplyBoundaryCondition(Rbc);
        rbc2=Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z;
        rrbc=sqrt(rbc2);
  
        Rbd.x=posD.x-posB.x;
        Rbd.y=posD.y-posB.y;
        Rbd.z=posD.z-posB.z;
        Rbd=ApplyBoundaryCondition(Rbd);
        rbd2=Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z;
  
        Rac.x=posC.x-posA.x;
        Rac.y=posC.y-posA.y;
        Rac.z=posC.z-posA.z;
        Rac=ApplyBoundaryCondition(Rac);
        rac2=Rac.x*Rac.x+Rac.y*Rac.y+Rac.z*Rac.z;
  
        Rad.x=posD.x-posA.x;
        Rad.y=posD.y-posA.y;
        Rad.z=posD.z-posA.z;
        Rad=ApplyBoundaryCondition(Rad);
        rad2=Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z;
  
        switch(Framework[CurrentSystem].InversionBendType[f1][i])
        {
          case HARMONIC_INVERSION:
          case HARMONIC_COSINE_INVERSION:
          case PLANAR_INVERSION:
            // w is a vector perpendicular to the B-C-D plane
            // c=w.w=(Rbc x Rbd).(Rbc x Rbd)= r_bc^2 r_bd^2 - (r_cb . r_bd)^2
            dot=Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z;
            c=rbc2*rbd2-SQR(dot);
            //c=(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z)*(Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z)-SQR(Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z);
            break;
          case MM3_INVERSION:
          case HARMONIC_INVERSION2:
          case HARMONIC_COSINE_INVERSION2:
          case PLANAR_INVERSION2:
            // w is a vector perpendicular to the A-C-D plane
            // c=w.w=(Rcd x Rad).(Rcd x Rad)=r_cd^2 r_ad^2 - (r_da . r_cd)^2
            dot=Rac.x*Rad.x+Rac.y*Rad.y+Rac.z*Rad.z;
            c=rac2*rad2-SQR(dot);
            //c=(Rcd.x*Rcd.x+Rcd.y*Rcd.y+Rcd.z*Rcd.z)*(Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z)-SQR(Rad.x*Rcd.x+Rad.y*Rcd.y+Rad.z*Rcd.z);
            break;
          default:
            printf("Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
            exit(0);
            break;
        }
  
        e=Rab.x*(Rbc.y*Rbd.z-Rbc.z*Rbd.y)+Rab.y*(Rbc.z*Rbd.x-Rbc.x*Rbd.z)+Rab.z*(Rbc.x*Rbd.y-Rbc.y*Rbd.x);
        CosChi=sqrt((Rab.x*Rab.x+Rab.y*Rab.y+Rab.z*Rab.z)-SQR(e)/c)/rrab;
  
        // Ensure CosChi is between -1 and 1.
        CosChi=SIGN(MIN2(fabs(CosChi),(REAL)1.0),CosChi);
  
        parms=Framework[CurrentSystem].InversionBendArguments[f1][i];
  
        switch(Framework[CurrentSystem].InversionBendType[f1][i])
        {
          case HARMONIC_INVERSION:
          case HARMONIC_INVERSION2:
            // (1/2)*p_0*(chi-p_1)^2
            // ===============================================
            // p_0/k_B [K]
            // p_1     [degrees]
            Chi=acos(CosChi);
            energy=0.5*parms[0]*SQR(Chi-parms[1]);
            dedcos=-SIGN(1.0,e)*(parms[0]*(Chi-parms[1])/sqrt(c*(rab2-e*e/c)));
            break;
          case HARMONIC_COSINE_INVERSION:
          case HARMONIC_COSINE_INVERSION2:
            // (1/2)*p_0*(cos(phi)-cos(p_1))^2
            // ===============================================
            // p_0/k_B [K]
            // p_1     [degrees]
            Chi=acos(CosChi);
            energy=0.5*parms[0]*SQR(CosChi-parms[1]);
            dedcos=SIGN(1.0,e)*parms[0]*(CosChi-parms[1])*sin(Chi)/sqrt(c*(rab2-e*e/c));
            break;
          case PLANAR_INVERSION:
          case PLANAR_INVERSION2:
            // (1/2)*p_0*(1-cos(phi))
            // ===============================================
            // p_0/k_B [K]
            Chi=acos(CosChi);
            energy=parms[0]*(1.0-CosChi);
            dedcos=-SIGN(1.0,e)*parms[0]*sin(Chi)/sqrt(c*(rab2-e*e/c));
            break;
          case MM3_INVERSION:
            // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
            // =================================================================================================
            // p_0/k_B [mdyne A/rad^2]
            // p_1     [degrees]
            Chi=acos(CosChi);
            temp=RAD2DEG*(Chi-parms[1]);
            temp2=SQR(temp);
            energy=parms[0]*temp2*(1.0-0.014*temp+5.6e-5*temp2-7.0e-7*temp*temp2+2.2e-8*SQR(temp2));
            dedcos=-SIGN(1.0,e)*parms[0]*temp*RAD2DEG*(2.0-3.0*0.014*temp+4.0*5.6e-5*temp2-5.0*7.0e-7*temp*temp2+6.0*2.2e-8*SQR(temp2))/sqrt(c*(rab2-e*e/c));
            break;
          default:
            printf("Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
            exit(0);
            break;
        }
  
        // energy
        *Energy=energy;
  
        switch(Framework[CurrentSystem].InversionBendType[f1][i])
        {
          case HARMONIC_COSINE_INVERSION:
          case PLANAR_INVERSION:
          case HARMONIC_INVERSION:
            term=e/c;
            dccdia.x=0.0;
            dccdia.y=0.0;
            dccdia.z=0.0;
            dccdic.x=(Rbc.x*rbd2-Rbd.x*dot)*term;
            dccdic.y=(Rbc.y*rbd2-Rbd.y*dot)*term;
            dccdic.z=(Rbc.z*rbd2-Rbd.z*dot)*term;
            dccdid.x=(Rbd.x*rbc2-Rbc.x*dot)*term;
            dccdid.y=(Rbd.y*rbc2-Rbc.y*dot)*term;
            dccdid.z=(Rbd.z*rbc2-Rbc.z*dot)*term;
            break;
          case HARMONIC_INVERSION2:
          case HARMONIC_COSINE_INVERSION2:
          case PLANAR_INVERSION2:
          case MM3_INVERSION:
            term=e/c;
            dccdic.x=(Rac.x*rad2-Rad.x*dot)*term;
            dccdic.y=(Rac.y*rad2-Rad.y*dot)*term;
            dccdic.z=(Rac.z*rad2-Rad.z*dot)*term;
            dccdid.x=(Rad.x*rac2-Rac.x*dot)*term;
            dccdid.y=(Rad.y*rac2-Rac.y*dot)*term;
            dccdid.z=(Rad.z*rac2-Rac.z*dot)*term;
            dccdia.x=-(dccdic.x+dccdid.x);
            dccdia.y=-(dccdic.y+dccdid.y);
            dccdia.z=-(dccdic.z+dccdid.z);
            break;
          default:
            printf("Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
            exit(0);
            break;
        }
  
  
        term=e/rab2;
        deedia.x=Rbd.y*Rbc.z-Rbd.z*Rbc.y+Rab.x*term;
        deedia.y=Rbd.z*Rbc.x-Rbd.x*Rbc.z+Rab.y*term;
        deedia.z=Rbd.x*Rbc.y-Rbd.y*Rbc.x+Rab.z*term;
        deedic.x=Rab.y*Rbd.z-Rab.z*Rbd.y;
        deedic.y=Rab.z*Rbd.x-Rab.x*Rbd.z;
        deedic.z=Rab.x*Rbd.y-Rab.y*Rbd.x;
        deedid.x=Rbc.y*Rab.z-Rbc.z*Rab.y;
        deedid.y=Rbc.z*Rab.x-Rbc.x*Rab.z;
        deedid.z=Rbc.x*Rab.y-Rbc.y*Rab.x;
  
        fa->x=-dedcos*(dccdia.x+deedia.x);
        fa->y=-dedcos*(dccdia.y+deedia.y);
        fa->z=-dedcos*(dccdia.z+deedia.z);
        fc->x=-dedcos*(dccdic.x+deedic.x);
        fc->y=-dedcos*(dccdic.y+deedic.y);
        fc->z=-dedcos*(dccdic.z+deedic.z);
        fd->x=-dedcos*(dccdid.x+deedid.x);
        fd->y=-dedcos*(dccdid.y+deedid.y);
        fd->z=-dedcos*(dccdid.z+deedid.z);
  
  fb->x=-(fa->x+fc->x+fd->x);
  fb->y=-(fa->y+fc->y+fd->y);
  fb->z=-(fa->z+fc->z+fd->z);

        strain_derivative->ax=Rab.x*fa->x+Rbc.x*fc->x+Rbd.x*fd->x;
        strain_derivative->ay=Rab.x*fa->y+Rbc.x*fc->y+Rbd.x*fd->y;
        strain_derivative->az=Rab.x*fa->z+Rbc.x*fc->z+Rbd.x*fd->z;
  
        strain_derivative->bx=Rab.y*fa->x+Rbc.y*fc->x+Rbd.y*fd->x;
        strain_derivative->by=Rab.y*fa->y+Rbc.y*fc->y+Rbd.y*fd->y;
        strain_derivative->bz=Rab.y*fa->z+Rbc.y*fc->z+Rbd.y*fd->z;
  
        strain_derivative->cx=Rab.z*fa->x+Rbc.z*fc->x+Rbd.z*fd->x;
        strain_derivative->cy=Rab.z*fa->y+Rbc.z*fc->y+Rbd.z*fd->y;
        strain_derivative->cz=Rab.z*fa->z+Rbc.z*fc->z+Rbd.z*fd->z;

}


//void CalculateFrameworkInversionBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX3x3 *StrainDerivativeTensor,REAL_MATRIX HessianMatrix,REAL_MATRIX CrossTerm)

void ComputeFrameworkInversionBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivativeTensor,int ComputeGradient,int ComputeHessian)
{
  int i,f1;
  int A,B,C,D;
  const REAL Delta=1e-7;
  VECTOR posA0,posB0,posC0,posD0;
  VECTOR posA,posB,posC,posD;
  VECTOR fa0,fb0,fc0,fd0,fa[4],fb[4],fc[4],fd[4];
  int index_i,index_j,index_k,index_l;
  REAL U;
  REAL_MATRIX3x3 strain_derivative;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfInversionBends[f1];i++)
    {
      A=Framework[CurrentSystem].InversionBends[f1][i].A;
      B=Framework[CurrentSystem].InversionBends[f1][i].B;
      C=Framework[CurrentSystem].InversionBends[f1][i].C;
      D=Framework[CurrentSystem].InversionBends[f1][i].D;

      index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
      index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;
      index_k=Framework[CurrentSystem].Atoms[f1][C].HessianIndex;
      index_l=Framework[CurrentSystem].Atoms[f1][D].HessianIndex;

      posA0=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB0=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC0=Framework[CurrentSystem].Atoms[f1][C].Position;
      posD0=Framework[CurrentSystem].Atoms[f1][D].Position;

      CalculateFrameworkInversionBendForces(f1,i,posA0,posB0,posC0,posD0,&U,&fa0,&fb0,&fc0,&fd0,&strain_derivative);

      *Energy+=U;

      if((index_i<0)&&(index_j<0)&&(index_k<0)&&(index_l<0)) continue;

      StrainDerivativeTensor->ax+=strain_derivative.ax;
      StrainDerivativeTensor->ay+=strain_derivative.ay;
      StrainDerivativeTensor->az+=strain_derivative.az;

      StrainDerivativeTensor->bx+=strain_derivative.bx;
      StrainDerivativeTensor->by+=strain_derivative.by;
      StrainDerivativeTensor->bz+=strain_derivative.bz;

      StrainDerivativeTensor->cx+=strain_derivative.cx;
      StrainDerivativeTensor->cy+=strain_derivative.cy;
      StrainDerivativeTensor->cz+=strain_derivative.cz;
      
      // add contribution to the first derivatives
      if(ComputeGradient)
      {
        if(index_i>=0)
        {
          Gradient[index_i]-=fa0.x;
          Gradient[index_i+1]-=fa0.y;
          Gradient[index_i+2]-=fa0.z;
        }

        if(index_j>=0)
        {
          Gradient[index_j]-=fb0.x;
          Gradient[index_j+1]-=fb0.y;
          Gradient[index_j+2]-=fb0.z;
        }

        if(index_k>=0)
        {
          Gradient[index_k]-=fc0.x;
          Gradient[index_k+1]-=fc0.y;
          Gradient[index_k+2]-=fc0.z;
        }

        if(index_l>=0)
        {
          Gradient[index_l]-=fd0.x;
          Gradient[index_l+1]-=fd0.y;
          Gradient[index_l+2]-=fd0.z;
        }
      }

      if(ComputeHessian)
      {
      // Atom A

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.x+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.x+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.x-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.x-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_i>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_i][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_i>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_i][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_i>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_i][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_i>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_i][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.y+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.y+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.y-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.y-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_i>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_i+1][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i+1][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i+1][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_i>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_i+1][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i+1][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i+1][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_i>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_i+1][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i+1][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i+1][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_i>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_i+1][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i+1][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i+1][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.z+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.z+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.z-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posA.z-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_i>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_i+2][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i+2][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i+2][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_i>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_i+2][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i+2][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i+2][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_i>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_i+2][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i+2][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i+2][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_i>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_i+2][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_i+2][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_i+2][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      // Atom B

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.x+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.x+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.x-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.x-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_j>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_j][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_j>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_j][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_j>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_j][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_j>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_j][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.y+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.y+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.y-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.y-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_j>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_j+1][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j+1][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j+1][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_j>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_j+1][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j+1][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j+1][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_j>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_j+1][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j+1][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j+1][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_j>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_j+1][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j+1][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j+1][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.z+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.z+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.z-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posB.z-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_j>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_j+2][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j+2][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j+2][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_j>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_j+2][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j+2][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j+2][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_j>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_j+2][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j+2][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j+2][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_j>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_j+2][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_j+2][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_j+2][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      // Atom C

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.x+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.x+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.x-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.x-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_k>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_k][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_k>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_k][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_k>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_k][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_k>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_k][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.y+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.y+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.y-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.y-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_k>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_k+1][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k+1][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k+1][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_k>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_k+1][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k+1][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k+1][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_k>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_k+1][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k+1][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k+1][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_k>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_k+1][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k+1][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k+1][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.z+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.z+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.z-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posC.z-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_k>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_k+2][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k+2][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k+2][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_k>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_k+2][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k+2][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k+2][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_k>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_k+2][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k+2][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k+2][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_k>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_k+2][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_k+2][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_k+2][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      // Atom D

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.x+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.x+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.x-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.x-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_l>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_l][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_l>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_l][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_l>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_l][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_l>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_l][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.y+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.y+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.y-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.y-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_l>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_l+1][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l+1][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l+1][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_l>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_l+1][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l+1][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l+1][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_l>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_l+1][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l+1][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l+1][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_l>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_l+1][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l+1][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l+1][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.z+=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[0],&fb[0],&fc[0],&fd[0],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.z+=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[1],&fb[1],&fc[1],&fd[1],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.z-=0.5*Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[2],&fb[2],&fc[2],&fd[2],&strain_derivative);

      posA=posA0; posB=posB0; posC=posC0; posD=posD0;
      posD.z-=Delta;
      CalculateFrameworkInversionBendForces(f1,i,posA,posB,posC,posD,&U,&fa[3],&fb[3],&fc[3],&fd[3],&strain_derivative);

      if((index_l>=0)&&(index_i>=0))
      {
        HessianMatrix.element[index_l+2][index_i]-=(-fa[0].x+8.0*fa[1].x-8.0*fa[2].x+fa[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l+2][index_i+1]-=(-fa[0].y+8.0*fa[1].y-8.0*fa[2].y+fa[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l+2][index_i+2]-=(-fa[0].z+8.0*fa[1].z-8.0*fa[2].z+fa[3].z)/(6.0*Delta);
      }

      if((index_l>=0)&&(index_j>=0))
      {
        HessianMatrix.element[index_l+2][index_j]-=(-fb[0].x+8.0*fb[1].x-8.0*fb[2].x+fb[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l+2][index_j+1]-=(-fb[0].y+8.0*fb[1].y-8.0*fb[2].y+fb[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l+2][index_j+2]-=(-fb[0].z+8.0*fb[1].z-8.0*fb[2].z+fb[3].z)/(6.0*Delta);
      }

      if((index_l>=0)&&(index_k>=0))
      {
        HessianMatrix.element[index_l+2][index_k]-=(-fc[0].x+8.0*fc[1].x-8.0*fc[2].x+fc[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l+2][index_k+1]-=(-fc[0].y+8.0*fc[1].y-8.0*fc[2].y+fc[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l+2][index_k+2]-=(-fc[0].z+8.0*fc[1].z-8.0*fc[2].z+fc[3].z)/(6.0*Delta);
      }

      if((index_l>=0)&&(index_l>=0))
      {
        HessianMatrix.element[index_l+2][index_l]-=(-fd[0].x+8.0*fd[1].x-8.0*fd[2].x+fd[3].x)/(6.0*Delta);
        HessianMatrix.element[index_l+2][index_l+1]-=(-fd[0].y+8.0*fd[1].y-8.0*fd[2].y+fd[3].y)/(6.0*Delta);
        HessianMatrix.element[index_l+2][index_l+2]-=(-fd[0].z+8.0*fd[1].z-8.0*fd[2].z+fd[3].z)/(6.0*Delta);
      }
      }
    }
  }
}

