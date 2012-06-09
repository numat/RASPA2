/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'inter_hessian.c' is part of RASPA.

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
#include <sys/stat.h>
#include "constants.h"
#include "utils.h"
#include "simulation.h"
#include "potentials.h"
#include "output.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "grids.h"
#include "ewald.h"
#include "inter_energy.h"
#include "inter_force.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "spectra.h"
#include "rigid.h"
#include "minimization.h"

// Hessian: Center of mass - Center of mass
// ========================================
static inline void HessianAtomicPositionPosition(COMPLEX phase_factor,COMPLEX_MATRIX HessianMatrix,int index_i,int index_j,
                                                 REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor)
{
  REAL_MATRIX3x3 Hessian;

  // the first and second term of Eq. S54 of Ref. Dubbeldam, Krishna, Snurr, 2009
  // ============================================================================
  // I==J: f_2 r_ij,alpha r_ij,beta + f_1 delta(alpha,beta)
  // I!=J: -f_2 r_ij,alpha r_ij,beta - f_1 delta(alpha,beta)

  Hessian.ax=f2*dr.x*dr.x+f1; Hessian.bx=f2*dr.y*dr.x;    Hessian.cx=f2*dr.z*dr.x;
  Hessian.ay=f2*dr.x*dr.y;    Hessian.by=f2*dr.y*dr.y+f1; Hessian.cy=f2*dr.z*dr.y;
  Hessian.az=f2*dr.x*dr.z;    Hessian.bz=f2*dr.y*dr.z;    Hessian.cz=f2*dr.z*dr.z+f1;

  // case [I,I]: Center of mass - Center of mass
  if(index_i>=0)
  { 
    switch(Dimension)
    {
      case 3:
        HessianMatrix.element[index_i][index_i+2].re+=ReplicaFactor*Hessian.az;
        HessianMatrix.element[index_i+1][index_i+2].re+=ReplicaFactor*Hessian.bz;
        HessianMatrix.element[index_i+2][index_i+2].re+=ReplicaFactor*Hessian.cz;
      case 2:
        HessianMatrix.element[index_i][index_i+1].re+=ReplicaFactor*Hessian.ay;
        HessianMatrix.element[index_i+1][index_i+1].re+=ReplicaFactor*Hessian.by;
      case 1:
        HessianMatrix.element[index_i][index_i].re+=ReplicaFactor*Hessian.ax;
        break;
    }
  }

  // case [J,J]: Center of mass - Center of mass
  if(index_j>=0)
  { 
    switch(Dimension)
    {
      case 3:
        HessianMatrix.element[index_j][index_j+2].re+=ReplicaFactor*Hessian.az;
        HessianMatrix.element[index_j+1][index_j+2].re+=ReplicaFactor*Hessian.bz;
        HessianMatrix.element[index_j+2][index_j+2].re+=ReplicaFactor*Hessian.cz;
      case 2:
        HessianMatrix.element[index_j][index_j+1].re+=ReplicaFactor*Hessian.ay;
        HessianMatrix.element[index_j+1][index_j+1].re+=ReplicaFactor*Hessian.by;
      case 1:
        HessianMatrix.element[index_j][index_j].re+=ReplicaFactor*Hessian.ax;
        break;
    }
  }

  // case [I,J]: Center of mass - Center of mass
  if((index_i>=0)&&(index_j>=0))
  {
    switch(Dimension)
    {
      case 3:
        HessianMatrix.element[index_i][index_j+2].re-=phase_factor.re*Hessian.az;
        HessianMatrix.element[index_i][index_j+2].im-=phase_factor.im*Hessian.az;
        HessianMatrix.element[index_i+1][index_j+2].re-=phase_factor.re*Hessian.bz;
        HessianMatrix.element[index_i+1][index_j+2].im-=phase_factor.im*Hessian.bz;
        HessianMatrix.element[index_i+2][index_j].re-=phase_factor.re*Hessian.az;
        HessianMatrix.element[index_i+2][index_j].im-=phase_factor.im*Hessian.az;
        HessianMatrix.element[index_i+2][index_j+1].re-=phase_factor.re*Hessian.bz;
        HessianMatrix.element[index_i+2][index_j+1].im-=phase_factor.im*Hessian.bz;
        HessianMatrix.element[index_i+2][index_j+2].re-=phase_factor.re*Hessian.cz;
        HessianMatrix.element[index_i+2][index_j+2].im-=phase_factor.im*Hessian.cz;
      case 2:
        HessianMatrix.element[index_i][index_j+1].re-=phase_factor.re*Hessian.ay;
        HessianMatrix.element[index_i][index_j+1].im-=phase_factor.im*Hessian.ay;
        HessianMatrix.element[index_i+1][index_j].re-=phase_factor.re*Hessian.ay;
        HessianMatrix.element[index_i+1][index_j].im-=phase_factor.im*Hessian.ay;
        HessianMatrix.element[index_i+1][index_j+1].re-=phase_factor.re*Hessian.by;
        HessianMatrix.element[index_i+1][index_j+1].im-=phase_factor.im*Hessian.by;
      case 1:
        HessianMatrix.element[index_i][index_j].re-=phase_factor.re*Hessian.ax;
        HessianMatrix.element[index_i][index_j].im-=phase_factor.im*Hessian.ax;
        break;
    }
  }
}


// Hessian: Center of mass - Orientation
// =====================================
static inline void HessianCenterOfMassOrientation(COMPLEX phase_factor,COMPLEX_MATRIX HessianMatrix,int index_i,int index_i2,int index_j,int index_j2,
              int index1,int index2,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor)
{
  REAL_MATRIX3x3 Hessian;
  VECTOR veci1,veci2,veci3;
  VECTOR vecj1,vecj2,vecj3;

  // the first and second term of Eq. S56 of Ref. Dubbeldam, Krishna, Snurr, 2009
  // ============================================================================
  // I==J: f_2 (r_ij . R^I_beta r_i^0) r_ij,alpha + f1 [R^I_beta r_i^0]_alpha
  // I!=J: -f_2 (r_ij . R^J_beta r_j^0) r_ij,alpha - f1 [R^J_beta r_j^0]_alpha

  // case [I,I]: Center of mass I - Orientation J
  if((index_i>=0)&&(index_i2>=0))
  {
    veci1=DVecX[index1];
    veci2=DVecY[index1];
    veci3=DVecZ[index1];

    vecj1=DVecX[index2];
    vecj2=DVecY[index2];
    vecj3=DVecZ[index2];

    Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
    Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
    Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

    Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
    Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
    Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

    Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
    Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
    Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

    switch(Dimension)
    {
      case 3:
        HessianMatrix.element[index_i][index_i2+2].re+=ReplicaFactor*Hessian.az;    // xz
        HessianMatrix.element[index_i+1][index_i2+2].re+=ReplicaFactor*Hessian.bz;  // yz
        HessianMatrix.element[index_i+2][index_i2].re+=ReplicaFactor*Hessian.cx;    // zx
        HessianMatrix.element[index_i+2][index_i2+1].re+=ReplicaFactor*Hessian.cy;  // zy
        HessianMatrix.element[index_i+2][index_i2+2].re+=ReplicaFactor*Hessian.cz;  // zz
      case 2:
        HessianMatrix.element[index_i][index_i2+1].re+=ReplicaFactor*Hessian.ay;    // xy
        HessianMatrix.element[index_i+1][index_i2].re+=ReplicaFactor*Hessian.bx;    // yx
        HessianMatrix.element[index_i+1][index_i2+1].re+=ReplicaFactor*Hessian.by;  // yy
      case 1:
        HessianMatrix.element[index_i][index_i2].re+=ReplicaFactor*Hessian.ax;      // xx
        break;
    }
  }

  // case [I,J]: Orientation I - Center of mass J
  if((index_i2>=0)&&(index_j>=0))
  {
    veci1=DVecX[index1];
    veci2=DVecY[index1];
    veci3=DVecZ[index1];

    vecj1=DVecX[index2];
    vecj2=DVecY[index2];
    vecj3=DVecZ[index2];

    Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.x+f1*veci1.x;
    Hessian.ay=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.x+f1*veci2.x;
    Hessian.az=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.x+f1*veci3.x;

    Hessian.bx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.y+f1*veci1.y;
    Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.y+f1*veci2.y;
    Hessian.bz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.y+f1*veci3.y;

    Hessian.cx=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*dr.z+f1*veci1.z;
    Hessian.cy=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*dr.z+f1*veci2.z;
    Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*dr.z+f1*veci3.z;

    switch(Dimension)
    {
      case 3:
        HessianMatrix.element[index_i2+2][index_j].re-=phase_factor.re*Hessian.az;    // zx
        HessianMatrix.element[index_i2+2][index_j].im-=phase_factor.im*Hessian.az;    // zx
        HessianMatrix.element[index_i2+2][index_j+1].re-=phase_factor.re*Hessian.bz;  // zy
        HessianMatrix.element[index_i2+2][index_j+1].im-=phase_factor.im*Hessian.bz;  // zy
        HessianMatrix.element[index_i2][index_j+2].re-=phase_factor.re*Hessian.cx;    // xz
        HessianMatrix.element[index_i2][index_j+2].im-=phase_factor.im*Hessian.cx;    // xz
        HessianMatrix.element[index_i2+1][index_j+2].re-=phase_factor.re*Hessian.cy;  // yz
        HessianMatrix.element[index_i2+1][index_j+2].im-=phase_factor.im*Hessian.cy;  // yz
        HessianMatrix.element[index_i2+2][index_j+2].re-=phase_factor.re*Hessian.cz;  // zz
        HessianMatrix.element[index_i2+2][index_j+2].im-=phase_factor.im*Hessian.cz;  // zz
      case 2:
        HessianMatrix.element[index_i2+1][index_j].re-=phase_factor.re*Hessian.ay;    // yx
        HessianMatrix.element[index_i2+1][index_j].im-=phase_factor.im*Hessian.ay;    // yx
        HessianMatrix.element[index_i2][index_j+1].re-=phase_factor.re*Hessian.bx;    // xy
        HessianMatrix.element[index_i2][index_j+1].im-=phase_factor.im*Hessian.bx;    // xy
        HessianMatrix.element[index_i2+1][index_j+1].re-=phase_factor.re*Hessian.by;  // yy
        HessianMatrix.element[index_i2+1][index_j+1].im-=phase_factor.im*Hessian.by;  // yy
      case 1:
        HessianMatrix.element[index_i2][index_j].re-=phase_factor.re*Hessian.ax;      // xx
        HessianMatrix.element[index_i2][index_j].im-=phase_factor.im*Hessian.ax;      // xx
        break;
    }
  }

  // case [J,J]: Center of mass J - Orientation J
  if((index_j>=0)&&(index_j2>=0))
  {
    veci1=DVecX[index1];
    veci2=DVecY[index1];
    veci3=DVecZ[index1];

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

    switch(Dimension)
    {
      case 3:
        HessianMatrix.element[index_j][index_j2+2].re+=ReplicaFactor*Hessian.az;
        HessianMatrix.element[index_j+1][index_j2+2].re+=ReplicaFactor*Hessian.bz;
        HessianMatrix.element[index_j+2][index_j2].re+=ReplicaFactor*Hessian.cx;
        HessianMatrix.element[index_j+2][index_j2+1].re+=ReplicaFactor*Hessian.cy;
        HessianMatrix.element[index_j+2][index_j2+2].re+=ReplicaFactor*Hessian.cz;
      case 2:
        HessianMatrix.element[index_j][index_j2+1].re+=ReplicaFactor*Hessian.ay;
        HessianMatrix.element[index_j+1][index_j2].re+=ReplicaFactor*Hessian.bx;
        HessianMatrix.element[index_j+1][index_j2+1].re+=ReplicaFactor*Hessian.by;
      case 1:
        HessianMatrix.element[index_j][index_j2].re+=ReplicaFactor*Hessian.ax;
        break;
    }
  }

  // case [I,J]: Center of mass I - Orientation J
  if((index_i>=0)&&(index_j2>=0))
  {
    veci1=DVecX[index1];
    veci2=DVecY[index1];
    veci3=DVecZ[index1];

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

    switch(Dimension)
    {
      case 3:
        HessianMatrix.element[index_i][index_j2+2].re-=phase_factor.re*Hessian.az;    // xz
        HessianMatrix.element[index_i][index_j2+2].im-=phase_factor.im*Hessian.az;    // xz
        HessianMatrix.element[index_i+1][index_j2+2].re-=phase_factor.re*Hessian.bz;  // yz
        HessianMatrix.element[index_i+1][index_j2+2].im-=phase_factor.im*Hessian.bz;  // yz
        HessianMatrix.element[index_i+2][index_j2].re-=phase_factor.re*Hessian.cx;    // zx
        HessianMatrix.element[index_i+2][index_j2].im-=phase_factor.im*Hessian.cx;    // zx
        HessianMatrix.element[index_i+2][index_j2+1].re-=phase_factor.re*Hessian.cy;  // zy
        HessianMatrix.element[index_i+2][index_j2+1].im-=phase_factor.im*Hessian.cy;  // zy
        HessianMatrix.element[index_i+2][index_j2+2].re-=phase_factor.re*Hessian.cz;  // zz
        HessianMatrix.element[index_i+2][index_j2+2].im-=phase_factor.im*Hessian.cz;  // zz
      case 2:
        HessianMatrix.element[index_i][index_j2+1].re-=phase_factor.re*Hessian.ay;    // xy
        HessianMatrix.element[index_i][index_j2+1].im-=phase_factor.im*Hessian.ay;    // xy
        HessianMatrix.element[index_i+1][index_j2].re-=phase_factor.re*Hessian.bx;    // yx
        HessianMatrix.element[index_i+1][index_j2].im-=phase_factor.im*Hessian.bx;    // yx
        HessianMatrix.element[index_i+1][index_j2+1].re-=phase_factor.re*Hessian.by;  // yy
        HessianMatrix.element[index_i+1][index_j2+1].im-=phase_factor.im*Hessian.by;  // yy
      case 1:
        HessianMatrix.element[index_i][index_j2].re-=phase_factor.re*Hessian.ax;      // xx
        HessianMatrix.element[index_i][index_j2].im-=phase_factor.im*Hessian.ax;      // xx
        break;
    }
  }
}


// Hessian: Orientation - Orientation
// ==================================
static inline void HessianOrientationOrientation(COMPLEX phase_factor,COMPLEX_MATRIX HessianMatrix,int index_i,int index_i2,int index_j,int index_j2,
             int index1,int index2,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor)
{
  REAL_MATRIX3x3 Hessian;
  VECTOR veci1,veci2,veci3;
  VECTOR vecj1,vecj2,vecj3;
  VECTOR DDvecIAX,DDvecIBY,DDvecICZ,DDvecIAY,DDvecIAZ,DDvecIBZ;
  VECTOR DDvecJAX,DDvecJBY,DDvecJCZ,DDvecJAY,DDvecJAZ,DDvecJBZ;

  // the first, second and third terms of Eq. S55 of Ref. Dubbeldam, Krishna, Snurr, 2009
  // ====================================================================================
  // I==J: f_2 (r_ij . R^I_alpha r_i^0) (r_ij . R^I_beta r_i^0)+
  //       f_1 (R^I_alpha r_i^0).(R^I_beta r_i^0)+
  //       f_1 (r_ij . R^I_alpha,beta r_i^0)
  // I!=J: -f_2 (r_ij . R^I_alpha r_i^0) (r_ij . R^J_beta r_j^0)-
  //       f_1 (R^I_alpha r_i^0).(R^J_beta r_j^0)

  // case [I,I]: Orientation I - Orientation I
  if(index_i2>=0)
  {
    veci1=DVecX[index1]; vecj1=DVecX[index2];
    veci2=DVecY[index1]; vecj2=DVecY[index2];
    veci3=DVecZ[index1]; vecj3=DVecZ[index2];

    DDvecIAX=DDVecAX[index1]; DDvecJAX=DDVecAX[index2];
    DDvecIBY=DDVecBY[index1]; DDvecJBY=DDVecBY[index2];
    DDvecICZ=DDVecCZ[index1]; DDvecJCZ=DDVecCZ[index2];
    DDvecIAY=DDVecAY[index1]; DDvecJAY=DDVecAY[index2];
    DDvecIAZ=DDVecAZ[index1]; DDvecJAZ=DDVecAZ[index2];
    DDvecIBZ=DDVecBZ[index1]; DDvecJBZ=DDVecBZ[index2];

    Hessian.ax=f1*(dr.x*DDvecIAX.x+dr.y*DDvecIAX.y+dr.z*DDvecIAX.z)+
               f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)+
               f1*(veci1.x*veci1.x+veci1.y*veci1.y+veci1.z*veci1.z);

    Hessian.by=f1*(dr.x*DDvecIBY.x+dr.y*DDvecIBY.y+dr.z*DDvecIBY.z)+
               f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)+
               f1*(veci2.x*veci2.x+veci2.y*veci2.y+veci2.z*veci2.z);

    Hessian.cz=f1*(dr.x*DDvecICZ.x+dr.y*DDvecICZ.y+dr.z*DDvecICZ.z)+
               f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)+
               f1*(veci3.x*veci3.x+veci3.y*veci3.y+veci3.z*veci3.z);


    Hessian.ay=f1*(dr.x*DDvecIAY.x+dr.y*DDvecIAY.y+dr.z*DDvecIAY.z)+
               f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)+
               f1*(veci1.x*veci2.x+veci1.y*veci2.y+veci1.z*veci2.z);

    Hessian.az=f1*(dr.x*DDvecIAZ.x+dr.y*DDvecIAZ.y+dr.z*DDvecIAZ.z)+
               f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)+
               f1*(veci1.x*veci3.x+veci1.y*veci3.y+veci1.z*veci3.z);

    Hessian.bz=f1*(dr.x*DDvecIBZ.x+dr.y*DDvecIBZ.y+dr.z*DDvecIBZ.z)+
               f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)+
               f1*(veci2.x*veci3.x+veci2.y*veci3.y+veci2.z*veci3.z);

    switch(Dimension)
    {
      case 3:
        HessianMatrix.element[index_i2+2][index_i2+2].re+=ReplicaFactor*Hessian.cz;
        HessianMatrix.element[index_i2][index_i2+2].re+=ReplicaFactor*Hessian.az;
        HessianMatrix.element[index_i2+1][index_i2+2].re+=ReplicaFactor*Hessian.bz;
      case 2:
        HessianMatrix.element[index_i2+1][index_i2+1].re+=ReplicaFactor*Hessian.by;
        HessianMatrix.element[index_i2][index_i2+1].re+=ReplicaFactor*Hessian.ay;
      case 1:
        HessianMatrix.element[index_i2][index_i2].re+=ReplicaFactor*Hessian.ax;
        break;
    }
  }

  // case [J,J]: Orientation J - Orientation J
  if(index_j2>=0)
  {
    veci1=DVecX[index1]; vecj1=DVecX[index2];
    veci2=DVecY[index1]; vecj2=DVecY[index2];
    veci3=DVecZ[index1]; vecj3=DVecZ[index2];

    DDvecIAX=DDVecAX[index1]; DDvecJAX=DDVecAX[index2];
    DDvecIBY=DDVecBY[index1]; DDvecJBY=DDVecBY[index2];
    DDvecICZ=DDVecCZ[index1]; DDvecJCZ=DDVecCZ[index2];
    DDvecIAY=DDVecAY[index1]; DDvecJAY=DDVecAY[index2];
    DDvecIAZ=DDVecAZ[index1]; DDvecJAZ=DDVecAZ[index2];
    DDvecIBZ=DDVecBZ[index1]; DDvecJBZ=DDVecBZ[index2];

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

    switch(Dimension)
    {
      case 3:
        HessianMatrix.element[index_j2+2][index_j2+2].re+=ReplicaFactor*Hessian.cz;
        HessianMatrix.element[index_j2][index_j2+2].re+=ReplicaFactor*Hessian.az;
        HessianMatrix.element[index_j2+1][index_j2+2].re+=ReplicaFactor*Hessian.bz;
      case 2:
        HessianMatrix.element[index_j2+1][index_j2+1].re+=ReplicaFactor*Hessian.by;
        HessianMatrix.element[index_j2][index_j2+1].re+=ReplicaFactor*Hessian.ay;
      case 1:
        HessianMatrix.element[index_j2][index_j2].re+=ReplicaFactor*Hessian.ax;
        break;
    }
  }

  // case [I,J]: Orientation I - Orientation J
  if((index_i2>=0)&&(index_j2>=0))
  {
    veci1=DVecX[index1]; vecj1=DVecX[index2];
    veci2=DVecY[index1]; vecj2=DVecY[index2];
    veci3=DVecZ[index1]; vecj3=DVecZ[index2];

    DDvecIAX=DDVecAX[index1]; DDvecJAX=DDVecAX[index2];
    DDvecIBY=DDVecBY[index1]; DDvecJBY=DDVecBY[index2];
    DDvecICZ=DDVecCZ[index1]; DDvecJCZ=DDVecCZ[index2];
    DDvecIAY=DDVecAY[index1]; DDvecJAY=DDVecAY[index2];
    DDvecIAZ=DDVecAZ[index1]; DDvecJAZ=DDVecAZ[index2];
    DDvecIBZ=DDVecBZ[index1]; DDvecJBZ=DDVecBZ[index2];

    Hessian.ax=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)+
               f1*(veci1.x*vecj1.x+veci1.y*vecj1.y+veci1.z*vecj1.z);
    Hessian.by=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)+
               f1*(veci2.x*vecj2.x+veci2.y*vecj2.y+veci2.z*vecj2.z);
    Hessian.cz=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
               f1*(veci3.x*vecj3.x+veci3.y*vecj3.y+veci3.z*vecj3.z);
    Hessian.ay=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)+
               f1*(veci1.x*vecj2.x+veci1.y*vecj2.y+veci1.z*vecj2.z);
    Hessian.bx=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)+
               f1*(veci2.x*vecj1.x+veci2.y*vecj1.y+veci2.z*vecj1.z);
    Hessian.az=f2*(dr.x*veci1.x+dr.y*veci1.y+dr.z*veci1.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
               f1*(veci1.x*vecj3.x+veci1.y*vecj3.y+veci1.z*vecj3.z);
    Hessian.cx=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*(dr.x*vecj1.x+dr.y*vecj1.y+dr.z*vecj1.z)+
               f1*(veci3.x*vecj1.x+veci3.y*vecj1.y+veci3.z*vecj1.z);
    Hessian.bz=f2*(dr.x*veci2.x+dr.y*veci2.y+dr.z*veci2.z)*(dr.x*vecj3.x+dr.y*vecj3.y+dr.z*vecj3.z)+
               f1*(veci2.x*vecj3.x+veci2.y*vecj3.y+veci2.z*vecj3.z);
    Hessian.cy=f2*(dr.x*veci3.x+dr.y*veci3.y+dr.z*veci3.z)*(dr.x*vecj2.x+dr.y*vecj2.y+dr.z*vecj2.z)+
               f1*(veci3.x*vecj2.x+veci3.y*vecj2.y+veci3.z*vecj2.z);

    switch(Dimension)
    {
      case 3:
        HessianMatrix.element[index_i2+2][index_j2+2].re-=phase_factor.re*Hessian.cz;
        HessianMatrix.element[index_i2+2][index_j2+2].im-=phase_factor.im*Hessian.cz;
        HessianMatrix.element[index_i2][index_j2+2].re-=phase_factor.re*Hessian.az;
        HessianMatrix.element[index_i2][index_j2+2].im-=phase_factor.im*Hessian.az;
        HessianMatrix.element[index_i2+2][index_j2].re-=phase_factor.re*Hessian.cx;
        HessianMatrix.element[index_i2+2][index_j2].im-=phase_factor.im*Hessian.cx;
        HessianMatrix.element[index_i2+1][index_j2+2].re-=phase_factor.re*Hessian.bz;
        HessianMatrix.element[index_i2+1][index_j2+2].im-=phase_factor.im*Hessian.bz;
        HessianMatrix.element[index_i2+2][index_j2+1].re-=phase_factor.re*Hessian.cy;
        HessianMatrix.element[index_i2+2][index_j2+1].im-=phase_factor.im*Hessian.cy;
      case 2:
        HessianMatrix.element[index_i2+1][index_j2+1].re-=phase_factor.re*Hessian.by;
        HessianMatrix.element[index_i2+1][index_j2+1].im-=phase_factor.im*Hessian.by;
        HessianMatrix.element[index_i2][index_j2+1].re-=phase_factor.re*Hessian.ay;
        HessianMatrix.element[index_i2][index_j2+1].im-=phase_factor.im*Hessian.ay;
        HessianMatrix.element[index_i2+1][index_j2].re-=phase_factor.re*Hessian.bx;
        HessianMatrix.element[index_i2+1][index_j2].im-=phase_factor.im*Hessian.bx;
      case 1:
        HessianMatrix.element[index_i2][index_j2].re-=phase_factor.re*Hessian.ax;
        HessianMatrix.element[index_i2][index_j2].im-=phase_factor.im*Hessian.ax;
        break;
    }
  }
}



static inline void GradientStrain(REAL *Gradient,REAL f1,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      switch(Dimension)
      {
        case 3:
          Gradient[n]+=f1*dr.z*dr.z;
        case 2:
          Gradient[n]+=f1*dr.y*dr.y;
        case 1:
          Gradient[n]+=f1*dr.x*dr.x;
          break;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=f1*dr.z*dr.z/Dimension;
            case 2:
              Gradient[n+1]+=f1*dr.y*dr.y/Dimension;
            case 1:
              Gradient[n]+=f1*dr.x*dr.x/Dimension;
              break;
          }
          break;
        case ANISOTROPIC:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=f1*dr.z*dr.z;
            case 2:
              Gradient[n+1]+=f1*dr.y*dr.y;
            case 1:
              Gradient[n]+=f1*dr.x*dr.x;
              break;
          }
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=f1*dr.x*dr.z;
              Gradient[n+4]+=f1*dr.y*dr.z;
              Gradient[n+5]+=f1*dr.z*dr.z;
            case 2:
              Gradient[n+1]+=f1*dr.x*dr.y;
              Gradient[n+Dimension]+=f1*dr.y*dr.y;
            case 1:
             Gradient[n]+=f1*dr.x*dr.x;
             break;
          }
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]+=f1*dr.x*dr.x;
              Gradient[n+1]+=f1*dr.y*dr.y;
              Gradient[n+2]+=f1*dr.y*dr.z;
              Gradient[n+3]+=f1*dr.z*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]+=f1*dr.x*dr.x;
              Gradient[n+1]+=f1*dr.x*dr.y;
              Gradient[n+2]+=f1*dr.y*dr.y;
              Gradient[n+3]+=f1*dr.z*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
            default:
              Gradient[n]+=f1*dr.x*dr.x;
              Gradient[n+1]+=f1*dr.x*dr.z;
              Gradient[n+2]+=f1*dr.y*dr.y;
              Gradient[n+3]+=f1*dr.z*dr.z;
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

static inline void GradientStrainI(REAL *Gradient,REAL f1,VECTOR dr,VECTOR posA,VECTOR comA)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      switch(Dimension)
      {
        case 3:
          Gradient[n]-=f1*(posA.z-comA.z)*dr.z;
        case 2:
          Gradient[n]-=f1*(posA.y-comA.y)*dr.y;
        case 1:
          Gradient[n]-=f1*(posA.x-comA.x)*dr.x;
          break;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]-=(posA.z-comA.z)*f1*dr.z;
            case 2:
              Gradient[n+1]-=(posA.y-comA.y)*f1*dr.y;
            case 1:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              break;
          }
          break;
        case REGULAR:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]-=0.5*f1*((posA.z-comA.z)*dr.x+(posA.x-comA.x)*dr.z);
              Gradient[n+4]-=0.5*f1*((posA.z-comA.z)*dr.y+(posA.y-comA.y)*dr.z);
              Gradient[n+5]-=f1*(posA.z-comA.z)*dr.z;
            case 2:
              Gradient[n+1]-=0.5*f1*((posA.y-comA.y)*dr.x+(posA.x-comA.x)*dr.y);
              Gradient[n+3]-=f1*(posA.y-comA.y)*dr.y;
            case 1:
              Gradient[n]-=f1*(posA.x-comA.x)*dr.x;
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]-=f1*(posA.z-comA.z)*dr.x;
              Gradient[n+4]-=f1*(posA.z-comA.z)*dr.y;
              Gradient[n+5]-=f1*(posA.z-comA.z)*dr.z;
            case 2:
              Gradient[n+1]-=f1*(posA.y-comA.y)*dr.x;
              Gradient[n+Dimension]-=f1*(posA.y-comA.y)*dr.y;
            case 1:
              Gradient[n]-=f1*(posA.x-comA.x)*dr.x;
              break;
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+2]-=0.5*f1*((posA.z-comA.z)*dr.y+(posA.y-comA.y)*dr.z);
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=0.5*f1*((posA.y-comA.y)*dr.x+(posA.x-comA.x)*dr.y);
              Gradient[n+2]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=0.5*f1*((posA.z-comA.z)*dr.x+(posA.x-comA.x)*dr.z);
              Gradient[n+2]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            default:
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+2]-=f1*(posA.z-comA.z)*dr.y;
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=f1*(posA.z-comA.z)*dr.x;
              Gradient[n+2]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]-=(posA.x-comA.x)*f1*dr.x;
              Gradient[n+1]-=f1*(posA.y-comA.y)*dr.x;
              Gradient[n+2]-=(posA.y-comA.y)*f1*dr.y;
              Gradient[n+3]-=(posA.z-comA.z)*f1*dr.z;
              break;
            default:
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
  REAL temp1,temp2,temp3;

  n=NumberOfCoordinatesMinimizationVariables;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      switch(Dimension)
      {
        case 3:
          Gradient[n]+=f1*(posB.z-comB.z)*dr.z;
        case 2:
          Gradient[n]+=f1*(posB.y-comB.y)*dr.y;
        case 1:
          Gradient[n]+=f1*(posB.x-comB.x)*dr.x;
          break;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=(posB.z-comB.z)*f1*dr.z;
            case 2:
              Gradient[n+1]+=(posB.y-comB.y)*f1*dr.y;
            case 1:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              break;
          }
          break;
        case REGULAR:
          temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
          temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
          temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=temp2;
              Gradient[n+4]+=temp3;
              Gradient[n+5]+=(posB.z-comB.z)*f1*dr.z;
            case 2:
              Gradient[n+1]+=temp1;
              Gradient[n+3]+=(posB.y-comB.y)*f1*dr.y;
            case 1:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              break;
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          temp1=(posB.y-comB.y)*f1*dr.x;
          temp2=(posB.z-comB.z)*f1*dr.x;
          temp3=(posB.z-comB.z)*f1*dr.y;
          switch(Dimension)
          {
            case 3:
              Gradient[n+2]+=temp2;
              Gradient[n+4]+=temp3;
              Gradient[n+5]+=(posB.z-comB.z)*f1*dr.z;
            case 2:
              Gradient[n+1]+=temp1;
              Gradient[n+Dimension]+=(posB.y-comB.y)*f1*dr.y;
            case 1:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              break;
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+2]+=temp3;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=temp1;
              Gradient[n+2]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=temp2;
              Gradient[n+2]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            default:
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+2]+=(posB.z-comB.z)*f1*dr.y;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=(posB.z-comB.z)*f1*dr.x;
              Gradient[n+2]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n]+=(posB.x-comB.x)*f1*dr.x;
              Gradient[n+1]+=(posB.y-comB.y)*f1*dr.x;
              Gradient[n+2]+=(posB.y-comB.y)*f1*dr.y;
              Gradient[n+3]+=(posB.z-comB.z)*f1*dr.z;
              break;
            default:
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

void ComputeInterVDWMolecularPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,
                                     REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int I,J,i,j,ig,jg,ia,ja;
  int typeA,typeB;
  int TypeMolA,TypeMolB;
  REAL rr;
  REAL energy,f1,f2;
  VECTOR posA,posB,dr;
  int index1,index2;
  int index_i,index_j;
  int index_i2,index_j2;
  VECTOR comA,comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3,start;
  REAL ReplicaFactor;
  int RigidI,RigidJ;
  REAL dot_product;
  COMPLEX phase_factor;
  VECTOR dr12,pos1,pos2;

  f1=f2=0.0;
  index1=index2=0;
  // first loop over adsorbate molecules
  for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
  {
    TypeMolA=Adsorbates[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      RigidI=Components[TypeMolA].Groups[ig].Rigid;
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(RigidI)
        {
          index_i=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndex;
          index_i2=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndexOrientation;
          comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

          index1=Adsorbates[CurrentSystem][I].Atoms[i].HessianAtomIndex;
        }
        else
        {
          index_i=Adsorbates[CurrentSystem][I].Atoms[i].HessianIndex;
          index_i2=-1;
          comA=Adsorbates[CurrentSystem][I].Atoms[i].Position;
        }

        typeA=Adsorbates[CurrentSystem][I].Atoms[i].Type;
        posA=Adsorbates[CurrentSystem][I].Atoms[i].Position;

        // second loop over adsorbates
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if(ncell==0) start=I+1;
              else start=0;

              for(J=start;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
              {
                TypeMolB=Adsorbates[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];
 
                    if(RigidJ)
                    {
                      index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                      index2=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=-1;
                      comB=Adsorbates[CurrentSystem][J].Atoms[j].Position;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

                    typeB=Adsorbates[CurrentSystem][J].Atoms[j].Type;
                    posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;

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

                      PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2);

                      *Energy+=ReplicaFactor*energy;

                      if((index_i<0)&&(index_i2<0)&&(index_j<0)&&(index_j2<0)) continue;

                      StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
                      StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
                      StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

                      StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
                      StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
                      StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

                      StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
                      StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
                      StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=ReplicaFactor*(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=ReplicaFactor*temp1;
                        StrainDerivative->az-=ReplicaFactor*temp2;
                        StrainDerivative->bx-=ReplicaFactor*temp1;
                        StrainDerivative->by-=ReplicaFactor*(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=ReplicaFactor*temp3;
                        StrainDerivative->cx-=ReplicaFactor*temp2;
                        StrainDerivative->cy-=ReplicaFactor*temp3;
                        StrainDerivative->cz-=ReplicaFactor*(posA.z-comA.z)*f1*dr.z;
                      }

                      if(RigidJ)
                      {
                        temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                        temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                        temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                        StrainDerivative->ax+=ReplicaFactor*(posB.x-comB.x)*f1*dr.x;
                        StrainDerivative->ay+=ReplicaFactor*temp1;
                        StrainDerivative->az+=ReplicaFactor*temp2;
                        StrainDerivative->bx+=ReplicaFactor*temp1;
                        StrainDerivative->by+=ReplicaFactor*(posB.y-comB.y)*f1*dr.y;
                        StrainDerivative->bz+=ReplicaFactor*temp3;
                        StrainDerivative->cx+=ReplicaFactor*temp2;
                        StrainDerivative->cy+=ReplicaFactor*temp3;
                        StrainDerivative->cz+=ReplicaFactor*(posB.z-comB.z)*f1*dr.z;
                      }

                      // add contribution to the first derivatives
                      if(ComputeGradient)
                      { 
                        if(index_i>=0)
                        {
                          Gradient[index_i]+=f1*dr.x;
                          Gradient[index_i+1]+=f1*dr.y;
                          Gradient[index_i+2]+=f1*dr.z;
                        }

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2>=0)
                          {
                            Gradient[index_i2]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                            Gradient[index_i2+1]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                            Gradient[index_i2+2]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                          }
                        }

                        GradientStrain(Gradient,ReplicaFactor*f1,dr);
  
                        if(ncell==0)
                        {
                          if(index_j>=0)
                          {
                            Gradient[index_j]-=f1*dr.x;
                            Gradient[index_j+1]-=f1*dr.y;
                            Gradient[index_j+2]-=f1*dr.z;
                          }

                          if(RigidJ)
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
                      }

                      if(ComputeHessian)
                      {
                        if(RigidI) pos1=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;
                        else pos1=Adsorbates[CurrentSystem][I].Atoms[i].Position;
                        if(RigidJ) pos2=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                        else pos2=Adsorbates[CurrentSystem][J].Atoms[j].Position;

                        pos2.x+=ReplicaShift[ncell].x;
                        pos2.y+=ReplicaShift[ncell].y;
                        pos2.z+=ReplicaShift[ncell].z;

                        // form e(-k.r)
                        dr12.x=pos1.x-pos2.x;
                        dr12.y=pos1.y-pos2.y;
                        dr12.z=pos1.z-pos2.z;
                        dr12=ApplyReplicaBoundaryCondition(dr12);

                        dot_product=k.x*dr12.x+k.y*dr12.y+k.z*dr12.z;
                        phase_factor.re=cos(dot_product);
                        phase_factor.im=-sin(dot_product);

                        HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                           f1,f2,dr,ReplicaFactor);
                        HessianOrientationOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                           f1,f2,dr,ReplicaFactor);
                      }
                    }
                  }
                }
              }
              for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
              {
                TypeMolB=Cations[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];

                    if(RigidJ)
                    {
                      index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                      index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=-1;
                      comB=Cations[CurrentSystem][J].Atoms[j].Position;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

                    typeB=Cations[CurrentSystem][J].Atoms[j].Type;
                    posB=Cations[CurrentSystem][J].Atoms[j].Position;

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

                      PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2);

                      *Energy+=ReplicaFactor*energy;

                      if((index_i<0)&&(index_i2<0)&&(index_j<0)&&(index_j2<0)) continue;

                      StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
                      StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
                      StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

                      StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
                      StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
                      StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

                      StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
                      StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
                      StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=ReplicaFactor*(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=ReplicaFactor*temp1;
                        StrainDerivative->az-=ReplicaFactor*temp2;
                        StrainDerivative->bx-=ReplicaFactor*temp1;
                        StrainDerivative->by-=ReplicaFactor*(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=ReplicaFactor*temp3;
                        StrainDerivative->cx-=ReplicaFactor*temp2;
                        StrainDerivative->cy-=ReplicaFactor*temp3;
                        StrainDerivative->cz-=ReplicaFactor*(posA.z-comA.z)*f1*dr.z;
                      }


                      if(RigidJ)
                      {
                        temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                        temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                        temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                        StrainDerivative->ax+=ReplicaFactor*(posB.x-comB.x)*f1*dr.x;
                        StrainDerivative->ay+=ReplicaFactor*temp1;
                        StrainDerivative->az+=ReplicaFactor*temp2;
                        StrainDerivative->bx+=ReplicaFactor*temp1;
                        StrainDerivative->by+=ReplicaFactor*(posB.y-comB.y)*f1*dr.y;
                        StrainDerivative->bz+=ReplicaFactor*temp3;
                        StrainDerivative->cx+=ReplicaFactor*temp2;
                        StrainDerivative->cy+=ReplicaFactor*temp3;
                        StrainDerivative->cz+=ReplicaFactor*(posB.z-comB.z)*f1*dr.z;
                      }

                      // add contribution to the first derivatives
                      if(ComputeGradient)
                      {
                        if(index_i>=0)
                        {
                          Gradient[index_i]+=f1*dr.x;
                          Gradient[index_i+1]+=f1*dr.y;
                          Gradient[index_i+2]+=f1*dr.z;
                        }

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2>=0)
                          {
                            Gradient[index_i2]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                            Gradient[index_i2+1]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                            Gradient[index_i2+2]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                          }
                        }

                        GradientStrain(Gradient,ReplicaFactor*f1,dr);

                        if(ncell==0)
                        {
                          if(index_j>=0)
                          {
                            Gradient[index_j]-=f1*dr.x;
                            Gradient[index_j+1]-=f1*dr.y;
                            Gradient[index_j+2]-=f1*dr.z;
                          }

                          if(RigidJ)
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
                      }

                      if(ComputeHessian)
                      {
                        if(RigidI) pos1=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;
                        else pos1=Adsorbates[CurrentSystem][I].Atoms[i].Position;
                        if(RigidJ) pos2=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                        else pos2=Cations[CurrentSystem][J].Atoms[j].Position;

                        pos2.x+=ReplicaShift[ncell].x;
                        pos2.y+=ReplicaShift[ncell].y;
                        pos2.z+=ReplicaShift[ncell].z;

                        // form e(-k.r)
                        dr12.x=pos1.x-pos2.x;
                        dr12.y=pos1.y-pos2.y;
                        dr12.z=pos1.z-pos2.z;
                        dr12=ApplyReplicaBoundaryCondition(dr12);

                        dot_product=k.x*dr12.x+k.y*dr12.y+k.z*dr12.z;
                        phase_factor.re=cos(dot_product);
                        phase_factor.im=-sin(dot_product);

                        HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);
                        HessianOrientationOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);
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


  // Cation-Cation
  for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
  {
    TypeMolA=Cations[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      RigidI=Components[TypeMolA].Groups[ig].Rigid;
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(RigidI)
        {
          index_i=Cations[CurrentSystem][I].Groups[ig].HessianIndex;
          index_i2=Cations[CurrentSystem][I].Groups[ig].HessianIndexOrientation;
          comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

          index1=Cations[CurrentSystem][I].Atoms[i].HessianAtomIndex;
        }
        else
        {
          index_i=Cations[CurrentSystem][I].Atoms[i].HessianIndex;
          index_i2=-1;
          comA=Cations[CurrentSystem][I].Atoms[i].Position;
        }

        typeA=Cations[CurrentSystem][I].Atoms[i].Type;
        posA=Cations[CurrentSystem][I].Atoms[i].Position;

        // second loop over adsorbates
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if(ncell==0) start=I+1;
              else start=0;

              for(J=start;J<NumberOfCationMolecules[CurrentSystem];J++)
              {
                TypeMolB=Cations[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];
 
                    if(RigidJ)
                    {
                      index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                      index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=-1;
                      comB=Cations[CurrentSystem][J].Atoms[j].Position;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

                    typeB=Cations[CurrentSystem][J].Atoms[j].Type;
                    posB=Cations[CurrentSystem][J].Atoms[j].Position;

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

                      PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2);

                      *Energy+=ReplicaFactor*energy;

                      if((index_i<0)&&(index_i2<0)&&(index_j<0)&&(index_j2<0)) continue;

                      StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
                      StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
                      StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

                      StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
                      StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
                      StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

                      StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
                      StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
                      StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=ReplicaFactor*(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=ReplicaFactor*temp1;
                        StrainDerivative->az-=ReplicaFactor*temp2;
                        StrainDerivative->bx-=ReplicaFactor*temp1;
                        StrainDerivative->by-=ReplicaFactor*(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=ReplicaFactor*temp3;
                        StrainDerivative->cx-=ReplicaFactor*temp2;
                        StrainDerivative->cy-=ReplicaFactor*temp3;
                        StrainDerivative->cz-=ReplicaFactor*(posA.z-comA.z)*f1*dr.z;
                      }


                      if(RigidJ)
                      {
                        temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                        temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                        temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                        StrainDerivative->ax+=ReplicaFactor*(posB.x-comB.x)*f1*dr.x;
                        StrainDerivative->ay+=ReplicaFactor*temp1;
                        StrainDerivative->az+=ReplicaFactor*temp2;
                        StrainDerivative->bx+=ReplicaFactor*temp1;
                        StrainDerivative->by+=ReplicaFactor*(posB.y-comB.y)*f1*dr.y;
                        StrainDerivative->bz+=ReplicaFactor*temp3;
                        StrainDerivative->cx+=ReplicaFactor*temp2;
                        StrainDerivative->cy+=ReplicaFactor*temp3;
                        StrainDerivative->cz+=ReplicaFactor*(posB.z-comB.z)*f1*dr.z;
                      }

                      // add contribution to the first derivatives
                      if(ComputeGradient)
                      { 
                        if(index_i>=0)
                        {
                          Gradient[index_i]+=f1*dr.x;
                          Gradient[index_i+1]+=f1*dr.y;
                          Gradient[index_i+2]+=f1*dr.z;
                        }

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2>=0)
                          {
                            Gradient[index_i2]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                            Gradient[index_i2+1]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                            Gradient[index_i2+2]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                          }
                        }

                        GradientStrain(Gradient,ReplicaFactor*f1,dr);
  
                        if(ncell==0)
                        {
                          if(index_j>=0)
                          {
                            Gradient[index_j]-=f1*dr.x;
                            Gradient[index_j+1]-=f1*dr.y;
                            Gradient[index_j+2]-=f1*dr.z;
                          }

                          if(RigidJ)
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
                      }

                      if(ComputeHessian)
                      {
                        if(RigidI) pos1=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;
                        else pos1=Cations[CurrentSystem][I].Atoms[i].Position;
                        if(RigidJ) pos2=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                        else pos2=Cations[CurrentSystem][J].Atoms[j].Position;

                        pos2.x+=ReplicaShift[ncell].x;
                        pos2.y+=ReplicaShift[ncell].y;
                        pos2.z+=ReplicaShift[ncell].z;

                        // form e(-k.r)
                        dr12.x=pos1.x-pos2.x;
                        dr12.y=pos1.y-pos2.y;
                        dr12.z=pos1.z-pos2.z;
                        dr12=ApplyReplicaBoundaryCondition(dr12);

                        dot_product=k.x*dr12.x+k.y*dr12.y+k.z*dr12.z;
                        phase_factor.re=cos(dot_product);
                        phase_factor.im=-sin(dot_product);

                        HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);
                        HessianOrientationOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);
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
}

void ComputeInterChargeChargeMolecularPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,
                                              REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int I,J,i,j,ig,jg,ia,ja;
  int typeA,typeB;
  int TypeMolA,TypeMolB;
  REAL ChargeA,ChargeB,rr;
  REAL U,f1,f2;
  VECTOR posA,posB,dr;
  int index_i,index_j;
  int index_i2,index_j2;
  VECTOR comA,comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3,start;
  REAL ReplicaFactor;
  int index1,index2;
  int RigidI,RigidJ;
  REAL SwitchingValue,SwitchingValueDerivative,SwitchingValueSecondDerivative;
  REAL TranslationValue,TranslationValueDerivative,TranslationValueSecondDerivative;
  COMPLEX phase_factor;
  REAL dot_product;
  VECTOR dr12,pos1,pos2;

  f1=f2=0.0;
  index1=index2=0;
  // first loop over adsorbate molecules
  for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
  {
    TypeMolA=Adsorbates[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      RigidI=Components[TypeMolA].Groups[ig].Rigid;
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(RigidI)
        {
          index_i=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndex;
          index_i2=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndexOrientation;
          comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

          index1=Adsorbates[CurrentSystem][I].Atoms[i].HessianAtomIndex;
        }
        else
        {
          index_i=Adsorbates[CurrentSystem][I].Atoms[i].HessianIndex;
          index_i2=-1;
          comA=Adsorbates[CurrentSystem][I].Atoms[i].Position;
        }

        typeA=Adsorbates[CurrentSystem][I].Atoms[i].Type;
        posA=Adsorbates[CurrentSystem][I].Atoms[i].Position;
        ChargeA=Adsorbates[CurrentSystem][I].Atoms[i].Charge;

        // second loop over adsorbates
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if(ncell==0) start=I+1;
              else start=0;

              for(J=start;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
              {
                TypeMolB=Adsorbates[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];
 
                    if(RigidJ)
                    {
                      index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                      index2=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=-1;
                      comB=Adsorbates[CurrentSystem][J].Atoms[j].Position;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

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
                      if(ncell==0) ReplicaFactor=1.0;
                      else ReplicaFactor=0.5;

                      switch(ChargeMethod)
                      {
                        case NONE:
                          U=f1=f2=0.0;
                          break;
                        case TRUNCATED_COULOMB:
                          *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r);
                          f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                          f2=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                          break;
                        case SHIFTED_COULOMB:
                          *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r-InverseCutOffChargeCharge);
                          f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                          f2=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                          break;
                        case SMOOTHED_COULOMB:
                          U=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                          f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/rr;
                          f2=2.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                          if(rr>CutOffChargeChargeSwitchSquared)
                          {
                            SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                           SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
                            SwitchingValueDerivative=5.0*SwitchingChargeChargeFactors5[5]*rr*rr+4.0*SwitchingChargeChargeFactors5[4]*rr*r+3.0*SwitchingChargeChargeFactors5[3]*rr+
                                                     2.0*SwitchingChargeChargeFactors5[2]*r+SwitchingChargeChargeFactors5[1];
                            SwitchingValueSecondDerivative=20.0*SwitchingChargeChargeFactors5[5]*(rr*r)+12.0*SwitchingChargeChargeFactors5[4]*(rr)+
                                                           6.0*SwitchingChargeChargeFactors5[3]*r+2.0*SwitchingChargeChargeFactors5[2];


                            TranslationValue=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                                            (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                             SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                             SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                            TranslationValueDerivative=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                                                      (7.0*SwitchingChargeChargeFactors7[7]*rr*rr*rr+6.0*SwitchingChargeChargeFactors7[6]*rr*rr*r+
                                                       5.0*SwitchingChargeChargeFactors7[5]*rr*rr+4.0*SwitchingChargeChargeFactors7[4]*rr*r+3.0*SwitchingChargeChargeFactors7[3]*rr+
                                                       2.0*SwitchingChargeChargeFactors7[2]*r+SwitchingChargeChargeFactors7[1]);
                            TranslationValueSecondDerivative=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                                                      (42.0*SwitchingChargeChargeFactors7[7]*rr*rr*r+30.0*SwitchingChargeChargeFactors7[6]*rr*rr+
                                                       20.0*SwitchingChargeChargeFactors7[5]*rr*r+12.0*SwitchingChargeChargeFactors7[4]*rr+6.0*SwitchingChargeChargeFactors7[3]*r+
                                                       2.0*SwitchingChargeChargeFactors7[2]);
                            f2=U*SwitchingValueSecondDerivative+2.0*f1*SwitchingValueDerivative+f2*SwitchingValue+TranslationValueSecondDerivative;
                            f1=U*SwitchingValueDerivative+f1*SwitchingValue+TranslationValueDerivative;
                            U=U*SwitchingValue+TranslationValue;
                          }
                          *Energy+=U;
                          f2=(f2*r-f1)/CUBE(r);
                          f1/=r;
                          break;
                        case EWALD:
                        default:
                          *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*erfc(Alpha[CurrentSystem]*r)*
                                                    ChargeA*ChargeB/r;

                          f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                              (erfc(Alpha[CurrentSystem]*r)+
                              2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                              (r*rr);

                          f2=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                              (((3.0*erfc(Alpha[CurrentSystem]*r)/(r*rr))+
                              (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                              (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));
                          break;
                      }

                      if((index_i<0)&&(index_i2<0)&&(index_j<0)&&(index_j2<0)) continue;

                      StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
                      StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
                      StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

                      StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
                      StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
                      StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

                      StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
                      StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
                      StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=ReplicaFactor*(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=ReplicaFactor*temp1;
                        StrainDerivative->az-=ReplicaFactor*temp2;
                        StrainDerivative->bx-=ReplicaFactor*temp1;
                        StrainDerivative->by-=ReplicaFactor*(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=ReplicaFactor*temp3;
                        StrainDerivative->cx-=ReplicaFactor*temp2;
                        StrainDerivative->cy-=ReplicaFactor*temp3;
                        StrainDerivative->cz-=ReplicaFactor*(posA.z-comA.z)*f1*dr.z;
                      }


                      if(RigidJ)
                      {
                        temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                        temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                        temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                        StrainDerivative->ax+=ReplicaFactor*(posB.x-comB.x)*f1*dr.x;
                        StrainDerivative->ay+=ReplicaFactor*temp1;
                        StrainDerivative->az+=ReplicaFactor*temp2;
                        StrainDerivative->bx+=ReplicaFactor*temp1;
                        StrainDerivative->by+=ReplicaFactor*(posB.y-comB.y)*f1*dr.y;
                        StrainDerivative->bz+=ReplicaFactor*temp3;
                        StrainDerivative->cx+=ReplicaFactor*temp2;
                        StrainDerivative->cy+=ReplicaFactor*temp3;
                        StrainDerivative->cz+=ReplicaFactor*(posB.z-comB.z)*f1*dr.z;
                      }

                      // add contribution to the first derivatives
                      if(ComputeGradient)
                      {
                        if(index_i>=0)
                        {
                          Gradient[index_i]+=f1*dr.x;
                          Gradient[index_i+1]+=f1*dr.y;
                          Gradient[index_i+2]+=f1*dr.z;
                        }

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2>=0)
                          {
                            Gradient[index_i2]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                            Gradient[index_i2+1]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                            Gradient[index_i2+2]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                          }
                        }

                        GradientStrain(Gradient,ReplicaFactor*f1,dr);
  
                        if(ncell==0)
                        {
                          if(index_j>=0)
                          {
                            Gradient[index_j]-=f1*dr.x;
                            Gradient[index_j+1]-=f1*dr.y;
                            Gradient[index_j+2]-=f1*dr.z;
                          }

                          if(RigidJ)
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
                      }

                      if(ComputeHessian)
                      {
                        if(RigidI) pos1=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;
                        else pos1=Adsorbates[CurrentSystem][I].Atoms[i].Position;

                        if(RigidJ) pos2=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;
                        else pos2=Adsorbates[CurrentSystem][J].Atoms[j].Position;

                        pos2.x+=ReplicaShift[ncell].x;
                        pos2.y+=ReplicaShift[ncell].y;
                        pos2.z+=ReplicaShift[ncell].z;

                        // form e(-k.r)
                        dr12.x=pos1.x-pos2.x;
                        dr12.y=pos1.y-pos2.y;
                        dr12.z=pos1.z-pos2.z;
                        dr12=ApplyReplicaBoundaryCondition(dr12);

                        dot_product=k.x*dr12.x+k.y*dr12.y+k.z*dr12.z;
                        phase_factor.re=cos(dot_product);
                        phase_factor.im=-sin(dot_product);

                        HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);
                        HessianOrientationOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);

                      }
                    }
                  }
                }
              }
              for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
              {
                TypeMolB=Cations[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];

                    if(RigidJ)
                    {
                      index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                      index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=-1;
                      comB=Cations[CurrentSystem][J].Atoms[j].Position;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

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
                      if(ncell==0) ReplicaFactor=1.0;
                      else ReplicaFactor=0.5;

                      switch(ChargeMethod)
                      {
                        case NONE:
                          f1=f2=0.0;
                          break;
                        case SHIFTED_COULOMB:
                          *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r-InverseCutOffChargeCharge);
                          f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                          f2=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                          break;
                        case TRUNCATED_COULOMB:
                          *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r);
                          f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                          f2=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                          break;
                        case EWALD:
                        default:
                          *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*erfc(Alpha[CurrentSystem]*r)*
                                                    ChargeA*ChargeB/r;

                          f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                              (erfc(Alpha[CurrentSystem]*r)+
                              2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                              (r*rr);

                          f2=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                              (((3.0*erfc(Alpha[CurrentSystem]*r)/(r*rr))+
                              (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                              (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));
                          break;
                      }

                      if((index_i<0)&&(index_i2<0)&&(index_j<0)&&(index_j2<0)) continue;

                      StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
                      StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
                      StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

                      StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
                      StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
                      StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

                      StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
                      StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
                      StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=ReplicaFactor*(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=ReplicaFactor*temp1;
                        StrainDerivative->az-=ReplicaFactor*temp2;
                        StrainDerivative->bx-=ReplicaFactor*temp1;
                        StrainDerivative->by-=ReplicaFactor*(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=ReplicaFactor*temp3;
                        StrainDerivative->cx-=ReplicaFactor*temp2;
                        StrainDerivative->cy-=ReplicaFactor*temp3;
                        StrainDerivative->cz-=ReplicaFactor*(posA.z-comA.z)*f1*dr.z;
                      }


                      if(RigidJ)
                      {
                        temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                        temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                        temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                        StrainDerivative->ax+=ReplicaFactor*(posB.x-comB.x)*f1*dr.x;
                        StrainDerivative->ay+=ReplicaFactor*temp1;
                        StrainDerivative->az+=ReplicaFactor*temp2;
                        StrainDerivative->bx+=ReplicaFactor*temp1;
                        StrainDerivative->by+=ReplicaFactor*(posB.y-comB.y)*f1*dr.y;
                        StrainDerivative->bz+=ReplicaFactor*temp3;
                        StrainDerivative->cx+=ReplicaFactor*temp2;
                        StrainDerivative->cy+=ReplicaFactor*temp3;
                        StrainDerivative->cz+=ReplicaFactor*(posB.z-comB.z)*f1*dr.z;
                      }

                      // add contribution to the first derivatives
                      if(ComputeGradient)
                      {
                        if(index_i>=0)
                        {
                          Gradient[index_i]+=f1*dr.x;
                          Gradient[index_i+1]+=f1*dr.y;
                          Gradient[index_i+2]+=f1*dr.z;
                        }

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2>=0)
                          {
                            Gradient[index_i2]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                            Gradient[index_i2+1]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                            Gradient[index_i2+2]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                          }
                        }

                        GradientStrain(Gradient,ReplicaFactor*f1,dr);

                        if(ncell==0)
                        {
                          if(index_j>=0)
                          {
                            Gradient[index_j]-=f1*dr.x;
                            Gradient[index_j+1]-=f1*dr.y;
                            Gradient[index_j+2]-=f1*dr.z;
                          }

                          if(RigidJ)
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
                      }

                      if(ComputeHessian)
                      {
                        HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);
                        HessianOrientationOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);

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


  // Cation-Cation
  for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
  {
    TypeMolA=Cations[CurrentSystem][I].Type;
    for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
    {
      RigidI=Components[TypeMolA].Groups[ig].Rigid;
      for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
      {
        i=Components[TypeMolA].Groups[ig].Atoms[ia];

        if(RigidI)
        {
          index_i=Cations[CurrentSystem][I].Groups[ig].HessianIndex;
          index_i2=Cations[CurrentSystem][I].Groups[ig].HessianIndexOrientation;
          comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

          index1=Cations[CurrentSystem][I].Atoms[i].HessianAtomIndex;
        }
        else
        {
          index_i=Cations[CurrentSystem][I].Atoms[i].HessianIndex;
          index_i2=-1;
          comA=Cations[CurrentSystem][I].Atoms[i].Position;
        }

        typeA=Cations[CurrentSystem][I].Atoms[i].Type;
        posA=Cations[CurrentSystem][I].Atoms[i].Position;
        ChargeA=Cations[CurrentSystem][I].Atoms[i].Charge;

        // second loop over adsorbates
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if(ncell==0) start=I+1;
              else start=0;

              for(J=start;J<NumberOfCationMolecules[CurrentSystem];J++)
              {
                TypeMolB=Cations[CurrentSystem][J].Type;
                for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                {
                  RigidJ=Components[TypeMolB].Groups[jg].Rigid;
                  for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                  {
                    j=Components[TypeMolB].Groups[jg].Atoms[ja];
 
                    if(RigidJ)
                    {
                      index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                      index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                      comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                      index2=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                    }
                    else
                    {
                      index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                      index_j2=-1;
                      comB=Cations[CurrentSystem][J].Atoms[j].Position;
                    }

                    comB.x+=ReplicaShift[ncell].x;
                    comB.y+=ReplicaShift[ncell].y;
                    comB.z+=ReplicaShift[ncell].z;

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
                      if(ncell==0) ReplicaFactor=1.0;
                      else ReplicaFactor=0.5;

                      switch(ChargeMethod)
                      {
                        case NONE:
                          f1=f2=0.0;
                          break;
                        case SHIFTED_COULOMB:
                          *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r-InverseCutOffChargeCharge);
                          f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                          f2=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                          break;
                        case TRUNCATED_COULOMB:
                          *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*(1.0/r);
                          f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(rr*r);
                          f2=3.0*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/(SQR(rr)*r);
                          break;
                        case EWALD:
                        default:
                          *Energy+=ReplicaFactor*COULOMBIC_CONVERSION_FACTOR*erfc(Alpha[CurrentSystem]*r)*
                                                    ChargeA*ChargeB/r;

                          f1=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                              (erfc(Alpha[CurrentSystem]*r)+
                              2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                              (r*rr);

                          f2=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                              (((3.0*erfc(Alpha[CurrentSystem]*r)/(r*rr))+
                              (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                              (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));
                          break;
                      }

                      StrainDerivative->ax+=ReplicaFactor*f1*dr.x*dr.x;
                      StrainDerivative->ay+=ReplicaFactor*f1*dr.x*dr.y;
                      StrainDerivative->az+=ReplicaFactor*f1*dr.x*dr.z;

                      StrainDerivative->bx+=ReplicaFactor*f1*dr.y*dr.x;
                      StrainDerivative->by+=ReplicaFactor*f1*dr.y*dr.y;
                      StrainDerivative->bz+=ReplicaFactor*f1*dr.y*dr.z;

                      StrainDerivative->cx+=ReplicaFactor*f1*dr.z*dr.x;
                      StrainDerivative->cy+=ReplicaFactor*f1*dr.z*dr.y;
                      StrainDerivative->cz+=ReplicaFactor*f1*dr.z*dr.z;

                      if(RigidI)
                      {
                        temp1=0.5*((posA.y-comA.y)*f1*dr.x+(posA.x-comA.x)*f1*dr.y);
                        temp2=0.5*((posA.z-comA.z)*f1*dr.x+(posA.x-comA.x)*f1*dr.z);
                        temp3=0.5*((posA.z-comA.z)*f1*dr.y+(posA.y-comA.y)*f1*dr.z);

                        StrainDerivative->ax-=ReplicaFactor*(posA.x-comA.x)*f1*dr.x;
                        StrainDerivative->ay-=ReplicaFactor*temp1;
                        StrainDerivative->az-=ReplicaFactor*temp2;
                        StrainDerivative->bx-=ReplicaFactor*temp1;
                        StrainDerivative->by-=ReplicaFactor*(posA.y-comA.y)*f1*dr.y;
                        StrainDerivative->bz-=ReplicaFactor*temp3;
                        StrainDerivative->cx-=ReplicaFactor*temp2;
                        StrainDerivative->cy-=ReplicaFactor*temp3;
                        StrainDerivative->cz-=ReplicaFactor*(posA.z-comA.z)*f1*dr.z;
                      }


                      if(RigidJ)
                      {
                        temp1=0.5*((posB.y-comB.y)*f1*dr.x+(posB.x-comB.x)*f1*dr.y);
                        temp2=0.5*((posB.z-comB.z)*f1*dr.x+(posB.x-comB.x)*f1*dr.z);
                        temp3=0.5*((posB.z-comB.z)*f1*dr.y+(posB.y-comB.y)*f1*dr.z);
                        StrainDerivative->ax+=ReplicaFactor*(posB.x-comB.x)*f1*dr.x;
                        StrainDerivative->ay+=ReplicaFactor*temp1;
                        StrainDerivative->az+=ReplicaFactor*temp2;
                        StrainDerivative->bx+=ReplicaFactor*temp1;
                        StrainDerivative->by+=ReplicaFactor*(posB.y-comB.y)*f1*dr.y;
                        StrainDerivative->bz+=ReplicaFactor*temp3;
                        StrainDerivative->cx+=ReplicaFactor*temp2;
                        StrainDerivative->cy+=ReplicaFactor*temp3;
                        StrainDerivative->cz+=ReplicaFactor*(posB.z-comB.z)*f1*dr.z;
                      }

                      // add contribution to the first derivatives
                      if(ComputeGradient)
                      {
                        if(index_i>=0)
                        {
                          Gradient[index_i]+=f1*dr.x;
                          Gradient[index_i+1]+=f1*dr.y;
                          Gradient[index_i+2]+=f1*dr.z;
                        }

                        if(RigidI)
                        {
                          GradientStrainI(Gradient,f1,dr,posA,comA);

                          // add contribution to the first derivatives
                          if(index_i2>=0)
                          {
                            Gradient[index_i2]+=f1*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                            Gradient[index_i2+1]+=f1*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                            Gradient[index_i2+2]+=f1*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
                          }
                        }

                        GradientStrain(Gradient,ReplicaFactor*f1,dr);
  
                        if(ncell==0)
                        {
                          if(index_j>=0)
                          {
                            Gradient[index_j]-=f1*dr.x;
                            Gradient[index_j+1]-=f1*dr.y;
                            Gradient[index_j+2]-=f1*dr.z;
                          }

                          if(RigidJ)
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
                      }

                      if(ComputeHessian)
                      {
                        HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,f1,f2,dr,ReplicaFactor);

                        HessianCenterOfMassOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);
                        HessianOrientationOrientation(phase_factor,HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,
                                                       f1,f2,dr,ReplicaFactor);
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
}
