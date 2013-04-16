/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_phonon.c' is part of RASPA-2.0

 *************************************************************************************************************/

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


// Hessian: Center of mass - Center of mass
// ========================================
static inline void HessianAtomicPositionPosition(COMPLEX phase_factor,COMPLEX_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor)
{
  REAL_MATRIX3x3 Hessian;

  Hessian.ax=f2*dr.x*dr.x+f1; Hessian.bx=f2*dr.y*dr.x;    Hessian.cx=f2*dr.z*dr.x;
  Hessian.ay=f2*dr.x*dr.y;    Hessian.by=f2*dr.y*dr.y+f1; Hessian.cy=f2*dr.z*dr.y;
  Hessian.az=f2*dr.x*dr.z;    Hessian.bz=f2*dr.y*dr.z;    Hessian.cz=f2*dr.z*dr.z+f1;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  { 
    if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x].re+=ReplicaFactor*Hessian.ax;
    if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y].re+=ReplicaFactor*Hessian.ay;
    if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z].re+=ReplicaFactor*Hessian.az;
    if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y].re+=ReplicaFactor*Hessian.by;
    if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z].re+=ReplicaFactor*Hessian.bz;
    if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z].re+=ReplicaFactor*Hessian.cz;
  }

  if((index_j.x>=0)&&(index_j.x>=0))  HessianMatrix.element[index_j.x][index_j.x].re+=ReplicaFactor*Hessian.ax;
  if((index_j.x>=0)&&(index_j.y>=0))  HessianMatrix.element[index_j.x][index_j.y].re+=ReplicaFactor*Hessian.ay;
  if((index_j.x>=0)&&(index_j.z>=0))  HessianMatrix.element[index_j.x][index_j.z].re+=ReplicaFactor*Hessian.az;
  if((index_j.y>=0)&&(index_j.y>=0))  HessianMatrix.element[index_j.y][index_j.y].re+=ReplicaFactor*Hessian.by;
  if((index_j.y>=0)&&(index_j.z>=0))  HessianMatrix.element[index_j.y][index_j.z].re+=ReplicaFactor*Hessian.bz;
  if((index_j.z>=0)&&(index_j.z>=0))  HessianMatrix.element[index_j.z][index_j.z].re+=ReplicaFactor*Hessian.cz;

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x].re-=phase_factor.re*Hessian.ax;
    if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x].im-=phase_factor.im*Hessian.ax;
    if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y].re-=phase_factor.re*Hessian.ay;
    if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y].im-=phase_factor.im*Hessian.ay;
    if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z].re-=phase_factor.re*Hessian.az;
    if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z].im-=phase_factor.im*Hessian.az;
    if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x].re-=phase_factor.re*Hessian.ay;
    if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x].im-=phase_factor.im*Hessian.ay;
    if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y].re-=phase_factor.re*Hessian.by;
    if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y].im-=phase_factor.im*Hessian.by;
    if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z].re-=phase_factor.re*Hessian.bz;
    if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z].im-=phase_factor.im*Hessian.bz;
    if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x].re-=phase_factor.re*Hessian.az;
    if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x].im-=phase_factor.im*Hessian.az;
    if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y].re-=phase_factor.re*Hessian.bz;
    if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y].im-=phase_factor.im*Hessian.bz;
    if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z].re-=phase_factor.re*Hessian.cz;
    if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z].im-=phase_factor.im*Hessian.cz;
  }
}


// Hessian: Center of mass - Orientation
// =====================================
static inline void HessianCenterOfMassOrientationJ(COMPLEX_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,INT_VECTOR3 index_j2,int index2,REAL f1,REAL f2,VECTOR dr)
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

  if((index_j.x>=0)&&(index_j2.x>=0))  HessianMatrix.element[index_j.x][index_j2.x].re+=Hessian.ax;
  if((index_j.x>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_j.x][index_j2.y].re+=Hessian.ay;
  if((index_j.x>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_j.x][index_j2.z].re+=Hessian.az;

  if((index_j.y>=0)&&(index_j2.x>=0))  HessianMatrix.element[index_j.y][index_j2.x].re+=Hessian.bx;
  if((index_j.y>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_j.y][index_j2.y].re+=Hessian.by;
  if((index_j.y>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_j.y][index_j2.z].re+=Hessian.bz;

  if((index_j.z>=0)&&(index_j2.x>=0))  HessianMatrix.element[index_j.z][index_j2.x].re+=Hessian.cx;
  if((index_j.z>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_j.z][index_j2.y].re+=Hessian.cy;
  if((index_j.z>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_j.z][index_j2.z].re+=Hessian.cz;

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
    if((index_i.x>=0)&&(index_j2.x>=0))  HessianMatrix.element[index_i.x][index_j2.x].re-=Hessian.ax;
    if((index_i.x>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_i.x][index_j2.y].re-=Hessian.ay;
    if((index_i.x>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_i.x][index_j2.z].re-=Hessian.az;

    if((index_i.y>=0)&&(index_j2.x>=0))  HessianMatrix.element[index_i.y][index_j2.x].re-=Hessian.bx;
    if((index_i.y>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_i.y][index_j2.y].re-=Hessian.by;
    if((index_i.y>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_i.y][index_j2.z].re-=Hessian.bz;

    if((index_i.z>=0)&&(index_j2.x>=0))  HessianMatrix.element[index_i.z][index_j2.x].re-=Hessian.cx;
    if((index_i.z>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_i.z][index_j2.y].re-=Hessian.cy;
    if((index_i.z>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_i.z][index_j2.z].re-=Hessian.cz;
  }
}

// Hessian: Orientation - Orientation
// ==================================
static inline void HessianOrientationOrientationJ(COMPLEX_MATRIX HessianMatrix,INT_VECTOR3 index_j2,int index2,REAL f1,REAL f2,VECTOR dr)
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

  if((index_j2.x>=0)&&(index_j2.x>=0))  HessianMatrix.element[index_j2.x][index_j2.x].re+=Hessian.ax;
  if((index_j2.y>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_j2.y][index_j2.y].re+=Hessian.by;
  if((index_j2.z>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_j2.z][index_j2.z].re+=Hessian.cz;
  if((index_j2.x>=0)&&(index_j2.y>=0))  HessianMatrix.element[index_j2.x][index_j2.y].re+=Hessian.ay;
  if((index_j2.x>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_j2.x][index_j2.z].re+=Hessian.az;
  if((index_j2.y>=0)&&(index_j2.z>=0))  HessianMatrix.element[index_j2.y][index_j2.z].re+=Hessian.bz;
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

void ComputeFrameworkAdsorbateVDWPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,
                                         REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr;
  REAL energy,f1,f2;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR pos,comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;
  COMPLEX phase_factor;
  REAL scalingB;

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
                    index_j2=UNDEFINED_INT_VECTOR3;
                    index2=-1;
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
                    scalingB=Adsorbates[CurrentSystem][J].Atoms[j].CFVDWScalingParameter;

                    PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2,scalingB);

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
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;
                      }
  
                      if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                      if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                      if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                        // add contribution to the first derivatives
                        if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                        if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                        if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,f1,f2,dr,1.0);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
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


void ComputeFrameworkAdsorbateChargeChargePhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,
                                                  REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr,U;
  REAL f1,f2;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;
  COMPLEX phase_factor;
  REAL scalingB;

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
                    index_j2=UNDEFINED_INT_VECTOR3;
                    index2=-1;
                  }

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

                  if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  {
                    scalingB=Adsorbates[CurrentSystem][J].Atoms[j].CFChargeScalingParameter;
                    ChargeB=scalingB*Adsorbates[CurrentSystem][J].Atoms[j].Charge;

                    PotentialSecondDerivativeCoulombic(ChargeA,ChargeB,rr,&U,&f1,&f2);

                    *Energy+=U;

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
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;
                      }

                      if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                      if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                      if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;
  
                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                         // add contribution to the first derivatives
                        if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                        if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                        if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,f1,f2,dr,1.0);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
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

void ComputeFrameworkCationVDWPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,
                                      REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr;
  REAL energy,f1,f2;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;
  COMPLEX phase_factor;
  REAL scalingB;

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
                    index_j2=UNDEFINED_INT_VECTOR3;
                    index2=-1;
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
                    scalingB=Cations[CurrentSystem][J].Atoms[j].CFVDWScalingParameter;

                    PotentialSecondDerivative(typeA,typeB,rr,&energy,&f1,&f2,scalingB);

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
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;
                      }

                      if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                      if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                      if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                        // add contribution to the first derivatives
                        if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                        if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                        if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,f1,f2,dr,1.0);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
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

void ComputeFrameworkCationChargeChargePhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,
                                               REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int J,i,j,jg,ja,fr;
  int typeA,typeB;
  int TypeMolB;
  REAL rr,U;
  REAL f1,f2;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j,index_j2;
  int index2;
  REAL ChargeA,ChargeB;
  VECTOR comB;
  REAL r,temp1,temp2,temp3;
  int ncell,k1,k2,k3;
  COMPLEX phase_factor;
  REAL scalingB;

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
                    index_j2=UNDEFINED_INT_VECTOR3;
                    index2=-1;
                  }

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

                  if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  {
                    scalingB=Adsorbates[CurrentSystem][J].Atoms[j].CFChargeScalingParameter;
                    ChargeB=scalingB*Cations[CurrentSystem][J].Atoms[j].Charge;

                    PotentialSecondDerivativeCoulombic(ChargeA,ChargeB,rr,&U,&f1,&f2);

                    *Energy+=U;

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
                        if(index_i.x>=0) Gradient[index_i.x]+=f1*dr.x;
                        if(index_i.y>=0) Gradient[index_i.y]+=f1*dr.y;
                        if(index_i.z>=0) Gradient[index_i.z]+=f1*dr.z;
                      }

                      if(index_j.x>=0) Gradient[index_j.x]-=f1*dr.x;
                      if(index_j.y>=0) Gradient[index_j.y]-=f1*dr.y;
                      if(index_j.z>=0) Gradient[index_j.z]-=f1*dr.z;

                      GradientStrain(Gradient,f1,dr);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        GradientStrainJ(Gradient,f1,dr,posB,comB);

                        // add contribution to the first derivatives
                        if(index_j2.x>=0) Gradient[index_j2.x]-=f1*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                        if(index_j2.y>=0) Gradient[index_j2.y]-=f1*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                        if(index_j2.z>=0) Gradient[index_j2.z]-=f1*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
                      }
                    }

                    if(ComputeHessian)
                    {
                      HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,f1,f2,dr,1.0);

                      if(Components[TypeMolB].Groups[jg].Rigid)
                      {
                        HessianCenterOfMassOrientationJ(HessianMatrix,index_i,index_j,index_j2,index2,f1,f2,dr);
                        HessianOrientationOrientationJ(HessianMatrix,index_j2,index2,f1,f2,dr);
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


void ComputeFrameworkIntraVDWPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,
                                     REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,j,typeA,typeB,f1,f2,start;
  REAL energy,DF,DDF;
  REAL rr,r;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j;
  int ncell,k1,k2,k3,indj;
  REAL ReplicaFactor;
  REAL dot_product;
  COMPLEX phase_factor;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
      posA=Framework[CurrentSystem].Atoms[f1][i].Position;
      index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            if(ncell==0) start=i+1;
            else start=0;
            for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
            {
              indj=j+ncell*Framework[CurrentSystem].NumberOfAtoms[f1];
              if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][indj],0))
              {
                typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
                posB=Framework[CurrentSystem].Atoms[f1][j].Position;
                index_j=Framework[CurrentSystem].Atoms[f1][j].HessianIndex;

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

                  PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,1.0);

                  // add contribution to the energy
                  *Energy+=ReplicaFactor*energy;

                  if(ComputeGradient)
                  {
                    // add contribution to the first derivatives
                    if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
                    if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
                    if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

                    if(ncell==0)
                    {
                      if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
                      if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
                      if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;
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
                    // form e(-k.r)
                    dot_product=k.x*dr.x+k.y*dr.y+k.z*dr.z;
                    phase_factor.re=cos(dot_product);
                    phase_factor.im=-sin(dot_product);

                    // add contribution to the second derivatives (Hessian matrix)
                    HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,DF,DDF,dr,ReplicaFactor);
                  }
                }
              }
            }
            ncell++;
          }
    }
  }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
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

                    PotentialSecondDerivative(typeA,typeB,rr,&energy,&DF,&DDF,1.0);

                    // add contribution to the energy
                    *Energy+=ReplicaFactor*energy;

                    if(ComputeGradient)
                    {
                      // add contribution to the first derivatives
                      if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
                      if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
                      if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

                      if(ncell==0)
                      {
                        if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
                        if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
                        if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;
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
                      dr.x=posA.x-posB.x;
                      dr.y=posA.y-posB.y;
                      dr.z=posA.z-posB.z;
                      dr=ApplyReplicaBoundaryCondition(dr);

                      // form e(-k.r)
                      dot_product=k.x*dr.x+k.y*dr.y+k.z*dr.z;
                      phase_factor.re=cos(dot_product);
                      phase_factor.im=-sin(dot_product);

                      // add contribution to the second derivatives (Hessian matrix)
                      HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,DF,DDF,dr,ReplicaFactor);
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

void ComputeFrameworkIntraChargeChargePhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,
                                              REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,j,typeA,typeB,f1,f2,start;
  REAL ChargeA,ChargeB,DF,DDF;
  REAL rr,r,U;
  VECTOR posA,posB,dr;
  INT_VECTOR3 index_i,index_j;
  int ncell,k1,k2,k3,indj;
  REAL ReplicaFactor;
  COMPLEX phase_factor;
  REAL dot_product;

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
                if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][indj],1))
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

                  if(rr<CutOffChargeChargeSquared[CurrentSystem])
                  {
                    if(ncell==0) ReplicaFactor=1.0;
                    else ReplicaFactor=0.5;

                    PotentialSecondDerivativeCoulombic(ChargeA,ChargeB,rr,&U,&DF,&DDF);

                    *Energy+=ReplicaFactor*U;

                    if(ComputeGradient)
                    {
                      // add contribution to the first derivatives
                      if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
                      if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
                      if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

                      if(ncell==0)
                      {
                        if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
                        if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
                        if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;
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
                      dr.x=posA.x-posB.x;
                      dr.y=posA.y-posB.y;
                      dr.z=posA.z-posB.z;
                      dr=ApplyReplicaBoundaryCondition(dr);

                      // form e(-k.r)
                      dot_product=k.x*dr.x+k.y*dr.y+k.z*dr.z;
                      phase_factor.re=cos(dot_product);
                      phase_factor.im=-sin(dot_product);

                      // add contribution to the second derivatives (Hessian matrix)
                      HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,DF,DDF,dr,ReplicaFactor);
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



void ComputeFrameworkBondPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
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
  INT_VECTOR3 index_i,index_j; // indices of the Hessian
  REAL *parms;  // pointer to potential parameter
  COMPLEX phase_factor;
  REAL dot_product;

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

      if(ComputeGradient)
      {
        // add contribution to the first derivatives
        if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
        if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
        if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

        if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
        if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
        if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

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
        // form e(-k.r)
        dot_product=k.x*dr.x+k.y*dr.y+k.z*dr.z;
        phase_factor.re=cos(dot_product);
        phase_factor.im=-sin(dot_product);

        // add contribution to the second derivatives (Hessian matrix)
        HessianAtomicPositionPosition(phase_factor,HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
      }
    }
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



void ComputeFrameworkBendPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,f1,A,B,C;
  REAL *parms,U;
  REAL CosTheta,Theta,SinTheta,temp,temp2;
  REAL rab,rbc,rac;
  POINT posA,posB,posC;
  VECTOR Rab,Rbc,Rac;
  INT_VECTOR3 index_i,index_j,index_k;
  int index1,index2,index3;
  REAL DTDX,DF,DDF;
  REAL_MATRIX3x3 D2I,D2K,D2IK;
  VECTOR dtA,dtB,dtC;
  REAL_MATRIX3x3 S;
  VECTOR vec_u,vec_v;
  REAL u,v,dot_product;
  COMPLEX phase_factor_ab,phase_factor_bc,phase_factor_ac;

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

      index1=Framework[CurrentSystem].Atoms[f1][A].HessianAtomIndex;
      index2=Framework[CurrentSystem].Atoms[f1][B].HessianAtomIndex;
      index3=Framework[CurrentSystem].Atoms[f1][C].HessianAtomIndex;

      posA=Framework[CurrentSystem].Atoms[f1][A].Position;
      posB=Framework[CurrentSystem].Atoms[f1][B].Position;
      posC=Framework[CurrentSystem].Atoms[f1][C].Position;

      Rab.x=posA.x-posB.x;
      Rab.y=posA.y-posB.y;
      Rab.z=posA.z-posB.z;
      Rab=ApplyBoundaryCondition(Rab);
      // form e(-k.r)
      dot_product=k.x*Rab.x+k.y*Rab.y+k.z*Rab.z;
      phase_factor_ab.re=cos(dot_product);
      phase_factor_ab.im=-sin(dot_product);
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
      // form e(-k.r)
      dot_product=k.x*Rbc.x+k.y*Rbc.y+k.z*Rbc.z;
      phase_factor_bc.re=cos(dot_product);
      phase_factor_bc.im=-sin(dot_product);
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
      // form e(-k.r)
      dot_product=k.x*Rac.x+k.y*Rac.y+k.z*Rac.z;
      phase_factor_ac.re=cos(dot_product);
      phase_factor_ac.im=-sin(dot_product);
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
        if(index_i.x>=0) Gradient[index_i.x]+=DF*dtA.x;
        if(index_i.y>=0) Gradient[index_i.y]+=DF*dtA.y;
        if(index_i.z>=0) Gradient[index_i.z]+=DF*dtA.z;

        if(index_j.x>=0) Gradient[index_j.x]+=DF*dtB.x;
        if(index_j.y>=0) Gradient[index_j.y]+=DF*dtB.y;
        if(index_j.z>=0) Gradient[index_j.z]+=DF*dtB.z;

        if(index_k.x>=0) Gradient[index_k.x]+=DF*dtC.x;
        if(index_k.y>=0) Gradient[index_k.y]+=DF*dtC.y;
        if(index_k.z>=0) Gradient[index_k.z]+=DF*dtC.z;

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
        if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x].re+=DDF*dtA.x*dtA.x+D2I.ax;
        if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y].re+=DDF*dtA.x*dtA.y+D2I.ay;
        if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y].re+=DDF*dtA.y*dtA.y+D2I.by;
        if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z].re+=DDF*dtA.x*dtA.z+D2I.az;
        if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z].re+=DDF*dtA.y*dtA.z+D2I.bz;
        if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z].re+=DDF*dtA.z*dtA.z+D2I.cz;

        // Calculate BB-block of the Hessian
        if((index_j.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_j.x][index_j.x].re+=DDF*dtB.x*dtB.x+(D2I.ax+D2K.ax+2.0*D2IK.ax);
        if((index_j.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.x][index_j.y].re+=DDF*dtB.x*dtB.y+(D2I.ay+D2K.ay+D2IK.ay+D2IK.bx);
        if((index_j.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.y][index_j.y].re+=DDF*dtB.y*dtB.y+(D2I.by+D2K.by+2.0*D2IK.by);
        if((index_j.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.x][index_j.z].re+=DDF*dtB.x*dtB.z+(D2I.az+D2K.az+D2IK.az+D2IK.cx);
        if((index_j.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.y][index_j.z].re+=DDF*dtB.y*dtB.z+(D2I.bz+D2K.bz+D2IK.bz+D2IK.cy);
        if((index_j.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.z][index_j.z].re+=DDF*dtB.z*dtB.z+(D2I.cz+D2K.cz+2.0*D2IK.cz);

        // Calculate the CC-block of the Hessian
        if((index_k.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_k.x][index_k.x].re+=DDF*dtC.x*dtC.x+D2K.ax;
        if((index_k.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_k.x][index_k.y].re+=DDF*dtC.x*dtC.y+D2K.ay;
        if((index_k.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_k.y][index_k.y].re+=DDF*dtC.y*dtC.y+D2K.by;
        if((index_k.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.x][index_k.z].re+=DDF*dtC.x*dtC.z+D2K.az;
        if((index_k.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.y][index_k.z].re+=DDF*dtC.y*dtC.z+D2K.bz;
        if((index_k.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_k.z][index_k.z].re+=DDF*dtC.z*dtC.z+D2K.cz;

        // Calculate the AB-block of the Hessian
        if(index1<index2)
        {
          if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x].re+=phase_factor_ab.re*(DDF*dtA.x*dtB.x-D2I.ax-D2IK.ax);
          if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x].im+=phase_factor_ab.im*(DDF*dtA.x*dtB.x-D2I.ax-D2IK.ax);
          if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y].re+=phase_factor_ab.re*(DDF*dtA.x*dtB.y-D2I.ay-D2IK.ay);
          if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y].im+=phase_factor_ab.im*(DDF*dtA.x*dtB.y-D2I.ay-D2IK.ay);
          if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z].re+=phase_factor_ab.re*(DDF*dtA.x*dtB.z-D2I.az-D2IK.az);
          if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z].im+=phase_factor_ab.im*(DDF*dtA.x*dtB.z-D2I.az-D2IK.az);
          if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x].re+=phase_factor_ab.re*(DDF*dtA.y*dtB.x-D2I.ay-D2IK.bx);
          if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x].im+=phase_factor_ab.im*(DDF*dtA.y*dtB.x-D2I.ay-D2IK.bx);
          if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y].re+=phase_factor_ab.re*(DDF*dtA.y*dtB.y-D2I.by-D2IK.by);
          if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y].im+=phase_factor_ab.im*(DDF*dtA.y*dtB.y-D2I.by-D2IK.by);
          if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z].re+=phase_factor_ab.re*(DDF*dtA.y*dtB.z-D2I.bz-D2IK.bz);
          if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z].im+=phase_factor_ab.im*(DDF*dtA.y*dtB.z-D2I.bz-D2IK.bz);
          if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x].re+=phase_factor_ab.re*(DDF*dtA.z*dtB.x-D2I.az-D2IK.cx);
          if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x].im+=phase_factor_ab.im*(DDF*dtA.z*dtB.x-D2I.az-D2IK.cx);
          if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y].re+=phase_factor_ab.re*(DDF*dtA.z*dtB.y-D2I.bz-D2IK.cy);
          if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y].im+=phase_factor_ab.im*(DDF*dtA.z*dtB.y-D2I.bz-D2IK.cy);
          if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z].re+=phase_factor_ab.re*(DDF*dtA.z*dtB.z-D2I.cz-D2IK.cz);
          if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z].im+=phase_factor_ab.im*(DDF*dtA.z*dtB.z-D2I.cz-D2IK.cz);
        }
        else
        {
          if((index_j.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.x][index_i.x].re+=phase_factor_ab.re*(DDF*dtA.x*dtB.x-D2I.ax-D2IK.ax);
          if((index_j.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.x][index_i.x].im+=phase_factor_ab.im*(DDF*dtA.x*dtB.x-D2I.ax-D2IK.ax);
          if((index_j.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.x][index_i.y].re+=phase_factor_ab.re*(DDF*dtA.y*dtB.x-D2I.ay-D2IK.bx);
          if((index_j.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.x][index_i.y].im+=phase_factor_ab.im*(DDF*dtA.y*dtB.x-D2I.ay-D2IK.bx);
          if((index_j.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.x][index_i.z].re+=phase_factor_ab.re*(DDF*dtA.z*dtB.x-D2I.az-D2IK.cx);
          if((index_j.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.x][index_i.z].im+=phase_factor_ab.im*(DDF*dtA.z*dtB.x-D2I.az-D2IK.cx);
          if((index_j.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.y][index_i.x].re+=phase_factor_ab.re*(DDF*dtA.x*dtB.y-D2I.ay-D2IK.ay);
          if((index_j.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.y][index_i.x].im+=phase_factor_ab.im*(DDF*dtA.x*dtB.y-D2I.ay-D2IK.ay);
          if((index_j.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.y][index_i.y].re+=phase_factor_ab.re*(DDF*dtA.y*dtB.y-D2I.by-D2IK.by);
          if((index_j.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.y][index_i.y].im+=phase_factor_ab.im*(DDF*dtA.y*dtB.y-D2I.by-D2IK.by);
          if((index_j.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.y][index_i.z].re+=phase_factor_ab.re*(DDF*dtA.z*dtB.y-D2I.bz-D2IK.cy);
          if((index_j.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.y][index_i.z].im+=phase_factor_ab.im*(DDF*dtA.z*dtB.y-D2I.bz-D2IK.cy);
          if((index_j.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.z][index_i.x].re+=phase_factor_ab.re*(DDF*dtA.x*dtB.z-D2I.az-D2IK.az);
          if((index_j.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_j.z][index_i.x].im+=phase_factor_ab.im*(DDF*dtA.x*dtB.z-D2I.az-D2IK.az);
          if((index_j.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.z][index_i.y].re+=phase_factor_ab.re*(DDF*dtA.y*dtB.z-D2I.bz-D2IK.bz);
          if((index_j.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_j.z][index_i.y].im+=phase_factor_ab.im*(DDF*dtA.y*dtB.z-D2I.bz-D2IK.bz);
          if((index_j.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.z][index_i.z].re+=phase_factor_ab.re*(DDF*dtA.z*dtB.z-D2I.cz-D2IK.cz);
          if((index_j.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_j.z][index_i.z].im+=phase_factor_ab.im*(DDF*dtA.z*dtB.z-D2I.cz-D2IK.cz);
        }

        // Calculate the AC-block of the Hessian
        if(index1<index3)
        {
          if((index_i.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.x][index_k.x].re+=phase_factor_ac.re*(DDF*dtA.x*dtC.x+D2IK.ax);
          if((index_i.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.x][index_k.x].im+=phase_factor_ac.im*(DDF*dtA.x*dtC.x+D2IK.ax);
          if((index_i.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.x][index_k.y].re+=phase_factor_ac.re*(DDF*dtA.x*dtC.y+D2IK.ay);
          if((index_i.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.x][index_k.y].im+=phase_factor_ac.im*(DDF*dtA.x*dtC.y+D2IK.ay);
          if((index_i.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.x][index_k.z].re+=phase_factor_ac.re*(DDF*dtA.x*dtC.z+D2IK.az);
          if((index_i.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.x][index_k.z].im+=phase_factor_ac.im*(DDF*dtA.x*dtC.z+D2IK.az);
          if((index_i.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.y][index_k.x].re+=phase_factor_ac.re*(DDF*dtA.y*dtC.x+D2IK.bx);
          if((index_i.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.y][index_k.x].im+=phase_factor_ac.im*(DDF*dtA.y*dtC.x+D2IK.bx);
          if((index_i.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.y][index_k.y].re+=phase_factor_ac.re*(DDF*dtA.y*dtC.y+D2IK.by);
          if((index_i.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.y][index_k.y].im+=phase_factor_ac.im*(DDF*dtA.y*dtC.y+D2IK.by);
          if((index_i.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.y][index_k.z].re+=phase_factor_ac.re*(DDF*dtA.y*dtC.z+D2IK.bz);
          if((index_i.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.y][index_k.z].im+=phase_factor_ac.im*(DDF*dtA.y*dtC.z+D2IK.bz);
          if((index_i.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.z][index_k.x].re+=phase_factor_ac.re*(DDF*dtA.z*dtC.x+D2IK.cx);
          if((index_i.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_i.z][index_k.x].im+=phase_factor_ac.im*(DDF*dtA.z*dtC.x+D2IK.cx);
          if((index_i.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.z][index_k.y].re+=phase_factor_ac.re*(DDF*dtA.z*dtC.y+D2IK.cy);
          if((index_i.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_i.z][index_k.y].im+=phase_factor_ac.im*(DDF*dtA.z*dtC.y+D2IK.cy);
          if((index_i.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.z][index_k.z].re+=phase_factor_ac.re*(DDF*dtA.z*dtC.z+D2IK.cz);
          if((index_i.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_i.z][index_k.z].im+=phase_factor_ac.im*(DDF*dtA.z*dtC.z+D2IK.cz);
        }
        else 
        {
          if((index_k.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.x][index_i.x].re+=phase_factor_ac.re*(DDF*dtA.x*dtC.x+D2IK.ax);
          if((index_k.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.x][index_i.x].im+=phase_factor_ac.im*(DDF*dtA.x*dtC.x+D2IK.ax);
          if((index_k.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.x][index_i.y].re+=phase_factor_ac.re*(DDF*dtA.y*dtC.x+D2IK.bx);
          if((index_k.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.x][index_i.y].im+=phase_factor_ac.im*(DDF*dtA.y*dtC.x+D2IK.bx);
          if((index_k.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.x][index_i.z].re+=phase_factor_ac.re*(DDF*dtA.z*dtC.x+D2IK.cx);
          if((index_k.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.x][index_i.z].im+=phase_factor_ac.im*(DDF*dtA.z*dtC.x+D2IK.cx);
          if((index_k.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.y][index_i.x].re+=phase_factor_ac.re*(DDF*dtA.x*dtC.y+D2IK.ay);
          if((index_k.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.y][index_i.x].im+=phase_factor_ac.im*(DDF*dtA.x*dtC.y+D2IK.ay);
          if((index_k.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.y][index_i.y].re+=phase_factor_ac.re*(DDF*dtA.y*dtC.y+D2IK.by);
          if((index_k.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.y][index_i.y].im+=phase_factor_ac.im*(DDF*dtA.y*dtC.y+D2IK.by);
          if((index_k.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.y][index_i.z].re+=phase_factor_ac.re*(DDF*dtA.z*dtC.y+D2IK.cy);
          if((index_k.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.y][index_i.z].im+=phase_factor_ac.im*(DDF*dtA.z*dtC.y+D2IK.cy);
          if((index_k.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.z][index_i.x].re+=phase_factor_ac.re*(DDF*dtA.x*dtC.z+D2IK.az);
          if((index_k.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_k.z][index_i.x].im+=phase_factor_ac.im*(DDF*dtA.x*dtC.z+D2IK.az);
          if((index_k.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.z][index_i.y].re+=phase_factor_ac.re*(DDF*dtA.y*dtC.z+D2IK.bz);
          if((index_k.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_k.z][index_i.y].im+=phase_factor_ac.im*(DDF*dtA.y*dtC.z+D2IK.bz);
          if((index_k.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.z][index_i.z].re+=phase_factor_ac.re*(DDF*dtA.z*dtC.z+D2IK.cz);
          if((index_k.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_k.z][index_i.z].im+=phase_factor_ac.im*(DDF*dtA.z*dtC.z+D2IK.cz);
        }
        
        // Calculate the BC-block of the Hessian
        if(index3<index2)
        {
          if((index_k.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.x][index_j.x].re+=phase_factor_bc.re*(DDF*dtC.x*dtB.x-D2K.ax-D2IK.ax);
          if((index_k.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.x][index_j.x].im+=phase_factor_bc.im*(DDF*dtC.x*dtB.x-D2K.ax-D2IK.ax);
          if((index_k.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.x][index_j.y].re+=phase_factor_bc.re*(DDF*dtC.x*dtB.y-D2K.ay-D2IK.bx);
          if((index_k.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.x][index_j.y].im+=phase_factor_bc.im*(DDF*dtC.x*dtB.y-D2K.ay-D2IK.bx);
          if((index_k.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.x][index_j.z].re+=phase_factor_bc.re*(DDF*dtC.x*dtB.z-D2K.az-D2IK.cx);
          if((index_k.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.x][index_j.z].im+=phase_factor_bc.im*(DDF*dtC.x*dtB.z-D2K.az-D2IK.cx);
          if((index_k.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.y][index_j.x].re+=phase_factor_bc.re*(DDF*dtC.y*dtB.x-D2K.ay-D2IK.ay);
          if((index_k.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.y][index_j.x].im+=phase_factor_bc.im*(DDF*dtC.y*dtB.x-D2K.ay-D2IK.ay);
          if((index_k.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.y][index_j.y].re+=phase_factor_bc.re*(DDF*dtC.y*dtB.y-D2K.by-D2IK.by);
          if((index_k.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.y][index_j.y].im+=phase_factor_bc.im*(DDF*dtC.y*dtB.y-D2K.by-D2IK.by);
          if((index_k.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.y][index_j.z].re+=phase_factor_bc.re*(DDF*dtC.y*dtB.z-D2K.bz-D2IK.cy);
          if((index_k.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.y][index_j.z].im+=phase_factor_bc.im*(DDF*dtC.y*dtB.z-D2K.bz-D2IK.cy);
          if((index_k.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.z][index_j.x].re+=phase_factor_bc.re*(DDF*dtC.z*dtB.x-D2K.az-D2IK.az);
          if((index_k.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_k.z][index_j.x].im+=phase_factor_bc.im*(DDF*dtC.z*dtB.x-D2K.az-D2IK.az);
          if((index_k.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.z][index_j.y].re+=phase_factor_bc.re*(DDF*dtC.z*dtB.y-D2K.bz-D2IK.bz);
          if((index_k.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_k.z][index_j.y].im+=phase_factor_bc.im*(DDF*dtC.z*dtB.y-D2K.bz-D2IK.bz);
          if((index_k.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.z][index_j.z].re+=phase_factor_bc.re*(DDF*dtC.z*dtB.z-D2K.cz-D2IK.cz);
          if((index_k.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_k.z][index_j.z].im+=phase_factor_bc.im*(DDF*dtC.z*dtB.z-D2K.cz-D2IK.cz);
        }
        else 
        {
          if((index_j.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.x][index_k.x].re+=phase_factor_bc.re*(DDF*dtC.x*dtB.x-D2K.ax-D2IK.ax);
          if((index_j.x>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.x][index_k.x].im+=phase_factor_bc.im*(DDF*dtC.x*dtB.x-D2K.ax-D2IK.ax);
          if((index_j.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.x][index_k.y].re+=phase_factor_bc.re*(DDF*dtC.y*dtB.x-D2K.ay-D2IK.ay);
          if((index_j.x>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.x][index_k.y].im+=phase_factor_bc.im*(DDF*dtC.y*dtB.x-D2K.ay-D2IK.ay);
          if((index_j.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.x][index_k.z].re+=phase_factor_bc.re*(DDF*dtC.z*dtB.x-D2K.az-D2IK.az);
          if((index_j.x>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.x][index_k.z].im+=phase_factor_bc.im*(DDF*dtC.z*dtB.x-D2K.az-D2IK.az);
          if((index_j.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.y][index_k.x].re+=phase_factor_bc.re*(DDF*dtC.x*dtB.y-D2K.ay-D2IK.bx);
          if((index_j.y>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.y][index_k.x].im+=phase_factor_bc.im*(DDF*dtC.x*dtB.y-D2K.ay-D2IK.bx);
          if((index_j.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.y][index_k.y].re+=phase_factor_bc.re*(DDF*dtC.y*dtB.y-D2K.by-D2IK.by);
          if((index_j.y>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.y][index_k.y].im+=phase_factor_bc.im*(DDF*dtC.y*dtB.y-D2K.by-D2IK.by);
          if((index_j.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.y][index_k.z].re+=phase_factor_bc.re*(DDF*dtC.z*dtB.y-D2K.bz-D2IK.bz);
          if((index_j.y>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.y][index_k.z].im+=phase_factor_bc.im*(DDF*dtC.z*dtB.y-D2K.bz-D2IK.bz);
          if((index_j.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.z][index_k.x].re+=phase_factor_bc.re*(DDF*dtC.x*dtB.z-D2K.az-D2IK.cx);
          if((index_j.z>=0)&&(index_k.x>=0)) HessianMatrix.element[index_j.z][index_k.x].im+=phase_factor_bc.im*(DDF*dtC.x*dtB.z-D2K.az-D2IK.cx);
          if((index_j.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.z][index_k.y].re+=phase_factor_bc.re*(DDF*dtC.y*dtB.z-D2K.bz-D2IK.cy);
          if((index_j.z>=0)&&(index_k.y>=0)) HessianMatrix.element[index_j.z][index_k.y].im+=phase_factor_bc.im*(DDF*dtC.y*dtB.z-D2K.bz-D2IK.cy);
          if((index_j.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.z][index_k.z].re+=phase_factor_bc.re*(DDF*dtC.z*dtB.z-D2K.cz-D2IK.cz);
          if((index_j.z>=0)&&(index_k.z>=0)) HessianMatrix.element[index_j.z][index_k.z].im+=phase_factor_bc.im*(DDF*dtC.z*dtB.z-D2K.cz-D2IK.cz);
        }
      }
    }
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



void ComputeFrameworkTorsionPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
}

void ComputeFrameworkImproperTorsionPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian)
{
}
