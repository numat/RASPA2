/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework_energy.c' is part of RASPA.

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
#include "simulation.h"
#include "ewald.h"
#include "potentials.h"
#include "utils.h"
#include "output.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "spacegroup.h"

// the routines for VDW and charge-charge interaction are splitted. The VDW interactions can potentially shift the
// actual interaction site across a bond (see e.g. MM3).

// Calculate the Energy of a bead of type "typeA" at a position "posA".
// posA is defined in xyz (also for triclinic).
REAL CalculateFrameworkVDWEnergyAtPosition(POINT posA,int typeA)
{
  int i,j,typeB,f1;
  REAL rr;
  REAL UVDW;
  VECTOR posB,dr,s;
  int icell,icell0;
  int ncell;

  UVDW=0.0;

  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))   // grid-interpolating for rigid frameworks
    {
      UVDW=InterpolateVDWGrid(typeA,posA);
      return UVDW;
    }
    else if(UseCellLists[CurrentSystem]) // energy using cell-lists
    {
      // convert from xyz to abc
      s.x=InverseBox[CurrentSystem].ax*posA.x+InverseBox[CurrentSystem].bx*posA.y+InverseBox[CurrentSystem].cx*posA.z;
      s.y=InverseBox[CurrentSystem].ay*posA.x+InverseBox[CurrentSystem].by*posA.y+InverseBox[CurrentSystem].cy*posA.z;
      s.z=InverseBox[CurrentSystem].az*posA.x+InverseBox[CurrentSystem].bz*posA.y+InverseBox[CurrentSystem].cz*posA.z;

      // apply boundary condition
      s.x-=(REAL)NINT(s.x);
      s.y-=(REAL)NINT(s.y);
      s.z-=(REAL)NINT(s.z);

      // s between 0 and 1
      s.x+=0.5;
      s.y+=0.5;
      s.z+=0.5;

      // compute the corresponding cell-id
      icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
             ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
             ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        // loop over cells
        for(i=0;i<27;i++)
        {
          icell=CellListMap[CurrentSystem][icell0][i];

          j=Framework[CurrentSystem].CellListHead[f1][icell];

          while(j>=0)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
            posB=Framework[CurrentSystem].Atoms[f1][j].AnisotropicPosition;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
              UVDW+=PotentialValue(typeA,typeB,rr);

            j=Framework[CurrentSystem].CellList[f1][j];
          }
        }
      }
    }
    else if(UseReplicas[CurrentSystem]) // energy using replicas
    {
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][i].Type;
          posB=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
              UVDW+=PotentialValue(typeA,typeB,rr);
          }
        }
      }
    }
    else // regular energy calculation
    {
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][i].Type;
          posB=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffVDWSquared)
            UVDW+=PotentialValue(typeA,typeB,rr);
        }
      }
    }
  }
  return UVDW;
}


void CalculateFrameworkChargeEnergyAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UchargeBondDipole)
{
  int i,typeB,f1,k;
  int B1,B2;
  REAL r,rr,chargeA,chargeB;
  REAL Bt1;
  REAL cosB;
  REAL DipoleMagnitudeB,temp;
  VECTOR posB,posB1,posB2,dr,dipoleB;
  REAL SwitchingValue,TranslationValue;
  REAL energy;
  int ncell;

  (*UChargeCharge)=0.0;
  (*UchargeBondDipole)=0.0;

  if(ChargeMethod==NONE) return;

  if(!PseudoAtoms[typeA].HasCharges) return;

  chargeA=PseudoAtoms[typeA].Charge1;
  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
    {
      (*UChargeCharge)=InterpolateCoulombGrid(typeA,posA);
    }
    else if(UseReplicas[CurrentSystem])
    {
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][i].Type;
          posB=Framework[CurrentSystem].Atoms[f1][i].Position;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              chargeB=Framework[CurrentSystem].Atoms[f1][i].Charge;
 
              switch(ChargeMethod)
              {
                case NONE:
                  break;
                case TRUNCATED_COULOMB:
                  (*UChargeCharge)+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  break;
               case SHIFTED_COULOMB:
                 (*UChargeCharge)+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                 break;
               case SMOOTHED_COULOMB:
                 energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                 if(rr>CutOffChargeChargeSwitchSquared)
                 {
                   SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                  SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                   TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                   (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                    SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                    SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                   energy=energy*SwitchingValue+TranslationValue;
                 }
                 (*UChargeCharge)+=energy;
                 break;
               case WOLFS_METHOD_DAMPED_FG:
                 (*UChargeCharge)+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                         -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                          (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                          (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                 break;
               case EWALD:
                 (*UChargeCharge)+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                         (erfc(Alpha[CurrentSystem]*r)/r);
                 break;
               default:
                 printf("Unknown charge-method in 'CalculateFrameworkChargeEnergyAtPosition'\n");
                 exit(0);
                 break;
              }
            }
          }
        }

        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyReplicaBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt1=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt1=1.0/(r*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                     SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                    SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkChargeEnergyAtPosition'\n");
                  exit(0);
                  break;
               }
 
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                (*UchargeBondDipole)+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
            }
          }
        }
      }
    }
    else
    {
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        {
          typeB=Framework[CurrentSystem].Atoms[f1][i].Type;
          posB=Framework[CurrentSystem].Atoms[f1][i].Position;
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared)
          {
            r=sqrt(rr);
            chargeB=Framework[CurrentSystem].Atoms[f1][i].Charge;

            switch(ChargeMethod)
            {
              case NONE:
                break;
              case TRUNCATED_COULOMB:
                (*UChargeCharge)+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                break;
              case SHIFTED_COULOMB:
                (*UChargeCharge)+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                break;
              case SMOOTHED_COULOMB:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                if(rr>CutOffChargeChargeSwitchSquared)
                {
                  SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                 SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                  TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                  (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                   SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                   SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                  energy=energy*SwitchingValue+TranslationValue;
                }
                (*UChargeCharge)+=energy;
                break;
              case WOLFS_METHOD_DAMPED_FG:
                (*UChargeCharge)+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                         -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                          (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                          (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                break;
              case EWALD:
                (*UChargeCharge)+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                             (erfc(Alpha[CurrentSystem]*r)/r);
                break;
              default:
                printf("Unknown charge-method in 'CalculateFrameworkChargeEnergyAtPosition'\n");
                exit(0);
                break;
            }
          }
        }

        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt1=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt1=1.0/(r*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkChargeEnergyAtPosition'\n");
                exit(0);
                break;
            }

            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            (*UchargeBondDipole)+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
          }
        }
      }
    }
  }
}

void CalculateFrameworkBondDipoleEnergyAtPosition(VECTOR posA1,VECTOR posA2,REAL DipoleMagnitudeA,REAL *UchargeBondDipole,REAL *UBondDipoleBondDipole)
{
  int k,f1;
  int B1,B2;
  int TypeB;
  REAL r;
  VECTOR posA,posB,posB1,posB2,dr;
  REAL DipoleMagnitudeB,chargeB;
  REAL ri2,rk2,rr,length,temp;
  VECTOR dipoleA,dipoleB;
  REAL cosAB,cosA,cosB;
  REAL Bt0,Bt1,Bt2;
  REAL SwitchingValue;
  int ncell;

  (*UchargeBondDipole)=0.0;
  (*UBondDipoleBondDipole)=0.0;

  dipoleA.x=posA2.x-posA1.x;
  dipoleA.y=posA2.y-posA1.y;
  dipoleA.z=posA2.z-posA1.z;
  posA.x=posA1.x+0.5*dipoleA.x;
  posA.y=posA1.y+0.5*dipoleA.y;
  posA.z=posA1.z+0.5*dipoleA.z;
  ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
  length=sqrt(ri2);
  temp=DipoleMagnitudeA/length;
  dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

  if(UseReplicas[CurrentSystem])
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
      {
        TypeB=Framework[CurrentSystem].Atoms[f1][k].Type;

        if(PseudoAtoms[TypeB].HasCharges)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].x);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt1=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt1=1.0/(r*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                     SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                    SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  break;
               default:
                 printf("Unknown charge-bonddipole method in 'CalculateFrameworkBondDipoleEnergyAtPosition'\n");
                 exit(0);
                 break;
              }

              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              (*UchargeBondDipole)-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
            }
          }
        }
      }
 
      for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
      {
        DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
        B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
        B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
        posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
        posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
        length=sqrt(rk2);
        temp=DipoleMagnitudeB/length;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
          dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
          dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
          dr=ApplyReplicaBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkBondDipoleEnergyAtPosition'\n");
                exit(0);
                break;
            }

            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            (*UBondDipoleBondDipole)+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
          }
        }
      }
    }
  }
  else
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
      {
        TypeB=Framework[CurrentSystem].Atoms[f1][k].Type;

        if(PseudoAtoms[TypeB].HasCharges)
        {
          posB=Framework[CurrentSystem].Atoms[f1][k].Position;
          chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt1=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt1=1.0/(r*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                   SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                  SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkBondDipoleEnergyAtPosition'\n");
                exit(0);
                break;
            }

            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            (*UchargeBondDipole)-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
          }
        }
      }
 
      for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
      {
        DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
        B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
        B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
        posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
        posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
        length=sqrt(rk2);
        temp=DipoleMagnitudeB/length;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffBondDipoleBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=Bt2=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              Bt2=3.0/(r*rr*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              Bt2=3.0/(r*rr*rr);
              if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
              {
                SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                               SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
              }
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                  4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
              break;
            default:
              printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkBondDipoleEnergyAtPosition'\n");
              exit(0);
              break;
          }

          cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
          cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          (*UBondDipoleBondDipole)+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
        }
      }
    }
  }
}


// Monte-Carlo full energy routines
// ================================

REAL CalculateFrameworkBondEnergy(int flag,int f2,int atom_id)
{
  int i,f1;
  REAL r,rr,r1,U;
  REAL temp,temp2,exp_term;
  REAL *parms,UHostBond;
  VECTOR posA,posB,dr;
  int A,B;

  UHostBond=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBonds[f1];i++)
      {
        A=Framework[CurrentSystem].Bonds[f1][i].A;
        B=Framework[CurrentSystem].Bonds[f1][i].B;
  
        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id)))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
  
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
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
              break;
            case CORE_SHELL_SPRING:
              U=0.5*parms[0]*SQR(r);
              break;
            case MORSE_BOND:
              // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
              // ===============================================
              // p_0/k_B [K]       force constant
              // p_1     [A^-1]    parameter
              // p_2     [A]       reference bond distance
              temp=exp(parms[1]*(parms[2]-r));
              U=parms[0]*(SQR(1.0-temp)-1.0);
              break;
            case LJ_12_6_BOND:
              // A/r_ij^12-B/r_ij^6
              // ===============================================
              // p_0/k_B [K A^12]
              // p_1/k_B [K A^6]
              temp=CUBE(1.0/rr);
              U=parms[0]*SQR(temp)-parms[1]*temp;
              break;
            case LENNARD_JONES_BOND:
              // 4*p_0*((p_1/r)^12-(p_1/r)^6)
              // ===============================================
              // p_0/k_B [K]
              // p_1     [A]
              temp=CUBE(parms[1]/rr);
              U=4.0*parms[0]*(temp*(temp-1.0));
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
              break;
            case RESTRAINED_HARMONIC_BOND:
              // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
              // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
              // ===============================================
              // p_0/k_B [K/A^2]
              // p_1     [A]
              // p_2     [A]
              r1=r-parms[1];
              U=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
                    +parms[0]*parms[2]*MAX2(fabs(r1)-parms[2],(REAL)0.0);
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
              break;
            case MM3_BOND:
              // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
              // =================================================================
              // p_0     [mdyne/A molecule]
              // p_1     [A]
              temp=r-parms[1];
              temp2=SQR(temp);
              U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
              break;
            case MEASURE_BOND:
              U=0.0;
              break;
            case RIGID_BOND:
              U=0.0;
              break;
            case FIXED_BOND:
              U=0.0;
              break;
            default:
              printf("Undefined Bond potential in routine 'CalculateFrameworkBondEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          // add contribution to the Adsorbate stretch energy
          UHostBond+=U;
        }
      }
    }
  }
  return UHostBond;
}

REAL CalculateFrameworkUreyBradleyEnergy(int flag,int f2,int atom_id)
{
  int i,f1;
  REAL r,rr,r1,U;
  REAL temp,temp2,exp_term;
  VECTOR dr,posA,posC;
  REAL *parms,UHostUreyBradley;
  int A,C;

  UHostUreyBradley=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfUreyBradleys[f1];i++)
      {
        A=Framework[CurrentSystem].UreyBradleys[f1][i].A;
        C=Framework[CurrentSystem].UreyBradleys[f1][i].C;
  
        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||C==atom_id)))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posC=Framework[CurrentSystem].Atoms[f1][C].Position;
  
          dr.x=posA.x-posC.x;
          dr.y=posA.y-posC.y;
          dr.z=posA.z-posC.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          r=sqrt(rr);
  
          parms=(REAL*)&Framework[CurrentSystem].UreyBradleyArguments[f1][i];
  
          switch(Framework[CurrentSystem].UreyBradleyType[f1][i])
          {
            case HARMONIC_UREYBRADLEY:
              // 0.5*p0*SQR(r-p1);
              // ===============================================
              // p_0/k_B [K/A^2]   force constant
              // p_1     [A]       reference bond distance
              U=0.5*parms[0]*SQR(r-parms[1]);
              break;
            case MORSE_UREYBRADLEY:
              // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
              // ===============================================
              // p_0/k_B [K]       force constant
              // p_1     [A^-1]    parameter
              // p_2     [A]       reference bond distance
              temp=exp(parms[1]*(parms[2]-r));
              U=parms[0]*(SQR(1.0-temp)-1.0);
              break;
            case LJ_12_6_UREYBRADLEY:
              // A/r_ij^12-B/r_ij^6
              // ===============================================
              // p_0/k_B [K A^12]
              // p_1/k_B [K A^6]
              temp=CUBE(1.0/rr);
              U=parms[0]*SQR(temp)-parms[1]*temp;
              break;
            case LENNARD_JONES_UREYBRADLEY:
              // 4*p_0*((p_1/r)^12-(p_1/r)^6)
              // ===============================================
              // p_0/k_B [K]
              // p_1     [A]
              temp=CUBE(parms[1]/rr);
              U=4.0*parms[0]*(temp*(temp-1.0));
              break;
            case BUCKINGHAM_UREYBRADLEY:
              // p_0*exp(-p_1 r)-p_2/r^6
              // ===============================================
              // p_0/k_B [K]
              // p_1     [A^-1]
              // p_2/k_B [K A^6]
              temp=parms[2]*CUBE(1.0/rr);
              exp_term=parms[0]*exp(-parms[1]*r);
              U=-temp+exp_term;
              break;
            case RESTRAINED_HARMONIC_UREYBRADLEY:
              // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
              // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
              // ===============================================
              // p_0/k_B [K/A^2]
              // p_1     [A]
              // p_2     [A]
              r1=r-parms[1];
              U=0.5*parms[0]*SQR(MIN2(fabs(r1),parms[2]))
                    +parms[0]*parms[2]*MAX2(fabs(r1)-parms[2],(REAL)0.0);
              break;
            case QUARTIC_UREYBRADLEY:
              // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
              // ===========================================================
              // p_0/k_B [K/A^2]
              // p_1     [A]
              // p_2/k_B [K/A^3]
              // p_3/k_B [K/A^4]
              temp=r-parms[1];
              temp2=SQR(r-parms[1]);
              U=0.5*parms[0]*temp2+(1.0/3.0)*parms[2]*temp*temp2+0.25*parms[3]*SQR(temp2);
              break;
            case CFF_QUARTIC_UREYBRADLEY:
              // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
              // ===============================================
              // p_0/k_B [K/A^2]
              // p_1     [A]
              // p_2/k_B [K/A^3]
              // p_3/k_B [K/A^4]
              temp=r-parms[1];
              temp2=SQR(r-parms[1]);
              U=parms[0]*temp2+parms[2]*temp*temp2+parms[3]*SQR(temp2);
              break;
            case MM3_UREYBRADLEY:
              // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
              // =================================================================
              // p_0     [mdyne/A molecule]
              // p_1     [A]
              temp=r-parms[1];
              temp2=SQR(temp);
              U=parms[0]*temp2*(1.0-2.55*temp+(7.0/12.0)*SQR(2.55)*temp2);
              break;
            case RIGID_UREYBRADLEY:
              U=0.0;
              break;
            case FIXED_UREYBRADLEY:
              U=0.0;
              break;
            default:
              printf("Undefined Urey-Bradley potential in routine 'CalculateFrameworkUreyBradleyEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          // add contribution to the Adsorbate Urey-Bradley energy
          UHostUreyBradley+=U;
        }
      }
    }
  }
  return UHostUreyBradley;
}

REAL CalculateFrameworkBendEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  REAL *parms,U,temp,temp2;
  REAL CosTheta,Theta,SinTheta;
  REAL rab,rbc,rac,DTDX,UHostBend;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rac;
  VECTOR Rad,Rbd,Rcd,t,ip,ap,cp;
  REAL delta,rt2,rap2,rcp2;

  UHostBend=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBends[f1];i++)
      {
        A=Framework[CurrentSystem].Bends[f1][i].A;
        B=Framework[CurrentSystem].Bends[f1][i].B;
        C=Framework[CurrentSystem].Bends[f1][i].C;
        D=Framework[CurrentSystem].Bends[f1][i].D;
  
        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&((A==atom_id)||(B==atom_id)||(C==atom_id)||(D==atom_id))))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          posC=Framework[CurrentSystem].Atoms[f1][C].Position;
  
          switch(Framework[CurrentSystem].BendType[f1][i])
          {
            case MM3_IN_PLANE_BEND:
              D=Framework[CurrentSystem].Bends[f1][i].D;
              posD=Framework[CurrentSystem].Atoms[f1][D].Position;
              Rad.x=posA.x-posD.x;
              Rad.y=posA.y-posD.y;
              Rad.z=posA.z-posD.z;
              Rad=ApplyBoundaryCondition(Rad);
  
              Rbd.x=posB.x-posD.x;
              Rbd.y=posB.y-posD.y;
              Rbd.z=posB.z-posD.z;
              Rbd=ApplyBoundaryCondition(Rbd);
  
              Rcd.x=posC.x-posD.x;
              Rcd.y=posC.y-posD.y;
              Rcd.z=posC.z-posD.z;
              Rcd=ApplyBoundaryCondition(Rcd);
  
              t.x=Rad.y*Rcd.z-Rad.z*Rcd.y;
              t.y=Rad.z*Rcd.x-Rad.x*Rcd.z;
              t.z=Rad.x*Rcd.y-Rad.y*Rcd.x;
              rt2=t.x*t.x+t.y*t.y+t.z*t.z;
              delta=-(t.x*Rbd.x+t.y*Rbd.y+t.z*Rbd.z)/rt2;
  
              ip.x=posB.x+t.x*delta;
              ip.y=posB.y+t.y*delta;
              ip.z=posB.z+t.z*delta;
              ap.x=posA.x-ip.x;
              ap.y=posA.y-ip.y;
              ap.z=posA.z-ip.z;
              cp.x=posC.x-ip.x;
              cp.y=posC.y-ip.y;
              cp.z=posC.z-ip.z;
              ap=ApplyBoundaryCondition(ap);
              cp=ApplyBoundaryCondition(cp);
  
              rap2=ap.x*ap.x+ap.y*ap.y+ap.z*ap.z;
              rcp2=cp.x*cp.x+cp.y*cp.y+cp.z*cp.z;
  
              CosTheta=(ap.x*cp.x+ap.y*cp.y+ap.z*cp.z)/sqrt(rap2*rcp2);
            break;
            default:
              Rab.x=posA.x-posB.x;
              Rab.y=posA.y-posB.y;
              Rab.z=posA.z-posB.z;
              Rab=ApplyBoundaryCondition(Rab);
              rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
              Rab.x/=rab;
              Rab.y/=rab;
              Rab.z/=rab;
  
              Rbc.x=posC.x-posB.x;
              Rbc.y=posC.y-posB.y;
              Rbc.z=posC.z-posB.z;
              Rbc=ApplyBoundaryCondition(Rbc);
              rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
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
              break;
          }
  
          CosTheta=MIN2(1.0,MAX2(-1.0,CosTheta));
          Theta=acos(CosTheta);
          SinTheta=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta)));
          DTDX=-1.0/sqrt(1.0-SQR(CosTheta));
  
          parms=Framework[CurrentSystem].BendArguments[f1][i];
  
          switch(Framework[CurrentSystem].BendType[f1][i])
          {
            case HARMONIC_BEND:
              // (1/2)p_0*(theta-p_1)^2
              // ===============================================
              // p_0/k_B [K/rad^2]
              // p_1     [degrees]
              U=0.5*parms[0]*SQR(Theta-parms[1]);
              break;
            case CORE_SHELL_BEND:
              // (1/2)p_0*(theta-p_1)^2
              // ===============================================
              // p_0/k_B [K/rad^2]
              // p_1     [degrees]
              U=0.5*parms[0]*SQR(Theta-parms[1]);
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
              break;
            case HARMONIC_COSINE_BEND:
              // (1/2)*p_0*(cos(theta)-cos(p_1))^2
              // ===============================================
              // p_0/k_B [K]
              // p_1     [degrees]
              U=0.5*parms[0]*SQR(CosTheta-parms[1]);
              break;
            case COSINE_BEND:
              // p_0*(1+cos(p_1*theta-p_2))
              // ===============================================
              // p_0/k_B [K]
              // p_1     [-]
              // p_2     [degrees]
              temp=parms[1]*Theta-parms[2];
              U=parms[0]*(1.0+cos(temp));
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
              break;
            case FIXED_BEND:
              U=0.0;
              break;
            case MEASURE_BEND:
              U=0.0;
              break;
            default:
              printf("Undefined Bend potential in routine 'CalculateFrameworkBendEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          // add contribution to the energy
          UHostBend+=U;
        }
      }
    }
  }
  return UHostBend;
}

REAL CalculateFrameworkInversionBendEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  REAL *parms;
  REAL c,e;
  REAL rrab,rab2;
  REAL CosChi,Chi,energy;
  REAL temp,temp2,UHostInversionBend;
  VECTOR Rab,Rbc,Rbd,Rcd,Rad;
  POINT posA,posB,posC,posD;

  UHostInversionBend=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfInversionBends[f1];i++)
      {
        A=Framework[CurrentSystem].InversionBends[f1][i].A;
        B=Framework[CurrentSystem].InversionBends[f1][i].B;
        C=Framework[CurrentSystem].InversionBends[f1][i].C;
        D=Framework[CurrentSystem].InversionBends[f1][i].D;
  
        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          posC=Framework[CurrentSystem].Atoms[f1][C].Position;
          posD=Framework[CurrentSystem].Atoms[f1][D].Position;
  
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
  
          Rbd.x=posD.x-posB.x;
          Rbd.y=posD.y-posB.y;
          Rbd.z=posD.z-posB.z;
          Rbd=ApplyBoundaryCondition(Rbd);
  
          Rcd.x=posD.x-posC.x;
          Rcd.y=posD.y-posC.y;
          Rcd.z=posD.z-posC.z;
          Rcd=ApplyBoundaryCondition(Rcd);
  
          Rad.x=posD.x-posA.x;
          Rad.y=posD.y-posA.y;
          Rad.z=posD.z-posA.z;
          Rad=ApplyBoundaryCondition(Rad);
  
          switch(Framework[CurrentSystem].InversionBendType[f1][i])
          {
            case HARMONIC_INVERSION:
            case HARMONIC_COSINE_INVERSION:
            case PLANAR_INVERSION:
              // w is a vector perpendicular to the B-C-D plane
              // c=w.w=(Rbc x Rbd).(Rbc x Rbd)= r_bc^2 r_bd^2 - (r_cb . r_bd)^2
              c=(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z)*(Rbd.x*Rbd.x+Rbd.y*Rbd.y+Rbd.z*Rbd.z)-SQR(Rbc.x*Rbd.x+Rbc.y*Rbd.y+Rbc.z*Rbd.z);
              break;
            case HARMONIC_INVERSION2:
            case HARMONIC_COSINE_INVERSION2:
            case PLANAR_INVERSION2:
            case MM3_INVERSION:
              // w is a vector perpendicular to the A-C-D plane
              // c=w.w=(Rcd x Rad).(Rcd x Rad)=r_cd^2 r_ad^2 - (r_da . r_cd)^2
              c=(Rcd.x*Rcd.x+Rcd.y*Rcd.y+Rcd.z*Rcd.z)*(Rad.x*Rad.x+Rad.y*Rad.y+Rad.z*Rad.z)-SQR(Rad.x*Rcd.x+Rad.y*Rcd.y+Rad.z*Rcd.z);
              break;
            default:
              printf("Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          e=Rab.x*(Rbd.y*Rbc.z-Rbd.z*Rbc.y)+Rab.y*(Rbd.z*Rbc.x-Rbd.x*Rbc.z)+Rab.z*(Rbd.x*Rbc.y-Rbd.y*Rbc.x);
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
              break;
            case HARMONIC_COSINE_INVERSION:
            case HARMONIC_COSINE_INVERSION2:
              // (1/2)*p_0*(cos(phi)-cos(p_1))^2
              // ===============================================
              // p_0/k_B [K]
              // p_1     [degrees]
              energy=0.5*parms[0]*SQR(CosChi-parms[1]);
              break;
            case PLANAR_INVERSION:
            case PLANAR_INVERSION2:
              // (1/2)*p_0*(1-cos(phi))
              // ===============================================
              // p_0/k_B [K]
              energy=parms[0]*(1.0-CosChi);
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
              break;
            default:
              printf("Undefined Inversion-Bend potential in routine 'CalculateFrameworkInversionBendEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          // energy
          UHostInversionBend+=energy;
        }
      }
    }
  }
  return UHostInversionBend;
}

REAL CalculateFrameworkTorsionEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL rbc,U;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi,UHostTorsion;
  VECTOR Pb,Pc;
  REAL *parms;

  UHostTorsion=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[f1];i++)
      {
        A=Framework[CurrentSystem].Torsions[f1][i].A;
        B=Framework[CurrentSystem].Torsions[f1][i].B;
        C=Framework[CurrentSystem].Torsions[f1][i].C;
        D=Framework[CurrentSystem].Torsions[f1][i].D;
  
        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          posC=Framework[CurrentSystem].Atoms[f1][C].Position;
          posD=Framework[CurrentSystem].Atoms[f1][D].Position;
  
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
          r=MAX2((REAL)1.0e-8,sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z)));
          dr.x/=r; dr.y/=r; dr.z/=r;
  
          ds.x=Ddc.x-dot_cd*Dcb.x;
          ds.y=Ddc.y-dot_cd*Dcb.y;
          ds.z=Ddc.z-dot_cd*Dcb.z;
          s=MAX2((REAL)1.0e-8,sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z)));
          ds.x/=s; ds.y/=s; ds.z/=s;
  
          // compute Cos(Phi)
          // Phi is defined in protein convention Phi(trans)=Pi
          CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;
  
          // Ensure CosPhi is between -1 and 1.
          CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
          CosPhi2=SQR(CosPhi);
  
          parms=Framework[CurrentSystem].TorsionArguments[f1][i];
  
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
              break;
            case HARMONIC_COSINE_DIHEDRAL:
              // (1/2)*p_0*(cos(phi)-cos(p_1))^2
              // ===============================================
              // p_0/k_B [K]
              // p_1     [degrees]
              U=0.5*parms[0]*SQR(CosPhi-cos(parms[1]));
              break;
            case THREE_COSINE_DIHEDRAL:
              // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
              // ========================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
              break;
            case MM3_DIHEDRAL:
              // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
              // ========================================================================
              // p_0     [kcal/mol]
              // p_1     [kcal/mol]
              // p_2     [kcal/mol]
              U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
              break;
            case CFF_DIHEDRAL:
              // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
              // ======================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
              break;
            case CFF_DIHEDRAL2:
              // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
              // ======================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
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
              break;
            case TRAPPE_DIHEDRAL:
              // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
              // =============================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              // p_3/k_B [K]
              U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
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
              SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
              U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
              break;
            case OPLS_DIHEDRAL:
              // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
              // =================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              // p_3/k_B [K]
              U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
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
              break;
            default:
              printf("Undefined Torsion potential in routine 'CalculateFrameworkTorsionEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          // energy
          UHostTorsion+=U;
        }
      }
    }
  }
  return UHostTorsion;
}

REAL CalculateFrameworkImproperTorsionEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL rbc,U;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,sign;
  REAL CosPhi,Phi,CosPhi2,SinPhi,UHostImproperTorsion;
  VECTOR Pb,Pc;
  REAL *parms;

  UHostImproperTorsion=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfImproperTorsions[f1];i++)
      {
        A=Framework[CurrentSystem].ImproperTorsions[f1][i].A;
        B=Framework[CurrentSystem].ImproperTorsions[f1][i].B;
        C=Framework[CurrentSystem].ImproperTorsions[f1][i].C;
        D=Framework[CurrentSystem].ImproperTorsions[f1][i].D;
  
        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          posC=Framework[CurrentSystem].Atoms[f1][C].Position;
          posD=Framework[CurrentSystem].Atoms[f1][D].Position;
  
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
  
          parms=Framework[CurrentSystem].ImproperTorsionArguments[f1][i];
  
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
              break;
            case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
              // (1/2)*p_0*(cos(phi)-cos(p_1))^2
              // ===============================================
              // p_0/k_B [K]
              // p_1     [degrees]
              U=0.5*parms[0]*SQR(CosPhi-cos(parms[1]));
              break;
            case THREE_COSINE_IMPROPER_DIHEDRAL:
              // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
              // ========================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
              break;
            case MM3_IMPROPER_DIHEDRAL:
              // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
              // ========================================================================
              // p_0     [kcal/mol]
              // p_1     [kcal/mol]
              // p_2     [kcal/mol]
              U=0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2);
              break;
            case CFF_IMPROPER_DIHEDRAL:
              // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
              // ======================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2);
              break;
            case CFF_IMPROPER_DIHEDRAL2:
              // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
              // ======================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi));
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
              break;
            case TRAPPE_IMPROPER_DIHEDRAL:
              // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
              // =============================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              // p_3/k_B [K]
              U=parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi));
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
              SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
              U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]));
              break;
            case OPLS_IMPROPER_DIHEDRAL:
              // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
              // =================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              // p_3/k_B [K]
              U=0.5*(parms[0]+(1.0+CosPhi)*(parms[1]+parms[3]-2.0*(CosPhi-1.0)*(parms[2]-2.0*parms[3]*CosPhi)));
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
              break;
            case FIXED_IMPROPER_DIHEDRAL:
              U=0.0;
              break;
            default:
              printf("Undefined Improper-Torsion potential in routine 'CalculateFrameworkImproperTorsionEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          // energy
          UHostImproperTorsion+=U;
        }
      }
    }
  }
  return UHostImproperTorsion;
}

REAL CalculateFrameworkBondBondEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,f1;
  REAL *parms;
  REAL energy;
  REAL rab,rbc,UHostBondBond;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc;

  UHostBondBond=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondBonds[f1];i++)
      {
        A=Framework[CurrentSystem].BondBonds[f1][i].A;
        B=Framework[CurrentSystem].BondBonds[f1][i].B;
        C=Framework[CurrentSystem].BondBonds[f1][i].C;
  
        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id)))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          posC=Framework[CurrentSystem].Atoms[f1][C].Position;
  
          Rab.x=posA.x-posB.x;
          Rab.y=posA.y-posB.y;
          Rab.z=posA.z-posB.z;
          Rab=ApplyBoundaryCondition(Rab);
          rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
  
          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;
          Rbc=ApplyBoundaryCondition(Rbc);
          rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  
          parms=Framework[CurrentSystem].BondBondArguments[f1][i];
  
          switch(Framework[CurrentSystem].BondBondType[f1][i])
          {
            case CVFF_BOND_BOND_CROSS:
            case CFF_BOND_BOND_CROSS:
              // p_0*(rab-p_1)*(rbc-p_2)
              // =======================
              // p_0/k_B [K/A^2]
              // p_1     [A]
              // p_2     [A]
              energy=parms[0]*(rab-parms[1])*(rbc-parms[2]);
            break;
            default:
              printf("Undefined Bond-Bond potential in routine 'CalculateFrameworkBondBondEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          // add contribution to the energy
          UHostBondBond+=energy;
        }
      }
    }
  }
  return UHostBondBond;
}

REAL CalculateFrameworkBondBendEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,f1;
  REAL *parms,U;
  REAL cost,theta,sint;
  REAL rab,rbc,rac,UHostBondBend;
  VECTOR posA,posB,posC;
  VECTOR Rab,Rbc,Rac;

  UHostBondBend=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondBends[f1];i++)
      {
        A=Framework[CurrentSystem].BondBends[f1][i].A;
        B=Framework[CurrentSystem].BondBends[f1][i].B;
        C=Framework[CurrentSystem].BondBends[f1][i].C;
  
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id)))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          posC=Framework[CurrentSystem].Atoms[f1][C].Position;
  
          Rab.x=posA.x-posB.x;
          Rab.y=posA.y-posB.y;
          Rab.z=posA.z-posB.z;
          Rab=ApplyBoundaryCondition(Rab);
          rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
          Rab.x/=rab;
          Rab.y/=rab;
          Rab.z/=rab;
  
          Rbc.x=posC.x-posB.x;
          Rbc.y=posC.y-posB.y;
          Rbc.z=posC.z-posB.z;
          Rbc=ApplyBoundaryCondition(Rbc);
          rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
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
  
          cost=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
          cost=SIGN(MIN2(fabs(cost),(REAL)1.0),cost);
          theta=acos(cost);
          sint=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(cost)));
  
          parms=Framework[CurrentSystem].BondBendArguments[f1][i];
  
          switch(Framework[CurrentSystem].BondBendType[f1][i])
          {
            case CVFF_BOND_BEND_CROSS:
            case CFF_BOND_BEND_CROSS:
              // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
              // =========================================
              // p_0     [degrees]
              // p_1/k_B [K/A/rad]
              // p_2     [A]
              // p_3/k_B [K/A/rad]
              // p_4     [A]
              U=(theta-parms[0])*(parms[1]*(rab-parms[2])+parms[3]*(rbc-parms[4]));
              break;
            case MM3_BOND_BEND_CROSS:
              // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
              // =====================================
              // p_0     [mdyne/rad]
              // p_1     [A]
              // p_2     [A]
              // p_3     [degrees]
              U=parms[0]*((rab-parms[1])+(rbc-parms[2]))*RAD2DEG*(theta-parms[3]);
              break;
            case TRUNCATED_HARMONIC:
              // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
              // ================================================================
              // p_0/k_B [K/rad^2]
              // p_1     [degrees]
              // p_2     [A]
              U=0.5*parms[0]*SQR(theta-parms[1])*exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[2],8));
              break;
            case SCREENED_HARMONIC:
              // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
              // ===============================================
              // p_0/k_B [K/rad^2]
              // p_1     [degrees]
              // p_2     [A]
              // p_3     [A]
              U=0.5*parms[0]*SQR(theta-parms[1])*exp(-(rab/parms[2]+rbc/parms[3]));
              break;
            case SCREENED_VESSAL:
              // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
              // ============================================================================
              // p_0/k_B [K/rad^2]
              // p_1     [degrees]
              // p_2     [A]
              // p_3     [A]
              U=(parms[0]/(8.0*SQR(theta-M_PI)))*SQR(SQR(parms[1]-M_PI)-SQR(theta-M_PI))
                    *exp(-(rab/parms[2]+rbc/parms[3]));
              break;
            case TRUNCATED_VESSAL:
              // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
              //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
              // ============================================================================
              // p_0/k_B [K/rad^(4+p_2)]
              // p_1     [degrees]
              // p_2     [-]
              // p_3     [A]
              U=parms[0]*(pow(theta,parms[2])*SQR(theta-parms[1])*SQR(theta+parms[1]-2.0*M_PI)
                    -0.5*parms[2]*pow((REAL)M_PI,(parms[2]-1.0))*SQR(theta-parms[1])*pow(M_PI-parms[1],(REAL)3.0))
                    *exp(-(pow(rab,8)+pow(rbc,8))/pow(parms[3],8));
              break;
            default:
              printf("Undefined Bond-Bend potential in routine 'CalculateFrameworkBondBendEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          // energy
          UHostBondBend+=U;
        }
      }
    }
  }
  return UHostBondBend;
}

// Bend/Bend cross term for a centered atom B
// the first angle is A-B-C
// the second angle is A-B-D
REAL CalculateFrameworkBendBendEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL U,rab,rbc,rbd;
  VECTOR Dab,Dbc,Dbd;
  REAL dot_abc,dot_abd,UHostBendBend;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL *parms;

  UHostBendBend=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBendBends[f1];i++)
      {
        A=Framework[CurrentSystem].BendBends[f1][i].A;
        B=Framework[CurrentSystem].BendBends[f1][i].B;
        C=Framework[CurrentSystem].BendBends[f1][i].C;
        D=Framework[CurrentSystem].BendBends[f1][i].D;
  
        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          posC=Framework[CurrentSystem].Atoms[f1][C].Position;
          posD=Framework[CurrentSystem].Atoms[f1][D].Position;
  
          Dab.x=posA.x-posB.x;
          Dab.y=posA.y-posB.y;
          Dab.z=posA.z-posB.z;
          Dab=ApplyBoundaryCondition(Dab);
          rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));
          Dab.x/=rab; Dab.y/=rab; Dab.z/=rab;
  
          Dbc.x=posC.x-posB.x;
          Dbc.y=posC.y-posB.y;
          Dbc.z=posC.z-posB.z;
          Dbc=ApplyBoundaryCondition(Dbc);
          rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
          Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;
  
          Dbd.x=posD.x-posB.x;
          Dbd.y=posD.y-posB.y;
          Dbd.z=posD.z-posB.z;
          Dbd=ApplyBoundaryCondition(Dbd);
          rbd=sqrt(SQR(Dbd.x)+SQR(Dbd.y)+SQR(Dbd.z));
          Dbd.x/=rbd; Dbd.y/=rbd; Dbd.z/=rbd;
  
          dot_abc=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
          CosTheta1=dot_abc;
          CosTheta1=SIGN(MIN2(fabs(CosTheta1),(REAL)1.0),CosTheta1);
          Theta1=acos(CosTheta1);
          SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));
  
          dot_abd=Dab.x*Dbd.x+Dab.y*Dbd.y+Dab.z*Dbd.z;
          CosTheta2=dot_abd;
          CosTheta2=SIGN(MIN2(fabs(CosTheta2),(REAL)1.0),CosTheta2);
          Theta2=acos(CosTheta2);
          SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));
  
          parms=(REAL*)&Framework[CurrentSystem].BendBendArguments[f1][i];
  
          switch(Framework[CurrentSystem].BendBendType[f1][i])
          {
            case CVFF_BEND_BEND_CROSS:
            case CFF_BEND_BEND_CROSS:
              // p_0*(Theta1-p_1)*(Theta2-p_2)
              // ===================================
              // p_0/k_B [K/rad^2)]
              // p_1     [degrees]
              // p_2     [degrees]
              U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2]);
              break;
            case MM3_BEND_BEND_CROSS:
              // -p_0*(Theta1-p_1)*(Theta2-p_2)
              // ===================================
              // p_0     [mdyne A/rad^2]
              // p_1     [degrees]
              // p_2     [degrees]
              U=-parms[0]*SQR(RAD2DEG)*(Theta1-parms[1])*(Theta2-parms[2]);
              break;
            default:
              printf("Undefined Bend-Bend potential in routine 'CalculateFrameworkBendBendEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
          // energy
          UHostBendBend+=U;
        }
      }
    }
  }
  return UHostBendBend;
}

REAL CalculateFrameworkBondTorsionEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL U,rab,rbc,rcd;
  VECTOR Dab,Dcb,Ddc,dr,ds;
  REAL dot_ab,dot_cd,r,s,temp;
  REAL CosPhi,CosPhi2;
  REAL *parms,UHostBondTorsion;

  UHostBondTorsion=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondTorsions[f1];i++)
      {
        A=Framework[CurrentSystem].BondTorsions[f1][i].A;
        B=Framework[CurrentSystem].BondTorsions[f1][i].B;
        C=Framework[CurrentSystem].BondTorsions[f1][i].C;
        D=Framework[CurrentSystem].BondTorsions[f1][i].D;
  
        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          posC=Framework[CurrentSystem].Atoms[f1][C].Position;
          posD=Framework[CurrentSystem].Atoms[f1][D].Position;
  
          Dab.x=posA.x-posB.x;
          Dab.y=posA.y-posB.y;
          Dab.z=posA.z-posB.z;
          Dab=ApplyBoundaryCondition(Dab);
          rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));
  
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
          rcd=sqrt(SQR(Ddc.x)+SQR(Ddc.y)+SQR(Ddc.z));
  
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
  
          parms=(REAL*)&Framework[CurrentSystem].BondTorsionArguments[f1][i];
  
          switch(Framework[CurrentSystem].BondTorsionType[f1][i])
          {
            case MM3_BOND_TORSION_CROSS:
              // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
              // =====================================================================================
              // p_0     [kcal/A mole]
              // p_1     [kcal/A mole]
              // p_2     [kcal/A mole]
              // p_3     [A]
              temp=(rbc-parms[3]);
              U=parms[0]*temp*CosPhi+parms[1]*temp*(2.0*CosPhi2-1.0)+parms[2]*temp*(4.0*CosPhi2*CosPhi-3.0*CosPhi);
              break;
            default:
              printf("Undefined Bond-Torsion potential in routine 'CalculateFrameworkBondTorsionEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          // energy
          UHostBondTorsion+=U;
        }
      }
    }
  }
  return UHostBondTorsion;
}

REAL CalculateFrameworkBendTorsionEnergy(int flag,int f2,int atom_id)
{
  int i,A,B,C,D,f1;
  POINT posA,posB,posC,posD;
  REAL U,rab,rbc,rcd;
  VECTOR Dab,Dbc,Dcd,dr,ds,Pb,Pc;
  REAL dot_ab,dot_cd,r,s;
  REAL CosPhi,CosPhi2;
  REAL CosTheta1,CosTheta2,Theta1,Theta2,SinTheta1,SinTheta2;
  REAL sign,Phi,SinPhi;
  REAL *parms,UHostBendTorsion;

  UHostBendTorsion=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBendTorsions[f1];i++)
      {
        A=Framework[CurrentSystem].BendTorsions[f1][i].A;
        B=Framework[CurrentSystem].BendTorsions[f1][i].B;
        C=Framework[CurrentSystem].BendTorsions[f1][i].C;
        D=Framework[CurrentSystem].BendTorsions[f1][i].D;
  
        // if flag is 'false', only execute when atom 'atom_id' is involved
        // else compute total energy
        if(flag||((f1==f2)&&(A==atom_id||B==atom_id||C==atom_id||D==atom_id)))
        {
          posA=Framework[CurrentSystem].Atoms[f1][A].Position;
          posB=Framework[CurrentSystem].Atoms[f1][B].Position;
          posC=Framework[CurrentSystem].Atoms[f1][C].Position;
          posD=Framework[CurrentSystem].Atoms[f1][D].Position;
  
          Dab.x=posA.x-posB.x;
          Dab.y=posA.y-posB.y;
          Dab.z=posA.z-posB.z;
          Dab=ApplyBoundaryCondition(Dab);
          rab=sqrt(SQR(Dab.x)+SQR(Dab.y)+SQR(Dab.z));
  
          Dbc.x=posC.x-posB.x;
          Dbc.y=posC.y-posB.y;
          Dbc.z=posC.z-posB.z;
          Dbc=ApplyBoundaryCondition(Dbc);
          rbc=sqrt(SQR(Dbc.x)+SQR(Dbc.y)+SQR(Dbc.z));
          Dbc.x/=rbc; Dbc.y/=rbc; Dbc.z/=rbc;
  
          Dcd.x=posD.x-posC.x;
          Dcd.y=posD.y-posC.y;
          Dcd.z=posD.z-posC.z;
          Dcd=ApplyBoundaryCondition(Dcd);
          rcd=sqrt(SQR(Dcd.x)+SQR(Dcd.y)+SQR(Dcd.z));
  
          dot_ab=Dab.x*Dbc.x+Dab.y*Dbc.y+Dab.z*Dbc.z;
          CosTheta1=dot_ab/rab;
          CosTheta1=MAX2(MIN2(CosTheta1,(REAL)1.0),-1.0);
          Theta1=acos(CosTheta1);
          SinTheta1=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta1)));
  
          dot_cd=Dcd.x*Dbc.x+Dcd.y*Dbc.y+Dcd.z*Dbc.z;
          CosTheta2=-dot_cd/rcd;
          CosTheta2=MAX2(MIN2(CosTheta2,(REAL)1.0),-1.0);
          Theta2=acos(CosTheta2);
          SinTheta2=MAX2((REAL)1.0e-8,sqrt(1.0-SQR(CosTheta2)));
  
          dr.x=Dab.x-dot_ab*Dbc.x;
          dr.y=Dab.y-dot_ab*Dbc.y;
          dr.z=Dab.z-dot_ab*Dbc.z;
          r=MAX2((REAL)1.0e-8,sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z)));
          dr.x/=r; dr.y/=r; dr.z/=r;
  
          ds.x=Dcd.x-dot_cd*Dbc.x;
          ds.y=Dcd.y-dot_cd*Dbc.y;
          ds.z=Dcd.z-dot_cd*Dbc.z;
          s=MAX2((REAL)1.0e-8,sqrt(SQR(ds.x)+SQR(ds.y)+SQR(ds.z)));
          ds.x/=s; ds.y/=s; ds.z/=s;
  
          // compute Cos(Phi)
          // Phi is defined in protein convention Phi(trans)=Pi
          CosPhi=dr.x*ds.x+dr.y*ds.y+dr.z*ds.z;
  
          // Ensure CosPhi is between -1 and 1.
          CosPhi=SIGN(MIN2(fabs(CosPhi),(REAL)1.0),CosPhi);
          CosPhi2=SQR(CosPhi);
  
          parms=(REAL*)&Framework[CurrentSystem].BendTorsionArguments[f1][i];
  
          switch(Framework[CurrentSystem].BendTorsionType[f1][i])
          {
            case CVFF_BEND_TORSION_CROSS:
            case CFF_BEND_TORSION_CROSS:
              // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
              // =====================================================================================
              // p_0/k_B [K/rad^3]
              // p_1     [degrees]
              // p_2     [degrees]
              U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi;
              break;
            case SMOOTHED_DIHEDRAL:
              // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2)]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K/rad^2]
              // p_1     [-]
              // p_2     [degrees]
              Pb.x=Dab.z*Dbc.y-Dab.y*Dbc.z;
              Pb.y=Dab.x*Dbc.z-Dab.z*Dbc.x;
              Pb.z=Dab.y*Dbc.x-Dab.x*Dbc.y;
              Pc.x=Dbc.y*Dcd.z-Dbc.z*Dcd.y;
              Pc.y=Dbc.z*Dcd.x-Dbc.x*Dcd.z;
              Pc.z=Dbc.x*Dcd.y-Dbc.y*Dcd.x;
              sign=(Dbc.x*(Pc.z*Pb.y-Pc.y*Pb.z)+Dbc.y*(Pb.z*Pc.x-Pb.x*Pc.z)
                    +Dbc.z*(Pc.y*Pb.x-Pc.x*Pb.y));
              Phi=SIGN(acos(CosPhi),sign);
              SinPhi=sin(Phi);
              SinPhi=SIGN(MAX2((REAL)1.0e-8,fabs(SinPhi)),SinPhi);  // remove singularity
              U=parms[0]*(1.0+cos(parms[1]*Phi-parms[2]))*Smoothing(Theta1)*Smoothing(Theta2);
              break;
            case SMOOTHED_THREE_COSINE_DIHEDRAL:
              // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                     Smoothing(Theta1)*Smoothing(Theta2);
              break;
            case NICHOLAS_DIHEDRAL:
              // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=(0.5*parms[0]*(1.0+CosPhi)+parms[1]*(1.0-CosPhi2)+0.5*parms[2]*(1.0-3.0*CosPhi+4.0*CosPhi*CosPhi2))*
                     Smoothing(Theta1);
              break;
            case SMOOTHED_CFF_DIHEDRAL:
              // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=(parms[0]*(1.0-CosPhi)+2.0*parms[1]*(1.0-CosPhi2)+parms[2]*(1.0+3.0*CosPhi-4.0*CosPhi*CosPhi2))*Smoothing(Theta1)*Smoothing(Theta2);
              break;
            case SMOOTHED_CFF_DIHEDRAL2:
              // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K]
              // p_1/k_B [K]
              // p_2/k_B [K]
              U=(parms[0]*(1.0+CosPhi)+parms[2]+CosPhi*(-3.0*parms[2]+2.0*CosPhi*(parms[1]+2.0*parms[2]*CosPhi)))*Smoothing(Theta1)*Smoothing(Theta2);
              break;
            case SMOOTHED_CFF_BEND_TORSION_CROSS:
              // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
              // ======================================================================================
              // p_0/k_B [K/rad^3]
              // p_1     [degrees]
              // p_2     [degrees]
              U=parms[0]*(Theta1-parms[1])*(Theta2-parms[2])*CosPhi*Smoothing(Theta1)*Smoothing(Theta2);
              break;
            default:
              printf("Undefined Bend-Torsion potential in routine 'CalculateFrameworkBendTorsionEnergy' ('framework_energy.c')\n");
              exit(0);
              break;
          }
  
          // energy
          UHostBendTorsion+=U;
        }
      }
    }
  }
  return UHostBendTorsion;
}

int CalculateFrameworkIntraVDWEnergy(void)
{
  int i,j,typeA,typeB;
  REAL rr;
  VECTOR posA,posB,dr;
  int f1,f2;

  // Framework-Framework energy
  UHostHostVDW[CurrentSystem]=0.0;

  if(!InternalFrameworkLennardJonesInteractions) return 0;

  // contributions from intra-framework
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

        for(j=i+1;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],0))
          {
            typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
            posB=Framework[CurrentSystem].Atoms[f1][j].AnisotropicPosition;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
              UHostHostVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);
          }
        }
      }
    }
  }

  // contributions from interactions between the frameworks
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1+1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
        {
          typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
          posB=Framework[CurrentSystem].Atoms[f2][j].AnisotropicPosition;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffVDWSquared)
            UHostHostVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraChargeChargeEnergy(void)
{
  int i,j,typeA,typeB;
  REAL chargeA,chargeB;
  REAL r,rr,energy;
  VECTOR posA,posB,dr;
  REAL SwitchingValue,TranslationValue;
  int f1,f2;

  // Framework-Framework energy
  UHostHostChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        chargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;

        for(j=i+1;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],1))
          {
            typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
            posB=Framework[CurrentSystem].Atoms[f1][j].Position;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

              switch(ChargeMethod)
              {
                case NONE:
                  break;
                case TRUNCATED_COULOMB:
                  UHostHostChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  break;
                case SHIFTED_COULOMB:
                  UHostHostChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(rr>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                   SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                      SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                      SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                    energy=energy*SwitchingValue+TranslationValue;
                  }
                  UHostHostChargeChargeReal[CurrentSystem]+=energy;
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  UHostHostChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                            (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                            (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  break;
                case EWALD:
                  UHostHostChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                   (erfc(Alpha[CurrentSystem]*r)/r);
                  break;
                default:
                  printf("Unknown charge-charge method in 'CalculateFrameworkIntraChargeChargeEnergy'\n");
                  exit(0);
                  break;
              }
            }
          }
        }
      }
    }
  }


  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1+1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        chargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;

        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
        {
          typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
          posB=Framework[CurrentSystem].Atoms[f2][j].Position;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared)
          {
            r=sqrt(rr);
            chargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;

            switch(ChargeMethod)
            {
              case NONE:
                break;
              case TRUNCATED_COULOMB:
                UHostHostChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                break;
              case SHIFTED_COULOMB:
                UHostHostChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                break;
              case SMOOTHED_COULOMB:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                if(rr>CutOffChargeChargeSwitchSquared)
                {
                  SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                 SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                  TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                   (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                    SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                    SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                  energy=energy*SwitchingValue+TranslationValue;
                }
                UHostHostChargeChargeReal[CurrentSystem]+=energy;
                break;
              case WOLFS_METHOD_DAMPED_FG:
                UHostHostChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                         -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                          (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                          (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                break;
              case EWALD:
                UHostHostChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                 (erfc(Alpha[CurrentSystem]*r)/r);
                break;
              default:
                printf("Unknown charge-charge method in 'CalculateFrameworkIntraChargeChargeEnergy'\n");
                exit(0);
                break;
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraChargeBondDipoleEnergy(void)
{
  int i,j,f1,f2;
  int A1,A2;
  int Type;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL r,rr,ri2,cosA,energy,temp,length,chargeB;
  VECTOR dipoleA;
  REAL Bt1;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;
  REAL SwitchingValue;

  Bt1=0.0;
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  UHostHostChargeBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
      {
        DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
        A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
        A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
        posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
        posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        dipoleA=ApplyBoundaryCondition(dipoleA);
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
        length=sqrt(ri2);
        temp=DipoleMagnitudeA/length;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          // if framework are different, or if they are the same but not excluded within the framework
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][j][i],2))
          {
            Type=Framework[CurrentSystem].Atoms[f1][j].Type;
            if(PseudoAtoms[Type].HasCharges)
            {
              posB=Framework[CurrentSystem].Atoms[f1][j].Position;
              chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeBondDipoleSquared)
              {
                r=sqrt(rr);
                SwitchingValue=1.0;
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt1=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                  case SHIFTED_COULOMB:
                    Bt1=1.0/(r*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt1=1.0/(r*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                    break;
                  case EWALD:
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    break;
                  default:
                    printf("Unknown charge-bonddipole method in 'CalculateFrameworkIntraChargeBondDipoleEnergy'\n");
                    exit(0);
                    break;
                }

                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
                UHostHostChargeBondDipoleReal[CurrentSystem]-=energy;
              }
            }
          }
        }
      }
    }
  }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=0;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      if(f1!=f2)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
        {
          DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
          A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
          A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
          posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
          posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
          dipoleA.x=posA2.x-posA1.x;
          dipoleA.y=posA2.y-posA1.y;
          dipoleA.z=posA2.z-posA1.z;
          dipoleA=ApplyBoundaryCondition(dipoleA);
          posA.x=posA1.x+0.5*dipoleA.x;
          posA.y=posA1.y+0.5*dipoleA.y;
          posA.z=posA1.z+0.5*dipoleA.z;
          ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
          length=sqrt(ri2);
          temp=DipoleMagnitudeA/length;
          dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

          for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
          {
            Type=Framework[CurrentSystem].Atoms[f2][j].Type;
            if(PseudoAtoms[Type].HasCharges)
            {
              posB=Framework[CurrentSystem].Atoms[f2][j].Position;
              chargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeBondDipoleSquared)
              {
                r=sqrt(rr);
                SwitchingValue=1.0;
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt1=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                  case SHIFTED_COULOMB:
                    Bt1=1.0/(r*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt1=1.0/(r*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                    break;
                  case EWALD:
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    break;
                  default:
                    printf("Unknown charge-bonddipole method in 'CalculateFrameworkIntraChargeBondDipoleEnergy'\n");
                    exit(0);
                    break;
                }

                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
                UHostHostChargeBondDipoleReal[CurrentSystem]-=energy;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraBondDipoleBondDipoleEnergy(void)
{
  int i,j,f1,f2;
  int A1,A2,B1,B2;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,rr,ri2,rk2,cosAB,cosA,cosB,energy,temp,length;
  VECTOR dipoleA,dipoleB;
  REAL Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL SwitchingValue;
  REAL_MATRIX3x3 v;

  Bt1=Bt2=0.0;  
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  UHostHostBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
      {
        DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
        A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
        A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
        posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
        posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        dipoleA=ApplyBoundaryCondition(dipoleA);
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
        length=sqrt(ri2);
        temp=DipoleMagnitudeA/length;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        for(j=i+1;j<Framework[CurrentSystem].NumberOfBondDipoles[f1];j++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][j];
          B1=Framework[CurrentSystem].BondDipoles[f1][j].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][j].B;

          // if framework are different, or if they are the same but not excluded within the framework
          if(!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],3))
          {
            posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
            posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

            dipoleB.x=posB2.x-posB1.x;
            dipoleB.y=posB2.y-posB1.y;
            dipoleB.z=posB2.z-posB1.z;
            dipoleB=ApplyBoundaryCondition(dipoleB);
            posB.x=posB1.x+0.5*dipoleB.x;
            posB.y=posB1.y+0.5*dipoleB.y;
            posB.z=posB1.z+0.5*dipoleB.z;
            rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
            length=sqrt(rk2);
            temp=DipoleMagnitudeB/length;
            dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt1=Bt2=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt1=1.0/(rr*r);
                  Bt2=3.0/(SQR(rr)*r);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                  }
                  break;
                case EWALD:
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkIntraBondDipoleBondDipoleEnergy'\n");
                  exit(0);
                  break;
              }

              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              UHostHostBondDipoleBondDipoleReal[CurrentSystem]+=energy;
            }
          }
        }
      }
    }
  }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1+1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
      {
        DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
        A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
        A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
        posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
        posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        dipoleA=ApplyBoundaryCondition(dipoleA);
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
        length=sqrt(ri2);
        temp=DipoleMagnitudeA/length;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        for(j=0;j<Framework[CurrentSystem].NumberOfBondDipoles[f2];j++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f2][j];
          B1=Framework[CurrentSystem].BondDipoles[f2][j].A;
          B2=Framework[CurrentSystem].BondDipoles[f2][j].B;

          // if framework are different, or if they are the same but not excluded within the framework
          posB1=Framework[CurrentSystem].Atoms[f2][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f2][B2].Position;

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt1=Bt2=0.0;  
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt1=1.0/(rr*r);
                Bt2=3.0/(SQR(rr)*r);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkIntraBondDipoleBondDipoleEnergy'\n");
                exit(0);
                break;
            }

            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            UHostHostBondDipoleBondDipoleReal[CurrentSystem]+=energy;
          }
        }
      }
    }
  }
  return 0;
}

void CalculateFrameworkAdsorbateVDWEnergy(void)
{
  int i,j,k,l,f1;
  int typeA,typeB,type;
  REAL rr;
  VECTOR posA,posB,dr;
  VECTOR s;
  int icell0,icell;

  UHostAdsorbateVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
      {
        UHostAdsorbateVDW[CurrentSystem]+=InterpolateVDWGrid(typeA,posA);
      }
      else 
      {
        if(UseCellLists[CurrentSystem])
        {
          // convert from xyz to abc
          s.x=InverseBox[CurrentSystem].ax*posA.x+InverseBox[CurrentSystem].bx*posA.y+InverseBox[CurrentSystem].cx*posA.z;
          s.y=InverseBox[CurrentSystem].ay*posA.x+InverseBox[CurrentSystem].by*posA.y+InverseBox[CurrentSystem].cy*posA.z;
          s.z=InverseBox[CurrentSystem].az*posA.x+InverseBox[CurrentSystem].bz*posA.y+InverseBox[CurrentSystem].cz*posA.z;

          // apply boundary condition
          s.x-=(REAL)NINT(s.x);
          s.y-=(REAL)NINT(s.y);
          s.z-=(REAL)NINT(s.z);

          // s between 0 and 1
          s.x+=0.5;
          s.y+=0.5;
          s.z+=0.5;

          // compute the corresponding cell-id
          icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
                 ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
                 ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            // loop over cells
            for(l=0;l<27;l++)
            {
              icell=CellListMap[CurrentSystem][icell0][l];

              k=Framework[CurrentSystem].CellListHead[f1][icell];

              while(k>=0)
              {
                posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
                typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                  UHostAdsorbateVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);

                k=Framework[CurrentSystem].CellList[f1][k];
              }
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
                UHostAdsorbateVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);
            }
          }
        }
      }
    }
  }
}

void CalculateFrameworkCationVDWEnergy(void)
{
  int i,j,k,l,f1;
  int typeA,typeB,type;
  REAL rr;
  VECTOR posA,posB,dr;
  VECTOR s;
  int icell0,icell;

  UHostCationVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
      {
        UHostCationVDW[CurrentSystem]+=InterpolateVDWGrid(typeA,posA);
      }
      else
      {
        if(UseCellLists[CurrentSystem])
        {
          // convert from xyz to abc
          s.x=InverseBox[CurrentSystem].ax*posA.x+InverseBox[CurrentSystem].bx*posA.y+InverseBox[CurrentSystem].cx*posA.z;
          s.y=InverseBox[CurrentSystem].ay*posA.x+InverseBox[CurrentSystem].by*posA.y+InverseBox[CurrentSystem].cy*posA.z;
          s.z=InverseBox[CurrentSystem].az*posA.x+InverseBox[CurrentSystem].bz*posA.y+InverseBox[CurrentSystem].cz*posA.z;

          // apply boundary condition
          s.x-=(REAL)NINT(s.x);
          s.y-=(REAL)NINT(s.y);
          s.z-=(REAL)NINT(s.z);

          // s between 0 and 1
          s.x+=0.5;
          s.y+=0.5;
          s.z+=0.5;

          // compute the corresponding cell-id
          icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
                 ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
                 ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            // loop over cells
            for(l=0;l<27;l++)
            {
              icell=CellListMap[CurrentSystem][icell0][l];

              k=Framework[CurrentSystem].CellListHead[f1][icell];

              while(k>=0)
              {
                posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
                typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                  UHostCationVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);

                k=Framework[CurrentSystem].CellList[f1][k];
              }
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
            {
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
                UHostCationVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);
            }
          }
        }
      }
    }
  }
}

int CalculateFrameworkAdsorbateChargeChargeEnergy(void)
{
  int i,j,k;
  int typeA,typeB;
  REAL r,rr,energy;
  REAL chargeA,chargeB;
  VECTOR posA,posB,dr;
  REAL SwitchingValue,TranslationValue;

  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
      {
        UHostAdsorbateChargeChargeReal[CurrentSystem]+=InterpolateCoulombGrid(typeA,posA);
      }
      else
      {
        for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
          {
            posB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;
            typeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              chargeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Charge;

              switch(ChargeMethod)
              {
                case NONE:
                  break;
                case TRUNCATED_COULOMB:
                  UHostAdsorbateChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  break;
                case SHIFTED_COULOMB:
                  UHostAdsorbateChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(rr>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                   SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                      SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                      SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                    energy=energy*SwitchingValue+TranslationValue;
                  }
                  UHostAdsorbateChargeChargeReal[CurrentSystem]+=energy;
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  UHostAdsorbateChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                            (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                            (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  break;
                case EWALD:
                  UHostAdsorbateChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (erfc(Alpha[CurrentSystem]*r)/r);
                  break;
                default:
                  printf("Unknown charge-charge method in 'CalculateFrameworkAdsorbateChargeChargeEnergy'\n");
                  exit(0);
                  break;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationChargeChargeEnergy(void)
{
  int i,j,k;
  int typeA,typeB;
  REAL r,rr,energy;
  REAL chargeA,chargeB;
  VECTOR posA,posB,dr;
  REAL SwitchingValue,TranslationValue;

  UHostCationChargeChargeReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;
      posA=Cations[CurrentSystem][i].Atoms[j].Position;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
      {
        UHostCationChargeChargeReal[CurrentSystem]+=InterpolateCoulombGrid(typeA,posA);
      }
      else
      {
        for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
          {
            posB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;
            typeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              chargeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Charge;

              switch(ChargeMethod)
              {
                case NONE:
                  break;
                case TRUNCATED_COULOMB:
                  UHostCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  break;
                case SHIFTED_COULOMB:
                  UHostCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(rr>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                   SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                      SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                      SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                    energy=energy*SwitchingValue+TranslationValue;
                  }
                  UHostCationChargeChargeReal[CurrentSystem]+=energy;
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  UHostCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                            (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                            (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  break;
                case EWALD:
                  UHostCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (erfc(Alpha[CurrentSystem]*r)/r);
                  break;
                default:
                  printf("Unknown charge-charge method in 'CalculateFrameworkCationChargeChargeEnergy'\n");
                  exit(0);
                  break;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateChargeBondDipoleEnergy(void)
{
  int i,j,k,f1;
  int A1,A2;
  int Type,TypeA;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL r,rr,ri2,cosA,energy,temp,length,chargeB;
  VECTOR dipoleA;
  REAL Bt0,Bt1;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;
  REAL SwitchingValue;

  UHostAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  Bt0=Bt1=0.0;  
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
      A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
      A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
      posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
      {
        for(k=0;k<Adsorbates[CurrentSystem][j].NumberOfAtoms;k++)
        {
          Type=Adsorbates[CurrentSystem][j].Atoms[k].Type;
          if(PseudoAtoms[Type].HasCharges)
          {
            posB=Adsorbates[CurrentSystem][j].Atoms[k].Position;
            chargeB=Adsorbates[CurrentSystem][j].Atoms[k].Charge;
  
            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkAdsorbateChargeBondDipoleEnergy'\n");
                  exit(0);
                  break;
              }

              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
              UHostAdsorbateChargeBondDipoleReal[CurrentSystem]-=energy;
            }
          }
        }
      }
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeA=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      posA1=Adsorbates[CurrentSystem][i].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][i].Atoms[A2].Position;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
        {
          posB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;
          Type=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;
          chargeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkAdsorbateChargeBondDipoleEnergy'\n");
                exit(0);
                break;
            }

            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
            UHostAdsorbateChargeBondDipoleReal[CurrentSystem]-=energy;
          }
        }
      }
    }
  }

  return 0;
}

int CalculateFrameworkCationChargeBondDipoleEnergy(void)
{
  int i,j,k,f1;
  int A1,A2;
  int Type,TypeA;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL r,rr,ri2,cosA,energy,temp,length,chargeB;
  VECTOR dipoleA;
  REAL Bt0,Bt1;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;
  REAL SwitchingValue;

  UHostCationChargeBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  Bt0=Bt1=0.0;  
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
      A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
      A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
      posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
      {
        for(k=0;k<Cations[CurrentSystem][j].NumberOfAtoms;k++)
        {
          Type=Cations[CurrentSystem][j].Atoms[k].Type;
          if(PseudoAtoms[Type].HasCharges)
          {
            posB=Cations[CurrentSystem][j].Atoms[k].Position;
            chargeB=Cations[CurrentSystem][j].Atoms[k].Charge;
  
            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkCationChargeBondDipoleEnergy'\n");
                  exit(0);
                  break;
              }

              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
              UHostCationChargeBondDipoleReal[CurrentSystem]-=energy;
            }
          }
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    TypeA=Cations[CurrentSystem][i].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      posA1=Cations[CurrentSystem][i].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][i].Atoms[A2].Position;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
        {
          posB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;
          Type=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;
          chargeB=Framework[CurrentSystem].Atoms[CurrentFramework][k].Charge;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkCationChargeBondDipoleEnergy'\n");
                exit(0);
                break;
            }

            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
            UHostCationChargeBondDipoleReal[CurrentSystem]-=energy;
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergy(void)
{
  int i,k,l,f1;
  int A1,A2,B1,B2;
  int TypeB;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,r2,r3,r4,r5,ri2,rk2,cosA,cosB,cosAB,energy,temp,length;
  VECTOR dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL SwitchingValue;

  UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  Bt0=Bt1=Bt2=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
      A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
      A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
      posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        TypeB=Adsorbates[CurrentSystem][k].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
          posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(r2<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(r2);
            r3=r2*r;
            r4=r2*r2;
            r5=r2*r3;
            SwitchingValue=1.0;

            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/r;
                Bt1=1.0/r3;
                Bt2=3.0/r5;
                break;
              case SMOOTHED_COULOMB:
                Bt0=1.0/r;
                Bt1=1.0/r3;
                Bt2=3.0/r5;
                if(r2>CutOffBondDipoleBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*r5+SwitchingBondDipoleBondDipoleFactors5[4]*r4+SwitchingBondDipoleBondDipoleFactors5[3]*r3+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*r2+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*r2)+
                    erfc(Alpha[CurrentSystem]*r)/r3;
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*r4)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*r2)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/r5;
                break;
              default:
                printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergy'\n");
                exit(0);
                break;
            }

            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=energy;
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationBondDipoleBondDipoleEnergy(void)
{
  int i,k,l,f1;
  int A1,A2,B1,B2;
  int TypeB;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,r2,r3,r4,r5,ri2,rk2,cosA,cosB,cosAB,energy,temp,length;
  VECTOR dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL SwitchingValue;

  UHostCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  if(ChargeMethod==NONE) return 0;

  Bt0=Bt1=Bt2=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
      A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
      A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
      posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        TypeB=Cations[CurrentSystem][k].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
          posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(r2<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(r2);
            r3=r2*r;
            r4=r2*r2;
            r5=r2*r3;
            SwitchingValue=1.0;

            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/r;
                Bt1=1.0/r3;
                Bt2=3.0/r5;
                break;
              case SMOOTHED_COULOMB:
                Bt0=1.0/r;
                Bt1=1.0/r3;
                Bt2=3.0/r5;
                if(r2>CutOffBondDipoleBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*r5+SwitchingBondDipoleBondDipoleFactors5[4]*r4+SwitchingBondDipoleBondDipoleFactors5[3]*r3+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*r2+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*r2)+
                    erfc(Alpha[CurrentSystem]*r)/r3;
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*r4)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*r2)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/r5;
                break;
              default:
                printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkCationBondDipoleBondDipoleEnergy'\n");
                exit(0);
                break;
            }

            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            UHostCationBondDipoleBondDipoleReal[CurrentSystem]+=energy;
          }
        }
      }
    }
  }
  return 0;
}

void CalculateFrameworkShiftEnergyDifferenceAdsorbateVDW(void)
{
  int i,j,m,typeA,typeB;
  POINT posB,posA_new,posA_old;
  REAL rr;
  VECTOR dr;

  // Framework-Adsorbate energy
  UAdsorbateVDWDelta[CurrentSystem]=0.0;
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][i].AnisotropicPosition;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferenceAnisotropicPosition;

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Adsorbates[CurrentSystem][m].Atoms[j].Type;
        posB=Adsorbates[CurrentSystem][m].Atoms[j].AnisotropicPosition;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<CutOffVDWSquared)
          UAdsorbateVDWDelta[CurrentSystem]+=PotentialValue(typeA,typeB,rr);

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<CutOffVDWSquared)
          UAdsorbateVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);

      }
    }
  }
}

void CalculateFrameworkShiftEnergyDifferenceAdsorbateCharge(void)
{
  int i,j,m,typeA,typeB;
  int A1,A2,B1,B2,TypeMolB;
  POINT posB,posA_new,posA_old;
  REAL rr,r,chargeA,chargeB,Bt0,Bt1,Bt2,cosA,cosB,cosAB,energy;
  VECTOR posA1,posA2,dipoleA_new,dipoleA_old;
  VECTOR posB1,posB2,dipoleB;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,temp;
  REAL SwitchingValue,TranslationValue;
  VECTOR dr;

  // Framework-Adsorbate energy
  UAdsorbateChargeChargeRealDelta[CurrentSystem]=0.0;
  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;


  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferencePosition;
    chargeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge;

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Adsorbates[CurrentSystem][m].Atoms[j].Type;
        posB=Adsorbates[CurrentSystem][m].Atoms[j].Position;
        chargeB=Adsorbates[CurrentSystem][m].Atoms[j].Charge;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeChargeSquared)
        {
          r=sqrt(rr);
          switch(ChargeMethod)
          {
            case NONE:
              break;
            case TRUNCATED_COULOMB:
              UAdsorbateChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
              break;
            case SHIFTED_COULOMB:
              UAdsorbateChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
              break;
            case SMOOTHED_COULOMB:
              energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
              if(rr>CutOffChargeChargeSwitchSquared)
              {
                SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                               SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                 (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                  SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                  SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                energy=energy*SwitchingValue+TranslationValue;
              }
              UAdsorbateChargeChargeRealDelta[CurrentSystem]+=energy;
              break;
            case WOLFS_METHOD_DAMPED_FG:
              UAdsorbateChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                        -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                         (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                        (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
               break;
            case EWALD:
              UAdsorbateChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (erfc(Alpha[CurrentSystem]*r)/r);
              break;
            default:
              printf("Unknown charge-charge method in 'CalculateFrameworkShiftEnergyDifferenceAdsorbate'\n");
              exit(0);
              break;
          }
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeChargeSquared)
        {
          r=sqrt(rr);
          switch(ChargeMethod)
          {
            case NONE:
              break;
            case TRUNCATED_COULOMB:
              UAdsorbateChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
              break;
            case SHIFTED_COULOMB:
              UAdsorbateChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
              break;
            case SMOOTHED_COULOMB:
              energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
              if(rr>CutOffChargeChargeSwitchSquared)
              {
                SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                               SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                 (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                  SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                  SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                energy=energy*SwitchingValue+TranslationValue;
              }
              UAdsorbateChargeChargeRealDelta[CurrentSystem]-=energy;
              break;
            case WOLFS_METHOD_DAMPED_FG:
              UAdsorbateChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                        -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                         (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                        (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
               break;
            case EWALD:
              UAdsorbateChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (erfc(Alpha[CurrentSystem]*r)/r);
              break;
            default:
              printf("Unknown charge-charge method in 'CalculateFrameworkShiftEnergyDifferenceAdsorbate'\n");
              exit(0);
              break;
          }
        }
      }

      TypeMolB=Adsorbates[CurrentSystem][m].Type;
      for(j=0;j<Components[TypeMolB].NumberOfBondDipoles;j++)
      {
        B1=Components[TypeMolB].BondDipoles[j].A;
        B2=Components[TypeMolB].BondDipoles[j].B;
        posB1=Adsorbates[CurrentSystem][m].Atoms[B1].Position;
        posB2=Adsorbates[CurrentSystem][m].Atoms[B2].Position;
        DipoleMagnitudeB=Components[TypeMolB].BondDipoleMagnitude[j];
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt1=1.0/(r*rr);
              if(rr>CutOffChargeBondDipoleSwitchSquared)
                SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                               SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              break;
            default:
              printf("Unknown charge-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceAdsorbate'\n");
              exit(0);
              break;
          }

          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt1=1.0/(r*rr);
              if(rr>CutOffChargeBondDipoleSwitchSquared)
                SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                               SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              break;
            default:
              printf("Unknown charge-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceAdsorbate'\n");
              exit(0);
              break;
          }
    
          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
        }
      }
    }
  }

  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
    dipoleA_new.x=posA2.x-posA1.x;
    dipoleA_new.y=posA2.y-posA1.y;
    dipoleA_new.z=posA2.z-posA1.z;
    dipoleA_new=ApplyBoundaryCondition(dipoleA_new);
    posA_new.x=posA1.x+0.5*dipoleA_new.x;
    posA_new.y=posA1.y+0.5*dipoleA_new.y;
    posA_new.z=posA1.z+0.5*dipoleA_new.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z));
    dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
    dipoleA_old.x=posA2.x-posA1.x;
    dipoleA_old.y=posA2.y-posA1.y;
    dipoleA_old.z=posA2.z-posA1.z;
    dipoleA_old=ApplyBoundaryCondition(dipoleA_old);
    posA_old.x=posA1.x+0.5*dipoleA_old.x;
    posA_old.y=posA1.y+0.5*dipoleA_old.y;
    posA_old.z=posA1.z+0.5*dipoleA_old.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z));
    dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      TypeMolB=Adsorbates[CurrentSystem][m].Type;
      for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Adsorbates[CurrentSystem][m].Atoms[j].Type;
        posB=Adsorbates[CurrentSystem][m].Atoms[j].Position;
        chargeB=Adsorbates[CurrentSystem][m].Atoms[j].Charge;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt1=1.0/(r*rr);
              if(rr>CutOffChargeBondDipoleSwitchSquared)
                SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                               SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              break;
            default:
              printf("Unknown charge-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceAdsorbate'\n");
              exit(0);
              break;
          }

          cosA=dipoleA_new.x*dr.x+dipoleA_new.y*dr.y+dipoleA_new.z*dr.z;
          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt1=1.0/(r*rr);
              if(rr>CutOffChargeBondDipoleSwitchSquared)
                SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                               SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              break;
            default:
              printf("Unknown charge-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceAdsorbate'\n");
              exit(0);
              break;
          }

          cosA=dipoleA_old.x*dr.x+dipoleA_old.y*dr.y+dipoleA_old.z*dr.z;
          UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
        }
      }

      for(j=0;j<Components[TypeMolB].NumberOfBondDipoles;j++)
      {
        B1=Components[TypeMolB].BondDipoles[j].A;
        B2=Components[TypeMolB].BondDipoles[j].B;
        posB1=Adsorbates[CurrentSystem][m].Atoms[B1].Position;
        posB2=Adsorbates[CurrentSystem][m].Atoms[B2].Position;
        DipoleMagnitudeB=Components[TypeMolB].BondDipoleMagnitude[j];
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffBondDipoleBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=Bt2=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/r;
              Bt1=1.0/(rr*r);
              Bt2=3.0/(SQR(rr)*r);
              break;
            case SMOOTHED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              Bt2=3.0/(r*rr*rr);
              if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
              {
                SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                               SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
              }
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                  4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
              break;
            default:
              printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceAdsorbate'\n");
              exit(0);
              break;
          }

          cosAB=dipoleA_new.x*dipoleB.x+dipoleA_new.y*dipoleB.y+dipoleA_new.z*dipoleB.z;
          cosA=dipoleA_new.x*dr.x+dipoleA_new.y*dr.y+dipoleA_new.z*dr.z;
          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffBondDipoleBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=Bt2=0.0;
              break;
           case TRUNCATED_COULOMB:
           case SHIFTED_COULOMB:
              Bt0=1.0/r;
              Bt1=1.0/(rr*r);
              Bt2=3.0/(SQR(rr)*r);
              break;
            case SMOOTHED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              Bt2=3.0/(r*rr*rr);
              if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
              {
                SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                               SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
              }
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                  4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
              break;
            default:
              printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceAdsorbate'\n");
              exit(0);
              break;
          }

          cosAB=dipoleA_old.x*dipoleB.x+dipoleA_old.y*dipoleB.y+dipoleA_old.z*dipoleB.z;
          cosA=dipoleA_old.x*dr.x+dipoleA_old.y*dr.y+dipoleA_old.z*dr.z;
          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
        }
      }
    }
  }
}

void CalculateFrameworkShiftEnergyDifferenceCationVDW(void)
{
  int i,j,m,typeA,typeB;
  POINT posB,posA_new,posA_old;
  REAL rr;
  VECTOR dr;

  // Framework-Cation energy
  UCationVDWDelta[CurrentSystem]=0.0;
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][i].AnisotropicPosition;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferenceAnisotropicPosition;

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Cations[CurrentSystem][m].Atoms[j].Type;
        posB=Cations[CurrentSystem][m].Atoms[j].AnisotropicPosition;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<CutOffVDWSquared)
          UCationVDWDelta[CurrentSystem]+=PotentialValue(typeA,typeB,rr);

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<CutOffVDWSquared)
          UCationVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);
      }
    }
  }
}


void CalculateFrameworkShiftEnergyDifferenceCationCharge(void)
{
  int i,j,m,typeA,typeB;
  int A1,A2,B1,B2,TypeMolB;
  POINT posB,posA_new,posA_old;
  REAL rr,r,chargeA,chargeB,Bt0,Bt1,Bt2,cosA,cosB,cosAB,energy;
  VECTOR posA1,posA2,dipoleA_new,dipoleA_old;
  VECTOR posB1,posB2,dipoleB;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,temp;
  REAL SwitchingValue,TranslationValue;
  VECTOR dr;

  // Framework-Cation energy
  UCationChargeChargeRealDelta[CurrentSystem]=0.0;
  UCationChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;


  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferencePosition;
    chargeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Cations[CurrentSystem][m].Atoms[j].Type;
        posB=Cations[CurrentSystem][m].Atoms[j].Position;
        chargeB=Cations[CurrentSystem][m].Atoms[j].Charge;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeChargeSquared)
        {
          r=sqrt(rr);
          switch(ChargeMethod)
          {
            case NONE:
              break;
            case TRUNCATED_COULOMB:
              UCationChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
              break;
            case SHIFTED_COULOMB:
              UCationChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
              break;
            case SMOOTHED_COULOMB:
              energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
              if(rr>CutOffChargeChargeSwitchSquared)
              {
                SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                               SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                 (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                  SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                  SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                energy=energy*SwitchingValue+TranslationValue;
              }
              UCationChargeChargeRealDelta[CurrentSystem]+=energy;
              break;
            case WOLFS_METHOD_DAMPED_FG:
              UCationChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                        -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                         (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                        (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
               break;
            case EWALD:
              UCationChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (erfc(Alpha[CurrentSystem]*r)/r);
              break;
            default:
              printf("Unknown charge-charge method in 'CalculateFrameworkShiftEnergyDifferenceCation'\n");
              exit(0);
              break;
          }
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeChargeSquared)
        {
          r=sqrt(rr);
          switch(ChargeMethod)
          {
            case NONE:
              break;
            case TRUNCATED_COULOMB:
              UCationChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
              break;
            case SHIFTED_COULOMB:
              UCationChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
              break;
            case SMOOTHED_COULOMB:
              energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
              if(rr>CutOffChargeChargeSwitchSquared)
              {
                SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                               SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                 (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                  SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                  SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                energy=energy*SwitchingValue+TranslationValue;
              }
              UCationChargeChargeRealDelta[CurrentSystem]-=energy;
              break;
            case WOLFS_METHOD_DAMPED_FG:
              UCationChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                        -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                         (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                        (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
               break;
            case EWALD:
              UCationChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (erfc(Alpha[CurrentSystem]*r)/r);
              break;
            default:
              printf("Unknown charge-charge method in 'CalculateFrameworkShiftEnergyDifferenceCation'\n");
              exit(0);
              break;
          }
        }
      }

      TypeMolB=Cations[CurrentSystem][m].Type;
      for(j=0;j<Components[TypeMolB].NumberOfBondDipoles;j++)
      {
        B1=Components[TypeMolB].BondDipoles[j].A;
        B2=Components[TypeMolB].BondDipoles[j].B;
        posB1=Cations[CurrentSystem][m].Atoms[B1].Position;
        posB2=Cations[CurrentSystem][m].Atoms[B2].Position;
        DipoleMagnitudeB=Components[TypeMolB].BondDipoleMagnitude[j];
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt1=1.0/(r*rr);
              if(rr>CutOffChargeBondDipoleSwitchSquared)
                SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                               SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              break;
            default:
              printf("Unknown charge-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceCation'\n");
              exit(0);
              break;
          }

          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          UCationChargeBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt1=1.0/(r*rr);
              if(rr>CutOffChargeBondDipoleSwitchSquared)
                SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                               SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              break;
            default:
              printf("Unknown charge-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceCation'\n");
              exit(0);
              break;
          }
    
          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          UCationChargeBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
        }
      }
    }
  }

  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
    dipoleA_new.x=posA2.x-posA1.x;
    dipoleA_new.y=posA2.y-posA1.y;
    dipoleA_new.z=posA2.z-posA1.z;
    dipoleA_new=ApplyBoundaryCondition(dipoleA_new);
    posA_new.x=posA1.x+0.5*dipoleA_new.x;
    posA_new.y=posA1.y+0.5*dipoleA_new.y;
    posA_new.z=posA1.z+0.5*dipoleA_new.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z));
    dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
    dipoleA_old.x=posA2.x-posA1.x;
    dipoleA_old.y=posA2.y-posA1.y;
    dipoleA_old.z=posA2.z-posA1.z;
    dipoleA_old=ApplyBoundaryCondition(dipoleA_old);
    posA_old.x=posA1.x+0.5*dipoleA_old.x;
    posA_old.y=posA1.y+0.5*dipoleA_old.y;
    posA_old.z=posA1.z+0.5*dipoleA_old.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z));
    dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      TypeMolB=Cations[CurrentSystem][m].Type;
      for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
      {
        typeB=Cations[CurrentSystem][m].Atoms[j].Type;
        posB=Cations[CurrentSystem][m].Atoms[j].Position;
        chargeB=Cations[CurrentSystem][m].Atoms[j].Charge;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt1=1.0/(r*rr);
              if(rr>CutOffChargeBondDipoleSwitchSquared)
                SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                               SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              break;
            default:
              printf("Unknown charge-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceCation'\n");
              exit(0);
              break;
          }

          cosA=dipoleA_new.x*dr.x+dipoleA_new.y*dr.y+dipoleA_new.z*dr.z;
          UCationChargeBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt1=1.0/(r*rr);
              if(rr>CutOffChargeBondDipoleSwitchSquared)
                SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                               SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              break;
            default:
              printf("Unknown charge-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceCation'\n");
              exit(0);
              break;
          }

          cosA=dipoleA_old.x*dr.x+dipoleA_old.y*dr.y+dipoleA_old.z*dr.z;
          UCationChargeBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
        }
      }

      for(j=0;j<Components[TypeMolB].NumberOfBondDipoles;j++)
      {
        B1=Components[TypeMolB].BondDipoles[j].A;
        B2=Components[TypeMolB].BondDipoles[j].B;
        posB1=Cations[CurrentSystem][m].Atoms[B1].Position;
        posB2=Cations[CurrentSystem][m].Atoms[B2].Position;
        DipoleMagnitudeB=Components[TypeMolB].BondDipoleMagnitude[j];
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffBondDipoleBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=Bt2=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/r;
              Bt1=1.0/(rr*r);
              Bt2=3.0/(SQR(rr)*r);
              break;
            case SMOOTHED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              Bt2=3.0/(r*rr*rr);
              if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
              {
                SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                               SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
              }
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                  4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
              break;
            default:
              printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceCation'\n");
              exit(0);
              break;
          }

          cosAB=dipoleA_new.x*dipoleB.x+dipoleA_new.y*dipoleB.y+dipoleA_new.z*dipoleB.z;
          cosA=dipoleA_new.x*dr.x+dipoleA_new.y*dr.y+dipoleA_new.z*dr.z;
          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffBondDipoleBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt0=Bt1=Bt2=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt0=1.0/r;
              Bt1=1.0/(rr*r);
              Bt2=3.0/(SQR(rr)*r);
              break;
            case SMOOTHED_COULOMB:
              Bt0=1.0/(r);
              Bt1=1.0/(r*rr);
              Bt2=3.0/(r*rr*rr);
              if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
              {
                SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                               SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
              }
              break;
            case EWALD:
              Bt0=erfc(Alpha[CurrentSystem]*r)/r;
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                  4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
              break;
            default:
              printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkShiftEnergyDifferenceCation'\n");
              exit(0);
              break;
          }

          cosAB=dipoleA_old.x*dipoleB.x+dipoleA_old.y*dipoleB.y+dipoleA_old.z*dipoleB.z;
          cosA=dipoleA_old.x*dr.x+dipoleA_old.y*dr.y+dipoleA_old.z*dr.z;
          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          UCationBondDipoleBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
        }
      }
    }
  }
}


REAL CalculateEnergyDifferenceFrameworkMoveVDW(int atom_id,VECTOR posA,int typeA)
{
  int j,f1,typeB;
  REAL rr;
  VECTOR dr,posB;
  REAL UHostVDW;

  UHostVDW=0.0;
  if(!InternalFrameworkLennardJonesInteractions) return 0.0;

  // loop over all frameworks
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
    {
      if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][atom_id][j],0))
      {
        posB=Framework[CurrentSystem].Atoms[f1][j].AnisotropicPosition;
        typeB=Framework[CurrentSystem].Atoms[f1][j].Type;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<CutOffVDWSquared)
          UHostVDW+=PotentialValue(typeA,typeB,rr);
      }
    }
  }
  return UHostVDW;
}


/*
void CalculateEnergyDifferenceFrameworkMoveVDW(int atom_id)
{
  int i,j,f1,A,typeA,typeB;
  REAL rr;
  VECTOR dr,posA_new,posA_old,posB;
  int Atoms[20];
  int NumberOfAtoms;

  UHostVDWDelta[CurrentSystem]=0.0;
  if(!InternalFrameworkLennardJonesInteractions) return;

  NumberOfAtoms=1;
  Atoms[0]=atom_id;

  // add neighbors for anisotropic sites
  for(i=0;i<Framework[CurrentSystem].Connectivity[CurrentFramework][atom_id];i++)
    Atoms[NumberOfAtoms++]=Framework[CurrentSystem].Neighbours[CurrentFramework][atom_id][i];

  for(i=0;i<NumberOfAtoms;i++)
  {
    atom_id=Atoms[i];
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].AnisotropicPosition;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].ReferenceAnisotropicPosition;
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Type;

    // loop over all frameworks
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
      {
        if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][atom_id][j],0))
        {
          posB=Framework[CurrentSystem].Atoms[f1][j].AnisotropicPosition;
          typeB=Framework[CurrentSystem].Atoms[f1][j].Type;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
            UHostVDWDelta[CurrentSystem]+=PotentialValue(typeA,typeB,rr);

          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
            UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);
        }
      }
    }
  }
}
*/


void CalculateEnergyDifferenceFrameworkMoveCharge(int atom_id)
{
  int i,j,f1,typeA,typeB;
  int A1,A2,B1,B2;
  REAL rr,r,chargeA,chargeB;
  VECTOR dr,posA_new,posA_old,posB,posA1,posA2,posB1,posB2,posB_new,posB_old;
  VECTOR dipoleB,dipoleA_new,dipoleA_old,dipoleB_new,dipoleB_old;
  REAL Bt1,Bt2,cosA,cosB,cosAB,temp,energy;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL SwitchingValue,TranslationValue;

  UHostChargeChargeRealDelta[CurrentSystem]=0.0;
  UHostChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UHostBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Position;
  posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].ReferencePosition;
  typeA=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Type;
  chargeA=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Charge;

  // loop over all frameworks
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
    {
      if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][atom_id][j],1))
      {
        posB=Framework[CurrentSystem].Atoms[f1][j].Position;
        typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
        chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeChargeSquared)
        {
          r=sqrt(rr);
          switch(ChargeMethod)
          {
            case NONE:
              break;
            case TRUNCATED_COULOMB:
              UHostChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
              break;
            case SHIFTED_COULOMB:
              UHostChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
              break;
            case SMOOTHED_COULOMB:
              energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
              if(rr>CutOffChargeChargeSwitchSquared)
              {
                SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                               SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                 (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                  SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                  SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                energy=energy*SwitchingValue+TranslationValue;
              }
              UHostChargeChargeRealDelta[CurrentSystem]+=energy;
              break;
            case WOLFS_METHOD_DAMPED_FG:
              UHostChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                        -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                         (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                        (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
               break;
            case EWALD:
              UHostChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                        (erfc(Alpha[CurrentSystem]*r)/r);
              break;
            default:
              printf("Unknown charge-charge method in 'CalculateEnergyDifferenceFrameworkMove'\n");
              exit(0);
              break;
          }
        }
        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeChargeSquared)
        {
          r=sqrt(rr);
          switch(ChargeMethod)
          {
            case NONE:
              break;
            case TRUNCATED_COULOMB:
              UHostChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
              break;
            case SHIFTED_COULOMB:
              UHostChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
              break;
            case SMOOTHED_COULOMB:
              energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
              if(rr>CutOffChargeChargeSwitchSquared)
              {
                SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                               SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                 (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                  SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                  SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                energy=energy*SwitchingValue+TranslationValue;
              }
              UHostChargeChargeRealDelta[CurrentSystem]-=energy;
              break;
            case WOLFS_METHOD_DAMPED_FG:
              UHostChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                        -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                         (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                        (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
               break;
            case EWALD:
              UHostChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                        (erfc(Alpha[CurrentSystem]*r)/r);
              break;
            default:
              printf("Unknown charge-charge method in 'CalculateEnergyDifferenceFrameworkMove'\n");
              exit(0);
              break;
          }
        }
      }
    }

    // compute the charge-bond-dipole interactions
    for(j=0;j<Framework[CurrentSystem].NumberOfBondDipoles[f1];j++)
    {
      if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][atom_id][j],2))
      {
        DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][j];
        B1=Framework[CurrentSystem].BondDipoles[f1][j].A;
        B2=Framework[CurrentSystem].BondDipoles[f1][j].B;
        posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
        posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA_new.x-posB.x;
        dr.y=posA_new.y-posB.y;
        dr.z=posA_new.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt1=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt1=1.0/(r*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt1=1.0/(r*rr);
              if(rr>CutOffChargeBondDipoleSwitchSquared)
                SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                               SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
              break;
            case EWALD:
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              break;
            default:
              printf("Unknown charge-bonddipole method in 'CalculateEnergyDifferenceFrameworkMove'\n");
              exit(0);
              break;
          }

          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          UHostChargeBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
        }

        dr.x=posA_old.x-posB.x;
        dr.y=posA_old.y-posB.y;
        dr.z=posA_old.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<CutOffChargeBondDipoleSquared)
        {
          r=sqrt(rr);
          SwitchingValue=1.0;
          switch(ChargeMethod)
          {
            case NONE:
              Bt1=0.0;
              break;
            case TRUNCATED_COULOMB:
            case SHIFTED_COULOMB:
              Bt1=1.0/(r*rr);
              break;
            case SMOOTHED_COULOMB:
              Bt1=1.0/(r*rr);
              if(rr>CutOffChargeBondDipoleSwitchSquared)
                SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                               SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
              break;
            case EWALD:
              Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                  erfc(Alpha[CurrentSystem]*r)/(rr*r);
              break;
            default:
              printf("Unknown charge-bonddipole method in 'CalculateEnergyDifferenceFrameworkMove'\n");
              exit(0);
              break;
          }

          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          UHostChargeBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
        }
      }
    }
  }

  // loop over all the atoms of the currently selected framework
  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    // check if the current-atom is part of the current bond-dipole i
    if((atom_id==A1)||(atom_id==A2))
    {
      posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
      dipoleA_new.x=posA2.x-posA1.x;
      dipoleA_new.y=posA2.y-posA1.y;
      dipoleA_new.z=posA2.z-posA1.z;
      dipoleA_new=ApplyBoundaryCondition(dipoleA_new);
      posA_new.x=posA1.x+0.5*dipoleA_new.x;
      posA_new.y=posA1.y+0.5*dipoleA_new.y;
      posA_new.z=posA1.z+0.5*dipoleA_new.z;
      temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z));
      dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;

      posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
      posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
      dipoleA_old.x=posA2.x-posA1.x;
      dipoleA_old.y=posA2.y-posA1.y;
      dipoleA_old.z=posA2.z-posA1.z;
      dipoleA_old=ApplyBoundaryCondition(dipoleA_old);
      posA_old.x=posA1.x+0.5*dipoleA_old.x;
      posA_old.y=posA1.y+0.5*dipoleA_old.y;
      posA_old.z=posA1.z+0.5*dipoleA_old.z;
      temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z));
      dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;

      // loop over all frameworks
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][j][i],2))
          {
            posB_new=Framework[CurrentSystem].Atoms[f1][j].Position;
            posB_old=Framework[CurrentSystem].Atoms[f1][j].ReferencePosition;
            typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
            chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

            dr.x=posA_new.x-posB_new.x;
            dr.y=posA_new.y-posB_new.y;
            dr.z=posA_new.z-posB_new.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt1=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt1=1.0/(r*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateEnergyDifferenceFrameworkMove'\n");
                  exit(0);
                  break;
              }

              cosA=dipoleA_new.x*dr.x+dipoleA_new.y*dr.y+dipoleA_new.z*dr.z;
              UHostChargeBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
            }

            dr.x=posA_old.x-posB_old.x;
            dr.y=posA_old.y-posB_old.y;
            dr.z=posA_old.z-posB_old.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt1=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt1=1.0/(r*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateEnergyDifferenceFrameworkMove'\n");
                  exit(0);
                  break;
              }
  
              cosA=dipoleA_old.x*dr.x+dipoleA_old.y*dr.y+dipoleA_old.z*dr.z;
              UHostChargeBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
            }
          }
        }

        // compute the bond-dipole-bond-dipole interactions
        for(j=0;j<Framework[CurrentSystem].NumberOfBondDipoles[f1];j++)
        {
          if(f1!=CurrentFramework?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][j],3))
          {
            DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][j];
            B1=Framework[CurrentSystem].BondDipoles[f1][j].A;
            B2=Framework[CurrentSystem].BondDipoles[f1][j].B;

            posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
            posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
            dipoleB_new.x=posB2.x-posB1.x;
            dipoleB_new.y=posB2.y-posB1.y;
            dipoleB_new.z=posB2.z-posB1.z;
            dipoleB_new=ApplyBoundaryCondition(dipoleB_new);
            posB_new.x=posB1.x+0.5*dipoleB_new.x;
            posB_new.y=posB1.y+0.5*dipoleB_new.y;
            posB_new.z=posB1.z+0.5*dipoleB_new.z;
            temp=DipoleMagnitudeB/sqrt(SQR(dipoleB_new.x)+SQR(dipoleB_new.y)+SQR(dipoleB_new.z));
            dipoleB_new.x*=temp; dipoleB_new.y*=temp; dipoleB_new.z*=temp;

            posB1=Framework[CurrentSystem].Atoms[f1][B1].ReferencePosition;
            posB2=Framework[CurrentSystem].Atoms[f1][B2].ReferencePosition;
            dipoleB_old.x=posB2.x-posB1.x;
            dipoleB_old.y=posB2.y-posB1.y;
            dipoleB_old.z=posB2.z-posB1.z;
            dipoleB_old=ApplyBoundaryCondition(dipoleB_old);
            posB_old.x=posB1.x+0.5*dipoleB_old.x;
            posB_old.y=posB1.y+0.5*dipoleB_old.y;
            posB_old.z=posB1.z+0.5*dipoleB_old.z;
            temp=DipoleMagnitudeB/sqrt(SQR(dipoleB_old.x)+SQR(dipoleB_old.y)+SQR(dipoleB_old.z));
            dipoleB_old.x*=temp; dipoleB_old.y*=temp; dipoleB_old.z*=temp;

            dr.x=posA_new.x-posB_new.x;
            dr.y=posA_new.y-posB_new.y;
            dr.z=posA_new.z-posB_new.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                  }
                  break;
                case EWALD:
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown bonddipole-bonddipole method in 'CalculateEnergyDifferenceFrameworkMove'\n");
                  exit(0);
                  break;
              }

              cosAB=dipoleA_new.x*dipoleB_new.x+dipoleA_new.y*dipoleB_new.y+dipoleA_new.z*dipoleB_new.z;
              cosA=dipoleA_new.x*dr.x+dipoleA_new.y*dr.y+dipoleA_new.z*dr.z;
              cosB=dipoleB_new.x*dr.x+dipoleB_new.y*dr.y+dipoleB_new.z*dr.z;
              UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            }

            dr.x=posA_old.x-posB_old.x;
            dr.y=posA_old.y-posB_old.y;
            dr.z=posA_old.z-posB_old.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                  }
                  break;
                case EWALD:
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown bonddipole-bonddipole method in 'CalculateEnergyDifferenceFrameworkMove'\n");
                  exit(0);
                  break;
              }

              cosAB=dipoleA_old.x*dipoleB_old.x+dipoleA_old.y*dipoleB_old.y+dipoleA_old.z*dipoleB_old.z;
              cosA=dipoleA_old.x*dr.x+dipoleA_old.y*dr.y+dipoleA_old.z*dr.z;
              cosB=dipoleB_old.x*dr.x+dipoleB_old.y*dr.y+dipoleB_old.z*dr.z;
              UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            }
          }
        }
      }
    }
  }
}

void CalculateFrameworkEnergyDifferenceShiftedFramework(void)
{
  int i,j,f1,typeA,typeB;
  int A1,A2,B1,B2;
  REAL rr,r,chargeA,chargeB,energy;
  VECTOR dr,posA_new,posA_old,posB,posA1,posA2,posB1,posB2;
  VECTOR dipoleB,dipoleA_new,dipoleA_old;
  REAL Bt1,Bt2,cosA,cosB,cosAB,temp;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL SwitchingValue,TranslationValue;

  // Framework-Adsorbate energy
  UHostVDWDelta[CurrentSystem]=0.0;
  UHostChargeChargeRealDelta[CurrentSystem]=0.0;
  UHostChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UHostBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  // loop over all the atoms of the currently selected framework
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    posA_new=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
    posA_old=Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferencePosition;
    typeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    chargeA=Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge;

    // loop over all frameworks
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      // except the current one (no self-interaction, simple translation of the whole framework)
      if(f1!=CurrentFramework)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][j].Position;
          typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
            if(InternalFrameworkLennardJonesInteractions)
              UHostVDWDelta[CurrentSystem]+=PotentialValue(typeA,typeB,rr);

          if(rr<CutOffChargeChargeSquared)
          {
            r=sqrt(rr);
            switch(ChargeMethod)
            {
              case NONE:
                break;
              case TRUNCATED_COULOMB:
                UHostChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                break;
              case SHIFTED_COULOMB:
                UHostChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                break;
              case SMOOTHED_COULOMB:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                if(rr>CutOffChargeChargeSwitchSquared)
                {
                  SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                 SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                  TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                   (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                    SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                    SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                  energy=energy*SwitchingValue+TranslationValue;
                }
                UHostChargeChargeRealDelta[CurrentSystem]+=energy;
                break;
              case WOLFS_METHOD_DAMPED_FG:
                UHostChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                          -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                           (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                          (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                break;
              case EWALD:
                UHostChargeChargeRealDelta[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                          (erfc(Alpha[CurrentSystem]*r)/r);
                break;
              default:
                printf("Unknown charge-charge method in 'CalculateFrameworkEnergyDifferenceShiftedFramework'\n");
                exit(0);
                break;
            }
          }
          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
            if(InternalFrameworkLennardJonesInteractions)
              UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);

          if(rr<CutOffChargeChargeSquared)
          {
            r=sqrt(rr);
            switch(ChargeMethod)
            {
              case NONE:
                break;
              case TRUNCATED_COULOMB:
                UHostChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                break;
              case SHIFTED_COULOMB:
                UHostChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                break;
              case SMOOTHED_COULOMB:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                if(rr>CutOffChargeChargeSwitchSquared)
                {
                  SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                 SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                  TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                   (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                    SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                    SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                  energy=energy*SwitchingValue+TranslationValue;
                }
                UHostChargeChargeRealDelta[CurrentSystem]-=energy;
                break;
              case WOLFS_METHOD_DAMPED_FG:
                UHostChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                          -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                           (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                          (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                break;
              case EWALD:
                UHostChargeChargeRealDelta[CurrentSystem]-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                          (erfc(Alpha[CurrentSystem]*r)/r);
                break;
              default:
                printf("Unknown charge-charge method in 'CalculateFrameworkEnergyDifferenceShiftedFramework'\n");
                exit(0);
                break;
            }
          }
        }

        // compute the charge-bond-dipole interactions
        for(j=0;j<Framework[CurrentSystem].NumberOfBondDipoles[f1];j++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][j];
          B1=Framework[CurrentSystem].BondDipoles[f1][j].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][j].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt1=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt1=1.0/(r*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkEnergyDifferenceShiftedFramework'\n");
                exit(0);
                break;
            }

            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            UHostChargeBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
          }

          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt1=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt1=1.0/(r*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkEnergyDifferenceShiftedFramework'\n");
                exit(0);
                break;
            }

            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            UHostChargeBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
          }
        }
      }
    }
  }

  // loop over all the atoms of the currently selected framework
  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
    dipoleA_new.x=posA2.x-posA1.x;
    dipoleA_new.y=posA2.y-posA1.y;
    dipoleA_new.z=posA2.z-posA1.z;
    dipoleA_new=ApplyBoundaryCondition(dipoleA_new);
    posA_new.x=posA1.x+0.5*dipoleA_new.x;
    posA_new.y=posA1.y+0.5*dipoleA_new.y;
    posA_new.z=posA1.z+0.5*dipoleA_new.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z));
    dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
    dipoleA_old.x=posA2.x-posA1.x;
    dipoleA_old.y=posA2.y-posA1.y;
    dipoleA_old.z=posA2.z-posA1.z;
    dipoleA_old=ApplyBoundaryCondition(dipoleA_old);
    posA_old.x=posA1.x+0.5*dipoleA_old.x;
    posA_old.y=posA1.y+0.5*dipoleA_old.y;
    posA_old.z=posA1.z+0.5*dipoleA_old.z;
    temp=DipoleMagnitudeA/sqrt(SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z));
    dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;


    // loop over all frameworks
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      // except the current one (no self-interaction, simple translation of the whole framework)
      if(f1!=CurrentFramework)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          posB=Framework[CurrentSystem].Atoms[f1][j].Position;
          typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
          chargeB=Framework[CurrentSystem].Atoms[f1][j].Charge;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt1=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt1=1.0/(r*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkEnergyDifferenceShiftedFramework'\n");
                exit(0);
                break;
            }

            cosA=dipoleA_new.x*dr.x+dipoleA_new.y*dr.y+dipoleA_new.z*dr.z;
            UHostChargeBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
          }

          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt1=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt1=1.0/(r*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkEnergyDifferenceShiftedFramework'\n");
                exit(0);
                break;
            }

            cosA=dipoleA_old.x*dr.x+dipoleA_old.y*dr.y+dipoleA_old.z*dr.z;
            UHostChargeBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
          }
        }

        // compute the bond-dipole-bond-dipole interactions
        for(j=0;j<Framework[CurrentSystem].NumberOfBondDipoles[f1];j++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][j];
          B1=Framework[CurrentSystem].BondDipoles[f1][j].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][j].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkEnergyDifferenceShiftedFramework'\n");
                exit(0);
                break;
            }

            cosAB=dipoleA_new.x*dipoleB.x+dipoleA_new.y*dipoleB.y+dipoleA_new.z*dipoleB.z;
            cosA=dipoleA_new.x*dr.x+dipoleA_new.y*dr.y+dipoleA_new.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
          }

          dr.x=posA_old.x-posB.x;
          dr.y=posA_old.y-posB.y;
          dr.z=posA_old.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkEnergyDifferenceShiftedFramework'\n");
                exit(0);
                break;
            }

            cosAB=dipoleA_old.x*dipoleB.x+dipoleA_old.y*dipoleB.y+dipoleA_old.z*dipoleB.z;
            cosA=dipoleA_old.x*dr.x+dipoleA_old.y*dr.y+dipoleA_old.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
          }
        }
      }
    }
  }
}

int CalculateFrameworkIntraReplicaVDWEnergy(void)
{
  int i,j,typeA,typeB,start;
  REAL energy,rr;
  VECTOR posA,posB,dr;
  int f1,f2,ncell,index_j;

  // Framework-Framework energy
  UHostHostVDW[CurrentSystem]=0.0;

  if(!InternalFrameworkLennardJonesInteractions) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition;

        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          if((f1==f2)&&(ncell==0)) start=i+1;
          else start=0;
          for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
          {
            index_j=j+ncell*Framework[CurrentSystem].NumberOfAtoms[f2];
            if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][index_j],0))
            {
              typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
              posB=Framework[CurrentSystem].Atoms[f2][j].AnisotropicPosition;

              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr);

                if(ncell==0)
                  UHostHostVDW[CurrentSystem]+=energy;
                else
                  UHostHostVDW[CurrentSystem]+=0.5*energy;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkIntraReplicaChargeChargeEnergy(void)
{
  int i,j,typeA,typeB,start;
  REAL chargeA,chargeB;
  REAL r,rr,energy;
  VECTOR posA,posB,dr;
  int f1,f2,ncell,index_j;
  REAL SwitchingValue,TranslationValue;

  // Framework-Framework energy
  UHostHostChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(f2=f1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;
        chargeA=Framework[CurrentSystem].Atoms[f1][i].Charge;

        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          if((f1==f2)&&(ncell==0)) start=i+1;
          else start=0;
          for(j=start;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
          {
            index_j=j+ncell*Framework[CurrentSystem].NumberOfAtoms[f2];
            if(f1!=f2?TRUE:!BITVAL(Framework[CurrentSystem].ExclusionMatrix[f1][i][index_j],0))
            {
              typeB=Framework[CurrentSystem].Atoms[f2][j].Type;
              posB=Framework[CurrentSystem].Atoms[f2][j].Position;
              chargeB=Framework[CurrentSystem].Atoms[f2][j].Charge;

              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared)
              {
                r=sqrt(rr);
                switch(ChargeMethod)
                {
                  case NONE:
                    energy=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                    break;
                  case SHIFTED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case WOLFS_METHOD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(rr>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                     SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                      (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                       SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                       SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                      energy=energy*SwitchingValue+TranslationValue;
                    }
                    break;
                  case WOLFS_METHOD_DAMPED_FG:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                           (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                           (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                    break;
                  case EWALD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                               (erfc(Alpha[CurrentSystem]*r)/r);
                    break;
                  default:
                    printf("Unknown charge-charge method in 'CalculateFrameworkIntraReplicaChargeChargeEnergy'\n");
                    exit(0);
                    break;
                }
                if(ncell==0)
                  UHostHostChargeChargeReal[CurrentSystem]+=energy;
                else
                  UHostHostChargeChargeReal[CurrentSystem]+=0.5*energy;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}


void CalculateFrameworkAdsorbateReplicaVDWEnergy(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr;
  VECTOR posA,posB,dr;
  int ncell;

  UHostAdsorbateVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
      {
        UHostAdsorbateVDW[CurrentSystem]+=InterpolateVDWGrid(typeA,posA);
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
                UHostAdsorbateVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);
            }
          }
        }
      }
    }
  }
}

void CalculateFrameworkCationReplicaVDWEnergy(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr;
  VECTOR posA,posB,dr;
  int ncell;

  UHostCationVDW[CurrentSystem]=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
      {
        UHostCationVDW[CurrentSystem]+=InterpolateVDWGrid(typeA,posA);
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
                UHostCationVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);
            }
          }
        }
      }
    }
  }
}

int CalculateFrameworkAdsorbateReplicaChargeChargeEnergy(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,r,chargeA,chargeB;
  REAL energy;
  VECTOR posA,posB,dr;
  int ncell;
  REAL SwitchingValue;
  REAL TranslationValue;

  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
      chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
      {
        UHostAdsorbateChargeChargeReal[CurrentSystem]+=InterpolateCoulombGrid(typeA,posA);
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared)
              {
                r=sqrt(rr);
                switch(ChargeMethod)
                {
                  case NONE:
                    energy=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                    break;
                  case SHIFTED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case WOLFS_METHOD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(rr>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                     SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                      (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                       SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                       SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                      energy=energy*SwitchingValue+TranslationValue;
                    }
                    break;
                  case WOLFS_METHOD_DAMPED_FG:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                           (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                           (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                    break;
                  case EWALD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                              (erfc(Alpha[CurrentSystem]*r)/r);
                    break;
                  default:
                    printf("Unknown charge-charge method in 'CalculateFrameworkAdsorbateReplicaChargeChargeEnergy'\n");
                    exit(0);
                    break;
                }

                // energy
                UHostAdsorbateChargeChargeReal[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationReplicaChargeChargeEnergy(void)
{
  int i,j,k,f1;
  int typeA,typeB,type;
  REAL rr,r,chargeA,chargeB;
  REAL energy;
  VECTOR posA,posB,dr;
  int ncell;
  REAL SwitchingValue,TranslationValue;

  UHostCationChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[j].Type;
      posA=Cations[CurrentSystem][i].Atoms[j].Position;
      chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
      {
        UHostCationChargeChargeReal[CurrentSystem]+=InterpolateCoulombGrid(typeA,posA);
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared)
              {
                r=sqrt(rr);
                switch(ChargeMethod)
                {
                  case NONE:
                    energy=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                    break;
                  case SHIFTED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case WOLFS_METHOD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(rr>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                     SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                      (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                       SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                       SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                      energy=energy*SwitchingValue+TranslationValue;
                    }
                    break;
                  case WOLFS_METHOD_DAMPED_FG:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                           (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                           (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                    break;
                  case EWALD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                              (erfc(Alpha[CurrentSystem]*r)/r);
                    break;
                  default:
                    printf("Unknown charge-charge method in 'CalculateFrameworkCationReplicaChargeChargeEnergy'\n");
                    exit(0);
                    break;
                }

                // energy
                UHostCationChargeChargeReal[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

// Monte Carlo routines
// ====================
int CalculateFrameworkAdsorbateVDWEnergyDifference(int m)
{
  int i,j,k,nr_atoms,typeA,typeB,f1;
  POINT posA_new,posA_old,posB;
  REAL rr,energy;
  VECTOR dr,s;
  int TypeMolA;
  int icell0,icell;
  int ncell;

  // Framework-Adsorbate energy
  OVERLAP=FALSE;
  UHostVDWDelta[CurrentSystem]=0.0;
  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    // grid interpolation; transfer to the correct coordinates
    nr_atoms=Adsorbates[CurrentSystem][m].NumberOfAtoms;
    TypeMolA=Adsorbates[CurrentSystem][m].Type;
    for(j=0;j<nr_atoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][m].Atoms[j].Type;
      posA_old=Adsorbates[CurrentSystem][m].Atoms[j].AnisotropicPosition;
      posA_new=TrialAnisotropicPosition[CurrentSystem][j];

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
      {
        UHostVDWDelta[CurrentSystem]+=InterpolateVDWGrid(typeA,posA_new);
        UHostVDWDelta[CurrentSystem]-=InterpolateVDWGrid(typeA,posA_old);
      }
      else if(UseCellLists[CurrentSystem])
      {
        // convert from xyz to abc
        s.x=InverseBox[CurrentSystem].ax*posA_new.x+InverseBox[CurrentSystem].bx*posA_new.y+InverseBox[CurrentSystem].cx*posA_new.z;
        s.y=InverseBox[CurrentSystem].ay*posA_new.x+InverseBox[CurrentSystem].by*posA_new.y+InverseBox[CurrentSystem].cy*posA_new.z;
        s.z=InverseBox[CurrentSystem].az*posA_new.x+InverseBox[CurrentSystem].bz*posA_new.y+InverseBox[CurrentSystem].cz*posA_new.z;

        // apply boundary condition
        s.x-=(REAL)NINT(s.x);
        s.y-=(REAL)NINT(s.y);
        s.z-=(REAL)NINT(s.z);

        // s between 0 and 1
        s.x+=0.5;
        s.y+=0.5;
        s.z+=0.5;

        // compute the corresponding cell-id
        icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
               ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
               ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          // loop over cells
          for(i=0;i<27;i++)
          {
            icell=CellListMap[CurrentSystem][icell0][i];

            k=Framework[CurrentSystem].CellListHead[f1][icell];

            while(k>=0)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UHostVDWDelta[CurrentSystem]+=energy;
              }

              k=Framework[CurrentSystem].CellList[f1][k];
            }
          }
        }

        // convert from xyz to abc
        s.x=InverseBox[CurrentSystem].ax*posA_old.x+InverseBox[CurrentSystem].bx*posA_old.y+InverseBox[CurrentSystem].cx*posA_old.z;
        s.y=InverseBox[CurrentSystem].ay*posA_old.x+InverseBox[CurrentSystem].by*posA_old.y+InverseBox[CurrentSystem].cy*posA_old.z;
        s.z=InverseBox[CurrentSystem].az*posA_old.x+InverseBox[CurrentSystem].bz*posA_old.y+InverseBox[CurrentSystem].cz*posA_old.z;

        // apply boundary condition
        s.x-=(REAL)NINT(s.x);
        s.y-=(REAL)NINT(s.y);
        s.z-=(REAL)NINT(s.z);

        // s between 0 and 1
        s.x+=0.5;
        s.y+=0.5;
        s.z+=0.5;

        // compute the corresponding cell-id
        icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
               ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
               ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          // loop over cells
          for(i=0;i<27;i++)
          {
            icell=CellListMap[CurrentSystem][icell0][i];

            k=Framework[CurrentSystem].CellListHead[f1][icell];

            while(k>=0)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);

              k=Framework[CurrentSystem].CellList[f1][k];
            }
          }
        }
      }
      else if(UseReplicas[CurrentSystem])
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UHostVDWDelta[CurrentSystem]+=energy;
              }

              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);
            }
          }
        }
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

            dr.x=posA_new.x-posB.x;
            dr.y=posA_new.y-posB.y;
            dr.z=posA_new.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              energy=PotentialValue(typeA,typeB,rr);
              if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
              UHostVDWDelta[CurrentSystem]+=energy;
            }

            dr.x=posA_old.x-posB.x;
            dr.y=posA_old.y-posB.y;
            dr.z=posA_old.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationVDWEnergyDifference(int m)
{
  int i,j,k,nr_atoms,typeA,typeB,f1;
  POINT posA_new,posA_old,posB;
  REAL rr,energy;
  VECTOR dr,s;
  int TypeMolA;
  int icell0,icell;
  int ncell;

  // Framework-Cation energy
  OVERLAP=FALSE;
  UHostVDWDelta[CurrentSystem]=0.0;
  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    // grid interpolation; transfer to the correct coordinates
    nr_atoms=Cations[CurrentSystem][m].NumberOfAtoms;
    TypeMolA=Cations[CurrentSystem][m].Type;
    for(j=0;j<nr_atoms;j++)
    {
      typeA=Cations[CurrentSystem][m].Atoms[j].Type;
      posA_old=Cations[CurrentSystem][m].Atoms[j].AnisotropicPosition;
      posA_new=TrialAnisotropicPosition[CurrentSystem][j];

      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(VDWGrid[typeA]))
      {
        UHostVDWDelta[CurrentSystem]+=InterpolateVDWGrid(typeA,posA_new);
        UHostVDWDelta[CurrentSystem]-=InterpolateVDWGrid(typeA,posA_old);
      }
      else if(UseCellLists[CurrentSystem])
      {
        // convert from xyz to abc
        s.x=InverseBox[CurrentSystem].ax*posA_new.x+InverseBox[CurrentSystem].bx*posA_new.y+InverseBox[CurrentSystem].cx*posA_new.z;
        s.y=InverseBox[CurrentSystem].ay*posA_new.x+InverseBox[CurrentSystem].by*posA_new.y+InverseBox[CurrentSystem].cy*posA_new.z;
        s.z=InverseBox[CurrentSystem].az*posA_new.x+InverseBox[CurrentSystem].bz*posA_new.y+InverseBox[CurrentSystem].cz*posA_new.z;

        // apply boundary condition
        s.x-=(REAL)NINT(s.x);
        s.y-=(REAL)NINT(s.y);
        s.z-=(REAL)NINT(s.z);

        // s between 0 and 1
        s.x+=0.5;
        s.y+=0.5;
        s.z+=0.5;

        // compute the corresponding cell-id
        icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
               ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
               ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          // loop over cells
          for(i=0;i<27;i++)
          {
            icell=CellListMap[CurrentSystem][icell0][i];

            k=Framework[CurrentSystem].CellListHead[f1][icell];

            while(k>=0)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UHostVDWDelta[CurrentSystem]+=energy;
              }

              k=Framework[CurrentSystem].CellList[f1][k];
            }
          }
        }

        // convert from xyz to abc
        s.x=InverseBox[CurrentSystem].ax*posA_old.x+InverseBox[CurrentSystem].bx*posA_old.y+InverseBox[CurrentSystem].cx*posA_old.z;
        s.y=InverseBox[CurrentSystem].ay*posA_old.x+InverseBox[CurrentSystem].by*posA_old.y+InverseBox[CurrentSystem].cy*posA_old.z;
        s.z=InverseBox[CurrentSystem].az*posA_old.x+InverseBox[CurrentSystem].bz*posA_old.y+InverseBox[CurrentSystem].cz*posA_old.z;

        // apply boundary condition
        s.x-=(REAL)NINT(s.x);
        s.y-=(REAL)NINT(s.y);
        s.z-=(REAL)NINT(s.z);

        // s between 0 and 1
        s.x+=0.5;
        s.y+=0.5;
        s.z+=0.5;

        // compute the corresponding cell-id
        icell0=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
               ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
               ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          // loop over cells
          for(i=0;i<27;i++)
          {
            icell=CellListMap[CurrentSystem][icell0][i];

            k=Framework[CurrentSystem].CellListHead[f1][icell];

            while(k>=0)
            {
              typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
              posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);

              k=Framework[CurrentSystem].CellList[f1][k];
            }
          }
        }
      }
      else if(UseReplicas[CurrentSystem])
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UHostVDWDelta[CurrentSystem]+=energy;
              }

              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
              UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);
            }
          }
        }
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            posB=Framework[CurrentSystem].Atoms[f1][k].AnisotropicPosition;

            dr.x=posA_new.x-posB.x;
            dr.y=posA_new.y-posB.y;
            dr.z=posA_new.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              energy=PotentialValue(typeA,typeB,rr);
              if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
              UHostVDWDelta[CurrentSystem]+=energy;
            }

            dr.x=posA_old.x-posB.x;
            dr.y=posA_old.y-posB.y;
            dr.z=posA_old.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UHostVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);
          }
        }
      }
    }
  }
  return 0.0;
}

void CalculateFrameworkAdsorbateChargeChargeEnergyDifference(int m)
{
  int j,k,nr_atoms,typeA,typeB,f1;
  POINT posA,posB;
  REAL rr,r,chargeA,chargeB,energy;
  REAL SwitchingValue,TranslationValue;
  VECTOR dr;
  int ncell;

  // Framework-Adsorbate energy
  UHostChargeChargeRealDelta[CurrentSystem]=0.0;
  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    // grid interpolation; transfer to the correct coordinates
    nr_atoms=Adsorbates[CurrentSystem][m].NumberOfAtoms;
    for(j=0;j<nr_atoms;j++)
    {
      typeA=Adsorbates[CurrentSystem][m].Atoms[j].Type;
      posA=Adsorbates[CurrentSystem][m].Atoms[j].Position;
      chargeA=Adsorbates[CurrentSystem][m].Atoms[j].Charge;
      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
      {
        if(ChargeMethod!=NONE)
        {
          UHostChargeChargeRealDelta[CurrentSystem]+=InterpolateCoulombGrid(typeA,TrialPosition[CurrentSystem][j]);
          UHostChargeChargeRealDelta[CurrentSystem]-=InterpolateCoulombGrid(typeA,posA);
        }
      }
      else if(UseReplicas[CurrentSystem])
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=TrialPosition[CurrentSystem][j].x-(posB.x+ReplicaShift[ncell].x);
              dr.y=TrialPosition[CurrentSystem][j].y-(posB.y+ReplicaShift[ncell].y);
              dr.z=TrialPosition[CurrentSystem][j].z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared)
              {
                r=sqrt(rr);
                switch(ChargeMethod)
                {
                  case NONE:
                    energy=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                    break;
                  case SHIFTED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(rr>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                     SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                      SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                      SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                      energy=energy*SwitchingValue+TranslationValue;
                    }
                    break;
                  case WOLFS_METHOD_DAMPED_FG:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                          (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                           (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                    break;
                  case EWALD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                 (erfc(Alpha[CurrentSystem]*r)/r);
                    break;
                  default:
                    printf("Unknown charge-charge method in 'CalculateFrameworkAdsorbateChargeChargeEnergyDifference'\n");
                    exit(0);
                    break;
                }
                UHostChargeChargeRealDelta[CurrentSystem]+=energy;
              }

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared)
              {
                r=sqrt(rr);
                switch(ChargeMethod)
                {
                  case NONE:
                    energy=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                    break;
                  case SHIFTED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(rr>CutOffChargeChargeSwitchSquared)
                    {
                       SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                      SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                       TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                        SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                        SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                       energy=energy*SwitchingValue+TranslationValue;
                    }
                    break;
                  case WOLFS_METHOD_DAMPED_FG:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                           (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                           (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                    break;
                  case EWALD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                (erfc(Alpha[CurrentSystem]*r)/r);
                    break;
                  default:
                    printf("Unknown charge-charge method in 'CalculateFrameworkAdsorbateChargeChargeEnergyDifference'\n");
                    exit(0);
                    break;
                }
                UHostChargeChargeRealDelta[CurrentSystem]-=energy;
              }
            }
          }
        }
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;
            dr.x=TrialPosition[CurrentSystem][j].x-posB.x;
            dr.y=TrialPosition[CurrentSystem][j].y-posB.y;
            dr.z=TrialPosition[CurrentSystem][j].z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              switch(ChargeMethod)
              {
                case NONE:
                  energy=0.0;
                  break;
                case TRUNCATED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  break;
                case SHIFTED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(rr>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                   SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
 
                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                    (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                     SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                     SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                    energy=energy*SwitchingValue+TranslationValue;
                  }
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                         -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                         (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                         (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  break;
                case EWALD:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (erfc(Alpha[CurrentSystem]*r)/r);
                  break;
                default:
                  printf("Unknown charge-charge method in 'CalculateFrameworkAdsorbateChargeChargeEnergyDifference'\n");
                  exit(0);
                  break;
               }
               UHostChargeChargeRealDelta[CurrentSystem]+=energy;
            }

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              switch(ChargeMethod)
              {
                case NONE:
                  energy=0.0;
                  break;
                case TRUNCATED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  break;
                case SHIFTED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(rr>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                   SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                    (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                     SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                     SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                    energy=energy*SwitchingValue+TranslationValue;
                  }
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                         -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                         (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                         (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  break;
                case EWALD:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (erfc(Alpha[CurrentSystem]*r)/r);
                  break;
                default:
                  printf("Unknown charge-charge method in 'CalculateFrameworkAdsorbateChargeChargeEnergyDifference'\n");
                  exit(0);
                  break;
              }
              UHostChargeChargeRealDelta[CurrentSystem]-=energy;
            }
          }
        }
      }
    }
  }
}

void CalculateFrameworkCationChargeChargeEnergyDifference(int m)
{
  int j,k,nr_atoms,typeA,typeB,f1;
  POINT posA,posB;
  REAL rr,r,chargeA,chargeB,energy;
  REAL SwitchingValue,TranslationValue;
  VECTOR dr;
  int ncell;

  // Framework-Cation energy
  UHostChargeChargeRealDelta[CurrentSystem]=0.0;
  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    // grid interpolation; transfer to the correct coordinates
    nr_atoms=Cations[CurrentSystem][m].NumberOfAtoms;
    for(j=0;j<nr_atoms;j++)
    {
      typeA=Cations[CurrentSystem][m].Atoms[j].Type;
      posA=Cations[CurrentSystem][m].Atoms[j].Position;
      chargeA=Cations[CurrentSystem][m].Atoms[j].Charge;
      if((Framework[CurrentSystem].FrameworkModel==GRID)&&(CoulombGrid))
      {
        if(ChargeMethod!=NONE)
        {
          UHostChargeChargeRealDelta[CurrentSystem]+=InterpolateCoulombGrid(typeA,TrialPosition[CurrentSystem][j]);
          UHostChargeChargeRealDelta[CurrentSystem]-=InterpolateCoulombGrid(typeA,posA);
        }
      }
      else if(UseReplicas[CurrentSystem])
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;

            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=TrialPosition[CurrentSystem][j].x-(posB.x+ReplicaShift[ncell].x);
              dr.y=TrialPosition[CurrentSystem][j].y-(posB.y+ReplicaShift[ncell].y);
              dr.z=TrialPosition[CurrentSystem][j].z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared)
              {
                r=sqrt(rr);
                switch(ChargeMethod)
                {
                  case NONE:
                    energy=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                    break;
                  case SHIFTED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(rr>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                     SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                      (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                       SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                       SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                      energy=energy*SwitchingValue+TranslationValue;
                    }
                    break;
                  case WOLFS_METHOD_DAMPED_FG:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                          (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                           (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                    break;
                  case EWALD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                 (erfc(Alpha[CurrentSystem]*r)/r);
                    break;
                  default:
                    printf("Unknown charge-charge method in 'CalculateFrameworkAdsorbateChargeChargeEnergyDifference'\n");
                    exit(0);
                    break;
                }
                UHostChargeChargeRealDelta[CurrentSystem]+=energy;
              }

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared)
              {
                r=sqrt(rr);
                switch(ChargeMethod)
                {
                  case NONE:
                    energy=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                    break;
                  case SHIFTED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(rr>CutOffChargeChargeSwitchSquared)
                    {
                       SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                      SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                       TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                        SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                        SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                       energy=energy*SwitchingValue+TranslationValue;
                    }
                    break;
                  case WOLFS_METHOD_DAMPED_FG:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                           (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                           (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                    break;
                  case EWALD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                 (erfc(Alpha[CurrentSystem]*r)/r);
                    break;
                  default:
                    printf("Unknown charge-charge method in 'CalculateFrameworkAdsorbateChargeChargeEnergyDifference'\n");
                    exit(0);
                    break;
                }
                UHostChargeChargeRealDelta[CurrentSystem]-=energy;
              }
            }
          }
        }
      }
      else
      {
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            typeB=Framework[CurrentSystem].Atoms[f1][k].Type;
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;
            dr.x=TrialPosition[CurrentSystem][j].x-posB.x;
            dr.y=TrialPosition[CurrentSystem][j].y-posB.y;
            dr.z=TrialPosition[CurrentSystem][j].z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              switch(ChargeMethod)
              {
                case NONE:
                  energy=0.0;
                  break;
                case TRUNCATED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  break;
                case SHIFTED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(rr>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                   SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
 
                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                    (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                     SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                     SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                    energy=energy*SwitchingValue+TranslationValue;
                  }
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                         -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                         (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                         (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  break;
                case EWALD:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (erfc(Alpha[CurrentSystem]*r)/r);
                  break;
                default:
                  printf("Unknown charge-charge method in 'CalculateFrameworkCationsChargeChargeEnergyDifference'\n");
                  exit(0);
                  break;
               }
               UHostChargeChargeRealDelta[CurrentSystem]+=energy;
            }

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              switch(ChargeMethod)
              {
                case NONE:
                  energy=0.0;
                  break;
                case TRUNCATED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  break;
                case SHIFTED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(rr>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                   SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                    (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                     SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                     SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                    energy=energy*SwitchingValue+TranslationValue;
                  }
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                         -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                         (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                         (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  break;
                case EWALD:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (erfc(Alpha[CurrentSystem]*r)/r);
                  break;
                default:
                  printf("Unknown charge-charge method in 'CalculateFrameworkCationsChargeChargeEnergyDifference'\n");
                  exit(0);
                  break;
              }
              UHostChargeChargeRealDelta[CurrentSystem]-=energy;
            }
          }
        }
      }
    }
  }
}


int CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference(int m)
{
  int j,k,f1;
  int A1,A2,B1,B2;
  int Type,TypeA;
  VECTOR new_posA,old_posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,rr,ri2,rk2,cosA,cosB,energy,temp,length,chargeB,chargeA;
  VECTOR new_dipoleA,old_dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL SwitchingValue;
  int ncell;

  UHostChargeBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  Bt0=Bt1=Bt2=0.0;  
  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][m].Atoms[j].Type;
      chargeA=Adsorbates[CurrentSystem][m].Atoms[j].Charge;
  
      new_posA=TrialPosition[CurrentSystem][j];
      old_posA=Adsorbates[CurrentSystem][m].Atoms[j].Position;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
  
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;


          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=new_posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=new_posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=new_posA.z-(posB.z+ReplicaShift[ncell].z);;
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
    
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
              UHostChargeBondDipoleRealDelta[CurrentSystem]+=energy;
            }
    
            dr.x=old_posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=old_posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=old_posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
    
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
              UHostChargeBondDipoleRealDelta[CurrentSystem]-=energy;
            }
          }
        }
      }
    }
  
    TypeA=Adsorbates[CurrentSystem][m].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
  
      posA1=TrialPosition[CurrentSystem][A1];
      posA2=TrialPosition[CurrentSystem][A2];
      new_dipoleA.x=posA2.x-posA1.x;
      new_dipoleA.y=posA2.y-posA1.y;
      new_dipoleA.z=posA2.z-posA1.z;
      new_posA.x=posA1.x+0.5*new_dipoleA.x;
      new_posA.y=posA1.y+0.5*new_dipoleA.y;
      new_posA.z=posA1.z+0.5*new_dipoleA.z;
      ri2=SQR(new_dipoleA.x)+SQR(new_dipoleA.y)+SQR(new_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      new_dipoleA.x*=temp; new_dipoleA.y*=temp; new_dipoleA.z*=temp;
  
      posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
      old_dipoleA.x=posA2.x-posA1.x;
      old_dipoleA.y=posA2.y-posA1.y;
      old_dipoleA.z=posA2.z-posA1.z;
      old_posA.x=posA1.x+0.5*old_dipoleA.x;
      old_posA.y=posA1.y+0.5*old_dipoleA.y;
      old_posA.z=posA1.z+0.5*old_dipoleA.z;
      ri2=SQR(old_dipoleA.x)+SQR(old_dipoleA.y)+SQR(old_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      old_dipoleA.x*=temp; old_dipoleA.y*=temp; old_dipoleA.z*=temp;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          Type=Framework[CurrentSystem].Atoms[f1][k].Type;
  
          if(PseudoAtoms[Type].HasCharges)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;
  
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=new_posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=new_posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=new_posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
              if(rr<CutOffChargeBondDipoleSquared)
              {
                r=sqrt(rr);
                SwitchingValue=1.0;
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt0=Bt1=Bt2=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt1=1.0/(r*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                    break;
                  case EWALD:
                    Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                    break;
                  default:
                    printf("Unknown charge-bonddipole method in 'CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference'\n");
                    exit(0);
                    break;
                }
    
                cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
                UHostChargeBondDipoleRealDelta[CurrentSystem]-=energy;
              }
    
              dr.x=old_posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=old_posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=old_posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
     
              if(rr<CutOffChargeBondDipoleSquared)
              {
                r=sqrt(rr);
                SwitchingValue=1.0;
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt0=Bt1=Bt2=0.0;  
                    break;
                  case TRUNCATED_COULOMB:
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt1=1.0/(r*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                    break;
                  case EWALD:
                    Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                    break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
    
                }
    
                cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
                UHostChargeBondDipoleRealDelta[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][m].Atoms[j].Type;
      chargeA=Adsorbates[CurrentSystem][m].Atoms[j].Charge;
  
      new_posA=TrialPosition[CurrentSystem][j];
      old_posA=Adsorbates[CurrentSystem][m].Atoms[j].Position;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
  
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
  
          dr.x=new_posA.x-posB.x;
          dr.y=new_posA.y-posB.y;
          dr.z=new_posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference'\n");
                exit(0);
                break;
            }
  
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
            UHostChargeBondDipoleRealDelta[CurrentSystem]+=energy;
          }
  
          dr.x=old_posA.x-posB.x;
          dr.y=old_posA.y-posB.y;
          dr.z=old_posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference'\n");
                exit(0);
                break;
            }
  
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
            UHostChargeBondDipoleRealDelta[CurrentSystem]-=energy;
          }
        }
      }
    }
  
    TypeA=Adsorbates[CurrentSystem][m].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
  
      posA1=TrialPosition[CurrentSystem][A1];
      posA2=TrialPosition[CurrentSystem][A2];
      new_dipoleA.x=posA2.x-posA1.x;
      new_dipoleA.y=posA2.y-posA1.y;
      new_dipoleA.z=posA2.z-posA1.z;
      new_posA.x=posA1.x+0.5*new_dipoleA.x;
      new_posA.y=posA1.y+0.5*new_dipoleA.y;
      new_posA.z=posA1.z+0.5*new_dipoleA.z;
      ri2=SQR(new_dipoleA.x)+SQR(new_dipoleA.y)+SQR(new_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      new_dipoleA.x*=temp; new_dipoleA.y*=temp; new_dipoleA.z*=temp;
  
      posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
      old_dipoleA.x=posA2.x-posA1.x;
      old_dipoleA.y=posA2.y-posA1.y;
      old_dipoleA.z=posA2.z-posA1.z;
      old_posA.x=posA1.x+0.5*old_dipoleA.x;
      old_posA.y=posA1.y+0.5*old_dipoleA.y;
      old_posA.z=posA1.z+0.5*old_dipoleA.z;
      ri2=SQR(old_dipoleA.x)+SQR(old_dipoleA.y)+SQR(old_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      old_dipoleA.x*=temp; old_dipoleA.y*=temp; old_dipoleA.z*=temp;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          Type=Framework[CurrentSystem].Atoms[f1][k].Type;
  
          if(PseudoAtoms[Type].HasCharges)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;
  
            dr.x=new_posA.x-posB.x;
            dr.y=new_posA.y-posB.y;
            dr.z=new_posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
  
              cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
              UHostChargeBondDipoleRealDelta[CurrentSystem]-=energy;
            }
  
    
            dr.x=old_posA.x-posB.x;
            dr.y=old_posA.y-posB.y;
            dr.z=old_posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
   
            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkAdsorbateChargeBondDipoleEnergyDifference'\n");
                exit(0);
                break;
  
              }
  
              cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
              UHostChargeBondDipoleRealDelta[CurrentSystem]+=energy;
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkCationChargeBondDipoleEnergyDifference(int m)
{
  int j,k,f1;
  int A1,A2,B1,B2;
  int Type,TypeA;
  VECTOR new_posA,old_posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,rr,ri2,rk2,cosA,cosB,energy,temp,length,chargeB,chargeA;
  VECTOR new_dipoleA,old_dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL SwitchingValue;
  int ncell;

  UHostChargeBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  Bt0=Bt1=Bt2=0.0;  
  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][m].Atoms[j].Type;
      chargeA=Cations[CurrentSystem][m].Atoms[j].Charge;
  
      new_posA=TrialPosition[CurrentSystem][j];
      old_posA=Cations[CurrentSystem][m].Atoms[j].Position;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
  
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
  
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=new_posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=new_posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=new_posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkCationChargeBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
    
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
              UHostChargeBondDipoleRealDelta[CurrentSystem]+=energy;
            }
    
            dr.x=old_posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=old_posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=old_posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkCationChargeBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
    
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
              UHostChargeBondDipoleRealDelta[CurrentSystem]-=energy;
            }
          }
        }
      }
    }
  
    TypeA=Cations[CurrentSystem][m].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
  
      posA1=TrialPosition[CurrentSystem][A1];
      posA2=TrialPosition[CurrentSystem][A2];
      new_dipoleA.x=posA2.x-posA1.x;
      new_dipoleA.y=posA2.y-posA1.y;
      new_dipoleA.z=posA2.z-posA1.z;
      new_posA.x=posA1.x+0.5*new_dipoleA.x;
      new_posA.y=posA1.y+0.5*new_dipoleA.y;
      new_posA.z=posA1.z+0.5*new_dipoleA.z;
      ri2=SQR(new_dipoleA.x)+SQR(new_dipoleA.y)+SQR(new_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      new_dipoleA.x*=temp; new_dipoleA.y*=temp; new_dipoleA.z*=temp;
  
      posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
      old_dipoleA.x=posA2.x-posA1.x;
      old_dipoleA.y=posA2.y-posA1.y;
      old_dipoleA.z=posA2.z-posA1.z;
      old_posA.x=posA1.x+0.5*old_dipoleA.x;
      old_posA.y=posA1.y+0.5*old_dipoleA.y;
      old_posA.z=posA1.z+0.5*old_dipoleA.z;
      ri2=SQR(old_dipoleA.x)+SQR(old_dipoleA.y)+SQR(old_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      old_dipoleA.x*=temp; old_dipoleA.y*=temp; old_dipoleA.z*=temp;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          Type=Framework[CurrentSystem].Atoms[f1][k].Type;
  
          if(PseudoAtoms[Type].HasCharges)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;
  
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=new_posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=new_posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=new_posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
              if(rr<CutOffChargeBondDipoleSquared)
              {
                r=sqrt(rr);
                SwitchingValue=1.0;
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt0=Bt1=Bt2=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt1=1.0/(r*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                    break;
                  case EWALD:
                    Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                    break;
                  default:
                    printf("Unknown charge-bonddipole method in 'CalculateFrameworkCationChargeBondDipoleEnergyDifference'\n");
                    exit(0);
                    break;
                }
    
                cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
                UHostChargeBondDipoleRealDelta[CurrentSystem]-=energy;
              }
    
      
              dr.x=old_posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=old_posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=old_posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
     
              if(rr<CutOffChargeBondDipoleSquared)
              {
                r=sqrt(rr);
                SwitchingValue=1.0;
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt0=Bt1=Bt2=0.0;  
                    break;
                  case TRUNCATED_COULOMB:
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt1=1.0/(r*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                    break;
                  case EWALD:
                    Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                    break;
                  default:
                    printf("Unknown charge-bonddipole method in 'CalculateFrameworkCationChargeBondDipoleEnergyDifference'\n");
                    exit(0);
                    break;
                }
    
                cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
                UHostChargeBondDipoleRealDelta[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][m].Atoms[j].Type;
      chargeA=Cations[CurrentSystem][m].Atoms[j].Charge;
  
      new_posA=TrialPosition[CurrentSystem][j];
      old_posA=Cations[CurrentSystem][m].Atoms[j].Position;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
  
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
  
          dr.x=new_posA.x-posB.x;
          dr.y=new_posA.y-posB.y;
          dr.z=new_posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkCationChargeBondDipoleEnergyDifference'\n");
                exit(0);
                break;
            }
  
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
            UHostChargeBondDipoleRealDelta[CurrentSystem]+=energy;
          }
  
          dr.x=old_posA.x-posB.x;
          dr.y=old_posA.y-posB.y;
          dr.z=old_posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
          if(rr<CutOffChargeBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown charge-bonddipole method in 'CalculateFrameworkCationChargeBondDipoleEnergyDifference'\n");
                exit(0);
                break;
            }
  
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
            UHostChargeBondDipoleRealDelta[CurrentSystem]-=energy;
          }
        }
      }
    }
  
    TypeA=Cations[CurrentSystem][m].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
  
      posA1=TrialPosition[CurrentSystem][A1];
      posA2=TrialPosition[CurrentSystem][A2];
      new_dipoleA.x=posA2.x-posA1.x;
      new_dipoleA.y=posA2.y-posA1.y;
      new_dipoleA.z=posA2.z-posA1.z;
      new_posA.x=posA1.x+0.5*new_dipoleA.x;
      new_posA.y=posA1.y+0.5*new_dipoleA.y;
      new_posA.z=posA1.z+0.5*new_dipoleA.z;
      ri2=SQR(new_dipoleA.x)+SQR(new_dipoleA.y)+SQR(new_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      new_dipoleA.x*=temp; new_dipoleA.y*=temp; new_dipoleA.z*=temp;
  
      posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
      old_dipoleA.x=posA2.x-posA1.x;
      old_dipoleA.y=posA2.y-posA1.y;
      old_dipoleA.z=posA2.z-posA1.z;
      old_posA.x=posA1.x+0.5*old_dipoleA.x;
      old_posA.y=posA1.y+0.5*old_dipoleA.y;
      old_posA.z=posA1.z+0.5*old_dipoleA.z;
      ri2=SQR(old_dipoleA.x)+SQR(old_dipoleA.y)+SQR(old_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      old_dipoleA.x*=temp; old_dipoleA.y*=temp; old_dipoleA.z*=temp;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
        {
          Type=Framework[CurrentSystem].Atoms[f1][k].Type;
  
          if(PseudoAtoms[Type].HasCharges)
          {
            posB=Framework[CurrentSystem].Atoms[f1][k].Position;
            chargeB=Framework[CurrentSystem].Atoms[f1][k].Charge;
  
            dr.x=new_posA.x-posB.x;
            dr.y=new_posA.y-posB.y;
            dr.z=new_posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkCationChargeBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
  
              cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
              UHostChargeBondDipoleRealDelta[CurrentSystem]-=energy;
            }
  
    
            dr.x=old_posA.x-posB.x;
            dr.y=old_posA.y-posB.y;
            dr.z=old_posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
   
            if(rr<CutOffChargeBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown charge-bonddipole method in 'CalculateFrameworkCationChargeBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
  
              cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeB*cosA);
              UHostChargeBondDipoleRealDelta[CurrentSystem]+=energy;
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference(int m)
{
  int j,k,f1;
  int A1,A2,B1,B2;
  int TypeA;
  REAL r,energy;
  VECTOR new_posA,old_posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL ri2,rk2,rr,length,temp;
  VECTOR new_dipoleA,old_dipoleA,dipoleB;
  REAL cosAB,cosA,cosB;
  REAL Bt0,Bt1,Bt2;
  REAL SwitchingValue;
  int ncell;

  OVERLAP=FALSE;
  UHostBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  TypeA=Adsorbates[CurrentSystem][m].Type;

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
  
      posA1=TrialPosition[CurrentSystem][A1];
      posA2=TrialPosition[CurrentSystem][A2];
      new_dipoleA.x=posA2.x-posA1.x;
      new_dipoleA.y=posA2.y-posA1.y;
      new_dipoleA.z=posA2.z-posA1.z;
      new_posA.x=posA1.x+0.5*new_dipoleA.x;
      new_posA.y=posA1.y+0.5*new_dipoleA.y;
      new_posA.z=posA1.z+0.5*new_dipoleA.z;
      ri2=SQR(new_dipoleA.x)+SQR(new_dipoleA.y)+SQR(new_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      new_dipoleA.x*=temp; new_dipoleA.y*=temp; new_dipoleA.z*=temp;
  
  
      posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
      old_dipoleA.x=posA2.x-posA1.x;
      old_dipoleA.y=posA2.y-posA1.y;
      old_dipoleA.z=posA2.z-posA1.z;
      old_posA.x=posA1.x+0.5*old_dipoleA.x;
      old_posA.y=posA1.y+0.5*old_dipoleA.y;
      old_posA.z=posA1.z+0.5*old_dipoleA.z;
      ri2=SQR(old_dipoleA.x)+SQR(old_dipoleA.y)+SQR(old_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      old_dipoleA.x*=temp; old_dipoleA.y*=temp; old_dipoleA.z*=temp;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
  
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=new_posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=new_posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=new_posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                  }
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
    
              cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
              cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              if(energy>=EnergyOverlapCriteria)
                return OVERLAP=TRUE; 
              UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
            }
    
            dr.x=old_posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=old_posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=old_posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                  }
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
    
              cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
              cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
  
      posA1=TrialPosition[CurrentSystem][A1];
      posA2=TrialPosition[CurrentSystem][A2];
      new_dipoleA.x=posA2.x-posA1.x;
      new_dipoleA.y=posA2.y-posA1.y;
      new_dipoleA.z=posA2.z-posA1.z;
      new_posA.x=posA1.x+0.5*new_dipoleA.x;
      new_posA.y=posA1.y+0.5*new_dipoleA.y;
      new_posA.z=posA1.z+0.5*new_dipoleA.z;
      ri2=SQR(new_dipoleA.x)+SQR(new_dipoleA.y)+SQR(new_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      new_dipoleA.x*=temp; new_dipoleA.y*=temp; new_dipoleA.z*=temp;
  
  
      posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
      old_dipoleA.x=posA2.x-posA1.x;
      old_dipoleA.y=posA2.y-posA1.y;
      old_dipoleA.z=posA2.z-posA1.z;
      old_posA.x=posA1.x+0.5*old_dipoleA.x;
      old_posA.y=posA1.y+0.5*old_dipoleA.y;
      old_posA.z=posA1.z+0.5*old_dipoleA.z;
      ri2=SQR(old_dipoleA.x)+SQR(old_dipoleA.y)+SQR(old_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      old_dipoleA.x*=temp; old_dipoleA.y*=temp; old_dipoleA.z*=temp;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
  
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
  
          dr.x=new_posA.x-posB.x;
          dr.y=new_posA.y-posB.y;
          dr.z=new_posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
          if(rr<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference'\n");
                exit(0);
                break;
            }
  
            cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
            cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            if(energy>=EnergyOverlapCriteria)
              return OVERLAP=TRUE; 
            UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
          }
  
          dr.x=old_posA.x-posB.x;
          dr.y=old_posA.y-posB.y;
          dr.z=old_posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
          if(rr<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkAdsorbateBondDipoleBondDipoleEnergyDifference'\n");
                exit(0);
                break;
            }
  
            cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
            cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
          }
        }
      }
    }
  }
  return 0;
}


int CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference(int m)
{
  int j,k,f1;
  int A1,A2,B1,B2;
  int TypeA;
  REAL r,energy;
  VECTOR new_posA,old_posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL ri2,rk2,rr,length,temp;
  VECTOR new_dipoleA,old_dipoleA,dipoleB;
  REAL cosAB,cosA,cosB;
  REAL Bt0,Bt1,Bt2;
  REAL SwitchingValue;
  int ncell;

  UHostBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  TypeA=Cations[CurrentSystem][m].Type;

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
  
      posA1=TrialPosition[CurrentSystem][A1];
      posA2=TrialPosition[CurrentSystem][A2];
      new_dipoleA.x=posA2.x-posA1.x;
      new_dipoleA.y=posA2.y-posA1.y;
      new_dipoleA.z=posA2.z-posA1.z;
      new_posA.x=posA1.x+0.5*new_dipoleA.x;
      new_posA.y=posA1.y+0.5*new_dipoleA.y;
      new_posA.z=posA1.z+0.5*new_dipoleA.z;
      ri2=SQR(new_dipoleA.x)+SQR(new_dipoleA.y)+SQR(new_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      new_dipoleA.x*=temp; new_dipoleA.y*=temp; new_dipoleA.z*=temp;
  
  
      posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
      old_dipoleA.x=posA2.x-posA1.x;
      old_dipoleA.y=posA2.y-posA1.y;
      old_dipoleA.z=posA2.z-posA1.z;
      old_posA.x=posA1.x+0.5*old_dipoleA.x;
      old_posA.y=posA1.y+0.5*old_dipoleA.y;
      old_posA.z=posA1.z+0.5*old_dipoleA.z;
      ri2=SQR(old_dipoleA.x)+SQR(old_dipoleA.y)+SQR(old_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      old_dipoleA.x*=temp; old_dipoleA.y*=temp; old_dipoleA.z*=temp;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
  
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
  
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=new_posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=new_posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=new_posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                  }
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
    
              cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
              cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
            }
    
            dr.x=old_posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=old_posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=old_posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              r=sqrt(rr);
              SwitchingValue=1.0;
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                  }
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  break;
                default:
                  printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference'\n");
                  exit(0);
                  break;
              }
    
              cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
              cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
  
      posA1=TrialPosition[CurrentSystem][A1];
      posA2=TrialPosition[CurrentSystem][A2];
      new_dipoleA.x=posA2.x-posA1.x;
      new_dipoleA.y=posA2.y-posA1.y;
      new_dipoleA.z=posA2.z-posA1.z;
      new_posA.x=posA1.x+0.5*new_dipoleA.x;
      new_posA.y=posA1.y+0.5*new_dipoleA.y;
      new_posA.z=posA1.z+0.5*new_dipoleA.z;
      ri2=SQR(new_dipoleA.x)+SQR(new_dipoleA.y)+SQR(new_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      new_dipoleA.x*=temp; new_dipoleA.y*=temp; new_dipoleA.z*=temp;
  
  
      posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
      old_dipoleA.x=posA2.x-posA1.x;
      old_dipoleA.y=posA2.y-posA1.y;
      old_dipoleA.z=posA2.z-posA1.z;
      old_posA.x=posA1.x+0.5*old_dipoleA.x;
      old_posA.y=posA1.y+0.5*old_dipoleA.y;
      old_posA.z=posA1.z+0.5*old_dipoleA.z;
      ri2=SQR(old_dipoleA.x)+SQR(old_dipoleA.y)+SQR(old_dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      old_dipoleA.x*=temp; old_dipoleA.y*=temp; old_dipoleA.z*=temp;
  
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(k=0;k<Framework[CurrentSystem].NumberOfBondDipoles[f1];k++)
        {
          DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][k];
          B1=Framework[CurrentSystem].BondDipoles[f1][k].A;
          B2=Framework[CurrentSystem].BondDipoles[f1][k].B;
          posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
          posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
  
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          dipoleB=ApplyBoundaryCondition(dipoleB);
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
  
          dr.x=new_posA.x-posB.x;
          dr.y=new_posA.y-posB.y;
          dr.z=new_posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
          if(rr<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference'\n");
                exit(0);
                break;
            }
  
            cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
            cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            UHostBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
          }
  
          dr.x=old_posA.x-posB.x;
          dr.y=old_posA.y-posB.y;
          dr.z=old_posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
          if(rr<CutOffBondDipoleBondDipoleSquared)
          {
            r=sqrt(rr);
            SwitchingValue=1.0;
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                break;
              default:
                printf("Unknown bonddipole-bonddipole method in 'CalculateFrameworkCationBondDipoleBondDipoleEnergyDifference'\n");
                exit(0);
                break;
            }
  
            cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
            cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            UHostBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
          }
        }
      }
    }
  }
  return 0;
}

REAL CalculateFrameworkElectrostaticPotential(POINT posA)
{
  int i,typeB,f;
  REAL r,rr,chargeB;
  VECTOR posB,dr;
  REAL ElectrostaticPotential,energy;
  REAL SwitchingValue,TranslationValue;

  ElectrostaticPotential=0.0;

  if(Framework[CurrentSystem].FrameworkModel!=NONE)
  {
    if(ChargeMethod!=NONE)
    {
      for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
      {
        for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f];i++)
        {
          typeB=Framework[CurrentSystem].Atoms[f][i].Type;
          posB=Framework[CurrentSystem].Atoms[f][i].Position;
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared)
          {
            r=sqrt(rr);
            chargeB=Framework[CurrentSystem].Atoms[f][i].Charge;

            switch(ChargeMethod)
            {
              case NONE:
                break;
              case TRUNCATED_COULOMB:
                ElectrostaticPotential+=COULOMBIC_CONVERSION_FACTOR*chargeB/r;
                break;
              case SHIFTED_COULOMB:
                ElectrostaticPotential+=COULOMBIC_CONVERSION_FACTOR*chargeB*(1.0/r-InverseCutOffChargeCharge);
                break;
              case SMOOTHED_COULOMB:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                if(rr>CutOffChargeChargeSwitchSquared)
                {
                  SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                 SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                  TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeB*
                                   (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                    SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                    SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                  energy=energy*SwitchingValue+TranslationValue;
                }
                ElectrostaticPotential+=energy;
                break;
              case WOLFS_METHOD_DAMPED_FG:
                ElectrostaticPotential+=COULOMBIC_CONVERSION_FACTOR*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                       -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                       (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                       (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                break;
              case EWALD:
                ElectrostaticPotential+=COULOMBIC_CONVERSION_FACTOR*chargeB*
                             (erfc(Alpha[CurrentSystem]*r)/r);
                break;
              default:
                printf("Unknown charge-charge method in 'CalculateFrameworkElectrostaticPotential'\n");
                exit(0);
                break;
            }
          }
        }
      }
    }
  }
  return ElectrostaticPotential;
}

REAL CalculateFrameworkVDWEnergyCorrection(VECTOR* Positions,VECTOR *AnisotropicPositions)
{
  int k;
  int typeA;
  REAL UVDWDelta;
  VECTOR posA;

  UVDWDelta=0.0;
  for(k=0;k<Components[CurrentComponent].NumberOfAtoms;k++)
  {
    typeA=Components[CurrentComponent].Type[k];

    if(PseudoAtoms[typeA].AnisotropicCorrection)
    {
      posA=AnisotropicPositions[k];
      UVDWDelta+=CalculateFrameworkVDWEnergyAtPosition(posA,typeA);

      posA=Positions[k];
      UVDWDelta-=CalculateFrameworkVDWEnergyAtPosition(posA,typeA);
    }
  }
  return UVDWDelta;
}

