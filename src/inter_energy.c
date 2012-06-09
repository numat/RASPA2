/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'inter_energy.c' is part of RASPA.

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "potentials.h"
#include "molecule.h"
#include "framework_energy.h"
#include "simulation.h"
#include "cbmc.h"
#include "ewald.h"
#include "utils.h"
#include "mc_moves.h"
#include "inter_energy.h"


REAL CalculateInterVDWEnergyAdsorbateAtPosition(POINT posA,int typeA,int exclude)
{
  int i,j,typeB;
  POINT posB;
  REAL r2,Uvdw;
  VECTOR dr;
  int TypeMolB;
  int ncell;

  if(OmitInterMolecularInteractions) return 0.0;

  Uvdw=0.0; 

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        TypeMolB=Adsorbates[CurrentSystem][j].Type;
        for(i=0;i<Components[Adsorbates[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Adsorbates[CurrentSystem][j].Atoms[i].AnisotropicPosition;
          typeB=Adsorbates[CurrentSystem][j].Atoms[i].Type;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
  
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
   
            if(r2<CutOffVDWSquared)
              Uvdw+=PotentialValue(typeA,typeB,r2);
          } 
        }
      }
    }
  }
  else
  {
    for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        TypeMolB=Adsorbates[CurrentSystem][j].Type;
        for(i=0;i<Components[Adsorbates[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Adsorbates[CurrentSystem][j].Atoms[i].AnisotropicPosition;
          typeB=Adsorbates[CurrentSystem][j].Atoms[i].Type;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);

          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(r2<CutOffVDWSquared)
            Uvdw+=PotentialValue(typeA,typeB,r2);
        }
      }
    }
  }

  return Uvdw;
}

int CalculateInterChargeEnergyAdsorbateAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,int exclude)
{
  int i,j,l,typeB,B1,B2,TypeB;
  POINT posB;
  REAL r,rr;
  REAL chargeA,chargeB;
  REAL Bt1;
  VECTOR dr,dipoleB;
  VECTOR posB1,posB2;
  REAL DipoleMagnitudeB,temp,cosB;
  REAL SwitchingValue,TranslationValue,energy;
  int ncell;

  *UChargeCharge=0.0;
  *UChargeBondDipole=0.0;

  if(OmitInterMolecularInteractions) return 0;

  chargeA=PseudoAtoms[typeA].Charge1;

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        for(i=0;i<Components[Adsorbates[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Adsorbates[CurrentSystem][j].Atoms[i].Position;
          typeB=Adsorbates[CurrentSystem][j].Atoms[i].Type;
          chargeB=Adsorbates[CurrentSystem][j].Atoms[i].Charge;

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
                case WOLFS_METHOD:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                         -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                          (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                          (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  break;
                case EWALD:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r);
                  break;
                default:
                  printf("Unknown charge-method in 'CalculateInterChargeEnergyAdsorbateAtPosition'\n");
                  exit(0);
                  break;
              }
              (*UChargeCharge)+=energy;
            }
          }
        }

        TypeB=Adsorbates[CurrentSystem][j].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Adsorbates[CurrentSystem][j].Atoms[B1].Position;
          posB2=Adsorbates[CurrentSystem][j].Atoms[B2].Position;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
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
                  Bt1=1.0/(r*rr);
                  break;
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
                  printf("Unknown charge-method in 'CalculateInterChargeEnergyAdsorbateAtPosition'\n");
                  exit(0);
                  break;
              }
 
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              (*UChargeBondDipole)+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        for(i=0;i<Components[Adsorbates[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Adsorbates[CurrentSystem][j].Atoms[i].Position;
          typeB=Adsorbates[CurrentSystem][j].Atoms[i].Type;
          chargeB=Adsorbates[CurrentSystem][j].Atoms[i].Charge;
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
              case WOLFS_METHOD:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                break;
              case WOLFS_METHOD_DAMPED_FG:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                             -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                              (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                              (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                break;
              case EWALD:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r);
                break;
              default:
                printf("Unknown charge-method in 'CalculateInterChargeEnergyAdsorbateAtPosition'\n");
                exit(0);
                break;
            }
            (*UChargeCharge)+=energy;
          }
        }

        TypeB=Adsorbates[CurrentSystem][j].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Adsorbates[CurrentSystem][j].Atoms[B1].Position;
          posB2=Adsorbates[CurrentSystem][j].Atoms[B2].Position;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
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
                Bt1=1.0/(r*rr);
                break;
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
                printf("Unknown charge-method in 'CalculateInterChargeEnergyAdsorbateAtPosition'\n");
                exit(0);
                break;
            }

            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            (*UChargeBondDipole)+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
          }
        }
      }
    }
  }
  return 0;
}

REAL CalculateInterVDWEnergyCationAtPosition(POINT posA,int typeA,int exclude)
{
  int i,j,typeB;
  POINT posB;
  REAL r2,Uvdw;
  VECTOR dr;
  int TypeMolB;
  int ncell;

  if(OmitInterMolecularInteractions) return 0.0;

  Uvdw=0.0; 
  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        TypeMolB=Cations[CurrentSystem][j].Type;
        for(i=0;i<Components[Cations[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Cations[CurrentSystem][j].Atoms[i].AnisotropicPosition;
          typeB=Cations[CurrentSystem][j].Atoms[i].Type;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
            if(r2<CutOffVDWSquared)
              Uvdw+=PotentialValue(typeA,typeB,r2);
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        TypeMolB=Cations[CurrentSystem][j].Type;
        for(i=0;i<Components[Cations[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Cations[CurrentSystem][j].Atoms[i].AnisotropicPosition;
          typeB=Cations[CurrentSystem][j].Atoms[i].Type;
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
  
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
          if(r2<CutOffVDWSquared)
            Uvdw+=PotentialValue(typeA,typeB,r2);
        }
      }
    }
  }
  return Uvdw;
}

int CalculateInterChargeEnergyCationAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,int exclude)
{
  int i,j,l,typeB,B1,B2,TypeB;
  POINT posB;
  REAL r,rr;
  REAL chargeA,chargeB;
  REAL Bt1;
  VECTOR dr,dipoleB;
  VECTOR posB1,posB2;
  REAL DipoleMagnitudeB,temp,cosB;
  REAL SwitchingValue,TranslationValue,energy;
  int ncell;

  *UChargeCharge=0.0;
  *UChargeBondDipole=0.0;

  if(OmitInterMolecularInteractions) return 0;

  chargeA=PseudoAtoms[typeA].Charge1;

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        for(i=0;i<Components[Cations[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Cations[CurrentSystem][j].Atoms[i].Position;
          typeB=Cations[CurrentSystem][j].Atoms[i].Type;
          chargeB=Cations[CurrentSystem][j].Atoms[i].Charge;

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
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r);
                  break;
                default:
                  printf("Unknown charge-method in 'CalculateInterChargeEnergyCationAtPosition'\n");
                  exit(0);
                  break;
              }
              (*UChargeCharge)+=energy;
            }
          }
        }
  
        TypeB=Cations[CurrentSystem][j].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Cations[CurrentSystem][j].Atoms[B1].Position;
          posB2=Cations[CurrentSystem][j].Atoms[B2].Position;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
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
                  Bt1=1.0/(r*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt1=1.0/(r*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt1=1.0/(r*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                  }
                  break;
                case EWALD:
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  break;
                default:
                  printf("Unknown charge-method in 'CalculateInterChargeEnergyCationAtPosition'\n");
                  exit(0);
                  break;
              }
    
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              (*UChargeBondDipole)+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        for(i=0;i<Components[Cations[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Cations[CurrentSystem][j].Atoms[i].Position;
          typeB=Cations[CurrentSystem][j].Atoms[i].Type;
          chargeB=Cations[CurrentSystem][j].Atoms[i].Charge;
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
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r);
                break;
              default:
                printf("Unknown charge-method in 'CalculateInterChargeEnergyCationAtPosition'\n");
                exit(0);
                break;
            }
            (*UChargeCharge)+=energy;
          }
        }
  
        TypeB=Cations[CurrentSystem][j].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Cations[CurrentSystem][j].Atoms[B1].Position;
          posB2=Cations[CurrentSystem][j].Atoms[B2].Position;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
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
                Bt1=1.0/(r*rr);
                break;
              case SHIFTED_COULOMB:
                Bt1=1.0/(r*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt1=1.0/(r*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                }
                break;
              case EWALD:
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                break;
              default:
                printf("Unknown charge-method in 'CalculateInterChargeEnergyCationAtPosition'\n");
                exit(0);
                break;
            }
  
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            (*UChargeBondDipole)+=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*chargeA*cosB);
          }
        }
      }
    }
  }
  return 0;
}

int CalculateInterBondDipoleEnergyAdsorbateAtPosition(POINT posA1,POINT posA2,REAL DipoleMagnitudeA,REAL *UChargeBondDipole,REAL *UBondDipoleBondDipole,int exclude)
{
  int i,j,l;
  int B1,B2;
  int TypeB;
  REAL r;
  VECTOR posA,posB,posB1,posB2,dr;
  REAL energy,DipoleMagnitudeB,ChargeB;
  REAL ri2,rr,length,temp;
  VECTOR dipoleA,dipoleB;
  REAL cosAB,cosA,cosB;
  REAL Bt0,Bt1,Bt2,SwitchingValue;
  int ncell;

 (*UChargeBondDipole)=0.0;
 (*UBondDipoleBondDipole)=0.0;

  if(OmitInterMolecularInteractions) return 0;

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

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        for(i=0;i<Components[Adsorbates[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Adsorbates[CurrentSystem][j].Atoms[i].Position;
          TypeB=Adsorbates[CurrentSystem][j].Atoms[i].Type;
          ChargeB=Adsorbates[CurrentSystem][j].Atoms[i].Charge;

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
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterBondDipoleEnergyAdsorbateAtPosition'\n");
                  exit(0);
                  break;
              }
    
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
              (*UChargeBondDipole)-=energy;
            }
          }
        }
  
        TypeB=Adsorbates[CurrentSystem][j].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Adsorbates[CurrentSystem][j].Atoms[B1].Position;
          posB2=Adsorbates[CurrentSystem][j].Atoms[B2].Position;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
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
                  printf("Unknown charge-method in 'CalculateInterBondDipoleEnergyAdsorbateAtPosition'\n");
                  exit(0);
                  break;
              }
    
              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              (*UBondDipoleBondDipole)+=energy;
    
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        for(i=0;i<Components[Adsorbates[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Adsorbates[CurrentSystem][j].Atoms[i].Position;
          TypeB=Adsorbates[CurrentSystem][j].Atoms[i].Type;
          ChargeB=Adsorbates[CurrentSystem][j].Atoms[i].Charge;
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
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                printf("Unknown charge-method in 'CalculateInterBondDipoleEnergyAdsorbateAtPosition'\n");
                exit(0);
                break;
            }
  
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
            (*UChargeBondDipole)-=energy;
          }
        }
  
        TypeB=Adsorbates[CurrentSystem][j].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Adsorbates[CurrentSystem][j].Atoms[B1].Position;
          posB2=Adsorbates[CurrentSystem][j].Atoms[B2].Position;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
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
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
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
                printf("Unknown charge-method in 'CalculateInterBondDipoleEnergyAdsorbateAtPosition'\n");
                exit(0);
                break;
            }
  
            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            (*UBondDipoleBondDipole)+=energy;
  
          }
        }
      }
    }
  }
  return 0;
}

int CalculateInterBondDipoleEnergyCationAtPosition(POINT posA1,POINT posA2,REAL DipoleMagnitudeA,REAL *UChargeBondDipole,REAL *UBondDipoleBondDipole,int exclude)
{
  int i,j,l;
  int B1,B2;
  int TypeB,typeB;
  REAL r;
  VECTOR posA,posB,posB1,posB2,dr;
  REAL energy,DipoleMagnitudeB,ChargeB;
  REAL ri2,rr,length,temp;
  VECTOR dipoleA,dipoleB;
  REAL cosAB,cosA,cosB;
  REAL Bt0,Bt1,Bt2,SwitchingValue;
  int ncell;

 (*UChargeBondDipole)=0.0;
 (*UBondDipoleBondDipole)=0.0;

  if(OmitInterMolecularInteractions) return 0;

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

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        for(i=0;i<Components[Cations[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Cations[CurrentSystem][j].Atoms[i].Position;
          typeB=Cations[CurrentSystem][j].Atoms[i].Type;
          ChargeB=Cations[CurrentSystem][j].Atoms[i].Charge;

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
                  Bt0=Bt1=Bt2=0.0;
                  break;
                case TRUNCATED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterBondDipoleEnergyCationAtPosition'\n");
                  exit(0);
                  break;
              }
    
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
              (*UChargeBondDipole)-=energy;
            }
          }
        }
  
        TypeB=Cations[CurrentSystem][j].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Cations[CurrentSystem][j].Atoms[B1].Position;
          posB2=Cations[CurrentSystem][j].Atoms[B2].Position;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
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
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterBondDipoleEnergyCationAtPosition'\n");
                  exit(0);
                  break;
              }
    
              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              (*UBondDipoleBondDipole)+=energy;
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
    {
      if(j!=exclude)
      {
        for(i=0;i<Components[Cations[CurrentSystem][j].Type].NumberOfAtoms;i++)
        {
          posB=Cations[CurrentSystem][j].Atoms[i].Position;
          typeB=Cations[CurrentSystem][j].Atoms[i].Type;
          ChargeB=Cations[CurrentSystem][j].Atoms[i].Charge;
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
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SMOOTHED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                if(rr>CutOffChargeBondDipoleSwitchSquared)
                {
                  SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                printf("Unknown charge-method in 'CalculateInterBondDipoleEnergyCationAtPosition'\n");
                exit(0);
                break;
            }
  
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
            (*UChargeBondDipole)-=energy;
          }
        }
  
        TypeB=Cations[CurrentSystem][j].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          posB1=Cations[CurrentSystem][j].Atoms[B1].Position;
          posB2=Cations[CurrentSystem][j].Atoms[B2].Position;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
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
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                printf("Unknown charge-method in 'CalculateInterBondDipoleEnergyCationAtPosition'\n");
                exit(0);
                break;
            }
  
            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            (*UBondDipoleBondDipole)+=energy;
          }
        }
      }
    }
  }
  return 0;
}

int CalculateTotalInterVDWEnergy(void)
{
  int i,j,k,l;
  int typeA,typeB,TypeMolA,TypeMolB;
  REAL rr,energy;
  VECTOR posA,posB,dr;
  int ncell,start;

  UAdsorbateAdsorbateVDW[CurrentSystem]=0.0;
  UAdsorbateCationVDW[CurrentSystem]=0.0;
  UCationCationVDW[CurrentSystem]=0.0;

  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      TypeMolA=Adsorbates[CurrentSystem][i].Type;
      for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
        posA=Adsorbates[CurrentSystem][i].Atoms[k].AnisotropicPosition;
  
        // loop over adsorbant molecules
        if(!OmitAdsorbateAdsorbateVDWInteractions)
        {
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            if(ncell==0) start=i+1;
            else start=0;
            for(j=start;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
            {
              TypeMolB=Adsorbates[CurrentSystem][j].Type;
              for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
              {
                typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
                posB=Adsorbates[CurrentSystem][j].Atoms[l].AnisotropicPosition;

                dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                {
                  energy=PotentialValue(typeA,typeB,rr);

                  // energy
                  if(ncell==0)
                    UAdsorbateAdsorbateVDW[CurrentSystem]+=energy;
                  else
                    UAdsorbateAdsorbateVDW[CurrentSystem]+=0.5*energy;
                }
              }
            }
          }
        }
  
        if(!OmitAdsorbateCationVDWInteractions)
        {
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
            {
              TypeMolB=Cations[CurrentSystem][j].Type;
              for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
              {
                posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
                typeB=Cations[CurrentSystem][j].Atoms[l].Type;

                dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                {
                  energy=PotentialValue(typeA,typeB,rr);

                  // energy
                  UAdsorbateCationVDW[CurrentSystem]+=energy;
                }
              }
            }
          }
        }
      }
    }
  
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      TypeMolA=Cations[CurrentSystem][i].Type;
      for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Cations[CurrentSystem][i].Atoms[k].Type;
        posA=Cations[CurrentSystem][i].Atoms[k].AnisotropicPosition;
  
        // loop over cation molecules
        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          if(ncell==0) start=i+1;
          else start=0;
          for(j=start;j<NumberOfCationMolecules[CurrentSystem];j++)
          {
            TypeMolB=Cations[CurrentSystem][j].Type;
            for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
            {
              posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
              typeB=Cations[CurrentSystem][j].Atoms[l].Type;

              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
              {
                typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                energy=PotentialValue(typeA,typeB,rr);

                // energy
                if(ncell==0)
                  UCationCationVDW[CurrentSystem]+=energy;
                else
                  UCationCationVDW[CurrentSystem]+=0.5*energy;
              }
            }
          }
        }
      }
    }
  }
  else
  { 
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      TypeMolA=Adsorbates[CurrentSystem][i].Type;
      for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
        posA=Adsorbates[CurrentSystem][i].Atoms[k].AnisotropicPosition;
  
        // loop over adsorbant molecules
        if(!OmitAdsorbateAdsorbateVDWInteractions)
        {
          for(j=i+1;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
          {
            TypeMolB=Adsorbates[CurrentSystem][j].Type;
            for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
            {
              typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
              posB=Adsorbates[CurrentSystem][j].Atoms[l].AnisotropicPosition;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(rr<CutOffVDWSquared)
                UAdsorbateAdsorbateVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);
            }
          }
        }
  
        if(!OmitAdsorbateCationVDWInteractions)
        {
          for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
          {
            TypeMolB=Cations[CurrentSystem][j].Type;
            for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
            {
              posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
              typeB=Cations[CurrentSystem][j].Atoms[l].Type;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UAdsorbateCationVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);
            }
          }
        }
      }
    }
  
    if(!OmitCationCationVDWInteractions)
    {
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        TypeMolA=Cations[CurrentSystem][i].Type;
        for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
        {
          typeA=Cations[CurrentSystem][i].Atoms[k].Type;
          posA=Cations[CurrentSystem][i].Atoms[k].AnisotropicPosition;
  
          // loop over cation molecules
          for(j=i+1;j<NumberOfCationMolecules[CurrentSystem];j++)
          {
            TypeMolB=Cations[CurrentSystem][j].Type;
            for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
            {
              posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
              typeB=Cations[CurrentSystem][j].Atoms[l].Type;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UCationCationVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr);
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateTotalInterChargeChargeCoulombEnergy(void)
{
  int i,j,k,l;
  int typeA,typeB;
  REAL r,r2;
  REAL chargeA,chargeB;
  REAL energy,SwitchingValue,TranslationValue;
  VECTOR posA,posB,dr;
  REAL NetChargeB,UWolfCorrection;
  int ncell,start;

  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeReal[CurrentSystem]=0.0;
  UCationCationChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
        posA=Adsorbates[CurrentSystem][i].Atoms[k].Position;
        chargeA=Adsorbates[CurrentSystem][i].Atoms[k].Charge;
  
        if(!OmitAdsorbateAdsorbateCoulombInteractions)
        {
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            if(ncell==0) start=i+1;
            else start=0;
            for(j=start;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
            {
              for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
              {
                posB=Adsorbates[CurrentSystem][j].Atoms[l].Position;

                dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(r2<CutOffChargeChargeSquared)
                {
                  r=sqrt(r2);
                  typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
                  chargeB=Adsorbates[CurrentSystem][j].Atoms[l].Charge;

                  switch(ChargeMethod)
                  {
                    case NONE:
                      energy=0.0;
                      break;
                    case SHIFTED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                      break;
                    case TRUNCATED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                      break;
                    case SMOOTHED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                      if(r2>CutOffChargeChargeSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                       SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                        TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                        (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                         SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                         SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                      printf("Unknown charge-method in 'CalculateTotalInterReplicaChargeChargeCoulombEnergy'\n");
                      exit(0);
                      break;
                  }

                  if(ncell==0)
                    UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=energy;
                  else
                    UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=0.5*energy;
                }
              }
            }
          }
        }
 
        if(!OmitAdsorbateCationCoulombInteractions)
        {
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++) 
          {
            for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
            {
              for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
              {
                posB=Cations[CurrentSystem][j].Atoms[l].Position;

                dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(r2<CutOffChargeChargeSquared)
                {
                  typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                  chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;

                  r=sqrt(r2);

                  switch(ChargeMethod)
                  {
                    case NONE:
                      energy=0.0;
                      break;
                    case SHIFTED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                      break;
                    case TRUNCATED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                      break;
                    case SMOOTHED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                      if(r2>CutOffChargeChargeSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                       SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                        TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                        (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                         SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                         SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                      printf("Unknown charge-method in 'CalculateTotalInterReplicaChargeChargeCoulombEnergy'\n");
                      exit(0);
                      break;
                  }

                  UAdsorbateCationChargeChargeReal[CurrentSystem]+=energy;
                }
              }
            }
          }
        }
      }
    }
  
    if(!OmitCationCationCoulombInteractions)
    {
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
        {
          typeA=Cations[CurrentSystem][i].Atoms[k].Type;
          posA=Cations[CurrentSystem][i].Atoms[k].Position;
          chargeA=Cations[CurrentSystem][i].Atoms[k].Charge;
  
          // loop over cation molecules
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            if(ncell==0) start=i+1;
            else start=0;
            for(j=start;j<NumberOfCationMolecules[CurrentSystem];j++)
            {
              for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
              {
                posB=Cations[CurrentSystem][j].Atoms[l].Position;

                dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(r2<CutOffChargeChargeSquared)
                {
                  typeB=Cations[CurrentSystem][j].Atoms[l].Type;

                  r=sqrt(r2);
                  chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;

                  switch(ChargeMethod)
                  {
                    case NONE:
                      energy=0.0;
                      break;
                    case SHIFTED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                      break;
                    case TRUNCATED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                      break;
                    case SMOOTHED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                      if(r2>CutOffChargeChargeSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                       SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

                        TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                        (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                         SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                         SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                      printf("Unknown charge-method in 'CalculateTotalInterReplicaChargeChargeCoulombEnergy'\n");
                      exit(0);
                      break;
                  }

                  if(ncell==0)
                    UCationCationChargeChargeReal[CurrentSystem]+=energy;
                  else
                    UCationCationChargeChargeReal[CurrentSystem]+=0.5*energy;
                }
              }
            }
          }
        }
      }
    }
  }
  else
  {
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
        posA=Adsorbates[CurrentSystem][i].Atoms[k].Position;
        chargeA=Adsorbates[CurrentSystem][i].Atoms[k].Charge;
  
        // loop over adsorbant molecules
        if(!OmitAdsorbateAdsorbateCoulombInteractions)
        {
          for(j=i+1;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
          {
            for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
            {
              posB=Adsorbates[CurrentSystem][j].Atoms[l].Position;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(r2<CutOffChargeChargeSquared)
              {
                typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
                chargeB=Adsorbates[CurrentSystem][j].Atoms[l].Charge;
  
                r=sqrt(r2);
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
                    if(r2>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                     SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
   
                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                        SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                        SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                      energy=energy*SwitchingValue+TranslationValue;
                    }
                    break;
                  case WOLFS_METHOD:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                    break;
                  case WOLFS_METHOD_DAMPED:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
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
                    printf("Unknown charge-method in 'CalculateTotalInterChargeChargeCoulombEnergy'\n");
                    exit(0);
                    break;
                }
                UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=energy;
              }
            }
          }
        }
  
        if(!OmitAdsorbateCationCoulombInteractions)
        {
          for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
          {
            for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
            {
              posB=Cations[CurrentSystem][j].Atoms[l].Position;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
              {
                typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;
  
                r=sqrt(r2);
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
                    if(r2>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                     SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                        SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                        SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                    printf("Unknown charge-method in 'CalculateTotalInterChargeChargeCoulombEnergy'\n");
                    exit(0);
                    break;
  
                }
                UAdsorbateCationChargeChargeReal[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  
    if(!OmitCationCationCoulombInteractions)
    {
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
        {
          typeA=Cations[CurrentSystem][i].Atoms[k].Type;
          posA=Cations[CurrentSystem][i].Atoms[k].Position;
          chargeA=Cations[CurrentSystem][i].Atoms[k].Charge;
  
          // loop over cation molecules
          for(j=i+1;j<NumberOfCationMolecules[CurrentSystem];j++)
          {
            for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
            {
              posB=Cations[CurrentSystem][j].Atoms[l].Position;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
              {
                typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;
  
                r=sqrt(r2);
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
                    if(r2>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                     SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                        SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                        SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                    printf("Unknown charge-method in 'CalculateTotalInterChargeChargeCoulombEnergy'\n");
                    exit(0);
                    break;
  
                }
                UCationCationChargeChargeReal[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }

  switch(ChargeMethod)
  {
    case WOLFS_METHOD:
      UWolfCorrection=0.0;
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=0.5*COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*InverseCutOffChargeCharge;
      }
      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UWolfCorrection;

      UWolfCorrection=0.0;
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Cations[CurrentSystem][i].Atoms[j].Type;
          chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=0.5*COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*InverseCutOffChargeCharge;
      }
      UCationCationChargeChargeReal[CurrentSystem]+=UWolfCorrection;
      break;
    case WOLFS_METHOD_DAMPED:
    case WOLFS_METHOD_DAMPED_FG:
      UWolfCorrection=0.0;
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*
             (0.5*erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
              Alpha[CurrentSystem]*M_1_SQRTPI);
      }
      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UWolfCorrection;

      UWolfCorrection=0.0;
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Cations[CurrentSystem][i].Atoms[j].Type;
          chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*
             (0.5*erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
              Alpha[CurrentSystem]*M_1_SQRTPI);
      }
      UCationCationChargeChargeReal[CurrentSystem]+=UWolfCorrection;
      break;
    default:
      break;
  }

  return 0;
}

int CalculateTotalInterChargeBondDipoleCoulombEnergy(void)
{
  int i,j,k,l;
  int A1,A2;
  int Type,TypeA;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL r,rr,ri2,cosA,energy,temp,length,ChargeB;
  VECTOR dipoleA;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,SwitchingValue;
  int start,ncell;

  UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleReal[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  Bt0=Bt1=Bt2=0.0;  
  if(UseReplicas[CurrentSystem])
  {
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
  
        if(!OmitAdsorbateAdsorbateCoulombInteractions)
        {
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            if(ncell==0) start=i+1;
            else start=0;
            for(k=start;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
            {
              if(i!=k)
              {
                for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
                {
                  Type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
                  if(PseudoAtoms[Type].HasCharges)
                  {
                    posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                    ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
    
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
                          Bt0=Bt1=Bt2=0.0;  
                          break;
                        case TRUNCATED_COULOMB:
                          Bt0=1.0/(r);
                          Bt1=1.0/(r*rr);
                          Bt2=3.0/(r*rr*rr);
                          break;
                        case SHIFTED_COULOMB:
                          Bt0=1.0/(r);
                          Bt1=1.0/(r*rr);
                          Bt2=3.0/(r*rr*rr);
                          break;
                        case SMOOTHED_COULOMB:
                          Bt0=1.0/(r);
                          Bt1=1.0/(r*rr);
                          Bt2=3.0/(r*rr*rr);
                          if(rr>CutOffChargeBondDipoleSwitchSquared)
                          {
                            SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                           SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                          printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombEnergy'\n");
                          exit(0);
                          break;
                      }
    
                      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                      energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                      if(ncell==0)
                        UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]-=energy;
                      else
                        UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]-=0.5*energy;
                    }
                  }
                }
              }
            }
          }
        }
 
        if(!OmitAdsorbateCationCoulombInteractions)
        {
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
            {
              for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
              {
                Type=Cations[CurrentSystem][k].Atoms[l].Type;
                if(PseudoAtoms[Type].HasCharges)
                {
                  posB=Cations[CurrentSystem][k].Atoms[l].Position;
                  ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;
      
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
                        Bt0=Bt1=Bt2=0.0;  
                        break;
                      case TRUNCATED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SHIFTED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SMOOTHED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        if(rr>CutOffChargeBondDipoleSwitchSquared)
                        {
                           SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                          SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                        printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombEnergy'\n");
                        exit(0);
                        break;
                    }
    
                    cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                    energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                    UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=energy;
                  }
                }
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
 
        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          if(ncell==0) start=i+1;
          else start=0;
          for(k=start;k<NumberOfCationMolecules[CurrentSystem];k++)
          {
            if(i!=k)
            {
              if(!OmitCationCationCoulombInteractions)
              {
                for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
                {
                  Type=Cations[CurrentSystem][k].Atoms[l].Type;
                  if(PseudoAtoms[Type].HasCharges)
                  {
                    posB=Cations[CurrentSystem][k].Atoms[l].Position;
                    ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;
      
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
                          Bt0=Bt1=Bt2=0.0;  
                          break;
                        case TRUNCATED_COULOMB:
                          Bt0=1.0/(r);
                          Bt1=1.0/(r*rr);
                          Bt2=3.0/(r*rr*rr);
                          break;
                        case SHIFTED_COULOMB:
                          Bt0=1.0/(r);
                          Bt1=1.0/(r*rr);
                          Bt2=3.0/(r*rr*rr);
                          break;
                        case SMOOTHED_COULOMB:
                          Bt0=1.0/(r);
                          Bt1=1.0/(r*rr);
                          Bt2=3.0/(r*rr*rr);
                          if(rr>CutOffChargeBondDipoleSwitchSquared)
                          {
                            SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                           SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                          printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombEnergy'\n");
                          exit(0);
                          break;
                      }
    
                      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                      energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                      if(ncell==0)
                        UCationCationChargeBondDipoleReal[CurrentSystem]-=energy;
                      else
                        UCationCationChargeBondDipoleReal[CurrentSystem]-=0.5*energy;
                    }
                  }
                }
              }
            }
          }
        }
  
        if(!OmitAdsorbateCationCoulombInteractions)
        {
          for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
          {
            for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
            {
              Type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[Type].HasCharges)
              {
                posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
    
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
                        Bt0=Bt1=Bt2=0.0;  
                        break;
                      case TRUNCATED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SHIFTED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SMOOTHED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        if(rr>CutOffChargeBondDipoleSwitchSquared)
                        {
                          SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                         SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                        printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombEnergy'\n");
                        exit(0);
                        break;
                    }
    
                    cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                    energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                    UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=energy;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  else
  {
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
  
        if(!OmitAdsorbateAdsorbateCoulombInteractions)
        {
          for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
          {
            if(i!=k)
            {
              for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
              {
                Type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
                if(PseudoAtoms[Type].HasCharges)
                {
                  posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                  ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
  
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
                        Bt0=Bt1=Bt2=0.0;  
                        break;
                      case TRUNCATED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SHIFTED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SMOOTHED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        if(rr>CutOffChargeBondDipoleSwitchSquared)
                        {
                          SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                         SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                        printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombEnergy'\n");
                        exit(0);
                        break;
                    }
  
                    cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                    energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                    UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]-=energy;
                  }
                }
              }
            }
          }
        }
  
        if(!OmitAdsorbateCationCoulombInteractions)
        {
          for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
          {
            for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
            {
              Type=Cations[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[Type].HasCharges)
              {
                posB=Cations[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;
    
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
                      Bt0=Bt1=Bt2=0.0;  
                      break;
                    case TRUNCATED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SMOOTHED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      if(rr>CutOffChargeBondDipoleSwitchSquared)
                      {
                         SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                        SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombEnergy'\n");
                      exit(0);
                      break;
                  }
  
                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                  UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=energy;
                }
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
  
        if(!OmitCationCationCoulombInteractions)
        {
          for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
          {
            if(i!=k)
            {
              for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
              {
                Type=Cations[CurrentSystem][k].Atoms[l].Type;
                if(PseudoAtoms[Type].HasCharges)
                {
                  posB=Cations[CurrentSystem][k].Atoms[l].Position;
                  ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;
    
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
                        Bt0=Bt1=Bt2=0.0;  
                        break;
                      case TRUNCATED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SHIFTED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SMOOTHED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        if(rr>CutOffChargeBondDipoleSwitchSquared)
                        {
                          SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                         SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                        printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombEnergy'\n");
                        exit(0);
                        break;
                    }
  
                    cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                    energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                    UCationCationChargeBondDipoleReal[CurrentSystem]-=energy;
                  }
                }
              }
            }
          }
        }
  
        if(!OmitAdsorbateCationCoulombInteractions)
        {
          for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
          {
            for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
            {
              Type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[Type].HasCharges)
              {
                posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
    
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
                      Bt0=Bt1=Bt2=0.0;  
                      break;
                    case TRUNCATED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SMOOTHED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      if(rr>CutOffChargeBondDipoleSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombEnergy'\n");
                      exit(0);
                      break;
                  }
  
                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                  UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=energy;
                }
              }
            }
          }
        }
      }
    }
  }
  return 0;
}


int CalculateTotalInterBondDipoleBondDipoleCoulombEnergy(void)
{
  int i,j,k,l;
  int A1,A2,B1,B2;
  int TypeA,TypeB;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,rr,ri2,rk2,cosAB,cosA,cosB,energy,temp,length;
  VECTOR dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2,SwitchingValue;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  int ncell,start;

  Bt0=Bt1=Bt2=0.0;  
  UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
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
  
        if(!OmitAdsorbateAdsorbateCoulombInteractions)
        {
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            if(ncell==0) start=i+1;
            else start=0;
            for(k=start;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
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
                        SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateTotalInterBondDipoleBondDipoleCoulombEnergy'\n");
                      exit(0);
                      break;
                  }
    
                  cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
                  if(ncell==0)
                    UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=energy;
                  else
                    UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=0.5*energy;
                }
              }
            }
          }
        }
  
        if(!OmitAdsorbateCationCoulombInteractions)
        {
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
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
                        SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateTotalInterBondDipoleBondDipoleCoulombEnergy'\n");
                      exit(0);
                      break;
                  }
    
                  cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
    
                  UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=energy;
                }
              }
            }
          }
        }
      }
    }
  
    if(!OmitCationCationCoulombInteractions)
    {
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
 
          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            if(ncell==0) start=i+1;
            else start=0;
            for(k=start;k<NumberOfCationMolecules[CurrentSystem];k++)
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
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
                        SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateTotalInterBondDipoleBondDipoleCoulombEnergy'\n");
                      exit(0);
                      break;
                  }
    
                  cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
   
                  if(ncell==0) 
                    UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=energy;
                  else
                    UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=0.5*energy;
                }
              }
            }
          }
        }
      }
    }
  }
  else
  {
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
  
        if(!OmitAdsorbateAdsorbateCoulombInteractions)
        {
          for(k=i+1;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
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
                      SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateTotalInterBondDipoleBondDipoleCoulombEnergy'\n");
                    exit(0);
                    break;
                }
  
                cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
                UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=energy;
              }
            }
          }
        }
  
        if(!OmitAdsorbateCationCoulombInteractions)
        {
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
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
                      SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateTotalInterBondDipoleBondDipoleCoulombEnergy'\n");
                    exit(0);
                    break;
                }
  
                cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
  
                UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  
    if(!OmitCationCationCoulombInteractions)
    {
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
  
          for(k=i+1;k<NumberOfCationMolecules[CurrentSystem];k++)
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
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
                      SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateTotalInterBondDipoleBondDipoleCoulombEnergy'\n");
                    exit(0);
                    break;
                }
  
                cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
  
                UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=energy;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

// HIER
int CalculateInterVDWEnergyDifferenceAdsorbate(int m)
{
  int i,j,k;
  int typeA,typeB;
  REAL rr,energy;
  VECTOR posA_new,posA_old,posB,dr;
  int TypeMolA,TypeMolB;
  int ncell;

  OVERLAP=FALSE;
  UAdsorbateVDWDelta[CurrentSystem]=0.0;
  UCationVDWDelta[CurrentSystem]=0.0;
  if(OmitInterMolecularInteractions) return 0;

  TypeMolA=Adsorbates[CurrentSystem][m].Type;
  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Adsorbates[CurrentSystem][m].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][m].Atoms[k].Type;
  
      posA_old=Adsorbates[CurrentSystem][m].Atoms[k].AnisotropicPosition;
      posA_new=TrialAnisotropicPosition[CurrentSystem][k];
  
      if(!OmitAdsorbateAdsorbateVDWInteractions)
      {
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            TypeMolB=Adsorbates[CurrentSystem][i].Type;
            for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
              typeB=Adsorbates[CurrentSystem][i].Atoms[j].Type;
 
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
                  UAdsorbateVDWDelta[CurrentSystem]+=energy;
                }
    
                dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                  UAdsorbateVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);
              }
            }
          }
        }
      }
  
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        TypeMolB=Cations[CurrentSystem][i].Type;
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          posB=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
          typeB=Cations[CurrentSystem][i].Atoms[j].Type;
  
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
              UCationVDWDelta[CurrentSystem]+=energy;
            }
    
            dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UCationVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);
            }
        }
      }
    }
  }
  else
  {
    for(k=0;k<Adsorbates[CurrentSystem][m].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][m].Atoms[k].Type;
  
      posA_old=Adsorbates[CurrentSystem][m].Atoms[k].AnisotropicPosition;
      posA_new=TrialAnisotropicPosition[CurrentSystem][k];
  
      if(!OmitAdsorbateAdsorbateVDWInteractions)
      {
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            TypeMolB=Adsorbates[CurrentSystem][i].Type;
            for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
              typeB=Adsorbates[CurrentSystem][i].Atoms[j].Type;
  
              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UAdsorbateVDWDelta[CurrentSystem]+=energy;
              }
  
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
  
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        TypeMolB=Cations[CurrentSystem][i].Type;
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          posB=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
          typeB=Cations[CurrentSystem][i].Atoms[j].Type;
  
          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
          {
            energy=PotentialValue(typeA,typeB,rr);
            if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
            UCationVDWDelta[CurrentSystem]+=energy;
          }
  
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
  return 0;
}

int CalculateInterVDWEnergyDifferenceCation(int m)
{
  int i,j,k;
  int typeA,typeB;
  REAL rr,energy;
  VECTOR posA_new,posA_old,posB,dr;
  int TypeMolA,TypeMolB;
  int ncell;

  OVERLAP=FALSE;
  UAdsorbateVDWDelta[CurrentSystem]=0.0;
  UCationVDWDelta[CurrentSystem]=0.0;
  if(OmitInterMolecularInteractions) return 0;

  TypeMolA=Cations[CurrentSystem][m].Type;
  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Cations[CurrentSystem][m].NumberOfAtoms;k++)
    {
      typeA=Cations[CurrentSystem][m].Atoms[k].Type;
      posA_old=Cations[CurrentSystem][m].Atoms[k].AnisotropicPosition;
      posA_new=TrialAnisotropicPosition[CurrentSystem][k];
  
      if(!OmitCationCationVDWInteractions)
      {
        for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            TypeMolB=Cations[CurrentSystem][i].Type;
            for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
              typeB=Cations[CurrentSystem][i].Atoms[j].Type;
 
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
                  UCationVDWDelta[CurrentSystem]+=energy;
                }
    
                dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                  UCationVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);
              }
            }
          }
        }
      }
  
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        TypeMolB=Adsorbates[CurrentSystem][i].Type;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          posB=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
          typeB=Adsorbates[CurrentSystem][i].Atoms[j].Type;
  
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
              UAdsorbateVDWDelta[CurrentSystem]+=energy;
            }
    
            dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UAdsorbateVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr);
          }
        }
      }
    }
  }
  else
  {
    for(k=0;k<Cations[CurrentSystem][m].NumberOfAtoms;k++)
    {
      typeA=Cations[CurrentSystem][m].Atoms[k].Type;
      posA_old=Cations[CurrentSystem][m].Atoms[k].AnisotropicPosition;
      posA_new=TrialAnisotropicPosition[CurrentSystem][k];
  
      if(!OmitCationCationVDWInteractions)
      {
        for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            TypeMolB=Cations[CurrentSystem][i].Type;
            for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
              typeB=Cations[CurrentSystem][i].Atoms[j].Type;
  
              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UCationVDWDelta[CurrentSystem]+=energy;
              }
  
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
  
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        TypeMolB=Adsorbates[CurrentSystem][i].Type;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          posB=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
          typeB=Adsorbates[CurrentSystem][i].Atoms[j].Type;
  
          dr.x=posA_new.x-posB.x;
          dr.y=posA_new.y-posB.y;
          dr.z=posA_new.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
          {
            energy=PotentialValue(typeA,typeB,rr);
            if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
            UAdsorbateVDWDelta[CurrentSystem]+=energy;
          }
  
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
  return 0;
}

int CalculateInterChargeChargeEnergyDifferenceAdsorbate(int m)
{
  int i,j,k;
  int typeA,typeB;
  REAL r,r2,chargeA,chargeB,energy;
  VECTOR posA,posB,dr;
  REAL SwitchingValue,TranslationValue;
  int ncell;

  UAdsorbateChargeChargeRealDelta[CurrentSystem]=0.0;
  UCationChargeChargeRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Adsorbates[CurrentSystem][m].NumberOfAtoms;k++)
    {
      posA=Adsorbates[CurrentSystem][m].Atoms[k].Position;
      typeA=Adsorbates[CurrentSystem][m].Atoms[k].Type;
      chargeA=Adsorbates[CurrentSystem][m].Atoms[k].Charge;
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Adsorbates[CurrentSystem][i].Atoms[j].Position;
              typeB=Adsorbates[CurrentSystem][i].Atoms[j].Type;
              chargeB=Adsorbates[CurrentSystem][i].Atoms[j].Charge;

              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=TrialPosition[CurrentSystem][k].x-(posB.x+ReplicaShift[ncell].x);
                dr.y=TrialPosition[CurrentSystem][k].y-(posB.y+ReplicaShift[ncell].y);
                dr.z=TrialPosition[CurrentSystem][k].z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(r2<CutOffChargeChargeSquared)
                {
                  r=sqrt(r2);
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
                    case WOLFS_METHOD_DAMPED:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                  -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                      break;
                    case SMOOTHED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                      if(r2>CutOffChargeChargeSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                       SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
    
                        TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                         (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                          SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                          SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                      printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceAdsorbate'\n");
                      exit(0);
                      break;
    
                  }
                  UAdsorbateChargeChargeRealDelta[CurrentSystem]+=energy;
                }
    
                dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(r2<CutOffChargeChargeSquared)
                {
                  r=sqrt(r2);
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
                    case WOLFS_METHOD_DAMPED:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                  -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                      break;
                    case SMOOTHED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                      if(r2>CutOffChargeChargeSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                       SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
    
                        TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                         (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                          SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                          SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                      printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceAdsorbate'\n");
                      exit(0);
                      break;
                  }
                  UAdsorbateChargeChargeRealDelta[CurrentSystem]-=energy;
                }
              }
            }
          }
        }
      }
  
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          posB=Cations[CurrentSystem][i].Atoms[j].Position;
          typeB=Cations[CurrentSystem][i].Atoms[j].Type;
          chargeB=Cations[CurrentSystem][i].Atoms[j].Charge;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=TrialPosition[CurrentSystem][k].x-(posB.x+ReplicaShift[ncell].x);
            dr.y=TrialPosition[CurrentSystem][k].y-(posB.y+ReplicaShift[ncell].y);
            dr.z=TrialPosition[CurrentSystem][k].z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(r2<CutOffChargeChargeSquared)
            {
              r=sqrt(r2);
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
                case WOLFS_METHOD_DAMPED:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                              -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(r2>CutOffChargeChargeSwitchSquared)
                  {
                   SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                   SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
    
                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                      (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                      SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                     SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                  printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceAdsorbate'\n");
                  exit(0);
                  break;
    
              }
              UCationChargeChargeRealDelta[CurrentSystem]+=energy;
            }
    
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(r2<CutOffChargeChargeSquared)
            {
              r=sqrt(r2);
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
                case WOLFS_METHOD_DAMPED:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                              -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(r2>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                   SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
    
                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                       SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                       SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                  printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceAdsorbate'\n");
                  exit(0);
                  break;
    
              }
              UCationChargeChargeRealDelta[CurrentSystem]-=energy;
            }
          }
        }
      }
    }
  }
  else
  {
    for(k=0;k<Adsorbates[CurrentSystem][m].NumberOfAtoms;k++)
    {
      posA=Adsorbates[CurrentSystem][m].Atoms[k].Position;
      typeA=Adsorbates[CurrentSystem][m].Atoms[k].Type;
      chargeA=Adsorbates[CurrentSystem][m].Atoms[k].Charge;
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Adsorbates[CurrentSystem][i].Atoms[j].Position;
              typeB=Adsorbates[CurrentSystem][i].Atoms[j].Type;
              chargeB=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
              dr.x=TrialPosition[CurrentSystem][k].x-posB.x;
              dr.y=TrialPosition[CurrentSystem][k].y-posB.y;
              dr.z=TrialPosition[CurrentSystem][k].z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
              {
                r=sqrt(r2);
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
                  case WOLFS_METHOD_DAMPED:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(r2>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                     SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                        SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                        SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*erfc(Alpha[CurrentSystem]*r)/r;
                    break;
                  default:
                    printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceAdsorbate'\n");
                    exit(0);
                    break;
  
                }
                UAdsorbateChargeChargeRealDelta[CurrentSystem]+=energy;
              }
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
              {
                r=sqrt(r2);
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
                  case WOLFS_METHOD_DAMPED:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(r2>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                     SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                        SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                        SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*erfc(Alpha[CurrentSystem]*r)/r;
                    break;
                  default:
                    printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceAdsorbate'\n");
                    exit(0);
                    break;
                }
                UAdsorbateChargeChargeRealDelta[CurrentSystem]-=energy;
              }
            }
          }
        }
      }
  
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          posB=Cations[CurrentSystem][i].Atoms[j].Position;
          typeB=Cations[CurrentSystem][i].Atoms[j].Type;
          chargeB=Cations[CurrentSystem][i].Atoms[j].Charge;
          dr.x=TrialPosition[CurrentSystem][k].x-posB.x;
          dr.y=TrialPosition[CurrentSystem][k].y-posB.y;
          dr.z=TrialPosition[CurrentSystem][k].z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(r2<CutOffChargeChargeSquared)
          {
            r=sqrt(r2);
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
              case WOLFS_METHOD_DAMPED:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                            -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                break;
              case SMOOTHED_COULOMB:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                if(r2>CutOffChargeChargeSwitchSquared)
                {
                 SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                 SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                  TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                    (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                    SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                   SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceAdsorbate'\n");
                exit(0);
                break;
  
            }
            UCationChargeChargeRealDelta[CurrentSystem]+=energy;
          }
  
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(r2<CutOffChargeChargeSquared)
          {
            r=sqrt(r2);
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
              case WOLFS_METHOD_DAMPED:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                            -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                break;
              case SMOOTHED_COULOMB:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                if(r2>CutOffChargeChargeSwitchSquared)
                {
                  SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                 SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                  TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                   (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                     SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                     SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceAdsorbate'\n");
                exit(0);
                break;
  
            }
            UCationChargeChargeRealDelta[CurrentSystem]-=energy;
          }
        }
      }
    }
  }
  return 0;
}

int CalculateInterChargeChargeEnergyDifferenceCation(int m)
{
  int i,j,k;
  int typeA,typeB;
  REAL r,r2,chargeA,chargeB,energy;
  VECTOR posA,posB,dr;
  REAL SwitchingValue,TranslationValue;
  int ncell;

  UAdsorbateChargeChargeRealDelta[CurrentSystem]=0.0;
  UCationChargeChargeRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Cations[CurrentSystem][m].NumberOfAtoms;k++)
    {
      posA=Cations[CurrentSystem][m].Atoms[k].Position;
      typeA=Cations[CurrentSystem][m].Atoms[k].Type;
      chargeA=Cations[CurrentSystem][m].Atoms[k].Charge;
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Cations[CurrentSystem][i].Atoms[j].Position;
              typeB=Cations[CurrentSystem][i].Atoms[j].Type;
              chargeB=Cations[CurrentSystem][i].Atoms[j].Charge;

              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              {
                dr.x=TrialPosition[CurrentSystem][k].x-(posB.x+ReplicaShift[ncell].x);
                dr.y=TrialPosition[CurrentSystem][k].y-(posB.y+ReplicaShift[ncell].y);
                dr.z=TrialPosition[CurrentSystem][k].z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(r2<CutOffChargeChargeSquared)
                {
                  r=sqrt(r2);
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
                    case WOLFS_METHOD_DAMPED:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                  -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                      break;
                    case SMOOTHED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                      if(r2>CutOffChargeChargeSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                       SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
    
                        TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                         (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                          SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                          SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                      printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceCation'\n");
                      exit(0);
                      break;
                  }
                  UCationChargeChargeRealDelta[CurrentSystem]+=energy;
                }

                dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(r2<CutOffChargeChargeSquared)
                {
                  r=sqrt(r2);
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
                    case WOLFS_METHOD_DAMPED:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                  -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                      break;
                    case SMOOTHED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                      if(r2>CutOffChargeChargeSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                       SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
    
                        TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                         (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                          SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                          SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                      printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceCation'\n");
                      exit(0);
                      break;
                  }
                  UCationChargeChargeRealDelta[CurrentSystem]-=energy;
                }
              }
            }
          }
        }
      }
  
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          posB=Adsorbates[CurrentSystem][i].Atoms[j].Position;
          typeB=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          chargeB=Adsorbates[CurrentSystem][i].Atoms[j].Charge;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=TrialPosition[CurrentSystem][k].x-(posB.x+ReplicaShift[ncell].x);
            dr.y=TrialPosition[CurrentSystem][k].y-(posB.y+ReplicaShift[ncell].y);
            dr.z=TrialPosition[CurrentSystem][k].z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(r2<CutOffChargeChargeSquared)
            {
              r=sqrt(r2);
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
                case WOLFS_METHOD_DAMPED:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                              -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(r2>CutOffChargeChargeSwitchSquared)
                  {
                   SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                  SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
    
                   TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                    (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                     SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                     SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                  printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceCation'\n");
                  exit(0);
                  break;
              }
              UAdsorbateChargeChargeRealDelta[CurrentSystem]+=energy;
            }
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(r2<CutOffChargeChargeSquared)
            {
              r=sqrt(r2);
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
                case WOLFS_METHOD_DAMPED:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                              -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                  break;
                case SMOOTHED_COULOMB:
                  energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  if(r2>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                   SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
    
                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                     (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                      SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                      SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                  printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceCation'\n");
                  exit(0);
                  break;
              } 
              UAdsorbateChargeChargeRealDelta[CurrentSystem]-=energy;
            }
          }
        }
      }
    }
  }
  else
  {
    for(k=0;k<Cations[CurrentSystem][m].NumberOfAtoms;k++)
    {
      posA=Cations[CurrentSystem][m].Atoms[k].Position;
      typeA=Cations[CurrentSystem][m].Atoms[k].Type;
      chargeA=Cations[CurrentSystem][m].Atoms[k].Charge;
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Cations[CurrentSystem][i].Atoms[j].Position;
              typeB=Cations[CurrentSystem][i].Atoms[j].Type;
              chargeB=Cations[CurrentSystem][i].Atoms[j].Charge;
              dr.x=TrialPosition[CurrentSystem][k].x-posB.x;
              dr.y=TrialPosition[CurrentSystem][k].y-posB.y;
              dr.z=TrialPosition[CurrentSystem][k].z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
              {
                r=sqrt(r2);
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
                  case WOLFS_METHOD_DAMPED:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(r2>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                     SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                        SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                        SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                    printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceCation'\n");
                    exit(0);
                    break;
                }
                UCationChargeChargeRealDelta[CurrentSystem]+=energy;
              }
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
              {
                r=sqrt(r2);
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
                  case WOLFS_METHOD_DAMPED:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                    break;
                  case SMOOTHED_COULOMB:
                    energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                    if(r2>CutOffChargeChargeSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                     SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                      TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                        SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                        SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                    printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceCation'\n");
                    exit(0);
                    break;
                }
                UCationChargeChargeRealDelta[CurrentSystem]-=energy;
              }
            }
          }
        }
      }
  
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          posB=Adsorbates[CurrentSystem][i].Atoms[j].Position;
          typeB=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          chargeB=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
          dr.x=TrialPosition[CurrentSystem][k].x-posB.x;
          dr.y=TrialPosition[CurrentSystem][k].y-posB.y;
          dr.z=TrialPosition[CurrentSystem][k].z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(r2<CutOffChargeChargeSquared)
          {
            r=sqrt(r2);
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
              case WOLFS_METHOD_DAMPED:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                            -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                break;
              case SMOOTHED_COULOMB:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                if(r2>CutOffChargeChargeSwitchSquared)
                {
                 SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                 TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                  (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                   SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                   SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceCation'\n");
                exit(0);
                break;
            }
            UAdsorbateChargeChargeRealDelta[CurrentSystem]+=energy;
          }
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(r2<CutOffChargeChargeSquared)
          {
            r=sqrt(r2);
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
              case WOLFS_METHOD_DAMPED:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                            -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
                break;
              case SMOOTHED_COULOMB:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                if(r2>CutOffChargeChargeSwitchSquared)
                {
                  SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                                 SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
  
                  TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                   (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                    SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                    SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
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
                printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceCation'\n");
                exit(0);
                break;
            } 
            UAdsorbateChargeChargeRealDelta[CurrentSystem]-=energy;
          }
        }
      }
    }
  }
  return 0;
}


int CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(int m)
{
  int j,k,l;
  int A1,A2,B1,B2;
  int Type,TypeA,TypeB;
  VECTOR new_posA,old_posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,rr,ri2,rk2,cosA,cosB,energy,temp,length,ChargeB,ChargeA;
  VECTOR new_dipoleA,old_dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL SwitchingValue;
  int ncell;

  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UCationChargeBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  Bt0=Bt1=Bt2=0.0;  
  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][m].Atoms[j].Type;
      ChargeA=Adsorbates[CurrentSystem][m].Atoms[j].Charge;
  
      new_posA=TrialPosition[CurrentSystem][j];
      old_posA=Adsorbates[CurrentSystem][m].Atoms[j].Position;
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            TypeB=Adsorbates[CurrentSystem][k].Type;
            for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
            {
              B1=Components[TypeB].BondDipoles[l].A;
              B2=Components[TypeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
  
              posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
              posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      if(rr>CutOffChargeBondDipoleSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                      exit(0);
                      break;
                  }
    
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=energy;
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      if(rr>CutOffChargeBondDipoleSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                      exit(0);
                      break;
                  }
    
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=energy;
                }
              } 
            }
          }
        }
  
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          TypeB=Cations[CurrentSystem][k].Type;
          for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
          {
            B1=Components[TypeB].BondDipoles[l].A;
            B2=Components[TypeB].BondDipoles[l].B;
            DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
  
            posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
            posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                    exit(0);
                    break;
                }
    
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                UCationChargeBondDipoleRealDelta[CurrentSystem]+=energy;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                    exit(0);
                    break;
                }
    
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                UCationChargeBondDipoleRealDelta[CurrentSystem]-=energy;
              }
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
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(m!=k)
          {
            for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
            {
              Type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[Type].HasCharges)
              {
                posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
  
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
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SHIFTED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SMOOTHED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        if(rr>CutOffChargeBondDipoleSwitchSquared)
                        {
                          SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                         SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                        printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                        exit(0);
                        break;
                    }
    
                    cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                    energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                    UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=energy;
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
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SHIFTED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SMOOTHED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        if(rr>CutOffChargeBondDipoleSwitchSquared)
                        {
                          SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                         SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                        printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                        exit(0);
                        break;
                    }
    
                    cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                    energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                    UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=energy;
                  }
                } 
              }
            }
          }
        }
      }
  
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
        {
          Type=Cations[CurrentSystem][k].Atoms[l].Type;
          if(PseudoAtoms[Type].HasCharges)
          {
            posB=Cations[CurrentSystem][k].Atoms[l].Position;
            ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;
  
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                    exit(0);
                    break;
                }
    
                cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                UCationChargeBondDipoleRealDelta[CurrentSystem]-=energy;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                    exit(0);
                    break;
                }
    
                cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                UCationChargeBondDipoleRealDelta[CurrentSystem]+=energy;
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
      ChargeA=Adsorbates[CurrentSystem][m].Atoms[j].Charge;
  
      new_posA=TrialPosition[CurrentSystem][j];
      old_posA=Adsorbates[CurrentSystem][m].Atoms[j].Position;
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            TypeB=Adsorbates[CurrentSystem][k].Type;
            for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
            {
              B1=Components[TypeB].BondDipoles[l].A;
              B2=Components[TypeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
  
              posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
              posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                    exit(0);
                    break;
                }
  
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=energy;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                    exit(0);
                    break;
                }
  
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=energy;
              }
  
            }
          }
        }
  
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          TypeB=Cations[CurrentSystem][k].Type;
          for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
          {
            B1=Components[TypeB].BondDipoles[l].A;
            B2=Components[TypeB].BondDipoles[l].B;
            DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
  
            posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
            posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                  exit(0);
                  break;
              }
  
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
              UCationChargeBondDipoleRealDelta[CurrentSystem]+=energy;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                  exit(0);
                  break;
              }
  
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
              UCationChargeBondDipoleRealDelta[CurrentSystem]-=energy;
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
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(m!=k)
          {
            for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
            {
              Type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[Type].HasCharges)
              {
                posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
  
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SMOOTHED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      if(rr>CutOffChargeBondDipoleSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                      exit(0);
                      break;
                  }
  
                  cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=energy;
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SMOOTHED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      if(rr>CutOffChargeBondDipoleSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                      exit(0);
                      break;
                  }
  
                  cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=energy;
                }
  
              }
            }
          }
        }
      }
  
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
        {
          Type=Cations[CurrentSystem][k].Atoms[l].Type;
          if(PseudoAtoms[Type].HasCharges)
          {
            posB=Cations[CurrentSystem][k].Atoms[l].Position;
            ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;
  
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                  exit(0);
                  break;
              }
  
              cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
              UCationChargeBondDipoleRealDelta[CurrentSystem]-=energy;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  if(rr>CutOffChargeBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate'\n");
                  exit(0);
                  break;
              }
  
              cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
              UCationChargeBondDipoleRealDelta[CurrentSystem]+=energy;
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(int m)
{
  int j,k,l;
  int A1,A2,B1,B2;
  int TypeA,TypeB;
  REAL r,energy;
  VECTOR new_posA,old_posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL ri2,rk2,rr,length,temp;
  VECTOR new_dipoleA,old_dipoleA,dipoleB;
  REAL cosAB,cosA,cosB;
  REAL Bt0,Bt1,Bt2,SwitchingValue;
  int ncell;

  UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;
  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

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
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            TypeB=Adsorbates[CurrentSystem][k].Type;
            for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
            {
              B1=Components[TypeB].BondDipoles[l].A;
              B2=Components[TypeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
              posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
              posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
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
                        SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate'\n");
                      exit(0);
                      break;
                  }
    
                  cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
                  cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
                  UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
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
                        SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate'\n");
                      exit(0);
                      break;
                  }
    
                  cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
                  cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
                  UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
                }
              } 
            }
          }
        }
      }
  
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        TypeB=Cations[CurrentSystem][k].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
          posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
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
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate'\n");
                  exit(0);
                  break;
              }
    
              cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
              cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
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
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate'\n");
                  exit(0);
                  break;
              }
    
              cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
              cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              UCationBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
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
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            TypeB=Adsorbates[CurrentSystem][k].Type;
            for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
            {
              B1=Components[TypeB].BondDipoles[l].A;
              B2=Components[TypeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
              posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
              posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
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
                      SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate'\n");
                    exit(0);
                    break;
                }
  
                cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
                cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
                UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
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
                      SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate'\n");
                    exit(0);
                    break;
                }
  
                cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
                cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
                UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
              }
  
            }
          }
        }
      }
  
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        TypeB=Cations[CurrentSystem][k].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
          posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
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
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
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
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate'\n");
                exit(0);
                break;
            }
  
            cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
            cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
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
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
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
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate'\n");
                exit(0);
                break;
            }
  
            cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
            cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            UCationBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
          }
        }
      }
    }
  }
  return 0;
}


int CalculateInterChargeBondDipoleEnergyDifferenceCation(int m)
{
  int j,k,l;
  int A1,A2,B1,B2;
  int Type,TypeA,TypeB;
  VECTOR new_posA,old_posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,rr,ri2,rk2,cosA,cosB,energy,temp,length,ChargeB,ChargeA;
  VECTOR new_dipoleA,old_dipoleA,dipoleB;
  REAL Bt0,Bt1,Bt2;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL SwitchingValue;
  int ncell;

  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UCationChargeBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  Bt0=Bt1=Bt2=0.0;  
  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][m].Atoms[j].Type;
      ChargeA=Cations[CurrentSystem][m].Atoms[j].Charge;
  
      new_posA=TrialPosition[CurrentSystem][j];
      old_posA=Cations[CurrentSystem][m].Atoms[j].Position;
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            TypeB=Cations[CurrentSystem][k].Type;
            for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
            {
              B1=Components[TypeB].BondDipoles[l].A;
              B2=Components[TypeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
  
              posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
              posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SMOOTHED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
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
                      printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                      exit(0);
                      break;
                  }
    
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                  UCationChargeBondDipoleRealDelta[CurrentSystem]+=energy;
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SMOOTHED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
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
                      printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                      exit(0);
                      break;
                  }
    
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                  UCationChargeBondDipoleRealDelta[CurrentSystem]-=energy;
                }
              } 
            }
          }
        }
  
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          TypeB=Adsorbates[CurrentSystem][k].Type;
          for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
          {
            B1=Components[TypeB].BondDipoles[l].A;
            B2=Components[TypeB].BondDipoles[l].B;
            DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
  
            posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
            posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                    exit(0);
                    break;
                }
    
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=energy;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                    exit(0);
                    break;
                }
    
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=energy;
              }
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
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(m!=k)
          {
            for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
            {
              Type=Cations[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[Type].HasCharges)
              {
                posB=Cations[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;

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
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SHIFTED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SMOOTHED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
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
                        printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                        exit(0);
                        break;
                    }
    
                    cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                    energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                    UCationChargeBondDipoleRealDelta[CurrentSystem]-=energy;
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
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SHIFTED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        break;
                      case SMOOTHED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
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
                        printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                        exit(0);
                        break;
                    }
    
                    cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                    energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                    UCationChargeBondDipoleRealDelta[CurrentSystem]+=energy;
                  }
                } 
              }
            }
          }
        }
      }
  
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
        {
          Type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
          if(PseudoAtoms[Type].HasCharges)
          {
            posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
            ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
 
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                    exit(0);
                    break;
                }
    
                cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=energy;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                    exit(0);
                    break;
                }
    
                cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=energy;
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
      ChargeA=Cations[CurrentSystem][m].Atoms[j].Charge;
  
      new_posA=TrialPosition[CurrentSystem][j];
      old_posA=Cations[CurrentSystem][m].Atoms[j].Position;
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            TypeB=Cations[CurrentSystem][k].Type;
            for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
            {
              B1=Components[TypeB].BondDipoles[l].A;
              B2=Components[TypeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
  
              posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
              posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                    exit(0);
                    break;
                }
  
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                UCationChargeBondDipoleRealDelta[CurrentSystem]+=energy;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
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
                    printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                    exit(0);
                    break;
                }
  
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
                UCationChargeBondDipoleRealDelta[CurrentSystem]-=energy;
              }
  
            }
          }
        }
  
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          TypeB=Adsorbates[CurrentSystem][k].Type;
          for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
          {
            B1=Components[TypeB].BondDipoles[l].A;
            B2=Components[TypeB].BondDipoles[l].B;
            DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
  
            posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
            posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
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
                  printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                  exit(0);
                  break;
              }
  
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
              UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=energy;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
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
                  printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                  exit(0);
                  break;
              }
  
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
              UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=energy;
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
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(m!=k)
          {
            for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
            {
              Type=Cations[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[Type].HasCharges)
              {
                posB=Cations[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;
  
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SMOOTHED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
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
                      printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                      exit(0);
                      break;
                  }
  
                  cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                  UCationChargeBondDipoleRealDelta[CurrentSystem]-=energy;
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SMOOTHED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
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
                      printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                      exit(0);
                      break;
                  }
  
                  cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
                  UCationChargeBondDipoleRealDelta[CurrentSystem]+=energy;
                }
  
              }
            }
          }
        }
      }
  
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
        {
          Type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
          if(PseudoAtoms[Type].HasCharges)
          {
            posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
            ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
  
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
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
                  printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                  exit(0);
                  break;
              }
  
              cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
              UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=energy;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
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
                  printf("Unknown charge-method in 'CalculateInterChargeBondDipoleEnergyDifferenceCation'\n");
                  exit(0);
                  break;
              }
  
              cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);
              UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=energy;
            }
          }
        }
      }
    }
  }
  return 0;
}


int CalculateInterBondDipoleBondDipoleEnergyDifferenceCation(int m)
{
  int j,k,l;
  int A1,A2,B1,B2;
  int TypeA,TypeB;
  REAL r,energy;
  VECTOR new_posA,old_posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL ri2,rk2,rr,length,temp;
  VECTOR new_dipoleA,old_dipoleA,dipoleB;
  REAL cosAB,cosA,cosB;
  REAL Bt0,Bt1,Bt2,SwitchingValue;
  int ncell;

  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;
  UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

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
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            TypeB=Cations[CurrentSystem][k].Type;
            for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
            {
              B1=Components[TypeB].BondDipoles[l].A;
              B2=Components[TypeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
              posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
              posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
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
                        SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceCation'\n");
                      exit(0);
                      break;
                  }
    
                  cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
                  cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
                  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
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
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
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
                        SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                      printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceCation'\n");
                      exit(0);
                      break;
                  }
    
                  cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
                  cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
                  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
                }
              }
            }
          }
        }
      }
  
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        TypeB=Adsorbates[CurrentSystem][k].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
          posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
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
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceCation'\n");
                  exit(0);
                  break;
              }
    
              cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
              cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
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
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  break;
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
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                  printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceCation'\n");
                  exit(0);
                  break;
              }
    
              cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
              cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
              UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
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
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            TypeB=Cations[CurrentSystem][k].Type;
            for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
            {
              B1=Components[TypeB].BondDipoles[l].A;
              B2=Components[TypeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
              posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
              posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
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
                      SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceCation'\n");
                    exit(0);
                    break;
                }
  
                cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
                cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
                UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
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
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
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
                      SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                    printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceCation'\n");
                    exit(0);
                    break;
                }
  
                cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
                cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
                UCationBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
              }
            }
          }
        }
      }
  
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        TypeB=Adsorbates[CurrentSystem][k].Type;
        for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
        {
          B1=Components[TypeB].BondDipoles[l].A;
          B2=Components[TypeB].BondDipoles[l].B;
          DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
          posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
          posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
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
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
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
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceCation'\n");
                exit(0);
                break;
            }
  
            cosAB=new_dipoleA.x*dipoleB.x+new_dipoleA.y*dipoleB.y+new_dipoleA.z*dipoleB.z;
            cosA=new_dipoleA.x*dr.x+new_dipoleA.y*dr.y+new_dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+=energy;
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
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
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
                  SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                 SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
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
                printf("Unknown charge-method in 'CalculateInterBondDipoleBondDipoleEnergyDifferenceCation'\n");
                exit(0);
                break;
            }
  
            cosAB=old_dipoleA.x*dipoleB.x+old_dipoleA.y*dipoleB.y+old_dipoleA.z*dipoleB.z;
            cosA=old_dipoleA.x*dr.x+old_dipoleA.y*dr.y+old_dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
            energy=SwitchingValue*COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
            UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]-=energy;
          }
        }
      }
    }
  }
  return 0;
}

REAL CalculateInterChargeElectrostaticPotential(POINT posA)
{
  int i,j,typeB;
  POINT posB;
  REAL r,r2;
  REAL chargeB;
  VECTOR dr;
  REAL UElectrostaticPotential,energy;
  REAL SwitchingValue,TranslationValue;

  UElectrostaticPotential=0;

  if(OmitInterMolecularInteractions) return 0;

  for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
  {
    for(i=0;i<Components[Adsorbates[CurrentSystem][j].Type].NumberOfAtoms;i++)
    {
      posB=Adsorbates[CurrentSystem][j].Atoms[i].Position;
      typeB=Adsorbates[CurrentSystem][j].Atoms[i].Type;
      chargeB=Adsorbates[CurrentSystem][j].Atoms[i].Charge;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);

      r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      if(r2<CutOffChargeChargeSquared)
      {
        r=sqrt(r2);
        switch(ChargeMethod)
        {
          case NONE:
            energy=0.0;
            break;
          case TRUNCATED_COULOMB:
            energy=COULOMBIC_CONVERSION_FACTOR*chargeB/r;
            break;
          case SHIFTED_COULOMB:
            energy=COULOMBIC_CONVERSION_FACTOR*chargeB*(1.0/r-InverseCutOffChargeCharge);
            break;
          case SMOOTHED_COULOMB:
            energy=COULOMBIC_CONVERSION_FACTOR*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
            if(r2>CutOffChargeChargeSwitchSquared)
            {
               SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                              SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

               TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeB*
                                (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                 SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                 SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
               energy=energy*SwitchingValue+TranslationValue;
             }
             break;
          case WOLFS_METHOD_DAMPED_FG:
            energy=COULOMBIC_CONVERSION_FACTOR*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                         -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                          (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                          (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
            break;
          case EWALD:
            energy=COULOMBIC_CONVERSION_FACTOR*chargeB*
                                 (erfc(Alpha[CurrentSystem]*r)/r);
            break;
          default:
            printf("Unknown charge-method in 'CalculateInterChargeElectrostaticPotential'\n");
            exit(0);
            break;
        }
        UElectrostaticPotential+=energy;
      }
    }
  }

  for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
  {
    for(i=0;i<Components[Cations[CurrentSystem][j].Type].NumberOfAtoms;i++)
    {
      posB=Cations[CurrentSystem][j].Atoms[i].Position;
      typeB=Cations[CurrentSystem][j].Atoms[i].Type;
      chargeB=Cations[CurrentSystem][j].Atoms[i].Charge;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);

      r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      if(r2<CutOffChargeChargeSquared)
      {
        r=sqrt(r2);
        switch(ChargeMethod)
        {
          case NONE:
            energy=0.0;
            break;
          case TRUNCATED_COULOMB:
            energy=COULOMBIC_CONVERSION_FACTOR*chargeB/r;
            break;
          case SHIFTED_COULOMB:
            energy=COULOMBIC_CONVERSION_FACTOR*chargeB*(1.0/r-InverseCutOffChargeCharge);
            break;
          case SMOOTHED_COULOMB:
            energy=COULOMBIC_CONVERSION_FACTOR*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
            if(r2>CutOffChargeChargeSwitchSquared)
            {
               SwitchingValue=SwitchingChargeChargeFactors5[5]*(r2*r2*r)+SwitchingChargeChargeFactors5[4]*(r2*r2)+SwitchingChargeChargeFactors5[3]*(r2*r)+
                              SwitchingChargeChargeFactors5[2]*r2+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];

               TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeB*
                                (SwitchingChargeChargeFactors7[7]*(r2*r2*r2*r)+SwitchingChargeChargeFactors7[6]*(r2*r2*r2)+
                                 SwitchingChargeChargeFactors7[5]*(r2*r2*r)+SwitchingChargeChargeFactors7[4]*(r2*r2)+SwitchingChargeChargeFactors7[3]*(r2*r)+
                                 SwitchingChargeChargeFactors7[2]*r2+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
               energy=energy*SwitchingValue+TranslationValue;
            }
            break;
          case WOLFS_METHOD_DAMPED_FG:
            energy=COULOMBIC_CONVERSION_FACTOR*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                         -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                          (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                          (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
            break;
          case EWALD:
            energy=COULOMBIC_CONVERSION_FACTOR*chargeB*
                                 (erfc(Alpha[CurrentSystem]*r)/r);
            break;
          default:
            printf("Unknown charge-method in 'CalculateInterChargeElectrostaticPotential'\n");
            exit(0);
            break;
        }
        UElectrostaticPotential+=energy;
      }
    }
  }

  return UElectrostaticPotential;
}

REAL CalculateInterVDWEnergyCorrectionAdsorbate(VECTOR* Positions,VECTOR *AnisotropicPositions,int exclude)
{
  int i,j,k;
  int typeA,typeB;
  REAL rr,UVDWDelta;
  VECTOR posA,posB,dr;
  int TypeMolB;

  UVDWDelta=0.0; 
  for(k=0;k<Components[CurrentComponent].NumberOfAtoms;k++)
  {
    typeA=Components[CurrentComponent].Type[k];

    if(PseudoAtoms[typeA].AnisotropicCorrection)
    {
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        if(i!=exclude)
        {
          TypeMolB=Adsorbates[CurrentSystem][i].Type;
          for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
          {
            posB=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
            typeB=Adsorbates[CurrentSystem][i].Atoms[j].Type;

            posA=AnisotropicPositions[k];
            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UVDWDelta+=PotentialValue(typeA,typeB,rr);

            posA=Positions[k];
            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UVDWDelta-=PotentialValue(typeA,typeB,rr);
          }
        }
      }
    }
  }
  return UVDWDelta;
}

REAL CalculateInterVDWEnergyCorrectionCation(VECTOR* Positions,VECTOR *AnisotropicPositions,int exclude)
{
  int i,j,k;
  int typeA,typeB;
  REAL rr,UVDWDelta;
  VECTOR posA,posB,dr;
  int TypeMolB;

  UVDWDelta=0.0; 
  for(k=0;k<Components[CurrentComponent].NumberOfAtoms;k++)
  {
    typeA=Components[CurrentComponent].Type[k];

    if(PseudoAtoms[typeA].AnisotropicCorrection)
    {
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        if(i!=exclude)
        {
          TypeMolB=Cations[CurrentSystem][i].Type;
          for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
          {
            posB=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
            typeB=Cations[CurrentSystem][i].Atoms[j].Type;

            posA=AnisotropicPositions[k];
            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UVDWDelta+=PotentialValue(typeA,typeB,rr);

            posA=Positions[k];
            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UVDWDelta-=PotentialValue(typeA,typeB,rr);
          }
        }
      }
    }
  }
  return UVDWDelta;
}
