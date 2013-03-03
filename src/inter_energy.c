/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'inter_energy.c' is part of RASPA-2.0

 *************************************************************************************************************/

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

REAL CalculateInterVDWSelfEnergyCorrectionNew(void)
{
  int j,k;
  int typeA,typeB;
  REAL UVDWCorrection;
  VECTOR posA,posB,dr;
  REAL scalingA,scalingB;
  REAL rr,energy;
  int ncell;

  UVDWCorrection=0.0;
  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Components[CurrentComponent].NumberOfAtoms;k++)
    {
      typeA=Components[CurrentComponent].Type[k];
      posA=TrialAnisotropicPosition[CurrentSystem][k];
      scalingA=CFVDWScaling[k];

      for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
      {
        typeB=Components[CurrentComponent].Type[j];
        posB=TrialAnisotropicPosition[CurrentSystem][j];
        scalingB=CFVDWScaling[j];

        for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
          dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
          dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
          dr=ApplyReplicaBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
          {
            energy=PotentialValue(typeA,typeB,rr,scalingA*scalingB);
            if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
              UVDWCorrection+=0.5*energy;
          }
        }
      }
    }
  }
  return UVDWCorrection;
}

REAL CalculateInterVDWSelfEnergyCorrectionAdsorbateOld(int mol)
{
  int j,k;
  int typeA,typeB;
  REAL UVDWCorrection;
  VECTOR posA,posB,dr;
  REAL scalingA,scalingB;
  REAL rr;
  int ncell;

  UVDWCorrection=0.0;
  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Adsorbates[CurrentSystem][mol].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][mol].Atoms[k].Type;
      posA=Adsorbates[CurrentSystem][mol].Atoms[k].AnisotropicPosition;
      scalingA=Adsorbates[CurrentSystem][mol].Atoms[k].CFVDWScalingParameter;

      for(j=0;j<Adsorbates[CurrentSystem][mol].NumberOfAtoms;j++)
      {
        typeB=Adsorbates[CurrentSystem][mol].Atoms[j].Type;
        posB=Adsorbates[CurrentSystem][mol].Atoms[j].AnisotropicPosition;
        scalingB=Adsorbates[CurrentSystem][mol].Atoms[j].CFVDWScalingParameter;

        for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
          dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
          dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
          dr=ApplyReplicaBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
            UVDWCorrection+=0.5*PotentialValue(typeA,typeB,rr,scalingA*scalingB);
        }
      }
    }
  }
  return UVDWCorrection;
}

REAL CalculateInterVDWSelfEnergyCorrectionCationOld(int mol)
{
  int j,k;
  int typeA,typeB;
  REAL UVDWCorrection;
  VECTOR posA,posB,dr;
  REAL scalingA,scalingB;
  REAL rr;
  int ncell;

  UVDWCorrection=0.0;
  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Cations[CurrentSystem][mol].NumberOfAtoms;k++)
    {
      typeA=Cations[CurrentSystem][mol].Atoms[k].Type;
      posA=Cations[CurrentSystem][mol].Atoms[k].AnisotropicPosition;
      scalingA=Cations[CurrentSystem][mol].Atoms[k].CFVDWScalingParameter;

      for(j=0;j<Cations[CurrentSystem][mol].NumberOfAtoms;j++)
      {
        typeB=Cations[CurrentSystem][mol].Atoms[j].Type;
        posB=Cations[CurrentSystem][mol].Atoms[j].AnisotropicPosition;
        scalingB=Cations[CurrentSystem][mol].Atoms[j].CFVDWScalingParameter;

        for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
          dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
          dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
          dr=ApplyReplicaBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if(rr<CutOffVDWSquared)
            UVDWCorrection+=0.5*PotentialValue(typeA,typeB,rr,scalingA*scalingB);
        }
      }
    }
  }
  return UVDWCorrection;
}

REAL CalculateInterChargeChargeSelfEnergyCorrectionNew(void)
{
  int j,k;
  int typeA,typeB;
  REAL UChargeChargeCorrection,chargeA,chargeB;
  REAL SwitchingValue,TranslationValue;
  REAL scalingA,scalingB;
  VECTOR posA,posB,dr;
  REAL rr,r,energy;
  int ncell;

  UChargeChargeCorrection=0.0;
  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Components[CurrentComponent].NumberOfAtoms;k++)
    {
      typeA=Components[CurrentComponent].Type[k];
      posA=TrialAnisotropicPosition[CurrentSystem][k];
      scalingA=CFChargeScaling[k];
      chargeA=scalingA*PseudoAtoms[typeA].Charge1;

      for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
      {
        typeB=Components[CurrentComponent].Type[j];
        posB=TrialAnisotropicPosition[CurrentSystem][j];
        scalingB=CFChargeScaling[j];
        chargeB=scalingB*PseudoAtoms[typeB].Charge1;

        for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
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
              case WOLFS_METHOD_DAMPED:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                            -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
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
                printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceAdsorbate'\n");
                exit(0);
                break;

            }
            UChargeChargeCorrection+=0.5*energy;
          }

        }
      }
    }
  }
  return UChargeChargeCorrection;
}

REAL CalculateInterChargeChargeSelfEnergyCorrectionAdsorbateOld(int mol)
{
  int j,k;
  int typeA,typeB;
  REAL UChargeChargeCorrection,chargeA,chargeB;
  REAL SwitchingValue,TranslationValue;
  REAL scalingA,scalingB;
  VECTOR posA,posB,dr;
  REAL rr,r,energy;
  int ncell;

  UChargeChargeCorrection=0.0;
  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Adsorbates[CurrentSystem][mol].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][mol].Atoms[k].Type;
      posA=Adsorbates[CurrentSystem][mol].Atoms[k].AnisotropicPosition;
      scalingA=Adsorbates[CurrentSystem][mol].Atoms[k].CFChargeScalingParameter;
      chargeA=scalingA*Adsorbates[CurrentSystem][mol].Atoms[k].Charge;

      for(j=0;j<Adsorbates[CurrentSystem][mol].NumberOfAtoms;j++)
      {
        typeB=Adsorbates[CurrentSystem][mol].Atoms[j].Type;
        posB=Adsorbates[CurrentSystem][mol].Atoms[j].AnisotropicPosition;
        scalingB=Adsorbates[CurrentSystem][mol].Atoms[j].CFChargeScalingParameter;
        chargeB=scalingB*Adsorbates[CurrentSystem][mol].Atoms[j].Charge;

        for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
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
              case WOLFS_METHOD_DAMPED:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                            -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
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
                printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceAdsorbate'\n");
                exit(0);
                break;

            }
            UChargeChargeCorrection+=0.5*energy;
          }

        }
      }
    }
  }
  return UChargeChargeCorrection;
}

REAL CalculateInterChargeChargeSelfEnergyCorrectionCationOld(int mol)
{
  int j,k;
  int typeA,typeB;
  REAL UChargeChargeCorrection,chargeA,chargeB;
  REAL SwitchingValue,TranslationValue;
  REAL scalingA,scalingB;
  VECTOR posA,posB,dr;
  REAL rr,r,energy;
  int ncell;

  UChargeChargeCorrection=0.0;
  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Cations[CurrentSystem][mol].NumberOfAtoms;k++)
    {
      typeA=Cations[CurrentSystem][mol].Atoms[k].Type;
      posA=Cations[CurrentSystem][mol].Atoms[k].AnisotropicPosition;
      scalingA=Cations[CurrentSystem][mol].Atoms[k].CFChargeScalingParameter;
      chargeA=scalingA*Cations[CurrentSystem][mol].Atoms[k].Charge;

      for(j=0;j<Cations[CurrentSystem][mol].NumberOfAtoms;j++)
      {
        typeB=Cations[CurrentSystem][mol].Atoms[j].Type;
        posB=Cations[CurrentSystem][mol].Atoms[j].AnisotropicPosition;
        scalingB=Cations[CurrentSystem][mol].Atoms[j].CFChargeScalingParameter;
        chargeB=scalingB*Cations[CurrentSystem][mol].Atoms[j].Charge;

        for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
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
              case WOLFS_METHOD_DAMPED:
                energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                            -erfc(Alpha[CurrentSystem]*CutOffChargeChargeSquared)*InverseCutOffChargeCharge);
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
                printf("Unknown charge-method in 'CalculateInterChargeChargeEnergyDifferenceCation'\n");
                exit(0);
                break;

            }
            UChargeChargeCorrection+=0.5*energy;
          }

        }
      }
    }
  }
  return UChargeChargeCorrection;
}


REAL CalculateInterVDWEnergyAdsorbateAtPosition(POINT posA,int typeA,int exclude,REAL scaling)
{
  int i,j,typeB;
  POINT posB;
  REAL r2,Uvdw;
  VECTOR dr;
  int TypeMolB;
  int ncell;
  REAL scalingB;

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
          scalingB=Adsorbates[CurrentSystem][j].Atoms[i].CFVDWScalingParameter;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
  
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
   
            if(r2<CutOffVDWSquared)
              Uvdw+=PotentialValue(typeA,typeB,r2,scaling*scalingB);
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
          scalingB=Adsorbates[CurrentSystem][j].Atoms[i].CFVDWScalingParameter;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(r2<CutOffVDWSquared)
            Uvdw+=PotentialValue(typeA,typeB,r2,scaling*scalingB);
        }
      }
    }
  }

  return Uvdw;
}

int CalculateInterChargeEnergyAdsorbateAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,int exclude,REAL scaling)
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
  REAL scalingB;
  int ncell;

  *UChargeCharge=0.0;
  *UChargeBondDipole=0.0;

  if(OmitInterMolecularInteractions) return 0;

  chargeA=scaling*PseudoAtoms[typeA].Charge1;

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
          scalingB=Adsorbates[CurrentSystem][j].Atoms[i].CFChargeScalingParameter;
          chargeB=scalingB*Adsorbates[CurrentSystem][j].Atoms[i].Charge;

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
          scalingB=Adsorbates[CurrentSystem][j].Atoms[i].CFChargeScalingParameter;
          chargeB=scalingB*Adsorbates[CurrentSystem][j].Atoms[i].Charge;
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

REAL CalculateInterVDWEnergyCationAtPosition(POINT posA,int typeA,int exclude,REAL scaling)
{
  int i,j,typeB;
  POINT posB;
  REAL r2,Uvdw;
  VECTOR dr;
  int TypeMolB;
  int ncell;
  REAL scalingB;

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
          scalingB=Cations[CurrentSystem][j].Atoms[i].CFVDWScalingParameter;

          for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
          {
            dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
            dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
            dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
            dr=ApplyReplicaBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
            if(r2<CutOffVDWSquared)
              Uvdw+=PotentialValue(typeA,typeB,r2,scaling*scalingB);
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
          scalingB=Cations[CurrentSystem][j].Atoms[i].CFVDWScalingParameter;
          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
  
          dr=ApplyBoundaryCondition(dr);
          r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
          if(r2<CutOffVDWSquared)
            Uvdw+=PotentialValue(typeA,typeB,r2,scaling*scalingB);
        }
      }
    }
  }
  return Uvdw;
}

int CalculateInterChargeEnergyCationAtPosition(POINT posA,int typeA,REAL *UChargeCharge,REAL *UChargeBondDipole,int exclude,REAL scaling)
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
  REAL scalingB;
  int ncell;

  *UChargeCharge=0.0;
  *UChargeBondDipole=0.0;

  if(OmitInterMolecularInteractions) return 0;

  chargeA=scaling*PseudoAtoms[typeA].Charge1;

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
          scalingB=Cations[CurrentSystem][j].Atoms[i].CFChargeScalingParameter;
          chargeB=scalingB*Cations[CurrentSystem][j].Atoms[i].Charge;

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
          scalingB=Cations[CurrentSystem][j].Atoms[i].CFChargeScalingParameter;
          chargeB=scalingB*Cations[CurrentSystem][j].Atoms[i].Charge;
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
  int typeA,typeB,TypeMolB;
  REAL rr,energy;
  VECTOR posA,posB,dr;
  int ncell,start;
  REAL scalingA,scalingB;

  UAdsorbateAdsorbateVDW[CurrentSystem]=0.0;
  UAdsorbateCationVDW[CurrentSystem]=0.0;
  UCationCationVDW[CurrentSystem]=0.0;

  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
        posA=Adsorbates[CurrentSystem][i].Atoms[k].AnisotropicPosition;
        scalingA=Adsorbates[CurrentSystem][i].Atoms[k].CFVDWScalingParameter;
  
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
                posB=Adsorbates[CurrentSystem][j].Atoms[l].AnisotropicPosition;

                dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffVDWSquared)
                {
                  typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
                  scalingB=Adsorbates[CurrentSystem][j].Atoms[l].CFVDWScalingParameter;
                  energy=PotentialValue(typeA,typeB,rr,scalingA*scalingB);

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

                dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                {
                  typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                  scalingB=Cations[CurrentSystem][j].Atoms[l].CFVDWScalingParameter;
                  energy=PotentialValue(typeA,typeB,rr,scalingA*scalingB);

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
      for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Cations[CurrentSystem][i].Atoms[k].Type;
        posA=Cations[CurrentSystem][i].Atoms[k].AnisotropicPosition;
        scalingA=Cations[CurrentSystem][i].Atoms[k].CFVDWScalingParameter;
  
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

              dr.x=posA.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffVDWSquared)
              {
                typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                scalingB=Cations[CurrentSystem][j].Atoms[l].CFVDWScalingParameter;
                energy=PotentialValue(typeA,typeB,rr,scalingA*scalingB);

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
      for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
        posA=Adsorbates[CurrentSystem][i].Atoms[k].AnisotropicPosition;
        scalingA=Adsorbates[CurrentSystem][i].Atoms[k].CFVDWScalingParameter;
  
        // loop over adsorbant molecules
        if(!OmitAdsorbateAdsorbateVDWInteractions)
        {
          for(j=i+1;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
          {
            TypeMolB=Adsorbates[CurrentSystem][j].Type;
            for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
            {
              posB=Adsorbates[CurrentSystem][j].Atoms[l].AnisotropicPosition;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(rr<CutOffVDWSquared)
              {
                typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
                scalingB=Adsorbates[CurrentSystem][j].Atoms[l].CFVDWScalingParameter;
                UAdsorbateAdsorbateVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr,scalingA*scalingB);
              }
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
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
              {
                typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                scalingB=Cations[CurrentSystem][j].Atoms[l].CFVDWScalingParameter;
                UAdsorbateCationVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr,scalingA*scalingB);
              }
            }
          }
        }
      }
    }
  
    if(!OmitCationCationVDWInteractions)
    {
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
        {
          typeA=Cations[CurrentSystem][i].Atoms[k].Type;
          posA=Cations[CurrentSystem][i].Atoms[k].AnisotropicPosition;
          scalingA=Cations[CurrentSystem][i].Atoms[k].CFVDWScalingParameter;
  
          // loop over cation molecules
          for(j=i+1;j<NumberOfCationMolecules[CurrentSystem];j++)
          {
            TypeMolB=Cations[CurrentSystem][j].Type;
            for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
            {
              posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
              {
                typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                scalingB=Cations[CurrentSystem][j].Atoms[l].CFVDWScalingParameter;
                UCationCationVDW[CurrentSystem]+=PotentialValue(typeA,typeB,rr,scalingA*scalingB);
              }
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
  REAL scalingA,scalingB;
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
        scalingA=Adsorbates[CurrentSystem][i].Atoms[k].CFChargeScalingParameter;
        chargeA=scalingA*Adsorbates[CurrentSystem][i].Atoms[k].Charge;
  
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
                  scalingB=Adsorbates[CurrentSystem][j].Atoms[l].CFChargeScalingParameter;
                  chargeB=scalingB*Adsorbates[CurrentSystem][j].Atoms[l].Charge;

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
                  scalingB=Cations[CurrentSystem][j].Atoms[l].CFChargeScalingParameter;
                  chargeB=scalingB*Cations[CurrentSystem][j].Atoms[l].Charge;

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
          scalingA=Cations[CurrentSystem][i].Atoms[k].CFChargeScalingParameter;
          chargeA=scalingA*Cations[CurrentSystem][i].Atoms[k].Charge;
  
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
                  scalingB=Cations[CurrentSystem][j].Atoms[l].CFChargeScalingParameter;
                  chargeB=scalingB*Cations[CurrentSystem][j].Atoms[l].Charge;

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
        scalingA=Adsorbates[CurrentSystem][i].Atoms[k].CFChargeScalingParameter;
        chargeA=scalingA*Adsorbates[CurrentSystem][i].Atoms[k].Charge;
  
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
                scalingB=Adsorbates[CurrentSystem][j].Atoms[l].CFChargeScalingParameter;
                chargeB=scalingB*Adsorbates[CurrentSystem][j].Atoms[l].Charge;
  
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
                scalingB=Cations[CurrentSystem][j].Atoms[l].CFChargeScalingParameter;
                chargeB=scalingB*Cations[CurrentSystem][j].Atoms[l].Charge;
  
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
          scalingA=Cations[CurrentSystem][i].Atoms[k].CFChargeScalingParameter;
          chargeA=scalingA*Cations[CurrentSystem][i].Atoms[k].Charge;
  
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
                scalingB=Cations[CurrentSystem][j].Atoms[l].CFChargeScalingParameter;
                chargeB=scalingB*Cations[CurrentSystem][j].Atoms[l].Charge;
  
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

/*********************************************************************************************************
 * Name       | CalculateInterVDWEnergyDifferenceAdsorbate                                               *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes inter-molecular energy differences for a single adsorbate molecule:             *
 *            | 1) New=TRUE, Old=FAlse, computes the VDW energy of the new molecule (positions in        *
 *            |    'TrialAnisotropicPosition' and the molecule is of type 'comp'                         *
 *            | 2) New=FALSE, Old=TRUE, computes the VDW energy of molecule 'm' of type 'comp'           *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in VDW energy of the molecule 'm' with    *
 *            |    the new positions in 'TrialAnisotropicPosition' and the old positions                 *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the VDW energy of the new positions stored in             *
 *            |      'TrialAnisotropicPosition'                                                          *
 *            | False: whether or not to compute the VDW energy of the old positions of molecule 'm'     *
 * Note       |                                                                                          *
 * Used       | int TranslationMoveAdsorbate(void)                                                       *
 *            | int RandomTranslationMoveAdsorbate(void)                                                 *
 *            | int RotationMoveAdsorbate(void)                                                          *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CFCBSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateInterVDWEnergyDifferenceAdsorbate(int m,int comp,int New,int Old)
{
  int i,j,k;
  int typeA,typeB;
  REAL rr,energy,scalingA_new,scalingA_old,scalingB;
  VECTOR posA_new,posA_old,posB,dr;
  int TypeMolB;
  int ncell;

  OVERLAP=FALSE;
  UAdsorbateVDWDelta[CurrentSystem]=0.0;
  UCationVDWDelta[CurrentSystem]=0.0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Components[comp].NumberOfAtoms;k++)
    {
      typeA=Components[comp].Type[k];

      if(New)
      {  
        posA_new=TrialAnisotropicPosition[CurrentSystem][k];
        scalingA_new=CFVDWScaling[k];
      }

      if(Old)
      {
        posA_old=Adsorbates[CurrentSystem][m].Atoms[k].AnisotropicPosition;
        scalingA_old=Adsorbates[CurrentSystem][m].Atoms[k].CFVDWScalingParameter;
      }

      if(!OmitAdsorbateAdsorbateVDWInteractions)
      {
        // self-interaction
        for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
        {
          typeB=Adsorbates[CurrentSystem][m].Atoms[j].Type;
          scalingB=Adsorbates[CurrentSystem][m].Atoms[j].CFVDWScalingParameter;
          if(New)
          {
            for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              posB=TrialAnisotropicPosition[CurrentSystem][j];
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr,scalingA_new*scalingB);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                  UAdsorbateVDWDelta[CurrentSystem]+=0.5*energy;
              }
            }
          }

          if(Old)
          {
            for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              posB=Adsorbates[CurrentSystem][m].Atoms[j].AnisotropicPosition;
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UAdsorbateVDWDelta[CurrentSystem]-=0.5*PotentialValue(typeA,typeB,rr,scalingA_old*scalingB);
            }
          }
        }

        // loop over all other adsorbates, except i==m
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            TypeMolB=Adsorbates[CurrentSystem][i].Type;
            for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
            {
              typeB=Adsorbates[CurrentSystem][i].Atoms[j].Type;
              posB=Adsorbates[CurrentSystem][i].Atoms[j].AnisotropicPosition;
              scalingB=Adsorbates[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

              if(New)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                { 
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  if(rr<CutOffVDWSquared)
                  {
                    energy=PotentialValue(typeA,typeB,rr,scalingA_new*scalingB);
                    if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                      UAdsorbateVDWDelta[CurrentSystem]+=energy;
                  }
                }
              }

              if(Old)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                { 
                  dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  if(rr<CutOffVDWSquared)
                    UAdsorbateVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scalingA_old*scalingB);
                }
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
          scalingB=Cations[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;
  
          if(New)
          {
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            { 
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr,scalingA_new*scalingB);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UCationVDWDelta[CurrentSystem]+=energy;
              }
            }
          }
   
          if(Old)
          { 
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            { 
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UCationVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scalingA_old*scalingB);
            }
          }
        }
      }
    }
  }
  else
  {
    for(k=0;k<Components[comp].NumberOfAtoms;k++)
    {
      typeA=Components[comp].Type[k];

      if(New)
      {
        posA_new=TrialAnisotropicPosition[CurrentSystem][k];
        scalingA_new=CFVDWScaling[k];
      }
      if(Old)
      {  
        posA_old=Adsorbates[CurrentSystem][m].Atoms[k].AnisotropicPosition;
        scalingA_old=Adsorbates[CurrentSystem][m].Atoms[k].CFVDWScalingParameter;
      }
  
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
              scalingB=Adsorbates[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;

              if(New)
              {  
                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                {
                  energy=PotentialValue(typeA,typeB,rr,scalingA_new*scalingB);
                  if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                  UAdsorbateVDWDelta[CurrentSystem]+=energy;
                }
              }
 
              if(Old) 
              {
                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                  UAdsorbateVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scalingA_old*scalingB);
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
          scalingB=Cations[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;
 
          if(New)
          { 
            dr.x=posA_new.x-posB.x;
            dr.y=posA_new.y-posB.y;
            dr.z=posA_new.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
            {
              energy=PotentialValue(typeA,typeB,rr,scalingA_new*scalingB);
              if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
              UCationVDWDelta[CurrentSystem]+=energy;
            }
          }
 
          if(Old) 
          {
            dr.x=posA_old.x-posB.x;
            dr.y=posA_old.y-posB.y;
            dr.z=posA_old.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UCationVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scalingA_old*scalingB);
          }
        }
      }
    }
  }
  return 0;
}

/*********************************************************************************************************
 * Name       | CalculateInterVDWEnergyDifferenceCation                                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes inter-molecular energy differences for a single cation molecule:                *
 *            | 1) New=TRUE, Old=FAlse, computes the VDW energy of the new molecule (positions in        *
 *            |    'TrialAnisotropicPosition' and the molecule is of type 'comp'                         *
 *            | 2) New=FALSE, Old=TRUE, computes the VDW energy of molecule 'm' of type 'comp'           *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in VDW energy of the molecule 'm' with    *
 *            |    the new positions in 'TrialAnisotropicPosition' and the old positions                 *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the VDW energy of the new positions stored in             *
 *            |      'TrialAnisotropicPosition'                                                          *
 *            | False: whether or not to compute the VDW energy of the old positions of molecule 'm'     *
 * Note       |                                                                                          *
 * Used       | int TranslationMoveCation(void)                                                          *
 *            | int RandomTranslationMoveCation(void)                                                    *
 *            | int RotationMoveCation(void)                                                             *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CFCBSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateInterVDWEnergyDifferenceCation(int m,int comp,int New,int Old)
{
  int i,j,k;
  int typeA,typeB;
  REAL rr,energy,scalingA_new,scalingA_old,scalingB;
  VECTOR posA_new,posA_old,posB,dr;
  int TypeMolB;
  int ncell;

  OVERLAP=FALSE;
  UAdsorbateVDWDelta[CurrentSystem]=0.0;
  UCationVDWDelta[CurrentSystem]=0.0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Components[comp].NumberOfAtoms;k++)
    {
      typeA=Components[comp].Type[k];

      if(New)
      {
        posA_new=TrialAnisotropicPosition[CurrentSystem][k];
        scalingA_new=CFVDWScaling[k];
      }

      if(Old)
      {
        posA_old=Cations[CurrentSystem][m].Atoms[k].AnisotropicPosition;
        scalingA_old=Cations[CurrentSystem][m].Atoms[k].CFVDWScalingParameter;
      }

      if(!OmitCationCationVDWInteractions)
      {
       // self-interaction
        for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
        {
          typeB=Cations[CurrentSystem][m].Atoms[j].Type;
          scalingB=Cations[CurrentSystem][m].Atoms[j].CFVDWScalingParameter;
          if(New)
          {
            for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              posB=TrialAnisotropicPosition[CurrentSystem][j];
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr,scalingA_new*scalingB);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                  UCationVDWDelta[CurrentSystem]+=0.5*energy;
              }
            }
          }

          if(Old)
          {
            for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              posB=Cations[CurrentSystem][m].Atoms[j].AnisotropicPosition;
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UCationVDWDelta[CurrentSystem]-=0.5*PotentialValue(typeA,typeB,rr,scalingA_old*scalingB);
            }
          }
        }

        // loop over all other cations, except i==m
        for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            TypeMolB=Cations[CurrentSystem][i].Type;
            for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Cations[CurrentSystem][i].Atoms[j].AnisotropicPosition;
              typeB=Cations[CurrentSystem][i].Atoms[j].Type;
              scalingB=Cations[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;
 
              if(New)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                {
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  if(rr<CutOffVDWSquared)
                  {
                    energy=PotentialValue(typeA,typeB,rr,scalingA_new*scalingB);
                    if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                    UCationVDWDelta[CurrentSystem]+=energy;
                  }
                }
              }

              if(Old)
              {    
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                {
                  dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  if(rr<CutOffVDWSquared)
                    UCationVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scalingA_old*scalingB);
                }
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
          scalingB=Adsorbates[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;
  
          if(New)
          {
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
              {
                energy=PotentialValue(typeA,typeB,rr,scalingA_new*scalingB);
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UAdsorbateVDWDelta[CurrentSystem]+=energy;
              }
            }
          }

          if(Old)
          {    
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<CutOffVDWSquared)
                UAdsorbateVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scalingA_old*scalingB);
            }
          }
        }
      }
    }
  }
  else
  {
    for(k=0;k<Components[comp].NumberOfAtoms;k++)
    {
      typeA=Components[comp].Type[k];

      if(New)
      {
        posA_new=TrialAnisotropicPosition[CurrentSystem][k];
        scalingA_new=CFVDWScaling[k];
      }

      if(Old)
      {
        posA_old=Cations[CurrentSystem][m].Atoms[k].AnisotropicPosition;
        scalingA_old=Cations[CurrentSystem][m].Atoms[k].CFVDWScalingParameter;
      }
  
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
              scalingB=Cations[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;
 
              if(New)
              { 
                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                {
                  energy=PotentialValue(typeA,typeB,rr,scalingA_new*scalingB);
                  if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                  UCationVDWDelta[CurrentSystem]+=energy;
                }
              }
 
              if(Old)
              { 
                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                  UCationVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scalingA_old*scalingB);
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
          scalingB=Adsorbates[CurrentSystem][i].Atoms[j].CFVDWScalingParameter;
 
          if(New)
          { 
            dr.x=posA_new.x-posB.x;
            dr.y=posA_new.y-posB.y;
            dr.z=posA_new.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
            {
              energy=PotentialValue(typeA,typeB,rr,scalingA_new*scalingB);
              if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
              UAdsorbateVDWDelta[CurrentSystem]+=energy;
            }
          }
  
          if(Old)
          {
            dr.x=posA_old.x-posB.x;
            dr.y=posA_old.y-posB.y;
            dr.z=posA_old.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UAdsorbateVDWDelta[CurrentSystem]-=PotentialValue(typeA,typeB,rr,scalingA_old*scalingB);
          }
        }
      }
    }
  }
  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateInterChargeChargeEnergyDifferenceAdsorbate                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes inter-molecular energy differences for a single adsorbate molecule:             *
 *            | 1) New=TRUE, Old=FAlse, computes the charge-charge energy of the new molecule (positions *
 *            |    in 'TrialPosition' and the molecule is of type 'comp'                                 *
 *            | 2) New=FALSE, Old=TRUE, computes the charge-charge energy of molecule 'm' of type 'comp' *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in charge-charge energy of the molecule   *
 *            |    'm' with the new positions in 'TrialPosition' and the old positions                   *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the charge-charge energy of the new positions stored in   *
 *            |      'TrialPosition'                                                                     *
 *            | False: whether or not to compute the charge-charge energy of the old positions of        *
 *            |        molecule 'm'                                                                      * 
 * Note       |                                                                                          *
 * Used       | int TranslationMoveAdsorbate(void)                                                       *
 *            | int RandomTranslationMoveAdsorbate(void)                                                 *
 *            | int RotationMoveAdsorbate(void)                                                          *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CFCBSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateInterChargeChargeEnergyDifferenceAdsorbate(int m,int comp,int New,int Old)
{
  int i,j,k;
  REAL r2,chargeB,energy;
  REAL chargeA_new,chargeA_old;
  VECTOR posA_new,posA_old,posB,dr;
  int ncell;

  UAdsorbateChargeChargeRealDelta[CurrentSystem]=0.0;
  UCationChargeChargeRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Components[comp].NumberOfAtoms;k++)
    {
      if(New)
      {
        posA_new=TrialPosition[CurrentSystem][k];
        chargeA_new=CFChargeScaling[k]*PseudoAtoms[Components[comp].Type[k]].Charge1;
      }
      if(Old)
      {
        posA_old=Adsorbates[CurrentSystem][m].Atoms[k].Position;
        chargeA_old=Adsorbates[CurrentSystem][m].Atoms[k].CFChargeScalingParameter*Adsorbates[CurrentSystem][m].Atoms[k].Charge;
      }
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        // self-interaction
        for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
        {
          chargeB=Adsorbates[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Adsorbates[CurrentSystem][m].Atoms[j].Charge;

          if(New)
          {
            for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              posB=TrialPosition[CurrentSystem][j];
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
              {
                energy=0.5*PotentialValueCoulombic(chargeA_new,chargeB,sqrt(r2));
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UAdsorbateChargeChargeRealDelta[CurrentSystem]+=energy;
              }
            }
          }

          if(Old)
          {
            for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              posB=Adsorbates[CurrentSystem][m].Atoms[j].Position;
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
                UAdsorbateChargeChargeRealDelta[CurrentSystem]-=0.5*PotentialValueCoulombic(chargeA_old,chargeB,sqrt(r2));
            }
          }
        }

        // loop over all other adsorbates, except i==m
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Adsorbates[CurrentSystem][i].Atoms[j].Position;
              chargeB=Adsorbates[CurrentSystem][i].Atoms[j].CFChargeScalingParameter*Adsorbates[CurrentSystem][i].Atoms[j].Charge;

              if(New)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                {
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  if(r2<CutOffChargeChargeSquared)
                  {
                    energy=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(r2));
                    if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                    UAdsorbateChargeChargeRealDelta[CurrentSystem]+=energy;
                  }
                }
              }
    
              if(Old)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                {
                  dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  if(r2<CutOffChargeChargeSquared)
                    UAdsorbateChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(r2));
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
          chargeB=Cations[CurrentSystem][i].Atoms[j].CFChargeScalingParameter*Cations[CurrentSystem][i].Atoms[j].Charge;

          if(New)
          {
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
              {
                energy=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(r2));
                if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                UCationChargeChargeRealDelta[CurrentSystem]+=energy;
              }
            }
          }
    
          if(Old)
          {
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
                UCationChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(r2));
            }
          }
        }
      }
    }
  }
  else
  {
    for(k=0;k<Components[comp].NumberOfAtoms;k++)
    {
      if(New)
      {
        posA_new=TrialPosition[CurrentSystem][k];
        chargeA_new=CFChargeScaling[k]*PseudoAtoms[Components[comp].Type[k]].Charge1;
      }
    
      if(Old) 
      {
        posA_old=Adsorbates[CurrentSystem][m].Atoms[k].Position;
        chargeA_old=Adsorbates[CurrentSystem][m].Atoms[k].CFChargeScalingParameter*Adsorbates[CurrentSystem][m].Atoms[k].Charge;
      }
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Adsorbates[CurrentSystem][i].Atoms[j].Position;
              chargeB=Adsorbates[CurrentSystem][i].Atoms[j].CFChargeScalingParameter*Adsorbates[CurrentSystem][i].Atoms[j].Charge;

              if(New)
              {
                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(r2<CutOffChargeChargeSquared)
                {
                  energy=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(r2));
                  if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
                  UAdsorbateChargeChargeRealDelta[CurrentSystem]+=energy;
                }
              }
 
              if(Old)
              { 
                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(r2<CutOffChargeChargeSquared)
                  UAdsorbateChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(r2));
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
          chargeB=Cations[CurrentSystem][i].Atoms[j].CFChargeScalingParameter*Cations[CurrentSystem][i].Atoms[j].Charge;

          if(New)
          {
            dr.x=posA_new.x-posB.x;
            dr.y=posA_new.y-posB.y;
            dr.z=posA_new.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(r2<CutOffChargeChargeSquared)
            {
              energy=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(r2));
              if(energy>=EnergyOverlapCriteria) return OVERLAP=TRUE;
              UCationChargeChargeRealDelta[CurrentSystem]+=energy;
            }
          }
 
          if(Old)
          { 
            dr.x=posA_old.x-posB.x;
            dr.y=posA_old.y-posB.y;
            dr.z=posA_old.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(r2<CutOffChargeChargeSquared)
              UCationChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(r2));
          }
        }
      }
    }
  }
  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateInterChargeChargeEnergyDifferenceCation                                         *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes inter-molecular energy differences for a single cation molecule:                *
 *            | 1) New=TRUE, Old=FAlse, computes the charge-charge energy of the new molecule (positions *
 *            |    in 'TrialPosition' and the molecule is of type 'comp'                                 *
 *            | 2) New=FALSE, Old=TRUE, computes the charge-charge energy of molecule 'm' of type 'comp' *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in charge-charge energy of the molecule   *
 *            |    'm' with the new positions in 'TrialPosition' and the old positions                   *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the charge-charge energy of the new positions stored in   *
 *            |      'TrialPosition'                                                                     *
 *            | False: whether or not to compute the charge-charge energy of the old positions of        *
 *            |        molecule 'm'                                                                      * 
 * Note       |                                                                                          *
 * Used       | int TranslationMoveCation(void)                                                          *
 *            | int RandomTranslationMoveCation(void)                                                    *
 *            | int RotationMoveCation(void)                                                             *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CFCBSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateInterChargeChargeEnergyDifferenceCation(int m,int comp,int New,int Old)
{
  int i,j,k;
  REAL r2,chargeB;
  REAL chargeA_new,chargeA_old;
  VECTOR posA_new,posA_old,posB,dr;
  int ncell;

  UAdsorbateChargeChargeRealDelta[CurrentSystem]=0.0;
  UCationChargeChargeRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(k=0;k<Components[comp].NumberOfAtoms;k++)
    {
      if(New)
      {
        posA_new=TrialPosition[CurrentSystem][k];
        chargeA_new=CFChargeScaling[k]*PseudoAtoms[Components[comp].Type[k]].Charge1;
      }
      if(Old)
      {
        posA_old=Cations[CurrentSystem][m].Atoms[k].Position;
        chargeA_old=Cations[CurrentSystem][m].Atoms[k].CFChargeScalingParameter*Cations[CurrentSystem][m].Atoms[k].Charge;
      }
  
      if(!OmitCationCationCoulombInteractions)
      {
        // self-interaction
        for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
        {
          chargeB=Cations[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Cations[CurrentSystem][m].Atoms[j].Charge;

          if(New)
          {
            for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              posB=TrialPosition[CurrentSystem][j];
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
                UCationChargeChargeRealDelta[CurrentSystem]+=0.5*PotentialValueCoulombic(chargeA_new,chargeB,sqrt(r2));
            }
          }

          if(Old)
          {
            for(ncell=1;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              posB=Cations[CurrentSystem][m].Atoms[j].Position;
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
                UCationChargeChargeRealDelta[CurrentSystem]-=0.5*PotentialValueCoulombic(chargeA_old,chargeB,sqrt(r2));
            }
          }
        }

        // loop over all other cations, except i==m
        for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Cations[CurrentSystem][i].Atoms[j].Position;
              chargeB=Cations[CurrentSystem][i].Atoms[j].CFChargeScalingParameter*Cations[CurrentSystem][i].Atoms[j].Charge;

              if(New)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                {
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  if(r2<CutOffChargeChargeSquared)
                    UCationChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(r2));
                }
              }
    
              if(Old)
              {
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                {
                  dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  if(r2<CutOffChargeChargeSquared)
                    UCationChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(r2));
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
          chargeB=Adsorbates[CurrentSystem][i].Atoms[j].CFChargeScalingParameter*Adsorbates[CurrentSystem][i].Atoms[j].Charge;

          if(New)
          {
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
                UAdsorbateChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(r2));
            }
          }
    
          if(Old)
          {
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            {
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(r2<CutOffChargeChargeSquared)
                UAdsorbateChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(r2));
            }
          }
        }
      }
    }
  }
  else
  {
    for(k=0;k<Components[comp].NumberOfAtoms;k++)
    {
      if(New)
      {
        posA_new=TrialPosition[CurrentSystem][k];
        chargeA_new=CFChargeScaling[k]*PseudoAtoms[Components[comp].Type[k]].Charge1;
      }
    
      if(Old) 
      {
        posA_old=Cations[CurrentSystem][m].Atoms[k].Position;
        chargeA_old=Cations[CurrentSystem][m].Atoms[k].CFChargeScalingParameter*Cations[CurrentSystem][m].Atoms[k].Charge;
      }
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
        {
          if(i!=m)
          {
            for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
            {
              posB=Cations[CurrentSystem][i].Atoms[j].Position;
              chargeB=Cations[CurrentSystem][i].Atoms[j].CFChargeScalingParameter*Cations[CurrentSystem][i].Atoms[j].Charge;

              if(New)
              {
                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(r2<CutOffChargeChargeSquared)
                  UCationChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(r2));
              }
 
              if(Old)
              { 
                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(r2<CutOffChargeChargeSquared)
                  UCationChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(r2));
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
          chargeB=Adsorbates[CurrentSystem][i].Atoms[j].CFChargeScalingParameter*Adsorbates[CurrentSystem][i].Atoms[j].Charge;

          if(New)
          {
            dr.x=posA_new.x-posB.x;
            dr.y=posA_new.y-posB.y;
            dr.z=posA_new.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(r2<CutOffChargeChargeSquared)
              UAdsorbateChargeChargeRealDelta[CurrentSystem]+=PotentialValueCoulombic(chargeA_new,chargeB,sqrt(r2));
          }
 
          if(Old)
          { 
            dr.x=posA_old.x-posB.x;
            dr.y=posA_old.y-posB.y;
            dr.z=posA_old.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(r2<CutOffChargeChargeSquared)
              UAdsorbateChargeChargeRealDelta[CurrentSystem]-=PotentialValueCoulombic(chargeA_old,chargeB,sqrt(r2));
          }
        }
      }
    }
  }
  return 0;
}

/*********************************************************************************************************
 * Name       | CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes inter-molecular energy differences for a single adsorbate molecule:             *
 *            | 1) New=TRUE, Old=FAlse, computes the charge-bonddipole energy of the new molecule        *
 *            |    (positions in 'TrialPosition' and the molecule is of type 'comp'                      *
 *            | 2) New=FALSE, Old=TRUE, computes the charge-bonddipole energy of molecule 'm' of type    *
 *            |    'comp'                                                                                *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in charge-bonddipole energy of the        *
 *            |    molecule 'm' with the new positions in 'TrialPosition' and the old positions          *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the charge-bonddipole energy of the new positions stored  *
 *            |      in 'TrialPosition'                                                                  *
 *            | False: whether or not to compute the charge-bonddipole energy of the old positions of    *
 *            |        molecule 'm'                                                                      * 
 * Note       |                                                                                          *
 * Used       | int TranslationMoveAdsorbate(void)                                                       *
 *            | int RandomTranslationMoveAdsorbate(void)                                                 *
 *            | int RotationMoveAdsorbate(void)                                                          *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CFCBSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateInterChargeBondDipoleEnergyDifferenceAdsorbate(int m,int comp,int New,int Old)
{
  int j,k,l;
  int A1,A2,B1,B2;
  int type,typeB;
  VECTOR posA_new,posA_old,posB,posA1,posA2,posB1,posB2,dr;
  REAL rr,ri2,rk2,temp,length,chargeB,chargeA_new,chargeA_old;
  VECTOR dipoleA_new,dipoleA_old,dipoleB;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL scalingA1,scalingA2,scalingB1,scalingB2;
  int ncell;

  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UCationChargeBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Components[comp].NumberOfAtoms;j++)
    {
      if(New)
      {
        posA_new=TrialPosition[CurrentSystem][j];
        chargeA_new=CFChargeScaling[j]*PseudoAtoms[Components[comp].Type[j]].Charge1;
      }

      if(Old)
      {
        posA_old=Adsorbates[CurrentSystem][m].Atoms[j].Position;
        chargeA_old=Adsorbates[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Adsorbates[CurrentSystem][m].Atoms[j].Charge;
      }
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            typeB=Adsorbates[CurrentSystem][k].Type;
            for(l=0;l<Components[typeB].NumberOfBondDipoles;l++)
            {
              B1=Components[typeB].BondDipoles[l].A;
              B2=Components[typeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[typeB].BondDipoleMagnitude[l];
  
              posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
              posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
              scalingB1=Adsorbates[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
              scalingB2=Adsorbates[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
              dipoleB.x=posB2.x-posB1.x;
              dipoleB.y=posB2.y-posB1.y;
              dipoleB.z=posB2.z-posB1.z;
              posB.x=posB1.x+0.5*dipoleB.x;
              posB.y=posB1.y+0.5*dipoleB.y;
              posB.z=posB1.z+0.5*dipoleB.z;
              rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
              length=sqrt(rk2);
              temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
              dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

              if(New)
              { 
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                { 
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                  if(rr<CutOffChargeBondDipoleSquared)
                    UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
                }
              }
    
              if(Old)
              { 
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                { 
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                  if(rr<CutOffChargeBondDipoleSquared)
                    UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
                }
              } 
            }
          }
        }
      }
  
      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        typeB=Cations[CurrentSystem][k].Type;
        for(l=0;l<Components[typeB].NumberOfBondDipoles;l++)
        {
          B1=Components[typeB].BondDipoles[l].A;
          B2=Components[typeB].BondDipoles[l].B;
          DipoleMagnitudeB=Components[typeB].BondDipoleMagnitude[l];
  
          posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
          posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
          scalingB1=Cations[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
          scalingB2=Cations[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
 
          if(New)
          { 
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            { 
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
              if(rr<CutOffChargeBondDipoleSquared)
                UCationChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
            }
          }
    
          if(Old)
          { 
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            { 
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
              if(rr<CutOffChargeBondDipoleSquared)
                UCationChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
            }
          }
        }
      }
    }
  
    for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
    {
      A1=Components[comp].BondDipoles[j].A;
      A2=Components[comp].BondDipoles[j].B;
      DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];
 
      if(New)
      { 
        posA1=TrialPosition[CurrentSystem][A1];
        posA2=TrialPosition[CurrentSystem][A2];
        scalingA1=CFChargeScaling[A1];
        scalingA2=CFChargeScaling[A2];
        dipoleA_new.x=posA2.x-posA1.x;
        dipoleA_new.y=posA2.y-posA1.y;
        dipoleA_new.z=posA2.z-posA1.z;
        posA_new.x=posA1.x+0.5*dipoleA_new.x;
        posA_new.y=posA1.y+0.5*dipoleA_new.y;
        posA_new.z=posA1.z+0.5*dipoleA_new.z;
        ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
      }
 
      if(Old)
      { 
        posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
        posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
        scalingA1=Adsorbates[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
        scalingA2=Adsorbates[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
        dipoleA_old.x=posA2.x-posA1.x;
        dipoleA_old.y=posA2.y-posA1.y;
        dipoleA_old.z=posA2.z-posA1.z;
        posA_old.x=posA1.x+0.5*dipoleA_old.x;
        posA_old.y=posA1.y+0.5*dipoleA_old.y;
        posA_old.z=posA1.z+0.5*dipoleA_old.z;
        ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
      }
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(m!=k)
          {
            for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
            {
              type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[type].HasCharges)
              {
                posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                chargeB=Adsorbates[CurrentSystem][k].Atoms[l].CFChargeScalingParameter*Adsorbates[CurrentSystem][k].Atoms[l].Charge;
 
                if(New) 
                {
                  for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                  { 
                    dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                    dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                    dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                    dr=ApplyReplicaBoundaryCondition(dr);
                    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                    if(rr<CutOffChargeBondDipoleSquared)
                      UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
                  }
                }
      
                if(Old) 
                {
                  for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                  { 
                    dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                    dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                    dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                    dr=ApplyReplicaBoundaryCondition(dr);
                    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
     
                    if(rr<CutOffChargeBondDipoleSquared)
                      UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
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
          type=Cations[CurrentSystem][k].Atoms[l].Type;
          if(PseudoAtoms[type].HasCharges)
          {
            posB=Cations[CurrentSystem][k].Atoms[l].Position;
            chargeB=Cations[CurrentSystem][k].Atoms[l].CFChargeScalingParameter*Cations[CurrentSystem][k].Atoms[l].Charge;
 
            if(New)
            { 
              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              { 
                dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                if(rr<CutOffChargeBondDipoleSquared)
                  UCationChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
              }
            }
      
            if(Old)
            { 
              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              { 
                dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                if(rr<CutOffChargeBondDipoleSquared)
                  UCationChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
              }
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<Components[comp].NumberOfAtoms;j++)
    {
      if(New)
      {
        posA_new=TrialPosition[CurrentSystem][j];
        chargeA_new=CFChargeScaling[j]*PseudoAtoms[Components[comp].Type[j]].Charge1;
      }

      if(Old)
      {
        posA_old=Adsorbates[CurrentSystem][m].Atoms[j].Position;
        chargeA_old=Adsorbates[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Adsorbates[CurrentSystem][m].Atoms[j].Charge;
      }
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            typeB=Adsorbates[CurrentSystem][k].Type;
            for(l=0;l<Components[typeB].NumberOfBondDipoles;l++)
            {
              B1=Components[typeB].BondDipoles[l].A;
              B2=Components[typeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[typeB].BondDipoleMagnitude[l];
  
              posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
              posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
              scalingB1=Adsorbates[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
              scalingB2=Adsorbates[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
              dipoleB.x=posB2.x-posB1.x;
              dipoleB.y=posB2.y-posB1.y;
              dipoleB.z=posB2.z-posB1.z;
              posB.x=posB1.x+0.5*dipoleB.x;
              posB.y=posB1.y+0.5*dipoleB.y;
              posB.z=posB1.z+0.5*dipoleB.z;
              rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
              length=sqrt(rk2);
              temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
              dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
 
              if(New)
              { 
                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
                if(rr<CutOffChargeBondDipoleSquared)
                  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
              }
 
              if(Old)
              { 
                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
                if(rr<CutOffChargeBondDipoleSquared)
                  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
              }
            }
          }
        }
  
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          typeB=Cations[CurrentSystem][k].Type;
          for(l=0;l<Components[typeB].NumberOfBondDipoles;l++)
          {
            B1=Components[typeB].BondDipoles[l].A;
            B2=Components[typeB].BondDipoles[l].B;
            DipoleMagnitudeB=Components[typeB].BondDipoleMagnitude[l];
  
            posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
            posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
            scalingB1=Cations[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
            scalingB2=Cations[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
            dipoleB.x=posB2.x-posB1.x;
            dipoleB.y=posB2.y-posB1.y;
            dipoleB.z=posB2.z-posB1.z;
            posB.x=posB1.x+0.5*dipoleB.x;
            posB.y=posB1.y+0.5*dipoleB.y;
            posB.z=posB1.z+0.5*dipoleB.z;
            rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
            length=sqrt(rk2);
            temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
            dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
 
            if(New)
            { 
              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(rr<CutOffChargeBondDipoleSquared)
                UCationChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
            }
 
            if(Old)
            { 
              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(rr<CutOffChargeBondDipoleSquared)
                UCationChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
            }
          }
        }
      }
    }
  
    for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
    {
      A1=Components[comp].BondDipoles[j].A;
      A2=Components[comp].BondDipoles[j].B;
      DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];
 
      if(New)
      { 
        posA1=TrialPosition[CurrentSystem][A1];
        posA2=TrialPosition[CurrentSystem][A2];
        scalingA1=CFChargeScaling[A1];
        scalingA2=CFChargeScaling[A2];
        dipoleA_new.x=posA2.x-posA1.x;
        dipoleA_new.y=posA2.y-posA1.y;
        dipoleA_new.z=posA2.z-posA1.z;
        posA_new.x=posA1.x+0.5*dipoleA_new.x;
        posA_new.y=posA1.y+0.5*dipoleA_new.y;
        posA_new.z=posA1.z+0.5*dipoleA_new.z;
        ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
      }
 
      if(Old)
      { 
        posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
        posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
        scalingA1=Adsorbates[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
        scalingA2=Adsorbates[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
        dipoleA_old.x=posA2.x-posA1.x;
        dipoleA_old.y=posA2.y-posA1.y;
        dipoleA_old.z=posA2.z-posA1.z;
        posA_old.x=posA1.x+0.5*dipoleA_old.x;
        posA_old.y=posA1.y+0.5*dipoleA_old.y;
        posA_old.z=posA1.z+0.5*dipoleA_old.z;
        ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
      }
  
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(m!=k)
          {
            for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
            {
              type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[type].HasCharges)
              {
                posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                chargeB=Adsorbates[CurrentSystem][k].Atoms[l].CFChargeScalingParameter*Adsorbates[CurrentSystem][k].Atoms[l].Charge;
 
                if(New)
                { 
                  dr.x=posA_new.x-posB.x;
                  dr.y=posA_new.y-posB.y;
                  dr.z=posA_new.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
                  if(rr<CutOffChargeBondDipoleSquared)
                    UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
                }
   
                if(Old)
                { 
                  dr.x=posA_old.x-posB.x;
                  dr.y=posA_old.y-posB.y;
                  dr.z=posA_old.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
   
                  if(rr<CutOffChargeBondDipoleSquared)
                    UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
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
          type=Cations[CurrentSystem][k].Atoms[l].Type;
          if(PseudoAtoms[type].HasCharges)
          {
            posB=Cations[CurrentSystem][k].Atoms[l].Position;
            chargeB=Cations[CurrentSystem][k].Atoms[l].CFChargeScalingParameter*Cations[CurrentSystem][k].Atoms[l].Charge;
 
            if(New)
            { 
              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(rr<CutOffChargeBondDipoleSquared)
                UCationChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
            }
   
            if(Old)
            { 
              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(rr<CutOffChargeBondDipoleSquared)
                UCationChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
            }
          }
        }
      }
    }
  }
  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes inter-molecular energy differences for a single adsorbate molecule:             *
 *            | 1) New=TRUE, Old=FAlse, computes the bonddipole-bonddipole energy of the new molecule    *
 *            |    (positions in 'TrialPosition' and the molecule is of type 'comp'                      *
 *            | 2) New=FALSE, Old=TRUE, computes the bonddipole-bonddipole energy of molecule 'm' of     *
 *            |    type 'comp'                                                                           *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in bonddipole-bonddipole energy of the    *
 *            |    molecule 'm' with the new positions in 'TrialPosition' and the old positions          *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the bonddipole-bonddipole energy of the new positions     *
 *            |      stored in 'TrialPosition'                                                           *
 *            | False: whether or not to compute the bonddipole-bonddipole energy of the old positions   *
 *            |        of molecule 'm'                                                                   * 
 * Note       |                                                                                          *
 * Used       | int TranslationMoveAdsorbate(void)                                                       *
 *            | int RandomTranslationMoveAdsorbate(void)                                                 *
 *            | int RotationMoveAdsorbate(void)                                                          *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CFCBSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateInterBondDipoleBondDipoleEnergyDifferenceAdsorbate(int m,int comp,int New,int Old)
{
  int j,k,l;
  int A1,A2,B1,B2;
  int TypeB;
  REAL r;
  VECTOR posA_new,posA_old,posB,posA1,posA2,posB1,posB2,dr;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL ri2,rk2,rr,length,temp;
  VECTOR dipoleA_new,dipoleA_old,dipoleB;
  REAL scalingA1,scalingA2,scalingB1,scalingB2;
  int ncell;

  UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;
  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
    {
      A1=Components[comp].BondDipoles[j].A;
      A2=Components[comp].BondDipoles[j].B;
      DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];
 
      if(New)
      { 
        posA1=TrialPosition[CurrentSystem][A1];
        posA2=TrialPosition[CurrentSystem][A2];
        scalingA1=CFChargeScaling[A1];
        scalingA2=CFChargeScaling[A2];
        dipoleA_new.x=posA2.x-posA1.x;
        dipoleA_new.y=posA2.y-posA1.y;
        dipoleA_new.z=posA2.z-posA1.z;
        posA_new.x=posA1.x+0.5*dipoleA_new.x;
        posA_new.y=posA1.y+0.5*dipoleA_new.y;
        posA_new.z=posA1.z+0.5*dipoleA_new.z;
        ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
      }
 
      if(Old)
      { 
        posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
        posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
        scalingA1=Adsorbates[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
        scalingA2=Adsorbates[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
        dipoleA_old.x=posA2.x-posA1.x;
        dipoleA_old.y=posA2.y-posA1.y;
        dipoleA_old.z=posA2.z-posA1.z;
        posA_old.x=posA1.x+0.5*dipoleA_old.x;
        posA_old.y=posA1.y+0.5*dipoleA_old.y;
        posA_old.z=posA1.z+0.5*dipoleA_old.z;
        ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
      }
  
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
              scalingB1=Adsorbates[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
              scalingB2=Adsorbates[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
              dipoleB.x=posB2.x-posB1.x;
              dipoleB.y=posB2.y-posB1.y;
              dipoleB.z=posB2.z-posB1.z;
              posB.x=posB1.x+0.5*dipoleB.x;
              posB.y=posB1.y+0.5*dipoleB.y;
              posB.z=posB1.z+0.5*dipoleB.z;
              rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
              length=sqrt(rk2);
              temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
              dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

              if(New)
              { 
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                { 
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                  if(rr<CutOffBondDipoleBondDipoleSquared)
                    UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));
                }
              }
    
              if(Old)
              { 
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                { 
                  dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                  if(rr<CutOffBondDipoleBondDipoleSquared)
                    UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
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
          scalingB1=Cations[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
          scalingB2=Cations[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          if(New)
          { 
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++) 
            {
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
              if(rr<CutOffBondDipoleBondDipoleSquared)
                UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));
            }
          }
    
          if(New)
          { 
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++) 
            {
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
              if(rr<CutOffBondDipoleBondDipoleSquared)
                UCationBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
    {
      A1=Components[comp].BondDipoles[j].A;
      A2=Components[comp].BondDipoles[j].B;
      DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];
 
      if(New)
      { 
        posA1=TrialPosition[CurrentSystem][A1];
        posA2=TrialPosition[CurrentSystem][A2];
        scalingA1=CFChargeScaling[A1];
        scalingA2=CFChargeScaling[A2];
        dipoleA_new.x=posA2.x-posA1.x;
        dipoleA_new.y=posA2.y-posA1.y;
        dipoleA_new.z=posA2.z-posA1.z;
        posA_new.x=posA1.x+0.5*dipoleA_new.x;
        posA_new.y=posA1.y+0.5*dipoleA_new.y;
        posA_new.z=posA1.z+0.5*dipoleA_new.z;
        ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
      } 
 
      if(Old)
      { 
        posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
        posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
        scalingA1=Adsorbates[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
        scalingA2=Adsorbates[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
        dipoleA_old.x=posA2.x-posA1.x;
        dipoleA_old.y=posA2.y-posA1.y;
        dipoleA_old.z=posA2.z-posA1.z;
        posA_old.x=posA1.x+0.5*dipoleA_old.x;
        posA_old.y=posA1.y+0.5*dipoleA_old.y;
        posA_old.z=posA1.z+0.5*dipoleA_old.z;
        ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
      }
  
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
              scalingB1=Adsorbates[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
              scalingB2=Adsorbates[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
              dipoleB.x=posB2.x-posB1.x;
              dipoleB.y=posB2.y-posB1.y;
              dipoleB.z=posB2.z-posB1.z;
              posB.x=posB1.x+0.5*dipoleB.x;
              posB.y=posB1.y+0.5*dipoleB.y;
              posB.z=posB1.z+0.5*dipoleB.z;
              rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
              length=sqrt(rk2);
              temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
              dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
 
              if(New)
              { 
                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
                if(rr<CutOffBondDipoleBondDipoleSquared)
                  UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));
              }
 
              if(Old)
              { 
                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
                if(rr<CutOffBondDipoleBondDipoleSquared)
                  UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
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
          scalingB1=Cations[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
          scalingB2=Cations[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          if(New)
          {  
            dr.x=posA_new.x-posB.x;
            dr.y=posA_new.y-posB.y;
            dr.z=posA_new.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
            if(rr<CutOffBondDipoleBondDipoleSquared)
              UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));
          }
 
          if(Old)
          { 
            dr.x=posA_old.x-posB.x;
            dr.y=posA_old.y-posB.y;
            dr.z=posA_old.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
            if(rr<CutOffBondDipoleBondDipoleSquared)
              UCationBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
          }
        }
      }
    }
  }
  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateInterChargeBondDipoleEnergyDifferenceCation                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes inter-molecular energy differences for a single cation molecule:                *
 *            | 1) New=TRUE, Old=FAlse, computes the charge-bonddipole energy of the new molecule        *
 *            |    (positions in 'TrialPosition' and the molecule is of type 'comp'                      *
 *            | 2) New=FALSE, Old=TRUE, computes the charge-bonddipole energy of molecule 'm' of type    *
 *            |    'comp'                                                                                *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in charge-bonddipole energy of the        *
 *            |    molecule 'm' with the new positions in 'TrialPosition' and the old positions          *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the charge-bonddipole energy of the new positions stored  *
 *            |      in 'TrialPosition'                                                                  *
 *            | False: whether or not to compute the charge-bonddipole energy of the old positions of    *
 *            |        molecule 'm'                                                                      * 
 * Note       |                                                                                          *
 * Used       | int TranslationMoveCation(void)                                                          *
 *            | int RandomTranslationMoveCation(void)                                                    *
 *            | int RotationMoveCation(void)                                                             *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CFCBSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateInterChargeBondDipoleEnergyDifferenceCation(int m,int comp,int New,int Old)
{
  int j,k,l;
  int A1,A2,B1,B2;
  int type,typeB;
  VECTOR posA_new,posA_old,posB,posA1,posA2,posB1,posB2,dr;
  REAL rr,ri2,rk2,temp,length,chargeB,chargeA_new,chargeA_old;
  VECTOR dipoleA_new,dipoleA_old,dipoleB;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL scalingA1,scalingA2,scalingB1,scalingB2;
  int ncell;

  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]=0.0;
  UCationChargeBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Components[comp].NumberOfAtoms;j++)
    {
      if(New)
      {
        posA_new=TrialPosition[CurrentSystem][j];
        chargeA_new=CFChargeScaling[j]*PseudoAtoms[Components[comp].Type[j]].Charge1;
      }

      if(Old)
      {
        posA_old=Cations[CurrentSystem][m].Atoms[j].Position;
        chargeA_old=Cations[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Cations[CurrentSystem][m].Atoms[j].Charge;
      }
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            typeB=Cations[CurrentSystem][k].Type;
            for(l=0;l<Components[typeB].NumberOfBondDipoles;l++)
            {
              B1=Components[typeB].BondDipoles[l].A;
              B2=Components[typeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[typeB].BondDipoleMagnitude[l];
  
              posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
              posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
              scalingB1=Cations[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
              scalingB2=Cations[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
              dipoleB.x=posB2.x-posB1.x;
              dipoleB.y=posB2.y-posB1.y;
              dipoleB.z=posB2.z-posB1.z;
              posB.x=posB1.x+0.5*dipoleB.x;
              posB.y=posB1.y+0.5*dipoleB.y;
              posB.z=posB1.z+0.5*dipoleB.z;
              rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
              length=sqrt(rk2);
              temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
              dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

              if(New)
              { 
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                { 
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                  if(rr<CutOffChargeBondDipoleSquared)
                    UCationChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
                }
              }
    
              if(Old)
              { 
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                { 
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                  if(rr<CutOffChargeBondDipoleSquared)
                    UCationChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
                }
              } 
            }
          }
        }
      }
  
      for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
      {
        typeB=Adsorbates[CurrentSystem][k].Type;
        for(l=0;l<Components[typeB].NumberOfBondDipoles;l++)
        {
          B1=Components[typeB].BondDipoles[l].A;
          B2=Components[typeB].BondDipoles[l].B;
          DipoleMagnitudeB=Components[typeB].BondDipoleMagnitude[l];
  
          posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
          posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
          scalingB1=Adsorbates[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
          scalingB2=Adsorbates[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
 
          if(New)
          { 
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            { 
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
              if(rr<CutOffChargeBondDipoleSquared)
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
            }
          }
    
          if(Old)
          { 
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
            { 
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
              if(rr<CutOffChargeBondDipoleSquared)
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
            }
          }
        }
      }
    }
  
    for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
    {
      A1=Components[comp].BondDipoles[j].A;
      A2=Components[comp].BondDipoles[j].B;
      DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];
 
      if(New)
      { 
        posA1=TrialPosition[CurrentSystem][A1];
        posA2=TrialPosition[CurrentSystem][A2];
        scalingA1=CFChargeScaling[A1];
        scalingA2=CFChargeScaling[A2];
        dipoleA_new.x=posA2.x-posA1.x;
        dipoleA_new.y=posA2.y-posA1.y;
        dipoleA_new.z=posA2.z-posA1.z;
        posA_new.x=posA1.x+0.5*dipoleA_new.x;
        posA_new.y=posA1.y+0.5*dipoleA_new.y;
        posA_new.z=posA1.z+0.5*dipoleA_new.z;
        ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
      }
 
      if(Old)
      { 
        posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
        posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
        scalingA1=Cations[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
        scalingA2=Cations[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
        dipoleA_old.x=posA2.x-posA1.x;
        dipoleA_old.y=posA2.y-posA1.y;
        dipoleA_old.z=posA2.z-posA1.z;
        posA_old.x=posA1.x+0.5*dipoleA_old.x;
        posA_old.y=posA1.y+0.5*dipoleA_old.y;
        posA_old.z=posA1.z+0.5*dipoleA_old.z;
        ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
      }
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(m!=k)
          {
            for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
            {
              type=Cations[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[type].HasCharges)
              {
                posB=Cations[CurrentSystem][k].Atoms[l].Position;
                chargeB=Cations[CurrentSystem][k].Atoms[l].CFChargeScalingParameter*Cations[CurrentSystem][k].Atoms[l].Charge;
 
                if(New) 
                {
                  for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                  { 
                    dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                    dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                    dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                    dr=ApplyReplicaBoundaryCondition(dr);
                    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                    if(rr<CutOffChargeBondDipoleSquared)
                      UCationChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
                  }
                }
      
                if(Old) 
                {
                  for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                  { 
                    dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                    dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                    dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                    dr=ApplyReplicaBoundaryCondition(dr);
                    rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
     
                    if(rr<CutOffChargeBondDipoleSquared)
                      UCationChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
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
          type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
          if(PseudoAtoms[type].HasCharges)
          {
            posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
            chargeB=Adsorbates[CurrentSystem][k].Atoms[l].CFChargeScalingParameter*Adsorbates[CurrentSystem][k].Atoms[l].Charge;
 
            if(New)
            { 
              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              { 
                dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                if(rr<CutOffChargeBondDipoleSquared)
                  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
              }
            }
      
            if(Old)
            { 
              for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
              { 
                dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                if(rr<CutOffChargeBondDipoleSquared)
                  UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
              }
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<Components[comp].NumberOfAtoms;j++)
    {
      if(New)
      {
        posA_new=TrialPosition[CurrentSystem][j];
        chargeA_new=CFChargeScaling[j]*PseudoAtoms[Components[comp].Type[j]].Charge1;
      }

      if(Old)
      {
        posA_old=Cations[CurrentSystem][m].Atoms[j].Position;
        chargeA_old=Cations[CurrentSystem][m].Atoms[j].CFChargeScalingParameter*Cations[CurrentSystem][m].Atoms[j].Charge;
      }
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(k!=m)
          {
            typeB=Cations[CurrentSystem][k].Type;
            for(l=0;l<Components[typeB].NumberOfBondDipoles;l++)
            {
              B1=Components[typeB].BondDipoles[l].A;
              B2=Components[typeB].BondDipoles[l].B;
              DipoleMagnitudeB=Components[typeB].BondDipoleMagnitude[l];
  
              posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
              posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
              scalingB1=Cations[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
              scalingB2=Cations[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
              dipoleB.x=posB2.x-posB1.x;
              dipoleB.y=posB2.y-posB1.y;
              dipoleB.z=posB2.z-posB1.z;
              posB.x=posB1.x+0.5*dipoleB.x;
              posB.y=posB1.y+0.5*dipoleB.y;
              posB.z=posB1.z+0.5*dipoleB.z;
              rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
              length=sqrt(rk2);
              temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
              dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
 
              if(New)
              { 
                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
                if(rr<CutOffChargeBondDipoleSquared)
                  UCationChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
              }
 
              if(Old)
              { 
                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
                if(rr<CutOffChargeBondDipoleSquared)
                  UCationChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
              }
            }
          }
        }
  
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          typeB=Adsorbates[CurrentSystem][k].Type;
          for(l=0;l<Components[typeB].NumberOfBondDipoles;l++)
          {
            B1=Components[typeB].BondDipoles[l].A;
            B2=Components[typeB].BondDipoles[l].B;
            DipoleMagnitudeB=Components[typeB].BondDipoleMagnitude[l];
  
            posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
            posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
            scalingB1=Adsorbates[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
            scalingB2=Adsorbates[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
            dipoleB.x=posB2.x-posB1.x;
            dipoleB.y=posB2.y-posB1.y;
            dipoleB.z=posB2.z-posB1.z;
            posB.x=posB1.x+0.5*dipoleB.x;
            posB.y=posB1.y+0.5*dipoleB.y;
            posB.z=posB1.z+0.5*dipoleB.z;
            rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
            length=sqrt(rk2);
            temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
            dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
 
            if(New)
            { 
              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(rr<CutOffChargeBondDipoleSquared)
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeA_new,dipoleB,dr,sqrt(rr));
            }
 
            if(Old)
            { 
              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(rr<CutOffChargeBondDipoleSquared)
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeA_old,dipoleB,dr,sqrt(rr));
            }
          }
        }
      }
    }
  
    for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
    {
      A1=Components[comp].BondDipoles[j].A;
      A2=Components[comp].BondDipoles[j].B;
      DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];
 
      if(New)
      { 
        posA1=TrialPosition[CurrentSystem][A1];
        posA2=TrialPosition[CurrentSystem][A2];
        scalingA1=CFChargeScaling[A1];
        scalingA2=CFChargeScaling[A2];
        dipoleA_new.x=posA2.x-posA1.x;
        dipoleA_new.y=posA2.y-posA1.y;
        dipoleA_new.z=posA2.z-posA1.z;
        posA_new.x=posA1.x+0.5*dipoleA_new.x;
        posA_new.y=posA1.y+0.5*dipoleA_new.y;
        posA_new.z=posA1.z+0.5*dipoleA_new.z;
        ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
      }
 
      if(Old)
      { 
        posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
        posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
        scalingA1=Cations[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
        scalingA2=Cations[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
        dipoleA_old.x=posA2.x-posA1.x;
        dipoleA_old.y=posA2.y-posA1.y;
        dipoleA_old.z=posA2.z-posA1.z;
        posA_old.x=posA1.x+0.5*dipoleA_old.x;
        posA_old.y=posA1.y+0.5*dipoleA_old.y;
        posA_old.z=posA1.z+0.5*dipoleA_old.z;
        ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
      }
  
      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(m!=k)
          {
            for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
            {
              type=Cations[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[type].HasCharges)
              {
                posB=Cations[CurrentSystem][k].Atoms[l].Position;
                chargeB=Cations[CurrentSystem][k].Atoms[l].CFChargeScalingParameter*Cations[CurrentSystem][k].Atoms[l].Charge;
 
                if(New)
                { 
                  dr.x=posA_new.x-posB.x;
                  dr.y=posA_new.y-posB.y;
                  dr.z=posA_new.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
                  if(rr<CutOffChargeBondDipoleSquared)
                    UCationChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
                }
   
                if(Old)
                { 
                  dr.x=posA_old.x-posB.x;
                  dr.y=posA_old.y-posB.y;
                  dr.z=posA_old.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
   
                  if(rr<CutOffChargeBondDipoleSquared)
                    UCationChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
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
          type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
          if(PseudoAtoms[type].HasCharges)
          {
            posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
            chargeB=Adsorbates[CurrentSystem][k].Atoms[l].CFChargeScalingParameter*Adsorbates[CurrentSystem][k].Atoms[l].Charge;
 
            if(New)
            { 
              dr.x=posA_new.x-posB.x;
              dr.y=posA_new.y-posB.y;
              dr.z=posA_new.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(rr<CutOffChargeBondDipoleSquared)
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]-=PotentialValueChargeBondDipole(chargeB,dipoleA_new,dr,sqrt(rr));
            }
   
            if(Old)
            { 
              dr.x=posA_old.x-posB.x;
              dr.y=posA_old.y-posB.y;
              dr.z=posA_old.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
              if(rr<CutOffChargeBondDipoleSquared)
                UAdsorbateChargeBondDipoleRealDelta[CurrentSystem]+=PotentialValueChargeBondDipole(chargeB,dipoleA_old,dr,sqrt(rr));
            }
          }
        }
      }
    }
  }
  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateInterBondDipoleBondDipoleEnergyDifferenceCation                                 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes inter-molecular energy differences for a single cation molecule:                *
 *            | 1) New=TRUE, Old=FAlse, computes the bonddipole-bonddipole energy of the new molecule    *
 *            |    (positions in 'TrialPosition' and the molecule is of type 'comp'                      *
 *            | 2) New=FALSE, Old=TRUE, computes the bonddipole-bonddipole energy of molecule 'm' of     *
 *            |    type 'comp'                                                                           *
 *            | 3) New=TRUE, Old=TRUE, computes the difference in bonddipole-bonddipole energy of the    *
 *            |    molecule 'm' with the new positions in 'TrialPosition' and the old positions          *
 * Parameters | m: the molecule-id                                                                       *
 *            | comp: the component-type                                                                 *
 *            | New: whether or not to compute the bonddipole-bonddipole energy of the new positions     *
 *            |      stored in 'TrialPosition'                                                           *
 *            | False: whether or not to compute the bonddipole-bonddipole energy of the old positions   *
 *            |        of molecule 'm'                                                                   * 
 * Note       |                                                                                          *
 * Used       | int TranslationMoveCation(void)                                                          *
 *            | int RandomTranslationMoveCation(void)                                                    *
 *            | int RotationMoveCation(void)                                                             *
 *            | int CFSwapLambaMove(void)                                                                *
 *            | int CFCBSwapLambaMove(void)                                                              *
 *********************************************************************************************************/

int CalculateInterBondDipoleBondDipoleEnergyDifferenceCation(int m,int comp,int New,int Old)
{
  int j,k,l;
  int A1,A2,B1,B2;
  int TypeB;
  REAL r;
  VECTOR posA_new,posA_old,posB,posA1,posA2,posB1,posB2,dr;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL ri2,rk2,rr,length,temp;
  VECTOR dipoleA_new,dipoleA_old,dipoleB;
  REAL scalingA1,scalingA2,scalingB1,scalingB2;
  int ncell;

  UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;
  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  if(UseReplicas[CurrentSystem])
  {
    for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
    {
      A1=Components[comp].BondDipoles[j].A;
      A2=Components[comp].BondDipoles[j].B;
      DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];
 
      if(New)
      { 
        posA1=TrialPosition[CurrentSystem][A1];
        posA2=TrialPosition[CurrentSystem][A2];
        scalingA1=CFChargeScaling[A1];
        scalingA2=CFChargeScaling[A2];
        dipoleA_new.x=posA2.x-posA1.x;
        dipoleA_new.y=posA2.y-posA1.y;
        dipoleA_new.z=posA2.z-posA1.z;
        posA_new.x=posA1.x+0.5*dipoleA_new.x;
        posA_new.y=posA1.y+0.5*dipoleA_new.y;
        posA_new.z=posA1.z+0.5*dipoleA_new.z;
        ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
      }
 
      if(Old)
      { 
        posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
        posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
        scalingA1=Cations[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
        scalingA2=Cations[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
        dipoleA_old.x=posA2.x-posA1.x;
        dipoleA_old.y=posA2.y-posA1.y;
        dipoleA_old.z=posA2.z-posA1.z;
        posA_old.x=posA1.x+0.5*dipoleA_old.x;
        posA_old.y=posA1.y+0.5*dipoleA_old.y;
        posA_old.z=posA1.z+0.5*dipoleA_old.z;
        ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
      }
  
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
              scalingB1=Cations[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
              scalingB2=Cations[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
              dipoleB.x=posB2.x-posB1.x;
              dipoleB.y=posB2.y-posB1.y;
              dipoleB.z=posB2.z-posB1.z;
              posB.x=posB1.x+0.5*dipoleB.x;
              posB.y=posB1.y+0.5*dipoleB.y;
              posB.z=posB1.z+0.5*dipoleB.z;
              rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
              length=sqrt(rk2);
              temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
              dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

              if(New)
              { 
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                { 
                  dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                  if(rr<CutOffBondDipoleBondDipoleSquared)
                    UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));
                }
              }
    
              if(Old)
              { 
                for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
                { 
                  dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
                  dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
                  dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
                  if(rr<CutOffBondDipoleBondDipoleSquared)
                    UCationBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
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
          scalingB1=Adsorbates[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
          scalingB2=Adsorbates[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          if(New)
          { 
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++) 
            {
              dr.x=posA_new.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_new.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_new.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
              if(rr<CutOffBondDipoleBondDipoleSquared)
                UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));
            }
          }
    
          if(New)
          { 
            for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++) 
            {
              dr.x=posA_old.x-(posB.x+ReplicaShift[ncell].x);
              dr.y=posA_old.y-(posB.y+ReplicaShift[ncell].y);
              dr.z=posA_old.z-(posB.z+ReplicaShift[ncell].z);
              dr=ApplyReplicaBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
    
              if(rr<CutOffBondDipoleBondDipoleSquared)
                UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
            }
          }
        }
      }
    }
  }
  else
  {
    for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
    {
      A1=Components[comp].BondDipoles[j].A;
      A2=Components[comp].BondDipoles[j].B;
      DipoleMagnitudeA=Components[comp].BondDipoleMagnitude[j];
 
      if(New)
      { 
        posA1=TrialPosition[CurrentSystem][A1];
        posA2=TrialPosition[CurrentSystem][A2];
        scalingA1=CFChargeScaling[A1];
        scalingA2=CFChargeScaling[A2];
        dipoleA_new.x=posA2.x-posA1.x;
        dipoleA_new.y=posA2.y-posA1.y;
        dipoleA_new.z=posA2.z-posA1.z;
        posA_new.x=posA1.x+0.5*dipoleA_new.x;
        posA_new.y=posA1.y+0.5*dipoleA_new.y;
        posA_new.z=posA1.z+0.5*dipoleA_new.z;
        ri2=SQR(dipoleA_new.x)+SQR(dipoleA_new.y)+SQR(dipoleA_new.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_new.x*=temp; dipoleA_new.y*=temp; dipoleA_new.z*=temp;
      } 
 
      if(Old)
      { 
        posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
        posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
        scalingA1=Cations[CurrentSystem][m].Atoms[A1].CFChargeScalingParameter;
        scalingA2=Cations[CurrentSystem][m].Atoms[A2].CFChargeScalingParameter;
        dipoleA_old.x=posA2.x-posA1.x;
        dipoleA_old.y=posA2.y-posA1.y;
        dipoleA_old.z=posA2.z-posA1.z;
        posA_old.x=posA1.x+0.5*dipoleA_old.x;
        posA_old.y=posA1.y+0.5*dipoleA_old.y;
        posA_old.z=posA1.z+0.5*dipoleA_old.z;
        ri2=SQR(dipoleA_old.x)+SQR(dipoleA_old.y)+SQR(dipoleA_old.z);
        length=sqrt(ri2);
        temp=0.5*(scalingA1+scalingA2)*DipoleMagnitudeA/length;
        dipoleA_old.x*=temp; dipoleA_old.y*=temp; dipoleA_old.z*=temp;
      }
  
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
              scalingB1=Cations[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
              scalingB2=Cations[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
              dipoleB.x=posB2.x-posB1.x;
              dipoleB.y=posB2.y-posB1.y;
              dipoleB.z=posB2.z-posB1.z;
              posB.x=posB1.x+0.5*dipoleB.x;
              posB.y=posB1.y+0.5*dipoleB.y;
              posB.z=posB1.z+0.5*dipoleB.z;
              rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
              length=sqrt(rk2);
              temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
              dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;
 
              if(New)
              { 
                dr.x=posA_new.x-posB.x;
                dr.y=posA_new.y-posB.y;
                dr.z=posA_new.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
                if(rr<CutOffBondDipoleBondDipoleSquared)
                  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));
              }
 
              if(Old)
              { 
                dr.x=posA_old.x-posB.x;
                dr.y=posA_old.y-posB.y;
                dr.z=posA_old.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
                if(rr<CutOffBondDipoleBondDipoleSquared)
                  UCationBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
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
          scalingB1=Adsorbates[CurrentSystem][k].Atoms[B1].CFChargeScalingParameter;
          scalingB2=Adsorbates[CurrentSystem][k].Atoms[B2].CFChargeScalingParameter;
          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
          length=sqrt(rk2);
          temp=0.5*(scalingB1+scalingB2)*DipoleMagnitudeB/length;
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          if(New)
          {  
            dr.x=posA_new.x-posB.x;
            dr.y=posA_new.y-posB.y;
            dr.z=posA_new.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
            if(rr<CutOffBondDipoleBondDipoleSquared)
              UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]+=PotentialValueBondDipoleBondDipole(dipoleA_new,dipoleB,dr,sqrt(r));
          }
 
          if(Old)
          { 
            dr.x=posA_old.x-posB.x;
            dr.y=posA_old.y-posB.y;
            dr.z=posA_old.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
  
            if(rr<CutOffBondDipoleBondDipoleSquared)
              UAdsorbateBondDipoleBondDipoleRealDelta[CurrentSystem]-=PotentialValueBondDipoleBondDipole(dipoleA_old,dipoleB,dr,sqrt(r));
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
              UVDWDelta+=PotentialValue(typeA,typeB,rr,1.0);

            posA=Positions[k];
            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UVDWDelta-=PotentialValue(typeA,typeB,rr,1.0);
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
              UVDWDelta+=PotentialValue(typeA,typeB,rr,1.0);

            posA=Positions[k];
            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffVDWSquared)
              UVDWDelta-=PotentialValue(typeA,typeB,rr,1.0);
          }
        }
      }
    }
  }
  return UVDWDelta;
}
