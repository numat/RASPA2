/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'cbmc.c' is part of RASPA.

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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cbmc.h"
#include "constants.h"
#include "simulation.h"
#include "utils.h"
#include "molecule.h"
#include "potentials.h"
#include "internal_energy.h"
#include "internal_force.h"
#include "inter_energy.h"
#include "inter_force.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "framework.h"
#include "ewald.h"
#include "utils.h"
#include "mc_moves.h"

int BiasingMethod;

static REAL RosenbluthNew;             // the new Rosenbluth weight (grow)
static REAL RosenbluthOld;             // the old Rosenbluth weight (retrace)
static REAL *RosenbluthTorsion;
static REAL RosenBluthFactorFirstBead;

static int CurrentBead;       // the current bead (new beads are grown from this one)
static int CurrentGroup;      // the current group

static int NumberOfPreviousBeads;  // the number of previous beads
static int PreviousBead;           // the first previous bead
static int PreviousGroup;          // the previous group

static int NumberOfBeadsToBePlaced;              // the numer of beads to be placed
static int *BeadsToBePlaced;

int NumberOfBeadsAlreadyPlaced;              // the number of beads already placed
int *BeadsAlreadyPlaced; // a list of atom-ids

static int NumberOfBranches;                                      // the number of branches
static int *NumberOfBranchAtoms;              // the number of atoms per branch
static int **BranchAtoms; // a list of atom-ids per branch

VECTOR FirstBeadPosition;
VECTOR **NewPosition;
VECTOR **NewVelocity;
VECTOR **NewForce;
VECTOR **TrialPositions;

int OVERLAP;
static REAL *Bfac;
static int *Trover;

int NumberOfTrialPositions;
int NumberOfTrialPositionsForTheFirstBead;
int NumberOfTrialMovesPerOpenBead;
int NumberOfTrialPositionsTorsion;
REAL TargetAccRatioSmallMCScheme;
REAL EnergyOverlapCriteria;
REAL MinimumRosenbluthFactor;

int MaxNumberOfBeads;
int MaxNumberOfBonds;
int MaxNumberOfBondDipoles;
int MaxNumberOfBends;
int MaxNumberOfBendBends;
int MaxNumberOfInversionBends;
int MaxNumberOfUreyBradleys;
int MaxNumberOfTorsions;
int MaxNumberOfImproperTorsions;
int MaxNumberOfOutOfPlanes;
int MaxNumberOfBondBonds;
int MaxNumberOfBondBends;
int MaxNumberOfBondTorsions;
int MaxNumberOfBendTorsions;
int MaxNumberOfIntraVDW;
int MaxNumberOfIntraChargeCharge;
int MaxNumberOfIntraChargeBondDipole;
int MaxNumberOfIntraBondDipoleBondDipole;

int MaxNumberOfPermanentDipoles;

static int NumberOfBonds;            // the number of bonds around the current-bead
static int *Bonds;                   // a list of bond-ids

static int NumberOfBondDipoles;      // the number of bond-dipoles
static int *BondDipoles;             // a list of bond-dipole-ids

static int NumberOfBends;            // the number of bonds around the current-bead
static int *Bends;                   // a list of bend-ids

static int NumberOfBendBends;        // the number of bonds around the current-bead
static int *BendBends;               // a list of bend-ids

static int NumberOfUreyBradleys;     // the number of Uray-Bradleys around the current-bead
static int *UreyBradleys;            // a list of Urey-Bradley-ids

static int NumberOfInversionBends;   // the number of inversion-bends around the current-bead
static int *InversionBends;          // a list of inversion-bend-ids

static int NumberOfTorsions;         // the number of torsions around the current-bead
static int *Torsions;                // a list of torsion-ids

static int NumberOfImproperTorsions; // the number of torsions around the current-bead
static int *ImproperTorsions;        // a list of torsion-ids

static int NumberOfBondBonds;        // the number of torsions around the current-bead
static int *BondBonds;               // a list of torsion-ids

static int NumberOfBondBends;        // the number of torsions around the current-bead
static int *BondBends;               // a list of torsion-ids

static int NumberOfBondTorsions;     // the number of torsions around the current-bead
static int *BondTorsions;            // a list of torsion-ids

static int NumberOfBendTorsions;     // the number of torsions around the current-bead
static int *BendTorsions;            // a list of torsion-ids

static int NumberOfIntraVDW;         // the number of intra-vdws around the current-bead
static int *VDW;                     // a list of VDW-ids

static int NumberOfIntraChargeCharge;  // the number of intra-ChargeCharge around the current-bead
static int *IntraChargeCharge;         // a list if intra-ChargeCharge ids

static int NumberOfIntraChargeBondDipole;  // the number of intra-Charge-Bondipoles around the current-bead
static int *IntraChargeBondDipole;         // a list if intra-Charge-BondDipoles ids

static int NumberOfIntraBondDipoleBondDipole;  // the number of intra-BondDipoles around the current-bead
static int *IntraBondDipoleBondDipole;         // a list if intra-BondDipoles ids

static int NumberOfIntraPermanentDipoles;
static int *IntraPermanentDipoles;

// grow energies
REAL *UBondNew;
REAL *UBendNew;
REAL *UBendBendNew;
REAL *UUreyBradleyNew;
REAL *UInversionBendNew;
REAL *UTorsionNew;
REAL *UImproperTorsionNew;
REAL *UBondBondNew;
REAL *UBondBendNew;
REAL *UBondTorsionNew;
REAL *UBendTorsionNew;
REAL *UIntraVDWNew;
REAL *UIntraChargeChargeNew;
REAL *UIntraChargeBondDipoleNew;
REAL *UIntraBondDipoleBondDipoleNew;

REAL *UAdsorbateVDWNew;
REAL *UCationVDWNew;
REAL *UHostVDWNew;
REAL *UAdsorbateChargeChargeNew;
REAL *UCationChargeChargeNew;
REAL *UHostChargeChargeNew;
REAL *UAdsorbateChargeBondDipoleNew;
REAL *UCationChargeBondDipoleNew;
REAL *UHostChargeBondDipoleNew;
REAL *UAdsorbateChargePermanentDipoleNew;
REAL *UCationChargePermanentDipoleNew;
REAL *UHostChargePermanentDipoleNew;
REAL *UAdsorbateBondDipoleBondDipoleNew;
REAL *UCationBondDipoleBondDipoleNew;
REAL *UHostBondDipoleBondDipoleNew;
REAL *UAdsorbateBondDipolePermanentDipoleNew;
REAL *UCationBondDipolePermanentDipoleNew;
REAL *UHostBondDipolePermanentDipoleNew;
REAL *UAdsorbatePermanentDipolePermanentDipoleNew;
REAL *UCationPermanentDipolePermanentDipoleNew;
REAL *UHostPermanentDipolePermanentDipoleNew;

// retrace energies
REAL *UBondOld;
REAL *UBendOld;
REAL *UBendBendOld;
REAL *UUreyBradleyOld;
REAL *UInversionBendOld;
REAL *UTorsionOld;
REAL *UImproperTorsionOld;
REAL *UBondBondOld;
REAL *UBondBendOld;
REAL *UBondTorsionOld;
REAL *UBendTorsionOld;
REAL *UIntraVDWOld;
REAL *UIntraChargeChargeOld;
REAL *UIntraChargeBondDipoleOld;
REAL *UIntraBondDipoleBondDipoleOld;

REAL *UAdsorbateVDWOld;
REAL *UCationVDWOld;
REAL *UHostVDWOld;
REAL *UAdsorbateChargeChargeOld;
REAL *UCationChargeChargeOld;
REAL *UHostChargeChargeOld;
REAL *UAdsorbateChargeBondDipoleOld;
REAL *UCationChargeBondDipoleOld;
REAL *UHostChargeBondDipoleOld;
REAL *UAdsorbateChargePermanentDipoleOld;
REAL *UCationChargePermanentDipoleOld;
REAL *UHostChargePermanentDipoleOld;
REAL *UAdsorbateBondDipoleBondDipoleOld;
REAL *UCationBondDipoleBondDipoleOld;
REAL *UHostBondDipoleBondDipoleOld;
REAL *UAdsorbateBondDipolePermanentDipoleOld;
REAL *UCationBondDipolePermanentDipoleOld;
REAL *UHostBondDipolePermanentDipoleOld;
REAL *UAdsorbatePermanentDipolePermanentDipoleOld;
REAL *UCationPermanentDipolePermanentDipoleOld;
REAL *UHostPermanentDipolePermanentDipoleOld;

// trial energies
static REAL *UBondTrial;
static REAL *UBendTrial;
static REAL *UBendBendTrial;
static REAL *UUreyBradleyTrial;
static REAL *UInversionBendTrial;
static REAL *UTorsionTrial;
static REAL *UImproperTorsionTrial;
static REAL *UBondBondTrial;
static REAL *UBondBendTrial;
static REAL *UBondTorsionTrial;
static REAL *UBendTorsionTrial;
static REAL *UIntraVDWTrial;
static REAL *UIntraChargeChargeTrial;
static REAL *UIntraChargeBondDipoleTrial;
static REAL *UIntraBondDipoleBondDipoleTrial;

static REAL *UHostVDWTrial;
static REAL *UAdsorbateVDWTrial;
static REAL *UCationVDWTrial;
static REAL *UHostChargeChargeTrial;
static REAL *UAdsorbateChargeChargeTrial;
static REAL *UCationChargeChargeTrial;
static REAL *UHostChargeBondDipoleTrial;
static REAL *UAdsorbateChargeBondDipoleTrial;
static REAL *UCationChargeBondDipoleTrial;
static REAL *UHostChargePermanentDipoleTrial;
static REAL *UAdsorbateChargePermanentDipoleTrial;
static REAL *UCationChargePermanentDipoleTrial;
static REAL *UHostBondDipoleBondDipoleTrial;
static REAL *UAdsorbateBondDipoleBondDipoleTrial;
static REAL *UCationBondDipoleBondDipoleTrial;
static REAL *UHostBondDipolePermanentDipoleTrial;
static REAL *UAdsorbateBondDipolePermanentDipoleTrial;
static REAL *UCationBondDipolePermanentDipoleTrial;
static REAL *UHostPermanentDipolePermanentDipoleTrial;
static REAL *UAdsorbatePermanentDipolePermanentDipoleTrial;
static REAL *UCationPermanentDipolePermanentDipoleTrial;

static REAL EnergyHostVDWFirstBead;
static REAL EnergyAdsorbateVDWFirstBead;
static REAL EnergyCationVDWFirstBead;
static REAL EnergyHostChargeChargeFirstBead;
static REAL EnergyAdsorbateChargeChargeFirstBead;
static REAL EnergyCationChargeChargeFirstBead;
static REAL EnergyHostChargeBondDipoleFirstBead;
static REAL EnergyAdsorbateChargeBondDipoleFirstBead;
static REAL EnergyCationChargeBondDipoleFirstBead;
static REAL EnergyHostChargePermanentDipoleFirstBead;
static REAL EnergyAdsorbateChargePermanentDipoleFirstBead;
static REAL EnergyCationChargePermanentDipoleFirstBead;
static REAL EnergyHostBondDipoleBondDipoleFirstBead;
static REAL EnergyAdsorbateBondDipoleBondDipoleFirstBead;
static REAL EnergyCationBondDipoleBondDipoleFirstBead;
static REAL EnergyHostBondDipolePermanentDipoleFirstBead;
static REAL EnergyAdsorbateBondDipolePermanentDipoleFirstBead;
static REAL EnergyCationBondDipolePermanentDipoleFirstBead;
static REAL EnergyHostPermanentDipolePermanentDipoleFirstBead;
static REAL EnergyAdsorbatePermanentDipolePermanentDipoleFirstBead;
static REAL EnergyCationPermanentDipolePermanentDipoleFirstBead;

static REAL *ShiftedBoltzmannFactors;
static VECTOR *Trial;
static VECTOR *TrialAnisotropicPositionRetrace;

static REAL *EnergiesHostVDW;
static REAL *EnergiesAdsorbateVDW;
static REAL *EnergiesCationVDW;
static REAL *EnergiesHostChargeCharge;
static REAL *EnergiesAdsorbateChargeCharge;
static REAL *EnergiesCationChargeCharge;
static REAL *EnergiesHostChargeBondDipole;
static REAL *EnergiesAdsorbateChargeBondDipole;
static REAL *EnergiesCationChargeBondDipole;
static REAL *EnergiesHostChargePermanentDipole;
static REAL *EnergiesAdsorbateChargePermanentDipole;
static REAL *EnergiesCationChargePermanentDipole;
static REAL *EnergiesHostBondDipoleBondDipole;
static REAL *EnergiesAdsorbateBondDipoleBondDipole;
static REAL *EnergiesCationBondDipoleBondDipole;
static REAL *EnergiesHostBondDipolePermanentDipole;
static REAL *EnergiesAdsorbateBondDipolePermanentDipole;
static REAL *EnergiesCationBondDipolePermanentDipole;
static REAL *EnergiesHostPermanentDipolePermanentDipole;
static REAL *EnergiesAdsorbatePermanentDipolePermanentDipole;
static REAL *EnergiesCationPermanentDipolePermanentDipole;

static VECTOR* store;
static REAL* angle;
static VECTOR* perpen;
static VECTOR* cord;
static REAL* bond_length;

static REAL *enbend;
static REAL *eobend;
static REAL *enbendbend;
static REAL *eobendbend;
static REAL *enbond;
static REAL *eobond;
static REAL *enureybradley;
static REAL *eoureybradley;
static REAL *eninversionbend;
static REAL *eoinversionbend;
static REAL *enimpropertorsion;
static REAL *eoimpropertorsion;

static int *lplace;
static int *lchose;
static int *beadn;
static int *PossibleCurrentBeads;
static int **MoleculeTodoConnectivity;
static int **MoleculeConnectivity;

static int HandleFirstBead(int Code);
static void Interactions(void);
static int GenerateTrialOrientationsSimpleSphere(int Old);
static int GenerateTrialOrientationsMCScheme(int Old);
static int ComputeExternalEnergies(void);

static int Rosen(void);
static int RosenOld(void);

static REAL ComputeSumRosenbluthWeight(REAL *Bfac,int *Trover,int NumberOfPositions);
static REAL ComputeNormalizedRosenbluthWeight(REAL *Bfac,int *Trover,int NumberOfPositions);
static int SelectTrialPosition(REAL *Bfac,int *Trover,int NumberOfPositions);


/*********************************************************************************************************
 * Name       | CalculateAnisotropicTrialPositions                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Calculate the anisotropic positions for the trial positions.                             *
 * Parameters | -                                                                                        *
 * Used in    | 'GrowMolecule' and 'RetraceMolecule'.                                                    *
 *********************************************************************************************************/

void CalculateAnisotropicTrialPositions(int TypeMolA,VECTOR *TrialPosition,VECTOR *TrialAnisotropicPosition)
{
  int k,A,B;
  int typeA;
  VECTOR vec,posA,posL,posR;
  VECTOR drAL,drAR;
  REAL ral,rar;
  REAL length,d;

  for(k=0;k<Components[TypeMolA].NumberOfAtoms;k++)
  {
    typeA=Components[TypeMolA].Type[k];
    posA=TrialPosition[k];

    if(PseudoAtoms[typeA].AnisotropicCorrection)
    {
      d=PseudoAtoms[typeA].AnisotropicDisplacement;
      switch(Components[TypeMolA].Connectivity[k])
      {
        case 0:
          break;
        case 1:
          A=Components[TypeMolA].ConnectivityList[k][0];
          posL=TrialPosition[A];
          vec.x=posA.x-posL.x;
          vec.y=posA.y-posL.y;
          vec.z=posA.z-posL.z;
          if(PseudoAtoms[typeA].AnisotropicType==ABSOLUTE)
          {
            length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
            vec.x/=length;
            vec.y/=length;
            vec.z/=length;
          }
          vec.x*=d;
          vec.y*=d;
          vec.z*=d;
          posA.x+=vec.x;
          posA.y+=vec.y;
          posA.z+=vec.z;
          break;
        case 2:
          switch(Components[TypeMolA].AnisotropicType)
          {
            case ANISOTROPIC_BISECTION:
              A=Components[TypeMolA].ConnectivityList[k][0];
              B=Components[TypeMolA].ConnectivityList[k][1];
              posL=TrialPosition[A];
              posR=TrialPosition[B];
              drAL.x=posA.x-posL.x;
              drAL.y=posA.y-posL.y;
              drAL.z=posA.z-posL.z;
              ral=sqrt(SQR(drAL.x)+SQR(drAL.y)+SQR(drAL.z));
              drAR.x=posA.x-posR.x;
              drAR.y=posA.y-posR.y;
              drAR.z=posA.z-posR.z;
              rar=sqrt(SQR(drAR.x)+SQR(drAR.y)+SQR(drAR.z));
              vec.x=(rar*drAL.x+ral*drAR.x)/(ral+rar);
              vec.y=(rar*drAL.y+ral*drAR.y)/(ral+rar);
              vec.z=(rar*drAL.z+ral*drAR.z)/(ral+rar);
              length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
              vec.x/=length;
              vec.y/=length;
              vec.z/=length;
              vec.x*=d;
              vec.y*=d;
              vec.z*=d;
              posA.x+=vec.x;
              posA.y+=vec.y;
              posA.z+=vec.z;
              break;
            default:
            case ANISOTROPIC_MID_POINT:
              A=Components[TypeMolA].ConnectivityList[k][0];
              B=Components[TypeMolA].ConnectivityList[k][1];
              posL=TrialPosition[A];
              posR=TrialPosition[B];
              vec.x=posA.x-0.5*(posL.x+posR.x);
              vec.y=posA.y-0.5*(posL.y+posR.y);
              vec.z=posA.z-0.5*(posL.z+posR.z);
              length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
              vec.x/=length;
              vec.y/=length;
              vec.z/=length;
              vec.x*=d;
              vec.y*=d;
              vec.z*=d;
              posA.x+=vec.x;
              posA.y+=vec.y;
              posA.z+=vec.z;
              break;
          }
          break;
        default:
          printf("ERROR: undefined anisotropic atom with connecvity: %d, in routine 'CalculateAnisotropicTrialPositions' (cbmc.c)\n",
                 Components[TypeMolA].Connectivity[k]);
          exit(0);
          break;
      }
    }
    TrialAnisotropicPosition[k]=posA;
  }
}


/*********************************************************************************************************
 * Name       | CalculateAnisotropicTrialPositions                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Calculate the sum of the Boltzmann factors                                               *
 * Parameters | -                                                                                        *
 * Used in    | 'GrowMolecule' and 'RetraceMolecule'.                                                    *
 *********************************************************************************************************/

REAL ComputeSumRosenbluthWeight(REAL *Bfac,int *Trover,int NumberOfPositions)
{
  int i;
  REAL SumBoltzmannFactors;

  SumBoltzmannFactors=0.0;
  for(i=0;i<NumberOfPositions;i++)
    if(!Trover[i])
      SumBoltzmannFactors+=exp(Bfac[i]);
  return SumBoltzmannFactors;

}

/*********************************************************************************************************
 * Name       | ComputeNormalizedRosenbluthWeight                                                        *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Calculate the normalized Rosenbluth weight                                               *
 * Parameters | -                                                                                        *
 * Used in    | 'GenerateTrialOrientationsMCScheme', 'Rosen', 'RosenOld'                                 *
 *********************************************************************************************************/

REAL ComputeNormalizedRosenbluthWeight(REAL *Bfac,int *Trover,int NumberOfPositions)
{
  int i;
  REAL SumBoltzmannFactors;

  SumBoltzmannFactors=0.0;
  for(i=0;i<NumberOfPositions;i++)
    if(!Trover[i])
      SumBoltzmannFactors+=exp(Bfac[i]);
  return SumBoltzmannFactors/(REAL)NumberOfPositions;

}

/*********************************************************************************************************
 * Name       | SelectTrialPosition                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Selects a trial position with the proper weight                                          *
 * Parameters | -                                                                                        *
 * Used in    | 'HandleFirstBead', 'GenerateTrialOrientationsMCScheme', 'Rosen'                          *
 *********************************************************************************************************/

int SelectTrialPosition(REAL *Bfac,int *Trover,int NumberOfPositions)
{
  int i;
  int selected,success;
  REAL SumShiftedBoltzmannFactors;
  REAL ws,cumw,largest_value;

  // find largest value of the energy
  largest_value=0.0;
  success=FALSE;
  for(i=0;i<NumberOfPositions;i++)
  {
    if(!Trover[i])
    {
      if(!success)
      {
        success=TRUE;
        largest_value=Bfac[i];
      }
      else
        if(Bfac[i]>largest_value) largest_value=Bfac[i];
    }
  }

  // shift the energies down
  SumShiftedBoltzmannFactors=0.0;
  for(i=0;i<NumberOfPositions;i++)
  {
    if(!Trover[i])
    {
      ShiftedBoltzmannFactors[i]=exp(Bfac[i]-largest_value);
      SumShiftedBoltzmannFactors+=ShiftedBoltzmannFactors[i];
    }
    else
      ShiftedBoltzmannFactors[i]=0.0;
  }

  // select the Boltzmann factor
  selected=0;
  cumw=ShiftedBoltzmannFactors[0];

  ws=RandomNumber()*SumShiftedBoltzmannFactors;
  while(cumw<ws)
    cumw+=ShiftedBoltzmannFactors[++selected];

  return selected;
}


enum {CBMC_INSERTION,CBMC_PARTIAL_INSERTION,CBMC_DELETION,CBMC_RETRACE_REINSERTION};

//     1. New Configuration
//        Generate NUMBER_OF_TRIAL_POSITIONS_FOR_THE_FIRST_BEAD random positions, 
//        one is chosen according to its Boltzmann weight;
//        used by Swap()/SwapAdd()/Regrow()/Initia()
//     2. Conventional Stuff: calculating weight for molecule with chosen firstbead-position
//        Calculate the Boltzmann factor in (FirstBeadPosition); used by PartialReinsertion()/IdentityChangeMove()/GibbsIdentityChangeMove()
//     3. Old Configuration
//        Generate NUMBER_OF_TRIAL_POSITIONS_FOR_THE_FIRST_BEAD-1 random positions; 
//        the first position is (FirstBeadPosition);
//        used by Swap() and Swaprem()
//     4. Old Configuration
//        Use the (NUMBER_OF_TRIAL_POSITIONS_FOR_THE_FIRST_BEAD-1) beads that were generated 
//        using the growth of a trial chain.
//        The first trial position is the old configuration. This is only used in the subroutine Reinsertion();
//        See Esselink/Loyens/Smit, Phys.Rev.E.,1995,51,2.
int HandleFirstBead(int Code)
{
  int i,type,start;
  int NumberOfFirstPositions;
  REAL EnergyHostVDW,EnergyAdsorbateVDW,EnergyCationVDW;
  REAL EnergyHostChargeCharge,EnergyAdsorbateChargeCharge,EnergyCationChargeCharge;
  REAL EnergyHostChargeBondDipole,EnergyAdsorbateChargeBondDipole,EnergyCationChargeBondDipole;
  POINT posA,s;
  static double ERR;

  OVERLAP=TRUE;
  RosenBluthFactorFirstBead=0.0;
  if((Code==2)||(Code==4)) NumberOfFirstPositions=1;
  else NumberOfFirstPositions=NumberOfTrialPositionsForTheFirstBead;

  for(i=0;i<NumberOfFirstPositions;i++)
  {
    if((i==0)&&(Code!=1))
      Trial[i]=FirstBeadPosition;
    else
    {
      FirstBeadPosition.x=FirstBeadPosition.y=FirstBeadPosition.z=0.0;
      if(Components[CurrentComponent].RestrictMoves)
      {
        // generate random fractional point between 0.0 and 1.0 and retry when it is not inside the specified range
        s.x=s.y=s.z=0.0;
        do
        {
          switch(Dimension)
          {
            case 3:
              s.z=RandomNumber();
            case 2:
              s.y=RandomNumber();
            case 1:
              s.x=RandomNumber();
              break;
          }

          // convert to a the Cartesian position in the simulation box
          FirstBeadPosition.x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
          FirstBeadPosition.y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
          FirstBeadPosition.z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;
        } while(!ValidCartesianPoint(CurrentComponent,FirstBeadPosition));
      }
      else
      {
        // generate random fractional point between 0 and 1
        s.x=s.y=s.z=0.0;
        switch(Dimension)
        {
          case 3:
            s.z=RandomNumber();
          case 2:
            s.y=RandomNumber();
          case 1:
            s.x=RandomNumber();
            break;
        }
        FirstBeadPosition=ConvertFromABCtoXYZ(s);
      }
      Trial[i]=FirstBeadPosition;
    }

    posA=Trial[i];
    start=Components[CurrentComponent].StartingBead;
    type=Components[CurrentComponent].Type[start];

    EnergyHostVDW=EnergyAdsorbateVDW=EnergyCationVDW=0.0;
    EnergyHostChargeCharge=EnergyAdsorbateChargeCharge=EnergyCationChargeCharge=0.0;
    EnergyHostChargeBondDipole=EnergyAdsorbateChargeBondDipole=EnergyCationChargeBondDipole=0.0;

    // calculate energies
    EnergyHostVDW=CalculateFrameworkVDWEnergyAtPosition(posA,type);
    CalculateFrameworkChargeEnergyAtPosition(posA,type,&EnergyHostChargeCharge,&EnergyHostChargeBondDipole);

    // compute VDW energy with adsorbates if no omit of adsorbate-adsorbate or the current molecule is a cation
    if(Components[CurrentComponent].ExtraFrameworkMolecule||(!OmitAdsorbateAdsorbateVDWInteractions))
      EnergyAdsorbateVDW=CalculateInterVDWEnergyAdsorbateAtPosition(posA,type,CurrentAdsorbateMolecule);

    // compute Coulomb energy with adsorbates if no omit of adsorbate-adsorbate or the current molecule is a cation
    if(Components[CurrentComponent].ExtraFrameworkMolecule||(!OmitAdsorbateAdsorbateCoulombInteractions))
      CalculateInterChargeEnergyAdsorbateAtPosition(posA,type,&EnergyAdsorbateChargeCharge,&EnergyAdsorbateChargeBondDipole,CurrentAdsorbateMolecule);

    // compute VDW energy with cations if no omit of cation-cation or the current molecule is an adsorbate
    if((!Components[CurrentComponent].ExtraFrameworkMolecule)||(!OmitCationCationVDWInteractions))
      EnergyCationVDW=CalculateInterVDWEnergyCationAtPosition(posA,type,CurrentCationMolecule);

    // compute Coulomb energy with cations if no omit of cation-cation or the current molecule is an adsorbate
    if((!Components[CurrentComponent].ExtraFrameworkMolecule)||(!OmitCationCationCoulombInteractions))
      CalculateInterChargeEnergyCationAtPosition(posA,type,&EnergyCationChargeCharge,&EnergyCationChargeBondDipole,CurrentCationMolecule);

    if((EnergyHostVDW>=EnergyOverlapCriteria)||(EnergyHostChargeCharge>=EnergyOverlapCriteria))
    {
      Trover[i]=TRUE;
      OVERLAP=OVERLAP&TRUE;
    }
    else
    {
      Trover[i]=FALSE;
      OVERLAP=OVERLAP&FALSE;

      EnergiesHostVDW[i]=EnergyHostVDW;
      EnergiesAdsorbateVDW[i]=EnergyAdsorbateVDW;
      EnergiesCationVDW[i]=EnergyCationVDW;

      EnergiesHostChargeCharge[i]=EnergyHostChargeCharge;
      EnergiesAdsorbateChargeCharge[i]=EnergyAdsorbateChargeCharge;
      EnergiesCationChargeCharge[i]=EnergyCationChargeCharge;

      EnergiesHostChargeBondDipole[i]=EnergyHostChargeBondDipole;
      EnergiesAdsorbateChargeBondDipole[i]=EnergyAdsorbateChargeBondDipole;
      EnergiesCationChargeBondDipole[i]=EnergyCationChargeBondDipole;

      EnergiesHostBondDipoleBondDipole[i]=0.0;
      EnergiesAdsorbateBondDipoleBondDipole[i]=0.0;
      EnergiesCationBondDipoleBondDipole[i]=0.0;

      if(BiasingMethod==LJ_BIASING)
        Bfac[i]=-Beta[CurrentSystem]*(EnergyHostVDW+EnergyAdsorbateVDW+EnergyCationVDW);
      else
      {
        Bfac[i]=-Beta[CurrentSystem]*(EnergyHostVDW+EnergyCationVDW+EnergyAdsorbateVDW+
                       EnergyHostChargeCharge+EnergyAdsorbateChargeCharge+EnergyCationChargeCharge+
                       EnergyHostChargeBondDipole+EnergyAdsorbateChargeBondDipole+EnergyCationChargeBondDipole);
      }
    }
  }
  if (OVERLAP) return 0;  // return when all trial-positions overlap

  i=0;
  RosenBluthFactorFirstBead=ComputeSumRosenbluthWeight(Bfac,Trover,NumberOfFirstPositions);
  if(Code==1)
  {
    i=SelectTrialPosition(Bfac,Trover,NumberOfFirstPositions);
    ERR=RosenBluthFactorFirstBead-exp(Bfac[i]);
  }
  else if(Code==4)
  {
    i=0;
    RosenBluthFactorFirstBead+=ERR;
  }

  // update positions and energies
  FirstBeadPosition=Trial[i];

  EnergyHostVDWFirstBead=EnergiesHostVDW[i];
  EnergyAdsorbateVDWFirstBead=EnergiesAdsorbateVDW[i];
  EnergyCationVDWFirstBead=EnergiesCationVDW[i];

  EnergyHostChargeChargeFirstBead=EnergiesHostChargeCharge[i];
  EnergyAdsorbateChargeChargeFirstBead=EnergiesAdsorbateChargeCharge[i];
  EnergyCationChargeChargeFirstBead=EnergiesCationChargeCharge[i];

  EnergyHostChargeBondDipoleFirstBead=EnergiesHostChargeBondDipole[i];
  EnergyAdsorbateChargeBondDipoleFirstBead=EnergiesAdsorbateChargeBondDipole[i];
  EnergyCationChargeBondDipoleFirstBead=EnergiesCationChargeBondDipole[i];

  EnergyHostBondDipoleBondDipoleFirstBead=0.0;
  EnergyAdsorbateBondDipoleBondDipoleFirstBead=0.0;
  EnergyCationBondDipoleBondDipoleFirstBead=0.0;

  if (Code!=2) RosenBluthFactorFirstBead/=(REAL)NumberOfTrialPositionsForTheFirstBead;

  // check if the Rosenbluth factor is reasonable; only for a new configuration
  if ((Code==1)&&(RosenBluthFactorFirstBead<MinimumRosenbluthFactor)) OVERLAP=TRUE;

  return 0;
}

int GetBondNumber(int A, int B)
{
  int i;

  for(i=0;i<Components[CurrentComponent].NumberOfBonds;i++)
    if(((Components[CurrentComponent].Bonds[i].A==A)&&(Components[CurrentComponent].Bonds[i].B==B))||
       ((Components[CurrentComponent].Bonds[i].B==A)&&(Components[CurrentComponent].Bonds[i].A==B)))
       return i;
  
  return -1;
}

REAL GetReferenceBondLength(int A,int B)
{
  int i;
  REAL *parms;

  for(i=0;i<Components[CurrentComponent].NumberOfBonds;i++)
    if(((Components[CurrentComponent].Bonds[i].A==A)&&(Components[CurrentComponent].Bonds[i].B==B))||
       ((Components[CurrentComponent].Bonds[i].B==A)&&(Components[CurrentComponent].Bonds[i].A==B)))
         break;

   if(i<0) return 1.5;

  parms=Components[CurrentComponent].BondArguments[i];
  switch(Components[CurrentComponent].BondType[i])
  {
    case HARMONIC_BOND:
      // 0.5*p0*SQR(r-p1);
      // ===============================================
      // p_0/k_B [K/A^2]   force constant
      // p_1     [A]       reference bond distance
      return parms[1];
    case MORSE_BOND:
      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
      // ===============================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference bond distance
      return parms[2];
    case LJ_12_6_BOND:
      // A/r_ij^12-B/r_ij^6
      // ===============================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
    case LENNARD_JONES_BOND:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A]
    case RESTRAINED_HARMONIC_BOND:
      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2     [A]
      return 1.5;
    case QUARTIC_BOND:
      // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
      // ===========================================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2/k_B [K/A^3]
      // p_3/k_B [K/A^4]
      return parms[1];
    case CFF_QUARTIC_BOND:
      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
      // ===============================================
      // p_0/k_B [K/A^2]
      // p_1     [A]
      // p_2/k_B [K/A^3]
      // p_3/k_B [K/A^4]
      return parms[1];
    case BUCKINGHAM_BOND:
      // p_0*exp(-p_1 r)-p_2/r^6
      // ===============================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2     [K A^6]
      return 1.5;
    case MM3_BOND:
      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
      // =================================================================
      // p_0     [mdyne/A molecule]
      // p_1     [A]
      return parms[1];
    case RIGID_BOND:
    case FIXED_BOND:
      return parms[0];
    default:
      printf("Unkown bond-length\n");
      return 1.5;
  }
}

// Calculates which interactions are needed for the next bead
void Interactions(void)
{
  int i,j,k;
  int A,B,C,D;
  int index;

  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    lplace[i]=FALSE;
    lchose[i]=FALSE;
  }

  for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
    lchose[BeadsAlreadyPlaced[i]]=TRUE;

  for(i=0;i<NumberOfBeadsToBePlaced;i++)
  {
    lplace[BeadsToBePlaced[i]]=TRUE;
    lchose[BeadsToBePlaced[i]]=TRUE;
  }

  // find all bonds from CurrentBead to a bead-to-place
  // skipping the ones which belong to a rigid-structure
  NumberOfBonds=0;
  for(i=0;i<NumberOfBeadsToBePlaced;i++)
  {
    A=MIN2(BeadsToBePlaced[i],CurrentBead);
    B=MAX2(BeadsToBePlaced[i],CurrentBead);
    for(j=0;j<Components[CurrentComponent].NumberOfBonds;j++)
      if(((A==(Components[CurrentComponent].Bonds[j].A))&&
          (B==(Components[CurrentComponent].Bonds[j].B)))||
         ((A==(Components[CurrentComponent].Bonds[j].B))&&
          (B==(Components[CurrentComponent].Bonds[j].A))))
         {
           if(Components[CurrentComponent].BondType[j]!=RIGID_BOND)
             Bonds[NumberOfBonds++]=j;
         }
  }

  NumberOfBondDipoles=0;
  for(i=0;i<NumberOfBeadsToBePlaced;i++)
  {
    A=MIN2(BeadsToBePlaced[i],CurrentBead);
    B=MAX2(BeadsToBePlaced[i],CurrentBead);
    for(k=0;k<Components[CurrentComponent].NumberOfBondDipoles;k++)
      if(((A==(Components[CurrentComponent].BondDipoles[k].A))&&
          (B==(Components[CurrentComponent].BondDipoles[k].B)))||
         ((A==(Components[CurrentComponent].BondDipoles[k].B))&&
          (B==(Components[CurrentComponent].BondDipoles[k].A))))
             BondDipoles[NumberOfBondDipoles++]=k;
  }
  for(i=0;i<NumberOfBeadsToBePlaced;i++)
  {
    A=BeadsToBePlaced[i];
    for(j=i+1;j<NumberOfBeadsToBePlaced;j++)
    {
      B=BeadsToBePlaced[j];
      for(k=0;k<Components[CurrentComponent].NumberOfBondDipoles;k++)
        if(((A==(Components[CurrentComponent].BondDipoles[k].A))&&
            (B==(Components[CurrentComponent].BondDipoles[k].B)))||
           ((A==(Components[CurrentComponent].BondDipoles[k].B))&&
            (B==(Components[CurrentComponent].BondDipoles[k].A))))
               BondDipoles[NumberOfBondDipoles++]=k;
    }
  }

  NumberOfBends=0;
  for(j=0;j<Components[CurrentComponent].NumberOfBends;j++)
  {
    if(lchose[Components[CurrentComponent].Bends[j].A]&&
       lchose[Components[CurrentComponent].Bends[j].B]&&
       lchose[Components[CurrentComponent].Bends[j].C]&&
      (lplace[Components[CurrentComponent].Bends[j].A]||
       lplace[Components[CurrentComponent].Bends[j].B]||
       lplace[Components[CurrentComponent].Bends[j].C]))
      Bends[NumberOfBends++]=j;
  }

  NumberOfBendBends=0;
  for(j=0;j<Components[CurrentComponent].NumberOfBendBends;j++)
  {
    if(lchose[Components[CurrentComponent].BendBends[j].A]&&
       lchose[Components[CurrentComponent].BendBends[j].B]&&
       lchose[Components[CurrentComponent].BendBends[j].C]&&
       lchose[Components[CurrentComponent].BendBends[j].D]&&
      (lplace[Components[CurrentComponent].BendBends[j].A]||
       lplace[Components[CurrentComponent].BendBends[j].B]||
       lplace[Components[CurrentComponent].BendBends[j].C]||
       lplace[Components[CurrentComponent].BendBends[j].D]))
      BendBends[NumberOfBendBends++]=j;
  }

  NumberOfUreyBradleys=0;
  for(j=0;j<Components[CurrentComponent].NumberOfUreyBradleys;j++)
  {
    if(lchose[Components[CurrentComponent].UreyBradleys[j].A]&&
       lchose[Components[CurrentComponent].UreyBradleys[j].B]&&
       lchose[Components[CurrentComponent].UreyBradleys[j].C]&&
      (lplace[Components[CurrentComponent].UreyBradleys[j].A]||
       lplace[Components[CurrentComponent].UreyBradleys[j].B]||
       lplace[Components[CurrentComponent].UreyBradleys[j].C]))
      UreyBradleys[NumberOfUreyBradleys++]=j;
  }

  NumberOfInversionBends=0;
  for(j=0;j<Components[CurrentComponent].NumberOfInversionBends;j++)
  {
    if(lchose[Components[CurrentComponent].InversionBends[j].A]&&
       lchose[Components[CurrentComponent].InversionBends[j].B]&&
       lchose[Components[CurrentComponent].InversionBends[j].C]&&
       lchose[Components[CurrentComponent].InversionBends[j].D]&&
      (lplace[Components[CurrentComponent].InversionBends[j].A]||
       lplace[Components[CurrentComponent].InversionBends[j].B]||
       lplace[Components[CurrentComponent].InversionBends[j].C]||
       lplace[Components[CurrentComponent].InversionBends[j].D]))
      InversionBends[NumberOfInversionBends++]=j;
  }

  NumberOfTorsions=0;
  for(j=0;j<Components[CurrentComponent].NumberOfTorsions;j++)
  {
    if(lchose[Components[CurrentComponent].Torsions[j].A]&&
       lchose[Components[CurrentComponent].Torsions[j].B]&&
       lchose[Components[CurrentComponent].Torsions[j].C]&&
       lchose[Components[CurrentComponent].Torsions[j].D]&&
      (lplace[Components[CurrentComponent].Torsions[j].A]||
       lplace[Components[CurrentComponent].Torsions[j].B]||
       lplace[Components[CurrentComponent].Torsions[j].C]||
       lplace[Components[CurrentComponent].Torsions[j].D]))
      Torsions[NumberOfTorsions++]=j;
  }

  NumberOfImproperTorsions=0;
  for(j=0;j<Components[CurrentComponent].NumberOfImproperTorsions;j++)
  {
    if(lchose[Components[CurrentComponent].ImproperTorsions[j].A]&&
       lchose[Components[CurrentComponent].ImproperTorsions[j].B]&&
       lchose[Components[CurrentComponent].ImproperTorsions[j].C]&&
       lchose[Components[CurrentComponent].ImproperTorsions[j].D]&&
      (lplace[Components[CurrentComponent].ImproperTorsions[j].A]||
       lplace[Components[CurrentComponent].ImproperTorsions[j].B]||
       lplace[Components[CurrentComponent].ImproperTorsions[j].C]||
       lplace[Components[CurrentComponent].ImproperTorsions[j].D]))
      ImproperTorsions[NumberOfImproperTorsions++]=j;
  }

  NumberOfBondBonds=0;
  for(j=0;j<Components[CurrentComponent].NumberOfBondBonds;j++)
  {
    if(lchose[Components[CurrentComponent].BondBonds[j].A]&&
       lchose[Components[CurrentComponent].BondBonds[j].B]&&
       lchose[Components[CurrentComponent].BondBonds[j].C]&&
      (lplace[Components[CurrentComponent].BondBonds[j].A]||
       lplace[Components[CurrentComponent].BondBonds[j].B]||
       lplace[Components[CurrentComponent].BondBonds[j].C]))
      BondBonds[NumberOfBondBonds++]=j;
  }

  NumberOfBondBends=0;
  for(j=0;j<Components[CurrentComponent].NumberOfBondBends;j++)
  {
    if(lchose[Components[CurrentComponent].BondBends[j].A]&&
       lchose[Components[CurrentComponent].BondBends[j].B]&&
       lchose[Components[CurrentComponent].BondBends[j].C]&&
      (lplace[Components[CurrentComponent].BondBends[j].A]||
       lplace[Components[CurrentComponent].BondBends[j].B]||
       lplace[Components[CurrentComponent].BondBends[j].C]))
      BondBends[NumberOfBondBends++]=j;
  }

  NumberOfBondTorsions=0;
  for(j=0;j<Components[CurrentComponent].NumberOfBondTorsions;j++)
  {
    if(lchose[Components[CurrentComponent].BondTorsions[j].A]&&
       lchose[Components[CurrentComponent].BondTorsions[j].B]&&
       lchose[Components[CurrentComponent].BondTorsions[j].C]&&
       lchose[Components[CurrentComponent].BondTorsions[j].D]&&
      (lplace[Components[CurrentComponent].BondTorsions[j].A]||
       lplace[Components[CurrentComponent].BondTorsions[j].B]||
       lplace[Components[CurrentComponent].BondTorsions[j].C]||
       lplace[Components[CurrentComponent].BondTorsions[j].D]))
      BondTorsions[NumberOfBondTorsions++]=j;
  }


  NumberOfBendTorsions=0;
  for(j=0;j<Components[CurrentComponent].NumberOfBendTorsions;j++)
  {
    if(lchose[Components[CurrentComponent].BendTorsions[j].A]&&
       lchose[Components[CurrentComponent].BendTorsions[j].B]&&
       lchose[Components[CurrentComponent].BendTorsions[j].C]&&
       lchose[Components[CurrentComponent].BendTorsions[j].D]&&
      (lplace[Components[CurrentComponent].BendTorsions[j].A]||
       lplace[Components[CurrentComponent].BendTorsions[j].B]||
       lplace[Components[CurrentComponent].BendTorsions[j].C]||
       lplace[Components[CurrentComponent].BendTorsions[j].D]))
      BendTorsions[NumberOfBendTorsions++]=j;
  }

  NumberOfIntraVDW=0;
  for(j=0;j<Components[CurrentComponent].NumberOfIntraVDW;j++)
  {
    if (lchose[Components[CurrentComponent].IntraVDW[j].A]&&
        lchose[Components[CurrentComponent].IntraVDW[j].B]&&
       (lplace[Components[CurrentComponent].IntraVDW[j].A]||
        lplace[Components[CurrentComponent].IntraVDW[j].B]))
      VDW[NumberOfIntraVDW++]=j;
  }

  NumberOfIntraChargeCharge=0;
  for(j=0;j<Components[CurrentComponent].NumberOfIntraChargeCharge;j++)
  {
    if (lchose[Components[CurrentComponent].IntraChargeCharge[j].A]&&
        lchose[Components[CurrentComponent].IntraChargeCharge[j].B]&&
       (lplace[Components[CurrentComponent].IntraChargeCharge[j].A]||
        lplace[Components[CurrentComponent].IntraChargeCharge[j].B]))
      IntraChargeCharge[NumberOfIntraChargeCharge++]=j;
  }

  NumberOfIntraChargeBondDipole=0;
  for(j=0;j<Components[CurrentComponent].NumberOfIntraChargeBondDipole;j++)
  {
    A=Components[CurrentComponent].IntraChargeBondDipole[j].A;
    index=Components[CurrentComponent].IntraChargeBondDipole[j].B;
    B=Components[CurrentComponent].BondDipoles[index].A;
    C=Components[CurrentComponent].BondDipoles[index].B;

    if (lchose[A]&&lchose[B]&&lchose[C]&&
       (lplace[A]||lplace[B]||lplace[C]))
      IntraChargeBondDipole[NumberOfIntraChargeBondDipole++]=j;
  }

  NumberOfIntraBondDipoleBondDipole=0;
  for(j=0;j<Components[CurrentComponent].NumberOfIntraBondDipoleBondDipole;j++)
  {
    index=Components[CurrentComponent].IntraBondDipoleBondDipole[j].A;
    A=Components[CurrentComponent].BondDipoles[index].A;
    B=Components[CurrentComponent].BondDipoles[index].B;
    index=Components[CurrentComponent].IntraBondDipoleBondDipole[j].B;
    C=Components[CurrentComponent].BondDipoles[index].A;
    D=Components[CurrentComponent].BondDipoles[index].B;

    if (lchose[A]&&lchose[B]&&lchose[C]&&lchose[D]&&
       (lplace[A]||lplace[B]||lplace[C]||lplace[D]))
      IntraBondDipoleBondDipole[NumberOfIntraBondDipoleBondDipole++]=j;
  }
}

int ComputeExternalEnergies(void)
{
  int i,j,Itrial,A,B,type,typeA,typeB;
  REAL rr,r;
  REAL EnergiesIntra,EnergiesIntraChargeCharge,EnergiesIntraChargeBondDipole,EnergiesIntraBondDipoleBondDipole;
  REAL EnergyHostVDW,EnergyAdsorbateVDW,EnergyCationVDW;
  REAL EnergyHostChargeCharge,EnergyAdsorbateChargeCharge,EnergyCationChargeCharge;
  REAL EnergyHostChargeBondDipole,EnergyAdsorbateChargeBondDipole,EnergyCationChargeBondDipole;
  REAL EnergyHostBondDipoleBondDipole,EnergyAdsorbateBondDipoleBondDipole,EnergyCationBondDipoleBondDipole;
  REAL EnergiesHostVDW,EnergiesAdsorbateVDW,EnergiesCationVDW;
  REAL EnergiesHostChargeCharge,EnergiesAdsorbateChargeCharge,EnergiesCationChargeCharge;
  REAL EnergiesHostChargeBondDipole,EnergiesAdsorbateChargeBondDipole,EnergiesCationChargeBondDipole;
  REAL EnergiesHostBondDipoleBondDipole,EnergiesAdsorbateBondDipoleBondDipole,EnergiesCationBondDipoleBondDipole;
  REAL DipoleMagnitudeA,DipoleMagnitudeB,cosA,cosB,cosAB;
  REAL Bt1,Bt2,temp,ChargeA,ChargeB;
  REAL Scaling,Magnitude;
  int index,A1,A2,B1,B2;
  VECTOR posA1,posA2,posB1,posB2,dipoleA,dipoleB,dr;

  int TRIAL_OVERLAP;
  POINT posA,posB,posAVDW;


  TRIAL_OVERLAP=FALSE;
  OVERLAP=TRUE;
      
  for(Itrial=0;Itrial<NumberOfTrialPositions;Itrial++)
  {
    EnergiesCationVDW=0.0;
    EnergiesAdsorbateVDW=0.0;
    EnergiesHostVDW=0.0;

    EnergiesHostChargeCharge=0.0;
    EnergiesAdsorbateChargeCharge=0.0;
    EnergiesCationChargeCharge=0.0;

    EnergiesHostChargeBondDipole=0.0;
    EnergiesAdsorbateChargeBondDipole=0.0;
    EnergiesCationChargeBondDipole=0.0;

    EnergiesHostBondDipoleBondDipole=0.0;
    EnergiesAdsorbateBondDipoleBondDipole=0.0;
    EnergiesCationBondDipoleBondDipole=0.0;

    for(j=0;j<NumberOfBeadsToBePlaced;j++)
    {
      TRIAL_OVERLAP=FALSE;
      if (BeadsToBePlaced[j]<Components[CurrentComponent].NumberOfAtoms)
      {
        posAVDW=posA=TrialPositions[Itrial][BeadsToBePlaced[j]];
        type=Components[CurrentComponent].Type[BeadsToBePlaced[j]];

        EnergyHostVDW=EnergyHostChargeCharge=EnergyHostChargeBondDipole=EnergyHostBondDipoleBondDipole=0.0;
        EnergyAdsorbateVDW=EnergyAdsorbateChargeCharge=EnergyAdsorbateChargeBondDipole=EnergyAdsorbateBondDipoleBondDipole=0.0;
        EnergyCationVDW=EnergyCationChargeCharge=EnergyCationChargeBondDipole=EnergyCationBondDipoleBondDipole=0.0;

        // Calculate External energy
        EnergyHostVDW=CalculateFrameworkVDWEnergyAtPosition(posAVDW,type);
        if(EnergyHostVDW>=EnergyOverlapCriteria)
        {
          TRIAL_OVERLAP=TRUE;
          break;
        }
        CalculateFrameworkChargeEnergyAtPosition(posA,type,&EnergyHostChargeCharge,&EnergyHostChargeBondDipole);

        EnergiesHostVDW+=EnergyHostVDW;
        EnergiesHostChargeCharge+=EnergyHostChargeCharge;
        EnergiesHostChargeBondDipole+=EnergyHostChargeBondDipole;
        EnergiesHostBondDipoleBondDipole+=EnergyHostBondDipoleBondDipole;

        // compute VDW energy with adsorbates if no omit of adsorbate-adsorbate or the current molecule is a cation
        if(Components[CurrentComponent].ExtraFrameworkMolecule||(!OmitAdsorbateAdsorbateVDWInteractions))
          EnergyAdsorbateVDW=CalculateInterVDWEnergyAdsorbateAtPosition(posAVDW,type,CurrentAdsorbateMolecule);

        // compute Coulomb energy with adsorbates if no omit of adsorbate-adsorbate or the current molecule is a cation
        if(Components[CurrentComponent].ExtraFrameworkMolecule||(!OmitAdsorbateAdsorbateCoulombInteractions))
          CalculateInterChargeEnergyAdsorbateAtPosition(posA,type,&EnergyAdsorbateChargeCharge,&EnergyAdsorbateChargeBondDipole,CurrentAdsorbateMolecule);

        // compute VDW energy with cations if no omit of cation-cation or the current molecule is an adsorbate
        if((!Components[CurrentComponent].ExtraFrameworkMolecule)||(!OmitCationCationVDWInteractions))
          EnergyCationVDW=CalculateInterVDWEnergyCationAtPosition(posAVDW,type,CurrentCationMolecule);

        // compute Coulomb energy with cations if no omit of cation-cation or the current molecule is an adsorbate
        if((!Components[CurrentComponent].ExtraFrameworkMolecule)||(!OmitCationCationCoulombInteractions))
          CalculateInterChargeEnergyCationAtPosition(posA,type,&EnergyCationChargeCharge,&EnergyCationChargeBondDipole,CurrentCationMolecule);

        EnergiesAdsorbateVDW+=EnergyAdsorbateVDW;
        EnergiesAdsorbateChargeCharge+=EnergyAdsorbateChargeCharge;
        EnergiesAdsorbateChargeBondDipole+=EnergyAdsorbateChargeBondDipole;
        EnergiesAdsorbateBondDipoleBondDipole+=EnergyAdsorbateBondDipoleBondDipole;

        EnergiesCationVDW+=EnergyCationVDW;
        EnergiesCationChargeCharge+=EnergyCationChargeCharge;
        EnergiesCationChargeBondDipole+=EnergyCationChargeBondDipole;
        EnergiesCationBondDipoleBondDipole+=EnergyCationBondDipoleBondDipole;
      }
    }
    for(j=0;j<NumberOfBondDipoles;j++)
    {
      A=Components[CurrentComponent].BondDipoles[BondDipoles[j]].A;
      posA=TrialPositions[Itrial][A];
      B=Components[CurrentComponent].BondDipoles[BondDipoles[j]].B;
      posB=TrialPositions[Itrial][B];
      Magnitude=Components[CurrentComponent].BondDipoleMagnitude[BondDipoles[j]];

      CalculateFrameworkBondDipoleEnergyAtPosition(posA,posB,Magnitude,&EnergyHostChargeBondDipole,&EnergyHostBondDipoleBondDipole);
      if(EnergyHostChargeBondDipole>=EnergyOverlapCriteria||EnergyHostBondDipoleBondDipole>=EnergyOverlapCriteria)
      {
         TRIAL_OVERLAP=TRUE;
         break;
      }
      EnergiesHostChargeBondDipole+=EnergyHostChargeBondDipole;
      EnergiesHostBondDipoleBondDipole+=EnergyHostBondDipoleBondDipole;

      // compute Coulomb energy with adsorbates if no omit of adsorbate-adsorbate or the current molecule is a cation
      if(Components[CurrentComponent].ExtraFrameworkMolecule||(!OmitAdsorbateAdsorbateCoulombInteractions))
        CalculateInterBondDipoleEnergyAdsorbateAtPosition(posA,posB,Magnitude,&EnergyAdsorbateChargeBondDipole,&EnergyAdsorbateBondDipoleBondDipole,CurrentAdsorbateMolecule);
      if(EnergyAdsorbateChargeBondDipole>=EnergyOverlapCriteria||EnergyAdsorbateBondDipoleBondDipole>=EnergyOverlapCriteria)
      {
         TRIAL_OVERLAP=TRUE;
         break;
      }
      EnergiesAdsorbateChargeBondDipole+=EnergyAdsorbateChargeBondDipole;
      EnergiesAdsorbateBondDipoleBondDipole+=EnergyAdsorbateBondDipoleBondDipole;

      // compute Coulomb energy with cations if no omit of cation-cation or the current molecule is an adsorbate
      if((!Components[CurrentComponent].ExtraFrameworkMolecule)||(!OmitCationCationCoulombInteractions))
        CalculateInterBondDipoleEnergyCationAtPosition(posA,posB,Magnitude,&EnergyCationChargeBondDipole,&EnergyCationBondDipoleBondDipole,CurrentCationMolecule);
      if(EnergyCationChargeBondDipole>=EnergyOverlapCriteria||EnergyCationBondDipoleBondDipole>=EnergyOverlapCriteria)
      {
         TRIAL_OVERLAP=TRUE;
         break;
      }
      EnergiesCationChargeBondDipole+=EnergyCationChargeBondDipole;
      EnergiesCationBondDipoleBondDipole+=EnergyCationBondDipoleBondDipole;
    }

    EnergiesIntra=0.0;
    EnergiesIntraChargeCharge=0.0;
    EnergiesIntraChargeBondDipole=0.0;
    EnergiesIntraBondDipoleBondDipole=0.0;

    if(TRIAL_OVERLAP==FALSE)
    {
      for(i=0;i<NumberOfIntraVDW;i++)
      {
        A=Components[CurrentComponent].IntraVDW[VDW[i]].A;
        B=Components[CurrentComponent].IntraVDW[VDW[i]].B;
        Scaling=Components[CurrentComponent].IntraVDWScaling[VDW[i]];
        rr=SQR(TrialPositions[Itrial][A].x-TrialPositions[Itrial][B].x)+
           SQR(TrialPositions[Itrial][A].y-TrialPositions[Itrial][B].y)+
           SQR(TrialPositions[Itrial][A].z-TrialPositions[Itrial][B].z);
        if(rr<CutOffVDWSquared)
        {
          typeA=Components[CurrentComponent].Type[A];
          typeB=Components[CurrentComponent].Type[B];
          EnergiesIntra+=Scaling*PotentialValue(typeA,typeB,rr);
        }
      }

      if(ChargeMethod!=NONE)
      {
        for(i=0;i<NumberOfIntraChargeCharge;i++)
        {
          A=Components[CurrentComponent].IntraChargeCharge[IntraChargeCharge[i]].A;
          B=Components[CurrentComponent].IntraChargeCharge[IntraChargeCharge[i]].B;
          Scaling=Components[CurrentComponent].IntraChargeChargeScaling[IntraChargeCharge[i]];
          r=sqrt(SQR(TrialPositions[Itrial][A].x-TrialPositions[Itrial][B].x)+
                 SQR(TrialPositions[Itrial][A].y-TrialPositions[Itrial][B].y)+
                 SQR(TrialPositions[Itrial][A].z-TrialPositions[Itrial][B].z));
          typeA=Components[CurrentComponent].Type[A];
          typeB=Components[CurrentComponent].Type[B];
          ChargeA=Components[CurrentComponent].Charge[A];
          ChargeB=Components[CurrentComponent].Charge[B];
          switch(ChargeMethod)
          {
            case NONE:
              break;
            case SHIFTED_COULOMB:
              EnergiesIntraChargeCharge+=Scaling*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/r;
              break;
            case TRUNCATED_COULOMB:
            case EWALD:
            default:
              EnergiesIntraChargeCharge+=Scaling*COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB/r;
              break;
          }
        }

        for(i=0;i<NumberOfIntraChargeBondDipole;i++)
        {
          index=IntraChargeBondDipole[i];

          A=Components[CurrentComponent].IntraChargeBondDipole[index].A;
          typeA=Components[CurrentComponent].Type[A];
          ChargeA=Components[CurrentComponent].Charge[A];
          posA=TrialPositions[Itrial][A];

          B=Components[CurrentComponent].IntraChargeBondDipole[index].B;
          DipoleMagnitudeB=Components[CurrentComponent].BondDipoleMagnitude[B];
          B1=Components[CurrentComponent].BondDipoles[B].A;
          B2=Components[CurrentComponent].BondDipoles[B].B;
          posB1=TrialPositions[Itrial][B1];
          posB2=TrialPositions[Itrial][B2];

          dipoleB.x=posB2.x-posB1.x;
          dipoleB.y=posB2.y-posB1.y;
          dipoleB.z=posB2.z-posB1.z;
          posB.x=posB1.x+0.5*dipoleB.x;
          posB.y=posB1.y+0.5*dipoleB.y;
          posB.z=posB1.z+0.5*dipoleB.z;
          temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
          dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

          dr.x=posB.x-posA.x;
          dr.y=posB.y-posA.y;
          dr.z=posB.z-posA.z;
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          r=sqrt(rr);

          Bt1=1.0/(r*rr);
          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          EnergiesIntraChargeBondDipole-=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
        }

        for(i=0;i<NumberOfIntraBondDipoleBondDipole;i++)
        {
          index=IntraBondDipoleBondDipole[i];

          A=Components[CurrentComponent].IntraBondDipoleBondDipole[index].A;
          DipoleMagnitudeA=Components[CurrentComponent].BondDipoleMagnitude[A];
          A1=Components[CurrentComponent].BondDipoles[A].A;
          A2=Components[CurrentComponent].BondDipoles[A].B;
          posA1=TrialPositions[Itrial][A1];
          posA2=TrialPositions[Itrial][A2];

          dipoleA.x=posA2.x-posA1.x;
          dipoleA.y=posA2.y-posA1.y;
          dipoleA.z=posA2.z-posA1.z;
          posA.x=posA1.x+0.5*dipoleA.x;
          posA.y=posA1.y+0.5*dipoleA.y;
          posA.z=posA1.z+0.5*dipoleA.z;
          temp=DipoleMagnitudeA/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
          dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

          B=Components[CurrentComponent].IntraBondDipoleBondDipole[index].B;
          DipoleMagnitudeB=Components[CurrentComponent].BondDipoleMagnitude[B];
          B1=Components[CurrentComponent].BondDipoles[B].A;
          B2=Components[CurrentComponent].BondDipoles[B].B;
          posB1=TrialPositions[Itrial][B1];
          posB2=TrialPositions[Itrial][B2];

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
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          r=sqrt(rr);

          Bt1=1.0/(r*rr);
          Bt2=3.0/(r*rr*rr);

          cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
          cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
          cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
          EnergiesIntraBondDipoleBondDipole+=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
        }

      }

      UHostVDWTrial[Itrial]=EnergiesHostVDW;
      UAdsorbateVDWTrial[Itrial]=EnergiesAdsorbateVDW;
      UCationVDWTrial[Itrial]=EnergiesCationVDW;

      UHostChargeChargeTrial[Itrial]=EnergiesHostChargeCharge;
      UAdsorbateChargeChargeTrial[Itrial]=EnergiesAdsorbateChargeCharge;
      UCationChargeChargeTrial[Itrial]=EnergiesCationChargeCharge;

      UHostChargeBondDipoleTrial[Itrial]=EnergiesHostChargeBondDipole;
      UAdsorbateChargeBondDipoleTrial[Itrial]=EnergiesAdsorbateChargeBondDipole;
      UCationChargeBondDipoleTrial[Itrial]=EnergiesCationChargeBondDipole;

      UHostBondDipoleBondDipoleTrial[Itrial]=EnergiesHostBondDipoleBondDipole;
      UAdsorbateBondDipoleBondDipoleTrial[Itrial]=EnergiesAdsorbateBondDipoleBondDipole;
      UCationBondDipoleBondDipoleTrial[Itrial]=EnergiesCationBondDipoleBondDipole;

      UIntraVDWTrial[Itrial]=EnergiesIntra;
      UIntraChargeChargeTrial[Itrial]=EnergiesIntraChargeCharge;
      UIntraChargeBondDipoleTrial[Itrial]=EnergiesIntraChargeBondDipole;
      UIntraBondDipoleBondDipoleTrial[Itrial]=EnergiesIntraBondDipoleBondDipole;

      if(BiasingMethod==LJ_BIASING)
      {
         // the biasing is done using only the VDW energies
         Bfac[Itrial]=-Beta[CurrentSystem]*(EnergiesHostVDW+EnergiesAdsorbateVDW+EnergiesCationVDW+EnergiesIntra);
      }
      else
      {
         // the biasing is done using all external energies
         Bfac[Itrial]=-Beta[CurrentSystem]*
                         (EnergiesHostVDW+EnergiesCationVDW+EnergiesAdsorbateVDW+
                          EnergiesHostChargeCharge+EnergiesHostChargeBondDipole+EnergiesHostBondDipoleBondDipole+
                          EnergiesAdsorbateChargeCharge+EnergiesAdsorbateChargeBondDipole+EnergiesAdsorbateBondDipoleBondDipole+
                          EnergiesCationChargeCharge+EnergiesCationChargeBondDipole+EnergiesCationBondDipoleBondDipole+
                          EnergiesIntra+EnergiesIntraChargeCharge+EnergiesIntraChargeBondDipole+EnergiesIntraBondDipoleBondDipole);
      }

      Trover[Itrial]=FALSE;
      OVERLAP=FALSE;
    }
    else
    {
      Bfac[Itrial]=0.0;
      Trover[Itrial]=TRUE;
    } 
  } 
  return 0;
}

int GenerateTrialOrientationsSimpleSphere(int Old)
{
  int j,iu;
  int A,B;
  VECTOR dr,vec;
  int atom_nr;

  for(iu=0;iu<NumberOfTrialPositions;iu++)
  {
    UBondTrial[iu]=0.0;
    UUreyBradleyTrial[iu]=0.0;
    UBendTrial[iu]=0.0;
    UBendBendTrial[iu]=0.0;
    UInversionBendTrial[iu]=0.0;
    UTorsionTrial[iu]=0.0;
    UImproperTorsionTrial[iu]=0.0;
    UBondBondTrial[iu]=0.0;
    UBondBendTrial[iu]=0.0;
    UBondTorsionTrial[iu]=0.0;
    UBendTorsionTrial[iu]=0.0;
  }

  if(Old)
  {
    // retracing: use the old bond-lengths
    for(j=0;j<NumberOfBonds;j++)
    {
      A=Components[CurrentComponent].Bonds[Bonds[j]].A;
      B=Components[CurrentComponent].Bonds[Bonds[j]].B;
      dr.x=TrialPositions[0][A].x-TrialPositions[0][B].x;
      dr.y=TrialPositions[0][A].y-TrialPositions[0][B].y;
      dr.z=TrialPositions[0][A].z-TrialPositions[0][B].z;
      bond_length[j]=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    }

    if(Components[CurrentComponent].Groups[CurrentGroup].Rigid)
    {
      for(j=0;j<Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms;j++)
      {
        atom_nr=Components[CurrentComponent].Groups[CurrentGroup].Atoms[j];
        if(Components[CurrentComponent].ExtraFrameworkMolecule)
          TrialPositions[0][atom_nr]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[atom_nr].Position;
        else
          TrialPositions[0][atom_nr]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[atom_nr].Position;
      }
    }
  }
  else
  {
    // growing: generate new bond-lengths
    for(j=0;j<NumberOfBonds;j++)
      bond_length[j]=GenerateBondlength(Bonds[j]);
  }

  for(iu=0;iu<NumberOfTrialPositions;iu++)
  {
    if(!(Old&&iu==0))    // generate 'k' new trial positions for the second bead except for iu=0 when retracing
    {
      vec=RandomNumberOnUnitSphere();

      if(Components[CurrentComponent].Groups[CurrentGroup].Rigid)
      {
        // case 1: a rigid part of the molecule has to be grown after the first bead.
        //         As a bias random orientations of the rigid part are generated
        for(j=0;j<Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms;j++)
        {
          atom_nr=Components[CurrentComponent].Groups[CurrentGroup].Atoms[j];
          cord[j].x=Components[CurrentComponent].Positions[atom_nr].x-Components[CurrentComponent].Positions[CurrentBead].x;
          cord[j].y=Components[CurrentComponent].Positions[atom_nr].y-Components[CurrentComponent].Positions[CurrentBead].y;
          cord[j].z=Components[CurrentComponent].Positions[atom_nr].z-Components[CurrentComponent].Positions[CurrentBead].z;
        }

        RotationAroundXYZAxis(vec,cord,Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms,(2.0*RandomNumber()-1.0)*M_PI);

        for(j=0;j<Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms;j++)
        {
          atom_nr=Components[CurrentComponent].Groups[CurrentGroup].Atoms[j];
          TrialPositions[iu][atom_nr].x=TrialPositions[0][CurrentBead].x+cord[j].x;
          TrialPositions[iu][atom_nr].y=TrialPositions[0][CurrentBead].y+cord[j].y;
          TrialPositions[iu][atom_nr].z=TrialPositions[0][CurrentBead].z+cord[j].z;
        } 
      }
      else
      {
        // case 2: a flexible part of the molecule has to be grown after the first bead.
        //         As a bias random orientations of the second bead are generated on a sphere
        TrialPositions[iu][BranchAtoms[0][0]].x=TrialPositions[0][CurrentBead].x+vec.x*bond_length[0];
        TrialPositions[iu][BranchAtoms[0][0]].y=TrialPositions[0][CurrentBead].y+vec.y*bond_length[0];
        TrialPositions[iu][BranchAtoms[0][0]].z=TrialPositions[0][CurrentBead].z+vec.z*bond_length[0];
      }
    }
  }

  // compute the bond-stretching energies for non-rigid units,
  // (for rigid-units this energy is zero)
  if(!Components[CurrentComponent].Groups[CurrentGroup].Rigid)
  {
    for(iu=0;iu<NumberOfTrialPositions;iu++)
      UBondTrial[iu]=CalculateBondEnergy(Bonds[0],iu);
  }
  return 0;
}

int GenerateTrialOrientationsMCScheme(int Old)
{
  int i,j,k,jj,iu,choise,jjj,trialmove,jjopen,trialchi,A,B,C,D;
  VECTOR vec,prev_vec,dr;
  VECTOR va,vb,vc;
  REAL BondLengthNew,BondLengthOld;
  REAL newbond,newbend,newbendbend,newureybradley,energy,normalize,tmp,
       newtrs,SinusNew,SinusOld,newinversionbend,newimpropertorsion;
  VECTOR com,com2;
  REAL_MATRIX3x3 rot;
  int sel,atom_nr;
  REAL bondlength;

  // intialize the trial-energies
  for(iu=0;iu<NumberOfTrialPositions;iu++)
  {
    UBondTrial[iu]=0.0;
    UUreyBradleyTrial[iu]=0.0;
    UBendTrial[iu]=0.0;
    UBendBendTrial[iu]=0.0;
    UInversionBendTrial[iu]=0.0;
    UTorsionTrial[iu]=0.0;
    UImproperTorsionTrial[iu]=0.0;
    UBondBondTrial[iu]=0.0;
    UBondBendTrial[iu]=0.0;
    UBondTorsionTrial[iu]=0.0;
    UBendTorsionTrial[iu]=0.0;
  }

  if(Old)
  {
    // use old bond length
    for(j=0;j<NumberOfBonds;j++)
    {
      A=Components[CurrentComponent].Bonds[Bonds[j]].A;
      B=Components[CurrentComponent].Bonds[Bonds[j]].B;
      dr.x=TrialPositions[0][A].x-TrialPositions[0][B].x;
      dr.y=TrialPositions[0][A].y-TrialPositions[0][B].y;
      dr.z=TrialPositions[0][A].z-TrialPositions[0][B].z;
      bond_length[j]=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    }
  }
  else
  {
    // generate new bondlengths
    for(j=0;j<NumberOfBonds;j++)
      bond_length[j]=GenerateBondlength(Bonds[j]);
  }

  // calculate the normalized vector prev_vec
  com.x=0.0;
  com.y=0.0;
  com.z=0.0;
  if(NumberOfPreviousBeads>1)
  {
    for(j=0;j<Components[CurrentComponent].Groups[PreviousGroup].NumberOfGroupAtoms;j++)
    {
      atom_nr=Components[CurrentComponent].Groups[PreviousGroup].Atoms[j];
      com.x+=TrialPositions[0][atom_nr].x;
      com.y+=TrialPositions[0][atom_nr].y;
      com.z+=TrialPositions[0][atom_nr].z;
    }
    com.x/=Components[CurrentComponent].Groups[PreviousGroup].NumberOfGroupAtoms;
    com.y/=Components[CurrentComponent].Groups[PreviousGroup].NumberOfGroupAtoms;
    com.z/=Components[CurrentComponent].Groups[PreviousGroup].NumberOfGroupAtoms;
    prev_vec.x=(com.x-TrialPositions[0][CurrentBead].x);
    prev_vec.y=(com.y-TrialPositions[0][CurrentBead].y);
    prev_vec.z=(com.z-TrialPositions[0][CurrentBead].z);
  }
  else
  {
    prev_vec.x=TrialPositions[0][PreviousBead].x-TrialPositions[0][CurrentBead].x;
    prev_vec.y=TrialPositions[0][PreviousBead].y-TrialPositions[0][CurrentBead].y;
    prev_vec.z=TrialPositions[0][PreviousBead].z-TrialPositions[0][CurrentBead].z;
  }
  normalize=sqrt(SQR(prev_vec.x)+SQR(prev_vec.y)+SQR(prev_vec.z));
  prev_vec.x/=normalize;
  prev_vec.y/=normalize;
  prev_vec.z/=normalize;

  // calculate an initial position for the first bead
  if(!Old)
  {
    // if no already grown configuration exists create one from scratch
    if(!Components[CurrentComponent].LMCMOL)
    {

      // for the branches place the first atoms of the branches at the proper bond-lengths and angles
      for(jj=0;jj<NumberOfBranches;jj++)
      {
        vec=RandomNumberOnCone(prev_vec,120*DEG2RAD);

        // get the reference bond-length for this bond
        bondlength=GetReferenceBondLength(CurrentBead,BranchAtoms[jj][0]);

        TrialPositions[0][BranchAtoms[jj][0]].x=TrialPositions[0][CurrentBead].x+bondlength*vec.x;
        TrialPositions[0][BranchAtoms[jj][0]].y=TrialPositions[0][CurrentBead].y+bondlength*vec.y;
        TrialPositions[0][BranchAtoms[jj][0]].z=TrialPositions[0][CurrentBead].z+bondlength*vec.z;
      }

      // if one of the new branches is rigid, while the previous group is flexible, then place all of the atoms of the current group
      if(Components[CurrentComponent].Groups[CurrentGroup].Rigid&&(!Components[CurrentComponent].Groups[PreviousGroup].Rigid))
      {
        va.x=Components[CurrentComponent].Positions[CurrentBead].x;
        va.y=Components[CurrentComponent].Positions[CurrentBead].y;
        va.z=Components[CurrentComponent].Positions[CurrentBead].z;

        com.x=prev_vec.x+RandomNumber()*0.1;
        com.y=prev_vec.y+RandomNumber()*0.1;
        com.z=prev_vec.z+RandomNumber()*0.1;
        CalculateRotationMatrix(va,com,&rot);

        vc.x=rot.ax*va.x+rot.bx*va.y+rot.cx*va.z;
        vc.y=rot.ay*va.x+rot.by*va.y+rot.cy*va.z;
        vc.z=rot.az*va.x+rot.bz*va.y+rot.cz*va.z;

        for(j=0;j<Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms;j++)
        {
          atom_nr=Components[CurrentComponent].Groups[CurrentGroup].Atoms[j];
          va.x=Components[CurrentComponent].Positions[atom_nr].x;
          va.y=Components[CurrentComponent].Positions[atom_nr].y;
          va.z=Components[CurrentComponent].Positions[atom_nr].z;

          vb.x=(rot.ax*va.x+rot.bx*va.y+rot.cx*va.z)-vc.x;
          vb.y=(rot.ay*va.x+rot.by*va.y+rot.cy*va.z)-vc.y;
          vb.z=(rot.az*va.x+rot.bz*va.y+rot.cz*va.z)-vc.z;

          TrialPositions[0][atom_nr].x=TrialPositions[0][CurrentBead].x+vb.x;
          TrialPositions[0][atom_nr].y=TrialPositions[0][CurrentBead].y+vb.y;
          TrialPositions[0][atom_nr].z=TrialPositions[0][CurrentBead].z+vb.z;
        }
      }
    }
    else     // calculate initial positions from a previously stored configuration
    {
      // calculate rotation matrix for "matching" the positions using rotation around the prev_vec
      if(NumberOfPreviousBeads>1)
      {
        com.x=0.0;
        com.y=0.0;
        com.z=0.0;
        for(j=0;j<Components[CurrentComponent].Groups[PreviousGroup].NumberOfGroupAtoms;j++)
        {
          atom_nr=Components[CurrentComponent].Groups[PreviousGroup].Atoms[j];
          com.x+=Components[CurrentComponent].RMCMOL[atom_nr].x;
          com.y+=Components[CurrentComponent].RMCMOL[atom_nr].y;
          com.z+=Components[CurrentComponent].RMCMOL[atom_nr].z;
        }
        com.x/=Components[CurrentComponent].Groups[PreviousGroup].NumberOfGroupAtoms;
        com.y/=Components[CurrentComponent].Groups[PreviousGroup].NumberOfGroupAtoms;
        com.z/=Components[CurrentComponent].Groups[PreviousGroup].NumberOfGroupAtoms;
        va.x=com.x-Components[CurrentComponent].RMCMOL[CurrentBead].x;
        va.y=com.y-Components[CurrentComponent].RMCMOL[CurrentBead].y;
        va.z=com.z-Components[CurrentComponent].RMCMOL[CurrentBead].z;
      }
      else
      {
        va.x=Components[CurrentComponent].RMCMOL[CurrentBead].x-
             Components[CurrentComponent].RMCMOL[PreviousBead].x;
        va.y=Components[CurrentComponent].RMCMOL[CurrentBead].y-
             Components[CurrentComponent].RMCMOL[PreviousBead].y;
        va.z=Components[CurrentComponent].RMCMOL[CurrentBead].z-
             Components[CurrentComponent].RMCMOL[PreviousBead].z;
      }

      com2.x=-prev_vec.x;
      com2.y=-prev_vec.y;
      com2.z=-prev_vec.z;

      CalculateRotationMatrix(va,com2,&rot);

      vc.x=rot.ax*va.x+rot.bx*va.y+rot.cx*va.z;
      vc.y=rot.ay*va.x+rot.by*va.y+rot.cy*va.z;
      vc.z=rot.az*va.x+rot.bz*va.y+rot.cz*va.z;

      for(i=0;i<NumberOfBeadsToBePlaced;i++)
      {
        if(NumberOfPreviousBeads>1)
        {
          va.x=com.x-Components[CurrentComponent].RMCMOL[BeadsToBePlaced[i]].x;
          va.y=com.y-Components[CurrentComponent].RMCMOL[BeadsToBePlaced[i]].y;
          va.z=com.z-Components[CurrentComponent].RMCMOL[BeadsToBePlaced[i]].z;
        }
        else
        {
          va.x=Components[CurrentComponent].RMCMOL[BeadsToBePlaced[i]].x-
               Components[CurrentComponent].RMCMOL[PreviousBead].x;
          va.y=Components[CurrentComponent].RMCMOL[BeadsToBePlaced[i]].y-
               Components[CurrentComponent].RMCMOL[PreviousBead].y;
          va.z=Components[CurrentComponent].RMCMOL[BeadsToBePlaced[i]].z-
               Components[CurrentComponent].RMCMOL[PreviousBead].z;
        }

        vb.x=(rot.ax*va.x+rot.bx*va.y+rot.cx*va.z)-vc.x;
        vb.y=(rot.ay*va.x+rot.by*va.y+rot.cy*va.z)-vc.y;
        vb.z=(rot.az*va.x+rot.bz*va.y+rot.cz*va.z)-vc.z;

        TrialPositions[0][BeadsToBePlaced[i]].x=TrialPositions[0][CurrentBead].x+vb.x;
        TrialPositions[0][BeadsToBePlaced[i]].y=TrialPositions[0][CurrentBead].y+vb.y;
        TrialPositions[0][BeadsToBePlaced[i]].z=TrialPositions[0][CurrentBead].z+vb.z;
      }
    }
  }

  for(j=0;j<NumberOfBonds;j++)
    eobond[j]=enbond[j]=CalculateBondEnergy(Bonds[j],0);

  for(j=0;j<NumberOfBends;j++)
    eobend[j]=enbend[j]=CalculateBendEnergy(Bends[j],0);

  for(j=0;j<NumberOfBendBends;j++)
    eobendbend[j]=enbendbend[j]=CalculateBendBendEnergy(BendBends[j],0);

  for(j=0;j<NumberOfUreyBradleys;j++)
    eoureybradley[j]=enureybradley[j]=CalculateUreyBradleyEnergy(UreyBradleys[j],0);

  for(j=0;j<NumberOfInversionBends;j++)
    eoinversionbend[j]=eninversionbend[j]=CalculateInversionBendEnergy(InversionBends[j],0);

  for(j=0;j<NumberOfImproperTorsions;j++)
    eoimpropertorsion[j]=enimpropertorsion[j]=CalculateImproperTorsionEnergy(ImproperTorsions[j],0);

  if(Components[CurrentComponent].LMCMOL)
    trialmove=2*NumberOfTrialMovesPerOpenBead*NumberOfBranches;
  else trialmove=4*NumberOfTrialMovesPerOpenBead*NumberOfBranches;

  trialchi=(int)(trialmove/4.0);
  // after trialchi steps, an attempt is made to correct the
  // stereochemistry of the molecule...

  if(!Old)
  {
    // When there is more than one bead open and no chirality: decide with a 
    // random number to change the identity of two beads
    // This is needed because we usually start from an previous stored configuration of this 
    // component type leading otherwise to an unaltered arrangement of the branches.

    if((!Components[CurrentComponent].Chirality[CurrentBead])&&(NumberOfBonds>=2)&&(RandomNumber()<0.5))
    {
      A=BranchAtoms[0][0];
      C=BranchAtoms[1][0];

      cord[A].x=TrialPositions[0][A].x-TrialPositions[0][CurrentBead].x;
      cord[A].y=TrialPositions[0][A].y-TrialPositions[0][CurrentBead].y;
      cord[A].z=TrialPositions[0][A].z-TrialPositions[0][CurrentBead].z;

      cord[C].x=TrialPositions[0][C].x-TrialPositions[0][CurrentBead].x;
      cord[C].y=TrialPositions[0][C].y-TrialPositions[0][CurrentBead].y;
      cord[C].z=TrialPositions[0][C].z-TrialPositions[0][CurrentBead].z;

      tmp=sqrt(SQR(cord[A].x)+SQR(cord[A].y)+SQR(cord[A].z))/sqrt(SQR(cord[C].x)+SQR(cord[C].y)+SQR(cord[C].z));
      TrialPositions[0][A].x=TrialPositions[0][CurrentBead].x+tmp*cord[C].x;
      TrialPositions[0][A].y=TrialPositions[0][CurrentBead].y+tmp*cord[C].y;
      TrialPositions[0][A].z=TrialPositions[0][CurrentBead].z+tmp*cord[C].z;

      tmp=1.0/tmp;
      TrialPositions[0][C].x=TrialPositions[0][CurrentBead].x+tmp*cord[A].x;
      TrialPositions[0][C].y=TrialPositions[0][CurrentBead].y+tmp*cord[A].y;
      TrialPositions[0][C].z=TrialPositions[0][CurrentBead].z+tmp*cord[A].z;

      for(j=0;j<NumberOfBonds;j++)
        eobond[j]=enbond[j]=CalculateBondEnergy(Bonds[j],0);

      for(j=0;j<NumberOfBends;j++)
        eobend[j]=enbend[j]=CalculateBendEnergy(Bends[j],0);

      for(j=0;j<NumberOfBendBends;j++)
        eobendbend[j]=enbendbend[j]=CalculateBendBendEnergy(BendBends[j],0);

      for(j=0;j<NumberOfUreyBradleys;j++)
        eoureybradley[j]=enureybradley[j]=CalculateUreyBradleyEnergy(UreyBradleys[j],0);

      for(j=0;j<NumberOfInversionBends;j++)
         eoinversionbend[j]=eninversionbend[j]=CalculateInversionBendEnergy(InversionBends[j],0);

      for(j=0;j<NumberOfImproperTorsions;j++)
         eoimpropertorsion[j]=enimpropertorsion[j]=CalculateImproperTorsionEnergy(ImproperTorsions[j],0);
    }

    for(jjj=0;jjj<trialmove;jjj++)
    {
      // select a trial move at random
      //
      // trialmove 1: change cone angle
      // trialmove 2: change position on the cone
      // trialmove 3: rotate everything around the cone so
      //              that there are no barriers
      //
      // probabilities of selecting trial moves:
      //
      //     only one bond-bending : tr1=0.9, tr2=0.0, tr3=0.1
      //     else                  : tr1=0.4, tr2=0.4, tr3=0.2
      //
      //     when IU=1, then the first 20 trialmoves TR2=1.0
      //     this is to equilibrate the positions on the cone very fast,
      //     because the cone-angles are chosen as in the equilibrium value

      if(NumberOfBends==1)
      {
        tmp=RandomNumber();
        if(tmp<0.3)
          choise=1;
        else
          choise=2;
      }
      else 
      {
        if(jjj<20)
          choise=3;
        else
        {
          tmp=RandomNumber();
          if(tmp<0.2) choise=1;
          else if(tmp<0.6) choise=2;
          else choise=3;
        }
      }

      for(j=0;j<NumberOfBonds;j++)
        enbond[j]=eobond[j];

      for(j=0;j<NumberOfBends;j++)
        enbend[j]=eobend[j];

      for(j=0;j<NumberOfBendBends;j++)
        enbendbend[j]=eobendbend[j];

      for(j=0;j<NumberOfUreyBradleys;j++)
        enureybradley[j]=eoureybradley[j];

      for(j=0;j<NumberOfInversionBends;j++)
        eninversionbend[j]=eoinversionbend[j];

      for(j=0;j<NumberOfImproperTorsions;j++)
        enimpropertorsion[j]=eoimpropertorsion[j];

      newbond=0.0;
      newbend=0.0;
      newbendbend=0.0;
      newureybradley=0.0;
      newinversionbend=0.0;
      newimpropertorsion=0.0;
      newtrs=0.0;

      switch(choise)
      {
        case 1:    // change bondlength of real bonds
          if(NumberOfBonds>0)
          {
            jj=(int)(RandomNumber()*(REAL)NumberOfBonds);

            if(Components[CurrentComponent].BondType[Bonds[jj]]!=FIXED_BOND)
            {
              jjopen=BranchAtoms[jj][0];
 
              store[0]=TrialPositions[0][jjopen];

              Components[CurrentComponent].ACNBOND[CurrentSystem][Bonds[jj]]+=1.0;

              vec.x=TrialPositions[0][jjopen].x-TrialPositions[0][CurrentBead].x;
              vec.y=TrialPositions[0][jjopen].y-TrialPositions[0][CurrentBead].y;
              vec.z=TrialPositions[0][jjopen].z-TrialPositions[0][CurrentBead].z;
              BondLengthOld=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));

              BondLengthNew=BondLengthOld+(2.0*RandomNumber()-1.0)*
                   Components[CurrentComponent].RMAX0[CurrentSystem][jjopen];
              TrialPositions[0][jjopen].x=TrialPositions[0][CurrentBead].x+(BondLengthNew/BondLengthOld)*vec.x;
              TrialPositions[0][jjopen].y=TrialPositions[0][CurrentBead].y+(BondLengthNew/BondLengthOld)*vec.y;
              TrialPositions[0][jjopen].z=TrialPositions[0][CurrentBead].z+(BondLengthNew/BondLengthOld)*vec.z;
          
              enbond[jj]=CalculateBondEnergy(Bonds[jj],0);
              newbond+=enbond[jj]-eobond[jj];

              for(j=0;j<NumberOfUreyBradleys;j++)
              {
                enureybradley[j]=CalculateUreyBradleyEnergy(UreyBradleys[j],0);
                newureybradley+=enureybradley[j]-eoureybradley[j];
              }

              // this move can change the inversion bend energy for 'PLANAR' inversion
              for(j=0;j<NumberOfInversionBends;j++)
              {
                eninversionbend[j]=CalculateInversionBendEnergy(InversionBends[j],0);
                newinversionbend+=eninversionbend[j]-eoinversionbend[j];
              }
              for(j=0;j<NumberOfImproperTorsions;j++)
              {
                enimpropertorsion[j]=CalculateImproperTorsionEnergy(ImproperTorsions[j],0);
                newimpropertorsion+=enimpropertorsion[j]-eoimpropertorsion[j];
              }


              if(RandomNumber()<SQR(BondLengthNew/BondLengthOld)*exp(-Beta[CurrentSystem]*(newbond+newureybradley+newinversionbend+newimpropertorsion)))
              {
                Components[CurrentComponent].ACSBOND[CurrentSystem][Bonds[jj]]+=1.0;
                eobond[jj]=enbond[jj];
                eoureybradley[jj]=enureybradley[jj];
                for(j=0;j<NumberOfInversionBends;j++)
                  eoinversionbend[j]=eninversionbend[j];
                for(j=0;j<NumberOfImproperTorsions;j++)
                  eoimpropertorsion[j]=enimpropertorsion[j];
              }
              else
                TrialPositions[0][jjopen]=store[0];
            }
          }
          break;
        case 2:    // change bend angle of one chosen bead
          if(NumberOfBranches>0)
          {
            jj=(int)(RandomNumber()*(REAL)NumberOfBranches);

            for(i=0;i<NumberOfBranchAtoms[jj];i++)
            {
              jjopen=BranchAtoms[jj][i];
              store[i]=TrialPositions[0][jjopen];
              cord[i].x=TrialPositions[0][jjopen].x-TrialPositions[0][CurrentBead].x;
              cord[i].y=TrialPositions[0][jjopen].y-TrialPositions[0][CurrentBead].y;
              cord[i].z=TrialPositions[0][jjopen].z-TrialPositions[0][CurrentBead].z;
            }

            if(NumberOfBranchAtoms[jj]==1)
            {
              SinusOld=sin(CalculateAngle(PreviousBead,CurrentBead,BranchAtoms[jj][0],0));
              dr.x=TrialPositions[0][CurrentBead].x-TrialPositions[0][BranchAtoms[jj][0]].x;
              dr.y=TrialPositions[0][CurrentBead].y-TrialPositions[0][BranchAtoms[jj][0]].y;
              dr.z=TrialPositions[0][CurrentBead].z-TrialPositions[0][BranchAtoms[jj][0]].z;
              perpen[jj]=Perpendicular(dr,prev_vec);
            }
            else
            {
              com.x=com.y=com.z=0.0;
              for(i=0;i<Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms;i++)
              {
                atom_nr=Components[CurrentComponent].Groups[CurrentGroup].Atoms[i];
                com.x+=TrialPositions[0][atom_nr].x;
                com.y+=TrialPositions[0][atom_nr].y;
                com.z+=TrialPositions[0][atom_nr].z;
              }
              com.x/=Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms;
              com.y/=Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms;
              com.z/=Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms;

              dr.x=com.x-TrialPositions[0][CurrentBead].x;
              dr.y=com.y-TrialPositions[0][CurrentBead].y;
              dr.z=com.z-TrialPositions[0][CurrentBead].z;
              perpen[jj]=Perpendicular(dr,prev_vec);

              SinusOld=sin(CalculateAngle2(PreviousBead,CurrentBead,CurrentGroup,0));
            }

            Components[CurrentComponent].ACNBEND[CurrentSystem][CurrentBead]+=1.0;

            RotationAroundXYZAxis(perpen[jj],cord,NumberOfBranchAtoms[jj],(2.0*RandomNumber()-1.0)*
                 Components[CurrentComponent].RMAX1[CurrentSystem][CurrentBead]);

            for(i=0;i<NumberOfBranchAtoms[jj];i++)
            {
              jjopen=BranchAtoms[jj][i];
              TrialPositions[0][jjopen].x=TrialPositions[0][CurrentBead].x+cord[i].x;
              TrialPositions[0][jjopen].y=TrialPositions[0][CurrentBead].y+cord[i].y;
              TrialPositions[0][jjopen].z=TrialPositions[0][CurrentBead].z+cord[i].z;
            }

            if(NumberOfBranchAtoms[jj]==1)
              SinusNew=sin(CalculateAngle(PreviousBead,CurrentBead,BranchAtoms[jj][0],0));
            else
            {
              SinusNew=sin(CalculateAngle2(PreviousBead,CurrentBead,CurrentGroup,0));
            }

            for(j=0;j<NumberOfBends;j++)
            {
              enbend[j]=CalculateBendEnergy(Bends[j],0);
              newbend+=enbend[j]-eobend[j];
            }

            for(j=0;j<NumberOfBendBends;j++)
            {
              enbendbend[j]=CalculateBendBendEnergy(BendBends[j],0);
              newbendbend+=enbendbend[j]-eobendbend[j];
            }

            for(j=0;j<NumberOfUreyBradleys;j++)
            {
              enureybradley[j]=CalculateUreyBradleyEnergy(UreyBradleys[j],0);
              newureybradley+=enureybradley[j]-eoureybradley[j];
            }

            for(j=0;j<NumberOfInversionBends;j++)
            {
              eninversionbend[j]=CalculateInversionBendEnergy(InversionBends[j],0);
              newinversionbend+=eninversionbend[j]-eoinversionbend[j];
            }

            for(j=0;j<NumberOfImproperTorsions;j++)
            {
              enimpropertorsion[j]=CalculateImproperTorsionEnergy(ImproperTorsions[j],0);
              newimpropertorsion+=enimpropertorsion[j]-eoimpropertorsion[j];
            }

            if(RandomNumber()<((SinusNew/SinusOld)*exp(-Beta[CurrentSystem]*(newbend+newureybradley+newinversionbend+newimpropertorsion))))
            {
              Components[CurrentComponent].ACSBEND[CurrentSystem][CurrentBead]+=1.0;

              for(j=0;j<NumberOfBonds;j++)
                eobond[jj]=enbond[jj];
              for(j=0;j<NumberOfBends;j++)
                eobend[j]=enbend[j];
              for(j=0;j<NumberOfBendBends;j++)
                eobendbend[j]=enbendbend[j];
              for(j=0;j<NumberOfUreyBradleys;j++)
                eoureybradley[j]=enureybradley[j];
              for(j=0;j<NumberOfInversionBends;j++)
                eoinversionbend[j]=eninversionbend[j];
              for(j=0;j<NumberOfImproperTorsions;j++)
                eoimpropertorsion[j]=enimpropertorsion[j];
            }
            else
            {
              for(i=0;i<NumberOfBranchAtoms[jj];i++)
              {
                jjopen=BranchAtoms[jj][i];
                TrialPositions[0][jjopen]=store[i];
              }
            }
          }
          break;
        case 3:        // rotate one position on the cone while keeping the cone angle fixed
          if(NumberOfBranches>0)
          {
            jj=(int)(RandomNumber()*(REAL)NumberOfBranches);

            for(i=0;i<NumberOfBranchAtoms[jj];i++)
            {
              jjopen=BranchAtoms[jj][i];
              store[i]=TrialPositions[0][jjopen];
              cord[i].x=TrialPositions[0][jjopen].x-TrialPositions[0][CurrentBead].x;
              cord[i].y=TrialPositions[0][jjopen].y-TrialPositions[0][CurrentBead].y;
              cord[i].z=TrialPositions[0][jjopen].z-TrialPositions[0][CurrentBead].z;
            }

            Components[CurrentComponent].ACNTRS[CurrentSystem][CurrentBead]+=1.0;

            RotationAroundXYZAxis(prev_vec,cord,NumberOfBranchAtoms[jj],(2.0*RandomNumber()-1.0)*
                 Components[CurrentComponent].RMAX2[CurrentSystem][CurrentBead]);

            for(i=0;i<NumberOfBranchAtoms[jj];i++)
            {
              jjopen=BranchAtoms[jj][i];
              TrialPositions[0][jjopen].x=TrialPositions[0][CurrentBead].x+cord[i].x;
              TrialPositions[0][jjopen].y=TrialPositions[0][CurrentBead].y+cord[i].y;
              TrialPositions[0][jjopen].z=TrialPositions[0][CurrentBead].z+cord[i].z;
            }

            for(j=0;j<NumberOfBends;j++)
            {
              enbend[j]=CalculateBendEnergy(Bends[j],0);
              newbend+=enbend[j]-eobend[j];
            }

            for(j=0;j<NumberOfBendBends;j++)
            {
              enbendbend[j]=CalculateBendBendEnergy(BendBends[j],0);
              newbendbend+=enbendbend[j]-eobendbend[j];
            }

            for(j=0;j<NumberOfUreyBradleys;j++)
            {
              enureybradley[j]=CalculateUreyBradleyEnergy(UreyBradleys[j],0);
              newureybradley+=enureybradley[j]-eoureybradley[j];
            }

            for(j=0;j<NumberOfInversionBends;j++)
            {
              eninversionbend[j]=CalculateInversionBendEnergy(InversionBends[j],0);
              newinversionbend+=eninversionbend[j]-eoinversionbend[j];
            }

            for(j=0;j<NumberOfImproperTorsions;j++)
            {
              enimpropertorsion[j]=CalculateImproperTorsionEnergy(ImproperTorsions[j],0);
              newimpropertorsion+=enimpropertorsion[j]-eoimpropertorsion[j];
            }

            if(RandomNumber()<exp(-Beta[CurrentSystem]*(newbend+newureybradley+newinversionbend+newimpropertorsion)))
            {
              Components[CurrentComponent].ACSTRS[CurrentSystem][CurrentBead]+=1.0;

              for(j=0;j<NumberOfBends;j++)
                eobend[j]=enbend[j];
              for(j=0;j<NumberOfBendBends;j++)
                eobendbend[j]=enbendbend[j];
              for(j=0;j<NumberOfUreyBradleys;j++)
                eoureybradley[j]=enureybradley[j];
              for(j=0;j<NumberOfInversionBends;j++)
                eoinversionbend[j]=eninversionbend[j];
              for(j=0;j<NumberOfImproperTorsions;j++)
                eoimpropertorsion[j]=enimpropertorsion[j];

              for(i=0;i<NumberOfBranchAtoms[jj];i++)
              {
                jjopen=BranchAtoms[jj][i];
                dr.x=TrialPositions[0][CurrentBead].x-TrialPositions[0][jjopen].x;
                dr.y=TrialPositions[0][CurrentBead].y-TrialPositions[0][jjopen].y;
                dr.z=TrialPositions[0][CurrentBead].z-TrialPositions[0][jjopen].z;
                perpen[jjopen]=Perpendicular(dr,prev_vec);
              }
            }
            else
            {
              for(i=0;i<NumberOfBranchAtoms[jj];i++)
              {
                jjopen=BranchAtoms[jj][i];
                TrialPositions[0][jjopen]=store[i];
              }
            }
          }
          break;
      }

      // check if the stereochemistry of the molecule is correct
      if((Components[CurrentComponent].Chirality[CurrentBead])&&
         (NumberOfBonds>=2)&&((jjj==trialchi)||(jjj==trialmove-1)))
      {
        A=Components[CurrentComponent].ChiralA[CurrentBead];
        B=Components[CurrentComponent].ChiralB[CurrentBead];
        C=Components[CurrentComponent].ChiralC[CurrentBead];
        D=Components[CurrentComponent].ChiralD[CurrentBead];

        cord[A].x=TrialPositions[0][A].x-TrialPositions[0][CurrentBead].x;
        cord[A].y=TrialPositions[0][A].y-TrialPositions[0][CurrentBead].y;
        cord[A].z=TrialPositions[0][A].z-TrialPositions[0][CurrentBead].z;

        cord[C].x=TrialPositions[0][C].x-TrialPositions[0][CurrentBead].x;
        cord[C].y=TrialPositions[0][C].y-TrialPositions[0][CurrentBead].y;
        cord[C].z=TrialPositions[0][C].z-TrialPositions[0][CurrentBead].z;

        cord[D].x=TrialPositions[0][D].x-TrialPositions[0][CurrentBead].x;
        cord[D].y=TrialPositions[0][D].y-TrialPositions[0][CurrentBead].y;
        cord[D].z=TrialPositions[0][D].z-TrialPositions[0][CurrentBead].z;

        // calculate (D x C) . A    if > 0 then R else S
        tmp=cord[D].x*cord[C].y*cord[A].z+cord[D].y*cord[C].z*cord[A].x+
            cord[D].z*cord[C].x*cord[A].y-cord[D].x*cord[C].z*cord[A].y-
            cord[D].y*cord[C].x*cord[A].z-cord[D].z*cord[C].y*cord[A].x;

        if(((tmp<0.0)&&(Components[CurrentComponent].ChiralityType[CurrentBead]==R_CHIRAL))||
           ((tmp>0.0)&&(Components[CurrentComponent].ChiralityType[CurrentBead]==S_CHIRAL)))
        {
          // the stereochemistry of the molecule is incorrect
          // swap the positions of two beads to be placed

          A=BranchAtoms[0][0];
          C=BranchAtoms[1][0];

          cord[A].x=TrialPositions[0][A].x-TrialPositions[0][CurrentBead].x;
          cord[A].y=TrialPositions[0][A].y-TrialPositions[0][CurrentBead].y;
          cord[A].z=TrialPositions[0][A].z-TrialPositions[0][CurrentBead].z;

          cord[C].x=TrialPositions[0][C].x-TrialPositions[0][CurrentBead].x;
          cord[C].y=TrialPositions[0][C].y-TrialPositions[0][CurrentBead].y;
          cord[C].z=TrialPositions[0][C].z-TrialPositions[0][CurrentBead].z;

          tmp=sqrt(SQR(cord[A].x)+SQR(cord[A].y)+SQR(cord[A].z))/sqrt(SQR(cord[C].x)+SQR(cord[C].y)+SQR(cord[C].z));
          TrialPositions[0][A].x=TrialPositions[0][CurrentBead].x+tmp*cord[C].x;
          TrialPositions[0][A].y=TrialPositions[0][CurrentBead].y+tmp*cord[C].y;
          TrialPositions[0][A].z=TrialPositions[0][CurrentBead].z+tmp*cord[C].z;

          tmp=1.0/tmp;
          TrialPositions[0][C].x=TrialPositions[0][CurrentBead].x+tmp*cord[A].x;
          TrialPositions[0][C].y=TrialPositions[0][CurrentBead].y+tmp*cord[A].y;
          TrialPositions[0][C].z=TrialPositions[0][CurrentBead].z+tmp*cord[A].z;

          for(j=0;j<NumberOfBonds;j++)
            eobond[j]=enbond[j]=CalculateBondEnergy(Bonds[j],0);

          for(j=0;j<NumberOfBends;j++)
            eobend[j]=enbend[j]=CalculateBendEnergy(Bends[j],0);

          for(j=0;j<NumberOfBendBends;j++)
            eobendbend[j]=enbendbend[j]=CalculateBendBendEnergy(BendBends[j],0);

          for(j=0;j<NumberOfUreyBradleys;j++)
            eoureybradley[j]=enureybradley[j]=CalculateUreyBradleyEnergy(UreyBradleys[j],0);

          for(j=0;j<NumberOfInversionBends;j++)
            eoinversionbend[j]=eninversionbend[j]=CalculateInversionBendEnergy(InversionBends[j],0);

          for(j=0;j<NumberOfImproperTorsions;j++)
            eoimpropertorsion[j]=enimpropertorsion[j]=CalculateImproperTorsionEnergy(ImproperTorsions[j],0);

          // do an additional 50 steps
          trialmove+=50;
        }
      }
    }
  }

  for(j=0;j<NumberOfBonds;j++)
  {
    UBondTrial[0]+=(energy=CalculateBondEnergy(Bonds[j],0));
    if(fabs(energy-eobond[j])>0.0001)
    {
      printf("Error: no energy conservation in internal MC scheme (Bond %f)\n",
              (double)fabs(energy-eobond[j]));
      exit(2);
    }
  }

  for(j=0;j<NumberOfUreyBradleys;j++)
  {
    UUreyBradleyTrial[0]+=(energy=CalculateUreyBradleyEnergy(UreyBradleys[j],0));
    if(fabs(energy-eoureybradley[j])>0.0001)
    {
      printf("Error: no energy conservation in internal MC scheme (UreyBradley %f)\n",
              (double)fabs(energy-eoureybradley[j]));
      exit(2);
    }
  }

  for(j=0;j<NumberOfBends;j++)
  {
    UBendTrial[0]+=(energy=CalculateBendEnergy(Bends[j],0));
    if(fabs(energy-eobend[j])>0.0001)
    {
      printf("Error: no energy conservation in internal MC scheme (Bend %f)\n",
              (double)fabs(energy-eobend[j]));
      exit(2);
    }
  }

  for(j=0;j<NumberOfBendBends;j++)
  {
    UBendBendTrial[0]+=(energy=CalculateBendBendEnergy(BendBends[j],0));
    if(fabs(energy-eobendbend[j])>0.0001)
    {
      printf("Error: no energy conservation in internal MC scheme (Bend-Bend %f)\n",
              (double)fabs(energy-eobendbend[j]));
      exit(2);
    }
  }

  for(j=0;j<NumberOfInversionBends;j++)
  {
    UInversionBendTrial[0]+=(energy=CalculateInversionBendEnergy(InversionBends[j],0));
    if(fabs(energy-eoinversionbend[j])>0.0001)
    {
      printf("Error: no energy conservation in internal MC scheme (Inversion-bend %f)\n",
              (double)fabs(energy-eoinversionbend[j]));
      exit(2);
    }
  }

  for(j=0;j<NumberOfImproperTorsions;j++)
  {
    UImproperTorsionTrial[0]+=(energy=CalculateImproperTorsionEnergy(ImproperTorsions[j],0));
    if(fabs(energy-eoimpropertorsion[j])>0.0001)
    {
      printf("Error: no energy conservation in internal MC scheme (Improper torsion %f)\n",
              (double)fabs(energy-eoimpropertorsion[j]));
      exit(2);
    }
  }

  // use rotation around the prev_vec to generate torsion-biasing
  // this rotation will leave all bending angles etc unchanged
  if(NumberOfPreviousBeads>1)
  {
    for(iu=0;iu<NumberOfTrialPositions-1;iu++)
      for(jj=0;jj<NumberOfBeadsToBePlaced;jj++)
        TrialPositions[iu+1][BeadsToBePlaced[jj]]=TrialPositions[iu][BeadsToBePlaced[jj]];
    for(iu=0;iu<NumberOfTrialPositions;iu++)
      RosenbluthTorsion[iu]=1.0;
  }
  else
  {
    for(jj=0;jj<NumberOfBeadsToBePlaced;jj++)
      store[jj]=TrialPositions[0][BeadsToBePlaced[jj]];

    for(iu=0;iu<NumberOfTrialPositions;iu++)
    {
      for(k=0;k<NumberOfTrialPositionsTorsion;k++)
      {
        for(jj=0;jj<NumberOfBeadsToBePlaced;jj++)
        {
          TrialPositions[iu][BeadsToBePlaced[jj]]=store[jj];
          cord[jj].x=TrialPositions[iu][BeadsToBePlaced[jj]].x-TrialPositions[iu][CurrentBead].x;
          cord[jj].y=TrialPositions[iu][BeadsToBePlaced[jj]].y-TrialPositions[iu][CurrentBead].y;
          cord[jj].z=TrialPositions[iu][BeadsToBePlaced[jj]].z-TrialPositions[iu][CurrentBead].z;
        }
        if((iu==0)&&(k==0)&&Old)
          angle[k]=0.0;
        else
        {
          angle[k]=(2.0*RandomNumber()-1.0)*M_PI;
          RotationAroundXYZAxis(prev_vec,cord,NumberOfBeadsToBePlaced,angle[k]);
        }

        for(jj=0;jj<NumberOfBeadsToBePlaced;jj++)
        {
          TrialPositions[iu][BeadsToBePlaced[jj]].x=cord[jj].x+TrialPositions[iu][CurrentBead].x;
          TrialPositions[iu][BeadsToBePlaced[jj]].y=cord[jj].y+TrialPositions[iu][CurrentBead].y;
          TrialPositions[iu][BeadsToBePlaced[jj]].z=cord[jj].z+TrialPositions[iu][CurrentBead].z;
        }

        newtrs=0.0;
        for(j=0;j<NumberOfTorsions;j++)
          newtrs+=CalculateTorsionEnergy(Torsions[j],iu);

        for(j=0;j<NumberOfBondTorsions;j++)
          newtrs+=CalculateBondTorsionEnergy(BondTorsions[j],iu);

        for(j=0;j<NumberOfBendTorsions;j++)
          newtrs+=CalculateBendTorsionEnergy(BendTorsions[j],iu);

        Bfac[k]=-Beta[CurrentSystem]*newtrs;
        Trover[k]=FALSE;
      }

      RosenbluthTorsion[iu]=ComputeNormalizedRosenbluthWeight(Bfac,Trover,NumberOfTrialPositionsTorsion);

      if((iu==0)&&Old)
        sel=0;
      else
      {
        if(NumberOfTorsions>0)
          sel=SelectTrialPosition(Bfac,Trover,NumberOfTrialPositionsTorsion);
        else
          sel=(int)(RandomNumber()*(REAL)NumberOfTrialPositionsTorsion);
      }

      for(jj=0;jj<NumberOfBeadsToBePlaced;jj++)
      {
        TrialPositions[iu][BeadsToBePlaced[jj]]=store[jj];
        cord[jj].x=TrialPositions[iu][BeadsToBePlaced[jj]].x-TrialPositions[iu][CurrentBead].x;
        cord[jj].y=TrialPositions[iu][BeadsToBePlaced[jj]].y-TrialPositions[iu][CurrentBead].y;
        cord[jj].z=TrialPositions[iu][BeadsToBePlaced[jj]].z-TrialPositions[iu][CurrentBead].z;
      }

      RotationAroundXYZAxis(prev_vec,cord,NumberOfBeadsToBePlaced,angle[sel]);

      for(jj=0;jj<NumberOfBeadsToBePlaced;jj++)
      {
        TrialPositions[iu][BeadsToBePlaced[jj]].x=cord[jj].x+TrialPositions[iu][CurrentBead].x;
        TrialPositions[iu][BeadsToBePlaced[jj]].y=cord[jj].y+TrialPositions[iu][CurrentBead].y;
        TrialPositions[iu][BeadsToBePlaced[jj]].z=cord[jj].z+TrialPositions[iu][CurrentBead].z;
      }
    }
  }


  // recompute internal bond, bend, Urey-Bradley, and torsion energy
  for(iu=0;iu<NumberOfTrialPositions;iu++)
  {
    UBondTrial[iu]=0.0;
    for(j=0;j<NumberOfBonds;j++)
      UBondTrial[iu]+=CalculateBondEnergy(Bonds[j],iu);

    UUreyBradleyTrial[iu]=0.0;
    for(j=0;j<NumberOfUreyBradleys;j++)
      UUreyBradleyTrial[iu]+=CalculateUreyBradleyEnergy(UreyBradleys[j],iu);

    UBendTrial[iu]=0.0;
    for(j=0;j<NumberOfBends;j++)
      UBendTrial[iu]+=CalculateBendEnergy(Bends[j],iu);

    UBendBendTrial[iu]=0.0;
    for(j=0;j<NumberOfBendBends;j++)
      UBendBendTrial[iu]+=CalculateBendBendEnergy(BendBends[j],iu);

    UInversionBendTrial[iu]=0.0;
    for(j=0;j<NumberOfInversionBends;j++)
      UInversionBendTrial[iu]+=CalculateInversionBendEnergy(InversionBends[j],iu);

    UTorsionTrial[iu]=0.0;
    for(j=0;j<NumberOfTorsions;j++)
      UTorsionTrial[iu]+=CalculateTorsionEnergy(Torsions[j],iu);

    UImproperTorsionTrial[iu]=0.0;
    for(j=0;j<NumberOfImproperTorsions;j++)
      UImproperTorsionTrial[iu]+=CalculateImproperTorsionEnergy(ImproperTorsions[j],iu);

    UBondBondTrial[iu]=0.0;
    for(j=0;j<NumberOfBondBonds;j++)
      UBondBondTrial[iu]+=CalculateBondBondEnergy(BondBonds[j],iu);

    UBondBendTrial[iu]=0.0;
    for(j=0;j<NumberOfBondBends;j++)
      UBondBendTrial[iu]+=CalculateBondBendEnergy(BondBends[j],iu);

    UBondTorsionTrial[iu]=0.0;
    for(j=0;j<NumberOfBondTorsions;j++)
      UBondTorsionTrial[iu]+=CalculateBondTorsionEnergy(BondTorsions[j],iu);

    UBendTorsionTrial[iu]=0.0;
    for(j=0;j<NumberOfBendTorsions;j++)
      UBendTorsionTrial[iu]+=CalculateBendTorsionEnergy(BendTorsions[j],iu);
  }

  return 0;
}

void CheckConfigMoves(void)
{
  int i,j,l,n,iu,comp,totalb,ip;

  for(comp=0;comp<NumberOfComponents;comp++)
  {
    for(n=0;n<Components[comp].NumberOfConfigMoves;n++)
    {
      for(iu=0;iu<Components[comp].NumberOfAtoms;iu++)
        for(l=0;l<Components[comp].NumberOfAtoms;l++)
          MoleculeTodoConnectivity[iu][l]=FALSE;

      for(l=0;l<Components[comp].NumberOfBonds;l++)
      {
        MoleculeTodoConnectivity[Components[comp].Bonds[l].A][Components[comp].Bonds[l].B]=TRUE;
        MoleculeTodoConnectivity[Components[comp].Bonds[l].B][Components[comp].Bonds[l].A]=TRUE;
      }

      for(iu=0;iu<Components[comp].NumberOfAtoms;iu++)
        for(l=0;l<Components[comp].NumberOfAtoms;l++)
          MoleculeConnectivity[iu][l]=MoleculeTodoConnectivity[iu][l];

      NumberOfBeadsAlreadyPlaced=Components[comp].NumberOfUnchangedAtomsConfig[n];
      for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
        BeadsAlreadyPlaced[i]=Components[comp].UnchangedAtomsConfig[n][i];

      for(iu=0;iu<NumberOfBeadsAlreadyPlaced;iu++)
        for(l=0;l<Components[comp].NumberOfAtoms;l++)
          MoleculeTodoConnectivity[BeadsAlreadyPlaced[iu]][l]=FALSE;

      do
      {
        totalb=0;
        for(iu=0;iu<NumberOfBeadsAlreadyPlaced;iu++)
        {
          l=BeadsAlreadyPlaced[iu];
          for(j=0;j<Components[comp].NumberOfAtoms;j++)
          {
            if(MoleculeTodoConnectivity[j][l])
            {
              beadn[totalb]=j;
              PossibleCurrentBeads[totalb]=l;
              totalb++;
            }
          }
        }
        // start always with the first one: fixed growth path
        CurrentBead=PossibleCurrentBeads[0];
        PreviousBead=-2;

        NumberOfPreviousBeads=0;
        NumberOfBeadsToBePlaced=0;
        for(iu=0;iu<Components[comp].NumberOfAtoms;iu++)
        {
          if(MoleculeConnectivity[iu][CurrentBead])
          {
            if(MoleculeTodoConnectivity[iu][CurrentBead])
              BeadsToBePlaced[NumberOfBeadsToBePlaced++]=iu;
            else
            {
              NumberOfPreviousBeads++;
              PreviousBead=iu;
            }
          }
        }

        
        if((Components[comp].Connectivity[CurrentBead]>1)&&
           ((Components[comp].Connectivity[CurrentBead]-NumberOfBeadsToBePlaced)>1))
        {
          printf("Error in config move %d of component: %d (%s)\n",n,comp,Components[comp].Name);
          printf("All branches need to be grown at the same time, this config move keeps two branches fixed\n");
          exit(0);
        }

        for(iu=0;iu<NumberOfBeadsToBePlaced;iu++)
        {
          j=BeadsToBePlaced[iu];
          BeadsAlreadyPlaced[NumberOfBeadsAlreadyPlaced++]=j;

          for(ip=0;ip<Components[comp].NumberOfAtoms;ip++)
             MoleculeTodoConnectivity[j][ip]=FALSE;
        }
      } while(NumberOfBeadsAlreadyPlaced!=Components[comp].NumberOfAtoms);
    }
  }
}

// determines wich atoms are bonded to eachother for
// a) the whole molecule: 'MoleculeConnectivity'-matrix
// b) the part of the molecule that needs to be regrown: 'MoleculeTodoConnectivity'-matrix
void SetConnectivityMatrix(void)
{
  int i,j;

  // set the bond-connectivity to all FALSE at first
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
      MoleculeConnectivity[i][j]=FALSE;

  // loop over the bonds and set the appropriate elements (i,j) to TRUE when i is bonded to j
  for(j=0;j<Components[CurrentComponent].NumberOfBonds;j++)
  {
    MoleculeConnectivity[Components[CurrentComponent].Bonds[j].A][Components[CurrentComponent].Bonds[j].B]=TRUE;
    MoleculeConnectivity[Components[CurrentComponent].Bonds[j].B][Components[CurrentComponent].Bonds[j].A]=TRUE;
  }

  // copy 'MoleculeTodoConnectivity' to 'MoleculeConnectivity'
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
      MoleculeTodoConnectivity[i][j]=MoleculeConnectivity[i][j];

  // set elements to FALSE for beads that are already placed (see CBMC-move which only regrows a part of a molecule)
  // the matrix can now be used to determine future atoms that still have to be placed
  for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
    for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
      MoleculeTodoConnectivity[BeadsAlreadyPlaced[i]][j]=FALSE;
}

void SetGrowingStatus(void)
{
  int i,j,k;
  int atom_nr,nchoice;
  int IsGroupAttached[1000];

  // set 'lplace' all FALSE
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    lplace[i]=FALSE;

  // set element 'lplace' to TRUE i already placed
  for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
    lplace[BeadsAlreadyPlaced[i]]=TRUE;

  nchoice=0;
  for(i=0;i<NumberOfBeadsAlreadyPlaced;i++)
  {
    k=BeadsAlreadyPlaced[i];
    for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
    {
      if(MoleculeTodoConnectivity[j][k])
      {
        beadn[nchoice]=j;
        PossibleCurrentBeads[nchoice]=k;
        nchoice++;
      }
    }
  }

  if(nchoice==0)
  {
    printf("Error in CBMC growing scheme.. No atoms can be grown, check the connectivity of your molecule.\n");
    exit(2);
  }

  // start always with the first one possible atom, i.e. a fixed growth path in both 'grow' and 'retrace'
  CurrentBead=PossibleCurrentBeads[0];
  PreviousBead=-2;

  NumberOfBranches=0;
  for(i=0;i<MaxNumberOfBeads;i++)
  {
    NumberOfBranchAtoms[i]=0;
    IsGroupAttached[i]=FALSE;
  }

  // determine the beads to be placed and the number of previous beads
  // the number of previous beads is usually 1, but can be >1 coming of a rigid unit and growing a flexible bond
  NumberOfPreviousBeads=0;
  NumberOfBeadsToBePlaced=0;
  for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
  {
    if(MoleculeConnectivity[i][CurrentBead])
    {
      if(MoleculeTodoConnectivity[i][CurrentBead])
        IsGroupAttached[Components[CurrentComponent].group[i]]=TRUE;
      else
      {
        NumberOfPreviousBeads++;
        PreviousBead=i;
      }
    }
  }

  // determine the current 'group' to which the 'currentbead' belongs
  CurrentGroup=Components[CurrentComponent].group[CurrentBead];


  // determine number of branches to be grown
  NumberOfBranches=0;
  for(i=0;i<MaxNumberOfBeads;i++)
  {
    if(IsGroupAttached[i])
    {
      // 'i' is group number that is a branch

      // if the group of the current-bead is rigid and the attached group is rigid, then add all of the atoms
      if((Components[CurrentComponent].Groups[CurrentGroup].Rigid)&&(Components[CurrentComponent].Groups[i].Rigid))
      {
        // if the group is rigid, add all the elements
        for(j=0;j<Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;j++)
        {
          // get the atom number of the group-elements
          atom_nr=Components[CurrentComponent].Groups[i].Atoms[j];

          // if the atom in the group is not the 'CurrentBead' and not already placed, add then to the 'BeadsTobePlaced'-array
          if(atom_nr!=CurrentBead)
          {
            BeadsToBePlaced[NumberOfBeadsToBePlaced++]=atom_nr;
            BranchAtoms[NumberOfBranches][NumberOfBranchAtoms[NumberOfBranches]++]=atom_nr;
          }
        }
        NumberOfBranches++;
      }
      else
      {
        // if the group is flexible, add only the atom connected to the current bead
        for(j=0;j<Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;j++)
        {
          // get the atom number of the group-elements
          atom_nr=Components[CurrentComponent].Groups[i].Atoms[j];

          // if the atom in the group is not the 'CurrentBead' and not already placed, add then to the 'BeadsTobePlaced'-array
          if(MoleculeTodoConnectivity[atom_nr][CurrentBead])
          {
            BeadsToBePlaced[NumberOfBeadsToBePlaced++]=atom_nr;
            BranchAtoms[NumberOfBranches][NumberOfBranchAtoms[NumberOfBranches]++]=atom_nr;
            NumberOfBranches++;
          }
        }
      }
    }
  }

  if(NumberOfPreviousBeads==0)
  {
    if(Components[CurrentComponent].Groups[CurrentGroup].Rigid)
    {
      NumberOfBeadsToBePlaced=0;
      for(i=0;i<Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms;i++)
        if(Components[CurrentComponent].Groups[CurrentGroup].Atoms[i]!=CurrentBead)
          BeadsToBePlaced[NumberOfBeadsToBePlaced++]=Components[CurrentComponent].Groups[CurrentGroup].Atoms[i];
    }
    else
    {
      NumberOfBeadsToBePlaced=1;
      BeadsToBePlaced[0]=beadn[0];
    }
  }
  else
  {
    if(PreviousBead==-2)
    {
      printf("Error: no prev and PreviousBead=-2\n");
      exit(2);
    }
  }

  // if there is a 'previousbead' determine the 'group' to which it belongs
  if(PreviousBead>=0)
    PreviousGroup=Components[CurrentComponent].group[PreviousBead];
  else
    PreviousGroup=-1;
}

int Rosen(void)
{
  int iu,ip,iwalk,j;
  REAL weight;

  SetConnectivityMatrix();

  // loop until all atoms of the molecule are placed
  do
  {
    SetGrowingStatus();

    Interactions();

    // fill trialpositions of the beads that are already grown
    for(iu=0;iu<NumberOfTrialPositions;iu++)
      for(j=0;j<NumberOfBeadsAlreadyPlaced;j++)
        TrialPositions[iu][BeadsAlreadyPlaced[j]]=NewPosition[CurrentSystem][BeadsAlreadyPlaced[j]];

    for(j=0;j<NumberOfTrialPositionsTorsion;j++)
      RosenbluthTorsion[j]=1.0;

    if(NumberOfPreviousBeads==0)
      GenerateTrialOrientationsSimpleSphere(FALSE);
    else
      GenerateTrialOrientationsMCScheme(FALSE);

    ComputeExternalEnergies();

    if(OVERLAP) return 2;

    weight=ComputeNormalizedRosenbluthWeight(Bfac,Trover,NumberOfTrialPositions);
    iwalk=SelectTrialPosition(Bfac,Trover,NumberOfTrialPositions);

    RosenbluthNew*=weight*RosenbluthTorsion[iwalk];

    if(RosenbluthNew<MinimumRosenbluthFactor)
    {
      OVERLAP=TRUE;
      return 2;
    }

    UBondNew[CurrentSystem]+=UBondTrial[iwalk];
    UBendNew[CurrentSystem]+=UBendTrial[iwalk];
    UBendBendNew[CurrentSystem]+=UBendBendTrial[iwalk];
    UInversionBendNew[CurrentSystem]+=UInversionBendTrial[iwalk];
    UUreyBradleyNew[CurrentSystem]+=UUreyBradleyTrial[iwalk];
    UTorsionNew[CurrentSystem]+=UTorsionTrial[iwalk];
    UImproperTorsionNew[CurrentSystem]+=UImproperTorsionTrial[iwalk];
    UBondBondNew[CurrentSystem]+=UBondBondTrial[iwalk];
    UBondBendNew[CurrentSystem]+=UBondBendTrial[iwalk];
    UBondTorsionNew[CurrentSystem]+=UBondTorsionTrial[iwalk];
    UBendTorsionNew[CurrentSystem]+=UBendTorsionTrial[iwalk];
    UIntraVDWNew[CurrentSystem]+=UIntraVDWTrial[iwalk];
    UIntraChargeChargeNew[CurrentSystem]+=UIntraChargeChargeTrial[iwalk];
    UIntraChargeBondDipoleNew[CurrentSystem]+=UIntraChargeBondDipoleTrial[iwalk];
    UIntraBondDipoleBondDipoleNew[CurrentSystem]+=UIntraBondDipoleBondDipoleTrial[iwalk];

    UCationVDWNew[CurrentSystem]+=UCationVDWTrial[iwalk];
    UAdsorbateVDWNew[CurrentSystem]+=UAdsorbateVDWTrial[iwalk];
    UHostVDWNew[CurrentSystem]+=UHostVDWTrial[iwalk];
    UHostChargeChargeNew[CurrentSystem]+=UHostChargeChargeTrial[iwalk];
    UAdsorbateChargeChargeNew[CurrentSystem]+=UAdsorbateChargeChargeTrial[iwalk];
    UCationChargeChargeNew[CurrentSystem]+=UCationChargeChargeTrial[iwalk];
    UHostChargeBondDipoleNew[CurrentSystem]+=UHostChargeBondDipoleTrial[iwalk];
    UAdsorbateChargeBondDipoleNew[CurrentSystem]+=UAdsorbateChargeBondDipoleTrial[iwalk];
    UCationChargeBondDipoleNew[CurrentSystem]+=UCationChargeBondDipoleTrial[iwalk];
    UHostBondDipoleBondDipoleNew[CurrentSystem]+=UHostBondDipoleBondDipoleTrial[iwalk];
    UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleTrial[iwalk];
    UCationBondDipoleBondDipoleNew[CurrentSystem]+=UCationBondDipoleBondDipoleTrial[iwalk];

    for(iu=0;iu<NumberOfBeadsToBePlaced;iu++)
    {
      j=BeadsToBePlaced[iu];

      // update the selected trial-position to the 'NewPosition'
      NewPosition[CurrentSystem][j]=TrialPositions[iwalk][j];

      // update the 'BeadsAlreadyPlaced'
      BeadsAlreadyPlaced[NumberOfBeadsAlreadyPlaced++]=j;

      // update the connectivity of the remainder of the molecule
      for(ip=0;ip<Components[CurrentComponent].NumberOfAtoms;ip++)
        MoleculeTodoConnectivity[j][ip]=FALSE;
    }
  } while(NumberOfBeadsAlreadyPlaced!=Components[CurrentComponent].NumberOfAtoms);

  return 0;
}

int RosenOld(void)
{
  int iu,ip,j,k;
  REAL weight;

  SetConnectivityMatrix();

  // loop until all atoms of the molecule are placed
  do
  {
    SetGrowingStatus();

    Interactions();

    // fill trialpositions of the beads that are already grown
    if(Components[CurrentComponent].ExtraFrameworkMolecule)
    {
      for(iu=0;iu<NumberOfTrialPositions;iu++)
        for(k=0;k<Components[CurrentComponent].NumberOfAtoms;k++)
          TrialPositions[iu][k]=Cations[CurrentSystem][CurrentCationMolecule].Atoms[k].Position;
    }
    else
    {
      for(iu=0;iu<NumberOfTrialPositions;iu++)
        for(k=0;k<Components[CurrentComponent].NumberOfAtoms;k++)
          TrialPositions[iu][k]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[k].Position;
    }

    for(j=0;j<NumberOfTrialPositionsTorsion;j++)
      RosenbluthTorsion[j]=1.0;

    if(NumberOfPreviousBeads==0)
      GenerateTrialOrientationsSimpleSphere(TRUE);
    else
      GenerateTrialOrientationsMCScheme(TRUE);

    ComputeExternalEnergies();

    if (OVERLAP) return 2;

    weight=ComputeNormalizedRosenbluthWeight(Bfac,Trover,NumberOfTrialPositions);
    RosenbluthOld*=weight*RosenbluthTorsion[0];

    UBondOld[CurrentSystem]+=UBondTrial[0];
    UUreyBradleyOld[CurrentSystem]+=UUreyBradleyTrial[0];
    UBendOld[CurrentSystem]+=UBendTrial[0];
    UBendBendOld[CurrentSystem]+=UBendBendTrial[0];
    UInversionBendOld[CurrentSystem]+=UInversionBendTrial[0];
    UTorsionOld[CurrentSystem]+=UTorsionTrial[0];
    UImproperTorsionOld[CurrentSystem]+=UImproperTorsionTrial[0];
    UBondBondOld[CurrentSystem]+=UBondBondTrial[0];
    UBondBendOld[CurrentSystem]+=UBondBendTrial[0];
    UBondTorsionOld[CurrentSystem]+=UBondTorsionTrial[0];
    UBendTorsionOld[CurrentSystem]+=UBendTorsionTrial[0];
    UIntraVDWOld[CurrentSystem]+=UIntraVDWTrial[0];
    UIntraChargeChargeOld[CurrentSystem]+=UIntraChargeChargeTrial[0];
    UIntraChargeBondDipoleOld[CurrentSystem]+=UIntraChargeBondDipoleTrial[0];
    UIntraBondDipoleBondDipoleOld[CurrentSystem]+=UIntraBondDipoleBondDipoleTrial[0];

    UHostVDWOld[CurrentSystem]+=UHostVDWTrial[0];
    UAdsorbateVDWOld[CurrentSystem]+=UAdsorbateVDWTrial[0];
    UCationVDWOld[CurrentSystem]+=UCationVDWTrial[0];
    UHostChargeChargeOld[CurrentSystem]+=UHostChargeChargeTrial[0];
    UAdsorbateChargeChargeOld[CurrentSystem]+=UAdsorbateChargeChargeTrial[0];
    UCationChargeChargeOld[CurrentSystem]+=UCationChargeChargeTrial[0];
    UHostChargeBondDipoleOld[CurrentSystem]+=UHostChargeBondDipoleTrial[0];
    UAdsorbateChargeBondDipoleOld[CurrentSystem]+=UAdsorbateChargeBondDipoleTrial[0];
    UCationChargeBondDipoleOld[CurrentSystem]+=UCationChargeBondDipoleTrial[0];
    UHostBondDipoleBondDipoleOld[CurrentSystem]+=UHostBondDipoleBondDipoleTrial[0];
    UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleTrial[0];
    UCationBondDipoleBondDipoleOld[CurrentSystem]+=UCationBondDipoleBondDipoleTrial[0];

    for(iu=0;iu<NumberOfBeadsToBePlaced;iu++)
    {
      j=BeadsToBePlaced[iu];
      BeadsAlreadyPlaced[NumberOfBeadsAlreadyPlaced++]=j;

      for(ip=0;ip<Components[CurrentComponent].NumberOfAtoms;ip++)
        MoleculeTodoConnectivity[j][ip]=FALSE;
    }
  }
  while(NumberOfBeadsAlreadyPlaced!=Components[CurrentComponent].NumberOfAtoms);

  return 0;
}

REAL RetraceMolecule(int Iicode)
{
  int start,k,j;
  REAL UVDWCorrectionAdsorbate,UVDWCorrectionCation;
  REAL UVDWCorrectionFramework;

  UBondOld[CurrentSystem]=0.0;
  UBendOld[CurrentSystem]=0.0;
  UBendBendOld[CurrentSystem]=0.0;
  UInversionBendOld[CurrentSystem]=0.0;
  UUreyBradleyOld[CurrentSystem]=0.0;
  UTorsionOld[CurrentSystem]=0.0;
  UImproperTorsionOld[CurrentSystem]=0.0;
  UBondBondOld[CurrentSystem]=0.0;
  UBondBendOld[CurrentSystem]=0.0;
  UBondTorsionOld[CurrentSystem]=0.0;
  UBendTorsionOld[CurrentSystem]=0.0;
  UIntraVDWOld[CurrentSystem]=0.0;
  UIntraChargeChargeOld[CurrentSystem]=0.0;
  UIntraChargeBondDipoleOld[CurrentSystem]=0.0;
  UIntraBondDipoleBondDipoleOld[CurrentSystem]=0.0;

  UHostVDWOld[CurrentSystem]=0.0;
  UAdsorbateVDWOld[CurrentSystem]=0.0;
  UCationVDWOld[CurrentSystem]=0.0;
  UHostChargeChargeOld[CurrentSystem]=0.0;
  UAdsorbateChargeChargeOld[CurrentSystem]=0.0;
  UCationChargeChargeOld[CurrentSystem]=0.0;
  UHostChargeBondDipoleOld[CurrentSystem]=0.0;
  UAdsorbateChargeBondDipoleOld[CurrentSystem]=0.0;
  UCationChargeBondDipoleOld[CurrentSystem]=0.0;
  UHostBondDipoleBondDipoleOld[CurrentSystem]=0.0;
  UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]=0.0;
  UCationBondDipoleBondDipoleOld[CurrentSystem]=0.0;

  RosenbluthOld=1.0;

  OVERLAP=FALSE;

  for(k=0;k<NumberOfTrialPositions;k++)
  {
    for(j=0;j<MaxNumberOfBeads;j++)
    {
      TrialPositions[k][j].x=0.0;
      TrialPositions[k][j].y=0.0;
      TrialPositions[k][j].z=0.0;
    }
  }

  start=Components[CurrentComponent].StartingBead;
  if (NumberOfBeadsAlreadyPlaced==0)
  {
    if(Components[CurrentComponent].ExtraFrameworkMolecule)
      FirstBeadPosition=Cations[CurrentSystem][CurrentCationMolecule].Atoms[start].Position;
    else
      FirstBeadPosition=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[start].Position;
    HandleFirstBead(Iicode);

    RosenbluthOld=RosenBluthFactorFirstBead;

    UCationVDWOld[CurrentSystem]=EnergyCationVDWFirstBead;
    UAdsorbateVDWOld[CurrentSystem]=EnergyAdsorbateVDWFirstBead;
    UHostVDWOld[CurrentSystem]=EnergyHostVDWFirstBead;

    UCationChargeChargeOld[CurrentSystem]=EnergyCationChargeChargeFirstBead;
    UAdsorbateChargeChargeOld[CurrentSystem]=EnergyAdsorbateChargeChargeFirstBead;
    UHostChargeChargeOld[CurrentSystem]=EnergyHostChargeChargeFirstBead;

    UCationChargeBondDipoleOld[CurrentSystem]=EnergyCationChargeBondDipoleFirstBead;
    UAdsorbateChargeBondDipoleOld[CurrentSystem]=EnergyAdsorbateChargeBondDipoleFirstBead;
    UHostChargeBondDipoleOld[CurrentSystem]=EnergyHostChargeBondDipoleFirstBead;

    UCationBondDipoleBondDipoleOld[CurrentSystem]=0.0;
    UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]=0.0;
    UHostBondDipoleBondDipoleOld[CurrentSystem]=0.0;

    NumberOfBeadsAlreadyPlaced=1;
    BeadsAlreadyPlaced[0]=start;
  }

  if(Components[CurrentComponent].NumberOfAtoms!=1)
    RosenOld();

  // calculate anisotropic sites
  CalculateAnisotropicTrialPositions(CurrentComponent,TrialPositions[0],TrialAnisotropicPositionRetrace);
  UVDWCorrectionFramework=CalculateFrameworkVDWEnergyCorrection(TrialPositions[0],TrialAnisotropicPositionRetrace);
  UVDWCorrectionAdsorbate=CalculateInterVDWEnergyCorrectionAdsorbate(TrialPositions[0],TrialAnisotropicPositionRetrace,CurrentAdsorbateMolecule);
  UVDWCorrectionCation=CalculateInterVDWEnergyCorrectionCation(TrialPositions[0],TrialAnisotropicPositionRetrace,CurrentCationMolecule);
  RosenbluthOld*=exp(-Beta[CurrentSystem]*(UVDWCorrectionFramework+UVDWCorrectionAdsorbate+UVDWCorrectionCation));
  UHostVDWOld[CurrentSystem]+=UVDWCorrectionFramework;
  UAdsorbateVDWOld[CurrentSystem]+=UVDWCorrectionAdsorbate;
  UCationVDWOld[CurrentSystem]+=UVDWCorrectionCation;

  // correction for certain cross-terms that can not be handled during the growth
  RosenbluthOld*=exp(-Beta[CurrentSystem]*(UBondBondOld[CurrentSystem]+UBondBendOld[CurrentSystem]));

  // biasing using only LJ requires a correction to the Rosenbluth weight aftwards
  if(BiasingMethod==LJ_BIASING)
  {
    RosenbluthOld*=exp(-Beta[CurrentSystem]*
            (UHostChargeChargeOld[CurrentSystem]+UAdsorbateChargeChargeOld[CurrentSystem]+UCationChargeChargeOld[CurrentSystem]+
             UHostChargeBondDipoleOld[CurrentSystem]+UAdsorbateChargeBondDipoleOld[CurrentSystem]+UCationChargeBondDipoleOld[CurrentSystem]+
             UHostBondDipoleBondDipoleOld[CurrentSystem]+UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]+UCationBondDipoleBondDipoleOld[CurrentSystem]+
             UIntraChargeChargeOld[CurrentSystem]+UIntraChargeBondDipoleOld[CurrentSystem]+UIntraBondDipoleBondDipoleOld[CurrentSystem]));
  }


  return RosenbluthOld;
}

REAL GrowMolecule(int Iicode)
{
  int j,start;
  REAL UVDWCorrectionAdsorbate,UVDWCorrectionCation;
  REAL UVDWCorrectionFramework;

  UBondNew[CurrentSystem]=0.0;
  UBendNew[CurrentSystem]=0.0;
  UBendBendNew[CurrentSystem]=0.0;
  UInversionBendNew[CurrentSystem]=0.0;
  UUreyBradleyNew[CurrentSystem]=0.0;
  UTorsionNew[CurrentSystem]=0.0;
  UImproperTorsionNew[CurrentSystem]=0.0;
  UBondBondNew[CurrentSystem]=0.0;
  UBondBendNew[CurrentSystem]=0.0;
  UBondTorsionNew[CurrentSystem]=0.0;
  UBendTorsionNew[CurrentSystem]=0.0;
  UIntraVDWNew[CurrentSystem]=0.0;
  UIntraChargeChargeNew[CurrentSystem]=0.0;
  UIntraChargeBondDipoleNew[CurrentSystem]=0.0;
  UIntraBondDipoleBondDipoleNew[CurrentSystem]=0.0;

  UHostVDWNew[CurrentSystem]=0.0;
  UAdsorbateVDWNew[CurrentSystem]=0.0;
  UCationVDWNew[CurrentSystem]=0.0;
  UHostChargeChargeNew[CurrentSystem]=0.0;
  UAdsorbateChargeChargeNew[CurrentSystem]=0.0;
  UCationChargeChargeNew[CurrentSystem]=0.0;
  UHostChargeBondDipoleNew[CurrentSystem]=0.0;
  UAdsorbateChargeBondDipoleNew[CurrentSystem]=0.0;
  UCationChargeBondDipoleNew[CurrentSystem]=0.0;
  UHostBondDipoleBondDipoleNew[CurrentSystem]=0.0;
  UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]=0.0;
  UCationBondDipoleBondDipoleNew[CurrentSystem]=0.0;

  RosenbluthNew=1.0;
  OVERLAP=FALSE;

  if(NumberOfBeadsAlreadyPlaced==0)
  {
    start=Components[CurrentComponent].StartingBead;
    if(Iicode==2)
      FirstBeadPosition=NewPosition[CurrentSystem][start];
    HandleFirstBead(Iicode);
    RosenbluthNew=RosenBluthFactorFirstBead;
    if(OVERLAP) return 0;


    UCationVDWNew[CurrentSystem]=EnergyCationVDWFirstBead;
    UAdsorbateVDWNew[CurrentSystem]=EnergyAdsorbateVDWFirstBead;
    UHostVDWNew[CurrentSystem]=EnergyHostVDWFirstBead;

    UCationChargeChargeNew[CurrentSystem]=EnergyCationChargeChargeFirstBead;
    UAdsorbateChargeChargeNew[CurrentSystem]=EnergyAdsorbateChargeChargeFirstBead;
    UHostChargeChargeNew[CurrentSystem]=EnergyHostChargeChargeFirstBead;

    UCationChargeBondDipoleNew[CurrentSystem]=EnergyCationChargeBondDipoleFirstBead;
    UAdsorbateChargeBondDipoleNew[CurrentSystem]=EnergyAdsorbateChargeBondDipoleFirstBead;
    UHostChargeBondDipoleNew[CurrentSystem]=EnergyHostChargeBondDipoleFirstBead;

    UCationBondDipoleBondDipoleNew[CurrentSystem]=0.0;
    UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]=0.0;
    UHostBondDipoleBondDipoleNew[CurrentSystem]=0.0;

    NumberOfBeadsAlreadyPlaced=1;
    BeadsAlreadyPlaced[0]=start;
    NewPosition[CurrentSystem][BeadsAlreadyPlaced[0]]=FirstBeadPosition;
  }

  if(Components[CurrentComponent].NumberOfAtoms!=1)
    Rosen();

  if (OVERLAP) return 0;

  // copy coordinates for small MC-scheme
  for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
  {
    Components[CurrentComponent].RMCMOL[j]=NewPosition[CurrentSystem][j];
    TrialPosition[CurrentSystem][j]=NewPosition[CurrentSystem][j];
  }

  // calculate anisotropic sites
  CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);
  UVDWCorrectionFramework=CalculateFrameworkVDWEnergyCorrection(TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);
  UVDWCorrectionAdsorbate=CalculateInterVDWEnergyCorrectionAdsorbate(TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem],CurrentAdsorbateMolecule);
  UVDWCorrectionCation=CalculateInterVDWEnergyCorrectionCation(TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem],CurrentCationMolecule);
  RosenbluthNew*=exp(-Beta[CurrentSystem]*(UVDWCorrectionFramework+UVDWCorrectionAdsorbate+UVDWCorrectionCation));
  UHostVDWNew[CurrentSystem]+=UVDWCorrectionFramework;
  UAdsorbateVDWNew[CurrentSystem]+=UVDWCorrectionAdsorbate;
  UCationVDWNew[CurrentSystem]+=UVDWCorrectionCation;
  

  // the old config can be used as a starting point
  Components[CurrentComponent].LMCMOL=TRUE;

  // correction for certain cross-terms that can not be handled during the growth
  RosenbluthNew*=exp(-Beta[CurrentSystem]*(UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]));

  // biasing using only LJ requires a correction to the Rosenbluth weight aftwards
  if(BiasingMethod==LJ_BIASING)
  {
    RosenbluthNew*=exp(-Beta[CurrentSystem]*
            (UHostChargeChargeNew[CurrentSystem]+UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+
             UHostChargeBondDipoleNew[CurrentSystem]+UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+
             UHostBondDipoleBondDipoleNew[CurrentSystem]+UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+
             UIntraChargeChargeNew[CurrentSystem]+UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]));
  }

  return RosenbluthNew;
}

void MakeInitialAdsorbate(void)
{
  int j,k;
  int StartingBead;
  POINT s;

  do
  {
    for(k=0;k<NumberOfTrialPositions;k++)
    {
      for(j=0;j<MaxNumberOfBeads;j++)
      {
        TrialPositions[k][j].x=0.0;
        TrialPositions[k][j].y=0.0;
        TrialPositions[k][j].z=0.0;
      }
    }
    CurrentAdsorbateMolecule=NumberOfAdsorbateMolecules[CurrentSystem];
    CurrentCationMolecule=NumberOfCationMolecules[CurrentSystem];
    NumberOfBeadsAlreadyPlaced=0;
    StartingBead=Components[CurrentComponent].StartingBead;

    if(Components[CurrentComponent].RestrictMoves)
    {
      do
      {
        s.x=s.y=s.z=0.0;
        switch(Dimension)
        {
          case 3:
            s.z=RandomNumber();
          case 2:
            s.y=RandomNumber();
          case 1:
            s.x=RandomNumber();
            break;
        }
        NewPosition[CurrentSystem][StartingBead].x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
        NewPosition[CurrentSystem][StartingBead].y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
        NewPosition[CurrentSystem][StartingBead].z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;
      } while(!ValidCartesianPoint(CurrentComponent,NewPosition[CurrentSystem][StartingBead]));
      GrowMolecule(2);
    }
    else GrowMolecule(1);
  }
  while(OVERLAP==TRUE||BlockedPocket(NewPosition[CurrentSystem][Components[CurrentComponent].StartingBead]));
  InsertAdsorbateMolecule();
}

void MakeInitialAdsorbates(int n,int type)
{
  int i;

  CurrentComponent=type;
  for(i=0;i<n;i++)
    MakeInitialAdsorbate(); 
}

void MakeInitialCation(void)
{
  int j,k;
  int StartingBead;
  POINT s;

  do
  {
    for(k=0;k<NumberOfTrialPositions;k++)
    {
      for(j=0;j<MaxNumberOfBeads;j++)
      {
        TrialPositions[k][j].x=0.0;
        TrialPositions[k][j].y=0.0;
        TrialPositions[k][j].z=0.0;
      }
    }
    CurrentCationMolecule=NumberOfCationMolecules[CurrentSystem];
    CurrentAdsorbateMolecule=NumberOfAdsorbateMolecules[CurrentSystem];
    NumberOfBeadsAlreadyPlaced=0;
    StartingBead=Components[CurrentComponent].StartingBead;

    if(Components[CurrentComponent].RestrictMoves)
    {
      do
      {
        s.x=s.y=s.z=0.0;
        switch(Dimension)
        {
          case 3:
            s.z=RandomNumber();
          case 2:
            s.y=RandomNumber();
          case 1:
            s.x=RandomNumber();
            break;
        }
        // convert to a the Cartesian position in the simulation box
        NewPosition[CurrentSystem][StartingBead].x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
        NewPosition[CurrentSystem][StartingBead].y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
        NewPosition[CurrentSystem][StartingBead].z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;
      }
      while(!ValidCartesianPoint(CurrentComponent,NewPosition[CurrentSystem][StartingBead]));
      GrowMolecule(2);
    }
    else GrowMolecule(1);
  }
  while(OVERLAP==TRUE||BlockedPocket(NewPosition[CurrentSystem][Components[CurrentComponent].StartingBead]));

  InsertCationMolecule();
}

void MakeInitialCations(int n,int type)
{
  int i;

  CurrentComponent=type;
  for(i=0;i<n;i++)
    MakeInitialCation(); 
}

void RescaleMaximumRotationAnglesSmallMC(void)
{
  int i,j;
  REAL rm;

  for(i=0;i<NumberOfComponents;i++)
  {
    for(j=0;j<Components[i].NumberOfAtoms;j++)
    {
      if (Components[i].ACNBOND[CurrentSystem][j]>5000.0)
      {
        rm=(1.0/TargetAccRatioSmallMCScheme)*Components[i].ACSBOND[CurrentSystem][j]/
           Components[i].ACNBOND[CurrentSystem][j];
        if (rm>2.0) rm=2.0;
        else if (rm<0.5) rm=0.5;
        Components[i].RMAX0[CurrentSystem][j]=rm*Components[i].RMAX0[CurrentSystem][j];
        Components[i].acn0[CurrentSystem][j]+=Components[i].ACNBOND[CurrentSystem][j];
        Components[i].acs0[CurrentSystem][j]+=Components[i].ACSBOND[CurrentSystem][j];
        Components[i].ACNBOND[CurrentSystem][j]=0.0;
        Components[i].ACSBOND[CurrentSystem][j]=0.0;
        if(Components[i].RMAX0[CurrentSystem][j]>0.5) Components[i].RMAX0[CurrentSystem][j]=0.5;
        if(Components[i].RMAX0[CurrentSystem][j]<0.01) Components[i].RMAX0[CurrentSystem][j]=0.01;
      }
    }

    for(j=0;j<Components[i].NumberOfAtoms;j++)
    {
      if (Components[i].ACNBEND[CurrentSystem][j]>5000.0)
      {
        rm=(1.0/TargetAccRatioSmallMCScheme)*Components[i].ACSBEND[CurrentSystem][j]/
           Components[i].ACNBEND[CurrentSystem][j];
        if (rm>2.0) rm=2.0;
        else if (rm<0.5) rm=0.5;
        Components[i].RMAX1[CurrentSystem][j]=rm*Components[i].RMAX1[CurrentSystem][j];
        Components[i].acn1[CurrentSystem][j]+=Components[i].ACNBEND[CurrentSystem][j];
        Components[i].acs1[CurrentSystem][j]+=Components[i].ACSBEND[CurrentSystem][j];
        Components[i].ACNBEND[CurrentSystem][j]=0.0;
        Components[i].ACSBEND[CurrentSystem][j]=0.0;
        if(Components[i].RMAX1[CurrentSystem][j]>M_PI) Components[i].RMAX1[CurrentSystem][j]=M_PI;
      }
    }

    for(j=0;j<Components[i].NumberOfAtoms;j++)
    {
      if (Components[i].ACNTRS[CurrentSystem][j]>5000.0)
      {
        rm=(1.0/TargetAccRatioSmallMCScheme)*Components[i].ACSTRS[CurrentSystem][j]/
           Components[i].ACNTRS[CurrentSystem][j];
        if (rm>2.0) rm = 2.0;
        else if(rm<0.5) rm=0.5;
        Components[i].RMAX2[CurrentSystem][j]=rm*Components[i].RMAX2[CurrentSystem][j];
        Components[i].acn2[CurrentSystem][j]+=Components[i].ACNTRS[CurrentSystem][j];
        Components[i].acs2[CurrentSystem][j]+=Components[i].ACSTRS[CurrentSystem][j];
        Components[i].ACNTRS[CurrentSystem][j]=0.0;
        Components[i].ACSTRS[CurrentSystem][j]=0.0;
        if(Components[i].RMAX2[CurrentSystem][j]>M_PI) Components[i].RMAX2[CurrentSystem][j]=M_PI;
      }
    }
  }
}

void InitializeSmallMCStatisticsAllSystems(void)
{
  int i,j,k;

  for(i=0;i<NumberOfComponents;i++)
    for(j=0;j<Components[i].NumberOfAtoms;j++)
      for(k=0;k<NumberOfSystems;k++)
      {
        Components[i].acn0[k][j]=0.0;
        Components[i].acs0[k][j]=0.0;
        Components[i].acn1[k][j]=0.0;
        Components[i].acs1[k][j]=0.0;
        Components[i].acn2[k][j]=0.0;
        Components[i].acs2[k][j]=0.0;
        Components[i].ACNALL[k][j]=0.0;
        Components[i].ACSALL[k][j]=0.0; 
        Components[i].ACNBOND[k][j]=0.0;
        Components[i].ACSBOND[k][j]=0.0;
        Components[i].ACNBEND[k][j]=0.0;
        Components[i].ACSBEND[k][j]=0.0;
        Components[i].ACNTRS[k][j]=0.0;
        Components[i].ACSTRS[k][j]=0.0;
      }
}

void PrintSmallMCAddStatistics(FILE *FilePtr)
{
  int i,j;

  fprintf(FilePtr,"Performance of the small-mc scheme\n");
  fprintf(FilePtr,"==================================\n");
  fprintf(FilePtr,"\n");
  for(i=0;i<NumberOfComponents;i++)
  {
    fprintf(FilePtr,"Component %d [%s]\n",i,Components[i].Name);
    fprintf(FilePtr,"----------------------------------------------\n");
    for(j=0;j<Components[i].NumberOfAtoms;j++)
    {
      if(Components[i].acn0[CurrentSystem][j]>0.5)
        Components[i].acn0[CurrentSystem][j]=100.0*Components[i].acs0[CurrentSystem][j]/
                                             Components[i].acn0[CurrentSystem][j];
      else
        Components[i].acn0[CurrentSystem][j]=0.0;

      if(Components[i].acn1[CurrentSystem][j]>0.5)
        Components[i].acn1[CurrentSystem][j]=100.0*Components[i].acs1[CurrentSystem][j]/
                                             Components[i].acn1[CurrentSystem][j];
      else
        Components[i].acn1[CurrentSystem][j]=0.0;

      if(Components[i].acn2[CurrentSystem][j]>0.5)
        Components[i].acn2[CurrentSystem][j]=100.0*Components[i].acs2[CurrentSystem][j]/
                                             Components[i].acn2[CurrentSystem][j];
      else
        Components[i].acn2[CurrentSystem][j]=0.0;

      if(Components[i].ACNALL[CurrentSystem][j]>0.5)
        Components[i].ACNALL[CurrentSystem][j]=100.0*Components[i].ACSALL[CurrentSystem][j]/
                                               Components[i].ACNALL[CurrentSystem][j];
      else
        Components[i].ACNALL[CurrentSystem][j]=0.0;

      fprintf(FilePtr,"Bead: %d\n",j);
      fprintf(FilePtr,"\tmaximum length change                         : %lf\n",(double)Components[i].RMAX0[CurrentSystem][j]);
      fprintf(FilePtr,"\tmaximum rotation angle for cone angle         : %lf\n",(double)Components[i].RMAX1[CurrentSystem][j]);
      fprintf(FilePtr,"\tmaximum rotation angle for position on a cone : %lf\n",(double)Components[i].RMAX2[CurrentSystem][j]);
      fprintf(FilePtr,"\tlength change acceptence ratio                : %lf [%%]\n",(double)Components[i].acn0[CurrentSystem][j]);
      fprintf(FilePtr,"\trotation cone angle acceptence ratio          : %lf [%%]\n",(double)Components[i].acn1[CurrentSystem][j]);
      fprintf(FilePtr,"\trotation angle cone position acceptence ratio : %lf [%%]\n",(double)Components[i].acn2[CurrentSystem][j]);
      fprintf(FilePtr,"\toverall acceptence ratio                      : %lf [%%]\n",(double)Components[i].ACNALL[CurrentSystem][j]);
      fprintf(FilePtr,"\n");
    }
  }
  fprintf(FilePtr,"\n\n");
}

void WriteRestartCBMC(FILE *FilePtr)
{
  fwrite(&BiasingMethod,sizeof(BiasingMethod),1,FilePtr);
  fwrite(&NumberOfTrialPositions,sizeof(NumberOfTrialPositions),1,FilePtr);
  fwrite(&NumberOfTrialPositionsForTheFirstBead,sizeof(NumberOfTrialPositionsForTheFirstBead),1,FilePtr);
  fwrite(&NumberOfTrialMovesPerOpenBead,sizeof(NumberOfTrialMovesPerOpenBead),1,FilePtr);
  fwrite(&NumberOfTrialPositionsTorsion,sizeof(NumberOfTrialPositionsTorsion),1,FilePtr);
  fwrite(&TargetAccRatioSmallMCScheme,sizeof(TargetAccRatioSmallMCScheme),1,FilePtr);
  fwrite(&EnergyOverlapCriteria,sizeof(EnergyOverlapCriteria),1,FilePtr);
  fwrite(&MinimumRosenbluthFactor,sizeof(MinimumRosenbluthFactor),1,FilePtr);

  fwrite(&MaxNumberOfBeads,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfBonds,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfBondDipoles,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfBends,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfBendBends,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfInversionBends,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfUreyBradleys,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfTorsions,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfImproperTorsions,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfBondBonds,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfBondBends,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfBondTorsions,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfBendTorsions,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfIntraVDW,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfIntraChargeCharge,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfIntraChargeBondDipole,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfIntraBondDipoleBondDipole,sizeof(int),1,FilePtr);
}

void AllocateCBMCMemory(void)
{
  int i;
  int MaxNumberOfIntraMolecularInteractions,MaxTrial;

  MaxTrial=MAX3(NumberOfTrialPositions,NumberOfTrialPositionsTorsion,NumberOfTrialPositionsForTheFirstBead);

  RosenbluthTorsion=(REAL*)calloc(NumberOfTrialPositionsTorsion,sizeof(REAL));

  BeadsToBePlaced=(int*)calloc(MaxNumberOfBeads,sizeof(int));
  BeadsAlreadyPlaced=(int*)calloc(MaxNumberOfBeads,sizeof(int));

  NumberOfBranchAtoms=(int*)calloc(MaxNumberOfBeads,sizeof(int));
  BranchAtoms=(int**)calloc(MaxNumberOfBeads,sizeof(int*));
  for(i=0;i<MaxNumberOfBeads;i++)
    BranchAtoms[i]=(int*)calloc(MaxNumberOfBeads,sizeof(int));

  NewPosition=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  NewVelocity=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  NewForce=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));

  for(i=0;i<NumberOfSystems;i++)
  {
    NewPosition[i]=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));
    NewVelocity[i]=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));
    NewForce[i]=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));
  }

  TrialAnisotropicPositionRetrace=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));

  TrialPositions=(VECTOR**)calloc(NumberOfTrialPositions,sizeof(VECTOR*));
  for(i=0;i<NumberOfTrialPositions;i++)
    TrialPositions[i]=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));

  Bfac=(REAL*)calloc(MaxTrial,sizeof(REAL));
  Trover=(int*)calloc(MaxTrial,sizeof(REAL));
  ShiftedBoltzmannFactors=(REAL*)calloc(MaxTrial,sizeof(REAL));

  Trial=(VECTOR*)calloc(MaxTrial,sizeof(VECTOR));

  EnergiesHostVDW=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesAdsorbateVDW=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesCationVDW=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesHostChargeCharge=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesAdsorbateChargeCharge=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesCationChargeCharge=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesHostChargeBondDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesAdsorbateChargeBondDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesCationChargeBondDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesHostChargePermanentDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesAdsorbateChargePermanentDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesCationChargePermanentDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesHostBondDipoleBondDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesAdsorbateBondDipoleBondDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesCationBondDipoleBondDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesHostBondDipolePermanentDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesAdsorbateBondDipolePermanentDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesCationBondDipolePermanentDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesHostPermanentDipolePermanentDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesAdsorbatePermanentDipolePermanentDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));
  EnergiesCationPermanentDipolePermanentDipole=(REAL*)calloc(MaxTrial,sizeof(REAL));

  lplace=(int*)calloc(MaxNumberOfBeads,sizeof(int));
  lchose=(int*)calloc(MaxNumberOfBeads,sizeof(int));

  bond_length=(REAL*)calloc(MaxNumberOfBeads,sizeof(REAL));
  cord=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));
  store=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));
  angle=(REAL*)calloc(MaxTrial,sizeof(REAL));
  perpen=(VECTOR*)calloc(MaxNumberOfBeads,sizeof(VECTOR));

  beadn=(int*)calloc(MaxNumberOfBeads,sizeof(int));
  PossibleCurrentBeads=(int*)calloc(MaxNumberOfBeads,sizeof(int));
  MoleculeTodoConnectivity=(int**)calloc(MaxNumberOfBeads,sizeof(int*));
  MoleculeConnectivity=(int**)calloc(MaxNumberOfBeads,sizeof(int*));
  for(i=0;i<MaxNumberOfBeads;i++)
  {
    MoleculeTodoConnectivity[i]=(int*)calloc(MaxNumberOfBeads,sizeof(int));
    MoleculeConnectivity[i]=(int*)calloc(MaxNumberOfBeads,sizeof(int));
  }

  MaxNumberOfIntraMolecularInteractions=0;
  if(MaxNumberOfBonds>MaxNumberOfIntraMolecularInteractions) MaxNumberOfIntraMolecularInteractions=MaxNumberOfBonds;
  if(MaxNumberOfBondDipoles>MaxNumberOfIntraMolecularInteractions) MaxNumberOfIntraMolecularInteractions=MaxNumberOfBondDipoles;
  if(MaxNumberOfBends>MaxNumberOfIntraMolecularInteractions) MaxNumberOfIntraMolecularInteractions=MaxNumberOfBends;
  if(MaxNumberOfBendBends>MaxNumberOfIntraMolecularInteractions) MaxNumberOfIntraMolecularInteractions=MaxNumberOfBendBends;
  if(MaxNumberOfUreyBradleys>MaxNumberOfIntraMolecularInteractions) MaxNumberOfIntraMolecularInteractions=MaxNumberOfUreyBradleys;
  if(MaxNumberOfInversionBends>MaxNumberOfIntraMolecularInteractions) MaxNumberOfIntraMolecularInteractions=MaxNumberOfInversionBends;
  if(MaxNumberOfBondBonds>MaxNumberOfIntraMolecularInteractions) MaxNumberOfIntraMolecularInteractions=MaxNumberOfBondBonds;
  if(MaxNumberOfBondBends>MaxNumberOfIntraMolecularInteractions) MaxNumberOfIntraMolecularInteractions=MaxNumberOfBondBends;
  if(MaxNumberOfBondTorsions>MaxNumberOfIntraMolecularInteractions) MaxNumberOfIntraMolecularInteractions=MaxNumberOfBondTorsions;
  if(MaxNumberOfBendTorsions>MaxNumberOfIntraMolecularInteractions) MaxNumberOfIntraMolecularInteractions=MaxNumberOfBendTorsions;

  enbend=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  eobend=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  enbendbend=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  eobendbend=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  enbond=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  eobond=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  enureybradley=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  eoureybradley=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  eninversionbend=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  eoinversionbend=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  enimpropertorsion=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));
  eoimpropertorsion=(REAL*)calloc(MaxNumberOfIntraMolecularInteractions,sizeof(REAL));

  Bonds=(int*)calloc(MaxNumberOfBonds,sizeof(int));
  BondDipoles=(int*)calloc(MaxNumberOfBondDipoles,sizeof(int));
  Bends=(int*)calloc(MaxNumberOfBends,sizeof(int));
  BendBends=(int*)calloc(MaxNumberOfBendBends,sizeof(int));
  UreyBradleys=(int*)calloc(MaxNumberOfUreyBradleys,sizeof(int));
  InversionBends=(int*)calloc(MaxNumberOfInversionBends,sizeof(int));
  Torsions=(int*)calloc(MaxNumberOfTorsions,sizeof(int));
  ImproperTorsions=(int*)calloc(MaxNumberOfImproperTorsions,sizeof(int));
  BondBonds=(int*)calloc(MaxNumberOfBondBonds,sizeof(int));
  BondBends=(int*)calloc(MaxNumberOfBondBends,sizeof(int));
  BondTorsions=(int*)calloc(MaxNumberOfBondTorsions,sizeof(int));
  BendTorsions=(int*)calloc(MaxNumberOfBendTorsions,sizeof(int));
  VDW=(int*)calloc(MaxNumberOfIntraVDW,sizeof(int));
  IntraChargeCharge=(int*)calloc(MaxNumberOfIntraChargeCharge,sizeof(int));
  IntraChargeBondDipole=(int*)calloc(MaxNumberOfIntraChargeBondDipole,sizeof(int));
  IntraBondDipoleBondDipole=(int*)calloc(MaxNumberOfIntraBondDipoleBondDipole,sizeof(int));

  // grow energies
  UBondNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBendNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBendBendNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UUreyBradleyNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UInversionBendNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UTorsionNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UImproperTorsionNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBondBondNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBondBendNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBondTorsionNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBendTorsionNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UIntraVDWNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UIntraChargeChargeNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UIntraChargeBondDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UIntraBondDipoleBondDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UAdsorbateVDWNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationVDWNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostVDWNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateChargeChargeNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationChargeChargeNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostChargeChargeNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateChargeBondDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationChargeBondDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostChargeBondDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateChargePermanentDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationChargePermanentDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostChargePermanentDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBondDipoleBondDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBondDipoleBondDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBondDipoleBondDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBondDipolePermanentDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBondDipolePermanentDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBondDipolePermanentDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbatePermanentDipolePermanentDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationPermanentDipolePermanentDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostPermanentDipolePermanentDipoleNew=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // retrace energies
  UBondOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBendOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBendBendOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UUreyBradleyOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UInversionBendOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UTorsionOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UImproperTorsionOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBondBondOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBondBendOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBondTorsionOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBendTorsionOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UIntraVDWOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UIntraChargeChargeOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UIntraChargeBondDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UIntraBondDipoleBondDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UAdsorbateVDWOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationVDWOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostVDWOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateChargeChargeOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationChargeChargeOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostChargeChargeOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateChargeBondDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationChargeBondDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostChargeBondDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateChargePermanentDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationChargePermanentDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostChargePermanentDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBondDipoleBondDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBondDipoleBondDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBondDipoleBondDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateBondDipolePermanentDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationBondDipolePermanentDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostBondDipolePermanentDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbatePermanentDipolePermanentDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationPermanentDipolePermanentDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostPermanentDipolePermanentDipoleOld=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // trial energies
  UBondTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UBendTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UBendBendTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UUreyBradleyTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UInversionBendTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UTorsionTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UImproperTorsionTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UBondBondTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UBondBendTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UBondTorsionTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UBendTorsionTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UIntraVDWTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UIntraChargeChargeTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UIntraChargeBondDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UIntraBondDipoleBondDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));

  UHostVDWTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UAdsorbateVDWTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UCationVDWTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UHostChargeChargeTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UAdsorbateChargeChargeTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UCationChargeChargeTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UHostChargeBondDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UAdsorbateChargeBondDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UCationChargeBondDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UHostChargePermanentDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UAdsorbateChargePermanentDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UCationChargePermanentDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UHostBondDipoleBondDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UAdsorbateBondDipoleBondDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UCationBondDipoleBondDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UHostBondDipolePermanentDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UAdsorbateBondDipolePermanentDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UCationBondDipolePermanentDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UHostPermanentDipolePermanentDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UAdsorbatePermanentDipolePermanentDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
  UCationPermanentDipolePermanentDipoleTrial=(REAL*)calloc(NumberOfTrialPositions,sizeof(REAL));
}

void ReadRestartCBMC(FILE *FilePtr)
{
  fread(&BiasingMethod,sizeof(BiasingMethod),1,FilePtr);
  fread(&NumberOfTrialPositions,sizeof(NumberOfTrialPositions),1,FilePtr);
  fread(&NumberOfTrialPositionsForTheFirstBead,sizeof(NumberOfTrialPositionsForTheFirstBead),1,FilePtr);
  fread(&NumberOfTrialMovesPerOpenBead,sizeof(NumberOfTrialMovesPerOpenBead),1,FilePtr);
  fread(&NumberOfTrialPositionsTorsion,sizeof(NumberOfTrialPositionsTorsion),1,FilePtr);
  fread(&TargetAccRatioSmallMCScheme,sizeof(TargetAccRatioSmallMCScheme),1,FilePtr);
  fread(&EnergyOverlapCriteria,sizeof(EnergyOverlapCriteria),1,FilePtr);
  fread(&MinimumRosenbluthFactor,sizeof(MinimumRosenbluthFactor),1,FilePtr);

  fread(&MaxNumberOfBeads,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfBonds,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfBondDipoles,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfBends,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfBendBends,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfInversionBends,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfUreyBradleys,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfTorsions,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfImproperTorsions,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfBondBonds,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfBondBends,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfBondTorsions,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfBendTorsions,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfIntraVDW,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfIntraChargeCharge,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfIntraChargeBondDipole,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfIntraBondDipoleBondDipole,sizeof(int),1,FilePtr);

  AllocateCBMCMemory();
}
