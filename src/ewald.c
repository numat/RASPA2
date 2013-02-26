/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'ewald.c' is part of RASPA.

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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "simulation.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "inter_hessian.h"
#include "utils.h"
#include "ewald.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "potentials.h"
#include "spectra.h"
#include "rigid.h"
#include "minimization.h"


// Some care is on order when dealing with flexible/rigid frameworks and net-charges.
// The whole system can be 'flexible' or 'rigid'. It is considered 'rigid' when all the
// frameworks in the system are 'rigid' and the FrameworkShiftMoveProbability is zero.
// Otherwise the system is considered as 'flexible'. For the flexible system the 
// individual frameworks can be either 'rigid' or 'flexible'. For a flexible system with
// several rigid frameworks, the framework distances can change but not their internal
// structure. Flexible frameworks do change internally.
// For a rigid system one can speed up computation in the Ewald-summation by assuming
// that no framework particle is allowed to move. Therefore their Ewald-sums can be
// precomputed.
// Even when the system is charge-neutral we need to pay attention to net-charges. For
// example when we have an aluminium subsituted framework with mobile cations. The cations
// counter the negative charge of the framework. The interaction energies are corrected
// for the net-charges as well as any net-charge of the whole system. Hence, a periodic
// system is allowed to have net-charges even when using the Ewald summation.

// interactions of A and B
// A-A:  Ewald + (netcharge A)^2xUIon - (alpha/sqrt(pi))* sum_i^A q_i^2
// B-B:  Ewald + (netcharge B)^2xUIon - (alpha/sqrt(pi))* sum_i^B q_i^2
// A-B:  Ewald + (netcharge A)x(netcharge B)xUIon


// the current net charges per system
REAL *NetChargeSystem;
REAL *NetChargeFramework;
REAL *NetChargeCations;
REAL *NetChargeAdsorbates;

// the squered cutoff in Fourier space
REAL *ReciprocalCutOffSquared;

// the net-charge difference for MC-moves
REAL NetChargeAdsorbateDelta;
REAL NetChargeCationDelta;

// the energy difference for MC-moves
REAL *UHostHostChargeChargeFourierDelta;
REAL *UAdsorbateAdsorbateChargeChargeFourierDelta;
REAL *UCationCationChargeChargeFourierDelta;
REAL *UAdsorbateCationChargeChargeFourierDelta;
REAL *UHostCationChargeChargeFourierDelta;
REAL *UHostAdsorbateChargeChargeFourierDelta;

REAL *UHostHostChargeBondDipoleFourierDelta;
REAL *UAdsorbateAdsorbateChargeBondDipoleFourierDelta;
REAL *UCationCationChargeBondDipoleFourierDelta;
REAL *UAdsorbateCationChargeBondDipoleFourierDelta;
REAL *UHostCationChargeBondDipoleFourierDelta;
REAL *UHostAdsorbateChargeBondDipoleFourierDelta;

REAL *UHostHostBondDipoleBondDipoleFourierDelta;
REAL *UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta;
REAL *UCationCationBondDipoleBondDipoleFourierDelta;
REAL *UAdsorbateCationBondDipoleBondDipoleFourierDelta;
REAL *UHostCationBondDipoleBondDipoleFourierDelta;
REAL *UHostAdsorbateBondDipoleBondDipoleFourierDelta;

int *NumberOfKVectors;              // number of wave-vectors for a system
static int MaxNumberOfWaveVectors;  // the maximum number of wave-vectors (over all systems)
static int MaxKvecX;                // the maximum number of wavevectors in the x-direction
static int MaxKvecY;                // the maximum number of wavevectors in the y-direction
static int MaxKvecZ;                // the maximum number of wavevectors in the z-direction
REAL *Alpha;                        // the Ewald alpha-parameter
INT_VECTOR3 *kvec;                  // the number of wavevectors in x,y,z per system
REAL EwaldPrecision;                // the relative precision of the Ewald computation
REAL DielectricConstantOfTheMedium; // the dielectric constant of the medium
int OmitEwaldFourier;               // option to skip the Fourier-part of the Ewald-summation

// stored structure factors of the rigid atoms
static COMPLEX **StoreRigidChargeFramework;
static COMPLEX **StoreRigidChargeAdsorbates;
static COMPLEX **StoreRigidChargeCations;

static COMPLEX **StoreRigidBondDipolesFramework;
static COMPLEX **StoreRigidBondDipolesAdsorbates;
static COMPLEX **StoreRigidBondDipolesCations;

// stored total structure factors
static COMPLEX **StoreTotalChargeFramework;
static COMPLEX **StoreTotalChargeAdsorbates;
static COMPLEX **StoreTotalChargeCations;

static COMPLEX **StoreTotalBondDipolesFramework;
static COMPLEX **StoreTotalBondDipolesAdsorbates;
static COMPLEX **StoreTotalBondDipolesCations;

// stored new structure factors for trial MC-moves
static COMPLEX **NewTotalChargeFramework;
static COMPLEX **NewTotalChargeAdsorbates;
static COMPLEX **NewTotalChargeCations;

static COMPLEX **NewTotalBondDipolesFramework;
static COMPLEX **NewTotalBondDipolesAdsorbates;
static COMPLEX **NewTotalBondDipolesCations;

// storage for charges
static COMPLEX *Eikx;                // exp(-ik.x) for k=0,..,kx (or k=-kx,..,kx when no symmetry is used, e.g. phonon disperion)
static COMPLEX *Eiky;                // exp(-ik.y) for k=-ky,..,ky
static COMPLEX *Eikz;                // exp(-ik.z) for k=-kz,..,kz
static COMPLEX *Eikr;                // exp(-ik.r) for i=1,..,N
static COMPLEX *Eikr_xy;             // temporary storage for exp(-ik.x)*exp(-ik.y) for i=1,..,N

// storage for bond-dipoles
static COMPLEX *Eikx_bd;             // exp(-ik.x) for k=0,..,kx (or k=-kx,..,kx when no symmetry is used, e.g. phonon disperion)
static COMPLEX *Eiky_bd;             // exp(-ik.y) for k=-ky,..,ky
static COMPLEX *Eikz_bd;             // exp(-ik.z) for k=-kz,..,kz
static COMPLEX *Eikr_bd;             // exp(-ik.r) for i=1,..,N
static COMPLEX *Eikr_xy_bd;          // temporary storage for exp(-ik.x)*exp(-ik.y) for i=1,..,N

static VECTOR **KVectors;            // the k-vectors (from 0,...,nvec-1) for each system
static REAL **KFactor;               // precomputed pre-factor COULOMBIC_CONVERSION_FACTOR*(4*Pi/V)*exp(-k^2/(4*SQR(alpha)))/k^2

static VECTOR **StoredKVectors;      // the k-vectors (from 0,...,nvec-1) for each system
static REAL **StoredKFactor;         // precomputed pre-factor

static VECTOR *Positions;            // the fractional charge positions times 2*Pi
static VECTOR **AtomVector;          // list of pointers to atom forces or electric field
static VECTOR *BondDipolePositions;  // the fractional bond-dipole positions times 2*Pi
static REAL *Charge;                 // the charges of the atoms
static REAL *Polarization;           // the polarization of the atoms
static VECTOR **BondDipoleForcesA;   // list of pointers to atom forces for the first atom of a bond-dipole
static VECTOR **BondDipoleForcesB;   // list of pointers to atom forces for the second atom of a bond-dipole
static VECTOR  *DipoleVector;        // the dipole-vector
static REAL *BondLength;             // the bond-length of a bond-dipole
static REAL *BondDipoleMagnitude;    // the magnitude of the bond-dipole

// data per framework, used in 'PrecomputeFixedEwaldContributions()'
static REAL *NetChargeFrameworks;                      // the net charge of a framework
static int *NumberOfCoulombicSitesFrameworks;          // the number of charge-sites of a framework
static int *NumberOfBondDipoleSitesFrameworks;         // the number of bond-dipoles of a framework
static COMPLEX *SumFrameworks;                         // the Ewald charge structure factor of a framework
static COMPLEX *SumBondDipoleFrameworks;               // the Ewald bond-dipole structure factor of a framework
static REAL *UChargeChargeFrameworkRigid;              // the charge-charge energy of the fixed framework atoms
static REAL *UChargeBondDipoleFrameworkRigid;          // the charge-bonddipole energy of the fixed framework atoms
static REAL *UBondDipoleBondDipoleFrameworkRigid;      // the bonddipole-bonddipole energy of the fixed framework atoms


/*********************************************************************************************************
 * Name       | InitializeEwald                                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initializes the alpha parameter and number of wavevectors.                               *
 * Parameters | -                                                                                        *
 * Note       | Either alpha, and the amount of k-vectors is explicitly given in the input-file or it is *
 *            | computed from the relative precision.                                                    *
 *********************************************************************************************************/

void InitializeEwald(REAL precision,int Automatic)
{
  int i;
  REAL eps,tol,tol1;
  int ii,jj,kk,nvec,kx,ky,kz;
  REAL ksqr,recip_cutoff;

  if(Automatic)
  {
    // compute the alpha-parameter and max k-vectors from the relative precision
    eps=MIN2(fabs(precision),(REAL)0.5);
    tol=sqrt(fabs(log(eps*CutOffChargeCharge)));

    for(i=0;i<NumberOfSystems;i++)
    {
      Alpha[i]=sqrt(fabs(log(eps*CutOffChargeCharge*tol)))/CutOffChargeCharge;
      tol1=sqrt(-log(eps*CutOffChargeCharge*SQR(2.0*tol*Alpha[i])));

      kvec[i].x=(int)NINT((REAL)0.25+BoxProperties[i].ax*Alpha[i]*tol1/M_PI);
      kvec[i].y=(int)NINT((REAL)0.25+BoxProperties[i].ay*Alpha[i]*tol1/M_PI);
      kvec[i].z=(int)NINT((REAL)0.25+BoxProperties[i].az*Alpha[i]*tol1/M_PI);

      NumberOfKVectors[i]=((kvec[i].x+1)*(2*kvec[i].y+1)*(2*kvec[i].z+1));
      if(ReciprocalCutOffSquared[i]<0.0)
        ReciprocalCutOffSquared[i]=SQR(1.05*MAX3(kvec[i].x,kvec[i].y,kvec[i].z));
    }
  }
  else
  {
    // the alpha-parameter and max k-vectors are set by hand
    for(i=0;i<NumberOfSystems;i++)
    {
      NumberOfKVectors[i]=((kvec[i].x+1)*(2*kvec[i].y+1)*(2*kvec[i].z+1));
      if(ReciprocalCutOffSquared[i]<0.0)
        ReciprocalCutOffSquared[i]=SQR(1.05*MAX3(kvec[i].x,kvec[i].y,kvec[i].z));
    }
  }

  for(i=0;i<NumberOfSystems;i++)
  {
    kx=kvec[i].x;
    ky=kvec[i].y;
    kz=kvec[i].z;
    recip_cutoff=ReciprocalCutOffSquared[i];

    nvec=0;
    for(ii=0;ii<=kx;ii++)
    {
      for(jj=-ky;jj<=ky;jj++)
      {
        for(kk=-kz;kk<=kz;kk++)
        {
          ksqr=SQR(ii)+SQR(jj)+SQR(kk);
          if((ksqr!=0)&&(ksqr<recip_cutoff))
            nvec++;
        }
      }
    }
    NumberOfKVectors[i]=nvec;
  }
}

/*********************************************************************************************************
 * Name       | AllocateEwaldMemory                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Allocates the memory for the Ewald-routines.                                             *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

// Used to allocate the array, and subsequently let another pointer point to the "center" of the array to allow negative indices.
// The pointers here are used for the re-allocation of memory.
static COMPLEX * Eikx0,* Eiky0,* Eikz0;
static COMPLEX * Eikx_bd0,* Eiky_bd0,* Eikz_bd0;

int AllocateEwaldMemory(void)
{
  int i,max;

  UChargeChargeFrameworkRigid=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UChargeBondDipoleFrameworkRigid=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UBondDipoleBondDipoleFrameworkRigid=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHostChargeChargeFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateChargeChargeFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationChargeChargeFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationChargeChargeFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationChargeChargeFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateChargeChargeFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHostChargeBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateChargeBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationChargeBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationChargeBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationChargeBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateChargeBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  UHostHostBondDipoleBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UCationCationBondDipoleBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UAdsorbateCationBondDipoleBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostCationBondDipoleBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  UHostAdsorbateBondDipoleBondDipoleFourierDelta=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  NetChargeSystem=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  NetChargeFramework=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  NetChargeCations=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  NetChargeAdsorbates=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  // find maximum number of frameworks for all systems
  max=0;
  for(i=0;i<NumberOfSystems;i++)
   if(Framework[i].NumberOfFrameworks>max) max=Framework[i].NumberOfFrameworks;

  NumberOfCoulombicSitesFrameworks=(int*)calloc(max,sizeof(int));
  NumberOfBondDipoleSitesFrameworks=(int*)calloc(max,sizeof(int));
  SumFrameworks=(COMPLEX*)calloc(max,sizeof(COMPLEX));
  SumBondDipoleFrameworks=(COMPLEX*)calloc(max,sizeof(COMPLEX));
  NetChargeFrameworks=(REAL*)calloc(max,sizeof(REAL));

  MaxNumberOfWaveVectors=0;
  MaxKvecX=MaxKvecY=MaxKvecZ=0;
  for(i=0;i<NumberOfSystems;i++)
  {
    if(NumberOfKVectors[i]>MaxNumberOfWaveVectors) MaxNumberOfWaveVectors=NumberOfKVectors[i];
    if(kvec[i].x>MaxKvecX) MaxKvecX=kvec[i].x;
    if(kvec[i].y>MaxKvecY) MaxKvecY=kvec[i].y;
    if(kvec[i].z>MaxKvecZ) MaxKvecZ=kvec[i].z;
  }

  // start with 1024 more sites to allow for added particles
  // if these 1024 are used up the memory will be reallocated
  MaxNumberOfCoulombicSites+=1024;
  MaxNumberOfBondDipoleSites+=1024;

  DipoleVector=(VECTOR*)calloc(MAX2(MaxNumberOfBondDipoleSites,MaxNumberOfCoulombicSites),sizeof(VECTOR));      // at least 'MaxNumberOfCoulombicSites' for induced dipoles (polarization)
  BondLength=(REAL*)calloc(MaxNumberOfBondDipoleSites,sizeof(REAL));
  BondDipoleMagnitude=(REAL*)calloc(MaxNumberOfBondDipoleSites,sizeof(REAL));

  Eikr=(COMPLEX*)calloc(MaxNumberOfCoulombicSites,sizeof(COMPLEX));
  Eikr_xy=(COMPLEX*)calloc(MaxNumberOfCoulombicSites,sizeof(COMPLEX));
  Eikx0=(COMPLEX*)calloc((2*MaxKvecX+1)*MaxNumberOfCoulombicSites,sizeof(COMPLEX));
  Eiky0=(COMPLEX*)calloc((2*MaxKvecY+1)*MaxNumberOfCoulombicSites,sizeof(COMPLEX));
  Eikz0=(COMPLEX*)calloc((2*MaxKvecZ+1)*MaxNumberOfCoulombicSites,sizeof(COMPLEX));

  Eikr_bd=(COMPLEX*)calloc(MaxNumberOfBondDipoleSites,sizeof(COMPLEX));
  Eikr_xy_bd=(COMPLEX*)calloc(MaxNumberOfBondDipoleSites,sizeof(COMPLEX));
  Eikx_bd0=(COMPLEX*)calloc((2*MaxKvecX+1)*MaxNumberOfBondDipoleSites,sizeof(COMPLEX));
  Eiky_bd0=(COMPLEX*)calloc((2*MaxKvecY+1)*MaxNumberOfBondDipoleSites,sizeof(COMPLEX));
  Eikz_bd0=(COMPLEX*)calloc((2*MaxKvecZ+1)*MaxNumberOfBondDipoleSites,sizeof(COMPLEX));

  // allow for negative indices
  Eikx=&Eikx0[MaxKvecX*MaxNumberOfCoulombicSites];
  Eiky=&Eiky0[MaxKvecY*MaxNumberOfCoulombicSites];
  Eikz=&Eikz0[MaxKvecZ*MaxNumberOfCoulombicSites];

  Eikx_bd=&Eikx_bd0[MaxKvecX*MaxNumberOfBondDipoleSites];
  Eiky_bd=&Eiky_bd0[MaxKvecY*MaxNumberOfBondDipoleSites];
  Eikz_bd=&Eikz_bd0[MaxKvecZ*MaxNumberOfBondDipoleSites];

  Positions=(VECTOR*)calloc(MaxNumberOfCoulombicSites,sizeof(VECTOR));
  AtomVector=(VECTOR**)calloc(MaxNumberOfCoulombicSites,sizeof(VECTOR*));
  Charge=(REAL*)calloc(MaxNumberOfCoulombicSites,sizeof(REAL));
  Polarization=(REAL*)calloc(MaxNumberOfCoulombicSites,sizeof(REAL));

  BondDipolePositions=(VECTOR*)calloc(MaxNumberOfBondDipoleSites,sizeof(VECTOR));
  BondDipoleForcesA=(VECTOR**)calloc(MaxNumberOfBondDipoleSites,sizeof(VECTOR*));
  BondDipoleForcesB=(VECTOR**)calloc(MaxNumberOfBondDipoleSites,sizeof(VECTOR*));

  KVectors=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  KFactor=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  StoreRigidChargeFramework=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));
  StoreRigidChargeCations=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));
  StoreRigidChargeAdsorbates=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));

  StoreRigidBondDipolesFramework=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));
  StoreRigidBondDipolesCations=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));
  StoreRigidBondDipolesAdsorbates=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));

  StoreTotalChargeFramework=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));
  StoreTotalChargeCations=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));
  StoreTotalChargeAdsorbates=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));

  StoreTotalBondDipolesFramework=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));
  StoreTotalBondDipolesCations=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));
  StoreTotalBondDipolesAdsorbates=(COMPLEX**)calloc(NumberOfSystems,sizeof(COMPLEX*));

  NewTotalChargeFramework=(COMPLEX**)calloc(2,sizeof(COMPLEX*));
  NewTotalChargeAdsorbates=(COMPLEX**)calloc(2,sizeof(COMPLEX*));
  NewTotalChargeCations=(COMPLEX**)calloc(2,sizeof(COMPLEX*));

  NewTotalBondDipolesFramework=(COMPLEX**)calloc(2,sizeof(COMPLEX*));
  NewTotalBondDipolesAdsorbates=(COMPLEX**)calloc(2,sizeof(COMPLEX*));
  NewTotalBondDipolesCations=(COMPLEX**)calloc(2,sizeof(COMPLEX*));

  StoredKVectors=(VECTOR**)calloc(2,sizeof(VECTOR*));
  StoredKFactor=(REAL**)calloc(2,sizeof(REAL*));

  for(i=0;i<NumberOfSystems;i++)
  {
    KVectors[i]=(VECTOR*)calloc(MaxNumberOfWaveVectors,sizeof(VECTOR));
    KFactor[i]=(REAL*)calloc(MaxNumberOfWaveVectors,sizeof(REAL));

    StoreRigidChargeFramework[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
    StoreRigidChargeCations[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
    StoreRigidChargeAdsorbates[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));

    StoreRigidBondDipolesFramework[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
    StoreRigidBondDipolesCations[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
    StoreRigidBondDipolesAdsorbates[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));

    StoreTotalChargeFramework[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
    StoreTotalChargeCations[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
    StoreTotalChargeAdsorbates[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));

    StoreTotalBondDipolesFramework[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
    StoreTotalBondDipolesCations[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
    StoreTotalBondDipolesAdsorbates[i]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
  }

  // temporary wavevectors for system 'A'
  NewTotalChargeFramework[0]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
  NewTotalChargeAdsorbates[0]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
  NewTotalChargeCations[0]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));

  // temporary wavevectors for system 'B'
  NewTotalChargeFramework[1]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
  NewTotalChargeAdsorbates[1]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
  NewTotalChargeCations[1]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));

  // temporary wavevectors for system 'A'
  NewTotalBondDipolesFramework[0]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
  NewTotalBondDipolesAdsorbates[0]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
  NewTotalBondDipolesCations[0]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));

  // temporary wavevectors for system 'B'
  NewTotalBondDipolesFramework[1]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
  NewTotalBondDipolesAdsorbates[1]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));
  NewTotalBondDipolesCations[1]=(COMPLEX*)calloc(MaxNumberOfWaveVectors,sizeof(COMPLEX));

  // temporary wavevectors for system 'A'
  StoredKVectors[0]=(VECTOR*)calloc(MaxNumberOfWaveVectors,sizeof(VECTOR));
  StoredKFactor[0]=(REAL*)calloc(MaxNumberOfWaveVectors,sizeof(REAL));

  // temporary wavevectors for system 'B'
  StoredKVectors[1]=(VECTOR*)calloc(MaxNumberOfWaveVectors,sizeof(VECTOR));
  StoredKFactor[1]=(REAL*)calloc(MaxNumberOfWaveVectors,sizeof(REAL));

  return 0;
}

/*********************************************************************************************************
 * Name       | ReallocateEwaldChargeMemory,ReallocateEwaldBondDipoleMemory                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Allocates more memory for the Ewald-routines.                                            *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

int ReallocateEwaldChargeMemory(void)
{
  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  printf("Realloc MaxNumberOfCoulombicSites: %d\n",MaxNumberOfCoulombicSites);

  Positions=(VECTOR*)realloc(Positions,MaxNumberOfCoulombicSites*sizeof(VECTOR));
  AtomVector=(VECTOR**)realloc(AtomVector,MaxNumberOfCoulombicSites*sizeof(VECTOR*));
  Charge=(REAL*)realloc(Charge,MaxNumberOfCoulombicSites*sizeof(REAL));

  Eikr=(COMPLEX*)realloc(Eikr,MaxNumberOfCoulombicSites*sizeof(COMPLEX));
  Eikr_xy=(COMPLEX*)realloc(Eikr_xy,MaxNumberOfCoulombicSites*sizeof(COMPLEX));
  Eikx0=(COMPLEX*)realloc(Eikx0,(2*MaxKvecX+1)*MaxNumberOfCoulombicSites*sizeof(COMPLEX));
  Eiky0=(COMPLEX*)realloc(Eiky0,(2*MaxKvecY+1)*MaxNumberOfCoulombicSites*sizeof(COMPLEX));
  Eikz0=(COMPLEX*)realloc(Eikz0,(2*MaxKvecZ+1)*MaxNumberOfCoulombicSites*sizeof(COMPLEX));

  // allow for negative indices
  Eikx=&Eikx0[MaxKvecX*MaxNumberOfCoulombicSites];
  Eiky=&Eiky0[MaxKvecY*MaxNumberOfCoulombicSites];
  Eikz=&Eikz0[MaxKvecZ*MaxNumberOfCoulombicSites];

  Positions=(VECTOR*)realloc(Positions,MaxNumberOfCoulombicSites*sizeof(VECTOR));
  AtomVector=(VECTOR**)realloc(AtomVector,MaxNumberOfCoulombicSites*sizeof(VECTOR*));
  Charge=(REAL*)realloc(Charge,MaxNumberOfCoulombicSites*sizeof(REAL));

  return 0;
}

int ReallocateEwaldBondDipoleMemory(void)
{
  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  printf("Realloc MaxNumberOfBondDipoleSites: %d\n",MaxNumberOfBondDipoleSites);

  DipoleVector=(VECTOR*)realloc(DipoleVector,MaxNumberOfBondDipoleSites*sizeof(VECTOR));
  BondLength=(REAL*)realloc(BondLength,MaxNumberOfBondDipoleSites*sizeof(REAL));
  BondDipoleMagnitude=(REAL*)realloc(BondDipoleMagnitude,MaxNumberOfBondDipoleSites*sizeof(REAL));

  Eikr_bd=(COMPLEX*)realloc(Eikr_bd,MaxNumberOfBondDipoleSites*sizeof(COMPLEX));
  Eikr_xy_bd=(COMPLEX*)realloc(Eikr_xy_bd,MaxNumberOfBondDipoleSites*sizeof(COMPLEX));
  Eikx_bd0=(COMPLEX*)realloc(Eikx_bd0,(2*MaxKvecX+1)*MaxNumberOfBondDipoleSites*sizeof(COMPLEX));
  Eiky_bd0=(COMPLEX*)realloc(Eiky_bd0,(2*MaxKvecY+1)*MaxNumberOfBondDipoleSites*sizeof(COMPLEX));
  Eikz_bd0=(COMPLEX*)realloc(Eikz_bd0,(2*MaxKvecZ+1)*MaxNumberOfBondDipoleSites*sizeof(COMPLEX));

  // allow for negative indices
  Eikx_bd=&Eikx_bd0[MaxKvecX*MaxNumberOfBondDipoleSites];
  Eiky_bd=&Eiky_bd0[MaxKvecY*MaxNumberOfBondDipoleSites];
  Eikz_bd=&Eikz_bd0[MaxKvecZ*MaxNumberOfBondDipoleSites];

  BondDipolePositions=(VECTOR*)realloc(BondDipolePositions,MaxNumberOfBondDipoleSites*sizeof(VECTOR));
  BondDipoleForcesA=(VECTOR**)realloc(BondDipoleForcesA,MaxNumberOfBondDipoleSites*sizeof(VECTOR*));
  BondDipoleForcesB=(VECTOR**)realloc(BondDipoleForcesB,MaxNumberOfBondDipoleSites*sizeof(VECTOR*));

  return 0;
}


/*********************************************************************************************************
 * Name       | SetupKVectors                                                                            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Precomputed rk.x, rk.y, and rk.z for each wavevector.                                    *
 * Parameters | -                                                                                        *
 * Note       | Needs to be recomputed for box-changes.                                                  *
 *********************************************************************************************************/

void SetupKVectors(void)
{
  int ii,jj,kk;
  int kx,ky,kz;
  VECTOR rk1,rk2,rk;
  REAL ksqr,rksqr;
  int nvec;
  REAL_MATRIX3x3 inv_box;
  VECTOR *kvecs;
  REAL *kfactor;
  REAL volume,alpha,recip_cutoff;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return;
  if(OmitEwaldFourier) return;

  // setup k-vectors
  nvec=0;
  kx=kvec[CurrentSystem].x;
  ky=kvec[CurrentSystem].y;
  kz=kvec[CurrentSystem].z;
  inv_box=InverseBox[CurrentSystem];
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  volume=Volume[CurrentSystem];
  alpha=Alpha[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  for(ii=0;ii<=kx;ii++)
  {
    rk2.x=2.0*M_PI*ii*inv_box.ax;
    rk2.y=2.0*M_PI*ii*inv_box.bx;
    rk2.z=2.0*M_PI*ii*inv_box.cx;

    for(jj=-ky;jj<=ky;jj++)
    {
      rk1.x=rk2.x+2.0*M_PI*jj*inv_box.ay;
      rk1.y=rk2.y+2.0*M_PI*jj*inv_box.by;
      rk1.z=rk2.z+2.0*M_PI*jj*inv_box.cy;

      for(kk=-kz;kk<=kz;kk++)
      {
        rk.x=rk1.x+2.0*M_PI*kk*inv_box.az;
        rk.y=rk1.y+2.0*M_PI*kk*inv_box.bz;
        rk.z=rk1.z+2.0*M_PI*kk*inv_box.cz;

        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          kvecs[nvec].x=rk.x;
          kvecs[nvec].y=rk.y;
          kvecs[nvec].z=rk.z;

          // In the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used.
          // Therefore the factor is 4*Pi instead of 2*Pi, except for kx=0.
          rksqr=SQR(rk.x)+SQR(rk.y)+SQR(rk.z);
          if(ii==0)
            kfactor[nvec]=COULOMBIC_CONVERSION_FACTOR*(2.0*M_PI/volume)*exp((-0.25/SQR(alpha))*rksqr)/rksqr;
          else
            kfactor[nvec]=COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/volume)*exp((-0.25/SQR(alpha))*rksqr)/rksqr;

          // next wavevector
          nvec++;
        }
      }
    }
  }
}

/*********************************************************************************************************
 * Name       | EwaldEnergyIon                                                                           *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of a single +1 charge in a truly periodic system.              *
 *            | bonddipole-bonddipole interactions for fixed atoms in a truly periodic system.           *                  
 * Parameters | -                                                                                        *
 * Note       | Largely based on the Fiche F.22 of Allen&Tildesly with a few modifications:              *
 *            | 1) exp(-ik.kx)*exp(-ik.ky) computed outside the 'kz' loop                                *
 *            | 2) energy are splitted in framework, adsorbate, cations groups                           *
 *            | 3) net-charge is taken into account                                                      *
 *            | 4) routine computes both charges and bond-dipoles                                        *
 *            | The Ewald sums of a single charge is needed to correct for non-zero net-charges.         *
 *********************************************************************************************************/

int EwaldEnergyIon(void)
{
  int j,ii,jj,kk;
  int nvec;
  REAL alpha,energy,ksqr,recip_cutoff;
  VECTOR pos;
  int index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z;
  REAL *kfactor;
  VECTOR *kvecs;
  COMPLEX sum;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  energy=0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kfactor=KFactor[CurrentSystem];
  kvecs=KVectors[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  Eikx[0].re=1.0; Eikx[0].im=0.0;
  Eiky[0].re=1.0; Eiky[0].im=0.0;
  Eikz[0].re=1.0; Eikz[0].im=0.0;

  pos.x=0.0; pos.y=0.0; pos.z=0.0;

  index_i=MaxNumberOfCoulombicSites;
  Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
  Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
  Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

  index_i=-MaxNumberOfCoulombicSites;
  Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
  Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
  Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
  {
    Eikx[j*MaxNumberOfCoulombicSites].re=Eikx[(j-1)*MaxNumberOfCoulombicSites].re*Eikx[MaxNumberOfCoulombicSites].re-
                                         Eikx[(j-1)*MaxNumberOfCoulombicSites].im*Eikx[MaxNumberOfCoulombicSites].im;
    Eikx[j*MaxNumberOfCoulombicSites].im=Eikx[(j-1)*MaxNumberOfCoulombicSites].im*Eikx[MaxNumberOfCoulombicSites].re+
                                         Eikx[(j-1)*MaxNumberOfCoulombicSites].re*Eikx[MaxNumberOfCoulombicSites].im;
  }

  for(j=2;j<=kmax_y;j++)
  {
    Eiky[j*MaxNumberOfCoulombicSites].re=Eiky[(j-1)*MaxNumberOfCoulombicSites].re*Eiky[MaxNumberOfCoulombicSites].re-
                                         Eiky[(j-1)*MaxNumberOfCoulombicSites].im*Eiky[MaxNumberOfCoulombicSites].im;
    Eiky[j*MaxNumberOfCoulombicSites].im=Eiky[(j-1)*MaxNumberOfCoulombicSites].im*Eiky[MaxNumberOfCoulombicSites].re+
                                         Eiky[(j-1)*MaxNumberOfCoulombicSites].re*Eiky[MaxNumberOfCoulombicSites].im;
    Eiky[-j*MaxNumberOfCoulombicSites].re=Eiky[j*MaxNumberOfCoulombicSites].re;
    Eiky[-j*MaxNumberOfCoulombicSites].im=-Eiky[j*MaxNumberOfCoulombicSites].im;
  }

  for(j=2;j<=kmax_z;j++)
  {
    Eikz[j*MaxNumberOfCoulombicSites].re=Eikz[(j-1)*MaxNumberOfCoulombicSites].re*Eikz[MaxNumberOfCoulombicSites].re-
                                         Eikz[(j-1)*MaxNumberOfCoulombicSites].im*Eikz[MaxNumberOfCoulombicSites].im;
    Eikz[j*MaxNumberOfCoulombicSites].im=Eikz[(j-1)*MaxNumberOfCoulombicSites].im*Eikz[MaxNumberOfCoulombicSites].re+
                                         Eikz[(j-1)*MaxNumberOfCoulombicSites].re*Eikz[MaxNumberOfCoulombicSites].im;
    Eikz[-j*MaxNumberOfCoulombicSites].re=Eikz[j*MaxNumberOfCoulombicSites].re;
    Eikz[-j*MaxNumberOfCoulombicSites].im=-Eikz[j*MaxNumberOfCoulombicSites].im;
  }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
      // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
      index_i=ii*MaxNumberOfCoulombicSites;
      index_j=jj*MaxNumberOfCoulombicSites;
      Eikr_xy[0].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
      Eikr_xy[0].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
          index_k=kk*MaxNumberOfCoulombicSites;
          sum.re=Eikr_xy[0].re*Eikz[index_k].re-Eikr_xy[0].im*Eikz[index_k].im;
          sum.im=Eikr_xy[0].im*Eikz[index_k].re+Eikr_xy[0].re*Eikz[index_k].im;

          energy+=kfactor[nvec]*((SQR(sum.re)+SQR(sum.im)));

          // next wavevector
          nvec++;
        }
      }
    }
  }
  UIon[CurrentSystem]=-(energy-COULOMBIC_CONVERSION_FACTOR*alpha/sqrt(M_PI));
  return 0;
}

/*********************************************************************************************************
 * Name       | PrecomputeFixedEwaldContributions                                                        *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Pre-computes the Fourier part of the charge-charge, charge-bonddipole, and               *
 *            | bonddipole-bonddipole interactions for fixed atoms in a truly periodic system.           *                  
 * Parameters | -                                                                                        *
 * Note       | Largely based on the Fiche F.22 of Allen&Tildesly with a few modifications:              *
 *            | 1) exp(-ik.kx)*exp(-ik.ky) computed outside the 'kz' loop                                *
 *            | 2) energy are splitted in framework, adsorbate, cations groups                           *
 *            | 3) net-charge is taken into account                                                      *
 *            | 4) routine computes both charges and bond-dipoles                                        *
 *            | The routine only loops over fixed atoms and bond-dipole and stores the contributions per *
 *            | wavevector. The routine is called once, from 'input.c'.                                  *
 *            | This routine is used in ewald-energy and force, and Hessian.                             *
 *********************************************************************************************************/

int PrecomputeFixedEwaldContributions(void)
{
  int i,j,f1,ii,jj,kk;
  int nvec,type_mol,type;
  REAL temp,alpha,charge;
  VECTOR pos,rk,posA,posB;
  int nr_molecules,nr_atoms,nr_frameworks;
  int index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_framework,nr_of_coulombic_sites_adsorbate,nr_of_coulombic_sites_cation;
  COMPLEX sum_framework,sum_adsorbate,sum_cation;
  int nr_of_bonddipoles,nr_of_bonddipole_sites,nr_of_bonddipole_sites_framework,nr_of_bonddipole_sites_adsorbate,nr_of_bonddipole_site_cation;
  COMPLEX sum_bonddipole_framework,sum_bonddipole_adsorbate,sum_bonddipole_cation;
  VECTOR dipole;
  VECTOR *kvecs;
  ATOM *atom_pointer;
  PAIR pair;
  REAL *kfactor,ksqr,recip_cutoff;
  int considered_charged;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  UChargeChargeFrameworkRigid[CurrentSystem]=0.0;
  UChargeBondDipoleFrameworkRigid[CurrentSystem]=0.0;
  UBondDipoleBondDipoleFrameworkRigid[CurrentSystem]=0.0;

  // put charge, bond-dipoles, and positions into appropriate arrays
  // ===============================================================

  nr_of_coulombic_sites=0;
  nr_frameworks=Framework[CurrentSystem].NumberOfFrameworks;
  for(f1=0;f1<nr_frameworks;f1++)
  {
    NumberOfCoulombicSitesFrameworks[f1]=0;
    NetChargeFrameworks[f1]=0.0;
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    for(i=0;i<nr_atoms;i++)
    {
      type=Framework[CurrentSystem].Atoms[f1][i].Type;
      charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged&&(Framework[CurrentSystem].Atoms[f1][i].Fixed.x||Framework[CurrentSystem].Atoms[f1][i].Fixed.y||Framework[CurrentSystem].Atoms[f1][i].Fixed.z))
      {
        Charge[nr_of_coulombic_sites]=charge;
        NetChargeFrameworks[f1]+=charge;
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        NumberOfCoulombicSitesFrameworks[f1]++;
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_framework=nr_of_coulombic_sites;
  nr_molecules=NumberOfAdsorbateMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Adsorbates[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        charge=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        if(considered_charged&&(Adsorbates[CurrentSystem][i].Atoms[j].Fixed.x||Adsorbates[CurrentSystem][i].Atoms[j].Fixed.y||Adsorbates[CurrentSystem][i].Atoms[j].Fixed.z))
        {
          Charge[nr_of_coulombic_sites]=charge;
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][i].Atoms[j].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_coulombic_sites_adsorbate=nr_of_coulombic_sites;
  nr_molecules=NumberOfCationMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Cations[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Cations[CurrentSystem][i].Atoms[j].Type;
        charge=Cations[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        if(considered_charged&&(Cations[CurrentSystem][i].Atoms[j].Fixed.x||Cations[CurrentSystem][i].Atoms[j].Fixed.y||Cations[CurrentSystem][i].Atoms[j].Fixed.z))
        {
          Charge[nr_of_coulombic_sites]=Cations[CurrentSystem][i].Atoms[j].Charge;
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][i].Atoms[j].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_coulombic_sites_cation=nr_of_coulombic_sites;

  nr_of_bonddipole_sites=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    NumberOfBondDipoleSitesFrameworks[f1]=0;
    nr_of_bonddipoles=Framework[CurrentSystem].NumberOfBondDipoles[f1];
    for(i=0;i<nr_of_bonddipoles;i++)
    {
      pair=Framework[CurrentSystem].BondDipoles[f1][i];
      atom_pointer=Framework[CurrentSystem].Atoms[f1];

      if((atom_pointer[pair.A].Fixed.x)&&(atom_pointer[pair.B].Fixed.x))
      {
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        dipole=ApplyBoundaryCondition(dipole);
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Framework[CurrentSystem].BondDipoleMagnitude[f1][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        NumberOfBondDipoleSitesFrameworks[f1]++;
        nr_of_bonddipole_sites++;
      }
    }
  }
  nr_of_bonddipole_sites_framework=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Adsorbates[CurrentSystem][i].Atoms;
      if((atom_pointer[pair.A].Fixed.x)&&(atom_pointer[pair.B].Fixed.x))
      {
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        nr_of_bonddipole_sites++;
      }
    }
  }
  nr_of_bonddipole_sites_adsorbate=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Cations[CurrentSystem][i].Atoms;
      if((atom_pointer[pair.A].Fixed.x)&&(atom_pointer[pair.B].Fixed.x))
      {
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        nr_of_bonddipole_sites++;
      }
    }
  }
  nr_of_bonddipole_site_cation=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          rk=kvecs[nvec];

          // loop over all the framework atoms
          for(i=0,f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            SumFrameworks[f1].re=0.0;
            SumFrameworks[f1].im=0.0;
            for(j=0;j<NumberOfCoulombicSitesFrameworks[f1];j++)
            {
              // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
              index_k=kk*MaxNumberOfCoulombicSites+i;
              Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
              Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

              // add contribution to the sum
              temp=Charge[i];
              SumFrameworks[f1].re+=temp*Eikr[i].re;
              SumFrameworks[f1].im+=temp*Eikr[i].im;
              i++;
            }
          }

          // loop over all the adsorbate atoms
          sum_adsorbate.re=0.0;
          sum_adsorbate.im=0.0;
          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_adsorbate.re+=temp*Eikr[i].re;
            sum_adsorbate.im+=temp*Eikr[i].im;
          }

          // loop over all the cation atoms
          sum_cation.re=0.0;
          sum_cation.im=0.0;
          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_cation.re+=temp*Eikr[i].re;
            sum_cation.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_framework.re=0.0;
          sum_bonddipole_framework.im=0.0;
          for(i=0;i<nr_of_bonddipole_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_framework.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_framework.im+=temp*Eikr_bd[i].im;
          }

          // loop over all the framework atoms
          for(i=0,f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            SumBondDipoleFrameworks[f1].re=0.0;
            SumBondDipoleFrameworks[f1].im=0.0;
            for(j=0;j<NumberOfBondDipoleSitesFrameworks[f1];j++)
            {
              // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
              index_k=kk*MaxNumberOfBondDipoleSites+i;
              Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
              Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

              dipole=DipoleVector[i];
              temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
              SumBondDipoleFrameworks[f1].re+=temp*Eikr_bd[i].re;
              SumBondDipoleFrameworks[f1].im+=temp*Eikr_bd[i].im;
              i++;
            }
          }

          sum_bonddipole_adsorbate.re=0.0;
          sum_bonddipole_adsorbate.im=0.0;
          for(;i<nr_of_bonddipole_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_adsorbate.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_adsorbate.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_cation.re=0.0;
          sum_bonddipole_cation.im=0.0;
          for(;i<nr_of_bonddipole_site_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_cation.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_cation.im+=temp*Eikr_bd[i].im;
          }

          // sum the individual frameworks
          sum_framework.re=sum_framework.im=0.0;
          for(i=0;i<Framework[CurrentSystem].NumberOfFrameworks;i++)
          {
            sum_framework.re+=SumFrameworks[i].re;
            sum_framework.im+=SumFrameworks[i].im;
          }

          sum_bonddipole_framework.re=sum_bonddipole_framework.im=0.0;
          for(i=0;i<Framework[CurrentSystem].NumberOfFrameworks;i++)
          {
            sum_bonddipole_framework.re+=SumBondDipoleFrameworks[i].re;
            sum_bonddipole_framework.im+=SumBondDipoleFrameworks[i].im;
          }

          // Store the Fourier-energies of rigid frameworks
          // These will be subtracted in 'EwaldFourierEnergy()' and 'EwaldFourierForce()'
          temp=kfactor[nvec];
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            if(Framework[CurrentSystem].FrameworkModels[f1]!=FLEXIBLE)
            {
              UChargeChargeFrameworkRigid[CurrentSystem]+=temp*(SQR(SumFrameworks[f1].re)+SQR(SumFrameworks[f1].im));
              UChargeBondDipoleFrameworkRigid[CurrentSystem]+=2.0*temp*(SumFrameworks[f1].im*SumBondDipoleFrameworks[f1].re-
                                                                           SumFrameworks[f1].re*SumBondDipoleFrameworks[f1].im);
              UBondDipoleBondDipoleFrameworkRigid[CurrentSystem]+=temp*(SQR(SumBondDipoleFrameworks[f1].re)+SQR(SumBondDipoleFrameworks[f1].im));
            }
          }

          // Store the sums 
          StoreRigidChargeFramework[CurrentSystem][nvec]=sum_framework;
          StoreRigidChargeAdsorbates[CurrentSystem][nvec]=sum_adsorbate;
          StoreRigidChargeCations[CurrentSystem][nvec]=sum_cation;

          StoreRigidBondDipolesFramework[CurrentSystem][nvec]=sum_bonddipole_framework;
          StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec]=sum_bonddipole_adsorbate;
          StoreRigidBondDipolesCations[CurrentSystem][nvec]=sum_bonddipole_cation;

          // next wavevector
          nvec++;
        }
      }
    }
  }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    if(Framework[CurrentSystem].FrameworkModels[f1]!=FLEXIBLE) UChargeChargeFrameworkRigid[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeFrameworks[f1]);

  return 0;
}

/*********************************************************************************************************
 * Name       | PrecomputeTotalEwaldContributions                                                        *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Pre-computes the total Fourier sums of the charge-charge, charge-bonddipole, and         *
 *            | bonddipole-bonddipole interactions for fixed atoms in a truly periodic system.           *                  
 * Parameters | -                                                                                        *
 * Note       | Largely based on the Fiche F.22 of Allen&Tildesly with a few modifications:              *
 *            | 1) exp(-ik.kx)*exp(-ik.ky) computed outside the 'kz' loop                                *
 *            | 2) energy are splitted in framework, adsorbate, cations groups                           *
 *            | 3) net-charge is taken into account                                                      *
 *            | 4) routine computes both charges and bond-dipoles                                        *
 *            | The routine loops over all atoms and bond-dipole and stores the contributions per        *
 *            | wavevector. It also stores the net-charges of the framework, adsorbates and cations.     *
 *            | The routine is called once, from 'input.c'.                                              *
 *            | It is used for Monte-Carlo moves. Here, only the difference in energy is computed and    *
 *            | the Ewald sum runs over the moving atoms. To be able to compute the true difference in   *
 *            | the sums of all the other atoms have to be stored.                                       *
 *********************************************************************************************************/

int PrecomputeTotalEwaldContributions(void)
{
  int i,j,f1,ii,jj,kk;
  int nvec,type_mol,type;
  REAL temp,alpha,ksqr,recip_cutoff,charge;
  VECTOR pos,rk,posA,posB;
  int nr_molecules,nr_atoms,nr_frameworks;
  int index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_framework,nr_of_coulombic_sites_adsorbate,nr_of_coulombic_sites_cation;
  COMPLEX sum_framework,sum_adsorbate,sum_cation;
  int nr_of_bonddipoles,nr_of_bonddipole_sites,nr_of_bonddipole_sites_framework,nr_of_bonddipole_sites_adsorbate,nr_of_bonddipole_site_cation;
  COMPLEX sum_bonddipole_framework,sum_bonddipole_adsorbate,sum_bonddipole_cation;
  VECTOR dipole;
  VECTOR *kvecs;
  ATOM *atom_pointer;
  PAIR pair;
  int considered_charged;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  NetChargeFramework[CurrentSystem]=0.0;
  NetChargeAdsorbates[CurrentSystem]=0.0;
  NetChargeCations[CurrentSystem]=0.0;

  // put charge, bond-dipoles, and positions into appropriate arrays
  // ===============================================================

  nr_of_coulombic_sites=0;
  nr_frameworks=Framework[CurrentSystem].NumberOfFrameworks;
  for(f1=0;f1<nr_frameworks;f1++)
  {
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    for(i=0;i<nr_atoms;i++)
    {
      type=Framework[CurrentSystem].Atoms[f1][i].Type;
      charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged)
      {
        Charge[nr_of_coulombic_sites]=charge;
        NetChargeFramework[CurrentSystem]+=charge;
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_framework=nr_of_coulombic_sites;
  nr_molecules=NumberOfAdsorbateMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Adsorbates[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        charge=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        if(considered_charged)
        {
          Charge[nr_of_coulombic_sites]=charge;
          NetChargeAdsorbates[CurrentSystem]+=charge;
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][i].Atoms[j].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_coulombic_sites_adsorbate=nr_of_coulombic_sites;
  nr_molecules=NumberOfCationMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Cations[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Cations[CurrentSystem][i].Atoms[j].Type;
        charge=Cations[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        if(considered_charged)
        {
          Charge[nr_of_coulombic_sites]=charge;
          NetChargeCations[CurrentSystem]+=charge;
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][i].Atoms[j].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_coulombic_sites_cation=nr_of_coulombic_sites;

  nr_of_bonddipole_sites=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    nr_of_bonddipoles=Framework[CurrentSystem].NumberOfBondDipoles[f1];
    for(i=0;i<nr_of_bonddipoles;i++)
    {
      pair=Framework[CurrentSystem].BondDipoles[f1][i];
      atom_pointer=Framework[CurrentSystem].Atoms[f1];

      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      dipole=ApplyBoundaryCondition(dipole);
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Framework[CurrentSystem].BondDipoleMagnitude[f1][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI; 
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI; 
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_framework=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Adsorbates[CurrentSystem][i].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_adsorbate=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Cations[CurrentSystem][i].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_site_cation=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          rk=kvecs[nvec];

          // loop over all the framework atoms
          sum_framework.re=0.0;
          sum_framework.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_framework.re+=temp*Eikr[i].re;
            sum_framework.im+=temp*Eikr[i].im;
          }

          // loop over all the adsorbate atoms
          sum_adsorbate.re=0.0;
          sum_adsorbate.im=0.0;
          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_adsorbate.re+=temp*Eikr[i].re;
            sum_adsorbate.im+=temp*Eikr[i].im;
          }

          // loop over all the cation atoms
          sum_cation.re=0.0;
          sum_cation.im=0.0;
          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_cation.re+=temp*Eikr[i].re;
            sum_cation.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_framework.re=0.0;
          sum_bonddipole_framework.im=0.0;
          for(i=0;i<nr_of_bonddipole_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_framework.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_framework.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_adsorbate.re=0.0;
          sum_bonddipole_adsorbate.im=0.0;
          for(;i<nr_of_bonddipole_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_adsorbate.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_adsorbate.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_cation.re=0.0;
          sum_bonddipole_cation.im=0.0;
          for(;i<nr_of_bonddipole_site_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_cation.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_cation.im+=temp*Eikr_bd[i].im;
          }

          // Store the sums 
          StoreTotalChargeFramework[CurrentSystem][nvec]=sum_framework;
          StoreTotalChargeAdsorbates[CurrentSystem][nvec]=sum_adsorbate;
          StoreTotalChargeCations[CurrentSystem][nvec]=sum_cation;

          StoreTotalBondDipolesFramework[CurrentSystem][nvec]=sum_bonddipole_framework;
          StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec]=sum_bonddipole_adsorbate;
          StoreTotalBondDipolesCations[CurrentSystem][nvec]=sum_bonddipole_cation;

          // next wavevector
          nvec++;
        }
      }
    }
  }

  NetChargeSystem[CurrentSystem]=NetChargeFramework[CurrentSystem]+NetChargeAdsorbates[CurrentSystem]+NetChargeCations[CurrentSystem];

  return 0;
}


/*********************************************************************************************************
 * Name       | EwaldFourierEnergy                                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of the charge-charge, charge-bonddipole, and                   *
 *            | bonddipole-bonddipole interactions in a truly periodic system.                           *                  
 * Parameters | -                                                                                        *
 * Note       | Largely based on the Fiche F.22 of Allen&Tildesly with a few modifications:              *
 *            | 1) exp(-ik.kx)*exp(-ik.ky) computed outside the 'kz' loop                                *
 *            | 2) energy are splitted in framework, adsorbate, cations groups                           *
 *            | 3) net-charge is taken into account                                                      *
 *            | 4) routine computes both charges and bond-dipoles                                        *
 *            | The total structure factors are stored.                                                  *
 *********************************************************************************************************/

int EwaldFourierEnergy(void)
{
  int i,j,m,f1,ii,jj,kk,A,B;
  int nvec,type_mol,type,typeA,typeB;
  REAL temp,alpha,r,rr,chargeA,chargeB,charge;
  VECTOR pos,rk,dr,posA,posB,posA1,posA2,posB1,posB2;
  int nr_molecules,nr_atoms,nr_frameworks,nr_of_excluded_pairs;
  int index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_framework,nr_of_coulombic_sites_adsorbate,nr_of_coulombic_sites_cation;
  COMPLEX sum_framework,sum_adsorbate,sum_cation;
  REAL energy_framework_framework,energy_adsorbate_adsorbate,energy_cation_cation;
  REAL energy_framework_adsorbate,energy_framework_cation,energy_adsorbate_cation;
  REAL energy_framework_self,energy_adsorbate_self,energy_cation_self;
  int nr_of_bonddipoles,nr_of_bonddipole_sites,nr_of_bonddipole_sites_framework,nr_of_bonddipole_sites_adsorbate,nr_of_bonddipole_site_cation;
  COMPLEX sum_bonddipole_framework,sum_bonddipole_adsorbate,sum_bonddipole_cation;
  VECTOR dipole,dipoleA,dipoleB;
  REAL energy_framework_framework_c_bd,energy_framework_adsorbate_c_bd,energy_framework_cation_c_bd;
  REAL energy_adsorbate_adsorbate_c_bd,energy_cation_cation_c_bd,energy_adsorbate_cation_c_bd;
  REAL energy_framework_framework_bd,energy_framework_adsorbate_bd,energy_framework_cation_bd;
  REAL energy_adsorbate_adsorbate_bd,energy_cation_cation_bd,energy_adsorbate_cation_bd;
  REAL energy_framework_self_bd,energy_adsorbate_self_bd,energy_cation_self_bd;
  REAL energy_framework_excluded,energy_adsorbate_excluded,energy_cation_excluded;
  REAL energy_framework_excluded_c_bd,energy_adsorbate_excluded_c_bd,energy_cation_excluded_c_bd;
  REAL energy_framework_excluded_bd,energy_adsorbate_excluded_bd,energy_cation_excluded_bd;
  REAL Bt0,Bt1,Bt2,cosA,cosB,cosAB,scaling,scalingA,scalingB;
  ATOM *atom_pointer;
  VECTOR *kvecs;
  REAL *kfactor,ksqr,recip_cutoff;
  PAIR pair;
  int considered_charged;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  // initialize energies
  UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UCationCationChargeChargeFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UHostHostChargeChargeFourier[CurrentSystem]=0.0;

  UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostHostChargeBondDipoleFourier[CurrentSystem]=0.0;

  UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=0.0;

  // initialize energies
  energy_framework_framework=energy_adsorbate_adsorbate=energy_cation_cation=0.0;
  energy_framework_adsorbate=energy_framework_cation=energy_adsorbate_cation=0.0;
  energy_framework_self=energy_adsorbate_self=energy_cation_self=0.0;

  energy_framework_framework_c_bd=energy_adsorbate_adsorbate_c_bd=energy_cation_cation_c_bd=0.0;
  energy_framework_adsorbate_c_bd=energy_framework_cation_c_bd=energy_adsorbate_cation_c_bd=0.0;

  energy_framework_framework_bd=energy_adsorbate_adsorbate_bd=energy_cation_cation_bd=0.0;
  energy_framework_adsorbate_bd=energy_framework_cation_bd=energy_adsorbate_cation_bd=0.0;
  energy_framework_self_bd=energy_adsorbate_self_bd=energy_cation_self_bd=0.0;

  energy_framework_excluded=energy_adsorbate_excluded=energy_cation_excluded=0.0;
  energy_framework_excluded_c_bd=energy_adsorbate_excluded_c_bd=energy_cation_excluded_c_bd=0.0;
  energy_framework_excluded_bd=energy_adsorbate_excluded_bd=energy_cation_excluded_bd=0.0;

  // put charge, bond-dipoles, and positions into appropriate arrays
  // ===============================================================

  nr_of_coulombic_sites=0;
  nr_frameworks=Framework[CurrentSystem].NumberOfFrameworks;
  for(f1=0;f1<nr_frameworks;f1++)
  {
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    for(i=0;i<nr_atoms;i++)
    {
      type=Framework[CurrentSystem].Atoms[f1][i].Type;
      charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged)
      {
        Charge[nr_of_coulombic_sites]=charge;
        if(!(Framework[CurrentSystem].Atoms[f1][i].Fixed.x&&Framework[CurrentSystem].Atoms[f1][i].Fixed.y&&Framework[CurrentSystem].Atoms[f1][i].Fixed.z))
        {
          energy_framework_self+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_coulombic_sites_framework=nr_of_coulombic_sites;
  nr_molecules=NumberOfAdsorbateMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Adsorbates[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        scaling=Adsorbates[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
        charge=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        charge*=scaling;
        if(considered_charged)
        {
          Charge[nr_of_coulombic_sites]=charge;
          energy_adsorbate_self+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
          if(!(Adsorbates[CurrentSystem][i].Atoms[j].Fixed.x&&Adsorbates[CurrentSystem][i].Atoms[j].Fixed.y&&Adsorbates[CurrentSystem][i].Atoms[j].Fixed.z))
          {
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][i].Atoms[j].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }
  nr_of_coulombic_sites_adsorbate=nr_of_coulombic_sites;
  nr_molecules=NumberOfCationMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Cations[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Cations[CurrentSystem][i].Atoms[j].Type;
        scaling=Cations[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
        charge=Cations[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        charge*=scaling;
        if(considered_charged)
        {
          Charge[nr_of_coulombic_sites]=charge;
          energy_cation_self+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
          if(!(Cations[CurrentSystem][i].Atoms[j].Fixed.x&&Cations[CurrentSystem][i].Atoms[j].Fixed.y&&Cations[CurrentSystem][i].Atoms[j].Fixed.z))
          {
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][i].Atoms[j].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }
  nr_of_coulombic_sites_cation=nr_of_coulombic_sites;

  nr_of_bonddipole_sites=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    nr_of_bonddipoles=Framework[CurrentSystem].NumberOfBondDipoles[f1];
    for(i=0;i<nr_of_bonddipoles;i++)
    {
      pair=Framework[CurrentSystem].BondDipoles[f1][i];
      atom_pointer=Framework[CurrentSystem].Atoms[f1];
      if((!atom_pointer[pair.A].Fixed.x)||(!atom_pointer[pair.B].Fixed.x))
      {
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        dipole=ApplyBoundaryCondition(dipole);
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Framework[CurrentSystem].BondDipoleMagnitude[f1][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;

        energy_framework_self_bd+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        nr_of_bonddipole_sites++;
      }
    }
  }
  nr_of_bonddipole_sites_framework=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Adsorbates[CurrentSystem][i].Atoms;
      if((!atom_pointer[pair.A].Fixed.x)||(!atom_pointer[pair.B].Fixed.x))
      {
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        energy_adsorbate_self_bd+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        nr_of_bonddipole_sites++;
      }
    }
  }
  nr_of_bonddipole_sites_adsorbate=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Cations[CurrentSystem][i].Atoms;
      if((!atom_pointer[pair.A].Fixed.x)||(!atom_pointer[pair.B].Fixed.x))
      {
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        energy_cation_self_bd+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        nr_of_bonddipole_sites++;
      }
    }
  }
  nr_of_bonddipole_site_cation=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          rk=kvecs[nvec];

          // loop over all the framework atoms
          sum_framework.re=0.0;
          sum_framework.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_framework.re+=temp*Eikr[i].re;
            sum_framework.im+=temp*Eikr[i].im;
          }

          // loop over all the adsorbate atoms
          sum_adsorbate.re=0.0;
          sum_adsorbate.im=0.0;
          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_adsorbate.re+=temp*Eikr[i].re;
            sum_adsorbate.im+=temp*Eikr[i].im;
          }

          // loop over all the cation atoms
          sum_cation.re=0.0;
          sum_cation.im=0.0;
          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_cation.re+=temp*Eikr[i].re;
            sum_cation.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_framework.re=0.0;
          sum_bonddipole_framework.im=0.0;
          for(i=0;i<nr_of_bonddipole_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_framework.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_framework.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_adsorbate.re=0.0;
          sum_bonddipole_adsorbate.im=0.0;
          for(;i<nr_of_bonddipole_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_adsorbate.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_adsorbate.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_cation.re=0.0;
          sum_bonddipole_cation.im=0.0;
          for(;i<nr_of_bonddipole_site_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_cation.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_cation.im+=temp*Eikr_bd[i].im;
          }

          // add the pre-computed contributions of fixed atoms
          sum_framework.re+=StoreRigidChargeFramework[CurrentSystem][nvec].re;
          sum_framework.im+=StoreRigidChargeFramework[CurrentSystem][nvec].im;
          sum_adsorbate.re+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].re;
          sum_adsorbate.im+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].im;
          sum_cation.re+=StoreRigidChargeCations[CurrentSystem][nvec].re;
          sum_cation.im+=StoreRigidChargeCations[CurrentSystem][nvec].im;

          sum_bonddipole_framework.re+=StoreRigidBondDipolesFramework[CurrentSystem][nvec].re;
          sum_bonddipole_framework.im+=StoreRigidBondDipolesFramework[CurrentSystem][nvec].im;
          sum_bonddipole_adsorbate.re+=StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec].re;
          sum_bonddipole_adsorbate.im+=StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec].im;
          sum_bonddipole_cation.re+=StoreRigidBondDipolesCations[CurrentSystem][nvec].re;
          sum_bonddipole_cation.im+=StoreRigidBondDipolesCations[CurrentSystem][nvec].im;

          // Store the total structure factors
          StoreTotalChargeFramework[CurrentSystem][nvec]=sum_framework;
          StoreTotalChargeAdsorbates[CurrentSystem][nvec]=sum_adsorbate;
          StoreTotalChargeCations[CurrentSystem][nvec]=sum_cation;

          StoreTotalBondDipolesFramework[CurrentSystem][nvec]=sum_bonddipole_framework;
          StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec]=sum_bonddipole_adsorbate;
          StoreTotalBondDipolesCations[CurrentSystem][nvec]=sum_bonddipole_cation;

          // precomputed wavevector dependent pre-factor
          temp=kfactor[nvec];

          // charge-charge energies
          energy_framework_framework+=temp*(SQR(sum_framework.re)+SQR(sum_framework.im));
          energy_adsorbate_adsorbate+=temp*(SQR(sum_adsorbate.re)+SQR(sum_adsorbate.im));
          energy_cation_cation+=temp*(SQR(sum_cation.re)+SQR(sum_cation.im));
          energy_framework_adsorbate+=temp*(sum_framework.re*sum_adsorbate.re+sum_framework.im*sum_adsorbate.im);
          energy_framework_cation+=temp*(sum_framework.re*sum_cation.re+sum_framework.im*sum_cation.im);
          energy_adsorbate_cation+=temp*(sum_adsorbate.re*sum_cation.re+sum_adsorbate.im*sum_cation.im);

          // charge-bond-dipole energies
          energy_framework_framework_c_bd+=2.0*temp*(sum_framework.im*sum_bonddipole_framework.re-sum_framework.re*sum_bonddipole_framework.im);
          energy_adsorbate_adsorbate_c_bd+=2.0*temp*(sum_adsorbate.im*sum_bonddipole_adsorbate.re-sum_adsorbate.re*sum_bonddipole_adsorbate.im);
          energy_cation_cation_c_bd+=2.0*temp*(sum_cation.im*sum_bonddipole_cation.re-sum_cation.re*sum_bonddipole_cation.im);
          energy_framework_adsorbate_c_bd+=temp*(sum_adsorbate.im*sum_bonddipole_framework.re+sum_framework.im*sum_bonddipole_adsorbate.re-
                                                (sum_adsorbate.re*sum_bonddipole_framework.im+sum_framework.re*sum_bonddipole_adsorbate.im));
          energy_framework_cation_c_bd+=temp*(sum_cation.im*sum_bonddipole_framework.re+sum_framework.im*sum_bonddipole_cation.re-
                                             (sum_cation.re*sum_bonddipole_framework.im+sum_framework.re*sum_bonddipole_cation.im));
          energy_adsorbate_cation_c_bd+=temp*(sum_cation.im*sum_bonddipole_adsorbate.re+sum_adsorbate.im*sum_bonddipole_cation.re-
                                             (sum_cation.re*sum_bonddipole_adsorbate.im+sum_adsorbate.re*sum_bonddipole_cation.im));

          // bonddipole-bond-dipole energies
          energy_framework_framework_bd+=temp*(SQR(sum_bonddipole_framework.re)+SQR(sum_bonddipole_framework.im));
          energy_adsorbate_adsorbate_bd+=temp*(SQR(sum_bonddipole_adsorbate.re)+SQR(sum_bonddipole_adsorbate.im));
          energy_cation_cation_bd+=temp*(SQR(sum_bonddipole_cation.re)+SQR(sum_bonddipole_cation.im));
          energy_framework_adsorbate_bd+=temp*(sum_bonddipole_framework.re*sum_bonddipole_adsorbate.re+sum_bonddipole_framework.im*sum_bonddipole_adsorbate.im);
          energy_framework_cation_bd+=temp*(sum_bonddipole_framework.re*sum_bonddipole_cation.re+sum_bonddipole_framework.im*sum_bonddipole_cation.im);
          energy_adsorbate_cation_bd+=temp*(sum_bonddipole_adsorbate.re*sum_bonddipole_cation.re+sum_bonddipole_adsorbate.im*sum_bonddipole_cation.im);

          // next wavevector
          nvec++;
        }
      }
    }
  }

  // Exclusion pairs for the framework
  // =================================
  // The list of exclusion pairs is pre-compute at the start of the simulation based on e.g. 1-2, 1-3 exclusion definitions.
  // Note: there is no cutoff used here

  energy_framework_excluded=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        atom_pointer=Framework[CurrentSystem].Atoms[f1];
        pair=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        chargeA=atom_pointer[pair.A].Charge;
        chargeB=atom_pointer[pair.B].Charge;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        if(r>1e-4)
          energy_framework_excluded+=erf(alpha*r)*chargeA*chargeB/r;
        else // if r->0 compute limiting value to avoid divergence when shell overlaps with core in core-shell models
          energy_framework_excluded+=chargeA*chargeB*alpha*(2.0/sqrt(M_PI));
      }
    }
  }
  energy_framework_excluded*=COULOMBIC_CONVERSION_FACTOR;

  energy_framework_excluded_c_bd=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeBondDipole[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        // the index of the bonddipole is the second index 'B'
        A=Framework[CurrentSystem].ExcludedIntraChargeBondDipole[f1][i].A;
        B=Framework[CurrentSystem].ExcludedIntraChargeBondDipole[f1][i].B;

        atom_pointer=Framework[CurrentSystem].Atoms[f1];

        typeA=atom_pointer[A].Type;
        chargeA=atom_pointer[A].Charge;
        posA=atom_pointer[A].Position;

        pair=Framework[CurrentSystem].BondDipoles[f1][B];
        posB1=atom_pointer[pair.A].Position;
        posB2=atom_pointer[pair.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=Framework[CurrentSystem].BondDipoleMagnitude[f1][B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posB.x-posA.x;
        dr.y=posB.y-posA.y;
        dr.z=posB.z-posA.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=erf(alpha*r);
        Bt0=-temp/r;
        Bt1=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr-temp/(rr*r);

        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_framework_excluded_c_bd+=Bt1*chargeA*cosB;
      }
    }
  }
  energy_framework_excluded_c_bd*=COULOMBIC_CONVERSION_FACTOR;

  energy_framework_excluded_bd=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraBondDipoleBondDipole[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Framework[CurrentSystem].ExcludedIntraBondDipoleBondDipole[f1][i].A;
        B=Framework[CurrentSystem].ExcludedIntraBondDipoleBondDipole[f1][i].B;

        atom_pointer=Framework[CurrentSystem].Atoms[f1];

        pair=Framework[CurrentSystem].BondDipoles[f1][A];
        posA1=atom_pointer[pair.A].Position;
        posA2=atom_pointer[pair.B].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        dipoleA=ApplyBoundaryCondition(dipoleA);
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        temp=Framework[CurrentSystem].BondDipoleMagnitude[f1][A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        pair=Framework[CurrentSystem].BondDipoles[f1][B];
        posB1=atom_pointer[pair.A].Position;
        posB2=atom_pointer[pair.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        temp=Framework[CurrentSystem].BondDipoleMagnitude[f1][B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        Bt0=-erf(alpha*r)/r;
        Bt1=2.0*alpha*exp(-SQR(-alpha*r))/(sqrt(M_PI)*rr)+
            -erf(alpha*r)/(rr*r);
        Bt2=6.0*alpha*exp(-SQR(alpha*r))/(sqrt(M_PI)*rr*rr)+
            4.0*alpha*SQR(alpha)*exp(-SQR(alpha*r))/(sqrt(M_PI)*rr)+
            -3.0*erf(alpha*r)/(rr*rr*r);

        cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
        cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_framework_excluded_bd-=(Bt1*cosAB-Bt2*cosA*cosB);
      }
    }
  }
  energy_framework_excluded_bd*=COULOMBIC_CONVERSION_FACTOR;

  // Exclusion pairs for the adsorbates
  // ==================================

  energy_adsorbate_excluded=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    atom_pointer=Adsorbates[CurrentSystem][m].Atoms;
    type=Adsorbates[CurrentSystem][m].Type;
    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      pair=Components[type].ExcludedIntraChargeCharge[i];
      typeA=atom_pointer[pair.A].Type;
      typeB=atom_pointer[pair.B].Type;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      scalingA=atom_pointer[pair.A].CFChargeScalingParameter;
      scalingB=atom_pointer[pair.B].CFChargeScalingParameter;
      chargeA=scalingA*atom_pointer[pair.A].Charge;
      chargeB=scalingB*atom_pointer[pair.B].Charge;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      energy_adsorbate_excluded+=erf(Alpha[CurrentSystem]*r)*chargeA*chargeB/r;
    }
  }
  energy_adsorbate_excluded*=COULOMBIC_CONVERSION_FACTOR;

  energy_adsorbate_excluded_c_bd=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    type=Adsorbates[CurrentSystem][m].Type;
    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraChargeBondDipole[i].A;
      B=Components[type].ExcludedIntraChargeBondDipole[i].B;

      atom_pointer=Adsorbates[CurrentSystem][m].Atoms;

      chargeA=Components[type].Charge[A];
      posA=atom_pointer[A].Position;

      pair=Components[type].BondDipoles[B];
      posB1=atom_pointer[pair.A].Position;
      posB2=atom_pointer[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(Alpha[CurrentSystem]*r)/r;
      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)-erf(Alpha[CurrentSystem]*r)/(rr*r);

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_adsorbate_excluded_c_bd+=Bt1*chargeA*cosB;
    }
  }
  energy_adsorbate_excluded_c_bd*=COULOMBIC_CONVERSION_FACTOR;

  energy_adsorbate_excluded_bd=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    type=Adsorbates[CurrentSystem][m].Type;
    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[type].ExcludedIntraBondDipoleBondDipole[i].B;

      atom_pointer=Adsorbates[CurrentSystem][m].Atoms;

      pair=Components[type].BondDipoles[A];
      posA1=atom_pointer[pair.A].Position;
      posA2=atom_pointer[pair.B].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=Components[type].BondDipoleMagnitude[A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      pair=Components[type].BondDipoles[B];
      posB1=atom_pointer[pair.A].Position;
      posB2=atom_pointer[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(Alpha[CurrentSystem]*r)/r;
      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -erf(Alpha[CurrentSystem]*r)/(rr*r);
      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -3.0*erf(Alpha[CurrentSystem]*r)/(rr*rr*r);

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_adsorbate_excluded_bd-=Bt1*cosAB-Bt2*cosA*cosB;
    }
  }
  energy_adsorbate_excluded_bd*=COULOMBIC_CONVERSION_FACTOR;

  // Exclusion pairs for the cations
  // ===============================

  energy_cation_excluded=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    atom_pointer=Cations[CurrentSystem][m].Atoms;
    type=Cations[CurrentSystem][m].Type;
    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      pair=Components[type].ExcludedIntraChargeCharge[i];
      typeA=atom_pointer[pair.A].Type;
      typeB=atom_pointer[pair.B].Type;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      scalingA=atom_pointer[pair.A].CFChargeScalingParameter;
      scalingB=atom_pointer[pair.B].CFChargeScalingParameter;
      chargeA=scalingA*atom_pointer[pair.A].Charge;
      chargeB=scalingB*atom_pointer[pair.B].Charge;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      energy_cation_excluded+=erf(Alpha[CurrentSystem]*r)*chargeA*chargeB/r;
    }
  }
  energy_cation_excluded*=COULOMBIC_CONVERSION_FACTOR;

  energy_cation_excluded_c_bd=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    type=Cations[CurrentSystem][m].Type;
    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraChargeBondDipole[i].A;
      B=Components[type].ExcludedIntraChargeBondDipole[i].B;

      atom_pointer=Cations[CurrentSystem][m].Atoms;

      chargeA=Components[type].Charge[A];
      posA=atom_pointer[A].Position;

      pair=Components[type].BondDipoles[B];
      posB1=atom_pointer[pair.A].Position;
      posB2=atom_pointer[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(Alpha[CurrentSystem]*r)/r;
      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)-erf(Alpha[CurrentSystem]*r)/(rr*r);

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_cation_excluded_c_bd+=Bt1*chargeA*cosB;
    }
  }
  energy_cation_excluded_c_bd*=COULOMBIC_CONVERSION_FACTOR;

  energy_cation_excluded_bd=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    type=Cations[CurrentSystem][m].Type;
    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[type].ExcludedIntraBondDipoleBondDipole[i].B;

      atom_pointer=Cations[CurrentSystem][m].Atoms;

      pair=Components[type].BondDipoles[A];
      posA1=atom_pointer[pair.A].Position;
      posA2=atom_pointer[pair.B].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=Components[type].BondDipoleMagnitude[A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      pair=Components[type].BondDipoles[B];
      posB1=atom_pointer[pair.A].Position;
      posB2=atom_pointer[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(Alpha[CurrentSystem]*r)/r;
      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -erf(Alpha[CurrentSystem]*r)/(rr*r);
      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -3.0*erf(Alpha[CurrentSystem]*r)/(rr*rr*r);

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_cation_excluded_bd-=Bt1*cosAB-Bt2*cosA*cosB;
    }
  }
  energy_cation_excluded_bd*=COULOMBIC_CONVERSION_FACTOR;


  // compute energies for host-host, host-adsorbate, host-cation, adsorbate-adsorbate, cation-cation, and adsorbate-cation
  // =====================================================================================================================
  // Note: total energy is Fourier minus self minus exclusion energy
  // Note: the net-charge is taken into account

  // Host-host energy (corrected for net-charge)
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    UHostHostChargeChargeFourier[CurrentSystem]=energy_framework_framework-energy_framework_self-energy_framework_excluded;
    UHostHostChargeBondDipoleFourier[CurrentSystem]=energy_framework_framework_c_bd-energy_framework_excluded_c_bd;
    UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=energy_framework_framework_bd-energy_framework_self_bd-energy_framework_excluded_bd;
    UHostHostChargeChargeFourier[CurrentSystem]-=UChargeChargeFrameworkRigid[CurrentSystem];
    UHostHostChargeBondDipoleFourier[CurrentSystem]-=UChargeBondDipoleFrameworkRigid[CurrentSystem];
    UHostHostBondDipoleBondDipoleFourier[CurrentSystem]-=UBondDipoleBondDipoleFrameworkRigid[CurrentSystem];
    UHostHostChargeChargeFourier[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeFramework[CurrentSystem]);
  }

  // Host-adsorbate energy (corrected for net-charge)
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=2.0*energy_framework_adsorbate;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=2.0*energy_framework_adsorbate_c_bd;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=2.0*energy_framework_adsorbate_bd;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*NetChargeAdsorbates[CurrentSystem];

  // Host-cation energy (corrected for net-charge)
  UHostCationChargeChargeFourier[CurrentSystem]=2.0*energy_framework_cation;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=2.0*energy_framework_cation_c_bd;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=2.0*energy_framework_cation_bd;
  UHostCationChargeChargeFourier[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*NetChargeCations[CurrentSystem];

  // Adsorbate-adsorbate energy (corrected for net-charge)
  if(!OmitAdsorbateAdsorbateCoulombInteractions)
  {
    UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=energy_adsorbate_adsorbate-energy_adsorbate_self-energy_adsorbate_excluded;
    UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=energy_adsorbate_adsorbate_c_bd-energy_adsorbate_excluded_c_bd;
    UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=energy_adsorbate_adsorbate_bd-energy_adsorbate_self_bd-energy_adsorbate_excluded_bd;
    UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeAdsorbates[CurrentSystem]);
  }

  // Cation-cation energy (corrected for net-charge)
  if(!OmitCationCationCoulombInteractions)
  {
    UCationCationChargeChargeFourier[CurrentSystem]=energy_cation_cation-energy_cation_self-energy_cation_excluded;
    UCationCationChargeBondDipoleFourier[CurrentSystem]=energy_cation_cation_c_bd-energy_cation_excluded_c_bd;
    UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=energy_cation_cation_bd-energy_cation_self_bd-energy_cation_excluded_bd;
    UCationCationChargeChargeFourier[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeCations[CurrentSystem]);
  }

  // Adsorbate-cation energy (corrected for net-charge)
  if(!OmitAdsorbateCationCoulombInteractions)
  {
    UAdsorbateCationChargeChargeFourier[CurrentSystem]=2.0*energy_adsorbate_cation;
    UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=2.0*energy_adsorbate_cation_c_bd;
    UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=2.0*energy_adsorbate_cation_bd;
    UAdsorbateCationChargeChargeFourier[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeAdsorbates[CurrentSystem]*NetChargeCations[CurrentSystem];
  }

  return 0;
}

/*********************************************************************************************************
 * Name       | EwaldFourierForce                                                                        *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of the charge-charge, charge-bonddipole, and                   *
 *            | bonddipole-bonddipole energies, forces and stress n a truly periodic system.             *                  
 * Parameters | -                                                                                        *
 * Note       | Largely based on the Fiche F.22 of Allen&Tildesly with a few modifications:              *
 *            | 1) exp(-ik.kx)*exp(-ik.ky) computed outside the 'kz' loop                                *
 *            | 2) energy are splitted in framework, adsorbate, cations groups                           *
 *            | 3) net-charge is taken into account                                                      *
 *            | 4) routine computes both charges and bond-dipoles                                        *
 *            | The total structure factors are not stored.                                              *
 *********************************************************************************************************/

int EwaldFourierForce(void)
{
  int i,j,m,f1,ii,jj,kk,A,B;
  int nvec,type_mol,type,typeA,typeB;
  REAL temp,alpha,r,rr,chargeA,chargeB,charge;
  VECTOR pos,rk,dr,posA,posB,posA1,posA2,posB1,posB2;
  int nr_molecules,nr_atoms,nr_frameworks,nr_of_excluded_pairs;
  int index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_framework,nr_of_coulombic_sites_adsorbate,nr_of_coulombic_sites_cation;
  COMPLEX sum,sum_framework,sum_adsorbate,sum_cation,temp_sum,temp_sum_bonddipole;
  REAL energy_framework_framework,energy_adsorbate_adsorbate,energy_cation_cation;
  REAL energy_framework_adsorbate,energy_framework_cation,energy_adsorbate_cation;
  REAL energy_framework_self,energy_adsorbate_self,energy_cation_self;
  int nr_of_bonddipoles,nr_of_bonddipole_sites,nr_of_bonddipole_sites_framework,nr_of_bonddipole_sites_adsorbate,nr_of_bonddipole_site_cation;
  COMPLEX sum_bonddipole,sum_bonddipole_framework,sum_bonddipole_adsorbate,sum_bonddipole_cation;
  VECTOR dipole,dipoleA,dipoleB;
  REAL energy_framework_framework_c_bd,energy_framework_adsorbate_c_bd,energy_framework_cation_c_bd;
  REAL energy_adsorbate_adsorbate_c_bd,energy_cation_cation_c_bd,energy_adsorbate_cation_c_bd;
  REAL energy_framework_framework_bd,energy_framework_adsorbate_bd,energy_framework_cation_bd;
  REAL energy_adsorbate_adsorbate_bd,energy_cation_cation_bd,energy_adsorbate_cation_bd;
  REAL energy_framework_self_bd,energy_adsorbate_self_bd,energy_cation_self_bd;
  REAL energy_framework_excluded,energy_adsorbate_excluded,energy_cation_excluded;
  REAL energy_framework_excluded_c_bd,energy_adsorbate_excluded_c_bd,energy_cation_excluded_c_bd;
  REAL energy_framework_excluded_bd,energy_adsorbate_excluded_bd,energy_cation_excluded_bd;
  REAL Bt0,Bt1,Bt2,Bt3,cosA,cosB,cosAB,fac,fac1,fac2,fac3,fac4,fac5,fac6,current_energy;
  REAL DipoleMagnitudeA,length,dot_product,lengthA,lengthB,dipole_magnitudeA,dipole_magnitudeB;
  VECTOR fa1,fa2,fb1,fb2,term,termA,termB;
  REAL scaling,scalingA,scalingB;
  ATOM *atom_pointer;
  REAL *kfactor,ksqr,recip_cutoff;
  PAIR pair,pairA,pairB;
  VECTOR *kvecs;
  REAL_MATRIX3x3 stress,v;
  int considered_charged;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  stress.ax=stress.bx=stress.cx=0.0;
  stress.ay=stress.by=stress.cy=0.0;
  stress.az=stress.bz=stress.cz=0.0;

  // initialize energies
  UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UCationCationChargeChargeFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UHostHostChargeChargeFourier[CurrentSystem]=0.0;

  UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostHostChargeBondDipoleFourier[CurrentSystem]=0.0;

  UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=0.0;

  // initialize energies
  energy_framework_framework=energy_adsorbate_adsorbate=energy_cation_cation=0.0;
  energy_framework_adsorbate=energy_framework_cation=energy_adsorbate_cation=0.0;
  energy_framework_self=energy_adsorbate_self=energy_cation_self=0.0;

  energy_framework_framework_c_bd=energy_adsorbate_adsorbate_c_bd=energy_cation_cation_c_bd=0.0;
  energy_framework_adsorbate_c_bd=energy_framework_cation_c_bd=energy_adsorbate_cation_c_bd=0.0;

  energy_framework_framework_bd=energy_adsorbate_adsorbate_bd=energy_cation_cation_bd=0.0;
  energy_framework_adsorbate_bd=energy_framework_cation_bd=energy_adsorbate_cation_bd=0.0;
  energy_framework_self_bd=energy_adsorbate_self_bd=energy_cation_self_bd=0.0;

  energy_framework_excluded=energy_adsorbate_excluded=energy_cation_excluded=0.0;
  energy_framework_excluded_c_bd=energy_adsorbate_excluded_c_bd=energy_cation_excluded_c_bd=0.0;
  energy_framework_excluded_bd=energy_adsorbate_excluded_bd=energy_cation_excluded_bd=0.0;

  // put charge, bond-dipoles, and positions into appropriate arrays
  // ===============================================================

  nr_of_coulombic_sites=0;
  nr_frameworks=Framework[CurrentSystem].NumberOfFrameworks;
  for(f1=0;f1<nr_frameworks;f1++)
  {
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    for(i=0;i<nr_atoms;i++)
    {
      type=Framework[CurrentSystem].Atoms[f1][i].Type;
      charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged)
      {
        Charge[nr_of_coulombic_sites]=charge;
        if(!(Framework[CurrentSystem].Atoms[f1][i].Fixed.x&&Framework[CurrentSystem].Atoms[f1][i].Fixed.y&&Framework[CurrentSystem].Atoms[f1][i].Fixed.z))
        {
          energy_framework_self+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          AtomVector[nr_of_coulombic_sites]=&Framework[CurrentSystem].Atoms[f1][i].Force;
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_coulombic_sites_framework=nr_of_coulombic_sites;
  nr_molecules=NumberOfAdsorbateMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Adsorbates[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        scaling=Adsorbates[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
        charge=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        charge*=scaling;
        if(considered_charged)
        {
          Charge[nr_of_coulombic_sites]=charge;
          energy_adsorbate_self+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
          if(!(Adsorbates[CurrentSystem][i].Atoms[j].Fixed.x&&Adsorbates[CurrentSystem][i].Atoms[j].Fixed.y&&Adsorbates[CurrentSystem][i].Atoms[j].Fixed.z))
          {
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][i].Atoms[j].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            AtomVector[nr_of_coulombic_sites]=&Adsorbates[CurrentSystem][i].Atoms[j].Force;
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }
  nr_of_coulombic_sites_adsorbate=nr_of_coulombic_sites;
  nr_molecules=NumberOfCationMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Cations[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Cations[CurrentSystem][i].Atoms[j].Type;
        scaling=Cations[CurrentSystem][i].Atoms[j].CFChargeScalingParameter;
        charge=Cations[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        charge*=scaling;
        if(considered_charged)
        {
          Charge[nr_of_coulombic_sites]=charge;
          energy_cation_self+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
          if(!(Cations[CurrentSystem][i].Atoms[j].Fixed.x&&Cations[CurrentSystem][i].Atoms[j].Fixed.y&&Cations[CurrentSystem][i].Atoms[j].Fixed.z))
          {
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][i].Atoms[j].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            AtomVector[nr_of_coulombic_sites]=&Cations[CurrentSystem][i].Atoms[j].Force;
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }
  nr_of_coulombic_sites_cation=nr_of_coulombic_sites;

  nr_of_bonddipole_sites=0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      nr_of_bonddipoles=Framework[CurrentSystem].NumberOfBondDipoles[f1];
      for(i=0;i<nr_of_bonddipoles;i++)
      {
        pair=Framework[CurrentSystem].BondDipoles[f1][i];
        atom_pointer=Framework[CurrentSystem].Atoms[f1];
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        dipole=ApplyBoundaryCondition(dipole);
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Framework[CurrentSystem].BondDipoleMagnitude[f1][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        BondDipoleMagnitude[nr_of_bonddipole_sites]=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
        BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));  
        energy_framework_self_bd+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].Force;
        BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].Force;
        nr_of_bonddipole_sites++;
      }
    }
  }
  nr_of_bonddipole_sites_framework=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Adsorbates[CurrentSystem][i].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[j];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));  
      energy_adsorbate_self_bd+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].Force;
      BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].Force;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_adsorbate=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Cations[CurrentSystem][i].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[j];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));  
      energy_cation_self_bd+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].Force;
      BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].Force;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_site_cation=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          rk=kvecs[nvec];

          // loop over all the framework atoms
          sum_framework.re=0.0;
          sum_framework.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_framework.re+=temp*Eikr[i].re;
            sum_framework.im+=temp*Eikr[i].im;
          }

          // loop over all the adsorbate atoms
          sum_adsorbate.re=0.0;
          sum_adsorbate.im=0.0;
          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_adsorbate.re+=temp*Eikr[i].re;
            sum_adsorbate.im+=temp*Eikr[i].im;
          }

          // loop over all the cation atoms
          sum_cation.re=0.0;
          sum_cation.im=0.0;
          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_cation.re+=temp*Eikr[i].re;
            sum_cation.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_framework.re=0.0;
          sum_bonddipole_framework.im=0.0;

          for(i=0;i<nr_of_bonddipole_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_framework.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_framework.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_adsorbate.re=0.0;
          sum_bonddipole_adsorbate.im=0.0;
          for(;i<nr_of_bonddipole_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_adsorbate.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_adsorbate.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_cation.re=0.0;
          sum_bonddipole_cation.im=0.0;
          for(;i<nr_of_bonddipole_site_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_cation.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_cation.im+=temp*Eikr_bd[i].im;
          }

          // add the pre-computed contributions of fixed atoms
          sum_framework.re+=StoreRigidChargeFramework[CurrentSystem][nvec].re;
          sum_framework.im+=StoreRigidChargeFramework[CurrentSystem][nvec].im;
          sum_adsorbate.re+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].re;
          sum_adsorbate.im+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].im;
          sum_cation.re+=StoreRigidChargeCations[CurrentSystem][nvec].re;
          sum_cation.im+=StoreRigidChargeCations[CurrentSystem][nvec].im;

          sum_bonddipole_framework.re+=StoreRigidBondDipolesFramework[CurrentSystem][nvec].re;
          sum_bonddipole_framework.im+=StoreRigidBondDipolesFramework[CurrentSystem][nvec].im;
          sum_bonddipole_adsorbate.re+=StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec].re;
          sum_bonddipole_adsorbate.im+=StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec].im;
          sum_bonddipole_cation.re+=StoreRigidBondDipolesCations[CurrentSystem][nvec].re;
          sum_bonddipole_cation.im+=StoreRigidBondDipolesCations[CurrentSystem][nvec].im;


          // precomputed wavevector dependent pre-factor
          temp=kfactor[nvec];

          // charge-charge, charge-bonddipole, bonddipole-donddipole energies
          energy_framework_framework+=temp*(SQR(sum_framework.re)+SQR(sum_framework.im));
          energy_framework_framework_c_bd+=2.0*temp*(sum_framework.im*sum_bonddipole_framework.re-sum_framework.re*sum_bonddipole_framework.im);
          energy_framework_framework_bd+=temp*(SQR(sum_bonddipole_framework.re)+SQR(sum_bonddipole_framework.im));

          energy_adsorbate_adsorbate+=temp*(SQR(sum_adsorbate.re)+SQR(sum_adsorbate.im));
          energy_adsorbate_adsorbate_c_bd+=2.0*temp*(sum_adsorbate.im*sum_bonddipole_adsorbate.re-sum_adsorbate.re*sum_bonddipole_adsorbate.im);
          energy_adsorbate_adsorbate_bd+=temp*(SQR(sum_bonddipole_adsorbate.re)+SQR(sum_bonddipole_adsorbate.im));

          energy_cation_cation+=temp*(SQR(sum_cation.re)+SQR(sum_cation.im));
          energy_cation_cation_c_bd+=2.0*temp*(sum_cation.im*sum_bonddipole_cation.re-sum_cation.re*sum_bonddipole_cation.im);
          energy_cation_cation_bd+=temp*(SQR(sum_bonddipole_cation.re)+SQR(sum_bonddipole_cation.im));

          energy_framework_adsorbate+=temp*(sum_framework.re*sum_adsorbate.re+sum_framework.im*sum_adsorbate.im);
          energy_framework_adsorbate_c_bd+=temp*(sum_adsorbate.im*sum_bonddipole_framework.re+sum_framework.im*sum_bonddipole_adsorbate.re-
                                                (sum_adsorbate.re*sum_bonddipole_framework.im+sum_framework.re*sum_bonddipole_adsorbate.im));
          energy_framework_adsorbate_bd+=temp*(sum_bonddipole_framework.re*sum_bonddipole_adsorbate.re+sum_bonddipole_framework.im*sum_bonddipole_adsorbate.im);

          energy_framework_cation+=temp*(sum_framework.re*sum_cation.re+sum_framework.im*sum_cation.im);
          energy_framework_cation_c_bd+=temp*(sum_cation.im*sum_bonddipole_framework.re+sum_framework.im*sum_bonddipole_cation.re-
                                             (sum_cation.re*sum_bonddipole_framework.im+sum_framework.re*sum_bonddipole_cation.im));
          energy_framework_cation_bd+=temp*(sum_bonddipole_framework.re*sum_bonddipole_cation.re+sum_bonddipole_framework.im*sum_bonddipole_cation.im);

          energy_adsorbate_cation+=temp*(sum_adsorbate.re*sum_cation.re+sum_adsorbate.im*sum_cation.im);
          energy_adsorbate_cation_c_bd+=temp*(sum_cation.im*sum_bonddipole_adsorbate.re+sum_adsorbate.im*sum_bonddipole_cation.re-
                                             (sum_cation.re*sum_bonddipole_adsorbate.im+sum_adsorbate.re*sum_bonddipole_cation.im));
          energy_adsorbate_cation_bd+=temp*(sum_bonddipole_adsorbate.re*sum_bonddipole_cation.re+sum_bonddipole_adsorbate.im*sum_bonddipole_cation.im);


          // get total sums
          sum.re=sum_framework.re+sum_adsorbate.re+sum_cation.re;
          sum.im=sum_framework.im+sum_adsorbate.im+sum_cation.im;
          sum_bonddipole.re=sum_bonddipole_framework.re+sum_bonddipole_adsorbate.re+sum_bonddipole_cation.re;
          sum_bonddipole.im=sum_bonddipole_framework.im+sum_bonddipole_adsorbate.im+sum_bonddipole_cation.im;

          // forces on atoms from other charges and bond-dipoles
          for(i=0;i<nr_of_coulombic_sites_framework;i++)
          {
            fac=2.0*temp*Charge[i]*((Eikr[i].im*(sum.re-sum_cation.re)-Eikr[i].re*(sum.im-sum_cation.im))+
                (-Eikr[i].re*sum_bonddipole.re-Eikr[i].im*sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitAdsorbateAdsorbateCoulombInteractions)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }
          if(OmitAdsorbateCationCoulombInteractions)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }

          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            fac=2.0*temp*Charge[i]*((Eikr[i].im*temp_sum.re-Eikr[i].re*temp_sum.im)+
                (-Eikr[i].re*temp_sum_bonddipole.re-Eikr[i].im*temp_sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitCationCationCoulombInteractions)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }
          if(OmitAdsorbateCationCoulombInteractions)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }

          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            fac=2.0*temp*Charge[i]*((Eikr[i].im*temp_sum.re-Eikr[i].re*temp_sum.im)+
                (-Eikr[i].re*temp_sum_bonddipole.re-Eikr[i].im*temp_sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          // stress tensor
          current_energy=temp*(SQR(sum.re)+SQR(sum.im)+
                               2.0*(sum.im*sum_bonddipole.re-sum.re*sum_bonddipole.im)+
                               SQR(sum_bonddipole.re)+SQR(sum_bonddipole.im));

          // take omitted interactions into account
          if(Framework[CurrentSystem].FrameworkModel!=FLEXIBLE) 
             current_energy-=temp*(SQR(sum_framework.re)+SQR(sum_framework.im)+
                               2.0*(sum_framework.im*sum_bonddipole_framework.re-sum_framework.re*sum_bonddipole_framework.im)+
                               SQR(sum_bonddipole_framework.re)+SQR(sum_bonddipole_framework.im));
          if(OmitAdsorbateAdsorbateCoulombInteractions)
            current_energy-=temp*(SQR(sum_adsorbate.re)+SQR(sum_adsorbate.im)+
                                  2.0*(sum_adsorbate.im*sum_bonddipole_adsorbate.re-sum_adsorbate.re*sum_bonddipole_adsorbate.im)+
                                  SQR(sum_bonddipole_adsorbate.re)+SQR(sum_bonddipole_adsorbate.im));
          if(OmitCationCationCoulombInteractions)
            current_energy-=temp*(SQR(sum_cation.re)+SQR(sum_cation.im)+
                                  2.0*(sum_cation.im*sum_bonddipole_cation.re-sum_cation.re*sum_bonddipole_cation.im)+
                                  SQR(sum_bonddipole_cation.re)+SQR(sum_bonddipole_cation.im));
          if(OmitAdsorbateCationCoulombInteractions)
            current_energy-=2.0*temp*(sum_adsorbate.re*sum_cation.re+sum_adsorbate.im*sum_cation.im+
                                      sum_adsorbate.im*sum_bonddipole_cation.re-sum_adsorbate.re*sum_bonddipole_cation.im+
                                      sum_cation.im*sum_bonddipole_adsorbate.re-sum_cation.re*sum_bonddipole_adsorbate.im+
                                      sum_bonddipole_adsorbate.re*sum_bonddipole_cation.re+sum_bonddipole_adsorbate.im*sum_bonddipole_cation.im);
          fac=2.0*(1.0/(SQR(rk.x)+SQR(rk.y)+SQR(rk.z))+(0.25/SQR(alpha)))*current_energy;
          stress.ax+=current_energy-fac*rk.x*rk.x;
          stress.ay+=-fac*rk.x*rk.y;
          stress.az+=-fac*rk.x*rk.z;

          stress.bx+=-fac*rk.y*rk.x;
          stress.by+=current_energy-fac*rk.y*rk.y;
          stress.bz+=-fac*rk.y*rk.z;

          stress.cx+=-fac*rk.z*rk.x;
          stress.cy+=-fac*rk.z*rk.y;
          stress.cz+=current_energy-fac*rk.z*rk.z;

          // forces on bond-dipoles from charges and other bond-dipoles
          for(i=0;i<nr_of_bonddipole_sites_framework;i++)
          {
            dipole=DipoleVector[i];
            DipoleMagnitudeA=BondDipoleMagnitude[i];
            length=BondLength[i];
            dot_product=(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);

            // charge-bonddipole contribution
            fac1=temp*sum.re*dot_product*Eikr_bd[i].re;
            fac2=2.0*temp*sum.re*Eikr_bd[i].im*DipoleMagnitudeA/length;
            fac3=2.0*temp*sum.re*Eikr_bd[i].im*dot_product/(length*DipoleMagnitudeA);

            fac4=temp*sum.im*dot_product*Eikr_bd[i].im;
            fac5=2.0*temp*sum.im*Eikr_bd[i].re*DipoleMagnitudeA/length;
            fac6=2.0*temp*sum.im*Eikr_bd[i].re*dot_product/(length*DipoleMagnitudeA);

            fa1.x=(fac1-fac2)*rk.x+fac3*dipole.x+(fac4+fac5)*rk.x-fac6*dipole.x;
            fa1.y=(fac1-fac2)*rk.y+fac3*dipole.y+(fac4+fac5)*rk.y-fac6*dipole.y;
            fa1.z=(fac1-fac2)*rk.z+fac3*dipole.z+(fac4+fac5)*rk.z-fac6*dipole.z;

            fa2.x=(fac1+fac2)*rk.x-fac3*dipole.x+(fac4-fac5)*rk.x+fac6*dipole.x;
            fa2.y=(fac1+fac2)*rk.y-fac3*dipole.y+(fac4-fac5)*rk.y+fac6*dipole.y;
            fa2.z=(fac1+fac2)*rk.z-fac3*dipole.z+(fac4-fac5)*rk.z+fac6*dipole.z;

            // bonddipole-bonddipole contribution
            fac1=temp*dot_product*(sum_bonddipole.re*Eikr_bd[i].im-sum_bonddipole.im*Eikr_bd[i].re);
            fac2=2.0*temp*(sum_bonddipole.re*Eikr_bd[i].re+sum_bonddipole.im*Eikr_bd[i].im)*DipoleMagnitudeA/length;
            fac3=2.0*temp*(sum_bonddipole.re*Eikr_bd[i].re+sum_bonddipole.im*Eikr_bd[i].im)*dot_product/
                 (length*DipoleMagnitudeA);

            fa1.x+=(fac1+fac2)*rk.x-fac3*dipole.x;
            fa1.y+=(fac1+fac2)*rk.y-fac3*dipole.y;
            fa1.z+=(fac1+fac2)*rk.z-fac3*dipole.z;

            fa2.x+=(fac1-fac2)*rk.x+fac3*dipole.x;
            fa2.y+=(fac1-fac2)*rk.y+fac3*dipole.y;
            fa2.z+=(fac1-fac2)*rk.z+fac3*dipole.z;

            BondDipoleForcesA[i]->x+=fa1.x;
            BondDipoleForcesA[i]->y+=fa1.y;
            BondDipoleForcesA[i]->z+=fa1.z;

            BondDipoleForcesB[i]->x+=fa2.x;
            BondDipoleForcesB[i]->y+=fa2.y;
            BondDipoleForcesB[i]->z+=fa2.z;

            // bond-dipole contribution to the stress
            fac=temp*(sum.im*Eikr_bd[i].re-sum.re*Eikr_bd[i].im
                      +sum_bonddipole.im*Eikr_bd[i].im+sum_bonddipole.re*Eikr_bd[i].re);
            stress.ax+=2.0*fac*dipole.x*rk.x;
            stress.ay+=fac*(dipole.x*rk.y+dipole.y*rk.x);
            stress.az+=fac*(dipole.x*rk.z+dipole.z*rk.x);
            stress.bx+=fac*(dipole.y*rk.x+dipole.x*rk.y);
            stress.by+=2.0*fac*dipole.y*rk.y;
            stress.bz+=fac*(dipole.y*rk.z+dipole.z*rk.y);
            stress.cx+=fac*(dipole.z*rk.x+dipole.x*rk.z);
            stress.cy+=fac*(dipole.z*rk.y+dipole.y*rk.z);
            stress.cz+=2.0*fac*dipole.z*rk.z;

            // convert forces on atoms to 'molecular stress'
            // it is equivalent to: v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x); etc.
            fac=0.5*(length/DipoleMagnitudeA);
            stress.ax+=fac*(fa2.x-fa1.x)*dipole.x;
            stress.bx+=fac*(fa2.y-fa1.y)*dipole.x;
            stress.cx+=fac*(fa2.z-fa1.z)*dipole.x;
            stress.ay+=fac*(fa2.x-fa1.x)*dipole.y;
            stress.by+=fac*(fa2.y-fa1.y)*dipole.y;
            stress.cy+=fac*(fa2.z-fa1.z)*dipole.y;
            stress.az+=fac*(fa2.x-fa1.x)*dipole.z;
            stress.bz+=fac*(fa2.y-fa1.y)*dipole.z;
            stress.cz+=fac*(fa2.z-fa1.z)*dipole.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitAdsorbateAdsorbateCoulombInteractions)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }
          if(OmitAdsorbateCationCoulombInteractions)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }

          for(;i<nr_of_bonddipole_sites_adsorbate;i++)
          {
            dipole=DipoleVector[i];
            DipoleMagnitudeA=BondDipoleMagnitude[i];
            length=BondLength[i];
            dot_product=(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);

            // charge-bonddipole contribution
            fac1=temp*temp_sum.re*dot_product*Eikr_bd[i].re;
            fac2=2.0*temp*temp_sum.re*Eikr_bd[i].im*DipoleMagnitudeA/length;
            fac3=2.0*temp*temp_sum.re*Eikr_bd[i].im*dot_product/(length*DipoleMagnitudeA);

            fac4=temp*temp_sum.im*dot_product*Eikr_bd[i].im;
            fac5=2.0*temp*temp_sum.im*Eikr_bd[i].re*DipoleMagnitudeA/length;
            fac6=2.0*temp*temp_sum.im*Eikr_bd[i].re*dot_product/(length*DipoleMagnitudeA);

            fa1.x=(fac1-fac2)*rk.x+fac3*dipole.x+(fac4+fac5)*rk.x-fac6*dipole.x;
            fa1.y=(fac1-fac2)*rk.y+fac3*dipole.y+(fac4+fac5)*rk.y-fac6*dipole.y;
            fa1.z=(fac1-fac2)*rk.z+fac3*dipole.z+(fac4+fac5)*rk.z-fac6*dipole.z;

            fa2.x=(fac1+fac2)*rk.x-fac3*dipole.x+(fac4-fac5)*rk.x+fac6*dipole.x;
            fa2.y=(fac1+fac2)*rk.y-fac3*dipole.y+(fac4-fac5)*rk.y+fac6*dipole.y;
            fa2.z=(fac1+fac2)*rk.z-fac3*dipole.z+(fac4-fac5)*rk.z+fac6*dipole.z;

            // bonddipole-bonddipole contribution
            fac1=temp*dot_product*(temp_sum_bonddipole.re*Eikr_bd[i].im-temp_sum_bonddipole.im*Eikr_bd[i].re);
            fac2=2.0*temp*(temp_sum_bonddipole.re*Eikr_bd[i].re+temp_sum_bonddipole.im*Eikr_bd[i].im)*DipoleMagnitudeA/length;
            fac3=2.0*temp*(temp_sum_bonddipole.re*Eikr_bd[i].re+temp_sum_bonddipole.im*Eikr_bd[i].im)*dot_product/
                 (length*DipoleMagnitudeA);

            fa1.x+=(fac1+fac2)*rk.x-fac3*dipole.x;
            fa1.y+=(fac1+fac2)*rk.y-fac3*dipole.y;
            fa1.z+=(fac1+fac2)*rk.z-fac3*dipole.z;

            fa2.x+=(fac1-fac2)*rk.x+fac3*dipole.x;
            fa2.y+=(fac1-fac2)*rk.y+fac3*dipole.y;
            fa2.z+=(fac1-fac2)*rk.z+fac3*dipole.z;

            BondDipoleForcesA[i]->x+=fa1.x;
            BondDipoleForcesA[i]->y+=fa1.y;
            BondDipoleForcesA[i]->z+=fa1.z;

            BondDipoleForcesB[i]->x+=fa2.x;
            BondDipoleForcesB[i]->y+=fa2.y;
            BondDipoleForcesB[i]->z+=fa2.z;

            // bond-dipole contribution to the stress
            fac=temp*(temp_sum.im*Eikr_bd[i].re-temp_sum.re*Eikr_bd[i].im
                      +temp_sum_bonddipole.im*Eikr_bd[i].im+temp_sum_bonddipole.re*Eikr_bd[i].re);
            stress.ax+=2.0*fac*dipole.x*rk.x;
            stress.ay+=fac*(dipole.x*rk.y+dipole.y*rk.x);
            stress.az+=fac*(dipole.x*rk.z+dipole.z*rk.x);
            stress.bx+=fac*(dipole.y*rk.x+dipole.x*rk.y);
            stress.by+=2.0*fac*dipole.y*rk.y;
            stress.bz+=fac*(dipole.y*rk.z+dipole.z*rk.y);
            stress.cx+=fac*(dipole.z*rk.x+dipole.x*rk.z);
            stress.cy+=fac*(dipole.z*rk.y+dipole.y*rk.z);
            stress.cz+=2.0*fac*dipole.z*rk.z;

            // convert forces on atoms to 'molecular stress'
            // it is equivalent to: v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x); etc.
            fac=0.5*(length/DipoleMagnitudeA);
            stress.ax+=fac*(fa2.x-fa1.x)*dipole.x;
            stress.bx+=fac*(fa2.y-fa1.y)*dipole.x;
            stress.cx+=fac*(fa2.z-fa1.z)*dipole.x;
            stress.ay+=fac*(fa2.x-fa1.x)*dipole.y;
            stress.by+=fac*(fa2.y-fa1.y)*dipole.y;
            stress.cy+=fac*(fa2.z-fa1.z)*dipole.y;
            stress.az+=fac*(fa2.x-fa1.x)*dipole.z;
            stress.bz+=fac*(fa2.y-fa1.y)*dipole.z;
            stress.cz+=fac*(fa2.z-fa1.z)*dipole.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitCationCationCoulombInteractions)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }
          if(OmitAdsorbateCationCoulombInteractions)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }

          for(;i<nr_of_bonddipole_site_cation;i++)
          {
            dipole=DipoleVector[i];
            DipoleMagnitudeA=BondDipoleMagnitude[i];
            length=BondLength[i];
            dot_product=(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);

            // charge-bonddipole contribution
            fac1=temp*temp_sum.re*dot_product*Eikr_bd[i].re;
            fac2=2.0*temp*temp_sum.re*Eikr_bd[i].im*DipoleMagnitudeA/length;
            fac3=2.0*temp*temp_sum.re*Eikr_bd[i].im*dot_product/(length*DipoleMagnitudeA);

            fac4=temp*temp_sum.im*dot_product*Eikr_bd[i].im;
            fac5=2.0*temp*temp_sum.im*Eikr_bd[i].re*DipoleMagnitudeA/length;
            fac6=2.0*temp*temp_sum.im*Eikr_bd[i].re*dot_product/(length*DipoleMagnitudeA);

            fa1.x=(fac1-fac2)*rk.x+fac3*dipole.x+(fac4+fac5)*rk.x-fac6*dipole.x;
            fa1.y=(fac1-fac2)*rk.y+fac3*dipole.y+(fac4+fac5)*rk.y-fac6*dipole.y;
            fa1.z=(fac1-fac2)*rk.z+fac3*dipole.z+(fac4+fac5)*rk.z-fac6*dipole.z;

            fa2.x=(fac1+fac2)*rk.x-fac3*dipole.x+(fac4-fac5)*rk.x+fac6*dipole.x;
            fa2.y=(fac1+fac2)*rk.y-fac3*dipole.y+(fac4-fac5)*rk.y+fac6*dipole.y;
            fa2.z=(fac1+fac2)*rk.z-fac3*dipole.z+(fac4-fac5)*rk.z+fac6*dipole.z;

            // bonddipole-bonddipole contribution
            fac1=temp*dot_product*(temp_sum_bonddipole.re*Eikr_bd[i].im-temp_sum_bonddipole.im*Eikr_bd[i].re);
            fac2=2.0*temp*(temp_sum_bonddipole.re*Eikr_bd[i].re+temp_sum_bonddipole.im*Eikr_bd[i].im)*DipoleMagnitudeA/length;
            fac3=2.0*temp*(temp_sum_bonddipole.re*Eikr_bd[i].re+temp_sum_bonddipole.im*Eikr_bd[i].im)*dot_product/
                 (length*DipoleMagnitudeA);

            fa1.x+=(fac1+fac2)*rk.x-fac3*dipole.x;
            fa1.y+=(fac1+fac2)*rk.y-fac3*dipole.y;
            fa1.z+=(fac1+fac2)*rk.z-fac3*dipole.z;

            fa2.x+=(fac1-fac2)*rk.x+fac3*dipole.x;
            fa2.y+=(fac1-fac2)*rk.y+fac3*dipole.y;
            fa2.z+=(fac1-fac2)*rk.z+fac3*dipole.z;

            BondDipoleForcesA[i]->x+=fa1.x;
            BondDipoleForcesA[i]->y+=fa1.y;
            BondDipoleForcesA[i]->z+=fa1.z;

            BondDipoleForcesB[i]->x+=fa2.x;
            BondDipoleForcesB[i]->y+=fa2.y;
            BondDipoleForcesB[i]->z+=fa2.z;

            // bond-dipole contribution to the stress
            fac=temp*(temp_sum.im*Eikr_bd[i].re-temp_sum.re*Eikr_bd[i].im
                      +temp_sum_bonddipole.im*Eikr_bd[i].im+temp_sum_bonddipole.re*Eikr_bd[i].re);
            stress.ax+=2.0*fac*dipole.x*rk.x;
            stress.ay+=fac*(dipole.x*rk.y+dipole.y*rk.x);
            stress.az+=fac*(dipole.x*rk.z+dipole.z*rk.x);
            stress.bx+=fac*(dipole.y*rk.x+dipole.x*rk.y);
            stress.by+=2.0*fac*dipole.y*rk.y;
            stress.bz+=fac*(dipole.y*rk.z+dipole.z*rk.y);
            stress.cx+=fac*(dipole.z*rk.x+dipole.x*rk.z);
            stress.cy+=fac*(dipole.z*rk.y+dipole.y*rk.z);
            stress.cz+=2.0*fac*dipole.z*rk.z;

            // convert forces on atoms to 'molecular stress'
            // it is equivalent to: v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x); etc.
            fac=0.5*(length/DipoleMagnitudeA);
            stress.ax+=fac*(fa2.x-fa1.x)*dipole.x;
            stress.bx+=fac*(fa2.y-fa1.y)*dipole.x;
            stress.cx+=fac*(fa2.z-fa1.z)*dipole.x;
            stress.ay+=fac*(fa2.x-fa1.x)*dipole.y;
            stress.by+=fac*(fa2.y-fa1.y)*dipole.y;
            stress.cy+=fac*(fa2.z-fa1.z)*dipole.y;
            stress.az+=fac*(fa2.x-fa1.x)*dipole.z;
            stress.bz+=fac*(fa2.y-fa1.y)*dipole.z;
            stress.cz+=fac*(fa2.z-fa1.z)*dipole.z;
          }

          // next wave-vector
          nvec++;
        }
      }
    }
  }

  // Exclusion pairs for the framework
  // =================================
  // The list of exclusion pairs is pre-compute at the start of the simulation based on e.g. 1-2, 1-3 exclusion definitions.
  // Note: there is no cutoff used here

  energy_framework_excluded=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        atom_pointer=Framework[CurrentSystem].Atoms[f1];
        pair=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        chargeA=atom_pointer[pair.A].Charge;
        chargeB=atom_pointer[pair.B].Charge;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;

        if(r>1e-4)
        {
          energy_framework_excluded-=chargeA*chargeB*Bt0;

          temp=chargeA*chargeB*Bt1;

          atom_pointer[pair.A].Force.x+=temp*dr.x;
          atom_pointer[pair.A].Force.y+=temp*dr.y;
          atom_pointer[pair.A].Force.z+=temp*dr.z;

          atom_pointer[pair.B].Force.x-=temp*dr.x;
          atom_pointer[pair.B].Force.y-=temp*dr.y;
          atom_pointer[pair.B].Force.z-=temp*dr.z;

          stress.ax+=temp*dr.x*dr.x;
          stress.bx+=temp*dr.y*dr.x;
          stress.cx+=temp*dr.z*dr.x;

          stress.ay+=temp*dr.x*dr.y;
          stress.by+=temp*dr.y*dr.y;
          stress.cy+=temp*dr.z*dr.y;

          stress.az+=temp*dr.x*dr.z;
          stress.bz+=temp*dr.y*dr.z;
          stress.cz+=temp*dr.z*dr.z;

        }
        else // if r->0 compute limiting value to avoid divergence when shell overlaps with core in core-shell models
          energy_framework_excluded+=chargeA*chargeB*alpha*(2.0/sqrt(M_PI));
      }
    }
  }

  energy_framework_excluded_c_bd=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeBondDipole[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        // the index of the bonddipole is the second index 'B'
        A=Framework[CurrentSystem].ExcludedIntraChargeBondDipole[f1][i].A;
        B=Framework[CurrentSystem].ExcludedIntraChargeBondDipole[f1][i].B;

        atom_pointer=Framework[CurrentSystem].Atoms[f1];

        typeA=atom_pointer[A].Type;
        chargeA=atom_pointer[A].Charge;
        posA=atom_pointer[A].Position;

        pairB=Framework[CurrentSystem].BondDipoles[f1][B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        dipole_magnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][B];
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posB.x-posA.x;
        dr.y=posB.y-posA.y;
        dr.z=posB.z-posA.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;

        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_framework_excluded_c_bd+=Bt1*chargeA*cosB;

        fa1.x=chargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
        fa1.y=chargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
        fa1.z=chargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

        atom_pointer[A].Force.x+=fa1.x;
        atom_pointer[A].Force.y+=fa1.y;
        atom_pointer[A].Force.z+=fa1.z;

        fac1=chargeA*Bt1*dipole_magnitudeB/lengthB;
        fac2=chargeA*Bt1*cosB/(dipole_magnitudeB*lengthB);

        fb1.x=0.5*fa1.x+fac1*dr.x-fac2*dipoleB.x;
        fb1.y=0.5*fa1.y+fac1*dr.y-fac2*dipoleB.y;
        fb1.z=0.5*fa1.z+fac1*dr.z-fac2*dipoleB.z;

        fb2.x=0.5*fa1.x-fac1*dr.x+fac2*dipoleB.x;
        fb2.y=0.5*fa1.y-fac1*dr.y+fac2*dipoleB.y;
        fb2.z=0.5*fa1.z-fac1*dr.z+fac2*dipoleB.z;

        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=fa1.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=fa1.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=fa1.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=fa1.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=fa1.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=fa1.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=fa1.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=fa1.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=fa1.z*dr.z+v.cz;
      }
    }
  }

  energy_framework_excluded_bd=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraBondDipoleBondDipole[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Framework[CurrentSystem].ExcludedIntraBondDipoleBondDipole[f1][i].A;
        B=Framework[CurrentSystem].ExcludedIntraBondDipoleBondDipole[f1][i].B;

        atom_pointer=Framework[CurrentSystem].Atoms[f1];

        pairA=Framework[CurrentSystem].BondDipoles[f1][A];
        posA1=atom_pointer[pairA.A].Position;
        posA2=atom_pointer[pairA.B].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        dipoleA=ApplyBoundaryCondition(dipoleA);
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        dipole_magnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][A];
        lengthA=sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
        temp=dipole_magnitudeA/lengthA;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        pairB=Framework[CurrentSystem].BondDipoles[f1][B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        dipole_magnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][B];
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;
        temp*=2.0*SQR(alpha);
        Bt3=temp+(5.0/rr)*Bt2;

        cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
        cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_framework_excluded_bd-=(Bt1*cosAB-Bt2*cosA*cosB);

        term.x=(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
        term.y=(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
        term.z=(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

        termA.x=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.x/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.x-Bt1*dipoleB.x)*dipole_magnitudeA/lengthA;
        termA.y=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.y/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.y-Bt1*dipoleB.y)*dipole_magnitudeA/lengthA;
        termA.z=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.z/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.z-Bt1*dipoleB.z)*dipole_magnitudeA/lengthA;

        termB.x=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.x/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.x-Bt1*dipoleA.x)*dipole_magnitudeB/lengthB;
        termB.y=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.y/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.y-Bt1*dipoleA.y)*dipole_magnitudeB/lengthB;
        termB.z=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.z/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.z-Bt1*dipoleA.z)*dipole_magnitudeB/lengthB;

        fa1.x=0.5*term.x+termA.x;
        fa1.y=0.5*term.y+termA.y;
        fa1.z=0.5*term.z+termA.z;
        fa2.x=0.5*term.x-termA.x;
        fa2.y=0.5*term.y-termA.y;
        fa2.z=0.5*term.z-termA.z;

        fb1.x=-0.5*term.x+termB.x;
        fb1.y=-0.5*term.y+termB.y;
        fb1.z=-0.5*term.z+termB.z;
        fb2.x=-0.5*term.x-termB.x;
        fb2.y=-0.5*term.y-termB.y;
        fb2.z=-0.5*term.z-termB.z;

        atom_pointer[pairA.A].Force.x-=fa1.x;
        atom_pointer[pairA.A].Force.y-=fa1.y;
        atom_pointer[pairA.A].Force.z-=fa1.z;

        atom_pointer[pairA.B].Force.x-=fa2.x;
        atom_pointer[pairA.B].Force.y-=fa2.y;
        atom_pointer[pairA.B].Force.z-=fa2.z;

        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
      }
    }
  }

  // Exclusion pairs for the adsorbates
  // ==================================

  // if adsorbate-adsorbate interactions are omitted, then the intra-molecular interactions are already fully elimated
  if(!OmitAdsorbateAdsorbateCoulombInteractions)
  {
    energy_adsorbate_excluded=0.0;
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      atom_pointer=Adsorbates[CurrentSystem][m].Atoms;
      type=Adsorbates[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        pair=Components[type].ExcludedIntraChargeCharge[i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        scalingA=atom_pointer[pair.A].CFChargeScalingParameter;
        scalingB=atom_pointer[pair.B].CFChargeScalingParameter;
        chargeA=scalingA*atom_pointer[pair.A].Charge;
        chargeB=scalingB*atom_pointer[pair.B].Charge;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;

        energy_adsorbate_excluded-=chargeA*chargeB*Bt0;

        temp=chargeA*chargeB*Bt1;

        atom_pointer[pair.A].Force.x+=temp*dr.x;
        atom_pointer[pair.A].Force.y+=temp*dr.y;
        atom_pointer[pair.A].Force.z+=temp*dr.z;

        atom_pointer[pair.B].Force.x-=temp*dr.x;
        atom_pointer[pair.B].Force.y-=temp*dr.y;
        atom_pointer[pair.B].Force.z-=temp*dr.z;

        stress.ax+=temp*dr.x*dr.x;
        stress.bx+=temp*dr.y*dr.x;
        stress.cx+=temp*dr.z*dr.x;

        stress.ay+=temp*dr.x*dr.y;
        stress.by+=temp*dr.y*dr.y;
        stress.cy+=temp*dr.z*dr.y;

        stress.az+=temp*dr.x*dr.z;
        stress.bz+=temp*dr.y*dr.z;
        stress.cz+=temp*dr.z*dr.z;
      }
    }

    energy_adsorbate_excluded_c_bd=0.0;
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      type=Adsorbates[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeBondDipole;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Components[type].ExcludedIntraChargeBondDipole[i].A;
        B=Components[type].ExcludedIntraChargeBondDipole[i].B;

        atom_pointer=Adsorbates[CurrentSystem][m].Atoms;

        chargeA=Components[type].Charge[A];
        posA=atom_pointer[A].Position;

        pairB=Components[type].BondDipoles[B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipole_magnitudeB=Components[type].BondDipoleMagnitude[B];
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posB.x-posA.x;
        dr.y=posB.y-posA.y;
        dr.z=posB.z-posA.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;

        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_adsorbate_excluded_c_bd+=Bt1*chargeA*cosB;

        fa1.x=chargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
        fa1.y=chargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
        fa1.z=chargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

        atom_pointer[A].Force.x+=fa1.x;
        atom_pointer[A].Force.y+=fa1.y;
        atom_pointer[A].Force.z+=fa1.z;

        fac1=chargeA*Bt1*dipole_magnitudeB/lengthB;
        fac2=chargeA*Bt1*cosB/(dipole_magnitudeB*lengthB);

        fb1.x=0.5*fa1.x+fac1*dr.x-fac2*dipoleB.x;
        fb1.y=0.5*fa1.y+fac1*dr.y-fac2*dipoleB.y;
        fb1.z=0.5*fa1.z+fac1*dr.z-fac2*dipoleB.z;

        fb2.x=0.5*fa1.x-fac1*dr.x+fac2*dipoleB.x;
        fb2.y=0.5*fa1.y-fac1*dr.y+fac2*dipoleB.y;
        fb2.z=0.5*fa1.z-fac1*dr.z+fac2*dipoleB.z;

        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=fa1.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=fa1.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=fa1.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=fa1.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=fa1.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=fa1.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=fa1.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=fa1.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=fa1.z*dr.z+v.cz;
      }
    }

    energy_adsorbate_excluded_bd=0.0;
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      type=Adsorbates[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraBondDipoleBondDipole;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Components[type].ExcludedIntraBondDipoleBondDipole[i].A;
        B=Components[type].ExcludedIntraBondDipoleBondDipole[i].B;

        atom_pointer=Adsorbates[CurrentSystem][m].Atoms;

        pairA=Components[type].BondDipoles[A];
        posA1=atom_pointer[pairA.A].Position;
        posA2=atom_pointer[pairA.B].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        dipole_magnitudeA=Components[type].BondDipoleMagnitude[A];
        lengthA=sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
        temp=dipole_magnitudeA/lengthA;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        pairB=Components[type].BondDipoles[B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipole_magnitudeB=Components[type].BondDipoleMagnitude[B];
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;
        temp*=2.0*SQR(alpha);
        Bt3=temp+(5.0/rr)*Bt2;

        cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
        cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_adsorbate_excluded_bd-=Bt1*cosAB-Bt2*cosA*cosB;

        term.x=(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
        term.y=(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
        term.z=(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

        termA.x=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.x/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.x-Bt1*dipoleB.x)*dipole_magnitudeA/lengthA;
        termA.y=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.y/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.y-Bt1*dipoleB.y)*dipole_magnitudeA/lengthA;
        termA.z=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.z/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.z-Bt1*dipoleB.z)*dipole_magnitudeA/lengthA;

        termB.x=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.x/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.x-Bt1*dipoleA.x)*dipole_magnitudeB/lengthB;
        termB.y=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.y/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.y-Bt1*dipoleA.y)*dipole_magnitudeB/lengthB;
        termB.z=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.z/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.z-Bt1*dipoleA.z)*dipole_magnitudeB/lengthB;

        fa1.x=0.5*term.x+termA.x;
        fa1.y=0.5*term.y+termA.y;
        fa1.z=0.5*term.z+termA.z;
        fa2.x=0.5*term.x-termA.x;
        fa2.y=0.5*term.y-termA.y;
        fa2.z=0.5*term.z-termA.z;

        fb1.x=-0.5*term.x+termB.x;
        fb1.y=-0.5*term.y+termB.y;
        fb1.z=-0.5*term.z+termB.z;
        fb2.x=-0.5*term.x-termB.x;
        fb2.y=-0.5*term.y-termB.y;
        fb2.z=-0.5*term.z-termB.z;

        atom_pointer[pairA.A].Force.x-=fa1.x;
        atom_pointer[pairA.A].Force.y-=fa1.y;
        atom_pointer[pairA.A].Force.z-=fa1.z;

        atom_pointer[pairA.B].Force.x-=fa2.x;
        atom_pointer[pairA.B].Force.y-=fa2.y;
        atom_pointer[pairA.B].Force.z-=fa2.z;

        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
      }
    }
  }

  // Exclusion pairs for the cations
  // ===============================

  // if cation-cation interactions are omitted, then the intra-molecular interactions are already fully elimated
  if(!OmitCationCationCoulombInteractions)
  {
    energy_cation_excluded=0.0;
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      atom_pointer=Cations[CurrentSystem][m].Atoms;
      type=Cations[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        pair=Components[type].ExcludedIntraChargeCharge[i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        scalingA=atom_pointer[pair.A].CFChargeScalingParameter;
        scalingB=atom_pointer[pair.B].CFChargeScalingParameter;
        chargeA=scalingA*atom_pointer[pair.A].Charge;
        chargeB=scalingB*atom_pointer[pair.B].Charge;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;

        energy_cation_excluded-=chargeA*chargeB*Bt0;

        temp=chargeA*chargeB*Bt1;

        atom_pointer[pair.A].Force.x+=temp*dr.x;
        atom_pointer[pair.A].Force.y+=temp*dr.y;
        atom_pointer[pair.A].Force.z+=temp*dr.z;

        atom_pointer[pair.B].Force.x-=temp*dr.x;
        atom_pointer[pair.B].Force.y-=temp*dr.y;
        atom_pointer[pair.B].Force.z-=temp*dr.z;

        stress.ax+=temp*dr.x*dr.x;
        stress.bx+=temp*dr.y*dr.x;
        stress.cx+=temp*dr.z*dr.x;

        stress.ay+=temp*dr.x*dr.y;
        stress.by+=temp*dr.y*dr.y;
        stress.cy+=temp*dr.z*dr.y;

        stress.az+=temp*dr.x*dr.z;
        stress.bz+=temp*dr.y*dr.z;
        stress.cz+=temp*dr.z*dr.z;
      }
    }

    energy_cation_excluded_c_bd=0.0;
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      type=Cations[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeBondDipole;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Components[type].ExcludedIntraChargeBondDipole[i].A;
        B=Components[type].ExcludedIntraChargeBondDipole[i].B;

        atom_pointer=Cations[CurrentSystem][m].Atoms;

        chargeA=atom_pointer[A].Charge;
        posA=atom_pointer[A].Position;

        pairB=Components[type].BondDipoles[B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipole_magnitudeB=Components[type].BondDipoleMagnitude[B];
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posB.x-posA.x;
        dr.y=posB.y-posA.y;
        dr.z=posB.z-posA.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;

        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_cation_excluded_c_bd+=Bt1*chargeA*cosB;

        term.x=chargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
        term.y=chargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
        term.z=chargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

        fa1.x=term.x;
        fa1.y=term.y;
        fa1.z=term.z;
        atom_pointer[A].Force.x+=fa1.x;
        atom_pointer[A].Force.y+=fa1.y;
        atom_pointer[A].Force.z+=fa1.z;

        fac1=chargeA*Bt1*dipole_magnitudeB/lengthB;
        fac2=chargeA*Bt1*cosB/(dipole_magnitudeB*lengthB);

        fb1.x=0.5*term.x+fac1*dr.x-fac2*dipoleB.x;
        fb1.y=0.5*term.y+fac1*dr.y-fac2*dipoleB.y;
        fb1.z=0.5*term.z+fac1*dr.z-fac2*dipoleB.z;

        fb2.x=0.5*term.x-fac1*dr.x+fac2*dipoleB.x;
        fb2.y=0.5*term.y-fac1*dr.y+fac2*dipoleB.y;
        fb2.z=0.5*term.z-fac1*dr.z+fac2*dipoleB.z;


        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=fa1.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=fa1.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=fa1.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=fa1.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=fa1.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=fa1.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=fa1.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=fa1.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=fa1.z*dr.z+v.cz;
      }
    }

    energy_cation_excluded_bd=0.0;
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      type=Cations[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraBondDipoleBondDipole;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Components[type].ExcludedIntraBondDipoleBondDipole[i].A;
        B=Components[type].ExcludedIntraBondDipoleBondDipole[i].B;

        atom_pointer=Cations[CurrentSystem][m].Atoms;

        pairA=Components[type].BondDipoles[A];
        posA1=atom_pointer[pairA.A].Position;
        posA2=atom_pointer[pairA.B].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        dipole_magnitudeA=Components[type].BondDipoleMagnitude[A];
        lengthA=sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
        temp=dipole_magnitudeA/lengthA;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        pairB=Components[type].BondDipoles[B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipole_magnitudeB=Components[type].BondDipoleMagnitude[B];
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;
        temp*=2.0*SQR(alpha);
        Bt3=temp+(5.0/rr)*Bt2;

        cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
        cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_cation_excluded_bd-=Bt1*cosAB-Bt2*cosA*cosB;

        term.x=(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
        term.y=(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
        term.z=(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

        termA.x=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.x/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.x-Bt1*dipoleB.x)*dipole_magnitudeA/lengthA;
        termA.y=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.y/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.y-Bt1*dipoleB.y)*dipole_magnitudeA/lengthA;
        termA.z=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.z/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.z-Bt1*dipoleB.z)*dipole_magnitudeA/lengthA;

        termB.x=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.x/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.x-Bt1*dipoleA.x)*dipole_magnitudeB/lengthB;
        termB.y=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.y/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.y-Bt1*dipoleA.y)*dipole_magnitudeB/lengthB;
        termB.z=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.z/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.z-Bt1*dipoleA.z)*dipole_magnitudeB/lengthB;

        fa1.x=0.5*term.x+termA.x;
        fa1.y=0.5*term.y+termA.y;
        fa1.z=0.5*term.z+termA.z;
        fa2.x=0.5*term.x-termA.x;
        fa2.y=0.5*term.y-termA.y;
        fa2.z=0.5*term.z-termA.z;

        fb1.x=-0.5*term.x+termB.x;
        fb1.y=-0.5*term.y+termB.y;
        fb1.z=-0.5*term.z+termB.z;
        fb2.x=-0.5*term.x-termB.x;
        fb2.y=-0.5*term.y-termB.y;
        fb2.z=-0.5*term.z-termB.z;

        atom_pointer[pairA.A].Force.x-=fa1.x;
        atom_pointer[pairA.A].Force.y-=fa1.y;
        atom_pointer[pairA.A].Force.z-=fa1.z;

        atom_pointer[pairA.B].Force.x-=fa2.x;
        atom_pointer[pairA.B].Force.y-=fa2.y;
        atom_pointer[pairA.B].Force.z-=fa2.z;

        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
      }
    }
  }

  StrainDerivativeTensor[CurrentSystem].ax-=stress.ax;
  StrainDerivativeTensor[CurrentSystem].ay-=stress.ay;
  StrainDerivativeTensor[CurrentSystem].az-=stress.az;

  StrainDerivativeTensor[CurrentSystem].bx-=stress.bx;
  StrainDerivativeTensor[CurrentSystem].by-=stress.by;
  StrainDerivativeTensor[CurrentSystem].bz-=stress.bz;

  StrainDerivativeTensor[CurrentSystem].cx-=stress.cx;
  StrainDerivativeTensor[CurrentSystem].cy-=stress.cy;
  StrainDerivativeTensor[CurrentSystem].cz-=stress.cz;


  // compute energies for host-host, host-adsorbate, host-cation, adsorbate-adsorbate, cation-cation, and adsorbate-cation
  // =====================================================================================================================
  // Note: total energy is Fourier minus self minus exclusion energy
  // Note: the net-charge is taken into account

  // Host-host energy
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    UHostHostChargeChargeFourier[CurrentSystem]=energy_framework_framework-energy_framework_self-energy_framework_excluded;
    UHostHostChargeBondDipoleFourier[CurrentSystem]=energy_framework_framework_c_bd-energy_framework_excluded_c_bd;
    UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=energy_framework_framework_bd-energy_framework_self_bd-energy_framework_excluded_bd;
    UHostHostChargeChargeFourier[CurrentSystem]-=UChargeChargeFrameworkRigid[CurrentSystem];
    UHostHostChargeBondDipoleFourier[CurrentSystem]-=UChargeBondDipoleFrameworkRigid[CurrentSystem];
    UHostHostBondDipoleBondDipoleFourier[CurrentSystem]-=UBondDipoleBondDipoleFrameworkRigid[CurrentSystem];
    UHostHostChargeChargeFourier[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeFramework[CurrentSystem]);
  }


  // Host-adsorbate energy
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=2.0*energy_framework_adsorbate;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=2.0*energy_framework_adsorbate_c_bd;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=2.0*energy_framework_adsorbate_bd;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*NetChargeAdsorbates[CurrentSystem];

  // Host-cation energy
  UHostCationChargeChargeFourier[CurrentSystem]=2.0*energy_framework_cation;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=2.0*energy_framework_cation_c_bd;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=2.0*energy_framework_cation_bd;
  UHostCationChargeChargeFourier[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*NetChargeCations[CurrentSystem];

  // Adsorbate-adsorbate energy
  if(!OmitAdsorbateAdsorbateCoulombInteractions)
  {
    UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=energy_adsorbate_adsorbate-energy_adsorbate_self-energy_adsorbate_excluded;
    UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=energy_adsorbate_adsorbate_c_bd-energy_adsorbate_excluded_c_bd;
    UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=energy_adsorbate_adsorbate_bd-energy_adsorbate_self_bd-energy_adsorbate_excluded_bd;
    UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeAdsorbates[CurrentSystem]);
  }

  // Cation-cation energy
  if(!OmitCationCationCoulombInteractions)
  {
    UCationCationChargeChargeFourier[CurrentSystem]=energy_cation_cation-energy_cation_self-energy_cation_excluded;
    UCationCationChargeBondDipoleFourier[CurrentSystem]=energy_cation_cation_c_bd-energy_cation_excluded_c_bd;
    UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=energy_cation_cation_bd-energy_cation_self_bd-energy_cation_excluded_bd;
    UCationCationChargeChargeFourier[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeCations[CurrentSystem]);
  }

  // Adsorbate-cation energy
  if(!OmitAdsorbateCationCoulombInteractions)
  {
    UAdsorbateCationChargeChargeFourier[CurrentSystem]=2.0*energy_adsorbate_cation;
    UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=2.0*energy_adsorbate_cation_c_bd;
    UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=2.0*energy_adsorbate_cation_bd;
    UAdsorbateCationChargeChargeFourier[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeAdsorbates[CurrentSystem]*NetChargeCations[CurrentSystem];
  }

  return 0;
}

/*********************************************************************************************************
 * Name       | EwaldFourierBornTerm                                                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of the charge-charge, charge-bonddipole, and                   *
 *            | bonddipole-bonddipole energies, forces and stress n a truly periodic system.             *                  
 * Parameters | -                                                                                        *
 * Note       | Largely based on the Fiche F.22 of Allen&Tildesly with a few modifications:              *
 *            | 1) exp(-ik.kx)*exp(-ik.ky) computed outside the 'kz' loop                                *
 *            | 2) energy are splitted in framework, adsorbate, cations groups                           *
 *            | 3) net-charge is taken into account                                                      *
 *            | 4) routine computes both charges and bond-dipoles                                        *
 *            | The total structure factors are not stored.                                              *
 *********************************************************************************************************/

int EwaldFourierBornTerm(void)
{
  int i,j,m,f1,ii,jj,kk,A,B;
  int nvec,type_mol,type,typeA,typeB;
  REAL temp,alpha,r,rr,chargeA,charge,DF,DDF;
  VECTOR pos,rk,dr,posA,posB,posA1,posA2,posB1,posB2;
  int nr_molecules,nr_atoms,nr_frameworks,nr_of_excluded_pairs;
  int index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_framework,nr_of_coulombic_sites_adsorbate,nr_of_coulombic_sites_cation;
  COMPLEX sum,sum_framework,sum_adsorbate,sum_cation,temp_sum,temp_sum_bonddipole;
  REAL energy_framework_framework,energy_adsorbate_adsorbate,energy_cation_cation;
  REAL energy_framework_adsorbate,energy_framework_cation,energy_adsorbate_cation;
  REAL energy_framework_self,energy_adsorbate_self,energy_cation_self;
  int nr_of_bonddipoles,nr_of_bonddipole_sites,nr_of_bonddipole_sites_framework,nr_of_bonddipole_sites_adsorbate,nr_of_bonddipole_site_cation;
  COMPLEX sum_bonddipole,sum_bonddipole_framework,sum_bonddipole_adsorbate,sum_bonddipole_cation;
  VECTOR dipole,dipoleA,dipoleB;
  REAL energy_framework_framework_c_bd,energy_framework_adsorbate_c_bd,energy_framework_cation_c_bd;
  REAL energy_adsorbate_adsorbate_c_bd,energy_cation_cation_c_bd,energy_adsorbate_cation_c_bd;
  REAL energy_framework_framework_bd,energy_framework_adsorbate_bd,energy_framework_cation_bd;
  REAL energy_adsorbate_adsorbate_bd,energy_cation_cation_bd,energy_adsorbate_cation_bd;
  REAL energy_framework_self_bd,energy_adsorbate_self_bd,energy_cation_self_bd;
  REAL energy_framework_excluded,energy_adsorbate_excluded,energy_cation_excluded;
  REAL energy_framework_excluded_c_bd,energy_adsorbate_excluded_c_bd,energy_cation_excluded_c_bd;
  REAL energy_framework_excluded_bd,energy_adsorbate_excluded_bd,energy_cation_excluded_bd;
  REAL Bt0,Bt1,Bt2,Bt3,cosA,cosB,cosAB,fac,fac1,fac2,fac3,fac4,fac5,fac6,current_energy;
  REAL DipoleMagnitudeA,length,dot_product,lengthA,lengthB,dipole_magnitudeA,dipole_magnitudeB;
  VECTOR fa1,fa2,fb1,fb2,term,termA,termB;
  ATOM *atom_pointer;
  REAL *kfactor,inverse_lambda_squared,rksq,recip_cutoff,ksqr;
  PAIR pair,pairA,pairB;
  REAL chargeB;
  VECTOR *kvecs;
  REAL_MATRIX3x3 stress,v,S,Theta;
  int considered_charged;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  stress.ax=stress.bx=stress.cx=0.0;
  stress.ay=stress.by=stress.cy=0.0;
  stress.az=stress.bz=stress.cz=0.0;

  // initialize energies
  UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UCationCationChargeChargeFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UHostHostChargeChargeFourier[CurrentSystem]=0.0;

  UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostHostChargeBondDipoleFourier[CurrentSystem]=0.0;

  UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=0.0;

  // initialize energies
  energy_framework_framework=energy_adsorbate_adsorbate=energy_cation_cation=0.0;
  energy_framework_adsorbate=energy_framework_cation=energy_adsorbate_cation=0.0;
  energy_framework_self=energy_adsorbate_self=energy_cation_self=0.0;

  energy_framework_framework_c_bd=energy_adsorbate_adsorbate_c_bd=energy_cation_cation_c_bd=0.0;
  energy_framework_adsorbate_c_bd=energy_framework_cation_c_bd=energy_adsorbate_cation_c_bd=0.0;

  energy_framework_framework_bd=energy_adsorbate_adsorbate_bd=energy_cation_cation_bd=0.0;
  energy_framework_adsorbate_bd=energy_framework_cation_bd=energy_adsorbate_cation_bd=0.0;
  energy_framework_self_bd=energy_adsorbate_self_bd=energy_cation_self_bd=0.0;

  energy_framework_excluded=energy_adsorbate_excluded=energy_cation_excluded=0.0;
  energy_framework_excluded_c_bd=energy_adsorbate_excluded_c_bd=energy_cation_excluded_c_bd=0.0;
  energy_framework_excluded_bd=energy_adsorbate_excluded_bd=energy_cation_excluded_bd=0.0;


  // put charge, bond-dipoles, and positions into appropriate arrays
  // ===============================================================

  nr_of_coulombic_sites=0;
  nr_frameworks=Framework[CurrentSystem].NumberOfFrameworks;
  for(f1=0;f1<nr_frameworks;f1++)
  {
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    for(i=0;i<nr_atoms;i++)
    {
      type=Framework[CurrentSystem].Atoms[f1][i].Type;
      charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged)
      {
        Charge[nr_of_coulombic_sites]=charge;
        if(!(Framework[CurrentSystem].Atoms[f1][i].Fixed.x&&Framework[CurrentSystem].Atoms[f1][i].Fixed.y&&Framework[CurrentSystem].Atoms[f1][i].Fixed.x))
        {
          energy_framework_self+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          AtomVector[nr_of_coulombic_sites]=&Framework[CurrentSystem].Atoms[f1][i].Force;
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_coulombic_sites_framework=nr_of_coulombic_sites;
  nr_molecules=NumberOfAdsorbateMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Adsorbates[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        charge=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        if(considered_charged)
        {
          Charge[nr_of_coulombic_sites]=charge;
          energy_adsorbate_self+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
          if(!(Adsorbates[CurrentSystem][i].Atoms[j].Fixed.x&&Adsorbates[CurrentSystem][i].Atoms[j].Fixed.y&&Adsorbates[CurrentSystem][i].Atoms[j].Fixed.z))
          {
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][i].Atoms[j].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            AtomVector[nr_of_coulombic_sites]=&Adsorbates[CurrentSystem][i].Atoms[j].Force;
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }
  nr_of_coulombic_sites_adsorbate=nr_of_coulombic_sites;
  nr_molecules=NumberOfCationMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Cations[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Cations[CurrentSystem][i].Atoms[j].Type;
        charge=Cations[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        if(considered_charged)
        {
          Charge[nr_of_coulombic_sites]=charge;
          energy_cation_self+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
          if(!(Cations[CurrentSystem][i].Atoms[j].Fixed.x&&Cations[CurrentSystem][i].Atoms[j].Fixed.y&&Cations[CurrentSystem][i].Atoms[j].Fixed.z))
          {
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][i].Atoms[j].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            AtomVector[nr_of_coulombic_sites]=&Cations[CurrentSystem][i].Atoms[j].Force;
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }
  nr_of_coulombic_sites_cation=nr_of_coulombic_sites;

  nr_of_bonddipole_sites=0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      nr_of_bonddipoles=Framework[CurrentSystem].NumberOfBondDipoles[f1];
      for(i=0;i<nr_of_bonddipoles;i++)
      {
        pair=Framework[CurrentSystem].BondDipoles[f1][i];
        atom_pointer=Framework[CurrentSystem].Atoms[f1];
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        dipole=ApplyBoundaryCondition(dipole);
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Framework[CurrentSystem].BondDipoleMagnitude[f1][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        BondDipoleMagnitude[nr_of_bonddipole_sites]=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
        BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));  
        energy_framework_self_bd+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].Force;
        BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].Force;
        nr_of_bonddipole_sites++;
      }
    }
  }
  nr_of_bonddipole_sites_framework=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Adsorbates[CurrentSystem][i].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[j];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));  
      energy_adsorbate_self_bd+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].Force;
      BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].Force;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_adsorbate=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Cations[CurrentSystem][i].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[j];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));  
      energy_cation_self_bd+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].Force;
      BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].Force;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_site_cation=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          // loop over all the framework atoms
          sum_framework.re=0.0;
          sum_framework.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_framework.re+=temp*Eikr[i].re;
            sum_framework.im+=temp*Eikr[i].im;
          }

          // loop over all the adsorbate atoms
          sum_adsorbate.re=0.0;
          sum_adsorbate.im=0.0;
          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_adsorbate.re+=temp*Eikr[i].re;
            sum_adsorbate.im+=temp*Eikr[i].im;
          }

          // loop over all the cation atoms
          sum_cation.re=0.0;
          sum_cation.im=0.0;
          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_cation.re+=temp*Eikr[i].re;
            sum_cation.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_framework.re=0.0;
          sum_bonddipole_framework.im=0.0;

          rk=kvecs[nvec];

          for(i=0;i<nr_of_bonddipole_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_framework.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_framework.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_adsorbate.re=0.0;
          sum_bonddipole_adsorbate.im=0.0;
          for(;i<nr_of_bonddipole_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_adsorbate.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_adsorbate.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_cation.re=0.0;
          sum_bonddipole_cation.im=0.0;
          for(;i<nr_of_bonddipole_site_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_cation.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_cation.im+=temp*Eikr_bd[i].im;
          }

          // add the pre-computed contributions of fixed atoms
          sum_framework.re+=StoreRigidChargeFramework[CurrentSystem][nvec].re;
          sum_framework.im+=StoreRigidChargeFramework[CurrentSystem][nvec].im;
          sum_adsorbate.re+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].re;
          sum_adsorbate.im+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].im;
          sum_cation.re+=StoreRigidChargeCations[CurrentSystem][nvec].re;
          sum_cation.im+=StoreRigidChargeCations[CurrentSystem][nvec].im;

          sum_bonddipole_framework.re+=StoreRigidBondDipolesFramework[CurrentSystem][nvec].re;
          sum_bonddipole_framework.im+=StoreRigidBondDipolesFramework[CurrentSystem][nvec].im;
          sum_bonddipole_adsorbate.re+=StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec].re;
          sum_bonddipole_adsorbate.im+=StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec].im;
          sum_bonddipole_cation.re+=StoreRigidBondDipolesCations[CurrentSystem][nvec].re;
          sum_bonddipole_cation.im+=StoreRigidBondDipolesCations[CurrentSystem][nvec].im;

          // precomputed wavevector dependent pre-factor
          temp=kfactor[nvec];

          // charge-charge, charge-bonddipole, bonddipole-donddipole energies
          energy_framework_framework+=temp*(SQR(sum_framework.re)+SQR(sum_framework.im));
          energy_framework_framework_c_bd+=2.0*temp*(sum_framework.im*sum_bonddipole_framework.re-sum_framework.re*sum_bonddipole_framework.im);
          energy_framework_framework_bd+=temp*(SQR(sum_bonddipole_framework.re)+SQR(sum_bonddipole_framework.im));

          energy_adsorbate_adsorbate+=temp*(SQR(sum_adsorbate.re)+SQR(sum_adsorbate.im));
          energy_adsorbate_adsorbate_c_bd+=2.0*temp*(sum_adsorbate.im*sum_bonddipole_adsorbate.re-sum_adsorbate.re*sum_bonddipole_adsorbate.im);
          energy_adsorbate_adsorbate_bd+=temp*(SQR(sum_bonddipole_adsorbate.re)+SQR(sum_bonddipole_adsorbate.im));

          energy_cation_cation+=temp*(SQR(sum_cation.re)+SQR(sum_cation.im));
          energy_cation_cation_c_bd+=2.0*temp*(sum_cation.im*sum_bonddipole_cation.re-sum_cation.re*sum_bonddipole_cation.im);
          energy_cation_cation_bd+=temp*(SQR(sum_bonddipole_cation.re)+SQR(sum_bonddipole_cation.im));

          energy_framework_adsorbate+=temp*(sum_framework.re*sum_adsorbate.re+sum_framework.im*sum_adsorbate.im);
          energy_framework_adsorbate_c_bd+=temp*(sum_adsorbate.im*sum_bonddipole_framework.re+sum_framework.im*sum_bonddipole_adsorbate.re-
                                                (sum_adsorbate.re*sum_bonddipole_framework.im+sum_framework.re*sum_bonddipole_adsorbate.im));
          energy_framework_adsorbate_bd+=temp*(sum_bonddipole_framework.re*sum_bonddipole_adsorbate.re+sum_bonddipole_framework.im*sum_bonddipole_adsorbate.im);

          energy_framework_cation+=temp*(sum_framework.re*sum_cation.re+sum_framework.im*sum_cation.im);
          energy_framework_cation_c_bd+=temp*(sum_cation.im*sum_bonddipole_framework.re+sum_framework.im*sum_bonddipole_cation.re-
                                             (sum_cation.re*sum_bonddipole_framework.im+sum_framework.re*sum_bonddipole_cation.im));
          energy_framework_cation_bd+=temp*(sum_bonddipole_framework.re*sum_bonddipole_cation.re+sum_bonddipole_framework.im*sum_bonddipole_cation.im);

          energy_adsorbate_cation+=temp*(sum_adsorbate.re*sum_cation.re+sum_adsorbate.im*sum_cation.im);
          energy_adsorbate_cation_c_bd+=temp*(sum_cation.im*sum_bonddipole_adsorbate.re+sum_adsorbate.im*sum_bonddipole_cation.re-
                                             (sum_cation.re*sum_bonddipole_adsorbate.im+sum_adsorbate.re*sum_bonddipole_cation.im));
          energy_adsorbate_cation_bd+=temp*(sum_bonddipole_adsorbate.re*sum_bonddipole_cation.re+sum_bonddipole_adsorbate.im*sum_bonddipole_cation.im);


          // get total sums
          sum.re=sum_framework.re+sum_adsorbate.re+sum_cation.re;
          sum.im=sum_framework.im+sum_adsorbate.im+sum_cation.im;
          sum_bonddipole.re=sum_bonddipole_framework.re+sum_bonddipole_adsorbate.re+sum_bonddipole_cation.re;
          sum_bonddipole.im=sum_bonddipole_framework.im+sum_bonddipole_adsorbate.im+sum_bonddipole_cation.im;

          // forces on atoms from other charges and bond-dipoles
          for(i=0;i<nr_of_coulombic_sites_framework;i++)
          {
            fac=2.0*temp*Charge[i]*((Eikr[i].im*(sum.re-sum_cation.re)-Eikr[i].re*(sum.im-sum_cation.im))+
                (-Eikr[i].re*sum_bonddipole.re-Eikr[i].im*sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitAdsorbateAdsorbateCoulombInteractions)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }
          if(OmitAdsorbateCationCoulombInteractions)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }

          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            fac=2.0*temp*Charge[i]*((Eikr[i].im*temp_sum.re-Eikr[i].re*temp_sum.im)+
                (-Eikr[i].re*temp_sum_bonddipole.re-Eikr[i].im*temp_sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitCationCationCoulombInteractions)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }
          if(OmitAdsorbateCationCoulombInteractions)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }

          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            fac=2.0*temp*Charge[i]*((Eikr[i].im*temp_sum.re-Eikr[i].re*temp_sum.im)+
                (-Eikr[i].re*temp_sum_bonddipole.re-Eikr[i].im*temp_sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }



          // stress tensor
          current_energy=temp*(SQR(sum.re)+SQR(sum.im)+
                               2.0*(sum.im*sum_bonddipole.re-sum.re*sum_bonddipole.im)+
                               SQR(sum_bonddipole.re)+SQR(sum_bonddipole.im));

          // take omitted interactions into account
          if(Framework[CurrentSystem].FrameworkModel!=FLEXIBLE) 
             current_energy-=temp*(SQR(sum_framework.re)+SQR(sum_framework.im)+
                               2.0*(sum_framework.im*sum_bonddipole_framework.re-sum_framework.re*sum_bonddipole_framework.im)+
                               SQR(sum_bonddipole_framework.re)+SQR(sum_bonddipole_framework.im));
          if(OmitAdsorbateAdsorbateCoulombInteractions)
            current_energy-=temp*(SQR(sum_adsorbate.re)+SQR(sum_adsorbate.im)+
                                  2.0*(sum_adsorbate.im*sum_bonddipole_adsorbate.re-sum_adsorbate.re*sum_bonddipole_adsorbate.im)+
                                  SQR(sum_bonddipole_adsorbate.re)+SQR(sum_bonddipole_adsorbate.im));
          if(OmitCationCationCoulombInteractions)
            current_energy-=temp*(SQR(sum_cation.re)+SQR(sum_cation.im)+
                                  2.0*(sum_cation.im*sum_bonddipole_cation.re-sum_cation.re*sum_bonddipole_cation.im)+
                                  SQR(sum_bonddipole_cation.re)+SQR(sum_bonddipole_cation.im));
          if(OmitAdsorbateCationCoulombInteractions)
            current_energy-=2.0*temp*(sum_adsorbate.re*sum_cation.re+sum_adsorbate.im*sum_cation.im+
                                      sum_adsorbate.im*sum_bonddipole_cation.re-sum_adsorbate.re*sum_bonddipole_cation.im+
                                      sum_cation.im*sum_bonddipole_adsorbate.re-sum_cation.re*sum_bonddipole_adsorbate.im+
                                      sum_bonddipole_adsorbate.re*sum_bonddipole_cation.re+sum_bonddipole_adsorbate.im*sum_bonddipole_cation.im);

          fac=2.0*(1.0/(SQR(rk.x)+SQR(rk.y)+SQR(rk.z))+(0.25/SQR(alpha)))*current_energy;
          S.ax=(current_energy-fac*rk.x*rk.x);
          S.ay=-fac*rk.x*rk.y;
          S.az=-fac*rk.x*rk.z;

          S.bx=-fac*rk.y*rk.x;
          S.by=(current_energy-fac*rk.y*rk.y);
          S.bz=-fac*rk.y*rk.z;

          S.cx=-fac*rk.z*rk.x;
          S.cy=-fac*rk.z*rk.y;
          S.cz=(current_energy-fac*rk.z*rk.z);

          S.ax=S.ay=S.az=0.0;
          S.bx=S.by=S.bz=0.0;
          S.cx=S.cy=S.cz=0.0;



          stress.ax+=current_energy-fac*rk.x*rk.x;
          stress.ay+=-fac*rk.x*rk.y;
          stress.az+=-fac*rk.x*rk.z;

          stress.bx+=-fac*rk.y*rk.x;
          stress.by+=current_energy-fac*rk.y*rk.y;
          stress.bz+=-fac*rk.y*rk.z;

          stress.cx+=-fac*rk.z*rk.x;
          stress.cy+=-fac*rk.z*rk.y;
          stress.cz+=current_energy-fac*rk.z*rk.z;

          rksq=SQR(rk.x)+SQR(rk.y)+SQR(rk.z);
          inverse_lambda_squared=0.25/(SQR(Alpha[CurrentSystem]))+1.0/rksq;

/*
(4.0*M_PI*COULOMBIC_CONVERSION_FACTOR/Volume[CurrentSystem])*energy_term*

COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/volume)*exp((-0.25/SQR(alpha))*rksqr)/rksqr
*/
          Theta.ax=1.0-2.0*rk.x*rk.x*inverse_lambda_squared;
          Theta.ay=-2.0*rk.x*rk.y*inverse_lambda_squared;
          Theta.az=-2.0*rk.x*rk.z*inverse_lambda_squared;

          Theta.bx=-2.0*rk.y*rk.x*inverse_lambda_squared;
          Theta.by=1.0-2.0*rk.y*rk.y*inverse_lambda_squared;
          Theta.bz=-2.0*rk.y*rk.z*inverse_lambda_squared;

          Theta.cx=-2.0*rk.z*rk.x*inverse_lambda_squared;
          Theta.cy=-2.0*rk.z*rk.y*inverse_lambda_squared;
          Theta.cz=1.0-2.0*rk.z*rk.z*inverse_lambda_squared;

/*
          BornTerm[CurrentSystem].xxxx+=current_energy*
            (Theta.ax*Theta.ax+2.0+4.0*rk.x*rk.x*rk.x*rk.x/(SQR(rksq))-8.0*rk.x*rk.x*inverse_lambda_squared-2.0*Theta.ax);
*/

          BornTerm[CurrentSystem].xxxx+=current_energy*
                 ((1.0-2.0*rk.x*rk.x*inverse_lambda_squared)*(1.0-2.0*rk.x*rk.x*inverse_lambda_squared)+(1.0+1.0)+
                 4.0*rk.x*rk.x*rk.x*rk.x/(SQR(rksq))
                 -2.0*(rk.x*rk.x+rk.x*rk.x+rk.x*rk.x+rk.x*rk.x)*inverse_lambda_squared-2.0*Theta.ax);
          BornTerm[CurrentSystem].xxyy+=current_energy*
                 ((1.0-2.0*rk.x*rk.x*inverse_lambda_squared)*(1.0-2.0*rk.y*rk.y*inverse_lambda_squared)+
                 4.0*rk.x*rk.x*rk.y*rk.y/(SQR(rksq)));
          BornTerm[CurrentSystem].xxzz+=current_energy*
                 ((1.0-2.0*rk.x*rk.x*inverse_lambda_squared)*(1.0-2.0*rk.z*rk.z*inverse_lambda_squared)+
                 4.0*rk.x*rk.x*rk.z*rk.z/(SQR(rksq)));
          BornTerm[CurrentSystem].xxyz+=current_energy*
                 ((1.0-2.0*rk.x*rk.x*inverse_lambda_squared)*(-2.0*rk.y*rk.z*inverse_lambda_squared)+
                 4.0*rk.x*rk.x*rk.y*rk.z/(SQR(rksq)));
          BornTerm[CurrentSystem].xxzx+=current_energy*
                 ((1.0-2.0*rk.x*rk.x*inverse_lambda_squared)*(-2.0*rk.z*rk.x*inverse_lambda_squared)+
                 4.0*rk.x*rk.x*rk.z*rk.x/(SQR(rksq))
                 -2.0*(rk.x*rk.z+rk.x*rk.z)*inverse_lambda_squared-Theta.az);
          BornTerm[CurrentSystem].xxxy+=current_energy*
                 ((1.0-2.0*rk.x*rk.x*inverse_lambda_squared)*(-2.0*rk.x*rk.y*inverse_lambda_squared)+
                 4.0*rk.x*rk.x*rk.x*rk.y/(SQR(rksq))
                 -2.0*(rk.x*rk.y+rk.x*rk.y)*inverse_lambda_squared-Theta.ay);

          BornTerm[CurrentSystem].yyyy+=current_energy*
                 ((1.0-2.0*rk.y*rk.y*inverse_lambda_squared)*(1.0-2.0*rk.y*rk.y*inverse_lambda_squared)+(1.0+1.0)+
                 4.0*rk.y*rk.y*rk.y*rk.y/(SQR(rksq))
                 -2.0*(rk.y*rk.y+rk.y*rk.y+rk.y*rk.y+rk.y*rk.y)*inverse_lambda_squared-2.0*Theta.by);
          BornTerm[CurrentSystem].yyzz+=current_energy*
                 ((1.0-2.0*rk.y*rk.y*inverse_lambda_squared)*(1.0-2.0*rk.z*rk.z*inverse_lambda_squared)+
                 4.0*rk.y*rk.y*rk.z*rk.z/(SQR(rksq)));
          BornTerm[CurrentSystem].yyyz+=current_energy*
                 ((1.0-2.0*rk.y*rk.y*inverse_lambda_squared)*(-2.0*rk.y*rk.z*inverse_lambda_squared)+
                 4.0*rk.y*rk.y*rk.y*rk.z/(SQR(rksq))
                 -2.0*(rk.y*rk.z+rk.y*rk.z)*inverse_lambda_squared-Theta.bz);
          BornTerm[CurrentSystem].yyzx+=current_energy*
                 ((1.0-2.0*rk.y*rk.y*inverse_lambda_squared)*(-2.0*rk.z*rk.x*inverse_lambda_squared)+
                 4.0*rk.y*rk.y*rk.z*rk.x/(SQR(rksq)));
          BornTerm[CurrentSystem].yyxy+=current_energy*
                 ((1.0-2.0*rk.y*rk.y*inverse_lambda_squared)*(-2.0*rk.x*rk.y*inverse_lambda_squared)+
                 4.0*rk.y*rk.y*rk.x*rk.y/(SQR(rksq))
                 -2.0*(rk.y*rk.x+rk.y*rk.x)*inverse_lambda_squared-Theta.bx);

          BornTerm[CurrentSystem].zzzz+=current_energy*
                 ((1.0-2.0*rk.z*rk.z*inverse_lambda_squared)*(1.0-2.0*rk.z*rk.z*inverse_lambda_squared)+(1.0+1.0)+
                 4.0*rk.z*rk.z*rk.z*rk.z/(SQR(rksq))
                 -2.0*(rk.z*rk.z+rk.z*rk.z+rk.z*rk.z+rk.z*rk.z)*inverse_lambda_squared-2.0*Theta.cz);
          BornTerm[CurrentSystem].zzyz+=current_energy*
                 ((1.0-2.0*rk.z*rk.z*inverse_lambda_squared)*(-2.0*rk.y*rk.z*inverse_lambda_squared)+
                 4.0*rk.z*rk.z*rk.y*rk.z/(SQR(rksq))
                 -2.0*(rk.z*rk.y+rk.z*rk.y)*inverse_lambda_squared-Theta.cy);
          BornTerm[CurrentSystem].zzzx+=current_energy*
                 ((1.0-2.0*rk.z*rk.z*inverse_lambda_squared)*(-2.0*rk.z*rk.x*inverse_lambda_squared)+
                 4.0*rk.z*rk.z*rk.z*rk.x/(SQR(rksq))
                 -2.0*(rk.z*rk.x+rk.z*rk.x)*inverse_lambda_squared-Theta.cx);
          BornTerm[CurrentSystem].zzxy+=current_energy*
                 ((1.0-2.0*rk.z*rk.z*inverse_lambda_squared)*(-2.0*rk.x*rk.y*inverse_lambda_squared)+
                 4.0*rk.z*rk.z*rk.x*rk.y/(SQR(rksq)));

          BornTerm[CurrentSystem].yzyz+=current_energy*
                 ((-2.0*rk.y*rk.z*inverse_lambda_squared)*(-2.0*rk.y*rk.z*inverse_lambda_squared)+1.0+
                 4.0*rk.y*rk.z*rk.y*rk.z/(SQR(rksq))
                 -2.0*(rk.z*rk.z+rk.y*rk.y)*inverse_lambda_squared-0.5*(Theta.by+Theta.cz));
          BornTerm[CurrentSystem].yzzx+=current_energy*
                 ((-2.0*rk.y*rk.z*inverse_lambda_squared)*(-2.0*rk.z*rk.x*inverse_lambda_squared)+
                 4.0*rk.y*rk.z*rk.z*rk.x/(SQR(rksq))
                 -2.0*(rk.y*rk.x)*inverse_lambda_squared-0.5*Theta.bx);
          BornTerm[CurrentSystem].yzxy+=current_energy*
                 ((-2.0*rk.y*rk.z*inverse_lambda_squared)*(-2.0*rk.x*rk.y*inverse_lambda_squared)+
                 4.0*rk.y*rk.z*rk.x*rk.y/(SQR(rksq))
                 -2.0*(rk.z*rk.x)*inverse_lambda_squared-0.5*Theta.az);

          BornTerm[CurrentSystem].zxzx+=current_energy*
                 ((-2.0*rk.z*rk.x*inverse_lambda_squared)*(-2.0*rk.z*rk.x*inverse_lambda_squared)+1.0+
                 4.0*rk.z*rk.x*rk.z*rk.x/(SQR(rksq))
                 -2.0*(rk.x*rk.x+rk.z*rk.z)*inverse_lambda_squared-0.5*(Theta.cz+Theta.ax));
          BornTerm[CurrentSystem].zxxy+=current_energy*
                 ((-2.0*rk.z*rk.x*inverse_lambda_squared)*(-2.0*rk.x*rk.y*inverse_lambda_squared)+
                 4.0*rk.z*rk.x*rk.x*rk.y/(SQR(rksq))
                 -2.0*(rk.z*rk.y)*inverse_lambda_squared-0.5*Theta.cy);

          BornTerm[CurrentSystem].xyxy+=current_energy*
                 ((-2.0*rk.x*rk.y*inverse_lambda_squared)*(-2.0*rk.x*rk.y*inverse_lambda_squared)+(1.0)+
                 4.0*rk.x*rk.y*rk.x*rk.y/(SQR(rksq))
                 -2.0*(rk.y*rk.y+rk.x*rk.x)*inverse_lambda_squared-0.5*(Theta.ax+Theta.by));


          // forces on bond-dipoles from charges and other bond-dipoles
          for(i=0;i<nr_of_bonddipole_sites_framework;i++)
          {
            dipole=DipoleVector[i];
            DipoleMagnitudeA=BondDipoleMagnitude[i];
            length=BondLength[i];
            dot_product=(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);

            // charge-bonddipole contribution
            fac1=temp*sum.re*dot_product*Eikr_bd[i].re;
            fac2=2.0*temp*sum.re*Eikr_bd[i].im*DipoleMagnitudeA/length;
            fac3=2.0*temp*sum.re*Eikr_bd[i].im*dot_product/(length*DipoleMagnitudeA);

            fac4=temp*sum.im*dot_product*Eikr_bd[i].im;
            fac5=2.0*temp*sum.im*Eikr_bd[i].re*DipoleMagnitudeA/length;
            fac6=2.0*temp*sum.im*Eikr_bd[i].re*dot_product/(length*DipoleMagnitudeA);

            fa1.x=(fac1-fac2)*rk.x+fac3*dipole.x+(fac4+fac5)*rk.x-fac6*dipole.x;
            fa1.y=(fac1-fac2)*rk.y+fac3*dipole.y+(fac4+fac5)*rk.y-fac6*dipole.y;
            fa1.z=(fac1-fac2)*rk.z+fac3*dipole.z+(fac4+fac5)*rk.z-fac6*dipole.z;

            fa2.x=(fac1+fac2)*rk.x-fac3*dipole.x+(fac4-fac5)*rk.x+fac6*dipole.x;
            fa2.y=(fac1+fac2)*rk.y-fac3*dipole.y+(fac4-fac5)*rk.y+fac6*dipole.y;
            fa2.z=(fac1+fac2)*rk.z-fac3*dipole.z+(fac4-fac5)*rk.z+fac6*dipole.z;

            // bonddipole-bonddipole contribution
            fac1=temp*dot_product*(sum_bonddipole.re*Eikr_bd[i].im-sum_bonddipole.im*Eikr_bd[i].re);
            fac2=2.0*temp*(sum_bonddipole.re*Eikr_bd[i].re+sum_bonddipole.im*Eikr_bd[i].im)*DipoleMagnitudeA/length;
            fac3=2.0*temp*(sum_bonddipole.re*Eikr_bd[i].re+sum_bonddipole.im*Eikr_bd[i].im)*dot_product/
                 (length*DipoleMagnitudeA);

            fa1.x+=(fac1+fac2)*rk.x-fac3*dipole.x;
            fa1.y+=(fac1+fac2)*rk.y-fac3*dipole.y;
            fa1.z+=(fac1+fac2)*rk.z-fac3*dipole.z;

            fa2.x+=(fac1-fac2)*rk.x+fac3*dipole.x;
            fa2.y+=(fac1-fac2)*rk.y+fac3*dipole.y;
            fa2.z+=(fac1-fac2)*rk.z+fac3*dipole.z;

            BondDipoleForcesA[i]->x+=fa1.x;
            BondDipoleForcesA[i]->y+=fa1.y;
            BondDipoleForcesA[i]->z+=fa1.z;

            BondDipoleForcesB[i]->x+=fa2.x;
            BondDipoleForcesB[i]->y+=fa2.y;
            BondDipoleForcesB[i]->z+=fa2.z;

            // bond-dipole contribution to the stress
            fac=temp*(sum.im*Eikr_bd[i].re-sum.re*Eikr_bd[i].im
                      +sum_bonddipole.im*Eikr_bd[i].im+sum_bonddipole.re*Eikr_bd[i].re);
            stress.ax+=2.0*fac*dipole.x*rk.x;
            stress.ay+=fac*(dipole.x*rk.y+dipole.y*rk.x);
            stress.az+=fac*(dipole.x*rk.z+dipole.z*rk.x);
            stress.bx+=fac*(dipole.y*rk.x+dipole.x*rk.y);
            stress.by+=2.0*fac*dipole.y*rk.y;
            stress.bz+=fac*(dipole.y*rk.z+dipole.z*rk.y);
            stress.cx+=fac*(dipole.z*rk.x+dipole.x*rk.z);
            stress.cy+=fac*(dipole.z*rk.y+dipole.y*rk.z);
            stress.cz+=2.0*fac*dipole.z*rk.z;

            // convert forces on atoms to 'molecular stress'
            // it is equivalent to: v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x); etc.
            fac=0.5*(length/DipoleMagnitudeA);
            stress.ax+=fac*(fa2.x-fa1.x)*dipole.x;
            stress.bx+=fac*(fa2.y-fa1.y)*dipole.x;
            stress.cx+=fac*(fa2.z-fa1.z)*dipole.x;
            stress.ay+=fac*(fa2.x-fa1.x)*dipole.y;
            stress.by+=fac*(fa2.y-fa1.y)*dipole.y;
            stress.cy+=fac*(fa2.z-fa1.z)*dipole.y;
            stress.az+=fac*(fa2.x-fa1.x)*dipole.z;
            stress.bz+=fac*(fa2.y-fa1.y)*dipole.z;
            stress.cz+=fac*(fa2.z-fa1.z)*dipole.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitAdsorbateAdsorbateCoulombInteractions)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }
          if(OmitAdsorbateCationCoulombInteractions)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }

          for(;i<nr_of_bonddipole_sites_adsorbate;i++)
          {
            dipole=DipoleVector[i];
            DipoleMagnitudeA=BondDipoleMagnitude[i];
            length=BondLength[i];
            dot_product=(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);

            // charge-bonddipole contribution
            fac1=temp*temp_sum.re*dot_product*Eikr_bd[i].re;
            fac2=2.0*temp*temp_sum.re*Eikr_bd[i].im*DipoleMagnitudeA/length;
            fac3=2.0*temp*temp_sum.re*Eikr_bd[i].im*dot_product/(length*DipoleMagnitudeA);

            fac4=temp*temp_sum.im*dot_product*Eikr_bd[i].im;
            fac5=2.0*temp*temp_sum.im*Eikr_bd[i].re*DipoleMagnitudeA/length;
            fac6=2.0*temp*temp_sum.im*Eikr_bd[i].re*dot_product/(length*DipoleMagnitudeA);

            fa1.x=(fac1-fac2)*rk.x+fac3*dipole.x+(fac4+fac5)*rk.x-fac6*dipole.x;
            fa1.y=(fac1-fac2)*rk.y+fac3*dipole.y+(fac4+fac5)*rk.y-fac6*dipole.y;
            fa1.z=(fac1-fac2)*rk.z+fac3*dipole.z+(fac4+fac5)*rk.z-fac6*dipole.z;

            fa2.x=(fac1+fac2)*rk.x-fac3*dipole.x+(fac4-fac5)*rk.x+fac6*dipole.x;
            fa2.y=(fac1+fac2)*rk.y-fac3*dipole.y+(fac4-fac5)*rk.y+fac6*dipole.y;
            fa2.z=(fac1+fac2)*rk.z-fac3*dipole.z+(fac4-fac5)*rk.z+fac6*dipole.z;

            // bonddipole-bonddipole contribution
            fac1=temp*dot_product*(temp_sum_bonddipole.re*Eikr_bd[i].im-temp_sum_bonddipole.im*Eikr_bd[i].re);
            fac2=2.0*temp*(temp_sum_bonddipole.re*Eikr_bd[i].re+temp_sum_bonddipole.im*Eikr_bd[i].im)*DipoleMagnitudeA/length;
            fac3=2.0*temp*(temp_sum_bonddipole.re*Eikr_bd[i].re+temp_sum_bonddipole.im*Eikr_bd[i].im)*dot_product/
                 (length*DipoleMagnitudeA);

            fa1.x+=(fac1+fac2)*rk.x-fac3*dipole.x;
            fa1.y+=(fac1+fac2)*rk.y-fac3*dipole.y;
            fa1.z+=(fac1+fac2)*rk.z-fac3*dipole.z;

            fa2.x+=(fac1-fac2)*rk.x+fac3*dipole.x;
            fa2.y+=(fac1-fac2)*rk.y+fac3*dipole.y;
            fa2.z+=(fac1-fac2)*rk.z+fac3*dipole.z;

            BondDipoleForcesA[i]->x+=fa1.x;
            BondDipoleForcesA[i]->y+=fa1.y;
            BondDipoleForcesA[i]->z+=fa1.z;

            BondDipoleForcesB[i]->x+=fa2.x;
            BondDipoleForcesB[i]->y+=fa2.y;
            BondDipoleForcesB[i]->z+=fa2.z;

            // bond-dipole contribution to the stress
            fac=temp*(temp_sum.im*Eikr_bd[i].re-temp_sum.re*Eikr_bd[i].im
                      +temp_sum_bonddipole.im*Eikr_bd[i].im+temp_sum_bonddipole.re*Eikr_bd[i].re);
            stress.ax+=2.0*fac*dipole.x*rk.x;
            stress.ay+=fac*(dipole.x*rk.y+dipole.y*rk.x);
            stress.az+=fac*(dipole.x*rk.z+dipole.z*rk.x);
            stress.bx+=fac*(dipole.y*rk.x+dipole.x*rk.y);
            stress.by+=2.0*fac*dipole.y*rk.y;
            stress.bz+=fac*(dipole.y*rk.z+dipole.z*rk.y);
            stress.cx+=fac*(dipole.z*rk.x+dipole.x*rk.z);
            stress.cy+=fac*(dipole.z*rk.y+dipole.y*rk.z);
            stress.cz+=2.0*fac*dipole.z*rk.z;

            // convert forces on atoms to 'molecular stress'
            // it is equivalent to: v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x); etc.
            fac=0.5*(length/DipoleMagnitudeA);
            stress.ax+=fac*(fa2.x-fa1.x)*dipole.x;
            stress.bx+=fac*(fa2.y-fa1.y)*dipole.x;
            stress.cx+=fac*(fa2.z-fa1.z)*dipole.x;
            stress.ay+=fac*(fa2.x-fa1.x)*dipole.y;
            stress.by+=fac*(fa2.y-fa1.y)*dipole.y;
            stress.cy+=fac*(fa2.z-fa1.z)*dipole.y;
            stress.az+=fac*(fa2.x-fa1.x)*dipole.z;
            stress.bz+=fac*(fa2.y-fa1.y)*dipole.z;
            stress.cz+=fac*(fa2.z-fa1.z)*dipole.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitCationCationCoulombInteractions)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }
          if(OmitAdsorbateCationCoulombInteractions)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }

          for(;i<nr_of_bonddipole_site_cation;i++)
          {
            dipole=DipoleVector[i];
            DipoleMagnitudeA=BondDipoleMagnitude[i];
            length=BondLength[i];
            dot_product=(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);

            // charge-bonddipole contribution
            fac1=temp*temp_sum.re*dot_product*Eikr_bd[i].re;
            fac2=2.0*temp*temp_sum.re*Eikr_bd[i].im*DipoleMagnitudeA/length;
            fac3=2.0*temp*temp_sum.re*Eikr_bd[i].im*dot_product/(length*DipoleMagnitudeA);

            fac4=temp*temp_sum.im*dot_product*Eikr_bd[i].im;
            fac5=2.0*temp*temp_sum.im*Eikr_bd[i].re*DipoleMagnitudeA/length;
            fac6=2.0*temp*temp_sum.im*Eikr_bd[i].re*dot_product/(length*DipoleMagnitudeA);

            fa1.x=(fac1-fac2)*rk.x+fac3*dipole.x+(fac4+fac5)*rk.x-fac6*dipole.x;
            fa1.y=(fac1-fac2)*rk.y+fac3*dipole.y+(fac4+fac5)*rk.y-fac6*dipole.y;
            fa1.z=(fac1-fac2)*rk.z+fac3*dipole.z+(fac4+fac5)*rk.z-fac6*dipole.z;

            fa2.x=(fac1+fac2)*rk.x-fac3*dipole.x+(fac4-fac5)*rk.x+fac6*dipole.x;
            fa2.y=(fac1+fac2)*rk.y-fac3*dipole.y+(fac4-fac5)*rk.y+fac6*dipole.y;
            fa2.z=(fac1+fac2)*rk.z-fac3*dipole.z+(fac4-fac5)*rk.z+fac6*dipole.z;

            // bonddipole-bonddipole contribution
            fac1=temp*dot_product*(temp_sum_bonddipole.re*Eikr_bd[i].im-temp_sum_bonddipole.im*Eikr_bd[i].re);
            fac2=2.0*temp*(temp_sum_bonddipole.re*Eikr_bd[i].re+temp_sum_bonddipole.im*Eikr_bd[i].im)*DipoleMagnitudeA/length;
            fac3=2.0*temp*(temp_sum_bonddipole.re*Eikr_bd[i].re+temp_sum_bonddipole.im*Eikr_bd[i].im)*dot_product/
                 (length*DipoleMagnitudeA);

            fa1.x+=(fac1+fac2)*rk.x-fac3*dipole.x;
            fa1.y+=(fac1+fac2)*rk.y-fac3*dipole.y;
            fa1.z+=(fac1+fac2)*rk.z-fac3*dipole.z;

            fa2.x+=(fac1-fac2)*rk.x+fac3*dipole.x;
            fa2.y+=(fac1-fac2)*rk.y+fac3*dipole.y;
            fa2.z+=(fac1-fac2)*rk.z+fac3*dipole.z;

            BondDipoleForcesA[i]->x+=fa1.x;
            BondDipoleForcesA[i]->y+=fa1.y;
            BondDipoleForcesA[i]->z+=fa1.z;

            BondDipoleForcesB[i]->x+=fa2.x;
            BondDipoleForcesB[i]->y+=fa2.y;
            BondDipoleForcesB[i]->z+=fa2.z;

            // bond-dipole contribution to the stress
            fac=temp*(temp_sum.im*Eikr_bd[i].re-temp_sum.re*Eikr_bd[i].im
                      +temp_sum_bonddipole.im*Eikr_bd[i].im+temp_sum_bonddipole.re*Eikr_bd[i].re);
            stress.ax+=2.0*fac*dipole.x*rk.x;
            stress.ay+=fac*(dipole.x*rk.y+dipole.y*rk.x);
            stress.az+=fac*(dipole.x*rk.z+dipole.z*rk.x);
            stress.bx+=fac*(dipole.y*rk.x+dipole.x*rk.y);
            stress.by+=2.0*fac*dipole.y*rk.y;
            stress.bz+=fac*(dipole.y*rk.z+dipole.z*rk.y);
            stress.cx+=fac*(dipole.z*rk.x+dipole.x*rk.z);
            stress.cy+=fac*(dipole.z*rk.y+dipole.y*rk.z);
            stress.cz+=2.0*fac*dipole.z*rk.z;

            // convert forces on atoms to 'molecular stress'
            // it is equivalent to: v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x); etc.
            fac=0.5*(length/DipoleMagnitudeA);
            stress.ax+=fac*(fa2.x-fa1.x)*dipole.x;
            stress.bx+=fac*(fa2.y-fa1.y)*dipole.x;
            stress.cx+=fac*(fa2.z-fa1.z)*dipole.x;
            stress.ay+=fac*(fa2.x-fa1.x)*dipole.y;
            stress.by+=fac*(fa2.y-fa1.y)*dipole.y;
            stress.cy+=fac*(fa2.z-fa1.z)*dipole.y;
            stress.az+=fac*(fa2.x-fa1.x)*dipole.z;
            stress.bz+=fac*(fa2.y-fa1.y)*dipole.z;
            stress.cz+=fac*(fa2.z-fa1.z)*dipole.z;
          }

          // next wave-vector
          nvec++;
        }
      }
    }
  }

  // Exclusion pairs for the framework
  // =================================
  // The list of exclusion pairs is pre-compute at the start of the simulation based on e.g. 1-2, 1-3 exclusion definitions.
  // Note: there is no cutoff used here

  energy_framework_excluded=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        atom_pointer=Framework[CurrentSystem].Atoms[f1];
        pair=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        chargeA=atom_pointer[pair.A].Charge;
        chargeB=atom_pointer[pair.B].Charge;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;

        if(r>1e-4)
        {
          energy_framework_excluded-=chargeA*chargeB*Bt0;

          DF=-chargeA*chargeB*Bt1;
          DDF=chargeA*chargeB*Bt2;

          temp=chargeA*chargeB*Bt1;

          atom_pointer[pair.A].Force.x+=temp*dr.x;
          atom_pointer[pair.A].Force.y+=temp*dr.y;
          atom_pointer[pair.A].Force.z+=temp*dr.z;

          atom_pointer[pair.B].Force.x-=temp*dr.x;
          atom_pointer[pair.B].Force.y-=temp*dr.y;
          atom_pointer[pair.B].Force.z-=temp*dr.z;

          stress.ax+=temp*dr.x*dr.x;
          stress.bx+=temp*dr.y*dr.x;
          stress.cx+=temp*dr.z*dr.x;

          stress.ay+=temp*dr.x*dr.y;
          stress.by+=temp*dr.y*dr.y;
          stress.cy+=temp*dr.z*dr.y;

          stress.az+=temp*dr.x*dr.z;
          stress.bz+=temp*dr.y*dr.z;
          stress.cz+=temp*dr.z*dr.z;

          // add contribution to the born term
          AddContributionToBornTerm(DDF,DF,dr);
        }
        else // if r->0 compute limiting value to avoid divergence when shell overlaps with core in core-shell models
        {
          energy_framework_excluded+=chargeA*chargeB*alpha*(2.0/sqrt(M_PI));
          DF=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*4.0*CUBE(Alpha[CurrentSystem])/(3.0*sqrt(M_PI));
          DDF=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*8.0*CUBE(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])/(5.0*sqrt(M_PI));
          AddContributionToBornTerm(DDF,DF,dr);
        }
      }
    }
  }

  energy_framework_excluded_c_bd=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeBondDipole[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        // the index of the bonddipole is the second index 'B'
        A=Framework[CurrentSystem].ExcludedIntraChargeBondDipole[f1][i].A;
        B=Framework[CurrentSystem].ExcludedIntraChargeBondDipole[f1][i].B;

        atom_pointer=Framework[CurrentSystem].Atoms[f1];

        typeA=atom_pointer[A].Type;
        chargeA=atom_pointer[A].Charge;
        posA=atom_pointer[A].Position;

        pairB=Framework[CurrentSystem].BondDipoles[f1][B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        dipole_magnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][B];
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posB.x-posA.x;
        dr.y=posB.y-posA.y;
        dr.z=posB.z-posA.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;

        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_framework_excluded_c_bd+=Bt1*chargeA*cosB;

        fa1.x=chargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
        fa1.y=chargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
        fa1.z=chargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

        atom_pointer[A].Force.x+=fa1.x;
        atom_pointer[A].Force.y+=fa1.y;
        atom_pointer[A].Force.z+=fa1.z;

        fac1=chargeA*Bt1*dipole_magnitudeB/lengthB;
        fac2=chargeA*Bt1*cosB/(dipole_magnitudeB*lengthB);

        fb1.x=0.5*fa1.x+fac1*dr.x-fac2*dipoleB.x;
        fb1.y=0.5*fa1.y+fac1*dr.y-fac2*dipoleB.y;
        fb1.z=0.5*fa1.z+fac1*dr.z-fac2*dipoleB.z;

        fb2.x=0.5*fa1.x-fac1*dr.x+fac2*dipoleB.x;
        fb2.y=0.5*fa1.y-fac1*dr.y+fac2*dipoleB.y;
        fb2.z=0.5*fa1.z-fac1*dr.z+fac2*dipoleB.z;

        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=fa1.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=fa1.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=fa1.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=fa1.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=fa1.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=fa1.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=fa1.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=fa1.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=fa1.z*dr.z+v.cz;
      }
    }
  }

  energy_framework_excluded_bd=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraBondDipoleBondDipole[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Framework[CurrentSystem].ExcludedIntraBondDipoleBondDipole[f1][i].A;
        B=Framework[CurrentSystem].ExcludedIntraBondDipoleBondDipole[f1][i].B;

        atom_pointer=Framework[CurrentSystem].Atoms[f1];

        pairA=Framework[CurrentSystem].BondDipoles[f1][A];
        posA1=atom_pointer[pairA.A].Position;
        posA2=atom_pointer[pairA.B].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        dipoleA=ApplyBoundaryCondition(dipoleA);
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        dipole_magnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][A];
        lengthA=sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
        temp=dipole_magnitudeA/lengthA;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        pairB=Framework[CurrentSystem].BondDipoles[f1][B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        dipole_magnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][B];
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;
        temp*=2.0*SQR(alpha);
        Bt3=temp+(5.0/rr)*Bt2;

        cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
        cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_framework_excluded_bd-=(Bt1*cosAB-Bt2*cosA*cosB);

        term.x=(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
        term.y=(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
        term.z=(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

        termA.x=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.x/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.x-Bt1*dipoleB.x)*dipole_magnitudeA/lengthA;
        termA.y=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.y/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.y-Bt1*dipoleB.y)*dipole_magnitudeA/lengthA;
        termA.z=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.z/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.z-Bt1*dipoleB.z)*dipole_magnitudeA/lengthA;

        termB.x=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.x/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.x-Bt1*dipoleA.x)*dipole_magnitudeB/lengthB;
        termB.y=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.y/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.y-Bt1*dipoleA.y)*dipole_magnitudeB/lengthB;
        termB.z=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.z/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.z-Bt1*dipoleA.z)*dipole_magnitudeB/lengthB;

        fa1.x=0.5*term.x+termA.x;
        fa1.y=0.5*term.y+termA.y;
        fa1.z=0.5*term.z+termA.z;
        fa2.x=0.5*term.x-termA.x;
        fa2.y=0.5*term.y-termA.y;
        fa2.z=0.5*term.z-termA.z;

        fb1.x=-0.5*term.x+termB.x;
        fb1.y=-0.5*term.y+termB.y;
        fb1.z=-0.5*term.z+termB.z;
        fb2.x=-0.5*term.x-termB.x;
        fb2.y=-0.5*term.y-termB.y;
        fb2.z=-0.5*term.z-termB.z;

        atom_pointer[pairA.A].Force.x-=fa1.x;
        atom_pointer[pairA.A].Force.y-=fa1.y;
        atom_pointer[pairA.A].Force.z-=fa1.z;

        atom_pointer[pairA.B].Force.x-=fa2.x;
        atom_pointer[pairA.B].Force.y-=fa2.y;
        atom_pointer[pairA.B].Force.z-=fa2.z;

        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
      }
    }
  }

  // Exclusion pairs for the adsorbates
  // ==================================

  // if adsorbate-adsorbate interactions are omitted, then the intra-molecular interactions are already fully elimated
  if(!OmitAdsorbateAdsorbateCoulombInteractions)
  {
    energy_adsorbate_excluded=0.0;
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      atom_pointer=Adsorbates[CurrentSystem][m].Atoms;
      type=Adsorbates[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        pair=Components[type].ExcludedIntraChargeCharge[i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        chargeA=atom_pointer[pair.A].Charge;
        chargeB=atom_pointer[pair.B].Charge;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;

        energy_adsorbate_excluded-=chargeA*chargeB*Bt0;

        temp=chargeA*chargeB*Bt1;

        DF=-chargeA*chargeB*Bt1;
        DDF=chargeA*chargeB*Bt2;

        atom_pointer[pair.A].Force.x+=temp*dr.x;
        atom_pointer[pair.A].Force.y+=temp*dr.y;
        atom_pointer[pair.A].Force.z+=temp*dr.z;

        atom_pointer[pair.B].Force.x-=temp*dr.x;
        atom_pointer[pair.B].Force.y-=temp*dr.y;
        atom_pointer[pair.B].Force.z-=temp*dr.z;

        stress.ax+=temp*dr.x*dr.x;
        stress.bx+=temp*dr.y*dr.x;
        stress.cx+=temp*dr.z*dr.x;

        stress.ay+=temp*dr.x*dr.y;
        stress.by+=temp*dr.y*dr.y;
        stress.cy+=temp*dr.z*dr.y;

        stress.az+=temp*dr.x*dr.z;
        stress.bz+=temp*dr.y*dr.z;
        stress.cz+=temp*dr.z*dr.z;

        AddContributionToBornTerm(DDF,DF,dr);
      }
    }

    energy_adsorbate_excluded_c_bd=0.0;
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      type=Adsorbates[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeBondDipole;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Components[type].ExcludedIntraChargeBondDipole[i].A;
        B=Components[type].ExcludedIntraChargeBondDipole[i].B;

        atom_pointer=Adsorbates[CurrentSystem][m].Atoms;

        chargeA=Components[type].Charge[A];
        posA=atom_pointer[A].Position;

        pairB=Components[type].BondDipoles[B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipole_magnitudeB=Components[type].BondDipoleMagnitude[B];
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posB.x-posA.x;
        dr.y=posB.y-posA.y;
        dr.z=posB.z-posA.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;

        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_adsorbate_excluded_c_bd+=Bt1*chargeA*cosB;

        fa1.x=chargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
        fa1.y=chargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
        fa1.z=chargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

        atom_pointer[A].Force.x+=fa1.x;
        atom_pointer[A].Force.y+=fa1.y;
        atom_pointer[A].Force.z+=fa1.z;

        fac1=chargeA*Bt1*dipole_magnitudeB/lengthB;
        fac2=chargeA*Bt1*cosB/(dipole_magnitudeB*lengthB);

        fb1.x=0.5*fa1.x+fac1*dr.x-fac2*dipoleB.x;
        fb1.y=0.5*fa1.y+fac1*dr.y-fac2*dipoleB.y;
        fb1.z=0.5*fa1.z+fac1*dr.z-fac2*dipoleB.z;

        fb2.x=0.5*fa1.x-fac1*dr.x+fac2*dipoleB.x;
        fb2.y=0.5*fa1.y-fac1*dr.y+fac2*dipoleB.y;
        fb2.z=0.5*fa1.z-fac1*dr.z+fac2*dipoleB.z;

        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=fa1.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=fa1.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=fa1.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=fa1.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=fa1.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=fa1.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=fa1.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=fa1.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=fa1.z*dr.z+v.cz;
      }
    }

    energy_adsorbate_excluded_bd=0.0;
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      type=Adsorbates[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraBondDipoleBondDipole;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Components[type].ExcludedIntraBondDipoleBondDipole[i].A;
        B=Components[type].ExcludedIntraBondDipoleBondDipole[i].B;

        atom_pointer=Adsorbates[CurrentSystem][m].Atoms;

        pairA=Components[type].BondDipoles[A];
        posA1=atom_pointer[pairA.A].Position;
        posA2=atom_pointer[pairA.B].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        dipole_magnitudeA=Components[type].BondDipoleMagnitude[A];
        lengthA=sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
        temp=dipole_magnitudeA/lengthA;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        pairB=Components[type].BondDipoles[B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipole_magnitudeB=Components[type].BondDipoleMagnitude[B];
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;
        temp*=2.0*SQR(alpha);
        Bt3=temp+(5.0/rr)*Bt2;

        cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
        cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_adsorbate_excluded_bd-=Bt1*cosAB-Bt2*cosA*cosB;

        term.x=(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
        term.y=(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
        term.z=(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

        termA.x=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.x/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.x-Bt1*dipoleB.x)*dipole_magnitudeA/lengthA;
        termA.y=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.y/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.y-Bt1*dipoleB.y)*dipole_magnitudeA/lengthA;
        termA.z=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.z/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.z-Bt1*dipoleB.z)*dipole_magnitudeA/lengthA;

        termB.x=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.x/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.x-Bt1*dipoleA.x)*dipole_magnitudeB/lengthB;
        termB.y=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.y/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.y-Bt1*dipoleA.y)*dipole_magnitudeB/lengthB;
        termB.z=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.z/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.z-Bt1*dipoleA.z)*dipole_magnitudeB/lengthB;

        fa1.x=0.5*term.x+termA.x;
        fa1.y=0.5*term.y+termA.y;
        fa1.z=0.5*term.z+termA.z;
        fa2.x=0.5*term.x-termA.x;
        fa2.y=0.5*term.y-termA.y;
        fa2.z=0.5*term.z-termA.z;

        fb1.x=-0.5*term.x+termB.x;
        fb1.y=-0.5*term.y+termB.y;
        fb1.z=-0.5*term.z+termB.z;
        fb2.x=-0.5*term.x-termB.x;
        fb2.y=-0.5*term.y-termB.y;
        fb2.z=-0.5*term.z-termB.z;

        atom_pointer[pairA.A].Force.x-=fa1.x;
        atom_pointer[pairA.A].Force.y-=fa1.y;
        atom_pointer[pairA.A].Force.z-=fa1.z;

        atom_pointer[pairA.B].Force.x-=fa2.x;
        atom_pointer[pairA.B].Force.y-=fa2.y;
        atom_pointer[pairA.B].Force.z-=fa2.z;

        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
      }
    }
  }

  // Exclusion pairs for the cations
  // ===============================

  // if cation-cation interactions are omitted, then the intra-molecular interactions are already fully elimated
  if(!OmitCationCationCoulombInteractions)
  {
    energy_cation_excluded=0.0;
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      atom_pointer=Cations[CurrentSystem][m].Atoms;
      type=Cations[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        pair=Components[type].ExcludedIntraChargeCharge[i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        chargeA=atom_pointer[pair.A].Charge;
        chargeB=atom_pointer[pair.B].Charge;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;

        energy_cation_excluded-=chargeA*chargeB*Bt0;

        temp=chargeA*chargeB*Bt1;

        DF=-chargeA*chargeB*Bt1;
        DDF=chargeA*chargeB*Bt2;

        atom_pointer[pair.A].Force.x+=temp*dr.x;
        atom_pointer[pair.A].Force.y+=temp*dr.y;
        atom_pointer[pair.A].Force.z+=temp*dr.z;

        atom_pointer[pair.B].Force.x-=temp*dr.x;
        atom_pointer[pair.B].Force.y-=temp*dr.y;
        atom_pointer[pair.B].Force.z-=temp*dr.z;

        stress.ax+=temp*dr.x*dr.x;
        stress.bx+=temp*dr.y*dr.x;
        stress.cx+=temp*dr.z*dr.x;

        stress.ay+=temp*dr.x*dr.y;
        stress.by+=temp*dr.y*dr.y;
        stress.cy+=temp*dr.z*dr.y;

        stress.az+=temp*dr.x*dr.z;
        stress.bz+=temp*dr.y*dr.z;
        stress.cz+=temp*dr.z*dr.z;

        AddContributionToBornTerm(DDF,DF,dr);

      }
    }

    energy_cation_excluded_c_bd=0.0;
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      type=Cations[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeBondDipole;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Components[type].ExcludedIntraChargeBondDipole[i].A;
        B=Components[type].ExcludedIntraChargeBondDipole[i].B;

        atom_pointer=Cations[CurrentSystem][m].Atoms;

        chargeA=Components[type].Charge[A];
        posA=atom_pointer[A].Position;

        pairB=Components[type].BondDipoles[B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipole_magnitudeB=Components[type].BondDipoleMagnitude[B];
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posB.x-posA.x;
        dr.y=posB.y-posA.y;
        dr.z=posB.z-posA.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;

        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_cation_excluded_c_bd+=Bt1*chargeA*cosB;

        term.x=chargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
        term.y=chargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
        term.z=chargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

        fa1.x=term.x;
        fa1.y=term.y;
        fa1.z=term.z;
        atom_pointer[A].Force.x+=fa1.x;
        atom_pointer[A].Force.y+=fa1.y;
        atom_pointer[A].Force.z+=fa1.z;

        fac1=chargeA*Bt1*dipole_magnitudeB/lengthB;
        fac2=chargeA*Bt1*cosB/(dipole_magnitudeB*lengthB);

        fb1.x=0.5*term.x+fac1*dr.x-fac2*dipoleB.x;
        fb1.y=0.5*term.y+fac1*dr.y-fac2*dipoleB.y;
        fb1.z=0.5*term.z+fac1*dr.z-fac2*dipoleB.z;

        fb2.x=0.5*term.x-fac1*dr.x+fac2*dipoleB.x;
        fb2.y=0.5*term.y-fac1*dr.y+fac2*dipoleB.y;
        fb2.z=0.5*term.z-fac1*dr.z+fac2*dipoleB.z;


        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=fa1.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=fa1.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=fa1.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=fa1.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=fa1.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=fa1.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=fa1.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=fa1.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=fa1.z*dr.z+v.cz;
      }
    }

    energy_cation_excluded_bd=0.0;
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      type=Cations[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraBondDipoleBondDipole;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        A=Components[type].ExcludedIntraBondDipoleBondDipole[i].A;
        B=Components[type].ExcludedIntraBondDipoleBondDipole[i].B;

        atom_pointer=Cations[CurrentSystem][m].Atoms;

        pairA=Components[type].BondDipoles[A];
        posA1=atom_pointer[pairA.A].Position;
        posA2=atom_pointer[pairA.B].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        dipole_magnitudeA=Components[type].BondDipoleMagnitude[A];
        lengthA=sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
        temp=dipole_magnitudeA/lengthA;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        pairB=Components[type].BondDipoles[B];
        posB1=atom_pointer[pairB.A].Position;
        posB2=atom_pointer[pairB.B].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        lengthB=sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
        dipole_magnitudeB=Components[type].BondDipoleMagnitude[B];
        temp=dipole_magnitudeB/lengthB;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;
        temp*=2.0*SQR(alpha);
        Bt2=temp+(3.0/rr)*Bt1;
        temp*=2.0*SQR(alpha);
        Bt3=temp+(5.0/rr)*Bt2;

        cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
        cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy_cation_excluded_bd-=Bt1*cosAB-Bt2*cosA*cosB;

        term.x=(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
        term.y=(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
        term.z=(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

        termA.x=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.x/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.x-Bt1*dipoleB.x)*dipole_magnitudeA/lengthA;
        termA.y=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.y/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.y-Bt1*dipoleB.y)*dipole_magnitudeA/lengthA;
        termA.z=(-Bt2*cosA*cosB+Bt1*cosAB)*dipoleA.z/(lengthA*dipole_magnitudeA)+
                (Bt2*cosB*dr.z-Bt1*dipoleB.z)*dipole_magnitudeA/lengthA;

        termB.x=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.x/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.x-Bt1*dipoleA.x)*dipole_magnitudeB/lengthB;
        termB.y=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.y/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.y-Bt1*dipoleA.y)*dipole_magnitudeB/lengthB;
        termB.z=(Bt1*cosAB-Bt2*cosA*cosB)*dipoleB.z/(lengthB*dipole_magnitudeB)+
                (Bt2*cosA*dr.z-Bt1*dipoleA.z)*dipole_magnitudeB/lengthB;

        fa1.x=0.5*term.x+termA.x;
        fa1.y=0.5*term.y+termA.y;
        fa1.z=0.5*term.z+termA.z;
        fa2.x=0.5*term.x-termA.x;
        fa2.y=0.5*term.y-termA.y;
        fa2.z=0.5*term.z-termA.z;

        fb1.x=-0.5*term.x+termB.x;
        fb1.y=-0.5*term.y+termB.y;
        fb1.z=-0.5*term.z+termB.z;
        fb2.x=-0.5*term.x-termB.x;
        fb2.y=-0.5*term.y-termB.y;
        fb2.z=-0.5*term.z-termB.z;

        atom_pointer[pairA.A].Force.x-=fa1.x;
        atom_pointer[pairA.A].Force.y-=fa1.y;
        atom_pointer[pairA.A].Force.z-=fa1.z;

        atom_pointer[pairA.B].Force.x-=fa2.x;
        atom_pointer[pairA.B].Force.y-=fa2.y;
        atom_pointer[pairA.B].Force.z-=fa2.z;

        atom_pointer[pairB.A].Force.x-=fb1.x;
        atom_pointer[pairB.A].Force.y-=fb1.y;
        atom_pointer[pairB.A].Force.z-=fb1.z;

        atom_pointer[pairB.B].Force.x-=fb2.x;
        atom_pointer[pairB.B].Force.y-=fb2.y;
        atom_pointer[pairB.B].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
      }
    }
  }

  StrainDerivativeTensor[CurrentSystem].ax-=stress.ax;
  StrainDerivativeTensor[CurrentSystem].ay-=stress.ay;
  StrainDerivativeTensor[CurrentSystem].az-=stress.az;

  StrainDerivativeTensor[CurrentSystem].bx-=stress.bx;
  StrainDerivativeTensor[CurrentSystem].by-=stress.by;
  StrainDerivativeTensor[CurrentSystem].bz-=stress.bz;

  StrainDerivativeTensor[CurrentSystem].cx-=stress.cx;
  StrainDerivativeTensor[CurrentSystem].cy-=stress.cy;
  StrainDerivativeTensor[CurrentSystem].cz-=stress.cz;

  // compute energies for host-host, host-adsorbate, host-cation, adsorbate-adsorbate, cation-cation, and adsorbate-cation
  // =====================================================================================================================
  // Note: total energy is Fourier minus self minus exclusion energy
  // Note: the net-charge is taken into account

  // Host-host energy
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    UHostHostChargeChargeFourier[CurrentSystem]=energy_framework_framework-energy_framework_self-energy_framework_excluded;
    UHostHostChargeBondDipoleFourier[CurrentSystem]=energy_framework_framework_c_bd-energy_framework_excluded_c_bd;
    UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=energy_framework_framework_bd-energy_framework_self_bd-energy_framework_excluded_bd;
    UHostHostChargeChargeFourier[CurrentSystem]-=UChargeChargeFrameworkRigid[CurrentSystem];
    UHostHostChargeBondDipoleFourier[CurrentSystem]-=UChargeBondDipoleFrameworkRigid[CurrentSystem];
    UHostHostBondDipoleBondDipoleFourier[CurrentSystem]-=UBondDipoleBondDipoleFrameworkRigid[CurrentSystem];
    UHostHostChargeChargeFourier[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeFramework[CurrentSystem]);
  }


  // Host-adsorbate energy (corrected for net-charge)
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=2.0*energy_framework_adsorbate;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=2.0*energy_framework_adsorbate_c_bd;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=2.0*energy_framework_adsorbate_bd;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*NetChargeAdsorbates[CurrentSystem];

  // Host-cation energy (corrected for net-charge)
  UHostCationChargeChargeFourier[CurrentSystem]=2.0*energy_framework_cation;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=2.0*energy_framework_cation_c_bd;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=2.0*energy_framework_cation_bd;
  UHostCationChargeChargeFourier[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*NetChargeCations[CurrentSystem];

  // Adsorbate-adsorbate energy (corrected for net-charge)
  if(!OmitAdsorbateAdsorbateCoulombInteractions)
  {
    UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=energy_adsorbate_adsorbate-energy_adsorbate_self-energy_adsorbate_excluded;
    UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=energy_adsorbate_adsorbate_c_bd-energy_adsorbate_excluded_c_bd;
    UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=energy_adsorbate_adsorbate_bd-energy_adsorbate_self_bd-energy_adsorbate_excluded_bd;
    UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeAdsorbates[CurrentSystem]);
  }

  // Cation-cation energy (corrected for net-charge)
  if(!OmitCationCationCoulombInteractions)
  {
    UCationCationChargeChargeFourier[CurrentSystem]=energy_cation_cation-energy_cation_self-energy_cation_excluded;
    UCationCationChargeBondDipoleFourier[CurrentSystem]=energy_cation_cation_c_bd-energy_cation_excluded_c_bd;
    UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=energy_cation_cation_bd-energy_cation_self_bd-energy_cation_excluded_bd;
    UCationCationChargeChargeFourier[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeCations[CurrentSystem]);
  }

  // Adsorbate-cation energy (corrected for net-charge)
  if(!OmitInterMolecularInteractions&&!OmitAdsorbateCationCoulombInteractions)
  {
    UAdsorbateCationChargeChargeFourier[CurrentSystem]=2.0*energy_adsorbate_cation;
    UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=2.0*energy_adsorbate_cation_c_bd;
    UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=2.0*energy_adsorbate_cation_bd;
    UAdsorbateCationChargeChargeFourier[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeAdsorbates[CurrentSystem]*NetChargeCations[CurrentSystem];
  }

  return 0;
}


/*********************************************************************************************************
 * Name       | CalculateEwaldFourierAdsorbate, CalculateEwaldFourierCation                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the change in Fourier sums of the charge-charge, charge-bonddipole, and         *
 *            | bonddipole-bonddipole energies in a truly periodic system.                               *                  
 * Parameters | 'New' boolean: whether this is a sum over a new molecule (insertion, Widom, re-insertion)*
 *            | 'Old' boolean: whether this is a sum over an old molecule.                               *
 *            | 'mol': the id of the chosen old molecule                                                 *
 *            | 'CurrentComponent' (global): the type of the new molecule.                               *
 * Note       | Largely based on the Fiche F.22 of Allen&Tildesly with a few modifications:              *
 *            | 1) exp(-ik.kx)*exp(-ik.ky) computed outside the 'kz' loop                                *
 *            | 2) energy are splitted in framework, adsorbate, cations groups                           *
 *            | 3) net-charge is taken into account                                                      *
 *            | 4) routine computes both charges and bond-dipoles                                        *
 *            | The modified sums are stored in 'NewTotalChargeAdsorbates' and                           *
 *            | 'NewTotalBondDipoleAdsorbates' and will be written to usually arrays on acceptance.      *
 *            | The routine is used in all mc-moves.                                                     *
 *********************************************************************************************************/

int CalculateEwaldFourierAdsorbate(int NewMolecule,int OldMolecule,int mol,int store)
{
  int i,j,ii,jj,kk;
  int A,B,nvec,nr_of_excluded_pairs;
  int kmax_x,kmax_y,kmax_z,index_i,index_j,index_k;
  int type_mol,nr_atoms,type,nr_of_bonddipoles;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_old,nr_of_coulombic_sites_new;
  int nr_of_bonddipole_sites,nr_of_bonddipole_sites_old,nr_of_bonddipole_sites_new;
  COMPLEX sum_old,sum_new,sum_adsorbates,sum_bonddipole_adsorbates;
  COMPLEX sum_bonddipole_old,sum_bonddipole_new;
  REAL fac,energy_charge_adsorbates,energy_charge_adsorbates_cations,energy_charge_framework_adsorbates;
  REAL energy_bonddipole_adsorbates,energy_bonddipole_adsorbates_cations,energy_bonddipole_framework_adsorbates;
  REAL alpha,chargeA,chargeB,charge,r,rr;
  REAL energy_self_new,energy_self_old;
  REAL net_charge_new,net_charge_old;
  REAL scaling,scalingA,scalingB;
  REAL cosA,cosB,cosAB,Bt0,Bt1,Bt2,temp;
  REAL energy_self_bd_old,energy_self_bd_new;
  REAL energy_charge_bonddipole_adsorbates,energy_charge_bonddipole_adsorbates_cations;
  REAL energy_charge_bonddipole_framework_adsorbates;
  REAL energy_excluded_new,energy_excluded_old;
  REAL energy_excluded_c_bd_new,energy_excluded_c_bd_old;
  REAL energy_excluded_bd_new,energy_excluded_bd_old;
  VECTOR pos,posA,posB,dr;
  VECTOR dipole,dipoleA,dipoleB,rk;
  VECTOR posA1,posA2,posB1,posB2;
  VECTOR *kvecs;
  REAL *kfactor,recip_cutoff,ksqr;
  PAIR pair;
  ATOM *atom_pointer;
  int considered_charged;

  // intialize differences in energy
  NetChargeAdsorbateDelta=0.0;
  UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]=0.0;
  UHostCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  UCationCationChargeChargeFourierDelta[CurrentSystem]=0.0;

  UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;

  UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  nvec=0;
  energy_self_new=energy_self_old=0.0;
  energy_self_bd_old=energy_self_bd_new=0.0;
  net_charge_new=net_charge_old=0.0;
  energy_charge_adsorbates=0.0;
  energy_charge_adsorbates_cations=0.0;
  energy_charge_framework_adsorbates=0.0;
  energy_charge_bonddipole_adsorbates=0.0;
  energy_charge_bonddipole_adsorbates_cations=0.0;
  energy_charge_bonddipole_framework_adsorbates=0.0;
  energy_bonddipole_adsorbates=0.0;
  energy_bonddipole_adsorbates_cations=0.0;
  energy_bonddipole_framework_adsorbates=0.0;
  fac=0.0;

  nr_of_coulombic_sites=nr_of_coulombic_sites_old=nr_of_coulombic_sites_new=0;
  nr_of_bonddipole_sites=nr_of_bonddipole_sites_old=nr_of_bonddipole_sites_new=0;

  if(OldMolecule)
  {
    nr_atoms=Adsorbates[CurrentSystem][mol].NumberOfAtoms;
    for(j=0;j<nr_atoms;j++)
    {
      type=Adsorbates[CurrentSystem][mol].Atoms[j].Type;
      scaling=Adsorbates[CurrentSystem][mol].Atoms[j].CFChargeScalingParameter;
      charge=Adsorbates[CurrentSystem][mol].Atoms[j].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      charge*=scaling;
      if(considered_charged)
      {
        Charge[nr_of_coulombic_sites]=charge;
        energy_self_old+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
        net_charge_old+=Charge[nr_of_coulombic_sites];
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][mol].Atoms[j].Position);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_old=nr_of_coulombic_sites;

  if(NewMolecule)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      type=Components[CurrentComponent].Type[i];
      if(PseudoAtoms[type].HasCharges)
      {
        scaling=CFChargeScaling[i];
        Charge[nr_of_coulombic_sites]=scaling*Components[CurrentComponent].Charge[i];
        energy_self_new+=COULOMBIC_CONVERSION_FACTOR*SQR(Charge[nr_of_coulombic_sites])*alpha/sqrt(M_PI);
        net_charge_new+=Charge[nr_of_coulombic_sites];
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(TrialPosition[CurrentSystem][i]);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_new=nr_of_coulombic_sites;

  if(OldMolecule)
  {
    type_mol=Adsorbates[CurrentSystem][mol].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Adsorbates[CurrentSystem][mol].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[j];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      energy_self_bd_old+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_old=nr_of_bonddipole_sites;

  if(NewMolecule)
  {
    type_mol=CurrentComponent;
    for(i=0;i<Components[type_mol].NumberOfBondDipoles;i++)
    {
      pair=Components[type_mol].BondDipoles[i];
      posA=TrialPosition[CurrentSystem][pair.A];
      posB=TrialPosition[CurrentSystem][pair.B];
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[i];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      energy_self_bd_new+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_new=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }


  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          sum_old.re=0.0;
          sum_old.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_old.re+=temp*Eikr[i].re;
            sum_old.im+=temp*Eikr[i].im;
          }
          sum_new.re=0.0;
          sum_new.im=0.0;
          for(;i<nr_of_coulombic_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_new.re+=temp*Eikr[i].re;
            sum_new.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_old.re=0.0;
          sum_bonddipole_old.im=0.0;
          rk=kvecs[nvec];
          for(i=0;i<nr_of_bonddipole_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_old.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_old.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_new.re=0.0;
          sum_bonddipole_new.im=0.0;
          for(;i<nr_of_bonddipole_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_new.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_new.im+=temp*Eikr_bd[i].im;
          }

          sum_adsorbates.re=StoreTotalChargeAdsorbates[CurrentSystem][nvec].re+(sum_new.re-sum_old.re);
          sum_adsorbates.im=StoreTotalChargeAdsorbates[CurrentSystem][nvec].im+(sum_new.im-sum_old.im);

          sum_bonddipole_adsorbates.re=StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re+(sum_bonddipole_new.re-sum_bonddipole_old.re);
          sum_bonddipole_adsorbates.im=StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im+(sum_bonddipole_new.im-sum_bonddipole_old.im);

          temp=kfactor[nvec];

          // compute energy differences using the stored total sums and the sum of the differences of the moving atoms
          if(Framework[CurrentSystem].FrameworkModel!=NONE)
          {
            energy_charge_framework_adsorbates+=temp*
              (StoreTotalChargeFramework[CurrentSystem][nvec].re*(sum_adsorbates.re-StoreTotalChargeAdsorbates[CurrentSystem][nvec].re)
              +StoreTotalChargeFramework[CurrentSystem][nvec].im*(sum_adsorbates.im-StoreTotalChargeAdsorbates[CurrentSystem][nvec].im));

            energy_charge_bonddipole_framework_adsorbates+=temp*
              (StoreTotalChargeFramework[CurrentSystem][nvec].im*(sum_bonddipole_adsorbates.re-StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re)
              +StoreTotalChargeFramework[CurrentSystem][nvec].re*(StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im-sum_bonddipole_adsorbates.im)
              +StoreTotalBondDipolesFramework[CurrentSystem][nvec].re*(sum_adsorbates.im-StoreTotalChargeAdsorbates[CurrentSystem][nvec].im)
              +StoreTotalBondDipolesFramework[CurrentSystem][nvec].im*(StoreTotalChargeAdsorbates[CurrentSystem][nvec].re-sum_adsorbates.re));

            energy_bonddipole_framework_adsorbates+=temp*
              (StoreTotalBondDipolesFramework[CurrentSystem][nvec].re*(sum_bonddipole_adsorbates.re-StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re)
              +StoreTotalBondDipolesFramework[CurrentSystem][nvec].im*(sum_bonddipole_adsorbates.im-StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im));
          }

          energy_charge_adsorbates+=temp*(SQR(sum_adsorbates.re)-SQR(StoreTotalChargeAdsorbates[CurrentSystem][nvec].re)+
                                          SQR(sum_adsorbates.im)-SQR(StoreTotalChargeAdsorbates[CurrentSystem][nvec].im));

          energy_charge_bonddipole_adsorbates+=2.0*temp*
                  (sum_adsorbates.im*sum_bonddipole_adsorbates.re-sum_adsorbates.re*sum_bonddipole_adsorbates.im
                  -StoreTotalChargeAdsorbates[CurrentSystem][nvec].im*StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re
                  +StoreTotalChargeAdsorbates[CurrentSystem][nvec].re*StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im);

          energy_bonddipole_adsorbates+=temp*(SQR(sum_bonddipole_adsorbates.re)-SQR(StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re)+
                                              SQR(sum_bonddipole_adsorbates.im)-SQR(StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im));

          if(NumberOfCationMolecules[CurrentSystem]>0)
          {
            energy_charge_adsorbates_cations+=temp*
              (StoreTotalChargeCations[CurrentSystem][nvec].re*(sum_adsorbates.re-StoreTotalChargeAdsorbates[CurrentSystem][nvec].re)
              +StoreTotalChargeCations[CurrentSystem][nvec].im*(sum_adsorbates.im-StoreTotalChargeAdsorbates[CurrentSystem][nvec].im));

            energy_charge_bonddipole_adsorbates_cations+=temp*
              (StoreTotalChargeCations[CurrentSystem][nvec].im*(sum_bonddipole_adsorbates.re-StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re)
              +StoreTotalChargeCations[CurrentSystem][nvec].re*(StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im-sum_bonddipole_adsorbates.im)
              +StoreTotalBondDipolesCations[CurrentSystem][nvec].re*(sum_adsorbates.im-StoreTotalChargeAdsorbates[CurrentSystem][nvec].im)
              +StoreTotalBondDipolesCations[CurrentSystem][nvec].im*(StoreTotalChargeAdsorbates[CurrentSystem][nvec].re-sum_adsorbates.re));

            energy_bonddipole_adsorbates_cations+=temp*
              (StoreTotalBondDipolesCations[CurrentSystem][nvec].re*(sum_bonddipole_adsorbates.re-StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re)
              +StoreTotalBondDipolesCations[CurrentSystem][nvec].im*(sum_bonddipole_adsorbates.im-StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im));
          }

          // store the new sums, these will be the current ones on acceptance of the mc-move
          NewTotalChargeAdsorbates[store][nvec]=sum_adsorbates;
          NewTotalBondDipolesAdsorbates[store][nvec]=sum_bonddipole_adsorbates;

          // next wave-vector
          nvec++;
        }
      }
    }
  }

  energy_excluded_new=0.0;
  energy_excluded_c_bd_new=0.0;
  energy_excluded_bd_new=0.0;
  if(NewMolecule)
  {
    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraChargeCharge[i].A;
      B=Components[CurrentComponent].ExcludedIntraChargeCharge[i].B;
      scalingA=CFChargeScaling[A];
      scalingB=CFChargeScaling[B];
      chargeA=scalingA*Components[CurrentComponent].Charge[A];
      chargeB=scalingB*Components[CurrentComponent].Charge[B];
      posA=TrialPosition[CurrentSystem][A];
      posB=TrialPosition[CurrentSystem][B];

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      Bt0=-erf(alpha*r)/r;
      energy_excluded_new-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*Bt0;
    }

    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraChargeBondDipole[i].A;
      B=Components[CurrentComponent].ExcludedIntraChargeBondDipole[i].B;

      chargeA=Components[CurrentComponent].Charge[A];
      posA=TrialPosition[CurrentSystem][A];

      pair=Components[CurrentComponent].BondDipoles[B];
      posB1=TrialPosition[CurrentSystem][pair.A];
      posB2=TrialPosition[CurrentSystem][pair.B];

      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(alpha*r)/r;
      Bt1=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr+Bt0/rr;

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_c_bd_new+=COULOMBIC_CONVERSION_FACTOR*Bt1*chargeA*cosB;
    }


    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[CurrentComponent].ExcludedIntraBondDipoleBondDipole[i].B;

      pair=Components[CurrentComponent].BondDipoles[A];
      posA1=TrialPosition[CurrentSystem][pair.A];
      posA2=TrialPosition[CurrentSystem][pair.B];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      pair=Components[CurrentComponent].BondDipoles[B];
      posB1=TrialPosition[CurrentSystem][pair.A];
      posB2=TrialPosition[CurrentSystem][pair.B];
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      temp=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
      Bt0=-erf(alpha*r)/r;
      Bt1=temp+Bt0/rr;
      temp*=2.0*SQR(alpha);
      Bt2=temp+(3.0/rr)*Bt1;

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_bd_new-=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
    }
  }

  energy_excluded_old=0.0;
  energy_excluded_c_bd_old=0.0;
  energy_excluded_bd_old=0.0;
  if(OldMolecule)
  {
    type=Adsorbates[CurrentSystem][mol].Type;

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      pair=Components[type].ExcludedIntraChargeCharge[i];
      scalingA=Adsorbates[CurrentSystem][mol].Atoms[pair.A].CFChargeScalingParameter;
      scalingB=Adsorbates[CurrentSystem][mol].Atoms[pair.B].CFChargeScalingParameter;
      chargeA=scalingA*Adsorbates[CurrentSystem][mol].Atoms[pair.A].Charge;
      chargeB=scalingB*Adsorbates[CurrentSystem][mol].Atoms[pair.B].Charge;
      posA=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posB=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

      Bt0=-erf(alpha*r)/r;
      energy_excluded_old-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*Bt0;
    }

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraChargeBondDipole[i].A;
      B=Components[type].ExcludedIntraChargeBondDipole[i].B;

      chargeA=Adsorbates[CurrentSystem][mol].Atoms[A].Charge;
      posA=Adsorbates[CurrentSystem][mol].Atoms[A].Position;

      pair=Components[type].BondDipoles[B];
      posB1=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posB2=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(alpha*r)/r;
      Bt1=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr+Bt0/rr;

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_c_bd_old+=COULOMBIC_CONVERSION_FACTOR*Bt1*chargeA*cosB;
    }

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[type].ExcludedIntraBondDipoleBondDipole[i].B;

      pair=Components[type].BondDipoles[A];
      posA1=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posA2=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=Components[type].BondDipoleMagnitude[A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      pair=Components[type].BondDipoles[B];
      posB1=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posB2=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      temp=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
      Bt0=-erf(alpha*r)/r;
      Bt1=temp+Bt0/rr;
      temp*=2.0*SQR(alpha);
      Bt2=temp+(3.0/rr)*Bt1;

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_bd_old-=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
    }
  }

  // set net-charge difference
  NetChargeAdsorbateDelta=net_charge_new-net_charge_old;

  // set energy differences
  if(!OmitAdsorbateAdsorbateCoulombInteractions)
  {
    UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]=energy_charge_adsorbates-(energy_excluded_new-energy_excluded_old)-(energy_self_new-energy_self_old);
    UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=energy_charge_bonddipole_adsorbates-(energy_excluded_c_bd_new-energy_excluded_c_bd_old);
    UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=energy_bonddipole_adsorbates-(energy_self_bd_new-energy_self_bd_old)-
                                     (energy_excluded_bd_new-energy_excluded_bd_old);
    UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeAdsorbates[CurrentSystem]+NetChargeAdsorbateDelta)-
                                                                UIon[CurrentSystem]*SQR(NetChargeAdsorbates[CurrentSystem]);
  }

  if(!OmitInterMolecularInteractions)
  {
    UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]=2.0*energy_charge_adsorbates_cations;
    UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]=2.0*energy_charge_bonddipole_adsorbates_cations;
    UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=2.0*energy_bonddipole_adsorbates_cations;
    UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeCations[CurrentSystem]*(NetChargeAdsorbates[CurrentSystem]+NetChargeAdsorbateDelta)-
                                                             2.0*UIon[CurrentSystem]*NetChargeCations[CurrentSystem]*(NetChargeAdsorbates[CurrentSystem]);
  }

  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]=2.0*energy_charge_framework_adsorbates;
  UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=2.0*energy_charge_bonddipole_framework_adsorbates;
  UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=2.0*energy_bonddipole_framework_adsorbates;
  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*(NetChargeAdsorbates[CurrentSystem]+NetChargeAdsorbateDelta)-
                                                         2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*(NetChargeAdsorbates[CurrentSystem]);

  return 0;
}

int CalculateEwaldFourierAdsorbate2(int NewMolecule,int OldMolecule,int mol,int store)
{
  int i,j,ii,jj,kk;
  int A,B,nvec,nr_of_excluded_pairs;
  int kmax_x,kmax_y,kmax_z,index_i,index_j,index_k;
  int type_mol,nr_atoms,type,nr_of_bonddipoles;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_old,nr_of_coulombic_sites_new;
  int nr_of_bonddipole_sites,nr_of_bonddipole_sites_old,nr_of_bonddipole_sites_new;
  COMPLEX sum_old,sum_new,sum_adsorbates,sum_bonddipole_adsorbates;
  COMPLEX sum_bonddipole_old,sum_bonddipole_new;
  REAL fac,energy_charge_adsorbates,energy_charge_adsorbates_cations,energy_charge_framework_adsorbates;
  REAL energy_bonddipole_adsorbates,energy_bonddipole_adsorbates_cations,energy_bonddipole_framework_adsorbates;
  REAL alpha,chargeA,chargeB,charge,r,rr;
  REAL energy_self_new,energy_self_old;
  REAL net_charge_new,net_charge_old;
  REAL scaling,scalingA,scalingB;
  REAL cosA,cosB,cosAB,Bt0,Bt1,Bt2,temp;
  REAL energy_self_bd_old,energy_self_bd_new;
  REAL energy_charge_bonddipole_adsorbates,energy_charge_bonddipole_adsorbates_cations;
  REAL energy_charge_bonddipole_framework_adsorbates;
  REAL energy_excluded_new,energy_excluded_old;
  REAL energy_excluded_c_bd_new,energy_excluded_c_bd_old;
  REAL energy_excluded_bd_new,energy_excluded_bd_old;
  VECTOR pos,posA,posB,dr;
  VECTOR dipole,dipoleA,dipoleB,rk;
  VECTOR posA1,posA2,posB1,posB2;
  VECTOR *kvecs;
  REAL *kfactor,recip_cutoff,ksqr;
  PAIR pair;
  ATOM *atom_pointer;
  int considered_charged;

  // intialize differences in energy
  NetChargeAdsorbateDelta=0.0;
  UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]=0.0;
  UHostCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  UCationCationChargeChargeFourierDelta[CurrentSystem]=0.0;

  UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;

  UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  nvec=0;
  energy_self_new=energy_self_old=0.0;
  energy_self_bd_old=energy_self_bd_new=0.0;
  net_charge_new=net_charge_old=0.0;
  energy_charge_adsorbates=0.0;
  energy_charge_adsorbates_cations=0.0;
  energy_charge_framework_adsorbates=0.0;
  energy_charge_bonddipole_adsorbates=0.0;
  energy_charge_bonddipole_adsorbates_cations=0.0;
  energy_charge_bonddipole_framework_adsorbates=0.0;
  energy_bonddipole_adsorbates=0.0;
  energy_bonddipole_adsorbates_cations=0.0;
  energy_bonddipole_framework_adsorbates=0.0;
  fac=0.0;

  nr_of_coulombic_sites=nr_of_coulombic_sites_old=nr_of_coulombic_sites_new=0;
  nr_of_bonddipole_sites=nr_of_bonddipole_sites_old=nr_of_bonddipole_sites_new=0;

  if(OldMolecule)
  {
    nr_atoms=Adsorbates[CurrentSystem][mol].NumberOfAtoms;
    for(j=0;j<nr_atoms;j++)
    {
      type=Adsorbates[CurrentSystem][mol].Atoms[j].Type;
      scaling=Adsorbates[CurrentSystem][mol].Atoms[j].CFChargeScalingParameter;
      charge=Adsorbates[CurrentSystem][mol].Atoms[j].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      charge*=scaling;
      if(considered_charged)
      {
        Charge[nr_of_coulombic_sites]=charge;
        energy_self_old+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
        net_charge_old+=Charge[nr_of_coulombic_sites];
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][mol].Atoms[j].Position);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_old=nr_of_coulombic_sites;

  if(NewMolecule)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      type=Components[CurrentComponent].Type[i];
      if(PseudoAtoms[type].HasCharges)
      {
        scaling=CFChargeScaling[i];
        Charge[nr_of_coulombic_sites]=scaling*Components[CurrentComponent].Charge[i];
        energy_self_new+=COULOMBIC_CONVERSION_FACTOR*SQR(Charge[nr_of_coulombic_sites])*alpha/sqrt(M_PI);
        net_charge_new+=Charge[nr_of_coulombic_sites];
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(TrialPosition[CurrentSystem][i]);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_new=nr_of_coulombic_sites;

  if(OldMolecule)
  {
    type_mol=Adsorbates[CurrentSystem][mol].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Adsorbates[CurrentSystem][mol].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[j];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      energy_self_bd_old+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_old=nr_of_bonddipole_sites;

  if(NewMolecule)
  {
    type_mol=CurrentComponent;
    for(i=0;i<Components[type_mol].NumberOfBondDipoles;i++)
    {
      pair=Components[type_mol].BondDipoles[i];
      posA=TrialPosition[CurrentSystem][pair.A];
      posB=TrialPosition[CurrentSystem][pair.B];
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[i];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      energy_self_bd_new+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_new=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }


  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          sum_old.re=0.0;
          sum_old.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_old.re+=temp*Eikr[i].re;
            sum_old.im+=temp*Eikr[i].im;
          }
          sum_new.re=0.0;
          sum_new.im=0.0;
          for(;i<nr_of_coulombic_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_new.re+=temp*Eikr[i].re;
            sum_new.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_old.re=0.0;
          sum_bonddipole_old.im=0.0;
          rk=kvecs[nvec];
          for(i=0;i<nr_of_bonddipole_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_old.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_old.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_new.re=0.0;
          sum_bonddipole_new.im=0.0;
          for(;i<nr_of_bonddipole_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_new.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_new.im+=temp*Eikr_bd[i].im;
          }

          sum_adsorbates.re=NewTotalChargeAdsorbates[CurrentSystem][nvec].re+(sum_new.re-sum_old.re);
          sum_adsorbates.im=NewTotalChargeAdsorbates[CurrentSystem][nvec].im+(sum_new.im-sum_old.im);

          sum_bonddipole_adsorbates.re=NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].re+(sum_bonddipole_new.re-sum_bonddipole_old.re);
          sum_bonddipole_adsorbates.im=NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].im+(sum_bonddipole_new.im-sum_bonddipole_old.im);

          temp=kfactor[nvec];

          // compute energy differences using the stored total sums and the sum of the differences of the moving atoms
          if(Framework[CurrentSystem].FrameworkModel!=NONE)
          {
            energy_charge_framework_adsorbates+=temp*
              (StoreTotalChargeFramework[CurrentSystem][nvec].re*(sum_adsorbates.re-NewTotalChargeAdsorbates[CurrentSystem][nvec].re)
              +StoreTotalChargeFramework[CurrentSystem][nvec].im*(sum_adsorbates.im-NewTotalChargeAdsorbates[CurrentSystem][nvec].im));

            energy_charge_bonddipole_framework_adsorbates+=temp*
              (StoreTotalChargeFramework[CurrentSystem][nvec].im*(sum_bonddipole_adsorbates.re-NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].re)
              +StoreTotalChargeFramework[CurrentSystem][nvec].re*(NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].im-sum_bonddipole_adsorbates.im)
              +NewTotalBondDipolesFramework[CurrentSystem][nvec].re*(sum_adsorbates.im-NewTotalChargeAdsorbates[CurrentSystem][nvec].im)
              +NewTotalBondDipolesFramework[CurrentSystem][nvec].im*(NewTotalChargeAdsorbates[CurrentSystem][nvec].re-sum_adsorbates.re));

            energy_bonddipole_framework_adsorbates+=temp*
              (NewTotalBondDipolesFramework[CurrentSystem][nvec].re*(sum_bonddipole_adsorbates.re-NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].re)
              +NewTotalBondDipolesFramework[CurrentSystem][nvec].im*(sum_bonddipole_adsorbates.im-NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].im));
          }

          energy_charge_adsorbates+=temp*(SQR(sum_adsorbates.re)-SQR(NewTotalChargeAdsorbates[CurrentSystem][nvec].re)+
                                          SQR(sum_adsorbates.im)-SQR(NewTotalChargeAdsorbates[CurrentSystem][nvec].im));

          energy_charge_bonddipole_adsorbates+=2.0*temp*
                  (sum_adsorbates.im*sum_bonddipole_adsorbates.re-sum_adsorbates.re*sum_bonddipole_adsorbates.im
                  -NewTotalChargeAdsorbates[CurrentSystem][nvec].im*NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].re
                  +NewTotalChargeAdsorbates[CurrentSystem][nvec].re*NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].im);

          energy_bonddipole_adsorbates+=temp*(SQR(sum_bonddipole_adsorbates.re)-SQR(NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].re)+
                                              SQR(sum_bonddipole_adsorbates.im)-SQR(NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].im));

          if(NumberOfCationMolecules[CurrentSystem]>0)
          {
            energy_charge_adsorbates_cations+=temp*
              (NewTotalChargeCations[CurrentSystem][nvec].re*(sum_adsorbates.re-NewTotalChargeAdsorbates[CurrentSystem][nvec].re)
              +NewTotalChargeCations[CurrentSystem][nvec].im*(sum_adsorbates.im-NewTotalChargeAdsorbates[CurrentSystem][nvec].im));

            energy_charge_bonddipole_adsorbates_cations+=temp*
              (NewTotalChargeCations[CurrentSystem][nvec].im*(sum_bonddipole_adsorbates.re-NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].re)
              +NewTotalChargeCations[CurrentSystem][nvec].re*(NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].im-sum_bonddipole_adsorbates.im)
              +NewTotalBondDipolesCations[CurrentSystem][nvec].re*(sum_adsorbates.im-NewTotalChargeAdsorbates[CurrentSystem][nvec].im)
              +NewTotalBondDipolesCations[CurrentSystem][nvec].im*(NewTotalChargeAdsorbates[CurrentSystem][nvec].re-sum_adsorbates.re));

            energy_bonddipole_adsorbates_cations+=temp*
              (NewTotalBondDipolesCations[CurrentSystem][nvec].re*(sum_bonddipole_adsorbates.re-NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].re)
              +NewTotalBondDipolesCations[CurrentSystem][nvec].im*(sum_bonddipole_adsorbates.im-NewTotalBondDipolesAdsorbates[CurrentSystem][nvec].im));
          }

          // store the new sums, these will be the current ones on acceptance of the mc-move
          NewTotalChargeAdsorbates[store][nvec]=sum_adsorbates;
          NewTotalBondDipolesAdsorbates[store][nvec]=sum_bonddipole_adsorbates;

          // next wave-vector
          nvec++;
        }
      }
    }
  }

  energy_excluded_new=0.0;
  energy_excluded_c_bd_new=0.0;
  energy_excluded_bd_new=0.0;
  if(NewMolecule)
  {
    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraChargeCharge[i].A;
      B=Components[CurrentComponent].ExcludedIntraChargeCharge[i].B;
      scalingA=CFChargeScaling[A];
      scalingB=CFChargeScaling[B];
      chargeA=scalingA*Components[CurrentComponent].Charge[A];
      chargeB=scalingB*Components[CurrentComponent].Charge[B];
      posA=TrialPosition[CurrentSystem][A];
      posB=TrialPosition[CurrentSystem][B];

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      Bt0=-erf(alpha*r)/r;
      energy_excluded_new-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*Bt0;
    }

    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraChargeBondDipole[i].A;
      B=Components[CurrentComponent].ExcludedIntraChargeBondDipole[i].B;

      chargeA=Components[CurrentComponent].Charge[A];
      posA=TrialPosition[CurrentSystem][A];

      pair=Components[CurrentComponent].BondDipoles[B];
      posB1=TrialPosition[CurrentSystem][pair.A];
      posB2=TrialPosition[CurrentSystem][pair.B];

      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(alpha*r)/r;
      Bt1=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr+Bt0/rr;

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_c_bd_new+=COULOMBIC_CONVERSION_FACTOR*Bt1*chargeA*cosB;
    }


    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[CurrentComponent].ExcludedIntraBondDipoleBondDipole[i].B;

      pair=Components[CurrentComponent].BondDipoles[A];
      posA1=TrialPosition[CurrentSystem][pair.A];
      posA2=TrialPosition[CurrentSystem][pair.B];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      pair=Components[CurrentComponent].BondDipoles[B];
      posB1=TrialPosition[CurrentSystem][pair.A];
      posB2=TrialPosition[CurrentSystem][pair.B];
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      temp=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
      Bt0=-erf(alpha*r)/r;
      Bt1=temp+Bt0/rr;
      temp*=2.0*SQR(alpha);
      Bt2=temp+(3.0/rr)*Bt1;

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_bd_new-=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
    }
  }

  energy_excluded_old=0.0;
  energy_excluded_c_bd_old=0.0;
  energy_excluded_bd_old=0.0;
  if(OldMolecule)
  {
    type=Adsorbates[CurrentSystem][mol].Type;

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      pair=Components[type].ExcludedIntraChargeCharge[i];
      scalingA=Adsorbates[CurrentSystem][mol].Atoms[pair.A].CFChargeScalingParameter;
      scalingB=Adsorbates[CurrentSystem][mol].Atoms[pair.B].CFChargeScalingParameter;
      chargeA=scalingA*Adsorbates[CurrentSystem][mol].Atoms[pair.A].Charge;
      chargeB=scalingB*Adsorbates[CurrentSystem][mol].Atoms[pair.B].Charge;
      posA=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posB=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

      Bt0=-erf(alpha*r)/r;
      energy_excluded_old-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*Bt0;
    }

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraChargeBondDipole[i].A;
      B=Components[type].ExcludedIntraChargeBondDipole[i].B;

      chargeA=Adsorbates[CurrentSystem][mol].Atoms[A].Charge;
      posA=Adsorbates[CurrentSystem][mol].Atoms[A].Position;

      pair=Components[type].BondDipoles[B];
      posB1=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posB2=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(alpha*r)/r;
      Bt1=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr+Bt0/rr;

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_c_bd_old+=COULOMBIC_CONVERSION_FACTOR*Bt1*chargeA*cosB;
    }

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[type].ExcludedIntraBondDipoleBondDipole[i].B;

      pair=Components[type].BondDipoles[A];
      posA1=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posA2=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=Components[type].BondDipoleMagnitude[A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      pair=Components[type].BondDipoles[B];
      posB1=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posB2=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      temp=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
      Bt0=-erf(alpha*r)/r;
      Bt1=temp+Bt0/rr;
      temp*=2.0*SQR(alpha);
      Bt2=temp+(3.0/rr)*Bt1;

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_bd_old-=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
    }
  }

  // set net-charge difference
  NetChargeAdsorbateDelta=net_charge_new-net_charge_old;

  // set energy differences
  if(!OmitAdsorbateAdsorbateCoulombInteractions)
  {
    UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]=energy_charge_adsorbates-(energy_excluded_new-energy_excluded_old)-(energy_self_new-energy_self_old);
    UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=energy_charge_bonddipole_adsorbates-(energy_excluded_c_bd_new-energy_excluded_c_bd_old);
    UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=energy_bonddipole_adsorbates-(energy_self_bd_new-energy_self_bd_old)-
                                     (energy_excluded_bd_new-energy_excluded_bd_old);
    UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeAdsorbates[CurrentSystem]+NetChargeAdsorbateDelta)-
                                                                UIon[CurrentSystem]*SQR(NetChargeAdsorbates[CurrentSystem]);
  }

  if(!OmitInterMolecularInteractions)
  {
    UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]=2.0*energy_charge_adsorbates_cations;
    UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]=2.0*energy_charge_bonddipole_adsorbates_cations;
    UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=2.0*energy_bonddipole_adsorbates_cations;
    UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeCations[CurrentSystem]*(NetChargeAdsorbates[CurrentSystem]+NetChargeAdsorbateDelta)-
                                                             2.0*UIon[CurrentSystem]*NetChargeCations[CurrentSystem]*(NetChargeAdsorbates[CurrentSystem]);
  }

  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]=2.0*energy_charge_framework_adsorbates;
  UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=2.0*energy_charge_bonddipole_framework_adsorbates;
  UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=2.0*energy_bonddipole_framework_adsorbates;
  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*(NetChargeAdsorbates[CurrentSystem]+NetChargeAdsorbateDelta)-
                                                         2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*(NetChargeAdsorbates[CurrentSystem]);

  return 0;
}

int CalculateEwaldFourierCation(int NewMolecule,int OldMolecule,int mol,int store)
{
  int i,j,ii,jj,kk;
  int A,B,nvec,nr_of_excluded_pairs;
  int kmax_x,kmax_y,kmax_z,index_i,index_j,index_k;
  int type_mol,nr_atoms,type,nr_of_bonddipoles;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_old,nr_of_coulombic_sites_new;
  int nr_of_bonddipole_sites,nr_of_bonddipole_sites_old,nr_of_bonddipole_sites_new;
  COMPLEX sum_old,sum_new,sum_cations,sum_bonddipole_cations;
  COMPLEX sum_bonddipole_old,sum_bonddipole_new;
  REAL fac,energy_charge_cations,energy_charge_adsorbates_cations,energy_charge_framework_cations;
  REAL energy_bonddipole_cations,energy_bonddipole_adsorbates_cations,energy_bonddipole_framework_cations;
  REAL alpha,chargeA,chargeB,charge,r,rr;
  REAL energy_self_new,energy_self_old;
  REAL net_charge_new,net_charge_old;
  REAL scaling,scalingA,scalingB;
  REAL cosA,cosB,cosAB,Bt0,Bt1,Bt2,temp;
  REAL energy_self_bd_old,energy_self_bd_new;
  REAL energy_charge_bonddipole_cations,energy_charge_bonddipole_adsorbates_cations;
  REAL energy_charge_bonddipole_framework_cations;
  REAL energy_excluded_new,energy_excluded_old;
  REAL energy_excluded_c_bd_new,energy_excluded_c_bd_old;
  REAL energy_excluded_bd_new,energy_excluded_bd_old;
  VECTOR pos,posA,posB,dr;
  VECTOR dipole,dipoleA,dipoleB,rk;
  VECTOR posA1,posA2,posB1,posB2;
  VECTOR *kvecs;
  REAL *kfactor,recip_cutoff,ksqr;
  PAIR pair;
  ATOM *atom_pointer;
  int considered_charged;

  // intialize differences in energy
  NetChargeCationDelta=0.0;
  UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]=0.0;
  UHostCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  UCationCationChargeChargeFourierDelta[CurrentSystem]=0.0;

  UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;

  UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  nvec=0;
  energy_self_new=energy_self_old=0.0;
  energy_self_bd_old=energy_self_bd_new=0.0;
  net_charge_new=net_charge_old=0.0;
  energy_charge_cations=0.0;
  energy_charge_adsorbates_cations=0.0;
  energy_charge_framework_cations=0.0;
  energy_charge_bonddipole_cations=0.0;
  energy_charge_bonddipole_adsorbates_cations=0.0;
  energy_charge_bonddipole_framework_cations=0.0;
  energy_bonddipole_cations=0.0;
  energy_bonddipole_adsorbates_cations=0.0;
  energy_bonddipole_framework_cations=0.0;
  fac=0.0;

  nr_of_coulombic_sites=nr_of_coulombic_sites_old=nr_of_coulombic_sites_new=0;
  nr_of_bonddipole_sites=nr_of_bonddipole_sites_old=nr_of_bonddipole_sites_new=0;

  if(OldMolecule)
  {
    nr_atoms=Cations[CurrentSystem][mol].NumberOfAtoms;
    for(j=0;j<nr_atoms;j++)
    {
      type=Cations[CurrentSystem][mol].Atoms[j].Type;
      charge=Cations[CurrentSystem][mol].Atoms[j].Charge;
      scaling=Cations[CurrentSystem][mol].Atoms[j].CFChargeScalingParameter;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      charge*=scaling;
      if(considered_charged)
      {
        Charge[nr_of_coulombic_sites]=charge;
        energy_self_old+=COULOMBIC_CONVERSION_FACTOR*SQR(charge)*alpha/sqrt(M_PI);
        net_charge_old+=Charge[nr_of_coulombic_sites];
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][mol].Atoms[j].Position);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_old=nr_of_coulombic_sites;

  if(NewMolecule)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      type=Components[CurrentComponent].Type[i];
      if(PseudoAtoms[type].HasCharges)
      {
        scaling=CFChargeScaling[i];
        Charge[nr_of_coulombic_sites]=scaling*Components[CurrentComponent].Charge[i];
        energy_self_new+=COULOMBIC_CONVERSION_FACTOR*SQR(Charge[nr_of_coulombic_sites])*alpha/sqrt(M_PI);
        net_charge_new+=Charge[nr_of_coulombic_sites];
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(TrialPosition[CurrentSystem][i]);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_new=nr_of_coulombic_sites;

  if(OldMolecule)
  {
    type_mol=Cations[CurrentSystem][mol].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Cations[CurrentSystem][mol].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[j];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      energy_self_bd_old+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_old=nr_of_bonddipole_sites;

  if(NewMolecule)
  {
    type_mol=CurrentComponent;
    for(i=0;i<Components[type_mol].NumberOfBondDipoles;i++)
    {
      pair=Components[type_mol].BondDipoles[i];
      posA=TrialPosition[CurrentSystem][pair.A];
      posB=TrialPosition[CurrentSystem][pair.B];
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[i];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      energy_self_bd_new+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_new=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }


  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          sum_old.re=0.0;
          sum_old.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_old.re+=temp*Eikr[i].re;
            sum_old.im+=temp*Eikr[i].im;
          }
          sum_new.re=0.0;
          sum_new.im=0.0;
          for(;i<nr_of_coulombic_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_new.re+=temp*Eikr[i].re;
            sum_new.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_old.re=0.0;
          sum_bonddipole_old.im=0.0;
          rk=kvecs[nvec];
          for(i=0;i<nr_of_bonddipole_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_old.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_old.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_new.re=0.0;
          sum_bonddipole_new.im=0.0;
          for(;i<nr_of_bonddipole_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_new.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_new.im+=temp*Eikr_bd[i].im;
          }

          sum_cations.re=StoreTotalChargeCations[CurrentSystem][nvec].re+(sum_new.re-sum_old.re);
          sum_cations.im=StoreTotalChargeCations[CurrentSystem][nvec].im+(sum_new.im-sum_old.im);

          sum_bonddipole_cations.re=StoreTotalBondDipolesCations[CurrentSystem][nvec].re+(sum_bonddipole_new.re-sum_bonddipole_old.re);
          sum_bonddipole_cations.im=StoreTotalBondDipolesCations[CurrentSystem][nvec].im+(sum_bonddipole_new.im-sum_bonddipole_old.im);

          temp=kfactor[nvec];

          // compute energy differences using the stored total sums and the sum of the differences of the moving atoms
          if(Framework[CurrentSystem].FrameworkModel!=NONE)
          {
            energy_charge_framework_cations+=temp*
              (StoreTotalChargeFramework[CurrentSystem][nvec].re*(sum_cations.re-StoreTotalChargeCations[CurrentSystem][nvec].re)
              +StoreTotalChargeFramework[CurrentSystem][nvec].im*(sum_cations.im-StoreTotalChargeCations[CurrentSystem][nvec].im));

            energy_charge_bonddipole_framework_cations+=temp*
              (StoreTotalChargeFramework[CurrentSystem][nvec].im*(sum_bonddipole_cations.re-StoreTotalBondDipolesCations[CurrentSystem][nvec].re)
              +StoreTotalChargeFramework[CurrentSystem][nvec].re*(StoreTotalBondDipolesCations[CurrentSystem][nvec].im-sum_bonddipole_cations.im)
              +StoreTotalBondDipolesFramework[CurrentSystem][nvec].re*(sum_cations.im-StoreTotalChargeCations[CurrentSystem][nvec].im)
              +StoreTotalBondDipolesFramework[CurrentSystem][nvec].im*(StoreTotalChargeCations[CurrentSystem][nvec].re-sum_cations.re));

            energy_bonddipole_framework_cations+=temp*
              (StoreTotalBondDipolesFramework[CurrentSystem][nvec].re*(sum_bonddipole_cations.re-StoreTotalBondDipolesCations[CurrentSystem][nvec].re)
              +StoreTotalBondDipolesFramework[CurrentSystem][nvec].im*(sum_bonddipole_cations.im-StoreTotalBondDipolesCations[CurrentSystem][nvec].im));
          }

          energy_charge_cations+=temp*(SQR(sum_cations.re)-SQR(StoreTotalChargeCations[CurrentSystem][nvec].re)+
                                       SQR(sum_cations.im)-SQR(StoreTotalChargeCations[CurrentSystem][nvec].im));

          energy_charge_bonddipole_cations+=2.0*temp*
                  (sum_cations.im*sum_bonddipole_cations.re-sum_cations.re*sum_bonddipole_cations.im
                  -StoreTotalChargeCations[CurrentSystem][nvec].im*StoreTotalBondDipolesCations[CurrentSystem][nvec].re
                  +StoreTotalChargeCations[CurrentSystem][nvec].re*StoreTotalBondDipolesCations[CurrentSystem][nvec].im);

          energy_bonddipole_cations+=temp*(SQR(sum_bonddipole_cations.re)-SQR(StoreTotalBondDipolesCations[CurrentSystem][nvec].re)+
                                           SQR(sum_bonddipole_cations.im)-SQR(StoreTotalBondDipolesCations[CurrentSystem][nvec].im));

          if(NumberOfAdsorbateMolecules[CurrentSystem]>0)
          {
            energy_charge_adsorbates_cations+=temp*
              (StoreTotalChargeAdsorbates[CurrentSystem][nvec].re*(sum_cations.re-StoreTotalChargeCations[CurrentSystem][nvec].re)
              +StoreTotalChargeAdsorbates[CurrentSystem][nvec].im*(sum_cations.im-StoreTotalChargeCations[CurrentSystem][nvec].im));

            energy_charge_bonddipole_adsorbates_cations+=temp*
              (StoreTotalChargeAdsorbates[CurrentSystem][nvec].im*(sum_bonddipole_cations.re-StoreTotalBondDipolesCations[CurrentSystem][nvec].re)
              +StoreTotalChargeAdsorbates[CurrentSystem][nvec].re*(StoreTotalBondDipolesCations[CurrentSystem][nvec].im-sum_bonddipole_cations.im)
              +StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re*(sum_cations.im-StoreTotalChargeCations[CurrentSystem][nvec].im)
              +StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im*(StoreTotalChargeCations[CurrentSystem][nvec].re-sum_cations.re));

            energy_bonddipole_adsorbates_cations+=temp*
              (StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re*(sum_bonddipole_cations.re-StoreTotalBondDipolesCations[CurrentSystem][nvec].re)
              +StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im*(sum_bonddipole_cations.im-StoreTotalBondDipolesCations[CurrentSystem][nvec].im));
          }

          // store the new sums, these will be the current ones on acceptance of the mc-move
          NewTotalChargeCations[store][nvec]=sum_cations;
          NewTotalBondDipolesCations[store][nvec]=sum_bonddipole_cations;

          // next wave-vector
          nvec++;
        }
      }
    }
  }

  energy_excluded_new=0.0;
  energy_excluded_c_bd_new=0.0;
  energy_excluded_bd_new=0.0;
  if(NewMolecule)
  {
    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraChargeCharge[i].A;
      B=Components[CurrentComponent].ExcludedIntraChargeCharge[i].B;
      scalingA=CFChargeScaling[A];
      scalingB=CFChargeScaling[B];
      chargeA=scalingA*Components[CurrentComponent].Charge[A];
      chargeB=scalingB*Components[CurrentComponent].Charge[B];
      posA=TrialPosition[CurrentSystem][A];
      posB=TrialPosition[CurrentSystem][B];

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      Bt0=-erf(alpha*r)/r;
      energy_excluded_new-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*Bt0;
    }

    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraChargeBondDipole[i].A;
      B=Components[CurrentComponent].ExcludedIntraChargeBondDipole[i].B;

      chargeA=Components[CurrentComponent].Charge[A];
      posA=TrialPosition[CurrentSystem][A];

      pair=Components[CurrentComponent].BondDipoles[B];
      posB1=TrialPosition[CurrentSystem][pair.A];
      posB2=TrialPosition[CurrentSystem][pair.B];

      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(alpha*r)/r;
      Bt1=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr+Bt0/rr;

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_c_bd_new+=COULOMBIC_CONVERSION_FACTOR*Bt1*chargeA*cosB;
    }


    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[CurrentComponent].ExcludedIntraBondDipoleBondDipole[i].B;

      pair=Components[CurrentComponent].BondDipoles[A];
      posA1=TrialPosition[CurrentSystem][pair.A];
      posA2=TrialPosition[CurrentSystem][pair.B];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      pair=Components[CurrentComponent].BondDipoles[B];
      posB1=TrialPosition[CurrentSystem][pair.A];
      posB2=TrialPosition[CurrentSystem][pair.B];
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      temp=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
      Bt0=-erf(alpha*r)/r;
      Bt1=temp+Bt0/rr;
      temp*=2.0*SQR(alpha);
      Bt2=temp+(3.0/rr)*Bt1;

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_bd_new-=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
    }
  }

  energy_excluded_old=0.0;
  energy_excluded_c_bd_old=0.0;
  energy_excluded_bd_old=0.0;
  if(OldMolecule)
  {
    type=Cations[CurrentSystem][mol].Type;

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      pair=Components[type].ExcludedIntraChargeCharge[i];
      scalingA=Cations[CurrentSystem][mol].Atoms[pair.A].CFChargeScalingParameter;
      scalingB=Cations[CurrentSystem][mol].Atoms[pair.B].CFChargeScalingParameter;
      chargeA=scalingA*Cations[CurrentSystem][mol].Atoms[pair.A].Charge;
      chargeB=scalingB*Cations[CurrentSystem][mol].Atoms[pair.B].Charge;
      posA=Cations[CurrentSystem][mol].Atoms[pair.A].Position;
      posB=Cations[CurrentSystem][mol].Atoms[pair.B].Position;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

      Bt0=-erf(alpha*r)/r;
      energy_excluded_old-=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*Bt0;
    }

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraChargeBondDipole[i].A;
      B=Components[type].ExcludedIntraChargeBondDipole[i].B;

      chargeA=Cations[CurrentSystem][mol].Atoms[A].Charge;
      posA=Cations[CurrentSystem][mol].Atoms[A].Position;

      pair=Components[type].BondDipoles[B];
      posB1=Cations[CurrentSystem][mol].Atoms[pair.A].Position;
      posB2=Cations[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(alpha*r)/r;
      Bt1=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr+Bt0/rr;

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_c_bd_old+=COULOMBIC_CONVERSION_FACTOR*Bt1*chargeA*cosB;
    }

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[type].ExcludedIntraBondDipoleBondDipole[i].B;

      pair=Components[type].BondDipoles[A];
      posA1=Cations[CurrentSystem][mol].Atoms[pair.A].Position;
      posA2=Cations[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=Components[type].BondDipoleMagnitude[A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      pair=Components[type].BondDipoles[B];
      posB1=Cations[CurrentSystem][mol].Atoms[pair.A].Position;
      posB2=Cations[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      temp=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
      Bt0=-erf(alpha*r)/r;
      Bt1=temp+Bt0/rr;
      temp*=2.0*SQR(alpha);
      Bt2=temp+(3.0/rr)*Bt1;

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy_excluded_bd_old-=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
    }
  }

  // set net-charge difference
  NetChargeCationDelta=net_charge_new-net_charge_old;

  // set energy differences
  if(!OmitCationCationCoulombInteractions)
  {
    UCationCationChargeChargeFourierDelta[CurrentSystem]=energy_charge_cations-(energy_excluded_new-energy_excluded_old)-(energy_self_new-energy_self_old);
    UCationCationChargeBondDipoleFourierDelta[CurrentSystem]=energy_charge_bonddipole_cations-(energy_excluded_c_bd_new-energy_excluded_c_bd_old);
    UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=energy_bonddipole_cations-(energy_self_bd_new-energy_self_bd_old)-
                                     (energy_excluded_bd_new-energy_excluded_bd_old);
    UCationCationChargeChargeFourierDelta[CurrentSystem]+=UIon[CurrentSystem]*SQR(NetChargeCations[CurrentSystem]+NetChargeCationDelta)-
                                                          UIon[CurrentSystem]*SQR(NetChargeCations[CurrentSystem]);
  }

  if(!OmitInterMolecularInteractions)
  {
    UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]=2.0*energy_charge_adsorbates_cations;
    UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]=2.0*energy_charge_bonddipole_adsorbates_cations;
    UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=2.0*energy_bonddipole_adsorbates_cations;
    UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeAdsorbates[CurrentSystem]*(NetChargeCations[CurrentSystem]+NetChargeCationDelta)-
                                                             2.0*UIon[CurrentSystem]*NetChargeAdsorbates[CurrentSystem]*(NetChargeCations[CurrentSystem]);
  }

  UHostCationChargeChargeFourierDelta[CurrentSystem]=2.0*energy_charge_framework_cations;
  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]=2.0*energy_charge_bonddipole_framework_cations;
  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=2.0*energy_bonddipole_framework_cations;
  UHostCationChargeChargeFourierDelta[CurrentSystem]+=2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*(NetChargeCations[CurrentSystem]+net_charge_new)-
                                                      2.0*UIon[CurrentSystem]*NetChargeFramework[CurrentSystem]*(NetChargeCations[CurrentSystem]+net_charge_old);

  return 0;
}

/*********************************************************************************************************
 * Name       | AcceptEwaldAdsorbateMove, AcceptEwaldCationMove                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Copies the new structure factors of the adsorbates or cations to the current ones.       *
 * Parameters |                                                                                          *
 * Note       |                                                                                          *
 *********************************************************************************************************/

void AcceptEwaldAdsorbateMove(int A)
{
  COMPLEX *temp_complex_pointer;

  SWAP(StoreTotalChargeAdsorbates[CurrentSystem],NewTotalChargeAdsorbates[A],temp_complex_pointer);
  SWAP(StoreTotalBondDipolesAdsorbates[CurrentSystem],NewTotalBondDipolesAdsorbates[A],temp_complex_pointer);
}

void AcceptEwaldCationMove(int A)
{
  COMPLEX *temp_complex_pointer;

  SWAP(StoreTotalChargeCations[CurrentSystem],NewTotalChargeCations[A],temp_complex_pointer);
  SWAP(StoreTotalBondDipolesCations[CurrentSystem],NewTotalBondDipolesCations[A],temp_complex_pointer);
}

/*********************************************************************************************************
 * Name       | SaveCurrentKVectors, RetrieveStoredKVectors                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Saves/restores the current k-vectors and pre-factors.                                    *
 * Parameters |                                                                                          *
 * Note       | A is either 0 or 1, B can be any system.                                                 *
 *********************************************************************************************************/

void SaveCurrentKVectors(int A,int B)
{
  memcpy(StoredKVectors[A],KVectors[B],MaxNumberOfWaveVectors*sizeof(VECTOR));
  memcpy(StoredKFactor[A],KFactor[B],MaxNumberOfWaveVectors*sizeof(REAL));
}

void RetrieveStoredKVectors(int A,int B)
{
  VECTOR *vector_pointer;
  REAL *real_pointer;

  SWAP(StoredKVectors[A],KVectors[B],vector_pointer);
  SWAP(StoredKFactor[A],KFactor[B],real_pointer);
}

/*********************************************************************************************************
 * Name       | SaveCurrentEwaldStructureFactors, RetrieveStoredEwaldStructureFactors                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Saves/restores the current structure factors.                                            *
 * Parameters |                                                                                          *
 * Note       | A is either 0 or 1, B can be any system.                                                 *
 *********************************************************************************************************/

void SaveCurrentEwaldStructureFactors(int A,int B)
{
  // memcpy copys the (non-overlapping) data at the second pointer to the first pointer
  memcpy(NewTotalChargeFramework[A],StoreTotalChargeFramework[B],MaxNumberOfWaveVectors*sizeof(COMPLEX));
  memcpy(NewTotalChargeAdsorbates[A],StoreTotalChargeAdsorbates[B],MaxNumberOfWaveVectors*sizeof(COMPLEX));
  memcpy(NewTotalChargeCations[A],StoreTotalChargeCations[B],MaxNumberOfWaveVectors*sizeof(COMPLEX));

  memcpy(NewTotalBondDipolesFramework[A],StoreTotalBondDipolesFramework[B],MaxNumberOfWaveVectors*sizeof(COMPLEX));
  memcpy(NewTotalBondDipolesAdsorbates[A],StoreTotalBondDipolesAdsorbates[B],MaxNumberOfWaveVectors*sizeof(COMPLEX));
  memcpy(NewTotalBondDipolesCations[A],StoreTotalBondDipolesCations[B],MaxNumberOfWaveVectors*sizeof(COMPLEX));
}

void RetrieveStoredEwaldStructureFactors(int A,int B)
{
  COMPLEX *complex_pointer;

  SWAP(NewTotalChargeFramework[A],StoreTotalChargeFramework[B],complex_pointer);
  SWAP(NewTotalChargeAdsorbates[A],StoreTotalChargeAdsorbates[B],complex_pointer);
  SWAP(NewTotalChargeCations[A],StoreTotalChargeCations[B],complex_pointer);

  SWAP(NewTotalBondDipolesFramework[A],StoreTotalBondDipolesFramework[B],complex_pointer);
  SWAP(NewTotalBondDipolesAdsorbates[A],StoreTotalBondDipolesAdsorbates[B],complex_pointer);
  SWAP(NewTotalBondDipolesCations[A],StoreTotalBondDipolesCations[B],complex_pointer);
}

/*********************************************************************************************************
 * Name       | SwapEwaldSystem                                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Swaps the stored Ewald information for two systems.                                      *
 * Parameters |                                                                                          *
 * Note       | A not equal B can be any two systems.                                                    *
 *            | Used in e.g. Parallel tempering.                                                         *
 *********************************************************************************************************/

void SwapEwaldSystem(int A,int B)
{
  VECTOR *vector_pointer;
  INT_VECTOR3 int_vector_temp;
  COMPLEX *complex_pointer;
  REAL *real_pointer,real_temp;
  int int_temp;

  if (A==B) return;

  SWAP(NumberOfKVectors[A],NumberOfKVectors[B],int_temp);
  SWAP(Alpha[A],Alpha[B],real_temp);
  SWAP(ReciprocalCutOffSquared[A],ReciprocalCutOffSquared[B],real_temp);
  SWAP(kvec[A],kvec[B],int_vector_temp);

  SWAP(KVectors[A],KVectors[B],vector_pointer);
  SWAP(KFactor[A],KFactor[B],real_pointer);

  SWAP(UChargeChargeFrameworkRigid[A],UChargeChargeFrameworkRigid[B],real_temp);
  SWAP(UChargeBondDipoleFrameworkRigid[A],UChargeBondDipoleFrameworkRigid[B],real_temp);
  SWAP(UBondDipoleBondDipoleFrameworkRigid[A],UBondDipoleBondDipoleFrameworkRigid[B],real_temp);

  SWAP(StoreRigidChargeFramework[A],StoreRigidChargeFramework[B],complex_pointer);
  SWAP(StoreRigidChargeCations[A],StoreRigidChargeCations[B],complex_pointer);
  SWAP(StoreRigidChargeAdsorbates[A],StoreRigidChargeAdsorbates[B],complex_pointer);

  SWAP(StoreRigidBondDipolesFramework[A],StoreRigidBondDipolesFramework[B],complex_pointer);
  SWAP(StoreRigidBondDipolesCations[A],StoreRigidBondDipolesCations[B],complex_pointer);
  SWAP(StoreRigidBondDipolesAdsorbates[A],StoreRigidBondDipolesAdsorbates[B],complex_pointer);

  SWAP(StoreTotalChargeFramework[A],StoreTotalChargeFramework[B],complex_pointer);
  SWAP(StoreTotalChargeCations[A],StoreTotalChargeCations[B],complex_pointer);
  SWAP(StoreTotalChargeAdsorbates[A],StoreTotalChargeAdsorbates[B],complex_pointer);

  SWAP(StoreTotalBondDipolesFramework[A],StoreTotalBondDipolesFramework[B],complex_pointer);
  SWAP(StoreTotalBondDipolesCations[A],StoreTotalBondDipolesCations[B],complex_pointer);
  SWAP(StoreTotalBondDipolesAdsorbates[A],StoreTotalBondDipolesAdsorbates[B],complex_pointer);
}

/*********************************************************************************************************
 * Name       | CalculateEwaldFourierFrameworkDisplacement                                               *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   |                                                                                          *
 * Parameters |                                                                                          *
 * Note       |                                                                                          *
 *********************************************************************************************************/

int CalculateEwaldFourierFrameworkDisplacement(void)
{
  int i,j;
  int nr_atoms,type;
  VECTOR pos;
  int nvec;
  REAL NetChargeNew,NetChargeOld;
  REAL USelfSumBondDipolesFrameworkNew,USelfSumBondDipolesFrameworkOld;
  REAL temp,charge;
  VECTOR posA1,posA2;
  int A1,A2;
  REAL energy_charge_framework,energy_charge_bonddipole_framework,energy_bonddipole_framework;
  REAL energy_charge_framework_adsorbates,energy_charge_bonddipole_framework_adsorbates,energy_bonddipole_framework_adsorbates;
  REAL energy_charge_framework_cations,energy_charge_bonddipole_framework_cations,energy_bonddipole_framework_cations;
  COMPLEX sum,sum_bonddipole;
  COMPLEX sum_new,sum_old;
  COMPLEX sum_bonddipole_new,sum_bonddipole_old;
  REAL alpha,*kfactor,recip_cutoff,ksqr;
  VECTOR *kvecs;
  int ii,jj,kk,index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_new,nr_of_coulombic_sites_old;
  int nr_of_bonddipole_sites,nr_of_bonddipole_sites_new,nr_of_bonddipole_sites_old;
  VECTOR rk,dipole;
  int considered_charged;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  energy_charge_framework=energy_charge_bonddipole_framework=energy_bonddipole_framework=0.0;
  energy_charge_framework_adsorbates=energy_charge_bonddipole_framework_adsorbates=energy_bonddipole_framework_adsorbates=0.0;
  energy_charge_framework_cations=energy_charge_bonddipole_framework_cations=energy_bonddipole_framework_cations=0.0;

  UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]=0.0;
  UHostCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  UCationCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  NetChargeAdsorbateDelta=0.0;

  USelfSumBondDipolesFrameworkNew=0.0;
  USelfSumBondDipolesFrameworkOld=0.0;

  nr_atoms=Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];

  nvec=0;
  NetChargeNew=NetChargeOld=0.0;

  nr_of_coulombic_sites=0;
  for(i=0;i<nr_atoms;i++)
  {
    type=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    charge=Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge;
    considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
    if(considered_charged)
    {
      Charge[nr_of_coulombic_sites]=charge;
      Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[CurrentFramework][i].Position);
      Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
      nr_of_coulombic_sites++;
    }
  }
  nr_of_coulombic_sites_new=nr_of_coulombic_sites;

  for(i=0;i<nr_atoms;i++)
  {
    type=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    charge=Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge;
    considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
    if(considered_charged)
    {
      Charge[nr_of_coulombic_sites]=charge;
      Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[CurrentFramework][i].ReferencePosition);
      Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
      nr_of_coulombic_sites++;
    }
  }
  nr_of_coulombic_sites_old=nr_of_coulombic_sites;

  nr_of_bonddipole_sites=0;
  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
    dipole.x=posA2.x-posA1.x;
    dipole.y=posA2.y-posA1.y;
    dipole.z=posA2.z-posA1.z;
    dipole=ApplyBoundaryCondition(dipole);
    pos.x=posA1.x+0.5*dipole.x;
    pos.y=posA1.y+0.5*dipole.y;
    pos.z=posA1.z+0.5*dipole.z;
    temp=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
    DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
    DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
    DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
    USelfSumBondDipolesFrameworkNew+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
    BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
    BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
    BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
    BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
    nr_of_bonddipole_sites++;
  }
  nr_of_bonddipole_sites_new=nr_of_bonddipole_sites;

  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
    posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
    dipole.x=posA2.x-posA1.x;
    dipole.y=posA2.y-posA1.y;
    dipole.z=posA2.z-posA1.z;
    dipole=ApplyBoundaryCondition(dipole);
    pos.x=posA1.x+0.5*dipole.x;
    pos.y=posA1.y+0.5*dipole.y;
    pos.z=posA1.z+0.5*dipole.z;
    temp=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
    DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
    DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
    DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
    USelfSumBondDipolesFrameworkNew+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
    BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
    BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
    BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
    BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
    nr_of_bonddipole_sites++;
  }
  nr_of_bonddipole_sites_old=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          sum_new.re=0.0;
          sum_new.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_new.re+=temp*Eikr[i].re;
            sum_new.im+=temp*Eikr[i].im;
          }

          sum_old.re=0.0;
          sum_old.im=0.0;
          for(;i<nr_of_coulombic_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_old.re+=temp*Eikr[i].re;
            sum_old.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_new.re=0.0;
          sum_bonddipole_new.im=0.0;
          rk=kvecs[nvec];
          for(i=0;i<nr_of_bonddipole_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_new.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_new.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_old.re=0.0;
          sum_bonddipole_old.im=0.0;
          for(;i<nr_of_bonddipole_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_old.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_old.im+=temp*Eikr_bd[i].im;
          }

          sum.re=StoreTotalChargeFramework[CurrentSystem][nvec].re+(sum_new.re-sum_old.re);
          sum.im=StoreTotalChargeFramework[CurrentSystem][nvec].im+(sum_new.im-sum_old.im);

          sum_bonddipole.re=StoreTotalBondDipolesFramework[CurrentSystem][nvec].re+(sum_bonddipole_new.re-sum_bonddipole_old.re);
          sum_bonddipole.im=StoreTotalBondDipolesFramework[CurrentSystem][nvec].im+(sum_bonddipole_new.im-sum_bonddipole_old.im);

          temp=kfactor[nvec];

          // energy differences
          energy_charge_framework+=temp*(SQR(sum.re)-SQR(StoreTotalChargeFramework[CurrentSystem][nvec].re)+
                                         SQR(sum.im)-SQR(StoreTotalChargeFramework[CurrentSystem][nvec].im));

          energy_charge_bonddipole_framework+=2.0*temp*
                  (sum.im*sum_bonddipole.re-sum.re*sum_bonddipole.im
                  -StoreTotalChargeFramework[CurrentSystem][nvec].im*StoreTotalBondDipolesFramework[CurrentSystem][nvec].re
                  +StoreTotalChargeFramework[CurrentSystem][nvec].re*StoreTotalBondDipolesFramework[CurrentSystem][nvec].im);

          energy_bonddipole_framework+=temp*(SQR(sum_bonddipole.re)-SQR(StoreTotalBondDipolesFramework[CurrentSystem][nvec].re)+
                                             SQR(sum_bonddipole.im)-SQR(StoreTotalBondDipolesFramework[CurrentSystem][nvec].im));

          if(NumberOfAdsorbateMolecules[CurrentSystem]>0)
          {
            energy_charge_framework_adsorbates+=temp*
              (StoreTotalChargeAdsorbates[CurrentSystem][nvec].re*(sum.re-StoreTotalChargeFramework[CurrentSystem][nvec].re)
              +StoreTotalChargeAdsorbates[CurrentSystem][nvec].im*(sum.im-StoreTotalChargeFramework[CurrentSystem][nvec].im));

            energy_charge_bonddipole_framework_adsorbates+=temp*
              (StoreTotalChargeAdsorbates[CurrentSystem][nvec].im*(sum_bonddipole.re-StoreTotalBondDipolesFramework[CurrentSystem][nvec].re)
              +StoreTotalChargeAdsorbates[CurrentSystem][nvec].re*(StoreTotalBondDipolesFramework[CurrentSystem][nvec].im-sum_bonddipole.im)
              +StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re*(sum.im-StoreTotalChargeFramework[CurrentSystem][nvec].im)
              +StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im*(StoreTotalChargeFramework[CurrentSystem][nvec].re-sum.re));

            energy_bonddipole_framework_adsorbates+=temp*
              (StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re*(sum_bonddipole.re-StoreTotalBondDipolesFramework[CurrentSystem][nvec].re)
              +StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im*(sum_bonddipole.im-StoreTotalBondDipolesFramework[CurrentSystem][nvec].im));
          }

          if(NumberOfCationMolecules[CurrentSystem]>0)
          {
            energy_charge_framework_cations+=temp*
              (StoreTotalChargeCations[CurrentSystem][nvec].re*(sum.re-StoreTotalChargeFramework[CurrentSystem][nvec].re)
              +StoreTotalChargeCations[CurrentSystem][nvec].im*(sum.im-StoreTotalChargeFramework[CurrentSystem][nvec].im));

            energy_charge_bonddipole_framework_cations+=temp*
              (StoreTotalChargeCations[CurrentSystem][nvec].im*(sum_bonddipole.re-StoreTotalBondDipolesFramework[CurrentSystem][nvec].re)
              +StoreTotalChargeCations[CurrentSystem][nvec].re*(StoreTotalBondDipolesFramework[CurrentSystem][nvec].im-sum_bonddipole.im)
              +StoreTotalBondDipolesCations[CurrentSystem][nvec].re*(sum.im-StoreTotalChargeFramework[CurrentSystem][nvec].im)
              +StoreTotalBondDipolesCations[CurrentSystem][nvec].im*(StoreTotalChargeFramework[CurrentSystem][nvec].re-sum.re));

            energy_bonddipole_framework_cations+=temp*
              (StoreTotalBondDipolesCations[CurrentSystem][nvec].re*(sum_bonddipole.re-StoreTotalBondDipolesFramework[CurrentSystem][nvec].re)
              +StoreTotalBondDipolesCations[CurrentSystem][nvec].im*(sum_bonddipole.im-StoreTotalBondDipolesFramework[CurrentSystem][nvec].im));
          }

          // Store Cksum of the Non-froozen particles
          NewTotalChargeFramework[0][nvec]=sum;
          NewTotalBondDipolesFramework[0][nvec]=sum_bonddipole;

          // if the framework is rigid, then also the structure factors of the fixed parts needs to be updated
          if(Framework[CurrentSystem].FrameworkModels[CurrentFramework]!=FLEXIBLE)
          {
            sum.re=StoreRigidChargeFramework[CurrentSystem][nvec].re+(sum_new.re-sum_old.re);
            sum.im=StoreRigidChargeFramework[CurrentSystem][nvec].im+(sum_new.im-sum_old.im);

            sum_bonddipole.re=StoreRigidBondDipolesFramework[CurrentSystem][nvec].re+(sum_bonddipole_new.re-sum_bonddipole_old.re);
            sum_bonddipole.im=StoreRigidBondDipolesFramework[CurrentSystem][nvec].im+(sum_bonddipole_new.im-sum_bonddipole_old.im);

            NewTotalChargeFramework[1][nvec]=sum;
            NewTotalBondDipolesFramework[1][nvec]=sum_bonddipole;
          }

          // next wavevector
          nvec++;
        }
      }
    }
  }

  UHostHostChargeChargeFourierDelta[CurrentSystem]=energy_charge_framework;
  UHostHostChargeBondDipoleFourierDelta[CurrentSystem]=energy_charge_bonddipole_framework;
  UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem]=energy_bonddipole_framework-(USelfSumBondDipolesFrameworkNew-USelfSumBondDipolesFrameworkOld);

  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]=2.0*energy_charge_framework_adsorbates;
  UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=2.0*energy_charge_bonddipole_framework_adsorbates;
  UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=2.0*energy_bonddipole_framework_adsorbates;

  UHostCationChargeChargeFourierDelta[CurrentSystem]=2.0*energy_charge_framework_cations;
  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]=2.0*energy_charge_bonddipole_framework_cations;
  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=2.0*energy_bonddipole_framework_cations;

  return 0;
}

/*********************************************************************************************************
 * Name       | AcceptEwaldFrameworkDiplacementMove, AcceptEwaldFrameworkDiplacementMoveRigid            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Updates the stored framework Ewald structure factors after a mc-move acceptance.         *
 * Parameters |                                                                                          *
 * Note       |                                                                                          *
 *********************************************************************************************************/

void AcceptEwaldFrameworkDiplacementMove(void)
{
  COMPLEX *complex_pointer;

  SWAP(StoreTotalChargeFramework[CurrentSystem],NewTotalChargeFramework[0],complex_pointer);
  SWAP(StoreTotalBondDipolesFramework[CurrentSystem],NewTotalBondDipolesFramework[0],complex_pointer);
}

void AcceptEwaldFrameworkDiplacementMoveRigid(void)
{
  COMPLEX *complex_pointer;

  SWAP(StoreRigidChargeFramework[CurrentSystem],NewTotalChargeFramework[1],complex_pointer);
  SWAP(StoreRigidBondDipolesFramework[CurrentSystem],NewTotalBondDipolesFramework[1],complex_pointer);
}



/*********************************************************************************************************
 * Name       | CalculateEwaldFourierFrameworkMoveDifference                                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier energy difference for a displaced framework atom.                   *
 * Parameters |                                                                                          *
 * Note       |                                                                                          *
 *********************************************************************************************************/

int CalculateEwaldFourierFrameworkAtomTranslate(int atom_id)
{
  int i,j;
  int Type;
  REAL fac;
  VECTOR Ss,pos;
  COMPLEX sum,sum_bonddipole;
  int nvec;
  REAL NetChargeNew,NetChargeOld;
  REAL USelfSumBondDipolesFrameworkNew,USelfSumBondDipolesFrameworkOld;
  REAL DipoleMagnitudeA,temp,r;
  VECTOR posA,posA1,posA2,dipoleA,dr,posB;
  int A1,A2,index,A,B;
  REAL EwaldIntraCorrectionChargeChargeNew,EwaldIntraCorrectionChargeChargeOld,chargeA,chargeB;
  REAL EwaldIntraCorrectionChargeBondDipoleNew,EwaldIntraCorrectionChargeBondDipoleOld;
  REAL EwaldIntraCorrectionBondDipoleBondDipoleNew,EwaldIntraCorrectionBondDipoleBondDipoleOld;
  REAL energy_charge_framework,energy_charge_bonddipole_framework,energy_bonddipole_framework;
  REAL energy_charge_framework_adsorbates,energy_charge_bonddipole_framework_adsorbates,energy_bonddipole_framework_adsorbates;
  REAL energy_charge_framework_cations,energy_charge_bonddipole_framework_cations,energy_bonddipole_framework_cations;
  REAL cosA,cosB,cosAB,rr,DipoleMagnitudeB,Bt1,Bt2,ChargeA;
  VECTOR posB1,posB2,dipoleB;
  int B1,B2,TypeA;

  int nr_of_bonddipole_sites_new,nr_of_bonddipole_sites_old,nr_of_bonddipole_sites,nr_of_coulombic_sites;
  int kmax_x,kmax_y,kmax_z;
  int ii,jj,kk,index_i,index_j,index_k;
  VECTOR rk,dipole;
  COMPLEX sum_new,sum_old;
  COMPLEX sum_bonddipole_new,sum_bonddipole_old;
  REAL alpha,*kfactor,recip_cutoff,ksqr;
  VECTOR *kvecs;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  energy_charge_framework=energy_charge_bonddipole_framework=energy_bonddipole_framework=0.0;
  energy_charge_framework_adsorbates=energy_charge_bonddipole_framework_adsorbates=energy_bonddipole_framework_adsorbates=0.0;
  energy_charge_framework_cations=energy_charge_bonddipole_framework_cations=energy_bonddipole_framework_cations=0.0;
  fac=0.0;
  Ss.x=Ss.y=Ss.z=0.0;

  UHostHostChargeChargeFourierDelta[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]=0.0;
  UHostCationChargeChargeFourierDelta[CurrentSystem]=0.0;
  UCationCationChargeChargeFourierDelta[CurrentSystem]=0.0;

  UHostHostChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourierDelta[CurrentSystem]=0.0;

  UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=0.0;

  NetChargeAdsorbateDelta=0.0;

  USelfSumBondDipolesFrameworkNew=0.0;
  USelfSumBondDipolesFrameworkOld=0.0;

  nvec=0;
  index=0;
  NetChargeNew=NetChargeOld=0.0;

  nr_of_bonddipole_sites=0;
  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    if((A1==atom_id)||(A2==atom_id))
    {
      posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
      dipole.x=posA2.x-posA1.x;
      dipole.y=posA2.y-posA1.y;
      dipole.z=posA2.z-posA1.z;
      dipole=ApplyBoundaryCondition(dipole);
      pos.x=posA1.x+0.5*dipole.x;
      pos.y=posA1.y+0.5*dipole.y;
      pos.z=posA1.z+0.5*dipole.z;
      temp=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      USelfSumBondDipolesFrameworkNew+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_new=nr_of_bonddipole_sites;

  nr_of_bonddipole_sites_old=0;
  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
  {
    DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;

    if((A1==atom_id)||(A2==atom_id))
    {
      posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
      posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
      dipole.x=posA2.x-posA1.x;
      dipole.y=posA2.y-posA1.y;
      dipole.z=posA2.z-posA1.z;
      dipole=ApplyBoundaryCondition(dipole);
      pos.x=posA1.x+0.5*dipole.x;
      pos.y=posA1.y+0.5*dipole.y;
      pos.z=posA1.z+0.5*dipole.z;
      temp=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      USelfSumBondDipolesFrameworkOld+=COULOMBIC_CONVERSION_FACTOR*2.0*CUBE(alpha)*SQR(temp)*(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z))/(3.0*sqrt(M_PI));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;

      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_old=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  nr_of_coulombic_sites=2;

  Eikx[0].re=1.0; Eikx[0].im=0.0;
  Eiky[0].re=1.0; Eiky[0].im=0.0;
  Eikz[0].re=1.0; Eikz[0].im=0.0;

  Type=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Type;
  Charge[0]=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Charge;
  pos=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Position);
  pos.x*=TWO_PI; pos.y*=TWO_PI; pos.z*=TWO_PI;

  index_i=MaxNumberOfCoulombicSites;
  Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
  Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
  Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

  index_i=-MaxNumberOfCoulombicSites;
  Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
  Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
  Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);

  Eikx[1].re=1.0; Eikx[1].im=0.0;
  Eiky[1].re=1.0; Eiky[1].im=0.0;
  Eikz[1].re=1.0; Eikz[1].im=0.0;

  Type=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Type;
  Charge[1]=Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].Charge;
  pos=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[CurrentFramework][atom_id].ReferencePosition);
  pos.x*=TWO_PI; pos.y*=TWO_PI; pos.z*=TWO_PI;

  index_i=MaxNumberOfCoulombicSites+1;
  Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
  Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
  Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

  index_i=-MaxNumberOfCoulombicSites+1;
  Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
  Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
  Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
          index_k=kk*MaxNumberOfCoulombicSites;
          temp=Charge[0];
          sum_new.re=temp*(Eikr_xy[0].re*Eikz[index_k].re-Eikr_xy[0].im*Eikz[index_k].im);
          sum_new.im=temp*(Eikr_xy[0].im*Eikz[index_k].re+Eikr_xy[0].re*Eikz[index_k].im);

          index_k=kk*MaxNumberOfCoulombicSites+1;
          temp=Charge[1];
          sum_old.re=temp*(Eikr_xy[1].re*Eikz[index_k].re-Eikr_xy[1].im*Eikz[index_k].im);
          sum_old.im=temp*(Eikr_xy[1].im*Eikz[index_k].re+Eikr_xy[1].re*Eikz[index_k].im);

          sum_bonddipole_new.re=0.0;
          sum_bonddipole_new.im=0.0;
          rk=kvecs[nvec];
          for(i=0;i<nr_of_bonddipole_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_new.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_new.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_old.re=0.0;
          sum_bonddipole_old.im=0.0;
          for(;i<nr_of_bonddipole_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_old.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_old.im+=temp*Eikr_bd[i].im;
          }

          sum.re=StoreTotalChargeFramework[CurrentSystem][nvec].re+(sum_new.re-sum_old.re);
          sum.im=StoreTotalChargeFramework[CurrentSystem][nvec].im+(sum_new.im-sum_old.im);

          sum_bonddipole.re=StoreTotalBondDipolesFramework[CurrentSystem][nvec].re+(sum_bonddipole_new.re-sum_bonddipole_old.re);
          sum_bonddipole.im=StoreTotalBondDipolesFramework[CurrentSystem][nvec].im+(sum_bonddipole_new.im-sum_bonddipole_old.im);

          temp=kfactor[nvec];

          // energy differences
          energy_charge_framework+=temp*(SQR(sum.re)-SQR(StoreTotalChargeFramework[CurrentSystem][nvec].re)+
                                         SQR(sum.im)-SQR(StoreTotalChargeFramework[CurrentSystem][nvec].im));

          energy_charge_bonddipole_framework+=2.0*temp*
                  (sum.im*sum_bonddipole.re-sum.re*sum_bonddipole.im
                  -StoreTotalChargeFramework[CurrentSystem][nvec].im*StoreTotalBondDipolesFramework[CurrentSystem][nvec].re
                  +StoreTotalChargeFramework[CurrentSystem][nvec].re*StoreTotalBondDipolesFramework[CurrentSystem][nvec].im);

          energy_bonddipole_framework+=temp*(SQR(sum_bonddipole.re)-SQR(StoreTotalBondDipolesFramework[CurrentSystem][nvec].re)+
                                             SQR(sum_bonddipole.im)-SQR(StoreTotalBondDipolesFramework[CurrentSystem][nvec].im));

          if(NumberOfAdsorbateMolecules[CurrentSystem]>0)
          {
            energy_charge_framework_adsorbates+=temp*
              (StoreTotalChargeAdsorbates[CurrentSystem][nvec].re*(sum.re-StoreTotalChargeFramework[CurrentSystem][nvec].re)
              +StoreTotalChargeAdsorbates[CurrentSystem][nvec].im*(sum.im-StoreTotalChargeFramework[CurrentSystem][nvec].im));

            energy_charge_bonddipole_framework_adsorbates+=temp*
              (StoreTotalChargeAdsorbates[CurrentSystem][nvec].im*(sum_bonddipole.re-StoreTotalBondDipolesFramework[CurrentSystem][nvec].re)
              +StoreTotalChargeAdsorbates[CurrentSystem][nvec].re*(StoreTotalBondDipolesFramework[CurrentSystem][nvec].im-sum_bonddipole.im)
              +StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re*(sum.im-StoreTotalChargeFramework[CurrentSystem][nvec].im)
              +StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im*(StoreTotalChargeFramework[CurrentSystem][nvec].re-sum.re));

            energy_bonddipole_framework_adsorbates+=temp*
              (StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re*(sum_bonddipole.re-StoreTotalBondDipolesFramework[CurrentSystem][nvec].re)
              +StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im*(sum_bonddipole.im-StoreTotalBondDipolesFramework[CurrentSystem][nvec].im));
          }

          if(NumberOfCationMolecules[CurrentSystem]>0)
          {
            energy_charge_framework_cations+=temp*
              (StoreTotalChargeCations[CurrentSystem][nvec].re*(sum.re-StoreTotalChargeFramework[CurrentSystem][nvec].re)
              +StoreTotalChargeCations[CurrentSystem][nvec].im*(sum.im-StoreTotalChargeFramework[CurrentSystem][nvec].im));

            energy_charge_bonddipole_framework_cations+=temp*
              (StoreTotalChargeCations[CurrentSystem][nvec].im*(sum_bonddipole.re-StoreTotalBondDipolesFramework[CurrentSystem][nvec].re)
              +StoreTotalChargeCations[CurrentSystem][nvec].re*(StoreTotalBondDipolesFramework[CurrentSystem][nvec].im-sum_bonddipole.im)
              +StoreTotalBondDipolesCations[CurrentSystem][nvec].re*(sum.im-StoreTotalChargeFramework[CurrentSystem][nvec].im)
              +StoreTotalBondDipolesCations[CurrentSystem][nvec].im*(StoreTotalChargeFramework[CurrentSystem][nvec].re-sum.re));

            energy_bonddipole_framework_cations+=temp*
              (StoreTotalBondDipolesCations[CurrentSystem][nvec].re*(sum_bonddipole.re-StoreTotalBondDipolesFramework[CurrentSystem][nvec].re)
              +StoreTotalBondDipolesCations[CurrentSystem][nvec].im*(sum_bonddipole.im-StoreTotalBondDipolesFramework[CurrentSystem][nvec].im));
          }

          // Store Cksum of the Non-froozen particles
          NewTotalChargeFramework[0][nvec]=sum;
          NewTotalBondDipolesFramework[0][nvec]=sum_bonddipole;

          // next wavevector
          nvec++;
        }
      }
    }
  }

  UHostHostChargeChargeFourierDelta[CurrentSystem]=energy_charge_framework;
  UHostHostChargeBondDipoleFourierDelta[CurrentSystem]=energy_charge_bonddipole_framework;
  UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem]=energy_bonddipole_framework
             -(USelfSumBondDipolesFrameworkNew-USelfSumBondDipolesFrameworkOld);

  UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]=2.0*energy_charge_framework_adsorbates;
  UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]=2.0*energy_charge_bonddipole_framework_adsorbates;
  UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]=2.0*energy_bonddipole_framework_adsorbates;

  UHostCationChargeChargeFourierDelta[CurrentSystem]=2.0*energy_charge_framework_cations;
  UHostCationChargeBondDipoleFourierDelta[CurrentSystem]=2.0*energy_charge_bonddipole_framework_cations;
  UHostCationBondDipoleBondDipoleFourierDelta[CurrentSystem]=2.0*energy_bonddipole_framework_cations;


  EwaldIntraCorrectionChargeChargeOld=0.0;
  EwaldIntraCorrectionChargeChargeNew=0.0;
  for(i=0;i<Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[CurrentFramework];i++)
  {
    A=Framework[CurrentSystem].ExcludedIntraChargeCharge[CurrentFramework][i].A;
    B=Framework[CurrentSystem].ExcludedIntraChargeCharge[CurrentFramework][i].B;

    if((A==atom_id)||(B==atom_id))
    {
      chargeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Charge;
      chargeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Charge;

      posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
      posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      EwaldIntraCorrectionChargeChargeNew+=COULOMBIC_CONVERSION_FACTOR*erf(Alpha[CurrentSystem]*r)*
            chargeA*chargeB/r;

      posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].ReferencePosition;
      posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].ReferencePosition;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      EwaldIntraCorrectionChargeChargeOld+=COULOMBIC_CONVERSION_FACTOR*erf(Alpha[CurrentSystem]*r)*
            chargeA*chargeB/r;
    }
  }
  UHostHostChargeChargeFourierDelta[CurrentSystem]-=EwaldIntraCorrectionChargeChargeNew-EwaldIntraCorrectionChargeChargeOld;

  EwaldIntraCorrectionChargeBondDipoleOld=0.0;
  EwaldIntraCorrectionChargeBondDipoleNew=0.0;
  for(i=0;i<Framework[CurrentSystem].NumberOfExcludedIntraChargeBondDipole[CurrentFramework];i++)
  {
    A=Framework[CurrentSystem].ExcludedIntraChargeBondDipole[CurrentFramework][i].A;
    B=Framework[CurrentSystem].ExcludedIntraChargeBondDipole[CurrentFramework][i].B;
    B1=Framework[CurrentSystem].BondDipoles[CurrentFramework][B].A;
    B2=Framework[CurrentSystem].BondDipoles[CurrentFramework][B].B;

    if((A==atom_id)||(B1==atom_id)||(B2==atom_id))
    {
      TypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
      ChargeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Charge;
      DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][B];

      posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
      posB1=Framework[CurrentSystem].Atoms[CurrentFramework][B1].Position;
      posB2=Framework[CurrentSystem].Atoms[CurrentFramework][B2].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      dipoleB=ApplyBoundaryCondition(dipoleB);
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -erf(Alpha[CurrentSystem]*r)/(rr*r);

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      EwaldIntraCorrectionChargeBondDipoleNew-=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);

      posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].ReferencePosition;
      posB1=Framework[CurrentSystem].Atoms[CurrentFramework][B1].ReferencePosition;
      posB2=Framework[CurrentSystem].Atoms[CurrentFramework][B2].ReferencePosition;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      dipoleB=ApplyBoundaryCondition(dipoleB);
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=DipoleMagnitudeB/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -erf(Alpha[CurrentSystem]*r)/(rr*r);

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      EwaldIntraCorrectionChargeBondDipoleOld-=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
    }
  }
  UHostHostChargeBondDipoleFourierDelta[CurrentSystem]+=EwaldIntraCorrectionChargeBondDipoleNew-EwaldIntraCorrectionChargeBondDipoleOld;

  EwaldIntraCorrectionBondDipoleBondDipoleOld=0.0;
  EwaldIntraCorrectionBondDipoleBondDipoleNew=0.0;
  for(i=0;i<Framework[CurrentSystem].NumberOfExcludedIntraBondDipoleBondDipole[CurrentFramework];i++)
  {
    A=Framework[CurrentSystem].ExcludedIntraBondDipoleBondDipole[CurrentFramework][i].A;
    B=Framework[CurrentSystem].ExcludedIntraBondDipoleBondDipole[CurrentFramework][i].B;

    A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][A].A;
    A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][A].B;
    B1=Framework[CurrentSystem].BondDipoles[CurrentFramework][B].A;
    B2=Framework[CurrentSystem].BondDipoles[CurrentFramework][B].B;
    if((A1==atom_id)||(A2==atom_id)||(B1==atom_id)||(B2==atom_id))
    {
      //NEW
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][A];
      posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=DipoleMagnitudeA/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][B];
      posB1=Framework[CurrentSystem].Atoms[CurrentFramework][B1].Position;
      posB2=Framework[CurrentSystem].Atoms[CurrentFramework][B2].Position;
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
      r=sqrt(rr);

      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -erf(Alpha[CurrentSystem]*r)/(rr*r);
      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -3.0*erf(Alpha[CurrentSystem]*r)/(rr*rr*r);

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      EwaldIntraCorrectionBondDipoleBondDipoleNew+=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);


      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][A];
      posA1=Framework[CurrentSystem].Atoms[CurrentFramework][A1].ReferencePosition;
      posA2=Framework[CurrentSystem].Atoms[CurrentFramework][A2].ReferencePosition;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=DipoleMagnitudeA/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][B];
      posB1=Framework[CurrentSystem].Atoms[CurrentFramework][B1].ReferencePosition;
      posB2=Framework[CurrentSystem].Atoms[CurrentFramework][B2].ReferencePosition;
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
      r=sqrt(rr);

      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -erf(Alpha[CurrentSystem]*r)/(rr*r);
      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -3.0*erf(Alpha[CurrentSystem]*r)/(rr*rr*r);

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      EwaldIntraCorrectionBondDipoleBondDipoleOld+=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);

    }
  }
  UHostHostBondDipoleBondDipoleFourierDelta[CurrentSystem]+=EwaldIntraCorrectionBondDipoleBondDipoleNew-EwaldIntraCorrectionBondDipoleBondDipoleOld;
  return 0;
}

void AcceptEwaldFrameworkMove(int A)
{
  COMPLEX *complex_pointer;

  SWAP(StoreTotalChargeFramework[CurrentSystem],NewTotalChargeFramework[A],complex_pointer);
  SWAP(StoreTotalBondDipolesFramework[CurrentSystem],NewTotalBondDipolesFramework[A],complex_pointer);
}


/*********************************************************************************************************
 * Name       | EwaldFourierElectrostaticPotentialAdsorbate,EwaldFourierElectrostaticPotentialCation     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of the electrostatic potential.                                *
 * Parameters |                                                                                          *
 * Limitation | Only charge-charge interaction is taken into account (no bond-dipoles yet).              *
 * Todo       | Add bond-dipoles to the static electric field.                                           *
 * Problem    | A bond-dipole exerts forces on both the atoms -> E!=F/mu                                 *
 * Used       | To compute the electric field numerically and validate the implementation.               *
 *********************************************************************************************************/

REAL EwaldFourierElectrostaticPotentialAdsorbate(int m,int l)
{
  int i,j;
  VECTOR pos,posA;
  int nvec;
  REAL UElectrostaticPotential;
  int kmax_x,kmax_y,kmax_z;
  int ii,jj,kk,index_i,index_j,index_k;
  COMPLEX sum,sum_charges,sum_bondipoles;
  REAL alpha,*kfactor,recip_cutoff,ksqr;
  VECTOR *kvecs;
  int nr_of_excluded_pairs;
  REAL r,rr,chargeA,chargeB;
  int typeA,typeB,type;
  VECTOR posB,dr;
  PAIR pair;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  UElectrostaticPotential=0.0;
  nvec=0;

  posA=Adsorbates[CurrentSystem][m].Atoms[l].Position;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  Eikx[0].re=1.0; Eikx[0].im=0.0;
  Eiky[0].re=1.0; Eiky[0].im=0.0;
  Eikz[0].re=1.0; Eikz[0].im=0.0;

  pos=ConvertFromXYZtoABC(posA);
  pos.x*=TWO_PI; pos.y*=TWO_PI; pos.z*=TWO_PI;

  index_i=MaxNumberOfCoulombicSites;
  Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
  Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
  Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

  index_i=-MaxNumberOfCoulombicSites;
  Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
  Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
  Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
  {
    Eikx[j*MaxNumberOfCoulombicSites].re=Eikx[(j-1)*MaxNumberOfCoulombicSites].re*Eikx[MaxNumberOfCoulombicSites].re-
                                         Eikx[(j-1)*MaxNumberOfCoulombicSites].im*Eikx[MaxNumberOfCoulombicSites].im;
    Eikx[j*MaxNumberOfCoulombicSites].im=Eikx[(j-1)*MaxNumberOfCoulombicSites].im*Eikx[MaxNumberOfCoulombicSites].re+
                                         Eikx[(j-1)*MaxNumberOfCoulombicSites].re*Eikx[MaxNumberOfCoulombicSites].im;
  }

  for(j=2;j<=kmax_y;j++)
  {
    Eiky[j*MaxNumberOfCoulombicSites].re=Eiky[(j-1)*MaxNumberOfCoulombicSites].re*Eiky[MaxNumberOfCoulombicSites].re-
                                         Eiky[(j-1)*MaxNumberOfCoulombicSites].im*Eiky[MaxNumberOfCoulombicSites].im;
    Eiky[j*MaxNumberOfCoulombicSites].im=Eiky[(j-1)*MaxNumberOfCoulombicSites].im*Eiky[MaxNumberOfCoulombicSites].re+
                                         Eiky[(j-1)*MaxNumberOfCoulombicSites].re*Eiky[MaxNumberOfCoulombicSites].im;
    Eiky[-j*MaxNumberOfCoulombicSites].re=Eiky[j*MaxNumberOfCoulombicSites].re;
    Eiky[-j*MaxNumberOfCoulombicSites].im=-Eiky[j*MaxNumberOfCoulombicSites].im;
  }

  for(j=2;j<=kmax_z;j++)
  {
    Eikz[j*MaxNumberOfCoulombicSites].re=Eikz[(j-1)*MaxNumberOfCoulombicSites].re*Eikz[MaxNumberOfCoulombicSites].re-
                                         Eikz[(j-1)*MaxNumberOfCoulombicSites].im*Eikz[MaxNumberOfCoulombicSites].im;
    Eikz[j*MaxNumberOfCoulombicSites].im=Eikz[(j-1)*MaxNumberOfCoulombicSites].im*Eikz[MaxNumberOfCoulombicSites].re+
                                         Eikz[(j-1)*MaxNumberOfCoulombicSites].re*Eikz[MaxNumberOfCoulombicSites].im;
    Eikz[-j*MaxNumberOfCoulombicSites].re=Eikz[j*MaxNumberOfCoulombicSites].re;
    Eikz[-j*MaxNumberOfCoulombicSites].im=-Eikz[j*MaxNumberOfCoulombicSites].im;
  }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
      // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
      index_i=ii*MaxNumberOfCoulombicSites;
      index_j=jj*MaxNumberOfCoulombicSites;
      Eikr_xy[0].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
      Eikr_xy[0].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
          index_k=kk*MaxNumberOfCoulombicSites;
          sum.re=(Eikr_xy[0].re*Eikz[index_k].re-Eikr_xy[0].im*Eikz[index_k].im);
          sum.im=(Eikr_xy[0].im*Eikz[index_k].re+Eikr_xy[0].re*Eikz[index_k].im);

          sum_charges.re=StoreTotalChargeFramework[CurrentSystem][nvec].re+StoreTotalChargeAdsorbates[CurrentSystem][nvec].re+StoreTotalChargeCations[CurrentSystem][nvec].re;
          sum_charges.im=StoreTotalChargeFramework[CurrentSystem][nvec].im+StoreTotalChargeAdsorbates[CurrentSystem][nvec].im+StoreTotalChargeCations[CurrentSystem][nvec].im;
          sum_bondipoles.re=StoreTotalBondDipolesFramework[CurrentSystem][nvec].re+StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re+StoreTotalBondDipolesCations[CurrentSystem][nvec].re;
          sum_bondipoles.im=StoreTotalBondDipolesFramework[CurrentSystem][nvec].im+StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im+StoreTotalBondDipolesCations[CurrentSystem][nvec].im;

          // electrostatic potential due to charges and dipoles
          UElectrostaticPotential+=2.0*kfactor[nvec]*
             (sum.re*sum_charges.re+sum.im*sum_charges.im+sum.im*sum_bondipoles.re-sum.re*sum_bondipoles.im);
          if(OmitAdsorbateAdsorbatePolarization)
            UElectrostaticPotential-=2.0*kfactor[nvec]*
               (sum.re*StoreTotalChargeAdsorbates[CurrentSystem][nvec].re+sum.im*StoreTotalChargeAdsorbates[CurrentSystem][nvec].im+
                sum.im*StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re-sum.re*StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im);
          if(OmitAdsorbateCationPolarization)
            UElectrostaticPotential-=2.0*kfactor[nvec]*
               (sum.re*StoreTotalChargeCations[CurrentSystem][nvec].re+sum.im*StoreTotalChargeCations[CurrentSystem][nvec].im+
                sum.im*StoreTotalBondDipolesCations[CurrentSystem][nvec].re-sum.re*StoreTotalBondDipolesCations[CurrentSystem][nvec].im);

          // next wavevector
          nvec++;
        }
      }
    }
  }

  // intra-molecular charge-charge exclusion
  if(!OmitAdsorbateAdsorbatePolarization)
  {
    type=Adsorbates[CurrentSystem][m].Type;
    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      pair=Components[type].ExcludedIntraChargeCharge[i];
      if((pair.A==l)||(pair.B==l))
      {
        typeA=Adsorbates[CurrentSystem][m].Atoms[pair.A].Type;
        typeB=Adsorbates[CurrentSystem][m].Atoms[pair.B].Type;
        posA=Adsorbates[CurrentSystem][m].Atoms[pair.A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[pair.B].Position;
        chargeA=Adsorbates[CurrentSystem][m].Atoms[pair.A].Charge;
        chargeB=Adsorbates[CurrentSystem][m].Atoms[pair.B].Charge;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        if(pair.A==l)
          UElectrostaticPotential-=COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)*chargeB/r;
        else if(pair.B==l)
          UElectrostaticPotential-=COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)*chargeA/r;
      }
    }
  }

  UElectrostaticPotential-=COULOMBIC_CONVERSION_FACTOR*((2.0*alpha)/sqrt(M_PI))*Adsorbates[CurrentSystem][m].Atoms[l].Charge;

  return UElectrostaticPotential;
}

REAL EwaldFourierElectrostaticPotentialCation(int m,int l)
{
  int i,j;
  VECTOR pos,posA;
  int nvec;
  REAL UElectrostaticPotential;
  int kmax_x,kmax_y,kmax_z;
  int ii,jj,kk,index_i,index_j,index_k;
  COMPLEX sum,sum_charges,sum_bondipoles;
  REAL alpha,*kfactor,recip_cutoff,ksqr;
  VECTOR *kvecs;
  int nr_of_excluded_pairs;
  REAL r,rr,chargeA,chargeB;
  int typeA,typeB,type;
  VECTOR posB,dr;
  PAIR pair;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  UElectrostaticPotential=0.0;
  nvec=0;

  posA=Cations[CurrentSystem][m].Atoms[l].Position;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  Eikx[0].re=1.0; Eikx[0].im=0.0;
  Eiky[0].re=1.0; Eiky[0].im=0.0;
  Eikz[0].re=1.0; Eikz[0].im=0.0;

  pos=ConvertFromXYZtoABC(posA);
  pos.x*=TWO_PI; pos.y*=TWO_PI; pos.z*=TWO_PI;

  index_i=MaxNumberOfCoulombicSites;
  Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
  Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
  Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

  index_i=-MaxNumberOfCoulombicSites;
  Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
  Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
  Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
  {
    Eikx[j*MaxNumberOfCoulombicSites].re=Eikx[(j-1)*MaxNumberOfCoulombicSites].re*Eikx[MaxNumberOfCoulombicSites].re-
                                         Eikx[(j-1)*MaxNumberOfCoulombicSites].im*Eikx[MaxNumberOfCoulombicSites].im;
    Eikx[j*MaxNumberOfCoulombicSites].im=Eikx[(j-1)*MaxNumberOfCoulombicSites].im*Eikx[MaxNumberOfCoulombicSites].re+
                                         Eikx[(j-1)*MaxNumberOfCoulombicSites].re*Eikx[MaxNumberOfCoulombicSites].im;
  }

  for(j=2;j<=kmax_y;j++)
  {
    Eiky[j*MaxNumberOfCoulombicSites].re=Eiky[(j-1)*MaxNumberOfCoulombicSites].re*Eiky[MaxNumberOfCoulombicSites].re-
                                         Eiky[(j-1)*MaxNumberOfCoulombicSites].im*Eiky[MaxNumberOfCoulombicSites].im;
    Eiky[j*MaxNumberOfCoulombicSites].im=Eiky[(j-1)*MaxNumberOfCoulombicSites].im*Eiky[MaxNumberOfCoulombicSites].re+
                                         Eiky[(j-1)*MaxNumberOfCoulombicSites].re*Eiky[MaxNumberOfCoulombicSites].im;
    Eiky[-j*MaxNumberOfCoulombicSites].re=Eiky[j*MaxNumberOfCoulombicSites].re;
    Eiky[-j*MaxNumberOfCoulombicSites].im=-Eiky[j*MaxNumberOfCoulombicSites].im;
  }

  for(j=2;j<=kmax_z;j++)
  {
    Eikz[j*MaxNumberOfCoulombicSites].re=Eikz[(j-1)*MaxNumberOfCoulombicSites].re*Eikz[MaxNumberOfCoulombicSites].re-
                                         Eikz[(j-1)*MaxNumberOfCoulombicSites].im*Eikz[MaxNumberOfCoulombicSites].im;
    Eikz[j*MaxNumberOfCoulombicSites].im=Eikz[(j-1)*MaxNumberOfCoulombicSites].im*Eikz[MaxNumberOfCoulombicSites].re+
                                         Eikz[(j-1)*MaxNumberOfCoulombicSites].re*Eikz[MaxNumberOfCoulombicSites].im;
    Eikz[-j*MaxNumberOfCoulombicSites].re=Eikz[j*MaxNumberOfCoulombicSites].re;
    Eikz[-j*MaxNumberOfCoulombicSites].im=-Eikz[j*MaxNumberOfCoulombicSites].im;
  }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
      // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
      index_i=ii*MaxNumberOfCoulombicSites;
      index_j=jj*MaxNumberOfCoulombicSites;
      Eikr_xy[0].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
      Eikr_xy[0].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
          index_k=kk*MaxNumberOfCoulombicSites;
          sum.re=(Eikr_xy[0].re*Eikz[index_k].re-Eikr_xy[0].im*Eikz[index_k].im);
          sum.im=(Eikr_xy[0].im*Eikz[index_k].re+Eikr_xy[0].re*Eikz[index_k].im);

          sum_charges.re=StoreTotalChargeFramework[CurrentSystem][nvec].re+StoreTotalChargeAdsorbates[CurrentSystem][nvec].re+StoreTotalChargeCations[CurrentSystem][nvec].re;
          sum_charges.im=StoreTotalChargeFramework[CurrentSystem][nvec].im+StoreTotalChargeAdsorbates[CurrentSystem][nvec].im+StoreTotalChargeCations[CurrentSystem][nvec].im;
          sum_bondipoles.re=StoreTotalBondDipolesFramework[CurrentSystem][nvec].re+StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re+StoreTotalBondDipolesCations[CurrentSystem][nvec].re;
          sum_bondipoles.im=StoreTotalBondDipolesFramework[CurrentSystem][nvec].im+StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im+StoreTotalBondDipolesCations[CurrentSystem][nvec].im;

          // electrostatic potential due to charges and dipoles
          UElectrostaticPotential+=2.0*kfactor[nvec]*
             (sum.re*sum_charges.re+sum.im*sum_charges.im+sum.im*sum_bondipoles.re-sum.re*sum_bondipoles.im);
          if(OmitCationCationPolarization)
            UElectrostaticPotential-=2.0*kfactor[nvec]*
               (sum.re*StoreTotalChargeCations[CurrentSystem][nvec].re+sum.im*StoreTotalChargeCations[CurrentSystem][nvec].im+
                sum.im*StoreTotalBondDipolesCations[CurrentSystem][nvec].re-sum.re*StoreTotalBondDipolesCations[CurrentSystem][nvec].im);
          if(OmitAdsorbateCationPolarization)
            UElectrostaticPotential-=2.0*kfactor[nvec]*
               (sum.re*StoreTotalChargeAdsorbates[CurrentSystem][nvec].re+sum.im*StoreTotalChargeAdsorbates[CurrentSystem][nvec].im+
                sum.im*StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].re-sum.re*StoreTotalBondDipolesAdsorbates[CurrentSystem][nvec].im);

          // next wavevector
          nvec++;
        }
      }
    }
  }

  // intra-molecular charge-charge exclusion
  if(!OmitCationCationPolarization)
  {
    type=Cations[CurrentSystem][m].Type;
    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      pair=Components[type].ExcludedIntraChargeCharge[i];
      if((pair.A==l)||(pair.B==l))
      {
        typeA=Cations[CurrentSystem][m].Atoms[pair.A].Type;
        typeB=Cations[CurrentSystem][m].Atoms[pair.B].Type;
        posA=Cations[CurrentSystem][m].Atoms[pair.A].Position;
        posB=Cations[CurrentSystem][m].Atoms[pair.B].Position;
        chargeA=Cations[CurrentSystem][m].Atoms[pair.A].Charge;
        chargeB=Cations[CurrentSystem][m].Atoms[pair.B].Charge;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        if(pair.A==l)
          UElectrostaticPotential-=COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)*chargeB/r;
        else if(pair.B==l)
          UElectrostaticPotential-=COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)*chargeA/r;
      }
    }
  }

  UElectrostaticPotential-=COULOMBIC_CONVERSION_FACTOR*((2.0*alpha)/sqrt(M_PI))*Cations[CurrentSystem][m].Atoms[l].Charge;

  return UElectrostaticPotential;
}



/*********************************************************************************************************
 * Name       | EwaldFourierStaticElectricField                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of the static electric field.                                  *
 * Parameters |                                                                                          *
 * Note       | The electric field of the dipole is the negative gradient of the electrostatic potential.*
 *            | The electric-field is computed on atoms that are ether charged or have polarization.     *
 * Limitation | Only charge-charge interaction is taken into account (no bond-dipoles yet).              *
 * Todo       | Add bond-dipoles to the static electric field.                                           *
 * Problem    | A bond-dipole exerts forces on both the atoms -> E!=F/mu                                 *
 *********************************************************************************************************/

int EwaldFourierStaticElectricField(void)
{
  int i,j,m,f1,ii,jj,kk;
  int nvec,type_mol,type,typeA,typeB;
  REAL temp,alpha,r,rr;
  VECTOR pos,rk,dr,posA,posB;
  int nr_molecules,nr_atoms,nr_frameworks,nr_of_excluded_pairs;
  int index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_framework,nr_of_coulombic_sites_adsorbate,nr_of_coulombic_sites_cation;
  COMPLEX sum,sum_framework,sum_adsorbate,sum_cation;
  int nr_of_bonddipoles,nr_of_bonddipole_sites,nr_of_bonddipole_sites_framework,nr_of_bonddipole_sites_adsorbate,nr_of_bonddipole_site_cation;
  COMPLEX sum_bonddipole,sum_bonddipole_framework,sum_bonddipole_adsorbate,sum_bonddipole_cation;
  VECTOR dipole;
  REAL Bt0,Bt1,fac;
  REAL chargeA,chargeB,charge;
  ATOM *atom_pointer;
  REAL *kfactor,recip_cutoff,ksqr;
  PAIR pair;
  VECTOR *kvecs;
  COMPLEX temp_sum_bonddipole,temp_sum;
  int considered_charged;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  // put charge, bond-dipoles, and positions into appropriate arrays
  // ===============================================================

  nr_of_coulombic_sites=0;
  nr_frameworks=Framework[CurrentSystem].NumberOfFrameworks;
  for(f1=0;f1<nr_frameworks;f1++)
  {
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    for(i=0;i<nr_atoms;i++)
    {
      type=Framework[CurrentSystem].Atoms[f1][i].Type;
      charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged)
      {
        Charge[nr_of_coulombic_sites]=charge;
        if(!(Framework[CurrentSystem].Atoms[f1][i].Fixed.x&&Framework[CurrentSystem].Atoms[f1][i].Fixed.y&&Framework[CurrentSystem].Atoms[f1][i].Fixed.z)||BackPolarization)
        {
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          AtomVector[nr_of_coulombic_sites]=&Framework[CurrentSystem].Atoms[f1][i].ElectricField;
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_coulombic_sites_framework=nr_of_coulombic_sites;
  nr_molecules=NumberOfAdsorbateMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Adsorbates[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        charge=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        if(considered_charged)
        {
          Charge[nr_of_coulombic_sites]=charge;
          if(!(Adsorbates[CurrentSystem][i].Atoms[j].Fixed.x&&Adsorbates[CurrentSystem][i].Atoms[j].Fixed.y&&Adsorbates[CurrentSystem][i].Atoms[j].Fixed.z)||BackPolarization)
          {
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][i].Atoms[j].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            AtomVector[nr_of_coulombic_sites]=&Adsorbates[CurrentSystem][i].Atoms[j].ElectricField;
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }
  nr_of_coulombic_sites_adsorbate=nr_of_coulombic_sites;
  nr_molecules=NumberOfCationMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Cations[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Cations[CurrentSystem][i].Atoms[j].Type;
        charge=Cations[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        if(considered_charged)
        {
          Charge[nr_of_coulombic_sites]=charge;
          if(!(Cations[CurrentSystem][i].Atoms[j].Fixed.x&&Cations[CurrentSystem][i].Atoms[j].Fixed.y&&Cations[CurrentSystem][i].Atoms[j].Fixed.z)||BackPolarization)
          {
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][i].Atoms[j].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            AtomVector[nr_of_coulombic_sites]=&Cations[CurrentSystem][i].Atoms[j].ElectricField;
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }
  nr_of_coulombic_sites_cation=nr_of_coulombic_sites;

  nr_of_bonddipole_sites=0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      nr_of_bonddipoles=Framework[CurrentSystem].NumberOfBondDipoles[f1];
      for(i=0;i<nr_of_bonddipoles;i++)
      {
        pair=Framework[CurrentSystem].BondDipoles[f1][i];
        atom_pointer=Framework[CurrentSystem].Atoms[f1];
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        dipole=ApplyBoundaryCondition(dipole);
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Framework[CurrentSystem].BondDipoleMagnitude[f1][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].ElectricField;
        BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].ElectricField;
        nr_of_bonddipole_sites++;
      }
    }
  }
  nr_of_bonddipole_sites_framework=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Adsorbates[CurrentSystem][i].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].ElectricField;
      BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].ElectricField;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_adsorbate=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Cations[CurrentSystem][i].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].ElectricField;
      BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].ElectricField;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_site_cation=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          // loop over all the framework atoms
          sum_framework.re=0.0;
          sum_framework.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_framework.re+=temp*Eikr[i].re;
            sum_framework.im+=temp*Eikr[i].im;
          }

          // loop over all the adsorbate atoms
          sum_adsorbate.re=0.0;
          sum_adsorbate.im=0.0;
          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_adsorbate.re+=temp*Eikr[i].re;
            sum_adsorbate.im+=temp*Eikr[i].im;
          }

          // loop over all the cation atoms
          sum_cation.re=0.0;
          sum_cation.im=0.0;
          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_cation.re+=temp*Eikr[i].re;
            sum_cation.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_framework.re=0.0;
          sum_bonddipole_framework.im=0.0;

          rk=kvecs[nvec];

          for(i=0;i<nr_of_bonddipole_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_framework.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_framework.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_adsorbate.re=0.0;
          sum_bonddipole_adsorbate.im=0.0;
          for(;i<nr_of_bonddipole_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_adsorbate.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_adsorbate.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_cation.re=0.0;
          sum_bonddipole_cation.im=0.0;
          for(;i<nr_of_bonddipole_site_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_cation.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_cation.im+=temp*Eikr_bd[i].im;
          }

          // add the pre-computed contributions of fixed atoms
          // if 'BackPolarization' then all the atoms are explicility taken into account
          if(!BackPolarization)
          {
            sum_framework.re+=StoreRigidChargeFramework[CurrentSystem][nvec].re;
            sum_framework.im+=StoreRigidChargeFramework[CurrentSystem][nvec].im;
            sum_adsorbate.re+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].re;
            sum_adsorbate.im+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].im;
            sum_cation.re+=StoreRigidChargeCations[CurrentSystem][nvec].re;
            sum_cation.im+=StoreRigidChargeCations[CurrentSystem][nvec].im;

            sum_bonddipole_framework.re+=StoreRigidBondDipolesFramework[CurrentSystem][nvec].re;
            sum_bonddipole_framework.im+=StoreRigidBondDipolesFramework[CurrentSystem][nvec].im;
            sum_bonddipole_adsorbate.re+=StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec].re;
            sum_bonddipole_adsorbate.im+=StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec].im;
            sum_bonddipole_cation.re+=StoreRigidBondDipolesCations[CurrentSystem][nvec].re;
            sum_bonddipole_cation.im+=StoreRigidBondDipolesCations[CurrentSystem][nvec].im;
          }


          // get total sums
          sum.re=sum_framework.re+sum_adsorbate.re+sum_cation.re;
          sum.im=sum_framework.im+sum_adsorbate.im+sum_cation.im;
          sum_bonddipole.re=sum_bonddipole_framework.re+sum_bonddipole_adsorbate.re+sum_bonddipole_cation.re;
          sum_bonddipole.im=sum_bonddipole_framework.im+sum_bonddipole_adsorbate.im+sum_bonddipole_cation.im;

          // precomputed wavevector dependent pre-factor
          temp=kfactor[nvec];


          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitIntraFrameworkPolarization)
          {
            temp_sum.re-=sum_framework.re;
            temp_sum.im-=sum_framework.im;
            temp_sum_bonddipole.re-=sum_bonddipole_framework.re;
            temp_sum_bonddipole.im-=sum_bonddipole_framework.im;
          }

          for(i=0;i<nr_of_coulombic_sites_framework;i++)
          {
            fac=2.0*temp*((Eikr[i].im*temp_sum.re-Eikr[i].re*temp_sum.im)-
                          (Eikr[i].re*temp_sum_bonddipole.re+Eikr[i].im*temp_sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitAdsorbateAdsorbatePolarization)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }
          if(OmitAdsorbateCationPolarization)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }

          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            fac=2.0*temp*((Eikr[i].im*temp_sum.re-Eikr[i].re*temp_sum.im)-
                          (Eikr[i].re*temp_sum_bonddipole.re+Eikr[i].im*temp_sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitCationCationPolarization)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }
          if(OmitAdsorbateCationPolarization)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }

          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            fac=2.0*temp*((Eikr[i].im*temp_sum.re-Eikr[i].re*temp_sum.im)-
                          (Eikr[i].re*temp_sum_bonddipole.re+Eikr[i].im*temp_sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          // next wave-vector
          nvec++;
        }
      }
    }
  }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      atom_pointer=Framework[CurrentSystem].Atoms[f1];
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        pair=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        chargeA=atom_pointer[pair.A].Charge;
        chargeB=atom_pointer[pair.B].Charge;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        if(r>1e-4)
        {
          temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
          Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
          Bt1=temp+Bt0/rr;
        }
        else
        {
          Bt1=-COULOMBIC_CONVERSION_FACTOR*4.0*CUBE(alpha)/(3.0*sqrt(M_PI));
        }

        temp=chargeB*Bt1;
        atom_pointer[pair.A].ElectricField.x+=temp*dr.x;
        atom_pointer[pair.A].ElectricField.y+=temp*dr.y;
        atom_pointer[pair.A].ElectricField.z+=temp*dr.z;

        temp=chargeA*Bt1;
        atom_pointer[pair.B].ElectricField.x-=temp*dr.x;
        atom_pointer[pair.B].ElectricField.y-=temp*dr.y;
        atom_pointer[pair.B].ElectricField.z-=temp*dr.z;
      }
    }
  }

  // subtract adsorbate intra-charge-charge
  if(!OmitAdsorbateAdsorbatePolarization)
  {
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      atom_pointer=Adsorbates[CurrentSystem][m].Atoms;
      type=Adsorbates[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        pair=Components[type].ExcludedIntraChargeCharge[i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        chargeA=atom_pointer[pair.A].Charge;
        chargeB=atom_pointer[pair.B].Charge;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;

        temp=chargeB*Bt1;
        atom_pointer[pair.A].ElectricField.x+=temp*dr.x;
        atom_pointer[pair.A].ElectricField.y+=temp*dr.y;
        atom_pointer[pair.A].ElectricField.z+=temp*dr.z;

        temp=chargeA*Bt1;
        atom_pointer[pair.B].ElectricField.x-=temp*dr.x;
        atom_pointer[pair.B].ElectricField.y-=temp*dr.y;
        atom_pointer[pair.B].ElectricField.z-=temp*dr.z;
      }
    }
  }

  // subtract cation intra-charge-charge
  if(!OmitCationCationPolarization)
  {
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      atom_pointer=Cations[CurrentSystem][m].Atoms;
      type=Cations[CurrentSystem][m].Type;
      nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        pair=Components[type].ExcludedIntraChargeCharge[i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        chargeA=atom_pointer[pair.A].Charge;
        chargeB=atom_pointer[pair.B].Charge;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;

        temp=chargeB*Bt1;
        atom_pointer[pair.A].ElectricField.x+=temp*dr.x;
        atom_pointer[pair.A].ElectricField.y+=temp*dr.y;
        atom_pointer[pair.A].ElectricField.z+=temp*dr.z;

        temp=chargeA*Bt1;
        atom_pointer[pair.B].ElectricField.x-=temp*dr.x;
        atom_pointer[pair.B].ElectricField.y-=temp*dr.y;
        atom_pointer[pair.B].ElectricField.z-=temp*dr.z;
      }
    }
  }

  return 0;
}


/*********************************************************************************************************
 * Name       | ComputeStaticElectricFieldEwald                                                          *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of the static electric field for MC-moves.                     *
 *            | 1) computes the total one.                                                               *
 *            | 2) computed the total one when a particle is added.                                      *
 *            | 3) computed the total one when a particle is deleted.                                    *
 * Parameters |                                                                                          *
 * Limitation | Only charge-charge interaction is taken into account (no bond-dipoles yet).              *
 * Todo       | Add bond-dipoles to the static electric field.                                           *
 * Problem    | A bond-dipole exerts forces on both the atoms -> E!=F/mu                                 *
 *********************************************************************************************************/

int ComputeStaticElectricFieldEwald(int NewMolecule,int excl_ads,int excl_cation)
{
  int i,j,m,f1,ii,jj,kk;
  int nvec,type_mol,type,typeA,typeB;
  REAL temp,alpha,r,rr;
  VECTOR pos,rk,dr,posA,posB;
  int nr_molecules,nr_atoms,nr_frameworks,nr_of_excluded_pairs;
  int index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_framework,nr_of_coulombic_sites_adsorbate,nr_of_coulombic_sites_cation;
  COMPLEX sum,sum_framework,sum_adsorbate,sum_cation;
  int nr_of_bonddipoles,nr_of_bonddipole_sites,nr_of_bonddipole_sites_framework,nr_of_bonddipole_sites_adsorbate,nr_of_bonddipole_site_cation;
  COMPLEX sum_bonddipole,sum_bonddipole_framework,sum_bonddipole_adsorbate,sum_bonddipole_cation;
  VECTOR dipole;
  REAL Bt0,Bt1,fac;
  ATOM *atom_pointer;
  REAL *kfactor,recip_cutoff,ksqr;
  PAIR pair;
  VECTOR *kvecs;
  REAL chargeA, chargeB,charge;
  int NumberOfExcludedPairs;
  COMPLEX temp_sum_bonddipole,temp_sum;
  int considered_charged;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  // put charge, bond-dipoles, and positions into appropriate arrays
  // ===============================================================

  nr_of_coulombic_sites=0;
  nr_frameworks=Framework[CurrentSystem].NumberOfFrameworks;
  for(f1=0;f1<nr_frameworks;f1++)
  {
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    for(i=0;i<nr_atoms;i++)
    {
      type=Framework[CurrentSystem].Atoms[f1][i].Type;
      charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged)
      {
        Charge[nr_of_coulombic_sites]=charge;
        if(!(Framework[CurrentSystem].Atoms[f1][i].Fixed.x&&Framework[CurrentSystem].Atoms[f1][i].Fixed.y&&Framework[CurrentSystem].Atoms[f1][i].Fixed.z)||BackPolarization)
        {
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          AtomVector[nr_of_coulombic_sites]=&Framework[CurrentSystem].Atoms[f1][i].ElectricField;
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_coulombic_sites_framework=nr_of_coulombic_sites;
  nr_molecules=NumberOfAdsorbateMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    if(i!=excl_ads)
    {
      type_mol=Adsorbates[CurrentSystem][i].Type;
      if(Components[type_mol].HasCharges)
      {
        nr_atoms=Adsorbates[CurrentSystem][i].NumberOfAtoms;
        for(j=0;j<nr_atoms;j++)
        {
          type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          charge=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
          if(considered_charged)
          {
            Charge[nr_of_coulombic_sites]=charge;
            if(!(Adsorbates[CurrentSystem][i].Atoms[j].Fixed.x&&Adsorbates[CurrentSystem][i].Atoms[j].Fixed.y&&Adsorbates[CurrentSystem][i].Atoms[j].Fixed.z)||BackPolarization)
            {
              Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][i].Atoms[j].Position);
              Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
              AtomVector[nr_of_coulombic_sites]=&Adsorbates[CurrentSystem][i].Atoms[j].ElectricField;
              nr_of_coulombic_sites++;
            }
          }
        }
      }
    }
  }
  if(NewMolecule&&(!Components[CurrentComponent].ExtraFrameworkMolecule))
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      type=Components[CurrentComponent].Type[i];
      if(PseudoAtoms[type].HasCharges)
      {
        Charge[nr_of_coulombic_sites]=Components[CurrentComponent].Charge[i];
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(TrialPosition[CurrentSystem][i]);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        AtomVector[nr_of_coulombic_sites]=&ElectricFieldAtTrialPosition[CurrentSystem][i];
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_adsorbate=nr_of_coulombic_sites;

  nr_molecules=NumberOfCationMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    if(i!=excl_cation)
    {
      type_mol=Cations[CurrentSystem][i].Type;
      if(Components[type_mol].HasCharges)
      {
        nr_atoms=Cations[CurrentSystem][i].NumberOfAtoms;
        for(j=0;j<nr_atoms;j++)
        {
          type=Cations[CurrentSystem][i].Atoms[j].Type;
          charge=Cations[CurrentSystem][i].Atoms[j].Charge;
          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
          if(considered_charged)
          {
            Charge[nr_of_coulombic_sites]=charge;
            if(!(Cations[CurrentSystem][i].Atoms[j].Fixed.x&&Cations[CurrentSystem][i].Atoms[j].Fixed.y&&Cations[CurrentSystem][i].Atoms[j].Fixed.z)||BackPolarization)
            {
              Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][i].Atoms[j].Position);
              Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
              AtomVector[nr_of_coulombic_sites]=&Cations[CurrentSystem][i].Atoms[j].ElectricField;
              nr_of_coulombic_sites++;
            }
          }
        }
      }
    }
  }
  if(NewMolecule&&Components[CurrentComponent].ExtraFrameworkMolecule)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      type=Components[CurrentComponent].Type[i];
      if(PseudoAtoms[type].HasCharges)
      {
        Charge[nr_of_coulombic_sites]=Components[CurrentComponent].Charge[i];
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(TrialPosition[CurrentSystem][i]);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        AtomVector[nr_of_coulombic_sites]=&ElectricFieldAtTrialPosition[CurrentSystem][i];
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_cation=nr_of_coulombic_sites;

  nr_of_bonddipole_sites=0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      nr_of_bonddipoles=Framework[CurrentSystem].NumberOfBondDipoles[f1];
      for(i=0;i<nr_of_bonddipoles;i++)
      {
        pair=Framework[CurrentSystem].BondDipoles[f1][i];
        atom_pointer=Framework[CurrentSystem].Atoms[f1];
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        dipole=ApplyBoundaryCondition(dipole);
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Framework[CurrentSystem].BondDipoleMagnitude[f1][i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI; 
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].ElectricField;
        BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].ElectricField;
        nr_of_bonddipole_sites++;
      }
    }
  }
  nr_of_bonddipole_sites_framework=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    if(i!=excl_ads)
    {
      type_mol=Adsorbates[CurrentSystem][i].Type;
      nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
      for(j=0;j<nr_of_bonddipoles;j++)
      {
        pair=Components[type_mol].BondDipoles[j];
        atom_pointer=Adsorbates[CurrentSystem][i].Atoms;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].ElectricField;
        BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].ElectricField;
        nr_of_bonddipole_sites++;
      }
    }
  }
  if(NewMolecule&&(!Components[CurrentComponent].ExtraFrameworkMolecule))
  {
    for(i=0;i<Components[CurrentComponent].NumberOfBondDipoles;i++)
    {
      pair=Components[CurrentComponent].BondDipoles[i];
      posA=TrialPosition[CurrentSystem][pair.A];
      posB=TrialPosition[CurrentSystem][pair.B];
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[CurrentComponent].BondDipoleMagnitude[i];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_adsorbate=nr_of_bonddipole_sites;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    if(i!=excl_cation)
    {
      type_mol=Cations[CurrentSystem][i].Type;
      nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
      for(j=0;j<nr_of_bonddipoles;j++)
      {
        pair=Components[type_mol].BondDipoles[j];
        atom_pointer=Cations[CurrentSystem][i].Atoms;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;
        dipole.x=posB.x-posA.x;
        dipole.y=posB.y-posA.y;
        dipole.z=posB.z-posA.z;
        pos.x=posA.x+0.5*dipole.x;
        pos.y=posA.y+0.5*dipole.y;
        pos.z=posA.z+0.5*dipole.z;
        temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
        DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
        DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
        DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
        BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
        BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
        BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
        BondDipoleForcesA[nr_of_bonddipole_sites]=&atom_pointer[pair.A].ElectricField;
        BondDipoleForcesB[nr_of_bonddipole_sites]=&atom_pointer[pair.B].ElectricField;
        nr_of_bonddipole_sites++;
      }
    }
  }
  if(NewMolecule&&Components[CurrentComponent].ExtraFrameworkMolecule)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfBondDipoles;i++)
    {
      pair=Components[CurrentComponent].BondDipoles[i];
      posA=TrialPosition[CurrentSystem][pair.A];
      posB=TrialPosition[CurrentSystem][pair.B];
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[CurrentComponent].BondDipoleMagnitude[i];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_site_cation=nr_of_bonddipole_sites;


  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          // loop over all the framework atoms
          sum_framework.re=0.0;
          sum_framework.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_framework.re+=temp*Eikr[i].re;
            sum_framework.im+=temp*Eikr[i].im;
          }

          // loop over all the adsorbate atoms
          sum_adsorbate.re=0.0;
          sum_adsorbate.im=0.0;
          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_adsorbate.re+=temp*Eikr[i].re;
            sum_adsorbate.im+=temp*Eikr[i].im;
          }

          // loop over all the cation atoms
          sum_cation.re=0.0;
          sum_cation.im=0.0;
          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_cation.re+=temp*Eikr[i].re;
            sum_cation.im+=temp*Eikr[i].im;
          }

          rk=kvecs[nvec];

          sum_bonddipole_framework.re=0.0;
          sum_bonddipole_framework.im=0.0;
          for(i=0;i<nr_of_bonddipole_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_framework.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_framework.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_adsorbate.re=0.0;
          sum_bonddipole_adsorbate.im=0.0;
          for(;i<nr_of_bonddipole_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_adsorbate.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_adsorbate.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_cation.re=0.0;
          sum_bonddipole_cation.im=0.0;
          for(;i<nr_of_bonddipole_site_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_cation.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_cation.im+=temp*Eikr_bd[i].im;
          }

          // add the pre-computed contributions of fixed atoms
          // if 'BackPolarization' then all the atoms are explicility taken into account
          if(!BackPolarization)
          {
            sum_framework.re+=StoreRigidChargeFramework[CurrentSystem][nvec].re;
            sum_framework.im+=StoreRigidChargeFramework[CurrentSystem][nvec].im;
            sum_adsorbate.re+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].re;
            sum_adsorbate.im+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].im;
            sum_cation.re+=StoreRigidChargeCations[CurrentSystem][nvec].re;
            sum_cation.im+=StoreRigidChargeCations[CurrentSystem][nvec].im;

            sum_bonddipole_framework.re+=StoreRigidBondDipolesFramework[CurrentSystem][nvec].re;
            sum_bonddipole_framework.im+=StoreRigidBondDipolesFramework[CurrentSystem][nvec].im;
            sum_bonddipole_adsorbate.re+=StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec].re;
            sum_bonddipole_adsorbate.im+=StoreRigidBondDipolesAdsorbates[CurrentSystem][nvec].im;
            sum_bonddipole_cation.re+=StoreRigidBondDipolesCations[CurrentSystem][nvec].re;
            sum_bonddipole_cation.im+=StoreRigidBondDipolesCations[CurrentSystem][nvec].im;
          }

          // get total sums
          sum.re=sum_framework.re+sum_adsorbate.re+sum_cation.re;
          sum.im=sum_framework.im+sum_adsorbate.im+sum_cation.im;
          sum_bonddipole.re=sum_bonddipole_framework.re+sum_bonddipole_adsorbate.re+sum_bonddipole_cation.re;
          sum_bonddipole.im=sum_bonddipole_framework.im+sum_bonddipole_adsorbate.im+sum_bonddipole_cation.im;

          // precomputed wavevector dependent pre-factor
          temp=kfactor[nvec];

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitIntraFrameworkPolarization)
          {
            temp_sum.re-=sum_framework.re;
            temp_sum.im-=sum_framework.im;
            temp_sum_bonddipole.re-=sum_bonddipole_framework.re;
            temp_sum_bonddipole.im-=sum_bonddipole_framework.im;
          }

          // forces on atoms from other charges and bond-dipoles
          for(i=0;i<nr_of_coulombic_sites_framework;i++)
          {
            fac=2.0*temp*((Eikr[i].im*temp_sum.re-Eikr[i].re*temp_sum.im)-
                          (Eikr[i].re*temp_sum_bonddipole.re+Eikr[i].im*temp_sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitAdsorbateAdsorbatePolarization)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }
          if(OmitAdsorbateCationPolarization)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }

          for(;i<nr_of_coulombic_sites_adsorbate;i++)
          {
            fac=2.0*temp*((Eikr[i].im*temp_sum.re-Eikr[i].re*temp_sum.im)-
                          (Eikr[i].re*temp_sum_bonddipole.re+Eikr[i].im*temp_sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          temp_sum.re=sum.re;
          temp_sum.im=sum.im;
          temp_sum_bonddipole.re=sum_bonddipole.re;
          temp_sum_bonddipole.im=sum_bonddipole.im;

          if(OmitCationCationPolarization)
          {
            temp_sum.re-=sum_cation.re;
            temp_sum.im-=sum_cation.im;
            temp_sum_bonddipole.re-=sum_bonddipole_cation.re;
            temp_sum_bonddipole.im-=sum_bonddipole_cation.im;
          }
          if(OmitAdsorbateCationPolarization)
          {
            temp_sum.re-=sum_adsorbate.re;
            temp_sum.im-=sum_adsorbate.im;
            temp_sum_bonddipole.re-=sum_bonddipole_adsorbate.re;
            temp_sum_bonddipole.im-=sum_bonddipole_adsorbate.im;
          }

          for(;i<nr_of_coulombic_sites_cation;i++)
          {
            fac=2.0*temp*((Eikr[i].im*temp_sum.re-Eikr[i].re*temp_sum.im)-
                          (Eikr[i].re*temp_sum_bonddipole.re+Eikr[i].im*temp_sum_bonddipole.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          // next wave-vector
          nvec++;
        }
      }
    }
  }

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      atom_pointer=Framework[CurrentSystem].Atoms[f1];
      nr_of_excluded_pairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[f1];
      for(i=0;i<nr_of_excluded_pairs;i++)
      {
        pair=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i];
        typeA=atom_pointer[pair.A].Type;
        typeB=atom_pointer[pair.B].Type;
        chargeA=atom_pointer[pair.A].Charge;
        chargeB=atom_pointer[pair.B].Charge;
        posA=atom_pointer[pair.A].Position;
        posB=atom_pointer[pair.B].Position;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        if(r>1e-4)
        {
          temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
          Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
          Bt1=temp+Bt0/rr;
        }
        else
        {
          Bt1=-COULOMBIC_CONVERSION_FACTOR*4.0*CUBE(alpha)/(3.0*sqrt(M_PI));
        }

        temp=chargeB*Bt1;
        atom_pointer[pair.A].ElectricField.x+=temp*dr.x;
        atom_pointer[pair.A].ElectricField.y+=temp*dr.y;
        atom_pointer[pair.A].ElectricField.z+=temp*dr.z;

        temp=chargeA*Bt1;
        atom_pointer[pair.B].ElectricField.x-=temp*dr.x;
        atom_pointer[pair.B].ElectricField.y-=temp*dr.y;
        atom_pointer[pair.B].ElectricField.z-=temp*dr.z;
      }
    }
  }

  // subtract adsorbate intra-charge-charge
  if(!OmitAdsorbateAdsorbatePolarization)
  {
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      if(m!=excl_ads)
      {
        atom_pointer=Adsorbates[CurrentSystem][m].Atoms;
        type=Adsorbates[CurrentSystem][m].Type;
        nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
        for(i=0;i<nr_of_excluded_pairs;i++)
        {
          pair=Components[type].ExcludedIntraChargeCharge[i];
          typeA=atom_pointer[pair.A].Type;
          typeB=atom_pointer[pair.B].Type;
          chargeA=atom_pointer[pair.A].Charge;
          chargeB=atom_pointer[pair.B].Charge;
          posA=atom_pointer[pair.A].Position;
          posB=atom_pointer[pair.B].Position;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          r=sqrt(rr);

          temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
          Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
          Bt1=temp+Bt0/rr;

          temp=chargeB*Bt1;
          atom_pointer[pair.A].ElectricField.x+=temp*dr.x;
          atom_pointer[pair.A].ElectricField.y+=temp*dr.y;
          atom_pointer[pair.A].ElectricField.z+=temp*dr.z;

          temp=chargeA*Bt1;
          atom_pointer[pair.B].ElectricField.x-=temp*dr.x;
          atom_pointer[pair.B].ElectricField.y-=temp*dr.y;
          atom_pointer[pair.B].ElectricField.z-=temp*dr.z;
        }
      }
    }
  }

  // subtract cation intra-charge-charge
  if(!OmitCationCationPolarization)
  {
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      if(m!=excl_cation)
      {
        atom_pointer=Cations[CurrentSystem][m].Atoms;
        type=Cations[CurrentSystem][m].Type;
        nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
        for(i=0;i<nr_of_excluded_pairs;i++)
        {
          pair=Components[type].ExcludedIntraChargeCharge[i];
          typeA=atom_pointer[pair.A].Type;
          typeB=atom_pointer[pair.B].Type;
          chargeA=atom_pointer[pair.A].Charge;
          chargeB=atom_pointer[pair.B].Charge;
          posA=atom_pointer[pair.A].Position;
          posB=atom_pointer[pair.B].Position;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          r=sqrt(rr);

          temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
          Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
          Bt1=temp+Bt0/rr;

          temp=chargeB*Bt1;
          atom_pointer[pair.A].ElectricField.x+=temp*dr.x;
          atom_pointer[pair.A].ElectricField.y+=temp*dr.y;
          atom_pointer[pair.A].ElectricField.z+=temp*dr.z;

          temp=chargeA*Bt1;
          atom_pointer[pair.B].ElectricField.x-=temp*dr.x;
          atom_pointer[pair.B].ElectricField.y-=temp*dr.y;
          atom_pointer[pair.B].ElectricField.z-=temp*dr.z;
        }
      }
    }
  }

  if(NewMolecule)
  {
    if((Components[CurrentComponent].ExtraFrameworkMolecule&&(!OmitCationCationPolarization))||
      ((!Components[CurrentComponent].ExtraFrameworkMolecule)&&(!OmitAdsorbateAdsorbatePolarization)))
    {
      NumberOfExcludedPairs=Components[CurrentComponent].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        pair=Components[CurrentComponent].ExcludedIntraChargeCharge[i];

        typeA=Components[CurrentSystem].Type[pair.A];
        typeB=Components[CurrentSystem].Type[pair.B];
        chargeA=Components[CurrentSystem].Charge[pair.A];
        chargeB=Components[CurrentSystem].Charge[pair.B];
        posA=TrialPosition[CurrentSystem][pair.A];
        posB=TrialPosition[CurrentSystem][pair.B];

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        temp=COULOMBIC_CONVERSION_FACTOR*(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
        Bt0=-COULOMBIC_CONVERSION_FACTOR*erf(alpha*r)/r;
        Bt1=temp+Bt0/rr;

        temp=chargeB*Bt1;
        ElectricFieldAtTrialPosition[CurrentSystem][pair.A].x+=temp*dr.x;
        ElectricFieldAtTrialPosition[CurrentSystem][pair.A].y+=temp*dr.y;
        ElectricFieldAtTrialPosition[CurrentSystem][pair.A].z+=temp*dr.z;

        temp=chargeA*Bt1;
        ElectricFieldAtTrialPosition[CurrentSystem][pair.B].x-=temp*dr.x;
        ElectricFieldAtTrialPosition[CurrentSystem][pair.B].y-=temp*dr.y;
        ElectricFieldAtTrialPosition[CurrentSystem][pair.B].z-=temp*dr.z;
      }
    }
  }

  return 0;
}

// TODO: write function that computes the diference in electric-field gradient for a molecule (much faster)
/*********************************************************************************************************
 * Name       | ComputeElectricFieldFromInducedDipolesEwald                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of the difference in polarization energy.                      *
 * Parameters |                                                                                          *
 * Note       | If no backpolarization is used then the electric field and difference in polarization    *
 *            | energy can be computed from the difference of the adsorbate atoms.                       *
 *********************************************************************************************************/
/*
int ComputeStaticElectricFieldEwaldAdsorbateDifference(int NewMolecule,int OldMolecule,int mol,int store)
{
  int i,j,ii,jj,kk;
  int A,B,nvec,nr_of_excluded_pairs;
  int kmax_x,kmax_y,kmax_z,index_i,index_j,index_k;
  int type_mol,nr_atoms,type,nr_of_bonddipoles;
  int nr_of_coulombic_sites,nr_of_coulombic_sites_old,nr_of_coulombic_sites_new;
  int nr_of_bonddipole_sites,nr_of_bonddipole_sites_old,nr_of_bonddipole_sites_new;
  COMPLEX sum_old,sum_new;
  COMPLEX sum_bonddipole_old,sum_bonddipole_new;
  REAL fac;
  REAL alpha,chargeA,chargeB,charge,r,rr;
  REAL cosA,cosB,cosAB,Bt0,Bt1,Bt2,temp;
  VECTOR pos,posA,posB,dr;
  VECTOR dipole,dipoleA,dipoleB,rk;
  VECTOR posA1,posA2,posB1,posB2;
  VECTOR *kvecs;
  REAL *kfactor,recip_cutoff,ksqr;
  PAIR pair;
  ATOM *atom_pointer;
  REAL polarization_energy_old,polarization_energy_new;
  COMPLEX sum_adsorbates_new,sum_adsorbates_old;
  int considered_charged;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return 0;
  if(OmitEwaldFourier) return 0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  nvec=0;
  polarization_energy_old=polarization_energy_new=0.0;
  fac=0.0;

  nr_of_coulombic_sites=nr_of_coulombic_sites_old=nr_of_coulombic_sites_new=0;
  nr_of_bonddipole_sites=nr_of_bonddipole_sites_old=nr_of_bonddipole_sites_new=0;

  if(OldMolecule)
  {
    nr_atoms=Adsorbates[CurrentSystem][mol].NumberOfAtoms;
    for(j=0;j<nr_atoms;j++)
    {
      type=Adsorbates[CurrentSystem][mol].Atoms[j].Type;
      charge=Adsorbates[CurrentSystem][mol].Atoms[j].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged)
      {
        Charge[nr_of_coulombic_sites]=charge;
        Polarization[nr_of_coulombic_sites]=PseudoAtoms[type].Polarization;
        AtomVector[nr_of_coulombic_sites]=&Adsorbates[CurrentSystem][mol].Atoms[j].ElectricField;
        AtomVector[nr_of_coulombic_sites]->x=0.0;
        AtomVector[nr_of_coulombic_sites]->y=0.0;
        AtomVector[nr_of_coulombic_sites]->z=0.0;
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][mol].Atoms[j].Position);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_old=nr_of_coulombic_sites;

  if(NewMolecule)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      type=Components[CurrentComponent].Type[i];
      if(PseudoAtoms[type].HasCharges)
      {
        Charge[nr_of_coulombic_sites]=Components[CurrentComponent].Type[i];
        Polarization[nr_of_coulombic_sites]=PseudoAtoms[type].Polarization;
        AtomVector[nr_of_coulombic_sites]=&ElectricFieldAtTrialPosition[CurrentSystem][i];
        AtomVector[nr_of_coulombic_sites]->x=0.0;
        AtomVector[nr_of_coulombic_sites]->y=0.0;
        AtomVector[nr_of_coulombic_sites]->z=0.0;
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(TrialPosition[CurrentSystem][i]);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_coulombic_sites_new=nr_of_coulombic_sites;

  if(OldMolecule)
  {
    type_mol=Adsorbates[CurrentSystem][mol].Type;
    nr_of_bonddipoles=Components[type_mol].NumberOfBondDipoles;
    for(j=0;j<nr_of_bonddipoles;j++)
    {
      pair=Components[type_mol].BondDipoles[j];
      atom_pointer=Adsorbates[CurrentSystem][mol].Atoms;
      posA=atom_pointer[pair.A].Position;
      posB=atom_pointer[pair.B].Position;
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[j]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[j];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_old=nr_of_bonddipole_sites;

  if(NewMolecule)
  {
    type_mol=CurrentComponent;
    for(i=0;i<Components[type_mol].NumberOfBondDipoles;i++)
    {
      pair=Components[type_mol].BondDipoles[i];
      posA=TrialPosition[CurrentSystem][pair.A];
      posB=TrialPosition[CurrentSystem][pair.B];
      dipole.x=posB.x-posA.x;
      dipole.y=posB.y-posA.y;
      dipole.z=posB.z-posA.z;
      pos.x=posA.x+0.5*dipole.x;
      pos.y=posA.y+0.5*dipole.y;
      pos.z=posA.z+0.5*dipole.z;
      temp=Components[type_mol].BondDipoleMagnitude[i]/sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      DipoleVector[nr_of_bonddipole_sites].x=temp*dipole.x;
      DipoleVector[nr_of_bonddipole_sites].y=temp*dipole.y;
      DipoleVector[nr_of_bonddipole_sites].z=temp*dipole.z;
      BondDipoleMagnitude[nr_of_bonddipole_sites]=Components[type_mol].BondDipoleMagnitude[i];
      BondLength[nr_of_bonddipole_sites]=sqrt(SQR(dipole.x)+SQR(dipole.y)+SQR(dipole.z));
      BondDipolePositions[nr_of_bonddipole_sites]=ConvertFromXYZtoABC(pos);
      BondDipolePositions[nr_of_bonddipole_sites].x*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].y*=TWO_PI;
      BondDipolePositions[nr_of_bonddipole_sites].z*=TWO_PI;
      nr_of_bonddipole_sites++;
    }
  }
  nr_of_bonddipole_sites_new=nr_of_bonddipole_sites;

  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }


  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_bonddipole_sites;i++)
  {
    Eikx_bd[i].re=1.0; Eikx_bd[i].im=0.0;
    Eiky_bd[i].re=1.0; Eiky_bd[i].im=0.0;
    Eikz_bd[i].re=1.0; Eikz_bd[i].im=0.0;

    pos=BondDipolePositions[i];

    index_i=MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x);  Eikx_bd[index_i].im=sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y);  Eiky_bd[index_i].im=sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z);  Eikz_bd[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfBondDipoleSites+i;
    Eikx_bd[index_i].re=cos(pos.x); Eikx_bd[index_i].im=-sin(pos.x);
    Eiky_bd[index_i].re=cos(pos.y); Eiky_bd[index_i].im=-sin(pos.y);
    Eikz_bd[index_i].re=cos(pos.z); Eikz_bd[index_i].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikx_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikx_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikx_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikx_bd[MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im=Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eiky_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eiky_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eiky_bd[MaxNumberOfBondDipoleSites+i].im;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eiky_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eiky_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eiky_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_bonddipole_sites;i++)
    {
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].re-
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im=Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].im*Eikz_bd[MaxNumberOfBondDipoleSites+i].re+
                                                 Eikz_bd[(j-1)*MaxNumberOfBondDipoleSites+i].re*Eikz_bd[MaxNumberOfBondDipoleSites+i].im;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].re=Eikz_bd[j*MaxNumberOfBondDipoleSites+i].re;
      Eikz_bd[-j*MaxNumberOfBondDipoleSites+i].im=-Eikz_bd[j*MaxNumberOfBondDipoleSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }
      for(i=0;i<nr_of_bonddipole_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfBondDipoleSites+i;
        index_j=jj*MaxNumberOfBondDipoleSites+i;
        Eikr_xy_bd[i].re=Eikx_bd[index_i].re*Eiky_bd[index_j].re-Eikx_bd[index_i].im*Eiky_bd[index_j].im;
        Eikr_xy_bd[i].im=Eikx_bd[index_i].im*Eiky_bd[index_j].re+Eikx_bd[index_i].re*Eiky_bd[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          temp=kfactor[nvec];
          rk=kvecs[nvec];

          sum_adsorbates_old.re=0.0;
          sum_adsorbates_old.im=0.0;
          for(i=0;i<nr_of_coulombic_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_adsorbates_old.re+=temp*Eikr[i].re;
            sum_adsorbates_old.im+=temp*Eikr[i].im;
          }
          sum_adsorbates_new.re=0.0;
          sum_adsorbates_new.im=0.0;
          for(;i<nr_of_coulombic_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            // add contribution to the sum
            temp=Charge[i];
            sum_adsorbates_new.re+=temp*Eikr[i].re;
            sum_adsorbates_new.im+=temp*Eikr[i].im;
          }

          sum_bonddipole_old.re=0.0;
          sum_bonddipole_old.im=0.0;
          for(i=0;i<nr_of_bonddipole_sites_old;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_old.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_old.im+=temp*Eikr_bd[i].im;
          }

          sum_bonddipole_new.re=0.0;
          sum_bonddipole_new.im=0.0;
          for(;i<nr_of_bonddipole_sites_new;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfBondDipoleSites+i;
            Eikr_bd[i].re=Eikr_xy_bd[i].re*Eikz_bd[index_k].re-Eikr_xy_bd[i].im*Eikz_bd[index_k].im;
            Eikr_bd[i].im=Eikr_xy_bd[i].im*Eikz_bd[index_k].re+Eikr_xy_bd[i].re*Eikz_bd[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_bonddipole_new.re+=temp*Eikr_bd[i].re;
            sum_bonddipole_new.im+=temp*Eikr_bd[i].im;
          }

          sum_new.re=StoreTotalChargeFramework[CurrentSystem][nvec].re+StoreTotalChargeAdsorbates[CurrentSystem][nvec].re+StoreTotalChargeCations[CurrentSystem][nvec].re+sum_adsorbates_new.re-sum_adsorbates_old.re;
          sum_new.im=StoreTotalChargeFramework[CurrentSystem][nvec].im+StoreTotalChargeAdsorbates[CurrentSystem][nvec].im+StoreTotalChargeCations[CurrentSystem][nvec].im+sum_adsorbates_new.im-sum_adsorbates_old.im;

          sum_old.re=StoreTotalChargeFramework[CurrentSystem][nvec].re+StoreTotalChargeAdsorbates[CurrentSystem][nvec].re+StoreTotalChargeCations[CurrentSystem][nvec].re;
          sum_old.im=StoreTotalChargeFramework[CurrentSystem][nvec].im+StoreTotalChargeAdsorbates[CurrentSystem][nvec].im+StoreTotalChargeCations[CurrentSystem][nvec].im;

          for(i=0;i<nr_of_coulombic_sites_old;i++)
          {
            fac=2.0*temp*((Eikr[i].im*sum_old.re-Eikr[i].re*sum_old.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }

          for(;i<nr_of_coulombic_sites_new;i++)
          {
            fac=2.0*temp*((Eikr[i].im*sum_new.re-Eikr[i].re*sum_new.im));
            AtomVector[i]->x+=fac*rk.x;
            AtomVector[i]->y+=fac*rk.y;
            AtomVector[i]->z+=fac*rk.z;
          }




          // store the new sums, these will be the current ones on acceptance of the mc-move
          //NewTotalChargeAdsorbates[store][nvec]=sum_adsorbates;
          //NewTotalBondDipolesAdsorbates[store][nvec]=sum_bonddipole_adsorbates;

          // next wave-vector
          nvec++;
        }
      }
    }
  }

  polarization_energy_old=0.0;
  for(i=0;i<nr_of_coulombic_sites_old;i++)
    polarization_energy_old-=0.5*Polarization[i]*(SQR(AtomVector[i]->x)+SQR(AtomVector[i]->y)+SQR(AtomVector[i]->z));

  polarization_energy_new=0.0;
  for(;i<nr_of_coulombic_sites_new;i++)
    polarization_energy_new-=0.5*Polarization[i]*(SQR(AtomVector[i]->x)+SQR(AtomVector[i]->y)+SQR(AtomVector[i]->z));

  printf("Difference in polarization energy: %18.10f %18.10f\n",polarization_energy_new,polarization_energy_old);

  if(NewMolecule)
  {
    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraChargeCharge[i].A;
      B=Components[CurrentComponent].ExcludedIntraChargeCharge[i].B;
      chargeA=Components[CurrentComponent].Charge[A];
      chargeB=Components[CurrentComponent].Charge[B];
      posA=TrialPosition[CurrentSystem][A];
      posB=TrialPosition[CurrentSystem][B];

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      Bt0=-erf(alpha*r)/r;
    }

    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraChargeBondDipole[i].A;
      B=Components[CurrentComponent].ExcludedIntraChargeBondDipole[i].B;

      chargeA=Components[CurrentComponent].Charge[A];
      posA=TrialPosition[CurrentSystem][A];

      pair=Components[CurrentComponent].BondDipoles[B];
      posB1=TrialPosition[CurrentSystem][pair.A];
      posB2=TrialPosition[CurrentSystem][pair.B];

      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(alpha*r)/r;
      Bt1=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr+Bt0/rr;

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
    }


    nr_of_excluded_pairs=Components[CurrentComponent].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[CurrentComponent].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[CurrentComponent].ExcludedIntraBondDipoleBondDipole[i].B;

      pair=Components[CurrentComponent].BondDipoles[A];
      posA1=TrialPosition[CurrentSystem][pair.A];
      posA2=TrialPosition[CurrentSystem][pair.B];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      pair=Components[CurrentComponent].BondDipoles[B];
      posB1=TrialPosition[CurrentSystem][pair.A];
      posB2=TrialPosition[CurrentSystem][pair.B];
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[CurrentComponent].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      temp=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
      Bt0=-erf(alpha*r)/r;
      Bt1=temp+Bt0/rr;
      temp*=2.0*SQR(alpha);
      Bt2=temp+(3.0/rr)*Bt1;

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
    }
  }

  if(OldMolecule)
  {
    type=Adsorbates[CurrentSystem][mol].Type;

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeCharge;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      pair=Components[type].ExcludedIntraChargeCharge[i];
      chargeA=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Charge;
      chargeB=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Charge;
      posA=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posB=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

      Bt0=-erf(alpha*r)/r;
    }

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraChargeBondDipole[i].A;
      B=Components[type].ExcludedIntraChargeBondDipole[i].B;

      chargeA=Adsorbates[CurrentSystem][mol].Atoms[A].Charge;
      posA=Adsorbates[CurrentSystem][mol].Atoms[A].Position;

      pair=Components[type].BondDipoles[B];
      posB1=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posB2=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      Bt0=-erf(alpha*r)/r;
      Bt1=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr+Bt0/rr;

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
    }

    nr_of_excluded_pairs=Components[type].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<nr_of_excluded_pairs;i++)
    {
      A=Components[type].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[type].ExcludedIntraBondDipoleBondDipole[i].B;

      pair=Components[type].BondDipoles[A];
      posA1=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posA2=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      temp=Components[type].BondDipoleMagnitude[A]/sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      pair=Components[type].BondDipoles[B];
      posB1=Adsorbates[CurrentSystem][mol].Atoms[pair.A].Position;
      posB2=Adsorbates[CurrentSystem][mol].Atoms[pair.B].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      temp=Components[type].BondDipoleMagnitude[B]/sqrt(SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z));
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      temp=(2.0/sqrt(M_PI))*alpha*exp(-SQR(-alpha*r))/rr;
      Bt0=-erf(alpha*r)/r;
      Bt1=temp+Bt0/rr;
      temp*=2.0*SQR(alpha);
      Bt2=temp+(3.0/rr)*Bt1;

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
    }
  }

  return 0;
}
*/

/*********************************************************************************************************
 * Name       | ComputeElectricFieldFromInducedDipolesEwald                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of the induced electric field.                                 *
 * Parameters |                                                                                          *
 * Options    | 'OmitIntraFrameworkPolarization', 'OmitAdsorbateAdsorbatePolarization',                  *
 *            | 'OmitCationCationPolarization', 'OmitAdsorbateCationPolarization'.                       *
 *********************************************************************************************************/

void ComputeElectricFieldFromInducedDipolesEwald(void)
{
  int i,j,f1,ii,jj,kk,index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z,nvec;
  int nr_of_coulombic_sites,nr_frameworks,nr_molecules,nr_atoms;
  int nr_of_induced_dipole_sites_framework,nr_of_induced_dipole_sites_adsorbate,nr_of_induced_dipole_sites_cation;
  int type,type_mol;
  COMPLEX sum_induced_dipole_framework,sum_induced_dipole_adsorbate,sum_induced_dipole_cation,sum_induced_dipole;
  COMPLEX temp_sum_induced_dipole;
  REAL alpha,fac,temp,fac_self,factor,charge;
  VECTOR *kvecs,pos,rk,dipole;
  REAL *kfactor,recip_cutoff,ksqr;
  int considered_charged;

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return;
  if(OmitEwaldFourier) return;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  nr_of_coulombic_sites=0;
  nr_frameworks=Framework[CurrentSystem].NumberOfFrameworks;
  for(f1=0;f1<nr_frameworks;f1++)
  {
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    for(i=0;i<nr_atoms;i++)
    {
      type=Framework[CurrentSystem].Atoms[f1][i].Type;
      charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged)
      {
        DipoleVector[nr_of_coulombic_sites]=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        AtomVector[nr_of_coulombic_sites]=&Framework[CurrentSystem].Atoms[f1][i].ElectricField;
        if(!OmitIntraFrameworkPolarization)
        {
          AtomVector[nr_of_coulombic_sites]->x+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].x/(3.0*Volume[CurrentSystem]);
          AtomVector[nr_of_coulombic_sites]->y+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].y/(3.0*Volume[CurrentSystem]);
          AtomVector[nr_of_coulombic_sites]->z+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].z/(3.0*Volume[CurrentSystem]);
        }
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_induced_dipole_sites_framework=nr_of_coulombic_sites;

  nr_molecules=NumberOfAdsorbateMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Adsorbates[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Adsorbates[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        charge=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        if(considered_charged)
        {
          DipoleVector[nr_of_coulombic_sites]=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][i].Atoms[j].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          AtomVector[nr_of_coulombic_sites]=&Adsorbates[CurrentSystem][i].Atoms[j].ElectricField;
          if(!OmitAdsorbateAdsorbatePolarization)
          {
            AtomVector[nr_of_coulombic_sites]->x+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].x/(3.0*Volume[CurrentSystem]);
            AtomVector[nr_of_coulombic_sites]->y+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].y/(3.0*Volume[CurrentSystem]);
            AtomVector[nr_of_coulombic_sites]->z+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].z/(3.0*Volume[CurrentSystem]);
          }
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_induced_dipole_sites_adsorbate=nr_of_coulombic_sites;

  nr_molecules=NumberOfCationMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    type_mol=Cations[CurrentSystem][i].Type;
    if(Components[type_mol].HasCharges)
    {
      nr_atoms=Cations[CurrentSystem][i].NumberOfAtoms;
      for(j=0;j<nr_atoms;j++)
      {
        type=Cations[CurrentSystem][i].Atoms[j].Type;
        charge=Cations[CurrentSystem][i].Atoms[j].Charge;
        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
        if(considered_charged)
        {
          DipoleVector[nr_of_coulombic_sites]=Cations[CurrentSystem][i].Atoms[j].InducedDipole;
          Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][i].Atoms[j].Position);
          Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
          AtomVector[nr_of_coulombic_sites]=&Cations[CurrentSystem][i].Atoms[j].ElectricField;
          if(!OmitCationCationPolarization)
          {
            AtomVector[nr_of_coulombic_sites]->x+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].x/(3.0*Volume[CurrentSystem]);
            AtomVector[nr_of_coulombic_sites]->y+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].y/(3.0*Volume[CurrentSystem]);
            AtomVector[nr_of_coulombic_sites]->z+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].z/(3.0*Volume[CurrentSystem]);
          }
          nr_of_coulombic_sites++;
        }
      }
    }
  }
  nr_of_induced_dipole_sites_cation=nr_of_coulombic_sites;


  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }


  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }


  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          rk=kvecs[nvec];

          sum_induced_dipole_framework.re=0.0;
          sum_induced_dipole_framework.im=0.0;
          for(i=0;i<nr_of_induced_dipole_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_induced_dipole_framework.re+=temp*Eikr[i].re;
            sum_induced_dipole_framework.im+=temp*Eikr[i].im;
          }

          sum_induced_dipole_adsorbate.re=0.0;
          sum_induced_dipole_adsorbate.im=0.0;
          for(;i<nr_of_induced_dipole_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_induced_dipole_adsorbate.re+=temp*Eikr[i].re;
            sum_induced_dipole_adsorbate.im+=temp*Eikr[i].im;
          }

          sum_induced_dipole_cation.re=0.0;
          sum_induced_dipole_cation.im=0.0;
          for(;i<nr_of_induced_dipole_sites_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_induced_dipole_cation.re+=temp*Eikr[i].re;
            sum_induced_dipole_cation.im+=temp*Eikr[i].im;
          }

          sum_induced_dipole.re=sum_induced_dipole_framework.re+sum_induced_dipole_adsorbate.re+sum_induced_dipole_cation.re;
          sum_induced_dipole.im=sum_induced_dipole_framework.im+sum_induced_dipole_adsorbate.im+sum_induced_dipole_cation.im;

          // precomputed wavevector dependent pre-factor
          factor=kfactor[nvec];

          // forces on atoms from other charges and bond-dipoles
          temp_sum_induced_dipole.re=sum_induced_dipole.re;
          temp_sum_induced_dipole.im=sum_induced_dipole.im;
          fac_self=1.0;

          if(OmitIntraFrameworkPolarization)
          {
            temp_sum_induced_dipole.re-=sum_induced_dipole_framework.re;
            temp_sum_induced_dipole.im-=sum_induced_dipole_framework.im;
            fac_self=0.0;
          }

          for(i=0;i<nr_of_induced_dipole_sites_framework;i++)
          {
            dipole=DipoleVector[i];
            temp=fac_self*(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);
            fac=2.0*factor*(Eikr[i].re*temp_sum_induced_dipole.re+Eikr[i].im*temp_sum_induced_dipole.im-temp);
            AtomVector[i]->x-=fac*rk.x;
            AtomVector[i]->y-=fac*rk.y;
            AtomVector[i]->z-=fac*rk.z;
          }

          temp_sum_induced_dipole.re=sum_induced_dipole.re;
          temp_sum_induced_dipole.im=sum_induced_dipole.im;
          fac_self=1.0;

          if(OmitAdsorbateAdsorbatePolarization)
          {
            temp_sum_induced_dipole.re-=sum_induced_dipole_adsorbate.re;
            temp_sum_induced_dipole.im-=sum_induced_dipole_adsorbate.im;
            fac_self=0.0;
          }
          if(OmitAdsorbateCationPolarization)
          {
            temp_sum_induced_dipole.re-=sum_induced_dipole_cation.re;
            temp_sum_induced_dipole.im-=sum_induced_dipole_cation.im;
          }

          for(;i<nr_of_induced_dipole_sites_adsorbate;i++)
          {
            dipole=DipoleVector[i];
            temp=fac_self*(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);
            fac=2.0*factor*(Eikr[i].re*temp_sum_induced_dipole.re+Eikr[i].im*temp_sum_induced_dipole.im-temp);
            AtomVector[i]->x-=fac*rk.x;
            AtomVector[i]->y-=fac*rk.y;
            AtomVector[i]->z-=fac*rk.z;
          }

          temp_sum_induced_dipole.re=sum_induced_dipole.re;
          temp_sum_induced_dipole.im=sum_induced_dipole.im;
          fac_self=1.0;

          if(OmitCationCationPolarization)
          {
            temp_sum_induced_dipole.re-=sum_induced_dipole_cation.re;
            temp_sum_induced_dipole.im-=sum_induced_dipole_cation.im;
            fac_self=0.0;
          }
          if(OmitAdsorbateCationPolarization)
          {
            temp_sum_induced_dipole.re-=sum_induced_dipole_adsorbate.re;
            temp_sum_induced_dipole.im-=sum_induced_dipole_adsorbate.im;
          }

          for(;i<nr_of_induced_dipole_sites_cation;i++)
          {
            dipole=DipoleVector[i];
            temp=fac_self*(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);
            fac=2.0*factor*(Eikr[i].re*temp_sum_induced_dipole.re+Eikr[i].im*temp_sum_induced_dipole.im-temp);
            AtomVector[i]->x-=fac*rk.x;
            AtomVector[i]->y-=fac*rk.y;
            AtomVector[i]->z-=fac*rk.z;
          }

          // next wavevector
          nvec++;
        }
      }
    }
  }
}

/*********************************************************************************************************
 * Name       | ComputeElectricFieldFromInducedDipolesEwaldMC                                            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of the induced electric field.                                 *
 * Parameters |                                                                                          *
 * Options    | 'OmitIntraFrameworkPolarization', 'OmitAdsorbateAdsorbatePolarization',                  *
 *            | 'OmitCationCationPolarization', 'OmitAdsorbateCationPolarization'.                       *
 *********************************************************************************************************/

void ComputeElectricFieldFromInducedDipolesEwaldMC(int NewMolecule,int excl_ads,int excl_cation)
{
  int i,j,f1,ii,jj,kk,index_i,index_j,index_k;
  int kmax_x,kmax_y,kmax_z,nvec;
  int nr_of_coulombic_sites,nr_frameworks,nr_molecules,nr_atoms;
  int nr_of_induced_dipole_sites_framework,nr_of_induced_dipole_sites_adsorbate,nr_of_induced_dipole_sites_cation;
  int type,type_mol;
  COMPLEX sum_induced_dipole_framework,sum_induced_dipole_adsorbate,sum_induced_dipole_cation,sum_induced_dipole;
  COMPLEX temp_sum_induced_dipole;
  REAL alpha,fac,temp,fac_self,factor,charge;
  VECTOR *kvecs,pos,rk,dipole;
  REAL *kfactor,recip_cutoff,ksqr;
  int considered_charged;
  

  // return immediately if the Ewald summation is not used
  if(ChargeMethod!=EWALD) return;
  if(OmitEwaldFourier) return;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  nr_of_coulombic_sites=0;
  nr_frameworks=Framework[CurrentSystem].NumberOfFrameworks;
  for(f1=0;f1<nr_frameworks;f1++)
  {
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    for(i=0;i<nr_atoms;i++)
    {
      type=Framework[CurrentSystem].Atoms[f1][i].Type;
      charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged)
      {
        DipoleVector[nr_of_coulombic_sites]=Framework[CurrentSystem].Atoms[f1][i].InducedDipole;
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        AtomVector[nr_of_coulombic_sites]=&Framework[CurrentSystem].Atoms[f1][i].ElectricField;
        if(!OmitIntraFrameworkPolarization)
        {
          AtomVector[nr_of_coulombic_sites]->x+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].x/(3.0*Volume[CurrentSystem]);
          AtomVector[nr_of_coulombic_sites]->y+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].y/(3.0*Volume[CurrentSystem]);
          AtomVector[nr_of_coulombic_sites]->z+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].z/(3.0*Volume[CurrentSystem]);
        }
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_induced_dipole_sites_framework=nr_of_coulombic_sites;

  nr_molecules=NumberOfAdsorbateMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    if(i!=excl_ads)
    {
      type_mol=Adsorbates[CurrentSystem][i].Type;
      if(Components[type_mol].HasCharges)
      {
        nr_atoms=Adsorbates[CurrentSystem][i].NumberOfAtoms;
        for(j=0;j<nr_atoms;j++)
        {
          type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          charge=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
          if(considered_charged)
          {
            DipoleVector[nr_of_coulombic_sites]=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][i].Atoms[j].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            AtomVector[nr_of_coulombic_sites]=&Adsorbates[CurrentSystem][i].Atoms[j].ElectricField;
            if(!OmitAdsorbateAdsorbatePolarization)
            {
              AtomVector[nr_of_coulombic_sites]->x+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].x/(3.0*Volume[CurrentSystem]);
              AtomVector[nr_of_coulombic_sites]->y+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].y/(3.0*Volume[CurrentSystem]);
              AtomVector[nr_of_coulombic_sites]->z+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].z/(3.0*Volume[CurrentSystem]);
            }
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }
  if(NewMolecule&&(!Components[CurrentComponent].ExtraFrameworkMolecule))
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      type=Components[CurrentComponent].Type[i];
      if(PseudoAtoms[type].HasCharges)
      {
        DipoleVector[nr_of_coulombic_sites]=InducedDipoleAtTrialPosition[CurrentSystem][i];
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(TrialPosition[CurrentSystem][i]);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        AtomVector[nr_of_coulombic_sites]=&ElectricFieldAtTrialPosition[CurrentSystem][i];
        if(!OmitAdsorbateAdsorbatePolarization)
        {
          AtomVector[nr_of_coulombic_sites]->x+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].x/(3.0*Volume[CurrentSystem]);
          AtomVector[nr_of_coulombic_sites]->y+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].y/(3.0*Volume[CurrentSystem]);
          AtomVector[nr_of_coulombic_sites]->z+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].z/(3.0*Volume[CurrentSystem]);
        }
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_induced_dipole_sites_adsorbate=nr_of_coulombic_sites;


  nr_molecules=NumberOfCationMolecules[CurrentSystem];
  for(i=0;i<nr_molecules;i++)
  {
    if(i!=excl_cation)
    {
      type_mol=Cations[CurrentSystem][i].Type;
      if(Components[type_mol].HasCharges)
      {
        nr_atoms=Cations[CurrentSystem][i].NumberOfAtoms;
        for(j=0;j<nr_atoms;j++)
        {
          type=Cations[CurrentSystem][i].Atoms[j].Type;
          charge=Cations[CurrentSystem][i].Atoms[j].Charge;
          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
          if(considered_charged)
          {
            DipoleVector[nr_of_coulombic_sites]=Cations[CurrentSystem][i].Atoms[j].InducedDipole;
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][i].Atoms[j].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            AtomVector[nr_of_coulombic_sites]=&Cations[CurrentSystem][i].Atoms[j].ElectricField;
            if(!OmitCationCationPolarization)
            {
              AtomVector[nr_of_coulombic_sites]->x+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].x/(3.0*Volume[CurrentSystem]);
              AtomVector[nr_of_coulombic_sites]->y+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].y/(3.0*Volume[CurrentSystem]);
              AtomVector[nr_of_coulombic_sites]->z+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].z/(3.0*Volume[CurrentSystem]);
            }
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }
  if(NewMolecule&&Components[CurrentComponent].ExtraFrameworkMolecule)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      type=Components[CurrentComponent].Type[i];
      if(PseudoAtoms[type].HasCharges)
      {
        DipoleVector[nr_of_coulombic_sites]=InducedDipoleAtTrialPosition[CurrentSystem][i];
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(TrialPosition[CurrentSystem][i]);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        AtomVector[nr_of_coulombic_sites]=&ElectricFieldAtTrialPosition[CurrentSystem][i];
        if(!OmitCationCationPolarization)
        {
          AtomVector[nr_of_coulombic_sites]->x+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].x/(3.0*Volume[CurrentSystem]);
          AtomVector[nr_of_coulombic_sites]->y+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].y/(3.0*Volume[CurrentSystem]);
          AtomVector[nr_of_coulombic_sites]->z+=COULOMBIC_CONVERSION_FACTOR*4.0*M_PI*DipoleVector[nr_of_coulombic_sites].z/(3.0*Volume[CurrentSystem]);
        }
        nr_of_coulombic_sites++;
      }
    }
  }
  nr_of_induced_dipole_sites_cation=nr_of_coulombic_sites;



  // Calculate the exp(ik.r) terms for all wavevectors
  // =================================================

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index_i=MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=sin(pos.z);

    index_i=-MaxNumberOfCoulombicSites+i;
    Eikx[index_i].re=cos(pos.x); Eikx[index_i].im=-sin(pos.x);
    Eiky[index_i].re=cos(pos.y); Eiky[index_i].im=-sin(pos.y);
    Eikz[index_i].re=cos(pos.z); Eikz[index_i].im=-sin(pos.z);
  }


  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }


  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_i=ii*MaxNumberOfCoulombicSites+i;
        index_j=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_i].re*Eiky[index_j].re-Eikx[index_i].im*Eiky[index_j].im;
        Eikr_xy[i].im=Eikx[index_i].im*Eiky[index_j].re+Eikx[index_i].re*Eiky[index_j].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          rk=kvecs[nvec];

          sum_induced_dipole_framework.re=0.0;
          sum_induced_dipole_framework.im=0.0;
          for(i=0;i<nr_of_induced_dipole_sites_framework;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_induced_dipole_framework.re+=temp*Eikr[i].re;
            sum_induced_dipole_framework.im+=temp*Eikr[i].im;
          }

          sum_induced_dipole_adsorbate.re=0.0;
          sum_induced_dipole_adsorbate.im=0.0;
          for(;i<nr_of_induced_dipole_sites_adsorbate;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_induced_dipole_adsorbate.re+=temp*Eikr[i].re;
            sum_induced_dipole_adsorbate.im+=temp*Eikr[i].im;
          }

          sum_induced_dipole_cation.re=0.0;
          sum_induced_dipole_cation.im=0.0;
          for(;i<nr_of_induced_dipole_sites_cation;i++)
          {
            // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
            index_k=kk*MaxNumberOfCoulombicSites+i;
            Eikr[i].re=Eikr_xy[i].re*Eikz[index_k].re-Eikr_xy[i].im*Eikz[index_k].im;
            Eikr[i].im=Eikr_xy[i].im*Eikz[index_k].re+Eikr_xy[i].re*Eikz[index_k].im;

            dipole=DipoleVector[i];
            temp=dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z;
            sum_induced_dipole_cation.re+=temp*Eikr[i].re;
            sum_induced_dipole_cation.im+=temp*Eikr[i].im;
          }

          sum_induced_dipole.re=sum_induced_dipole_framework.re+sum_induced_dipole_adsorbate.re+sum_induced_dipole_cation.re;
          sum_induced_dipole.im=sum_induced_dipole_framework.im+sum_induced_dipole_adsorbate.im+sum_induced_dipole_cation.im;

          // precomputed wavevector dependent pre-factor
          factor=kfactor[nvec];

          // forces on atoms from other charges and bond-dipoles
          temp_sum_induced_dipole.re=sum_induced_dipole.re;
          temp_sum_induced_dipole.im=sum_induced_dipole.im;
          fac_self=1.0;

          if(OmitIntraFrameworkPolarization)
          {
            temp_sum_induced_dipole.re-=sum_induced_dipole_framework.re;
            temp_sum_induced_dipole.im-=sum_induced_dipole_framework.im;
            fac_self=0.0;
          }

          for(i=0;i<nr_of_induced_dipole_sites_framework;i++)
          {
            dipole=DipoleVector[i];
            temp=fac_self*(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);
            fac=2.0*factor*(Eikr[i].re*temp_sum_induced_dipole.re+Eikr[i].im*temp_sum_induced_dipole.im-temp);
            AtomVector[i]->x-=fac*rk.x;
            AtomVector[i]->y-=fac*rk.y;
            AtomVector[i]->z-=fac*rk.z;
          }

          temp_sum_induced_dipole.re=sum_induced_dipole.re;
          temp_sum_induced_dipole.im=sum_induced_dipole.im;
          fac_self=1.0;
          if(OmitAdsorbateAdsorbatePolarization)
          {
            temp_sum_induced_dipole.re-=sum_induced_dipole_adsorbate.re;
            temp_sum_induced_dipole.im-=sum_induced_dipole_adsorbate.im;
            fac_self=0.0;
          }
          if(OmitAdsorbateCationPolarization)
          {
            temp_sum_induced_dipole.re-=sum_induced_dipole_cation.re;
            temp_sum_induced_dipole.im-=sum_induced_dipole_cation.im;
          }

          for(;i<nr_of_induced_dipole_sites_adsorbate;i++)
          {
            dipole=DipoleVector[i];
            temp=fac_self*(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);
            fac=2.0*factor*(Eikr[i].re*temp_sum_induced_dipole.re+Eikr[i].im*temp_sum_induced_dipole.im-temp);
            AtomVector[i]->x-=fac*rk.x;
            AtomVector[i]->y-=fac*rk.y;
            AtomVector[i]->z-=fac*rk.z;
          }

          temp_sum_induced_dipole.re=sum_induced_dipole.re;
          temp_sum_induced_dipole.im=sum_induced_dipole.im;
          fac_self=1.0;

          if(OmitCationCationPolarization)
          {
            temp_sum_induced_dipole.re-=sum_induced_dipole_cation.re;
            temp_sum_induced_dipole.im-=sum_induced_dipole_cation.im;
            fac_self=0.0;
          }
          if(OmitAdsorbateCationPolarization)
          {
            temp_sum_induced_dipole.re-=sum_induced_dipole_adsorbate.re;
            temp_sum_induced_dipole.im-=sum_induced_dipole_adsorbate.im;
          }

          for(;i<nr_of_induced_dipole_sites_cation;i++)
          {
            dipole=DipoleVector[i];
            temp=fac_self*(dipole.x*rk.x+dipole.y*rk.y+dipole.z*rk.z);
            fac=2.0*factor*(Eikr[i].re*temp_sum_induced_dipole.re+Eikr[i].im*temp_sum_induced_dipole.im-temp);
            AtomVector[i]->x-=fac*rk.x;
            AtomVector[i]->y-=fac*rk.y;
            AtomVector[i]->z-=fac*rk.z;
          }

          // next wavevector
          nvec++;
        }
      }
    }
  }

}



// TODO not working at all
void ComputeInducedDipolesForcesEwald(void)
{
}



// routines for the generalized Hessian matrix (energy minimization, spectra, phonon-dispersion)
// =============================================================================================



// Hessian: Strain - Center of mass (part I)
// =========================================
static inline void HessianAtomicPositionStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,REAL fac,REAL_MATRIX3x3 Theta,VECTOR Rk)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
      // ================================================================
      if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=fac*(Theta.ax*Rk.x+Theta.by*Rk.x+Theta.cz*Rk.x);    // xx x + yy x + zz x
      if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=fac*(Theta.ax*Rk.y+Theta.by*Rk.y+Theta.cz*Rk.y);    // xx y + yy y + zz y
      if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=fac*(Theta.ax*Rk.z+Theta.by*Rk.z+Theta.cz*Rk.z);    // xx z + yy z + zz z
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=fac*Theta.ax*Rk.x;  // xx x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=fac*Theta.by*Rk.x;  // yy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=fac*Theta.cz*Rk.x;  // xz x

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=fac*Theta.ax*Rk.y;  // xx y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=fac*Theta.by*Rk.y;  // yy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=fac*Theta.cz*Rk.y;  // zz y

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=fac*Theta.ax*Rk.z;  // xx z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=fac*Theta.by*Rk.z;  // yy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=fac*Theta.cz*Rk.z;  // zz z
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=fac*Theta.ax*Rk.x;  // xx x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=fac*Theta.ay*Rk.x;  // xy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=fac*Theta.az*Rk.x;  // xz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=fac*Theta.by*Rk.x;  // yy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]-=fac*Theta.bz*Rk.x;  // yz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]-=fac*Theta.cz*Rk.x;  // xz x

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=fac*Theta.ax*Rk.y;  // xx y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=fac*Theta.ay*Rk.y;  // xy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=fac*Theta.az*Rk.y;  // xz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=fac*Theta.by*Rk.y;  // yy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]-=fac*Theta.bz*Rk.y;  // yz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]-=fac*Theta.cz*Rk.y;  // yz y

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=fac*Theta.ax*Rk.z;  // xx z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=fac*Theta.ay*Rk.z;  // xy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=fac*Theta.az*Rk.z;  // xz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=fac*Theta.by*Rk.z;  // yy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]-=fac*Theta.bz*Rk.z;  // yz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]-=fac*Theta.cz*Rk.z;  // zz z
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=fac*Theta.ax*Rk.x;  // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=fac*Theta.by*Rk.x;  // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=fac*Theta.bz*Rk.x;  // yz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=fac*Theta.cz*Rk.x;  // xz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=fac*Theta.ax*Rk.y;  // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=fac*Theta.by*Rk.y;  // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=fac*Theta.bz*Rk.y;  // yz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=fac*Theta.cz*Rk.y;  // yz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=fac*Theta.ax*Rk.z;  // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=fac*Theta.by*Rk.z;  // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=fac*Theta.bz*Rk.z;  // yz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=fac*Theta.cz*Rk.z;  // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=fac*Theta.ax*Rk.x;  // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=fac*Theta.az*Rk.x;  // xz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=fac*Theta.by*Rk.x;  // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=fac*Theta.cz*Rk.x;  // xz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=fac*Theta.ax*Rk.y;  // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=fac*Theta.az*Rk.y;  // xz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=fac*Theta.by*Rk.y;  // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=fac*Theta.cz*Rk.y;  // yz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=fac*Theta.ax*Rk.z;  // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=fac*Theta.az*Rk.z;  // xz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=fac*Theta.by*Rk.z;  // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=fac*Theta.cz*Rk.z;  // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=fac*Theta.ax*Rk.x;  // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=fac*Theta.ay*Rk.x;  // xy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=fac*Theta.by*Rk.x;  // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=fac*Theta.cz*Rk.x;  // xz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=fac*Theta.ax*Rk.y;  // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=fac*Theta.ay*Rk.y;  // xy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=fac*Theta.by*Rk.y;  // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=fac*Theta.cz*Rk.y;  // yz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=fac*Theta.ax*Rk.z;  // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=fac*Theta.ay*Rk.z;  // xy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=fac*Theta.by*Rk.z;  // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=fac*Theta.cz*Rk.z;  // zz z
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

// Hessian: Strain - Center of mass (part II)
// ==========================================
static inline void HessianCenterOfMassStrainI(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,REAL f2_I,VECTOR posA,VECTOR comA,VECTOR Rk)
{
  int n;
  VECTOR dI;

  n=NumberOfCoordinatesMinimizationVariables;

  dI.x=posA.x-comA.x;
  dI.y=posA.y-comA.y;
  dI.z=posA.z-comA.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
      // =================================================================
      if(index_i.x>=0) HessianMatrix.element[index_i.x][n]-=f2_I*(Rk.x*dI.x*Rk.x+Rk.y*dI.y*Rk.x+Rk.z*dI.z*Rk.x); // xx x + yy x + zz x
      if(index_i.y>=0) HessianMatrix.element[index_i.y][n]-=f2_I*(Rk.x*dI.x*Rk.y+Rk.y*dI.y*Rk.y+Rk.z*dI.z*Rk.y); // xx y + yy y + zz y
      if(index_i.z>=0) HessianMatrix.element[index_i.z][n]-=f2_I*(Rk.x*dI.x*Rk.z+Rk.y*dI.y*Rk.z+Rk.z*dI.z*Rk.z); // xx z + yy z + zz z
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=f2_I*Rk.x*dI.x*Rk.x; // xx x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=f2_I*Rk.y*dI.y*Rk.x; // yy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=f2_I*Rk.z*dI.z*Rk.x; // zz x

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=f2_I*Rk.x*dI.x*Rk.y; // xx y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=f2_I*Rk.y*dI.y*Rk.y; // yy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=f2_I*Rk.z*dI.z*Rk.y; // zz y

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=f2_I*Rk.x*dI.x*Rk.z; // xx z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=f2_I*Rk.y*dI.y*Rk.z; // yy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=f2_I*Rk.z*dI.z*Rk.z; // zz z
          break;
        case REGULAR:
          // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=f2_I*Rk.x*dI.x*Rk.x;                 // xx x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.x; // xy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.x; // xz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=f2_I*Rk.y*dI.y*Rk.x;                 // yy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.x; // yz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]-=f2_I*Rk.z*dI.z*Rk.x;                 // zz x

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=f2_I*Rk.x*dI.x*Rk.y;                 // xx y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.y; // xy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.y; // xz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=f2_I*Rk.y*dI.y*Rk.y;                 // yy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.y; // yz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]-=f2_I*Rk.z*dI.z*Rk.y;                 // zz y

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=f2_I*Rk.x*dI.x*Rk.z;                 // xx z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.z; // xy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.z; // xz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=f2_I*Rk.y*dI.y*Rk.z;                 // yy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.z; // yz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]-=f2_I*Rk.z*dI.z*Rk.z;                 // zz z
          break;
        case REGULAR_UPPER_TRIANGLE:
          // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=f2_I*Rk.x*dI.x*Rk.x; // xx x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=f2_I*Rk.x*dI.y*Rk.x; // xy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=f2_I*Rk.x*dI.z*Rk.x; // xz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=f2_I*Rk.y*dI.y*Rk.x; // yy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]-=f2_I*Rk.y*dI.z*Rk.x; // yz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]-=f2_I*Rk.z*dI.z*Rk.x; // zz x

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=f2_I*Rk.x*dI.x*Rk.y; // xx y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=f2_I*Rk.x*dI.y*Rk.y; // xy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=f2_I*Rk.x*dI.z*Rk.y; // xz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=f2_I*Rk.y*dI.y*Rk.y; // yy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]-=f2_I*Rk.y*dI.z*Rk.y; // yz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]-=f2_I*Rk.z*dI.z*Rk.y; // zz y

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=f2_I*Rk.x*dI.x*Rk.z; // xx z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=f2_I*Rk.x*dI.y*Rk.z; // xy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=f2_I*Rk.x*dI.z*Rk.z; // xz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=f2_I*Rk.y*dI.y*Rk.z; // yy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]-=f2_I*Rk.y*dI.z*Rk.z; // yz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]-=f2_I*Rk.z*dI.z*Rk.z; // zz z
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=f2_I*Rk.x*dI.x*Rk.x;                 // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=f2_I*Rk.y*dI.y*Rk.x;                 // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.x; // yz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=f2_I*Rk.z*dI.z*Rk.x;                 // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=f2_I*Rk.x*dI.x*Rk.y;                 // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=f2_I*Rk.y*dI.y*Rk.y;                 // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.y; // yz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=f2_I*Rk.z*dI.z*Rk.y;                 // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=f2_I*Rk.x*dI.x*Rk.z;                 // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=f2_I*Rk.y*dI.y*Rk.z;                 // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.z; // yz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=f2_I*Rk.z*dI.z*Rk.z;                 // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=f2_I*Rk.x*dI.x*Rk.x;                 // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.x; // xz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=f2_I*Rk.y*dI.y*Rk.x;                 // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=f2_I*Rk.z*dI.z*Rk.x;                 // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=f2_I*Rk.x*dI.x*Rk.y;                 // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.y; // xz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=f2_I*Rk.y*dI.y*Rk.y;                 // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=f2_I*Rk.z*dI.z*Rk.y;                 // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=f2_I*Rk.x*dI.x*Rk.z;                 // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.z; // xz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=f2_I*Rk.y*dI.y*Rk.z;                 // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=f2_I*Rk.z*dI.z*Rk.z;                 // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=f2_I*Rk.x*dI.x*Rk.x;                 // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.x; // xy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=f2_I*Rk.y*dI.y*Rk.x;                 // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=f2_I*Rk.z*dI.z*Rk.x;                 // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=f2_I*Rk.x*dI.x*Rk.y;                 // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.y; // xy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=f2_I*Rk.y*dI.y*Rk.y;                 // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=f2_I*Rk.z*dI.z*Rk.y;                 // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=f2_I*Rk.x*dI.x*Rk.z;                 // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.z; // xy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=f2_I*Rk.y*dI.y*Rk.z;                 // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=f2_I*Rk.z*dI.z*Rk.z;                 // zz z
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=f2_I*Rk.x*dI.x*Rk.x; // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=f2_I*Rk.y*dI.y*Rk.x; // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=f2_I*Rk.y*dI.z*Rk.x; // yz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=f2_I*Rk.z*dI.z*Rk.x; // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=f2_I*Rk.x*dI.x*Rk.y; // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=f2_I*Rk.y*dI.y*Rk.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=f2_I*Rk.y*dI.z*Rk.y; // yz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=f2_I*Rk.z*dI.z*Rk.y; // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=f2_I*Rk.x*dI.x*Rk.z; // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=f2_I*Rk.y*dI.y*Rk.z; // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=f2_I*Rk.y*dI.z*Rk.z; // yz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=f2_I*Rk.z*dI.z*Rk.z; // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]-=f2_I*Rk.x*dI.x*Rk.x; // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]-=f2_I*Rk.x*dI.z*Rk.x; // xz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]-=f2_I*Rk.y*dI.y*Rk.x; // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]-=f2_I*Rk.z*dI.z*Rk.x; // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]-=f2_I*Rk.x*dI.x*Rk.y; // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]-=f2_I*Rk.x*dI.z*Rk.y; // xz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]-=f2_I*Rk.y*dI.y*Rk.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]-=f2_I*Rk.z*dI.z*Rk.y; // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]-=f2_I*Rk.x*dI.x*Rk.z; // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]-=f2_I*Rk.x*dI.z*Rk.z; // xz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]-=f2_I*Rk.y*dI.y*Rk.z; // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]-=f2_I*Rk.z*dI.z*Rk.z; // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              if(index_i.x>=0)  HessianMatrix.element[index_i.x][n  ]-=f2_I*Rk.x*dI.x*Rk.x; // xx x
              if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+1]-=f2_I*Rk.x*dI.y*Rk.x; // xy x
              if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+2]-=f2_I*Rk.y*dI.y*Rk.x; // yy x
              if(index_i.x>=0)  HessianMatrix.element[index_i.x][n+3]-=f2_I*Rk.z*dI.z*Rk.x; // zz x

              if(index_i.y>=0)  HessianMatrix.element[index_i.y][n  ]-=f2_I*Rk.x*dI.x*Rk.y; // xx y
              if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+1]-=f2_I*Rk.x*dI.y*Rk.y; // xy y
              if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+2]-=f2_I*Rk.y*dI.y*Rk.y; // yy y
              if(index_i.y>=0)  HessianMatrix.element[index_i.y][n+3]-=f2_I*Rk.z*dI.z*Rk.y; // zz y

              if(index_i.z>=0)  HessianMatrix.element[index_i.z][n  ]-=f2_I*Rk.x*dI.x*Rk.z; // xx z
              if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+1]-=f2_I*Rk.x*dI.y*Rk.z; // xy z
              if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+2]-=f2_I*Rk.y*dI.y*Rk.z; // yy z
              if(index_i.z>=0)  HessianMatrix.element[index_i.z][n+3]-=f2_I*Rk.z*dI.z*Rk.z; // zz z
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

// Hessian: Strain - Center of mass (part III)
// ===========================================
static inline void HessianCenterOfMassStrainJ(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,REAL f2_IJ,VECTOR posB,VECTOR comB,VECTOR Rk)
{
  int n;
  VECTOR dJ;

  n=NumberOfCoordinatesMinimizationVariables;

  dJ.x=posB.x-comB.x;
  dJ.y=posB.y-comB.y;
  dJ.z=posB.z-comB.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
      // ================================================================
      if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2_IJ*(Rk.x*dJ.x*Rk.x+Rk.y*dJ.y*Rk.x+Rk.z*dJ.z*Rk.x); // xx x + yy x + zz x
      if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2_IJ*(Rk.x*dJ.x*Rk.y+Rk.y*dJ.y*Rk.y+Rk.z*dJ.z*Rk.y); // xx y + yy y + zz y
      if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2_IJ*(Rk.x*dJ.x*Rk.z+Rk.y*dJ.y*Rk.z+Rk.z*dJ.z*Rk.z); // xx z + yy z + zz z
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x; // xx x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.x; // yy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2_IJ*Rk.z*dJ.z*Rk.x; // zz x

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y; // xx y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.y; // yy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2_IJ*Rk.z*dJ.z*Rk.y; // zz y

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z; // xx z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.z; // yy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2_IJ*Rk.z*dJ.z*Rk.z; // zz z
          break;
        case REGULAR:
          // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x;                 // xx x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.x; // xy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.x; // xz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.x;                 // yy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.x; // yz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.x;                 // zz x

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y;                 // xx y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.y; // xy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.y; // xz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.y;                 // yy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.y; // yz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.y;                 // zz y

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z;                 // xx z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.z; // xy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.z; // xz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.z;                 // yy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.z; // yz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.z;                 // zz z
          break;
        case REGULAR_UPPER_TRIANGLE:
          // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x; // xx x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.x; // xy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2_IJ*Rk.x*dJ.z*Rk.x; // xz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.x; // yy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=f2_IJ*Rk.y*dJ.z*Rk.x; // yz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.x; // zz x

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y; // xx y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.y; // xy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2_IJ*Rk.x*dJ.z*Rk.y; // xz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.y; // yy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=f2_IJ*Rk.y*dJ.z*Rk.y; // yz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.y; // zz y

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z; // xx z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.z; // xy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2_IJ*Rk.x*dJ.z*Rk.z; // xz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.z; // yy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=f2_IJ*Rk.y*dJ.z*Rk.z; // yz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.z; // zz z
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x;                 // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.x;                 // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.x; // yz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x;                 // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y;                 // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.y;                 // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.y; // yz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y;                 // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z;                 // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.z;                 // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.z; // yz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z;                 // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x;                 // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.x; // xz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.x;                 // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x;                 // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y;                 // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.y; // xz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.y;                 // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y;                 // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z;                 // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.z; // xz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.z;                 // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z;                 // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x;                 // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.x; // xy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.x;                 // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x;                 // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y;                 // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.y; // xy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.y;                 // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y;                 // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z;                 // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.z; // xy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.z;                 // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z;                 // zz z
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x; // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.x; // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2_IJ*Rk.y*dJ.z*Rk.x; // yz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x; // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y; // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2_IJ*Rk.y*dJ.z*Rk.y; // yz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y; // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z; // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.z; // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2_IJ*Rk.y*dJ.z*Rk.z; // yz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z; // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x; // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2_IJ*Rk.x*dJ.z*Rk.x; // xz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.x; // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x; // zz x
          
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y; // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2_IJ*Rk.x*dJ.z*Rk.y; // xz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y; // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z; // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2_IJ*Rk.x*dJ.z*Rk.z; // xz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.z; // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z; // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x; // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.x; // xy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.x; // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x; // zz x
          
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y; // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.y; // xy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y; // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z; // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.z; // xy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.z; // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z; // zz z
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


// Hessian: Strain - Orientation (part I)
// ======================================
static inline void HessianOrientationStrainI(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i2,int index1_rigid,
         REAL f1I,REAL f2I,VECTOR posA,VECTOR comA,VECTOR Rk,REAL_MATRIX3x3 Theta)
{
  int n;
  VECTOR dot_product;
  VECTOR veci1,veci2,veci3;
  REAL temp;
  VECTOR dI;

  n=NumberOfCoordinatesMinimizationVariables;

  dI.x=posA.x-comA.x;
  dI.y=posA.y-comA.y;
  dI.z=posA.z-comA.z;

  veci1=DVecX[index1_rigid];
  veci2=DVecY[index1_rigid];
  veci3=DVecZ[index1_rigid];

  dot_product.x=Rk.x*veci1.x+Rk.y*veci1.y+Rk.z*veci1.z;
  dot_product.y=Rk.x*veci2.x+Rk.y*veci2.y+Rk.z*veci2.z;
  dot_product.z=Rk.x*veci3.x+Rk.y*veci3.y+Rk.z*veci3.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      temp=f1I*Theta.ax+f2I*dI.x*Rk.x;
      if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]-=temp*dot_product.x+f1I*veci1.x*Rk.x; // xx x
      if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]-=temp*dot_product.y+f1I*veci2.x*Rk.x; // xx y
      if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]-=temp*dot_product.z+f1I*veci3.x*Rk.x; // xx z

      temp=f1I*Theta.by+f2I*dI.y*Rk.y;
      if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]-=temp*dot_product.x+f1I*veci1.y*Rk.y; // yy z
      if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]-=temp*dot_product.y+f1I*veci2.y*Rk.y; // yy y
      if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]-=temp*dot_product.z+f1I*veci3.y*Rk.y; // yy z

      temp=f1I*Theta.cz+f2I*dI.z*Rk.z;
      if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]-=temp*dot_product.x+f1I*veci1.z*Rk.z; // zz x
      if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]-=temp*dot_product.y+f1I*veci2.z*Rk.z; // zz y
      if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]-=temp*dot_product.z+f1I*veci3.z*Rk.z; // zz z
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          temp=f1I*Theta.ax+f2I*dI.x*Rk.x;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]-=temp*dot_product.x+f1I*veci1.x*Rk.x;   // xx x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]-=temp*dot_product.y+f1I*veci2.x*Rk.x;   // xx y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]-=temp*dot_product.z+f1I*veci3.x*Rk.x;   // xx z

          temp=f1I*Theta.by+f2I*dI.y*Rk.y;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]-=temp*dot_product.x+f1I*veci1.y*Rk.y; // yy z
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]-=temp*dot_product.y+f1I*veci2.y*Rk.y; // yy y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]-=temp*dot_product.z+f1I*veci3.y*Rk.y; // yy z

          temp=f1I*Theta.cz+f2I*dI.z*Rk.z;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]-=temp*dot_product.x+f1I*veci1.z*Rk.z; // zz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]-=temp*dot_product.y+f1I*veci2.z*Rk.z; // zz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]-=temp*dot_product.z+f1I*veci3.z*Rk.z; // zz z
          break;
        case REGULAR:
          // the first, second, and third term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ===================================================================================
          temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

          temp=f2I*0.5*(Rk.x*dI.y+Rk.y*dI.x)+f1I*Theta.ay;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]-=temp*dot_product.x+f1I*0.5*(Rk.x*veci1.y+Rk.y*veci1.x); // xy x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]-=temp*dot_product.y+f1I*0.5*(Rk.x*veci2.y+Rk.y*veci2.x); // xy y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]-=temp*dot_product.z+f1I*0.5*(Rk.x*veci3.y+Rk.y*veci3.x); // xy z

          temp=f2I*0.5*(Rk.x*dI.z+Rk.z*dI.x)+f1I*Theta.az;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]-=temp*dot_product.x+f1I*0.5*(veci1.z*Rk.x+veci1.x*Rk.z); // xz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]-=temp*dot_product.y+f1I*0.5*(veci2.z*Rk.x+veci2.x*Rk.z); // xz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]-=temp*dot_product.z+f1I*0.5*(veci3.z*Rk.x+veci3.x*Rk.z); // xz z

          temp=f2I*Rk.y*dI.y+f1I*Theta.by;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

          temp=f2I*0.5*(Rk.y*dI.z+Rk.z*dI.y)+f1I*Theta.bz;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+4]-=temp*dot_product.x+f1I*0.5*(Rk.y*veci1.z+Rk.z*veci1.y); // yz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+4]-=temp*dot_product.y+f1I*0.5*(Rk.y*veci2.z+Rk.z*veci2.y); // yz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+4]-=temp*dot_product.z+f1I*0.5*(Rk.y*veci3.z+Rk.z*veci3.y); // yz z

          temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+5]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+5]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+5]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
          break;
        case REGULAR_UPPER_TRIANGLE:
          // the first, second, and third term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ===================================================================================
          temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

          temp=f2I*Rk.x*dI.y+f1I*Theta.ay;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]-=temp*dot_product.x+f1I*Rk.x*veci1.y; // xy x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]-=temp*dot_product.y+f1I*Rk.x*veci2.y; // xy y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]-=temp*dot_product.z+f1I*Rk.x*veci3.y; // xy z

          temp=f2I*Rk.x*dI.z+f1I*Theta.az;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]-=temp*dot_product.x+f1I*Rk.x*veci1.z; // xz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]-=temp*dot_product.y+f1I*Rk.x*veci2.z; // xz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]-=temp*dot_product.z+f1I*Rk.x*veci3.z; // xz z

          temp=f2I*Rk.y*dI.y+f1I*Theta.by;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

          temp=f2I*Rk.y*dI.z+f1I*Theta.bz;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+4]-=temp*dot_product.x+f1I*Rk.y*veci1.z; // yz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+4]-=temp*dot_product.y+f1I*Rk.y*veci2.z; // yz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+4]-=temp*dot_product.z+f1I*Rk.y*veci3.z; // yz z

          temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+5]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+5]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+5]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first, second, and third term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===================================================================================
              temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
              if(index_i2.x>=0)  HessianMatrix.element[index_i2.x][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
              if(index_i2.y>=0)  HessianMatrix.element[index_i2.y][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
              if(index_i2.z>=0)  HessianMatrix.element[index_i2.z][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

              temp=f2I*Rk.y*dI.y+f1I*Theta.by;
              if(index_i2.x>=0)  HessianMatrix.element[index_i2.x][n+1]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
              if(index_i2.y>=0)  HessianMatrix.element[index_i2.y][n+1]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
              if(index_i2.z>=0)  HessianMatrix.element[index_i2.z][n+1]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

              temp=f2I*0.5*(Rk.y*dI.z+Rk.z*dI.y)+f1I*Theta.bz;
              if(index_i2.x>=0)  HessianMatrix.element[index_i2.x][n+2]-=temp*dot_product.x+f1I*0.5*(Rk.y*veci1.z+Rk.z*veci1.y); // yz x
              if(index_i2.y>=0)  HessianMatrix.element[index_i2.y][n+2]-=temp*dot_product.y+f1I*0.5*(Rk.y*veci2.z+Rk.z*veci2.y); // yz y
              if(index_i2.z>=0)  HessianMatrix.element[index_i2.z][n+2]-=temp*dot_product.z+f1I*0.5*(Rk.y*veci3.z+Rk.z*veci3.y); // yz z

              temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
              if(index_i2.x>=0)  HessianMatrix.element[index_i2.x][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
              if(index_i2.y>=0)  HessianMatrix.element[index_i2.y][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
              if(index_i2.z>=0)  HessianMatrix.element[index_i2.z][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the first, second, and third term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===================================================================================
              temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

              temp=f2I*0.5*(Rk.x*dI.z+Rk.z*dI.x)+f1I*Theta.az;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]-=temp*dot_product.x+f1I*0.5*(veci1.z*Rk.x+veci1.x*Rk.z); // xz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]-=temp*dot_product.y+f1I*0.5*(veci2.z*Rk.x+veci2.x*Rk.z); // xz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]-=temp*dot_product.z+f1I*0.5*(veci3.z*Rk.x+veci3.x*Rk.z); // xz z

              temp=f2I*Rk.y*dI.y+f1I*Theta.by;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

              temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first, second, and third term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ===================================================================================
              temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

              temp=f2I*0.5*(Rk.x*dI.y+Rk.y*dI.x)+f1I*Theta.ay;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]-=temp*dot_product.x+f1I*0.5*(Rk.x*veci1.y+Rk.y*veci1.x); // xy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]-=temp*dot_product.y+f1I*0.5*(Rk.x*veci2.y+Rk.y*veci2.x); // xy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]-=temp*dot_product.z+f1I*0.5*(Rk.x*veci3.y+Rk.y*veci3.x); // xy z

              temp=f2I*Rk.y*dI.y+f1I*Theta.by;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

              temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first and second term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

              temp=f2I*Rk.y*dI.y+f1I*Theta.by;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

              temp=f2I*Rk.y*dI.z+f1I*Theta.bz;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]-=temp*dot_product.x+f1I*Rk.y*veci1.z; // yz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]-=temp*dot_product.y+f1I*Rk.y*veci2.z; // yz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]-=temp*dot_product.z+f1I*Rk.y*veci3.z; // yz z

              temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the first and second term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

              temp=f2I*Rk.x*dI.z+f1I*Theta.az;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]-=temp*dot_product.x+f1I*Rk.x*veci1.z; // xz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]-=temp*dot_product.y+f1I*Rk.x*veci2.z; // xz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]-=temp*dot_product.z+f1I*Rk.x*veci3.z; // xz z

              temp=f2I*Rk.y*dI.y+f1I*Theta.by;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

              temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first and second term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

              temp=f2I*Rk.x*dI.y+f1I*Theta.ay;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]-=temp*dot_product.x+f1I*Rk.x*veci1.y; // xy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]-=temp*dot_product.y+f1I*Rk.x*veci2.y; // xy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]-=temp*dot_product.z+f1I*Rk.x*veci3.y; // xy z

              temp=f2I*Rk.y*dI.y+f1I*Theta.by;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

              temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
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

// Hessian: Strain - Orientation (part II)
// =======================================
static inline void HessianOrientationStrainJ(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i2,int index2_rigid,
                                             REAL f2_IJ,VECTOR posB,VECTOR comB,VECTOR Rk)
{
  int n;
  VECTOR dot_product;
  VECTOR vecj1,vecj2,vecj3;
  REAL temp;
  VECTOR dJ;

  n=NumberOfCoordinatesMinimizationVariables;

  dJ.x=posB.x-comB.x;
  dJ.y=posB.y-comB.y;
  dJ.z=posB.z-comB.z;

  vecj1=DVecX[index2_rigid];
  vecj2=DVecY[index2_rigid];
  vecj3=DVecZ[index2_rigid];

  dot_product.x=Rk.x*vecj1.x+Rk.y*vecj1.y+Rk.z*vecj1.z;
  dot_product.y=Rk.x*vecj2.x+Rk.y*vecj2.y+Rk.z*vecj2.z;
  dot_product.z=Rk.x*vecj3.x+Rk.y*vecj3.y+Rk.z*vecj3.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      temp=f2_IJ*Rk.x*dJ.x;
      if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=temp*dot_product.x; // xx x
      if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=temp*dot_product.y; // xx y
      if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=temp*dot_product.z; // xx z

      temp=f2_IJ*Rk.y*dJ.y;
      if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=temp*dot_product.x; // yy x
      if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=temp*dot_product.y; // yy y
      if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=temp*dot_product.z; // yy z

      temp=f2_IJ*Rk.z*dJ.z;
      if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n]+=temp*dot_product.x; // zz x
      if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n]+=temp*dot_product.y; // zz y
      if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n]+=temp*dot_product.z; // zz z
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          temp=f2_IJ*Rk.x*dJ.x;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]+=temp*dot_product.x; // xx x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]+=temp*dot_product.y; // xx y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]+=temp*dot_product.z; // xx z

          temp=f2_IJ*Rk.y*dJ.y;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=temp*dot_product.x; // yy x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=temp*dot_product.y; // yy y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=temp*dot_product.z; // yy z

          temp=f2_IJ*Rk.z*dJ.z;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=temp*dot_product.x; // zz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=temp*dot_product.y; // zz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=temp*dot_product.z; // zz z
          break;
        case REGULAR:
          // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          temp=f2_IJ*Rk.x*dJ.x;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]+=temp*dot_product.x; // xx x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]+=temp*dot_product.y; // xx y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]+=temp*dot_product.z; // xx z

          temp=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x);
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=temp*dot_product.x; // xy x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=temp*dot_product.y; // xy y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=temp*dot_product.z; // xy z

          temp=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x);
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=temp*dot_product.x; // xz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=temp*dot_product.y; // xz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=temp*dot_product.z; // xz z

          temp=f2_IJ*Rk.y*dJ.y;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=temp*dot_product.x; // yy x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=temp*dot_product.y; // yy y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=temp*dot_product.z; // yy z

          temp=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y);
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+4]+=temp*dot_product.x; // yz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+4]+=temp*dot_product.y; // yz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+4]+=temp*dot_product.z; // yz z

          temp=f2_IJ*Rk.z*dJ.z;
          if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+5]+=temp*dot_product.x; // zz x
          if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+5]+=temp*dot_product.y; // zz y
          if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+5]+=temp*dot_product.z; // zz z
          break;
        case REGULAR_UPPER_TRIANGLE:
          // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          temp=f2_IJ*dJ.x*Rk.x;
          if(index_i2.x>=0)  HessianMatrix.element[index_i2.x][n  ]+=temp*dot_product.x; // xx x
          if(index_i2.y>=0)  HessianMatrix.element[index_i2.y][n  ]+=temp*dot_product.y; // xx y
          if(index_i2.z>=0)  HessianMatrix.element[index_i2.z][n  ]+=temp*dot_product.z; // xx z

          temp=f2_IJ*dJ.y*Rk.x;
          if(index_i2.x>=0)  HessianMatrix.element[index_i2.x][n+1]+=temp*dot_product.x; // xy x
          if(index_i2.y>=0)  HessianMatrix.element[index_i2.y][n+1]+=temp*dot_product.y; // xy y
          if(index_i2.z>=0)  HessianMatrix.element[index_i2.z][n+1]+=temp*dot_product.z; // xy z

          temp=f2_IJ*dJ.z*Rk.x;
          if(index_i2.x>=0)  HessianMatrix.element[index_i2.x][n+2]+=temp*dot_product.x; // xz x
          if(index_i2.y>=0)  HessianMatrix.element[index_i2.y][n+2]+=temp*dot_product.y; // xz y
          if(index_i2.z>=0)  HessianMatrix.element[index_i2.z][n+2]+=temp*dot_product.z; // xz z

          temp=f2_IJ*dJ.y*Rk.y;
          if(index_i2.x>=0)  HessianMatrix.element[index_i2.x][n+3]+=temp*dot_product.x; // yy x
          if(index_i2.y>=0)  HessianMatrix.element[index_i2.y][n+3]+=temp*dot_product.y; // yy y
          if(index_i2.z>=0)  HessianMatrix.element[index_i2.z][n+3]+=temp*dot_product.z; // yy z

          temp=f2_IJ*dJ.z*Rk.y;
          if(index_i2.x>=0)  HessianMatrix.element[index_i2.x][n+4]+=temp*dot_product.x; // yz x
          if(index_i2.y>=0)  HessianMatrix.element[index_i2.y][n+4]+=temp*dot_product.y; // yz y
          if(index_i2.z>=0)  HessianMatrix.element[index_i2.z][n+4]+=temp*dot_product.z; // yz z

          temp=f2_IJ*dJ.z*Rk.z;
          if(index_i2.x>=0)  HessianMatrix.element[index_i2.x][n+5]+=temp*dot_product.x; // zz x
          if(index_i2.y>=0)  HessianMatrix.element[index_i2.y][n+5]+=temp*dot_product.y; // zz y
          if(index_i2.z>=0)  HessianMatrix.element[index_i2.z][n+5]+=temp*dot_product.z; // zz z
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp=f2_IJ*Rk.x*dJ.x;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]+=temp*dot_product.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]+=temp*dot_product.y; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]+=temp*dot_product.z; // xx z

              temp=f2_IJ*Rk.y*dJ.y;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=temp*dot_product.x; // yy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=temp*dot_product.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=temp*dot_product.z; // yy z

              temp=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y);
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=temp*dot_product.x; // yz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=temp*dot_product.y; // yz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=temp*dot_product.z; // yz z

              temp=f2_IJ*Rk.z*dJ.z;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=temp*dot_product.x; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=temp*dot_product.y; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=temp*dot_product.z; // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp=f2_IJ*Rk.x*dJ.x;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]+=temp*dot_product.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]+=temp*dot_product.y; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]+=temp*dot_product.z; // xx z

              temp=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x);
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=temp*dot_product.x; // xz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=temp*dot_product.y; // xz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=temp*dot_product.z; // xz z

              temp=f2_IJ*Rk.y*dJ.y;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=temp*dot_product.x; // yy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=temp*dot_product.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=temp*dot_product.z; // yy z

              temp=f2_IJ*Rk.z*dJ.z;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=temp*dot_product.x; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=temp*dot_product.y; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=temp*dot_product.z; // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp=f2_IJ*Rk.x*dJ.x;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]+=temp*dot_product.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]+=temp*dot_product.y; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]+=temp*dot_product.z; // xx z

              temp=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x);
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=temp*dot_product.x; // xy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=temp*dot_product.y; // xy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=temp*dot_product.z; // xy z

              temp=f2_IJ*Rk.y*dJ.y;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=temp*dot_product.x; // yy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=temp*dot_product.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=temp*dot_product.z; // yy z

              temp=f2_IJ*Rk.z*dJ.z;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=temp*dot_product.x; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=temp*dot_product.y; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=temp*dot_product.z; // zz z
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp=f2_IJ*Rk.x*dJ.x;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]+=temp*dot_product.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]+=temp*dot_product.y; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]+=temp*dot_product.z; // xx z

              temp=f2_IJ*Rk.y*dJ.y;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=temp*dot_product.x; // yy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=temp*dot_product.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=temp*dot_product.z; // yy z

              temp=f2_IJ*Rk.y*dJ.z;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=temp*dot_product.x; // yz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=temp*dot_product.y; // yz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=temp*dot_product.z; // yz z

              temp=f2_IJ*Rk.z*dJ.z;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=temp*dot_product.x; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=temp*dot_product.y; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=temp*dot_product.z; // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp=f2_IJ*Rk.x*dJ.x;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]+=temp*dot_product.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]+=temp*dot_product.y; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]+=temp*dot_product.z; // xx z

              temp=f2_IJ*Rk.x*dJ.z;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=temp*dot_product.x; // xz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=temp*dot_product.y; // xz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=temp*dot_product.z; // xz z

              temp=f2_IJ*Rk.y*dJ.y;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=temp*dot_product.x; // yy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=temp*dot_product.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=temp*dot_product.z; // yy z

              temp=f2_IJ*Rk.z*dJ.z;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=temp*dot_product.x; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=temp*dot_product.y; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=temp*dot_product.z; // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp=f2_IJ*Rk.x*dJ.x;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n  ]+=temp*dot_product.x; // xx x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n  ]+=temp*dot_product.y; // xx y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n  ]+=temp*dot_product.z; // xx z

              temp=f2_IJ*Rk.x*dJ.y;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+1]+=temp*dot_product.x; // xy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+1]+=temp*dot_product.y; // xy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+1]+=temp*dot_product.z; // xy z

              temp=f2_IJ*Rk.y*dJ.y;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+2]+=temp*dot_product.x; // yy x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+2]+=temp*dot_product.y; // yy y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+2]+=temp*dot_product.z; // yy z

              temp=f2_IJ*Rk.z*dJ.z;
              if(index_i2.x>=0) HessianMatrix.element[index_i2.x][n+3]+=temp*dot_product.x; // zz x
              if(index_i2.y>=0) HessianMatrix.element[index_i2.y][n+3]+=temp*dot_product.y; // zz y
              if(index_i2.z>=0) HessianMatrix.element[index_i2.z][n+3]+=temp*dot_product.z; // zz z
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


// Hessian: Strain - Strain (part I)
// =================================
static inline void HessianAtomicStrainStrainLocal(REAL_MATRIX HessianMatrix,REAL f,REAL InverseLamdaSquared,
                                             VECTOR Rk,REAL Rksq,REAL_MATRIX3x3 Theta)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT: 
      HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+(1.0+1.0)+
             4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-2.0*(4.0*Rk.x*Rk.x)*InverseLamdaSquared-2.0*Theta.ax);
      HessianMatrix.element[n][n]+=2.0*f*(Theta.ax*Theta.by+
             4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
      HessianMatrix.element[n][n]+=2.0*f*(Theta.ax*Theta.cz+
             4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));
      HessianMatrix.element[n][n]+=f*(Theta.by*Theta.by+(1.0+1.0)+
             4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-2.0*(4.0*Rk.y*Rk.y)*InverseLamdaSquared-2.0*Theta.by);
      HessianMatrix.element[n][n]+=2.0*f*(Theta.by*Theta.cz+
             4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));
      HessianMatrix.element[n][n]+=f*(Theta.cz*Theta.cz+(1.0+1.0)+
             4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-2.0*(4.0*Rk.z*Rk.z)*InverseLamdaSquared-2.0*Theta.cz);
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+(1.0+1.0)+
                 4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-2.0*(4.0*Rk.x*Rk.x)*InverseLamdaSquared-2.0*Theta.ax);
          HessianMatrix.element[n][n]+=2.0*f*(Theta.ax*Theta.by+
                 4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
          HessianMatrix.element[n][n]+=2.0*f*(Theta.ax*Theta.cz+
                 4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));
          HessianMatrix.element[n][n]+=f*(Theta.by*Theta.by+(1.0+1.0)+
                 4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-2.0*(4.0*Rk.y*Rk.y)*InverseLamdaSquared-2.0*Theta.by);
          HessianMatrix.element[n][n]+=2.0*f*(Theta.by*Theta.cz+
                 4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));
          HessianMatrix.element[n][n]+=f*(Theta.cz*Theta.cz+(1.0+1.0)+
                 4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-2.0*(4.0*Rk.z*Rk.z)*InverseLamdaSquared-2.0*Theta.cz);
          break;
        case ANISOTROPIC:
            HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+(1.0+1.0)+
                   4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-2.0*(4.0*Rk.x*Rk.x)*InverseLamdaSquared-2.0*Theta.ax);
            HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.by+
                   4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
            HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.cz+
                   4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));
            HessianMatrix.element[n+1][n+1]+=f*(Theta.by*Theta.by+(1.0+1.0)+
                   4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-2.0*(4.0*Rk.y*Rk.y)*InverseLamdaSquared-2.0*Theta.by);
            HessianMatrix.element[n+1][n+2]+=f*(Theta.by*Theta.cz+
                   4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));
            HessianMatrix.element[n+2][n+2]+=f*(Theta.cz*Theta.cz+(1.0+1.0)+
                   4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-2.0*(4.0*Rk.z*Rk.z)*InverseLamdaSquared-2.0*Theta.cz);
          break;
        case REGULAR:
            // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ================================================================
            HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
            HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.ay+4.0*Rk.x*Rk.x*Rk.x*Rk.y/(SQR(Rksq))-4.0*Rk.y*Rk.x*InverseLamdaSquared-Theta.ay);
            HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.az+4.0*Rk.x*Rk.x*Rk.x*Rk.z/(SQR(Rksq))-4.0*Rk.z*Rk.x*InverseLamdaSquared-Theta.az);
            HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
            HessianMatrix.element[n][n+4]+=f*(Theta.ax*Theta.bz+4.0*Rk.x*Rk.x*Rk.y*Rk.z/(SQR(Rksq)));
            HessianMatrix.element[n][n+5]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+1][n+1]+=f*(Theta.ay*Theta.ay+1.0+4.0*Rk.x*Rk.y*Rk.x*Rk.y/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.x*Rk.x)*InverseLamdaSquared-0.5*(Theta.ax+Theta.by));
            HessianMatrix.element[n+1][n+2]+=f*(Theta.ay*Theta.az+4.0*Rk.x*Rk.y*Rk.x*Rk.z/(SQR(Rksq))-2.0*Rk.y*Rk.z*InverseLamdaSquared-0.5*Theta.bz);
            HessianMatrix.element[n+1][n+3]+=f*(Theta.ay*Theta.by+4.0*Rk.x*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-4.0*Rk.x*Rk.y*InverseLamdaSquared-Theta.ay);
            HessianMatrix.element[n+1][n+4]+=f*(Theta.ay*Theta.bz+4.0*Rk.x*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-2.0*Rk.x*Rk.z*InverseLamdaSquared-0.5*Theta.az);
            HessianMatrix.element[n+1][n+5]+=f*(Theta.ay*Theta.cz+4.0*Rk.x*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+2][n+2]+=f*(Theta.az*Theta.az+1.0+4.0*Rk.x*Rk.z*Rk.x*Rk.z/(SQR(Rksq))-2.0*(Rk.z*Rk.z+Rk.x*Rk.x)*InverseLamdaSquared-0.5*(Theta.ax+Theta.cz));
            HessianMatrix.element[n+2][n+3]+=f*(Theta.az*Theta.by+4.0*Rk.x*Rk.z*Rk.y*Rk.y/(SQR(Rksq)));
            HessianMatrix.element[n+2][n+4]+=f*(Theta.az*Theta.bz+4.0*Rk.x*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*Rk.x*Rk.y*InverseLamdaSquared-0.5*Theta.ay);
            HessianMatrix.element[n+2][n+5]+=f*(Theta.az*Theta.cz+4.0*Rk.x*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.x*Rk.z*InverseLamdaSquared-Theta.az);

            HessianMatrix.element[n+3][n+3]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
            HessianMatrix.element[n+3][n+4]+=f*(Theta.by*Theta.bz+4.0*Rk.y*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);
            HessianMatrix.element[n+3][n+5]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+4][n+4]+=f*(Theta.bz*Theta.bz+1.0+4.0*Rk.y*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.z*Rk.z)*InverseLamdaSquared-0.5*(Theta.by+Theta.cz));
            HessianMatrix.element[n+4][n+5]+=f*(Theta.bz*Theta.cz+4.0*Rk.y*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);

            HessianMatrix.element[n+5][n+5]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
          break;
        case REGULAR_UPPER_TRIANGLE:
            // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ================================================================
            HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
            HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.ay+4.0*Rk.x*Rk.x*Rk.x*Rk.y/(SQR(Rksq))-4.0*Rk.y*Rk.x*InverseLamdaSquared-Theta.ay);
            HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.az+4.0*Rk.x*Rk.x*Rk.x*Rk.z/(SQR(Rksq))-4.0*Rk.z*Rk.x*InverseLamdaSquared-Theta.az);
            HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
            HessianMatrix.element[n][n+4]+=f*(Theta.ax*Theta.bz+4.0*Rk.x*Rk.x*Rk.y*Rk.z/(SQR(Rksq)));
            HessianMatrix.element[n][n+5]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+1][n+1]+=f*(Theta.ay*Theta.ay+1.0+4.0*Rk.x*Rk.y*Rk.x*Rk.y/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.x*Rk.x)*InverseLamdaSquared-Theta.by);
            HessianMatrix.element[n+1][n+2]+=f*(Theta.ay*Theta.az+4.0*Rk.x*Rk.y*Rk.x*Rk.z/(SQR(Rksq))-2.0*(Rk.y*Rk.z)*InverseLamdaSquared-Theta.bz);
            HessianMatrix.element[n+1][n+3]+=f*(Theta.ay*Theta.by+4.0*Rk.x*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-4.0*Rk.x*Rk.y*InverseLamdaSquared);
            HessianMatrix.element[n+1][n+4]+=f*(Theta.ay*Theta.bz+4.0*Rk.x*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-2.0*Rk.x*Rk.z*InverseLamdaSquared);
            HessianMatrix.element[n+1][n+5]+=f*(Theta.ay*Theta.cz+4.0*Rk.x*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+2][n+2]+=f*(Theta.az*Theta.az+1.0+4.0*Rk.x*Rk.z*Rk.x*Rk.z/(SQR(Rksq))-2.0*(Rk.z*Rk.z+Rk.x*Rk.x)*InverseLamdaSquared-Theta.cz);
            HessianMatrix.element[n+2][n+3]+=f*(Theta.az*Theta.by+4.0*Rk.x*Rk.z*Rk.y*Rk.y/(SQR(Rksq)));
            HessianMatrix.element[n+2][n+4]+=f*(Theta.az*Theta.bz+4.0*Rk.x*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*Rk.x*Rk.y*InverseLamdaSquared);
            HessianMatrix.element[n+2][n+5]+=f*(Theta.az*Theta.cz+4.0*Rk.x*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.x*Rk.z*InverseLamdaSquared);

            HessianMatrix.element[n+3][n+3]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
            HessianMatrix.element[n+3][n+4]+=f*(Theta.by*Theta.bz+4.0*Rk.y*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);
            HessianMatrix.element[n+3][n+5]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+4][n+4]+=f*(Theta.bz*Theta.bz+1.0+4.0*Rk.y*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.z*Rk.z)*InverseLamdaSquared-Theta.cz);
            HessianMatrix.element[n+4][n+5]+=f*(Theta.bz*Theta.cz+4.0*Rk.y*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared);

            HessianMatrix.element[n+5][n+5]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
            break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.bz+4.0*Rk.x*Rk.x*Rk.y*Rk.z/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+1][n+2]+=f*(Theta.by*Theta.bz+4.0*Rk.y*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);
              HessianMatrix.element[n+1][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+2][n+2]+=f*(Theta.bz*Theta.bz+1.0+4.0*Rk.y*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.z*Rk.z)*InverseLamdaSquared-0.5*(Theta.by+Theta.cz));
              HessianMatrix.element[n+2][n+3]+=f*(Theta.bz*Theta.cz+4.0*Rk.y*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.az+4.0*Rk.x*Rk.x*Rk.x*Rk.z/(SQR(Rksq))-4.0*Rk.z*Rk.x*InverseLamdaSquared-Theta.az);
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.az*Theta.az+1.0+4.0*Rk.x*Rk.z*Rk.x*Rk.z/(SQR(Rksq))-2.0*(Rk.z*Rk.z+Rk.x*Rk.x)*InverseLamdaSquared-0.5*(Theta.ax+Theta.cz));
              HessianMatrix.element[n+1][n+2]+=f*(Theta.az*Theta.by+4.0*Rk.x*Rk.z*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n+1][n+3]+=f*(Theta.az*Theta.cz+4.0*Rk.x*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.x*Rk.z*InverseLamdaSquared-Theta.az);

              HessianMatrix.element[n+2][n+2]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+2][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.ay+4.0*Rk.x*Rk.x*Rk.x*Rk.y/(SQR(Rksq))-4.0*Rk.y*Rk.x*InverseLamdaSquared-Theta.ay);
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.ay*Theta.ay+1.0+4.0*Rk.x*Rk.y*Rk.x*Rk.y/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.x*Rk.x)*InverseLamdaSquared-0.5*(Theta.ax+Theta.by));
              HessianMatrix.element[n+1][n+2]+=f*(Theta.ay*Theta.by+4.0*Rk.x*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-4.0*Rk.x*Rk.y*InverseLamdaSquared-Theta.ay);
              HessianMatrix.element[n+1][n+3]+=f*(Theta.ay*Theta.cz+4.0*Rk.x*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+2][n+2]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+2][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.bz+4.0*Rk.x*Rk.x*Rk.y*Rk.z/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+1][n+2]+=f*(Theta.by*Theta.bz+4.0*Rk.y*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);
              HessianMatrix.element[n+1][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+2][n+2]+=f*(Theta.bz*Theta.bz+1.0+4.0*Rk.y*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.z*Rk.z)*InverseLamdaSquared-Theta.cz);
              HessianMatrix.element[n+2][n+3]+=f*(Theta.bz*Theta.cz+4.0*Rk.y*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared);

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.az+4.0*Rk.x*Rk.x*Rk.x*Rk.z/(SQR(Rksq))-4.0*Rk.z*Rk.x*InverseLamdaSquared-Theta.az);
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.az*Theta.az+1.0+4.0*Rk.x*Rk.z*Rk.x*Rk.z/(SQR(Rksq))-2.0*(Rk.z*Rk.z+Rk.x*Rk.x)*InverseLamdaSquared-Theta.cz);
              HessianMatrix.element[n+1][n+2]+=f*(Theta.az*Theta.by+4.0*Rk.x*Rk.z*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n+1][n+3]+=f*(Theta.az*Theta.cz+4.0*Rk.x*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.x*Rk.z*InverseLamdaSquared);

              HessianMatrix.element[n+2][n+2]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+2][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.ay+4.0*Rk.x*Rk.x*Rk.x*Rk.y/(SQR(Rksq))-4.0*Rk.y*Rk.x*InverseLamdaSquared-Theta.ay);
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.ay*Theta.ay+1.0+4.0*Rk.x*Rk.y*Rk.x*Rk.y/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.x*Rk.x)*InverseLamdaSquared-Theta.by);
              HessianMatrix.element[n+1][n+2]+=f*(Theta.ay*Theta.by+4.0*Rk.x*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-4.0*Rk.x*Rk.y*InverseLamdaSquared);
              HessianMatrix.element[n+1][n+3]+=f*(Theta.ay*Theta.cz+4.0*Rk.x*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+2][n+2]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+2][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
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
static inline void HessianAtomicCorrectionStrainStrainI(REAL_MATRIX HessianMatrix,REAL f1_I,REAL f2_I,VECTOR posA,VECTOR comA,VECTOR Rk,REAL_MATRIX3x3 Theta)
{
  int n;
  REAL temp1;
  VECTOR dI;

  n=NumberOfCoordinatesMinimizationVariables;

  dI.x=posA.x-comA.x;
  dI.y=posA.y-comA.y;
  dI.z=posA.z-comA.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                         // xx xx
      HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x);   // xx yy
      HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x);   // xx zz
      HessianMatrix.element[n][n]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                     // yy yy
      HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y); // yy zz
      HessianMatrix.element[n][n]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                     // zz zz

      temp1=f2_I*dI.x*Rk.x;
      HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
      HessianMatrix.element[n][n]+=2.0*temp1*dI.y*Rk.y;                   // xx yy
      HessianMatrix.element[n][n]+=2.0*temp1*dI.z*Rk.z;                   // xx zz
      temp1=f2_I*dI.y*Rk.y;
      HessianMatrix.element[n][n]+=temp1*dI.y*Rk.y;                 // yy yy
      HessianMatrix.element[n][n]+=2.0*temp1*dI.z*Rk.z;                 // yy zz
      temp1=f2_I*(posA.z-comA.z)*Rk.z;
      HessianMatrix.element[n][n]+=temp1*dI.z*Rk.z;                 // zz zz

      HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
      HessianMatrix.element[n][n]+=f1_I*Rk.y*dI.y;                  // yy yy
      HessianMatrix.element[n][n]+=f1_I*Rk.z*dI.z;                  // zz zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                         // xx xx
          HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x);   // xx yy
          HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x);   // xx zz
          HessianMatrix.element[n][n]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                     // yy yy
          HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y); // yy zz
          HessianMatrix.element[n][n]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                     // zz zz

          temp1=f2_I*dI.x*Rk.x;
          HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
          HessianMatrix.element[n][n]+=2.0*temp1*dI.y*Rk.y;                   // xx yy
          HessianMatrix.element[n][n]+=2.0*temp1*dI.z*Rk.z;                   // xx zz
          temp1=f2_I*dI.y*Rk.y;
          HessianMatrix.element[n][n]+=temp1*dI.y*Rk.y;                 // yy yy
          HessianMatrix.element[n][n]+=2.0*temp1*dI.z*Rk.z;                 // yy zz
          temp1=f2_I*(posA.z-comA.z)*Rk.z;
          HessianMatrix.element[n][n]+=temp1*dI.z*Rk.z;                 // zz zz

          HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
          HessianMatrix.element[n][n]+=f1_I*Rk.y*dI.y;                  // yy yy
          HessianMatrix.element[n][n]+=f1_I*Rk.z*dI.z;                  // zz zz
          break;
        case ANISOTROPIC:
          HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                         // xx xx
          HessianMatrix.element[n][n+1]+=f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x;   // xx yy
          HessianMatrix.element[n][n+2]+=f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x;   // xx zz
          HessianMatrix.element[n+1][n+1]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                     // yy yy
          HessianMatrix.element[n+1][n+2]+=f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y; // yy zz
          HessianMatrix.element[n+2][n+2]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                     // zz zz

          temp1=f2_I*dI.x*Rk.x;
          HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
          HessianMatrix.element[n][n+1]+=temp1*dI.y*Rk.y;                   // xx yy
          HessianMatrix.element[n][n+2]+=temp1*dI.z*Rk.z;                   // xx zz
          temp1=f2_I*dI.y*Rk.y;
          HessianMatrix.element[n+1][n+1]+=temp1*dI.y*Rk.y;                 // yy yy
          HessianMatrix.element[n+1][n+2]+=temp1*dI.z*Rk.z;                 // yy zz
          temp1=f2_I*(posA.z-comA.z)*Rk.z;
          HessianMatrix.element[n+2][n+2]+=temp1*dI.z*Rk.z;                 // zz zz

          HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
          HessianMatrix.element[n+1][n+1]+=f1_I*Rk.y*dI.y;                  // yy yy
          HessianMatrix.element[n+2][n+2]+=f1_I*Rk.z*dI.z;                  // zz zz
          break;
        case REGULAR:
          // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                                       // xx xx
          HessianMatrix.element[n][n+1]+=0.5*f1_I*Theta.ax*(dI.y*Rk.x+dI.x*Rk.y)+f1_I*Theta.ay*dI.x*Rk.x; // xx xy
          HessianMatrix.element[n][n+2]+=0.5*f1_I*Theta.ax*(dI.z*Rk.x+dI.x*Rk.z)+f1_I*Theta.az*dI.x*Rk.x; // xx xz
          HessianMatrix.element[n][n+3]+=f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x;                 // xx yy
          HessianMatrix.element[n][n+4]+=0.5*f1_I*Theta.ax*(dI.z*Rk.y+dI.y*Rk.z)+f1_I*Theta.bz*dI.x*Rk.x; // xx yz
          HessianMatrix.element[n][n+5]+=f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x;                 // xx zz

          HessianMatrix.element[n+1][n+1]+=f1_I*Theta.ay*(dI.y*Rk.x+dI.x*Rk.y);                                      // xy xy
          HessianMatrix.element[n+1][n+2]+=0.5*f1_I*(Theta.ay*(dI.z*Rk.x+dI.x*Rk.z)+Theta.az*(dI.y*Rk.x+dI.x*Rk.y)); // xy xz
          HessianMatrix.element[n+1][n+3]+=f1_I*Theta.ay*dI.y*Rk.y+0.5*f1_I*Theta.by*(dI.y*Rk.x+dI.x*Rk.y);          // xy yy
          HessianMatrix.element[n+1][n+4]+=0.5*f1_I*(Theta.ay*(dI.z*Rk.y+dI.y*Rk.z)+Theta.bz*(dI.y*Rk.x+dI.x*Rk.y)); // xy yz
          HessianMatrix.element[n+1][n+5]+=f1_I*Theta.ay*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.y*Rk.x+dI.x*Rk.y);          // xy zz

          HessianMatrix.element[n+2][n+2]+=f1_I*Theta.az*(dI.z*Rk.x+dI.x*Rk.z);                                      // xz xz
          HessianMatrix.element[n+2][n+3]+=f1_I*Theta.az*dI.y*Rk.y+0.5*f1_I*Theta.by*(dI.z*Rk.x+dI.x*Rk.z);          // xz yy
          HessianMatrix.element[n+2][n+4]+=0.5*f1_I*(Theta.az*(dI.z*Rk.y+dI.y*Rk.z)+Theta.bz*(dI.z*Rk.x+dI.x*Rk.z)); // xz yz
          HessianMatrix.element[n+2][n+5]+=f1_I*Theta.az*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.z*Rk.x+dI.x*Rk.z);          // xz zz

          HessianMatrix.element[n+3][n+3]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                                     // yy yy
          HessianMatrix.element[n+3][n+4]+=0.5*f1_I*Theta.by*(dI.z*Rk.y+dI.y*Rk.z)+f1_I*Theta.bz*dI.y*Rk.y; // yy yz
          HessianMatrix.element[n+3][n+5]+=f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y;                 // yy zz

          HessianMatrix.element[n+4][n+4]+=f1_I*Theta.bz*(dI.z*Rk.y+dI.y*Rk.z);                             // yz yz
          HessianMatrix.element[n+4][n+5]+=f1_I*Theta.bz*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.z*Rk.y+dI.y*Rk.z); // yz zz

          HessianMatrix.element[n+5][n+5]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                                     // zz zz


          // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          temp1=f2_I*dI.x*Rk.x;
          HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
          HessianMatrix.element[n][n+1]+=temp1*0.5*(dI.y*Rk.x+dI.x*Rk.y);   // xx xy
          HessianMatrix.element[n][n+2]+=temp1*0.5*(dI.z*Rk.x+dI.x*Rk.z);   // xx xz
          HessianMatrix.element[n][n+3]+=temp1*dI.y*Rk.y;                   // xx yy
          HessianMatrix.element[n][n+4]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z);   // xx yz
          HessianMatrix.element[n][n+5]+=temp1*dI.z*Rk.z;                   // xx zz

          temp1=f2_I*0.5*(dI.x*Rk.y+dI.y*Rk.x);
          HessianMatrix.element[n+1][n+1]+=temp1*0.5*(dI.y*Rk.x+dI.x*Rk.y); // xy xy
          HessianMatrix.element[n+1][n+2]+=temp1*0.5*(dI.z*Rk.x+dI.x*Rk.z); // xy xz
          HessianMatrix.element[n+1][n+3]+=temp1*dI.y*Rk.y;                 // xy yy
          HessianMatrix.element[n+1][n+4]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // xy yz
          HessianMatrix.element[n+1][n+5]+=temp1*dI.z*Rk.z;                 // xy zz

          temp1=f2_I*0.5*(dI.x*Rk.z+dI.z*Rk.x);
          HessianMatrix.element[n+2][n+2]+=temp1*0.5*(dI.z*Rk.x+dI.x*Rk.z); // xz xz
          HessianMatrix.element[n+2][n+3]+=temp1*dI.y*Rk.y;                 // xz yy
          HessianMatrix.element[n+2][n+4]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // xz yz
          HessianMatrix.element[n+2][n+5]+=temp1*dI.z*Rk.z;                 // xz zz

          temp1=f2_I*dI.y*Rk.y;
          HessianMatrix.element[n+3][n+3]+=temp1*dI.y*Rk.y;                 // yy yy
          HessianMatrix.element[n+3][n+4]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // yy yz
          HessianMatrix.element[n+3][n+5]+=temp1*dI.z*Rk.z;                 // yy zz

          temp1=f2_I*0.5*(dI.y*Rk.z+dI.z*Rk.y);
          HessianMatrix.element[n+4][n+4]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // yz yz
          HessianMatrix.element[n+4][n+5]+=temp1*dI.z*Rk.z;                 // yz zz

          temp1=f2_I*(posA.z-comA.z)*Rk.z; 
          HessianMatrix.element[n+5][n+5]+=temp1*dI.z*Rk.z;                 // zz zz


          // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
          HessianMatrix.element[n][n+1]+=0.5*f1_I*Rk.x*dI.y;                // xx xy
          HessianMatrix.element[n][n+2]+=0.5*f1_I*Rk.x*dI.z;                // xx xz

          HessianMatrix.element[n+1][n+1]+=0.25*f1_I*(Rk.y*dI.y+Rk.x*dI.x); // xy xy
          HessianMatrix.element[n+1][n+2]+=0.25*f1_I*Rk.y*dI.z;             // xy xz
          HessianMatrix.element[n+1][n+3]+=0.5*f1_I*Rk.x*dI.y;              // xy yy
          HessianMatrix.element[n+1][n+4]+=0.25*f1_I*Rk.x*dI.z;             // xy yz

          HessianMatrix.element[n+2][n+2]+=0.25*f1_I*(Rk.z*dI.z+Rk.x*dI.x); // xz xz
          HessianMatrix.element[n+2][n+4]+=0.25*f1_I*Rk.x*dI.y;             // xz yz
          HessianMatrix.element[n+2][n+5]+=0.5*f1_I*Rk.x*dI.z;              // xz zz

          HessianMatrix.element[n+3][n+3]+=f1_I*Rk.y*dI.y;                  // yy yy
          HessianMatrix.element[n+3][n+4]+=0.5*f1_I*Rk.y*dI.z;              // yy yz

          HessianMatrix.element[n+4][n+4]+=0.25*f1_I*(Rk.z*dI.z+Rk.y*dI.y); // yz yz
          HessianMatrix.element[n+4][n+5]+=0.5*f1_I*Rk.y*dI.z;              // yz zz

          HessianMatrix.element[n+5][n+5]+=f1_I*Rk.z*dI.z;                  // zz zz
          break;
        case REGULAR_UPPER_TRIANGLE:
          // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          HessianMatrix.element[n][n  ]+=f1_I*(Theta.ax*dI.x*Rk.x+Rk.x*dI.x*Theta.ax); // xx xx
          HessianMatrix.element[n][n+1]+=f1_I*(Theta.ax*dI.y*Rk.x+Rk.x*dI.x*Theta.ay); // xx xy
          HessianMatrix.element[n][n+2]+=f1_I*(Theta.ax*dI.z*Rk.x+Rk.x*dI.x*Theta.az); // xx xz
          HessianMatrix.element[n][n+3]+=f1_I*(Theta.ax*dI.y*Rk.y+Rk.x*dI.x*Theta.by); // xx yy
          HessianMatrix.element[n][n+4]+=f1_I*(Theta.ax*dI.z*Rk.y+Rk.x*dI.x*Theta.bz); // xx yz
          HessianMatrix.element[n][n+5]+=f1_I*(Theta.ax*dI.z*Rk.z+Rk.x*dI.x*Theta.cz); // xx zz

          HessianMatrix.element[n+1][n+1]+=f1_I*(Theta.ay*Rk.x*dI.y+Rk.x*dI.y*Theta.ay); // xy xy
          HessianMatrix.element[n+1][n+2]+=f1_I*(Theta.ay*Rk.x*dI.z+Rk.x*dI.y*Theta.az); // xy xz
          HessianMatrix.element[n+1][n+3]+=f1_I*(Theta.ay*Rk.y*dI.y+Rk.x*dI.y*Theta.by); // xy yy
          HessianMatrix.element[n+1][n+4]+=f1_I*(Theta.ay*Rk.y*dI.z+Rk.x*dI.y*Theta.bz); // xy yz
          HessianMatrix.element[n+1][n+5]+=f1_I*(Theta.ay*Rk.z*dI.z+Rk.x*dI.y*Theta.cz); // xy zz

          HessianMatrix.element[n+2][n+2]+=f1_I*(Theta.az*Rk.x*dI.z+Rk.x*dI.z*Theta.az); // xz xz
          HessianMatrix.element[n+2][n+3]+=f1_I*(Theta.az*Rk.y*dI.y+Rk.x*dI.z*Theta.by); // xz yy
          HessianMatrix.element[n+2][n+4]+=f1_I*(Theta.az*Rk.y*dI.z+Rk.x*dI.z*Theta.bz); // xz yz
          HessianMatrix.element[n+2][n+5]+=f1_I*(Theta.az*Rk.z*dI.z+Rk.x*dI.z*Theta.cz); // xz zz

          HessianMatrix.element[n+3][n+3]+=f1_I*(Theta.by*Rk.y*dI.y+Rk.y*dI.y*Theta.by); // yy yy
          HessianMatrix.element[n+3][n+4]+=f1_I*(Theta.by*Rk.y*dI.z+Rk.y*dI.y*Theta.bz); // yy yz
          HessianMatrix.element[n+3][n+5]+=f1_I*(Theta.by*Rk.z*dI.z+Rk.y*dI.y*Theta.cz); // yy zz

          HessianMatrix.element[n+4][n+4]+=f1_I*(Theta.bz*Rk.y*dI.z+Rk.y*dI.z*Theta.bz); // yz yz
          HessianMatrix.element[n+4][n+5]+=f1_I*(Theta.bz*Rk.z*dI.z+Rk.y*dI.z*Theta.cz); // yz zz

          HessianMatrix.element[n+5][n+5]+=f1_I*(Theta.cz*Rk.z*dI.z+Rk.z*dI.z*Theta.cz); // zz zz


          // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          temp1=f2_I*Rk.x*dI.x;
          HessianMatrix.element[n][n  ]+=temp1*Rk.x*dI.x;   // xx xx
          HessianMatrix.element[n][n+1]+=temp1*Rk.x*dI.y;   // xx xy
          HessianMatrix.element[n][n+2]+=temp1*Rk.x*dI.z;   // xx xz
          HessianMatrix.element[n][n+3]+=temp1*Rk.y*dI.y;   // xx yy
          HessianMatrix.element[n][n+4]+=temp1*Rk.y*dI.z;   // xx yz
          HessianMatrix.element[n][n+5]+=temp1*Rk.z*dI.z;   // xx zz

          temp1=f2_I*Rk.x*dI.y;
          HessianMatrix.element[n+1][n+1]+=temp1*Rk.x*dI.y; // xy xy
          HessianMatrix.element[n+1][n+2]+=temp1*Rk.x*dI.z; // xy xz
          HessianMatrix.element[n+1][n+3]+=temp1*Rk.y*dI.y; // xy yy
          HessianMatrix.element[n+1][n+4]+=temp1*Rk.y*dI.z; // xy yz
          HessianMatrix.element[n+1][n+5]+=temp1*Rk.z*dI.z; // xy zz

          temp1=f2_I*Rk.x*dI.z;
          HessianMatrix.element[n+2][n+2]+=temp1*Rk.x*dI.z; // xy xz
          HessianMatrix.element[n+2][n+3]+=temp1*Rk.y*dI.y; // xy yy
          HessianMatrix.element[n+2][n+4]+=temp1*Rk.y*dI.z; // xy yz
          HessianMatrix.element[n+2][n+5]+=temp1*Rk.z*dI.z; // xy zz

          temp1=f2_I*Rk.y*dI.y;
          HessianMatrix.element[n+3][n+3]+=temp1*Rk.y*dI.y; // xy yy
          HessianMatrix.element[n+3][n+4]+=temp1*Rk.y*dI.z; // xy yz
          HessianMatrix.element[n+3][n+5]+=temp1*Rk.z*dI.z; // xy zz

          temp1=f2_I*Rk.y*dI.z;
          HessianMatrix.element[n+4][n+4]+=temp1*Rk.y*dI.z; // xy yz
          HessianMatrix.element[n+4][n+5]+=temp1*Rk.z*dI.z; // xy zz

          temp1=f2_I*Rk.z*dI.z;
          HessianMatrix.element[n+5][n+5]+=temp1*Rk.z*dI.z; // xy zz

          // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;     // xx xx
          HessianMatrix.element[n][n+1]+=f1_I*Rk.x*dI.y;   // xx xy
          HessianMatrix.element[n][n+2]+=f1_I*Rk.x*dI.z;   // xx xz

          HessianMatrix.element[n+1][n+3]+=f1_I*Rk.x*dI.y; // xy yy
          HessianMatrix.element[n+1][n+4]+=f1_I*Rk.x*dI.z; // xy yz

          HessianMatrix.element[n+2][n+5]+=f1_I*Rk.x*dI.z; // xz zz

          HessianMatrix.element[n+3][n+3]+=f1_I*Rk.y*dI.y; // yy yy
          HessianMatrix.element[n+3][n+4]+=f1_I*Rk.y*dI.z; // yy yz

          HessianMatrix.element[n+4][n+5]+=f1_I*Rk.y*dI.z; // yz zz

          HessianMatrix.element[n+5][n+5]+=f1_I*Rk.z*dI.z; // zz zz
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                                       // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x;                 // xx yy
              HessianMatrix.element[n][n+2]+=0.5*f1_I*Theta.ax*(dI.z*Rk.y+dI.y*Rk.z)+f1_I*Theta.bz*dI.x*Rk.x; // xx yz
              HessianMatrix.element[n][n+3]+=f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x;                 // xx zz

              HessianMatrix.element[n+1][n+1]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                                     // yy yy
              HessianMatrix.element[n+1][n+2]+=0.5*f1_I*Theta.by*(dI.z*Rk.y+dI.y*Rk.z)+f1_I*Theta.bz*dI.y*Rk.y; // yy yz
              HessianMatrix.element[n+1][n+3]+=f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y;                 // yy zz

              HessianMatrix.element[n+2][n+2]+=f1_I*Theta.bz*(dI.z*Rk.y+dI.y*Rk.z);                             // yz yz
              HessianMatrix.element[n+2][n+3]+=f1_I*Theta.bz*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.z*Rk.y+dI.y*Rk.z); // yz zz

              HessianMatrix.element[n+3][n+3]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                                     // zz zz


              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*dI.x*Rk.x;
              HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]+=temp1*dI.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+2]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z);   // xx yz
              HessianMatrix.element[n][n+3]+=temp1*dI.z*Rk.z;                   // xx zz

              temp1=f2_I*dI.y*Rk.y;
              HessianMatrix.element[n+1][n+1]+=temp1*dI.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+1][n+2]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // yy yz
              HessianMatrix.element[n+1][n+3]+=temp1*dI.z*Rk.z;                 // yy zz

              temp1=f2_I*0.5*(dI.y*Rk.z+dI.z*Rk.y);
              HessianMatrix.element[n+2][n+2]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // yz yz
              HessianMatrix.element[n+2][n+3]+=temp1*dI.z*Rk.z;                 // yz zz

              temp1=f2_I*(posA.z-comA.z)*Rk.z; 
              HessianMatrix.element[n+3][n+3]+=temp1*dI.z*Rk.z;                 // zz zz


              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx

              HessianMatrix.element[n+1][n+1]+=f1_I*Rk.y*dI.y;                  // yy yy
              HessianMatrix.element[n+1][n+2]+=0.5*f1_I*Rk.y*dI.z;              // yy yz

              HessianMatrix.element[n+2][n+2]+=0.25*f1_I*(Rk.z*dI.z+Rk.y*dI.y); // yz yz
              HessianMatrix.element[n+2][n+3]+=0.5*f1_I*Rk.y*dI.z;              // yz zz

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z;                  // zz zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                                       // xx xx
              HessianMatrix.element[n][n+1]+=0.5*f1_I*Theta.ax*(dI.z*Rk.x+dI.x*Rk.z)+f1_I*Theta.az*dI.x*Rk.x; // xx xz
              HessianMatrix.element[n][n+2]+=f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x;                 // xx yy
              HessianMatrix.element[n][n+3]+=f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x;                 // xx zz

              HessianMatrix.element[n+1][n+1]+=f1_I*Theta.az*(dI.z*Rk.x+dI.x*Rk.z);                                      // xz xz
              HessianMatrix.element[n+1][n+2]+=f1_I*Theta.az*dI.y*Rk.y+0.5*f1_I*Theta.by*(dI.z*Rk.x+dI.x*Rk.z);          // xz yy
              HessianMatrix.element[n+1][n+3]+=f1_I*Theta.az*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.z*Rk.x+dI.x*Rk.z);          // xz zz

              HessianMatrix.element[n+2][n+2]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                                     // yy yy
              HessianMatrix.element[n+2][n+3]+=f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y;                 // yy zz

              HessianMatrix.element[n+3][n+3]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                                     // zz zz


              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*dI.x*Rk.x;
              HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]+=temp1*0.5*(dI.z*Rk.x+dI.x*Rk.z);   // xx xz
              HessianMatrix.element[n][n+2]+=temp1*dI.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+3]+=temp1*dI.z*Rk.z;                   // xx zz

              temp1=f2_I*0.5*(dI.x*Rk.z+dI.z*Rk.x);
              HessianMatrix.element[n+1][n+1]+=temp1*0.5*(dI.z*Rk.x+dI.x*Rk.z); // xz xz
              HessianMatrix.element[n+1][n+2]+=temp1*dI.y*Rk.y;                 // xz yy
              HessianMatrix.element[n+1][n+3]+=temp1*dI.z*Rk.z;                 // xz zz

              temp1=f2_I*dI.y*Rk.y;
              HessianMatrix.element[n+2][n+2]+=temp1*dI.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+2][n+3]+=temp1*dI.z*Rk.z;                 // yy zz

              temp1=f2_I*(posA.z-comA.z)*Rk.z; 
              HessianMatrix.element[n+3][n+3]+=temp1*dI.z*Rk.z;                 // zz zz


              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
              HessianMatrix.element[n][n+1]+=0.5*f1_I*Rk.x*dI.z;                // xx xz

              HessianMatrix.element[n+1][n+1]+=0.25*f1_I*(Rk.z*dI.z+Rk.x*dI.x); // xz xz
              HessianMatrix.element[n+1][n+3]+=0.5*f1_I*Rk.x*dI.z;              // xz zz

              HessianMatrix.element[n+2][n+2]+=f1_I*Rk.y*dI.y;                  // yy yy

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z;                  // zz zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                                       // xx xx
              HessianMatrix.element[n][n+1]+=0.5*f1_I*Theta.ax*(dI.y*Rk.x+dI.x*Rk.y)+f1_I*Theta.ay*dI.x*Rk.x; // xx xy
              HessianMatrix.element[n][n+2]+=f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x;                 // xx yy
              HessianMatrix.element[n][n+3]+=f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x;                 // xx zz

              HessianMatrix.element[n+1][n+1]+=f1_I*Theta.ay*(dI.y*Rk.x+dI.x*Rk.y);                                      // xy xy
              HessianMatrix.element[n+1][n+2]+=f1_I*Theta.ay*dI.y*Rk.y+0.5*f1_I*Theta.by*(dI.y*Rk.x+dI.x*Rk.y);          // xy yy
              HessianMatrix.element[n+1][n+3]+=f1_I*Theta.ay*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.y*Rk.x+dI.x*Rk.y);          // xy zz

              HessianMatrix.element[n+2][n+2]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                                     // yy yy
              HessianMatrix.element[n+2][n+3]+=f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y;                 // yy zz

              HessianMatrix.element[n+3][n+3]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                                     // zz zz

              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*dI.x*Rk.x;
              HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]+=temp1*0.5*(dI.y*Rk.x+dI.x*Rk.y);   // xx xy
              HessianMatrix.element[n][n+2]+=temp1*dI.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+3]+=temp1*dI.z*Rk.z;                   // xx zz

              temp1=f2_I*0.5*(dI.x*Rk.y+dI.y*Rk.x);
              HessianMatrix.element[n+1][n+1]+=temp1*0.5*(dI.y*Rk.x+dI.x*Rk.y); // xy xy
              HessianMatrix.element[n+1][n+2]+=temp1*dI.y*Rk.y;                 // xy yy
              HessianMatrix.element[n+1][n+3]+=temp1*dI.z*Rk.z;                 // xy zz

              temp1=f2_I*dI.y*Rk.y;
              HessianMatrix.element[n+2][n+2]+=temp1*dI.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+2][n+3]+=temp1*dI.z*Rk.z;                 // yy zz

              temp1=f2_I*(posA.z-comA.z)*Rk.z; 
              HessianMatrix.element[n+3][n+3]+=temp1*dI.z*Rk.z;                 // zz zz


              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
              HessianMatrix.element[n][n+1]+=0.5*f1_I*Rk.x*dI.y;                // xx xy

              HessianMatrix.element[n+1][n+1]+=0.25*f1_I*(Rk.y*dI.y+Rk.x*dI.x); // xy xy
              HessianMatrix.element[n+1][n+2]+=0.5*f1_I*Rk.x*dI.y;              // xy yy

              HessianMatrix.element[n+2][n+2]+=f1_I*Rk.y*dI.y;                  // yy yy

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z;                  // zz zz
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n  ]+=f1_I*(Theta.ax*dI.x*Rk.x+Rk.x*dI.x*Theta.ax); // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*(Theta.ax*dI.y*Rk.y+Rk.x*dI.x*Theta.by); // xx yy
              HessianMatrix.element[n][n+2]+=f1_I*(Theta.ax*dI.z*Rk.y+Rk.x*dI.x*Theta.bz); // xx yz
              HessianMatrix.element[n][n+3]+=f1_I*(Theta.ax*dI.z*Rk.z+Rk.x*dI.x*Theta.cz); // xx zz

              HessianMatrix.element[n+1][n+1]+=f1_I*(Theta.by*Rk.y*dI.y+Rk.y*dI.y*Theta.by); // yy yy
              HessianMatrix.element[n+1][n+2]+=f1_I*(Theta.by*Rk.y*dI.z+Rk.y*dI.y*Theta.bz); // yy yz
              HessianMatrix.element[n+1][n+3]+=f1_I*(Theta.by*Rk.z*dI.z+Rk.y*dI.y*Theta.cz); // yy zz

              HessianMatrix.element[n+2][n+2]+=f1_I*(Theta.bz*Rk.y*dI.z+Rk.y*dI.z*Theta.bz); // yz yz
              HessianMatrix.element[n+2][n+3]+=f1_I*(Theta.bz*Rk.z*dI.z+Rk.y*dI.z*Theta.cz); // yz zz

              HessianMatrix.element[n+3][n+3]+=f1_I*(Theta.cz*Rk.z*dI.z+Rk.z*dI.z*Theta.cz); // zz zz


              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*Rk.x*dI.x;
              HessianMatrix.element[n][n  ]+=temp1*Rk.x*dI.x;   // xx xx
              HessianMatrix.element[n][n+1]+=temp1*Rk.y*dI.y;   // xx yy
              HessianMatrix.element[n][n+2]+=temp1*Rk.y*dI.z;   // xx yz
              HessianMatrix.element[n][n+3]+=temp1*Rk.z*dI.z;   // xx zz

              temp1=f2_I*Rk.y*dI.y;
              HessianMatrix.element[n+1][n+1]+=temp1*Rk.y*dI.y; // xy yy
              HessianMatrix.element[n+1][n+2]+=temp1*Rk.y*dI.z; // xy yz
              HessianMatrix.element[n+1][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.y*dI.z;
              HessianMatrix.element[n+2][n+2]+=temp1*Rk.y*dI.z; // xy yz
              HessianMatrix.element[n+2][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]+=temp1*Rk.z*dI.z; // xy zz

              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;     // xx xx

              HessianMatrix.element[n+1][n+1]+=f1_I*Rk.y*dI.y; // yy yy
              HessianMatrix.element[n+1][n+2]+=f1_I*Rk.y*dI.z; // yy yz

              HessianMatrix.element[n+2][n+3]+=f1_I*Rk.y*dI.z; // yz zz

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z; // zz zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n  ]+=f1_I*(Theta.ax*dI.x*Rk.x+Rk.x*dI.x*Theta.ax); // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*(Theta.ax*dI.z*Rk.x+Rk.x*dI.x*Theta.az); // xx xz
              HessianMatrix.element[n][n+2]+=f1_I*(Theta.ax*dI.y*Rk.y+Rk.x*dI.x*Theta.by); // xx yy
              HessianMatrix.element[n][n+3]+=f1_I*(Theta.ax*dI.z*Rk.z+Rk.x*dI.x*Theta.cz); // xx zz

              HessianMatrix.element[n+1][n+1]+=f1_I*(Theta.az*Rk.x*dI.z+Rk.x*dI.z*Theta.az); // xz xz
              HessianMatrix.element[n+1][n+2]+=f1_I*(Theta.az*Rk.y*dI.y+Rk.x*dI.z*Theta.by); // xz yy
              HessianMatrix.element[n+1][n+3]+=f1_I*(Theta.az*Rk.z*dI.z+Rk.x*dI.z*Theta.cz); // xz zz

              HessianMatrix.element[n+2][n+2]+=f1_I*(Theta.by*Rk.y*dI.y+Rk.y*dI.y*Theta.by); // yy yy
              HessianMatrix.element[n+2][n+3]+=f1_I*(Theta.by*Rk.z*dI.z+Rk.y*dI.y*Theta.cz); // yy zz

              HessianMatrix.element[n+3][n+3]+=f1_I*(Theta.cz*Rk.z*dI.z+Rk.z*dI.z*Theta.cz); // zz zz


              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*Rk.x*dI.x;
              HessianMatrix.element[n][n  ]+=temp1*Rk.x*dI.x;   // xx xx
              HessianMatrix.element[n][n+1]+=temp1*Rk.x*dI.z;   // xx xz
              HessianMatrix.element[n][n+2]+=temp1*Rk.y*dI.y;   // xx yy
              HessianMatrix.element[n][n+3]+=temp1*Rk.z*dI.z;   // xx zz

              temp1=f2_I*Rk.x*dI.z;
              HessianMatrix.element[n+1][n+1]+=temp1*Rk.x*dI.z; // xy xz
              HessianMatrix.element[n+1][n+2]+=temp1*Rk.y*dI.y; // xy yy
              HessianMatrix.element[n+1][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.y*dI.y;
              HessianMatrix.element[n+2][n+2]+=temp1*Rk.y*dI.y; // xy yy
              HessianMatrix.element[n+2][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]+=temp1*Rk.z*dI.z; // xy zz

              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;     // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*Rk.x*dI.z;   // xx xz

              HessianMatrix.element[n+1][n+3]+=f1_I*Rk.x*dI.z; // xz zz

              HessianMatrix.element[n+2][n+2]+=f1_I*Rk.y*dI.y; // yy yy

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z; // zz zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n  ]+=f1_I*(Theta.ax*dI.x*Rk.x+Rk.x*dI.x*Theta.ax); // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*(Theta.ax*dI.y*Rk.x+Rk.x*dI.x*Theta.ay); // xx xy
              HessianMatrix.element[n][n+2]+=f1_I*(Theta.ax*dI.y*Rk.y+Rk.x*dI.x*Theta.by); // xx yy
              HessianMatrix.element[n][n+3]+=f1_I*(Theta.ax*dI.z*Rk.z+Rk.x*dI.x*Theta.cz); // xx zz

              HessianMatrix.element[n+1][n+1]+=f1_I*(Theta.ay*Rk.x*dI.y+Rk.x*dI.y*Theta.ay); // xy xy
              HessianMatrix.element[n+1][n+2]+=f1_I*(Theta.ay*Rk.y*dI.y+Rk.x*dI.y*Theta.by); // xy yy
              HessianMatrix.element[n+1][n+3]+=f1_I*(Theta.ay*Rk.z*dI.z+Rk.x*dI.y*Theta.cz); // xy zz

              HessianMatrix.element[n+2][n+2]+=f1_I*(Theta.by*Rk.y*dI.y+Rk.y*dI.y*Theta.by); // yy yy
              HessianMatrix.element[n+2][n+3]+=f1_I*(Theta.by*Rk.z*dI.z+Rk.y*dI.y*Theta.cz); // yy zz

              HessianMatrix.element[n+3][n+3]+=f1_I*(Theta.cz*Rk.z*dI.z+Rk.z*dI.z*Theta.cz); // zz zz


              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*Rk.x*dI.x;
              HessianMatrix.element[n][n  ]+=temp1*Rk.x*dI.x;   // xx xx
              HessianMatrix.element[n][n+1]+=temp1*Rk.x*dI.y;   // xx xy
              HessianMatrix.element[n][n+2]+=temp1*Rk.y*dI.y;   // xx yy
              HessianMatrix.element[n][n+3]+=temp1*Rk.z*dI.z;   // xx zz

              temp1=f2_I*Rk.x*dI.y;
              HessianMatrix.element[n+1][n+1]+=temp1*Rk.x*dI.y; // xy xy
              HessianMatrix.element[n+1][n+2]+=temp1*Rk.y*dI.y; // xy yy
              HessianMatrix.element[n+1][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.y*dI.y;
              HessianMatrix.element[n+2][n+2]+=temp1*Rk.y*dI.y; // xy yy
              HessianMatrix.element[n+2][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]+=temp1*Rk.z*dI.z; // xy zz

              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;     // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*Rk.x*dI.y;   // xx xy

              HessianMatrix.element[n+1][n+2]+=f1_I*Rk.x*dI.y; // xy yy

              HessianMatrix.element[n+2][n+2]+=f1_I*Rk.y*dI.y; // yy yy

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z; // zz zz
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

// Hessian: Strain - Strain (part III)
// ===================================
static inline void HessianAtomicCorrectionStrainStrainJ(REAL_MATRIX HessianMatrix,REAL f2_IJ,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB,VECTOR Rk)
{
  int n;
  REAL temp1;
  VECTOR dI,dJ;

  n=NumberOfCoordinatesMinimizationVariables;

  dI.x=posA.x-comA.x;
  dI.y=posA.y-comA.y;
  dI.z=posA.z-comA.z;

  dJ.x=posB.x-comB.x;
  dJ.y=posB.y-comB.y;
  dJ.z=posB.z-comB.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      temp1=f2_IJ*dI.x*Rk.x; // xx
      HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
      HessianMatrix.element[n][n]-=2.0*temp1*dJ.y*Rk.y;                   // xx yy
      HessianMatrix.element[n][n]-=2.0*temp1*dJ.z*Rk.z;                   // xx zz
      temp1=f2_IJ*dI.y*Rk.y; // yy
      HessianMatrix.element[n][n]-=temp1*dJ.y*Rk.y;                 // yy yy
      HessianMatrix.element[n][n]-=2.0*temp1*dJ.z*Rk.z;                 // yy zz
      temp1=f2_IJ*dI.z*Rk.z; // zz
      HessianMatrix.element[n][n]-=temp1*dJ.z*Rk.z;                 // zz zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          temp1=f2_IJ*dI.x*Rk.x; // xx
          HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
          HessianMatrix.element[n][n+1]-=temp1*dJ.y*Rk.y;                   // xx yy
          HessianMatrix.element[n][n+2]-=temp1*dJ.z*Rk.z;                   // xx zz
          temp1=f2_IJ*dI.y*Rk.y; // yy
          HessianMatrix.element[n+1][n+1]-=temp1*dJ.y*Rk.y;                 // yy yy
          HessianMatrix.element[n+1][n+2]-=temp1*dJ.z*Rk.z;                 // yy zz
          temp1=f2_IJ*dI.z*Rk.z; // zz
          HessianMatrix.element[n+2][n+2]-=temp1*dJ.z*Rk.z;                 // zz zz
          break;
        case REGULAR:
          // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          temp1=f2_IJ*dI.x*Rk.x; // xx
          HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
          HessianMatrix.element[n][n+1]-=0.5*temp1*(dJ.y*Rk.x+dJ.x*Rk.y);   // xx xy
          HessianMatrix.element[n][n+2]-=0.5*temp1*(dJ.z*Rk.x+dJ.x*Rk.z);   // xx xz
          HessianMatrix.element[n][n+3]-=temp1*dJ.y*Rk.y;                   // xx yy
          HessianMatrix.element[n][n+4]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z);   // xx yz
          HessianMatrix.element[n][n+5]-=temp1*dJ.z*Rk.z;                   // xx zz

          temp1=0.5*f2_IJ*(dI.x*Rk.y+dI.y*Rk.x); // xy
          HessianMatrix.element[n+1][n+1]-=0.5*temp1*(dJ.y*Rk.x+dJ.x*Rk.y); // xy xy
          HessianMatrix.element[n+1][n+2]-=0.5*temp1*(dJ.z*Rk.x+dJ.x*Rk.z); // xy xz
          HessianMatrix.element[n+1][n+3]-=temp1*dJ.y*Rk.y;                 // xy yy
          HessianMatrix.element[n+1][n+4]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // xy yz
          HessianMatrix.element[n+1][n+5]-=temp1*dJ.z*Rk.z;                 // xy zz

          temp1=0.5*f2_IJ*(dI.x*Rk.z+dI.z*Rk.x); // xz
          HessianMatrix.element[n+2][n+2]-=0.5*temp1*(dJ.z*Rk.x+dJ.x*Rk.z); // xz xz
          HessianMatrix.element[n+2][n+3]-=temp1*dJ.y*Rk.y;                 // xz yy
          HessianMatrix.element[n+2][n+4]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // xz yz
          HessianMatrix.element[n+2][n+5]-=temp1*dJ.z*Rk.z;                 // xz zz

          temp1=f2_IJ*dI.y*Rk.y; // yy
          HessianMatrix.element[n+3][n+3]-=temp1*dJ.y*Rk.y;                 // yy yy
          HessianMatrix.element[n+3][n+4]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // yy yz
          HessianMatrix.element[n+3][n+5]-=temp1*dJ.z*Rk.z;                 // yy zz

          temp1=0.5*f2_IJ*(dI.y*Rk.z+dI.z*Rk.y); // yz
          HessianMatrix.element[n+4][n+4]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // yz yz
          HessianMatrix.element[n+4][n+5]-=temp1*dJ.z*Rk.z;                 // yz zz

          temp1=f2_IJ*dI.z*Rk.z; // zz
          HessianMatrix.element[n+5][n+5]-=temp1*dJ.z*Rk.z;                 // zz zz
          break;
        case REGULAR_UPPER_TRIANGLE:
          // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================

          temp1=f2_IJ*Rk.x*dI.x; 
          HessianMatrix.element[n][n  ]-=temp1*Rk.x*dJ.x;   // xx xx
          HessianMatrix.element[n][n+1]-=temp1*Rk.x*dJ.y;   // xx xy
          HessianMatrix.element[n][n+2]-=temp1*Rk.x*dJ.z;   // xx xz
          HessianMatrix.element[n][n+3]-=temp1*Rk.y*dJ.y;   // xx yy
          HessianMatrix.element[n][n+4]-=temp1*Rk.y*dJ.z;   // xx yz
          HessianMatrix.element[n][n+5]-=temp1*Rk.z*dJ.z;   // xx zz

          temp1=f2_IJ*Rk.x*dI.y;
          HessianMatrix.element[n+1][n+1]-=temp1*Rk.x*dJ.y; // xy xy
          HessianMatrix.element[n+1][n+2]-=temp1*Rk.x*dJ.z; // xy xz
          HessianMatrix.element[n+1][n+3]-=temp1*Rk.y*dJ.y; // xy yy
          HessianMatrix.element[n+1][n+4]-=temp1*Rk.y*dJ.z; // xy yz
          HessianMatrix.element[n+1][n+5]-=temp1*Rk.z*dJ.z; // xy zz

          temp1=f2_IJ*Rk.x*dI.z;
          HessianMatrix.element[n+2][n+2]-=temp1*Rk.x*dJ.z; // xz xz
          HessianMatrix.element[n+2][n+3]-=temp1*Rk.y*dJ.y; // xz yy
          HessianMatrix.element[n+2][n+4]-=temp1*Rk.y*dJ.z; // xz yz
          HessianMatrix.element[n+2][n+5]-=temp1*Rk.z*dJ.z; // xz zz

          temp1=f2_IJ*Rk.y*dI.y;
          HessianMatrix.element[n+3][n+3]-=temp1*Rk.y*dJ.y; // yy yy
          HessianMatrix.element[n+3][n+4]-=temp1*Rk.y*dJ.z; // yy yz
          HessianMatrix.element[n+3][n+5]-=temp1*Rk.z*dJ.z; // yy zz

          temp1=f2_IJ*Rk.y*dI.z;
          HessianMatrix.element[n+4][n+4]-=temp1*Rk.y*dJ.z; // yz yz
          HessianMatrix.element[n+4][n+5]-=temp1*Rk.z*dJ.z; // yz zz

          temp1=f2_IJ*Rk.z*dI.z;
          HessianMatrix.element[n+5][n+5]-=temp1*Rk.z*dJ.z; // zz zz
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_IJ*dI.x*Rk.x; // xx
              HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]-=temp1*dJ.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+2]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z);   // xx yz
              HessianMatrix.element[n][n+3]-=temp1*dJ.z*Rk.z;                   // xx zz

              temp1=f2_IJ*dI.y*Rk.y; // yy
              HessianMatrix.element[n+1][n+1]-=temp1*dJ.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+1][n+2]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // yy yz
              HessianMatrix.element[n+1][n+3]-=temp1*dJ.z*Rk.z;                 // yy zz

              temp1=0.5*f2_IJ*(dI.y*Rk.z+dI.z*Rk.y); // yz
              HessianMatrix.element[n+2][n+2]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // yz yz
              HessianMatrix.element[n+2][n+3]-=temp1*dJ.z*Rk.z;                 // yz zz

              temp1=f2_IJ*dI.z*Rk.z; // zz
              HessianMatrix.element[n+3][n+3]-=temp1*dJ.z*Rk.z;                 // zz zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_IJ*dI.x*Rk.x; // xx
              HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]-=0.5*temp1*(dJ.z*Rk.x+dJ.x*Rk.z);   // xx xz
              HessianMatrix.element[n][n+2]-=temp1*dJ.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+3]-=temp1*dJ.z*Rk.z;                   // xx zz

              temp1=0.5*f2_IJ*(dI.x*Rk.z+dI.z*Rk.x); // xz
              HessianMatrix.element[n+1][n+1]-=0.5*temp1*(dJ.z*Rk.x+dJ.x*Rk.z); // xz xz
              HessianMatrix.element[n+1][n+2]-=temp1*dJ.y*Rk.y;                 // xz yy
              HessianMatrix.element[n+1][n+3]-=temp1*dJ.z*Rk.z;                 // xz zz

              temp1=f2_IJ*dI.y*Rk.y; // yy
              HessianMatrix.element[n+2][n+2]-=temp1*dJ.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+2][n+3]-=temp1*dJ.z*Rk.z;                 // yy zz

              temp1=f2_IJ*dI.z*Rk.z; // zz
              HessianMatrix.element[n+3][n+3]-=temp1*dJ.z*Rk.z;                 // zz zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_IJ*dI.x*Rk.x; // xx
              HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]-=0.5*temp1*(dJ.y*Rk.x+dJ.x*Rk.y);   // xx xy
              HessianMatrix.element[n][n+2]-=temp1*dJ.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+3]-=temp1*dJ.z*Rk.z;                   // xx zz

              temp1=0.5*f2_IJ*(dI.x*Rk.y+dI.y*Rk.x); // xy
              HessianMatrix.element[n+1][n+1]-=0.5*temp1*(dJ.y*Rk.x+dJ.x*Rk.y); // xy xy
              HessianMatrix.element[n+1][n+2]-=temp1*dJ.y*Rk.y;                 // xy yy
              HessianMatrix.element[n+1][n+3]-=temp1*dJ.z*Rk.z;                 // xy zz

              temp1=f2_IJ*dI.y*Rk.y; // yy
              HessianMatrix.element[n+2][n+2]-=temp1*dJ.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+2][n+3]-=temp1*dJ.z*Rk.z;                 // yy zz

              temp1=f2_IJ*dI.z*Rk.z; // zz
              HessianMatrix.element[n+3][n+3]-=temp1*dJ.z*Rk.z;                 // zz zz
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================

              temp1=f2_IJ*Rk.x*dI.x; 
              HessianMatrix.element[n][n  ]-=temp1*Rk.x*dJ.x;   // xx xx
              HessianMatrix.element[n][n+1]-=temp1*Rk.y*dJ.y;   // xx yy
              HessianMatrix.element[n][n+2]-=temp1*Rk.y*dJ.z;   // xx yz
              HessianMatrix.element[n][n+3]-=temp1*Rk.z*dJ.z;   // xx zz

              temp1=f2_IJ*Rk.y*dI.y;
              HessianMatrix.element[n+1][n+1]-=temp1*Rk.y*dJ.y; // yy yy
              HessianMatrix.element[n+1][n+2]-=temp1*Rk.y*dJ.z; // yy yz
              HessianMatrix.element[n+1][n+3]-=temp1*Rk.z*dJ.z; // yy zz

              temp1=f2_IJ*Rk.y*dI.z;
              HessianMatrix.element[n+2][n+2]-=temp1*Rk.y*dJ.z; // yz yz
              HessianMatrix.element[n+2][n+3]-=temp1*Rk.z*dJ.z; // yz zz

              temp1=f2_IJ*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]-=temp1*Rk.z*dJ.z; // zz zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================

              temp1=f2_IJ*Rk.x*dI.x;
              HessianMatrix.element[n][n  ]-=temp1*Rk.x*dJ.x;   // xx xx
              HessianMatrix.element[n][n+1]-=temp1*Rk.x*dJ.z;   // xx xz
              HessianMatrix.element[n][n+2]-=temp1*Rk.y*dJ.y;   // xx yy
              HessianMatrix.element[n][n+3]-=temp1*Rk.z*dJ.z;   // xx zz

              temp1=f2_IJ*Rk.x*dI.z;
              HessianMatrix.element[n+1][n+1]-=temp1*Rk.x*dJ.z; // xz xz
              HessianMatrix.element[n+1][n+2]-=temp1*Rk.y*dJ.y; // xz yy
              HessianMatrix.element[n+1][n+3]-=temp1*Rk.z*dJ.z; // xz zz

              temp1=f2_IJ*Rk.y*dI.y;
              HessianMatrix.element[n+2][n+2]-=temp1*Rk.y*dJ.y; // yy yy
              HessianMatrix.element[n+2][n+3]-=temp1*Rk.z*dJ.z; // yy zz

              temp1=f2_IJ*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]-=temp1*Rk.z*dJ.z; // zz zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================

              temp1=f2_IJ*Rk.x*dI.x;
              HessianMatrix.element[n][n  ]-=temp1*Rk.x*dJ.x;   // xx xx
              HessianMatrix.element[n][n+1]-=temp1*Rk.x*dJ.y;   // xx xy
              HessianMatrix.element[n][n+2]-=temp1*Rk.y*dJ.y;   // xx yy
              HessianMatrix.element[n][n+3]-=temp1*Rk.z*dJ.z;   // xx zz

              temp1=f2_IJ*Rk.x*dI.y;
              HessianMatrix.element[n+1][n+1]-=temp1*Rk.x*dJ.y; // xy xy
              HessianMatrix.element[n+1][n+2]-=temp1*Rk.y*dJ.y; // xy yy
              HessianMatrix.element[n+1][n+3]-=temp1*Rk.z*dJ.z; // xy zz

              temp1=f2_IJ*Rk.y*dI.y;
              HessianMatrix.element[n+2][n+2]-=temp1*Rk.y*dJ.y; // yy yy
              HessianMatrix.element[n+2][n+3]-=temp1*Rk.z*dJ.z; // yy zz

              temp1=f2_IJ*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]-=temp1*Rk.z*dJ.z; // zz zz
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


static inline void GradientStrain(REAL fac,REAL *Gradient,REAL_MATRIX3x3 Theta)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]-=fac*(Theta.ax+Theta.by+Theta.cz); // xx+yy+zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n  ]-=fac*Theta.ax; // xx
          Gradient[n+1]-=fac*Theta.by; // yy
          Gradient[n+2]-=fac*Theta.cz; // zz
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n  ]-=fac*Theta.ax; // xx
          Gradient[n+1]-=fac*Theta.bx; // xy
          Gradient[n+2]-=fac*Theta.cx; // xz
          Gradient[n+3]-=fac*Theta.by; // yy
          Gradient[n+4]-=fac*Theta.cy; // yz
          Gradient[n+5]-=fac*Theta.cz; // zz
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]-=fac*Theta.ax; // xx
              Gradient[n+1]-=fac*Theta.by; // yy
              Gradient[n+2]-=fac*Theta.cy; // yz
              Gradient[n+3]-=fac*Theta.cz; // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]-=fac*Theta.ax; // xx
              Gradient[n+1]-=fac*Theta.cx; // xz
              Gradient[n+2]-=fac*Theta.by; // yy
              Gradient[n+3]-=fac*Theta.cz; // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]-=fac*Theta.ax; // xx
              Gradient[n+1]-=fac*Theta.bx; // xy
              Gradient[n+2]-=fac*Theta.by; // yy
              Gradient[n+3]-=fac*Theta.cz; // zz
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

static inline void GradientStrainILocal(REAL *Gradient,REAL fac,VECTOR Rk,VECTOR posA,VECTOR comA)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]-=fac*Rk.x*(posA.x-comA.x)+fac*Rk.y*(posA.y-comA.y)+fac*Rk.z*(posA.z-comA.z); // xx+yy+zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x); // xx
          Gradient[n+1]-=fac*Rk.y*(posA.y-comA.y); // yy
          Gradient[n+2]-=fac*Rk.z*(posA.z-comA.z); // zz
          break;
        case REGULAR:
          Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x);                            // xx
          Gradient[n+1]-=0.5*fac*(Rk.x*(posA.y-comA.y)+Rk.y*(posA.x-comA.x)); // xy
          Gradient[n+2]-=0.5*fac*(Rk.x*(posA.z-comA.z)+Rk.z*(posA.x-comA.x)); // xz
          Gradient[n+3]-=fac*Rk.y*(posA.y-comA.y);                            // yy
          Gradient[n+4]-=0.5*fac*(Rk.y*(posA.z-comA.z)+Rk.z*(posA.y-comA.y)); // yz
          Gradient[n+5]-=fac*Rk.z*(posA.z-comA.z);                            // zz
          break;
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x); // xx
          Gradient[n+1]-=fac*Rk.x*(posA.y-comA.y); // xy
          Gradient[n+2]-=fac*Rk.x*(posA.z-comA.z); // xz
          Gradient[n+3]-=fac*Rk.y*(posA.y-comA.y); // yy
          Gradient[n+4]-=fac*Rk.y*(posA.z-comA.z); // yz
          Gradient[n+5]-=fac*Rk.z*(posA.z-comA.z); // zz
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x);                            // xx
              Gradient[n+1]-=fac*Rk.y*(posA.y-comA.y);                            // yy
              Gradient[n+2]-=0.5*fac*(Rk.y*(posA.z-comA.z)+Rk.z*(posA.y-comA.y)); // yz
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z);                            // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x);                            // xx
              Gradient[n+1]-=0.5*fac*(Rk.x*(posA.z-comA.z)+Rk.z*(posA.x-comA.x)); // xz
              Gradient[n+2]-=fac*Rk.y*(posA.y-comA.y);                            // yy
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z);                            // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x);                            // xx
              Gradient[n+1]-=0.5*fac*(Rk.x*(posA.y-comA.y)+Rk.y*(posA.x-comA.x)); // xy
              Gradient[n+2]-=fac*Rk.y*(posA.y-comA.y);                            // yy
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z);                            // zz
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x); // xx
              Gradient[n+1]-=fac*Rk.y*(posA.y-comA.y); // yy
              Gradient[n+2]-=fac*Rk.y*(posA.z-comA.z); // yz
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z); // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x); // xx
              Gradient[n+1]-=fac*Rk.x*(posA.z-comA.z); // xz
              Gradient[n+2]-=fac*Rk.y*(posA.y-comA.y); // yy
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z); // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x); // xx
              Gradient[n+1]-=fac*Rk.x*(posA.y-comA.y); // xy
              Gradient[n+2]-=fac*Rk.y*(posA.y-comA.y); // yy
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z); // zz
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

static inline void HessianAtomicPositionStrainExcluded(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,REAL f1,REAL f2,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;
      if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;
      if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;

      if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;
      if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;
      if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.y*dr.y*dr.x;
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.z*dr.z*dr.x;

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n]+=f2*dr.x*dr.x*dr.y;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.z*dr.z*dr.y;

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n]+=f2*dr.x*dr.x*dr.z;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.y*dr.y*dr.z;
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.y*dr.y*dr.x;
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.z*dr.z*dr.x;

          if(index_j.y>=0) HessianMatrix.element[index_j.y][n]-=f2*dr.x*dr.x*dr.y;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.z*dr.z*dr.y;

          if(index_j.z>=0) HessianMatrix.element[index_j.z][n]-=f2*dr.x*dr.x*dr.z;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.y*dr.y*dr.z;
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;     // xy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.y*dr.y*dr.x;             // yy x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+4]+=f2*dr.y*dr.z*dr.x;             // yz x
          if(index_i.x>=0) HessianMatrix.element[index_i.x][n+5]+=f2*dr.z*dr.z*dr.x;             // zz x

          if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2*dr.x*dr.x*dr.y;             // xx y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;     // xy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.x*dr.z*dr.y;             // xz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+4]+=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
          if(index_i.y>=0) HessianMatrix.element[index_i.y][n+5]+=f2*dr.z*dr.z*dr.y;             // zz y

          if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2*dr.x*dr.x*dr.z;             // xx z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.y*dr.z;             // xy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.y*dr.y*dr.z;             // yy z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+4]+=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
          if(index_i.z>=0) HessianMatrix.element[index_i.z][n+5]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

          if(index_j.x>=0) HessianMatrix.element[index_j.x][n  ]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;     // xy x
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=f2*dr.y*dr.y*dr.x;             // yy x
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+4]-=f2*dr.y*dr.z*dr.x;             // yz x
          if(index_j.x>=0) HessianMatrix.element[index_j.x][n+5]-=f2*dr.z*dr.z*dr.x;             // zz x

          if(index_j.y>=0) HessianMatrix.element[index_j.y][n  ]-=f2*dr.x*dr.x*dr.y;             // xx y
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;     // xy y
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.x*dr.z*dr.y;             // xz y
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+4]-=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
          if(index_j.y>=0) HessianMatrix.element[index_j.y][n+5]-=f2*dr.z*dr.z*dr.y;             // zz y

          if(index_j.z>=0) HessianMatrix.element[index_j.z][n  ]-=f2*dr.x*dr.x*dr.z;             // xx z
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.y*dr.z;             // xy z
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.y*dr.y*dr.z;             // yy z
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+4]-=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
          if(index_j.z>=0) HessianMatrix.element[index_j.z][n+5]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
          break;                                   
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.z*dr.x;             // yz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2*dr.x*dr.x*dr.y;             // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2*dr.x*dr.x*dr.z;             // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n  ]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.y*dr.z*dr.x;             // yz x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n  ]-=f2*dr.x*dr.x*dr.y;             // xx y
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n  ]-=f2*dr.x*dr.x*dr.z;             // xx z
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2*dr.x*dr.x*dr.y;             // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.z*dr.y;             // xz y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2*dr.x*dr.x*dr.z;             // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n  ]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n  ]-=f2*dr.x*dr.x*dr.y;             // xx y
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.z*dr.y;             // xz y
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n  ]-=f2*dr.x*dr.x*dr.z;             // xx z
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n  ]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;     // xy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+2]+=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_i.x>=0) HessianMatrix.element[index_i.x][n+3]+=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_i.y>=0) HessianMatrix.element[index_i.y][n  ]+=f2*dr.x*dr.x*dr.y;             // xx y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;     // xy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_i.y>=0) HessianMatrix.element[index_i.y][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_i.z>=0) HessianMatrix.element[index_i.z][n  ]+=f2*dr.x*dr.x*dr.z;             // xx z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+1]+=f2*dr.x*dr.y*dr.z;             // xy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+2]+=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_i.z>=0) HessianMatrix.element[index_i.z][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z

              if(index_j.x>=0) HessianMatrix.element[index_j.x][n  ]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;     // xy x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+2]-=f2*dr.y*dr.y*dr.x;             // yy x
              if(index_j.x>=0) HessianMatrix.element[index_j.x][n+3]-=f2*dr.z*dr.z*dr.x;             // zz x

              if(index_j.y>=0) HessianMatrix.element[index_j.y][n  ]-=f2*dr.x*dr.x*dr.y;             // xx y
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;     // xy y
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+2]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
              if(index_j.y>=0) HessianMatrix.element[index_j.y][n+3]-=f2*dr.z*dr.z*dr.y;             // zz y

              if(index_j.z>=0) HessianMatrix.element[index_j.z][n  ]-=f2*dr.x*dr.x*dr.z;             // xx z
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+1]-=f2*dr.x*dr.y*dr.z;             // xy z
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+2]-=f2*dr.y*dr.y*dr.z;             // yy z
              if(index_j.z>=0) HessianMatrix.element[index_j.z][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
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

static inline void HessianAtomicStrainStrainLocalExcluded(REAL_MATRIX HessianMatrix,REAL f1,REAL f2,VECTOR dr)
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
          HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
          HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
          HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz
          HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
          HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz
          HessianMatrix.element[n+2][n+2]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
          break;
        case REGULAR:
          HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
          HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;    // xxxy
          HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;    // xxxz
          HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
          HessianMatrix.element[n][n+4]+=f2*dr.x*dr.x*dr.y*dr.z;                 // xxyz
          HessianMatrix.element[n][n+5]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

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
          HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
          HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;       // xxxy
          HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;       // xxxz
          HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
          HessianMatrix.element[n][n+4]+=f2*dr.x*dr.x*dr.y*dr.z;                    // xxyz
          HessianMatrix.element[n][n+5]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz


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
              HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.z;                 // xxyz
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;                  // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.z*dr.y*dr.z+0.5*f1*(dr.y*dr.y+dr.z*dr.z);  // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.z*dr.z*dr.z+f1*dr.y*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;    // xxxz
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.z*dr.x*dr.z+0.5*f1*(dr.x*dr.x+dr.z*dr.z);  // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.z*dr.y*dr.y;                               // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.z*dr.z*dr.z+f1*dr.x*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;    // xxxy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

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
              HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.z;                    // xxyz
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;     // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.z*dr.y*dr.z+f1*dr.z*dr.z;     // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.z*dr.z*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;       // xxxz
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.z*dr.x*dr.z+f1*dr.z*dr.z;     // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.z*dr.y*dr.y;                  // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.z*dr.z*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;       // xxxy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz


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

/*********************************************************************************************************
 * Name       | CalculateEwaldFourierDerivatives                                                         *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the Fourier part of energy, stress, and energy/strain  derivatives.             *
 * Parameters |                                                                                          *
 *********************************************************************************************************/

int CalculateEwaldFourierDerivatives(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainFirstDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,j,l,m,f1,f2;
  int type,typeA,typeB,Type;
  REAL Rksq,fac;
  REAL exp_term,energy_term;
  VECTOR pos,Rk;
  int nvec;
  REAL ChargeA,ChargeB,charge;
  int A,B,NumberOfExcludedPairs,TypeA,TypeB;
  REAL r,rr;
  POINT posA,posB,comA,comB;
  VECTOR dr;
  int index;
  INT_VECTOR3 index_i,index_j,index_i2,index_j2;
  int index_ic,index_jc;
  REAL_MATRIX3x3 Hessian,Stress;
  REAL USum,Uself_sum,DF,DDF;
  COMPLEX Cksum,CksumFramework;
  REAL InverseLamdaSquared;
  REAL Usum_framework;
  int TypeMolA,TypeMolB;
  int I,J,ia,ja,ig,jg;
  int grpA,grpB,RigidI,RigidJ;
  int index1,index2,index1_rigid,index2_rigid;
  REAL temp1,temp2,temp3;
  REAL_MATRIX3x3 S,Theta;
  REAL f,f1_I,f2_I,f2_IJ;
  VECTOR dot_product_i,dot_product_j;
  REAL dot_product_AX,dot_product_BY,dot_product_CZ;
  REAL dot_product_AY,dot_product_AZ,dot_product_BZ;

  int nr_of_coulombic_sites;
  int ii,jj,kk,index_k;
  int kmax_x,kmax_y,kmax_z;
  REAL *kfactor,alpha,temp;
  VECTOR *kvecs;
  REAL factor,recip_cutoff,ksqr;
  int considered_charged;

  if(ChargeMethod==NONE) return 0;

  alpha=Alpha[CurrentSystem];
  kmax_x=kvec[CurrentSystem].x;
  kmax_y=kvec[CurrentSystem].y;
  kmax_z=kvec[CurrentSystem].z;
  kvecs=KVectors[CurrentSystem];
  kfactor=KFactor[CurrentSystem];
  recip_cutoff=ReciprocalCutOffSquared[CurrentSystem];

  fac=0.0;
  USum=0.0;
  Uself_sum=0.0;
  Usum_framework=0.0;

  Stress.ax=Stress.bx=Stress.cx=0.0;
  Stress.ay=Stress.by=Stress.cy=0.0;
  Stress.az=Stress.bz=Stress.cz=0.0;

  nr_of_coulombic_sites=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      type=Framework[CurrentSystem].Atoms[f1][i].Type;
      charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
      considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
      if(considered_charged&&(!(Framework[CurrentSystem].Atoms[f1][i].Fixed.x&&Framework[CurrentSystem].Atoms[f1][i].Fixed.y&&Framework[CurrentSystem].Atoms[f1][i].Fixed.z)))
      {
        Charge[nr_of_coulombic_sites]=charge;
        Uself_sum+=SQR(Charge[nr_of_coulombic_sites]);
        Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[f1][i].Position);
        Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
        nr_of_coulombic_sites++;
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        if(!(Adsorbates[CurrentSystem][m].Groups[l].FixedCenterOfMass.x&&Adsorbates[CurrentSystem][m].Groups[l].FixedCenterOfMass.y&&
             Adsorbates[CurrentSystem][m].Groups[l].FixedCenterOfMass.z&&Adsorbates[CurrentSystem][m].Groups[l].FixedOrientation.x&&
             Adsorbates[CurrentSystem][m].Groups[l].FixedOrientation.y&&Adsorbates[CurrentSystem][m].Groups[l].FixedOrientation.z))
        {
          for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
          {
            A=Components[Type].Groups[l].Atoms[i];
            type=Adsorbates[CurrentSystem][m].Atoms[A].Type;
            charge=Adsorbates[CurrentSystem][m].Atoms[A].Charge;
            considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
            if(considered_charged)
            {
              Charge[nr_of_coulombic_sites]=charge;
              Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][m].Atoms[A].Position);
              Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
              nr_of_coulombic_sites++;
            }
          }
        }
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          type=Adsorbates[CurrentSystem][m].Atoms[A].Type;
          charge=Adsorbates[CurrentSystem][m].Atoms[A].Charge;
          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
          if(considered_charged&&(!(Adsorbates[CurrentSystem][m].Atoms[A].Fixed.x&&Adsorbates[CurrentSystem][m].Atoms[A].Fixed.y&&Adsorbates[CurrentSystem][m].Atoms[A].Fixed.z)))
          {
            Charge[nr_of_coulombic_sites]=charge;
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Adsorbates[CurrentSystem][m].Atoms[A].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        if(!(Cations[CurrentSystem][m].Groups[l].FixedCenterOfMass.x&&Cations[CurrentSystem][m].Groups[l].FixedCenterOfMass.y&&
             Cations[CurrentSystem][m].Groups[l].FixedCenterOfMass.z&&Cations[CurrentSystem][m].Groups[l].FixedOrientation.x&&
             Cations[CurrentSystem][m].Groups[l].FixedOrientation.y&&Cations[CurrentSystem][m].Groups[l].FixedOrientation.z))
        {
          for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
          {
            A=Components[Type].Groups[l].Atoms[i];
            type=Cations[CurrentSystem][m].Atoms[A].Type;
            charge=Cations[CurrentSystem][m].Atoms[A].Charge;
            considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
            if(considered_charged)
            {
              Charge[nr_of_coulombic_sites]=charge;
              Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][m].Atoms[A].Position);
              Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
              nr_of_coulombic_sites++;
            }
          }
        }
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          type=Cations[CurrentSystem][m].Atoms[A].Type;
          charge=Cations[CurrentSystem][m].Atoms[A].Charge;
          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
          if(considered_charged&&(!(Cations[CurrentSystem][m].Atoms[A].Fixed.x&&Cations[CurrentSystem][m].Atoms[A].Fixed.y&&Cations[CurrentSystem][m].Atoms[A].Fixed.z)))
          {
            Charge[nr_of_coulombic_sites]=charge;
            Positions[nr_of_coulombic_sites]=ConvertFromXYZtoABC(Cations[CurrentSystem][m].Atoms[A].Position);
            Positions[nr_of_coulombic_sites].x*=TWO_PI; Positions[nr_of_coulombic_sites].y*=TWO_PI; Positions[nr_of_coulombic_sites].z*=TWO_PI;
            nr_of_coulombic_sites++;
          }
        }
      }
    }
  }

  // calculate kx,ky,kz=-1,0,1 explicitly
  for(i=0;i<nr_of_coulombic_sites;i++)
  {
    Eikx[i].re=1.0; Eikx[i].im=0.0;
    Eiky[i].re=1.0; Eiky[i].im=0.0;
    Eikz[i].re=1.0; Eikz[i].im=0.0;

    pos=Positions[i];

    index=MaxNumberOfCoulombicSites+i;
    Eikx[index].re=cos(pos.x); Eikx[index].im=sin(pos.x);
    Eiky[index].re=cos(pos.y); Eiky[index].im=sin(pos.y);
    Eikz[index].re=cos(pos.z); Eikz[index].im=sin(pos.z);

    index=-MaxNumberOfCoulombicSites+i;
    Eikx[index].re=cos(pos.x); Eikx[index].im=-sin(pos.x);
    Eiky[index].re=cos(pos.y); Eiky[index].im=-sin(pos.y);
    Eikz[index].re=cos(pos.z); Eikz[index].im=-sin(pos.z);
  }

  // calculate remaining kx,jy,kz by a recurrence formula (to avoid using 'cos' and 'sin' explicitly)
  // in the x-direction symmetry is used (-kx=kx) and only positive wavevectors are used
  for(j=2;j<=kmax_x;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikx[j*MaxNumberOfCoulombicSites+i].re=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].re-
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].im;
      Eikx[j*MaxNumberOfCoulombicSites+i].im=Eikx[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikx[MaxNumberOfCoulombicSites+i].re+
                                             Eikx[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikx[MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_y;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eiky[j*MaxNumberOfCoulombicSites+i].re=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].re-
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[j*MaxNumberOfCoulombicSites+i].im=Eiky[(j-1)*MaxNumberOfCoulombicSites+i].im*Eiky[MaxNumberOfCoulombicSites+i].re+
                                             Eiky[(j-1)*MaxNumberOfCoulombicSites+i].re*Eiky[MaxNumberOfCoulombicSites+i].im;
      Eiky[-j*MaxNumberOfCoulombicSites+i].re=Eiky[j*MaxNumberOfCoulombicSites+i].re;
      Eiky[-j*MaxNumberOfCoulombicSites+i].im=-Eiky[j*MaxNumberOfCoulombicSites+i].im;
    }

  for(j=2;j<=kmax_z;j++)
    for(i=0;i<nr_of_coulombic_sites;i++)
    {
      Eikz[j*MaxNumberOfCoulombicSites+i].re=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].re-
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[j*MaxNumberOfCoulombicSites+i].im=Eikz[(j-1)*MaxNumberOfCoulombicSites+i].im*Eikz[MaxNumberOfCoulombicSites+i].re+
                                             Eikz[(j-1)*MaxNumberOfCoulombicSites+i].re*Eikz[MaxNumberOfCoulombicSites+i].im;
      Eikz[-j*MaxNumberOfCoulombicSites+i].re=Eikz[j*MaxNumberOfCoulombicSites+i].re;
      Eikz[-j*MaxNumberOfCoulombicSites+i].im=-Eikz[j*MaxNumberOfCoulombicSites+i].im;
    }

  // Main loop of the Fourier-routine
  // ================================

  nvec=0;
  for(ii=0;ii<=kmax_x;ii++)
  {
    for(jj=-kmax_y;jj<=kmax_y;jj++)
    {
      for(i=0;i<nr_of_coulombic_sites;i++)
      {
        // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
        // precompute exp(-ik.kx)*exp(-ik.ky) outside the 'kk' loop
        index_ic=ii*MaxNumberOfCoulombicSites+i;
        index_jc=jj*MaxNumberOfCoulombicSites+i;
        Eikr_xy[i].re=Eikx[index_ic].re*Eiky[index_jc].re-Eikx[index_ic].im*Eiky[index_jc].im;
        Eikr_xy[i].im=Eikx[index_ic].im*Eiky[index_jc].re+Eikx[index_ic].re*Eiky[index_jc].im;
      }

      for(kk=-kmax_z;kk<=kmax_z;kk++)
      {
        ksqr=SQR(ii)+SQR(jj)+SQR(kk);
        if((ksqr!=0)&&(ksqr<recip_cutoff)) // explicitly exclude |k|=0
        {
          Rk=kvecs[nvec];
          factor=kfactor[nvec];
          Rksq=SQR(Rk.x)+SQR(Rk.y)+SQR(Rk.z);

          Cksum.re=0.0;
          Cksum.im=0.0;

          CksumFramework.re=0.0;
          CksumFramework.im=0.0;

          index=0;
          if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
          {
            for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
            {
              for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
              {
                if(!Framework[CurrentSystem].Atoms[f1][i].Fixed.x)
                {
                  type=Framework[CurrentSystem].Atoms[f1][i].Type;
                  charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
                  considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                  if(considered_charged)
                  {
                    // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
                    index_k=kk*MaxNumberOfCoulombicSites+index;
                    temp=Charge[index];
                    Eikr[index].re=temp*(Eikr_xy[index].re*Eikz[index_k].re-Eikr_xy[index].im*Eikz[index_k].im);
                    Eikr[index].im=temp*(Eikr_xy[index].im*Eikz[index_k].re+Eikr_xy[index].re*Eikz[index_k].im);
                    Cksum.re+=Eikr[index].re;
                    Cksum.im+=Eikr[index].im;
                    index++;
                  }
                }
              }
            }
          }

          for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
          {
            TypeMolA=Adsorbates[CurrentSystem][I].Type;
            for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
            {
              if(Components[TypeMolA].Groups[ig].Rigid)
              {
                if(!((Adsorbates[CurrentSystem][I].Groups[ig].FixedCenterOfMass.x)&&(Adsorbates[CurrentSystem][I].Groups[ig].FixedOrientation.x)))
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    type=Adsorbates[CurrentSystem][I].Atoms[A].Type;
                    charge=Adsorbates[CurrentSystem][I].Atoms[A].Charge;
                    considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                    if(considered_charged)
                    {
                      // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
                      index_k=kk*MaxNumberOfCoulombicSites+index;
                      temp=Charge[index];
                      Eikr[index].re=temp*(Eikr_xy[index].re*Eikz[index_k].re-Eikr_xy[index].im*Eikz[index_k].im);
                      Eikr[index].im=temp*(Eikr_xy[index].im*Eikz[index_k].re+Eikr_xy[index].re*Eikz[index_k].im);
                      Cksum.re+=Eikr[index].re;
                      Cksum.im+=Eikr[index].im;
                      index++;
                    }
                  }
                }
              }
              else
              {
                for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                {
                  A=Components[TypeMolA].Groups[ig].Atoms[ia];
                  if(!Adsorbates[CurrentSystem][I].Atoms[A].Fixed.x)
                  {
                    type=Adsorbates[CurrentSystem][I].Atoms[A].Type;
                    charge=Adsorbates[CurrentSystem][I].Atoms[A].Charge;
                    considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                    if(considered_charged)
                    {
                      // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
                      index_k=kk*MaxNumberOfCoulombicSites+index;
                      temp=Charge[index];
                      Eikr[index].re=temp*(Eikr_xy[index].re*Eikz[index_k].re-Eikr_xy[index].im*Eikz[index_k].im);
                      Eikr[index].im=temp*(Eikr_xy[index].im*Eikz[index_k].re+Eikr_xy[index].re*Eikz[index_k].im);
                      Cksum.re+=Eikr[index].re;
                      Cksum.im+=Eikr[index].im;
                      index++;
                    }
                  }
                }
              }
            }
          }

          for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
          {
            TypeMolA=Cations[CurrentSystem][I].Type;
            for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
            {
              if(Components[TypeMolA].Groups[ig].Rigid)
              {
                if(!((Cations[CurrentSystem][I].Groups[ig].FixedCenterOfMass.x)&&(Cations[CurrentSystem][I].Groups[ig].FixedOrientation.x)))
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    type=Cations[CurrentSystem][I].Atoms[A].Type;
                    charge=Cations[CurrentSystem][I].Atoms[A].Charge;
                    considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                    if(considered_charged)
                    {
                      // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
                      index_k=kk*MaxNumberOfCoulombicSites+index;
                      temp=Charge[index];
                      Eikr[index].re=temp*(Eikr_xy[index].re*Eikz[index_k].re-Eikr_xy[index].im*Eikz[index_k].im);
                      Eikr[index].im=temp*(Eikr_xy[index].im*Eikz[index_k].re+Eikr_xy[index].re*Eikz[index_k].im);
                      Cksum.re+=Eikr[index].re;
                      Cksum.im+=Eikr[index].im;
                      index++;
                    }
                  }
                }
              }
              else
              {
                for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                {
                  A=Components[TypeMolA].Groups[ig].Atoms[ia];
                  if(!Cations[CurrentSystem][I].Atoms[A].Fixed.x)
                  {
                    type=Cations[CurrentSystem][I].Atoms[A].Type;
                    charge=Cations[CurrentSystem][I].Atoms[A].Charge;
                    considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                    if(considered_charged)
                    {
                      // exp(-ik.r)=exp(-ik.kx)*exp(-ik.ky)*exp(-ik.kz)
                      index_k=kk*MaxNumberOfCoulombicSites+index;
                      temp=Charge[index];
                      Eikr[index].re=temp*(Eikr_xy[index].re*Eikz[index_k].re-Eikr_xy[index].im*Eikz[index_k].im);
                      Eikr[index].im=temp*(Eikr_xy[index].im*Eikz[index_k].re+Eikr_xy[index].re*Eikz[index_k].im);
                      Cksum.re+=Eikr[index].re;
                      Cksum.im+=Eikr[index].im;
                      index++;
                    }
                  }
                }
              }
            }
          }

          Cksum.re+=StoreRigidChargeFramework[CurrentSystem][nvec].re;
          Cksum.im+=StoreRigidChargeFramework[CurrentSystem][nvec].im;

          Cksum.re+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].re;
          Cksum.im+=StoreRigidChargeAdsorbates[CurrentSystem][nvec].im;

          Cksum.re+=StoreRigidChargeCations[CurrentSystem][nvec].re;
          Cksum.im+=StoreRigidChargeCations[CurrentSystem][nvec].im;


          exp_term=exp((-0.25/SQR(Alpha[CurrentSystem]))*Rksq)/Rksq;
          energy_term=(SQR(Cksum.re)+SQR(Cksum.im))*exp_term; 

          f=factor*(SQR(Cksum.re)+SQR(Cksum.im));

          // energy sums
          USum+=f;

          if(Framework[CurrentSystem].FrameworkModel!=FLEXIBLE)
          {
            USum-=factor*
                  (SQR(StoreRigidChargeFramework[CurrentSystem][nvec].re)+SQR(StoreRigidChargeFramework[CurrentSystem][nvec].im));
            f-=factor*
                  (SQR(StoreRigidChargeFramework[CurrentSystem][nvec].re)+SQR(StoreRigidChargeFramework[CurrentSystem][nvec].im));
          }

          InverseLamdaSquared=0.25/(SQR(Alpha[CurrentSystem]))+1.0/Rksq;

          Theta.ax=1.0-2.0*Rk.x*Rk.x*InverseLamdaSquared;
          Theta.ay=-2.0*Rk.x*Rk.y*InverseLamdaSquared;
          Theta.az=-2.0*Rk.x*Rk.z*InverseLamdaSquared;

          Theta.bx=-2.0*Rk.y*Rk.x*InverseLamdaSquared;
          Theta.by=1.0-2.0*Rk.y*Rk.y*InverseLamdaSquared;
          Theta.bz=-2.0*Rk.y*Rk.z*InverseLamdaSquared;

          Theta.cx=-2.0*Rk.z*Rk.x*InverseLamdaSquared;
          Theta.cy=-2.0*Rk.z*Rk.y*InverseLamdaSquared;
          Theta.cz=1.0-2.0*Rk.z*Rk.z*InverseLamdaSquared;

          StrainFirstDerivative->ax-=f*Theta.ax;
          StrainFirstDerivative->bx-=f*Theta.bx;
          StrainFirstDerivative->cx-=f*Theta.cx;

          StrainFirstDerivative->ay-=f*Theta.ay;
          StrainFirstDerivative->by-=f*Theta.by;
          StrainFirstDerivative->cy-=f*Theta.cy;

          StrainFirstDerivative->az-=f*Theta.az;
          StrainFirstDerivative->bz-=f*Theta.bz;
          StrainFirstDerivative->cz-=f*Theta.cz;

          // Stress tensor
          index1=0;
          if(ComputeGradient)
          {
            GradientStrain(f,Gradient,Theta);

            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
                {
                  type=Framework[CurrentSystem].Atoms[f1][i].Type;
                  charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
                  considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                  if(considered_charged)
                  { 
                    index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

                    f1_I=2.0*factor*(-Eikr[index1].im*Cksum.re+Eikr[index1].re*Cksum.im);
                    if(index_i.x>=0) Gradient[index_i.x]+=f1_I*Rk.x;
                    if(index_i.y>=0) Gradient[index_i.y]+=f1_I*Rk.y;
                    if(index_i.z>=0) Gradient[index_i.z]+=f1_I*Rk.z;

                    if((index_i.x>=0)||(index_i.y>=0)||(index_i.z>=0)) index1++;
                  }
                }
              }
            }
          }

          for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
          {
            TypeMolA=Adsorbates[CurrentSystem][I].Type;
            for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
            {
              comA.x=comA.y=comA.z=0.0;
              if(Components[TypeMolA].Groups[ig].Rigid)
                comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

              for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
              {
                i=Components[TypeMolA].Groups[ig].Atoms[ia];

                type=Adsorbates[CurrentSystem][I].Atoms[i].Type;
                charge=Adsorbates[CurrentSystem][I].Atoms[i].Charge;
                considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                if(considered_charged)
                {
                  f1_I=2.0*factor*(-Eikr[index1].im*Cksum.re+Eikr[index1].re*Cksum.im);

                  posA=Adsorbates[CurrentSystem][I].Atoms[i].Position;

                  if(Components[TypeMolA].Groups[ig].Rigid)
                  {
                    index_i=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndex;
                    index_i2=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndexOrientation;

                    index1_rigid=Adsorbates[CurrentSystem][I].Atoms[i].HessianAtomIndex;
                  }
                  else
                  {
                    index_i=Adsorbates[CurrentSystem][I].Atoms[i].HessianIndex;
                    index_i2=UNDEFINED_INT_VECTOR3;
                    index1_rigid=-1;
                  }

                  if(ComputeGradient)
                  {
                    if(index_i.x>=0) Gradient[index_i.x]+=f1_I*Rk.x;
                    if(index_i.y>=0) Gradient[index_i.y]+=f1_I*Rk.y;
                    if(index_i.z>=0) Gradient[index_i.z]+=f1_I*Rk.z;
                  }

                  if(Components[TypeMolA].Groups[ig].Rigid)
                  {
                    if(ComputeGradient)
                    {
                      if(index_i2.x>=0) Gradient[index_i2.x]+=f1_I*(Rk.x*DVecX[index1_rigid].x+Rk.y*DVecX[index1_rigid].y+Rk.z*DVecX[index1_rigid].z);
                      if(index_i2.y>=0) Gradient[index_i2.y]+=f1_I*(Rk.x*DVecY[index1_rigid].x+Rk.y*DVecY[index1_rigid].y+Rk.z*DVecY[index1_rigid].z);
                      if(index_i2.z>=0) Gradient[index_i2.z]+=f1_I*(Rk.x*DVecZ[index1_rigid].x+Rk.y*DVecZ[index1_rigid].y+Rk.z*DVecZ[index1_rigid].z);

                      GradientStrainILocal(Gradient,f1_I,Rk,posA,comA);
                    }

                    pos=Components[TypeMolA].Positions[i];
                    temp1=0.5*f1_I*((posA.y-comA.y)*Rk.x+(posA.x-comA.x)*Rk.y);
                    temp2=0.5*f1_I*((posA.z-comA.z)*Rk.x+(posA.x-comA.x)*Rk.z);
                    temp3=0.5*f1_I*((posA.z-comA.z)*Rk.y+(posA.y-comA.y)*Rk.z);

                    StrainFirstDerivative->ax-=f1_I*(posA.x-comA.x)*Rk.x;
                    StrainFirstDerivative->bx-=temp1;
                    StrainFirstDerivative->cx-=temp2;
                    StrainFirstDerivative->ay-=temp1;
                    StrainFirstDerivative->by-=f1_I*(posA.y-comA.y)*Rk.y;
                    StrainFirstDerivative->cy-=temp3;
                    StrainFirstDerivative->az-=temp2;
                    StrainFirstDerivative->bz-=temp3;
                    StrainFirstDerivative->cz-=f1_I*(posA.z-comA.z)*Rk.z;
                  }

                  if((index_i.x>=0)||(index_i.y>=0)||(index_i.z>=0)||(index_i2.x>=0)||(index_i2.y>=0)||(index_i2.z>=0))
                    index1++;
                }
              }
            }
          }

          for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
          {
            TypeMolA=Cations[CurrentSystem][I].Type;
            for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
            {
              comA.x=comA.y=comA.z=0.0;
              if(Components[TypeMolA].Groups[ig].Rigid)
                comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

              for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
              {
                i=Components[TypeMolA].Groups[ig].Atoms[ia];

                type=Cations[CurrentSystem][I].Atoms[i].Type;
                charge=Cations[CurrentSystem][I].Atoms[i].Charge;
                considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                if(considered_charged)
                {
                  f1_I=2.0*factor*(-Eikr[index1].im*Cksum.re+Eikr[index1].re*Cksum.im);

                  posA=Cations[CurrentSystem][I].Atoms[i].Position;

                  if(Components[TypeMolA].Groups[ig].Rigid)
                  {
                    index_i=Cations[CurrentSystem][I].Groups[ig].HessianIndex;
                    index_i2=Cations[CurrentSystem][I].Groups[ig].HessianIndexOrientation;

                    index1_rigid=Cations[CurrentSystem][I].Atoms[i].HessianAtomIndex;
                  }
                  else
                  {
                    index_i=Cations[CurrentSystem][I].Atoms[i].HessianIndex;
                    index_i2=UNDEFINED_INT_VECTOR3;
                    index1_rigid=-1;
                  }

                  if(ComputeGradient)
                  {
                    if(index_i.x>=0) Gradient[index_i.x]+=f1_I*Rk.x;
                    if(index_i.y>=0) Gradient[index_i.y]+=f1_I*Rk.y;
                    if(index_i.z>=0) Gradient[index_i.z]+=f1_I*Rk.z;
                  }

                  if(Components[TypeMolA].Groups[ig].Rigid)
                  {
                    if(ComputeGradient)
                    {
                      if(index_i2.x>=0) Gradient[index_i2.x]+=f1_I*(Rk.x*DVecX[index1_rigid].x+Rk.y*DVecX[index1_rigid].y+Rk.z*DVecX[index1_rigid].z);
                      if(index_i2.y>=0) Gradient[index_i2.y]+=f1_I*(Rk.x*DVecY[index1_rigid].x+Rk.y*DVecY[index1_rigid].y+Rk.z*DVecY[index1_rigid].z);
                      if(index_i2.z>=0) Gradient[index_i2.z]+=f1_I*(Rk.x*DVecZ[index1_rigid].x+Rk.y*DVecZ[index1_rigid].y+Rk.z*DVecZ[index1_rigid].z);

                      GradientStrainILocal(Gradient,f1_I,Rk,posA,comA);
                    }

                    pos=Components[TypeMolA].Positions[i];
                    temp1=0.5*f1_I*((posA.y-comA.y)*Rk.x+(posA.x-comA.x)*Rk.y);
                    temp2=0.5*f1_I*((posA.z-comA.z)*Rk.x+(posA.x-comA.x)*Rk.z);
                    temp3=0.5*f1_I*((posA.z-comA.z)*Rk.y+(posA.y-comA.y)*Rk.z);

                    StrainFirstDerivative->ax-=f1_I*(posA.x-comA.x)*Rk.x;
                    StrainFirstDerivative->bx-=temp1;
                    StrainFirstDerivative->cx-=temp2;
                    StrainFirstDerivative->ay-=temp1;
                    StrainFirstDerivative->by-=f1_I*(posA.y-comA.y)*Rk.y;
                    StrainFirstDerivative->cy-=temp3;
                    StrainFirstDerivative->az-=temp2;
                    StrainFirstDerivative->bz-=temp3;
                    StrainFirstDerivative->cz-=f1_I*(posA.z-comA.z)*Rk.z;
                  }

                  if((index_i.x>=0)||(index_i.y>=0)||(index_i.z>=0)||(index_i2.x>=0)||(index_i2.y>=0)||(index_i2.z>=0))
                    index1++;
                }
              }
            }
          }


          if(ComputeHessian)
          {
            HessianAtomicStrainStrainLocal(HessianMatrix,f,InverseLamdaSquared,Rk,Rksq,Theta);

            index1=0;
            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
                {
                  type=Framework[CurrentSystem].Atoms[f1][i].Type;
                  charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
                  considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                  if(considered_charged)
                  {
                    index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;
                    index_i2=UNDEFINED_INT_VECTOR3;

                    f1_I=2.0*factor*(-Eikr[index1].im*Cksum.re+Eikr[index1].re*Cksum.im);
                    f2_I=2.0*factor*(-Eikr[index1].re*Cksum.re-Eikr[index1].im*Cksum.im);

                    if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x]+=f2_I*Rk.x*Rk.x;
                    if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y]+=f2_I*Rk.x*Rk.y;
                    if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z]+=f2_I*Rk.x*Rk.z;
                    if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y]+=f2_I*Rk.y*Rk.y;
                    if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z]+=f2_I*Rk.y*Rk.z;
                    if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z]+=f2_I*Rk.z*Rk.z;
                    if((index_i.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.y][index_i.x]+=f2_I*Rk.x*Rk.y;
                    if((index_i.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.z][index_i.x]+=f2_I*Rk.x*Rk.z;
                    if((index_i.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.z][index_i.y]+=f2_I*Rk.y*Rk.z;

                    HessianAtomicPositionStrain(HessianMatrix,index_i,f1_I,Theta,Rk);
 
                    index2=0;
                    for(f2=0;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
                    {
                      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
                      {
                        type=Framework[CurrentSystem].Atoms[f2][j].Type;
                        charge=Framework[CurrentSystem].Atoms[f2][j].Charge;
                        considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                        if(considered_charged)
                        {
                          index_j=Framework[CurrentSystem].Atoms[f2][j].HessianIndex;
                          index_j2=UNDEFINED_INT_VECTOR3;

                          f2_IJ=2.0*factor*(-Eikr[index2].re*Eikr[index1].re-Eikr[index2].im*Eikr[index1].im);
                          if((index_i.x>=0)&&(index_j.x>=0)&&(index_i.x<=index_j.x)) HessianMatrix.element[index_i.x][index_j.x]-=f2_IJ*Rk.x*Rk.x;
                          if((index_i.x>=0)&&(index_j.y>=0)&&(index_i.x<=index_j.y)) HessianMatrix.element[index_i.x][index_j.y]-=f2_IJ*Rk.x*Rk.y;
                          if((index_i.x>=0)&&(index_j.z>=0)&&(index_i.x<=index_j.z)) HessianMatrix.element[index_i.x][index_j.z]-=f2_IJ*Rk.x*Rk.z;
                          if((index_i.y>=0)&&(index_j.x>=0)&&(index_i.y<=index_j.x)) HessianMatrix.element[index_i.y][index_j.x]-=f2_IJ*Rk.y*Rk.x;
                          if((index_i.y>=0)&&(index_j.y>=0)&&(index_i.y<=index_j.y)) HessianMatrix.element[index_i.y][index_j.y]-=f2_IJ*Rk.y*Rk.y;
                          if((index_i.y>=0)&&(index_j.z>=0)&&(index_i.y<=index_j.z)) HessianMatrix.element[index_i.y][index_j.z]-=f2_IJ*Rk.y*Rk.z;
                          if((index_i.z>=0)&&(index_j.x>=0)&&(index_i.z<=index_j.x)) HessianMatrix.element[index_i.z][index_j.x]-=f2_IJ*Rk.z*Rk.x;
                          if((index_i.z>=0)&&(index_j.y>=0)&&(index_i.z<=index_j.y)) HessianMatrix.element[index_i.z][index_j.y]-=f2_IJ*Rk.z*Rk.y;
                          if((index_i.z>=0)&&(index_j.z>=0)&&(index_i.z<=index_j.z)) HessianMatrix.element[index_i.z][index_j.z]-=f2_IJ*Rk.z*Rk.z;

                          if((index_j.x>=0)||(index_j.y>=0)||(index_j.z>=0)||(index_j2.x>=0)||(index_j2.y>=0)||(index_j2.z>=0))
                            index2++;
                        }
                      }
                    }

                    for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
                    {
                      TypeMolB=Adsorbates[CurrentSystem][J].Type;
                      for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                      {
                        if(Components[TypeMolB].Groups[jg].Rigid)
                          comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                        for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                        {
                          j=Components[TypeMolB].Groups[jg].Atoms[ja];
                          type=Adsorbates[CurrentSystem][J].Atoms[j].Type;
                          charge=Adsorbates[CurrentSystem][J].Atoms[j].Charge;
                          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                          if(considered_charged)
                          {
                            f2_IJ=2.0*factor*(-Eikr[index2].re*Eikr[index1].re-Eikr[index2].im*Eikr[index1].im);

                            if(Components[TypeMolB].Groups[jg].Rigid)
                              pos=Components[TypeMolB].Positions[j];

                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                              index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                              index2_rigid=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;

                              dot_product_j.x=Rk.x*DVecX[index2_rigid].x+Rk.y*DVecX[index2_rigid].y+Rk.z*DVecX[index2_rigid].z;
                              dot_product_j.y=Rk.x*DVecY[index2_rigid].x+Rk.y*DVecY[index2_rigid].y+Rk.z*DVecY[index2_rigid].z;
                              dot_product_j.z=Rk.x*DVecZ[index2_rigid].x+Rk.y*DVecZ[index2_rigid].y+Rk.z*DVecZ[index2_rigid].z;
                            }
                            else
                            {
                              index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                              index_j2=UNDEFINED_INT_VECTOR3;
                              index2_rigid=-1;

                              dot_product_j.x=dot_product_j.y=dot_product_j.z=0.0;
                            }

                            posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;
                            comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                            if(Components[TypeMolB].Groups[jg].Rigid)
                              HessianCenterOfMassStrainJ(HessianMatrix,index_i,f2_IJ,posB,comB,Rk);

                            if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]-=f2_IJ*Rk.x*Rk.x;
                            if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]-=f2_IJ*Rk.x*Rk.y;
                            if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]-=f2_IJ*Rk.x*Rk.z;
                            if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]-=f2_IJ*Rk.y*Rk.x;
                            if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]-=f2_IJ*Rk.y*Rk.y;
                            if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]-=f2_IJ*Rk.y*Rk.z;
                            if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]-=f2_IJ*Rk.z*Rk.x;
                            if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]-=f2_IJ*Rk.z*Rk.y;
                            if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]-=f2_IJ*Rk.z*Rk.z;

                            // com of I with orientation of J
                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              if((index_i.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.x][index_j2.x]-=f2_IJ*Rk.x*dot_product_j.x;
                              if((index_i.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.x][index_j2.y]-=f2_IJ*Rk.x*dot_product_j.y;
                              if((index_i.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.x][index_j2.z]-=f2_IJ*Rk.x*dot_product_j.z;
                              if((index_i.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.y][index_j2.x]-=f2_IJ*Rk.y*dot_product_j.x;
                              if((index_i.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.y][index_j2.y]-=f2_IJ*Rk.y*dot_product_j.y;
                              if((index_i.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.y][index_j2.z]-=f2_IJ*Rk.y*dot_product_j.z;
                              if((index_i.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.z][index_j2.x]-=f2_IJ*Rk.z*dot_product_j.x;
                              if((index_i.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.z][index_j2.y]-=f2_IJ*Rk.z*dot_product_j.y;
                              if((index_i.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.z][index_j2.z]-=f2_IJ*Rk.z*dot_product_j.z;
                            }

                            if((index_j.x>=0)||(index_j.y>=0)||(index_j.z>=0)||(index_j2.x>=0)||(index_j2.y>=0)||(index_j2.z>=0))
                              index2++;
                          }
                        }
                      }
                    }

                    for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
                    {
                      TypeMolB=Cations[CurrentSystem][J].Type;
                      for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                      {
                        if(Components[TypeMolB].Groups[jg].Rigid)
                          comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                        for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                        {
                          j=Components[TypeMolB].Groups[jg].Atoms[ja];
                          type=Cations[CurrentSystem][J].Atoms[j].Type;
                          charge=Cations[CurrentSystem][J].Atoms[j].Charge;
                          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                          if(considered_charged)
                          {
                            f2_IJ=2.0*factor*(-Eikr[index2].re*Eikr[index1].re-Eikr[index2].im*Eikr[index1].im);

                            if(Components[TypeMolB].Groups[jg].Rigid)
                              pos=Components[TypeMolB].Positions[j];

                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                              index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                              index2_rigid=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;

                              dot_product_j.x=Rk.x*DVecX[index2_rigid].x+Rk.y*DVecX[index2_rigid].y+Rk.z*DVecX[index2_rigid].z;
                              dot_product_j.y=Rk.x*DVecY[index2_rigid].x+Rk.y*DVecY[index2_rigid].y+Rk.z*DVecY[index2_rigid].z;
                              dot_product_j.z=Rk.x*DVecZ[index2_rigid].x+Rk.y*DVecZ[index2_rigid].y+Rk.z*DVecZ[index2_rigid].z;
                            }
                            else
                            {
                              index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                              index_j2=UNDEFINED_INT_VECTOR3;
                              index2_rigid=-1;

                              dot_product_j.x=dot_product_j.y=dot_product_j.z=0.0;
                            }

                            posB=Cations[CurrentSystem][J].Atoms[j].Position;
                            comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                            if(Components[TypeMolB].Groups[jg].Rigid)
                              HessianCenterOfMassStrainJ(HessianMatrix,index_i,f2_IJ,posB,comB,Rk);

                            if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]-=f2_IJ*Rk.x*Rk.x;
                            if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]-=f2_IJ*Rk.x*Rk.y;
                            if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]-=f2_IJ*Rk.x*Rk.z;
                            if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]-=f2_IJ*Rk.y*Rk.x;
                            if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]-=f2_IJ*Rk.y*Rk.y;
                            if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]-=f2_IJ*Rk.y*Rk.z;
                            if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]-=f2_IJ*Rk.z*Rk.x;
                            if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]-=f2_IJ*Rk.z*Rk.y;
                            if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]-=f2_IJ*Rk.z*Rk.z;

                              // com of I with orientation of J
                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              if((index_i.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.x][index_j2.x]-=f2_IJ*Rk.x*dot_product_j.x;
                              if((index_i.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.x][index_j2.y]-=f2_IJ*Rk.x*dot_product_j.y;
                              if((index_i.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.x][index_j2.z]-=f2_IJ*Rk.x*dot_product_j.z;
                              if((index_i.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.y][index_j2.x]-=f2_IJ*Rk.y*dot_product_j.x;
                              if((index_i.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.y][index_j2.y]-=f2_IJ*Rk.y*dot_product_j.y;
                              if((index_i.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.y][index_j2.z]-=f2_IJ*Rk.y*dot_product_j.z;
                              if((index_i.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.z][index_j2.x]-=f2_IJ*Rk.z*dot_product_j.x;
                              if((index_i.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.z][index_j2.y]-=f2_IJ*Rk.z*dot_product_j.y;
                              if((index_i.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.z][index_j2.z]-=f2_IJ*Rk.z*dot_product_j.z;
                            }

                            if((index_j.x>=0)||(index_j.y>=0)||(index_j.z>=0)||(index_j2.x>=0)||(index_j2.y>=0)||(index_j2.z>=0))
                              index2++;
                          }
                        }
                      }
                    }
                    if((index_i.x>=0)||(index_i.y>=0)||(index_i.z>=0)||(index_i2.x>=0)||(index_i2.y>=0)||(index_i2.z>=0))
                      index1++;
                  }//added
                }
              }
            } // end flexible framework

            // AF, AA and AC
            for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
            {
              TypeMolA=Adsorbates[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                  comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

                for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                {
                  i=Components[TypeMolA].Groups[ig].Atoms[ia];
                  type=Adsorbates[CurrentSystem][I].Atoms[i].Type;
                  charge=Adsorbates[CurrentSystem][I].Atoms[i].Charge;
                  considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                  if(considered_charged)
                  {
                    f1_I=2.0*factor*(-Eikr[index1].im*Cksum.re+Eikr[index1].re*Cksum.im);
                    f2_I=2.0*factor*(-Eikr[index1].re*Cksum.re-Eikr[index1].im*Cksum.im);

                    if(Components[TypeMolA].Groups[ig].Rigid)
                    {
                      index_i=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndex;
                      index_i2=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndexOrientation;

                      index1_rigid=Adsorbates[CurrentSystem][I].Atoms[i].HessianAtomIndex;

                      dot_product_i.x=Rk.x*DVecX[index1_rigid].x+Rk.y*DVecX[index1_rigid].y+Rk.z*DVecX[index1_rigid].z;
                      dot_product_i.y=Rk.x*DVecY[index1_rigid].x+Rk.y*DVecY[index1_rigid].y+Rk.z*DVecY[index1_rigid].z;
                      dot_product_i.z=Rk.x*DVecZ[index1_rigid].x+Rk.y*DVecZ[index1_rigid].y+Rk.z*DVecZ[index1_rigid].z;

                      dot_product_AX=Rk.x*DDVecAX[index1_rigid].x+Rk.y*DDVecAX[index1_rigid].y+Rk.z*DDVecAX[index1_rigid].z;
                      dot_product_BY=Rk.x*DDVecBY[index1_rigid].x+Rk.y*DDVecBY[index1_rigid].y+Rk.z*DDVecBY[index1_rigid].z;
                      dot_product_CZ=Rk.x*DDVecCZ[index1_rigid].x+Rk.y*DDVecCZ[index1_rigid].y+Rk.z*DDVecCZ[index1_rigid].z;
                      dot_product_AY=Rk.x*DDVecAY[index1_rigid].x+Rk.y*DDVecAY[index1_rigid].y+Rk.z*DDVecAY[index1_rigid].z;
                      dot_product_AZ=Rk.x*DDVecAZ[index1_rigid].x+Rk.y*DDVecAZ[index1_rigid].y+Rk.z*DDVecAZ[index1_rigid].z;
                      dot_product_BZ=Rk.x*DDVecBZ[index1_rigid].x+Rk.y*DDVecBZ[index1_rigid].y+Rk.z*DDVecBZ[index1_rigid].z;
                    }
                    else
                    {
                      index_i=Adsorbates[CurrentSystem][I].Atoms[i].HessianIndex;
                      index_i2=UNDEFINED_INT_VECTOR3;
                      index1_rigid=-1;

                      dot_product_i.x=dot_product_i.y=dot_product_i.z=0.0;
                      dot_product_AX=dot_product_BY=dot_product_CZ=0.0;
                      dot_product_AY=dot_product_AZ=dot_product_BZ=0.0;
                    }

                    posA=Adsorbates[CurrentSystem][I].Atoms[i].Position;
                    comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

                    if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x]+=f2_I*Rk.x*Rk.x;
                    if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y]+=f2_I*Rk.x*Rk.y;
                    if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z]+=f2_I*Rk.x*Rk.z;
                    if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y]+=f2_I*Rk.y*Rk.y;
                    if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z]+=f2_I*Rk.y*Rk.z;
                    if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z]+=f2_I*Rk.z*Rk.z;

                    if((index_i.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.y][index_i.x]+=f2_I*Rk.x*Rk.y;
                    if((index_i.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.z][index_i.x]+=f2_I*Rk.x*Rk.z;
                    if((index_i.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.z][index_i.y]+=f2_I*Rk.y*Rk.z;
                    

                    // Crossterm: derivative of the energy with respect to strain and position
                    HessianAtomicPositionStrain(HessianMatrix,index_i,f1_I,Theta,Rk);


                    if(Components[TypeMolA].Groups[ig].Rigid)
                    {
                      pos=Components[TypeMolA].Positions[i];

                      HessianCenterOfMassStrainI(HessianMatrix,index_i,f2_I,posA,comA,Rk);

                      if((index_i2.x>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i2.x][index_i2.x]+=f2_I*SQR(dot_product_i.x);
                      if((index_i2.y>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i2.y][index_i2.y]+=f2_I*SQR(dot_product_i.y);
                      if((index_i2.z>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.z][index_i2.z]+=f2_I*SQR(dot_product_i.z);
                      if((index_i2.x>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i2.x][index_i2.y]+=f2_I*dot_product_i.x*dot_product_i.y;
                      if((index_i2.x>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.x][index_i2.z]+=f2_I*dot_product_i.x*dot_product_i.z;
                      if((index_i2.y>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.y][index_i2.z]+=f2_I*dot_product_i.y*dot_product_i.z;

                      if((index_i.x>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i.x][index_i2.x]+=f2_I*Rk.x*dot_product_i.x;
                      if((index_i.x>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i.x][index_i2.y]+=f2_I*Rk.x*dot_product_i.y;
                      if((index_i.x>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i.x][index_i2.z]+=f2_I*Rk.x*dot_product_i.z;
                      if((index_i.y>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i.y][index_i2.x]+=f2_I*Rk.y*dot_product_i.x;
                      if((index_i.y>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i.y][index_i2.y]+=f2_I*Rk.y*dot_product_i.y;
                      if((index_i.y>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i.y][index_i2.z]+=f2_I*Rk.y*dot_product_i.z;
                      if((index_i.z>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i.z][index_i2.x]+=f2_I*Rk.z*dot_product_i.x;
                      if((index_i.z>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i.z][index_i2.y]+=f2_I*Rk.z*dot_product_i.y;
                      if((index_i.z>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i.z][index_i2.z]+=f2_I*Rk.z*dot_product_i.z;

                      if((index_i2.x>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i2.x][index_i2.x]+=f1_I*dot_product_AX;
                      if((index_i2.y>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i2.y][index_i2.y]+=f1_I*dot_product_BY;
                      if((index_i2.z>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.z][index_i2.z]+=f1_I*dot_product_CZ;
                      if((index_i2.x>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i2.x][index_i2.y]+=f1_I*dot_product_AY;
                      if((index_i2.x>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.x][index_i2.z]+=f1_I*dot_product_AZ;
                      if((index_i2.y>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.y][index_i2.z]+=f1_I*dot_product_BZ;

                      HessianOrientationStrainI(HessianMatrix,index_i2,index1_rigid,f1_I,f2_I,posA,comA,Rk,Theta);

                      // correction term for first Born term
                      HessianAtomicCorrectionStrainStrainI(HessianMatrix,f1_I,f2_I,posA,comA,Rk,Theta);
                    }

                    index2=0;
                    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
                         for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
                         {
                           type=Framework[CurrentSystem].Atoms[f1][j].Type;
                           charge=Framework[CurrentSystem].Atoms[f1][j].Charge;
                           considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                           if(considered_charged)
                           {
                             index_j=Framework[CurrentSystem].Atoms[f1][j].HessianIndex;
                             index_j2=UNDEFINED_INT_VECTOR3;

                             if((index_j.x>=0)||(index_j.y>=0)||(index_j.z>=0)||(index_j2.x>=0)||(index_j2.y>=0)||(index_j2.z>=0))
                               index2++;
                           }
                         }

                    for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
                    {
                      TypeMolB=Adsorbates[CurrentSystem][J].Type;
 
                      for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                      {
                        if(Components[TypeMolB].Groups[jg].Rigid)
                          comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                        for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                        {
                          j=Components[TypeMolB].Groups[jg].Atoms[ja];
                          type=Adsorbates[CurrentSystem][J].Atoms[j].Type;
                          charge=Adsorbates[CurrentSystem][J].Atoms[j].Charge;
                          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                          if(considered_charged)
                          {
                            f2_IJ=2.0*factor*(-Eikr[index2].re*Eikr[index1].re-Eikr[index2].im*Eikr[index1].im);
                            if(Components[TypeMolB].Groups[jg].Rigid)
                              pos=Components[TypeMolB].Positions[j];


                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                              index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                              index2_rigid=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;

                              dot_product_j.x=Rk.x*DVecX[index2_rigid].x+Rk.y*DVecX[index2_rigid].y+Rk.z*DVecX[index2_rigid].z;
                              dot_product_j.y=Rk.x*DVecY[index2_rigid].x+Rk.y*DVecY[index2_rigid].y+Rk.z*DVecY[index2_rigid].z;
                              dot_product_j.z=Rk.x*DVecZ[index2_rigid].x+Rk.y*DVecZ[index2_rigid].y+Rk.z*DVecZ[index2_rigid].z;
                            }
                            else
                            {
                              index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                              index_j2=UNDEFINED_INT_VECTOR3;
                              index2_rigid=-1;

                              dot_product_j.x=dot_product_j.y=dot_product_j.z=0.0;
                            }

                            posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;

                            if(Components[TypeMolB].Groups[jg].Rigid)
                              HessianCenterOfMassStrainJ(HessianMatrix,index_i,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianOrientationStrainJ(HessianMatrix,index_i2,index1_rigid,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianAtomicCorrectionStrainStrainJ(HessianMatrix,f2_IJ,posA,comA,posB,comB,Rk);

                            if(I<=J)
                            {
                              if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]-=f2_IJ*Rk.x*Rk.x;
                              if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]-=f2_IJ*Rk.x*Rk.y;
                              if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]-=f2_IJ*Rk.x*Rk.z;
                              if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]-=f2_IJ*Rk.y*Rk.x;
                              if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]-=f2_IJ*Rk.y*Rk.y;
                              if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]-=f2_IJ*Rk.y*Rk.z;
                              if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]-=f2_IJ*Rk.z*Rk.x;
                              if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]-=f2_IJ*Rk.z*Rk.y;
                              if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]-=f2_IJ*Rk.z*Rk.z;

                              // com of I with orientation of J
                              if(Components[TypeMolB].Groups[jg].Rigid)
                              {
                                if((index_i.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.x][index_j2.x]-=f2_IJ*Rk.x*dot_product_j.x;
                                if((index_i.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.x][index_j2.y]-=f2_IJ*Rk.x*dot_product_j.y;
                                if((index_i.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.x][index_j2.z]-=f2_IJ*Rk.x*dot_product_j.z;
                                if((index_i.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.y][index_j2.x]-=f2_IJ*Rk.y*dot_product_j.x;
                                if((index_i.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.y][index_j2.y]-=f2_IJ*Rk.y*dot_product_j.y;
                                if((index_i.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.y][index_j2.z]-=f2_IJ*Rk.y*dot_product_j.z;
                                if((index_i.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.z][index_j2.x]-=f2_IJ*Rk.z*dot_product_j.x;
                                if((index_i.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.z][index_j2.y]-=f2_IJ*Rk.z*dot_product_j.y;
                                if((index_i.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.z][index_j2.z]-=f2_IJ*Rk.z*dot_product_j.z;
                              }

                              // orientation of I with com of J
                              if(Components[TypeMolA].Groups[ig].Rigid)
                              {
                                if((index_i2.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.x][index_j.x]-=f2_IJ*Rk.x*dot_product_i.x;
                                if((index_i2.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.y][index_j.x]-=f2_IJ*Rk.x*dot_product_i.y;
                                if((index_i2.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.z][index_j.x]-=f2_IJ*Rk.x*dot_product_i.z;
                                if((index_i2.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.x][index_j.y]-=f2_IJ*Rk.y*dot_product_i.x;
                                if((index_i2.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.y][index_j.y]-=f2_IJ*Rk.y*dot_product_i.y;
                                if((index_i2.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.z][index_j.y]-=f2_IJ*Rk.y*dot_product_i.z;
                                if((index_i2.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.x][index_j.z]-=f2_IJ*Rk.z*dot_product_i.x;
                                if((index_i2.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.y][index_j.z]-=f2_IJ*Rk.z*dot_product_i.y;
                                if((index_i2.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.z][index_j.z]-=f2_IJ*Rk.z*dot_product_i.z;
                              }

                              // orientation of I with orientation of J
                              if((Components[TypeMolA].Groups[ig].Rigid)&&(Components[TypeMolB].Groups[jg].Rigid))
                              {
                                if((index_i2.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.x][index_j2.x]-=f2_IJ*dot_product_i.x*dot_product_j.x;
                                if((index_i2.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.x][index_j2.y]-=f2_IJ*dot_product_i.x*dot_product_j.y;
                                if((index_i2.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.x][index_j2.z]-=f2_IJ*dot_product_i.x*dot_product_j.z;
                                if((index_i2.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.y][index_j2.x]-=f2_IJ*dot_product_i.y*dot_product_j.x;
                                if((index_i2.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.y][index_j2.y]-=f2_IJ*dot_product_i.y*dot_product_j.y;
                                if((index_i2.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.y][index_j2.z]-=f2_IJ*dot_product_i.y*dot_product_j.z;
                                if((index_i2.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.z][index_j2.x]-=f2_IJ*dot_product_i.z*dot_product_j.x;
                                if((index_i2.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.z][index_j2.y]-=f2_IJ*dot_product_i.z*dot_product_j.y;
                                if((index_i2.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.z][index_j2.z]-=f2_IJ*dot_product_i.z*dot_product_j.z;
                              }
                            }

                             if((index_j.x>=0)||(index_j.y>=0)||(index_j.z>=0)||(index_j2.x>=0)||(index_j2.y>=0)||(index_j2.z>=0))
                              index2++;
                          }
                        }
                      }
                    }

                    for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
                    {
                      TypeMolB=Cations[CurrentSystem][J].Type;

                      for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                      {
                        if(Components[TypeMolB].Groups[jg].Rigid)
                          comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                        for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                        {
                          j=Components[TypeMolB].Groups[jg].Atoms[ja];
                          type=Cations[CurrentSystem][J].Atoms[j].Type;
                          charge=Cations[CurrentSystem][J].Atoms[j].Charge;
                          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                          if(considered_charged)
                          {
                            f2_IJ=2.0*factor*(-Eikr[index2].re*Eikr[index1].re-Eikr[index2].im*Eikr[index1].im);

                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                              index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                              index2_rigid=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                              pos=Components[TypeMolB].Positions[j];

                              dot_product_j.x=Rk.x*DVecX[index2_rigid].x+Rk.y*DVecX[index2_rigid].y+Rk.z*DVecX[index2_rigid].z;
                              dot_product_j.y=Rk.x*DVecY[index2_rigid].x+Rk.y*DVecY[index2_rigid].y+Rk.z*DVecY[index2_rigid].z;
                              dot_product_j.z=Rk.x*DVecZ[index2_rigid].x+Rk.y*DVecZ[index2_rigid].y+Rk.z*DVecZ[index2_rigid].z;
                            }
                            else
                            {
                              index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                              index_j2=UNDEFINED_INT_VECTOR3;
                              index2_rigid=-1;

                              dot_product_j.x=dot_product_j.y=dot_product_j.z=0.0;
                            }

                            posB=Cations[CurrentSystem][J].Atoms[j].Position;

                            if(Components[TypeMolB].Groups[jg].Rigid)
                              HessianCenterOfMassStrainJ(HessianMatrix,index_i,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianOrientationStrainJ(HessianMatrix,index_i2,index1_rigid,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianAtomicCorrectionStrainStrainJ(HessianMatrix,f2_IJ,posA,comA,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid)
                              HessianCenterOfMassStrainJ(HessianMatrix,index_j,f2_IJ,posA,comA,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianOrientationStrainJ(HessianMatrix,index_j2,index2_rigid,f2_IJ,posA,comA,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianAtomicCorrectionStrainStrainJ(HessianMatrix,f2_IJ,posB,comB,posA,comA,Rk);

                            if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]-=f2_IJ*Rk.x*Rk.x;
                            if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]-=f2_IJ*Rk.x*Rk.y;
                            if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]-=f2_IJ*Rk.x*Rk.z;
                            if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]-=f2_IJ*Rk.y*Rk.x;
                            if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]-=f2_IJ*Rk.y*Rk.y;
                            if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]-=f2_IJ*Rk.y*Rk.z;
                            if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]-=f2_IJ*Rk.z*Rk.x;
                            if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]-=f2_IJ*Rk.z*Rk.y;
                            if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]-=f2_IJ*Rk.z*Rk.z;

                            // com of I with orientation of J
                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              if((index_i.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.x][index_j2.x]-=f2_IJ*Rk.x*dot_product_j.x;
                              if((index_i.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.x][index_j2.y]-=f2_IJ*Rk.x*dot_product_j.y;
                              if((index_i.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.x][index_j2.z]-=f2_IJ*Rk.x*dot_product_j.z;
                              if((index_i.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.y][index_j2.x]-=f2_IJ*Rk.y*dot_product_j.x;
                              if((index_i.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.y][index_j2.y]-=f2_IJ*Rk.y*dot_product_j.y;
                              if((index_i.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.y][index_j2.z]-=f2_IJ*Rk.y*dot_product_j.z;
                              if((index_i.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.z][index_j2.x]-=f2_IJ*Rk.z*dot_product_j.x;
                              if((index_i.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.z][index_j2.y]-=f2_IJ*Rk.z*dot_product_j.y;
                              if((index_i.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.z][index_j2.z]-=f2_IJ*Rk.z*dot_product_j.z;
                            }

                            // orientation of I with com of J
                            if(Components[TypeMolA].Groups[ig].Rigid)
                            {
                              if((index_i2.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.x][index_j.x]-=f2_IJ*Rk.x*dot_product_i.x;
                              if((index_i2.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.y][index_j.x]-=f2_IJ*Rk.x*dot_product_i.y;
                              if((index_i2.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.z][index_j.x]-=f2_IJ*Rk.x*dot_product_i.z;
                              if((index_i2.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.x][index_j.y]-=f2_IJ*Rk.y*dot_product_i.x;
                              if((index_i2.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.y][index_j.y]-=f2_IJ*Rk.y*dot_product_i.y;
                              if((index_i2.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.z][index_j.y]-=f2_IJ*Rk.y*dot_product_i.z;
                              if((index_i2.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.x][index_j.z]-=f2_IJ*Rk.z*dot_product_i.x;
                              if((index_i2.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.y][index_j.z]-=f2_IJ*Rk.z*dot_product_i.y;
                              if((index_i2.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.z][index_j.z]-=f2_IJ*Rk.z*dot_product_i.z;
                            }

                            // orientation of I with orientation of J
                            if((Components[TypeMolA].Groups[ig].Rigid)&&(Components[TypeMolB].Groups[jg].Rigid))
                            {
                              if((index_i2.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.x][index_j2.x]-=f2_IJ*dot_product_i.x*dot_product_j.x;
                              if((index_i2.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.x][index_j2.y]-=f2_IJ*dot_product_i.x*dot_product_j.y;
                              if((index_i2.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.x][index_j2.z]-=f2_IJ*dot_product_i.x*dot_product_j.z;
                              if((index_i2.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.y][index_j2.x]-=f2_IJ*dot_product_i.y*dot_product_j.x;
                              if((index_i2.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.y][index_j2.y]-=f2_IJ*dot_product_i.y*dot_product_j.y;
                              if((index_i2.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.y][index_j2.z]-=f2_IJ*dot_product_i.y*dot_product_j.z;
                              if((index_i2.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.z][index_j2.x]-=f2_IJ*dot_product_i.z*dot_product_j.x;
                              if((index_i2.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.z][index_j2.y]-=f2_IJ*dot_product_i.z*dot_product_j.y;
                              if((index_i2.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.z][index_j2.z]-=f2_IJ*dot_product_i.z*dot_product_j.z;
                            }

                            if((index_j.x>=0)||(index_j.y>=0)||(index_j.z>=0)||(index_j2.x>=0)||(index_j2.y>=0)||(index_j2.z>=0))
                              index2++;
                          }
                        }
                      }
                    }

                    if((index_i.x>=0)||(index_i.y>=0)||(index_i.z>=0)||(index_i2.x>=0)||(index_i2.y>=0)||(index_i2.z>=0))
                      index1++;
                  }
                }
              }
            }

            // CF, CA and CC
            for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
            {
              TypeMolA=Cations[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                  comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

                for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                {
                  i=Components[TypeMolA].Groups[ig].Atoms[ia];
                  type=Cations[CurrentSystem][I].Atoms[i].Type;
                  charge=Cations[CurrentSystem][I].Atoms[i].Charge;
                  considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                  if(considered_charged)
                  {
                    f1_I=2.0*factor*(-Eikr[index1].im*Cksum.re+Eikr[index1].re*Cksum.im);
                    f2_I=2.0*factor*(-Eikr[index1].re*Cksum.re-Eikr[index1].im*Cksum.im);
                    if(Components[TypeMolA].Groups[ig].Rigid)
                    {
                      index_i=Cations[CurrentSystem][I].Groups[ig].HessianIndex;
                      index_i2=Cations[CurrentSystem][I].Groups[ig].HessianIndexOrientation;
                      index1_rigid=Cations[CurrentSystem][I].Atoms[i].HessianAtomIndex;

                      dot_product_i.x=Rk.x*DVecX[index1_rigid].x+Rk.y*DVecX[index1_rigid].y+Rk.z*DVecX[index1_rigid].z;
                      dot_product_i.y=Rk.x*DVecY[index1_rigid].x+Rk.y*DVecY[index1_rigid].y+Rk.z*DVecY[index1_rigid].z;
                      dot_product_i.z=Rk.x*DVecZ[index1_rigid].x+Rk.y*DVecZ[index1_rigid].y+Rk.z*DVecZ[index1_rigid].z;

                      dot_product_AX=Rk.x*DDVecAX[index1_rigid].x+Rk.y*DDVecAX[index1_rigid].y+Rk.z*DDVecAX[index1_rigid].z;
                      dot_product_BY=Rk.x*DDVecBY[index1_rigid].x+Rk.y*DDVecBY[index1_rigid].y+Rk.z*DDVecBY[index1_rigid].z;
                      dot_product_CZ=Rk.x*DDVecCZ[index1_rigid].x+Rk.y*DDVecCZ[index1_rigid].y+Rk.z*DDVecCZ[index1_rigid].z;
                      dot_product_AY=Rk.x*DDVecAY[index1_rigid].x+Rk.y*DDVecAY[index1_rigid].y+Rk.z*DDVecAY[index1_rigid].z;
                      dot_product_AZ=Rk.x*DDVecAZ[index1_rigid].x+Rk.y*DDVecAZ[index1_rigid].y+Rk.z*DDVecAZ[index1_rigid].z;
                      dot_product_BZ=Rk.x*DDVecBZ[index1_rigid].x+Rk.y*DDVecBZ[index1_rigid].y+Rk.z*DDVecBZ[index1_rigid].z;
                    }
                    else
                    {
                      index_i=Cations[CurrentSystem][I].Atoms[i].HessianIndex;
                      index_i2=UNDEFINED_INT_VECTOR3;
                      index1_rigid=-1;

                      dot_product_i.x=dot_product_i.y=dot_product_i.z=0.0;
                      dot_product_AX=dot_product_BY=dot_product_CZ=0.0;
                      dot_product_AY=dot_product_AZ=dot_product_BZ=0.0;
                    }

                    posA=Cations[CurrentSystem][I].Atoms[i].Position;
                    comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

                    if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x]+=f2_I*Rk.x*Rk.x;
                    if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y]+=f2_I*Rk.x*Rk.y;
                    if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z]+=f2_I*Rk.x*Rk.z;
                    if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y]+=f2_I*Rk.y*Rk.y;
                    if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z]+=f2_I*Rk.y*Rk.z;
                    if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z]+=f2_I*Rk.z*Rk.z;

                    if((index_i.y>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.y][index_i.x]+=f2_I*Rk.x*Rk.y;
                    if((index_i.z>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.z][index_i.x]+=f2_I*Rk.x*Rk.z;
                    if((index_i.z>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.z][index_i.y]+=f2_I*Rk.y*Rk.z;

                    // Crossterm: derivative of the energy with respect to strain and position
                    HessianAtomicPositionStrain(HessianMatrix,index_i,f1_I,Theta,Rk);

                    if(Components[TypeMolA].Groups[ig].Rigid)
                    {
                      pos=Components[TypeMolA].Positions[i];

                      HessianCenterOfMassStrainI(HessianMatrix,index_i,f2_I,posA,comA,Rk);

                      if((index_i2.x>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i2.x][index_i2.x]+=f2_I*dot_product_i.x*dot_product_i.x;
                      if((index_i2.y>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i2.y][index_i2.y]+=f2_I*dot_product_i.y*dot_product_i.y;
                      if((index_i2.z>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.z][index_i2.z]+=f2_I*dot_product_i.z*dot_product_i.z;
                      if((index_i2.x>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i2.x][index_i2.y]+=f2_I*dot_product_i.x*dot_product_i.y;
                      if((index_i2.x>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.x][index_i2.z]+=f2_I*dot_product_i.x*dot_product_i.z;
                      if((index_i2.y>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.y][index_i2.z]+=f2_I*dot_product_i.y*dot_product_i.z;

                      if((index_i.x>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i.x][index_i2.x]+=f2_I*Rk.x*dot_product_i.x;
                      if((index_i.x>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i.x][index_i2.y]+=f2_I*Rk.x*dot_product_i.y;
                      if((index_i.x>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i.x][index_i2.z]+=f2_I*Rk.x*dot_product_i.z;
                      if((index_i.y>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i.y][index_i2.x]+=f2_I*Rk.y*dot_product_i.x;
                      if((index_i.y>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i.y][index_i2.y]+=f2_I*Rk.y*dot_product_i.y;
                      if((index_i.y>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i.y][index_i2.z]+=f2_I*Rk.y*dot_product_i.z;
                      if((index_i.z>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i.z][index_i2.x]+=f2_I*Rk.z*dot_product_i.x;
                      if((index_i.z>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i.z][index_i2.y]+=f2_I*Rk.z*dot_product_i.y;
                      if((index_i.z>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i.z][index_i2.z]+=f2_I*Rk.z*dot_product_i.z;

                      pos=Components[TypeMolA].Positions[i];
                      if((index_i2.x>=0)&&(index_i2.x>=0)) HessianMatrix.element[index_i2.x][index_i2.x]+=f1_I*dot_product_AX;
                      if((index_i2.y>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i2.y][index_i2.y]+=f1_I*dot_product_BY;
                      if((index_i2.z>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.z][index_i2.z]+=f1_I*dot_product_CZ;
                      if((index_i2.x>=0)&&(index_i2.y>=0)) HessianMatrix.element[index_i2.x][index_i2.y]+=f1_I*dot_product_AY;
                      if((index_i2.x>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.x][index_i2.z]+=f1_I*dot_product_AZ;
                      if((index_i2.y>=0)&&(index_i2.z>=0)) HessianMatrix.element[index_i2.y][index_i2.z]+=f1_I*dot_product_BZ;

                      // derivative of stress with respect to the orientation 
                      HessianOrientationStrainI(HessianMatrix,index_i2,index1_rigid,f1_I,f2_I,posA,comA,Rk,Theta);

                      // correction term for first Born term
                      HessianAtomicCorrectionStrainStrainI(HessianMatrix,f1_I,f2_I,posA,comA,Rk,Theta);
                    }

                    index2=0;
                    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
                         for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
                         {
                           index_j=Framework[CurrentSystem].Atoms[f1][j].HessianIndex;
                           index_j2=UNDEFINED_INT_VECTOR3;

                           if((index_j.x>=0)||(index_j.y>=0)||(index_j.z>=0)||(index_j2.x>=0)||(index_j2.y>=0)||(index_j2.z>=0))
                             index2++;
                         }

                    for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
                    {
                      TypeMolB=Adsorbates[CurrentSystem][J].Type;
 
                      for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                      {
                        for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                        {
                          j=Components[TypeMolB].Groups[jg].Atoms[ja];
                          type=Adsorbates[CurrentSystem][J].Atoms[j].Type;
                          charge=Adsorbates[CurrentSystem][J].Atoms[j].Charge;
                          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                          if(considered_charged)
                          {
                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                              index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                            }
                            else
                            {
                              index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                              index_j2=UNDEFINED_INT_VECTOR3;
                            }

                            // SKIPPING

                            if((index_j.x>=0)||(index_j.y>=0)||(index_j.z>=0)||(index_j2.x>=0)||(index_j2.y>=0)||(index_j2.z>=0))
                              index2++;
                          }
                        }
                      }
                    }

                    for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
                    {
                      TypeMolB=Cations[CurrentSystem][J].Type;

                      for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                      {
                        if(Components[TypeMolB].Groups[jg].Rigid)
                          comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                        for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                        {
                          j=Components[TypeMolB].Groups[jg].Atoms[ja];
                          type=Cations[CurrentSystem][J].Atoms[j].Type;
                          charge=Cations[CurrentSystem][J].Atoms[j].Charge;
                          considered_charged=(fabs(charge)>1e-10)||(PseudoAtoms[type].IsPolarizable);
                          if(considered_charged)
                          {
                            f2_IJ=2.0*factor*
                                  (-Eikr[index2].re*Eikr[index1].re-Eikr[index2].im*Eikr[index1].im);
                            if(Components[TypeMolB].Groups[jg].Rigid)
                              pos=Components[TypeMolB].Positions[j];


                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                              index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                              index2_rigid=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;

                              dot_product_j.x=Rk.x*DVecX[index2_rigid].x+Rk.y*DVecX[index2_rigid].y+Rk.z*DVecX[index2_rigid].z;
                              dot_product_j.y=Rk.x*DVecY[index2_rigid].x+Rk.y*DVecY[index2_rigid].y+Rk.z*DVecY[index2_rigid].z;
                              dot_product_j.z=Rk.x*DVecZ[index2_rigid].x+Rk.y*DVecZ[index2_rigid].y+Rk.z*DVecZ[index2_rigid].z;
                            }
                            else
                            {
                              index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                              index_j2=UNDEFINED_INT_VECTOR3;
                              index2_rigid=-1;

                              dot_product_j.x=dot_product_j.y=dot_product_j.z=0.0;
                            }

                            posB=Cations[CurrentSystem][J].Atoms[j].Position;

                            if(Components[TypeMolB].Groups[jg].Rigid)
                              HessianCenterOfMassStrainJ(HessianMatrix,index_i,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianOrientationStrainJ(HessianMatrix,index_i2,index1_rigid,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianAtomicCorrectionStrainStrainJ(HessianMatrix,f2_IJ,posA,comA,posB,comB,Rk);

                            if(I<=J)
                            {
                              if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.x][index_j.x]-=f2_IJ*Rk.x*Rk.x;
                              if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.x][index_j.y]-=f2_IJ*Rk.x*Rk.y;
                              if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.x][index_j.z]-=f2_IJ*Rk.x*Rk.z;
                              if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.y][index_j.x]-=f2_IJ*Rk.y*Rk.x;
                              if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.y][index_j.y]-=f2_IJ*Rk.y*Rk.y;
                              if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.y][index_j.z]-=f2_IJ*Rk.y*Rk.z;
                              if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i.z][index_j.x]-=f2_IJ*Rk.z*Rk.x;
                              if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i.z][index_j.y]-=f2_IJ*Rk.z*Rk.y;
                              if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i.z][index_j.z]-=f2_IJ*Rk.z*Rk.z;

                              // com of I with orientation of J
                              if(Components[TypeMolB].Groups[jg].Rigid)
                              {
                                if((index_i.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.x][index_j2.x]-=f2_IJ*Rk.x*dot_product_j.x;
                                if((index_i.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.x][index_j2.y]-=f2_IJ*Rk.x*dot_product_j.y;
                                if((index_i.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.x][index_j2.z]-=f2_IJ*Rk.x*dot_product_j.z;
                                if((index_i.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.y][index_j2.x]-=f2_IJ*Rk.y*dot_product_j.x;
                                if((index_i.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.y][index_j2.y]-=f2_IJ*Rk.y*dot_product_j.y;
                                if((index_i.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.y][index_j2.z]-=f2_IJ*Rk.y*dot_product_j.z;
                                if((index_i.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i.z][index_j2.x]-=f2_IJ*Rk.z*dot_product_j.x;
                                if((index_i.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i.z][index_j2.y]-=f2_IJ*Rk.z*dot_product_j.y;
                                if((index_i.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i.z][index_j2.z]-=f2_IJ*Rk.z*dot_product_j.z;
                              }

                              // orientation of I with com of J
                              if(Components[TypeMolA].Groups[ig].Rigid)
                              {
                                if((index_i2.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.x][index_j.x]-=f2_IJ*Rk.x*dot_product_i.x;
                                if((index_i2.y>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.y][index_j.x]-=f2_IJ*Rk.x*dot_product_i.y;
                                if((index_i2.z>=0)&&(index_j.x>=0)) HessianMatrix.element[index_i2.z][index_j.x]-=f2_IJ*Rk.x*dot_product_i.z;
                                if((index_i2.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.x][index_j.y]-=f2_IJ*Rk.y*dot_product_i.x;
                                if((index_i2.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.y][index_j.y]-=f2_IJ*Rk.y*dot_product_i.y;
                                if((index_i2.z>=0)&&(index_j.y>=0)) HessianMatrix.element[index_i2.z][index_j.y]-=f2_IJ*Rk.y*dot_product_i.z;
                                if((index_i2.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.x][index_j.z]-=f2_IJ*Rk.z*dot_product_i.x;
                                if((index_i2.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.y][index_j.z]-=f2_IJ*Rk.z*dot_product_i.y;
                                if((index_i2.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_i2.z][index_j.z]-=f2_IJ*Rk.z*dot_product_i.z;
                              }

                              // orientation of I with orientation of J
                              if((Components[TypeMolA].Groups[ig].Rigid)&&(Components[TypeMolB].Groups[jg].Rigid))
                              {
                                if((index_i2.x>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.x][index_j2.x]-=f2_IJ*dot_product_i.x*dot_product_j.x;
                                if((index_i2.x>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.x][index_j2.y]-=f2_IJ*dot_product_i.x*dot_product_j.y;
                                if((index_i2.x>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.x][index_j2.z]-=f2_IJ*dot_product_i.x*dot_product_j.z;
                                if((index_i2.y>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.y][index_j2.x]-=f2_IJ*dot_product_i.y*dot_product_j.x;
                                if((index_i2.y>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.y][index_j2.y]-=f2_IJ*dot_product_i.y*dot_product_j.y;
                                if((index_i2.y>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.y][index_j2.z]-=f2_IJ*dot_product_i.y*dot_product_j.z;
                                if((index_i2.z>=0)&&(index_j2.x>=0)) HessianMatrix.element[index_i2.z][index_j2.x]-=f2_IJ*dot_product_i.z*dot_product_j.x;
                                if((index_i2.z>=0)&&(index_j2.y>=0)) HessianMatrix.element[index_i2.z][index_j2.y]-=f2_IJ*dot_product_i.z*dot_product_j.y;
                                if((index_i2.z>=0)&&(index_j2.z>=0)) HessianMatrix.element[index_i2.z][index_j2.z]-=f2_IJ*dot_product_i.z*dot_product_j.z;
                              }
                            }
                            
                            if((index_j.x>=0)||(index_j.y>=0)||(index_j.z>=0)||(index_j2.x>=0)||(index_j2.y>=0)||(index_j2.z>=0))
                              index2++;
                          }
                        }
                      }
                    }
                    if((index_i.x>=0)||(index_i.y>=0)||(index_i.z>=0)||(index_i2.x>=0)||(index_i2.y>=0)||(index_i2.z>=0))
                      index1++;
                  }
                }
              }
            }
          }

          // next wavevector
          nvec++;
        }
      }
    }
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {

    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[f1];i++)
      {
        A=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i].A;
        B=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i].B;

        index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
        index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;

        typeA=Framework[CurrentSystem].Atoms[f1][A].Type;
        typeB=Framework[CurrentSystem].Atoms[f1][B].Type;

        ChargeA=Framework[CurrentSystem].Atoms[f1][A].Charge;
        ChargeB=Framework[CurrentSystem].Atoms[f1][B].Charge;

        posA=Framework[CurrentSystem].Atoms[f1][A].Position;
        posB=Framework[CurrentSystem].Atoms[f1][B].Position;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        // Coulomb subtraction has a problem at small distances, leading to large negative values in the Hessian
        // i.e. minimization of Coesite fails (a high density mineral)
        // if the distance is smaller than 1e-4 Angsgtrom, use the limit values
        if(r>1e-4)
        {
          // add contribution to the energy
          (*Energy)-=COULOMBIC_CONVERSION_FACTOR*erf(Alpha[CurrentSystem]*r)*
                                ChargeA*ChargeB/r;

          DF=-(COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
             (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-erf(Alpha[CurrentSystem]*r))/
             (r*rr));

          DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
              (((-3.0*erf(Alpha[CurrentSystem]*r)/(r*rr))+
              (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
             (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));
        }
        else
        {
          (*Energy)-=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*2.0*Alpha[CurrentSystem]/sqrt(M_PI);
          DF=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*4.0*CUBE(Alpha[CurrentSystem])/(3.0*sqrt(M_PI));
          DDF=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*8.0*CUBE(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])/(5.0*sqrt(M_PI));
        }

        //if((index_i<0)&&(index_j<0)) continue;

        S.ax=-dr.x*dr.x;
        S.bx=-dr.y*dr.x;
        S.cx=-dr.z*dr.x;

        S.ay=-dr.x*dr.y;
        S.by=-dr.y*dr.y;
        S.cy=-dr.z*dr.y;

        S.az=-dr.x*dr.z;
        S.bz=-dr.y*dr.z;
        S.cz=-dr.z*dr.z;

        // add contribution to the first derivatives
        if(ComputeGradient)
        {
          if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
          if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
          if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

          if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
          if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
          if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

          GradientStrain(DF,Gradient,S);
        }

        StrainFirstDerivative->ax-=DF*S.ax;
        StrainFirstDerivative->bx-=DF*S.bx;
        StrainFirstDerivative->cx-=DF*S.cx;

        StrainFirstDerivative->ay-=DF*S.ay;
        StrainFirstDerivative->by-=DF*S.by;
        StrainFirstDerivative->cy-=DF*S.cy;

        StrainFirstDerivative->az-=DF*S.az;
        StrainFirstDerivative->bz-=DF*S.bz;
        StrainFirstDerivative->cz-=DF*S.cz;

        if(ComputeHessian)
        {
          Hessian.ax=DDF*dr.x*dr.x+DF; Hessian.bx=DDF*dr.y*dr.x;    Hessian.cx=DDF*dr.z*dr.x;
          Hessian.ay=DDF*dr.x*dr.y;    Hessian.by=DDF*dr.y*dr.y+DF; Hessian.cy=DDF*dr.z*dr.y;
          Hessian.az=DDF*dr.x*dr.z;    Hessian.bz=DDF*dr.y*dr.z;    Hessian.cz=DDF*dr.z*dr.z+DF;

          if((index_i.x>=0)&&(index_i.x>=0)) HessianMatrix.element[index_i.x][index_i.x]+=Hessian.ax;
          if((index_i.x>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.x][index_i.y]+=Hessian.ay;
          if((index_i.x>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.x][index_i.z]+=Hessian.az;
          if((index_i.y>=0)&&(index_i.y>=0)) HessianMatrix.element[index_i.y][index_i.y]+=Hessian.by;
          if((index_i.y>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.y][index_i.z]+=Hessian.bz;
          if((index_i.z>=0)&&(index_i.z>=0)) HessianMatrix.element[index_i.z][index_i.z]+=Hessian.cz;

          if((index_j.x>=0)&&(index_j.x>=0)) HessianMatrix.element[index_j.x][index_j.x]+=Hessian.ax;
          if((index_j.x>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.x][index_j.y]+=Hessian.ay;
          if((index_j.x>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.x][index_j.z]+=Hessian.az;
          if((index_j.y>=0)&&(index_j.y>=0)) HessianMatrix.element[index_j.y][index_j.y]+=Hessian.by;
          if((index_j.y>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.y][index_j.z]+=Hessian.bz;
          if((index_j.z>=0)&&(index_j.z>=0)) HessianMatrix.element[index_j.z][index_j.z]+=Hessian.cz;

          if((index_i.x>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.x)][MAX2(index_i.x,index_j.x)]-=Hessian.ax;
          if((index_i.x>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.y)][MAX2(index_i.x,index_j.y)]-=Hessian.ay;
          if((index_i.x>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.x,index_j.z)][MAX2(index_i.x,index_j.z)]-=Hessian.az;
          if((index_i.y>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.x)][MAX2(index_i.y,index_j.x)]-=Hessian.ay;
          if((index_i.y>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.y)][MAX2(index_i.y,index_j.y)]-=Hessian.by;
          if((index_i.y>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.y,index_j.z)][MAX2(index_i.y,index_j.z)]-=Hessian.bz;
          if((index_i.z>=0)&&(index_j.x>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.x)][MAX2(index_i.z,index_j.x)]-=Hessian.az;
          if((index_i.z>=0)&&(index_j.y>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.y)][MAX2(index_i.z,index_j.y)]-=Hessian.bz;
          if((index_i.z>=0)&&(index_j.z>=0)) HessianMatrix.element[MIN2(index_i.z,index_j.z)][MAX2(index_i.z,index_j.z)]-=Hessian.cz;

          HessianAtomicPositionStrainExcluded(HessianMatrix,index_i,index_j,DF,DDF,dr);
          HessianAtomicStrainStrainLocalExcluded(HessianMatrix,DF,DDF,dr);
        }
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      NumberOfExcludedPairs=Components[Type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        A=MIN2(Components[Type].ExcludedIntraChargeCharge[i].A,Components[Type].ExcludedIntraChargeCharge[i].B);
        B=MAX2(Components[Type].ExcludedIntraChargeCharge[i].A,Components[Type].ExcludedIntraChargeCharge[i].B);

        TypeA=Adsorbates[CurrentSystem][m].Atoms[A].Type;
        ChargeA=Adsorbates[CurrentSystem][m].Atoms[A].Charge;
        considered_charged=(fabs(ChargeA)>1e-10)||(PseudoAtoms[TypeA].IsPolarizable);
        if(considered_charged)
        {
          TypeB=Adsorbates[CurrentSystem][m].Atoms[B].Type;
          ChargeB=Adsorbates[CurrentSystem][m].Atoms[B].Charge;
          considered_charged=(fabs(ChargeB)>1e-10)||(PseudoAtoms[TypeB].IsPolarizable);
          if(considered_charged)
          {
            comA=posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
            comB=posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

            index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
            index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
            index_i2=UNDEFINED_INT_VECTOR3;
            index_j2=UNDEFINED_INT_VECTOR3;
            index1=-1;
            index2=-1;

            grpA=Components[Type].group[A];
            RigidI=Components[Type].Groups[grpA].Rigid;
            if(RigidI)
            {
              comA=Adsorbates[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
              index_i=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndex;
              index_i2=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
              index1=Adsorbates[CurrentSystem][m].Atoms[A].HessianAtomIndex;
            }

            grpB=Components[Type].group[B];
            RigidJ=Components[Type].Groups[grpB].Rigid;
            if(RigidJ)
            {
              comB=Adsorbates[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
              index_j=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndex;
              index_j2=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
              index2=Adsorbates[CurrentSystem][m].Atoms[B].HessianAtomIndex;
            }

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            // add contribution to the energy
            (*Energy)-=COULOMBIC_CONVERSION_FACTOR*erf(Alpha[CurrentSystem]*r)*
                               ChargeA*ChargeB/r;

            DF=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
               (erf(Alpha[CurrentSystem]*r)-2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
               (r*rr);

            DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                 (((-3.0*erf(Alpha[CurrentSystem]*r)/(r*rr))+
                 (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                 (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));

            //if((index_i<0)&&(index_j<0)) continue;

            // skip the remainder of the loop if both atoms belong to the same rigid group
            if((grpA==grpB)&&(Components[Type].Groups[grpA].Rigid)) continue;

            S.ax=-dr.x*dr.x;
            S.bx=-dr.y*dr.x;
            S.cx=-dr.z*dr.x;

            S.ay=-dr.x*dr.y;
            S.by=-dr.y*dr.y;
            S.cy=-dr.z*dr.y;

            S.az=-dr.x*dr.z;
            S.bz=-dr.y*dr.z;
            S.cz=-dr.z*dr.z;

            // add contribution to the first derivatives
            if(ComputeGradient)
            {
              if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
              if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
              if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

              if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
              if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
              if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

              if(RigidI)
              {
                GradientStrainI(Gradient,DF,dr,posA,comA);

                if(index_i2.x>=0) Gradient[index_i2.x]+=DF*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                if(index_i2.y>=0) Gradient[index_i2.y]+=DF*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                if(index_i2.z>=0) Gradient[index_i2.z]+=DF*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
              }

              if(RigidJ)
              {
                GradientStrainI(Gradient,-DF,dr,posB,comB);

                if(index_j2.x>=0) Gradient[index_j2.x]-=DF*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                if(index_j2.y>=0) Gradient[index_j2.y]-=DF*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                if(index_j2.z>=0) Gradient[index_j2.z]-=DF*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
              }

              GradientStrain(DF,Gradient,S);
            }

            StrainFirstDerivative->ax-=DF*S.ax;
            StrainFirstDerivative->bx-=DF*S.bx;
            StrainFirstDerivative->cx-=DF*S.cx;

            StrainFirstDerivative->ay-=DF*S.ay;
            StrainFirstDerivative->by-=DF*S.by;
            StrainFirstDerivative->cy-=DF*S.cy;

            StrainFirstDerivative->az-=DF*S.az;
            StrainFirstDerivative->bz-=DF*S.bz;
            StrainFirstDerivative->cz-=DF*S.cz;

            if(ComputeHessian)
            {
              HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);
              HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);

              HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB,RigidI,RigidJ);
              HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,DF,DDF,posA,comA,posB,comB,dr);

              HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
              HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB);
            }
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      NumberOfExcludedPairs=Components[Type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        A=MIN2(Components[Type].ExcludedIntraChargeCharge[i].A,Components[Type].ExcludedIntraChargeCharge[i].B);
        B=MAX2(Components[Type].ExcludedIntraChargeCharge[i].A,Components[Type].ExcludedIntraChargeCharge[i].B);

        TypeA=Cations[CurrentSystem][m].Atoms[A].Type;
        ChargeA=Cations[CurrentSystem][m].Atoms[A].Charge;
        considered_charged=(fabs(ChargeA)>1e-10)||(PseudoAtoms[TypeA].IsPolarizable);
        if(considered_charged)
        {
          TypeB=Cations[CurrentSystem][m].Atoms[B].Type;
          ChargeB=Cations[CurrentSystem][m].Atoms[B].Charge;
          considered_charged=(fabs(ChargeB)>1e-10)||(PseudoAtoms[TypeB].IsPolarizable);
          if(considered_charged)
          {
            comA=posA=Cations[CurrentSystem][m].Atoms[A].Position;
            comB=posB=Cations[CurrentSystem][m].Atoms[B].Position;

            index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
            index_j=Cations[CurrentSystem][m].Atoms[B].HessianIndex;
            index_i2=UNDEFINED_INT_VECTOR3;
            index_j2=UNDEFINED_INT_VECTOR3;
            index1=-1;
            index2=-1;

            grpA=Components[Type].group[A];
            RigidI=Components[Type].Groups[grpA].Rigid;
            if(RigidI)
            {
              comA=Cations[CurrentSystem][m].Groups[grpA].CenterOfMassPosition;
              index_i=Cations[CurrentSystem][m].Groups[grpA].HessianIndex;
              index_i2=Cations[CurrentSystem][m].Groups[grpA].HessianIndexOrientation;
              index1=Cations[CurrentSystem][m].Atoms[A].HessianAtomIndex;
            }

            grpB=Components[Type].group[B];
            RigidJ=Components[Type].Groups[grpB].Rigid;
            if(RigidJ)
            {
              comB=Cations[CurrentSystem][m].Groups[grpB].CenterOfMassPosition;
              index_j=Cations[CurrentSystem][m].Groups[grpB].HessianIndex;
              index_j2=Cations[CurrentSystem][m].Groups[grpB].HessianIndexOrientation;
              index2=Cations[CurrentSystem][m].Atoms[B].HessianAtomIndex;
            }

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            // add contribution to the energy
            (*Energy)-=COULOMBIC_CONVERSION_FACTOR*erf(Alpha[CurrentSystem]*r)*
                               ChargeA*ChargeB/r;

            DF=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
               (erf(Alpha[CurrentSystem]*r)-2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
               (r*rr);

            DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                 (((-3.0*erf(Alpha[CurrentSystem]*r)/(r*rr))+
                 (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                 (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));

            //if((index_i<0)&&(index_j<0)) continue;

            // skip the remainder of the loop if both atoms belong to the same rigid group
            if((grpA==grpB)&&(Components[Type].Groups[grpA].Rigid)) continue;

            S.ax=-dr.x*dr.x;
            S.bx=-dr.y*dr.x;
            S.cx=-dr.z*dr.x;

            S.ay=-dr.x*dr.y;
            S.by=-dr.y*dr.y;
            S.cy=-dr.z*dr.y;

            S.az=-dr.x*dr.z;
            S.bz=-dr.y*dr.z;
            S.cz=-dr.z*dr.z;

            // add contribution to the first derivatives
            if(ComputeGradient)
            {
              if(index_i.x>=0) Gradient[index_i.x]+=DF*dr.x;
              if(index_i.y>=0) Gradient[index_i.y]+=DF*dr.y;
              if(index_i.z>=0) Gradient[index_i.z]+=DF*dr.z;

              if(index_j.x>=0) Gradient[index_j.x]-=DF*dr.x;
              if(index_j.y>=0) Gradient[index_j.y]-=DF*dr.y;
              if(index_j.z>=0) Gradient[index_j.z]-=DF*dr.z;

              if(RigidI)
              {
                GradientStrainI(Gradient,DF,dr,posA,comA);

                if(index_i2.x>=0) Gradient[index_i2.x]+=DF*(dr.x*DVecX[index1].x+dr.y*DVecX[index1].y+dr.z*DVecX[index1].z);
                if(index_i2.y>=0) Gradient[index_i2.y]+=DF*(dr.x*DVecY[index1].x+dr.y*DVecY[index1].y+dr.z*DVecY[index1].z);
                if(index_i2.z>=0) Gradient[index_i2.z]+=DF*(dr.x*DVecZ[index1].x+dr.y*DVecZ[index1].y+dr.z*DVecZ[index1].z);
              }

              if(RigidJ)
              {
                GradientStrainI(Gradient,-DF,dr,posB,comB);

                if(index_j2.x>=0) Gradient[index_j2.x]-=DF*(dr.x*DVecX[index2].x+dr.y*DVecX[index2].y+dr.z*DVecX[index2].z);
                if(index_j2.y>=0) Gradient[index_j2.y]-=DF*(dr.x*DVecY[index2].x+dr.y*DVecY[index2].y+dr.z*DVecY[index2].z);
                if(index_j2.z>=0) Gradient[index_j2.z]-=DF*(dr.x*DVecZ[index2].x+dr.y*DVecZ[index2].y+dr.z*DVecZ[index2].z);
              }

              GradientStrain(DF,Gradient,S);
            }

            StrainFirstDerivative->ax-=DF*S.ax;
            StrainFirstDerivative->bx-=DF*S.bx;
            StrainFirstDerivative->cx-=DF*S.cx;

            StrainFirstDerivative->ay-=DF*S.ay;
            StrainFirstDerivative->by-=DF*S.by;
            StrainFirstDerivative->cy-=DF*S.cy;

            StrainFirstDerivative->az-=DF*S.az;
            StrainFirstDerivative->bz-=DF*S.bz;
            StrainFirstDerivative->cz-=DF*S.cz;

            if(ComputeHessian)
            {
              HessianCenterOfMassOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);
              HessianOrientationOrientation(HessianMatrix,index_i,index_i2,index_j,index_j2,index1,index2,DF,DDF,dr,1.0);

              HessianCenterOfMassStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB,RigidI,RigidJ);
              HessianOrientationStrain(HessianMatrix,index_i2,index_j2,index1,index2,DF,DDF,posA,comA,posB,comB,dr);

              HessianAtomicPositionPosition(HessianMatrix,index_i,index_j,DF,DDF,dr,1.0);
              HessianAtomicStrainStrain(HessianMatrix,index_i,index_j,DF,DDF,dr,posA,comA,posB,comB);
            }
          }
        }
      }
    }
  }

  Uself_sum=0.0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        Uself_sum+=SQR(Framework[CurrentSystem].Atoms[f1][i].Charge);
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      Uself_sum+=SQR(Adsorbates[CurrentSystem][i].Atoms[j].Charge);

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      Uself_sum+=SQR(Cations[CurrentSystem][i].Atoms[j].Charge);

  Uself_sum*=Alpha[CurrentSystem]/sqrt(M_PI);

  *Energy+=USum-COULOMBIC_CONVERSION_FACTOR*Uself_sum;

  return 0;
}


// Remark 10 February 2011: Check and fix atoms without charge
int CalculateEwaldFourierPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainFirstDerivative,int ComputeGradient,int ComputeHessian)
{
  return 0;
}

void WriteRestartEwald(FILE *FilePtr)
{
  int i;

  fwrite(&OmitEwaldFourier,sizeof(int),1,FilePtr);
  fwrite(Alpha,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(kvec,sizeof(INT_VECTOR3),NumberOfSystems,FilePtr);
  fwrite(NumberOfKVectors,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(ReciprocalCutOffSquared,sizeof(REAL),NumberOfSystems,FilePtr);

  fwrite(&EwaldPrecision,sizeof(REAL),1,FilePtr);
  fwrite(&DielectricConstantOfTheMedium,sizeof(REAL),1,FilePtr);

  fwrite(&MaxNumberOfWaveVectors,sizeof(int),1,FilePtr);
  fwrite(&MaxKvecX,sizeof(int),1,FilePtr);
  fwrite(&MaxKvecY,sizeof(int),1,FilePtr);
  fwrite(&MaxKvecZ,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfCoulombicSites,sizeof(int),1,FilePtr);
  fwrite(&MaxNumberOfBondDipoleSites,sizeof(int),1,FilePtr);
  fwrite(NumberOfKVectors,sizeof(int),NumberOfSystems,FilePtr);

  fwrite(NetChargeSystem,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(NetChargeFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(NetChargeCations,sizeof(REAL),NumberOfSystems,FilePtr);
  fwrite(NetChargeAdsorbates,sizeof(REAL),NumberOfSystems,FilePtr);

  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    for(i=0;i<NumberOfSystems;i++)
    {
      fwrite(KVectors[i],sizeof(VECTOR),NumberOfKVectors[i],FilePtr);
      fwrite(KFactor[i],sizeof(REAL),NumberOfKVectors[i],FilePtr);

      fwrite(StoreRigidChargeFramework[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fwrite(StoreRigidChargeCations[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fwrite(StoreRigidChargeAdsorbates[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);

      fwrite(StoreRigidBondDipolesFramework[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fwrite(StoreRigidBondDipolesCations[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fwrite(StoreRigidBondDipolesAdsorbates[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);

      fwrite(StoreTotalChargeFramework[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fwrite(StoreTotalChargeCations[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fwrite(StoreTotalChargeAdsorbates[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);

      fwrite(StoreTotalBondDipolesFramework[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fwrite(StoreTotalBondDipolesCations[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fwrite(StoreTotalBondDipolesAdsorbates[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
    }
  }
}

void ReadRestartEwald(FILE *FilePtr)
{
  int i;

  Alpha=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  kvec=(INT_VECTOR3*)calloc(NumberOfSystems,sizeof(INT_VECTOR3));
  NumberOfKVectors=(int*)calloc(NumberOfSystems,sizeof(int));
  ReciprocalCutOffSquared=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  fread(&OmitEwaldFourier,sizeof(int),1,FilePtr);
  fread(Alpha,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(kvec,sizeof(INT_VECTOR3),NumberOfSystems,FilePtr);
  fread(NumberOfKVectors,sizeof(int),NumberOfSystems,FilePtr);
  fread(ReciprocalCutOffSquared,sizeof(REAL),NumberOfSystems,FilePtr);

  fread(&EwaldPrecision,sizeof(REAL),1,FilePtr);
  fread(&DielectricConstantOfTheMedium,sizeof(REAL),1,FilePtr);

  fread(&MaxNumberOfWaveVectors,sizeof(int),1,FilePtr);
  fread(&MaxKvecX,sizeof(int),1,FilePtr);
  fread(&MaxKvecY,sizeof(int),1,FilePtr);
  fread(&MaxKvecZ,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfCoulombicSites,sizeof(int),1,FilePtr);
  fread(&MaxNumberOfBondDipoleSites,sizeof(int),1,FilePtr);
  fread(NumberOfKVectors,sizeof(int),NumberOfSystems,FilePtr);

  AllocateEwaldMemory();

  fread(NetChargeSystem,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(NetChargeFramework,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(NetChargeCations,sizeof(REAL),NumberOfSystems,FilePtr);
  fread(NetChargeAdsorbates,sizeof(REAL),NumberOfSystems,FilePtr);


  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
  {
    for(i=0;i<NumberOfSystems;i++)
    {
      fread(KVectors[i],sizeof(VECTOR),NumberOfKVectors[i],FilePtr);
      fread(KFactor[i],sizeof(REAL),NumberOfKVectors[i],FilePtr);

      fread(StoreRigidChargeFramework[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fread(StoreRigidChargeCations[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fread(StoreRigidChargeAdsorbates[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);

      fread(StoreRigidBondDipolesFramework[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fread(StoreRigidBondDipolesCations[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fread(StoreRigidBondDipolesAdsorbates[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);

      fread(StoreTotalChargeFramework[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fread(StoreTotalChargeCations[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fread(StoreTotalChargeAdsorbates[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);

      fread(StoreTotalBondDipolesFramework[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fread(StoreTotalBondDipolesCations[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
      fread(StoreTotalBondDipolesAdsorbates[i],sizeof(COMPLEX),NumberOfKVectors[i],FilePtr);
    }
  }

  //CurrentSystem=0;
  //SetupKVectors();
  //EwaldEnergyIon();
  //PrecomputeFixedEwaldContributions();
  //PrecomputeTotalEwaldContributions();
}
