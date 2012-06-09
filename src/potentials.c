/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'potentials.c' is part of RASPA.

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

// 30 September 2010: Morse-potentials and Feynman-Hibbs potentials added by Rachel Getman and Decai Yu (Northwestern University)


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "potentials.h"
#include "simulation.h"
#include "molecule.h"
#include "cbmc.h"
#include "utils.h"
#include "inter_energy.h"
#include "inter_force.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "internal_energy.h"
#include "internal_force.h"
#include "ewald.h"
#include "grids.h"
#include "cubic_spline_1d.h"
#include "spacegroup.h"
#include "thermo_baro_stats.h"
#include "statistics.h"
#include "warnings.h"
#include "sample.h"

int GeneralMixingRule;
int IndividualMixingRules;
int IndividualInteractions;

POTENTIAL BondTypes[NR_BOND_TYPES]=
  {{2,"HARMONIC_BOND"},
   {1,"CORE_SHELL_SPRING"},
   {3,"MORSE_BOND"},
   {2,"LJ_12_6_BOND"},
   {2,"LENNARD_JONES_BOND"},
   {3,"BUCKINGHAM_BOND"},
   {3,"RESTRAINED_HARMONIC_BOND"},
   {4,"QUARTIC_BOND"},
   {4,"CFF_QUARTIC_BOND"},
   {2,"MM3_BOND"},
   {0,"RIGID_BOND"},
   {1,"FIXED_BOND"},
   {0,"MEASURE_BOND"}};

POTENTIAL UreyBradleyTypes[NR_UREYBRADLEY_TYPES]=
  {{2,"HARMONIC_UREYBRADLEY"},
   {3,"MORSE_UREYBRADLEY"},
   {2,"LJ_12_6_UREYBRADLEY"},
   {2,"LENNARD_JONES_UREYBRADLEY"},
   {3,"BUCKINGHAM_UREYBRADLEY"},
   {3,"RESTRAINED_HARMONIC_UREYBRADLEY"},
   {4,"QUARTIC_UREYBRADLEY"},
   {4,"CFF_QUARTIC_UREYBRADLEY"},
   {2,"MM3_UREYBRADLEY"},
   {0,"RIGID_UREYBRADLEY"},
   {1,"FIXED_UREYBRADLEY"},
   {0,"MEASURE_UREYBRADLEY"}};


POTENTIAL BendTypes[NR_BEND_TYPES]=
  {{2,"HARMONIC_BEND"},
   {2,"CORE_SHELL_BEND"},
   {4,"QUARTIC_BEND"},
   {4,"CFF_QUARTIC_BEND"},
   {2,"HARMONIC_COSINE_BEND"},
   {2,"COSINE_BEND"},
   {2,"MM3_BEND"},
   {2,"MM3_IN_PLANE_BEND"},
   {1,"FIXED_BEND"},
   {0,"MEASURE_BEND"}};

POTENTIAL InversionBendTypes[NR_INVERSION_BEND_TYPES]=
  {{2,"HARMONIC_INVERSION"},
   {2,"HARMONIC_COSINE_INVERSION"},
   {1,"PLANAR_INVERSION"},
   {2,"MM3_INVERSION"},
   {2,"HARMONIC_INVERSION2"},
   {2,"HARMONIC_COSINE_INVERSION2"},
   {1,"PLANAR_INVERSION2"},
   {1,"FIXED_INVERSION_BEND"}};

POTENTIAL TorsionTypes[NR_TORSION_TYPES]=
  {{3,"CVFF_DIHEDRAL"},
   {2,"HARMONIC_DIHEDRAL"},
   {2,"HARMONIC_COSINE_DIHEDRAL"},
   {3,"THREE_COSINE_DIHEDRAL"},
   {3,"CFF_DIHEDRAL"},
   {3,"CFF_DIHEDRAL2"},
   {6,"FOURIER_SERIES_DIHEDRAL"},
   {6,"FOURIER_SERIES_DIHEDRAL2"},
   {6,"SIX_COSINE_DIHEDRAL"},
   {4,"TRAPPE_DIHEDRAL"},
   {4,"OPLS_DIHEDRAL"},
   {3,"MM3_DIHEDRAL"},
   {1,"FIXED_DIHEDRAL"}};

POTENTIAL ImproperTorsionTypes[NR_IMPROPER_TORSION_TYPES]=
  {{3,"CVFF_IMPROPER_DIHEDRAL"},
   {2,"HARMONIC_IMPROPER_DIHEDRAL"},
   {2,"HARMONIC_COSINE_IMPROPER_DIHEDRAL"},
   {3,"THREE_COSINE_IMPROPER_DIHEDRAL"},
   {3,"CFF_IMPROPER_DIHEDRAL"},
   {3,"CFF_IMPROPER_DIHEDRAL2"},
   {6,"FOURIER_SERIES_IMPROPER_DIHEDRAL"},
   {6,"FOURIER_SERIES_IMPROPER_DIHEDRAL2"},
   {6,"SIX_COSINE_IMPROPER_DIHEDRAL"},
   {4,"TRAPPE_IMPROPER_DIHEDRAL"},
   {4,"OPLS_IMPROPER_DIHEDRAL"},
   {3,"MM3_IMPROPER_DIHEDRAL"},
   {1,"FIXED_IMPROPER_DIHEDRAL"}};

POTENTIAL OutOfPlaneTypes[NR_OUT_OF_PLANE_TYPES]=
   {{1,"FIXED_OUT_OF_PLANE_DISTANCE"}};

POTENTIAL BondBondTypes[NR_BOND_BOND_TYPES]=
  {{3,"CVFF_BOND_BOND_CROSS"},
   {3,"CFF_BOND_BOND_CROSS"}};

POTENTIAL BondBendTypes[NR_BOND_BEND_TYPES]=
  {{5,"CVFF_BOND_BEND_CROSS"},
   {5,"CFF_BOND_BEND_CROSS"},
   {4,"MM3_BOND_BEND_CROSS"},
   {3,"TRUNCATED_HARMONIC"},
   {4,"SCREENED_HARMONIC"},
   {4,"SCREENED_VESSAL"},
   {4,"TRUNCATED_VESSAL"}};

POTENTIAL BendBendTypes[NR_BEND_BEND_TYPES]=
  {{3,"CVFF_BEND_BEND_CROSS"},
   {3,"CFF_BEND_BEND_CROSS"},
   {3,"MM3_BEND_BEND_CROSS"}};

POTENTIAL BondTorsionTypes[NR_BOND_TORSION_TYPES]=
  {{4,"MM3_BOND_TORSION_CROSS"}};

POTENTIAL BendTorsionTypes[NR_BEND_TORSION_TYPES]=
  {{3,"SMOOTHED_DIHEDRAL"},
   {3,"SMOOTHED_THREE_COSINE_DIHEDRAL"},
   {3,"NICHOLAS_DIHEDRAL"},
   {3,"SMOOTHED_CFF_DIHEDRAL"},
   {3,"SMOOTHED_CFF_DIHEDRAL2"},
   {3,"CVFF_BEND_TORSION_CROSS"},
   {3,"CFF_BEND_TORSION_CROSS"},
   {3,"SMOOTHED_CFF_BEND_TORSION_CROSS"}};

REAL (**PotentialParms)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS];
int **PotentialType;

int TailCorrections;
int ShiftPotentials;
int **TailCorrection;
int **ShiftPotential;

REAL SwitchingVDWFactors3[4];
REAL SwitchingVDWFactors5[6];
REAL SwitchingVDWFactors7[8];
REAL SwitchingChargeChargeFactors3[4];
REAL SwitchingChargeChargeFactors5[6];
REAL SwitchingChargeChargeFactors7[8];
REAL SwitchingChargeBondDipoleFactors3[4];
REAL SwitchingChargeBondDipoleFactors5[6];
REAL SwitchingChargeBondDipoleFactors7[8];
REAL SwitchingBondDipoleBondDipoleFactors3[4];
REAL SwitchingBondDipoleBondDipoleFactors5[6];
REAL SwitchingBondDipoleBondDipoleFactors7[8];

int NumberOfPseudoAtomsWithoutVDWInteraction=0;
char (*PseudoAtomsWithoutVDWInteraction)[32]=NULL;


void ParseForceFieldSelfParameters(char *Arguments,int i,char *PotentialName)
{
  double arg1,arg2,arg3,arg4,arg5,arg6;

  // zero potential
  if((strcasecmp(PotentialName,"NONE")==0)||(strcasecmp(PotentialName,"ZERO_POTENTIAL")==0))
    PotentialType[i][i]=ZERO_POTENTIAL;

  // 4*p_0*((p_1/r)^12-(p_1/r)^6)
  // ======================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_3/k_B [K]    (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"LENNARD_JONES")==0)||(strcasecmp(PotentialName,"LENNARD-JONES")==0))
  {
    PotentialType[i][i]=LENNARD_JONES;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
  // ======================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  if((strcasecmp(PotentialName,"LENNARD_JONES_SMOOTHED3")==0)||(strcasecmp(PotentialName,"LENNARD-JONES-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=LENNARD_JONES_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
  // ======================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  if((strcasecmp(PotentialName,"LENNARD_JONES_SMOOTHED5")==0)||(strcasecmp(PotentialName,"LENNARD-JONES-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=LENNARD_JONES_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // 4*p_0*((p_1/r)^12-(p_1/r)^6)
  // ======================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_3/k_B [K]    (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"WCA")==0)||(strcasecmp(PotentialName,"WCA")==0))
  {
    PotentialType[i][i]=WCA;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
  // =============================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2     [u]    reduced mass in unified atomic mass units
  // p_3/k_B [K]    (non-zero for a shifted potential)
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES")==0))
  {
    PotentialType[i][i]=FEYNMAN_HIBBS_LENNARD_JONES;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
  // ====================================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2     [u]    reduced mass in unified atomic mass units
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
  // ====================================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2     [u]    reduced mass in unified atomic mass units
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY; //epsilon
    PotentialParms[i][i][1]=arg2;                  //sigma
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*(4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2)
  // ==============================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2/T   [A^2]  correction factor
  // p_3/k_B [K]    (non-zero for a shifted potential)
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES2")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES2")==0))
  {
    PotentialType[i][i]=FEYNMAN_HIBBS_LENNARD_JONES2;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*(4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2)}*S(r)
  // =====================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2/T   [A^2]  correction factor
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES2-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*(4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2)}*S(r)
  // =====================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2/T   [A^2]  correction factor
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES2-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]-[(p_1/rc)^12-(p_1/rc)^6]}+[12*(p_1/rc)^12-6*(p_1/rc)^6]*(r-rc)/rc
  // ===============================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  if((strcasecmp(PotentialName,"LENNARD_JONES_SHIFTED_FORCE")==0)||(strcasecmp(PotentialName,"LENNARD-JONES-SHIFTED-FORCE")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=LENNARD_JONES_SHIFTED_FORCE;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]+[6*(p_1/rc)^12-3*(p_1/rc)^6]}*r^2/rc^2+[7*(p_1/rc)^12+4*(p_1/rc)^6]
  // =================================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  if((strcasecmp(PotentialName,"LENNARD_JONES_SHIFTED_FORCE2")==0)||(strcasecmp(PotentialName,"LENNARD-JONES-SHIFTED-FORCE2")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=LENNARD_JONES_SHIFTED_FORCE2;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY; //epsilon
    PotentialParms[i][i][1]=arg2;                  //sigma
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // p_0/r^12-p_1/r^6
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  // p_2/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"12_6")==0)||(strcasecmp(PotentialName,"12-6")==0)||
     (strcasecmp(PotentialName,"POTENTIAL_12_6")==0)||(strcasecmp(PotentialName,"POTENTIAL-12-6")==0))
  {
    PotentialType[i][i]=POTENTIAL_12_6;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {p_0/r^12-p_1/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  if((strcasecmp(PotentialName,"12_6_SMOOTHED3")==0)||(strcasecmp(PotentialName,"12-6-SMOOTHED3")==0)||
     (strcasecmp(PotentialName,"POTENTIAL_12_6_SMOOTHED3")==0)||(strcasecmp(PotentialName,"POTENTIAL-12-6-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=POTENTIAL_12_6_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {p_0/r^12-p_1/r^6]}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  if((strcasecmp(PotentialName,"12_6_SMOOTHED5")==0)||(strcasecmp(PotentialName,"12-6-SMOOTHED5")==0)||
     (strcasecmp(PotentialName,"POTENTIAL_12_6_SMOOTHED5")==0)||(strcasecmp(PotentialName,"POTENTIAL-12-6-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=POTENTIAL_12_6_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // p_0/r^12+p_1/r^6+p_2/r^2+p_3
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  // p_2/k_B [K A^2]
  // p_3/k_B [K] 
  // p_4/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"12_6_2_0")==0)||(strcasecmp(PotentialName,"12-6-2-0")==0)||
     (strcasecmp(PotentialName,"POTENTIAL_12_6_2_0")==0)||(strcasecmp(PotentialName,"POTENTAL-12-6-2-0")==0))
  {
    PotentialType[i][i]=POTENTIAL_12_6_2_0;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  // p_2/k_B [K A^2]
  // p_3/k_B [K] 
  if((strcasecmp(PotentialName,"12_6_2_0_SMOOTHED3")==0)||(strcasecmp(PotentialName,"12-6-2-0-SMOOTHED3")==0)||
     (strcasecmp(PotentialName,"POTENTIAL_12_6_2_0_SMOOTHED3")==0)||(strcasecmp(PotentialName,"POTENTAL-12-6-2-0-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=POTENTIAL_12_6_2_0_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  // p_2/k_B [K A^2]
  // p_3/k_B [K] 
  if((strcasecmp(PotentialName,"12_6_2_0_SMOOTHED5")==0)||(strcasecmp(PotentialName,"12-6-2-0-SMOOTHED5")==0)||
     (strcasecmp(PotentialName,"POTENTIAL_12_6_2_0_SMOOTHED5")==0)||(strcasecmp(PotentialName,"POTENTAL-12-6-2-0-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=POTENTIAL_12_6_2_0_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  // p_3/k_B [K]       (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"MORSE")==0)
  {
    PotentialType[i][i]=MORSE;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  if((strcasecmp(PotentialName,"MORSE_SMOOTHED3")==0)||(strcasecmp(PotentialName,"MORSE-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MORSE_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  if((strcasecmp(PotentialName,"MORSE_SMOOTHED5")==0)||(strcasecmp(PotentialName,"MORSE-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MORSE_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  // p_3/k_B [K]       (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"MORSE2")==0)
  {
    PotentialType[i][i]=MORSE2;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // {p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]}*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  if(strcasecmp(PotentialName,"MORSE2_SMOOTHED3")==0)
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MORSE2_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // {p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]}*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  if(strcasecmp(PotentialName,"MORSE2_SMOOTHED5")==0)
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MORSE2_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))*(r/p_2-2^(1/6))])^2-1}
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A]       reference distance
  // p_2/k_B [K]       (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"MORSE3")==0)
  {
    PotentialType[i][i]=MORSE3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))*(r/p_2-2^(1/6))])^2-1}*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A]       reference distance
  if(strcasecmp(PotentialName,"MORSE3_SMOOTHED3")==0)
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MORSE3_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))*(r/p_2-2^(1/6))])^2-1}*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A]       reference distance
  if(strcasecmp(PotentialName,"MORSE3_SMOOTHED5")==0)
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MORSE3_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // p_0/r^9-p_1/r^6
  // ======================================================================================
  // p_0/k_B [K A^9]
  // p_1/k_B [K A^6]
  // p_2/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"CFF_9_6")==0)||(strcasecmp(PotentialName,"CFF-9-6")==0))
  {
    PotentialType[i][i]=CFF_9_6;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY; 
    PotentialParms[i][i][2]=(REAL)0.0;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {p_0/r^9-p_1/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^9]
  // p_1/k_B [K A^6]
  if((strcasecmp(PotentialName,"CFF_9_6_SMOOTHED3")==0)||(strcasecmp(PotentialName,"CFF-9-6-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=CFF_9_6_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=(REAL)0.0;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {p_0/r^9-p_1/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^9]
  // p_1/k_B [K A^6]
  if((strcasecmp(PotentialName,"CFF_9_6_SMOOTHED5")==0)||(strcasecmp(PotentialName,"CFF-9-6-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=CFF_9_6_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=(REAL)0.0;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // p_0*[2*(p_1/r)^9-3(p_1/r)^6]
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A]
  // p_2/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"CFF_EPS_SIGMA")==0)||(strcasecmp(PotentialName,"CFF-EPS-SIGMA")==0))
  {
    PotentialType[i][i]=CFF_EPS_SIGMA;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A]
  if((strcasecmp(PotentialName,"CFF_EPS_SIGMA_SMOOTHED3")==0)||(strcasecmp(PotentialName,"CFF-EPS-SIGMA-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=CFF_EPS_SIGMA_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A]
  if((strcasecmp(PotentialName,"CFF_EPS_SIGMA_SMOOTHED5")==0)||(strcasecmp(PotentialName,"CFF-EPS-SIGMA-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=CFF_EPS_SIGMA_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=(REAL)0.0;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // p_0*exp(-p_1*r)-p_2/r^6
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"BUCKINGHAM")==0)
  {
    PotentialType[i][i]=BUCKINGHAM;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  if((strcasecmp(PotentialName,"BUCKINGHAM_SMOOTHED3")==0)||(strcasecmp(PotentialName,"BUCKINGHAM-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=BUCKINGHAM_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  if((strcasecmp(PotentialName,"BUCKINGHAM_SMOOTHED5")==0)||(strcasecmp(PotentialName,"BUCKINGHAM-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=BUCKINGHAM_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=(REAL)0.0;
  }
  // if(r<p_3) 1e10 else p_0*exp(-p_1*r)-p_2/r^6 
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3     [A]
  // p_4/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"BUCKINGHAM2")==0)
  {
    PotentialType[i][i]=BUCKINGHAM2;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // if(r<p_3) 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3     [A]
  if((strcasecmp(PotentialName,"BUCKINGHAM2_SMOOTHED3")==0)||(strcasecmp(PotentialName,"BUCKINGHAM2-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=BUCKINGHAM2_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // if(r<p_3) 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3     [A]
  if((strcasecmp(PotentialName,"BUCKINGHAM2_SMOOTHED5")==0)||(strcasecmp(PotentialName,"BUCKINGHAM2-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=BUCKINGHAM2_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]  if P>=3.02
  // sqrt(p_0^i*p_0^j)*192.27*P^2                    if P<3.02
  // ======================================================================================
  // p_0     [kcal/mol]
  // p_1     [A]
  // p_3     [kcal/mol]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"MM3_HYDROGEN_VDW")==0)||(strcasecmp(PotentialName,"MM3-HYDROGEN-VDW")==0))
  {
    PotentialType[i][i]=MM3_HYDROGEN_VDW;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[i][i][0]=arg1*KCAL_PER_MOL_TO_ENERGY;
    PotentialParms[i][i][1]=2.0*arg2;    // the sum of the VDW radii
    PotentialParms[i][i][2]=(REAL)0.0;
    PotentialParms[i][i][3]=arg3;
  }
  // sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]  if P>=3.02
  // sqrt(p_0^i*p_0^j)*192.27*P^2                    if P<3.02
  // ======================================================================================
  // p_0     [kcal/mol]
  // p_1     [A]
  // p_3     [kcal/mol]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"MM3_VDW")==0)||(strcasecmp(PotentialName,"MM3-VDW")==0))
  {
    PotentialType[i][i]=MM3_VDW;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KCAL_PER_MOL_TO_ENERGY;
    PotentialParms[i][i][1]=2.0*arg2;    // the sum of the VDW radii
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
  // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                    if P<3.02
  // ======================================================================================
  // p_0     [kcal/mol]
  // p_1     [A]
  if((strcasecmp(PotentialName,"MM3_VDW_SMOOTHED3")==0)||(strcasecmp(PotentialName,"MM3-VDW-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MM3_VDW_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KCAL_PER_MOL_TO_ENERGY;
    PotentialParms[i][i][1]=2.0*arg2;    // the sum of the VDW radii
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
  // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                    if P<3.02
  // ======================================================================================
  // p_0     [kcal/mol]
  // p_1     [A]
  if((strcasecmp(PotentialName,"MM3_VDW_SMOOTHED5")==0)||(strcasecmp(PotentialName,"MM3-VDW-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MM3_VDW_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KCAL_PER_MOL_TO_ENERGY;
    PotentialParms[i][i][1]=2.0*arg2;    // the sum of the VDW radii
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K]
  // p_3     [A^-1]
  // p_4/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"MATSUOKA_CLEMENTI_YOSHIMINE")==0)||(strcasecmp(PotentialName,"MATSUOKA-CLEMENTI-YOSHIMINE")==0))
  {
    PotentialType[i][i]=MATSUOKA_CLEMENTI_YOSHIMINE;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K]
  // p_3     [A^-1]
  if((strcasecmp(PotentialName,"MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3")==0)||(strcasecmp(PotentialName,"MATSUOKA-CLEMENTI-YOSHIMINE-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K]
  // p_3     [A^-1]
  if((strcasecmp(PotentialName,"MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5")==0)||(strcasecmp(PotentialName,"MATSUOKA-CLEMENTI-YOSHIMINE-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10
  // ======================================================================================
  // p_0/k_B [K]
  // p_1/k_B [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  // p_5/k_B [K A^10]
  // p_6/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"GENERIC")==0)
  {
    PotentialType[i][i]=GENERIC;
    sscanf(Arguments,"%lf %lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5,&arg6);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=arg6*KELVIN_TO_ENERGY;
    PotentialParms[i][i][6]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1/k_B [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  // p_5/k_B [K A^10]
  if((strcasecmp(PotentialName,"GENERIC_SMOOTHED3")==0)||(strcasecmp(PotentialName,"GENERIC-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=GENERIC_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5,&arg6);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=arg6*KELVIN_TO_ENERGY;
    PotentialParms[i][i][6]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1/k_B [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  // p_5/k_B [K A^10]
  if((strcasecmp(PotentialName,"GENERIC_SMOOTHED5")==0)||(strcasecmp(PotentialName,"GENERIC-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=GENERIC_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5,&arg6);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=arg6*KELVIN_TO_ENERGY;
    PotentialParms[i][i][6]=(REAL)0.0;
  }
  // p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3/k_B [K A^8]
  // p_4/k_B [K A^10]
  // p_5/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"PELLENQ_NICHOLSON")==0)||(strcasecmp(PotentialName,"PELLENQ-NICHOLSON")==0))
  {
    PotentialType[i][i]=PELLENQ_NICHOLSON;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3/k_B [K A^8]
  // p_4/k_B [K A^10]
  if((strcasecmp(PotentialName,"PELLENQ_NICHOLSON_SMOOTHED3")==0)||(strcasecmp(PotentialName,"PELLENQ-NICHOLSON-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=PELLENQ_NICHOLSON_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3/k_B [K A^8]
  // p_4/k_B [K A^10]
  if((strcasecmp(PotentialName,"PELLENQ_NICHOLSON_SMOOTHED5")==0)||(strcasecmp(PotentialName,"PELLENQ-NICHOLSON-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=PELLENQ_NICHOLSON_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=(REAL)0.0;
  }
  // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^12]
  // p_5/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"HYDRATED_ION_WATER")==0)||(strcasecmp(PotentialName,"HYDRATED-ION-WATER")==0))
  {
    PotentialType[i][i]=HYDRATED_ION_WATER;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^12]
  if((strcasecmp(PotentialName,"HYDRATED_ION_WATER_SMOOTHED3")==0)||(strcasecmp(PotentialName,"HYDRATED-ION-WATER-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=HYDRATED_ION_WATER_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^12]
  if((strcasecmp(PotentialName,"HYDRATED_ION_WATER_SMOOTHED5")==0)||(strcasecmp(PotentialName,"HYDRATED-ION-WATER-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=HYDRATED_ION_WATER_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=(REAL)0.0;
  }
  // p_0/r^p_1-p_2/r^p_3
  // ======================================================================================
  // p_0/k_B [K A^p_1]
  // p_1     [-]
  // p_2/k_B [K A^p_3]
  // p_3     [-]
  // p_4/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"MIE")==0)
  {
    PotentialType[i][i]=MIE;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // {p_0/r^p_1-p_2/r^p_3}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^p_1]
  // p_1     [-]
  // p_2/k_B [K A^p_3]
  // p_3     [-]
  if((strcasecmp(PotentialName,"MIE_SMOOTHED3")==0)||(strcasecmp(PotentialName,"MIE-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MIE_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // {p_0/r^p_1-p_2/r^p_3}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^p_1]
  // p_1     [-]
  // p_2/k_B [K A^p_3]
  // p_3     [-]
  if((strcasecmp(PotentialName,"MIE_SMOOTHED5")==0)||(strcasecmp(PotentialName,"MIE-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=MIE_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][i][3]=arg4;
    PotentialParms[i][i][4]=(REAL)0.0;
  }
  // p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [-]
  // p_2     [A]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  // p_5/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"BORN_HUGGINS_MEYER")==0)||(strcasecmp(PotentialName,"BORN-HUGGINS-MEYER")==0))
  {
    PotentialType[i][i]=BORN_HUGGINS_MEYER;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=(REAL)0.0;
  }
  // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [-]
  // p_2     [A]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  if((strcasecmp(PotentialName,"BORN_HUGGINS_MEYER_SMOOTHED3")==0)||(strcasecmp(PotentialName,"BORN-HUGGINS-MEYER-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=BORN_HUGGINS_MEYER_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=(REAL)0.0;
  }
  // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [-]
  // p_2     [A]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  if((strcasecmp(PotentialName,"BORN_HUGGINS_MEYER_SMOOTHED5")==0)||(strcasecmp(PotentialName,"BORN-HUGGINS-MEYER-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=BORN_HUGGINS_MEYER_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2;
    PotentialParms[i][i][2]=arg3;
    PotentialParms[i][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][i][5]=(REAL)0.0;
  }
  // p_0/r^12-p_1/r^10
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^10]
  // p_2/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"HYDROGEN")==0)
  {
    PotentialType[i][i]=HYDROGEN;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {p_0/r^12-p_1/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^10]
  if((strcasecmp(PotentialName,"HYDROGEN_SMOOTHED3")==0)||(strcasecmp(PotentialName,"HYDROGEN-SMOOTHED3")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=HYDROGEN_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
  // {p_0/r^12-p_1/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^10]
  if((strcasecmp(PotentialName,"HYDROGEN_SMOOTHED5")==0)||(strcasecmp(PotentialName,"HYDROGEN-SMOOTHED5")==0))
  {
    TailCorrection[i][i]=FALSE;
    ShiftPotential[i][i]=FALSE;
    PotentialType[i][i]=HYDROGEN_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[i][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][i][2]=(REAL)0.0;
  }
}

void ParseForceFieldBinaryParameters(char *Arguments,int i,int j,char *PotentialName)
{
  double arg1,arg2,arg3,arg4,arg5,arg6;

  // zero potential
  if(strcasecmp(PotentialName,"None")==0)
  {
    PotentialType[i][j]=ZERO_POTENTIAL;
    PotentialType[j][i]=ZERO_POTENTIAL;
  }
  if(strncasecmp(PotentialName,"zero-potential",strlen("zero-potential"))==0)
  {
    PotentialType[i][j]=ZERO_POTENTIAL;
    PotentialType[j][i]=ZERO_POTENTIAL;
  }
  // 4*p_0*((p_1/r)^12-(p_1/r)^6)
  // ======================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_3/k_B [K]    (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"LENNARD_JONES")==0)||(strcasecmp(PotentialName,"LENNARD-JONES")==0))
  {
    PotentialType[i][j]=LENNARD_JONES;
    PotentialType[j][i]=LENNARD_JONES;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
  // ======================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  if((strcasecmp(PotentialName,"LENNARD_JONES_SMOOTHED3")==0)||(strcasecmp(PotentialName,"LENNARD-JONES-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=LENNARD_JONES_SMOOTHED3;
    PotentialType[j][i]=LENNARD_JONES_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
  // ======================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  if((strcasecmp(PotentialName,"LENNARD_JONES_SMOOTHED5")==0)||(strcasecmp(PotentialName,"LENNARD-JONES-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=LENNARD_JONES_SMOOTHED5;
    PotentialType[j][i]=LENNARD_JONES_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // 4*p_0*((p_1/r)^12-(p_1/r)^6)
  // ======================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_3/k_B [K]    (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"WCA")==0)||(strcasecmp(PotentialName,"WCA-JONES")==0))
  {
    PotentialType[i][j]=WCA;
    PotentialType[j][i]=WCA;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
  // =============================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2     [u]    reduced mass in unified atomic mass units
  // p_3/k_B [K]    (non-zero for a shifted potential)
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES")==0))
  {
    PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES;
    PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
  // ====================================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2     [u]    reduced mass in unified atomic mass units
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
    PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;   
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;     
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
  // ====================================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2     [u]    reduced mass in unified atomic mass units
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
    PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*(4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2)
  // ==============================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2/T   [A^2]  correction factor
  // p_3/k_B [K]    (non-zero for a shifted potential)
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES2")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES2")==0))
  {
    PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2;
    PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*(4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2)}*S(r)
  // =====================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2/T   [A^2]  correction factor
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES2-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
    PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*(4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2)}*S(r)
  // =====================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  // p_2/T   [A^2]  correction factor
  // T       [K]    the temperature
  if((strcasecmp(PotentialName,"FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5")==0)||(strcasecmp(PotentialName,"FEYNMAN-HIBBS-LENNARD-JONES2-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
    PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]-[(p_1/rc)^12-(p_1/rc)^6]}+[12*(p_1/rc)^12-6*(p_1/rc)^6]*(r-rc)/rc
  // ===============================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  if((strcasecmp(PotentialName,"LENNARD_JONES_SHIFTED_FORCE")==0)||(strcasecmp(PotentialName,"LENNARD-JONES-SHIFTED-FORCE")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=LENNARD_JONES_SHIFTED_FORCE;
    PotentialType[j][i]=LENNARD_JONES_SHIFTED_FORCE;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]+[6*(p_1/rc)^12-3*(p_1/rc)^6]}*r^2/rc^2+[7*(p_1/rc)^12+4*(p_1/rc)^6]
  // =================================================================================================
  // p_0/k_B [K]    strength parameter epsilon
  // p_1     [A]    size parameter sigma
  if((strcasecmp(PotentialName,"LENNARD_JONES_SHIFTED_FORCE2")==0)||(strcasecmp(PotentialName,"LENNARD_JONES_SHIFTED_FORCE2")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=LENNARD_JONES_SHIFTED_FORCE2;
    PotentialType[j][i]=LENNARD_JONES_SHIFTED_FORCE2;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // p_0/r^12-p_1/r^6
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  // p_2/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"12_6")==0)||(strcasecmp(PotentialName,"12-6")==0)||
     (strcasecmp(PotentialName,"POTENTIAL_12_6")==0)||(strcasecmp(PotentialName,"POTENTIAL-12-6")==0))
  {
    PotentialType[i][j]=POTENTIAL_12_6;
    PotentialType[j][i]=POTENTIAL_12_6;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {p_0/r^12-p_1/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  if((strcasecmp(PotentialName,"12_6_SMOOTHED3")==0)||(strcasecmp(PotentialName,"12-6-SMOOTHED3")==0)||
     (strcasecmp(PotentialName,"POTENTIAL_12_6_SMOOTHED3")==0)||(strcasecmp(PotentialName,"POTENTIAL-12-6-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=POTENTIAL_12_6_SMOOTHED3;
    PotentialType[j][i]=POTENTIAL_12_6_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {p_0/r^12-p_1/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  if((strcasecmp(PotentialName,"12_6_SMOOTHED5")==0)||(strcasecmp(PotentialName,"12-6-SMOOTHED5")==0)||
     (strcasecmp(PotentialName,"POTENTIAL_12_6_SMOOTHED5")==0)||(strcasecmp(PotentialName,"POTENTIAL-12-6-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=POTENTIAL_12_6_SMOOTHED5;
    PotentialType[j][i]=POTENTIAL_12_6_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // p_0/r^12+p_1/r^6+p_2/r^2+p_3
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  // p_2/k_B [K A^2]
  // p_3/k_B [K] 
  // p_4/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"12_6_2_0")==0)||(strcasecmp(PotentialName,"12-6-2-0")==0))
  {
    PotentialType[i][j]=POTENTIAL_12_6_2_0;
    PotentialType[j][i]=POTENTIAL_12_6_2_0;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  // p_2/k_B [K A^2]
  // p_3/k_B [K] 
  if((strcasecmp(PotentialName,"12_6_2_0_SMOOTHED3")==0)||(strcasecmp(PotentialName,"12-6-2-0-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=POTENTIAL_12_6_2_0_SMOOTHED3;
    PotentialType[j][i]=POTENTIAL_12_6_2_0_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^6]
  // p_2/k_B [K A^2]
  // p_3/k_B [K] 
  if((strcasecmp(PotentialName,"12_6_2_0_SMOOTHED5")==0)||(strcasecmp(PotentialName,"12-6-2-0-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=POTENTIAL_12_6_2_0_SMOOTHED5;
    PotentialType[j][i]=POTENTIAL_12_6_2_0_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  // p_3/k_B [K]       (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"MORSE")==0)
  {
    PotentialType[i][j]=MORSE;
    PotentialType[j][i]=MORSE;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][j][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  if(strcasecmp(PotentialName,"MORSE_SMOOTHED3")==0)
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MORSE_SMOOTHED3;
    PotentialType[j][i]=MORSE_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][j][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  if(strcasecmp(PotentialName,"MORSE_SMOOTHED5")==0)
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MORSE_SMOOTHED5;
    PotentialType[j][i]=MORSE_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][j][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  // p_3/k_B [K]       (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"MORSE2")==0)
  {
    PotentialType[i][j]=MORSE2;
    PotentialType[j][i]=MORSE2;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][j][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  if(strcasecmp(PotentialName,"MORSE2_SMOOTHED3")==0)
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MORSE2_SMOOTHED3;
    PotentialType[j][i]=MORSE2_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][j][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A^-1]    parameter
  // p_2     [A]       reference distance
  if(strcasecmp(PotentialName,"MORSE2_SMOOTHED5")==0)
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MORSE2_SMOOTHED5;
    PotentialType[j][i]=MORSE2_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][j][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))*(r/p_2-2^(1/6))])^2-1}
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A]       reference distance
  // p_2/k_B [K]       (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"MORSE3")==0)
  {
    PotentialType[i][j]=MORSE3;
    PotentialType[j][i]=MORSE3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][j][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))*(r/p_2-2^(1/6))])^2-1}*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A]       reference distance
  if(strcasecmp(PotentialName,"MORSE3_SMOOTHED3")==0)
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MORSE3_SMOOTHED3;
    PotentialType[j][i]=MORSE3_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][j][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))*(r/p_2-2^(1/6))])^2-1}*S(r)
  // =================================================================================
  // p_0/k_B [K]       force constant
  // p_1     [A]       reference distance
  if(strcasecmp(PotentialName,"MORSE3_SMOOTHED5")==0)
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MORSE3_SMOOTHED5;
    PotentialType[j][i]=MORSE3_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[i][j][0]=arg1/ENERGY_TO_KELVIN;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // p_0/r^9-p_1/r^6
  // ======================================================================================
  // p_0/k_B [K A^9]
  // p_1/k_B [K A^6]
  // p_2/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"CFF_9_6")==0)||(strcasecmp(PotentialName,"CFF-9-6")==0))
  {
    PotentialType[i][j]=CFF_9_6;
    PotentialType[j][i]=CFF_9_6;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {p_0/r^9-p_1/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^9]
  // p_1/k_B [K A^6]
  if((strcasecmp(PotentialName,"CFF_9_6_SMOOTHED3")==0)||(strcasecmp(PotentialName,"CFF-9-6-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=CFF_9_6_SMOOTHED3;
    PotentialType[j][i]=CFF_9_6_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;   
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {p_0/r^9-p_1/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^9]
  // p_1/k_B [K A^6]
  if((strcasecmp(PotentialName,"CFF_9_6_SMOOTHED5")==0)||(strcasecmp(PotentialName,"CFF-9-6-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=CFF_9_6_SMOOTHED5;
    PotentialType[j][i]=CFF_9_6_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // p_0*[2*(p_1/r)^9-3(p_1/r)^6]
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A]
  // p_2/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"CFF_EPS_SIGMA")==0)||(strcasecmp(PotentialName,"CFF-EPS-SIGMA")==0))
  {
    PotentialType[i][j]=CFF_EPS_SIGMA;
    PotentialType[j][i]=CFF_EPS_SIGMA;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A]
  if((strcasecmp(PotentialName,"CFF_EPS_SIGMA_SMOOTHED3")==0)||(strcasecmp(PotentialName,"CFF-EPS-SIGMA-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=CFF_EPS_SIGMA_SMOOTHED3;
    PotentialType[j][i]=CFF_EPS_SIGMA_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;   
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;     
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A]
  if((strcasecmp(PotentialName,"CFF_EPS_SIGMA_SMOOTHED5")==0)||(strcasecmp(PotentialName,"CFF-EPS-SIGMA-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=CFF_EPS_SIGMA_SMOOTHED5;
    PotentialType[j][i]=CFF_EPS_SIGMA_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // p_0*exp(-p_1*r)-p_2/r^6
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"BUCKINGHAM")==0)
  {
    PotentialType[i][j]=BUCKINGHAM;
    PotentialType[j][i]=BUCKINGHAM;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  if((strcasecmp(PotentialName,"BUCKINGHAM_SMOOTHED3")==0)||(strcasecmp(PotentialName,"BUCKINGHAM-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=BUCKINGHAM_SMOOTHED3;
    PotentialType[j][i]=BUCKINGHAM_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  if((strcasecmp(PotentialName,"BUCKINGHAM_SMOOTHED5")==0)||(strcasecmp(PotentialName,"BUCKINGHAM-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=BUCKINGHAM_SMOOTHED5;
    PotentialType[j][i]=BUCKINGHAM_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=(REAL)0.0;
    PotentialParms[i][j][3]=(REAL)0.0;
  }

  // if(r<p_3) 1e10 else p_0*exp(-p_1*r)-p_2/r^6
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3/k_B [K]  (non-zero for a shifted potential)
  // p_3     [A]
  if(strcasecmp(PotentialName,"BUCKINGHAM2")==0)
  {
    PotentialType[i][j]=BUCKINGHAM2;
    PotentialType[j][i]=BUCKINGHAM2;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4;
    PotentialParms[i][j][3]=arg4;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // if(r<p_3) 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3     [A]
  if((strcasecmp(PotentialName,"BUCKINGHAM2_SMOOTHED3")==0)||(strcasecmp(PotentialName,"BUCKINGHAM2-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=BUCKINGHAM2_SMOOTHED3;
    PotentialType[j][i]=BUCKINGHAM2_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4;
    PotentialParms[i][j][3]=arg4;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // if(r<p-3) 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3     [A]
  if((strcasecmp(PotentialName,"BUCKINGHAM2_SMOOTHED5")==0)||(strcasecmp(PotentialName,"BUCKINGHAM2-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=BUCKINGHAM2_SMOOTHED5;
    PotentialType[j][i]=BUCKINGHAM2_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4;
    PotentialParms[i][j][3]=arg4;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }

  // sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]  if P>=3.02
  // sqrt(p_0^i*p_0^j)*192.27*P^2                    if P<3.02
  // ======================================================================================
  // p_0     [kcal/mol]
  // p_1     [A]
  // p_3     [kcal/mol]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"MM3_VDW")==0)||(strcasecmp(PotentialName,"MM3-VDW")==0))
  {
    PotentialType[i][j]=MM3_VDW;
    PotentialType[j][i]=MM3_VDW;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KCAL_PER_MOL_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KCAL_PER_MOL_TO_ENERGY;
    PotentialParms[j][i][1]=arg2; // input should be the sum of the VDW radii
    PotentialParms[i][j][1]=arg2; 
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
  // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                    if P<3.02
  // ======================================================================================
  // p_0     [kcal/mol]
  // p_1     [A]
  if((strcasecmp(PotentialName,"MM3_VDW_SMOOTHED3")==0)||(strcasecmp(PotentialName,"MM3-VDW-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MM3_VDW_SMOOTHED3;
    PotentialType[j][i]=MM3_VDW_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KCAL_PER_MOL_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KCAL_PER_MOL_TO_ENERGY;
    PotentialParms[j][i][1]=arg2; // input should be the sum of the VDW radii
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
  // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                    if P<3.02
  // ======================================================================================
  // p_0     [kcal/mol]
  // p_1     [A]
  if((strcasecmp(PotentialName,"MM3_VDW_SMOOTHED5")==0)||(strcasecmp(PotentialName,"MM3-VDW-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MM3_VDW_SMOOTHED5;
    PotentialType[j][i]=MM3_VDW_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KCAL_PER_MOL_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KCAL_PER_MOL_TO_ENERGY;
    PotentialParms[j][i][1]=arg2; // input should be the sum of the VDW radii
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  } 
  // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K]
  // p_3     [A^-1]
  // p_4/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"MATSUOKA_CLEMENTI_YOSHIMINE")==0)
  {
    PotentialType[i][j]=MATSUOKA_CLEMENTI_YOSHIMINE;
    PotentialType[j][i]=MATSUOKA_CLEMENTI_YOSHIMINE;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4;
    PotentialParms[i][j][3]=arg4;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K]
  // p_3     [A^-1]
  if((strcasecmp(PotentialName,"MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3")==0)||(strcasecmp(PotentialName,"MATSUOKA-CLEMENTI-YOSHIMINE-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3;
    PotentialType[j][i]=MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4;
    PotentialParms[i][j][3]=arg4;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K]
  // p_3     [A^-1]
  if((strcasecmp(PotentialName,"MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5")==0)||(strcasecmp(PotentialName,"MATSUOKA-CLEMENTI-YOSHIMINE-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5;
    PotentialType[j][i]=MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4;
    PotentialParms[i][j][3]=arg4;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  // p_5/k_B [K A^10]
  // p_6/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"GENERIC")==0)
  {
    PotentialType[i][j]=GENERIC;
    PotentialType[j][i]=GENERIC;
    sscanf(Arguments,"%lf %lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5,&arg6);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][j][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[j][i][5]=arg6*KELVIN_TO_ENERGY;
    PotentialParms[i][j][5]=arg6*KELVIN_TO_ENERGY;
    PotentialParms[j][i][6]=(REAL)0.0;
    PotentialParms[i][j][6]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  // p_5/k_B [K A^10]
  if((strcasecmp(PotentialName,"GENERIC_SMOOTHED3")==0)||(strcasecmp(PotentialName,"GENERIC-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=GENERIC_SMOOTHED3;
    PotentialType[j][i]=GENERIC_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5,&arg6);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][j][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[j][i][5]=arg6*KELVIN_TO_ENERGY;
    PotentialParms[i][j][5]=arg6*KELVIN_TO_ENERGY;
    PotentialParms[j][i][6]=(REAL)0.0;
    PotentialParms[i][j][6]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  // p_5/k_B [K A^10]
  if((strcasecmp(PotentialName,"GENERIC_SMOOTHED5")==0)||(strcasecmp(PotentialName,"GENERIC-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=GENERIC_SMOOTHED5;
    PotentialType[j][i]=GENERIC_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5,&arg6);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][j][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[j][i][5]=arg6*KELVIN_TO_ENERGY;
    PotentialParms[i][j][5]=arg6*KELVIN_TO_ENERGY;
    PotentialParms[j][i][6]=(REAL)0.0;
    PotentialParms[i][j][6]=(REAL)0.0;
  }
  // p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3/k_B [K A^8]
  // p_4/k_B [K A^10]
  // p_5/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"PELLENQ_NICHOLSON")==0)||(strcasecmp(PotentialName,"PELLENQ-NICHOLSON")==0))
  {
    PotentialType[i][j]=PELLENQ_NICHOLSON;
    PotentialType[j][i]=PELLENQ_NICHOLSON;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[j][i][0]=(arg1*KELVIN_TO_ENERGY);
    PotentialParms[i][j][0]=(arg1*KELVIN_TO_ENERGY);
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(arg3*KELVIN_TO_ENERGY);
    PotentialParms[i][j][2]=(arg3*KELVIN_TO_ENERGY);
    PotentialParms[j][i][3]=(arg4*KELVIN_TO_ENERGY);
    PotentialParms[i][j][3]=(arg4*KELVIN_TO_ENERGY);
    PotentialParms[j][i][4]=(arg5*KELVIN_TO_ENERGY);
    PotentialParms[i][j][4]=(arg5*KELVIN_TO_ENERGY);
    PotentialParms[j][i][5]=(REAL)0.0;
    PotentialParms[i][j][5]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3/k_B [K A^8]
  // p_4/k_B [K A^10]
  if((strcasecmp(PotentialName,"PELLENQ_NICHOLSON_SMOOTHED3")==0)||(strcasecmp(PotentialName,"PELLENQ-NICHOLSON-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=PELLENQ_NICHOLSON_SMOOTHED3;
    PotentialType[j][i]=PELLENQ_NICHOLSON_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[j][i][0]=(arg1*KELVIN_TO_ENERGY);
    PotentialParms[i][j][0]=(arg1*KELVIN_TO_ENERGY);
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(arg3*KELVIN_TO_ENERGY);
    PotentialParms[i][j][2]=(arg3*KELVIN_TO_ENERGY);
    PotentialParms[j][i][3]=(arg4*KELVIN_TO_ENERGY);
    PotentialParms[i][j][3]=(arg4*KELVIN_TO_ENERGY);
    PotentialParms[j][i][4]=(arg5*KELVIN_TO_ENERGY);
    PotentialParms[i][j][4]=(arg5*KELVIN_TO_ENERGY);
    PotentialParms[j][i][5]=(REAL)0.0;
    PotentialParms[i][j][5]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^6]
  // p_3/k_B [K A^8]
  // p_4/k_B [K A^10]
  if((strcasecmp(PotentialName,"PELLENQ_NICHOLSON_SMOOTHED5")==0)||(strcasecmp(PotentialName,"PELLENQ-NICHOLSON-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=PELLENQ_NICHOLSON_SMOOTHED5;
    PotentialType[j][i]=PELLENQ_NICHOLSON_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[j][i][0]=(arg1*KELVIN_TO_ENERGY);
    PotentialParms[i][j][0]=(arg1*KELVIN_TO_ENERGY);
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=(arg3*KELVIN_TO_ENERGY);
    PotentialParms[i][j][2]=(arg3*KELVIN_TO_ENERGY);
    PotentialParms[j][i][3]=(arg4*KELVIN_TO_ENERGY);
    PotentialParms[i][j][3]=(arg4*KELVIN_TO_ENERGY);
    PotentialParms[j][i][4]=(arg5*KELVIN_TO_ENERGY);
    PotentialParms[i][j][4]=(arg5*KELVIN_TO_ENERGY);
    PotentialParms[j][i][5]=(REAL)0.0;
    PotentialParms[i][j][5]=(REAL)0.0;
  }
  // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^12]
  // p_5/k_B [K]  (non-zero for a shifted potential)
  if((strcasecmp(PotentialName,"HYDRATED_ION_WATER")==0)||(strcasecmp(PotentialName,"HYDRATED-ION-WATER")==0))
  {
    PotentialType[i][j]=HYDRATED_ION_WATER;
    PotentialType[j][i]=HYDRATED_ION_WATER;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][j][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[j][i][5]=(REAL)0.0;
    PotentialParms[i][j][5]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^12]
  if((strcasecmp(PotentialName,"HYDRATED_ION_WATER_SMOOTHED3")==0)||(strcasecmp(PotentialName,"HYDRATED-ION-WATER-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=HYDRATED_ION_WATER_SMOOTHED3;
    PotentialType[j][i]=HYDRATED_ION_WATER_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][j][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[j][i][5]=(REAL)0.0;
    PotentialParms[i][j][5]=(REAL)0.0;
  }
  // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [A^-1]
  // p_2/k_B [K A^4]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^12]
  if((strcasecmp(PotentialName,"HYDRATED_ION_WATER_SMOOTHED5")==0)||(strcasecmp(PotentialName,"HYDRATED-ION-WATER-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=HYDRATED_ION_WATER_SMOOTHED5;
    PotentialType[j][i]=HYDRATED_ION_WATER_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf",&arg1,&arg2,&arg3);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][j][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[j][i][5]=(REAL)0.0;
    PotentialParms[i][j][5]=(REAL)0.0;
  }
  // p_0*[p_1/r^p_2-p_1/r^p_3]
  // ======================================================================================
  // p_0/k_B [K]
  // p_1/k_B [K A^p_2]
  // p_2     [-]
  // p_3     [-]
  // p_4/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"MIE")==0)
  {
    PotentialType[i][j]=MIE;
    PotentialType[j][i]=MIE;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4;
    PotentialParms[i][j][3]=arg4;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // {p_0*[p_1/r^p_2-p_1/r^p_3]}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1/k_B [K A^p_2]
  // p_2     [-]
  // p_3     [-]
  if((strcasecmp(PotentialName,"MIE_SMOOTHED3")==0)||(strcasecmp(PotentialName,"MIE_SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MIE_SMOOTHED3;
    PotentialType[j][i]=MIE_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[i][j][2]=arg3*KELVIN_TO_ENERGY;
    PotentialParms[j][i][3]=arg4;
    PotentialParms[i][j][3]=arg4;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // {p_0*[p_1/r^p_2-p_1/r^p_3]}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1/k_B [K A^p_2]
  // p_2     [-]
  // p_3     [-]
  if((strcasecmp(PotentialName,"MIE_SMOOTHED5")==0)||(strcasecmp(PotentialName,"MIE_SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=MIE_SMOOTHED5;
    PotentialType[j][i]=MIE_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=arg4;
    PotentialParms[i][j][3]=arg4;
    PotentialParms[j][i][4]=(REAL)0.0;
    PotentialParms[i][j][4]=(REAL)0.0;
  }
  // p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [-]
  // p_2     [A]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  // p_5/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"BORN_HUGGINS_MEYER")==0)
  {
    PotentialType[i][j]=BORN_HUGGINS_MEYER;
    PotentialType[j][i]=BORN_HUGGINS_MEYER;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][j][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[j][i][5]=(REAL)0.0;
    PotentialParms[i][j][5]=(REAL)0.0;
  }
  // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [-]
  // p_2     [A]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  if((strcasecmp(PotentialName,"BORN_HUGGINS_MEYER_SMOOTHED3")==0)||(strcasecmp(PotentialName,"BORN-HUGGINS-MEYER-SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=BORN_HUGGINS_MEYER_SMOOTHED3;
    PotentialType[j][i]=BORN_HUGGINS_MEYER_SMOOTHED3;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][j][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[j][i][5]=(REAL)0.0;
    PotentialParms[i][j][5]=(REAL)0.0;
  }
  // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
  // ======================================================================================
  // p_0/k_B [K]
  // p_1     [-]
  // p_2     [A]
  // p_3/k_B [K A^6]
  // p_4/k_B [K A^8]
  if((strcasecmp(PotentialName,"BORN_HUGGINS_MEYER_SMOOTHED5")==0)||(strcasecmp(PotentialName,"BORN-HUGGINS-MEYER-SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=BORN_HUGGINS_MEYER_SMOOTHED5;
    PotentialType[j][i]=BORN_HUGGINS_MEYER_SMOOTHED5;
    sscanf(Arguments,"%lf %lf %lf %lf %lf",&arg1,&arg2,&arg3,&arg4,&arg5);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2;
    PotentialParms[i][j][1]=arg2;
    PotentialParms[j][i][2]=arg3;
    PotentialParms[i][j][2]=arg3;
    PotentialParms[j][i][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[i][j][3]=arg4*KELVIN_TO_ENERGY;
    PotentialParms[j][i][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[i][j][4]=arg5*KELVIN_TO_ENERGY;
    PotentialParms[j][i][5]=(REAL)0.0;
    PotentialParms[i][j][5]=(REAL)0.0;
  }
  // p_0/r^12-p_1/r^10
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^10]
  // p_2/k_B [K]  (non-zero for a shifted potential)
  if(strcasecmp(PotentialName,"HYDROGEN")==0)
  {
    PotentialType[i][j]=HYDROGEN;
    PotentialType[j][i]=HYDROGEN;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {p_0/r^12-p_1/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^10]
  if((strcasecmp(PotentialName,"HYDROGEN_SMOOTHED3")==0)||(strcasecmp(PotentialName,"HYDROGEN_SMOOTHED3")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=HYDROGEN_SMOOTHED3;
    PotentialType[j][i]=HYDROGEN_SMOOTHED3;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
  // {p_0/r^12-p_1/r^10}*S(r)
  // ======================================================================================
  // p_0/k_B [K A^12]
  // p_1/k_B [K A^10]
  if((strcasecmp(PotentialName,"HYDROGEN_SMOOTHED5")==0)||(strcasecmp(PotentialName,"HYDROGEN_SMOOTHED5")==0))
  {
    TailCorrection[i][j]=TailCorrection[j][i]=FALSE;
    ShiftPotential[i][j]=ShiftPotential[j][i]=FALSE;
    PotentialType[i][j]=HYDROGEN_SMOOTHED5;
    PotentialType[j][i]=HYDROGEN_SMOOTHED5;
    sscanf(Arguments,"%lf %lf",&arg1,&arg2);
    PotentialParms[j][i][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[i][j][0]=arg1*KELVIN_TO_ENERGY;
    PotentialParms[j][i][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[i][j][1]=arg2*KELVIN_TO_ENERGY;
    PotentialParms[j][i][2]=(REAL)0.0;
    PotentialParms[i][j][2]=(REAL)0.0;
  }
}

/*********************************************************************************************************
 * Name       | ReadForceFieldDefinitionsMixingRules                                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Reads the force field definitions using mixing rules ('force_field_mixing_rules.def')    *
 * Parameters | -                                                                                        *
 * Note       | This function is performed before 'ReadForceFieldDefinitions'. This allows the user to   *
 *            | first define generic force fields using 'ReadForceFieldDefinitionsMixingRules' before    *
 *            | overwriting specific rules with 'ReadForceFieldDefinitions'.                             *
 *********************************************************************************************************/

void ReadForceFieldDefinitionsMixingRules(void)
{
  int i,j,k,temp;
  char buffer[256],AtomTypeA[256],PotentialName[256],MixingRuleName[256],Arguments[256];
  int NumberOfInteractions;
  int PatternMatchA;
  FILE *FilePtr;

  // allocate here the force field memory
  PotentialParms=(REAL(**)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS])calloc(NumberOfPseudoAtoms,sizeof(REAL(*)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]));
  PotentialType=(int**)calloc(NumberOfPseudoAtoms,sizeof(int*));
  TailCorrection=(int**)calloc(NumberOfPseudoAtoms,sizeof(int*));
  ShiftPotential=(int**)calloc(NumberOfPseudoAtoms,sizeof(int*));
  MapPseudoAtom=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    PotentialParms[i]=(REAL(*)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS])calloc(NumberOfPseudoAtoms,sizeof(REAL[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]));
    PotentialType[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
    TailCorrection[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
    ShiftPotential[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  }


  GeneralMixingRule=NO_MIXING_RULE;
  for(i=0;i<NumberOfPseudoAtoms;i++)
    for(j=0;j<NumberOfPseudoAtoms;j++)
      PotentialType[i][j]=UNDEFINED_POTENTIAL;

  if(!(FilePtr=fopen("./force_field_mixing_rules.def","r")))
  {
    sprintf(buffer,"%s/share/raspa/forcefield/%s/force_field_mixing_rules.def",RASPA_DIRECTORY,ForceField);
    if(!(FilePtr=fopen(buffer,"r")))
    {
      printf("'force_field_mixing_rules.def' file not found and therefore not used\n");
      return;
    }
  }

  fscanf(FilePtr,"%*[^\n]");fscanf(FilePtr,"%*c");  //skip line
  fscanf(FilePtr,"%s\n",MixingRuleName);
  if(strncasecmp(MixingRuleName,"shifted",strlen("shifted"))==0)
    ShiftPotentials=TRUE;
  else if(strncasecmp(MixingRuleName,"truncated",strlen("truncated"))==0)
    ShiftPotentials=FALSE;
  else
    ShiftPotentials=TRUE;
   
  if(ShiftPotentials)
  {
    for(i=0;i<NumberOfPseudoAtoms;i++)
      for(j=0;j<NumberOfPseudoAtoms;j++)
        ShiftPotential[i][j]=TRUE;
    printf("Shift all potentials\n");
  }
  else
  {
    for(i=0;i<NumberOfPseudoAtoms;i++)
      for(j=0;j<NumberOfPseudoAtoms;j++)
        ShiftPotential[i][j]=FALSE;
  }

  fscanf(FilePtr,"%*[^\n]");fscanf(FilePtr,"%*c");  //skip line
  fscanf(FilePtr,"%s\n",MixingRuleName);
  if(strncasecmp(MixingRuleName,"no",strlen("no"))==0)
    TailCorrections=FALSE;
  else if(strncasecmp(MixingRuleName,"yes",strlen("yes"))==0)
    TailCorrections=TRUE;
  else
    TailCorrections=TRUE;

  if(TailCorrections)
  {
    for(i=0;i<NumberOfPseudoAtoms;i++)
      for(j=0;j<NumberOfPseudoAtoms;j++)
        TailCorrection[i][j]=TRUE;
  }
  else
  {
    for(i=0;i<NumberOfPseudoAtoms;i++)
      for(j=0;j<NumberOfPseudoAtoms;j++)
        TailCorrection[i][j]=FALSE;
  }

  fscanf(FilePtr,"%*[^\n]");fscanf(FilePtr,"%*c");
  fscanf(FilePtr,"%d\n",&temp);     // read NumberOfInteractions
  NumberOfInteractions=temp;

  fscanf(FilePtr,"%*[^\n]");fscanf(FilePtr,"%*c");  //skip line
  for(k=0;k<NumberOfInteractions;k++)
  {
    PatternMatchA=FALSE;
    strcpy(AtomTypeA,"");
    strcpy(PotentialName,"");
    strcpy(Arguments,"");
    fscanf(FilePtr,"%s %s%[^\n]",AtomTypeA,PotentialName,Arguments); fscanf(FilePtr,"%*c");

    // if the atom-type ends with '_' then it applies to every atom-type that starts with that string
    // e.g. 'O_' matches with 'O1', 'O2', etc.
    if(AtomTypeA[strlen(AtomTypeA)-1]=='_')
    {
      AtomTypeA[strlen(AtomTypeA)-1]='\0';
      PatternMatchA=TRUE;
    }

    // NOTE: could possibly also use the chemical-element for pattern-matching
    for(i=0;i<NumberOfPseudoAtoms;i++)
    {
      if(((PatternMatchA==FALSE)&&(strncasecmp(PseudoAtoms[i].Name,AtomTypeA,MAX2(strlen(PseudoAtoms[i].Name),strlen(AtomTypeA)))==0))||
        //((PatternMatchA==TRUE)&&(strncasecmp(PseudoAtoms[i].ChemicalElement,AtomTypeA,MAX2(strlen(PseudoAtoms[i].ChemicalElement),strlen(AtomTypeA)))==0)))
        ((PatternMatchA==TRUE)&&(strncasecmp(PseudoAtoms[i].Name,AtomTypeA,strlen(AtomTypeA))==0)))
      {
        ParseForceFieldSelfParameters(Arguments,i,PotentialName);
      }
    }
  }

  // LORENTZ-BERTHELOT mixing rules assumes the cross-element Lennard-Jones potential
  // are in the middle of the two corresponding pure Lennard-Jones potentials, e.g. arithmetic mean of sigma
  // JORGENSEN mixing rules uses the geometric (geometric mean of sigma and epsilon) mixing rules
  // The arithmetic mean gives marginally better equilibrium distances
  // for van der Waals interactions than the geometric combination rule (Halgren 1992)
  // The 6th-power rule (not available with all forcefields) yields even better results (Waldman and Hagler 1993).
  // With the Ewald method for LJ-potentials (Karasawa & Goddard 1989) (Ewald sums for periodic systems),
  // the geometric mean leads to faster convergence than the arithmetic mean

  fscanf(FilePtr,"%*[^\n]");fscanf(FilePtr,"%*c");  //skip line
  fscanf(FilePtr,"%s\n",buffer);     // read general mixing rule
  if((strcasecmp(buffer,"LORENTZ-BERTHELOT")==0)||
     (strcasecmp(buffer,"LORENTZ_BERTHELOT")==0)||
     (strcasecmp(buffer,"LORENTZBERTHELOT")==0))
  {
    GeneralMixingRule=LORENTZ_BERTHELOT;
    for(i=0;i<NumberOfPseudoAtoms;i++)
      for(j=i+1;j<NumberOfPseudoAtoms;j++)
      {
        if((PotentialType[i][i]==ZERO_POTENTIAL)||(PotentialType[j][j]==ZERO_POTENTIAL))
        {
          PotentialType[i][j]=ZERO_POTENTIAL;
          PotentialType[j][i]=ZERO_POTENTIAL;
        }
        if((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES))
        {
          PotentialType[i][j]=LENNARD_JONES;
          PotentialType[j][i]=LENNARD_JONES;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED3))||
           ((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED3))||
           ((PotentialType[i][i]==LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==LENNARD_JONES)))
        {
          PotentialType[i][j]=LENNARD_JONES_SMOOTHED3;
          PotentialType[j][i]=LENNARD_JONES_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }

        if(((PotentialType[i][i]==LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==LENNARD_JONES))||
           ((PotentialType[i][i]==LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED3)))
        {
          PotentialType[i][j]=LENNARD_JONES_SMOOTHED5;
          PotentialType[j][i]=LENNARD_JONES_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }

        if((PotentialType[i][i]==WCA)&&(PotentialType[j][j]==WCA))
        {
          PotentialType[i][j]=WCA;
          PotentialType[j][i]=WCA;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }

        if((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES)))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3)))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ji=(p_1^j+p_1^i)/2
        }

        if((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2)))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3)))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_1^ji=(p_1^j+p_1^i)/2
        }


        if((PotentialType[i][i]==LENNARD_JONES_SHIFTED_FORCE)&&(PotentialType[j][j]==LENNARD_JONES_SHIFTED_FORCE))
        {
          PotentialType[i][j]=LENNARD_JONES_SHIFTED_FORCE;
          PotentialType[j][i]=LENNARD_JONES_SHIFTED_FORCE;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }

        if((PotentialType[i][i]==LENNARD_JONES_SHIFTED_FORCE2)&&(PotentialType[j][j]==LENNARD_JONES_SHIFTED_FORCE2))
        {
          PotentialType[i][j]=LENNARD_JONES_SHIFTED_FORCE2;
          PotentialType[j][i]=LENNARD_JONES_SHIFTED_FORCE2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }

        if((PotentialType[i][i]==POTENTIAL_12_6)&&(PotentialType[j][j]==POTENTIAL_12_6))
        {
          PotentialType[i][j]=POTENTIAL_12_6;
          PotentialType[j][i]=POTENTIAL_12_6;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED3)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED3))||
           ((PotentialType[i][i]==POTENTIAL_12_6)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED3))||
           ((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED3)&&(PotentialType[j][j]==POTENTIAL_12_6)))
        {
          PotentialType[i][j]=POTENTIAL_12_6_SMOOTHED3;
          PotentialType[j][i]=POTENTIAL_12_6_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED5)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED5))||
           ((PotentialType[i][i]==POTENTIAL_12_6)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED5))||
           ((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED3)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED5))||
           ((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED5)&&(PotentialType[j][j]==POTENTIAL_12_6))||
           ((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED5)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED3)))
        {
          PotentialType[i][j]=POTENTIAL_12_6_SMOOTHED5;
          PotentialType[j][i]=POTENTIAL_12_6_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }

        if((PotentialType[i][i]==MORSE)&&(PotentialType[j][j]==MORSE))
        {
          PotentialType[i][j]=MORSE;
          PotentialType[j][i]=MORSE;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
        }
        if(((PotentialType[i][i]==MORSE_SMOOTHED3)&&(PotentialType[j][j]==MORSE_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE)&&(PotentialType[j][j]==MORSE_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE_SMOOTHED3)&&(PotentialType[j][j]==MORSE)))
        {
          PotentialType[i][j]=MORSE_SMOOTHED3;
          PotentialType[j][i]=MORSE_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
        }
        if(((PotentialType[i][i]==MORSE_SMOOTHED5)&&(PotentialType[j][j]==MORSE_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE)&&(PotentialType[j][j]==MORSE_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE_SMOOTHED3)&&(PotentialType[j][j]==MORSE_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE_SMOOTHED5)&&(PotentialType[j][j]==MORSE))||
           ((PotentialType[i][i]==MORSE_SMOOTHED5)&&(PotentialType[j][j]==MORSE_SMOOTHED3)))
        {
          PotentialType[i][j]=MORSE_SMOOTHED5;
          PotentialType[j][i]=MORSE_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
        }

        if((PotentialType[i][i]==MORSE2)&&(PotentialType[j][j]==MORSE2))
        {
          PotentialType[i][j]=MORSE2;
          PotentialType[j][i]=MORSE2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
        }
        if(((PotentialType[i][i]==MORSE2_SMOOTHED3)&&(PotentialType[j][j]==MORSE2_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE2)&&(PotentialType[j][j]==MORSE2_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE2_SMOOTHED3)&&(PotentialType[j][j]==MORSE2)))
        {
          PotentialType[i][j]=MORSE2_SMOOTHED3;
          PotentialType[j][i]=MORSE2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
        }
        if(((PotentialType[i][i]==MORSE2_SMOOTHED5)&&(PotentialType[j][j]==MORSE2_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE2)&&(PotentialType[j][j]==MORSE2_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE2_SMOOTHED3)&&(PotentialType[j][j]==MORSE2_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE2_SMOOTHED5)&&(PotentialType[j][j]==MORSE2))||
           ((PotentialType[i][i]==MORSE2_SMOOTHED5)&&(PotentialType[j][j]==MORSE2_SMOOTHED3)))
        {
          PotentialType[i][j]=MORSE2_SMOOTHED5;
          PotentialType[j][i]=MORSE2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
        }

        if((PotentialType[i][i]==MORSE3)&&(PotentialType[j][j]==MORSE3))
        {
          PotentialType[i][j]=MORSE3;
          PotentialType[j][i]=MORSE3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==MORSE3_SMOOTHED3)&&(PotentialType[j][j]==MORSE3_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE3)&&(PotentialType[j][j]==MORSE3_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE3_SMOOTHED3)&&(PotentialType[j][j]==MORSE3)))
        {
          PotentialType[i][j]=MORSE3_SMOOTHED3;
          PotentialType[j][i]=MORSE3_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==MORSE3_SMOOTHED5)&&(PotentialType[j][j]==MORSE3_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE3)&&(PotentialType[j][j]==MORSE3_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE3_SMOOTHED3)&&(PotentialType[j][j]==MORSE3_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE3_SMOOTHED5)&&(PotentialType[j][j]==MORSE3))||
           ((PotentialType[i][i]==MORSE3_SMOOTHED5)&&(PotentialType[j][j]==MORSE3_SMOOTHED3)))
        {
          PotentialType[i][j]=MORSE3_SMOOTHED5;
          PotentialType[j][i]=MORSE3_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }


        if((PotentialType[i][i]==CFF_9_6)&&(PotentialType[j][j]==CFF_9_6))
        {
          PotentialType[i][j]=CFF_9_6;
          PotentialType[j][i]=CFF_9_6;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==CFF_9_6_SMOOTHED3)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED3))||
           ((PotentialType[i][i]==CFF_9_6)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED3))||
           ((PotentialType[i][i]==CFF_9_6_SMOOTHED3)&&(PotentialType[j][j]==CFF_9_6)))
        {
          PotentialType[i][j]=CFF_9_6_SMOOTHED3;
          PotentialType[j][i]=CFF_9_6_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==CFF_9_6_SMOOTHED5)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_9_6)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_9_6_SMOOTHED3)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_9_6_SMOOTHED5)&&(PotentialType[j][j]==CFF_9_6))||
           ((PotentialType[i][i]==CFF_9_6_SMOOTHED5)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED3)))
        {
          PotentialType[i][j]=CFF_9_6_SMOOTHED5;
          PotentialType[j][i]=CFF_9_6_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }

        if((PotentialType[i][i]==CFF_EPS_SIGMA)&&(PotentialType[j][j]==CFF_EPS_SIGMA))
        {
          PotentialType[i][j]=CFF_EPS_SIGMA;
          PotentialType[j][i]=CFF_EPS_SIGMA;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED3)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED3))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED3))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED3)&&(PotentialType[j][j]==CFF_EPS_SIGMA)))
        {
          PotentialType[i][j]=CFF_EPS_SIGMA_SMOOTHED3;
          PotentialType[j][i]=CFF_EPS_SIGMA_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED5)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED3)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED5)&&(PotentialType[j][j]==CFF_EPS_SIGMA))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED5)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED3)))
        {
          PotentialType[i][j]=CFF_EPS_SIGMA_SMOOTHED5;
          PotentialType[j][i]=CFF_EPS_SIGMA_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }

        if((PotentialType[i][i]==BUCKINGHAM)&&(PotentialType[j][j]==BUCKINGHAM))
        {
          PotentialType[i][j]=BUCKINGHAM;
          PotentialType[j][i]=BUCKINGHAM;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }
        if(((PotentialType[i][i]==BUCKINGHAM_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED3))||
           ((PotentialType[i][i]==BUCKINGHAM)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED3))||
           ((PotentialType[i][i]==BUCKINGHAM_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM)))
        {
          PotentialType[i][j]=BUCKINGHAM_SMOOTHED3;
          PotentialType[j][i]=BUCKINGHAM_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }
        if(((PotentialType[i][i]==BUCKINGHAM_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM))||
           ((PotentialType[i][i]==BUCKINGHAM_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED3)))
        {
          PotentialType[i][j]=BUCKINGHAM_SMOOTHED5;
          PotentialType[j][i]=BUCKINGHAM_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }

        if((PotentialType[i][i]==BUCKINGHAM2)&&(PotentialType[j][j]==BUCKINGHAM2))
        {
          PotentialType[i][j]=BUCKINGHAM2;
          PotentialType[j][i]=BUCKINGHAM2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ij=0.5*(p_3^i+p_3^j)
          PotentialParms[j][i][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ji=0.5*(p_3^j+p_3^i)
        }
        if(((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED3))||
           ((PotentialType[i][i]==BUCKINGHAM2)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED3))||
           ((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM2)))
        {
          PotentialType[i][j]=BUCKINGHAM2_SMOOTHED3;
          PotentialType[j][i]=BUCKINGHAM2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ij=0.5*(p_3^i+p_3^j)
          PotentialParms[j][i][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ji=0.5*(p_3^j+p_3^i)
        }
        if(((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM2)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM2))||
           ((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED3)))
        {
          PotentialType[i][j]=BUCKINGHAM2_SMOOTHED5;
          PotentialType[j][i]=BUCKINGHAM2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ij=0.5*(p_3^i+p_3^j)
          PotentialParms[j][i][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ji=0.5*(p_3^j+p_3^i)
        }


        if((PotentialType[i][i]==MM3_VDW)&&(PotentialType[j][j]==MM3_VDW))
        {
          PotentialType[i][j]=MM3_VDW;
          PotentialType[j][i]=MM3_VDW;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==MM3_VDW_SMOOTHED3)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED3))||
           ((PotentialType[i][i]==MM3_VDW)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED3))||
           ((PotentialType[i][i]==MM3_VDW_SMOOTHED3)&&(PotentialType[j][j]==MM3_VDW)))
        {
          PotentialType[i][j]=MM3_VDW_SMOOTHED3;
          PotentialType[j][i]=MM3_VDW_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
        if(((PotentialType[i][i]==MM3_VDW_SMOOTHED5)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED5))||
           ((PotentialType[i][i]==MM3_VDW)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED5))||
           ((PotentialType[i][i]==MM3_VDW_SMOOTHED3)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED5))||
           ((PotentialType[i][i]==MM3_VDW_SMOOTHED5)&&(PotentialType[j][j]==MM3_VDW))||
           ((PotentialType[i][i]==MM3_VDW_SMOOTHED5)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED3)))
        {
          PotentialType[i][j]=MM3_VDW_SMOOTHED5;
          PotentialType[j][i]=MM3_VDW_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
        }
      }
  }
  else if(strncasecmp(buffer,"JORGENSEN",strlen("JORGENSEN"))==0)
  {
    GeneralMixingRule=JORGENSEN;
    for(i=0;i<NumberOfPseudoAtoms;i++)
      for(j=i+1;j<NumberOfPseudoAtoms;j++)
      {
        if((PotentialType[i][i]==ZERO_POTENTIAL)||(PotentialType[j][j]==ZERO_POTENTIAL))
        {
          PotentialType[i][j]=ZERO_POTENTIAL;
          PotentialType[j][i]=ZERO_POTENTIAL;
        }
        if((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES))
        {
          PotentialType[i][j]=LENNARD_JONES;
          PotentialType[j][i]=LENNARD_JONES;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED3))||
           ((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED3))||
           ((PotentialType[i][i]==LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==LENNARD_JONES)))
        {
          PotentialType[i][j]=LENNARD_JONES_SMOOTHED3;
          PotentialType[j][i]=LENNARD_JONES_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }

        if(((PotentialType[i][i]==LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==LENNARD_JONES))||
           ((PotentialType[i][i]==LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==LENNARD_JONES_SMOOTHED3)))
        {
          PotentialType[i][j]=LENNARD_JONES_SMOOTHED5;
          PotentialType[j][i]=LENNARD_JONES_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }

        if((PotentialType[i][i]==WCA)&&(PotentialType[j][j]==WCA))
        {
          PotentialType[i][j]=WCA;
          PotentialType[j][i]=WCA;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }

        if((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES)))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3)))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }

        if((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2)))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2))||
           ((PotentialType[i][i]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5)&&(PotentialType[j][j]==FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3)))
        {
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }


        if((PotentialType[i][i]==LENNARD_JONES_SHIFTED_FORCE)&&(PotentialType[j][j]==LENNARD_JONES_SHIFTED_FORCE))
        {
          PotentialType[i][j]=LENNARD_JONES_SHIFTED_FORCE;
          PotentialType[j][i]=LENNARD_JONES_SHIFTED_FORCE;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }

        if((PotentialType[i][i]==LENNARD_JONES_SHIFTED_FORCE2)&&(PotentialType[j][j]==LENNARD_JONES_SHIFTED_FORCE2))
        {
          PotentialType[i][j]=LENNARD_JONES_SHIFTED_FORCE2;
          PotentialType[j][i]=LENNARD_JONES_SHIFTED_FORCE2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }

        if((PotentialType[i][i]==POTENTIAL_12_6)&&(PotentialType[j][j]==POTENTIAL_12_6))
        {
          PotentialType[i][j]=POTENTIAL_12_6;
          PotentialType[j][i]=POTENTIAL_12_6;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED3)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED3))||
           ((PotentialType[i][i]==POTENTIAL_12_6)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED3))||
           ((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED3)&&(PotentialType[j][j]==POTENTIAL_12_6)))
        {
          PotentialType[i][j]=POTENTIAL_12_6_SMOOTHED3;
          PotentialType[j][i]=POTENTIAL_12_6_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED5)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED5))||
           ((PotentialType[i][i]==POTENTIAL_12_6)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED5))||
           ((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED3)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED5))||
           ((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED5)&&(PotentialType[j][j]==POTENTIAL_12_6))||
           ((PotentialType[i][i]==POTENTIAL_12_6_SMOOTHED5)&&(PotentialType[j][j]==POTENTIAL_12_6_SMOOTHED3)))
        {
          PotentialType[i][j]=POTENTIAL_12_6_SMOOTHED5;
          PotentialType[j][i]=POTENTIAL_12_6_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }

        if((PotentialType[i][i]==MORSE)&&(PotentialType[j][j]==MORSE))
        {
          PotentialType[i][j]=MORSE;
          PotentialType[j][i]=MORSE;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }
        if(((PotentialType[i][i]==MORSE_SMOOTHED3)&&(PotentialType[j][j]==MORSE_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE)&&(PotentialType[j][j]==MORSE_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE_SMOOTHED3)&&(PotentialType[j][j]==MORSE)))
        {
          PotentialType[i][j]=MORSE_SMOOTHED3;
          PotentialType[j][i]=MORSE_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }
        if(((PotentialType[i][i]==MORSE_SMOOTHED5)&&(PotentialType[j][j]==MORSE_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE)&&(PotentialType[j][j]==MORSE_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE_SMOOTHED3)&&(PotentialType[j][j]==MORSE_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE_SMOOTHED5)&&(PotentialType[j][j]==MORSE))||
           ((PotentialType[i][i]==MORSE_SMOOTHED5)&&(PotentialType[j][j]==MORSE_SMOOTHED3)))
        {
          PotentialType[i][j]=MORSE_SMOOTHED5;
          PotentialType[j][i]=MORSE_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }

        if((PotentialType[i][i]==MORSE2)&&(PotentialType[j][j]==MORSE2))
        {
          PotentialType[i][j]=MORSE2;
          PotentialType[j][i]=MORSE2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }
        if(((PotentialType[i][i]==MORSE2_SMOOTHED3)&&(PotentialType[j][j]==MORSE2_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE2)&&(PotentialType[j][j]==MORSE2_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE2_SMOOTHED3)&&(PotentialType[j][j]==MORSE2)))
        {
          PotentialType[i][j]=MORSE2_SMOOTHED3;
          PotentialType[j][i]=MORSE2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }
        if(((PotentialType[i][i]==MORSE2_SMOOTHED5)&&(PotentialType[j][j]==MORSE2_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE2)&&(PotentialType[j][j]==MORSE2_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE2_SMOOTHED3)&&(PotentialType[j][j]==MORSE2_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE2_SMOOTHED5)&&(PotentialType[j][j]==MORSE2))||
           ((PotentialType[i][i]==MORSE2_SMOOTHED5)&&(PotentialType[j][j]==MORSE2_SMOOTHED3)))
        {
          PotentialType[i][j]=MORSE2_SMOOTHED5;
          PotentialType[j][i]=MORSE2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }

        if((PotentialType[i][i]==MORSE3)&&(PotentialType[j][j]==MORSE3))
        {
          PotentialType[i][j]=MORSE3;
          PotentialType[j][i]=MORSE3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==MORSE3_SMOOTHED3)&&(PotentialType[j][j]==MORSE3_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE3)&&(PotentialType[j][j]==MORSE3_SMOOTHED3))||
           ((PotentialType[i][i]==MORSE3_SMOOTHED3)&&(PotentialType[j][j]==MORSE3)))
        {
          PotentialType[i][j]=MORSE3_SMOOTHED3;
          PotentialType[j][i]=MORSE3_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==MORSE3_SMOOTHED5)&&(PotentialType[j][j]==MORSE3_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE3)&&(PotentialType[j][j]==MORSE3_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE3_SMOOTHED3)&&(PotentialType[j][j]==MORSE3_SMOOTHED5))||
           ((PotentialType[i][i]==MORSE3_SMOOTHED5)&&(PotentialType[j][j]==MORSE3))||
           ((PotentialType[i][i]==MORSE3_SMOOTHED5)&&(PotentialType[j][j]==MORSE3_SMOOTHED3)))
        {
          PotentialType[i][j]=MORSE3_SMOOTHED5;
          PotentialType[j][i]=MORSE3_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }


        if((PotentialType[i][i]==CFF_9_6)&&(PotentialType[j][j]==CFF_9_6))
        {
          PotentialType[i][j]=CFF_9_6;
          PotentialType[j][i]=CFF_9_6;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==CFF_9_6_SMOOTHED3)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED3))||
           ((PotentialType[i][i]==CFF_9_6)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED3))||
           ((PotentialType[i][i]==CFF_9_6_SMOOTHED3)&&(PotentialType[j][j]==CFF_9_6)))
        {
          PotentialType[i][j]=CFF_9_6_SMOOTHED3;
          PotentialType[j][i]=CFF_9_6_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==CFF_9_6_SMOOTHED5)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_9_6)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_9_6_SMOOTHED3)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_9_6_SMOOTHED5)&&(PotentialType[j][j]==CFF_9_6))||
           ((PotentialType[i][i]==CFF_9_6_SMOOTHED5)&&(PotentialType[j][j]==CFF_9_6_SMOOTHED3)))
        {
          PotentialType[i][j]=CFF_9_6_SMOOTHED5;
          PotentialType[j][i]=CFF_9_6_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }

        if((PotentialType[i][i]==CFF_EPS_SIGMA)&&(PotentialType[j][j]==CFF_EPS_SIGMA))
        {
          PotentialType[i][j]=CFF_EPS_SIGMA;
          PotentialType[j][i]=CFF_EPS_SIGMA;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED3)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED3))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED3))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED3)&&(PotentialType[j][j]==CFF_EPS_SIGMA)))
        {
          PotentialType[i][j]=CFF_EPS_SIGMA_SMOOTHED3;
          PotentialType[j][i]=CFF_EPS_SIGMA_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED5)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED3)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED5))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED5)&&(PotentialType[j][j]==CFF_EPS_SIGMA))||
           ((PotentialType[i][i]==CFF_EPS_SIGMA_SMOOTHED5)&&(PotentialType[j][j]==CFF_EPS_SIGMA_SMOOTHED3)))
        {
          PotentialType[i][j]=CFF_EPS_SIGMA_SMOOTHED5;
          PotentialType[j][i]=CFF_EPS_SIGMA_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }

        if((PotentialType[i][i]==BUCKINGHAM)&&(PotentialType[j][j]==BUCKINGHAM))
        {
          PotentialType[i][j]=BUCKINGHAM;
          PotentialType[j][i]=BUCKINGHAM;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }
        if(((PotentialType[i][i]==BUCKINGHAM_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED3))||
           ((PotentialType[i][i]==BUCKINGHAM)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED3))||
           ((PotentialType[i][i]==BUCKINGHAM_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM)))
        {
          PotentialType[i][j]=BUCKINGHAM_SMOOTHED3;
          PotentialType[j][i]=BUCKINGHAM_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }
        if(((PotentialType[i][i]==BUCKINGHAM_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM))||
           ((PotentialType[i][i]==BUCKINGHAM_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM_SMOOTHED3)))
        {
          PotentialType[i][j]=BUCKINGHAM_SMOOTHED5;
          PotentialType[j][i]=BUCKINGHAM_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
        }

        if((PotentialType[i][i]==BUCKINGHAM2)&&(PotentialType[j][j]==BUCKINGHAM2))
        {
          PotentialType[i][j]=BUCKINGHAM2;
          PotentialType[j][i]=BUCKINGHAM2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ij=sqrt(p_3^i*p_3^j)
          PotentialParms[j][i][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ji=sqrt(p_3^j*p_3^i)
        }
        if(((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED3))||
           ((PotentialType[i][i]==BUCKINGHAM2)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED3))||
           ((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM2)))
        {
          PotentialType[i][j]=BUCKINGHAM2_SMOOTHED3;
          PotentialType[j][i]=BUCKINGHAM2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ij=sqrt(p_3^i*p_3^j)
          PotentialParms[j][i][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ji=sqrt(p_3^j*p_3^i)
        }
        if(((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM2)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED3)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED5))||
           ((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM2))||
           ((PotentialType[i][i]==BUCKINGHAM2_SMOOTHED5)&&(PotentialType[j][j]==BUCKINGHAM2_SMOOTHED3)))
        {
          PotentialType[i][j]=BUCKINGHAM2_SMOOTHED5;
          PotentialType[j][i]=BUCKINGHAM2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ij=sqrt(p_3^i*p_3^j)
          PotentialParms[j][i][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ji=sqrt(p_3^j*p_3^i)
        }


        if((PotentialType[i][i]==MM3_VDW)&&(PotentialType[j][j]==MM3_VDW))
        {
          PotentialType[i][j]=MM3_VDW;
          PotentialType[j][i]=MM3_VDW;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==MM3_VDW_SMOOTHED3)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED3))||
           ((PotentialType[i][i]==MM3_VDW)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED3))||
           ((PotentialType[i][i]==MM3_VDW_SMOOTHED3)&&(PotentialType[j][j]==MM3_VDW)))
        {
          PotentialType[i][j]=MM3_VDW_SMOOTHED3;
          PotentialType[j][i]=MM3_VDW_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
        if(((PotentialType[i][i]==MM3_VDW_SMOOTHED5)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED5))||
           ((PotentialType[i][i]==MM3_VDW)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED5))||
           ((PotentialType[i][i]==MM3_VDW_SMOOTHED3)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED5))||
           ((PotentialType[i][i]==MM3_VDW_SMOOTHED5)&&(PotentialType[j][j]==MM3_VDW))||
           ((PotentialType[i][i]==MM3_VDW_SMOOTHED5)&&(PotentialType[j][j]==MM3_VDW_SMOOTHED3)))
        {
          PotentialType[i][j]=MM3_VDW_SMOOTHED5;
          PotentialType[j][i]=MM3_VDW_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
        }
      }
  }
}

/*********************************************************************************************************
 * Name       | ReadForceFieldDefinitions                                                                *
* ------------------------------------------------------------------------------------------------------ *
 * Function   | Reads the force field definitions using individual rules ('force_field.def')             *
 * Parameters | -                                                                                        *
 * Note       | This function is performed after 'ReadForceFieldDefinitionsMixingRules'. This allows the *
 *            | user to first define generic force fields using 'ReadForceFieldDefinitionsMixingRules'   *
 *            | before overwriting specific rules with 'ReadForceFieldDefinitions'.                      *
 *********************************************************************************************************/

void ReadForceFieldDefinitions(void)
{
  int i,j,k,temp; 
  char buffer[256],AtomTypeA[256],AtomTypeB[256],PotentialName[256],Arguments[256];
  char MixingRuleName[256],Tail[256];
  int NumberOfInteractions;
  FILE *FilePtr;
  int PatternMatchA,PatternMatchB;

  if(!(FilePtr=fopen("./force_field.def","r")))
  {
    sprintf(buffer,"%s/share/raspa/forcefield/%s/force_field.def",RASPA_DIRECTORY,ForceField);
    if(!(FilePtr=fopen(buffer,"r")))
    {
      printf("'force_field.def' file not found and therefore not used\n");
      return;
    }
  }

  // read the truncated vs shifted to overwrite
  fscanf(FilePtr,"%*[^\n]");fscanf(FilePtr,"%*c");
  fscanf(FilePtr,"%d\n",&temp);     // read NumberOfInteractions
  NumberOfInteractions=temp;
  if(NumberOfInteractions>0)
  {
    fscanf(FilePtr,"%*[^\n]");fscanf(FilePtr,"%*c");  //skip line
      for(k=0;k<NumberOfInteractions;k++)
    {
      fscanf(FilePtr,"%s %s %s %s\n",AtomTypeA,AtomTypeB,MixingRuleName,Tail);
      i=ReturnPseudoAtomNumber(AtomTypeA);
      j=ReturnPseudoAtomNumber(AtomTypeB);

      if(strncasecmp(MixingRuleName,"truncated",strlen("truncated"))==0)
      {
        ShiftPotential[i][j]=FALSE;
        ShiftPotential[j][i]=FALSE;
      }
      else if(strncasecmp(MixingRuleName,"shifted",strlen("shifted"))==0)
      {
        ShiftPotential[i][j]=TRUE;
        ShiftPotential[j][i]=TRUE;
      }
      else
      {
        ShiftPotential[i][j]=TRUE;
        ShiftPotential[j][i]=TRUE;
      }

      if(strncasecmp(Tail,"no",strlen("no"))==0)
      {
        TailCorrection[i][j]=FALSE;
        TailCorrection[j][i]=FALSE;
      }
      else if(strncasecmp(Tail,"yes",strlen("yes"))==0)
      {
        TailCorrection[i][j]=TRUE;
        TailCorrection[j][i]=TRUE;
      }
      else
      {
        TailCorrection[i][j]=TRUE;
        TailCorrection[j][i]=TRUE;
      }
    }
  }

  fscanf(FilePtr,"%*[^\n]");fscanf(FilePtr,"%*c");
  fscanf(FilePtr,"%d\n",&temp);     // read NumberOfInteractions
  NumberOfInteractions=temp;
  IndividualInteractions=temp;
  if(NumberOfInteractions>0)
    fscanf(FilePtr,"%*[^\n]");fscanf(FilePtr,"%*c");  //skip line

  for(k=0;k<NumberOfInteractions;k++)
  {
    PatternMatchA=PatternMatchB=FALSE;
    strcpy(AtomTypeA,"");
    strcpy(AtomTypeB,"");
    strcpy(PotentialName,"");
    strcpy(Arguments,"");
    fscanf(FilePtr,"%s %s %s%[^\n]",AtomTypeA,AtomTypeB,PotentialName,Arguments); fscanf(FilePtr,"%*c");

    // if the atom-type ends with '_' then it applies to every atom-type that starts with that string
    // e.g. 'O_' matches with 'O1', 'O2', etc.
    if(AtomTypeA[strlen(AtomTypeA)-1]=='_')
    {
      AtomTypeA[strlen(AtomTypeA)-1]='\0';
      PatternMatchA=TRUE;
    }

    if(AtomTypeB[strlen(AtomTypeB)-1]=='_')
    {
      AtomTypeB[strlen(AtomTypeB)-1]='\0';
      PatternMatchB=TRUE;
    }

    for(i=0;i<NumberOfPseudoAtoms;i++)
    {
      if(((PatternMatchA==FALSE)&&(strncasecmp(PseudoAtoms[i].Name,AtomTypeA,MAX2(strlen(PseudoAtoms[i].Name),strlen(AtomTypeA)))==0))||
        //((PatternMatchA==TRUE)&&(strncasecmp(PseudoAtoms[i].ChemicalElement,AtomTypeA,MAX2(strlen(PseudoAtoms[i].ChemicalElement),strlen(AtomTypeA)))==0)))
        ((PatternMatchA==TRUE)&&(strncasecmp(PseudoAtoms[i].Name,AtomTypeA,strlen(AtomTypeA))==0)))
      {
        for(j=0;j<NumberOfPseudoAtoms;j++)
        {
          if(((PatternMatchB==FALSE)&&(strncasecmp(PseudoAtoms[j].Name,AtomTypeB,MAX2(strlen(PseudoAtoms[j].Name),strlen(AtomTypeB)))==0))||
            //((PatternMatchB==TRUE)&&(strncasecmp(PseudoAtoms[j].ChemicalElement,AtomTypeB,MAX2(strlen(PseudoAtoms[j].ChemicalElement),strlen(AtomTypeB)))==0)))
            ((PatternMatchB==TRUE)&&(strncasecmp(PseudoAtoms[j].Name,AtomTypeB,strlen(AtomTypeB))==0)))
          {
            ParseForceFieldBinaryParameters(Arguments,i,j,PotentialName);
          }
        }
      }
    }
  }


  fscanf(FilePtr,"%[^\n]",Arguments);fscanf(FilePtr,"%*c");
  fscanf(FilePtr,"%d\n",&temp);     // read NumberOfInteractions
  NumberOfInteractions=temp;
  IndividualMixingRules=temp;
  fscanf(FilePtr,"%*[^\n]");fscanf(FilePtr,"%*c");  //skip line

  for(k=0;k<NumberOfInteractions;k++)
  {
    fscanf(FilePtr,"%s %s %s\n",AtomTypeA,AtomTypeB,MixingRuleName);

    i=ReturnPseudoAtomNumber(AtomTypeA);
    j=ReturnPseudoAtomNumber(AtomTypeB);

    if((strcasecmp(MixingRuleName,"LORENTZ_BERTHELOT")==0)||
       (strcasecmp(MixingRuleName,"LORENTZ-BERTHELOT")==0))
    {
      switch(PotentialType[i][i])
      {
        case ZERO_POTENTIAL:
          PotentialType[i][j]=ZERO_POTENTIAL;
          PotentialType[j][i]=ZERO_POTENTIAL;
          break;
        case LENNARD_JONES:
          PotentialType[i][j]=LENNARD_JONES;
          PotentialType[j][i]=LENNARD_JONES;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case LENNARD_JONES_SMOOTHED3:
          PotentialType[i][j]=LENNARD_JONES_SMOOTHED3;
          PotentialType[j][i]=LENNARD_JONES_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case LENNARD_JONES_SMOOTHED5:
          PotentialType[i][j]=LENNARD_JONES_SMOOTHED5;
          PotentialType[j][i]=LENNARD_JONES_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case WCA:
          PotentialType[i][j]=WCA;
          PotentialType[j][i]=WCA;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES2:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case LENNARD_JONES_SHIFTED_FORCE:
          PotentialType[i][j]=LENNARD_JONES_SHIFTED_FORCE;
          PotentialType[j][i]=LENNARD_JONES_SHIFTED_FORCE;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case LENNARD_JONES_SHIFTED_FORCE2:
          PotentialType[i][j]=LENNARD_JONES_SHIFTED_FORCE2;
          PotentialType[j][i]=LENNARD_JONES_SHIFTED_FORCE2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case POTENTIAL_12_6:
          PotentialType[i][j]=POTENTIAL_12_6;
          PotentialType[j][i]=POTENTIAL_12_6;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case POTENTIAL_12_6_SMOOTHED3:
          PotentialType[i][j]=POTENTIAL_12_6_SMOOTHED3;
          PotentialType[j][i]=POTENTIAL_12_6_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case POTENTIAL_12_6_SMOOTHED5:
          PotentialType[i][j]=POTENTIAL_12_6_SMOOTHED5;
          PotentialType[j][i]=POTENTIAL_12_6_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case MORSE:
          PotentialType[i][j]=MORSE;
          PotentialType[j][i]=MORSE;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case MORSE_SMOOTHED3:
          PotentialType[i][j]=MORSE_SMOOTHED3;
          PotentialType[j][i]=MORSE_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case MORSE_SMOOTHED5:
          PotentialType[i][j]=MORSE_SMOOTHED5;
          PotentialType[j][i]=MORSE_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case MORSE2:
          PotentialType[i][j]=MORSE2;
          PotentialType[j][i]=MORSE2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case MORSE2_SMOOTHED3:
          PotentialType[i][j]=MORSE2_SMOOTHED3;
          PotentialType[j][i]=MORSE2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case MORSE2_SMOOTHED5:
          PotentialType[i][j]=MORSE2_SMOOTHED5;
          PotentialType[j][i]=MORSE2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          PotentialParms[i][j][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ij=(p_2^i+p_2^j)/2
          PotentialParms[j][i][2]=0.5*(PotentialParms[i][i][2]+PotentialParms[j][j][2]); // p_2^ji=(p_2^j+p_2^i)/2
          break;
        case MORSE3:
          PotentialType[i][j]=MORSE3;
          PotentialType[j][i]=MORSE3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case MORSE3_SMOOTHED3:
          PotentialType[i][j]=MORSE3_SMOOTHED3;
          PotentialType[j][i]=MORSE3_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case MORSE3_SMOOTHED5:
          PotentialType[i][j]=MORSE3_SMOOTHED5;
          PotentialType[j][i]=MORSE3_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case CFF_9_6:
          PotentialType[i][j]=CFF_9_6;
          PotentialType[j][i]=CFF_9_6;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case CFF_9_6_SMOOTHED3:
          PotentialType[i][j]=CFF_9_6_SMOOTHED3;
          PotentialType[j][i]=CFF_9_6_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case CFF_9_6_SMOOTHED5:
          PotentialType[i][j]=CFF_9_6_SMOOTHED5;
          PotentialType[j][i]=CFF_9_6_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case CFF_EPS_SIGMA:
          PotentialType[i][j]=CFF_EPS_SIGMA;
          PotentialType[j][i]=CFF_EPS_SIGMA;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case CFF_EPS_SIGMA_SMOOTHED3:
          PotentialType[i][j]=CFF_EPS_SIGMA_SMOOTHED3;
          PotentialType[j][i]=CFF_EPS_SIGMA_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case CFF_EPS_SIGMA_SMOOTHED5:
          PotentialType[i][j]=CFF_EPS_SIGMA_SMOOTHED5;
          PotentialType[j][i]=CFF_EPS_SIGMA_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
          break;
        case BUCKINGHAM:
          PotentialType[i][j]=BUCKINGHAM;
          PotentialType[j][i]=BUCKINGHAM;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case BUCKINGHAM_SMOOTHED3:
          PotentialType[i][j]=BUCKINGHAM_SMOOTHED3;
          PotentialType[j][i]=BUCKINGHAM_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case BUCKINGHAM_SMOOTHED5:
          PotentialType[i][j]=BUCKINGHAM_SMOOTHED5;
          PotentialType[j][i]=BUCKINGHAM_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case BUCKINGHAM2:
          PotentialType[i][j]=BUCKINGHAM2;
          PotentialType[j][i]=BUCKINGHAM2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ij=0.5*(p_3^i+p_3^j)
          PotentialParms[j][i][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ji=0.5*(p_3^j+p_3^i)
          break;
        case BUCKINGHAM2_SMOOTHED3:
          PotentialType[i][j]=BUCKINGHAM2_SMOOTHED3;
          PotentialType[j][i]=BUCKINGHAM2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ij=0.5*(p_3^i+p_3^j)
          PotentialParms[j][i][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ji=0.5*(p_3^j+p_3^i)
          break;
        case BUCKINGHAM2_SMOOTHED5:
          PotentialType[i][j]=BUCKINGHAM2_SMOOTHED5;
          PotentialType[j][i]=BUCKINGHAM2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ij=0.5*(p_3^i+p_3^j)
          PotentialParms[j][i][3]=0.5*(PotentialParms[i][i][3]+PotentialParms[j][j][3]); // p_3^ji=0.5*(p_3^j+p_3^i)
          break;
        case MM3_VDW:
          PotentialType[i][j]=MM3_VDW;
          PotentialType[j][i]=MM3_VDW;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          break;
        case MM3_VDW_SMOOTHED3:
          PotentialType[i][j]=MM3_VDW_SMOOTHED3;
          PotentialType[j][i]=MM3_VDW_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          break;
        case MM3_VDW_SMOOTHED5:
          PotentialType[i][j]=MM3_VDW_SMOOTHED5;
          PotentialType[j][i]=MM3_VDW_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=0.5*(p_1^i+p_1^j)
          PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=0.5*(p_1^j+p_1^i)
          break;
        default:
          printf("No mixing rule defined for these type of interaction %d-%d\n",i,j);
          exit(0);
          break;
      }
    }

    if(strcasecmp(MixingRuleName,"JORGENSEN")==0)
    {
      switch(PotentialType[i][i])
      {
       case ZERO_POTENTIAL:
          PotentialType[i][j]=ZERO_POTENTIAL;
          PotentialType[j][i]=ZERO_POTENTIAL;
          break;
        case LENNARD_JONES:
          PotentialType[i][j]=LENNARD_JONES;
          PotentialType[j][i]=LENNARD_JONES;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case LENNARD_JONES_SMOOTHED3:
          PotentialType[i][j]=LENNARD_JONES_SMOOTHED3;
          PotentialType[j][i]=LENNARD_JONES_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case LENNARD_JONES_SMOOTHED5:
          PotentialType[i][j]=LENNARD_JONES_SMOOTHED5;
          PotentialType[j][i]=LENNARD_JONES_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case WCA:
          PotentialType[i][j]=WCA;
          PotentialType[j][i]=WCA;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES2:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5:
          PotentialType[i][j]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5;
          PotentialType[j][i]=FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case LENNARD_JONES_SHIFTED_FORCE:
          PotentialType[i][j]=LENNARD_JONES_SHIFTED_FORCE;
          PotentialType[j][i]=LENNARD_JONES_SHIFTED_FORCE;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case LENNARD_JONES_SHIFTED_FORCE2:
          PotentialType[i][j]=LENNARD_JONES_SHIFTED_FORCE2;
          PotentialType[j][i]=LENNARD_JONES_SHIFTED_FORCE2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case POTENTIAL_12_6:
          PotentialType[i][j]=POTENTIAL_12_6;
          PotentialType[j][i]=POTENTIAL_12_6;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case POTENTIAL_12_6_SMOOTHED3:
          PotentialType[i][j]=POTENTIAL_12_6_SMOOTHED3;
          PotentialType[j][i]=POTENTIAL_12_6_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case POTENTIAL_12_6_SMOOTHED5:
          PotentialType[i][j]=POTENTIAL_12_6_SMOOTHED5;
          PotentialType[j][i]=POTENTIAL_12_6_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case MORSE:
          PotentialType[i][j]=MORSE;
          PotentialType[j][i]=MORSE;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case MORSE_SMOOTHED3:
          PotentialType[i][j]=MORSE_SMOOTHED3;
          PotentialType[j][i]=MORSE_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case MORSE_SMOOTHED5:
          PotentialType[i][j]=MORSE_SMOOTHED5;
          PotentialType[j][i]=MORSE_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case MORSE2:
          PotentialType[i][j]=MORSE2;
          PotentialType[j][i]=MORSE2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case MORSE2_SMOOTHED3:
          PotentialType[i][j]=MORSE2_SMOOTHED3;
          PotentialType[j][i]=MORSE2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case MORSE2_SMOOTHED5:
          PotentialType[i][j]=MORSE2_SMOOTHED5;
          PotentialType[j][i]=MORSE2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case MORSE3:
          PotentialType[i][j]=MORSE3;
          PotentialType[j][i]=MORSE3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case MORSE3_SMOOTHED3:
          PotentialType[i][j]=MORSE3_SMOOTHED3;
          PotentialType[j][i]=MORSE3_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case MORSE3_SMOOTHED5:
          PotentialType[i][j]=MORSE3_SMOOTHED5;
          PotentialType[j][i]=MORSE3_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case CFF_9_6:
          PotentialType[i][j]=CFF_9_6;
          PotentialType[j][i]=CFF_9_6;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case CFF_9_6_SMOOTHED3:
          PotentialType[i][j]=CFF_9_6_SMOOTHED3;
          PotentialType[j][i]=CFF_9_6_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case CFF_9_6_SMOOTHED5:
          PotentialType[i][j]=CFF_9_6_SMOOTHED5;
          PotentialType[j][i]=CFF_9_6_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case CFF_EPS_SIGMA:
          PotentialType[i][j]=CFF_EPS_SIGMA;
          PotentialType[j][i]=CFF_EPS_SIGMA;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case CFF_EPS_SIGMA_SMOOTHED3:
          PotentialType[i][j]=CFF_EPS_SIGMA_SMOOTHED3;
          PotentialType[j][i]=CFF_EPS_SIGMA_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case CFF_EPS_SIGMA_SMOOTHED5:
          PotentialType[i][j]=CFF_EPS_SIGMA_SMOOTHED5;
          PotentialType[j][i]=CFF_EPS_SIGMA_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case BUCKINGHAM:
          PotentialType[i][j]=BUCKINGHAM;
          PotentialType[j][i]=BUCKINGHAM;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case BUCKINGHAM_SMOOTHED3:
          PotentialType[i][j]=BUCKINGHAM_SMOOTHED3;
          PotentialType[j][i]=BUCKINGHAM_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case BUCKINGHAM_SMOOTHED5:
          PotentialType[i][j]=BUCKINGHAM_SMOOTHED5;
          PotentialType[j][i]=BUCKINGHAM_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          break;
        case BUCKINGHAM2:
          PotentialType[i][j]=BUCKINGHAM2;
          PotentialType[j][i]=BUCKINGHAM2;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ij=sqrt(p_3^i*p_3^j)
          PotentialParms[j][i][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ji=sqrt(p_3^j*p_3^i)
          break;
        case BUCKINGHAM2_SMOOTHED3:
          PotentialType[i][j]=BUCKINGHAM2_SMOOTHED3;
          PotentialType[j][i]=BUCKINGHAM2_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ij=sqrt(p_3^i*p_3^j)
          PotentialParms[j][i][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ji=sqrt(p_3^j*p_3^i)
          break;
        case BUCKINGHAM2_SMOOTHED5:
          PotentialType[i][j]=BUCKINGHAM2_SMOOTHED5;
          PotentialType[j][i]=BUCKINGHAM2_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          PotentialParms[i][j][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ij=sqrt(p_2^i*p_2^j)
          PotentialParms[j][i][2]=sqrt(PotentialParms[i][i][2]*PotentialParms[j][j][2]); // p_2^ji=sqrt(p_2^j*p_2^i)
          PotentialParms[i][j][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ij=sqrt(p_3^i*p_3^j)
          PotentialParms[j][i][3]=sqrt(PotentialParms[i][i][3]*PotentialParms[j][j][3]); // p_3^ji=sqrt(p_3^j*p_3^i)
          break;
        case MM3_VDW:
          PotentialType[i][j]=MM3_VDW;
          PotentialType[j][i]=MM3_VDW;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case MM3_VDW_SMOOTHED3:
          PotentialType[i][j]=MM3_VDW_SMOOTHED3;
          PotentialType[j][i]=MM3_VDW_SMOOTHED3;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        case MM3_VDW_SMOOTHED5:
          PotentialType[i][j]=MM3_VDW_SMOOTHED5;
          PotentialType[j][i]=MM3_VDW_SMOOTHED5;
          PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
          PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
          PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
          PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
          break;
        default:
          printf("No mixing rule defined for this type of interaction\n");
          exit(0);
          break;
      }
    }
  }
}

void ComputePotentialShifts(void)
{
  int i,j;
  REAL arg3,arg4,arg5,arg6,arg7;

  for(i=0;i<NumberOfPseudoAtoms;i++)
    for(j=i;j<NumberOfPseudoAtoms;j++)
    {
      switch(PotentialType[i][j])
      {
        case UNDEFINED_POTENTIAL:
        case ZERO_POTENTIAL:
          TailCorrection[i][j]=FALSE;
          break;
        case LENNARD_JONES:
          if(ShiftPotential[i][j])
            arg3=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg3=(REAL)0.0;
          PotentialParms[j][i][2]=arg3;
          PotentialParms[i][j][2]=arg3;
          break;
        case LENNARD_JONES_SMOOTHED3:
        case LENNARD_JONES_SMOOTHED5:
          PotentialParms[j][i][2]=0.0;
          PotentialParms[i][j][2]=0.0;
          break;
        case WCA:
          if(ShiftPotential[i][j])
            arg3=PotentialValue(i,j,SQR(pow(2.0,1.0/6.0)*PotentialParms[i][j][1]));
          else
            arg3=(REAL)0.0;
          PotentialParms[j][i][2]=arg3;
          PotentialParms[i][j][2]=arg3;
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES:
          if(ShiftPotential[i][j])
            arg4=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg4=(REAL)0.0;
          PotentialParms[j][i][3]=arg4;
          PotentialParms[i][j][3]=arg4;
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3:
        case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5:
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES2:
          if(ShiftPotential[i][j])
            arg4=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg4=(REAL)0.0;
          PotentialParms[j][i][3]=arg4;
          PotentialParms[i][j][3]=arg4;
          break;
        case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3:
        case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5:
          break;
        case LENNARD_JONES_SHIFTED_FORCE:
        case LENNARD_JONES_SHIFTED_FORCE2:
          break;
        case POTENTIAL_12_6:
          if(ShiftPotential[i][j])
            arg3=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg3=(REAL)0.0;
          PotentialParms[j][i][2]=arg3;
          PotentialParms[i][j][2]=arg3;
          break;
        case POTENTIAL_12_6_SMOOTHED3:
        case POTENTIAL_12_6_SMOOTHED5:
          break;
        case POTENTIAL_12_6_2_0:
          if(ShiftPotential[i][j])
            arg5=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg5=(REAL)0.0;
          PotentialParms[j][i][4]=arg5;
          PotentialParms[i][j][4]=arg5;
          break;
        case POTENTIAL_12_6_2_0_SMOOTHED3:
        case POTENTIAL_12_6_2_0_SMOOTHED5:
          break;
        case MORSE:
          if(ShiftPotential[i][j])
            arg4=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg4=(REAL)0.0;
          PotentialParms[j][i][3]=arg4;
          PotentialParms[i][j][3]=arg4;
          break;
        case MORSE_SMOOTHED3:
        case MORSE_SMOOTHED5:
          break;
        case MORSE2:
          if(ShiftPotential[i][j])
            arg4=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg4=(REAL)0.0;
          PotentialParms[j][i][3]=arg4;
          PotentialParms[i][j][3]=arg4;
          break;
        case MORSE2_SMOOTHED3:
        case MORSE2_SMOOTHED5:
          break;
        case MORSE3:
          if(ShiftPotential[i][j])
            arg3=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg3=(REAL)0.0;
          PotentialParms[j][i][2]=arg3;
          PotentialParms[i][j][2]=arg3;
          break;
        case MORSE3_SMOOTHED3:
        case MORSE3_SMOOTHED5:
          break;
        case CFF_9_6:
          if(ShiftPotential[i][j])
            arg3=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg3=(REAL)0.0;
          PotentialParms[j][i][2]=arg3;
          PotentialParms[i][j][2]=arg3;
          break;
        case CFF_9_6_SMOOTHED3:
        case CFF_9_6_SMOOTHED5:
          break;
        case CFF_EPS_SIGMA:
          if(ShiftPotential[i][j])
            arg3=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg3=(REAL)0.0;
          PotentialParms[j][i][2]=arg3;
          PotentialParms[i][j][2]=arg3;
          break;
        case CFF_EPS_SIGMA_SMOOTHED3:
        case CFF_EPS_SIGMA_SMOOTHED5:
          break;
        case BUCKINGHAM:
          if(ShiftPotential[i][j])
            arg4=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg4=(REAL)0.0;
          PotentialParms[j][i][3]=arg4;
          PotentialParms[i][j][3]=arg4;
          break;
        case BUCKINGHAM_SMOOTHED3:
        case BUCKINGHAM_SMOOTHED5:
          break;
        case BUCKINGHAM2:
          if(ShiftPotential[i][j])
            arg5=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg5=(REAL)0.0;
          PotentialParms[j][i][4]=arg5;
          PotentialParms[i][j][4]=arg5;
          break;
        case BUCKINGHAM2_SMOOTHED3:
        case BUCKINGHAM2_SMOOTHED5:
          break;
        case MM3_VDW:
        case MM3_HYDROGEN_VDW:
          if(ShiftPotential[i][j])
            arg3=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg3=(REAL)0.0;
          PotentialParms[j][i][2]=arg3;
          PotentialParms[i][j][2]=arg3;
          break;
        case MM3_VDW_SMOOTHED3:
        case MM3_VDW_SMOOTHED5:
          break;
        case MATSUOKA_CLEMENTI_YOSHIMINE:
          if(ShiftPotential[i][j])
            arg5=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg5=(REAL)0.0;
          PotentialParms[j][i][4]=arg5;
          PotentialParms[i][j][4]=arg5;
          break;
        case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3:
        case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5:
          break;
        case GENERIC:
          if(ShiftPotential[i][j])
            arg7=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg7=(REAL)0.0;
          PotentialParms[j][i][6]=arg7;
          PotentialParms[i][j][6]=arg7;
          break;
        case GENERIC_SMOOTHED3:
        case GENERIC_SMOOTHED5:
          break;
        case PELLENQ_NICHOLSON:
          if(ShiftPotential[i][j])
            arg6=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg6=(REAL)0.0;
          PotentialParms[j][i][5]=arg6;
          PotentialParms[i][j][5]=arg6;
          break;
        case PELLENQ_NICHOLSON_SMOOTHED3:
        case PELLENQ_NICHOLSON_SMOOTHED5:
          break;
        case HYDRATED_ION_WATER:
          if(ShiftPotential[i][j])
            arg6=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg6=(REAL)0.0;
          PotentialParms[j][i][5]=arg6;
          PotentialParms[i][j][5]=arg6;
          break;
        case HYDRATED_ION_WATER_SMOOTHED3:
        case HYDRATED_ION_WATER_SMOOTHED5:
          break;
        case MIE:
          if(ShiftPotential[i][j])
            arg5=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg5=(REAL)0.0;
          PotentialParms[j][i][4]=arg5;
          PotentialParms[i][j][4]=arg5;
          break;
        case MIE_SMOOTHED3:
        case MIE_SMOOTHED5:
          break;
        case BORN_HUGGINS_MEYER:
          if(ShiftPotential[i][j])
            arg6=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg6=(REAL)0.0;
          PotentialParms[j][i][5]=arg6;
          PotentialParms[i][j][5]=arg6;
          break;
        case BORN_HUGGINS_MEYER_SMOOTHED3:
        case BORN_HUGGINS_MEYER_SMOOTHED5:
          break;
        case HYDROGEN:
          if(ShiftPotential[i][j])
            arg3=PotentialValue(i,j,CutOffVDWSquared);
          else
            arg3=(REAL)0.0;
          PotentialParms[j][i][2]=arg3;
          PotentialParms[i][j][2]=arg3;
          break;
        case HYDROGEN_SMOOTHED3:
        case HYDROGEN_SMOOTHED5:
          break;
        default:
          break;
      }
    }
}

/*******************************************************************************************
 * Name       | ConvertFromABCtoXYZ                                                        *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the Cartesian position converted from a fractional position        *
 * Parameters | VECTOR t: the fractional position                                          *
 *******************************************************************************************/

VECTOR ConvertFromABCtoXYZ(VECTOR t)
{
  VECTOR dr;

  dr.x=Box[CurrentSystem].ax*t.x+Box[CurrentSystem].bx*t.y+Box[CurrentSystem].cx*t.z;
  dr.y=Box[CurrentSystem].ay*t.x+Box[CurrentSystem].by*t.y+Box[CurrentSystem].cy*t.z;
  dr.z=Box[CurrentSystem].az*t.x+Box[CurrentSystem].bz*t.y+Box[CurrentSystem].cz*t.z;
  return dr;
}


/*******************************************************************************************
 * Name       | ConvertFromXYZtoABC                                                        *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the fractional position converted from a Cartesian position        *
 * Parameters | VECTOR t: the Cartesian position                                           *
 *******************************************************************************************/

VECTOR ConvertFromXYZtoABC(VECTOR t)
{
  VECTOR s;

  s.x=InverseBox[CurrentSystem].ax*t.x+InverseBox[CurrentSystem].bx*t.y+InverseBox[CurrentSystem].cx*t.z;
  s.y=InverseBox[CurrentSystem].ay*t.x+InverseBox[CurrentSystem].by*t.y+InverseBox[CurrentSystem].cy*t.z;
  s.z=InverseBox[CurrentSystem].az*t.x+InverseBox[CurrentSystem].bz*t.y+InverseBox[CurrentSystem].cz*t.z;

  return s;
}


/*******************************************************************************************
 * Name       | ConvertFromABCtoXYZUnitCell                                                *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the Cartesian position converted from a fractional position        *
 * Parameters | VECTOR t: the fractional position                                          *
 *******************************************************************************************/

VECTOR ConvertFromABCtoXYZUnitCell(VECTOR t)
{
  VECTOR dr;

  dr.x=UnitCellBox[CurrentSystem].ax*t.x+UnitCellBox[CurrentSystem].bx*t.y+UnitCellBox[CurrentSystem].cx*t.z;
  dr.y=UnitCellBox[CurrentSystem].ay*t.x+UnitCellBox[CurrentSystem].by*t.y+UnitCellBox[CurrentSystem].cy*t.z;
  dr.z=UnitCellBox[CurrentSystem].az*t.x+UnitCellBox[CurrentSystem].bz*t.y+UnitCellBox[CurrentSystem].cz*t.z;
  return dr;
}


/*******************************************************************************************
 * Name       | ConvertFromXYZtoABCUnitCell                                                *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the fractional position converted from a Cartesian position        *
 * Parameters | VECTOR t: the Cartesian position                                           *
 *******************************************************************************************/

VECTOR ConvertFromXYZtoABCUnitCell(VECTOR t)
{
  VECTOR s;

  s.x=InverseUnitCellBox[CurrentSystem].ax*t.x+InverseUnitCellBox[CurrentSystem].bx*t.y+InverseUnitCellBox[CurrentSystem].cx*t.z;
  s.y=InverseUnitCellBox[CurrentSystem].ay*t.x+InverseUnitCellBox[CurrentSystem].by*t.y+InverseUnitCellBox[CurrentSystem].cy*t.z;
  s.z=InverseUnitCellBox[CurrentSystem].az*t.x+InverseUnitCellBox[CurrentSystem].bz*t.y+InverseUnitCellBox[CurrentSystem].cz*t.z;

  return s;
}

/*******************************************************************************************
 * Name       | ApplyBoundaryCondition                                                     *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the Cartesian distance vector corrected for periodic boundaries    *
 * Parameters | VECTOR dr: the Cartesian distance vector                                   *
 * Note       | RASPA does not restrict the positions to be within the simulation cell     *
 *******************************************************************************************/
/*
VECTOR ApplyBoundaryCondition(VECTOR dr)
{
  VECTOR s,t;

  switch(BoundaryCondition[CurrentSystem])
  {
    case FINITE:
      break;
    case RECTANGULAR:
    case CUBIC:
      dr.x-=Box[CurrentSystem].ax*(REAL)NINT(dr.x*InverseBox[CurrentSystem].ax);
      dr.y-=Box[CurrentSystem].by*(REAL)NINT(dr.y*InverseBox[CurrentSystem].by);
      dr.z-=Box[CurrentSystem].cz*(REAL)NINT(dr.z*InverseBox[CurrentSystem].cz);
      break;
    case TRICLINIC:
      // convert from xyz to abc
      s.x=InverseBox[CurrentSystem].ax*dr.x+InverseBox[CurrentSystem].bx*dr.y+InverseBox[CurrentSystem].cx*dr.z;
      s.y=InverseBox[CurrentSystem].ay*dr.x+InverseBox[CurrentSystem].by*dr.y+InverseBox[CurrentSystem].cy*dr.z;
      s.z=InverseBox[CurrentSystem].az*dr.x+InverseBox[CurrentSystem].bz*dr.y+InverseBox[CurrentSystem].cz*dr.z;

      // apply boundary condition
      t.x=s.x-(REAL)NINT(s.x);
      t.y=s.y-(REAL)NINT(s.y);
      t.z=s.z-(REAL)NINT(s.z);

      // convert from abc to xyz
      dr.x=Box[CurrentSystem].ax*t.x+Box[CurrentSystem].bx*t.y+Box[CurrentSystem].cx*t.z;
      dr.y=Box[CurrentSystem].ay*t.x+Box[CurrentSystem].by*t.y+Box[CurrentSystem].cy*t.z;
      dr.z=Box[CurrentSystem].az*t.x+Box[CurrentSystem].bz*t.y+Box[CurrentSystem].cz*t.z;
      break;
    default:
      fprintf(stderr,"Error: Unkown boundary condition....\n");
      exit(0);
      break;
  }
  return dr;
}
*/

/*******************************************************************************************
 * Name       | ApplyReplicaBoundaryCondition                                              *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the Cartesian distance vector corrected for periodic boundaries    *
 * Parameters | VECTOR dr: the Cartesian distance vector                                   *
 * Note       | This function operates on the periodic replica cell                        *
 *******************************************************************************************/

VECTOR ApplyReplicaBoundaryCondition(VECTOR dr)
{
  VECTOR s,t;

  switch(BoundaryCondition[CurrentSystem])
  {
    case RECTANGULAR:
    case CUBIC:
      dr.x-=ReplicaBox[CurrentSystem].ax*(REAL)NINT(dr.x*InverseReplicaBox[CurrentSystem].ax);
      dr.y-=ReplicaBox[CurrentSystem].by*(REAL)NINT(dr.y*InverseReplicaBox[CurrentSystem].by);
      dr.z-=ReplicaBox[CurrentSystem].cz*(REAL)NINT(dr.z*InverseReplicaBox[CurrentSystem].cz);
      break;
    case FINITE:
      break;
    case TRICLINIC:
      // convert from xyz to abc
      s.x=InverseReplicaBox[CurrentSystem].ax*dr.x+InverseReplicaBox[CurrentSystem].bx*dr.y+InverseReplicaBox[CurrentSystem].cx*dr.z;
      s.y=InverseReplicaBox[CurrentSystem].ay*dr.x+InverseReplicaBox[CurrentSystem].by*dr.y+InverseReplicaBox[CurrentSystem].cy*dr.z;
      s.z=InverseReplicaBox[CurrentSystem].az*dr.x+InverseReplicaBox[CurrentSystem].bz*dr.y+InverseReplicaBox[CurrentSystem].cz*dr.z;

      // apply boundary condition
      t.x=s.x-(REAL)NINT(s.x);
      t.y=s.y-(REAL)NINT(s.y);
      t.z=s.z-(REAL)NINT(s.z);

      // convert from abc to xyz
      dr.x=ReplicaBox[CurrentSystem].ax*t.x+ReplicaBox[CurrentSystem].bx*t.y+ReplicaBox[CurrentSystem].cx*t.z;
      dr.y=ReplicaBox[CurrentSystem].ay*t.x+ReplicaBox[CurrentSystem].by*t.y+ReplicaBox[CurrentSystem].cy*t.z;
      dr.z=ReplicaBox[CurrentSystem].az*t.x+ReplicaBox[CurrentSystem].bz*t.y+ReplicaBox[CurrentSystem].cz*t.z;
      break;
    default:
      fprintf(stderr,"Error: Unkown boundary condition....\n");
      exit(0);
      break;
  }
  return dr;
}

/*******************************************************************************************
 * Name       | ApplyBoundaryConditionUnitCell                                             *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the Cartesian distance vector corrected for periodic boundaries    *
 * Parameters | VECTOR dr: the Cartesian distance vector                                   *
 * Note       | This function operates on the periodic unit cell                           *
 *******************************************************************************************/

VECTOR ApplyBoundaryConditionUnitCell(VECTOR dr)
{
  VECTOR s,t;

  switch(BoundaryCondition[CurrentSystem])
  {
    case FINITE:
      break;
    case CUBIC:
      dr.x-=UnitCellSize[CurrentSystem].x*NINT(dr.x/UnitCellSize[CurrentSystem].x);
      dr.y-=UnitCellSize[CurrentSystem].y*NINT(dr.y/UnitCellSize[CurrentSystem].y);
      dr.z-=UnitCellSize[CurrentSystem].z*NINT(dr.z/UnitCellSize[CurrentSystem].z);
      break;
    case RECTANGULAR:
      dr.x-=UnitCellSize[CurrentSystem].x*NINT(dr.x/UnitCellSize[CurrentSystem].x);
      dr.y-=UnitCellSize[CurrentSystem].y*NINT(dr.y/UnitCellSize[CurrentSystem].y);
      dr.z-=UnitCellSize[CurrentSystem].z*NINT(dr.z/UnitCellSize[CurrentSystem].z);
      break;
    case TRICLINIC:
      // convert from xyz to abc
      s.x=InverseUnitCellBox[CurrentSystem].ax*dr.x+InverseUnitCellBox[CurrentSystem].bx*dr.y+InverseUnitCellBox[CurrentSystem].cx*dr.z;
      s.y=InverseUnitCellBox[CurrentSystem].ay*dr.x+InverseUnitCellBox[CurrentSystem].by*dr.y+InverseUnitCellBox[CurrentSystem].cy*dr.z;
      s.z=InverseUnitCellBox[CurrentSystem].az*dr.x+InverseUnitCellBox[CurrentSystem].bz*dr.y+InverseUnitCellBox[CurrentSystem].cz*dr.z;

      // apply boundary condition
      t.x=s.x-NINT(s.x);
      t.y=s.y-NINT(s.y);
      t.z=s.z-NINT(s.z);

      // convert from abc to xyz
      dr.x=UnitCellBox[CurrentSystem].ax*t.x+UnitCellBox[CurrentSystem].bx*t.y+UnitCellBox[CurrentSystem].cx*t.z;
      dr.y=UnitCellBox[CurrentSystem].ay*t.x+UnitCellBox[CurrentSystem].by*t.y+UnitCellBox[CurrentSystem].cy*t.z;
      dr.z=UnitCellBox[CurrentSystem].az*t.x+UnitCellBox[CurrentSystem].bz*t.y+UnitCellBox[CurrentSystem].cz*t.z;
      break;
    default:
      fprintf(stderr,"Error: Unkown boundary condition....\n");
      exit(0);
      break;
  }
  return dr;
}

/*******************************************************************************************
 * Name       | ComputeDampingCoefficients                                                 *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the damping coefficients for the Pellenq-potential                 *
 * Parameters | REAL r: the Cartesian distance                                             *
 *            | REAL b: a parameter                                                        *
 * Returns    | REAL *f6:  the coefficient of the r^-6 term                                *
 *            | REAL *f8:  the coefficient of the r^-8 term                                *
 *            | REAL *f10: the coefficient of the r^-10 term                               *
 *******************************************************************************************/

void ComputeDampingCoefficients(REAL r, REAL b,REAL *f6,REAL *f8,REAL *f10)
{
  REAL sum,val;

  val=exp(-b*r);

  sum=val+
      b*r*val+
      pow(b*r,2)*val/2.0+
      pow(b*r,3)*val/6.0+
      pow(b*r,4)*val/24.0+
      pow(b*r,5)*val/120.0+
      pow(b*r,6)*val/720.0;
  *f6=1.0-sum;

  sum+=pow(b*r,7)*val/5040.0+pow(b*r,8)*val/40320.0;
  *f8=1.0-sum;
 
  sum+=pow(b*r,9)*val/362880.0+pow(b*r,10)*val/3628800.0;
  *f10=1.0-sum;
}

/*******************************************************************************************
 * Name       | ComputeDampingCoefficientsDerivatives                                      *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the damping coefficient first derivative for the Pellenq-potential *
 * Parameters | REAL r: the Cartesian distance                                             *
 *            | REAL b: a parameter                                                        *
 * Returns    | REAL *f6:  the damping first derivative of the r^-6 term                   *
 *            | REAL *f8:  the damping first derivative of the r^-8 term                   *
 *            | REAL *f10: the damping first derivative of the r^-10 term                  *
 *******************************************************************************************/

void ComputeDampingCoefficientsDerivatives(REAL r, REAL b,REAL *f6,REAL *f8,REAL *f10)
{
  REAL val;

  val=exp(-b*r);

  *f6=(1.0/720.0)*pow(b,7)*val*pow(r,6);
  *f8=(1.0/40320.0)*pow(b,9)*val*pow(r,8);
  *f10=(1.0/3628800.0)*pow(b,11)*val*pow(r,10);
}

/*******************************************************************************************
 * Name       | ComputeDampingCoefficientsSecondDerivatives                                *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the damping second derivatives for the Pellenq-potential           *
 * Parameters | REAL r: the Cartesian distance                                             *
 *            | REAL b: a parameter                                                        *
 * Returns    | REAL *f6:  the damping second derivative of the r^-6 term                  *
 *            | REAL *f8:  the damping second derivative of the r^-8 term                  *
 *            | REAL *f10: the damping second derivative of the r^-10 term                 *
 *******************************************************************************************/

void ComputeDampingCoefficientsSecondDerivatives(REAL r, REAL b,REAL *f6,REAL *f8,REAL *f10)
{
  REAL val;

  val=exp(-b*r);

  *f6=(1.0/120.0)*pow(b,7)*val*pow(r,5)-(1.0/720.0)*pow(b,8)*val*pow(r,6);
  *f8=(1.0/5040.0)*pow(b,9)*val*pow(r,7)-(1.0/40320.0)*pow(b,10)*val*pow(r,8);
  *f10=(1.0/362880.0)*pow(b,11)*val*pow(r,9)-(1.0/3628800.0)*pow(b,12)*val*pow(r,10);
}

/*******************************************************************************************
 * Name       | PotentialValue                                                             *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the energy of a potential U[r]                                     *
 * Parameters | int typeA: the atom-type of the first atom of the pair interaction         *
 *            | int typeB: the atom-type of the second atom of the pair interaction        *
 *            | REAL rr:   the squared Cartesian distance between the atoms                *
 * Note       | The square is passed to avoid a sqrt when possible (e.g. Lennard-Jones)    *
 *******************************************************************************************/

REAL PotentialValue(int typeA,int typeB,REAL rr)
{
  REAL fcVal,r,U,rri3,rri3_2;
  REAL arg1,arg2,arg3,arg4,arg5,arg6,arg7;
  REAL ri6,ri9;
  REAL exp1,exp2,exp_term,P;
  REAL f6,f8,f10;
  REAL rri2,rri4,rri6,rri8,rri10,rri12,rri14,rri16;
  REAL term1,term2;
  REAL energy,SwitchingValue;

  switch(PotentialType[typeA][typeB])
  {
    case UNDEFINED_POTENTIAL:
    case ZERO_POTENTIAL:
      U=0.0;
      fcVal=0.0;
      break;
    case LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_3/k_B [K]    (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2];
      rri3=CUBE(arg2/rr);
      return 4.0*arg1*(rri3*(rri3-1.0))-arg3;
    case LENNARD_JONES_SMOOTHED3:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*4.0*arg1*(rri3*(rri3-1.0));
      }
      return 4.0*arg1*(rri3*(rri3-1.0));
    case LENNARD_JONES_SMOOTHED5:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*4.0*arg1*(rri3*(rri3-1.0));
      }
      return 4.0*arg1*(rri3*(rri3-1.0));
    case WCA:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_3/k_B [K]    (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      if(rr>SQR(pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1])) return 0.0;
      arg3=PotentialParms[typeA][typeB][2];
      rri3=CUBE(arg2/rr);
      return 4.0*arg1*(rri3*(rri3-1.0))-arg3;
    case FEYNMAN_HIBBS_LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      // =============================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // p_3/k_B [K]    (non-zero for a shifted potential)
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(arg2/rr);
      return 4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr))-arg4;
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ====================================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // p_3/k_B [K]    (non-zero for a shifted potential)
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      rri3=CUBE(arg2/rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      }
      return 4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ====================================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // p_3/k_B [K]    (non-zero for a shifted potential)
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      rri3=CUBE(arg2/rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      }
      return 4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
    case FEYNMAN_HIBBS_LENNARD_JONES2:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      // ============================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // p_3/k_B [K]    (non-zero for a shifted potential)
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(arg2/rr);
      return 4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr))-arg4;
    case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ===================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      rri3=CUBE(arg2/rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      }
      return 4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
    case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ===================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      rri3=CUBE(arg2/rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      }
      return 4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
    case LENNARD_JONES_SHIFTED_FORCE:
      // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]-[(p_1/rc)^12-(p_1/rc)^6]}+[12*(p_1/rc)^12-6*(p_1/rc)^6]*(r-rc)/rc
      // ===============================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      r=sqrt(rr);
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      rri3_2=CUBE(arg2/CutOffVDWSquared);
      return 4.0*arg1*(rri3*(rri3-1.0)-rri3_2*(rri3_2-1.0)+(12*SQR(rri3_2)-6.0*rri3_2)*(r-CutOffVDW)/CutOffVDW);
    case LENNARD_JONES_SHIFTED_FORCE2:
      // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]+[6*(p_1/rc)^12-3*(p_1/rc)^6]}*r^2/rc^2+[7*(p_1/rc)^12+4*(p_1/rc)^6]
      // =================================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      rri3_2=CUBE(arg2/CutOffVDWSquared);
      return 4.0*arg1*((rri3*(rri3-1.0))+6.0*(rri3_2*(rri3_2-0.5))*rr/CutOffVDWSquared-7.0*SQR(rri3_2)+4.0*rri3_2);
    case POTENTIAL_12_6:
      // p_0/r^12-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      rri3=CUBE(1.0/rr);
      return arg1*SQR(rri3)-arg2*rri3-arg3;
    case POTENTIAL_12_6_SMOOTHED3:
      // {p_0/r^12-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri3=CUBE(1.0/rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(arg1*SQR(rri3)-arg2*rri3);
      }
      return arg1*SQR(rri3)-arg2*rri3;
    case POTENTIAL_12_6_SMOOTHED5:
      // {p_0/r^12-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri3=CUBE(1.0/rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(arg1*SQR(rri3)-arg2*rri3);
      }
      return arg1*SQR(rri3)-arg2*rri3;
    case POTENTIAL_12_6_2_0:
      // p_0/r^12+p_1/r^6+p_2/r^2+p_3
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K] 
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      rri3=CUBE(1.0/rr);
      return arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4-arg5;
    case POTENTIAL_12_6_2_0_SMOOTHED3:
      // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K] 
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(1.0/rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4);
      }
      return arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4;
    case POTENTIAL_12_6_2_0_SMOOTHED5:
      // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K] 
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(1.0/rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4);
      }
      return arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4;
    case MORSE:
      // p_0*{(1.0-exp[-p_1*(r-p_2)])^2-1.0}  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      return arg1*(SQR(1.0-exp_term)-1.0)-arg4;
    case MORSE_SMOOTHED3:
      // p_0*[(1.0-exp(-p_1*(r-p_2)))^2-1.0]  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*arg1*(SQR(1.0-exp_term)-1.0);
      }
      return arg1*(SQR(1.0-exp_term)-1.0);
    case MORSE_SMOOTHED5:
      // p_0*[(1.0-exp(-p_1*(r-p_2)))^2-1.0]  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*arg1*(SQR(1.0-exp_term)-1.0);
      }
      return arg1*(SQR(1.0-exp_term)-1.0);
    case MORSE2:
      // p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp_term=exp(arg2*(1-r/arg3))-2*exp(arg2/2*(1-r/arg3));
      return arg1*exp_term-arg4;
    case MORSE2_SMOOTHED3:
      // p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]*S(r)
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(arg2*(1-r/arg3))-2*exp(arg2/2*(1-r/arg3));
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*arg1*exp_term;
      }
      return arg1*exp_term;
    case MORSE2_SMOOTHED5:
      // p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]*S(r)
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(arg2*(1-r/arg3))-2*exp(arg2/2*(1-r/arg3));
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*arg1*exp_term;
      }
      return arg1*exp_term;
    case MORSE3:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      arg4=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      return arg1*(SQR(1.0-exp_term)-1.0)-arg4;
    case MORSE3_SMOOTHED3:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*arg1*(SQR(1.0-exp_term)-1.0);
      }
      return arg1*(SQR(1.0-exp_term)-1.0);
    case MORSE3_SMOOTHED5:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*arg1*(SQR(1.0-exp_term)-1.0);
      }
      return arg1*(SQR(1.0-exp_term)-1.0);
    case CFF_9_6:
      // p_0/r^9-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      ri6=CUBE(1.0/rr);
      ri9=ri6*sqrt(ri6);
      return arg1*ri9-arg2*ri6-arg3;
    case CFF_9_6_SMOOTHED3:
      // {p_0/r^9-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      ri6=CUBE(1.0/rr);
      ri9=ri6*sqrt(ri6);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(arg1*ri9-arg2*ri6);
      }
      return arg1*ri9-arg2*ri6;
    case CFF_9_6_SMOOTHED5:
      // {p_0/r^9-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      ri6=CUBE(1.0/rr);
      ri9=ri6*sqrt(ri6);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(arg1*ri9-arg2*ri6);
      }
      return arg1*ri9-arg2*ri6;
    case CFF_EPS_SIGMA:
      // p_0*[2*(p_1/r)^9-3(p_1/r)^6]
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2];
      ri6=CUBE(arg2/rr);
      ri9=ri6*sqrt(ri6);
      return arg1*(2.0*ri9-3.0*ri6)-arg3;
    case CFF_EPS_SIGMA_SMOOTHED3:
      // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      ri6=CUBE(arg2/rr);
      ri9=ri6*sqrt(ri6);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*arg1*(2.0*ri9-3.0*ri6);
      }
      return arg1*(2.0*ri9-3.0*ri6);
    case CFF_EPS_SIGMA_SMOOTHED5:
      // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      ri6=CUBE(arg2/rr);
      ri9=ri6*sqrt(ri6);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*arg1*(2.0*ri9-3.0*ri6);
      }
      return arg1*(2.0*ri9-3.0*ri6);
    case BUCKINGHAM:
      // p_0*exp(-p_1*r)-p_2/r^6
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      return exp_term-rri3-arg4;
    case BUCKINGHAM_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(exp_term-rri3);
      }
      return exp_term-rri3;
    case BUCKINGHAM_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(exp_term-rri3);
      }
      return exp_term-rri3;
    case BUCKINGHAM2:
      // if(r<p_3) 1e10 else p_0*exp(-p_1*r)-p_2/r^6
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3     [A]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      if(r<arg4) return 1e10;
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      return exp_term-rri3-arg5;
    case BUCKINGHAM2_SMOOTHED3:
      // if(r<p_3 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      if(r<arg4) return 1e10;
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(exp_term-rri3);
      }
      return exp_term-rri3;
    case BUCKINGHAM2_SMOOTHED5:
      // if(r<p_3) 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      if(r<arg4) return 1e10;
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(exp_term-rri3);
      }
      return exp_term-rri3;
    case MM3_VDW:
    case MM3_HYDROGEN_VDW:
      // sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]  if P>=3.02
      // sqrt(p_0^i*p_0^j)*192.27*P^2                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      // p_3     [kcal/mol]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      P=arg2/r;
      if(P>3.02)
        return arg1*192.270*SQR(P)-arg3;
      else
      {
        rri3=arg1*2.25*CUBE(SQR(P));
        exp_term=arg1*1.84e5*exp(-12.0/P);
        return exp_term-rri3-arg3;
      }
      break;
    case MM3_VDW_SMOOTHED3:
      // {sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
      // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      r=sqrt(rr);
      P=arg2/r;
      if(P>3.02)
        energy=arg1*192.270*SQR(P);
      else
      {
        rri3=arg1*2.25*CUBE(SQR(P));
        exp_term=arg1*1.84e5*exp(-12.0/P);
        energy=exp_term-rri3;
      }
      if(rr>CutOffVDWSwitchSquared)
      {
         SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
         energy*=SwitchingValue;
      }
      return energy;
    case MM3_VDW_SMOOTHED5:
      // {sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
      // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      r=sqrt(rr);
      P=arg2/r;
      if(P>3.02)
        energy=arg1*192.270*SQR(P);
      else
      {
        rri3=arg1*2.25*CUBE(SQR(P));
        exp_term=arg1*1.84e5*exp(-12.0/P);
        energy=exp_term-rri3;
      }
      if(rr>CutOffVDWSwitchSquared)
      {
         SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
         energy*=SwitchingValue;
      }
      return energy;
    case MATSUOKA_CLEMENTI_YOSHIMINE:
      // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      exp1=arg1*exp(-arg2*r);
      exp2=arg3*exp(-arg4*r);
      return exp1+exp2-arg5;
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3:
      // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp1=arg1*exp(-arg2*r);
      exp2=arg3*exp(-arg4*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(exp1+exp2);
      }
      return exp1+exp2;
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5:
      // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp1=arg1*exp(-arg2*r);
      exp2=arg3*exp(-arg4*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(exp1+exp2);
      }
      return exp1+exp2;
    case GENERIC:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      // p_6/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      arg7=PotentialParms[typeA][typeB][6];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      exp_term=arg1*exp(-arg2*r);
      return exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10-arg7;
    case GENERIC_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      exp_term=arg1*exp(-arg2*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10);
      }
      return exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10;
    case GENERIC_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      exp_term=arg1*exp(-arg2*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10);
      }
      return exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10;
    case PELLENQ_NICHOLSON:
      // p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      return arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10-arg6;
    case PELLENQ_NICHOLSON_SMOOTHED3:
      // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10);
      }
      return arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10;
    case PELLENQ_NICHOLSON_SMOOTHED5:
      // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10);
      }
      return arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10;
    case HYDRATED_ION_WATER:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      exp_term=arg1*exp(-arg2*r);
      return exp_term-arg3*rri4-arg4*rri6-arg5*rri12-arg6;
    case HYDRATED_ION_WATER_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      exp_term=arg1*exp(-arg2*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(exp_term-arg3*rri4-arg4*rri6-arg5*rri12);
      }
      return exp_term-arg3*rri4-arg4*rri6-arg5*rri12;
    case HYDRATED_ION_WATER_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      exp_term=arg1*exp(-arg2*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(exp_term-arg3*rri4-arg4*rri6-arg5*rri12);
      }
      return exp_term-arg3*rri4-arg4*rri6-arg5*rri12;
    case MIE:
      // p_0/r^p_1-p_2/r^p_3
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      return (arg1/pow(r,arg2)-arg3/pow(r,arg4))-arg5;
    case MIE_SMOOTHED3:
      // {p_0/r^p_1-p_2/r^p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      term1=arg1/pow(r,arg2);
      term2=arg3/pow(r,arg4);
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(term1-term2);
      }
      return (term1-term2);
    case MIE_SMOOTHED5:
      // {p_0/r^p_1-p_2/r^p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      term1=arg1/pow(r,arg2);
      term2=arg3/pow(r,arg4);
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(term1-term2);
      }
      return (term1-term2);
    case BORN_HUGGINS_MEYER:
      // p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      exp_term=arg1*exp(arg2*(arg3-r));
      return exp_term-arg4*rri6-arg5*rri8-arg6;
    case BORN_HUGGINS_MEYER_SMOOTHED3:
      // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      exp_term=arg1*exp(arg2*(arg3-r));
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(exp_term-arg4*rri6-arg5*rri8);
      }
      return exp_term-arg4*rri6-arg5*rri8;
    case BORN_HUGGINS_MEYER_SMOOTHED5:
      // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      exp_term=arg1*exp(arg2*(arg3-r));
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(exp_term-arg4*rri6-arg5*rri8);
      }
      return exp_term-arg4*rri6-arg5*rri8;
    case HYDROGEN:
      // p_0/r^12-p_1/r^10
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      return arg1*rri12-arg2*rri10-arg3;
    case HYDROGEN_SMOOTHED3:
      // {p_0/r^12-p_1/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        return SwitchingValue*(arg1*rri12-arg2*rri10);
      }
      return arg1*rri12-arg2*rri10;
    case HYDROGEN_SMOOTHED5:
      // {p_0/r^12-p_1/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                        SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        return SwitchingValue*(arg1*rri12-arg2*rri10);
      }
      return arg1*rri12-arg2*rri10;
    case LENNARD_JONES_CONTINUOUS_FRACTIONAL:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_3/k_B [K]    (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2];
      arg4=CurrentLambda;
      rri3=CUBE(arg2/rr);
      term1=0.5*SQR(1.0-arg4)+pow(rr/arg2,3);
      return 4.0*arg1*arg4*(1.0/SQR(term1)-1.0/term1);
    default:
      fprintf(stderr,"Undefined Potential in 'Potential Value'\n");
      exit(0);
      break;
  }
  return 0;
}

/*******************************************************************************************
 * Name       | PotentialGradient                                                          *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the energy and first derivative of a potential U[r]                *
 * Parameters | int typeA: the atom-type of the first atom of the pair interaction         *
 *            | int typeB: the atom-type of the second atom of the pair interaction        *
 *            | REAL rr:   the squared Cartesian distance between the atoms                *
 * Returns    | REAL *energy:  the contribution to the energy                              *
 *            | REAL *factor1: the first derivative     D[U[r],r]/r                        *
 * Note       | The square is passed to avoid a sqrt when possible (e.g. Lennard-Jones)    *
 *******************************************************************************************/

void PotentialGradient(int typeA,int typeB,REAL rr,REAL *energy,REAL *force_factor)
{
  REAL fcVal,r,U,rri3,rri3_2;
  REAL arg1,arg2,arg3,arg4,arg5,arg6,arg7;
  REAL ri6,ri9;
  REAL exp1,exp2,exp_term,P;
  REAL f6,f8,f10,f6d,f8d,f10d;
  REAL rri2,rri4,rri6,rri8,rri10,rri12,rri14,rri16;
  REAL term1,term2;
  REAL SwitchingValue,SwitchingValueDerivative;

  U=0.0;
  fcVal=0.0;
  switch(PotentialType[typeA][typeB])
  {
    case UNDEFINED_POTENTIAL:
    case ZERO_POTENTIAL:
      U=0.0;
      fcVal=0.0;
      break;
    case LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_3/k_B [K]    (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0))-arg3;
      fcVal=48.0*arg1*(rri3*(0.5-rri3))/rr;
      break;
    case LENNARD_JONES_SMOOTHED3:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0));
      fcVal=48.0*arg1*(rri3*(0.5-rri3))/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case LENNARD_JONES_SMOOTHED5:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0));
      fcVal=48.0*arg1*(rri3*(0.5-rri3))/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case WCA:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_3/k_B [K]    (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      if(rr<=SQR(pow(2.0,1.0/6.0))*arg2)
      {
        arg3=PotentialParms[typeA][typeB][2];
        rri3=CUBE(arg2/rr);
        U=4.0*arg1*(rri3*(rri3-1.0))-arg3;
        fcVal=48.0*arg1*(rri3*(0.5-rri3))/rr;
      }
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      // =============================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [u]    reduced mass in unified atomic mass units
      // p_3/k_B [K]    (non-zero for a shifted potential)
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr))-arg4;
      fcVal=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ====================================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [u]    reduced mass in unified atomic mass units
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      fcVal=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ====================================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [u]    reduced mass in unified atomic mass units
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      fcVal=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES2:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      // ============================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // p_3/k_B [K]    (non-zero for a shifted potential)
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr))-arg4;
      fcVal=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ===================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      fcVal=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ===================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      fcVal=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case LENNARD_JONES_SHIFTED_FORCE:
      // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]-[(p_1/rc)^12-(p_1/rc)^6]}+[12*(p_1/rc)^12-6*(p_1/rc)^6]*(r-rc)/rc
      // ===============================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      r=sqrt(rr);
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      rri3_2=CUBE(arg2/CutOffVDWSquared);
      U=4.0*arg1*(rri3*(rri3-1.0)-rri3_2*(rri3_2-1.0)+(12*SQR(rri3_2)-6.0*rri3_2)*(r-CutOffVDW)/CutOffVDW);
      fcVal=48.0*arg1*((rri3_2*(rri3_2-0.5)/(r*CutOffVDW))-(rri3*(rri3-0.5)/rr));
      break;
    case LENNARD_JONES_SHIFTED_FORCE2:
      // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]+[6*(p_1/rc)^12-3*(p_1/rc)^6]}*r^2/rc^2+[7*(p_1/rc)^12+4*(p_1/rc)^6]
      // =================================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      rri3_2=CUBE(arg2/CutOffVDWSquared);
      U=4.0*arg1*((rri3*(rri3-1.0))+6.0*(rri3_2*(rri3_2-0.5))*rr/CutOffVDWSquared-7.0*SQR(rri3_2)+4.0*rri3_2);
      fcVal=-48.0*arg1*(rri3*(rri3-0.5)/rr)+48.0*arg1*(rri3_2*(rri3_2-0.5)/CutOffVDWSquared);
      break;
    case POTENTIAL_12_6:
      // p_0/r^12-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)-arg2*rri3-arg3;
      fcVal=-(12.0/rr)*(arg1*SQR(rri3)-0.5*arg2*rri3);
      break;
    case POTENTIAL_12_6_SMOOTHED3:
      // {p_0/r^12-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)-arg2*rri3;
      fcVal=-(12.0/rr)*(arg1*SQR(rri3)-0.5*arg2*rri3);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case POTENTIAL_12_6_SMOOTHED5:
      // {p_0/r^12-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)-arg2*rri3;
      fcVal=-(12.0/rr)*(arg1*SQR(rri3)-0.5*arg2*rri3);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case POTENTIAL_12_6_2_0:
      // p_0/r^12+p_1/r^6+p_2/r^2+p_3
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K] 
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4-arg5;
      fcVal=-(12.0*arg1*SQR(rri3)+6.0*arg2*rri3+2.0*arg3/rr)/rr;
      break;
    case POTENTIAL_12_6_2_0_SMOOTHED3:
      // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K] 
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4;
      fcVal=-(12.0*arg1*SQR(rri3)+6.0*arg2*rri3+2.0*arg3/rr)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case POTENTIAL_12_6_2_0_SMOOTHED5:
      // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K] 
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4;
      fcVal=-(12.0*arg1*SQR(rri3)+6.0*arg2*rri3+2.0*arg3/rr)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE:
      // p_0*[(1.0-exp(-p_1*(r-p_2)))^2-1.0]  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0)-arg4;
      fcVal=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      break;
    case MORSE_SMOOTHED3:
      // p_0*[(1.0-exp(-p_1*(r-p_2)))^2-1.0]  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0)-arg4;
      fcVal=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE_SMOOTHED5:
      // p_0*[(1.0-exp(-p_1*(r-p_2)))^2-1.0]  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0)-arg4;
      fcVal=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE2:
      // p_0*[exp{p_1*(1-r/p_2)}-2*exp{0.5*p_1*(1-r/p_2)}]
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      U=arg1*(exp(arg2*(1.0-r/arg3))-2.0*exp(0.5*arg2*(1.0-r/arg3)))-arg4;
      fcVal=arg1*arg2*(exp(0.5*arg2*(1.0-r/arg3))-exp(arg2*(1.0-r/arg3)))/(arg3*r);
      break;
    case MORSE2_SMOOTHED3:
      // p_0*[exp{p_1*(1-r/p_2)}-2*exp{0.5*p_1*(1-r/p_2)}]
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      U=arg1*(exp(arg2*(1.0-r/arg3))-2.0*exp(0.5*arg2*(1.0-r/arg3)))-arg4;
      fcVal=arg1*arg2*(exp(0.5*arg2*(1.0-r/arg3))-exp(arg2*(1.0-r/arg3)))/(arg3*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE2_SMOOTHED5:
      // p_0*[exp{p_1*(1-r/p_2)}-2*exp{0.5*p_1*(1-r/p_2)}]
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      U=arg1*(exp(arg2*(1.0-r/arg3))-2.0*exp(0.5*arg2*(1.0-r/arg3)))-arg4;
      fcVal=arg1*arg2*(exp(0.5*arg2*(1.0-r/arg3))-exp(arg2*(1.0-r/arg3)))/(arg3*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE3:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      arg4=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0)-arg4;
      fcVal=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      break;
    case MORSE3_SMOOTHED3:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      arg4=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0)-arg4;
      fcVal=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE3_SMOOTHED5:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      arg4=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0)-arg4;
      fcVal=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case CFF_9_6:
      // p_0/r^9-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      ri6=CUBE(1.0/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*ri9-arg2*ri6-arg3;
      fcVal=-(9.0*arg1*ri9-6.0*arg2*ri6)/rr;
      break;
    case CFF_9_6_SMOOTHED3:
      // {p_0/r^9-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      ri6=CUBE(1.0/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*ri9-arg2*ri6;
      fcVal=-(9.0*arg1*ri9-6.0*arg2*ri6)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case CFF_9_6_SMOOTHED5:
      // {p_0/r^9-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      ri6=CUBE(1.0/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*ri9-arg2*ri6;
      fcVal=-(9.0*arg1*ri9-6.0*arg2*ri6)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case CFF_EPS_SIGMA:
      // p_0*[2*(p_1/r)^9-3(p_1/r)^6]
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2];
      ri6=CUBE(arg2/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*(2.0*ri9-3.0*ri6)-arg3;
      fcVal=-18.0*arg1*(ri9-ri6)/rr;
      break;
    case CFF_EPS_SIGMA_SMOOTHED3:
      // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      ri6=CUBE(arg2/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*(2.0*ri9-3.0*ri6);
      fcVal=-18.0*arg1*(ri9-ri6)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case CFF_EPS_SIGMA_SMOOTHED5:
      // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      ri6=CUBE(arg2/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*(2.0*ri9-3.0*ri6);
      fcVal=-18.0*arg1*(ri9-ri6)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BUCKINGHAM:
      // p_0*exp(-p_1*r)-p_2/r^6
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-rri3-arg4;
      fcVal=-(arg2*exp_term/r-(6.0/rr)*rri3);
      break;
    case BUCKINGHAM_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-rri3;
      fcVal=-(arg2*exp_term/r-(6.0/rr)*rri3);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BUCKINGHAM_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-rri3;
      fcVal=-(arg2*exp_term/r-(6.0/rr)*rri3);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BUCKINGHAM2:
      // if(r<p_3) 1e10 else p_0*exp(-p_1*r)-p_2/r^6
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3     [A]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      // Note: the case r<p_3 can not occur here as it is supposed to be avoided by the high energy in 'PotentialValue'
      //       in Molecular Dynamics the usual Buckingham should suffice
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-rri3-arg5;
      fcVal=-(arg2*exp_term/r-(6.0/rr)*rri3);
      break;
    case BUCKINGHAM2_SMOOTHED3:
      // if(r<p_3) 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3     [A]
      // Note: the case r<p_3 can not occur here as it is supposed to be avoided by the high energy in 'PotentialValue'
      //       in Molecular Dynamics the usual Buckingham should suffice
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-rri3;
      fcVal=-(arg2*exp_term/r-(6.0/rr)*rri3);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BUCKINGHAM2_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3     [A]
      // Note: the case r<p_3 can not occur here as it is supposed to be avoided by the high energy in 'PotentialValue'
      //       in Molecular Dynamics the usual Buckingham should suffice
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-rri3;
      fcVal=-(arg2*exp_term/r-(6.0/rr)*rri3);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MM3_VDW:
    case MM3_HYDROGEN_VDW:
      // sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]  if P>=3.02
      // sqrt(p_0^i*p_0^j)*192.27*P^2                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      // p_3     [kcal/mol]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      P=arg2/r;
      if(P>3.02) 
      { 
        U=arg1*192.270*SQR(P)-arg3;
        fcVal=-arg1*SQR(arg2)*384.54/(rr*rr);
      }
      else
      {
        rri3=arg1*2.25*CUBE(SQR(P));
        exp_term=arg1*1.84e5*exp(-12.0/P);
        U=exp_term-rri3-arg3;
        fcVal=-(12.0*exp_term/(arg2*r)-(6.0/rr)*rri3);
      }
      break;
    case MM3_VDW_SMOOTHED3:
      // {sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
      // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      P=arg2/r;
      if(P>3.02)
      {
        U=arg1*192.270*SQR(P);
        fcVal=-arg1*SQR(arg2)*384.54/(rr*rr);
      }
      else
      {
        rri3=arg1*2.25*CUBE(SQR(P));
        exp_term=arg1*1.84e5*exp(-12.0/P);
        U=exp_term-rri3;
        fcVal=-(12.0*exp_term/(arg2*r)-(6.0/rr)*rri3);
      }
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MM3_VDW_SMOOTHED5:
      // {sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
      // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      P=arg2/r;
      if(P>3.02)
      {
        U=arg1*192.270*SQR(P);
        fcVal=-arg1*SQR(arg2)*384.54/(rr*rr);
      }
      else
      {
        rri3=arg1*2.25*CUBE(SQR(P));
        exp_term=arg1*1.84e5*exp(-12.0/P);
        U=exp_term-rri3;
        fcVal=-(12.0*exp_term/(arg2*r)-(6.0/rr)*rri3);
      }
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MATSUOKA_CLEMENTI_YOSHIMINE:
      // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      exp1=arg1*exp(-arg2*r);
      exp2=arg3*exp(-arg4*r);
      U=exp1+exp2-arg5;
      fcVal=-(arg2*exp1+arg4*exp2)/r;
      break;
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3:
      // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp1=arg1*exp(-arg2*r);
      exp2=arg3*exp(-arg4*r);
      U=exp1+exp2;
      fcVal=-(arg2*exp1+arg4*exp2)/r;
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5:
      // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp1=arg1*exp(-arg2*r);
      exp2=arg3*exp(-arg4*r);
      U=exp1+exp2;
      fcVal=-(arg2*exp1+arg4*exp2)/r;
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case GENERIC:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      // p_6/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      arg7=PotentialParms[typeA][typeB][6];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10-arg7;
      fcVal=-(arg2*exp_term/r-4.0*arg3*rri6-6.0*arg4*rri8-8.0*arg5*rri10-10.0*arg6*rri12);
      break;
    case GENERIC_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10;
      fcVal=-(arg2*exp_term/r-4.0*arg3*rri6-6.0*arg4*rri8-8.0*arg5*rri10-10.0*arg6*rri12);
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case GENERIC_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10;
      fcVal=-(arg2*exp_term/r-4.0*arg3*rri6-6.0*arg4*rri8-8.0*arg5*rri10-10.0*arg6*rri12);
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case PELLENQ_NICHOLSON:
      // p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      ComputeDampingCoefficientsDerivatives(r,arg2,&f6d,&f8d,&f10d);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      U=arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10-arg6;
      fcVal=-((arg1*arg2*exp(-arg2*r)+f6d*arg3*rri6+f8d*arg4*rri8+f10d*arg5*rri10)/r
             -(f6*(6.0*arg3*rri8)+f8*(8.0*arg4*rri10)+f10*(10.0*arg5*rri12)));
      break;
    case PELLENQ_NICHOLSON_SMOOTHED3:
      // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      ComputeDampingCoefficientsDerivatives(r,arg2,&f6d,&f8d,&f10d);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      U=arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10;
      fcVal=-((arg1*arg2*exp(-arg2*r)+f6d*arg3*rri6+f8d*arg4*rri8+f10d*arg5*rri10)/r
             -(f6*(6.0*arg3*rri8)+f8*(8.0*arg4*rri10)+f10*(10.0*arg5*rri12)));
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case PELLENQ_NICHOLSON_SMOOTHED5:
      // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      ComputeDampingCoefficientsDerivatives(r,arg2,&f6d,&f8d,&f10d);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      U=arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10;
      fcVal=-((arg1*arg2*exp(-arg2*r)+f6d*arg3*rri6+f8d*arg4*rri8+f10d*arg5*rri10)/r
             -(f6*(6.0*arg3*rri8)+f8*(8.0*arg4*rri10)+f10*(10.0*arg5*rri12)));
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case HYDRATED_ION_WATER:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri12-arg6;
      fcVal=-(arg2*exp_term/r-4.0*arg3*rri6-6.0*arg4*rri8-12.0*arg5*rri14);
      break;
    case HYDRATED_ION_WATER_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri12;
      fcVal=-(arg2*exp_term/r-4.0*arg3*rri6-6.0*arg4*rri8-12.0*arg5*rri14);
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case HYDRATED_ION_WATER_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri12;
      fcVal=-(arg2*exp_term/r-4.0*arg3*rri6-6.0*arg4*rri8-12.0*arg5*rri14);
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MIE:
      // p_0/r^p_1-p_2/r^p_3
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      term1=arg1/pow(r,arg2);
      term2=arg3/pow(r,arg4);
      U=(term1-term2)-arg5;
      fcVal=(arg4*term2-arg2*term1)/rr;
      break;
    case MIE_SMOOTHED3:
      // {p_0/r^p_1-p_2/r^p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      term1=arg1/pow(r,arg2);
      term2=arg3/pow(r,arg4);
      U=(term1-term2);
      fcVal=(arg4*term2-arg2*term1)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MIE_SMOOTHED5:
      // {p_0/r^p_1-p_2/r^p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      term1=arg1/pow(r,arg2);
      term2=arg3/pow(r,arg4);
      U=(term1-term2);
      fcVal=(arg4*term2-arg2*term1)/rr;
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BORN_HUGGINS_MEYER:
      // p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      exp_term=arg1*exp(arg2*(arg3-r));
      U=exp_term-arg4*rri6-arg5*rri8-arg6;
      fcVal=-(arg2*exp_term/r-6.0*arg4*rri8-8.0*arg5*rri10);
      break;
    case BORN_HUGGINS_MEYER_SMOOTHED3:
      // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      exp_term=arg1*exp(arg2*(arg3-r));
      U=exp_term-arg4*rri6-arg5*rri8;
      fcVal=-(arg2*exp_term/r-6.0*arg4*rri8-8.0*arg5*rri10);
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BORN_HUGGINS_MEYER_SMOOTHED5:
      // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      exp_term=arg1*exp(arg2*(arg3-r));
      U=exp_term-arg4*rri6-arg5*rri8;
      fcVal=-(arg2*exp_term/r-6.0*arg4*rri8-8.0*arg5*rri10);
      if(rr>CutOffVDWSwitchSquared)
      {
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case HYDROGEN:
      // p_0/r^12-p_1/r^10
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      U=arg1*rri12-arg2*rri10-arg3;
      fcVal=-(12.0*arg1*rri14-10.0*arg2*rri12);
      break;
    case HYDROGEN_SMOOTHED3:
      // {p_0/r^12-p_1/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      U=arg1*rri12-arg2*rri10;
      fcVal=-(12.0*arg1*rri14-10.0*arg2*rri12);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=(SwitchingVDWFactors3[3]*(rr*r)+SwitchingVDWFactors3[2]*rr+
                        SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0]);
        SwitchingValueDerivative=(3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1]);
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case HYDROGEN_SMOOTHED5:
      // {p_0/r^12-p_1/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2; 
      rri16=rri14*rri2;
      U=arg1*rri12-arg2*rri10;
      fcVal=-(12.0*arg1*rri14-10.0*arg2*rri12);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        SwitchingValue=SwitchingVDWFactors5[5]*(rr*rr*r)+SwitchingVDWFactors5[4]*(rr*rr)+SwitchingVDWFactors5[3]*(rr*r)+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*rr*rr+4.0*SwitchingVDWFactors5[4]*rr*r+3.0*SwitchingVDWFactors5[3]*rr+
                                 2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        fcVal=U*SwitchingValueDerivative/r+fcVal*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    default:
      fprintf(stderr,"Undefined Potential in 'Potential Force'\n");
      exit(0);
      break;
  }
  *energy=U;
  *force_factor=fcVal;
}



/*******************************************************************************************
 * Name       | PotentialSecondDerivative                                                  *
 * --------------------------------------------------------------------------------------- *
 * Function   | Returns the energy, first and second derivative of a potential U[r]        *
 * Parameters | int typeA: the atom-type of the first atom of the pair interaction         *
 *            | int typeB: the atom-type of the second atom of the pair interaction        *
 *            | REAL rr:   the squared Cartesian distance between the atoms                *
 * Returns    | REAL *energy:  the contribution to the energy                              *
 *            | REAL *factor1: the first derivative     D[U[r],r]/r                        *
 *            | REAL *factor2: the second derivative    D[D[U[r],r]/r,r]/r                 *
 * Note       | The square is passed to avoid a sqrt when possible (e.g. Lennard-Jones)    *
 *******************************************************************************************/

void PotentialSecondDerivative(int typeA,int typeB,REAL rr,REAL *energy,REAL *factor1,REAL *factor2)
{
  REAL fcVal1,fcVal2,r,r3,r4,r5,U,rri3,ri6,ri9,rri3_2;
  REAL arg1,arg2,arg3,arg4,arg5,arg6,arg7;
  REAL exp_term,P;
  REAL exp1,exp2;
  REAL f6,f8,f10,f6d,f8d,f10d,f6d2,f8d2,f10d2;
  REAL rri2,rri4,rri6,rri8,rri10,rri12,rri14,rri16;
  REAL term1,term2;
  REAL SwitchingValue,SwitchingValueDerivative,SwitchingValueDerivative2;

  U=0.0;
  fcVal1=0.0;
  fcVal2=0.0;
  switch(PotentialType[typeA][typeB])
  {
    case UNDEFINED_POTENTIAL:
    case ZERO_POTENTIAL:
      U=0.0;
      fcVal1=0.0;
      fcVal2=0.0;
      break;
    case LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_3/k_B [K]    (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0))-arg3;
      fcVal1=24.0*arg1*(rri3*(1.0-2.0*rri3))/rr;
      fcVal2=96.0*arg1*(rri3*(7.0*rri3-2.0))/SQR(rr);
      break;
    case LENNARD_JONES_SMOOTHED3:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0));
      fcVal1=24.0*arg1*(rri3*(1.0-2.0*rri3))/rr;
      fcVal2=96.0*arg1*(rri3*(7.0*rri3-2.0))/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case LENNARD_JONES_SMOOTHED5:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0));
      fcVal1=24.0*arg1*(rri3*(1.0-2.0*rri3))/rr;
      fcVal2=96.0*arg1*(rri3*(7.0*rri3-2.0))/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case WCA:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_3/k_B [K]    (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      if(rr<=SQR(pow(2.0,1.0/6.0))*arg2)
      {
        arg3=PotentialParms[typeA][typeB][2];
        rri3=CUBE(arg2/rr);
        U=4.0*arg1*(rri3*(rri3-1.0))-arg3;
        fcVal1=24.0*arg1*(rri3*(1.0-2.0*rri3))/rr;
        fcVal2=96.0*arg1*(rri3*(7.0*rri3-2.0))/SQR(rr);
      }
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      // =============================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [u]    reduced mass in unified atomic mass units
      // p_3/k_B [K]    (non-zero for a shifted potential)
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr))-arg4;
      fcVal1=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      fcVal2=96.0*arg1*((rri3*(7.0*rri3-2.0))+arg3*4.0*rri3*(308.0*rri3-25.0)/rr)/SQR(rr);
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ====================================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [u]    reduced mass in unified atomic mass units
      // p_3/k_B [K]    (non-zero for a shifted potential)
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      fcVal1=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      fcVal2=96.0*arg1*((rri3*(7.0*rri3-2.0))+arg3*4.0*rri3*(308.0*rri3-25.0)/rr)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ====================================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [u]    reduced mass in unified atomic mass units
      // p_3/k_B [K]    (non-zero for a shifted potential)
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      fcVal1=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      fcVal2=96.0*arg1*((rri3*(7.0*rri3-2.0))+arg3*4.0*rri3*(308.0*rri3-25.0)/rr)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES2:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      // ============================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // p_3/k_B [K]    (non-zero for a shifted potential)
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr))-arg4;
      fcVal1=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      fcVal2=96.0*arg1*((rri3*(7.0*rri3-2.0))+arg3*4.0*rri3*(308.0*rri3-25.0)/rr)/SQR(rr);
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ===================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      fcVal1=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      fcVal2=96.0*arg1*((rri3*(7.0*rri3-2.0))+arg3*4.0*rri3*(308.0*rri3-25.0)/rr)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5:
      // {4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2}*S(r)
      // ===================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0)+arg3*(rri3*(132.0*rri3-30.0)/rr));
      fcVal1=48.0*arg1*(rri3*(0.5-rri3)+arg3*rri3*(20.0-154.0*rri3)/rr)/rr;
      fcVal2=96.0*arg1*((rri3*(7.0*rri3-2.0))+arg3*4.0*rri3*(308.0*rri3-25.0)/rr)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case LENNARD_JONES_SHIFTED_FORCE:
      // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]-[(p_1/rc)^12-(p_1/rc)^6]}+[12*(p_1/rc)^12-6*(p_1/rc)^6]*(r-rc)/rc
      // ===============================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      r=sqrt(rr);
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      rri3_2=CUBE(arg2/CutOffVDWSquared);
      U=4.0*arg1*(rri3*(rri3-1.0)-rri3_2*(rri3_2-1.0)+(12*SQR(rri3_2)-6.0*rri3_2)*(r-CutOffVDW)/CutOffVDW);
      fcVal1=48.0*arg1*((rri3_2*(rri3_2-0.5)/(r*CutOffVDW))-(rri3*(rri3-0.5)/rr));
      fcVal2=4.0*arg1*(156*rri3*rri3-42.0*rri3)/SQR(rr)-fcVal1/rr;
      break;
    case LENNARD_JONES_SHIFTED_FORCE2:
      // 4*p_0*{[(p_1/r)^12-(p_1/r)^6]+[6*(p_1/rc)^12-3*(p_1/rc)^6]}*r^2/rc^2+[7*(p_1/rc)^12+4*(p_1/rc)^6]
      // =================================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      rri3=CUBE(arg2/rr);
      rri3_2=CUBE(arg2/CutOffVDWSquared);
      U=4.0*arg1*((rri3*(rri3-1.0))+6.0*(rri3_2*(rri3_2-0.5))*rr/CutOffVDWSquared-7.0*SQR(rri3_2)+4.0*rri3_2);
      fcVal1=-48.0*arg1*(rri3*(rri3-0.5)/rr)+48.0*arg1*(rri3_2*(rri3_2-0.5)/CutOffVDWSquared);
      fcVal2=96.0*arg1*(rri3*(7.0*rri3-2.0))/SQR(rr);
      break;
    case POTENTIAL_12_6:
      // p_0/r^12-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)-arg2*rri3-arg3;
      fcVal1=-(12.0/rr)*(arg1*SQR(rri3)-0.5*arg2*rri3);
      fcVal2=(24.0/SQR(rr))*(7.0*arg1*SQR(rri3)-2.0*arg2*rri3);
      break;
    case POTENTIAL_12_6_SMOOTHED3:
      // {p_0/r^12-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)-arg2*rri3;
      fcVal1=-(12.0/rr)*(arg1*SQR(rri3)-0.5*arg2*rri3);
      fcVal2=(24.0/SQR(rr))*(7.0*arg1*SQR(rri3)-2.0*arg2*rri3);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case POTENTIAL_12_6_SMOOTHED5:
      // {p_0/r^12-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)-arg2*rri3;
      fcVal1=-(12.0/rr)*(arg1*SQR(rri3)-0.5*arg2*rri3);
      fcVal2=(24.0/SQR(rr))*(7.0*arg1*SQR(rri3)-2.0*arg2*rri3);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case POTENTIAL_12_6_2_0:
      // p_0/r^12+p_1/r^6+p_2/r^2+p_3
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K] 
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4-arg5;
      fcVal1=-(12.0*arg1*SQR(rri3)+6.0*arg2*rri3+2.0*arg3/rr)/rr;
      fcVal2=8.0*(21.0*arg1*SQR(rri3)+6.0*arg2*rri3+arg3/rr)/SQR(rr);
      break;
    case POTENTIAL_12_6_2_0_SMOOTHED3:
      // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K] 
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4;
      fcVal1=-(12.0*arg1*SQR(rri3)+6.0*arg2*rri3+2.0*arg3/rr)/rr;
      fcVal2=8.0*(21.0*arg1*SQR(rri3)+6.0*arg2*rri3+arg3/rr)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case POTENTIAL_12_6_2_0_SMOOTHED5:
      // {p_0/r^12+p_1/r^6+p_2/r^2+p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      // p_2/k_B [K A^2]
      // p_3/k_B [K] 
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      rri3=CUBE(1.0/rr);
      U=arg1*SQR(rri3)+arg2*rri3+arg3/rr+arg4;
      fcVal1=-(12.0*arg1*SQR(rri3)+6.0*arg2*rri3+2.0*arg3/rr)/rr;
      fcVal2=8.0*(21.0*arg1*SQR(rri3)+6.0*arg2*rri3+arg3/rr)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE:
      // p_0*[(1.0-exp(-p_1*(r-p_2)))^2-1.0]  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0)-arg4;
      fcVal1=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      fcVal2=2.0*arg1*arg2*exp_term*(exp_term*(2.0*arg2*r+1.0)-arg2*r-1.0)/(rr*r);
      break;
    case MORSE_SMOOTHED3:
      // p_0*[(1.0-exp(-p_1*(r-p_2)))^2-1.0]*S(r)
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0);
      fcVal1=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      fcVal2=2.0*arg1*arg2*exp_term*(exp_term*(2.0*arg2*r+1.0)-arg2*r-1.0)/(rr*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE_SMOOTHED5:
      // p_0*[(1.0-exp(-p_1*(r-p_2)))^2-1.0]*S(r)
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0);
      fcVal1=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      fcVal2=2.0*arg1*arg2*exp_term*(exp_term*(2.0*arg2*r+1.0)-arg2*r-1.0)/(rr*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE2:
      // p_0*[exp{p_1*(1-r/p_2)}-2*exp{0.5*p_1*(1-r/p_2)}]
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      U=arg1*(exp(arg2*(1.0-r/arg3))-2.0*exp(0.5*arg2*(1.0-r/arg3)))-arg4;
      fcVal1=arg1*arg2*(exp(0.5*arg2*(1.0-r/arg3))-exp(arg2*(1.0-r/arg3)))/(arg3*r);
      fcVal2=arg1*arg2*exp(-arg2*r/arg3)*(exp(arg2*(0.5+0.5*r/arg3))*(-arg3-0.5*arg2*r)+exp(arg2)*(arg3+arg2*r))/(rr*r*SQR(arg3));
      break;
    case MORSE2_SMOOTHED3:
      // p_0*[exp{p_1*(1-r/p_2)}-2*exp{0.5*p_1*(1-r/p_2)}]
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      U=arg1*(exp(arg2*(1.0-r/arg3))-2.0*exp(0.5*arg2*(1.0-r/arg3)))-arg4;
      fcVal1=arg1*arg2*(exp(0.5*arg2*(1.0-r/arg3))-exp(arg2*(1.0-r/arg3)))/(arg3*r);
      fcVal2=arg1*arg2*exp(-arg2*r/arg3)*(exp(arg2*(0.5+0.5*r/arg3))*(-arg3-0.5*arg2*r)+exp(arg2)*(arg3+arg2*r))/(rr*r*SQR(arg3));
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE2_SMOOTHED5:
      // p_0*[exp{p_1*(1-r/p_2)}-2*exp{0.5*p_1*(1-r/p_2)}]
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      U=arg1*(exp(arg2*(1.0-r/arg3))-2.0*exp(0.5*arg2*(1.0-r/arg3)))-arg4;
      fcVal1=arg1*arg2*(exp(0.5*arg2*(1.0-r/arg3))-exp(arg2*(1.0-r/arg3)))/(arg3*r);
      fcVal2=arg1*arg2*exp(-arg2*r/arg3)*(exp(arg2*(0.5+0.5*r/arg3))*(-arg3-0.5*arg2*r)+exp(arg2)*(arg3+arg2*r))/(rr*r*SQR(arg3));
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE3:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      arg4=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0)-arg4;
      fcVal1=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      fcVal2=2.0*arg1*arg2*exp_term*(exp_term*(2.0*arg2*r+1.0)-arg2*r-1.0)/(rr*r);
      break;
    case MORSE3_SMOOTHED3:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}*S(r)
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      arg4=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0)-arg4;
      fcVal1=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      fcVal2=2.0*arg1*arg2*exp_term*(exp_term*(2.0*arg2*r+1.0)-arg2*r-1.0)/(rr*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MORSE3_SMOOTHED5:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}*S(r)
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      arg4=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      exp_term=exp(-arg2*(r-arg3));
      U=arg1*(SQR(1.0-exp_term)-1.0)-arg4;
      fcVal1=2.0*arg1*arg2*(1.0-exp_term)*exp_term/r;
      fcVal2=2.0*arg1*arg2*exp_term*(exp_term*(2.0*arg2*r+1.0)-arg2*r-1.0)/(rr*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case CFF_9_6:
      // p_0/r^9-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      ri6=CUBE(1.0/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*ri9-arg2*ri6-arg3;
      fcVal1=-(9.0*arg1*ri9-6.0*arg2*ri6)/rr;
      fcVal2=(99.0*arg1*ri9-48.0*arg2*ri6)/SQR(rr);
      break;
    case CFF_9_6_SMOOTHED3:
      // {p_0/r^9-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      ri6=CUBE(1.0/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*ri9-arg2*ri6;
      fcVal1=-(9.0*arg1*ri9-6.0*arg2*ri6)/rr;
      fcVal2=(99.0*arg1*ri9-48.0*arg2*ri6)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case CFF_9_6_SMOOTHED5:
      // {p_0/r^9-p_1/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      ri6=CUBE(1.0/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*ri9-arg2*ri6;
      fcVal1=-(9.0*arg1*ri9-6.0*arg2*ri6)/rr;
      fcVal2=(99.0*arg1*ri9-48.0*arg2*ri6)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case CFF_EPS_SIGMA:
      // p_0*[2*(p_1/r)^9-3(p_1/r)^6]
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2];
      ri6=CUBE(arg2/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*(2.0*ri9-3.0*ri6)-arg3;
      fcVal1=-18.0*arg1*(ri9-ri6)/rr;
      fcVal2=18.0*arg1*(11.0*ri9-8.0*ri6)/SQR(rr);
      break;
    case CFF_EPS_SIGMA_SMOOTHED3:
      // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      ri6=CUBE(arg2/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*(2.0*ri9-3.0*ri6);
      fcVal1=-18.0*arg1*(ri9-ri6)/rr;
      fcVal2=18.0*arg1*(11.0*ri9-8.0*ri6)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case CFF_EPS_SIGMA_SMOOTHED5:
      // {p_0*[2*(p_1/r)^9-3(p_1/r)^6]}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      ri6=CUBE(arg2/rr);
      ri9=ri6*sqrt(ri6);
      U=arg1*(2.0*ri9-3.0*ri6);
      fcVal1=-18.0*arg1*(ri9-ri6)/rr;
      fcVal2=18.0*arg1*(11.0*ri9-8.0*ri6)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BUCKINGHAM:
      // p_0*exp(-p_1*r)-p_2/r^6
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=-rri3+exp_term-arg4;
      fcVal1=-arg2*exp_term/r+(6.0/rr)*rri3;
      fcVal2=-48.0*rri3/(rr*rr)+arg2*exp_term*(1.0+arg2*r)/(rr*r);
      break;
    case BUCKINGHAM_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=-rri3+exp_term;
      fcVal1=-arg2*exp_term/r+(6.0/rr)*rri3;
      fcVal2=-48.0*rri3/(rr*rr)+arg2*exp_term*(1.0+arg2*r)/(rr*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BUCKINGHAM_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=-rri3+exp_term;
      fcVal1=-arg2*exp_term/r+(6.0/rr)*rri3;
      fcVal2=-48.0*rri3/(rr*rr)+arg2*exp_term*(1.0+arg2*r)/(rr*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BUCKINGHAM2:
      // if(r<p_3) 1e10 else p_0*exp(-p_1*r)-p_2/r^6
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3     [A]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      // Note: the case r<p_3 can not occur here as it is supposed to be avoided by the high energy in 'PotentialValue'
      //       in Molecular Dynamics the usual Buckingham should suffice
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=-rri3+exp_term-arg5;
      fcVal1=-arg2*exp_term/r+(6.0/rr)*rri3;
      fcVal2=-48.0*rri3/(rr*rr)+arg2*exp_term*(1.0+arg2*r)/(rr*r);
      break;
    case BUCKINGHAM2_SMOOTHED3:
      // if(r<p_3) 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3     [A]
      // Note: the case r<p_3 can not occur here as it is supposed to be avoided by the high energy in 'PotentialValue'
      //       in Molecular Dynamics the usual Buckingham should suffice
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=-rri3+exp_term;
      fcVal1=-arg2*exp_term/r+(6.0/rr)*rri3;
      fcVal2=-48.0*rri3/(rr*rr)+arg2*exp_term*(1.0+arg2*r)/(rr*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BUCKINGHAM2_SMOOTHED5:
      // if(r<p_3) 1e10 else {p_0*exp(-p_1*r)-p_2/r^6}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3     [A]
      // Note: the case r<p_3 can not occur here as it is supposed to be avoided by the high energy in 'PotentialValue'
      //       in Molecular Dynamics the usual Buckingham should suffice
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      rri3=arg3*CUBE(1.0/rr);
      exp_term=arg1*exp(-arg2*r);
      U=-rri3+exp_term;
      fcVal1=-arg2*exp_term/r+(6.0/rr)*rri3;
      fcVal2=-48.0*rri3/(rr*rr)+arg2*exp_term*(1.0+arg2*r)/(rr*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MM3_VDW:
    case MM3_HYDROGEN_VDW:
      // sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]  if P>=3.02
      // sqrt(p_0^i*p_0^j)*192.27*P^2                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      // p_3     [kcal/mol]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      P=arg2/r;
      if(P>3.02)
      {
        U=arg1*192.270*SQR(P)-arg3;
        fcVal1=-arg1*SQR(arg2)*384.54/SQR(rr);
        fcVal2=arg1*SQR(arg2)*1538.16/(SQR(rr)*rr);
      }
      else
      {
        rri3=arg1*2.25*CUBE(SQR(P));
        exp_term=arg1*1.84e5*exp(-12.0/P);
        U=exp_term-rri3-arg3;
        fcVal1=-(12.0*exp_term/(arg2*r)-(6.0/rr)*rri3);
        fcVal2=12.0*(12.0*r+arg2)*exp_term/(r*r*r*arg2*arg2)-48.0*rri3/SQR(rr);
      }
      break;
    case MM3_VDW_SMOOTHED3:
      // {sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
      // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      P=arg2/r;
      if(P>3.02)
      {
        U=arg1*192.270*SQR(P)-arg3;
        fcVal1=-arg1*SQR(arg2)*384.54/SQR(rr);
        fcVal2=arg1*SQR(arg2)*1538.16/(SQR(rr)*rr);
      }
      else
      {
        rri3=arg1*2.25*CUBE(SQR(P));
        exp_term=arg1*1.84e5*exp(-12.0/P);
        U=exp_term-rri3-arg3;
        fcVal1=-(12.0*exp_term/(arg2*r)-(6.0/rr)*rri3);
        fcVal2=12.0*(12.0*r+arg2)*exp_term/(r*r*r*arg2*arg2)-48.0*rri3/SQR(rr);
      }
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MM3_VDW_SMOOTHED5:
      // {sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]}*S(r)  if P>=3.02
      // {sqrt(p_0^i*p_0^j)*192.27*P^2}*S(r)                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      r=sqrt(rr);
      P=arg2/r;
      if(P>3.02)
      {
        U=arg1*192.270*SQR(P)-arg3;
        fcVal1=-arg1*SQR(arg2)*384.54/SQR(rr);
        fcVal2=arg1*SQR(arg2)*1538.16/(SQR(rr)*rr);
      }
      else
      {
        rri3=arg1*2.25*CUBE(SQR(P));
        exp_term=arg1*1.84e5*exp(-12.0/P);
        U=exp_term-rri3-arg3;
        fcVal1=-(12.0*exp_term/(arg2*r)-(6.0/rr)*rri3);
        fcVal2=12.0*(12.0*r+arg2)*exp_term/(r*r*r*arg2*arg2)-48.0*rri3/SQR(rr);
      }
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MATSUOKA_CLEMENTI_YOSHIMINE:
      // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      exp1=arg1*exp(-arg2*r);
      exp2=arg3*exp(-arg4*r);
      U=exp1+exp2-arg5;
      fcVal1=-(arg2*exp1+arg4*exp2)/r;
      fcVal2=arg2*exp1*(1.0+arg2*r)/(rr*r)+arg4*exp2*(1.0+arg4*r)/(rr*r);
      break;
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3:
      // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp1=arg1*exp(-arg2*r);
      exp2=arg3*exp(-arg4*r);
      U=exp1+exp2;
      fcVal1=-(arg2*exp1+arg4*exp2)/r;
      fcVal2=arg2*exp1*(1.0+arg2*r)/(rr*r)+arg4*exp2*(1.0+arg4*r)/(rr*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5:
      // {p_0*exp(-p_1*r)+p_2*exp(-p_3*r)}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      exp1=arg1*exp(-arg2*r);
      exp2=arg3*exp(-arg4*r);
      U=exp1+exp2;
      fcVal1=-(arg2*exp1+arg4*exp2)/r;
      fcVal2=arg2*exp1*(1.0+arg2*r)/(rr*r)+arg4*exp2*(1.0+arg4*r)/(rr*r);
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case GENERIC:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      // p_6/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      arg7=PotentialParms[typeA][typeB][6];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10-arg7;
      fcVal1=-arg2*exp_term/r+4.0*arg3*rri6+6.0*arg4*rri8+8.0*arg5*rri10+10.0*arg6*rri12;
      fcVal2=arg2*exp_term*(1.0+arg2*r)/(rr*r)-24.0*arg3*rri8-48.0*arg4*rri10-80.0*arg5*rri12-120*arg6*rri14;
      break;
    case GENERIC_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10;
      fcVal1=-arg2*exp_term/r+4.0*arg3*rri6+6.0*arg4*rri8+8.0*arg5*rri10+10.0*arg6*rri12;
      fcVal2=arg2*exp_term*(1.0+arg2*r)/(rr*r)-24.0*arg3*rri8-48.0*arg4*rri10-80.0*arg5*rri12-120*arg6*rri14;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case GENERIC_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri8-arg6*rri10;
      fcVal1=-arg2*exp_term/r+4.0*arg3*rri6+6.0*arg4*rri8+8.0*arg5*rri10+10.0*arg6*rri12;
      fcVal2=arg2*exp_term*(1.0+arg2*r)/(rr*r)-24.0*arg3*rri8-48.0*arg4*rri10-80.0*arg5*rri12-120*arg6*rri14;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case PELLENQ_NICHOLSON:
      // p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      ComputeDampingCoefficientsDerivatives(r,arg2,&f6d,&f8d,&f10d);
      ComputeDampingCoefficientsSecondDerivatives(r,arg2,&f6d2,&f8d2,&f10d2);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      U=arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10-arg6;
      fcVal1=f6*(6.0*arg3*rri8)+f8*(8.0*arg4*rri10)+f10*(10.0*arg5*rri12)-arg1*arg2*exp(-arg2*r)/r
             -(f6d*arg3*rri6+f8d*arg4*rri8+f10d*arg5*rri10)/r;
      fcVal2=-fcVal1/rr
             +arg1*SQR(arg2)*exp(-arg2*r)/rr-42.0*arg3*f6*rri10-72.0*arg4*f8*rri12-110.0*arg5*f10*rri14
             -arg3*f6d2*rri8-arg4*f8d2*rri10-arg5*f10d2*rri12
             +12.0*arg3*f6d*rri8/r+16.0*arg4*f8d*rri10/r+20.0*arg5*f10d*rri12/r;
      break;
    case PELLENQ_NICHOLSON_SMOOTHED3:
      // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      ComputeDampingCoefficientsDerivatives(r,arg2,&f6d,&f8d,&f10d);
      ComputeDampingCoefficientsSecondDerivatives(r,arg2,&f6d2,&f8d2,&f10d2);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      U=arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10;
      fcVal1=f6*(6.0*arg3*rri8)+f8*(8.0*arg4*rri10)+f10*(10.0*arg5*rri12)-arg1*arg2*exp(-arg2*r)/r
             -(f6d*arg3*rri6+f8d*arg4*rri8+f10d*arg5*rri10)/r;
      fcVal2=-fcVal1/rr
             +arg1*SQR(arg2)*exp(-arg2*r)/rr-42.0*arg3*f6*rri10-72.0*arg4*f8*rri12-110.0*arg5*f10*rri14
             -arg3*f6d2*rri8-arg4*f8d2*rri10-arg5*f10d2*rri12
             +12.0*arg3*f6d*rri8/r+16.0*arg4*f8d*rri10/r+20.0*arg5*f10d*rri12/r;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case PELLENQ_NICHOLSON_SMOOTHED5:
      // {p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      ComputeDampingCoefficients(r,arg2,&f6,&f8,&f10);
      ComputeDampingCoefficientsDerivatives(r,arg2,&f6d,&f8d,&f10d);
      ComputeDampingCoefficientsSecondDerivatives(r,arg2,&f6d2,&f8d2,&f10d2);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      U=arg1*exp(-arg2*r)-f6*arg3*rri6-f8*arg4*rri8-f10*arg5*rri10;
      fcVal1=f6*(6.0*arg3*rri8)+f8*(8.0*arg4*rri10)+f10*(10.0*arg5*rri12)-arg1*arg2*exp(-arg2*r)/r
             -(f6d*arg3*rri6+f8d*arg4*rri8+f10d*arg5*rri10)/r;
      fcVal2=-fcVal1/rr
             +arg1*SQR(arg2)*exp(-arg2*r)/rr-42.0*arg3*f6*rri10-72.0*arg4*f8*rri12-110.0*arg5*f10*rri14
             -arg3*f6d2*rri8-arg4*f8d2*rri10-arg5*f10d2*rri12
             +12.0*arg3*f6d*rri8/r+16.0*arg4*f8d*rri10/r+20.0*arg5*f10d*rri12/r;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case HYDRATED_ION_WATER:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri12-arg6;
      fcVal1=-arg2*exp_term/r+4.0*arg3*rri6+6.0*arg4*rri8+12.0*arg5*rri14;
      fcVal2=arg2*exp_term*(1.0+arg2*r)/(rr*r)-24.0*arg3*rri8-48.0*arg4*rri10-168.0*arg5*rri16;
      break;
    case HYDRATED_ION_WATER_SMOOTHED3:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri12;
      fcVal1=-arg2*exp_term/r+4.0*arg3*rri6+6.0*arg4*rri8+12.0*arg5*rri14;
      fcVal2=arg2*exp_term*(1.0+arg2*r)/(rr*r)-24.0*arg3*rri8-48.0*arg4*rri10-168.0*arg5*rri16;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case HYDRATED_ION_WATER_SMOOTHED5:
      // {p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=SQR(rri2);
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      exp_term=arg1*exp(-arg2*r);
      U=exp_term-arg3*rri4-arg4*rri6-arg5*rri12;
      fcVal1=-arg2*exp_term/r+4.0*arg3*rri6+6.0*arg4*rri8+12.0*arg5*rri14;
      fcVal2=arg2*exp_term*(1.0+arg2*r)/(rr*r)-24.0*arg3*rri8-48.0*arg4*rri10-168.0*arg5*rri16;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MIE:
      // p_0/r^p_1-p_2/r^p_3
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      // p_4/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      term1=arg1/pow(r,arg2);
      term2=arg3/pow(r,arg4);
      U=(term1-term2)-arg5;
      fcVal1=(arg4*term2-arg2*term1)/rr;
      fcVal2=(arg2*(2+arg2)*term1-arg4*(2+arg4)*term2)/SQR(rr);
      break;
    case MIE_SMOOTHED3:
      // {p_0/r^p_1-p_2/r^p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      term1=arg1/pow(r,arg2);
      term2=arg3/pow(r,arg4);
      U=(term1-term2);
      fcVal1=(arg4*term2-arg2*term1)/rr;
      fcVal2=(arg2*(2+arg2)*term1-arg4*(2+arg4)*term2)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case MIE_SMOOTHED5:
      // {p_0/r^p_1-p_2/r^p_3}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      r=sqrt(rr);
      term1=arg1/pow(r,arg2);
      term2=arg3/pow(r,arg4);
      U=(term1-term2);
      fcVal1=(arg4*term2-arg2*term1)/rr;
      fcVal2=(arg2*(2+arg2)*term1-arg4*(2+arg4)*term2)/SQR(rr);
      if(rr>CutOffVDWSwitchSquared)
      {
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BORN_HUGGINS_MEYER:
      // p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      exp_term=arg1*exp(arg2*(arg3-r));
      U=exp_term-arg4*rri6-arg5*rri8-arg6;
      fcVal1=6.0*arg4*rri8+8.0*arg5*rri10-arg2*exp_term/r;
      fcVal2=arg2*(1.0+arg2*r)*exp_term/(rr*r)-48.0*arg4*rri10-80*arg5*rri12;
      break;
    case BORN_HUGGINS_MEYER_SMOOTHED3:
      // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      exp_term=arg1*exp(arg2*(arg3-r));
      U=exp_term-arg4*rri6-arg5*rri8;
      fcVal1=6.0*arg4*rri8+8.0*arg5*rri10-arg2*exp_term/r;
      fcVal2=arg2*(1.0+arg2*r)*exp_term/(rr*r)-48.0*arg4*rri10-80*arg5*rri12;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case BORN_HUGGINS_MEYER_SMOOTHED5:
      // {p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8}*S(r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      r=sqrt(rr);
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      exp_term=arg1*exp(arg2*(arg3-r));
      U=exp_term-arg4*rri6-arg5*rri8;
      fcVal1=6.0*arg4*rri8+8.0*arg5*rri10-arg2*exp_term/r;
      fcVal2=arg2*(1.0+arg2*r)*exp_term/(rr*r)-48.0*arg4*rri10-80*arg5*rri12;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case HYDROGEN:
      // p_0/r^12-p_1/r^10
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      // p_2/k_B [K]  (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      U=arg1*rri12-arg2*rri10-arg3;
      fcVal1=10.0*arg2*rri12-12.0*arg1*rri14;
      fcVal2=168*arg1*rri16-120*arg2*rri14;
      break;
    case HYDROGEN_SMOOTHED3:
      // {p_0/r^12-p_1/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      U=arg1*rri12-arg2*rri10;
      fcVal1=10.0*arg2*rri12-12.0*arg1*rri14;
      fcVal2=168*arg1*rri16-120*arg2*rri14;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        SwitchingValue=SwitchingVDWFactors3[3]*r3+SwitchingVDWFactors3[2]*rr+SwitchingVDWFactors3[1]*r+SwitchingVDWFactors3[0];
        SwitchingValueDerivative=3.0*SwitchingVDWFactors3[3]*rr+2.0*SwitchingVDWFactors3[2]*r+SwitchingVDWFactors3[1];
        SwitchingValueDerivative2=6.0*SwitchingVDWFactors3[3]*r+2.0*SwitchingVDWFactors3[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case HYDROGEN_SMOOTHED5:
      // {p_0/r^12-p_1/r^10}*S(r)
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rri2=1.0/rr;
      rri4=rri2*rri2;
      rri6=rri4*rri2;
      rri8=rri6*rri2;
      rri10=rri8*rri2;
      rri12=rri10*rri2;
      rri14=rri12*rri2;
      rri16=rri14*rri2;
      U=arg1*rri12-arg2*rri10;
      fcVal1=10.0*arg2*rri12-12.0*arg1*rri14;
      fcVal2=168*arg1*rri16-120*arg2*rri14;
      if(rr>CutOffVDWSwitchSquared)
      {
        r=sqrt(rr);
        r3=rr*r;
        r4=rr*rr;
        r5=r3*rr;
        SwitchingValue=SwitchingVDWFactors5[5]*r5+SwitchingVDWFactors5[4]*r4+SwitchingVDWFactors5[3]*r3+
                       SwitchingVDWFactors5[2]*rr+SwitchingVDWFactors5[1]*r+SwitchingVDWFactors5[0];
        SwitchingValueDerivative=5.0*SwitchingVDWFactors5[5]*r4+4.0*SwitchingVDWFactors5[4]*r3+3.0*SwitchingVDWFactors5[3]*rr+
                                  2.0*SwitchingVDWFactors5[2]*r+SwitchingVDWFactors5[1];
        SwitchingValueDerivative2=20.0*SwitchingVDWFactors5[5]*r3+12.0*SwitchingVDWFactors5[4]*rr+6.0*SwitchingVDWFactors5[3]*r+2.0*SwitchingVDWFactors5[2];

        fcVal2=U*(SwitchingValueDerivative2-SwitchingValueDerivative/r)/rr+2.0*fcVal1*SwitchingValueDerivative/r+fcVal2*SwitchingValue;
        fcVal1=U*SwitchingValueDerivative/r+fcVal1*SwitchingValue;
        U*=SwitchingValue;
      }
      break;
    case LENNARD_JONES_CONTINUOUS_FRACTIONAL:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_3/k_B [K]    (non-zero for a shifted potential)
      r=sqrt(rr);
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2];
      arg4=CurrentLambda;
      rri3=CUBE(arg2/rr);
      term1=0.5*SQR(1.0-arg4);
      fcVal1=-(48.0*arg4*rr*rr*r/SQR(term1+rr*rr*rr))*((1.0/(term1+rr*rr*rr))-0.5);
      break;
    default:
      U=0.0;
      fcVal1=0.0;
      fcVal2=0.0;
      break;
  }
  *energy=U;
  *factor1=fcVal1;
  *factor2=fcVal2;
}

/********************************************************************************************
 * Name       | PotentialThirdDerivative                                                    *
 * ---------------------------------------------------------------------------------------- *
 * Function   | Returns the energy, first, second, and third derivative of a potential U[r] *
 * Parameters | int typeA: the atom-type of the first atom of the pair interaction          *
 *            | int typeB: the atom-type of the second atom of the pair interaction         *
 *            | REAL rr:   the squared Cartesian distance between the atoms                 *
 * Returns    | REAL *Energy:  the contribution to the energy                               *
 *            | REAL *Factor1: the first derivative     D[U[r],r]/r                         *
 *            | REAL *Factor2: the second derivative    D[D[U[r],r]/r,r]/r                  *
 *            | REAL *Factor3: the third derivative     D[D[D[U[r],r]/r,r]/r,r]/r           *
 * Note       | The square is passed to avoid a sqrt when possible (e.g. Lennard-Jones)     *
 ********************************************************************************************/

void PotentialThirdDerivative(int typeA,int typeB,REAL rr,REAL *energy,REAL *factor1,REAL *factor2,REAL *factor3)
{
  REAL fcVal1,fcVal2,fcVal3,U,rri3;
  REAL arg1,arg2,arg3;

  switch(PotentialType[typeA][typeB])
  {
    case UNDEFINED_POTENTIAL:
    case ZERO_POTENTIAL:
      U=0.0;
      fcVal1=0.0;
      fcVal2=0.0;
      fcVal3=0.0;
      break;
    case LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_3/k_B [K]    (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=SQR(PotentialParms[typeA][typeB][1]);
      arg3=PotentialParms[typeA][typeB][2];
      rri3=CUBE(arg2/rr);
      U=4.0*arg1*(rri3*(rri3-1.0))-arg3;
      fcVal1=24.0*arg1*(rri3*(1.0-2.0*rri3))/rr;
      fcVal2=96.0*arg1*(rri3*(7.0*rri3-2.0))/SQR(rr);
      fcVal3=384.0*arg1*(rri3*(5.0-28.0*rri3))/SQR(SQR(rr));
      break;
    default:
      U=0.0;
      fcVal1=0.0;
      fcVal2=0.0;
      fcVal3=0.0;
      break;
  }
  *energy=U;
  *factor1=fcVal1;
  *factor2=fcVal2;
  *factor3=fcVal3;
}



/************************************************************************************************
 * Name       | PotentialCorrection                                                             *
 * -------------------------------------------------------------------------------------------- *
 * Function   | Returns the attractive part of the potential U[r] integrated outside the cutoff *
 *            | int_{r_c}^{infty} r^2 U(r) dr                                                   *
 * Parameters | int typeA: the atom-type of the first atom of the pair interaction              *
 *            | int typeB: the atom-type of the second atom of the pair interaction             *
 *            | REAL r:    the cutoff-distance                                                  *
 * Notes      | Smoothed and shifted potentials do not need a correction                        *
 ************************************************************************************************/

REAL PotentialCorrection(int typeA,int typeB,REAL r)
{     
  REAL arg1,arg2,arg3,arg4,arg5,arg6;
  REAL rr,ri3,ri9;
  REAL term1,term2,term3,term4;

  switch(PotentialType[typeA][typeB])
  {             
    case UNDEFINED_POTENTIAL:
    case ZERO_POTENTIAL:
      return 0.0;
    case LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      term1=CUBE(arg2/r);
      term2=CUBE(term1);
      return (4.0/3.0)*arg1*CUBE(arg2)*((1.0/3.0)*term2-term1);
    case LENNARD_JONES_SMOOTHED3:
    case LENNARD_JONES_SMOOTHED5:
      return 0.0;
    case FEYNMAN_HIBBS_LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      // =============================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [u]    reduced mass in unified atomic mass units
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      term1=CUBE(arg2/r);
      term2=CUBE(term1);
      term3=SQR(term1);
      term4=SQR(term3);
      return (4.0/3.0)*arg1*CUBE(arg2)*(term2/3.0-term1)+arg3*24.0*arg1*r*(2.0*term4-term3);
    case WCA:
      return 0.0;
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3:
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5:
      return 0.0;
    case FEYNMAN_HIBBS_LENNARD_JONES2:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      // ============================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      term1=CUBE(arg2/r);
      term2=CUBE(term1);
      term3=SQR(term1);
      term4=SQR(term3);
      return (4.0/3.0)*arg1*CUBE(arg2)*(term2/3.0-term1)+arg3*24.0*arg1*r*(2.0*term4-term3);
    case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3:
    case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5:
      return 0.0;
    case LENNARD_JONES_SHIFTED_FORCE:
    case LENNARD_JONES_SHIFTED_FORCE2:
      return 0.0;
    case POTENTIAL_12_6:
      // p_0/r^12-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      rr=r*r;
      ri3=1.0/(rr*r);
      ri9=ri3*ri3*ri3;
      return (arg1*ri9-3.0*arg2*ri3)/9.0;
    case POTENTIAL_12_6_SMOOTHED3:
    case POTENTIAL_12_6_SMOOTHED5:
      return 0.0;
    case POTENTIAL_12_6_2_0:
      printf("The 12_6_2_0 potential has no converging tail-correction.\n");
      exit(0);
      break;
    case POTENTIAL_12_6_2_0_SMOOTHED3:
    case POTENTIAL_12_6_2_0_SMOOTHED5:
      break;
    case MORSE:
      // p_0*{(1.0-exp[-p_1*(r-p_2)])^2-1.0}  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      return (arg1*exp(arg2*(arg3 - 2*r))*(exp(arg2*arg3)*(1 + 2*arg2*r*(1 + arg2*r)) - 8*exp(arg2*r)*(2 + arg2*r*(2 + arg2*r))))/(4.*pow(arg2,3));
    case MORSE_SMOOTHED3:
    case MORSE_SMOOTHED5:
      return 0.0;
      break;
    case MORSE2:
      // p_0*[exp{p_1*(1-r/p_2)}-2*exp{(p_1/2)*(1-r/p_2)}]
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      return (arg1*arg3*(exp(arg2)*(2.0*SQR(arg3)+2.0*arg2*arg3*r+SQR(arg2)*SQR(r))- 
             4.0*exp((arg2*(arg3+r))/(2.0*arg3))*(8.0*SQR(arg3)+4.0*arg2*arg3*r+SQR(arg2)*SQR(r))))/
            (CUBE(arg2)*exp((arg2*r)/arg3));
    case MORSE2_SMOOTHED3:
    case MORSE2_SMOOTHED5:
      return 0.0;
      break;
    case MORSE3:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      return (arg1*exp(arg2*(arg3 - 2*r))*(exp(arg2*arg3)*(1 + 2*arg2*r*(1 + arg2*r)) - 8*exp(arg2*r)*(2 + arg2*r*(2 + arg2*r))))/(4.*pow(arg2,3));
    case MORSE3_SMOOTHED3:
    case MORSE3_SMOOTHED5:
      return 0.0;
    case CFF_9_6:
      // p_0/r^9-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      return (arg1-2.0*arg2*CUBE(r))/(6.0*pow(r,6));
    case CFF_9_6_SMOOTHED3:
    case CFF_9_6_SMOOTHED5:
      return 0.0;
    case CFF_EPS_SIGMA:
      // p_0*[2*(p_1/r)^9-3(p_1/r)^6]
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      term1=CUBE(arg2);
      term2=SQR(term1);
      return arg1*term2*(term1-3.0*CUBE(r))/(3.0*pow(r,6));
    case CFF_EPS_SIGMA_SMOOTHED3:
    case CFF_EPS_SIGMA_SMOOTHED5:
      return 0.0;
    case BUCKINGHAM:
    case BUCKINGHAM2:
      // p_0*exp(-p_1*r)-p_2/r^6
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      return arg1*exp(-arg2*r)*(2.0+arg2*r*(2.0+arg2*r))/CUBE(arg2)-arg3/(3.0*CUBE(r));
    case BUCKINGHAM_SMOOTHED3:
    case BUCKINGHAM_SMOOTHED5:
    case BUCKINGHAM2_SMOOTHED3:
    case BUCKINGHAM2_SMOOTHED5:
      break;
    case MM3_VDW:
    case MM3_HYDROGEN_VDW:
      // sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]  if P>=3.02
      // sqrt(p_0^i*p_0^j)*192.27*P^2                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      return (1.0/864.0)*(arg1*arg2*1.84e5*exp(-12.0*r/arg2)*(SQR(arg2)+12.0*arg2*r+72.0*SQR(r)))-arg1*2.25*pow(arg2,6)/(3.0*CUBE(r));
    case MM3_VDW_SMOOTHED3:
    case MM3_VDW_SMOOTHED5:
      return 0.0;
    case MATSUOKA_CLEMENTI_YOSHIMINE:
      // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      return arg1*exp(-arg2*r)*(2.0+arg2*r*(2.0+arg2*r))/CUBE(arg2)+arg3*exp(-arg4*r)*(2.0+arg4*r*(2.0+arg4*r))/CUBE(arg4);
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3:
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5:
      return 0.0;
    case GENERIC:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      return arg1*exp(-arg2*r)*(2.0+arg2*r*(2.0+arg2*r))/CUBE(arg2)-arg3/r-arg4/(3.0*CUBE(r))-arg5/(5.0*pow(r,5))-arg6/(7.0*pow(r,7));
    case GENERIC_SMOOTHED3:
    case GENERIC_SMOOTHED5:
      return 0.0;
    case PELLENQ_NICHOLSON:
      // p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      term1=arg1*exp(-arg2*r)*(2.0+arg2*r*(2.0+arg2*r))/CUBE(arg2);
      term2=(arg3*(-240.0+240.0*exp(arg2*r)-arg2*r*(240.0+arg2*r*(120.0+arg2*r*(38.0+arg2*r*(8.0+arg2*r))))))/(720.0*CUBE(r));
      term3=(arg4*(-8064.0+8064.0*exp(arg2*r)-arg2*r*(8064.0+arg2*r*(4032.0+arg2*r*(1344.0+arg2*r*(336.0+arg2*r*(66.0+arg2*r*(10.0+arg2*r))))))))/(40320.0*pow(r,5));
      term4=(arg5*(-518400.0+518400.0*exp(arg2*r)-arg2*r*(518400.0+arg2*r*(259200.0+ 
             arg2*r*(86400.0+arg2*r*(21600.0+arg2*r*(4320.0+arg2*r*(720.0+arg2*r*(102.0+arg2*r*(12.0+arg2*r))))))))))/(3.6288e6*pow(r,7));
      return term1-exp(-arg2*r)*(term2+term3+term4);
    case PELLENQ_NICHOLSON_SMOOTHED3:
    case PELLENQ_NICHOLSON_SMOOTHED5:
      return 0.0;
    case HYDRATED_ION_WATER:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      return arg1*exp(-arg2*r)*(2.0+arg2*r*(2.0+arg2*r))/CUBE(arg2)-arg3/r-arg4/(3.0*pow(r,3))-arg5/(9.0*pow(r,9));
    case HYDRATED_ION_WATER_SMOOTHED3:
    case HYDRATED_ION_WATER_SMOOTHED5:
      return 0.0;
    case MIE:
      // p_0/r^p_1-p_2/r^p_3
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      return arg1*pow(r,3-arg2)/(arg2-3)-arg3*pow(r,3-arg4)/(arg4-3);
    case MIE_SMOOTHED3:
    case MIE_SMOOTHED5:
      return 0.0;
    case BORN_HUGGINS_MEYER:
      // p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      return (arg1*exp(arg2*(arg3-r))*(2.0+arg2*r*(2.0+arg2*r)))/CUBE(arg2)-
             (3*arg5+5.0*arg4*r*r)/(15.0*pow(r,5));
    case BORN_HUGGINS_MEYER_SMOOTHED3:
    case BORN_HUGGINS_MEYER_SMOOTHED5:
      return 0.0;
    case HYDROGEN:
      // p_0/r^12-p_1/r^10
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      return (7.0*arg1-9.0*arg2*r*r)/(63.0*pow(r,9));
    case HYDROGEN_SMOOTHED3:
    case HYDROGEN_SMOOTHED5:
      return 0.0;
    default:
      printf("Undefined potential in routine 'PotentialCorrection' ('potential.c')\n");
      exit(0);
      break;
  }
  return 0;
}

/************************************************************************************************
 * Name       | PotentialCorrectionPressure                                                     *
 * -------------------------------------------------------------------------------------------- *
 * Function   | Returns the pressure tail-correction integrated outside the cutoff              *
 *            | Integrate[D[U[r],r]*r^3, {r, rc, Infinity}]                                     *
 * Parameters | int typeA: the atom-type of the first atom of the pair interaction              *
 *            | int typeB: the atom-type of the second atom of the pair interaction             *
 *            | REAL r:    the cutoff-distance                                                  *
 * Notes      | Smoothed and shifted potentials do not need a correction                        *
 ************************************************************************************************/


REAL PotentialCorrectionPressure(int typeA,int typeB,REAL r)
{
  REAL arg1,arg2,arg3,arg4,arg5,arg6;
  REAL term1,term2,term3,term4;

  switch(PotentialType[typeA][typeB])
  {
    case UNDEFINED_POTENTIAL:
    case ZERO_POTENTIAL:
      return 0.0;
    case LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
      // ======================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      term1=CUBE(arg2);
      term2=CUBE(arg2/r);
      return 8.0*arg1*term1*(term2-(2.0/3.0)*CUBE(term2));
    case WCA:
      return 0.0;
    case FEYNMAN_HIBBS_LENNARD_JONES:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(h_bar^2/(24 p_2 K_B T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      // =============================================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [u]    reduced mass in unified atomic mass units
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=FH_CONVERSION_FACTOR/(PotentialParms[typeA][typeB][2]*therm_baro_stats.ExternalTemperature[CurrentSystem]);
      term1=CUBE(arg2/r);
      term2=CUBE(term1);
      term3=SQR(term1);
      term4=SQR(term3);
      return 8.0*arg1*CUBE(arg2)*(term1-(2.0/3.0)*CUBE(term1))+arg3*96.0*arg1*r*(-7.0*term4+2.0*term3);
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED3:
    case FEYNMAN_HIBBS_LENNARD_JONES_SMOOTHED5:
      return 0.0;
    case FEYNMAN_HIBBS_LENNARD_JONES2:
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+(p_2/T)*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      // ============================================================================
      // p_0/k_B [K]    strength parameter epsilon
      // p_1     [A]    size parameter sigma
      // p_2     [A^2]  correction factor
      // T       [K]    the temperature
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2]/therm_baro_stats.ExternalTemperature[CurrentSystem];
      term1=CUBE(arg2/r);
      term2=CUBE(term1);
      term3=SQR(term1);
      term4=SQR(term3);
      return 8.0*arg1*CUBE(arg2)*(term1-(2.0/3.0)*CUBE(term1))+arg3*96.0*arg1*r*(-7.0*term4+2.0*term3);
    case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED3:
    case FEYNMAN_HIBBS_LENNARD_JONES2_SMOOTHED5:
      return 0.0;
    case LENNARD_JONES_SHIFTED_FORCE:
    case LENNARD_JONES_SHIFTED_FORCE2:
      return 0.0;
    case POTENTIAL_12_6:
      // p_0/r^12-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      return (-4.0*arg1)/(3.0*pow(r,9))+(2.0*arg2)/pow(r,3);
    case POTENTIAL_12_6_SMOOTHED3:
    case POTENTIAL_12_6_SMOOTHED5:
      return 0.0;
    case POTENTIAL_12_6_2_0:
      printf("The 12_6_2_0 potential has no converging tail-correction.\n");
      exit(0);
      break;
    case POTENTIAL_12_6_2_0_SMOOTHED3:
    case POTENTIAL_12_6_2_0_SMOOTHED5:
      break;
    case MORSE:
      // p_0*{(1.0-exp[-p_1*(r-p_2)])^2-1.0}  (also -p_0*[1.0-{1.0-exp(-p_1(r-p_2))}^2])
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      return (arg1*exp(arg2*(arg3 - 2*r))*(8*exp(arg2*r)*(6 + arg2*r*(6 + arg2*r*(3 + arg2*r))) - 
              exp(arg2*arg3)*(3 + 2*arg2*r*(3 + arg2*r*(3 + 2*arg2*r)))))/(4.*pow(arg2,3));

    case MORSE_SMOOTHED3:
    case MORSE_SMOOTHED5:
      return 0.0;
    case MORSE2:
      // p_0*(exp[p_1*(1-r/p_2)]-2*exp[(p_1/2)*(1-r/p_2)])
      // =================================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      return (-(arg1*exp(arg2)*(6*pow(arg3,3) + 6*arg2*pow(arg3,2)*r + 3*pow(arg2,2)*arg3*pow(r,2) + pow(arg2,3)*pow(r,3))) + 
             2*arg1*exp((arg2*(arg3 + r))/(2.*arg3))*(48*pow(arg3,3) + 24*arg2*pow(arg3,2)*r + 6*pow(arg2,2)*arg3*pow(r,2) + 
             pow(arg2,3)*pow(r,3)))/(pow(arg2,3)*exp((arg2*r)/arg3));
    case MORSE2_SMOOTHED3:
    case MORSE2_SMOOTHED5:
      return 0.0;
    case MORSE3:
      // p_0*{(1-exp[(-ln(2)/(2^(1/6)-1))(r/p_1-2^(1/6))])^2-1}
      // ===============================================================================
      // p_0/k_B [K]       force constant
      // p_1     [A^-1]    parameter
      // p_2     [A]       reference distance
      // p_3/k_B [K]       (non-zero for a shifted potential)
      arg1=PotentialParms[typeA][typeB][0];
      arg2=log(2.0)/(PotentialParms[typeA][typeB][1]*(pow(2.0,1.0/6.0)-1.0));
      arg3=pow(2.0,1.0/6.0)*PotentialParms[typeA][typeB][1];
      return (arg1*exp(arg2*(arg3 - 2*r))*(8*exp(arg2*r)*(6 + arg2*r*(6 + arg2*r*(3 + arg2*r))) - 
              exp(arg2*arg3)*(3 + 2*arg2*r*(3 + arg2*r*(3 + 2*arg2*r)))))/(4.*pow(arg2,3));
    case CFF_9_6:
      // p_0/r^9-p_1/r^6
      // ======================================================================================
      // p_0/k_B [K A^9]
      // p_1/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      term1=1.0/CUBE(r);
      return 2.0*arg2*term1-1.5*arg1*SQR(term1);
    case CFF_9_6_SMOOTHED3:
    case CFF_9_6_SMOOTHED5:
      return 0.0;
    case CFF_EPS_SIGMA:
      // p_0*[2*(p_1/r)^9-3(p_1/r)^6]
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      term1=CUBE(arg2);
      term2=3.0*arg1*SQR(term1);
      term3=1.0/CUBE(r);
      return 2.0*term2*term3-term2*term1*SQR(term3);
    case CFF_EPS_SIGMA_SMOOTHED3:
    case CFF_EPS_SIGMA_SMOOTHED5:
      return 0.0;
    case BUCKINGHAM:
    case BUCKINGHAM2:
      // p_0*exp(-p_1*r)-p_2/r^6
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      return (2.0*arg3)/pow(r,3)-(arg1*(6.0+arg2*r*(6.0+arg2*r*(3.0+arg2*r))))/(pow(arg2,3)*exp(arg2*r));
    case BUCKINGHAM_SMOOTHED3:
    case BUCKINGHAM_SMOOTHED5:
    case BUCKINGHAM2_SMOOTHED3:
    case BUCKINGHAM2_SMOOTHED5:
      break;
    case MM3_VDW:
    case MM3_HYDROGEN_VDW:
      // sqrt(p_0^i*p_0^j)*[1.84e5*exp(-12/P)-2.25*P^6]  if P>=3.02
      // sqrt(p_0^i*p_0^j)*192.27*P^2                    if P<3.02
      // ======================================================================================
      // p_0     [kcal/mol]
      // p_1     [A]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      return (2.0*arg1*pow(arg2,6)*2.25)/pow(r,3)-(1.84e5*arg1*(pow(arg2,3)+12.0*pow(arg2,2)*r+72.0*arg2*pow(r,2)+288.0*pow(r,3)))/
               (288.0*exp((12.0*r)/arg2));
    case MM3_VDW_SMOOTHED3:
    case MM3_VDW_SMOOTHED5:
      return 0.0;
    case MATSUOKA_CLEMENTI_YOSHIMINE:
      // p_0*exp(-p_1*r)+p_2*exp(-p_3*r)
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K]
      // p_3     [A^-1]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      return (exp((-arg2 - arg4)*r)*(-(arg1*pow(arg4,3)*exp(arg4*r)*(6 + arg2*r*(6 + arg2*r*(3 + arg2*r)))) - 
              pow(arg2,3)*arg3*exp(arg2*r)*(6 + arg4*r*(6 + arg4*r*(3 + arg4*r)))))/(pow(arg2,3)*pow(arg4,3));
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED3:
    case MATSUOKA_CLEMENTI_YOSHIMINE_SMOOTHED5:
      return 0.0;
    case GENERIC:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^8-p_5/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1/k_B [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      // p_5/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      arg6=PotentialParms[typeA][typeB][5];
      return -(arg1*(6.0+arg2*r*(6.0+arg2*r*(3.0+arg2*r)))*exp(-arg2*r))/CUBE(arg2)+(4.0*arg3)/r+(2.0*arg4)/pow(r,3)+(8.0*arg5)/(5.0*pow(r,5))+(10.0*arg6)/(7.0*pow(r,7));
    case GENERIC_SMOOTHED3:
    case GENERIC_SMOOTHED5:
      return 0.0;
    case PELLENQ_NICHOLSON:
      // p_0*exp(-p_1*r)-f_6*p_2/r^6-f_8*p_3/r^8-f_10*p_4/r^10
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^6]
      // p_3/k_B [K A^8]
      // p_4/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      term1=-arg1*(6.0+arg2*r*(6.0+arg2*r*(3.0+arg2*r)))*exp(-arg2*r)/CUBE(arg2);
      term2=(arg3*(1440.0-1440.0*exp(arg2*r)+arg2*r*(1440.0+arg2*r*(720.0+arg2*r*(234.0+arg2*r*(54.0+arg2*r*(9.0+arg2*r)))))))/
             (720.0*pow(r,3));
      term3=(arg4*(64512.0-64512.0*exp(arg2*r)+arg2*r*(64512.0+arg2*r*
           (32256.0+arg2*r*(10752.0+arg2*r*(2688.0+arg2*r*(534.0+arg2*r*(86.0+arg2*r*(11.0+arg2*r)))))))))/
            (40320.0*pow(r,5));
      term4=(arg5*(5184000.0-5184000.0*exp(arg2*r)+arg2*r*(5184000.0+arg2*r*(2592000.0+arg2*r*(864000.0+arg2*r*(216000.0+
             arg2*r*(43200.0+arg2*r*(7200.0+arg2*r*(1026.0+arg2*r*(126.0+arg2*r*(13.0+arg2*r)))))))))))/
             (3.6288e6*pow(r,7));
      return term1-exp(-arg2*r)*(term2+term3+term4);
    case PELLENQ_NICHOLSON_SMOOTHED3:
    case PELLENQ_NICHOLSON_SMOOTHED5:
      return 0.0;
    case HYDRATED_ION_WATER:
      // p_0*exp(-p_1*r)-p_2/r^4-p_3/r^6-p_4/r^12
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [A^-1]
      // p_2/k_B [K A^4]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^12]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      return (4.0*arg5)/(3.0*pow(r,9))+(2.0*arg4)/pow(r,3)+(4.0*arg3)/r-exp(-arg2*r)*(arg1*(6.0+arg2*r*(6.0+arg2*r*(3.0+arg2*r))))/CUBE(arg2);
    case HYDRATED_ION_WATER_SMOOTHED3:
    case HYDRATED_ION_WATER_SMOOTHED5:
      return 0.0;
    case MIE:
      // p_0/r^p_1-p_2/r^p_3
      // ======================================================================================
      // p_0/k_B [K A^p_1]
      // p_1     [-]
      // p_2/k_B [K A^p_3]
      // p_3     [-]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      return (arg3*arg4*pow(r,3.0-arg4))/(arg4-3.0)-(arg1*arg2*pow(r,3.0-arg2))/(arg2-3.0);
    case MIE_SMOOTHED3:
    case MIE_SMOOTHED5:
      return 0.0;
    case BORN_HUGGINS_MEYER:
      // p_0*exp(p_1*(p_2-r))-p_3/r^6-p_4/r^8
      // ======================================================================================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [A]
      // p_3/k_B [K A^6]
      // p_4/k_B [K A^8]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      arg3=PotentialParms[typeA][typeB][2];
      arg4=PotentialParms[typeA][typeB][3];
      arg5=PotentialParms[typeA][typeB][4];
      return (8.0*arg5)/(5.0*pow(r,5))+(2.0*arg4)/pow(r,3)-(arg1*exp(arg2*(arg3-r))*(6.0+arg2*r*(6.0+arg2*r*(3.0+arg2*r))))/CUBE(arg2);
    case BORN_HUGGINS_MEYER_SMOOTHED3:
    case BORN_HUGGINS_MEYER_SMOOTHED5:
      return 0.0;
    case HYDROGEN:
      // p_0/r^12-p_1/r^10
      // ======================================================================================
      // p_0/k_B [K A^12]
      // p_1/k_B [K A^10]
      arg1=PotentialParms[typeA][typeB][0];
      arg2=PotentialParms[typeA][typeB][1];
      return  (10.0*arg2)/(7.0*pow(r,7))-4.0*arg1/(3.0*pow(r,9));
    case HYDROGEN_SMOOTHED3:
    case HYDROGEN_SMOOTHED5:
      return 0.0;
    default:
      printf("Undefined potential in routine 'PotentialCorrectionPressure' ('potential.c')\n");
      exit(0);
      break;
  }
  return 0.0;
}

/*********************************************************************************************************
 * Name       | CalculateTailCorrection                                                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the tail-correction contribution for the energy, pressure and strain-derivative *
 * Parameters | -                                                                                        *
 * Note       | This function is used in computing the total energy                                      *
 *********************************************************************************************************/

void CalculateTailCorrection(void)
{
  int i,j;
  REAL energy,pressure;

  energy=pressure=0.0;
  for(i=0;i<NumberOfPseudoAtoms;i++)
    for(j=0;j<NumberOfPseudoAtoms;j++)
    {
      if(TailCorrection[i][j])
      {
        energy+=2.0*M_PI*NumberOfPseudoAtomsType[CurrentSystem][i]*NumberOfPseudoAtomsType[CurrentSystem][j]*PotentialCorrection(i,j,CutOffVDW);
        pressure+=(2.0/3.0)*M_PI*NumberOfPseudoAtomsType[CurrentSystem][i]*NumberOfPseudoAtomsType[CurrentSystem][j]*PotentialCorrectionPressure(i,j,CutOffVDW);
      }
      else if(!ShiftPotential[i][j])
      {
        // impulsive correction
        pressure-=(2.0/3.0)*M_PI*NumberOfPseudoAtomsType[CurrentSystem][i]*NumberOfPseudoAtomsType[CurrentSystem][j]*CUBE(CutOffVDW)*PotentialValue(i,j,CutOffVDWSquared);
      }
    }

  StrainDerivativeTailCorrection[CurrentSystem]=pressure/Volume[CurrentSystem];
  UTailCorrection[CurrentSystem]=energy/Volume[CurrentSystem];
}

REAL CalculatePressureTailCorrection(void)
{
  int i,j;
  REAL pressure;

  pressure=0.0;
  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    for(j=0;j<NumberOfPseudoAtoms;j++)
    {
      if(TailCorrection[i][j])
        pressure+=(2.0/3.0)*M_PI*NumberOfPseudoAtomsType[CurrentSystem][i]*NumberOfPseudoAtomsType[CurrentSystem][j]*PotentialCorrectionPressure(i,j,CutOffVDW);
      else if(!ShiftPotential[i][j])
        pressure-=(2.0/3.0)*M_PI*NumberOfPseudoAtomsType[CurrentSystem][i]*NumberOfPseudoAtomsType[CurrentSystem][j]*CUBE(CutOffVDW)*PotentialValue(i,j,CutOffVDWSquared);
    }
  }
  return pressure/Volume[CurrentSystem];
}


/*********************************************************************************************************
 * Name       | TailMolecularEnergyDifferenceAdd                                                         *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Returns the tail-correction energy difference when adding a molecule                     *
 * Parameters | -                                                                                        *
 * Note       | This function is used in the swap-addition Monte Carlo move                              *
 *********************************************************************************************************/

REAL TailMolecularEnergyDifferenceAdd(void)
{
  int i,j;
  int nr_atoms,type;
  REAL energy;

  for(i=0;i<NumberOfPseudoAtoms;i++)
    NumberOfPseudoAtomsTypeNew[i]=NumberOfPseudoAtomsType[CurrentSystem][i];

  nr_atoms=Components[CurrentComponent].NumberOfAtoms;
  for(i=0;i<nr_atoms;i++)
  {
    type=Components[CurrentComponent].Type[i];
    NumberOfPseudoAtomsTypeNew[type]++;
  }

  energy=0.0;
  for(i=0;i<NumberOfPseudoAtoms;i++)
    for(j=0;j<NumberOfPseudoAtoms;j++)
    {
      if(TailCorrection[i][j])
      {
        energy+=2.0*M_PI*(REAL)NumberOfPseudoAtomsTypeNew[i]*(REAL)NumberOfPseudoAtomsTypeNew[j]*
                PotentialCorrection(i,j,CutOffVDW);
      }
    }

  return energy/Volume[CurrentSystem]-UTailCorrection[CurrentSystem];
}

/*********************************************************************************************************
 * Name       | TailMolecularEnergyDifferenceRemove                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Returns the tail-correction energy difference when removing a molecule                   *
 * Parameters | -                                                                                        *
 * Note       | This function is used in the swap-remove Monte Carlo move                                *
 *********************************************************************************************************/

REAL TailMolecularEnergyDifferenceRemove(void)
{
  int i,j;
  int nr_atoms,type;
  REAL energy;

  for(i=0;i<NumberOfPseudoAtoms;i++)
    NumberOfPseudoAtomsTypeNew[i]=NumberOfPseudoAtomsType[CurrentSystem][i];

  nr_atoms=Components[CurrentComponent].NumberOfAtoms;
  for(i=0;i<nr_atoms;i++)
  {
    type=Components[CurrentComponent].Type[i];
    NumberOfPseudoAtomsTypeNew[type]--;
  }

  energy=0.0;
  for(i=0;i<NumberOfPseudoAtoms;i++)
    for(j=0;j<NumberOfPseudoAtoms;j++)
    {
      if(TailCorrection[i][j])
      {
        energy+=2.0*M_PI*NumberOfPseudoAtomsTypeNew[i]*NumberOfPseudoAtomsTypeNew[j]*
                PotentialCorrection(i,j,CutOffVDW);
      }
    }

  return UTailCorrection[CurrentSystem]-energy/Volume[CurrentSystem];
}

/*********************************************************************************************************
 * Name       | TailMolecularEnergyDifferenceAddRemove                                                   *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Returns the tail-correction energy difference when adding a molecule of a type A and     *
 *            | removing a molecule of type B                                                            *
 * Parameters | int ComponentToAdd:    the component-id of the molecule that is added                    *
 *            | int ComponentToRemove: the component-id of the molecule that is removed                  *
 * Note       | This function is used in e.g. the identity-switch Monte Carlo move                       *
 *********************************************************************************************************/

REAL TailMolecularEnergyDifferenceAddRemove(int ComponentToAdd,int ComponentToRemove)
{
  int i,j;
  int nr_atoms;
  REAL energy_new,energy_old;

  for(i=0;i<NumberOfPseudoAtoms;i++)
    NumberOfPseudoAtomsTypeNew[i]=NumberOfPseudoAtomsType[CurrentSystem][i];

  nr_atoms=Components[ComponentToAdd].NumberOfAtoms;
  for(i=0;i<nr_atoms;i++)
    NumberOfPseudoAtomsTypeNew[Components[ComponentToAdd].Type[i]]++;

  nr_atoms=Components[ComponentToRemove].NumberOfAtoms;
  for(i=0;i<nr_atoms;i++)
    NumberOfPseudoAtomsTypeNew[Components[ComponentToRemove].Type[i]]--;

  energy_new=0.0;
  energy_old=0.0;
  for(i=0;i<NumberOfPseudoAtoms;i++)
    for(j=0;j<NumberOfPseudoAtoms;j++)
    {
      if(TailCorrection[i][j])
      {
        energy_new+=2.0*M_PI*NumberOfPseudoAtomsTypeNew[i]*NumberOfPseudoAtomsTypeNew[j]*
                    PotentialCorrection(i,j,CutOffVDW);

        energy_old+=2.0*M_PI*NumberOfPseudoAtomsType[CurrentSystem][i]*NumberOfPseudoAtomsType[CurrentSystem][j]*
                    PotentialCorrection(i,j,CutOffVDW);
      }
    }
  return (energy_new-energy_old)/Volume[CurrentSystem];
}

/*********************************************************************************************************
 * Name       | BiasingPotential                                                                         *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Returns the biasing energy  for component 'i' at Cartesian position 'pos'                *
 * Parameters | int i:      the component-id                                                             *
 *            | VECTOR pos: the Cartesian position at which the bias is computed                         *
 *********************************************************************************************************/

REAL BiasingPotential(int i,VECTOR pos)
{
  REAL F;
  REAL Derivatives[5];
  VECTOR s;
  VECTOR dr;
  REAL q;

  F=0.0;
  if(Components[i].Biased!=NO_BIASING)
  {
    s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
    s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
    s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;

    // apply boundary condition
    s.x-=(REAL)NINT(s.x);
    s.y-=(REAL)NINT(s.y);
    s.z-=(REAL)NINT(s.z);

    // 's' the fractional position within the first unit cell
    if(s.x<0.0) s.x+=1.0;
    if(s.y<0.0) s.y+=1.0;
    if(s.z<0.0) s.z+=1.0;

    // map the fractional position to a reaction coordinate
    switch(Components[i].BiasingDirection)
    {
      case NO_MAPPING:
        return 0.0;
      case A_MAPPING:
        q=s.x;
        break;
      case B_MAPPING:
        q=s.y;
        break;
      case C_MAPPING:
        q=s.z;
        break;
      case MAP_AB_DIAGONAL:
        dr.x=M_SQRT1_2;
        dr.y=-M_SQRT1_2;
        dr.z=0.0;
        q=((1.0-s.x)*dr.x+(-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_2;
        break;
      case MAP_AC_DIAGONAL:
        dr.x=M_SQRT1_2;
        dr.y=0.0;
        dr.z=-M_SQRT1_2;
        q=((1.0-s.x)*dr.x+(-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_2;
        break;
      case MAP_BC_DIAGONAL:
        dr.x=M_SQRT1_2;
        dr.y=M_SQRT1_2;
        dr.z=0.0;
        q=(s.x*dr.x+s.y*dr.y+s.z*dr.z)*M_SQRT1_2;
        break;
      case MAP_O_AB_DIAGONAL:
        dr.x=M_SQRT1_2;
        dr.y=M_SQRT1_2;
        dr.z=0.0;
        q=(s.x*dr.x+s.y*dr.y+s.z*dr.z)*M_SQRT1_2;
        break;
      case MAP_O_AC_DIAGONAL:
        dr.x=M_SQRT1_2;
        dr.y=0.0;
        dr.z=M_SQRT1_2;
        q=(s.x*dr.x+s.y*dr.y+s.z*dr.z)*M_SQRT1_2;
        break;
      case MAP_O_BC_DIAGONAL:
        dr.x=0.0;
        dr.y=M_SQRT1_2;
        dr.z=M_SQRT1_2;
        q=(s.x*dr.x+s.y*dr.y+s.z*dr.z)*M_SQRT1_2;
        break;
      case MAP_A_BC_DIAGONAL:
        dr.x=M_SQRT1_3;
        dr.y=-M_SQRT1_3;
        dr.z=-M_SQRT1_3;
        q=((1.0-s.x)*dr.x+(-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_3;
        break;
      case MAP_B_AC_DIAGONAL:
        dr.x=-M_SQRT1_3;
        dr.y=M_SQRT1_3;
        dr.z=-M_SQRT1_3;
        q=((-s.x)*dr.x+(1.0-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_3;
        break;
      case MAP_C_AB_DIAGONAL:
        dr.x=-M_SQRT1_3;
        dr.y=-M_SQRT1_3;
        dr.z=M_SQRT1_3;
        q=((-s.x)*dr.x+(-s.y)*dr.y+(1.0-s.z)*dr.z)*M_SQRT1_3;
        break;
      case MAP_O_ABC_DIAGONAL:
        dr.x=-M_SQRT1_3;
        dr.y=-M_SQRT1_3;
        dr.z=-M_SQRT1_3;
        q=((-s.x)*dr.x+(-s.y)*dr.y+(-s.z)*dr.z)*M_SQRT1_3;
        break;
      default:
        printf("Undefined biasing-direction in routine 'BiasingPotential' (potentials.c)\n");
        exit(0);
        break;
    }
  }

  // evaluate the cubic spline at the 1-D fractional coordinate 'q'
  F=EvaluateCubicSpline(Components[i].BiasingFunction,
                    Components[i].BiasingFunction.n,
                    q,
                    Components[i].FreeEnergyXdata,
                    Derivatives);

  // return the biasing free energy
  return F;
}

REAL BiasingPotentialQ(int i,REAL q)
{
  REAL F;
  REAL Derivatives[5];

  F=EvaluateCubicSpline(Components[i].BiasingFunction,
                    Components[i].BiasingFunction.n,
                    q,
                    Components[i].FreeEnergyXdata,
                    Derivatives);
  return F;
}

REAL BiasingPotentialUmbrellaQ(int i,REAL q)
{
  REAL F,Factor;
  REAL Derivatives[5];

  Factor=Components[i].UmbrellaFactor;
  F=EvaluateCubicSpline(Components[i].BiasingFunction,
                    Components[i].BiasingFunction.n,
                    q,
                    Components[i].FreeEnergyXdata,
                    Derivatives);
  return Factor*F;
}

REAL BiasingPotentialDerivatives(int i,REAL q,REAL *Derivatives)
{
  REAL F;

  F=EvaluateCubicSpline(Components[i].BiasingFunction,
                    Components[i].BiasingFunction.n,
                    q,
                    Components[i].FreeEnergyXdata,
                    Derivatives);
  return F;
}

/*********************************************************************************************************
 * Name       | BiasingPotentialUmbrella                                                                 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Returns the biasing energy  for component 'i' at Cartesian position 'pos'                *
 * Parameters | int i:      the component-id                                                             *
 *            | VECTOR pos: the Cartesian position at which the bias is computed                         *
 * Note       | Used for Umbrella sampling in Monte Carlo simulations                                    *
 *********************************************************************************************************/

REAL BiasingPotentialUmbrella(int i, VECTOR pos)
{
  REAL F,Factor;

  F=BiasingPotential(i,pos);
  Factor=Components[i].UmbrellaFactor;
  return exp(Factor*F);
}

/*********************************************************************************************************
 * Name       | BiasingPotentialRuizMontero                                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Returns the biasing energy  for component 'i' at Cartesian position 'pos'                *
 * Parameters | int i:      the component-id                                                             *
 *            | VECTOR pos: the Cartesian position at which the bias is computed                         *
 * Note       | Used for Ruiz-Montero TST in Monte Carlo simulations                                     *
 *********************************************************************************************************/

REAL BiasingPotentialRuizMontero(int i,VECTOR pos)
{
  REAL F,Factor;

  F=BiasingPotential(i,pos);
  Factor=Components[i].RuizMonteroFactor;
  return exp(Factor*F);
}

void ComputeDummyInteractions(void)
{
  int i,j;
  int interaction;

  for(i=0;i<NumberOfPseudoAtoms;i++)
    PseudoAtoms[i].Interaction=TRUE;

  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    interaction=FALSE;
    for(j=0;j<NumberOfPseudoAtoms;j++)
      interaction|=(PotentialType[i][j]!=NONE);
    PseudoAtoms[i].Interaction=interaction;
  }
}

/*********************************************************************************************************
 * Name       | SwitchingFunctionBend                                                                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the switching funtion for bending near 180 degrees                              *
 * Parameters | REAL theta: the bending angle                                                            *
 * Note       | This function is used to e.g. avoid bending angles of 180 degrees in bend-torsion        *
 *            | zeolite potentials (Nicholas model)                                                      *
 *********************************************************************************************************/

REAL SwitchingFunctionBend(REAL theta)
{
  REAL on,off;

  on=170.0*DEG2RAD;
  off=180.0*DEG2RAD;

  if(theta<on) return 1.0;
  return SQR(off-theta)*(off+2.0*theta-3.0*on)/CUBE(off-on);
}

REAL SwitchingFunctionBendDerivative(REAL theta)
{
  REAL on,off;

  on=170.0*DEG2RAD;
  off=180.0*DEG2RAD;

  if(theta<on) return 0.0;
  return 6.0*(off-theta)*(on-theta)/CUBE(off-on);
}

// TitleTitle: COMPUTER-SIMULATION OF POLAR LIQUIDS 
// Author(s): ADAMS DJ, ADAMS EM, HILLS GJ 
// Source: MOLECULAR PHYSICS   Volume: 38   Issue: 2   Pages: 387-400   Published: 1979 

// Title: THE ROLE OF LONG RANGED FORCES IN DETERMINING THE STRUCTURE AND PROPERTIES OF LIQUID WATER 
// Author(s): ANDREA TA, SWOPE WC, ANDERSEN HC 
// Source: JOURNAL OF CHEMICAL PHYSICS   Volume: 79   Issue: 9   Pages: 4576-4584   Published: 1983 

// Title: NEW SPHERICAL-CUTOFF METHODS FOR LONG-RANGE FORCES IN MACROMOLECULAR SIMULATION 
// Author(s): STEINBACH PJ, BROOKS BR
// Source: JOURNAL OF COMPUTATIONAL CHEMISTRY   Volume: 15   Issue: 7   Pages: 667-683   Published: JUL 1994 


/*********************************************************************************************************
 * Name       | ComputeSwitchingFactors                                                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Computes the factors used in energy switching                                            *
 * Parameters | -                                                                                        *
 * Note       | Polynomial energy  switching over a window  is used for terms whose  energy is small     *
 *            | near the cutoff  distance. For monopole electrostatic  interactions, which are  quite    *
 *            | large in typical  cutoff ranges, a two  polynomial multiplicative-additive shifted       *
 *            | energy switch unique to TINKER is applied. The  TINKER  method  is  similar  in  spirit  *
 *            | to  the  force switching methods of Steinbach and Brooks, J. Comput. Chem., 15, 667-683  *
 *            | (1994).                                                                                  *
 *********************************************************************************************************/

void ComputeSwitchingFactors(void)
{
  REAL temp;
  REAL off,off2,off3,off4,off5,off6,off7;
  REAL cut,cut2,cut3,cut4,cut5,cut6,cut7;

  off=CutOffVDW;
  off2=SQR(off);
  off3=off2*off;
  off4=off3*off;
  off5=off4*off;
  off6=off5*off;
  off7=off6*off;

  cut=CutOffVDWSwitch;
  cut2=SQR(cut);
  cut3=cut2*cut;
  cut4=cut3*cut;
  cut5=cut4*cut;
  cut6=cut5*cut;
  cut7=cut6*cut;

  temp=pow(off-cut,3.0);
  SwitchingVDWFactors3[0]=off2*(off-3.0*cut)/temp;
  SwitchingVDWFactors3[1]=6.0*off*cut/temp;
  SwitchingVDWFactors3[2]=-3.0*(off+cut)/temp;
  SwitchingVDWFactors3[3]=2.0/temp;


  temp=pow(off-cut,5.0);
  SwitchingVDWFactors5[0]=off*off2*(off2-5.0*off*cut+10.0*cut2)/temp;
  SwitchingVDWFactors5[1]=-30.0*off2*cut2/temp;
  SwitchingVDWFactors5[2]=30.0*(off2*cut+off*cut2)/temp;
  SwitchingVDWFactors5[3]=-10.0*(off2+4.0*off*cut+cut2)/temp;
  SwitchingVDWFactors5[4]=15.0*(off+cut)/temp;
  SwitchingVDWFactors5[5]=-6.0/temp;


  temp=(9.3*cut*off/(off-cut))*
        (cut7-7.0*cut6*off+21.0*cut5*off2
         - 35.0*cut4*off3+35.0*cut3*off4
         - 21.0*cut2*off5+7.0*cut*off6-off7);
  SwitchingVDWFactors7[0]=cut3*off3*(-39.0*cut+64.0*off)/temp;
  SwitchingVDWFactors7[1]=cut2*off2*(117.0*cut2-100.0*cut*off-192.0*off2)/temp;
  SwitchingVDWFactors7[2]=cut*off*(-117.0*cut3-84.0*cut2*off+534.0*cut*off2+192.0*off3)/temp;
  SwitchingVDWFactors7[3]=(39.0*cut4+212.0*cut3*off-450.0*cut2*off2-612.0*cut*off3-64.0*off4)/temp;
  SwitchingVDWFactors7[4]=(-92.0*cut3+66.0*cut2*off+684.0*cut*off2+217.0*off3)/temp;
  SwitchingVDWFactors7[5]=(42.0*cut2-300.0*cut*off-267.0*off2)/temp;
  SwitchingVDWFactors7[6]=(36.0*cut+139.0*off)/temp;
  SwitchingVDWFactors7[7]=-25.0/temp;

  off=CutOffChargeCharge;
  off2=SQR(off);
  off3=off2*off;
  off4=off3*off;
  off5=off4*off;
  off6=off5*off;
  off7=off6*off;

  cut=CutOffChargeChargeSwitch;
  cut2=SQR(cut);
  cut3=cut2*cut;
  cut4=cut3*cut;
  cut5=cut4*cut;
  cut6=cut5*cut;
  cut7=cut6*cut;

  temp=pow(off-cut,3.0);
  SwitchingChargeChargeFactors3[0]=off2*(off-3.0*cut)/temp;
  SwitchingChargeChargeFactors3[1]=6.0*off*cut/temp;
  SwitchingChargeChargeFactors3[2]=-3.0*(off+cut)/temp;
  SwitchingChargeChargeFactors3[3]=2.0/temp;


  temp=pow(off-cut,5.0);
  SwitchingChargeChargeFactors5[0]=off*off2*(off2-5.0*off*cut+10.0*cut2)/temp;
  SwitchingChargeChargeFactors5[1]=-30.0*off2*cut2/temp;
  SwitchingChargeChargeFactors5[2]=30.0*(off2*cut+off*cut2)/temp;
  SwitchingChargeChargeFactors5[3]=-10.0*(off2+4.0*off*cut+cut2)/temp;
  SwitchingChargeChargeFactors5[4]=15.0*(off+cut)/temp;
  SwitchingChargeChargeFactors5[5]=-6.0/temp;

  temp=(9.3*cut*off /(off-cut))*
          (cut7-7.0*cut6*off+21.0*cut5*off2
           -35.0*cut4*off3+35.0*cut3*off4
           -21.0*cut2*off5+7.0*cut*off6-off7);
  SwitchingChargeChargeFactors7[0]=cut3*off3*(-39.0*cut+64.0*off)/temp;
  SwitchingChargeChargeFactors7[1]=cut2*off2*(117.0*cut2-100.0*cut*off-192.0*off2)/temp;
  SwitchingChargeChargeFactors7[2]=cut*off*(-117.0*cut3-84.0*cut2*off+534.0*cut*off2+192.0*off3)/temp;
  SwitchingChargeChargeFactors7[3]=(39.0*cut4+212.0*cut3*off-450.0*cut2*off2-612.0*cut*off3-64.0*off4)/temp;
  SwitchingChargeChargeFactors7[4]=(-92.0*cut3+66.0*cut2*off+684.0*cut*off2+217.0*off3)/temp;
  SwitchingChargeChargeFactors7[5]=(42.0*cut2-300.0*cut*off-267.0*off2)/temp;
  SwitchingChargeChargeFactors7[6]=(36.0*cut+139.0*off)/temp;
  SwitchingChargeChargeFactors7[7]= -25.0/temp;


  off=CutOffChargeBondDipole;
  off2=SQR(off);
  off3=off2*off;
  off4=off3*off;
  off5=off4*off;
  off6=off5*off;
  off7=off6*off;

  cut=CutOffChargeBondDipoleSwitch;
  cut2=SQR(cut);
  cut3=cut2*cut;
  cut4=cut3*cut;
  cut5=cut4*cut;
  cut6=cut5*cut;
  cut7=cut6*cut;

  temp=pow(off-cut,3.0);
  SwitchingChargeBondDipoleFactors3[0]=off2*(off-3.0*cut)/temp;
  SwitchingChargeBondDipoleFactors3[1]=6.0*off*cut/temp;
  SwitchingChargeBondDipoleFactors3[2]=-3.0*(off+cut)/temp;
  SwitchingChargeBondDipoleFactors3[3]=2.0/temp;


  temp=pow(off-cut,5.0);
  SwitchingChargeBondDipoleFactors5[0]=off*off2*(off2-5.0*off*cut+10.0*cut2)/temp;
  SwitchingChargeBondDipoleFactors5[1]=-30.0*off2*cut2/temp;
  SwitchingChargeBondDipoleFactors5[2]=30.0*(off2*cut+off*cut2)/temp;
  SwitchingChargeBondDipoleFactors5[3]=-10.0*(off2+4.0*off*cut+cut2)/temp;
  SwitchingChargeBondDipoleFactors5[4]=15.0*(off+cut)/temp;
  SwitchingChargeBondDipoleFactors5[5]=-6.0/temp;

  temp=(9.3*cut*off /(off-cut))*
          (cut7-7.0*cut6*off+21.0*cut5*off2
           -35.0*cut4*off3+35.0*cut3*off4
           -21.0*cut2*off5+7.0*cut*off6-off7);
  SwitchingChargeBondDipoleFactors7[0]=cut3*off3*(-39.0*cut+64.0*off)/temp;
  SwitchingChargeBondDipoleFactors7[1]=cut2*off2*(117.0*cut2-100.0*cut*off-192.0*off2)/temp;
  SwitchingChargeBondDipoleFactors7[2]=cut*off*(-117.0*cut3-84.0*cut2*off+534.0*cut*off2+192.0*off3)/temp;
  SwitchingChargeBondDipoleFactors7[3]=(39.0*cut4+212.0*cut3*off-450.0*cut2*off2-612.0*cut*off3-64.0*off4)/temp;
  SwitchingChargeBondDipoleFactors7[4]=(-92.0*cut3+66.0*cut2*off+684.0*cut*off2+217.0*off3)/temp;
  SwitchingChargeBondDipoleFactors7[5]=(42.0*cut2-300.0*cut*off-267.0*off2)/temp;
  SwitchingChargeBondDipoleFactors7[6]=(36.0*cut+139.0*off)/temp;
  SwitchingChargeBondDipoleFactors7[7]= -25.0/temp;

  off=CutOffBondDipoleBondDipole;
  off2=SQR(off);
  off3=off2*off;
  off4=off3*off;
  off5=off4*off;
  off6=off5*off;
  off7=off6*off;

  cut=CutOffBondDipoleBondDipoleSwitch;
  cut2=SQR(cut);
  cut3=cut2*cut;
  cut4=cut3*cut;
  cut5=cut4*cut;
  cut6=cut5*cut;
  cut7=cut6*cut;

  temp=pow(off-cut,3.0);
  SwitchingBondDipoleBondDipoleFactors3[0]=off2*(off-3.0*cut)/temp;
  SwitchingBondDipoleBondDipoleFactors3[1]=6.0*off*cut/temp;
  SwitchingBondDipoleBondDipoleFactors3[2]=-3.0*(off+cut)/temp;
  SwitchingBondDipoleBondDipoleFactors3[3]=2.0/temp;


  temp=pow(off-cut,5.0);
  SwitchingBondDipoleBondDipoleFactors5[0]=off*off2*(off2-5.0*off*cut+10.0*cut2)/temp;
  SwitchingBondDipoleBondDipoleFactors5[1]=-30.0*off2*cut2/temp;
  SwitchingBondDipoleBondDipoleFactors5[2]=30.0*(off2*cut+off*cut2)/temp;
  SwitchingBondDipoleBondDipoleFactors5[3]=-10.0*(off2+4.0*off*cut+cut2)/temp;
  SwitchingBondDipoleBondDipoleFactors5[4]=15.0*(off+cut)/temp;
  SwitchingBondDipoleBondDipoleFactors5[5]=-6.0/temp;

  temp=(9.3*cut*off /(off-cut))*
          (cut7-7.0*cut6*off+21.0*cut5*off2
           -35.0*cut4*off3+35.0*cut3*off4
           -21.0*cut2*off5+7.0*cut*off6-off7);
  SwitchingBondDipoleBondDipoleFactors7[0]=cut3*off3*(-39.0*cut+64.0*off)/temp;
  SwitchingBondDipoleBondDipoleFactors7[1]=cut2*off2*(117.0*cut2-100.0*cut*off-192.0*off2)/temp;
  SwitchingBondDipoleBondDipoleFactors7[2]=cut*off*(-117.0*cut3-84.0*cut2*off+534.0*cut*off2+192.0*off3)/temp;
  SwitchingBondDipoleBondDipoleFactors7[3]=(39.0*cut4+212.0*cut3*off-450.0*cut2*off2-612.0*cut*off3-64.0*off4)/temp;
  SwitchingBondDipoleBondDipoleFactors7[4]=(-92.0*cut3+66.0*cut2*off+684.0*cut*off2+217.0*off3)/temp;
  SwitchingBondDipoleBondDipoleFactors7[5]=(42.0*cut2-300.0*cut*off-267.0*off2)/temp;
  SwitchingBondDipoleBondDipoleFactors7[6]=(36.0*cut+139.0*off)/temp;
  SwitchingBondDipoleBondDipoleFactors7[7]= -25.0/temp;
}
