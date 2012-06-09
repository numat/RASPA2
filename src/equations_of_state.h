/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'equation_of_state.h' is part of RASPA.

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

#ifndef EQUATIONS_OF_STATE_H
#define EQUATIONS_OF_STATE_H

#include <stdio.h>

enum{VAN_DER_WAALS,REDLICH_KWONG,SOAVE_REDLICH_KWONG,PENG_ROBINSON,PENG_ROBINSON_GASEM};
enum{VAN_DER_WAALS_MIXING_RULES};
enum{VAPOUR_STABLE,LIQUID_STABLE,VAPOUR_LIQUID_STABLE,SUPER_CRITICAL_FLUID,LIQUID,VAPOUR};

extern REAL *HeliumVoidFraction;
extern REAL *ExcessVolume;

extern int EquationOfState;
extern int MultiComponentMixingRules;
extern REAL **BinaryInteractionParameter;

extern int *FluidState;
extern int **ComputeFugacityCoefficient;

void RescaleMolarFractions(void);
void ComputeGasPropertiesForAllSystems(void);

void WriteRestartEquationOfState(FILE *FilePtr);
void AllocateEquationOfStateMemory(void);
void ReadRestartEquationOfState(FILE *FilePtr);

#endif
