/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'output.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include "simulation.h"

extern FILE **OutputFilePtr;

void OpenOutputFile(void);
void CloseOutputFile(void);

void PrintPreSimulationStatus(void);
void PrintPreSimulationStatusCurrentSystem(int system);

void PrintPostSimulationStatus(void);

void PrintEnergyStatus(FILE *FilePtr);

void PrintEnergyDriftStatus(FILE *FilePtr);

void PrintRestartFile(void);

void WriteBinaryRestartFiles(void);

void ReadRestartOutput(FILE* FilePtr);
void WriteRestartOutput(FILE* FilePtr);

#endif
