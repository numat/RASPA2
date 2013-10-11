/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'input.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef INPUT_H
#define INPUT_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern int EwaldAutomatic;

char *ReadLine(char *buffer, size_t length, FILE *file);
char *LoadFile(char *);

int ReadInputFile(char *);
int ReadInput(char *);
void ReadRestartFile(void);
void ReadRestartFileOld(void);
void ReadBinaryRestartFiles(void);

#endif
