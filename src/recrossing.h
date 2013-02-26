/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'recrossing.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef RECROSSING_H
#define RECROSSING_H

extern VECTOR *BarrierPosition;
extern VECTOR *BarrierNormal;
extern VECTOR *BarrierAngle;
extern REAL *MaxBarrierDistance;
extern REAL *MaxBarrierTime;
extern int *NumberOfVelocities;
extern int *PutMoleculeOnBarrier;

void SetBarierNormal(void);
int BarrierRecrossing(void);
void AllocateRecrossingMemory(void);

#endif

