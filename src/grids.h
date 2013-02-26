/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'grids.c' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef GRIDS_H
#define GRIDS_H

extern int UseTabularGrid;
extern int NumberOfGrids;
extern int *GridTypeList;
extern char (*GridTypeListName)[256];

extern float *****VDWGrid;
extern float ****CoulombGrid;

extern INT_VECTOR3 NumberOfVDWGridPoints;

extern REAL SpacingVDWGrid;
extern REAL SpacingCoulombGrid;

extern int BlockEnergyGrids;
extern int BlockGridPockets;
extern int BlockGridPores;
extern REAL BlockEnergyGridOverlapCriteria;
extern int NumberOfGridSeeds;
extern VECTOR *GridSeeds;

void MakeASCIGrid(void);


VECTOR MapToUnitCell(VECTOR pos);
VECTOR MapZToBox(VECTOR pos);
void MakeRigidFrameworkList(void);
void RigidFrameworkGrid(VECTOR pos,int typeA,REAL *Uvdw,REAL *Ucoul);

void MakeGrid(void);
int WriteVDWGrid(int l);
void ReadVDWGrid(void);
int WriteCoulombGrid(void);
void ReadCoulombGrid(void);
REAL InterpolateVDWGrid(int typeA,VECTOR pos);
REAL InterpolateVDWForceGrid(int typeA,VECTOR pos,VECTOR *Force);
REAL InterpolateCoulombGrid(int typeA,VECTOR pos);
REAL InterpolateCoulombForceGrid(int typeA,VECTOR pos,VECTOR *Force);
void TestGrid(FILE *FilePtr);
void TestForceGrid(FILE *FilePtr);
INT_VECTOR3 ConvertXYZPositionToGridIndex(VECTOR pos);
VECTOR ConvertGridIndexToXYZIndex(INT_VECTOR3 GridIndex);
void BlockingVDWGrid(void);

void WriteRestartGrids(FILE *FilePtr);
void AllocateGridMemory(void);
void ReadRestartGrids(FILE *FilePtr);
#endif
