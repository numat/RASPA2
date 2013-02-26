/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2013 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'movies.h' is part of RASPA-2.0

 *************************************************************************************************************/

#ifndef MOVIES_H
#define MOVIES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

enum{PDB,CSSR,MOL,XYZ,VASP};
enum{VTK_UNIT_CELL,VTK_FULL_BOX};

extern REAL MovieScale;

extern int WriteVTKGrids;

extern VECTOR VTKFractionalFrameworkAtomsMin;
extern VECTOR VTKFractionalFrameworkAtomsMax;
extern VECTOR VTKFractionalFrameworkBondsMin;
extern VECTOR VTKFractionalFrameworkBondsMax;
extern VECTOR VTKFractionalAdsorbateComMin;
extern VECTOR VTKFractionalAdsorbateComMax;
extern VECTOR VTKFractionalCationComMin;
extern VECTOR VTKFractionalCationComMax;

extern int FreeEnergyAveragingTypeVTK;
extern int DensityAveragingTypeVTK;
extern int AverageDensityOverUnitCellsVTK;

extern int *Movies;
extern int *WriteMoviesEvery;

VECTOR ConvertPositionToVTKPosition(VECTOR pos);

int SamplePDBMovies(int Choice,int Subdir);

void WriteVTK(int system);

void WriteMolecule(int mol);
void WriteIonCssr(void);

void WriteSnapshotIonsCssr(void);
void WriteSnapshotIonsCssrUsingSymmetry(void);

void WriteAsymetricPositions(void);
void WriteDlpolyInputFiles(void);

void FreeEnergyProfile3D(void);

void AllocateMovieMemory(void);

void WriteRestartMovies(FILE *FilePtr);
void ReadRestartMovies(FILE *FilePtr);

void WriteSnapshotTinker(char *string);


#endif
