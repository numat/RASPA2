/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'phonon_dispersion.c' is part of RASPA.

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
//#include <Accelerate/Accelerate.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "constants.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "framework_hessian.h"
#include "framework_phonon.h"
#include "inter_force.h"
#include "internal_force.h"
#include "internal_hessian.h"
#include "internal_phonon.h"
#include "inter_hessian.h"
#include "inter_phonon.h"
#include "molecule.h"
#include "potentials.h"
#include "integration.h"
#include "matrix.h"
#include "thermo_baro_stats.h"
#include "ewald.h"
#include "numerical.h"
#include "statistics.h"
#include "output.h"
#include "movies.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "scattering_factors.h"
#include "rigid.h"
#include "spectra.h"
#include "minimization.h"

static int NumberOfPositionVariables;
static int NumberOfBoxVariables;
static int NumberOfVariables;

int NumberOfDispersionCurvePoints;
VECTOR WaveVector;

// Compute the energy, generalized gradient, generalized Hessian matrix
// This routine expects that state 'x' is expanded into the positions and cell-info
void EvaluateDerivativesPhonon(VECTOR k,int n,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX Hessian,REAL_MATRIX3x3 *StrainFirstDerivative,
                         int ComputeGradient,int ComputeHessian)
{
  REAL ExtPressure;

  *Energy=0.0;
  StrainFirstDerivative->ax=0.0; StrainFirstDerivative->bx=0.0; StrainFirstDerivative->cx=0.0;
  StrainFirstDerivative->ay=0.0; StrainFirstDerivative->by=0.0; StrainFirstDerivative->cy=0.0;
  StrainFirstDerivative->az=0.0; StrainFirstDerivative->bz=0.0; StrainFirstDerivative->cz=0.0;

  ExtPressure=therm_baro_stats.ExternalPressure[CurrentSystem][0];

  // compute the first and second derivatives of the rotation matrix
  PreComputeRotationDerivatives();

  CalculateAnisotropicSites();

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    ComputeFrameworkBondPhonon(k,Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
    //ComputeFrameworkBendPhonon(k,Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
    ComputeFrameworkIntraVDWPhonon(k,Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
    ComputeFrameworkIntraChargeChargePhonon(k,Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  }

  // compute the adsorbate-adsorbate, cation-cation and adsorbate-cations intermolecular contributions
  ComputeInterVDWMolecularPhonon(k,Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  ComputeInterChargeChargeMolecularPhonon(k,Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);
  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    CalculateEwaldFourierPhonon(k,Energy,Gradient,Hessian,StrainFirstDerivative,ComputeGradient,ComputeHessian);

}

void ComputeDerivativesPhonon(VECTOR k,int np,int nb,REAL *x,REAL* Energy,REAL *Gradient,COMPLEX_MATRIX Hessian,REAL_MATRIX3x3 *StrainFirstDerivative)
{
  int i,j;

  *Energy=0.0;
  for(i=0;i<np+nb;i++)
  {
    Gradient[i]=0.0;
    for(j=0;j<np+nb;j++)
      Hessian.element[i][j].re=Hessian.element[i][j].im=0.0;
  }

  StrainFirstDerivative->ax=0.0; StrainFirstDerivative->bx=0.0; StrainFirstDerivative->cx=0.0;
  StrainFirstDerivative->ay=0.0; StrainFirstDerivative->by=0.0; StrainFirstDerivative->cy=0.0;
  StrainFirstDerivative->az=0.0; StrainFirstDerivative->bz=0.0; StrainFirstDerivative->cz=0.0;

  // generate the positions of all the atoms from array 'x'
  CreatePositionsFromGeneralizedCoordinates(np,nb,x);

  // correct the positions for constraints
  ShakeInMinimization();

  CreateGeneralizedCoordinatesFromPositions(np,nb,x);

  EvaluateDerivativesPhonon(k,np,Energy,Gradient,Hessian,StrainFirstDerivative,TRUE,TRUE);

  for(i=0;i<np+nb;i++)
    for(j=i+1;j<np+nb;j++)
    {
      Hessian.element[j][i].re=Hessian.element[i][j].re;
      Hessian.element[j][i].im=-Hessian.element[i][j].im;
    }
}

void MassWeightPhononMatrix(int n,COMPLEX_MATRIX Hessian,REAL *Weights)
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
    {
      Hessian.element[i][j].re*=(Weights[i]*Weights[j]);
      Hessian.element[i][j].im*=(Weights[i]*Weights[j]);
    }
}

/*********************************************************************************************************
 * Name       | RemoveShellInteractionsFromHessian                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Remove the presence of the massless shell from the Hessian matrix                        *
 * Parameters | -                                                                                        *
 * Note       | The number of degrees of freedom is the same as for the rigid ion model. The corrected   *
 *            | Hessian matrix where the shell interactions have been removed can be written as:         *
 *            | H'= H_cc - H_cs H_ss^-1 H_sc                                                             *
 *            | cc is core-core, cs is core-shell, sc is shell-core, and ss is the shell-shell part      *
 *********************************************************************************************************/

void RemoveShellInteractionsFromPhonons(COMPLEX_MATRIX Hessian,COMPLEX_MATRIX *CorrectedHessian)
{
  int i,j,k,l;
  COMPLEX_MATRIX ss;

  if(ShellSize>0)
  {
    ss=CreateComplexMatrix(ShellSize,ShellSize);

    for(i=0;i<ShellSize;i++)
      for(j=0;j<ShellSize;j++)
      {
        ss.element[i][j].re=Hessian.element[i+ShellIndex][j+ShellIndex].re;
        ss.element[i][j].im=Hessian.element[i+ShellIndex][j+ShellIndex].im;
      }

    InverseComplexMatrix(ss);

    for(i=0;i<CoreSize;i++)
    {
      for(j=0;j<CoreSize;j++)
      {
        CorrectedHessian->element[i][j].re=Hessian.element[i][j].re;
        CorrectedHessian->element[i][j].im=Hessian.element[i][j].im;
        for(k=0;k<ShellSize;k++)
        {
          for (l=0;l<ShellSize;l++)
          {
            CorrectedHessian->element[i][j].re-=Hessian.element[i][ShellIndex+l].re*ss.element[l][k].re*Hessian.element[ShellIndex+k][j].re-
                                                Hessian.element[i][ShellIndex+l].im*ss.element[l][k].im*Hessian.element[ShellIndex+k][j].re-
                                                Hessian.element[i][ShellIndex+l].im*ss.element[l][k].re*Hessian.element[ShellIndex+k][j].im-
                                                Hessian.element[i][ShellIndex+l].re*ss.element[l][k].im*Hessian.element[ShellIndex+k][j].im;
            CorrectedHessian->element[i][j].im-=Hessian.element[i][ShellIndex+l].im*ss.element[l][k].re*Hessian.element[ShellIndex+k][j].re+
                                                Hessian.element[i][ShellIndex+l].re*ss.element[l][k].im*Hessian.element[ShellIndex+k][j].re+
                                                Hessian.element[i][ShellIndex+l].re*ss.element[l][k].re*Hessian.element[ShellIndex+k][j].im-
                                                Hessian.element[i][ShellIndex+l].im*ss.element[l][k].im*Hessian.element[ShellIndex+k][j].im;
          }
        }
      }
    }
    DeleteComplexMatrix(ss);
  }
  else
  {
    for(i=0;i<Hessian.m;i++)
      for(j=0;j<Hessian.n;j++)
      {
        CorrectedHessian->element[i][j].re=Hessian.element[i][j].re;
        CorrectedHessian->element[i][j].im=Hessian.element[i][j].im;
      }
  }
}


#ifdef HAVE_LAPACK

extern void zheev_(char *jobz, char *uplo, int *n, COMPLEX *a,
           int *lda, REAL *w, COMPLEX *work,int *lwork,double *rwork,
           int *info);

// eigenvector i: [0][i],[1][i],[2][i], etc.

void SolveEigenValuesAndVectorsPhonon(COMPLEX_MATRIX HessianMatrix,REAL* Frequencies)
{
  char jobz, uplo;
  int n,lda, lwork, info;
  double  *rwork;
  COMPLEX *work;

  info=0;
  jobz = 'V';
  uplo = 'U';
  n=HessianMatrix.n;
  lda = n; // The leading dimension of the matrix to be solved.

  work=(COMPLEX *)malloc(1*sizeof(COMPLEX));
  rwork=(double *)malloc(1*sizeof(double));
  lwork=-1;

  // Workspace size computation
  zheev_(&jobz, &uplo, &n, &HessianMatrix.element[0][0], &lda, Frequencies, work, &lwork, rwork,&info);
  lwork=(int)work[0].re;
  free(work);
  free(rwork);
  work=(COMPLEX *)malloc(lwork*sizeof(COMPLEX));
  rwork=(double *)malloc((n*lwork+1)*sizeof(double));

  // compute eigenvalues and vectors
  zheev_(&jobz, &uplo, &n, &HessianMatrix.element[0][0], &lda, Frequencies, work, &lwork, rwork,&info);
  free(work);
  free(rwork);
}

#else

void SolveEigenValuesAndVectorsPhonon(COMPLEX_MATRIX HessianMatrix,REAL* Frequencies)
{
}

#endif

int PhononDisperionCurves(void)
{
  int i,j;
  REAL Energy;
  REAL *Gradients;
  REAL *Weights;
  REAL *Positions;
  REAL *Charges;
  REAL *Frequencies;
  COMPLEX_MATRIX ReducedGeneralizedHessianMatrix;
  COMPLEX_MATRIX GeneralizedHessianMatrix;
  REAL_MATRIX3x3 StrainDerivativeTensor;
  REAL length,max_length;
  VECTOR k;
  FILE *FilePtr;
  char buffer[256];

  printf("Computing phonon dipersion curve\n");

  WaveVector.x=0.0; WaveVector.y=0.0; WaveVector.z=1.0;
  NumberOfDispersionCurvePoints=1000;

  CurrentSystem=0;
  StoredBox=Box[CurrentSystem];
  StoredReplicaBox=ReplicaBox[CurrentSystem];
  StoredInverseBox=InverseBox[CurrentSystem];

  MinimizationVariables=CARTESIAN;
  Ensemble[CurrentSystem]=NVE;

  // index the molecule-atoms into a single list
  NumberOfPositionVariables=OrderNumberOfMinimiationVariables();
  NumberOfBoxVariables=0;   // the box matrix is fixed
  NumberOfVariables=NumberOfPositionVariables+NumberOfBoxVariables;

  GeneralizedHessianMatrix=CreateComplexMatrix(CoreSize+ShellSize,CoreSize+ShellSize);
  ReducedGeneralizedHessianMatrix=CreateComplexMatrix(CoreSize,CoreSize);

  Gradients=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Weights=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Positions=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Charges=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Frequencies=(REAL*)calloc(NumberOfVariables,sizeof(REAL));

  AllocateMinimizationLocalMemory();

  for(i=0;i<NumberOfVariables;i++)
    Gradients[i]=0.0;

  SetWeights(NumberOfCoordinatesMinimizationVariables,Weights,Charges);
  SetStrainToZero(NumberOfCoordinatesMinimizationVariables,Positions);
  CreateGeneralizedCoordinatesFromPositions(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,Positions);

  k.x=4.0*M_PI*WaveVector.x;
  k.y=4.0*M_PI*WaveVector.y;
  k.z=4.0*M_PI*WaveVector.z;
  k=ConvertFromXYZtoABC(k);
  max_length=sqrt(SQR(k.x)+SQR(k.y)+SQR(k.z));

  // make the output directory
  mkdir("PhononDisperion",S_IRWXU);

  // make the system directory in the output directory
  sprintf(buffer,"PhononDisperion/System_%d",CurrentSystem);
  mkdir(buffer,S_IRWXU);

  sprintf(buffer,"PhononDisperion/System_%d/PhononDisperion.dat",CurrentSystem);
  FilePtr=fopen(buffer,"w");

  fprintf(FilePtr,"# phonon dispersion curve\n");
  fprintf(FilePtr,"# =======================\n");
  fprintf(FilePtr,"# %d points\n",NumberOfDispersionCurvePoints);
  fprintf(FilePtr,"# wave vector: %g %g %g\n",WaveVector.x,WaveVector.y,WaveVector.z);
  fprintf(FilePtr,"# column 1: |k|\n");
  fprintf(FilePtr,"# column 2: f(k) [cm^-1]\n");
  fprintf(FilePtr,"# column 3: f(k) [THz]\n\n");
  for(i=0;i<=NumberOfDispersionCurvePoints;i++)
  {
    k.x=4.0*M_PI*i*WaveVector.x/NumberOfDispersionCurvePoints;
    k.y=4.0*M_PI*i*WaveVector.y/NumberOfDispersionCurvePoints;
    k.z=4.0*M_PI*i*WaveVector.z/NumberOfDispersionCurvePoints;
    k=ConvertFromXYZtoABC(k);
    length=sqrt(SQR(k.x)+SQR(k.y)+SQR(k.z));

    // compute the generalized Hessian
    ComputeDerivativesPhonon(k,NumberOfPositionVariables,NumberOfBoxVariables,Positions,&Energy,Gradients,GeneralizedHessianMatrix,&StrainDerivativeTensor);

    //ProjectConstraintsFromHessianMatrix(NumberOfPositionVariables,NumberOfBoxVariables,Gradients,GeneralizedHessianMatrix);

    RemoveShellInteractionsFromPhonons(GeneralizedHessianMatrix,&ReducedGeneralizedHessianMatrix);
    MassWeightPhononMatrix(CoreSize,ReducedGeneralizedHessianMatrix,Weights);

    SolveEigenValuesAndVectorsPhonon(ReducedGeneralizedHessianMatrix,Frequencies);

    for(j=0;j<CoreSize;j++)
      fprintf(FilePtr,"%g %g %g\n",length/max_length,SIGN(sqrt(fabs(Frequencies[j]))*TO_WAVENUMBERS,Frequencies[j]),SIGN(sqrt(fabs(Frequencies[j]))*TO_THZ,Frequencies[j]));
    fprintf(FilePtr,"\n");
  }
  fclose(FilePtr);

  free(Frequencies); 
  free(Charges); 
  free(Positions); 
  free(Weights); 
  free(Gradients); 

  DeleteComplexMatrix(ReducedGeneralizedHessianMatrix);
  DeleteComplexMatrix(GeneralizedHessianMatrix);

  return 0;
}
