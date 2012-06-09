/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'Inter_force.c' is part of RASPA.

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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "constants.h"
#include "utils.h"
#include "simulation.h"
#include "potentials.h"
#include "output.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "grids.h"
#include "ewald.h"
#include "inter_energy.h"
#include "inter_force.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "spectra.h"
#include "minimization.h"


void CalculateTotalInterVDWForce(void)
{
  int i,j,k,l;
  int typeA,typeB,TypeMolA,TypeMolB;
  REAL rr;
  REAL energy,force_factor;
  VECTOR pos,posA,posB,dr,f;
  REAL ReductionA,ReductionB;
  VECTOR drA,drB,fa,fb;
  VECTOR posA1,posA2,posB1,posB2;
  REAL ra,rb;

  int ConnectedAtomA1,ConnectedAtomA2;
  int ConnectedAtomB1,ConnectedAtomB2;
  VECTOR v,w;
  REAL length_v,length_w;
  REAL dot_product;
  

  UAdsorbateAdsorbateVDW[CurrentSystem]=0.0;
  UAdsorbateCationVDW[CurrentSystem]=0.0;
  UCationCationVDW[CurrentSystem]=0.0;

  if(OmitInterMolecularInteractions) return;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeMolA=Adsorbates[CurrentSystem][i].Type;
    for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[k].AnisotropicPosition;

      ReductionA=0.0;
      ConnectedAtomA1=ConnectedAtomA2=-1;
      ra=drA.x=drA.y=drA.z=0.0;
      v.x=v.y=v.z=length_v=0.0;
      if(PseudoAtoms[typeA].AnisotropicCorrection)
      {
        switch(Components[TypeMolA].Connectivity[k])
        {
          case 0:
            break;
          case 1:
            ConnectedAtomA1=Components[TypeMolA].ConnectivityList[k][0];
            if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
              ReductionA=1.0+PseudoAtoms[typeA].AnisotropicDisplacement;
            else
            {
              posA1=Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Position;
              pos=Adsorbates[CurrentSystem][i].Atoms[k].Position;
              drA.x=posA1.x-pos.x;
              drA.y=posA1.y-pos.y;
              drA.z=posA1.z-pos.z;
              drA=ApplyBoundaryCondition(drA);
              ra=sqrt(SQR(drA.x)+SQR(drA.y)+SQR(drA.z));
              ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
            }
            break;
          case 2:
            switch(Components[TypeMolA].AnisotropicType)
            {
              case ANISOTROPIC_BISECTION:
                printf("ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateTotalInterVDWForce' (inter_force.c)\n");
                exit(0);
                break;
              case ANISOTROPIC_MID_POINT:
                ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
                ConnectedAtomA1=Components[TypeMolA].ConnectivityList[k][0];
                ConnectedAtomA2=Components[TypeMolA].ConnectivityList[k][1];
                pos=Adsorbates[CurrentSystem][i].Atoms[k].Position;
                posA1=Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Position;
                posA2=Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Position;
                v.x=pos.x-0.5*(posA1.x+posA2.x);
                v.y=pos.y-0.5*(posA1.y+posA2.y);
                v.z=pos.z-0.5*(posA1.z+posA2.z);
                length_v=sqrt(SQR(v.x)+SQR(v.y)+SQR(v.z));
                break;
              default:
                printf("ERROR!\n");
                exit(0);
                break; 
            }
            break;
          default:
            break;
        }
      }

      // loop over adsorbant molecules
      if(!OmitAdsorbateAdsorbateVDWInteractions)
      {
        for(j=i+1;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
        {
          TypeMolB=Adsorbates[CurrentSystem][j].Type;
          for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
          {
            typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
            posB=Adsorbates[CurrentSystem][j].Atoms[l].AnisotropicPosition;

            ReductionB=0.0;
            ConnectedAtomB1=ConnectedAtomB2=-1;
            rb=drB.x=drB.y=drB.z=0.0;
            w.x=w.y=w.z=length_w=0.0;
            if(PseudoAtoms[typeB].AnisotropicCorrection)
            {
              switch(Components[TypeMolB].Connectivity[l])
              {
                case 0:
                  break;
                case 1:
                  ConnectedAtomB1=Components[TypeMolB].ConnectivityList[l][0];
                  if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                    ReductionB=1.0+PseudoAtoms[typeB].AnisotropicDisplacement;
                  else
                  {
                    posB1=Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Position;
                    pos=Adsorbates[CurrentSystem][j].Atoms[l].Position;
                    drB.x=posB1.x-pos.x;
                    drB.y=posB1.y-pos.y;
                    drB.z=posB1.z-pos.z;
                    drB=ApplyBoundaryCondition(drB);
                    rb=sqrt(SQR(drB.x)+SQR(drB.y)+SQR(drB.z));
                    ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                  }
                  break;
                case 2:
                  switch(Components[TypeMolB].AnisotropicType)
                  {
                    case ANISOTROPIC_BISECTION:
                      printf("ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateTotalInterVDWForce' (inter_force.c)\n");
                      exit(0);
                      break;
                   case ANISOTROPIC_MID_POINT:
                      ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                      ConnectedAtomB1=Components[TypeMolB].ConnectivityList[l][0];
                      ConnectedAtomB2=Components[TypeMolB].ConnectivityList[l][1];
                      posB1=Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Position;
                      posB2=Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB2].Position;
                      pos=Adsorbates[CurrentSystem][j].Atoms[l].Position;
                      w.x=pos.x-0.5*(posB1.x+posB2.x);
                      w.y=pos.y-0.5*(posB1.y+posB2.y);
                      w.z=pos.z-0.5*(posB1.z+posB2.z);
                      length_w=sqrt(SQR(w.x)+SQR(w.y)+SQR(w.z));
                      break;
                  }
                  break;
                default:
                  break;
              }
            }

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              PotentialGradient(typeA,typeB,rr,&energy,&force_factor);

              UAdsorbateAdsorbateVDW[CurrentSystem]+=energy;

              // forces
              f.x=force_factor*dr.x;
              f.y=force_factor*dr.y;
              f.z=force_factor*dr.z;

              if(PseudoAtoms[typeA].AnisotropicCorrection)
              {
                switch(Components[TypeMolA].Connectivity[k])
                {
                  case 0:
                    break;
                  case 1:
                    if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
                    {
                      Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=ReductionA*f.x;
                      Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=ReductionA*f.y;
                      Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=ReductionA*f.z;

                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=(1.0-ReductionA)*f.x;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=(1.0-ReductionA)*f.y;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=(1.0-ReductionA)*f.z;
                    }
                    else
                    {
                      dot_product=dr.x*drA.x+dr.y*drA.y+dr.z*drA.z;

                      fa.x=(-ReductionA*drA.x*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.x)*force_factor;
                      fa.y=(-ReductionA*drA.y*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.y)*force_factor;
                      fa.z=(-ReductionA*drA.z*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.z)*force_factor;

                      fb.x=(ReductionA*drA.x*dot_product/CUBE(ra)-(ReductionA/ra)*dr.x)*force_factor;
                      fb.y=(ReductionA*drA.y*dot_product/CUBE(ra)-(ReductionA/ra)*dr.y)*force_factor;
                      fb.z=(ReductionA*drA.z*dot_product/CUBE(ra)-(ReductionA/ra)*dr.z)*force_factor;

                      Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=fa.x;
                      Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=fa.y;
                      Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=fa.z;

                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=fb.x;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=fb.y;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=fb.z;
                    }
                    break;
                  case 2:
                    dot_product=v.x*dr.x+v.y*dr.y+v.z*dr.z;

                    fa.x=0.5*(ReductionA*v.x*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.x)*force_factor;
                    fa.y=0.5*(ReductionA*v.y*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.y)*force_factor;
                    fa.z=0.5*(ReductionA*v.z*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.z)*force_factor;

                    fb.x=(-ReductionA*v.x*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.x)*force_factor;
                    fb.y=(-ReductionA*v.y*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.y)*force_factor;
                    fb.z=(-ReductionA*v.z*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.z)*force_factor;

                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=fa.x;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=fa.y;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=fa.z;

                    Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=fb.x;
                    Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=fb.y;
                    Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=fb.z;

                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.x-=fa.x;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.y-=fa.y;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.z-=fa.z;
                    break;
                  default:
                    printf("Not yet implemented in routine 'CalculateTotalInterVDWForce'\n");
                    exit(0);
                    break;
                }
              }
              else
              {
                Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=f.x;
                Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=f.y;
                Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=f.z;
              }

              if(PseudoAtoms[typeB].AnisotropicCorrection)
              {
                switch(Components[TypeMolB].Connectivity[l])
                {
                  case 0:
                    break;
                  case 1:
                    if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                    {
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.x+=ReductionB*f.x;
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.y+=ReductionB*f.y;
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.z+=ReductionB*f.z;

                      Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.x+=(1.0-ReductionB)*f.x;
                      Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.y+=(1.0-ReductionB)*f.y;
                      Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.z+=(1.0-ReductionB)*f.z;
                    }
                    else
                    {
                      dot_product=dr.x*drB.x+dr.y*drB.y+dr.z*drB.z;

                      fa.x=(-ReductionB*drB.x*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.x)*force_factor;
                      fa.y=(-ReductionB*drB.y*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.y)*force_factor;
                      fa.z=(-ReductionB*drB.z*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.z)*force_factor;

                      fb.x=(ReductionB*drB.x*dot_product/CUBE(rb)-(ReductionB/rb)*dr.x)*force_factor;
                      fb.y=(ReductionB*drB.y*dot_product/CUBE(rb)-(ReductionB/rb)*dr.y)*force_factor;
                      fb.z=(ReductionB*drB.z*dot_product/CUBE(rb)-(ReductionB/rb)*dr.z)*force_factor;

                      Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.x+=fb.x;
                      Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.y+=fb.y;
                      Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.z+=fb.z;

                      Adsorbates[CurrentSystem][j].Atoms[l].Force.x+=fa.x;
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.y+=fa.y;
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.z+=fa.z;
                    }
                    break;
                  case 2:
                    dot_product=w.x*dr.x+w.y*dr.y+w.z*dr.z;

                    fa.x=0.5*(ReductionB*w.x*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.x)*force_factor;
                    fa.y=0.5*(ReductionB*w.y*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.y)*force_factor;
                    fa.z=0.5*(ReductionB*w.z*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.z)*force_factor;

                    fb.x=(-ReductionB*w.x*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.x)*force_factor;
                    fb.y=(-ReductionB*w.y*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.y)*force_factor;
                    fb.z=(-ReductionB*w.z*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.z)*force_factor;

                    Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.x+=fa.x;
                    Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.y+=fa.y;
                    Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.z+=fa.z;

                    Adsorbates[CurrentSystem][j].Atoms[l].Force.x+=fb.x;
                    Adsorbates[CurrentSystem][j].Atoms[l].Force.y+=fb.y;
                    Adsorbates[CurrentSystem][j].Atoms[l].Force.z+=fb.z;

                    Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB2].Force.x+=fa.x;
                    Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB2].Force.y+=fa.y;
                    Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB2].Force.z+=fa.z;
                    break;
                  default:
                    printf("Not yet implemented in routine 'CalculateTotalInterVDWForce'\n");
                    exit(0);
                    break;
                }
              }
              else
              {
                Adsorbates[CurrentSystem][j].Atoms[l].Force.x+=f.x;
                Adsorbates[CurrentSystem][j].Atoms[l].Force.y+=f.y;
                Adsorbates[CurrentSystem][j].Atoms[l].Force.z+=f.z;
              }

              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
            }
          }
        }
      }

      if(!OmitAdsorbateCationVDWInteractions)
      {
        for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
        {
          TypeMolB=Cations[CurrentSystem][j].Type;
          for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
          {
            posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
            typeB=Cations[CurrentSystem][j].Atoms[l].Type;

            ReductionB=0.0;
            ConnectedAtomB1=ConnectedAtomB2=-1;
            rb=drB.x=drB.y=drB.z=0.0;
            w.x=w.y=w.z=length_w=0.0;
            if(PseudoAtoms[typeB].AnisotropicCorrection)
            {
              switch(Components[TypeMolB].Connectivity[l])
              {
                case 0:
                  break;
                case 1:
                  ConnectedAtomB1=Components[TypeMolB].ConnectivityList[l][0];
                  if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                    ReductionB=1.0+PseudoAtoms[typeB].AnisotropicDisplacement;
                  else
                  {
                    posB1=Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Position;
                    pos=Cations[CurrentSystem][j].Atoms[l].Position;
                    drB.x=posB1.x-pos.x;
                    drB.y=posB1.y-pos.y;
                    drB.z=posB1.z-pos.z;
                    drB=ApplyBoundaryCondition(drB);
                    rb=sqrt(SQR(drB.x)+SQR(drB.y)+SQR(drB.z));
                    ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                  }
                  break;
                case 2:
                  switch(Components[TypeMolB].AnisotropicType)
                  {
                    case ANISOTROPIC_BISECTION:
                      printf("ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateTotalInterVDWForce' (inter_force.c)\n");
                      exit(0);
                      break;
                   case ANISOTROPIC_MID_POINT:
                      ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                      ConnectedAtomB1=Components[TypeMolB].ConnectivityList[l][0];
                      ConnectedAtomB2=Components[TypeMolB].ConnectivityList[l][1];
                      posB1=Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Position;
                      posB2=Cations[CurrentSystem][j].Atoms[ConnectedAtomB2].Position;
                      pos=Cations[CurrentSystem][j].Atoms[l].Position;
                      w.x=pos.x-0.5*(posB1.x+posB2.x);
                      w.y=pos.y-0.5*(posB1.y+posB2.y);
                      w.z=pos.z-0.5*(posB1.z+posB2.z);
                      length_w=sqrt(SQR(w.x)+SQR(w.y)+SQR(w.z));
                      break;
                  }
                  break;
                default:
                  break;
              }
            }

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              PotentialGradient(typeA,typeB,rr,&energy,&force_factor);

              // energy
              UAdsorbateCationVDW[CurrentSystem]+=energy;

              // forces
              f.x=force_factor*dr.x;
              f.y=force_factor*dr.y;
              f.z=force_factor*dr.z;

              if(PseudoAtoms[typeA].AnisotropicCorrection)
              {
                switch(Components[TypeMolA].Connectivity[k])
                {
                  case 0:
                    break;
                  case 1:
                    if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
                    {
                      Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=ReductionA*f.x;
                      Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=ReductionA*f.y;
                      Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=ReductionA*f.z;

                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=(1.0-ReductionA)*f.x;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=(1.0-ReductionA)*f.y;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=(1.0-ReductionA)*f.z;
                    }
                    else
                    {
                      dot_product=dr.x*drA.x+dr.y*drA.y+dr.z*drA.z;

                      fa.x=(-ReductionA*drA.x*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.x)*force_factor;
                      fa.y=(-ReductionA*drA.y*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.y)*force_factor;
                      fa.z=(-ReductionA*drA.z*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.z)*force_factor;

                      fb.x=(ReductionA*drA.x*dot_product/CUBE(ra)-(ReductionA/ra)*dr.x)*force_factor;
                      fb.y=(ReductionA*drA.y*dot_product/CUBE(ra)-(ReductionA/ra)*dr.y)*force_factor;
                      fb.z=(ReductionA*drA.z*dot_product/CUBE(ra)-(ReductionA/ra)*dr.z)*force_factor;

                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=fb.x;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=fb.y;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=fb.z;

                      Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=fa.x;
                      Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=fa.y;
                      Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=fa.z;
                    }
                    break;
                  case 2:
                    dot_product=v.x*dr.x+v.y*dr.y+v.z*dr.z;

                    fa.x=0.5*(ReductionA*v.x*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.x)*force_factor;
                    fa.y=0.5*(ReductionA*v.y*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.y)*force_factor;
                    fa.z=0.5*(ReductionA*v.z*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.z)*force_factor;

                    fb.x=(-ReductionA*v.x*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.x)*force_factor;
                    fb.y=(-ReductionA*v.y*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.y)*force_factor;
                    fb.z=(-ReductionA*v.z*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.z)*force_factor;

                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=fa.x;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=fa.y;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=fa.z;

                    Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=fb.x;
                    Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=fb.y;
                    Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=fb.z;

                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.x-=fa.x;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.y-=fa.y;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.z-=fa.z;
                    break;
                  default:
                    printf("Not yet implemented in routine 'CalculateTotalInterVDWForce'\n");
                    exit(0);
                    break;
                }
              }
              else
              {
                Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=f.x;
                Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=f.y;
                Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=f.z;
              }

              if(PseudoAtoms[typeB].AnisotropicCorrection)
              {
                switch(Components[TypeMolB].Connectivity[l])
                {
                  case 0:
                    break;
                  case 1:
                    if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                    {
                      Cations[CurrentSystem][j].Atoms[l].Force.x+=ReductionB*f.x;
                      Cations[CurrentSystem][j].Atoms[l].Force.y+=ReductionB*f.y;
                      Cations[CurrentSystem][j].Atoms[l].Force.z+=ReductionB*f.z;

                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.x+=(1.0-ReductionB)*f.x;
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.y+=(1.0-ReductionB)*f.y;
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.z+=(1.0-ReductionB)*f.z;
                    }
                    else
                    {
                      dot_product=dr.x*drB.x+dr.y*drB.y+dr.z*drB.z;

                      fa.x=(-ReductionB*drB.x*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.x)*force_factor;
                      fa.y=(-ReductionB*drB.y*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.y)*force_factor;
                      fa.z=(-ReductionB*drB.z*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.z)*force_factor;

                      fb.x=(ReductionB*drB.x*dot_product/CUBE(rb)-(ReductionB/rb)*dr.x)*force_factor;
                      fb.y=(ReductionB*drB.y*dot_product/CUBE(rb)-(ReductionB/rb)*dr.y)*force_factor;
                      fb.z=(ReductionB*drB.z*dot_product/CUBE(rb)-(ReductionB/rb)*dr.z)*force_factor;

                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.x+=fb.x;
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.y+=fb.y;
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.z+=fb.z;

                      Cations[CurrentSystem][j].Atoms[l].Force.x+=fa.x;
                      Cations[CurrentSystem][j].Atoms[l].Force.y+=fa.y;
                      Cations[CurrentSystem][j].Atoms[l].Force.z+=fa.z;
                    }
                    break;
                  case 2:
                    dot_product=w.x*dr.x+w.y*dr.y+w.z*dr.z;

                    fa.x=0.5*(ReductionB*w.x*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.x)*force_factor;
                    fa.y=0.5*(ReductionB*w.y*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.y)*force_factor;
                    fa.z=0.5*(ReductionB*w.z*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.z)*force_factor;

                    fb.x=(-ReductionB*w.x*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.x)*force_factor;
                    fb.y=(-ReductionB*w.y*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.y)*force_factor;
                    fb.z=(-ReductionB*w.z*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.z)*force_factor;

                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.x+=fa.x;
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.y+=fa.y;
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.z+=fa.z;

                    Cations[CurrentSystem][j].Atoms[l].Force.x+=fb.x;
                    Cations[CurrentSystem][j].Atoms[l].Force.y+=fb.y;
                    Cations[CurrentSystem][j].Atoms[l].Force.z+=fb.z;

                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB2].Force.x+=fa.x;
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB2].Force.y+=fa.y;
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB2].Force.z+=fa.z;
                    break;
                  default:
                    printf("Not yet implemented in routine 'CalculateTotalInterVDWForce'\n");
                    exit(0);
                    break;
                }
              }
              else
              {
                Cations[CurrentSystem][j].Atoms[l].Force.x+=f.x;
                Cations[CurrentSystem][j].Atoms[l].Force.y+=f.y;
                Cations[CurrentSystem][j].Atoms[l].Force.z+=f.z;
              }

              // stress-tensor
              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
            }
          }
        }
      }
    }
  }

  if(!OmitCationCationVDWInteractions)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      TypeMolA=Cations[CurrentSystem][i].Type;
      for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Cations[CurrentSystem][i].Atoms[k].Type;
        posA=Cations[CurrentSystem][i].Atoms[k].AnisotropicPosition;

        ReductionA=0.0;
        ConnectedAtomA1=ConnectedAtomA2=-1;
        ra=drA.x=drA.y=drA.z=0.0;
        v.x=v.y=v.z=length_v=0.0;
        if(PseudoAtoms[typeA].AnisotropicCorrection)
        {
          switch(Components[TypeMolA].Connectivity[k])
          {
            case 0:
              break;
            case 1:
              ConnectedAtomA1=Components[TypeMolA].ConnectivityList[k][0];
              if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
                ReductionA=1.0+PseudoAtoms[typeA].AnisotropicDisplacement;
              else
              {
                posA1=Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Position;
                pos=Cations[CurrentSystem][i].Atoms[k].Position;
                drA.x=posA1.x-pos.x;
                drA.y=posA1.y-pos.y;
                drA.z=posA1.z-pos.z;
                drA=ApplyBoundaryCondition(drA);
                ra=sqrt(SQR(drA.x)+SQR(drA.y)+SQR(drA.z));
                ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
              }
              break;
            case 2:
              switch(Components[TypeMolA].AnisotropicType)
              {
                case ANISOTROPIC_BISECTION:
                  printf("ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateTotalInterVDWForce' (inter_force.c)\n");
                  exit(0);
                  break;
                case ANISOTROPIC_MID_POINT:
                  ReductionA=PseudoAtoms[typeA].AnisotropicDisplacement;
                  ConnectedAtomA1=Components[TypeMolA].ConnectivityList[k][0];
                  ConnectedAtomA2=Components[TypeMolA].ConnectivityList[k][1];
                  pos=Cations[CurrentSystem][i].Atoms[k].Position;
                  posA1=Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Position;
                  posA2=Cations[CurrentSystem][i].Atoms[ConnectedAtomA2].Position;
                  v.x=pos.x-0.5*(posA1.x+posA2.x);
                  v.y=pos.y-0.5*(posA1.y+posA2.y);
                  v.z=pos.z-0.5*(posA1.z+posA2.z);
                  length_v=sqrt(SQR(v.x)+SQR(v.y)+SQR(v.z));
                  break;
                default:
                  printf("ERROR!\n");
                  exit(0);
                  break;
              }
              break;
            default:
              break;
          }
        }

        // loop over cation molecules
        for(j=i+1;j<NumberOfCationMolecules[CurrentSystem];j++)
        {
          TypeMolB=Cations[CurrentSystem][j].Type;
          for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
          {
            posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
            typeB=Cations[CurrentSystem][j].Atoms[l].Type;

            ReductionB=0.0;
            ConnectedAtomB1=ConnectedAtomB2=-1;
            rb=drB.x=drB.y=drB.z=0.0;
            w.x=w.y=w.z=length_w=0.0;
            if(PseudoAtoms[typeB].AnisotropicCorrection)
            {
              switch(Components[TypeMolB].Connectivity[l])
              {
                case 0:
                  break;
                case 1:
                  ConnectedAtomB1=Components[TypeMolB].ConnectivityList[l][0];
                  if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                    ReductionB=1.0+PseudoAtoms[typeB].AnisotropicDisplacement;
                  else
                  {
                    posB1=Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Position;
                    pos=Cations[CurrentSystem][j].Atoms[l].Position;
                    drB.x=posB1.x-pos.x;
                    drB.y=posB1.y-pos.y;
                    drB.z=posB1.z-pos.z;
                    drB=ApplyBoundaryCondition(drB);
                    rb=sqrt(SQR(drB.x)+SQR(drB.y)+SQR(drB.z));
                    ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                  }
                  break;
                case 2:
                  switch(Components[TypeMolB].AnisotropicType)
                  {
                    case ANISOTROPIC_BISECTION:
                      printf("ERROR: 'ANISOTROPIC_BISECTION' not implemented yet in routine 'CalculateTotalInterVDWForce' (inter_force.c)\n");
                      exit(0);
                      break;
                   case ANISOTROPIC_MID_POINT:
                      ReductionB=PseudoAtoms[typeB].AnisotropicDisplacement;
                      ConnectedAtomB1=Components[TypeMolB].ConnectivityList[l][0];
                      ConnectedAtomB2=Components[TypeMolB].ConnectivityList[l][1];
                      posB1=Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Position;
                      posB2=Cations[CurrentSystem][j].Atoms[ConnectedAtomB2].Position;
                      pos=Cations[CurrentSystem][j].Atoms[l].Position;
                      w.x=pos.x-0.5*(posB1.x+posB2.x);
                      w.y=pos.y-0.5*(posB1.y+posB2.y);
                      w.z=pos.z-0.5*(posB1.z+posB2.z);
                      length_w=sqrt(SQR(w.x)+SQR(w.y)+SQR(w.z));
                      break;
                  }
                  break;
                default:
                  break;
              }
            }

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffVDWSquared)
            {
              PotentialGradient(typeA,typeB,rr,&energy,&force_factor);

              // energy
              UCationCationVDW[CurrentSystem]+=energy;

              // forces
              f.x=force_factor*dr.x;
              f.y=force_factor*dr.y;
              f.z=force_factor*dr.z;

              if(PseudoAtoms[typeA].AnisotropicCorrection)
              {
                switch(Components[TypeMolA].Connectivity[k])
                {
                  case 0:
                    break;
                  case 1:
                    if(PseudoAtoms[typeA].AnisotropicType==RELATIVE)
                    {
                      Cations[CurrentSystem][i].Atoms[k].Force.x-=ReductionA*f.x;
                      Cations[CurrentSystem][i].Atoms[k].Force.y-=ReductionA*f.y;
                      Cations[CurrentSystem][i].Atoms[k].Force.z-=ReductionA*f.z;

                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=(1.0-ReductionA)*f.x;
                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=(1.0-ReductionA)*f.y;
                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=(1.0-ReductionA)*f.z;
                    }
                    else
                    {
                      dot_product=dr.x*drA.x+dr.y*drA.y+dr.z*drA.z;

                      fa.x=(-ReductionA*drA.x*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.x)*force_factor;
                      fa.y=(-ReductionA*drA.y*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.y)*force_factor;
                      fa.z=(-ReductionA*drA.z*dot_product/CUBE(ra)+(1.0+ReductionA/ra)*dr.z)*force_factor;

                      fb.x=(ReductionA*drA.x*dot_product/CUBE(ra)-(ReductionA/ra)*dr.x)*force_factor;
                      fb.y=(ReductionA*drA.y*dot_product/CUBE(ra)-(ReductionA/ra)*dr.y)*force_factor;
                      fb.z=(ReductionA*drA.z*dot_product/CUBE(ra)-(ReductionA/ra)*dr.z)*force_factor;

                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=fb.x;
                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=fb.y;
                      Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=fb.z;

                      Cations[CurrentSystem][i].Atoms[k].Force.x-=fa.x;
                      Cations[CurrentSystem][i].Atoms[k].Force.y-=fa.y;
                      Cations[CurrentSystem][i].Atoms[k].Force.z-=fa.z;
                    }
                    break;
                  case 2:
                    dot_product=v.x*dr.x+v.y*dr.y+v.z*dr.z;

                    fa.x=0.5*(ReductionA*v.x*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.x)*force_factor;
                    fa.y=0.5*(ReductionA*v.y*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.y)*force_factor;
                    fa.z=0.5*(ReductionA*v.z*dot_product/CUBE(length_v)-(ReductionA/length_v)*dr.z)*force_factor;

                    fb.x=(-ReductionA*v.x*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.x)*force_factor;
                    fb.y=(-ReductionA*v.y*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.y)*force_factor;
                    fb.z=(-ReductionA*v.z*dot_product/CUBE(length_v)+(1.0+ReductionA/length_v)*dr.z)*force_factor;

                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.x-=fa.x;
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.y-=fa.y;
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA1].Force.z-=fa.z;

                    Cations[CurrentSystem][i].Atoms[k].Force.x-=fb.x;
                    Cations[CurrentSystem][i].Atoms[k].Force.y-=fb.y;
                    Cations[CurrentSystem][i].Atoms[k].Force.z-=fb.z;

                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.x-=fa.x;
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.y-=fa.y;
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA2].Force.z-=fa.z;
                    break;
                  default:
                    printf("Not yet implemented in routine 'CalculateTotalInterVDWForce'\n");
                    exit(0);
                    break;
                }
              }
              else
              {
                Cations[CurrentSystem][i].Atoms[k].Force.x-=f.x;
                Cations[CurrentSystem][i].Atoms[k].Force.y-=f.y;
                Cations[CurrentSystem][i].Atoms[k].Force.z-=f.z;
              }

              if(PseudoAtoms[typeB].AnisotropicCorrection)
              {
                switch(Components[TypeMolB].Connectivity[l])
                {
                  case 0:
                    break;
                  case 1:
                    if(PseudoAtoms[typeB].AnisotropicType==RELATIVE)
                    {
                      Cations[CurrentSystem][j].Atoms[l].Force.x+=ReductionB*f.x;
                      Cations[CurrentSystem][j].Atoms[l].Force.y+=ReductionB*f.y;
                      Cations[CurrentSystem][j].Atoms[l].Force.z+=ReductionB*f.z;

                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.x+=(1.0-ReductionB)*f.x;
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.y+=(1.0-ReductionB)*f.y;
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.z+=(1.0-ReductionB)*f.z;
                    }
                    else
                    {
                      dot_product=dr.x*drB.x+dr.y*drB.y+dr.z*drB.z;

                      fa.x=(-ReductionB*drB.x*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.x)*force_factor;
                      fa.y=(-ReductionB*drB.y*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.y)*force_factor;
                      fa.z=(-ReductionB*drB.z*dot_product/CUBE(rb)+(1.0+ReductionB/rb)*dr.z)*force_factor;

                      fb.x=(ReductionB*drB.x*dot_product/CUBE(rb)-(ReductionB/rb)*dr.x)*force_factor;
                      fb.y=(ReductionB*drB.y*dot_product/CUBE(rb)-(ReductionB/rb)*dr.y)*force_factor;
                      fb.z=(ReductionB*drB.z*dot_product/CUBE(rb)-(ReductionB/rb)*dr.z)*force_factor;

                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.x+=fb.x;
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.y+=fb.y;
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.z+=fb.z;

                      Cations[CurrentSystem][j].Atoms[l].Force.x+=fa.x;
                      Cations[CurrentSystem][j].Atoms[l].Force.y+=fa.y;
                      Cations[CurrentSystem][j].Atoms[l].Force.z+=fa.z;
                    }
                    break;
                  case 2:
                    dot_product=w.x*dr.x+w.y*dr.y+w.z*dr.z;

                    fa.x=0.5*(ReductionB*w.x*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.x)*force_factor;
                    fa.y=0.5*(ReductionB*w.y*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.y)*force_factor;
                    fa.z=0.5*(ReductionB*w.z*dot_product/CUBE(length_w)-(ReductionB/length_w)*dr.z)*force_factor;

                    fb.x=(-ReductionB*w.x*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.x)*force_factor;
                    fb.y=(-ReductionB*w.y*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.y)*force_factor;
                    fb.z=(-ReductionB*w.z*dot_product/CUBE(length_w)+(1.0+ReductionB/length_w)*dr.z)*force_factor;

                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.x+=fa.x;
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.y+=fa.y;
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB1].Force.z+=fa.z;

                    Cations[CurrentSystem][j].Atoms[l].Force.x+=fb.x;
                    Cations[CurrentSystem][j].Atoms[l].Force.y+=fb.y;
                    Cations[CurrentSystem][j].Atoms[l].Force.z+=fb.z;

                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB2].Force.x+=fa.x;
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB2].Force.y+=fa.y;
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB2].Force.z+=fa.z;
                    break;
                  default:
                    printf("Not yet implemented in routine 'CalculateTotalInterVDWForce'\n");
                    exit(0);
                    break;
                }
              }
              else
              {
                Cations[CurrentSystem][j].Atoms[l].Force.x+=f.x;
                Cations[CurrentSystem][j].Atoms[l].Force.y+=f.y;
                Cations[CurrentSystem][j].Atoms[l].Force.z+=f.z;
              }

              // stress-tensor
              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
            }
          }
        }
      }
    }
  }
}

int CalculateTotalInterChargeChargeCoulombForce(void)
{
  int i,j,k,l;
  int typeA,typeB;
  REAL r,rr;
  REAL chargeA,chargeB;
  REAL U,force_factor;
  VECTOR posA,posB,dr,f;
  REAL SwitchingValueDerivative,SwitchingValue;
  REAL TranslationValue,TranslationValueDerivative;
  REAL UWolfCorrection,NetChargeB;

  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeReal[CurrentSystem]=0.0;
  UCationCationChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[k].Position;
      chargeA=Adsorbates[CurrentSystem][i].Atoms[k].Charge;

      // loop over adsorbant molecules
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(j=i+1;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
        {
          for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
          {
            posB=Adsorbates[CurrentSystem][j].Atoms[l].Position;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
              chargeB=Adsorbates[CurrentSystem][j].Atoms[l].Charge;

              switch(ChargeMethod)
              {
                case NONE:
                  force_factor=0.0;
                  break;
                case TRUNCATED_COULOMB:
                  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                  break;
                case SHIFTED_COULOMB:
                  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                  break;
                case SMOOTHED_COULOMB:
                  U=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/rr;
                  if(rr>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                   SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
                    SwitchingValueDerivative=5.0*SwitchingChargeChargeFactors5[5]*rr*rr+4.0*SwitchingChargeChargeFactors5[4]*rr*r+3.0*SwitchingChargeChargeFactors5[3]*rr+
                                             2.0*SwitchingChargeChargeFactors5[2]*r+SwitchingChargeChargeFactors5[1];

                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                   (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                    SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                    SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                    TranslationValueDerivative=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                             (7.0*SwitchingChargeChargeFactors7[7]*rr*rr*rr+6.0*SwitchingChargeChargeFactors7[6]*rr*rr*r+
                                              5.0*SwitchingChargeChargeFactors7[5]*rr*rr+4.0*SwitchingChargeChargeFactors7[4]*rr*r+3.0*SwitchingChargeChargeFactors7[3]*rr+
                                              2.0*SwitchingChargeChargeFactors7[2]*r+SwitchingChargeChargeFactors7[1]);
                    force_factor=U*SwitchingValueDerivative+force_factor*SwitchingValue+TranslationValueDerivative;
                    U=U*SwitchingValue+TranslationValue;
                  }
                  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=U;
                  force_factor/=r;
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                            (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                            (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/rr+
                                 2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                  2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                  break;
                case EWALD:
                  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                          (erfc(Alpha[CurrentSystem]*r)/r);

                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                        (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                        (r*rr);
                  break;
                default:
                  printf("Unknown charge-method in 'CalculateTotalInterChargeChargeCoulombForce'\n");
                  exit(0);
                  break;
              }

              // forces
              f.x=force_factor*dr.x;
              f.y=force_factor*dr.y;
              f.z=force_factor*dr.z;

              Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=f.x;
              Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=f.y;
              Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=f.z;

              Adsorbates[CurrentSystem][j].Atoms[l].Force.x+=f.x;
              Adsorbates[CurrentSystem][j].Atoms[l].Force.y+=f.y;
              Adsorbates[CurrentSystem][j].Atoms[l].Force.z+=f.z;

              // the strain derivative
              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
            }
          }
        }
      }

      if(!OmitAdsorbateCationCoulombInteractions)
      {
        for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
        {
          for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
          {
            posB=Cations[CurrentSystem][j].Atoms[l].Position;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffChargeChargeSquared)
            {
              typeB=Cations[CurrentSystem][j].Atoms[l].Type;
              chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;

              r=sqrt(rr);

              switch(ChargeMethod)
              {
                case NONE:
                  force_factor=0.0;
                  break;
                case TRUNCATED_COULOMB:
                  UAdsorbateCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                  break;
                case SHIFTED_COULOMB:
                  UAdsorbateCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                  break;
                case SMOOTHED_COULOMB:
                  U=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/rr;
                  if(rr>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                   SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
                    SwitchingValueDerivative=5.0*SwitchingChargeChargeFactors5[5]*rr*rr+4.0*SwitchingChargeChargeFactors5[4]*rr*r+3.0*SwitchingChargeChargeFactors5[3]*rr+
                                             2.0*SwitchingChargeChargeFactors5[2]*r+SwitchingChargeChargeFactors5[1];

                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                   (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                    SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                    SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                    TranslationValueDerivative=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                             (7.0*SwitchingChargeChargeFactors7[7]*rr*rr*rr+6.0*SwitchingChargeChargeFactors7[6]*rr*rr*r+
                                              5.0*SwitchingChargeChargeFactors7[5]*rr*rr+4.0*SwitchingChargeChargeFactors7[4]*rr*r+3.0*SwitchingChargeChargeFactors7[3]*rr+
                                              2.0*SwitchingChargeChargeFactors7[2]*r+SwitchingChargeChargeFactors7[1]);
                    force_factor=U*SwitchingValueDerivative+force_factor*SwitchingValue+TranslationValueDerivative;
                    U=U*SwitchingValue+TranslationValue;
                  }
                  UAdsorbateCationChargeChargeReal[CurrentSystem]+=U;
                  force_factor/=r;
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  UAdsorbateCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                            (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                            (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/rr+
                                 2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                  2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                  break;
                case EWALD:
                  UAdsorbateCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (erfc(Alpha[CurrentSystem]*r)/r);

                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                         (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                         (r*rr);
                 break;
               default:
                 printf("Unknown charge-method in 'CalculateTotalInterChargeChargeCoulombForce'\n");
                 exit(0);
                 break;
              }

              // forces
              f.x=force_factor*dr.x;
              f.y=force_factor*dr.y;
              f.z=force_factor*dr.z;

              Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=f.x;
              Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=f.y;
              Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=f.z;

              Cations[CurrentSystem][j].Atoms[l].Force.x+=f.x;
              Cations[CurrentSystem][j].Atoms[l].Force.y+=f.y;
              Cations[CurrentSystem][j].Atoms[l].Force.z+=f.z;

              // stress-tensor
              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
            }
          }
        }
      }
    }
  }

  if(!OmitCationCationCoulombInteractions)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
      {
        typeA=Cations[CurrentSystem][i].Atoms[k].Type;
        posA=Cations[CurrentSystem][i].Atoms[k].Position;
        chargeA=Cations[CurrentSystem][i].Atoms[k].Charge;

        // loop over cation molecules
        for(j=i+1;j<NumberOfCationMolecules[CurrentSystem];j++)
        {
          for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
          {
            posB=Cations[CurrentSystem][j].Atoms[l].Position;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<CutOffChargeChargeSquared)
            {
              typeB=Cations[CurrentSystem][j].Atoms[l].Type;
              chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;

              r=sqrt(rr);

              switch(ChargeMethod)
              {
                case NONE:
                  force_factor=0.0;
                  break;
                case TRUNCATED_COULOMB:
                  UCationCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                  break;
                case SHIFTED_COULOMB:
                  UCationCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                  break;
                case SMOOTHED_COULOMB:
                  U=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-2.0/(CutOffChargeChargeSwitch+CutOffChargeCharge));
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/rr;
                  if(rr>CutOffChargeChargeSwitchSquared)
                  {
                    SwitchingValue=SwitchingChargeChargeFactors5[5]*(rr*rr*r)+SwitchingChargeChargeFactors5[4]*(rr*rr)+SwitchingChargeChargeFactors5[3]*(rr*r)+
                                   SwitchingChargeChargeFactors5[2]*rr+SwitchingChargeChargeFactors5[1]*r+SwitchingChargeChargeFactors5[0];
                    SwitchingValueDerivative=5.0*SwitchingChargeChargeFactors5[5]*rr*rr+4.0*SwitchingChargeChargeFactors5[4]*rr*r+3.0*SwitchingChargeChargeFactors5[3]*rr+
                                             2.0*SwitchingChargeChargeFactors5[2]*r+SwitchingChargeChargeFactors5[1];

                    TranslationValue=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                   (SwitchingChargeChargeFactors7[7]*(rr*rr*rr*r)+SwitchingChargeChargeFactors7[6]*(rr*rr*rr)+
                                    SwitchingChargeChargeFactors7[5]*(rr*rr*r)+SwitchingChargeChargeFactors7[4]*(rr*rr)+SwitchingChargeChargeFactors7[3]*(rr*r)+
                                    SwitchingChargeChargeFactors7[2]*rr+SwitchingChargeChargeFactors7[1]*r+SwitchingChargeChargeFactors7[0]);
                    TranslationValueDerivative=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                             (7.0*SwitchingChargeChargeFactors7[7]*rr*rr*rr+6.0*SwitchingChargeChargeFactors7[6]*rr*rr*r+
                                              5.0*SwitchingChargeChargeFactors7[5]*rr*rr+4.0*SwitchingChargeChargeFactors7[4]*rr*r+3.0*SwitchingChargeChargeFactors7[3]*rr+
                                              2.0*SwitchingChargeChargeFactors7[2]*r+SwitchingChargeChargeFactors7[1]);
                    force_factor=U*SwitchingValueDerivative+force_factor*SwitchingValue+TranslationValueDerivative;
                    U=U*SwitchingValue+TranslationValue;
                  }
                  UCationCationChargeChargeReal[CurrentSystem]+=U;
                  force_factor/=r;
                  break;
                case WOLFS_METHOD_DAMPED_FG:
                  UCationCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                           -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                            (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                            (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/rr+
                                 2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                  2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                  break;
                case EWALD:
                  UCationCationChargeChargeReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                                       (erfc(Alpha[CurrentSystem]*r)/r);
                  force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                        (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                        (r*rr);
                  break;
                default:
                  printf("Unknown charge-method in 'CalculateTotalInterChargeChargeCoulombForce'\n");
                  exit(0);
                  break;
              }

              // forces
              f.x=force_factor*dr.x;
              f.y=force_factor*dr.y;
              f.z=force_factor*dr.z;

              Cations[CurrentSystem][i].Atoms[k].Force.x-=f.x;
              Cations[CurrentSystem][i].Atoms[k].Force.y-=f.y;
              Cations[CurrentSystem][i].Atoms[k].Force.z-=f.z;

              Cations[CurrentSystem][j].Atoms[l].Force.x+=f.x;
              Cations[CurrentSystem][j].Atoms[l].Force.y+=f.y;
              Cations[CurrentSystem][j].Atoms[l].Force.z+=f.z;

              // stress-tensor
              StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
              StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
              StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

              StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
              StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
              StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

              StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
              StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
              StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
            }
          }
        }
      }
    }
  }

  switch(ChargeMethod)
  {
    case WOLFS_METHOD:
      UWolfCorrection=0.0;
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=0.5*COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*InverseCutOffChargeCharge;
      }
      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UWolfCorrection;

      UWolfCorrection=0.0;
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Cations[CurrentSystem][i].Atoms[j].Type;
          chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=0.5*COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*InverseCutOffChargeCharge;
      }
      UCationCationChargeChargeReal[CurrentSystem]+=UWolfCorrection;
      break;
    case WOLFS_METHOD_DAMPED:
    case WOLFS_METHOD_DAMPED_FG:
      UWolfCorrection=0.0;
      for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
          chargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*
             (0.5*erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
              Alpha[CurrentSystem]*M_1_SQRTPI);
      }
      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UWolfCorrection;

      UWolfCorrection=0.0;
      for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      {
        NetChargeB=0.0;
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          typeA=Cations[CurrentSystem][i].Atoms[j].Type;
          chargeA=Cations[CurrentSystem][i].Atoms[j].Charge;
          NetChargeB+=chargeA;
        }
        UWolfCorrection-=COULOMBIC_CONVERSION_FACTOR*SQR(NetChargeB)*
             (0.5*erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
              Alpha[CurrentSystem]*M_1_SQRTPI);
      }
      UCationCationChargeChargeReal[CurrentSystem]+=UWolfCorrection;
      break;
    default:
      break;
  }

  return 0;
}

int CalculateTotalInterChargeBondDipoleCoulombForce(void)
{
  int i,j,k,l;
  int A1,A2;
  int Type,TypeA;
  VECTOR posA,posB,posA1,posA2,dr;
  REAL r,rr,ri2,cosA,energy,temp,length,ChargeB;
  VECTOR dipoleA,fb1,fa1,fa2,term;
  REAL Bt0,Bt1,Bt2,fac1,fac2;
  REAL DipoleMagnitudeA;
  REAL_MATRIX3x3 v;
  VECTOR rabi;
  REAL SwitchingValueDerivative,SwitchingValue;

  UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleReal[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  Bt0=Bt1=Bt2=0.0;  
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeA=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      posA1=Adsorbates[CurrentSystem][i].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][i].Atoms[A2].Position;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      rabi=dipoleA;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      rabi.x/=length;
      rabi.y/=length;
      rabi.z/=length;

      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          if(i!=k)
          {
            for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
            {
              Type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[Type].HasCharges)
              {
                posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
  
                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                r=sqrt(rr);

                if(rr<CutOffChargeBondDipoleSquared)
                {
                  switch(ChargeMethod)
                  {
                    case NONE:
                      Bt0=Bt1=Bt2=0.0;  
                      break;
                    case TRUNCATED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SMOOTHED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case EWALD:
                      Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          erfc(Alpha[CurrentSystem]*r)/(rr*r);
                      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                      break;
                    default:
                      printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombForce'\n");
                      exit(0);
                      break;
                  }

                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                  energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);

                  term.x=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                  term.y=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                  term.z=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                  fb1.x=term.x;
                  fb1.y=term.y;
                  fb1.z=term.z;

                  fac1=COULOMBIC_CONVERSION_FACTOR*ChargeB*Bt1*DipoleMagnitudeA/length;
                  fac2=COULOMBIC_CONVERSION_FACTOR*ChargeB*Bt1*cosA/(DipoleMagnitudeA*length);

                  fa1.x=-0.5*term.x-fac1*dr.x+fac2*dipoleA.x;
                  fa1.y=-0.5*term.y-fac1*dr.y+fac2*dipoleA.y;
                  fa1.z=-0.5*term.z-fac1*dr.z+fac2*dipoleA.z;

                  fa2.x=-0.5*term.x+fac1*dr.x-fac2*dipoleA.x;
                  fa2.y=-0.5*term.y+fac1*dr.y-fac2*dipoleA.y;
                  fa2.z=-0.5*term.z+fac1*dr.z-fac2*dipoleA.z;

                  switch(ChargeMethod)
                  {
                    case SMOOTHED_COULOMB:
                      if(rr>CutOffChargeBondDipoleSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                        SwitchingValueDerivative=5.0*SwitchingChargeBondDipoleFactors5[5]*rr*rr+4.0*SwitchingChargeBondDipoleFactors5[4]*rr*r+3.0*SwitchingChargeBondDipoleFactors5[3]*rr+
                                                 2.0*SwitchingChargeBondDipoleFactors5[2]*r+SwitchingChargeBondDipoleFactors5[1];

                        SwitchingValueDerivative*=-energy/r;

                        energy*=SwitchingValue;

                        term.x=term.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                        term.y=term.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                        term.z=term.z*SwitchingValue+SwitchingValueDerivative*dr.z;

                        fa1.x=fa1.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                        fa1.y=fa1.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                        fa1.z=fa1.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;
                        fa2.x=fa2.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                        fa2.y=fa2.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                        fa2.z=fa2.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;

                        fb1.x=fb1.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                        fb1.y=fb1.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                        fb1.z=fb1.z*SwitchingValue+SwitchingValueDerivative*dr.z;
                      }
                      break;
                  }

                  UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]-=energy;

                  Adsorbates[CurrentSystem][k].Atoms[l].Force.x+=fb1.x;
                  Adsorbates[CurrentSystem][k].Atoms[l].Force.y+=fb1.y;
                  Adsorbates[CurrentSystem][k].Atoms[l].Force.z+=fb1.z;


                  Adsorbates[CurrentSystem][i].Atoms[A1].Force.x+=fa1.x;
                  Adsorbates[CurrentSystem][i].Atoms[A1].Force.y+=fa1.y;
                  Adsorbates[CurrentSystem][i].Atoms[A1].Force.z+=fa1.z;

                  Adsorbates[CurrentSystem][i].Atoms[A2].Force.x+=fa2.x;
                  Adsorbates[CurrentSystem][i].Atoms[A2].Force.y+=fa2.y;
                  Adsorbates[CurrentSystem][i].Atoms[A2].Force.z+=fa2.z;

                  // convert forces on atoms to molecular virial
                  v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x);
                  v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x);
                  v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x);

                  v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y);
                  v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y);
                  v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y);

                  v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z);
                  v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z);
                  v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z);

                  // the strain derivative
                  StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
                  StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
                  StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

                  StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
                  StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
                  StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

                  StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
                  StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
                  StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
                }
              }
            }
          }
        }
      }

      if(!OmitAdsorbateCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
          {
            Type=Cations[CurrentSystem][k].Atoms[l].Type;
            if(PseudoAtoms[Type].HasCharges)
            {
              posB=Cations[CurrentSystem][k].Atoms[l].Position;
              ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              r=sqrt(rr);

              if(rr<CutOffChargeBondDipoleSquared)
              {
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt0=Bt1=Bt2=0.0;  
                    break;
                  case TRUNCATED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case EWALD:
                    Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                    break;
                  default:
                    printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombForce'\n");
                    exit(0);
                    break;
                }

                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);

                term.x=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                term.y=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                term.z=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                fb1.x=term.x;
                fb1.y=term.y;
                fb1.z=term.z;

                fac1=COULOMBIC_CONVERSION_FACTOR*ChargeB*Bt1*DipoleMagnitudeA/length;
                fac2=COULOMBIC_CONVERSION_FACTOR*ChargeB*Bt1*cosA/(DipoleMagnitudeA*length);

                fa1.x=-0.5*term.x-fac1*dr.x+fac2*dipoleA.x;
                fa1.y=-0.5*term.y-fac1*dr.y+fac2*dipoleA.y;
                fa1.z=-0.5*term.z-fac1*dr.z+fac2*dipoleA.z;

                fa2.x=-0.5*term.x+fac1*dr.x-fac2*dipoleA.x;
                fa2.y=-0.5*term.y+fac1*dr.y-fac2*dipoleA.y;
                fa2.z=-0.5*term.z+fac1*dr.z-fac2*dipoleA.z;

                switch(ChargeMethod)
                {
                  case SMOOTHED_COULOMB:
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                      SwitchingValueDerivative=5.0*SwitchingChargeBondDipoleFactors5[5]*rr*rr+4.0*SwitchingChargeBondDipoleFactors5[4]*rr*r+3.0*SwitchingChargeBondDipoleFactors5[3]*rr+
                                               2.0*SwitchingChargeBondDipoleFactors5[2]*r+SwitchingChargeBondDipoleFactors5[1];

                      SwitchingValueDerivative*=-energy/r;

                      energy*=SwitchingValue;

                      term.x=term.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                      term.y=term.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                      term.z=term.z*SwitchingValue+SwitchingValueDerivative*dr.z;

                      fa1.x=fa1.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                      fa1.y=fa1.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                      fa1.z=fa1.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;
                      fa2.x=fa2.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                      fa2.y=fa2.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                      fa2.z=fa2.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;

                      fb1.x=fb1.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                      fb1.y=fb1.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                      fb1.z=fb1.z*SwitchingValue+SwitchingValueDerivative*dr.z;
                    }
                    break;
                }

                UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=energy;

                Cations[CurrentSystem][k].Atoms[l].Force.x+=fb1.x;
                Cations[CurrentSystem][k].Atoms[l].Force.y+=fb1.y;
                Cations[CurrentSystem][k].Atoms[l].Force.z+=fb1.z;

                Adsorbates[CurrentSystem][i].Atoms[A1].Force.x+=fa1.x;
                Adsorbates[CurrentSystem][i].Atoms[A1].Force.y+=fa1.y;
                Adsorbates[CurrentSystem][i].Atoms[A1].Force.z+=fa1.z;

                Adsorbates[CurrentSystem][i].Atoms[A2].Force.x+=fa2.x;
                Adsorbates[CurrentSystem][i].Atoms[A2].Force.y+=fa2.y;
                Adsorbates[CurrentSystem][i].Atoms[A2].Force.z+=fa2.z;

                // convert forces on atoms to molecular virial
                v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x);
                v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x);
                v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x);

                v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y);
                v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y);
                v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y);

                v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z);
                v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z);
                v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z);

                // the strain derivative
                StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
                StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
                StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

                StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
                StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
                StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

                StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
                StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
                StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
              }
            }
          }
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    TypeA=Cations[CurrentSystem][i].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      posA1=Cations[CurrentSystem][i].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][i].Atoms[A2].Position;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      if(!OmitCationCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          if(i!=k)
          {
            for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
            {
              Type=Cations[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[Type].HasCharges)
              {
                posB=Cations[CurrentSystem][k].Atoms[l].Position;
                ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;
  
                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                r=sqrt(rr);

                if(rr<CutOffChargeBondDipoleSquared)
                {
                  switch(ChargeMethod)
                  {
                    case NONE:
                      Bt0=Bt1=Bt2=0.0;  
                      break;
                    case TRUNCATED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SMOOTHED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case EWALD:
                      Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          erfc(Alpha[CurrentSystem]*r)/(rr*r);
                      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                      break;
                    default:
                      printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombForce'\n");
                      exit(0);
                      break;
                  }

                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                  energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);

                  term.x=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                  term.y=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                  term.z=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                  fb1.x=term.x;
                  fb1.y=term.y;
                  fb1.z=term.z;

                  fac1=COULOMBIC_CONVERSION_FACTOR*ChargeB*Bt1*DipoleMagnitudeA/length;
                  fac2=COULOMBIC_CONVERSION_FACTOR*ChargeB*Bt1*cosA/(DipoleMagnitudeA*length);

                  fa1.x=-0.5*term.x-fac1*dr.x+fac2*dipoleA.x;
                  fa1.y=-0.5*term.y-fac1*dr.y+fac2*dipoleA.y;
                  fa1.z=-0.5*term.z-fac1*dr.z+fac2*dipoleA.z;

                  fa2.x=-0.5*term.x+fac1*dr.x-fac2*dipoleA.x;
                  fa2.y=-0.5*term.y+fac1*dr.y-fac2*dipoleA.y;
                  fa2.z=-0.5*term.z+fac1*dr.z-fac2*dipoleA.z;

                  switch(ChargeMethod)
                  {
                    case SMOOTHED_COULOMB:
                      if(rr>CutOffChargeBondDipoleSwitchSquared)
                      {
                        SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                       SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                        SwitchingValueDerivative=5.0*SwitchingChargeBondDipoleFactors5[5]*rr*rr+4.0*SwitchingChargeBondDipoleFactors5[4]*rr*r+3.0*SwitchingChargeBondDipoleFactors5[3]*rr+
                                                 2.0*SwitchingChargeBondDipoleFactors5[2]*r+SwitchingChargeBondDipoleFactors5[1];

                        SwitchingValueDerivative*=-energy/r;

                        energy*=SwitchingValue;
  
                        term.x=term.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                        term.y=term.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                        term.z=term.z*SwitchingValue+SwitchingValueDerivative*dr.z;

                        fa1.x=fa1.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                        fa1.y=fa1.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                        fa1.z=fa1.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;
                        fa2.x=fa2.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                        fa2.y=fa2.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                        fa2.z=fa2.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;

                        fb1.x=fb1.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                        fb1.y=fb1.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                        fb1.z=fb1.z*SwitchingValue+SwitchingValueDerivative*dr.z;
                      }
                      break;
                  }

                  UCationCationChargeBondDipoleReal[CurrentSystem]-=energy;

                  Cations[CurrentSystem][k].Atoms[l].Force.x+=fb1.x;
                  Cations[CurrentSystem][k].Atoms[l].Force.y+=fb1.y;
                  Cations[CurrentSystem][k].Atoms[l].Force.z+=fb1.z;

                  Cations[CurrentSystem][i].Atoms[A1].Force.x+=fa1.x;
                  Cations[CurrentSystem][i].Atoms[A1].Force.y+=fa1.y;
                  Cations[CurrentSystem][i].Atoms[A1].Force.z+=fa1.z;

                  Cations[CurrentSystem][i].Atoms[A2].Force.x+=fa2.x;
                  Cations[CurrentSystem][i].Atoms[A2].Force.y+=fa2.y;
                  Cations[CurrentSystem][i].Atoms[A2].Force.z+=fa2.z;
 
                  // convert forces on atoms to molecular virial
                  v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x);
                  v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x);
                  v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x);

                  v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y);
                  v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y);
                  v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y);

                  v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z);
                  v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z);
                  v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z);

                  // the strain derivative
                  StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
                  StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
                  StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

                  StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
                  StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
                  StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

                  StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
                  StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
                  StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
                }
              }
            }
          }
        }
      }

      if(!OmitAdsorbateCationCoulombInteractions)
      {
        for(k=0;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
          {
            Type=Adsorbates[CurrentSystem][k].Atoms[l].Type;
            if(PseudoAtoms[Type].HasCharges)
            {
              posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
              ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
  
              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              r=sqrt(rr);

              if(rr<CutOffChargeBondDipoleSquared)
              {
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt0=Bt1=Bt2=0.0;  
                    break;
                  case TRUNCATED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SMOOTHED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case EWALD:
                    Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                    break;
                  default:
                    printf("Unknown charge-method in 'CalculateTotalInterChargeBondDipoleCoulombForce'\n");
                    exit(0);
                    break;
                }

                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeB*cosA);

                term.x=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                term.y=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                term.z=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                fb1.x=term.x;
                fb1.y=term.y;
                fb1.z=term.z;

                fac1=COULOMBIC_CONVERSION_FACTOR*ChargeB*Bt1*DipoleMagnitudeA/length;
                fac2=COULOMBIC_CONVERSION_FACTOR*ChargeB*Bt1*cosA/(DipoleMagnitudeA*length);

                fa1.x=-0.5*term.x-fac1*dr.x+fac2*dipoleA.x;
                fa1.y=-0.5*term.y-fac1*dr.y+fac2*dipoleA.y;
                fa1.z=-0.5*term.z-fac1*dr.z+fac2*dipoleA.z;

                fa2.x=-0.5*term.x+fac1*dr.x-fac2*dipoleA.x;
                fa2.y=-0.5*term.y+fac1*dr.y-fac2*dipoleA.y;
                fa2.z=-0.5*term.z+fac1*dr.z-fac2*dipoleA.z;

                switch(ChargeMethod)
                {
                  case SMOOTHED_COULOMB:
                    if(rr>CutOffChargeBondDipoleSwitchSquared)
                    {
                      SwitchingValue=SwitchingChargeBondDipoleFactors5[5]*(rr*rr*r)+SwitchingChargeBondDipoleFactors5[4]*(rr*rr)+SwitchingChargeBondDipoleFactors5[3]*(rr*r)+
                                     SwitchingChargeBondDipoleFactors5[2]*rr+SwitchingChargeBondDipoleFactors5[1]*r+SwitchingChargeBondDipoleFactors5[0];
                      SwitchingValueDerivative=5.0*SwitchingChargeBondDipoleFactors5[5]*rr*rr+4.0*SwitchingChargeBondDipoleFactors5[4]*rr*r+3.0*SwitchingChargeBondDipoleFactors5[3]*rr+
                                               2.0*SwitchingChargeBondDipoleFactors5[2]*r+SwitchingChargeBondDipoleFactors5[1];

                      SwitchingValueDerivative*=-energy/r;

                      energy*=SwitchingValue;

                      term.x=term.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                      term.y=term.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                      term.z=term.z*SwitchingValue+SwitchingValueDerivative*dr.z;

                      fa1.x=fa1.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                      fa1.y=fa1.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                      fa1.z=fa1.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;
                      fa2.x=fa2.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                      fa2.y=fa2.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                      fa2.z=fa2.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;

                      fb1.x=fb1.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                      fb1.y=fb1.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                      fb1.z=fb1.z*SwitchingValue+SwitchingValueDerivative*dr.z;
                    }
                    break;
                }

                UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=energy;

                Adsorbates[CurrentSystem][k].Atoms[l].Force.x+=fb1.x;
                Adsorbates[CurrentSystem][k].Atoms[l].Force.y+=fb1.y;
                Adsorbates[CurrentSystem][k].Atoms[l].Force.z+=fb1.z;

                Cations[CurrentSystem][i].Atoms[A1].Force.x+=fa1.x;
                Cations[CurrentSystem][i].Atoms[A1].Force.y+=fa1.y;
                Cations[CurrentSystem][i].Atoms[A1].Force.z+=fa1.z;

                Cations[CurrentSystem][i].Atoms[A2].Force.x+=fa2.x;
                Cations[CurrentSystem][i].Atoms[A2].Force.y+=fa2.y;
                Cations[CurrentSystem][i].Atoms[A2].Force.z+=fa2.z;

                // convert forces on atoms to molecular virial
                v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x);
                v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x);
                v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x);

                v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y);
                v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y);
                v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y);

                v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z);
                v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z);
                v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z);

                // the strain derivative
                StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
                StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
                StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

                StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
                StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
                StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;

                StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
                StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
                StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
              }
            }
          }
        }
      }
    }
  }

  return 0;
}

int CalculateTotalInterBondDipoleBondDipoleCoulombForce(void)
{
  int i,j,k,l;
  int A1,A2,B1,B2;
  int TypeA,TypeB;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr;
  REAL r,rr,ri2,rk2,cosAB,cosA,cosB,energy,temp,length;
  VECTOR dipoleA,dipoleB,fb1,fb2,fa1,fa2,termA,termB,term;
  REAL Bt0,Bt1,Bt2,Bt3;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL_MATRIX3x3 v;
  REAL SwitchingValue,SwitchingValueDerivative;

  Bt0=Bt1=Bt2=Bt3=0.0;  
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;
  if(OmitInterMolecularInteractions) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeA=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
    {
      A1=Components[TypeA].BondDipoles[j].A;
      A2=Components[TypeA].BondDipoles[j].B;
      posA1=Adsorbates[CurrentSystem][i].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][i].Atoms[A2].Position;
      DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        for(k=i+1;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          TypeB=Adsorbates[CurrentSystem][k].Type;
          for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
          {  
            B1=Components[TypeB].BondDipoles[l].A;
            B2=Components[TypeB].BondDipoles[l].B;
            posB1=Adsorbates[CurrentSystem][k].Atoms[B1].Position;
            posB2=Adsorbates[CurrentSystem][k].Atoms[B2].Position;
            DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
            dipoleB.x=posB2.x-posB1.x;
            dipoleB.y=posB2.y-posB1.y;
            dipoleB.z=posB2.z-posB1.z;
            posB.x=posB1.x+0.5*dipoleB.x;
            posB.y=posB1.y+0.5*dipoleB.y;
            posB.z=posB1.z+0.5*dipoleB.z;
            rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
            length=sqrt(rk2);
            temp=DipoleMagnitudeB/length;
            dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=Bt3=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                      20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                  break;
                default:
                  printf("Unknown charge-method in 'CalculateTotalInterBondDipoleBondDipoleCoulombForce'\n");
                  exit(0);
                  break;
              }

              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);

              term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
              term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
              term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

              termA.x=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt2*cosA*cosB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt1*DipoleMagnitudeA*dipoleB.x/sqrt(ri2)
                      +Bt2*DipoleMagnitudeA*cosB*dr.x/sqrt(ri2));

              termA.y=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt2*cosA*cosB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt1*DipoleMagnitudeA*dipoleB.y/sqrt(ri2)
                      +Bt2*DipoleMagnitudeA*cosB*dr.y/sqrt(ri2));

              termA.z=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt2*cosA*cosB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt1*DipoleMagnitudeA*dipoleB.z/sqrt(ri2)
                      +Bt2*DipoleMagnitudeA*cosB*dr.z/sqrt(ri2));

              termB.x=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt2*cosA*cosB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt1*DipoleMagnitudeB*dipoleA.x/sqrt(rk2)
                      +Bt2*DipoleMagnitudeB*cosA*dr.x/sqrt(rk2));

              termB.y=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt2*cosA*cosB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt1*DipoleMagnitudeB*dipoleA.y/sqrt(rk2)
                      +Bt2*DipoleMagnitudeB*cosA*dr.y/sqrt(rk2));

              termB.z=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt2*cosA*cosB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt1*DipoleMagnitudeB*dipoleA.z/sqrt(rk2)
                      +Bt2*DipoleMagnitudeB*cosA*dr.z/sqrt(rk2));

              fa1.x=-0.5*term.x-termA.x;
              fa1.y=-0.5*term.y-termA.y;
              fa1.z=-0.5*term.z-termA.z;
              fa2.x=-0.5*term.x+termA.x;
              fa2.y=-0.5*term.y+termA.y;
              fa2.z=-0.5*term.z+termA.z;

              fb1.x=0.5*term.x-termB.x;
              fb1.y=0.5*term.y-termB.y;
              fb1.z=0.5*term.z-termB.z;
              fb2.x=0.5*term.x+termB.x;
              fb2.y=0.5*term.y+termB.y;
              fb2.z=0.5*term.z+termB.z;

              switch(ChargeMethod)
              {
                case SMOOTHED_COULOMB:
                  if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                    SwitchingValueDerivative=5.0*SwitchingBondDipoleBondDipoleFactors5[5]*rr*rr+4.0*SwitchingBondDipoleBondDipoleFactors5[4]*rr*r+3.0*SwitchingBondDipoleBondDipoleFactors5[3]*rr+
                                             2.0*SwitchingBondDipoleBondDipoleFactors5[2]*r+SwitchingBondDipoleBondDipoleFactors5[1];

                    SwitchingValueDerivative*=energy/r;

                    energy*=SwitchingValue;

                    term.x=term.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                    term.y=term.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                    term.z=term.z*SwitchingValue+SwitchingValueDerivative*dr.z;

                    fa1.x=fa1.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                    fa1.y=fa1.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                    fa1.z=fa1.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;
                    fa2.x=fa2.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                    fa2.y=fa2.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                    fa2.z=fa2.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;
            
                    fb1.x=fb1.x*SwitchingValue+0.5*SwitchingValueDerivative*dr.x;
                    fb1.y=fb1.y*SwitchingValue+0.5*SwitchingValueDerivative*dr.y;
                    fb1.z=fb1.z*SwitchingValue+0.5*SwitchingValueDerivative*dr.z;
                    fb2.x=fb2.x*SwitchingValue+0.5*SwitchingValueDerivative*dr.x;
                    fb2.y=fb2.y*SwitchingValue+0.5*SwitchingValueDerivative*dr.y;
                    fb2.z=fb2.z*SwitchingValue+0.5*SwitchingValueDerivative*dr.z;
                    break;
                }
              }

              UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=energy;

              Adsorbates[CurrentSystem][i].Atoms[A1].Force.x+=fa1.x;
              Adsorbates[CurrentSystem][i].Atoms[A1].Force.y+=fa1.y;
              Adsorbates[CurrentSystem][i].Atoms[A1].Force.z+=fa1.z;

              Adsorbates[CurrentSystem][i].Atoms[A2].Force.x+=fa2.x;
              Adsorbates[CurrentSystem][i].Atoms[A2].Force.y+=fa2.y;
              Adsorbates[CurrentSystem][i].Atoms[A2].Force.z+=fa2.z;

              Adsorbates[CurrentSystem][k].Atoms[B1].Force.x+=fb1.x;
              Adsorbates[CurrentSystem][k].Atoms[B1].Force.y+=fb1.y;
              Adsorbates[CurrentSystem][k].Atoms[B1].Force.z+=fb1.z;

              Adsorbates[CurrentSystem][k].Atoms[B2].Force.x+=fb2.x;
              Adsorbates[CurrentSystem][k].Atoms[B2].Force.y+=fb2.y;
              Adsorbates[CurrentSystem][k].Atoms[B2].Force.z+=fb2.z;

              // convert forces on atoms to molecular virial
              // usually this conversion produces a torque on the center of mass and an asymmetric stress
              // by making it symmetric we regain the 'atomic' strain derivative
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
              StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
              StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
              StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;
 
              StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
              StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
              StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;
 
              StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
              StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
              StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
            }
          }
        }
      }

      if(!OmitAdsorbateCationCoulombInteractions)
      {
        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          TypeB=Cations[CurrentSystem][k].Type;
          for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
          {
            B1=Components[TypeB].BondDipoles[l].A;
            B2=Components[TypeB].BondDipoles[l].B;
            posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
            posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
            DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
            dipoleB.x=posB2.x-posB1.x;
            dipoleB.y=posB2.y-posB1.y;
            dipoleB.z=posB2.z-posB1.z;
            posB.x=posB1.x+0.5*dipoleB.x;
            posB.y=posB1.y+0.5*dipoleB.y;
            posB.z=posB1.z+0.5*dipoleB.z;
            rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
            length=sqrt(rk2);
            temp=DipoleMagnitudeB/length;
            dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=Bt3=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                      20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                  break;
                default:
                  printf("Unknown charge-method in 'CalculateTotalInterBondDipoleBondDipoleCoulombForce'\n");
                  exit(0);
                  break;
              }

              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);

              term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
              term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
              term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

              termA.x=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt2*cosA*cosB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt1*DipoleMagnitudeA*dipoleB.x/sqrt(ri2)
                      +Bt2*DipoleMagnitudeA*cosB*dr.x/sqrt(ri2));

              termA.y=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt2*cosA*cosB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt1*DipoleMagnitudeA*dipoleB.y/sqrt(ri2)
                      +Bt2*DipoleMagnitudeA*cosB*dr.y/sqrt(ri2));

              termA.z=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt2*cosA*cosB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt1*DipoleMagnitudeA*dipoleB.z/sqrt(ri2)
                      +Bt2*DipoleMagnitudeA*cosB*dr.z/sqrt(ri2));

              termB.x=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt2*cosA*cosB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt1*DipoleMagnitudeB*dipoleA.x/sqrt(rk2)
                      +Bt2*DipoleMagnitudeB*cosA*dr.x/sqrt(rk2));

              termB.y=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt2*cosA*cosB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt1*DipoleMagnitudeB*dipoleA.y/sqrt(rk2)
                      +Bt2*DipoleMagnitudeB*cosA*dr.y/sqrt(rk2));

              termB.z=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt2*cosA*cosB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt1*DipoleMagnitudeB*dipoleA.z/sqrt(rk2)
                      +Bt2*DipoleMagnitudeB*cosA*dr.z/sqrt(rk2));

              fa1.x=-0.5*term.x-termA.x;
              fa1.y=-0.5*term.y-termA.y;
              fa1.z=-0.5*term.z-termA.z;
              fa2.x=-0.5*term.x+termA.x;
              fa2.y=-0.5*term.y+termA.y;
              fa2.z=-0.5*term.z+termA.z;

              fb1.x=0.5*term.x-termB.x;
              fb1.y=0.5*term.y-termB.y;
              fb1.z=0.5*term.z-termB.z;
              fb2.x=0.5*term.x+termB.x;
              fb2.y=0.5*term.y+termB.y;
              fb2.z=0.5*term.z+termB.z;

              switch(ChargeMethod)
              {
                case SMOOTHED_COULOMB:
                  if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                    SwitchingValueDerivative=5.0*SwitchingBondDipoleBondDipoleFactors5[5]*rr*rr+4.0*SwitchingBondDipoleBondDipoleFactors5[4]*rr*r+3.0*SwitchingBondDipoleBondDipoleFactors5[3]*rr+
                                             2.0*SwitchingBondDipoleBondDipoleFactors5[2]*r+SwitchingBondDipoleBondDipoleFactors5[1];

                    SwitchingValueDerivative*=energy/r;

                    energy*=SwitchingValue;

                    term.x=term.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                    term.y=term.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                    term.z=term.z*SwitchingValue+SwitchingValueDerivative*dr.z;

                    fa1.x=fa1.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                    fa1.y=fa1.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                    fa1.z=fa1.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;
                    fa2.x=fa2.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                    fa2.y=fa2.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                    fa2.z=fa2.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;

                    fb1.x=fb1.x*SwitchingValue+0.5*SwitchingValueDerivative*dr.x;
                    fb1.y=fb1.y*SwitchingValue+0.5*SwitchingValueDerivative*dr.y;
                    fb1.z=fb1.z*SwitchingValue+0.5*SwitchingValueDerivative*dr.z;
                    fb2.x=fb2.x*SwitchingValue+0.5*SwitchingValueDerivative*dr.x;
                    fb2.y=fb2.y*SwitchingValue+0.5*SwitchingValueDerivative*dr.y;
                    fb2.z=fb2.z*SwitchingValue+0.5*SwitchingValueDerivative*dr.z;
                    break;
                }
              }

              UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=energy;

              Adsorbates[CurrentSystem][i].Atoms[A1].Force.x+=fa1.x;
              Adsorbates[CurrentSystem][i].Atoms[A1].Force.y+=fa1.y;
              Adsorbates[CurrentSystem][i].Atoms[A1].Force.z+=fa1.z;

              Adsorbates[CurrentSystem][i].Atoms[A2].Force.x+=fa2.x;
              Adsorbates[CurrentSystem][i].Atoms[A2].Force.y+=fa2.y;
              Adsorbates[CurrentSystem][i].Atoms[A2].Force.z+=fa2.z;

              Cations[CurrentSystem][k].Atoms[B1].Force.x+=fb1.x;
              Cations[CurrentSystem][k].Atoms[B1].Force.y+=fb1.y;
              Cations[CurrentSystem][k].Atoms[B1].Force.z+=fb1.z;

              Cations[CurrentSystem][k].Atoms[B2].Force.x+=fb2.x;
              Cations[CurrentSystem][k].Atoms[B2].Force.y+=fb2.y;
              Cations[CurrentSystem][k].Atoms[B2].Force.z+=fb2.z;

              // convert forces on atoms to molecular virial
              // usually this conversion produces a torque on the center of mass and an asymmetric stress
              // by making it symmetric we regain the 'atomic' strain derivative
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
              StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
              StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
              StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

              StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
              StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
              StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;
 
              StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
              StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
              StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
            }
          }
        }
      }
    }
  }

  if(!OmitCationCationCoulombInteractions)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      TypeA=Cations[CurrentSystem][i].Type;
      for(j=0;j<Components[TypeA].NumberOfBondDipoles;j++)
      {
        A1=Components[TypeA].BondDipoles[j].A;
        A2=Components[TypeA].BondDipoles[j].B;
        posA1=Cations[CurrentSystem][i].Atoms[A1].Position;
        posA2=Cations[CurrentSystem][i].Atoms[A2].Position;
        DipoleMagnitudeA=Components[TypeA].BondDipoleMagnitude[j];
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
        length=sqrt(ri2);
        temp=DipoleMagnitudeA/length;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        for(k=i+1;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          TypeB=Cations[CurrentSystem][k].Type;
          for(l=0;l<Components[TypeB].NumberOfBondDipoles;l++)
          {
            B1=Components[TypeB].BondDipoles[l].A;
            B2=Components[TypeB].BondDipoles[l].B;
            posB1=Cations[CurrentSystem][k].Atoms[B1].Position;
            posB2=Cations[CurrentSystem][k].Atoms[B2].Position;
            DipoleMagnitudeB=Components[TypeB].BondDipoleMagnitude[l];
            dipoleB.x=posB2.x-posB1.x;
            dipoleB.y=posB2.y-posB1.y;
            dipoleB.z=posB2.z-posB1.z;
            posB.x=posB1.x+0.5*dipoleB.x;
            posB.y=posB1.y+0.5*dipoleB.y;
            posB.z=posB1.z+0.5*dipoleB.z;
            rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
            length=sqrt(rk2);
            temp=DipoleMagnitudeB/length;
            dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            if(rr<CutOffBondDipoleBondDipoleSquared)
            {
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=Bt3=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case SMOOTHED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                      20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                  break;
                default:
                  printf("Unknown charge-method in 'CalculateTotalInterBondDipoleBondDipoleCoulombForce'\n");
                  exit(0);
                  break;
              }

              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);

              term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
              term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
              term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

              termA.x=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt2*cosA*cosB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt1*DipoleMagnitudeA*dipoleB.x/sqrt(ri2)
                      +Bt2*DipoleMagnitudeA*cosB*dr.x/sqrt(ri2));

              termA.y=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt2*cosA*cosB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt1*DipoleMagnitudeA*dipoleB.y/sqrt(ri2)
                      +Bt2*DipoleMagnitudeA*cosB*dr.y/sqrt(ri2));

              termA.z=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt2*cosA*cosB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
                      -Bt1*DipoleMagnitudeA*dipoleB.z/sqrt(ri2)
                      +Bt2*DipoleMagnitudeA*cosB*dr.z/sqrt(ri2));

              termB.x=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt2*cosA*cosB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt1*DipoleMagnitudeB*dipoleA.x/sqrt(rk2)
                      +Bt2*DipoleMagnitudeB*cosA*dr.x/sqrt(rk2));

              termB.y=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt2*cosA*cosB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt1*DipoleMagnitudeB*dipoleA.y/sqrt(rk2)
                      +Bt2*DipoleMagnitudeB*cosA*dr.y/sqrt(rk2));

              termB.z=COULOMBIC_CONVERSION_FACTOR*(
                      Bt1*cosAB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt2*cosA*cosB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
                      -Bt1*DipoleMagnitudeB*dipoleA.z/sqrt(rk2)
                      +Bt2*DipoleMagnitudeB*cosA*dr.z/sqrt(rk2));

              fa1.x=-0.5*term.x-termA.x;
              fa1.y=-0.5*term.y-termA.y;
              fa1.z=-0.5*term.z-termA.z;
              fa2.x=-0.5*term.x+termA.x;
              fa2.y=-0.5*term.y+termA.y;
              fa2.z=-0.5*term.z+termA.z;

              fb1.x=0.5*term.x-termB.x;
              fb1.y=0.5*term.y-termB.y;
              fb1.z=0.5*term.z-termB.z;
              fb2.x=0.5*term.x+termB.x;
              fb2.y=0.5*term.y+termB.y;
              fb2.z=0.5*term.z+termB.z;

              switch(ChargeMethod)
              {
                case SMOOTHED_COULOMB:
                  if(rr>CutOffBondDipoleBondDipoleSwitchSquared)
                  {
                    SwitchingValue=SwitchingBondDipoleBondDipoleFactors5[5]*(rr*rr*r)+SwitchingBondDipoleBondDipoleFactors5[4]*(rr*rr)+SwitchingBondDipoleBondDipoleFactors5[3]*(rr*r)+
                                   SwitchingBondDipoleBondDipoleFactors5[2]*rr+SwitchingBondDipoleBondDipoleFactors5[1]*r+SwitchingBondDipoleBondDipoleFactors5[0];
                    SwitchingValueDerivative=5.0*SwitchingBondDipoleBondDipoleFactors5[5]*rr*rr+4.0*SwitchingBondDipoleBondDipoleFactors5[4]*rr*r+3.0*SwitchingBondDipoleBondDipoleFactors5[3]*rr+
                                             2.0*SwitchingBondDipoleBondDipoleFactors5[2]*r+SwitchingBondDipoleBondDipoleFactors5[1];

                    SwitchingValueDerivative*=energy/r;

                    energy*=SwitchingValue;

                    term.x=term.x*SwitchingValue+SwitchingValueDerivative*dr.x;
                    term.y=term.y*SwitchingValue+SwitchingValueDerivative*dr.y;
                    term.z=term.z*SwitchingValue+SwitchingValueDerivative*dr.z;

                    fa1.x=fa1.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                    fa1.y=fa1.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                    fa1.z=fa1.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;
                    fa2.x=fa2.x*SwitchingValue-0.5*SwitchingValueDerivative*dr.x;
                    fa2.y=fa2.y*SwitchingValue-0.5*SwitchingValueDerivative*dr.y;
                    fa2.z=fa2.z*SwitchingValue-0.5*SwitchingValueDerivative*dr.z;

                    fb1.x=fb1.x*SwitchingValue+0.5*SwitchingValueDerivative*dr.x;
                    fb1.y=fb1.y*SwitchingValue+0.5*SwitchingValueDerivative*dr.y;
                    fb1.z=fb1.z*SwitchingValue+0.5*SwitchingValueDerivative*dr.z;
                    fb2.x=fb2.x*SwitchingValue+0.5*SwitchingValueDerivative*dr.x;
                    fb2.y=fb2.y*SwitchingValue+0.5*SwitchingValueDerivative*dr.y;
                    fb2.z=fb2.z*SwitchingValue+0.5*SwitchingValueDerivative*dr.z;
                    break;
                }
              }

              UCationCationBondDipoleBondDipoleReal[CurrentSystem]+=energy;

              Cations[CurrentSystem][i].Atoms[A1].Force.x+=fa1.x;
              Cations[CurrentSystem][i].Atoms[A1].Force.y+=fa1.y;
              Cations[CurrentSystem][i].Atoms[A1].Force.z+=fa1.z;

              Cations[CurrentSystem][i].Atoms[A2].Force.x+=fa2.x;
              Cations[CurrentSystem][i].Atoms[A2].Force.y+=fa2.y;
              Cations[CurrentSystem][i].Atoms[A2].Force.z+=fa2.z;

              Cations[CurrentSystem][k].Atoms[B1].Force.x+=fb1.x;
              Cations[CurrentSystem][k].Atoms[B1].Force.y+=fb1.y;
              Cations[CurrentSystem][k].Atoms[B1].Force.z+=fb1.z;

              Cations[CurrentSystem][k].Atoms[B2].Force.x+=fb2.x;
              Cations[CurrentSystem][k].Atoms[B2].Force.y+=fb2.y;
              Cations[CurrentSystem][k].Atoms[B2].Force.z+=fb2.z;

              // convert forces on atoms to molecular virial
              // usually this conversion produces a torque on the center of mass and an asymmetric stress
              // by making it symmetric we regain the 'atomic' strain derivative
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
              StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x-v.ax;
              StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x-v.bx;
              StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x-v.cx;

              StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y-v.ay;
              StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y-v.by;
              StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y-v.cy;
 
              StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z-v.az;
              StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z-v.bz;
              StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z-v.cz;
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateTotalInterChargeChargeCoulombElectricFieldMC(int New,int excl_ads,int excl_cation)
{
  int i,j,k,l;
  int typeA,typeB;
  REAL r,rr;
  REAL chargeA,chargeB;
  REAL force_factor,force_factor_A,force_factor_B;
  VECTOR posA,posB,dr;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
    {
      if(i!=excl_ads)
      { 
        typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
        posA=Adsorbates[CurrentSystem][i].Atoms[k].Position;
        chargeA=Adsorbates[CurrentSystem][i].Atoms[k].Charge;

        // loop over adsorbant molecules
        if(!OmitAdsorbateAdsorbatePolarization)
        {
          for(j=i+1;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
          {
            for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
            {
              if(j!=excl_ads)
              {
                posB=Adsorbates[CurrentSystem][j].Atoms[l].Position;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeChargeSquared)
                {
                  r=sqrt(rr);
                  typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
                  chargeB=Adsorbates[CurrentSystem][j].Atoms[l].Charge;

                  switch(ChargeMethod)
                  {
                    case NONE:
                      force_factor=0.0;
                      break;
                    case SHIFTED_COULOMB:
                    case TRUNCATED_COULOMB:
                      force_factor=-COULOMBIC_CONVERSION_FACTOR/(rr*r);
                      break;
                    case WOLFS_METHOD_DAMPED_FG:
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*(erfc(Alpha[CurrentSystem]*r)/rr+
                                     2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                    -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                      2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                      break;
                    case EWALD:
                    default:
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*
                            (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                            (r*rr);
                      break;
                  }
                  force_factor_A=force_factor*chargeB;
                  force_factor_B=force_factor*chargeA;

                  // electric field
                  Adsorbates[CurrentSystem][i].Atoms[k].ElectricField.x-=force_factor_A*dr.x;
                  Adsorbates[CurrentSystem][i].Atoms[k].ElectricField.y-=force_factor_A*dr.y;
                  Adsorbates[CurrentSystem][i].Atoms[k].ElectricField.z-=force_factor_A*dr.z;

                  Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.x+=force_factor_B*dr.x;
                  Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.y+=force_factor_B*dr.y;
                  Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.z+=force_factor_B*dr.z;
                }
              }
            }
          }
        }

        if(!OmitAdsorbateCationPolarization)
        {
          for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
          {
            for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
            {
              if(j!=excl_cation)
              {
                posB=Cations[CurrentSystem][j].Atoms[l].Position;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffChargeChargeSquared)
                {
                  typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                  chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;

                  r=sqrt(rr);

                  switch(ChargeMethod)
                  {
                    case NONE:
                      force_factor=0.0;
                      break;
                    case SHIFTED_COULOMB:
                    case TRUNCATED_COULOMB:
                      force_factor=-COULOMBIC_CONVERSION_FACTOR/(rr*r);
                      break;
                    case WOLFS_METHOD_DAMPED_FG:
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*(erfc(Alpha[CurrentSystem]*r)/rr+
                                     2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                    -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                      2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                      break;
                    case EWALD:
                    default:
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*
                            (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                            (r*rr);
                      break;
                  }
                  force_factor_A=force_factor*chargeB;
                  force_factor_B=force_factor*chargeA;

                  // electric field
                  Adsorbates[CurrentSystem][i].Atoms[k].ElectricField.x-=force_factor_A*dr.x;
                  Adsorbates[CurrentSystem][i].Atoms[k].ElectricField.y-=force_factor_A*dr.y;
                  Adsorbates[CurrentSystem][i].Atoms[k].ElectricField.z-=force_factor_A*dr.z;

                  Cations[CurrentSystem][j].Atoms[l].ElectricField.x+=force_factor_B*dr.x;
                  Cations[CurrentSystem][j].Atoms[l].ElectricField.y+=force_factor_B*dr.y;
                  Cations[CurrentSystem][j].Atoms[l].ElectricField.z+=force_factor_B*dr.z;
                }
              }
            }
          }
        }
      }
    }
  }

  if(!OmitCationCationPolarization)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
      {
        if(i!=excl_cation)
        {
          typeA=Cations[CurrentSystem][i].Atoms[k].Type;
          posA=Cations[CurrentSystem][i].Atoms[k].Position;
          chargeA=Cations[CurrentSystem][i].Atoms[k].Charge;

          // loop over cation molecules
          for(j=i+1;j<NumberOfCationMolecules[CurrentSystem];j++)
          {
            for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
            {
              if(j!=excl_cation)
              {
                posB=Cations[CurrentSystem][j].Atoms[l].Position;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffChargeChargeSquared)
                {
                  typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                  chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;

                  r=sqrt(rr);

                  switch(ChargeMethod)
                  {
                    case NONE:
                      force_factor=0.0;
                      break;
                    case SHIFTED_COULOMB:
                    case TRUNCATED_COULOMB:
                      force_factor=-COULOMBIC_CONVERSION_FACTOR/(rr*r);
                      break;
                    case WOLFS_METHOD_DAMPED_FG:
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*(erfc(Alpha[CurrentSystem]*r)/rr+
                                     2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                    -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                      2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                      break;
                    case EWALD:
                    default:
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*
                             (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                             (r*rr);
                      break;
                  }
                  force_factor_A=force_factor*chargeB;
                  force_factor_B=force_factor*chargeA;

                  // electric field
                  Cations[CurrentSystem][i].Atoms[k].ElectricField.x-=force_factor_A*dr.x;
                  Cations[CurrentSystem][i].Atoms[k].ElectricField.y-=force_factor_A*dr.y;
                  Cations[CurrentSystem][i].Atoms[k].ElectricField.z-=force_factor_A*dr.z;

                  Cations[CurrentSystem][j].Atoms[l].ElectricField.x+=force_factor_B*dr.x;
                  Cations[CurrentSystem][j].Atoms[l].ElectricField.y+=force_factor_B*dr.y;
                  Cations[CurrentSystem][j].Atoms[l].ElectricField.z+=force_factor_B*dr.z;
                }
              }
            }
          }
        }
      }
    }
  }

  if(New)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      posA=TrialPosition[CurrentSystem][i];
      typeA=Components[CurrentComponent].Type[i];
      chargeA=Components[CurrentComponent].Charge[i];

      // loop over adsorbate molecules
      if(((!Components[CurrentComponent].ExtraFrameworkMolecule)&&(!OmitAdsorbateAdsorbatePolarization))||
         (Components[CurrentComponent].ExtraFrameworkMolecule&&(!OmitAdsorbateCationPolarization)))
      {
        for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
        {
          for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
          {
            if(j!=excl_ads)
            {
              posB=Adsorbates[CurrentSystem][j].Atoms[l].Position;

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared)
              {
                r=sqrt(rr);
                typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
                chargeB=Adsorbates[CurrentSystem][j].Atoms[l].Charge;

                switch(ChargeMethod)
                {
                  case NONE:
                    force_factor=0.0;
                    break;
                  case SHIFTED_COULOMB:
                  case TRUNCATED_COULOMB:
                    force_factor=-COULOMBIC_CONVERSION_FACTOR/(rr*r);
                    break;
                  case WOLFS_METHOD_DAMPED_FG:
                    force_factor=-COULOMBIC_CONVERSION_FACTOR*(erfc(Alpha[CurrentSystem]*r)/rr+
                                   2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                  -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                    2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                    break;
                  case EWALD:
                  default:
                    force_factor=-COULOMBIC_CONVERSION_FACTOR*
                          (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                          (r*rr);
                    break;
                }
                force_factor_A=force_factor*chargeB;
                force_factor_B=force_factor*chargeA;

                // electric field
                ElectricFieldAtTrialPosition[CurrentSystem][i].x-=force_factor_A*dr.x;
                ElectricFieldAtTrialPosition[CurrentSystem][i].y-=force_factor_A*dr.y;
                ElectricFieldAtTrialPosition[CurrentSystem][i].z-=force_factor_A*dr.z;

                Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.x+=force_factor_B*dr.x;
                Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.y+=force_factor_B*dr.y;
                Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.z+=force_factor_B*dr.z;
              }
            }
          }
        }
      }

      // loop over cation molecules
      if(((Components[CurrentComponent].ExtraFrameworkMolecule)&&(!OmitCationCationPolarization))||
         (!Components[CurrentComponent].ExtraFrameworkMolecule&&(!OmitAdsorbateCationPolarization)))
      {
        for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
        {
          for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
          {
            if(j!=excl_cation)
            {
              posB=Cations[CurrentSystem][j].Atoms[l].Position;

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared)
              {
                r=sqrt(rr);
                typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;

                switch(ChargeMethod)
                {
                  case NONE:
                    force_factor=0.0;
                    break;
                  case SHIFTED_COULOMB:
                  case TRUNCATED_COULOMB:
                    force_factor=-COULOMBIC_CONVERSION_FACTOR/(rr*r);
                    break;
                  case WOLFS_METHOD_DAMPED_FG:
                    force_factor=-COULOMBIC_CONVERSION_FACTOR*(erfc(Alpha[CurrentSystem]*r)/rr+
                                   2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                  -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                    2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                    break;
                  case EWALD:
                  default:
                    force_factor=-COULOMBIC_CONVERSION_FACTOR*
                          (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                          (r*rr);
                    break;
                }
                force_factor_A=force_factor*chargeB;
                force_factor_B=force_factor*chargeA;

                // electric field
                ElectricFieldAtTrialPosition[CurrentSystem][i].x-=force_factor_A*dr.x;
                ElectricFieldAtTrialPosition[CurrentSystem][i].y-=force_factor_A*dr.y;
                ElectricFieldAtTrialPosition[CurrentSystem][i].z-=force_factor_A*dr.z;

                Cations[CurrentSystem][j].Atoms[l].ElectricField.x+=force_factor_B*dr.x;
                Cations[CurrentSystem][j].Atoms[l].ElectricField.y+=force_factor_B*dr.y;
                Cations[CurrentSystem][j].Atoms[l].ElectricField.z+=force_factor_B*dr.z;
              }
            }
          }
        }
      }
    }
  }

  return 0;
}

int CalculateInterElectricFieldFromInducedDipoles(void)
{
  int i,j,k,l;
  int TypeA,TypeB;
  VECTOR posA,posB,dr;
  REAL r,rr,cosA,cosB;
  VECTOR dipoleA,dipoleB,term;
  REAL Bt0,Bt1,Bt2,Bt3;
  REAL_MATRIX3x3 v;

  Bt0=Bt1=Bt2=Bt3=0.0;  
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeA=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
      dipoleA=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;

      if(!OmitAdsorbateAdsorbatePolarization)
      {
        for(k=i+1;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          TypeB=Adsorbates[CurrentSystem][k].Type;
          for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
          {
            posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
            dipoleB=Adsorbates[CurrentSystem][k].Atoms[l].InducedDipole;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            if(rr<CutOffChargeChargeSquared)
            {
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=Bt3=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case EWALD:
                default:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                      20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                  break;
              }

              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
              term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
              term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
              term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

              Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x+=term.x;
              Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y+=term.y;
              Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z+=term.z;

              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
              term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
              term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

              Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.x+=term.x;
              Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.y+=term.y;
              Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.z+=term.z;
            }
          }
        }
      }

      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        TypeB=Cations[CurrentSystem][k].Type;
        for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
        {
          posB=Cations[CurrentSystem][k].Atoms[l].Position;
          dipoleB=Cations[CurrentSystem][k].Atoms[l].InducedDipole;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          r=sqrt(rr);

          if(rr<CutOffChargeChargeSquared)
          {
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=Bt3=0.0;  
                break;
              case EWALD:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                    20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                break;
              case TRUNCATED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                Bt3=15.0/(r*rr*rr*rr);
                break;
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                Bt3=15.0/(r*rr*rr*rr);
                break;
            }

            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

            term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
            term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
            term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

            Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x+=term.x;
            Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y+=term.y;
            Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z+=term.z;

            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

            term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
            term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
            term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

            Cations[CurrentSystem][k].Atoms[l].ElectricField.x+=term.x;
            Cations[CurrentSystem][k].Atoms[l].ElectricField.y+=term.y;
            Cations[CurrentSystem][k].Atoms[l].ElectricField.z+=term.z;
          }
        }
      }
    }
  }

  if(!OmitCationCationPolarization)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      TypeA=Cations[CurrentSystem][i].Type;
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        posA=Cations[CurrentSystem][i].Atoms[j].Position;
        dipoleA=Cations[CurrentSystem][i].Atoms[j].InducedDipole;

        for(k=i+1;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          TypeB=Cations[CurrentSystem][k].Type;
          for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
          {
            posB=Cations[CurrentSystem][k].Atoms[l].Position;
            dipoleB=Cations[CurrentSystem][k].Atoms[l].InducedDipole;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            if(rr<CutOffChargeChargeSquared)
            {
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=Bt3=0.0;  
                  break;
                case EWALD:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                  Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                      20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                  break;
                case TRUNCATED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
              }

              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

              term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
              term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
              term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.z-Bt1*dipoleA.z);
    
              Cations[CurrentSystem][i].Atoms[j].ElectricField.x+=term.x;
              Cations[CurrentSystem][i].Atoms[j].ElectricField.y+=term.y;
              Cations[CurrentSystem][i].Atoms[j].ElectricField.z+=term.z;
  
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
    
              term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
              term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
              term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

              Cations[CurrentSystem][k].Atoms[l].ElectricField.x+=term.x;
              Cations[CurrentSystem][k].Atoms[l].ElectricField.y+=term.y;
              Cations[CurrentSystem][k].Atoms[l].ElectricField.z+=term.z;
            }
          }
        }
      }
    }
  }

  return 0;
}


// NEW implementation (work in progress)
int CalculateInterElectricFieldFromInducedDipoleMC(int New,int excl_ads,int excl_cation)
{
  int i,j,k,l;
  int TypeA,TypeB;
  VECTOR posA,posB,dr;
  REAL r,rr,cosA,cosB;
  VECTOR dipoleA,dipoleB,term;
  REAL Bt0,Bt1,Bt2,Bt3;
  REAL_MATRIX3x3 v;

  Bt0=Bt1=Bt2=Bt3=0.0;  
  v.ax=v.bx=v.cx=0.0;
  v.ay=v.by=v.cy=0.0;
  v.az=v.bz=v.cz=0.0;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    if(i!=excl_ads)
    {
      TypeA=Adsorbates[CurrentSystem][i].Type;
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
        dipoleA=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;

        if(!OmitAdsorbateAdsorbatePolarization)
        {
          for(k=i+1;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
          {
            if(k!=excl_ads)
            {
              TypeB=Adsorbates[CurrentSystem][k].Type;
              for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
              {
                posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                dipoleB=Adsorbates[CurrentSystem][k].Atoms[l].InducedDipole;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                r=sqrt(rr);

                if(rr<CutOffChargeChargeSquared)
                {
                  switch(ChargeMethod)
                  {
                    case NONE:
                      Bt0=Bt1=Bt2=Bt3=0.0;  
                      break;
                    case TRUNCATED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      Bt3=15.0/(r*rr*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      Bt3=15.0/(r*rr*rr*rr);
                      break;
                    case EWALD:
                    default:
                      Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          erfc(Alpha[CurrentSystem]*r)/(rr*r);
                      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                      Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                          20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                          8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                      break;
                  }

                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                  term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
                  term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
                  term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

                  Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x+=term.x;
                  Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y+=term.y;
                  Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z+=term.z;

                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                  term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                  term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                  term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                  Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.x+=term.x;
                  Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.y+=term.y;
                  Adsorbates[CurrentSystem][k].Atoms[l].ElectricField.z+=term.z;
                }
              }
            }
          }

          if(!OmitAdsorbateCationPolarization)
          {
            for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
            {
              if(k!=excl_cation)
              {
                TypeB=Cations[CurrentSystem][k].Type;
                for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
                {
                  posB=Cations[CurrentSystem][k].Atoms[l].Position;
                  dipoleB=Cations[CurrentSystem][k].Atoms[l].InducedDipole;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                  r=sqrt(rr);

                  if(rr<CutOffChargeChargeSquared)
                  {
                    switch(ChargeMethod)
                    {
                      case NONE:
                        Bt0=Bt1=Bt2=Bt3=0.0;  
                        break;
                      case EWALD:
                        Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                        Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                            erfc(Alpha[CurrentSystem]*r)/(rr*r);
                        Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                            4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                            3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                        Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                            20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                            8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                            15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                        break;
                      case TRUNCATED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        Bt3=15.0/(r*rr*rr*rr);
                        break;
                      case SHIFTED_COULOMB:
                        Bt0=1.0/(r);
                        Bt1=1.0/(r*rr);
                        Bt2=3.0/(r*rr*rr);
                        Bt3=15.0/(r*rr*rr*rr);
                        break;
                    }

                    cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

                    term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                    term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                    term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                    Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.x+=term.x;
                    Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.y+=term.y;
                    Adsorbates[CurrentSystem][i].Atoms[j].ElectricField.z+=term.z;

                    cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

                    term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
                    term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
                    term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

                    Cations[CurrentSystem][k].Atoms[l].ElectricField.x+=term.x;
                    Cations[CurrentSystem][k].Atoms[l].ElectricField.y+=term.y;
                    Cations[CurrentSystem][k].Atoms[l].ElectricField.z+=term.z;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if(!OmitCationCationPolarization)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      if(i!=excl_cation)
      {
        TypeA=Cations[CurrentSystem][i].Type;
        for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
        {
          posA=Cations[CurrentSystem][i].Atoms[j].Position;
          dipoleA=Cations[CurrentSystem][i].Atoms[j].InducedDipole;

          for(k=i+1;k<NumberOfCationMolecules[CurrentSystem];k++)
          {
            if(k!=excl_cation)
            {
              TypeB=Cations[CurrentSystem][k].Type;
              for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
              {
                posB=Cations[CurrentSystem][k].Atoms[l].Position;
                dipoleB=Cations[CurrentSystem][k].Atoms[l].InducedDipole;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                r=sqrt(rr);

                if(rr<CutOffChargeChargeSquared)
                {
                  switch(ChargeMethod)
                  {
                    case NONE:
                      Bt0=Bt1=Bt2=Bt3=0.0;  
                      break;
                    case EWALD:
                      Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          erfc(Alpha[CurrentSystem]*r)/(rr*r);
                      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                      Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                          20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                          8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                      break;
                    case TRUNCATED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      Bt3=15.0/(r*rr*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      Bt3=15.0/(r*rr*rr*rr);
                      break;
                  }

                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

                  term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                  term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                  term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.z-Bt1*dipoleA.z);
    
                  Cations[CurrentSystem][i].Atoms[j].ElectricField.x+=term.x;
                  Cations[CurrentSystem][i].Atoms[j].ElectricField.y+=term.y;
                  Cations[CurrentSystem][i].Atoms[j].ElectricField.z+=term.z;
  
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
    
                  term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
                  term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
                  term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

                  Cations[CurrentSystem][k].Atoms[l].ElectricField.x+=term.x;
                  Cations[CurrentSystem][k].Atoms[l].ElectricField.y+=term.y;
                  Cations[CurrentSystem][k].Atoms[l].ElectricField.z+=term.z;
                }
              }
            }
          }
        }
      }
    }
  }

  if(New)
  {
    for(i=0;i<Components[CurrentComponent].NumberOfAtoms;i++)
    {
      posA=TrialPosition[CurrentSystem][i];
      dipoleA=InducedDipoleAtTrialPosition[CurrentSystem][i];

      // loop over adsorbate molecules
      if(((!Components[CurrentComponent].ExtraFrameworkMolecule)&&(!OmitAdsorbateAdsorbatePolarization))||
         (Components[CurrentComponent].ExtraFrameworkMolecule&&(!OmitAdsorbateCationPolarization)))
      {
        for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
        {
          if(j!=excl_ads)
          {
            for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
            {
              posB=Adsorbates[CurrentSystem][j].Atoms[l].Position;
              dipoleB=Adsorbates[CurrentSystem][j].Atoms[l].InducedDipole;

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              r=sqrt(rr);

              if(rr<CutOffChargeChargeSquared)
              {
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt0=Bt1=Bt2=Bt3=0.0;
                    break;
                  case EWALD:
                    Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                    Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                        20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                    break;
                  case TRUNCATED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    Bt3=15.0/(r*rr*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    Bt3=15.0/(r*rr*rr*rr);
                    break;
                }

                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
                term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
                term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

                ElectricFieldAtTrialPosition[CurrentSystem][i].x+=term.x;
                ElectricFieldAtTrialPosition[CurrentSystem][i].y+=term.y;
                ElectricFieldAtTrialPosition[CurrentSystem][i].z+=term.z;

                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.x+=term.x;
                Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.y+=term.y;
                Adsorbates[CurrentSystem][j].Atoms[l].ElectricField.z+=term.z;
              }
            }
          }
        }
      }

      // loop over cation molecules
      if(((Components[CurrentComponent].ExtraFrameworkMolecule)&&(!OmitCationCationPolarization))||
         (!Components[CurrentComponent].ExtraFrameworkMolecule&&(!OmitAdsorbateCationPolarization)))
      {
        for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
        {
          for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
          {
            if(j!=excl_cation)
            {
              posB=Cations[CurrentSystem][j].Atoms[l].Position;
              dipoleB=Cations[CurrentSystem][j].Atoms[l].InducedDipole;

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              r=sqrt(rr);

              if(rr<CutOffChargeChargeSquared)
              {
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt0=Bt1=Bt2=Bt3=0.0;
                    break;
                  case EWALD:
                    Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                    Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                        20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                    break;
                  case TRUNCATED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    Bt3=15.0/(r*rr*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    Bt3=15.0/(r*rr*rr*rr);
                    break;
                }

                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
                term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
                term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
                term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

                ElectricFieldAtTrialPosition[CurrentSystem][i].x+=term.x;
                ElectricFieldAtTrialPosition[CurrentSystem][i].y+=term.y;
                ElectricFieldAtTrialPosition[CurrentSystem][i].z+=term.z;

                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
                term.x=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                term.y=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                term.z=COULOMBIC_CONVERSION_FACTOR*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                Cations[CurrentSystem][j].Atoms[l].ElectricField.x+=term.x;
                Cations[CurrentSystem][j].Atoms[l].ElectricField.y+=term.y;
                Cations[CurrentSystem][j].Atoms[l].ElectricField.z+=term.z;
              }
            }
          }
        }
      }
    }
  }

  return 0;
}


int CalculateInterChargeInducedDipoleForce(void)
{
  int i,j,k,l;
  int TypeA,TypeB;
  VECTOR posA,posB,dr;
  REAL r,rr,cosA,cosB,ChargeB,ChargeA;
  VECTOR dipoleA,dipoleB,term;
  REAL Bt0,Bt1,Bt2;

  if(ChargeMethod==NONE) return 0;

  Bt0=Bt1=Bt2=0.0;  

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      TypeA=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      //if(PseudoAtoms[TypeA].HasCharges)
      {
        posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
        dipoleA=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;
        ChargeA=Adsorbates[CurrentSystem][i].Atoms[j].Charge;

        if(!OmitAdsorbateAdsorbatePolarization)
        {
          for(k=i+1;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
          {
            for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
            {
              TypeB=Adsorbates[CurrentSystem][k].Atoms[l].Type;
              //if(PseudoAtoms[TypeB].HasCharges)
              {
                posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
                dipoleB=Adsorbates[CurrentSystem][k].Atoms[l].InducedDipole;
                ChargeB=Adsorbates[CurrentSystem][k].Atoms[l].Charge;
   
                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeChargeSquared)
                {
                  r=sqrt(rr);
                  switch(ChargeMethod)
                  {
                    case NONE:
                      Bt0=Bt1=Bt2=0.0;  
                      break;
                    case TRUNCATED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case EWALD:
                    default:
                      Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          erfc(Alpha[CurrentSystem]*r)/(rr*r);
                      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                      break;
                  }
  
                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

                  term.x=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                  term.y=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                  term.z=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                  Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=term.x;
                  Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=term.y;
                  Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=term.z;

                  Adsorbates[CurrentSystem][k].Atoms[l].Force.x+=term.x;
                  Adsorbates[CurrentSystem][k].Atoms[l].Force.y+=term.y;
                  Adsorbates[CurrentSystem][k].Atoms[l].Force.z+=term.z;

                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

                  term.x=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
                  term.y=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
                  term.z=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

                  Adsorbates[CurrentSystem][i].Atoms[j].Force.x+=term.x;
                  Adsorbates[CurrentSystem][i].Atoms[j].Force.y+=term.y;
                  Adsorbates[CurrentSystem][i].Atoms[j].Force.z+=term.z;

                  Adsorbates[CurrentSystem][k].Atoms[l].Force.x-=term.x;
                  Adsorbates[CurrentSystem][k].Atoms[l].Force.y-=term.y;
                  Adsorbates[CurrentSystem][k].Atoms[l].Force.z-=term.z;
                }
              }
            }
          }
        }

        for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
          {
            TypeB=Cations[CurrentSystem][k].Atoms[l].Type;
            if(PseudoAtoms[TypeB].HasCharges)
            {
              posB=Cations[CurrentSystem][k].Atoms[l].Position;
              dipoleB=Cations[CurrentSystem][k].Atoms[l].InducedDipole;
              ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;

              dr.x=posA.x-posB.x;
              dr.y=posA.y-posB.y;
              dr.z=posA.z-posB.z;
              dr=ApplyBoundaryCondition(dr);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

              if(rr<CutOffChargeChargeSquared)
              {
                r=sqrt(rr);
                switch(ChargeMethod)
                {
                  case NONE:
                    Bt0=Bt1=Bt2=0.0;
                    break;
                  case TRUNCATED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case SHIFTED_COULOMB:
                    Bt0=1.0/(r);
                    Bt1=1.0/(r*rr);
                    Bt2=3.0/(r*rr*rr);
                    break;
                  case EWALD:
                  default:
                    Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                    Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        erfc(Alpha[CurrentSystem]*r)/(rr*r);
                    Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                        4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                        3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                    break;
                }

                cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

                term.x=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                term.y=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                term.z=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=term.x;
                Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=term.y;
                Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=term.z;

                Cations[CurrentSystem][k].Atoms[l].Force.x+=term.x;
                Cations[CurrentSystem][k].Atoms[l].Force.y+=term.y;
                Cations[CurrentSystem][k].Atoms[l].Force.z+=term.z;

                cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

                term.x=-COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
                term.y=-COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
                term.z=-COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

                Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=term.x;
                Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=term.y;
                Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=term.z;

                Cations[CurrentSystem][k].Atoms[l].Force.x+=term.x;
                Cations[CurrentSystem][k].Atoms[l].Force.y+=term.y;
                Cations[CurrentSystem][k].Atoms[l].Force.z+=term.z;
              }
            }
          }
        }
      }
    }
  }

  if(!OmitCationCationPolarization)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        posA=Cations[CurrentSystem][i].Atoms[j].Position;
        dipoleA=Cations[CurrentSystem][i].Atoms[j].InducedDipole;
        TypeA=Cations[CurrentSystem][i].Atoms[j].Type;

        if(PseudoAtoms[TypeA].HasCharges)
        {
          posA=Cations[CurrentSystem][i].Atoms[j].Position;
          ChargeA=Cations[CurrentSystem][i].Atoms[j].Charge;

          for(k=i+1;k<NumberOfCationMolecules[CurrentSystem];k++)
          {
            for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
            {
              TypeB=Cations[CurrentSystem][k].Atoms[l].Type;
              if(PseudoAtoms[TypeB].HasCharges)
              {
                posB=Cations[CurrentSystem][k].Atoms[l].Position;
                dipoleB=Cations[CurrentSystem][k].Atoms[l].InducedDipole;
                ChargeB=Cations[CurrentSystem][k].Atoms[l].Charge;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                if(rr<CutOffChargeChargeSquared)
                {
                  r=sqrt(rr);
                  switch(ChargeMethod)
                  {
                    case NONE:
                      Bt0=Bt1=Bt2=0.0;
                      break;
                    case TRUNCATED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case SHIFTED_COULOMB:
                      Bt0=1.0/(r);
                      Bt1=1.0/(r*rr);
                      Bt2=3.0/(r*rr*rr);
                      break;
                    case EWALD:
                    default:
                      Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          erfc(Alpha[CurrentSystem]*r)/(rr*r);
                      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                          3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                      break;
                  }

                  cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

                  term.x=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
                  term.y=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
                  term.z=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

                  Cations[CurrentSystem][i].Atoms[j].Force.x-=term.x;
                  Cations[CurrentSystem][i].Atoms[j].Force.y-=term.y;
                  Cations[CurrentSystem][i].Atoms[j].Force.z-=term.z;

                  Cations[CurrentSystem][k].Atoms[l].Force.x+=term.x;
                  Cations[CurrentSystem][k].Atoms[l].Force.y+=term.y;
                  Cations[CurrentSystem][k].Atoms[l].Force.z+=term.z;
 
                  cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

                  term.x=-COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
                  term.y=-COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
                  term.z=-COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

                  Cations[CurrentSystem][i].Atoms[j].Force.x-=term.x;
                  Cations[CurrentSystem][i].Atoms[j].Force.y-=term.y;
                  Cations[CurrentSystem][i].Atoms[j].Force.z-=term.z;
 
                  Cations[CurrentSystem][k].Atoms[l].Force.x+=term.x;
                  Cations[CurrentSystem][k].Atoms[l].Force.y+=term.y;
                  Cations[CurrentSystem][k].Atoms[l].Force.z+=term.z;
                }
              }
            }
          }
        }
      }
    }
  }
  
  return 0;
}

int CalculateInterInducedDipoleInducedDipoleForce(void)
{

  int i,j,k,l;
  VECTOR posA,posB,dr;
  REAL r,rr,cosA,cosB,cosAB;
  VECTOR dipoleA,dipoleB,term;
  REAL Bt0,Bt1,Bt2,Bt3;

  if(ChargeMethod==NONE) return 0;

  if(!BackPolarization) return 0;

  Bt0=Bt1=Bt2=Bt3=0.0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      posA=Adsorbates[CurrentSystem][i].Atoms[j].Position;
      dipoleA=Adsorbates[CurrentSystem][i].Atoms[j].InducedDipole;

      if(!OmitAdsorbateAdsorbatePolarization)
      {
        for(k=i+1;k<NumberOfAdsorbateMolecules[CurrentSystem];k++)
        {
          for(l=0;l<Adsorbates[CurrentSystem][k].NumberOfAtoms;l++)
          {
            posB=Adsorbates[CurrentSystem][k].Atoms[l].Position;
            dipoleB=Adsorbates[CurrentSystem][k].Atoms[l].InducedDipole;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=Bt3=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case EWALD:
                default:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                    20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                break;
              }

              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

              //UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);

              term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
              term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
              term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

              Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=term.x;
              Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=term.y;
              Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=term.z;

              Adsorbates[CurrentSystem][k].Atoms[l].Force.x+=term.x;
              Adsorbates[CurrentSystem][k].Atoms[l].Force.y+=term.y;
              Adsorbates[CurrentSystem][k].Atoms[l].Force.z+=term.z;
            }
          }
        }
      }

      for(k=0;k<NumberOfCationMolecules[CurrentSystem];k++)
      {
        for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
        {
          posB=Cations[CurrentSystem][k].Atoms[l].Position;
          dipoleB=Cations[CurrentSystem][k].Atoms[l].InducedDipole;

          dr.x=posA.x-posB.x;
          dr.y=posA.y-posB.y;
          dr.z=posA.z-posB.z;
          dr=ApplyBoundaryCondition(dr);
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

          if(rr<CutOffChargeChargeSquared)
          {
            r=sqrt(rr);
            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=Bt3=0.0;  
                break;
              case TRUNCATED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                Bt3=15.0/(r*rr*rr*rr);
                break;
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                Bt3=15.0/(r*rr*rr*rr);
                break;
              case EWALD:
              default:
                Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    erfc(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                    20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                break;
            }

            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

            term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
            term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
            term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

            Adsorbates[CurrentSystem][i].Atoms[j].Force.x-=term.x;
            Adsorbates[CurrentSystem][i].Atoms[j].Force.y-=term.y;
            Adsorbates[CurrentSystem][i].Atoms[j].Force.z-=term.z;

            Cations[CurrentSystem][k].Atoms[l].Force.x+=term.x;
            Cations[CurrentSystem][k].Atoms[l].Force.y+=term.y;
            Cations[CurrentSystem][k].Atoms[l].Force.z+=term.z;

          }
        }
      }
    }
  }

  if(!OmitCationCationPolarization)
  {
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        posA=Cations[CurrentSystem][i].Atoms[j].Position;
        dipoleA=Cations[CurrentSystem][i].Atoms[j].InducedDipole;

        for(k=i+1;k<NumberOfCationMolecules[CurrentSystem];k++)
        {
          for(l=0;l<Cations[CurrentSystem][k].NumberOfAtoms;l++)
          {
            posB=Cations[CurrentSystem][k].Atoms[l].Position;
            dipoleB=Cations[CurrentSystem][k].Atoms[l].InducedDipole;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyBoundaryCondition(dr);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            if(rr<CutOffChargeChargeSquared)
            {
              r=sqrt(rr);
              switch(ChargeMethod)
              {
                case NONE:
                  Bt0=Bt1=Bt2=Bt3=0.0;  
                  break;
                case TRUNCATED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case SHIFTED_COULOMB:
                  Bt0=1.0/(r);
                  Bt1=1.0/(r*rr);
                  Bt2=3.0/(r*rr*rr);
                  Bt3=15.0/(r*rr*rr*rr);
                  break;
                case EWALD:
                default:
                  Bt0=erfc(Alpha[CurrentSystem]*r)/r;
                  Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      erfc(Alpha[CurrentSystem]*r)/(rr*r);
                  Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                      4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                      3.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*r);
                Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                    20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    15.0*erfc(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                break;
              }

              cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
              cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
              cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

              term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
              term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
              term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

              Cations[CurrentSystem][i].Atoms[j].Force.x-=term.x;
              Cations[CurrentSystem][i].Atoms[j].Force.y-=term.y;
              Cations[CurrentSystem][i].Atoms[j].Force.z-=term.z;

              Cations[CurrentSystem][k].Atoms[l].Force.x+=term.x;
              Cations[CurrentSystem][k].Atoms[l].Force.y+=term.y;
              Cations[CurrentSystem][k].Atoms[l].Force.z+=term.z;
            }
          }
        }
      }
    }
  }
  return 0;
}

int CalculateTotalInterReplicaVDWForce(void)
{
  int i,j,k,l;
  int typeA,typeB,TypeMolA,TypeMolB;
  REAL rr;
  REAL energy,force_factor;
  VECTOR posA,posB,dr,f;
  REAL ReductionA,ReductionB;
  int ConnectedAtomA,ConnectedAtomB;
  int ncell,k1,k2,k3,start;

  UAdsorbateAdsorbateVDW[CurrentSystem]=0.0;
  UAdsorbateCationVDW[CurrentSystem]=0.0;
  UCationCationVDW[CurrentSystem]=0.0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeMolA=Adsorbates[CurrentSystem][i].Type;
    for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[k].AnisotropicPosition;

      ReductionA=1.0;
      ConnectedAtomA=-1;

      // loop over adsorbant molecules
      if(!OmitAdsorbateAdsorbateVDWInteractions)
      {
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if(ncell==0) start=i+1;
              else start=0;
              for(j=start;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
              {
                TypeMolB=Adsorbates[CurrentSystem][j].Type;
                for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
                {
                  typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
                  posB=Adsorbates[CurrentSystem][j].Atoms[l].AnisotropicPosition;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  ReductionB=1.0;
                  ConnectedAtomB=-1;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffVDWSquared)
                  {
                    PotentialGradient(typeA,typeB,rr,&energy,&force_factor);

                    // energy
                    if(ncell==0)
                      UAdsorbateAdsorbateVDW[CurrentSystem]+=energy;
                    else
                      UAdsorbateAdsorbateVDW[CurrentSystem]+=0.5*energy;

                    // forces
                    f.x=force_factor*dr.x;
                    f.y=force_factor*dr.y;
                    f.z=force_factor*dr.z;

                    Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=ReductionA*f.x;
                    Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=ReductionA*f.y;
                    Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=ReductionA*f.z;

                    if(ConnectedAtomA>=0)
                    {
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.x-=(1.0-ReductionA)*f.x;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.y-=(1.0-ReductionA)*f.y;
                      Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.z-=(1.0-ReductionA)*f.z;
                    }

                    if(ncell==0)
                    {
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.x+=ReductionB*f.x;
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.y+=ReductionB*f.y;
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.z+=ReductionB*f.z;

                      if(ConnectedAtomB>=0)
                      {
                        Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB].Force.x+=(1.0-ReductionB)*f.x;
                        Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB].Force.y+=(1.0-ReductionB)*f.y;
                        Adsorbates[CurrentSystem][j].Atoms[ConnectedAtomB].Force.z+=(1.0-ReductionB)*f.z;
                      }

                      StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
                    }
                    else
                    {
                      StrainDerivativeTensor[CurrentSystem].ax+=0.5*f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx+=0.5*f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx+=0.5*f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay+=0.5*f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by+=0.5*f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy+=0.5*f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az+=0.5*f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz+=0.5*f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz+=0.5*f.z*dr.z;
                    }
                  }
                }
              }
              ncell++;
            }
      }

      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
            {
              TypeMolB=Cations[CurrentSystem][j].Type;
              for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
              {
                posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
                typeB=Cations[CurrentSystem][j].Atoms[l].Type;

                posB.x+=ReplicaShift[ncell].x;
                posB.y+=ReplicaShift[ncell].y;
                posB.z+=ReplicaShift[ncell].z;
 
                ReductionB=1.0;
                ConnectedAtomB=-1;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                {
                  PotentialGradient(typeA,typeB,rr,&energy,&force_factor);

                  // energy
                  UAdsorbateCationVDW[CurrentSystem]+=energy;

                  // forces
                  f.x=force_factor*dr.x;
                  f.y=force_factor*dr.y;
                  f.z=force_factor*dr.z;

                  Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=ReductionA*f.x;
                  Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=ReductionA*f.y;
                  Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=ReductionA*f.z;

                  if(ConnectedAtomA>=0)
                  {
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.x-=(1.0-ReductionA)*f.x;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.y-=(1.0-ReductionA)*f.y;
                    Adsorbates[CurrentSystem][i].Atoms[ConnectedAtomA].Force.z-=(1.0-ReductionA)*f.z;
                  }

                  Cations[CurrentSystem][j].Atoms[l].Force.x+=ReductionB*f.x;
                  Cations[CurrentSystem][j].Atoms[l].Force.y+=ReductionB*f.y;
                  Cations[CurrentSystem][j].Atoms[l].Force.z+=ReductionB*f.z;

                  if(ConnectedAtomB>=0)
                  {
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB].Force.x+=(1.0-ReductionB)*f.x;
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB].Force.y+=(1.0-ReductionB)*f.y;
                    Cations[CurrentSystem][j].Atoms[ConnectedAtomB].Force.z+=(1.0-ReductionB)*f.z;
                  }

                  // stress-tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
                }
              }
            }
            ncell++;
          }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    TypeMolA=Cations[CurrentSystem][i].Type;
    for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[k].Type;
      posA=Cations[CurrentSystem][i].Atoms[k].AnisotropicPosition;

      ReductionA=1.0;
      ConnectedAtomA=-1;

      // loop over cation molecules
      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            if(ncell==0) start=i+1;
            else start=0;
            for(j=start;j<NumberOfCationMolecules[CurrentSystem];j++)
            {
              TypeMolB=Cations[CurrentSystem][j].Type;
              for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
              {
                posB=Cations[CurrentSystem][j].Atoms[l].AnisotropicPosition;
                typeB=Cations[CurrentSystem][j].Atoms[l].Type;

                posB.x+=ReplicaShift[ncell].x;
                posB.y+=ReplicaShift[ncell].y;
                posB.z+=ReplicaShift[ncell].z;

                ReductionB=1.0;
                ConnectedAtomB=-1;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                {
                  typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                  PotentialGradient(typeA,typeB,rr,&energy,&force_factor);

                  // energy
                  if(ncell==0)
                    UCationCationVDW[CurrentSystem]+=energy;
                  else
                    UCationCationVDW[CurrentSystem]+=0.5*energy;

                  // forces
                  f.x=force_factor*dr.x;
                  f.y=force_factor*dr.y;
                  f.z=force_factor*dr.z;

                  Cations[CurrentSystem][i].Atoms[k].Force.x-=ReductionA*f.x;
                  Cations[CurrentSystem][i].Atoms[k].Force.y-=ReductionA*f.y;
                  Cations[CurrentSystem][i].Atoms[k].Force.z-=ReductionA*f.z;

                  if(ConnectedAtomA>=0)
                  {
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA].Force.x-=(1.0-ReductionA)*f.x;
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA].Force.y-=(1.0-ReductionA)*f.y;
                    Cations[CurrentSystem][i].Atoms[ConnectedAtomA].Force.z-=(1.0-ReductionA)*f.z;
                  }

                  if(ncell==0)
                  {
                    Cations[CurrentSystem][j].Atoms[l].Force.x+=ReductionB*f.x;
                    Cations[CurrentSystem][j].Atoms[l].Force.y+=ReductionB*f.y;
                    Cations[CurrentSystem][j].Atoms[l].Force.z+=ReductionB*f.z;

                    if(ConnectedAtomB>=0)
                    {
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB].Force.x+=(1.0-ReductionB)*f.x;
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB].Force.y+=(1.0-ReductionB)*f.y;
                      Cations[CurrentSystem][j].Atoms[ConnectedAtomB].Force.z+=(1.0-ReductionB)*f.z;
                    }

                    // stress-tensor
                    StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                    StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                    StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                    StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                    StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                    StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                    StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                    StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                    StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
                  }
                  else
                  {
                    // stress-tensor
                    StrainDerivativeTensor[CurrentSystem].ax+=0.5*f.x*dr.x;
                    StrainDerivativeTensor[CurrentSystem].bx+=0.5*f.y*dr.x;
                    StrainDerivativeTensor[CurrentSystem].cx+=0.5*f.z*dr.x;

                    StrainDerivativeTensor[CurrentSystem].ay+=0.5*f.x*dr.y;
                    StrainDerivativeTensor[CurrentSystem].by+=0.5*f.y*dr.y;
                    StrainDerivativeTensor[CurrentSystem].cy+=0.5*f.z*dr.y;

                    StrainDerivativeTensor[CurrentSystem].az+=0.5*f.x*dr.z;
                    StrainDerivativeTensor[CurrentSystem].bz+=0.5*f.y*dr.z;
                    StrainDerivativeTensor[CurrentSystem].cz+=0.5*f.z*dr.z;
                  }
                }
              }
            }
            ncell++;
          }
    }
  }
  return 0;
}

int CalculateTotalInterReplicaChargeChargeCoulombForce(void)
{
  int i,j,k,l;
  int typeA,typeB;
  REAL r,rr;
  REAL chargeA,chargeB;
  REAL force_factor,energy;
  VECTOR posA,posB,dr,f;
  int ncell,k1,k2,k3,start;

  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeReal[CurrentSystem]=0.0;
  UCationCationChargeChargeReal[CurrentSystem]=0.0;

  if(ChargeMethod==NONE) return 0;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[k].Position;
      chargeA=Adsorbates[CurrentSystem][i].Atoms[k].Charge;

      // loop over adsorbant molecules
      if(!OmitAdsorbateAdsorbateCoulombInteractions)
      {
        ncell=0;
        for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
          for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
            for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
            {
              if(ncell==0) start=i+1;
              else start=0;
              for(j=start;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
              {
                for(l=0;l<Adsorbates[CurrentSystem][j].NumberOfAtoms;l++)
                {
                  posB=Adsorbates[CurrentSystem][j].Atoms[l].Position;

                  posB.x+=ReplicaShift[ncell].x;
                  posB.y+=ReplicaShift[ncell].y;
                  posB.z+=ReplicaShift[ncell].z;

                  dr.x=posA.x-posB.x;
                  dr.y=posA.y-posB.y;
                  dr.z=posA.z-posB.z;
                  dr=ApplyReplicaBoundaryCondition(dr);
                  rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

                  if(rr<CutOffChargeChargeSquared)
                  {
                    r=sqrt(rr);
                    typeB=Adsorbates[CurrentSystem][j].Atoms[l].Type;
                    chargeB=Adsorbates[CurrentSystem][j].Atoms[l].Charge;

                    switch(ChargeMethod)
                    {
                      case NONE:
                        energy=0.0;
                        force_factor=0.0;
                        break;
                      case TRUNCATED_COULOMB:
                        energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                        force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                        break;
                      case SHIFTED_COULOMB:
                        energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                        force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                        break;
                      case WOLFS_METHOD_DAMPED_FG:
                        energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                 -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                                  (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                  (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                        force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/rr+
                                       2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                       -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                        2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                        break;
                      case EWALD:
                        energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                               (erfc(Alpha[CurrentSystem]*r)/r);
                        force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                               (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                               (r*rr);
                        break;
                      default:
                        printf("Unknown charge-method in 'CalculateTotalInterReplicaChargeChargeCoulombForce'\n");
                        exit(0);
                        break;
                    }

                    if(ncell==0)
                      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=energy;
                    else
                      UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=0.5*energy;

                    // forces
                    f.x=force_factor*dr.x;
                    f.y=force_factor*dr.y;
                    f.z=force_factor*dr.z;

                    Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=f.x;
                    Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=f.y;
                    Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=f.z;

                    if(ncell==0)
                    {
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.x+=f.x;
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.y+=f.y;
                      Adsorbates[CurrentSystem][j].Atoms[l].Force.z+=f.z;

                      // the strain derivative
                      StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
                    }
                    else
                    {

                      // the strain derivative
                      StrainDerivativeTensor[CurrentSystem].ax+=0.5*f.x*dr.x;
                      StrainDerivativeTensor[CurrentSystem].bx+=0.5*f.y*dr.x;
                      StrainDerivativeTensor[CurrentSystem].cx+=0.5*f.z*dr.x;

                      StrainDerivativeTensor[CurrentSystem].ay+=0.5*f.x*dr.y;
                      StrainDerivativeTensor[CurrentSystem].by+=0.5*f.y*dr.y;
                      StrainDerivativeTensor[CurrentSystem].cy+=0.5*f.z*dr.y;

                      StrainDerivativeTensor[CurrentSystem].az+=0.5*f.x*dr.z;
                      StrainDerivativeTensor[CurrentSystem].bz+=0.5*f.y*dr.z;
                      StrainDerivativeTensor[CurrentSystem].cz+=0.5*f.z*dr.z;
                    }
                  }
                }
              }
              ncell++;
            }
      }

      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
            {
              for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
              {
                posB=Cations[CurrentSystem][j].Atoms[l].Position;

                posB.x+=ReplicaShift[ncell].x;
                posB.y+=ReplicaShift[ncell].y;
                posB.z+=ReplicaShift[ncell].z;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                {
                  typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                  chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;

                  r=sqrt(rr);

                  switch(ChargeMethod)
                  {
                    case NONE:
                      energy=0.0;
                      force_factor=0.0;
                      break;
                    case TRUNCATED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                      break;
                    case SHIFTED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                      break;
                    case WOLFS_METHOD_DAMPED_FG:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                               -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                                (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/rr+
                                    2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                    -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                    2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                      break;
                    case EWALD:
                      energy=COULOMBIC_CONVERSION_FACTOR*erfc(Alpha[CurrentSystem]*r)*
                                                       chargeA*chargeB/r;

                      force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                           (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                           (r*rr);
                      break;
                    default:
                      printf("Unknown charge-method in 'CalculateTotalInterReplicaChargeChargeCoulombForce'\n");
                      exit(0);
                      break;
                  }

                  UAdsorbateCationChargeChargeReal[CurrentSystem]+=energy;

                  // forces
                  f.x=force_factor*dr.x;
                  f.y=force_factor*dr.y;
                  f.z=force_factor*dr.z;

                  Adsorbates[CurrentSystem][i].Atoms[k].Force.x-=f.x;
                  Adsorbates[CurrentSystem][i].Atoms[k].Force.y-=f.y;
                  Adsorbates[CurrentSystem][i].Atoms[k].Force.z-=f.z;

                  Cations[CurrentSystem][j].Atoms[l].Force.x+=f.x;
                  Cations[CurrentSystem][j].Atoms[l].Force.y+=f.y;
                  Cations[CurrentSystem][j].Atoms[l].Force.z+=f.z;

                  // stress-tensor
                  StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                  StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                  StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                  StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                  StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                  StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;
 
                  StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                  StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                  StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
                }
              }
            }
            ncell++;
        }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[k].Type;
      posA=Cations[CurrentSystem][i].Atoms[k].Position;
      chargeA=Cations[CurrentSystem][i].Atoms[k].Charge;

      // loop over cation molecules
      ncell=0;
      for(k1=0;k1<NumberOfReplicaCells[CurrentSystem].x;k1++)
        for(k2=0;k2<NumberOfReplicaCells[CurrentSystem].y;k2++)
          for(k3=0;k3<NumberOfReplicaCells[CurrentSystem].z;k3++)
          {
            if(ncell==0) start=i+1;
            else start=0;
            for(j=start;j<NumberOfCationMolecules[CurrentSystem];j++)
            {
              for(l=0;l<Cations[CurrentSystem][j].NumberOfAtoms;l++)
              {
                posB=Cations[CurrentSystem][j].Atoms[l].Position;

                posB.x+=ReplicaShift[ncell].x;
                posB.y+=ReplicaShift[ncell].y;
                posB.z+=ReplicaShift[ncell].z;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                dr=ApplyReplicaBoundaryCondition(dr);
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                if(rr<CutOffVDWSquared)
                {
                  typeB=Cations[CurrentSystem][j].Atoms[l].Type;
                  chargeB=Cations[CurrentSystem][j].Atoms[l].Charge;
 
                  r=sqrt(rr);

                  switch(ChargeMethod)
                  {
                    case NONE:
                      energy=0.0;
                      force_factor=0.0;
                      break;
                    case TRUNCATED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/r;
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                      break;
                    case SHIFTED_COULOMB:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(1.0/r-InverseCutOffChargeCharge);
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB/(rr*r);
                      break;
                    case WOLFS_METHOD_DAMPED_FG:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/r
                                -erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*InverseCutOffChargeCharge+
                                (r-CutOffChargeCharge)*(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                (2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge)));
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*(erfc(Alpha[CurrentSystem]*r)/rr+
                                    2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)*M_1_SQRTPI/r
                                   -(erfc(Alpha[CurrentSystem]*CutOffChargeCharge)*SQR(InverseCutOffChargeCharge)+
                                    2.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*CutOffChargeChargeSquared)*M_1_SQRTPI*InverseCutOffChargeCharge))/r;
                      break;
                    case EWALD:
                      energy=COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                             (erfc(Alpha[CurrentSystem]*r)/r);
                      force_factor=-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
                              (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
                              (r*rr);
                      break;
                    default:
                      printf("Unknown charge-method in 'CalculateTotalInterReplicaChargeChargeCoulombForce'\n");
                      exit(0);
                      break;
                  }

                  if(ncell==0)
                    UCationCationChargeChargeReal[CurrentSystem]+=energy;
                  else
                    UCationCationChargeChargeReal[CurrentSystem]+=0.5*energy;

                  // forces
                  f.x=force_factor*dr.x;
                  f.y=force_factor*dr.y;
                  f.z=force_factor*dr.z;

                  Cations[CurrentSystem][i].Atoms[k].Force.x-=f.x;
                  Cations[CurrentSystem][i].Atoms[k].Force.y-=f.y;
                  Cations[CurrentSystem][i].Atoms[k].Force.z-=f.z;

                  if(ncell==0)
                  {
                    Cations[CurrentSystem][j].Atoms[l].Force.x+=f.x;
                    Cations[CurrentSystem][j].Atoms[l].Force.y+=f.y;
                    Cations[CurrentSystem][j].Atoms[l].Force.z+=f.z;

                    // stress-tensor
                    StrainDerivativeTensor[CurrentSystem].ax+=f.x*dr.x;
                    StrainDerivativeTensor[CurrentSystem].bx+=f.y*dr.x;
                    StrainDerivativeTensor[CurrentSystem].cx+=f.z*dr.x;

                    StrainDerivativeTensor[CurrentSystem].ay+=f.x*dr.y;
                    StrainDerivativeTensor[CurrentSystem].by+=f.y*dr.y;
                    StrainDerivativeTensor[CurrentSystem].cy+=f.z*dr.y;

                    StrainDerivativeTensor[CurrentSystem].az+=f.x*dr.z;
                    StrainDerivativeTensor[CurrentSystem].bz+=f.y*dr.z;
                    StrainDerivativeTensor[CurrentSystem].cz+=f.z*dr.z;
                  }
                  else
                  {
                    // stress-tensor
                    StrainDerivativeTensor[CurrentSystem].ax+=0.5*f.x*dr.x;
                    StrainDerivativeTensor[CurrentSystem].bx+=0.5*f.y*dr.x;
                    StrainDerivativeTensor[CurrentSystem].cx+=0.5*f.z*dr.x;

                    StrainDerivativeTensor[CurrentSystem].ay+=0.5*f.x*dr.y;
                    StrainDerivativeTensor[CurrentSystem].by+=0.5*f.y*dr.y;
                    StrainDerivativeTensor[CurrentSystem].cy+=0.5*f.z*dr.y;

                    StrainDerivativeTensor[CurrentSystem].az+=0.5*f.x*dr.z;
                    StrainDerivativeTensor[CurrentSystem].bz+=0.5*f.y*dr.z;
                    StrainDerivativeTensor[CurrentSystem].cz+=0.5*f.z*dr.z;
                  }
                }
              }
            }
            ncell++;
          }
    }
  }
  return 0;
}
