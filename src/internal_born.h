/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'internal_born.h' is part of RASPA.

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

#ifndef INTERNAL_BORN_H
#define INTERNAL_BORN_H

void CalculateAdsorbateBondBornTerm(void);
void CalculateAdsorbateUreyBradleyBornTerm(void);
void CalculateAdsorbateBendBornTerm(void);
void CalculateAdsorbateTorsionBornTerm(void);
void CalculateAdsorbateImproperTorsionBornTerm(void);
void CalculateAdsorbateIntraVDWBornTerm(void);
void CalculateAdsorbateIntraCoulombBornTerm(void);

void CalculateCationBondBornTerm(void);
void CalculateCationUreyBradleyBornTerm(void);
void CalculateCationBendBornTerm(void);
void CalculateCationTorsionBornTerm(void);
void CalculateCationImproperTorsionBornTerm(void);
void CalculateCationIntraVDWBornTerm(void);
void CalculateCationIntraCoulombBornTerm(void);

#endif
