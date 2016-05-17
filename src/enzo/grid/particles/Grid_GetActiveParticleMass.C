/***********************************************************************
/
/  GRID CLASS (RETURN ACTIVE PARTICLE MASS)
/
/  written by: John Regan
/  date:       September, 2014
/  modified1:  
/
/
************************************************************************/
#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "fortran.def"
#include "CosmologyParameters.h"

#include "ActiveParticle.h"

void grid::GetActiveParticleMass(float * ActiveParticleMass) 
{
  int i = 0;
  for (i = 0; i < NumberOfActiveParticles; i++) {
    ActiveParticleMass[i] = ActiveParticles[i]->ReturnMass();
  }

  return;
}

void grid::GetActiveParticleFixedInSpace(int * ActiveParticleFixedInSpace) 
{
  int i = 0;
  for (i = 0; i < NumberOfActiveParticles; i++) {
    ActiveParticleFixedInSpace[i] = ActiveParticles[i]->ReturnFixedInSpace();
  }

  return;
}
