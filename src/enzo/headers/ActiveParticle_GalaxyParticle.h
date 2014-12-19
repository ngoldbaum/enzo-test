/***********************************************************************
/
/  GALAXY PARTICLE TYPE
/
/  written by: Stephen Skory
/  date:       August, 2012
/
/  PURPOSE:
/
************************************************************************/

#include "preincludes.h"
#include "hdf5.h"
#include "h5utilities.h"

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
#include "EventHooks.h"
#include "ActiveParticle.h"

class ActiveParticleType_GalaxyParticle : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_GalaxyParticle(void) : ActiveParticleType() {
    Radius = 0;
    initialized = 0;
  };
  ActiveParticleType_GalaxyParticle(ActiveParticleType_GalaxyParticle* part) :
    ActiveParticleType(static_cast<ActiveParticleType*>(part)) {
    Radius = part->Radius;
    initialized = part->initialized;
  };

  // Static members
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  template <class active_particle_class>
    static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				 int ThisLevel, bool CallEvolvePhotons,
				 int TotalStarParticleCountPrevious[],
				 int GalaxyParticleID);
  template <class active_particle_class>
    static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalStarParticleCountPrevious[],
				int GalaxyParticleID);
  template <class active_particle_class>
    static int DepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int GalaxyParticleID);
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, int TopGridDims[], int ActiveParticleID);
  // Galaxy Particle helper routines.
  static int SubtractMassFromGrid(int nParticles,
    ActiveParticleType** ParticleList, LevelHierarchyEntry *LevelArray[],
    FLOAT dx, int ThisLevel);
  static int InitializeParticleType(void);
  static int GalaxyParticleFeedback(int nParticles, ActiveParticleType** ParticleList,
		     FLOAT dx, LevelHierarchyEntry *LevelArray[], int ThisLevel,
		     FLOAT period[3]);
  
  static int GalaxyParticleGravity(int nParticles, ActiveParticleType** ParticleList,
		     FLOAT dx, LevelHierarchyEntry *LevelArray[], int ThisLevel,
		     FLOAT period[3]);
  
  static std::vector<ParticleAttributeHandler *> AttributeHandlers;

  // Need this to make active particle ID work correctly.
  int GetEnabledParticleID(int myid = -1) {
    static int ParticleID = -1;
    if (myid >= 0) {
      if (ParticleID != -1) ENZO_FAIL("Setting Particle ID Twice!");
      ParticleID = myid;
    }
    return ParticleID;
  };
  
  // Galaxy Particle specific stuff.
  float Radius;
  int initialized; // Unused right now...
  
  
};

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int AssignActiveParticlesToGrids(ActiveParticleType** ParticleList, int nParticles, 
				 LevelHierarchyEntry *LevelArray[]); 

template <class active_particle_class>
int ActiveParticleType_GalaxyParticle::BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							 int ThisLevel, bool CallEvolvePhotons,
							 int TotalStarParticleCountPrevious[],
							 int GalaxyParticleID)
{

  return SUCCESS;
}

template <class active_particle_class>
int ActiveParticleType_GalaxyParticle::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							   int ThisLevel, int TotalStarParticleCountPrevious[],
							   int GalaxyParticleID)
{

  /* Galaxy particles live on the maximum refinement level.
     If we are on a lower level, this does not concern us */

  if (ThisLevel == MaximumRefinementLevel)
    {

      /* Generate a list of all galaxy particles in the simulation box */
      int i,nParticles;
      ActiveParticleType** ParticleList = NULL;

      ParticleList = ActiveParticleFindAll(LevelArray, &nParticles, GalaxyParticleID);

      /* Return if there are no galaxy particles */
      
      if (nParticles == 0)
	return SUCCESS;

      /* Calculate CellWidth on maximum refinement level */

      // This assumes a cubic box and may not work for simulations with
      // MinimumMassForRefinementLevelExponent
      FLOAT dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
	(MetaData->TopGridDims[0]*POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));
	  // Find the period of the box.
	  // I'm going to be cheap here because galaxy particles are only ever
	  // meant to be run in 3D.
	  FLOAT period[3];
	  for (int dim = 0; dim < 3; dim++) {
	    period[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
	  }

     /* Apply feedback */
      if (GalaxyParticleFeedback(nParticles, ParticleList,
        dx, LevelArray, ThisLevel, period) == FAIL)
	ENZO_FAIL("Galaxy Particle Feedback failed. \n");

    delete [] ParticleList;

    }

  return SUCCESS;

}

template <class active_particle_class>
int ActiveParticleType_GalaxyParticle::DepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
							   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							   int ThisLevel, int GalaxyParticleID)
{

  if (ThisLevel == MaximumRefinementLevel)
    {

      /* Generate a list of all galaxy particles in the simulation box */
      int i,nParticles;
      ActiveParticleType** ParticleList = NULL;

      ParticleList = ActiveParticleFindAll(LevelArray, &nParticles, GalaxyParticleID);

      /* Return if there are no galaxy particles */
      
      if (nParticles == 0)
	   return SUCCESS;

      /* Calculate CellWidth on maximum refinement level */

      // This assumes a cubic box and may not work for simulations with
      // MinimumMassForRefinementLevelExponent
      FLOAT dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
	(MetaData->TopGridDims[0]*POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));
	  // Find the period of the box.
	  // I'm going to be cheap here because galaxy particles are only ever
	  // meant to be run in 3D.
	  FLOAT period[3];
	  for (int dim = 0; dim < 3; dim++) {
	    period[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
	  }

     /* Apply feedback */
      if (GalaxyParticleGravity(nParticles, ParticleList,
        dx, LevelArray, ThisLevel, period) == FAIL)
	ENZO_FAIL("Galaxy Particle Gravity failed. \n");

    delete [] ParticleList;

    }

  return SUCCESS;

}
