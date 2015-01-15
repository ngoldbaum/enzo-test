/***********************************************************************
/
/ Accreting Particle
/
************************************************************************/

#ifdef USE_MPI
#include "communicators.h"
#endif 

#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CommunicationUtilities.h"
#include "communication.h"
#include "phys_constants.h"
#include "FofLib.h"

/* Prototypes */

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

class ActiveParticleType_AccretingParticle : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_AccretingParticle(void) : ActiveParticleType() {
    AccretionRate = 0;
  };
  ActiveParticleType_AccretingParticle(ActiveParticleType_AccretingParticle* part) :
    ActiveParticleType(static_cast<ActiveParticleType*>(part)) {
    AccretionRate = part->AccretionRate;
  };
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  template <class active_particle_class>
    static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				 int ThisLevel, bool CallEvolvePhotons,
				 int TotalStarParticleCountPrevious[],
				 int AccretingParticleID);
  template <class active_particle_class>
    static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalStarParticleCountPrevious[],
				int AccretingParticleID);
  template <class active_particle_class>
    static int DepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int GalaxyParticleID) {return SUCCESS; };
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, int TopGridDims[], int ActiveParticleID);
  static int InitializeParticleType();
  ENABLED_PARTICLE_ID_ACCESSOR
  bool IsARadiationSource(FLOAT Time);
  
  // sink helper routines

  template <class active_particle_class>
    static active_particle_class** MergeAccretingParticles(int *nParticles, ActiveParticleType** ParticleList, 
							   FLOAT LinkingLength,int *ngroups, 
							   LevelHierarchyEntry *LevelArray[], int ThisLevel);

  static int Accrete(int nParticles, ActiveParticleType** ParticleList,
		     int AccretionRadius, FLOAT dx, 
		     LevelHierarchyEntry *LevelArray[], int ThisLevel);

  static int AccreteOntoAccretingParticle(grid* AccretionZone, 
					  ActiveParticleType_AccretingParticle* ThisParticle,
					  FLOAT AccretionRadius);

  static float OverflowFactor;
  static int AccretionRadius;   // in units of CellWidth on the maximum refinement level
  static int LinkingLength;     // Should be equal to AccretionRadius
  static int RadiationParticle;
  static double LuminosityPerSolarMass;
  static int RadiationSEDNumberOfBins;
  static float* RadiationEnergyBins;
  static float* RadiationSED;
  static float RadiationLifetime;
  float AccretionRate;
  static std::vector<ParticleAttributeHandler *> AttributeHandlers;
};

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

template <class active_particle_class>
active_particle_class** ActiveParticleType_AccretingParticle::MergeAccretingParticles
(int *nParticles, ActiveParticleType** ParticleList, FLOAT LinkingLength, 
 int *ngroups, LevelHierarchyEntry *LevelArray[], int ThisLevel)
{
  int i,j;
  int dim;
  int GroupNumberAssignment[*nParticles];
  FLOAT* tempPos = NULL;
  int *groupsize = NULL;
  int **grouplist = NULL;
  active_particle_class **MergedParticles = NULL;

  HierarchyEntry** LevelGrids = NULL;

  int NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &LevelGrids);

  /* Construct list of sink particle positions to pass to Foflist */
  FLOAT ParticleCoordinates[3*(*nParticles)];
  
  for (i=0; i<(*nParticles); i++) {
    tempPos = ParticleList[i]->ReturnPosition();
    for (dim=0; dim<3; dim++)
      ParticleCoordinates[3*i+dim] = tempPos[dim];
  }

  /* Find mergeable groups using an FOF search */

  *ngroups = FofList((*nParticles), ParticleCoordinates, LinkingLength, GroupNumberAssignment, &groupsize, &grouplist);
  
  MergedParticles = new active_particle_class*[*ngroups]();
  
  /* Merge the mergeable groups */

  for (i=0; i<*ngroups; i++) {
    MergedParticles[i] = static_cast<active_particle_class*>(ParticleList[grouplist[i][0]]);
    if (groupsize[i] != 1) {
      for (j=1; j<groupsize[i]; j++) {
	MergedParticles[i]->Merge(static_cast<active_particle_class*>(ParticleList[grouplist[i][j]]));
	if (ParticleList[grouplist[i][j]]->DisableParticle(LevelArray,MergedParticles[i]->
							   ReturnCurrentGrid()->ReturnProcessorNumber()) == FAIL)
	  ENZO_FAIL("MergeAccretingParticles: DisableParticle failed!\n");
	if (NumberOfProcessors > 1) {
	  delete ParticleList[grouplist[i][j]];
	  ParticleList[grouplist[i][j]] = NULL;
	}
      }
    }
  }

  delete [] groupsize;
  groupsize = NULL;
  for (i=0; i<*ngroups; i++)
    delete [] grouplist[i];
  delete [] grouplist;
  grouplist = NULL;

  /* Loop over the grids and check if any of the merged particles have
     moved. If so, disable the particle on the current grid and assign
     it to the new grid*/

  int NewGrid = -1;

  for (i = 0; i < *ngroups; i++) {
    if (MergedParticles[i]->ReturnCurrentGrid()->PointInGrid(MergedParticles[i]->ReturnPosition()) == false) {
      // Find the grid to transfer to 
      for (j = 0; j < NumberOfGrids; j++) {
	if (LevelGrids[j]->GridData->PointInGrid(MergedParticles[i]->ReturnPosition())) {
	  NewGrid = j;
	  break;
	}
      }
      if (NewGrid == -1)
	ENZO_FAIL("Cannot assign particle to grid after merging!\n");
      int OldProc = MergedParticles[i]->CurrentGrid->ReturnProcessorNumber();
      active_particle_class *temp = new active_particle_class(MergedParticles[i]);
      MergedParticles[i]->DisableParticle(LevelArray,LevelGrids[NewGrid]->GridData->ReturnProcessorNumber()); 
      if (LevelGrids[NewGrid]->GridData->AddActiveParticle(static_cast<ActiveParticleType*>(temp)) == FAIL)
      	ENZO_FAIL("Active particle grid assignment failed!\n");
      if (MyProcessorNumber == OldProc) {
	delete MergedParticles[i];
	MergedParticles[i] = new active_particle_class(temp);
      }
      else if (MyProcessorNumber != temp->CurrentGrid->ReturnProcessorNumber())
	delete temp;
      MergedParticles[i]->AssignCurrentGrid(LevelGrids[NewGrid]->GridData);
    }
  }

  delete [] LevelGrids;

  *nParticles = *ngroups;

  return MergedParticles;
}


int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
 
int AssignActiveParticlesToGrids(ActiveParticleType** ParticleList, int nParticles, 
				 LevelHierarchyEntry *LevelArray[]); 

template <class active_particle_class>
int ActiveParticleType_AccretingParticle::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
							   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
							   int ThisLevel, int TotalStarParticleCountPrevious[],
							   int AccretingParticleID)
{

  /* Accreting particles live on the maximum refinement level.  If we are on a lower level, this does not concern us */

  if (ThisLevel == MaximumRefinementLevel)
    {

      /* Generate a list of all sink particles in the simulation box */
      int i,nParticles,NumberOfMergedParticles;
      ActiveParticleType** ParticleList = NULL;

      ParticleList = ActiveParticleFindAll(LevelArray, &nParticles, AccretingParticleID);

      /* Return if there are no accreting particles */
      
      if (nParticles == 0)
	return SUCCESS;

      /* Calculate CellWidth on maximum refinement level */

      // This assumes a cubic box and may not work for simulations with MinimumMassForRefinementLevelExponent
      FLOAT dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
	(MetaData->TopGridDims[0]*POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));

      /* Do Merging   */

      active_particle_class **MergedParticles = NULL;
      
      /* Generate new merged list of sink particles */
      
      MergedParticles = MergeAccretingParticles
	<active_particle_class>(&nParticles, ParticleList, LinkingLength*dx,
				&NumberOfMergedParticles,LevelArray,ThisLevel);

      delete [] ParticleList;
   
      // Do merging twice to catch pathological cases where merging
      // leaves multiple sinks inside the same accretion zone.

      ParticleList = new ActiveParticleType*[NumberOfMergedParticles];

      for (i = 0; i<NumberOfMergedParticles; i++)
	ParticleList[i] = static_cast<ActiveParticleType*>(MergedParticles[i]);

      delete [] MergedParticles;

      MergedParticles = MergeAccretingParticles
	<active_particle_class>(&nParticles, ParticleList, LinkingLength*dx,
				&NumberOfMergedParticles,LevelArray,ThisLevel);

      delete [] ParticleList;

      if (debug)
	printf("Number of particles after merging: %"ISYM"\n",NumberOfMergedParticles);

      /* Assign local particles to grids */
 
      ParticleList = new ActiveParticleType*[NumberOfMergedParticles];

      // need to use a bit of redirection because C++ pointer arrays have
      // trouble with polymorphism
      for (i = 0; i<NumberOfMergedParticles; i++)
	ParticleList[i] = static_cast<ActiveParticleType*>(MergedParticles[i]);

      if (AssignActiveParticlesToGrids(ParticleList,NumberOfMergedParticles, LevelArray) == FAIL)
	return FAIL;

      delete [] ParticleList;

      for (i = 0; i<NumberOfMergedParticles; i++)
	if (MergedParticles[i]->ReturnCurrentGrid()->ReturnProcessorNumber() != MyProcessorNumber) {
	  delete MergedParticles[i];
	  MergedParticles[i] = NULL;
	}

      delete [] MergedParticles;
      
      /* Regenerate the global active particle list */
      
      ParticleList = ActiveParticleFindAll(LevelArray, &nParticles, AccretingParticleID);

      /* Do accretion */
      if (Accrete(nParticles,ParticleList,AccretionRadius,dx,LevelArray,ThisLevel) == FAIL)
	ENZO_FAIL("Accreting Particle accretion failed. \n");
      
      for (i = 0; i < nParticles; i++)
	if (ParticleList[i]->ReturnCurrentGrid()->ReturnProcessorNumber() != MyProcessorNumber) {
	  delete ParticleList[i];
	  ParticleList[i] = NULL;
	}

      // Accrete takes care of its own memory management - no need to free ParticleList.

      delete [] ParticleList;

    }

  return SUCCESS;
}

