/***********************************************************************

  

  written by: John Regan
  date:       March, 2014
  modified1: John Regan - ported to Enzo-3.0

  PURPOSE: Somewhat idealised treatment of radiation sources. 
           This routine creates massless radiation sources 
           at predefined points in the simulation. The sources
           emit radiation for a user defined time. 
          
 
  RETURNS:
     SUCCESS if the particle is successfully created.
     FAIL otherwise
***********************************************************************/

#include "ActiveParticle_RadiationParticle.h"
#include "phys_constants.h"
#include "CosmologyParameters.h"


#define APDEBUG 0
#define NUMPARAMS 4
#define RadiationSourceLifeTime 1e10   //i.e. infinite

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class RadiationParticleGrid : private grid {
  friend class ActiveParticleType_RadiationParticle;
};

//Initialize static variables
int ActiveParticleType_RadiationParticle::RadiationNumSources          = INT_UNDEFINED;
InitData *ActiveParticleType_RadiationParticle::Root                   = NULL;
int ActiveParticleType_RadiationParticle::RadiationSEDNumberOfBins     = INT_UNDEFINED;
float* ActiveParticleType_RadiationParticle::RadiationEnergyBins       = NULL;
float* ActiveParticleType_RadiationParticle::RadiationSED              = NULL;
float ActiveParticleType_RadiationParticle::RadiationLifetime          = FLOAT_UNDEFINED;
float ActiveParticleType_RadiationParticle::RadiationPhotonsPerSecond  = FLOAT_UNDEFINED;

//Redshift
float CR = 0.0;

int ActiveParticleType_RadiationParticle::InitializeParticleType() {

#ifdef NEW_CONFIG
  
  // Update the parameter config to include the local defaults.  Note
  // that this does not overwrite any values previously specified.
  Param.Update(config_radiation_particle_defaults);

  // Retrieve parameters from Param structure
  Param.GetScalar(RadiationNumSources, 
		  "Physics.ActiveParticles.RadiationParticle.RadiationNumSources");
  Param.GetScalar(RadiationSEDNumberOfBins, 
		  "Physics.ActiveParticles.RadiationParticle.RadiationSEDNumberOfBins");
  Param.GetScalar(RadiationPhotonsPerSecond, 
		 "Physics.ActiveParticles.RadiationParticle.RadiationPhotonsPerSecond");
  Param.GetArray(RadiationEnergyBins, 
		 "Physics.ActiveParticles.RadiationParticle.RadiationEnergyBins");
  Param.GetArray(RadiationSED, 
		 "Physics.ActiveParticles.RadiationParticle.RadiationSED");
  
#else

  RadiationNumSources = NumberOfRadiationParticles;
  RadiationSEDNumberOfBins = NumberOfEnergyBins;
  RadiationEnergyBins = new float[RadiationSEDNumberOfBins];
  RadiationSED = new float[RadiationSEDNumberOfBins];
  for(int i = 0; i < RadiationSEDNumberOfBins; i++)
    {
      RadiationEnergyBins[i] = RadiationEnergyInBin[i];
      RadiationSED[i] = RadiationBinSED[i];
    }
  // Example of what values may look like
  // RadiationEnergyBins[1] = 25.0;
  // RadiationSED[1] = 100.0/RadiationSEDNumberOfBins;
  // RadiationEnergyBins[2] = 55.0;
  // RadiationSED[2] = 100.0/RadiationSEDNumberOfBins;
  // RadiationEnergyBins[3] = 12.8;
  // RadiationSED[3] = 100.0/RadiationSEDNumberOfBins;
  // RadiationEnergyBins[4] = 1;
  // RadiationSED[4] = 100.0/RadiationSEDNumberOfBins;
  // RadiationEnergyBins[5] = 200.0;
  // RadiationSED[5] = 100.0/RadiationSEDNumberOfBins;
  RadiationLifetime = RadiationSourceLifeTime;
  RadiationPhotonsPerSecond = PhotonsPerSecond;

#endif
  /* Error check parameters */

  if (RadiativeTransfer == FALSE)
    ENZO_FAIL("RadiativeTransfer must be turned on with RadiationParticle.");
  if (MultiSpecies < 2)
    ENZO_FAIL("Multispecies with H2 must be turned on with RadiationParticle.");

  fprintf(stderr, "%s: Initialize RadiationParticle\n", __FUNCTION__);
  if(SUCCESS != ReadRadiationParameterFile())
    {
      fprintf(stderr, "ERROR: Failure to read Radiation Source File.\n");
      ENZO_FAIL("Failure in ActiveParticle_RadiationParticle.\n");
    }
      
  /* Add on the Particle Array Handlers */
  typedef ActiveParticleType_RadiationParticle ap;
  /* AttributeVector is a vector of ParticleAttributeHandler */
  AttributeVector &ah = ap::AttributeHandlers;
  /* Create the vector ah and populate with the basic members*/
  ActiveParticleType::SetupBaseParticleAttributes(ah);
  /* Add RadiationLifeTime to the vector ah via the push_back library call*/
  ah.push_back(new Handler<ap, float, &ap::RadiationLifeTime>("RadiationLifeTime"));
  return SUCCESS;
}

template <class active_particle_class>
int ActiveParticleType_RadiationParticle::BeforeEvolveLevel
(HierarchyEntry *Grids[], TopGridData *MetaData,
 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
 int ThisLevel, bool CallEvolvePhotons, 
 int TotalStarParticleCountPrevious[],
 int RadiationParticleID)
{
  
  FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  
  float creation_redshift = 0.0;
  float current_redshift = 0.0 , a = 0.0, dadt = 0.0;
  InitData *node = Root;
  ActiveParticleType **RadiationParticleList = NULL;
  int nParticles = 0, ipart = 0, j = 0;
  CosmologyComputeExpansionFactor(Time, &a, &dadt);
  current_redshift = (1 + InitialRedshift)/a - 1;
  CR = current_redshift;
  
  RadiationParticleList = 
	 ActiveParticleFindAll(LevelArray, &nParticles, RadiationParticleID);

  /* 
   * First lets see if its time to create this particle yet 
   * If a particle is ready mark it as alive
   */
  while(node->next != NULL)
    {
      creation_redshift = node->Redshift;
      // creation time is given w.r.t. the current time
      if (current_redshift > creation_redshift) 
	{
	  node->Create = false;
	}
      else
	{
	  if(nParticles >= NumberOfRadiationParticles)  { //Check that we don't duplicate
	    node->Create = false;
	  }
	  else {
	    node->Create = true;
	  }
	  /* Set up radiative transfer variables */
#ifdef TRANSFER
	  if (CallEvolvePhotons)
	    {
	      const double LConv = (double) TimeUnits / pow(LengthUnits,3);
	      RadiationParticleList = 
		ActiveParticleFindAll(LevelArray, &nParticles, RadiationParticleID);
	      RadiationSourceEntry* source;
	      ActiveParticleType_RadiationParticle *ThisParticle;
	      for (ipart = 0; ipart < nParticles; ipart++) {
		if (RadiationParticleList[ipart]->IsARadiationSource(Time)) {
		  ThisParticle =
		    static_cast<ActiveParticleType_RadiationParticle*>(RadiationParticleList[ipart]);
		  source = ThisParticle->RadiationSourceInitialize();
		  source->LifeTime   = RadiationLifetime;
		  source->Luminosity = RadiationPhotonsPerSecond * LConv;
		  source->EnergyBins = RadiationSEDNumberOfBins;
		  source->Energy = new float[RadiationSEDNumberOfBins];
		  source->SED = new float[RadiationSEDNumberOfBins];
		  source->Type = Isotropic; //Can be Isotropic, Beamed or Episodic
		  for (j = 0; j < RadiationSEDNumberOfBins; j++) {
		    source->Energy[j] = RadiationEnergyBins[j];
		    source->SED[j] = RadiationSED[j];
		  }
		  if (GlobalRadiationSources->NextSource != NULL)
		    GlobalRadiationSources->NextSource->PreviousSource = source;
		  GlobalRadiationSources->NextSource = source;


		}
	      }
	    }
#endif

	}
       node = node->next;
    }

  if(APDEBUG)
    {
       RadiationParticleList = 
	 ActiveParticleFindAll(LevelArray, &nParticles, RadiationParticleID);
       
       ActiveParticleType_RadiationParticle *ThisParticle;
       for (ipart = 0; ipart < nParticles; ipart++) {
	 if (RadiationParticleList[ipart]->IsARadiationSource(Time)) {
	   ThisParticle =
	     static_cast<ActiveParticleType_RadiationParticle*>(RadiationParticleList[ipart]);
		  
	   ThisParticle->PrintInfo();
	 }
       }
    }

  return SUCCESS;
}


int ActiveParticleType_RadiationParticle::EvaluateFormation(grid *thisgrid_orig, 
							    ActiveParticleFormationData &data)
{
  if (CheckForCreateParticle(Root) == false)  
    return SUCCESS;  //No particle to form yet
  RadiationParticleGrid *thisGrid = static_cast<RadiationParticleGrid *>(thisgrid_orig);
  int i = 0, j = 0, k = 0;
  int GhostZones = thisGrid->GridStartIndex[0];
  int mindex = 0;
  int *GridIndices = NULL;
  float creation_redshift = 0.0;
  FLOAT edge[3], cellwidth=0.0;
  InitData *node = Root, *pnode = Root, *curnode = NULL;
  FLOAT ppos[3];
  int GridDimension[3] = {thisGrid->GridDimension[0],
                          thisGrid->GridDimension[1],
                          thisGrid->GridDimension[2]};

  /* Covenience local */
  edge[0] = thisGrid->CellLeftEdge[0][0];
  edge[1] = thisGrid->CellLeftEdge[1][0];
  edge[2] = thisGrid->CellLeftEdge[2][0];
  cellwidth = thisGrid->CellWidth[0][0];
  while(node->next != NULL)
    {
      if(node->Create == false || node->Alive ==true) //Check which particles are ready to be formed
	{
	  node = node->next;
	  continue;
	}
      for(i = 0; i < 3; i++)
	  ppos[i]  = node->Position[i];
      creation_redshift = node->Redshift;
      /* Find the indices in the grid */
      GridIndices = GetGridIndices(ppos, edge, cellwidth);
		  
      /* Only if the indices are all good-looking, and if this grid is at the 
	 finest level in the hierarchy, then insert the particles; otherwise,
	 this particle should be created in one of the other grids! */
      
      if ( (GridIndices[0] >= GhostZones)               && 
	   (GridIndices[0] < (GridDimension[0] - GhostZones)) && 
	   (GridIndices[1] >= GhostZones)               && 
	   (GridIndices[1] < (GridDimension[1] - GhostZones)) &&
	   (GridIndices[2] >= GhostZones)               && 
	   (GridIndices[2] < (GridDimension[2] - GhostZones)) ) 
	{	      
	  mindex = (GridIndices[2] * GridDimension[1] + GridIndices[1]) * 
	    GridDimension[0] + GridIndices[0];
	  
	  if (thisGrid->BaryonField[thisGrid->NumberOfBaryonFields][mindex] == 0.0) 
	    {
	      
	      /*
	       * ====================================================================
	       * PARTICLE CREATION
	       * ====================================================================
	       */
	      ActiveParticleType_RadiationParticle *np = 
		new ActiveParticleType_RadiationParticle();
	      data.NewParticles[data.NumberOfNewParticles++] = np;
	      np->BirthTime = thisGrid->Time;
	      
	      //Mass
	      np->Mass = 0.0; //massless source particles
	      
	      //Positions
	      np->pos[0] = ppos[0]; 
	      np->pos[1] = ppos[1]; 
	      np->pos[2] = ppos[2]; 
	     
	      // Velocities are initialised to zero
	      np->vel[0] = 0.0;
	      np->vel[1] = 0.0;
	      np->vel[2] = 0.0;
	      np->RadiationLifeTime = RadiationSourceLifeTime;
	      np->type =  np->GetEnabledParticleID(); //PARTICLE_TYPE_RAD;
	      np->level = data.level;
	      np->GridID = data.GridID;
	      np->CurrentGrid = thisGrid;
	      fprintf(stderr, "%s: A radiation particle inserted at (%"PSYM",%"PSYM",%"PSYM") " \
		      "with v=(%f,%f,%f), m=%f, type=%ld, redshift = %f\n", __FUNCTION__,
		      np->pos[0], 
		      np->pos[1],
		      np->pos[2],
		      np->vel[0], 
		      np->vel[1],
		      np->vel[2],
		      np->Mass,
		      np->type,
		      CR);
	      np->Metallicity = 0.0;
	      node->Alive = true;
	      //node->Deleteme = true;
	    } // refinement field
	} // indices all OK
     
      curnode = node;
      node = node->next;
      if(curnode->Deleteme == true)
	{
	  /* Pop this entry from the LL */
	  if(curnode == Root)
	    {
	      Root = curnode->next;
	      pnode = Root;
	      delete curnode;
	    }
	  else
	    {
	      pnode->next = node->next;
	      delete node;
	    }
	  
	}
      else
	{
	  pnode = node;
	}
    } //end while
  /*
   * This is not ideal but if RadiativeTransfer is turned on 
   * midrun then the RT values are all zero rather than 
   * their default values. Here we set the RT values to 
   * defaults suitable for the radiation particle. 
   */
  
  SetRadiationDefaults();

  if (APDEBUG && data.NumberOfNewParticles > 0) {
    fprintf(stderr, "AP_RadiationParticle: Have created %"ISYM" new particles\n",
	    data.NumberOfNewParticles);
  }


  return SUCCESS;
}

// Radiation particle Feedback - done through Radiative Transfer
int ActiveParticleType_RadiationParticle::EvaluateFeedback
(grid *thisgrid_orig, ActiveParticleFormationData &supp_data)
{
  return SUCCESS;
}

int ActiveParticleType_RadiationParticle::SetFlaggingField(LevelHierarchyEntry *LevelArray[], 
							   int level, int TopGridDims[], 
							   int PopIIIParticleID)
{
  return SUCCESS;
}


void ActiveParticleType_RadiationParticle::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{

  /*
   * For every entry in the ActiveparticleFormationData struct, we have a bool
   * here. If a flag is set to true, some derived data is calculated and attached
   * to the ActiveParticleFormationData struct.
   *
   * DarkMatterDensity, CoolingTime, Temperature, MetalField, H2Fraction, and
   * CoolingRate return an array containing the quantity. Since it is sometimes expensive
   * to cache these fields, they should only be turned on if your formation or feedback
   * algorithm require it.
   */

  flags.DarkMatterDensity = false;
  flags.CoolingTime = false;
  flags.Temperature = false;
  flags.MetalField = false;
  flags.H2Fraction = false;
  flags.CoolingRate = false;
}


/*
 * Get the Cell indices given the positions and the cell bounds. 
 */
int* ActiveParticleType_RadiationParticle::GetGridIndices(double *position, FLOAT *GridLeftEdge, float CellWidth)
{
  static int ind[3];
  ind[0] = (int) ((position[0] - GridLeftEdge[0]) / CellWidth);
  ind[1] = (int) ((position[1] - GridLeftEdge[1]) / CellWidth);
  ind[2] = (int) ((position[2] - GridLeftEdge[2]) / CellWidth);
  return ind;
}


void ActiveParticleType_RadiationParticle::SetRadiationDefaults()
{
  if (RadiativeTransfer == TRUE)
    {
     
      if(RadiativeTransferInitialHEALPixLevel == 0)
	RadiativeTransferInitialHEALPixLevel = 2;
      if(RadiativeTransferPropagationSpeedFraction == 0) 
	RadiativeTransferPropagationSpeedFraction = 1.0;
      if(RadiativeTransferRaysPerCell == 0.0)
	RadiativeTransferRaysPerCell = 5.1;
      if(RadiativeTransferCoupledRateSolver == 0)
	RadiativeTransferCoupledRateSolver = 1;
      if(RadiativeTransferFluxBackgroundLimit == 0.0)
	RadiativeTransferFluxBackgroundLimit = 0.1;
      if(RadiativeTransferSplitPhotonRadius == 0)
	RadiativeTransferSplitPhotonRadius = FLOAT_UNDEFINED; // kpc
      if(RadiativeTransferPhotonMergeRadius == 0.0)
	RadiativeTransferPhotonMergeRadius = 10.0;
      if(RadiativeTransferTimestepVelocityLimit == 0.0)
	RadiativeTransferTimestepVelocityLimit = 100.0;
      if(RadiativeTransferAdaptiveTimestep == FALSE)
	RadiativeTransferAdaptiveTimestep = TRUE;
      if(RadiativeTransferPeriodicBoundary == 0)
	RadiativeTransferPeriodicBoundary = 1;
      
    }
  return;
}

bool ActiveParticleType_RadiationParticle::IsARadiationSource(FLOAT Time)
{
  return true;
}


int ActiveParticleType_RadiationParticle::ReadRadiationParameterFile()
{

  int i = 0;
  FILE *fd = NULL; 
  char line[MAX_LINE_LENGTH];
  
  InitData *node = NULL, *pnode = NULL, *cnode = NULL;
  Root = new InitData;
  Root->next = NULL;
  pnode = Root;
  //Create nodes for LL
  for(i = 0; i < RadiationNumSources; i++)
    {
      node = new InitData;
      node->Redshift = 0.0;
      node->Create = false;
      node->Alive = false;
      node->Deleteme = false;
      node->next = NULL;
      pnode->next = node;
      pnode = node; 
    }
  cnode = Root;
  i = 0;
  while(i < RadiationNumSources)
    {
      /* Open the file and read in the MBH masses and locations */
      if ((fd = fopen(RadiationSourcesFileName, "r")) == NULL) 
	{
	  
	  fprintf(stderr, "%s: Error opening file %s\n", __FUNCTION__, 
		  RadiationSourcesFileName);
	  delete Root;	      
	  return FAIL;     
	} 
      else 
	{
	  /* Now read the files line by line */
	  while (fgets(line, MAX_LINE_LENGTH, fd) != NULL) 
	    {
	      if (line[0] != '#' && strlen(line) > 1) 
		{  
		  /* order: Position[3], Creation redshift */
		  
		  if (sscanf(line, " %"PSYM"  %"PSYM"  %"PSYM"  %"FSYM, 
			     &(cnode->Position[0]), &(cnode->Position[1]), 
			     &(cnode->Position[2]), &(cnode->Redshift)) != NUMPARAMS) 
		    {
		      fprintf(stderr, "%s: Unrecognised line found in %s - ignoring\n", 
			      __FUNCTION__, RadiationSourcesFileName);
		      fprintf(stderr, "%s: line[0] = %c\t length = %ld\n", __FUNCTION__, line[0], 
			      strlen(line));
		      continue;
		   
		    }
		  else
		    {
		      fprintf(stdout, "Particle Positions = (%"PSYM", %"PSYM", %"PSYM")\n", 
			      cnode->Position[0], cnode->Position[1], cnode->Position[2]);
		      fprintf(stdout,"Particle will be created at z <= %f\n", cnode->Redshift);
		    }

		}
	    }
	  cnode = cnode->next;
	  if(cnode->next == NULL)
	    break;
	}
    }
		      



  return SUCCESS;
}

/* Check to see if we should form any radiation particle yet */
bool CheckForCreateParticle(InitData *Root)
{
   InitData *node = Root;
   
   while(node->next != NULL)
     {
       if(node->Create == true)
	 return true;
       node = node->next;
     }
   return false;

}



namespace {
  ActiveParticleType_info *RadiationParticleInfo = 
    register_ptype <ActiveParticleType_RadiationParticle> ("RadiationParticle");
}

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_RadiationParticle::AttributeHandlers;

