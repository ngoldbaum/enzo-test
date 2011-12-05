/*-*-C++-*-*/
/***********************************************************************
/
/  STAR PARTICLE STRUCTURE
/
/  written by: John Wise
/  date:       September, 2005
/  modified1:  John Wise
/  date:       March, 2009 (converted into a class)
/  modified2:  John Wise, Greg Bryan, Britton Smith, Cameron Hummels,
/              Matt Turk
/  date:       May, 2011 (converting from Star to ActiveParticle)
/
/  PURPOSE:
/
************************************************************************/
#ifndef __ACTIVE_PARTICLE_H
#define __ACTIVE_PARTICLE_H

#include "hdf5.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "StarBuffer.h"

struct ActiveParticleFormationData;
struct ActiveParticleFormationDataFlags;
class ParticleBufferHandler;

class ActiveParticleType
{
public:
  void static ConstructData(grid *_grid,
			    ActiveParticleFormationDataFlags &flags,
			    ActiveParticleFormationData &data);
  void static DestroyData(grid *_grid,
			  ActiveParticleFormationData &data);
  int static WriteDataset(int ndims, hsize_t *dims, char *name, hid_t group,
			  hid_t data_type, void *data);
  int static ReadDataset(int ndims, hsize_t *dims, char *name, hid_t group,
			 hid_t data_type, void *read_to);

  /* Several pure virtual functions */
  
  /* This should return the number of new star particles created, and should
   * create them. */
//  ActiveParticleType(){};
//  ~ActiveParticleType(){};
//  ActiveParticleType(ActiveParticleType*& part){};

  ActiveParticleType(void);
  ActiveParticleType(grid *_grid, ActiveParticleFormationData &data);
  ActiveParticleType(grid *_grid, int _id, int _level);
  ActiveParticleType(ActiveParticleType*& part);
  ~ActiveParticleType(void);

  void operator=(ActiveParticleType *a);

  template <class active_particle_class> active_particle_class *copy(void);

  int   ReturnID(void) { return Identifier; };
  double ReturnMass(void) { return Mass; };
  float ReturnBirthTime(void) { return BirthTime; };
  float ReturnDynamicalTime(void) { return DynamicalTime; };
  float ReturnMetallicity(void) { return Metallicity; };
  star_type ReturnType(void) {return type; };
  int   ReturnLevel(void) { return level; };

  void  ReduceLevel(void) { level--; };
  void  IncreaseLevel(void) { level++; };
  void  SetLevel(int i) { level = i; };
  void  SetGridID(int i) { GridID = i; };
  grid *ReturnCurrentGrid(void) { return CurrentGrid; };
  void  AssignCurrentGrid(grid *a) { this->CurrentGrid = a; };
  void  AddMass(double dM) { Mass += dM; };
  void  AdjustMassByFactor(double factor) { Mass *= factor; };

  FLOAT *ReturnPosition(void) { return pos; }
  float *ReturnVelocity(void) { return vel; }
  void   ConvertAllMassesToSolar(void);
  void   ConvertMassToSolar(void);
  void   Merge(ActiveParticleType *a);
  //bool  MergableMBH(ActiveParticleType *a);
  float Separation(ActiveParticleType *a);
  float Separation2(ActiveParticleType *a);
  float RelativeVelocity2(ActiveParticleType *a);
  void  UpdatePositionVelocity(void);
  void  MirrorToParticle(void);
  void  CopyFromParticle(grid *_grid, int _id, int _level);
  int   DisableParticle(LevelHierarchyEntry *LevelArray[]);
  int   SphereContained(LevelHierarchyEntry *LevelArray[], int level, 
			float Radius);
  void  PrintInfo(void);
  //void  ActivateNewStar(FLOAT Time, float Timestep);
  //void  DeleteCopyInGrid(void);
  //int   DeleteCopyInGridGlobal(LevelHierarchyEntry *LevelArray[]);

  /* Virtual and pure virtual functions in this base class */

  //virtual void FillBuffer(ParticleBufferHandler *buffer);
  //virtual ParticleBufferHandler *AllocateBuffer(void) = 0;
  virtual bool IsARadiationSource(FLOAT Time) { return FALSE; };
  virtual bool Mergable(ActiveParticleType *a);
  virtual int GetEnabledParticleID(int id = -1) = 0;

#ifdef TRANSFER
  RadiationSourceEntry* RadiationSourceInitialize(void);
#endif

protected:
  grid *CurrentGrid;
  FLOAT	pos[MAX_DIMENSION];
  float vel[MAX_DIMENSION];
  double Mass;		// Msun
  float BirthTime;
  float DynamicalTime;      
  float Metallicity;
  PINT Identifier;
  int level;
  int GridID;
  star_type type;
  
  bool Active;
  int EnabledParticleID;

private: /* Cannot be accessed by subclasses! */
  
  friend class grid;
  friend class ActiveParticleType_info;
  friend class ParticleBufferHandler;

};

struct ActiveParticleFormationData {
  int NumberOfNewParticles;
  int MaxNumberOfNewParticles;
  ActiveParticleType **NewParticles;
  /* This is where all the pointers that normally get passed into
   * formation routines gets placed. Things like fractional h2, dark
   * matter density, etc etc. Anything that's derived.  It's okay to
   * add to this.  */
  float *DarkMatterDensity;
  float *H2Fraction;
  float *CoolingTime;
  float *CoolingRate;
  float *Temperature;
  float *TotalMetals;
  float DensityUnits;
  float LengthUnits;
  float TemperatureUnits;
  float TimeUnits;
  float VelocityUnits;
  double MassUnits;
  int DensNum;
  int Vel1Num;
  int Vel2Num;
  int Vel3Num;
  int TENum;
  int GENum;
  int MetalNum;
  int MetalIaNum;
  int ColourNum;
  int level;
  int GridID;
};

const struct ActiveParticleFormationData data_default = {
  0,        // NumberOfNewParticles
  0,        // MaxNumberOfNewParticles
  NULL,     // NewParticles
  NULL,     // DarkMatterDensity
  NULL,     // H2Fraction
  NULL,     // CoolingTime
  NULL,     // CoolingRate
  NULL,     // Temperature
  NULL,     // TotalMetals
  0.0,      // DensityUnits
  0.0,      // LengthUnits
  0.0,      // TemperatureUnits
  0.0,      // TimeUnits
  0.0,      // VelocityUnits
  0.0,      //  MassUnits
  -1,       // DensNum
  -1,       // Vel1Num
  -1,       // Vel2Num
  -1,       // Vel3Num
  -1,       // TENum
  -1,       // GENum
  -1,       // MetalNum
  -1,       // MetalIaNum
  -1,       // ColourNum
  -1,       // level
  -1        // GridID
};


struct ActiveParticleFormationDataFlags {
  /* For every entry in the ActiveParticleFormationData struct, we
   * have a bool here. */
  bool DarkMatterDensity;
  bool H2Fraction;
  bool CoolingTime;
  bool CoolingRate;
  bool Temperature;
  bool UnitConversions;
  bool DataFieldNumbers;
  bool MetalField;
};

const struct ActiveParticleFormationDataFlags flags_default = {
  false,    // DarkMatterDensity
  false,    // H2Fraction
  false,    // CoolingTime
  false,    // CoolingRate
  false,    // Temperature
  false,    // UnitConversions
  false,    // DataFieldNumbers
  false     // MetalField
};


//! maps the name of a plug-in to a pointer of the factory pattern
class ActiveParticleType_info;
typedef std::map<std::string, ActiveParticleType_info *> ActiveParticleMap;

ActiveParticleMap &get_active_particle_types();

void EnableActiveParticleType(char *active_particle_type_name);

class ParticleBufferHandler
{
public:
  ParticleBufferHandler() {};
  ~ParticleBufferHandler() {};
  /*virtual void WriteBuffers(hid_t group);*/
  virtual void FillBuffer(ActiveParticleType *np);
  virtual void UnpackBuffer(ActiveParticleType *np);
  int return_proc(void) { return dest_processor; }
  int return_grid(void) { return GridID; }
  int return_type(void) { return type; }
  void set_grid(int num) { this->GridID = num; }
  void set_proc(int num) { this->dest_processor = num; }
  bool compare_grids(const ParticleBufferHandler *lhs, const ParticleBufferHandler *rhs) {
    if (lhs->GridID < rhs->GridID) return true;
    else return false;
  };

protected:
  FLOAT	pos[MAX_DIMENSION];
  float vel[MAX_DIMENSION];
  double Mass;		// Msun
  float BirthTime;
  float DynamicalTime;      
  float Metallicity;
  PINT Identifier;
  int dest_processor;
  int level;
  int GridID;
  int type;

};

/* Comparer functions for sorting particle buffers with std::sort */

struct cmp_ap_grid {
  bool operator()(ParticleBufferHandler* const& a, ParticleBufferHandler* const& b) const {
    if (a->return_grid() < b->return_grid()) return true;
    else return false;
  }
};

struct cmp_ap_proc {
  bool operator()(ParticleBufferHandler* const& a, ParticleBufferHandler* const& b) const {
    if (a->return_proc() < b->return_proc()) return true;
    else return false;
  }
};

class ActiveParticleType_info
{
public:
       
  /* We will add more functions to this as necessary */
  ActiveParticleType_info
  (std::string this_name,
   int (*ffunc)(grid *thisgrid_orig, ActiveParticleFormationData &data),
   void (*dfunc)(ActiveParticleFormationDataFlags &flags),
   ParticleBufferHandler *(*abfunc)(ActiveParticleType *np),
   void (*unfunc)(ActiveParticleType *np, ParticleBufferHandler **buffer, int place),
   int (*ifunc)(),
   int (*feedfunc)(grid *thisgrid_orig, ActiveParticleFormationData &data),
   int (*writefunc)(ActiveParticleType *these_particles, int n, int GridRank, hid_t group_id),
   int (*readfunc)(ActiveParticleType **particles_to_read, int *n, int GridRank, hid_t group_id),
   ActiveParticleType *particle
   ){
    this->formation_function = ffunc;
    this->describe_data_flags = dfunc;
    this->allocate_buffer = abfunc;
    this->unpack_buffer = unfunc;
    this->particle_instance = particle;
    this->initialize = ifunc;
    this->feedback_function = feedfunc;
    this->write_function = writefunc;
    this->read_function = readfunc;
    get_active_particle_types()[this_name] = this;
  }

  static int count(){return get_active_particle_types().size();}
  int GetEnabledParticleID(){return this->MyEnabledParticleID;}

  int Enable(){
    /* 0-indexed */
    this->MyEnabledParticleID = this->TotalEnabledParticleCount++;
    this->particle_instance->GetEnabledParticleID(this->MyEnabledParticleID);
    return this->MyEnabledParticleID;
  }

  int (*initialize)(void);
  int (*formation_function)(grid *thisgrid_orig, ActiveParticleFormationData &data);
  int (*feedback_function)(grid *thisgrid_orig, ActiveParticleFormationData &data);
  int (*write_function)(ActiveParticleType *these_particles, int n, int GridRank, hid_t group_id);
  int (*read_function)(ActiveParticleType **particles_to_read, int *n, int GridRank, hid_t group_id);
  void (*describe_data_flags)(ActiveParticleFormationDataFlags &flags);
  ParticleBufferHandler* (*allocate_buffer)(ActiveParticleType *np);
  void (*unpack_buffer)(ActiveParticleType *np, ParticleBufferHandler **buffer, int place);
  ActiveParticleType* particle_instance;
private:
  /* This is distinct from the global as a redundant error-checking
     pattern */
  static int TotalEnabledParticleCount;
  int MyEnabledParticleID; /* Defaults to 0 */
  int *EnabledParticleIDPointer;

};

template <class active_particle_class>
ActiveParticleType_info *register_ptype(std::string name)
{
  active_particle_class *pp = new active_particle_class();
  ActiveParticleType_info *pinfo = new ActiveParticleType_info
    (name,
     (&active_particle_class::EvaluateFormation),
     (&active_particle_class::DescribeSupplementalData),
     (&active_particle_class::AllocateBuffer),
     (&active_particle_class::UnpackBuffer),
     (&active_particle_class::InitializeParticleType),
     (&active_particle_class::EvaluateFeedback),
     (&active_particle_class::WriteToOutput),
     (&active_particle_class::ReadFromOutput),
     pp);
  return pinfo;
}

#define ENABLED_PARTICLE_ID_ACCESSOR					\
  int GetEnabledParticleID(int myid = -1) {				\
    static int ParticleID = -1;						\
    if (myid >= 0) {							\
      if (ParticleID != -1) ENZO_FAIL("Setting Particle ID Twice!");	\
      ParticleID = myid;						\
    }									\
    return ParticleID;							\
  };


#endif

