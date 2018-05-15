.. _active_particles:

Active Particles: Particles That do Non-local Feedback
======================================================

The particles described in :ref`star_particles` for the most part do feedback entirely local to a grid. 
This presents difficulty when one wants to do a more sophisticated form of feedback, for example a sink particle 
that can accrete gas from more than one grid, or a particle that can emit radiation. Active particles are a 
software framework for developing new particles and running simulations that include particles that run complex 
non-local feedback or formation algorithms. For a general introduction, see `Intro to Active Particles on figshare <https://figshare.com/articles/Enzo_Active_Particle_API/1254507>`

The following are Active Particle Types currently implemented in Enzo. 

Cen & Ostriker
--------------
Select this method by including this line in your .enzo file:
AppendActiveParticleType = CenOstriker

*Source: ActiveParticle_CenOstriker.C*

Identical in function to star particle Method 0, it reads the same parameters described here :ref:`StarParticleParameters`. Stochastic star formation is supported, and is enabled by setting the parameter .........

Accreting Particle
------------------
Select this method by including this line in your .enzo file: AppendActiveParticleType = AccretingParticle

*Source: ActiveParticle_AccretingParticle.C*

This routine implements the method of `Krumholz et al. (2004) <http://adsabs.harvard.edu/abs/2004ApJ...611..399K>` on top of the star particle Method 4 (sink particles). 

Parameters:
* OverflowFactor (default 1.01)
* LinkingLength (default 4)
* AccretionRadius (default 4)

Radiating Particle
------------------
Select this method by including this line in your .enzo file: AppendActiveParticleType = RadiationParticle

*Source: ActiveParticle_RadiationParticle.C*

Identical in function to star particle Method 5 (Radiative Stellar Clusters), it reads the same parameters described here :ref:`StarParticleParameters`. 

Skeleton Code
-------------
This method, which executes but doesn't do much, can be selected by including this line in your .enzo file: AppendActiveParticleType = Skeleton.

*Source: ActiveParticle_Skeleton.C*

This well commented code should be used for constructing new active particle types. 
