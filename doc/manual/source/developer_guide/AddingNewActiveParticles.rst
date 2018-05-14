.. _AddingNewActiveParticles:

Adding new Active Particles
===========================

One goal of the active particle infrastructure is to make it possible to add a
new particle types simply by adding two new files to the Enzo source code. This
document will walk you through this process. By the end of this tutorial you
should have a modified version of Enzo including a skeleton version of your
particle. From there you will be able to modify the skeleton particle to include
the physics you would like your particle to capture.

In this guide we will be creating a new active particle named
AwesomeParticle. You should come up with a descriptive name for your particle
and replace occurances of the AwesomeParticle name with the name of your
particle below. To start with, navigate to the ``enzo`` source directory and
make a copy of the `SkeletonParticle`:

    $ hg cp ActiveParticle_Skeleton.C ActiveParticle_Awesome.C
    $ hg cp ActiveParticle_Skeleton.h ActiveParticle_Awesome.h

Next you will need to make sure the new ``.C`` file you have created actually
gets built by Enzo's build system. To do this, edit ``Make.config.objects`` and
add an entry for ``ActiveParticle_Awesome.o`` in the appropriate alphabetically
sorted location::

    --- a/src/enzo/Make.config.objects
    +++ b/src/enzo/Make.config.objects
    @@ -26,6 +26,7 @@ OBJS_CONFIG_LIB = \
             ActiveParticleResetAccelerations.o \
             ActiveParticleRoutines.o \
             ActiveParticle_AccretingParticle.o \
    +        ActiveParticle_Awesome.o \
             ActiveParticle_CenOstriker.o \
             ActiveParticle_DisableParticle.o \
             ActiveParticle_GalaxyParticle.o \

Next, you should replace all occurances of the string ``Skeleton`` with
``Awesome``. You can either do that manually in your editor or you could use
``sed``:

  $ ls ActiveParticle_Awesome.[h,C] | xargs sed -i 's/Skeleton/Awesome/g'

If you use a Mac, the command is slightly different:

  $ ls ActiveParticle_Awesome.[h,C] | xargs sed -i '' -e 's/Skeleton/Awesome/g'

Now try building Enzo:
  
  $ make -j4

You should see in the stdout from the compilation process that
ActiveParticle_Awesome.C is getting compiled. Assuming everything was done
correctly, Enzo should build successfully. Congratulations! You've now added a
skeleton version of your particle. To make your particle actually do something
useful, you will need to modify the new particle you've added, filling in the
functions in the active particle API. This process is described below.


  
Available hooks
---------------

Rather than needing to modify the main evolution loop or any other place in the
Enzo source code where you would like to add special treatment for a particle,
the active particle framework provides hooks into various places in the Enzo
codebase where an activated active particle can modify the simulation.

Below we describe what these hooks are
