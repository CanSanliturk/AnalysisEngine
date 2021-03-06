# AnalysisEngine
3D Structural Analysis Engine

A finite element tool to solve structural models.

Current Capabilities
---------------------
* Tool can be used to create structural model including beam and truss elements, restraints and element releases.

* Static analysis can be executed on created structural model.

* Modal analysis can be executed on created structural model.

Future Development
--------------------
* Damping matrices will be introduced to members (and to overall structural) which uses Rayleigh damping.

* Newmark-Beta algorithm will be implemented for dynamic analysis.

* Material and geometrical nonlinearities will be considered.

* More element types will be added like Timoshenko Beam, shell element etc.

Solver
------------------
* Eigen and Armadillo linear algebra solvers are implemented. 

[![Build status](https://ci.appveyor.com/api/projects/status/6hinynxesc3hqfnk/branch/master?svg=true)](https://ci.appveyor.com/project/CanSanliturk/analysisengine/branch/master)
