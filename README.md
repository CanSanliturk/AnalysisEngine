# AnalysisEngine
3D Structural Analysis Engine

A finite element tool to solve structural models.

Current Capabilities
---------------------
* Tool can be used to create structural model including frame elements, truss elements, various shell elements, restraints and element releases.

* Bilinear membrane, incompatible membrane and drilling membrane are available for shell elements membrane action. For plate action, Mindlin's 4-Node Thick Plate formulation and 8-Node Thick Plate Formulation (Serendipity plate) are implemented. 

* Static analysis can be executed on created structural model.

* Modal analysis can be executed on created structural model.

* Implicit Newmark-Beta Algorithm can be called for dynamic analysis of structures with both time and space varying loads.

Future Development
--------------------
* Material and geometrical nonlinearities will be considered.

* Buckling analysis will be added.

Solver
------------------
* Eigen and Armadillo linear algebra solvers are implemented. 

[![Build status](https://ci.appveyor.com/api/projects/status/6hinynxesc3hqfnk/branch/master?svg=true)](https://ci.appveyor.com/project/CanSanliturk/analysisengine/branch/master)
