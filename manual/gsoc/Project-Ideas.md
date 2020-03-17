# Ideas for Google Summer of Code projects
In this page you'll find some ideas for projects to be developed in the context of the [Google Summer of Code](https://summerofcode.withgoogle.com/). If you are a student looking forward to apply, please remember that these are just suggestions and general ideas, you are free to propose you own project and **strongly encouraged** to discuss with the regular developers your ideas before submitting an application.

<!--    - [Couple MBDyn with OpenFOAM](#couple-mbdyn-with-openfoam)-->

## Ideas by category
  - **Modeling Capabilities**
    - [Implement new integration schemes](#implement-new-integration-schemes)
    - [Improve the cycloidal rotor module](#improve-the-cycloidal-rotor-module)
    - [Friction in joints](#friction-in-joints)
    - [Embedded Optimization](#embedded-optimization)
    - [Ground Vehicle Model Development](#ground-vehicle-model-development)
    - [Improve the tire model](#improve-the-tire-model)
    - [Implement a PID controller element](#PID-controller-element)
    - [Revamp the modal joint](#revamp-the-modal-joint)
  - **IPC/RT**
    - [Improve FMI support](#improve-fmi-support)
    - [Improve coupling with other software ](#improve-coupling-with-other-software)
  - **User Interface**
    - [Blendyn Development](#blendyn-development)
    - [Port Blendyn to Blender 2.8](#port-blendyn-to-blender-28)
    - [FreeCAD gui](#freecad-gui)
    - [Flexible element graphical interface](#flexible-element-graphical-interface)
    - [Improve the parser](#improve-the-parser)
  - **Miscellanea**
    - [Convergence criteria](#convergence-criteria)
    - [Libraries update](#libraries-update)
    - [Package update](#package-update)
    - [Configuration update](#configuration-update)
    - [Cascaded Analysis](#cascaded-analysis)
    - [Automatic Differentiation](#automatic-differentiation)


-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------

## Modeling Capabilities
<!--
### Couple MBDyn with OpenFOAM
The MBDyn software implements basic aerodynamic models which do not allow to properly model complex flow structures. A common approach is to rely on CFD simulations to validate the models used within MBDyn.

A more precise and faster approach would be to have the flow simulations embedded directly into the multibody dynamics model.

A preliminary model with a 1-DoF communication has been implemented and can serve as a starting point.

This project is open-ended, but some avenues for thought are:
  - implement translation, rotation, or both into the fluid structure interaction approach
  - allow for mesh deformation inside the flow solver, the multibody dynamics solver, or both
  - use OpenFOAM as the client and MBDyn as the server for the simulation or attempt a different appproach which would require more substantial code changes to both OpenFOAM and MBDyn

Requirements:
  - use MBDyn for the structural part
  - use OpenFOAM for the fluid part
  - provide one or more test cases that exemplify all of the features
  - provide detailed installation instructions for anyone starting from the official releases of MBDyn and OpenFOAM

Documentation:
If this project interests you, contact the MBDyn mailing list or the mentor to have the proper initial documentation

**Category**: [Modeling Capabilities](#modeling-capabilities)   
**Programming Languages**: C++   
**Keywords**: aerodynamics, fluid-structure interaction, OpenFOAM   
**Priority**: High   
**Difficulty**: Intermediate/High   
**Mentors**: Louis Gagnon, TBD
-->

-------------------------------------------------------------------------------------------

### Implement new integration schemes
The MBDyn software supports several multi-step integration schemes, which must be written for first-order differential equations.

New single-step, second-order accurate schemes with controllable algorithmic dissipation have been formulated, which can be rearranged for first-order differential equations, and could be implemented in MBDyn without excessive effort.

Requirements:
  - familiarity with numerical integration
  - ability to interpret the code (numerical integration is currently buried in the core code)

Documentation:
If this project interests you, contact the MBDyn mailing list or the mentor to have the proper initial documentation

**Category**: [Modeling Capabilities](#modeling-capabilities)   
**Programming Languages**: C++
**Keywords**: numerical integration, differential equations   
**Priority**: High   
**Difficulty**: Low/Intermediate   
**Mentors**: [Pierangelo Masarati](@10102934)

-------------------------------------------------------------------------------------------

### Improve the cycloidal rotor module
During GSoC2016, a double multiple streamtube inflow model for [cycloidal rotors](https://en.wikipedia.org/wiki/Cyclorotor) was implemented into MBDyn. The code can be found [here](https://www.mbdyn.org/?News&id=28). It is desired to improve the current implementation by, for example, mitigating the following limitations:

  - eliminate divergence issues by carefully checking and fixing the code;
  - perform validation tests with experimental results;
  - implement averaging and filtering.

For this project, the student should review the work done during 2016 and come up with his/her own suggestions.

**Category**: [Modeling Capabilities](#modeling-capabilities)   
**Programming Languages**: C++   
**Keywords**: aeronautics, rotorcraft   
**Priority**: Medium   
**Difficulty**: Intermediate   
**Mentors**: Louis Gagnon, Giuseppe Quaranta

-------------------------------------------------------------------------------------------

### Friction in Joints
It would be nice if every MBDyn joint could optionally add friction to its internal forces and moments.

A possible to-do list for the project would be:
  - improve the friction model and add it to all joints (suggestion: try a transition force from stick to slip that is higher than slip force);
  - implement contacts with non-flat surfaces.
  - start from the total joint

**Category**: [Modeling Capabilities](#modeling-capabilities)  
**Programming Languages**: C++   
**Keywords**: mechanical friction       
**Priority**: Low   
**Difficulty**: Intermediate   
**Mentors**: [Marco Morandini](@marco.morandini)

-------------------------------------------------------------------------------------------

### Embedded Optimization
One of the strenghts of MBDyn simulation capabilities lies in dealing with (possibily kinematically underdetermined), overactuated systems. One prime example are biomechanical systems, where the overactuation originates from the redundancy of extensors or flexors muscles acting upon the same joint ([example](http://dx.doi.org/10.1177/1464419313490680)).

Control forces and moments on overactuated systems can be found using optimization techniques that, at the moment, have to be performed outside of MBDyn. The project aims at developing an internal, as flexible as possible, static optimization module, that can be configured easily to accept user-defined, state-dependent cost functions to be minimized.

**Category**: [Modeling Capabilities](#modeling-capabilities)   
**Programming Languages**: C++   
**Keywords**: biomechanics, overactuated systems, optimization       
**Priority**: Medium    
**Difficulty**: Intermediate/High    
**Mentors**: [Andrea Zanoni](@10260632), Louis Gagnon

-------------------------------------------------------------------------------------------

### Ground vehicle model development
A [semi-trailer truck](https://www.mbdyn.org/userfiles/documents/examples/semitrailerPub.tar.bz2) model was developed in MBDyn and this was a first step into adapting the software to ground vehicle analysis. See an [animation](https://www.mbdyn.org/userfiles/documents/examples/animCamion.gif) of the truck in action. It would be nice to take some concepts used for such an analysis and implement them directly in the MBDyn code.

Some possible improvements include: - take the scripting tricks and joints used in the truck model and incorporate them directly into the MBDyn code - improve the MBDyn output to include car, truck, and general ground vehicle data - implement a set of general procedures to carry-on with a vehicle model (ie: SAE testing procedures for fuel consumption, stability, and maneuverability) - develop an output element anologous to the **aircraft instruments** element, that exposes relevant measures to potential vehicle control systems (e.g. vehicle sideslip angle, wheel angular velocities, tyre slips, yaw angular acceleration, etc.)

**Category**: [Modeling Capabilities](#modeling-capabilities)   
**Programming Languages**: C++   
**Keywords**: vehicle dynamics, mbdyn scripting, modules    
**Priority**: Low    
**Difficulty**: Intermediate    
**Mentors**: Louis Gagnon, [Andrea Zanoni](@10260632)

-------------------------------------------------------------------------------------------

### Improve the tire model
The wheel4 module intends to reproduce the dynamic behavior of tires using a rigid ring type of model. During the development of the module, various avenues were envisioned to accelerate the resolution procedure. The documentation and parameter definition of the model would also benefit from a simplification.

Tentative to-do list:
  - simplify the module documentation
  - simplify the mathematical implementation
  - implement a simpler parameter definition method
  - add camber variation effects
  - filter the road profile directly inside MBDyn instead of during preprocessing with other computing environments (e.g. Octave or Matlab)
  - implement a lateral tire section profile (currently, only the longitudinal profile is considered by the model)
  - simplify the ability to change the orientation of the tire

**Category**: [Modeling Capabilities](#modeling-capabilities)   
**Programming Languages**: C++   
**Keywords**: vehicle dynamics, modules, tire mechanics    
**Priority**: Low    
**Difficulty**: Intermediate/High    
**Mentors**: Louis Gagnon, [Andrea Zanoni](@10260632)

-------------------------------------------------------------------------------------------

### Implement a PID controller element
The current method used to implement a PID controller into a MBDyn model usually puts the user in front of two choices:

  - **a** rely on an external software to compute the response of the controller, or
  - **b** use a combination of abstract nodes, genel elements, and drives to:

    - read the user-requested system response and compare it with the actual system's response, eventually capping the maximum and minimum response
    - filter the system's response to remove dependency of the PID on the user-chosen timestep
    - use abstract nodes to store outputs from genel elements which perform filtering and/or integration
    - use abstract nodes to store outputs from drives which read the velocity of the concerned degree  of freedom

Option (a) requires that the user couples MBDyn to external software.
The contributions to the Jacobian are in this case also not taken into consideration by MBDyn and that may lead to divergence.
This also slows down the simulation and adequate software is sometimes not available to the user.

As an alternative, option (b) is proved to work but renders the input files extremely complicated to write and read.
It also limits implementation to elements which can return a velocity for the D part.
Finally, the complexity of the input files of option (b) also makes the model highly prone to error which are then difficult to identify.

A dedicated high-level element which could handle multiple degrees of freedom, filter inputs, and provide a response based on the sum of (filtered) contributions of proportional,
integrated, and derivative errors, each multiplied by their respective coefficients would provide an attractive alternative to (a) and (b).

Depending on the skillset of the student, the element could be coded with Python into the preprocessor or with C++ as a MBDyn module.

List of tasks:

  - implement a 1 degree-of-freedom (DOF) PID element
  - add option to filter the inputs read by it
  - add option to filter the individual contributions of P, I, and D
  - add option to cap the maximum and minimum output values of each P, I, and D contributions
  - enable proper element-output in both text and NetCDF format
  - extend the element to be able to handle multiple DOFs.

The ideal candidate will demonstrate resourcefulness and a capacity to resolve most technical problems independently.

**Category**: [Modeling Capabilities](#modeling-capabilities)   
**Programming Languages**: C++ or Python   
**Keywords**: system response, control, autopilot, Proportional-Integrative-Derivative    
**Priority**: Medium    
**Difficulty**: Intermediate    
**Mentors**: [Louis Gagnon](@louis.gagnon), [Andrea Zanoni](@10260632)    

-------------------------------------------------------------------------------------------
### Revamp the modal joint
MBDyn's modal joint implements a Component Mode Synthesis (CMS) deformable body. The CMS element interacts with the multibody model only in specific nodes (interface nodes). The full dynamics of the body is condensed into the superposition of the responses of its modes of vibrations, that it outputs in the form of the time histories of its modal coordinates, together with the ridig body motion of a specific node used to define the floating reference frame, the modal node.

The modal joint implementation can be improved in several different ways. For example:
  - break down the modal joint element into a set of elements/components:
    - modal database: FEM model; split `const` ( = instantiated once) and `local` (= duplicated) data
    - place `const` data in shared database, `local` data in modal joint; avoid unnecessary copy of memory
    - use standard containers for structured data, for example
      - `Mat3xN *pOffsetFEMNodes --> std::vector`
      - `Mat3xN *pOffsetMBNodes --> std::vector`
      - `Mat3xN *pRotMBNodes --> std::vector`;
  - define a modal dynamics as a "deformable body" element (the structural dynamics part of the existing element);
  - define a modal clamp constraint element (as it is now, but separate from dynamics element);
  - add a total modal joint: this makes it possible to add arbitrary elements connected to the modal element;
  - add support for “reference” definition based on modal database geometry;
  - implement/fix support of modal joint in initial assembly.

**Category**: [Modeling Capabilities](#modeling-capabilities)   
**Programming Languages**: C++   
**Keywords**: lineal algebra, physics, eigenanalysis    
**Priority**: Low    
**Difficulty**: Intermediate/High    
**Mentors**: Giuseppe Quaranta, [Andrea Zanoni](@10260632), Louis Gagnon

-------------------------------------------------------------------------------------------

## IPC/RT
### Improve FMI support
During the GSoC2016 the Functional Mock-up Interface (FMI) standard was implemented into MBDyn. The code can be seen [here](https://www.mbdyn.org/?News&id=28). The current implementation could be improved by, for example, mitigating the following limitations:
  - the module does not support the string drive;
  - a segmentation fault occurs if FMU with version 1.0 is used with simulation type not supported by the FMU
  - requires a constant timestep for cosimulation

For this project, the student should review the work done during 2016 and come up with his/her own suggestions.

**Category**: [IPC/RT](#ipc-rt)   
**Programming Languages**: C++   
**Keywords**: fmi, co-simulation    
**Priority**: Low   
**Difficulty**: High   
**Mentors**: [Marco Morandini](@marco.morandini), [Pierangelo Masarati](@10102934)

-------------------------------------------------------------------------------------------

### Improve coupling with other software
Over the years, various interfaces to popular modeling software have been developed, in some cases even duplicated. However, the implementation in some cases remains unreliable. A good coupling interface rework is needed:
  - improve MATLAB hooks;
  - streamline Octave hooks;
  - add/improve Scilab hooks;
  - add/improve MATLAB/Octave/Scilab and Simulink/Scicos co-simulation interface.

Large projects:
  - add support for coupling with [SU2 CFD software](https://su2code.github.io/);
  - turn MBDyn into a module for Simulink/Scicos (see contrib/SimulinkInterface/, need an implicit extension);
  - add package-wide support for WinSocks in communication entities, including SimulinkInterface and ScicosInterface (using MSYS/MinGW).

**Category**: [IPC/RT](#ipc-rt)   
**Programming Languages**: C++, Python    
**Keywords**: fmi, co-simulation, Matlab, Octave, SU2, Scilab     
**Priority**: Low    
**Difficulty**: Intermediate/High    
**Mentors**: [Pierangelo Masarati](@10102934), Giuseppe Quaranta

-------------------------------------------------------------------------------------------

## User Interface
### Blendyn development
MBDyn is a multibody dynamics solver which comes without any default graphical user interface for pre- and post-processing. There exist a few standalone post-processing tools based on [EasyAnim](https://www.mbdyn.org/?Software_Download___EasyAnim), [OpenDx](http://www.opendx.org/) and [Blender](https://www.blender.org/).

However, [Blendyn](https://github.com/zanoni-mbdyn/blendyn), based on [Blender](https://www.blender.org/), is the most up-to-date.

See some example [videos](https://youtu.be/x5n0OgskIMc?list=PLTtFbiep140gc-f-x14ltv0N7YZNzvioF) of its output and the [tutorials](https://github.com/zanoni-mbdyn/blendyn/wiki/Tutorials) to understand better what Blendyn is about.

It is simple to use and generates 3D animations that represent the exact model movement and joints. Blendyn has got a great push in the development in the 2017 edition of the GSoC by the work of [Reddy Janga](https://github.com/janga1997), but some desirable features are still missing and several other need completion. For example:
  - the automatic visualization of deformable elements (see [PR#26](https://github.com/zanoni-mbdyn/blendyn/pull/26);
  - the [Live Animation](https://github.com/zanoni-mbdyn/blendyn/projects/5) project, that aims to allow running MBDyn directly from the Blender interface and animate the results in realtime;
  - the internal forces and stress/strain fields visualization of deformable components during the animation;
  - complete writing for each element in Blendyn the updated input for the MBDyn entity, in the Blender text editor.

**Category**: [User Interface](#user-interface)   
**Programming Languages**: Python    
**Keywords**: Blender, UI, post-process    
**Priority**: High    
**Difficulty**: Low/Intermediate    
**Mentors**: [Andrea Zanoni](@10260632), Louis Gagnon    
**ENTRY TEST**: Complete step 1 of standard MBDyn GSoC [entry test](https://gitlab.polimi.it/Pub/mbdyn/wikis/Google-Summer-of-Code/Entry-Test), then use the Blender Python API (or console) to create a simple Blender model

-------------------------------------------------------------------------------------------

### Port Blendyn to Blender 2.8
MBDyn is a multibody dynamics solver which comes without any default graphical user interface for pre- and post-processing. There exist a few standalone post-processing tools based on [EasyAnim](https://www.mbdyn.org/?Software_Download___EasyAnim), [OpenDx](http://www.opendx.org/) and [Blender](https://www.blender.org/).

However, [Blendyn](https://github.com/zanoni-mbdyn/blendyn), based on [Blender](https://www.blender.org/), is the most up-to-date.

See some example [videos](https://youtu.be/x5n0OgskIMc?list=PLTtFbiep140gc-f-x14ltv0N7YZNzvioF) of its output and the [tutorials](https://github.com/zanoni-mbdyn/blendyn/wiki/Tutorials) to understand better what Blendyn is about.

A new minor version of [Blender](https://www.blender.org/), version 2.8, has been in development for some time and it appears to be ready to 
soon be promoted to stable. To ensure forward compatibility and to exploit the numerous exciting features of Blender 2.8, however, the 
[Blendyn](https://github.com/zanoni-mbdyn/blendyn) code has to be updated, since both the Python API and the Blender UI are changing 
with the new version.

**Category**: [User Interface](#user-interface)   
**Programming Languages**: Python    
**Keywords**: Blender, UI, post-process    
**Priority**: High 
**Difficulty**: Low/Intermediate    
**Mentors**: [Andrea Zanoni](@10260632), Louis Gagnon    
**ENTRY TEST**: Complete step 1 of standard MBDyn GSoC [entry test](https://gitlab.polimi.it/Pub/mbdyn/wikis/Google-Summer-of-Code/Entry-Test), then use the Blender Python API (or console) to create a simple Blender model

-------------------------------------------------------------------------------------------

### FreeCAD gui
MBDyn is a multibody dynamics solver which comes without any default graphical user interface for pre- and post-processing. There exist a few standalone post-processing tools based on [EasyAnim](https://www.mbdyn.org/?Software_Download___EasyAnim), [OpenDx](http://www.opendx.org/) and [Blender](https://www.blender.org/).

The most up-to-date and still under active development is [Blendyn](https://github.com/zanoni-mbdyn/blendyn), which is based on Blender. It is a tool that allows to create [nice-looking](https://github.com/zanoni-mbdyn/blendyn/wiki/Tutorial-3-Full-PUMA-Main-Rotor) animations as well as to quickly visualize the results of a simulation.

What Blendyn is not, however, is a well-integrated piece in a typical engineering design toolchain. It is not trivial to convert an engineering CAD model into a model suited to be animated based on the results of an MBDyn simulation in Blendyn. Also, modeling in Blender for engineering purposes is rather inconvenient, since that is not the main purpose of the software.

In order to make MBDyn more user-friendly for the average engineer, a new tool must be introduced. Now comes [FreeCAD](https://www.freecadweb.org/?lang=en) on the scene. Freecad is probably the most advanced free engineering CAD software available and it offers a fully-fledged [Python API](https://www.freecadweb.org/wiki/FreeCAD_Scripting_Basics#Python_scripting_in_FreeCAD). Since Freecad supports CAD assembly operations, it could be a natural platform for a graphical pre- and post-processing tool for MBDyn.

The scope of this project is to begin **from scratch** the development of such a tool. Everything is starting anew, so we will make all the software design choices together!

**Category**: [User Interface](#user-interface)   
**Programming Languages**: Python    
**Keywords**: Freecad, UI, post-process    
**Priority**: Medium     
**Difficulty**: Low/Intermediate    
**Mentors**: [Andrea Zanoni](@10260632)    
**ENTRY TEST**: Complete step 1 of standard MBDyn GSoC [entry test](https://gitlab.polimi.it/Pub/mbdyn/wikis/Google-Summer-of-Code/Entry-Test), then use the Freecad Python API to create a simple Freecad model.

-------------------------------------------------------------------------------------------

### Flexible element graphical interface
MBDyn's modal joint implements a Component Mode Synthesis (CMS) deformable body. The CMS element interacts with the multibody model only in specific nodes (interface nodes). The full dynamics of the body is condensed into the superposition of the responses of its modes of vibrations, that it outputs in the form of the time histories of its modal coordinates, together with the ridig body motion of a specific node used to define the floating reference frame, the modal node.

Currently, the post processing (i.e. the visualization of the superposition of the rigid and the deformable motion of the body) for the element is not implemented in any of the existing post-processing tools. The project aims at filling this gap, by adding support for the flexible superelement to [Blendyn](https://github.com/zanoni-mbdyn/blendyn).

**Category**: [User Interface](#user-interface)   
**Programming Languages**: Python    
**Keywords**: Blender, UI, post-process    
**Priority**: Medium    
**Difficulty**: Intermediate    
**Mentors**: [Andrea Zanoni](@10260632), [Marco Morandini](@marco.morandini)    
**ENTRY TEST**: Complete step 1 of standard MBDyn GSoC [entry test](https://gitlab.polimi.it/Pub/mbdyn/wikis/Google-Summer-of-Code/Entry-Test), visualize the results of your simulation using [Blendyn](https://github.com/zanoni-mbdyn/blendyn), simulate the [ssbeam](https://www.mbdyn.org/userfiles/documents/examples/ssbeam.tar.gz) MBDyn example, visualize the trajectory of at least one of the FEM nodes of the ssbeam beam model, using a software of your choice (e.g. Octave, Scilab, Python, gnuplot, ...): refer to section 8.12.32 of the MBDyn input manual to gather information about the modal element output.

-------------------------------------------------------------------------------------------

### Improve the parser
The current MBDyn parser has a clear and consistent language. However, it does not let the user define if-else conditions, for-while loops, arrays of variables, and use of variables in the strings of node and element drives definitions. This leads the users to often rely on external scripting (e.g. Matlab or Python) to generate their own input files. This supplementary step creates additional opportunities for bugs to appear within a model. Adding these features directly to the MBDyn parser would thus make the generation of complex multibody models much simpler and more focused.

Another functionality which would simplify the generation of input files would be to allow defining drives with multiple inputs. A string drive with such a feature could be: "(Var1 < Var2) * Var3" and could also allow access to the value of a drive at times different from the current time (offset from the current time). These improved drive options could be implemented in a similar fashion than the node and element plugin variables.

**Category**: [User Interface](#user-interface)   
**Programming Languages**: C++, MBDyn scripting    
**Keywords**: parser, UI     
**Priority**: Low    
**Difficulty**: Intermediate    
**Mentors**: [Pierangelo Masarati](@10102934), [Marco Morandini](@marco.morandini), Louis Gagnon

-------------------------------------------------------------------------------------------

## Miscellanea
### Convergence criteria
The convergence criteria determine if the result of the simulation during a timestep has been sufficiently accurate to advance to the next timestep. Currently implemented convergence criteria could be upgraded, for example, by:
  - flagging different equations and variables according to their physical dimensions;
  - computing different residual norms for different physical domains;
  - for each physical domain, computing a reference comparison value for the convergence check;
  - developing and test the new convergence test.

**Category**: [Miscellanea](#miscellanea)   
**Programming Languages**: C++    
**Keywords**: numerical integration    
**Priority**: High    
**Difficulty**: Intermediate    
**Mentors**: [Marco Morandini](@marco.morandini), [Pierangelo Masarati](@10102934)

-------------------------------------------------------------------------------------------

### Libraries Update
The internal libraries used by MBDyn have a somewhat complicated configuration and could be simplified. Some of them could be generally useful.
Some ideas:
  - make libraries installable
  - abstract libraries in a rational manner, to reduce/eliminate cross-dependencies
  - replace standard features from broadly available libraries (e.g. STL, Boost); an example is to use shared pointers for features like constitutive laws, drive callers and so forth when useful (native C++ 2011 shared pointers are experimentally in use in some portions of code).

**Category**: [Miscellanea](#miscellanea)     
**Programming Languages**: C++, GNU autotools      
**Keywords**: libraries, software packaging   
**Priority**: High   
**Difficulty**: Intermediate    
**Mentors**: [Pierangelo Masarati](@10102934), [Marco Morandini](@marco.morandini)   

-------------------------------------------------------------------------------------------

### Package Update
A number of modifications can be made to the MBDyn package to extend its usability and improve the user experience and learning curve:
Some ideas:
  - add a comprehensive test suite
  - make it more compliant to GNU style
  - implement CMake support
  - improve Windows compatibility

**Category**: [Miscellanea](#miscellanea)   
**Programming Languages**: C++, GNU autotools, CMake scripting    
**Keywords**: software packaging, GNU autotools, software testing, GNU, cmake, Microsoft Windows
**Priority**: Medium   
**Difficulty**: Beginner    
**Mentors**: [Pierangelo Masarati](@10102934), Louis Gagnon   

-------------------------------------------------------------------------------------------

### Configuration update
The configuration of MBDyn can benefit from some improvements which would make it easier to use by a greater number of users. The tasks involved are:
  - test for `_` appended to F77 symbols (recent autoconf does that), (or better: remove all F77 code)
  - audit all the suite for not-so-GNU systems
  - add check for make

**Category**: [Miscellanea](#miscellanea)   
**Programming Languages**: C++, GNU autotools, make    
**Keywords**: software packaging, GNU autotools, make
**Priority**: Medium   
**Difficulty**: Intermediate    
**Mentors**: [Pierangelo Masarati](@10102934)

-------------------------------------------------------------------------------------------

### Cascaded analysis
Add support for cascaded models and solutions with rules for each solution to inherit the initial configuration from the final solution of another solution; the typical use case is a static solution that prepares the initial configuration for a subsequent initial value problem.

**Category**: [Miscellanea](#miscellanea)   
**Programming Languages**: C++    
**Keywords**: simulation, model restart    
**Priority**: Medium/High       
**Difficulty**: Intermediate/High    
**Mentors**: [Pierangelo Masarati](@10102934), [Marco Morandini](@marco.morandini)

-------------------------------------------------------------------------------------------

### Automatic differentiation
Automatic differentiation is a technique for generating the Jacobian matrix automatically, needed for the nonlinear solver, given the information of the residual vector which is coded in C++ language. At the moment almost all elements in MBDyn are using a hand-written contribution to the Jacobian matrix. That approach is time consuming, error prone and difficult to maintain. In addition to that, there is no easy way to check for the correctness of the Jacobian matrix.

Automatic differentiation is also implemented in MBDyn but it is still less efficient in terms of CPU time and size of generated machine code.

This project aims at perfecting the automatic differentiation implementation of MBDyn. For example by
  1. Perform different benchmarks with open source template meta programming libraries for the matrix/vector domain like blitz++, armadillo and eigen3 and also for the scalar (automatic differentiation) domain like Adolc and CppAD. Check also for recent papers in the area of automatic differentiation and template meta programming (e.g. on www.autodiff.org). Create a report about the techniques those libraries are using (aliasing, lazy evaluation, common sub-expressions, order of memory access, code generation, or just in time code generation, tape mode or tapeless mode, sparse Jacobians). That report also should show how efficient those libraries are in comparison to our current implementation in terms of CPU time and code size. In addition to that it should show how hardware techniques like vectorization (sse/avx) and intrinsic functions of compilers are exploited;
  2. Use a different approach for evaluation of template meta expressions: Instead of scalar domain template meta programs and matrix/vector domain template meta programs, only one domain should be used in order to make it easier for the compiler to optimize the code in terms of CPU time and size of machine code. Make use of the most promising techniques from 1) for the implementation. This part is probably the most difficult task because large parts of the existing library have to be redesigned. However a prototype implementation, which covers the most critical operation, namely the matrix-matrix products, is available and performs even better than equivalent Fortran code.
  The following topics should be considered:
    - not all loops should be unrolled in order to reduce code size. In that way the instruction cache should be used less intensively and stall conditions of the CPU should be avoided
    - add support for matrices and vetors with dynamic size
    - automatically create temporary results for sub-expressions which are evaluated more than once
    - implement a simple “cost model” in order to decide whether to create a temporary sub-expression or not
    - use sequential memory access whenever possible
    - use explicit vectorization for all derivatives and also for matrix vector products if the compiler supports it
    - check for aliasing at compile time or at least at runtime: if runtime checking is required, checks should be done in debug mode only
    - check if the `__restrict__` keyword of C++ can help to improve performance
    - perform benchmarks with Fortran code generated by automatic differentiation by source code transformation at different optimization levels (`-O2`, `-O3`, `-Os`)

**Category**: [Miscellanea](#miscellanea)   
**Programming Languages**: C++    
**Keywords**: simulation, automatic differentiation, code optimization    
**Priority**: Low       
**Difficulty**: Intermediate/High    
**Mentors**: Unassigned
