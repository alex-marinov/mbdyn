(This is the final report for GSoC2020)
# Final report of a MBDyn GSoC2020 project: A user-defined run-time loadable module template for co-simulation with Chrono::Engine
by Runsen ZHANG, August, 2020
## Project overview
The goal of the project is to enable tight coupling co-simulation between [MBDyn](https://www.mbdyn.org) and another open-source software, [Chrono::Engine](http://projectchrono.org/) (C::E). For that, a user-defined runtime loadable module template is provided. The template includes a MBDyn-C::E interface runtime module, by which users can built their own C::E models in the opaque object, and link it to the MBDyn’s module. Several demos are aslo provided as references.

## What I have done
**Codes:**([contributions within GSoC2020](https://public.gitlab.polimi.it/DAER/mbdyn/-/tree/gsoc_chrono_interface/modules/module-chrono-interface))
- Created an interface for single coupling node case;
- Created the run-time loadable module named as "module-chrono-interface" for MBDyn-side codes;
- Finished codes for C::E-side codes;
- Extended the codes to cases for arbitrary number of nodes;

**Demos:**([demos for module-chrono-interface](https://public.gitlab.polimi.it/DAER/mbdyn/-/tree/gsoc_chrono_interface/modules/module-chrono-interface/examples))
- Provided three demos, "two-mass oscillator" case, "pendulum" case, and "multi-node" case, to test the module.

## Usage
To use this module, users need to preinstall Chrono::Engine and the develop version of MBDyn. After installation, Chrono::Engine models defined by users should be built with the interface for C::E-side, and be dynamically linked to MBDyn. The syntax in MBDyn is:

> User defined:
>
> Element_Number, ChronoInterface,
>
> coupling, \[none | tight\],
>
> iteration_number, \[tolerance, <real>\],
>
> motor type, \[setpoint | spline\],
>
> \[length scale, Real, mass scale, Real, \]
>
> nodes number, \<Int\>,
>
> \<nodes_ID\>, \[offset, \<Vec3\>\], \<CEBodies_ID\>, output, \[yes | no\]
>
> \[ground, \<C::E_ground_ID\>, \<Vec3\>, \<Mat3x3\>\],
>
> verbose, \[yes | no\];

## My strace for GSoC 2020
- [Weekly reports](https://public.gitlab.polimi.it/DAER/mbdyn/-/wikis/Google-Summer-of-Code/Students-Blogs)
- [Reports for the first and second evaluation](https://docs.google.com/document/d/13tWTA6T8iw5FbyF6oxnNy7YXeSlEV_LGjM01L2IpAqg/edit?usp=sharing)
- [Contributions with GSoC 2020](https://public.gitlab.polimi.it/DAER/mbdyn/-/tree/gsoc_chrono_interface/modules/module-chrono-interface)
- [List of commits](https://public.gitlab.polimi.it/DAER/mbdyn/-/commits/gsoc_chrono_interface)

## Future study
Although goals of GSoC have been achieved, algorithms for co-simulation between MBDyn and Chrono::Engine are still worth of further in-depth study, such as:
- Loose/explicit co-simulation scheme;
- Multirate co-simulation scheme;
- ...
