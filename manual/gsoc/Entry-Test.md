# Google Summer of Code students entry test

This entry test is divided into 4 steps: the further you get through it, the better your chances are of being selected for the GSoC.

_**Important note:** you must complete step 1 for your application to be considered.
Once you complete a step, notify your tentative mentors to obtain feedback._

_Some projects require specific entry tests, if your project does not, then follow these instructions._

## Step 1: first contact
  - download MBDyn, compile it, create a test case, and run it;
  - describe how you compiled MBDyn, including any difficulties you encountered and how you tackled them, the case you chose to run (provide a link to the actual case) and provide at least one figure showing the results as a plot of relevant data, such as positions, forces, 
velocities, etc...

Useful links: [MBDyn latest release and input manual](https://www.mbdyn.org/?Software_Download), [examples](https://www.mbdyn.org/?Documentation___Official_Documentation___Examples).

_Note that this step does not require programming, you just need to be able to write a consistent MBDyn input file_

_Should require roughly 2 hours_

## Step 2
  - implement a modification to the MBDyn code used by the test case you ran;
  - re-run your test case with the modified code;
  - describe the effect of your modification on the test case's results;
  - issue a merge request as described in the [Developers Guidelines](https://gitlab.polimi.it/Pub/mbdyn/blob/master/CONTRIBUTING.md).

_Should require roughly 4 hours_

## Step 3 create your own module
  - compile MBDyn with support for run-time loadable modules, following instructions reported in the related [FAQ](https://www.mbdyn.org/?Documentation___Official_Documentation___FAQ#HOW_CAN_I_BUILD_RUN_TIME_LOADABLE_MODULES);
  - develop a simple run-time loadable module that implements a custom instance of one of the following (in order of complexity)
    - a drive caller
    - a constitutive law
    - an element;
  - write a simple test case (an input file) that demonstrates the correct usage of your custom module;
  - document its usage and its implementation;

Some additional tips:
  - name the module `module-gsoc-yourname.cc`;
  - document the purpose and the implementation of the module, either in the [mailing list](https://gitlab.polimi.it/Pub/mbdyn/wikis/gsoc/faqs) or opening an [issue](https://gitlab.polimi.it/Pub/mbdyn/issues);
  - look at distributed modules (in the "modules/" folder) as useful working examples;
  - non-creative copy’n’paste is ** strongly discouraged**;
  - look at the mailing list archive to see past years' module related questions;

_Note that you could have to follow additional steps depending on your build target OS, read carefully the FAQs on the official website_

_The project proposed in the application does not need to be related to the module that is developed as the entry test_

_Depending on the complexity of the module, this step can require anything between one day and a few weeks of work_