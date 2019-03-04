# MBDyn<sup>&copy;</sup>
**MBDyn**<sup>&copy;</sup> is a multibody analysis code.   
http://www.mbdyn.org

Copyright<sup>&copy;</sup> 1996-2019

Pierangelo Masarati     <masarati@aero.polimi.it>    
Paolo Mantegazza        <mantegazza@aero.polimi.it>

Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano    
via La Masa, 34 - 20156 Milano, Italy   
http://www.aero.polimi.it



### developers:
  - see CONTRIBUTING.md in the root directory of the source tree



### distribution:
  - set `build/version` to the right number
  - \# set `MBDYN_DEVEL=no` in `configure.ac`
  - follow the above instructions
  - run `./configure` with the desired options
  - add version and release date to `NEWS`
  - add version and release notes to `CHANGES`
  - check `BUGS`
  - advance `MBDYN_REL_ENG` tag after running `./check_cvs.sh`
  - set `MBDYN_REL_ENG_<version>` tag
  - run `make dist` \&\& `make distcheck`
  - update https://www.mbdyn.org/ web page
  - announce on mailing lists



### installation:
  - clone the source tree from this repo or download a release 
    from the [website](https://www.mbdyn.org/?Software_Download)
  - (*if needed*) run `sh bootstrap.sh` to generate configure scripts
  - configure the package; `./configure` will suffice in most cases,
    unless you need some of the supported packages and they are not 
	in the standard place
  - `make` the package
  - run `make install`; the command `mbdyn` and some utilities will
	be installed in `$PREFIX/bin`, while some libraries will be put
	in `$PREFIX/lib` (`PREFIX` defaults to `/usr/local/mbdyn`).
	Those libraries are neither required to run mbdyn nor the utilities
	(unless you force the shared build).
	A bare-bone man page will be installed in `$PREFIX/share/man/man1`.