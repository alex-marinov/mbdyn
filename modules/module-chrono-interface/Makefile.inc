# $Header$
###############################################################################
#
# Add module dependencies here; e.g., if in dir 'module-X' there's
# a file 'myfunc.cc' write:
#
# MODULE_DEPENDENCIES=myfunc.o
#
# Then put this file in directory 'module-X' as 'Makefile.inc'
#
###############################################################################

# where I generate the shared lib including C_E models, and some funcs.
# Notice!!!
# you should do at lease one of the following (refer to the notice for the usage of MBDyn shared libraries):
# 1. add LIBPATH to the 'LD_LIBRARY_PATH' environment variable during execution (worked)
# ### export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/your/libraries/path
# ### or add to bash.profile to add the path permanently
# 2. add LIBPATH to the 'LD_RUN_PATH' environment variable during linking (I didn't try this)
# 3. use the '-Wl, -rpath -Wl, LIBPATH' linker flag (I didn't try this)
# 4. have your system administrator add LIBPATH to '/etc/ld.so.conf' (what I used)
# ### sudo vim /etc/ld.so.conf;
# ### sudo ldconfig
# 5. copy your shared libraries to /lib or /usr/lib (default path)

LIBPATH=/home/runsen/Documents/GSoC/C::E_models  #/home/runsen/Documents/GSoC/C::E_models/paper_example#


# the header for the funcs “mbdyn_ce.h”.
MODULE_INCLUDE = -I .
# link the shared library, such as “libCE_models.so”
MODULE_LINK= -L$(LIBPATH)  -lCE_models  #-lMassSpring
