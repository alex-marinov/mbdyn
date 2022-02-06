#!/bin/bash

# compile inverse dynamics embedded solver example
TGT=direct
SRCDIR="../.."
BLDDIR="../.."
CPPFLAGS="-I${BLDDIR}/include -I${SRCDIR}/include -I${SRCDIR}/libraries/libmbutil -I${SRCDIR}/libraries/libmbmath -I${SRCDIR}/libraries/libmbwrap -I${SRCDIR}/libraries/libmbc -I${SRCDIR}/mbdyn -I${SRCDIR}/mbdyn/base -I${SRCDIR}/mbdyn/hydr -I${SRCDIR}/mbdyn/aero -I${SRCDIR}/mbdyn/struct -I${SRCDIR}/mbdyn/elec -I${SRCDIR}/mbdyn/thermo"
OBJS="${BLDDIR}/mbdyn/base/.libs/libbase.a ${BLDDIR}/mbdyn/aero/.libs/libaero.a ${BLDDIR}/mbdyn/struct/.libs/libstruct.a ${BLDDIR}/mbdyn/elec/.libs/libelec.a ${BLDDIR}/mbdyn/thermo/.libs/libthermo.a ${BLDDIR}/mbdyn/hydr/.libs/libhydr.a ${BLDDIR}/mbdyn/base/.libs/libbase.a ${BLDDIR}/libraries/libmbwrap/.libs/libmbwrap.a ${BLDDIR}/libraries/libann/.libs/libmbann.a ${BLDDIR}/libraries/libmbc/.libs/libmbc_static.a"
LIBS="-L/usr/lib/x86_64-linux-gnu -lltdl -pthread -llapack -lblas -larpack -lginac -lnetcdf -lnetcdf_c++ -lsuitesparseconfig -lumfpack -lcholmod  -lm -lcamd -lm -lccolamd -lm -lcolamd -lm -lamd -lm -lklu -lamd -lcolamd -lbtf -lqrupdate -lspqr -lcholmod -lnetcdf_c++4 -lnetcdf"

g++ -g -O0 ${CPPFLAGS} -o ${TGT} ${TGT}.cc $OBJS $LIBS
