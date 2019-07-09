#!/bin/bash
#
# script to cross-build and pack up windows build output requires
# M Cross Environment (mxe.cc). The location of the mxe installation 
# (e.g. /opt/mxe) must be set in the "mxedir" variable in 
# win_package_config.sh
#

set -ex

# A POSIX variable
OPTIND=1  # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
replace_mxe_netcdf=true
use_netcdf=false
configure_enable_netcdf=no

usage="$(basename "$0") [-h] [-n] -- creates Renewnet Foundry release

where:
    -h  show this help text
    -N  enable netcdf (disabled by default)
    -n  don't replace netcdf makefile in MXE (replacement makefile builds netcdf 4.1.3 instead of default)"

while getopts "h?nN" opt; do
    case "$opt" in
    h|\?)
        echo "$usage"
        exit 0
        ;;
    n)  replace_mxe_netcdf=false
        echo "Replace MXE netcdf makefile: $replace_mxe_netcdf"
        ;;
    N)  use_netcdf=true
        configure_enable_netcdf=yes
        echo "Replace MXE netcdf makefile: $replace_mxe_netcdf"
        ;;
    esac
done

#host=x86_64-w64-mingw32.shared

# store current directory so we can restore it when we're done
restorepwd=$(pwd)

# get this file's containing directory
# see https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
# for explanation
thisfiledir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# get config variables: mxedir, install_dir and host
source ${thisfiledir}/win_package_config.sh

echo host contains: ${host}

# make sure all the required dependancies are in MXE
mxedeps="boost libltdl suitesparse lapack hdf5 hdf4"

if [ "${use_netcdf}" = true]; then

  # add netcdf to the mxe build
  mxedeps="${mxedeps} netcdf"

  if [ "${replace_mxe_netcdf}" = true ]; then
    # test if netcdf files are the same
    cmp --silent ${thisfiledir}/netcdf.mk ${mxedir}/src/netcdf.mk && cmpstatus=0

    #cmpstatus=$?

    if [ ${cmpstatus} -eq 0 ]; then
      echo "Not copying netcdf files as netcdf.mk is the same."
    else
      # copy makefile for netcdf 4.1.3 and associated patch file into mxe, will
      # trigger a recompilation of netcdf
      cp ${thisfiledir}/netcdf* ${mxedir}/src/
    fi

  else

    echo "Not replacing MXE netcdf makefile."

    # we are using the latest netcdf which has a separate C++ library
    mxedeps="${mxedeps} netcdf-cxx4"

  fi

fi

# build the dependancies using MXE
cd ${mxedir}
make ${mxedeps}

srcdir=${thisfiledir}/../..

cd ${srcdir}

mkdir -p ${install_dir}

# clear out any previous installation
rm -rf ${install_dir}/*

./bootstrap.sh

# The following prevents wine from automatically running .exe files
#
# https://askubuntu.com/questions/344088/how-to-ensure-wine-does-not-auto-run-exe-files
if hash wine 2>/dev/null; then
    if hash update-binfmts 2>/dev/null; then
        sudo update-binfmts --disable wine
    fi
fi

# _USE_MATH_DEFINES required for windows to define M_PI etc see:
#
# https://stackoverflow.com/questions/26065359/m-pi-flagged-as-undeclared-identifier
#
./configure --host=$host \
            --with-boost \
            --prefix=${install_dir} \
            --enable-install_test_progs=yes \
            --with-umfpack \
            --with-metis=no \
            --enable-netcdf=${configure_enable_netcdf} \
            `# --with-module=fabricate` \
            LINKFLAGS="$LINKFLAGS -L${mxedir}/usr/x86_64-w64-mingw32.static/lib/" \
            CXXFLAGS="$CXXFLAGS -std=c++11 -D_USE_MATH_DEFINES" \
            LIBS="$LIBS -lstdc++ -lgcc -lgfortran -lamd -lcurl -lsuitesparseconfig -lcholmod  -lportablexdr -lws2_32 -lnetcdf -lhdf5_hl -lhdf5 -lmfhdf -ldf -lz -lm -ljpeg  -lcurl -lgomp -lgnutls -lidn2 -lssh2 -lnettle -lmetis -lcolamd -lccolamd -lcamd " \
            CPPFLAGS="$CPPFLAGS -I${mxedir}/usr/x86_64-w64-mingw32.static/include/suitesparse -I${mxedir}/usr/x86_64-w64-mingw32.static/include/"

if hash wine 2>/dev/null; then
    if hash update-binfmts 2>/dev/null; then
        sudo update-binfmts --enable wine
    fi
fi

make -j3

make install

${thisfiledir}/copy_win_libs.sh

# copy test files
cp -r ${thisfiledir}/test ${install_dir}/share/

# restore working directory
cd ${restorepwd}



