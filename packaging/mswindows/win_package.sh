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
replace_mxe_netcdf=false
use_netcdf=true
configure_enable_netcdf=yes
use_ginac=false
build_ginac=true
make_jobs=1
configure_enable_ginac=no

usage="$(basename "$0") [-h] [-n] -- Builds MBDyn for Microsoft Windows

where:
    -h  show this help text
    -N  disable netcdf (disabled by default)
    -n  replace netcdf makefile in MXE (replacement makefile builds netcdf 4.1.3 instead of default)
    -G  enable GiNaC
    -F  don't download and build GiNaC, assume it's already available
    -j  Number of jobs to use for make and MXE make
    -b  branch of mbdyn to checkout"

while getopts "h?nNGj:F" opt; do
    case "$opt" in
    h|\?)
        echo "$usage"
        exit 0
        ;;
    n)  replace_mxe_netcdf=true
        echo "Replace MXE netcdf makefile: $replace_mxe_netcdf"
        ;;
    N)  use_netcdf=false
        configure_enable_netcdf=no
        echo "Use NetCDF: $use_netcdf"
        ;;
    G)  use_ginac=true
        echo "Use GiNaC: $use_ginac"
        ;;
    F)  build_ginac=false
        echo "Build GiNaC: $build_ginac"
        ;;
    j)  make_jobs=$OPTARG
        echo "Make Jobs: $make_jobs"
        ;;
    esac
done
shift $((OPTIND -1))

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

if [ "${use_netcdf}" = true ]; then

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

remove_mxe_libcln_mk=false
if [ "${use_ginac}" = true ]; then
  if [ ! -f "${mxedir}/src/libcln.mk" ]; then
    cp ${thisfiledir}/libcln.mk ${mxedir}/src/
    remove_mxe_libcln_mk=true
  fi
  mxedeps="${mxedeps} libcln"
fi

# build the dependancies using MXE
cd ${mxedir}
make ${mxedeps} JOBS=${make_jobs}

if [ "${remove_mxe_libcln_mk}" = true ]; then
  rm "${mxedir}/src/libcln.mk"
fi

if [ "${use_ginac}" = true ]; then

  configure_enable_ginac=yes

  ginactmpdir=/tmp

  cd ${ginactmpdir}

  ginac_git_dir=ginac-git-4-mbdyn

  ginacroot=${ginactmpdir}/${ginac_git_dir}

  ginac_install_dir="${mxedir}/usr/${host}"  # ${ginacroot}/install


  if [ "${build_ginac}" = true ]; then

    rm -rf ${ginacroot}

    cd ${ginactmpdir}

    git clone git://www.ginac.de/ginac.git ${ginac_git_dir}

    mkdir -p ${ginac_install_dir}

    cd ${ginacroot}

    # this line fixes the Makefile.am to avoid building the docs etc (which fails)
    sed -i 's/ginac check ginsh tools doc/ginac/g' Makefile.am

    autoreconf -i

    ./configure --host=x86_64-w64-mingw32.static --enable-static --disable-shared --prefix=${ginac_install_dir}

    make -j${make_jobs}

    make install

  fi

fi


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

cppflags="-I${mxedir}/usr/x86_64-w64-mingw32.static/include/suitesparse -I${mxedir}/usr/x86_64-w64-mingw32.static/include/"
if [ "${use_ginac}" = true ]; then
  cppflags="${cppflags} -I${ginacroot}/include"
fi

linkflags="-L${mxedir}/usr/x86_64-w64-mingw32.static/lib/"
if [ "${use_ginac}" = true ]; then
  linkflags="${linkflags} -L${ginacroot}/lib"
fi

libs="-lstdc++ -lgcc -lgfortran -lamd -lcurl -lsuitesparseconfig -lcholmod -lgomp -lgnutls -lnettle -lmetis -lcolamd -lccolamd -lcamd"
if [ "${use_ginac}" = true ]; then
  libs="${libs} -lginac"
fi
if [ "${use_netcdf}" = true ]; then
  libs="${libs} -lnetcdf-cxx4 $(/opt/mxe/usr/x86_64-w64-mingw32.static/bin/nc-config  --libs) -lportablexdr -ljpeg $(/opt/mxe/usr/x86_64-w64-mingw32.static/bin/curl-config  --static-libs)"
  #sed -i 's/netcdf_c++4/netcdf-cxx4/g' configure.ac
fi

sed -i 's/-Wno-error=format-truncation//g' configure.ac

./bootstrap.sh

# _USE_MATH_DEFINES required for windows to define M_PI etc see:
#
# https://stackoverflow.com/questions/26065359/m-pi-flagged-as-undeclared-identifier
#
./configure --host=$host \
            `# --enable-Werror=no` \
            --with-boost \
            --prefix=${install_dir} \
            --enable-install_test_progs=yes \
            --with-umfpack \
            --with-metis=no \
            --enable-netcdf=${configure_enable_netcdf} \
            --with-ginac=${configure_enable_ginac} \
            `# --with-module=fabricate` \
            LINKFLAGS="$LINKFLAGS ${linkflags}" \
            CXXFLAGS="$CXXFLAGS -std=c++11 -D_USE_MATH_DEFINES" \
            LIBS="$LIBS ${libs}" \
            CPPFLAGS="$CPPFLAGS ${cppflags}"

if hash wine 2>/dev/null; then
    if hash update-binfmts 2>/dev/null; then
        sudo update-binfmts --enable wine
    fi
fi

make clean

make -j${make_jobs}

make install

${thisfiledir}/copy_win_libs.sh

# copy test files
cp -r ${thisfiledir}/test ${install_dir}/share/

# restore working directory
cd ${restorepwd}



