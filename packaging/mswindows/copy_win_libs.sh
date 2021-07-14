#!/bin/bash
#
#
# copy_win_libs.sh
#

set -ex

# get this file's containing directory
# see https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
# for explanation
thisfiledir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# get config variables: mxedir, install_dir and host
source ${thisfiledir}/win_package_config.sh

srcdir=${thisfiledir}/../..

${mxedir}/tools/copydlldeps.sh -c -d ${install_dir}/bin -R $HOME/src/mbdyn/mbdyn/    \
--indir ${srcdir}/mbdyn/    \
--indir ${srcdir}/base/    \
--indir ${srcdir}/utils/    \
--indir ${srcdir}/libraries/libmbwrap/    \
--indir ${srcdir}/libraries/liby12    \
--indir ${srcdir}/libraries/libmbutil    \
--indir ${srcdir}/libraries/libmbmath

case "$host" in
  *static*)
    # don't copy dlls
    echo "Not copying dll's as it is a static build"
    ;;
  *shared*)
    # copy over other shared libs required by executeable
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libgcc_s_seh-1.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libgfortran-3.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libklu.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/liblapack.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libltdl-7.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libopenblas.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libstdc++-6.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libumfpack.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libwinpthread-1.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libdl.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libamd.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libbtf.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libcolamd.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libsuitesparseconfig.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libcholmod.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libnetcdf*.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libblas.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libcurl-4.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libhdf5-8.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libhdf5_hl-8.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libmfhdf-0.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libportablexdr-0.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libgomp-1.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libgnutls-30.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libidn2-0.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libssh2-1.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libnettle-6.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/zlib1.dll ${install_dir}/bin/

    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libcamd.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libccolamd.dll ${install_dir}/bin/
    cp ${mxedir}/usr/x86_64-w64-mingw32.shared/bin/libmetis.dll ${install_dir}/bin/

    # copy socket library to lib dir
    #cp ${mxedir}/usr/x86_64-w64-mingw32.static/lib/libws2_32.a ${install_dir}/lib/
    ;;
  *)
    ;;
esac





