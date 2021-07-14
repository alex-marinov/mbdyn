#!/bin/bash
#
#


#host=x86_64-w64-mingw32.shared
host=x86_64-w64-mingw32.static

# mxe installation directory
mxedir=/opt/mxe

# set to where you would like to put the windows build. If the
# host is x86_64-w64-mingw32.shared the following puts it in
# something like :
# /home/jbloggs/build/mbdyn/x86_64-w64-mingw32_shared
# note dot swapped for underscore (that's what last bit does)
install_dir=$HOME/build/mbdyn-official-git/"${host//./_}"
