#!/bin/sh

which autoconf > /dev/null
if [ $? -ne 0 ]; then 
  "Error: autoconf must be installed to run this script\n"; exit;
  exit 1;
fi


which pkg-config > /dev/null
if [ $? -ne 0 ]; then 
  "Error: pkg-config must be installed to run this script\n"; exit;
  exit 1;
fi


export CFLAGS="-O0"
export CXXFLAGS="-O0"
export PKG_CONFIG_PATH="`pwd`/osx_build/lib/pkgconfig"
export BOOTSTRAPFLAGS="  --with-libraries=filesystem,iostreams,program_options,system " && \
export B2FLAGS=" --toolset=darwin threading=single " && \
export CONFIGFLAGS=" --disable-shared --disable-dependency-tracking "
export PATCHBOOST=""
export PATCHARGS=""
export ANTDISTCMD="ant dmg"
mkdir -p build && \
cd build && \
../release_build.sh && \
cd .. 

