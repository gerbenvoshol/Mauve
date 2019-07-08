#!/bin/sh
which pkg-config > /dev/null
if [ $? -ne 0 ]; then 
  "Error: pkg-config must be installed to run this script\n"; exit;
  exit 1;
fi


mkdir -p build && \
cd build && \
export NO_COMPRESSION="true" && \
export ANTDISTCMD="ant dist" && \
export SVNPROTO="svn" && \
../release_build.sh && \
cd .. && \
echo "done building binaries!\n"

