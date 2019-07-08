#!/bin/bash
#
# program to build a snapshot of mauveAligner/progressiveMauve etc and post it
# to the website
#


PLATFORM=$1   # this should be something like "linux-x32"
SCRIPT=$2     # this should be something like "build_linux_x32.sh"

#
# bootstrap to ensure we've got the latest snapshot build script!
#

mkdir -p /tmp/mauve_snapshot
cd /tmp/mauve_snapshot
svn co https://svn.code.sf.net/p/mauve/code/build_scripts/
cd build_scripts
details/build_snapshot.sh $PLATFORM $SCRIPT

