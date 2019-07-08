#!/bin/sh
#
# script to create a mauve source release
# this should be run from within a development dir that contains 
# the appropriate source code working copies

VERSION="2.2.0"
VERSIONDASH="2-2-0"


# step 1.  Create source distributions for all except mauve java 
mkdir mauve_$VERSION
cd libGenome && make dist && cp *.tar.gz ../mauve_$VERSION && cd .. && \
cd muscle && make dist && cp *.tar.gz ../mauve_$VERSION && cd .. && \
cd libMems && make dist && cp *.tar.gz ../mauve_$VERSION && cd .. && \
cd mauveAligner && make dist && cp *.tar.gz ../mauve_$VERSION && cd .. && \
# step 2.  Create mauve java source dist
echo "\n\nCreating mauve java source distribution\n\n" && \
cd mauve && svn export . ../mauve-$VERSION && cd .. \
tar czf mauve-$VERSION.tar.gz mauve-$VERSION && \
mv mauve-$VERSION.tar.gz mauve_$VERSION && \
rm -rf mauve-$VERSION && \
# step 3.  Tag the subversion repositories
echo "\n\nTagging repositories...\n\n" && \
cd libGenome && svn copy . https://mauve.svn.sourceforge.net/svnroot/mauve/libGenome/tags/mauve-$VERSIONDASH-release -m "Tagging mauve $VERSION source release" && cd .. && \
cd muscle && svn copy . https://mauve.svn.sourceforge.net/svnroot/mauve/muscle/tags/mauve-$VERSIONDASH-release -m "Tagging mauve $VERSION source release" && cd .. && \
cd libMems && svn copy . https://mauve.svn.sourceforge.net/svnroot/mauve/libMems/tags/mauve-$VERSIONDASH-release -m "Tagging mauve $VERSION source release" && cd .. && \
cd mauveAligner && svn copy . https://mauve.svn.sourceforge.net/svnroot/mauve/mauveAligner/tags/mauve-$VERSIONDASH-release -m "Tagging mauve $VERSION source release" && cd .. && \
 cd mauve && svn copy . https://mauve.svn.sourceforge.net/svnroot/mauve/mauve/tags/mauve-$VERSIONDASH-release -m "Tagging mauve $VERSION source release" && cd .. && \
# step 4. Copy source release to server
echo "\n\nCopying source release to web server\n\n"&& \
scp -r mauve_$VERSION koadman@gel.ahabs.wisc.edu:/srv/www/vhosts/gel/mauve/source

