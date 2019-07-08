#!/bin/sh

BPATH=`pwd`
export PATH=$PATH:$BPATH/bin
export PKG_CONFIG_PATH=$BPATH/lib/pkgconfig

if [ -z "$SVNPROTO" ]; then SVNPROTO="http"
fi

curl -LO http://downloads.sourceforge.net/project/boost/boost/1.57.0/boost_1_57_0.tar.bz2 2> /dev/null && \
svn co $SVNPROTO://svn.code.sf.net/p/mauve/code/libGenome/trunk libGenome && \
svn co $SVNPROTO://svn.code.sf.net/p/mauve/code/libMems/trunk libMems && \
svn co $SVNPROTO://svn.code.sf.net/p/mauve/code/muscle/trunk muscle && \
svn co $SVNPROTO://svn.code.sf.net/p/mauve/code/sgEvolver/trunk sgEvolver && \
svn co $SVNPROTO://svn.code.sf.net/p/mauve/code/repeatoire/trunk repeatoire && \
svn co $SVNPROTO://svn.code.sf.net/p/mauve/code/mauveAligner/trunk mauveAligner || \
exit 21

tar xjf boost_1_57_0.tar.bz2 && \
cd boost_1_57_0 && \
$PATCHBOOST $PATCHARGS && \
./bootstrap.sh --with-libraries=filesystem,iostreams,program_options,system --prefix=$BPATH && \
./b2 -a -j4 $B2FLAGS release link=static  install  || \
exit 22
cd ..


cd libGenome
./autogen.sh && ./configure $CONFIGFLAGS --prefix=$BPATH && make clean && make -j2 install && make dist || \
exit 23
cd ..


cd muscle
./autogen.sh && ./configure $CONFIGFLAGS --prefix=$BPATH && make clean && make -j2 ; make && make install && make dist || \
exit 24
cd ..

cd libMems
./autogen.sh && ./configure $CONFIGFLAGS --prefix=$BPATH --with-boost=$BPATH && make clean && make -j2 install && make dist || \
exit 25
cd ..

cd mauveAligner
./autogen.sh && ./configure $CONFIGFLAGS --prefix=$BPATH && make clean && make && find src -maxdepth 1 -perm +111 -exec strip {} \; && make dist || \
exit 26
cd ..


svn co $SVNPROTO://svn.code.sf.net/p/mauve/code/mauve/trunk mauve && \
cd mauve && \
if [ -n "$MAUVE_VERSION_OVERRIDE" ]; then
  perl -p -i -e "s/property.+name\=\"release\.version\".+value\=.+\"/property name\=\"release\.version\" value\=\"$MAUVE_VERSION_OVERRIDE\"/g" build.xml
fi

cp ../mauveAligner/src/mauveStatic linux-x64/mauveAligner && \
cp ../mauveAligner/src/mauveStatic osx/mauveAligner && \
cp ../mauveAligner/src/progressiveMauveStatic linux-x64/progressiveMauve && \
cp ../mauveAligner/src/progressiveMauveStatic osx/progressiveMauve && \
$ANTDISTCMD 
cd ..

cd sgEvolver
./autogen.sh && ./configure $CONFIGFLAGS --prefix=$BPATH && make clean && make && find src -maxdepth 1 -perm +111 -exec strip {} \; && make dist || \
exit 27
cd ..

cd repeatoire
./autogen.sh && ./configure $CONFIGFLAGS --prefix=$BPATH && make clean && make && find src -maxdepth 1 -perm +111 -exec strip {} \; && make dist || \
exit 28
cd ..

