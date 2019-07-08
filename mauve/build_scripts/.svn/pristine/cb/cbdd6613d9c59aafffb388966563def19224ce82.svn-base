#!/bin/bash
#
# program to build a snapshot of mauveAligner/progressiveMauve etc and post it
# to the website
#

DEPLOYHOST=darlini8@darlinglab.org

PLATFORM=$1   # this should be something like "linux-x64"
SCRIPT=$2     # this should be something like "build_linux_x64.sh"
DATE=`date +%F`
YEAR=`date +%Y`
export MAUVE_VERSION_OVERRIDE="snapshot_$DATE"

./$SCRIPT > build.log 2> build.err
SUCCESS=$?
SNAPDIR="$YEAR/$DATE/$PLATFORM"
mkdir -p $SNAPDIR
cp build.log $SNAPDIR/build.log.txt
cp build.err $SNAPDIR/build.err.txt
bzip2 $SNAPDIR/build.log.txt
bzip2 $SNAPDIR/build.err.txt

scp -r $YEAR $DEPLOYHOST:www/mauve/snapshots/
ssh $DEPLOYHOST bunzip2 -f www/mauve/snapshots/$SNAPDIR/build.log.txt.bz2
ssh $DEPLOYHOST bunzip2 -f www/mauve/snapshots/$SNAPDIR/build.err.txt.bz2

find build/mauveAligner/src -perm +111 -exec cp {} $SNAPDIR \;
cp build/mauveAligner/src/mauveStatic $SNAPDIR/mauveAligner
cp build/mauveAligner/src/progressiveMauveStatic $SNAPDIR/progressiveMauve
cp build/sgEvolver/src/scoreAlignment2 $SNAPDIR/
cp build/repeatoire/src/repeatoire $SNAPDIR/
cp build/sgEvolver/src/sgEvolver $SNAPDIR/
cp build/mauve/dist/* $SNAPDIR/
find $SNAPDIR -perm +111 -exec bzip2 {} \;
find $SNAPDIR -name "*.bz2" -exec chmod 755 {} \;

if [ -d "build/mauve/dist" ] 
then
	cd build/mauve/dist
	export COMPLETEPACKAGE=`ls *.dmg *.tar.gz`
	cd ../../..
else
	export COMPLETEPACKAGE=""
fi

scp $SNAPDIR/mauveAligner.bz2 $DEPLOYHOST:www/mauve/snapshots/$SNAPDIR/
scp $SNAPDIR/progressiveMauve.bz2 $DEPLOYHOST:www/mauve/snapshots/$SNAPDIR/
scp $SNAPDIR/repeatoire.bz2 $DEPLOYHOST:www/mauve/snapshots/$SNAPDIR/
scp -r $YEAR $DEPLOYHOST:www/mauve/snapshots/

if [ $SUCCESS -ne 0 ] 
then
	cat > $PLATFORM.snaptab <<DONE
<tr> 
	<td> 
		$PLATFORM 
	</td> 
	<td bgcolor="#ff6666">failure</td> 
	<td><a href="$SNAPDIR/$COMPLETEPACKAGE">GUI Installer</a></td>
	<td><a href="$SNAPDIR/mauveAligner.bz2">mauveAligner</a></td> 
	<td><a href="$SNAPDIR/progressiveMauve.bz2">progressiveMauve</a></td> 
	<td><a href="$SNAPDIR/repeatoire.bz2">repeatoire</a></td> 
	<td><a href="$SNAPDIR">All other programs</a></td> 
	<td><a href="$SNAPDIR/build.log.txt">build.log</a></td>
	<td><a href="$SNAPDIR/build.err.txt">build.err</a></td>
</tr> 
DONE
else
	cat > $PLATFORM.snaptab <<DONE
<tr> 
	<td> 
		$PLATFORM 
	</td> 
	<td bgcolor="#66ff66">success</td> 
	<td><a href="$SNAPDIR/$COMPLETEPACKAGE">GUI Installer</a></td>
	<td><a href="$SNAPDIR/mauveAligner.bz2">mauveAligner</a></td> 
	<td><a href="$SNAPDIR/progressiveMauve.bz2">progressiveMauve</a></td> 
	<td><a href="$SNAPDIR/repeatoire.bz2">repeatoire</a></td> 
	<td><a href="$SNAPDIR">All other programs</a></td> 
	<td><a href="$SNAPDIR/build.log.txt">build.log</a></td>
	<td><a href="$SNAPDIR/build.err.txt">build.err</a></td>
</tr> 
DONE
fi

scp $PLATFORM.snaptab $DEPLOYHOST:www/mauve/snapshots
if [ $SUCCESS -eq 0 ]; then
	scp $PLATFORM.snaptab $DEPLOYHOST:www/mauve/snapshots/$PLATFORM.success
fi

# send the error log in an e-mail if the build failed
#if [ $SUCCESS -ne 0 ]; then
#	cat > errormail.txt << DONE
#From: daemon@mauve.mooo.com
#To: aarondarling@ucdavis.edu
#Subject: mauve build failed

#Error log below.
#Output log http://darlinglab.org/mauve/snapshots/$SNAPDIR/build.err.txt
#DONE
#cat build.err >> errormail.txt
#echo "." >> errormail.txt
#sendmail -f daemon@mauve.mooo.com -t < errormail.txt
#fi

rm -rf /tmp/mauve_snapshot 
