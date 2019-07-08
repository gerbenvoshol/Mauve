::
:: Windows build script for mauve software
:: (c)2008,2015 aaron darling
::
:: REQUIREMENTS
:: Subversion command-line, available from http://subversion.tigris.org
:: Visual Studio with 64-bit support
:: Boost compiled for 32- and 64-bit
:: NSIS -- Nullsoft scriptable install system
::
:: CONFIGURATION
:: Set paths to your subversion and ant jar
:: Set paths to your boost installation
:: Set the DEFAULTINCLUDE and DEFAULTLIB to whatever your visual studio paths are...get these from "set" on the VC command prompt
:: Set the JAVA_HOME to your JDK.  You MUST use progra~1 to avoid spaces in the path.
:: Note the first time it is run, you will need to accept SSL certificates manually
::

@ECHO ON

:: SET devenv="C:\Program Files\Microsoft Visual Studio 8\Common7\IDE\VCExpress.exe"
SET devenv="C:\Progra~2\Microsoft Visual Studio 8\"
SET svn="C:\Program Files\SlikSvn\bin\svn.exe"
SET boostinc32="C:\boost\include\boost-1_37"
SET boostinc64="C:\boost\include64\boost-1_37"
SET boostlib32="C:\boost\lib"
SET boostlib64="C:\boost\lib64"
SET JAVA_HOME="C:\Progra~1\Java\jdk1.7.0_71"
SET antbat="C:\Program Files (x86)\eclipse\plugins\org.apache.ant_1.7.0.v200803061910\bin\ant.bat"
SET puttydir="C:\Program Files (x86)\PuTTY\"
SET deploywww="darlini8@darlinglab.org:www/mauve"


%svn% co http://svn.code.sf.net/p//mauve/code/libGenome/trunk libGenome 
%svn% co http://svn.code.sf.net/p//mauve/code/muscle/trunk muscle
%svn% co http://svn.code.sf.net/p//mauve/code/libMems/trunk libMems 
%svn% co http://svn.code.sf.net/p//mauve/code/mauveAligner/trunk mauveAligner 
%svn% co http://svn.code.sf.net/p//mauve/code/sgEvolver/trunk sgEvolver 
%svn% co http://svn.code.sf.net/p//mauve/code/mauve/trunk mauve 

del errlog.txt

::
:: build 64-bit
::
call %devenv%\VC\vcvarsall.bat amd64
@ECHO ON
set INCLUDEBACKUP=%INCLUDE%
set LIBBACKUP=%LIB%
set INCLUDE="%INCLUDE%;%boostinc64%;C:\build_temp\libGenome\;C:\build_temp\libMems\;C:\build_temp\muscle\;C:\build_temp\muscle\libMuscle\"
set LIB="%LIB%;%boostlib64%;C:\build_temp\libGenome\lib;C:\build_temp\libMems\lib;C:\build_temp\libMuscle\lib"
start /wait "" %devenv%\Common7\IDE\devenv.exe mauveAligner\projects\everything.sln /rebuild "Release OpenMP|x64" /UseEnv /Out errlog.txt 
set LIB=%LIBBACKUP%


::
:: build 32-bit
::
set INCLUDE=%INCLUDEBACKUP%
call %devenv%\VC\vcvarsall.bat x86
@ECHO ON
set INCLUDE="%INCLUDE%;%boostinc32%;C:\build_temp\libGenome\;C:\build_temp\libMems\;C:\build_temp\muscle\;C:\build_temp\muscle\libMuscle\"
set LIB="%LIB%;%boostlib32%;C:\build_temp\libGenome\lib;C:\build_temp\libMems\lib;C:\build_temp\libMuscle\lib"
:: start /wait "" %devenv%Common7\IDE\devenv.exe /clean "Release OpenMP|Win32" mauveAligner\projects\everything.sln /Out errlog.txt 
start /wait "" %devenv%Common7\IDE\devenv.exe mauveAligner\projects\everything.sln /rebuild "Release OpenMP|Win32" /UseEnv /Out errlog.txt 
set INCLUDE=%INCLUDEBACKUP%
set LIB=%LIBBACKUP%


if not exist "mauveAligner\projects\Release OpenMP\mauveAligner.exe" goto postbuild
if not exist "mauveAligner\projects\Release OpenMP\progressiveMauve.exe" goto postbuild
if not exist "mauveAligner\projects\x64\Release OpenMP\mauveAligner.exe" goto postbuild
if not exist "mauveAligner\projects\x64\Release OpenMP\progressiveMauve.exe" goto postbuild

copy "mauveAligner\projects\Release OpenMP\mauveAligner.exe" mauve\win32\mauveAligner.exe
copy "mauveAligner\projects\Release OpenMP\progressiveMauve.exe" mauve\win32\progressiveMauve.exe

copy "mauveAligner\projects\x64\Release OpenMP\mauveAligner.exe" mauve\win64\mauveAligner.exe
copy "mauveAligner\projects\x64\Release OpenMP\progressiveMauve.exe" mauve\win64\progressiveMauve.exe

set MAUVE_SNAPSHOT="any_value_means_yes"
cd mauve
call %antbat% nsicompile
cd ..


:postbuild


::
:: The following gets the date and was ripped from ss64.com
::
SETLOCAL

:: This will return date into environment vars
:: Works on any NT/2K/XP machine independent of regional date settings
:: 20 March 2002

FOR /f "tokens=1-4 delims=/-. " %%G IN ('date /t') DO (call :s_fixdate %%G %%H %%I %%J)
goto :s_print_the_date
   
:s_fixdate
if "%1:~0,1%" GTR "9" shift
FOR /f "skip=1 tokens=2-4 delims=(-)" %%G IN ('echo.^|date') DO (
    set %%G=%1&set %%H=%2&set %%I=%3)
goto :eof

:s_print_the_date
ENDLOCAL&SET mm=%mm%&SET dd=%dd%&SET yy=%yy%
:eof

set SNAPDIR=%yy%\%yy%-%mm%-%dd%\windows
set SNAPDIRUNIX=%yy%/%yy%-%mm%-%dd%/windows
mkdir %yy%
mkdir %yy%\%yy%-%mm%-%dd%
mkdir %yy%\%yy%-%mm%-%dd%\windows


if not exist "mauveAligner\projects\Release OpenMP\mauveAligner.exe" goto buildFailed
if not exist "mauveAligner\projects\Release OpenMP\progressiveMauve.exe" goto buildFailed
if not exist "mauveAligner\projects\x64\Release OpenMP\mauveAligner.exe" goto buildFailed
if not exist "mauveAligner\projects\x64\Release OpenMP\progressiveMauve.exe" goto buildFailed
if not exist "mauve\dist\mauve_installer_%yy%%mm%%dd%.exe" goto buildFailed


copy mauve\dist\mauve_installer_%yy%%mm%%dd%.exe %SNAPDIR%
copy "mauveAligner\projects\Release OpenMP\*.exe" %SNAPDIR%
copy "sgEvolver\projects\Release OpenMP\sgEvolver.exe" %SNAPDIR%
copy "sgEvolver\projects\Release OpenMP\scoreAlignment2.exe" %SNAPDIR%

del windows.snaptab
echo ^<tr^> > windows.snaptab
echo ^<td^>Windows^</td^> >> windows.snaptab
echo ^<td bgcolor=#66ff66^>success^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/mauve_installer_%yy%%mm%%dd%.exe^>GUI Installer^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/mauveAligner.exe^>mauveAligner^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/progressiveMauve.exe^>progressiveMauve^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/repeatoire.exe^>repeatoire^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/^>Other programs^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/errlog.txt^>build log^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/errlog.txt^>error log^</a^>^</td^> >> windows.snaptab
echo ^</tr^> >> windows.snaptab

echo "I am the eggman"
echo "I am the walrus"

%puttydir%pscp -r %SNAPDIR% %deploywww%/snapshots/%yy%/%yy%-%mm%-%dd%

echo "coo coo "

%puttydir%pscp -r mauve\win32\Microsoft.VC80.OpenMP %deploywww%/snapshots/%yy%/%yy%-%mm%-%dd%/windows

echo " ka choo"

%puttydir%pscp windows.snaptab %deploywww%/snapshots/windows.success

goto copyEnd


:buildFailed


del windows.snaptab
echo ^<tr^> > windows.snaptab
echo ^<td^>Windows^</td^> >> windows.snaptab
echo ^<td bgcolor=#ff6666^>failure^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/^>GUI Installer^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/^>mauveAligner^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/^>progressiveMauve^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/^>repeatoire^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/^>Other programs^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/errlog.txt^>build log^</a^>^</td^> >> windows.snaptab
echo ^<td^>^<a href=%SNAPDIRUNIX%/errlog.txt^>error log^</a^>^</td^> >> windows.snaptab
echo ^</tr^> >> windows.snaptab

%puttydir%pscp -r %SNAPDIR% %deploywww%/snapshots/%yy%/%yy%-%mm%-%dd%


:copyEnd

%puttydir%pscp errlog.txt %deploywww%/snapshots/%yy%/%yy%-%mm%-%dd%/windows/

%puttydir%pscp windows.snaptab %deploywww%/snapshots/
