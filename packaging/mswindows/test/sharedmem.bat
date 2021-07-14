rem Set to MBDyn executable path
rem set BLDPATH=/home/masarati/Lavoro/mbdyn/mbdyn-devel
rem PATH="$BLDPATH/mbdyn:$BLDPATH/utils:$PATH"

rem expect this batch file to be in the same directory as sharedmem.mbd
rem see https://stackoverflow.com/questions/17063947/get-current-batchfile-directory
rem for more information on the commnd below (uses batch parameters)
set thisfilepath=%~dp0

set mbdyndir=G:\mbdyn\windows\bin

set outputdir=%thisfilepath%\sharedmem_output_windows
rem create a directory for the output
if not exist "%outputdir%" mkdir "%outputdir%"

echo thisfilepath is %thisfilepath%
echo mbdyndir is %mbdyndir%

echo "*** Starting MBDyn... "
rem note character escaping required using ^ because we're using 'start'. See
rem https://superuser.com/questions/338277/windows-cmd-batch-start-and-output-redirection
rem for more details
start call "%mbdyndir%\mbdyn.exe" -f "%thisfilepath%\sharedmem.mbd" -o "%outputdir%\sharedmem_output" ^> "%outputdir%\sharedmem_output.txt"
rem "mbdyn running %thisfilepath%\sharedmem.mbd"
rem ^2^>^&^1
rem "G:\mbdyn\windows\bin\mbdyn.exe" -PPP -f "sharedmem.mbd" -o "sharedmem_output" > "sharedmem_output.txt" 2>&1
echo "sleeping 2 seconds..."
timeout 2

echo "*** Starting test..."

rem usage: testsocket [options]
	rem -a		use accelerations
	rem -c [random:]<c>	number of iterations
	rem -f {fx,fy,fz,mx,my,mz} reference node force/moment
	rem -H <url>	URL (local://path | inet://host:port)
	rem -l		labels
	rem -i <filename>	input file
	rem -n		only forces, no moments
	rem -N <nodes>	nodes number
	rem -o <filename>	 output file
	rem -p {f0x,f0y,f0z,m0x,m0y,m0z,...}	nodal forces (need -N first)
	rem -r		use reference node data
	rem -R {mat|theta|euler123}	orientation format
	rem -s <sleeptime>	sleep time between tries
	rem -t <timeout>	how long to wait for connection
	rem -v		verbose
	rem -x		data_and_next

rem "G:\mbdyn\windows\bin\test_strext_sharedmemory_lib_cxx.exe" -c 3 -N 2 -p 0.,0.,1.,0.,0.,0.,0.,0.,1.,0.,0.,0.  -R mat -s 0 -x -H shared_mem_test_shmName -v
"%mbdyndir%\test_strext_sharedmemory.exe" ^
	-c 3 ^
	-N 2 ^
	-p 0.,0.,1.,0.,0.,0.,0.,0.,1.,0.,0.,0. ^
	-R mat ^
	-s 0 ^
	-x ^
	-H shared_mem_test_shmName

echo "*** Done"


