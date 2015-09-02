Workflow for creation of database for body waves with all waveforms etc.

1) [DONE] Make origin table
	- use makeCAdb.make_origin python tool
	- event selection criteria in this script.
	
2) [DONE] Make site, sitechan, sensor, instrument tables
	- need to first request station inventory, which will include response information
	- use make_instrument(inv,dbdir,dbnam,respdir) where inv = read_inventory(invfile)
	- use make_site(invfile)
	- use make_sitechan_sensor(invfile)
	- for the Cascadia dataset have used existing sac pole-zero files, and in fact made sitechan table with references to these files, but can write own pz files with the tool makeCAdb.write_sac_response_file(...)
	
3) Request data for each event
	- evdata_DOWNLOAD_IRIS.m 
	- produces different mat files for each station in each event directory
	- have to specify a range of orids, uses origin table and site table
	- saves in .raw structure
	
The _datinfo.mat file in each event directory keeps track of all stations, their chans, and any data processing steps.

Steps 4 and 5 can be run in either order, they are designed to be commutative :)
	
4) Rotate channels
	- evdata_ROTATE.m 
	- for the Cascadia dataset I used Helen's orientations from receiver functions.
	- where possible corrects OBS orientations so all stations have data in N,E,Z format (left-handed)
	- where no orientation is available, original data left in data.dat and data.chans - the channel names will make it clear it has not been rotated.
	- the raw data + chans are still there in data.raw
	
5) Correct instrument response
	- evdata_RMRESP.m 
	- seeks out pole-zero file for each station and uses it to remove instrument response from data channels
	- the raw data is still there in data.raw
	
4/5a) Reset
	- evdata_RESET.m 
	- if something went wrong/an error is discovered, this will reset the data and channels to the raw values from IRIS, as well as resetting the _datinfo file