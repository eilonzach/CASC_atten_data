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
	- produces different mat files for each station in each event directory
	- use matlab get_eventsdata_IRIS.m 
	- have to specify a range of orids, uses origin table and site table
	- saves in .raw structure
	
4) Rotate channels
	- use matlab rotate_eventsdata.m
	- for the Cascadia dataset I used Helen's orientations from receiver functions.
	- saves raw channels in data.raw_prerot, but where possible corrects OBS orientations so all stations have data in N,E,Z format (left-handed)
	- where no orientation is available, original data left in data.dat and data.chans - the channel names will make it clear it has not been rotated.
	
5) Correct instrument response
	- use matlab rmresp_eventsdata.m
	- seeks out pole-zero file for each station and uses it to remove instrument response from data channels
	- the raw data is still there in data.raw_prerot