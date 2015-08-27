#!/usr/bin/env python
"""This module contains the various subroutines
Usage import """
#######################################################################################

#####################  IMPORT STANDARD MODULES   ######################################   

import sys,os, site, signal, shutil
import numpy as np #for numerical analysis
from scipy import integrate  #integration routine
from subprocess import call #for calling unix commands
from datetime import date  #to give a timestamp to output
import time
import pdb	#for the debugger pdb.set_trace()
from matplotlib import pyplot as plt
from math import pi,exp,log,sqrt
from ConfigParser import SafeConfigParser
from obspy import UTCDateTime


####################       IMPORT OBSPY MODULES     ######################################

#  Initialise client object
# from obspy.arclink.client import Client
from obspy.fdsn import Client
client = Client("IRIS")
from obspy.fdsn.header import URL_MAPPINGS

print "Available Clients are:"
for key in sorted(URL_MAPPINGS.keys()):
	print("{0:<7} {1}".format(key,  URL_MAPPINGS[key]))


# setup ObsPy-Antelope
signal.signal(signal.SIGINT, signal.SIG_DFL)
sys.path.append(os.environ['ANTELOPE'] + "/data/python")
site.addsitedir('/opt/antelope/local/lib/python' + sys.version[:3])
# import antelope.datascope as datascope
from antelope.datascope import dbtmp, dbopen, destroying, closing
from antelope.datascope import DbfindEnd, DbaddvError

####################     IMPORT OWN MODULES     ######################################

from myFUNCTIONS import tup2str

#####################  MODULES   ######################################################   


def make_origin():
	"""This subroutine contains the workflow to select a bunch of event parameters and then use the IRIS request tools to build a catalogue of events that satisfy the parameters. Then the catalogue is exported to the Antelope origin table format."""
	print "\n\n THE EVENT PARAMETERS ARE DEFINED IN THE MAKE_ORIGIN SUBROUTINE. EDIT THERE\n\n"
# 	DEFINE EVENT PARAMETERS (or leave empty for unconstrained)
	dbdir = "/Users/zeilon/Work/CASCADIA/CAdb/" # needs final slash
	dbnam= "cascattendb"

	deplims = [    0 ,1000 ] # [mindep, maxdep] ... in km

	maglims = [  6.5 , 10 ]  # [minmag, maxmag]

	starttime = "2006-01-01" # "yyyy-mm-dd"
	endtime   = "2010-12-31" # "yyyy-mm-dd"

	# geographic search params - leave either blank to ignore. Will fail if both filled in.

	# 						SEARCH IN BOX
	latlims = [   ] # [minlat, maxlat]
	lonlims = [   ] # [minlon, maxlon]
	# 	   --------------------------OR--------------------------
	# 					SEARCH BY DISTANCE
	refpt  = [45, -120] # [ref_lat, ref_lon] 
	gcarclims = [ 30, 135 ]  # [min_gcarc, max_gcarc] ... in degrees
#  ================== NO NEED TO ALTER ANYTHING BELOW THIS LINE ================== # 
			
# 	parse conditions to right format
	starttime =  UTCDateTime(starttime)
	endtime = UTCDateTime(endtime)
	
	magmin = np.float(maglims[0])
	magmax = np.float(maglims[1])
	depmin = np.float(deplims[0])
	depmax = np.float(deplims[1])

	# if distance search prefs empty, must be a BOX search
	if (not refpt) | (not gcarclims): 
		print "No distance limits"
		refla = None
		reflo = None
		gcmin = None
		gcmax = None
	else:    
		refla = np.float(refpt[0])
		reflo = np.float(refpt[1])
		gcmin = np.float(gcarclims[0])
		gcmax = np.float(gcarclims[1])

	# if box search prefs empty, must be a DISTANCE search
	if (not latlims) | (not lonlims):
		print "No box limits"
		latmin = None
		latmax = None
		lonmin = None
		lonmax = None
	else:    
		latmin = np.float(latlims[0])
		latmax = np.float(latlims[1])
		lonmin = np.float(lonlims[0])
		lonmax = np.float(lonlims[1])

	# Get the catalogue
	cat = client.get_events(starttime=starttime,endtime=endtime,\
		minmagnitude=magmin,maxmagnitude=magmax,\
		mindepth=depmin,maxdepth=depmax,\
		minlatitude=latmin,maxlatitude=latmax,\
		minlongitude=lonmin,maxlongitude=lonmax,\
		latitude=refla,longitude=reflo,minradius=gcmin,maxradius=gcmax,\
		orderby="time-asc",catalog="GCMT")
	print(cat) 

	nevts = len(cat)
	
# 	REMOVE ANY OVERLAPPING EVENTS
# 	Remove any events from the catalog that are < 2000s after or 1000s before another event

	evtimes = np.zeros([nevts,1])

	for ie in range(0,nevts):
		evt =  cat.events[ie]
		elat = evt.origins[0].latitude
		elon = evt.origins[0].longitude
		edep = evt.origins[0].depth/1000
		evtime = evt.origins[0].time
		emag = evt.magnitudes[0].mag
		emagt = evt.magnitudes[0].magnitude_type
		print "%s=%.2f event at [%6.2f, %7.2f] and %5.1f km depth on %s " % \
			(emagt, emag, elat, elon, edep, evtime)
		evtimes[ie] = evtime
	kill = []
	for ie in range(0,nevts):
		if ( ( (evtimes[ie]-evtimes)<2000 ) & ( (evtimes[ie]-evtimes)>0 ) ).any():
			print("One recently before:")
			print(UTCDateTime(evtimes[ie]))
			print(UTCDateTime(evtimes[ie+1]))
			kill.append(ie)
		elif ( ( (evtimes[ie]-evtimes)>-1000 ) & ( (evtimes[ie]-evtimes)<0 ) ).any():
			print("One soon after:")
			print(UTCDateTime(evtimes[ie]))
			print(UTCDateTime(evtimes[ie-1]))
			kill.append(ie)
		
	evts = cat.events[:]
	for ik in range(0,len(kill)):
		evt_to_kill = evts[kill[ik]]
		print "Removing event %s" % (UTCDateTime(evt_to_kill.origins[0].time))
		cat.events.remove(evt_to_kill)

	print "%.0f events removed becauss of overlapping data... plotting" % nevts-len(cat)
	nevts = len(cat)
# 	cat.plot();
	
	ifwrite = raw_input("Write %.0f events to origin table? [y/n]" % nevts)
	if ifwrite=="y":
	# 	What if origin table already exists
		origintabfile = dbdir+dbnam+".origin"
		if os.path.isfile(origintabfile):
			kok = 0;
			while kok==0:
				yn = raw_input("Origin table already seems to exist. Overwrite/append/cancel? [o/a/c] ")
				if   yn=="o":
					print("Ok,deleting and overwriting...")
					os.remove(origintabfile)
					kok=1;
				elif yn=="a":
					print("Ok, appending...")
					kok=1;
				elif yn=="c":
					raise SystemExit("Ok, not deleting... aborting")
					exit
					kok=1;

	# 	Add rows to origin table
		with closing(dbopen(dbdir+dbnam,'r+')) as db:
			dbor = db.lookup(table = 'origin')
			nrecs = dbor.record_count
			dbor.record=dbor.record_count # set current record to last row to append
			for ie in range(0,nevts):
				evt =  cat.events[ie]
				elat = evt.origins[0].latitude
				elon = evt.origins[0].longitude
				edep = evt.origins[0].depth/1000
				evtime=evt.origins[0].time
				emag = evt.magnitudes[0].mag
				emagt = evt.magnitudes[0].magnitude_type
				print "Event",nrecs+ie+1,"  ",elat,elon,edep,emag,utc2epoch(evtime)
				flds = ('orid','lat','lon','depth','time','ms','nass','ndef' )
				vals = ( nrecs+ie+1,elat,elon,edep,utc2epoch(evtime),emag,0,0)
				fldvals = zip(flds, vals)
				try:
					dbaddrec = dbor.addv(*fldvals)		
				except DbaddvError:
					print "Error with adding this row..."
					try:
						ifrep = dbor.find('orid=='+nrecs+ie+1)
						print "Skipping repeated orid %s" % sta
					except DbfindEnd:
						pdb.set_trace()
						raise SystemExit("Something wrong - won't add orid,but not already in there")
		
	return 

def request_inv():
	"""This subroutine contains the workflow to define a bunch of sta+chan parameters and then use the IRIS request tools to obtain all stations satisfying these parameters in an obspy inventory."""
	print "\n\n THE STATION + CHANNEL PARAMETERS ARE DEFINED IN THE MAKE_INV SUBROUTINE. EDIT THERE\n\n"

# 	DEFINE SITE PARAMETERS (or leave empty for unconstrained)	
	net_dos = "_CASCADIA" # or '_CASCADIA_OBS' for only obs, or '*' for all
	
	chan_dos = "B*,H*"
	
#  ================== NO NEED TO ALTER ANYTHING BELOW THIS LINE ================== # 
			
	# Get the inventory of stations
	inv = client.get_stations(network=net_dos,channel=chan_dos,level="response")

	print "Saving IRIS inventory with response files as Casc_inventory.xml" 
	inv.write(dbdir+'Casc_inventory.xml', format='STATIONXML')

	return inv
	
def make_site(invfile):
	"""This subroutine reads an obspy inventory and creates site table."""
	
	dbdir = "/Users/zeilon/Work/CASCADIA/CAdb/" # needs final slash
	dbnam= "cascattendb"

	
	print("Reading inventory file %s") % invfile
	from obspy import read_inventory
	inv = read_inventory(invfile)
	
	print "Writing .site table"
	flds_site = ('sta','ondate','offdate','lat','lon','elev','staname','statype','refsta')
	# will store network in the 'refsta' field.	 
			
	ifwrite = raw_input("Write %.0f stations to site table? [y/n] " % count_stas(inv))
	if ifwrite=="y":
		# 	What if site table already exists
		sitetabfile = dbdir+dbnam+".site"
		if os.path.isfile(sitetabfile):
			kok = 0;
			while kok==0:
				yn = raw_input("Site table already seems to exist. Overwrite/append/cancel? [o/a/c] ")
				if yn=="o":
					print("Ok,deleting and overwriting...")
					os.remove(sitetabfile)
					kok=1;
				elif yn=="a":
					print("Ok, appending...")
					kok=1;
				elif yn=="c":
					raise SystemExit("Ok, not deleting... aborting")
					exit
					kok=1;

	
		with closing(dbopen(dbdir+dbnam,'r+')) as db:
			dbsi = db.lookup(table = 'site')
	
			ista = 0
			for inn in range(0,len(inv.networks)):
				network = inv.networks[inn].code.encode()
				print "\n ====Network %s" % network
				for iss in range(0,len(inv.networks[inn].stations)):
					
					station = inv.networks[inn].stations[iss]
					sta = station.code.encode()
					on_date = station.start_date.year*1000 + station.start_date.julday
					off_date = station.end_date.year*1000 + station.end_date.julday
					
					if network == '7D':
						statype = 'OBS'
					else:
						statype = 'LAND'
					
					vals_site = (sta,on_date,off_date,\
							station.latitude,station.longitude,station.elevation,\
							station.site.name,statype,network)
					fldvals_site = zip(flds_site, vals_site)
					print "--------Station %5s %5.2f %7.2f  --  %s"\
									% (sta,station.latitude,station.longitude, station.site.name)
					try:
						dbaddrec = dbsi.addv(*fldvals_site)		
					except DbaddvError:
						print "Error with adding this row..."
						try:
							print "Maybe OBS/LAND station repeat?"
							ifrep = dbsi.find('sta=="'+sta+'"')
							dbsi.record=ifrep
							statype_existing = tup2str(dbsi.getv('statype'))
							if statype_existing!=statype: # same staname, different statype
								print("Different statypes: ammending old sta %s to %s"\
										 % (sta,sta+statype_existing[0]))
								# amend old row
								dbsi.putv(('sta',sta+statype_existing[0]))
								# add new row with amended sta name
								print("Ammending new sta %s to %s"\
										 % (sta,sta+statype[0]))
								vals_site = (sta+statype[0],on_date,off_date,\
									station.latitude,station.longitude,\
									station.elevation,station.site.name,statype,network)
								fldvals_site = zip(flds_site, vals_site)
								try:
									dbaddrec = dbsi.addv(*fldvals_site)	
								except DbaddvError:
									print("Some other error")
									pdb.set_trace()
									raise SystemExit("Something wrong - won't add sta,\
							 						but not already in there")
							else: # same staname and statype: skip
								print "Skipping repeated station %s" % sta
						except DbfindEnd:
							pdb.set_trace()
							raise SystemExit("Something wrong - won't add sta,\
							 						but not already in there")
							 
					
					
	return

def make_sitechan_sensor(invfile):
	"""This subroutine reads an obspy inventory and creates sitechan and sensor tables."""
	
	dbdir = "/Users/zeilon/Work/CASCADIA/CAdb/" # needs final slash
	dbnam= "cascattendb"
	respdir = "/Users/zeilon/Work/CASCADIA/CAdb/response/" # needs final slash
	datemin = UTCDateTime(2011,1,1)

	print("Reading inventory file %s") % invfile
	from obspy import read_inventory
	inv = read_inventory(invfile)

	print "Writing .sitechan + .sensor tables"
	flds_sch = ('sta','chan','ondate','offdate','chanid','hang','vang','descrip')
	flds_sen = ('sta','chan','time','endtime','inid','chanid')

	# What if sitechan table already exists?
	sitechantabfile = dbdir+dbnam+".sitechan"
	if os.path.isfile(sitechantabfile):
		kok = 0;
		while kok==0:
			yn = raw_input("sitechan table already seems to exist. Overwrite/append/cancel? [o/a/c] ")
			if yn=="o":
				print("Ok,deleting and overwriting...")
				os.remove(sitechantabfile)
				kok=1;
			elif yn=="a":
				print("Ok, appending...")
				kok=1;
			elif yn=="c":
				raise SystemExit("Ok, not deleting... aborting")
				exit
				kok=1;
	# What if sensor table already exists?
	sensortabfile = dbdir+dbnam+".sensor"
	if os.path.isfile(sensortabfile):
		kok = 0;
		while kok==0:
			yn = raw_input("sensor table already seems to exist. Overwrite/append/cancel? [o/a/c] ")
			if yn=="o":
				print("Ok,deleting and overwriting...")
				os.remove(sensortabfile)
				kok=1;
			elif yn=="a":
				print("Ok, appending...")
				kok=1;
			elif yn=="c":
				raise SystemExit("Ok, not deleting... aborting")
				exit
				kok=1;

	with closing(dbopen(dbdir+dbnam,'r+')) as db:
		dbsch = db.lookup(table = 'sitechan')
		dbsen = db.lookup(table = 'sensor')
		dbins = db.lookup(table = 'instrument')

		for inn in range(0,len(inv.networks)):
			network = inv[inn]
			for iss in range(0, len(network.stations)):
				station = inv[inn][iss]
				for icc in range(0,len(station.channels)):
					channel = inv[inn][iss][icc]
					if (channel.code[1]=='H' or channel.code[1:3]=='DH') and channel.start_date>datemin:
						insnam = channel.sensor.description
						try:
							dbins.record = dbins.find('insname == "%s"' %insnam)
							inid = dbins.getv('inid')
						except DbfindEnd:
							print("No such instrument, sorry =======================================")
							inid = []
						sens_on_date = utc2epoch(station.start_date)
						sens_off_date = utc2epoch(station.end_date)
				
						chan_on_date = channel.start_date.year*1000 + channel.start_date.julday
						chan_off_date = channel.end_date.year*1000 + channel.end_date.julday
						chanid = dbsch.nextid('chanid')

						# flds_sch = ('sta','chan','ondate','offdate','chanid','hang','vang','descrip')
						# flds_sen = ('sta','chan','time','endtime','inid','chanid')
						vals_sch = (station.code,channel.code,chan_on_date,chan_off_date,chanid,\
									float(channel.azimuth),float(channel.dip)+90,\
									find_respfile(respdir,station.code,channel.code))
						vals_sen = (station.code,channel.code,sens_on_date,sens_off_date,inid[0],chanid)
						fldvals_sch = zip(flds_sch, vals_sch)
						fldvals_sen = zip(flds_sen, vals_sen)

						try:
							print("chanid=%.0f, inid=%.0f, %s, %s, %s, %s" %\
									  (chanid,inid[0],station.code,channel.code,chan_on_date,chan_off_date))
							dbaddrec = dbsch.addv(*fldvals_sch)
						except DbaddvError:
							chanid = chanid - 1
							print "Error with adding this row to sitechan..."
							try:
								ifrep = dbsch.find('sta=="%s" && chan=="%s"' % (station.code,channel.code))
								print "Skipping repeated station+chan %s, %s" % (station.code,channel.code)
							except DbfindEnd:
								pdb.set_trace()
								raise SystemExit("Something wrong - won't add sta,\
														but not already in there")
						try:
							dbaddrec = dbsen.addv(*fldvals_sen)
						except DbaddvError:
							chanid = chanid - 1
							print "Error with adding this row to sensor..."
							try:
								ifrep = dbsen.find('sta=="%s" && chan=="%s"' % (station.code,channel.code))
								print "Skipping repeated station+chan %s, %s" % (station.code,channel.code)
							except DbfindEnd:
								pdb.set_trace()
								raise SystemExit("Something wrong - won't add sta,\
														but not already in there")
	return

def make_instrument(inv,dbdir,dbnam,respdir):
    flds = ('inid','insname','instype','samprate','dir','dfile')
    
    datemin = UTCDateTime(2011,1,1)
    
    #  What if instrument table already exists?
    instabfile = dbdir+dbnam+".instrument"
    if os.path.isfile(instabfile):
        kok = 0;
        while kok==0:
            yn = raw_input("Instrument table already seems to exist. Overwrite/append/cancel? [o/a/c] ")
            if yn=="o":
                print("Ok,deleting and overwriting...")
                os.remove(instabfile)
                kok=1;
            elif yn=="a":
                print("Ok, appending...")
                kok=1;
            elif yn=="c":
                raise SystemExit("Ok, not deleting... aborting")
                exit
                kok=1;
    
    with closing(dbopen(dbdir+dbnam,'r+')) as db:
		dbins = db.lookup(table='instrument')
		inid = dbins.record_count
		insnams_all = []
		for inn in range(0,len(inv.networks)):
			network = inv[inn]
			for iss in range(0, len(network.stations)):
				station = inv[inn][iss]
				for icc in range(0,len(station.channels)):
					channel = inv[inn][iss][icc]
					if (channel.code[1]=='H' or channel.code[1:3]=='DH') and channel.start_date>datemin:
						insnam = channel.sensor.description
						if not list(set(insnams_all) & set([insnam])):
							inid=inid+1
							insnams_all.append(insnam)
							respfile = find_respfile(respdir,station.code,channel.code)
							# if respfile is not found write own respfile with write_sac_response_file)
							if not respfile
								respfile = write_sac_response_file(channel,\
								 "junk_respfile",network.code,station.code)

							cprespfile = "SAC_PZs__CAdb_instr"+str(inid).zfill(2)							
							shutil.copy(respdir+respfile,respdir+cprespfile)

							if network.code == '7D':
								instype = 'O_'+channel.code[0:2]+"x"
							else:
								instype =	 'L_'+channel.code[0:2]+"x"

							print("Adding instrument "+str(inid)+", " +insnam)
							vals = (inid,insnam,instype,channel.sample_rate,respdir,cprespfile)
							fldvals = zip(flds, vals)
							dbaddrec = dbins.addv(*fldvals)
	
def count_stas(inv):
	"""Counts the number of stations in an inventory by looping through networks"""
	nstas = 0
	for inn in range(0,len(inv.networks)):
		nstas = nstas + len(inv.networks[inn].stations)
	
	return nstas 
		
def find_respfile(respdir,sta,chan):
    respfile=[]
    for irf in os.listdir(respdir):
        if irf.split('_')[3]==sta and irf.split('_')[4]==chan:
            respfile = irf
            print("Station %s channel %s has response file %s" % (sta,chan,irf) )
            break
    if not respfile:
       print("No instrument response file found for station %s channel %s" % (sta,chan))
    return respfile

def utc2epoch(utctime):
    epochtime = utctime - UTCDateTime(1970,1,1,0,0,0)
    return epochtime
    
def write_sac_response_file(channel_obj,respfile,network=None,station=None,comment=None):
	"""writes a sac-format instrument/channel response file from an IRIS-requested channel
	 object of the sort included in inventories that are requested at the 'response' level
	 i.e. channel_obj = inventory[inn][iss][icc]
	 respfile is the full path to the output file (inc. its name)
	 network and station are strings that give the network code and station code
	 comment is for optional comments (a string)"""
	resp = channel.response
	pzs = resp.get_paz()
	
	sensitivity = 1;
	for istg in range(0,len(resp.response_stages)): 
		sensitivity = sensitivity*resp.response_stages[istg].stage_gain; 
	A0 = resp.response_stages[0].normalization_factor
	C = A0*sensitivity
	
	if not comment: comment='N/A'
	
	with open(respfile,'w') as f:
		f.write('* **********************************\n')
		f.write('* NETWORK   (KNETWK): %s\n' % network.code)
		f.write('* STATION    (KSTNM): %s\n' % station.code)
		f.write('* LOCATION   (KSTNM): \n'                 )
		f.write('* CHANNEL   (KCMPNM): %s\n' % channel.code)
		f.write('* CREATED           : %s\n' % str(UTCDateTime()).split('.')[0])
		f.write('* START             : %s\n' % str(channel.start_date).split('.')[0])
		f.write('* END               : %s\n' % str(channel.end_date).split('.')[0])
		f.write('* DESCRIPTION       : %s\n' % station.site.name)
		f.write('* LATITUDE          : %s\n' % channel.latitude)
		f.write('* LONGITUDE         : %s\n' % channel.longitude)
		f.write('* ELEVATION         : %s\n' % channel.elevation)
		f.write('* DEPTH             : %s\n' % channel.depth)
		f.write('* DIP               : %s\n' % str(float(channel.dip)+90)) # to satisfy convention of dip from vertical +ive
		f.write('* AZIMIUTH          : %s\n' % channel.azimuth)
		f.write('* SAMPLE RATE       : %s\n' % channel.sample_rate)
		f.write('* INPUT UNIT        : %s\n' % resp.response_stages[0].input_units)
		f.write('* OUTPUT UNIT       : %s\n' % resp.response_stages[-1].output_units)
		f.write('* INSTTYPE          : %s\n' % channel.sensor.description)
		f.write('* INSTGAIN          : %.6e (%s)\n' \
					   % (resp.response_stages[0].stage_gain,resp.response_stages[0].input_units))
		f.write('* COMMENT           : %s\n' % comment)
		f.write('* SENSITIVITY       : %.6e (%s)\n' \
					   % (sensitivity,resp.response_stages[0].input_units))
		f.write('* A0                : %.6e\n' % resp.response_stages[0].normalization_factor)
		f.write('* **********************************\n')
		
		pzs.zeros.insert(0,0j) # add extra zero that always seems to be there in the SAC examples
		f.write('ZEROS   %.0f\n' % len(pzs.zeros))
		for izz in range(0,len(pzs.zeros)):
			zz = pzs.zeros[izz]; 
			mpzr = '+'; mpzi = '+'
			if np.real(zz)<0: mpzr = ''
			if np.imag(zz)<0: mpzi = ''
			f.write('        %s%.6e   %s%.6e\n' % (mpzr,np.real(zz),mpzi,np.imag(zz)))

		f.write('POLES   %.0f\n' % len(pzs.poles))
		for ipp in range(0,len(pzs.poles)):
			pp = pzs.poles[ipp]; 
			mppr = '+'; mppi = '+'
			if np.real(pp)<0: mppr = ''
			if np.imag(pp)<0: mppi = ''
			f.write('        %s%.6e   %s%.6e\n' % (mppr,np.real(pp),mppi,np.imag(pp)))


		f.write('CONSTANT        +%.6e\n' % C)

	return respfile

def make_instrument_old(inv):
    dbdir = "/Users/zeilon/Work/CASCADIA/CAdb/" # needs final slash
    dbnam = "cascattendb"
    
    respdir = "/Users/zeilon/Work/CASCADIA/CAdb/response/" # needs final slash
    
    flds = ('inid','insname','instype','samprate','dir','dfile')
    
    ignorechans = ['HKO']
    datemin = UTCDateTime(2011,1,1)
    
    #  What if site table already exists

    instabfile = dbdir+dbnam+".instrument"
    if os.path.isfile(instabfile):
        kok = 0;
        while kok==0:
            yn = raw_input("Instrument table already seems to exist. Overwrite/append/cancel? [o/a/c] ")
            if yn=="o":
                print("Ok,deleting and overwriting...")
                os.remove(instabfile)
                kok=1;
            elif yn=="a":
                print("Ok, appending...")
                kok=1;
            elif yn=="c":
                raise SystemExit("Ok, not deleting... aborting")
                exit
                kok=1;
    
    with closing(dbopen(dbdir+dbnam,'r+')) as db:
        dbins = db.lookup(table='instrument')
        
        inid = dbins.record_count

        for inn in range(0,len(inv.networks)):
            network = inv[inn]
            for iss in range(0, len(network.stations)):
                station = inv[inn][iss]
                for icc in range(0,len(station.channels)):
                    channel = inv[inn][iss][icc]
#                     if list(set(ignorechans) & set([channel.code])) \
#                                         or channel.code[1]=='X' \
#                                         or channel.code[1]=='Y' \
#                                         or channel.code[1]=='L' \
#                                         or channel.code[1]=='N' \
#                                         or channel.code[1]=='B':
#                         print("Ignoring channel %s..." % channel.code)
#                     else:
                    if (channel.code[1]=='H' or channel.code[1:3]=='DH') and channel.start_date>datemin:
                        if network.code == '7D':
                            instype = 'O_'+channel.code
                        else:
                            instype = 'L_'+channel.code

                        insnam = channel.sensor.description
                        inid=inid+1
                        respfile = find_respfile(respdir,station.code,channel.code)

                        print("Adding instrument "+str(inid)+", " +insnam+" for channel "+channel.code)
                        vals = (inid,insnam,instype,channel.sample_rate,respdir,respfile)
                        fldvals = zip(flds, vals)
                        dbaddrec = dbins.addv(*fldvals)
		
def make_sitechan_sensor_old(invfile):
	"""This subroutine reads an obspy inventory and creates sitechan and sensor tables."""
	
	dbdir = "/Users/zeilon/Work/CASCADIA/CAdb/" # needs final slash
	dbnam= "cascattendb"
	respdir = "/Users/zeilon/Work/CASCADIA/CAdb/response/" # needs final slash

	print("Reading inventory file %s") % invfile
	from obspy import read_inventory
	inv = read_inventory(invfile)

	#Check if instrument table exists
	with closing(dbopen(dbdir+dbnam,'r')) as db:
		dbins = db.lookup(table='instrument')
		insrecs = dbins.record_count
	if not (insrecs > 0):
		print("No instrument table yet... writing one")
		make_instrument(inv,dbdir,dbnam,respdir)
	
	print "Writing .sitechan + .sensor tables"
	flds_sitechan = ('sta','chan','ondate','offdate','chanid','hang','vang','description')
	flds_sensor = ('sta','chan','time','endtime','inid','chanid')
			
	ifwrite = raw_input("Write %.0f stations to origin table? [y/n] " % count_stas(inv))
	if ifwrite=="y":
		# 	What if site table already exists
		sitetabfile = dbdir+dbnam+".site"
		if os.path.isfile(sitetabfile):
			kok = 0;
			while kok==0:
				yn = raw_input("Site table already seems to exist. Overwrite/append/cancel? [o/a/c] ")
				if yn=="o":
					print("Ok,deleting and overwriting...")
					os.remove(sitetabfile)
					kok=1;
				elif yn=="a":
					print("Ok, appending...")
					kok=1;
				elif yn=="c":
					raise SystemExit("Ok, not deleting... aborting")
					exit
					kok=1;

	
		with closing(dbopen(dbdir+dbnam,'r+')) as db:
			dbsi = db.lookup(table = 'site')
	
			ista = 0
			for inn in range(0,len(inv.networks)):
				network = inv.networks[inn].code.encode()
				print "\n ====Network %s" % network
				for iss in range(0,len(inv.networks[inn].stations)):
					station = inv.networks[inn].stations[iss]
					sta = station.code.encode()
					on_date = station.start_date.year*1000 + station.start_date.julday
					off_date = station.end_date.year*1000 + station.end_date.julday
					vals_site = (sta,on_date,off_date,\
							station.latitude,station.longitude,station.elevation)
					fldvals_site = zip(flds_site, vals_site)
					print "--------Station %5s %5.2f %7.2f "\
									% (sta,station.latitude,station.longitude)
					try:
						dbaddrec = dbsi.addv(*fldvals_site)		
					except DbaddvError:
						print "Error with adding this row..."
						try:
							ifrep = dbsi.find('sta=="'+sta+'"')
							print "Skipping repeated station %s" % sta
						except DbfindEnd:
							pdb.set_trace()
							raise SystemExit("Something wrong - won't add sta,\
							 						but not already in there")
			
	return

#########################  MAIN   ######################################################   	
def usage():
	"Displays the correct usage of makeAntelope_origin.py"
	print "Usage:  %s arg1 arg2 arg3" % os.path.basename(sys.argv[0])

def main():
	if len(sys.argv) == 9:
		# Give name of aspect output file
		ifile=sys.argv[1]
		# Give name of translated output file
		ofile=sys.argv[2]
		#model="M11"
		model=sys.argv[3]
		#freq=1
		freq=float(sys.argv[4])
		#freq=1
		N=int(sys.argv[5])
		RFmodule=sys.argv[6]
		cfile=sys.argv[7]
		gsopt=float(sys.argv[8])
		set_constants_file(cfile)
		convaspect(ifile,ofile,model,freq,N,RFmodule,gsopt)
	else:	
		usage()
		sys.exit(2)
		
	return

if __name__ == "__main__":
    main()