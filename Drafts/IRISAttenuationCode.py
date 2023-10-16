#obspy
from obspy import UTCDateTime
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")
#rtergpy (Andy's code)
from rtergpy.run import defaults, event, etime2name
from rtergpy.waveforms import getwaves, get_respinv
#attenuation (Hiba's code)
from Attenuation.AttenuationFunctions import processANSS, filtering, maxamp_calc, organize_data, maxamp_plot

# Processing and Reading information about event stored in ANSS_data.txt
Defaults = defaults()
Event = event()

Defaults.src='IRIS'
Defaults.network='??'
Defaults.chan='BHZ'
Defaults.stationrange=[1.,10.]
Event.ecount='00'
Event.iter='IRIS'
Event.newData = True   # False means use already downloaded data
edateold=""
ANSS = processANSS() 
for index, EQ in ANSS.iterrows(): #organizing data to use details
    eloc = [EQ.Latitude,EQ.Longitude,EQ.Depth] 
    MagType = [EQ.Mtype]
    MagValue = [EQ.Mag]
    Magnitude = [MagType, MagValue]
    year,mo,dy = EQ.Date.split('-')
    hh,mn,sec = EQ.Time.split(':')
    etime=(UTCDateTime(int(year),int(mo),int(dy),int(hh),int(mn),float(sec)))
    if EQ.Date == edateold:
        Event.ecount=str(int(Event.ecount)+1).zfill(2)
    else:
        Event.ecount='00'
    edateold=EQ.Date
    Event.eventname=etime2name(etime,ecount=Event.ecount)
    Event.origin=[eloc,etime]
    print("\n\n"+Event.eventname+" ===============================")
    try:
        st, df = [], []
        st, df = getwaves(Defaults=Defaults,Event=Event)
    except:
        print("ERROR: running on "+Event.eventname+" failed!!!!\n\n") 
    
eventID = (str(etime)+Event.iter+Event.ecount) # creates unique event name
#inventory = get_respinv(network=Defaults.network,eloc,etime,rads,chan=Defaults.chan,src=Defaults.src) # make an inventory incase needed.
stp = filtering(st) # filter stream for instrument response and taper
df_freq = maxamp_calc(stp, eventID, EQ, Defaults) #calculate max amps and dist for each tr at freq bands
organize_data(df_freq, EQ, etime, eloc, eventID) #creates dataframe of all results
maxamp_plot(eventID,df_freq) #plots maxamps vs distance