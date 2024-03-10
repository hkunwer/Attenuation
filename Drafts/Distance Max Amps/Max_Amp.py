#obspy
from obspy import UTCDateTime
from obspy.taup import TauPyModel
import pandas as pd 
model = TauPyModel(model="iasp91")
#rtergpy (Andy's code)
from rtergpy_reference.run import defaults, event, etime2name
from rtergpy_reference.waveforms import getwaves
#attenuation (Hiba's code)
from AttenuationFunctionsTesting import check_existing_dataframe, processANSStxt, filteringst, maxamp_calc_freq_bands, max_and_normalized_max_plots

# Processing and Reading information about event stored in ANSS_data.txt
Defaults = defaults()
Event = event()
Defaults.src='RASPISHAKE'
Defaults.network='AM'
Defaults.chan= 'EHZ'
Defaults.stationrange=[0.5,10.0]
Event.ecount='00'
Event.newData = False  # False means use already downloaded data
edateold=""

#ANSS = processANSStxt() #If you have raw ANSS list, process with this
ANSS = pd.read_csv('ANSS_processed_data.csv')
check_existing_dataframe() # Check if pkl file exists to add new results too, create new one if not accessible.

for index, EQ in ANSS.iterrows(): #organizing data to use details
    rads = Defaults.stationrange
    eloc = [EQ.Latitude,EQ.Longitude,EQ.Depth] 
    MagType = [EQ.Mtype]
    MagValue = [EQ.Mag]
    Magnitude = [MagType, MagValue]
    year,mo,dy = EQ.Date.split('-')
    hh,mn,sec = EQ.Time.split(':')
    etime=(UTCDateTime(int(year),int(mo),int(dy),int(hh),int(mn),float(sec)))
    #print(etime)
    
    if EQ.Date == edateold:
        Event.ecount=str(int(Event.ecount)+1).zfill(2)
    else:
        Event.ecount='00'
        
    if Defaults.src == 'RASPISHAKE':
        Event.iter = 'RS'

    else: 
        Event.iter = 'IRIS'
        
    edateold = EQ.Date
    Event.eventname=etime2name(etime,ecount=Event.ecount)+Event.iter
    Event.origin=[eloc,etime]
    #print(Event.origin)
    
    print("\n\n"+Event.eventname+" ===============================")
    
    st, df = [], []
    try:
        st, df = getwaves(Defaults, Event)
        print("Completed getwaves")
        stp = filteringst(st) 
        maxamp_calc_freq_bands(stp, EQ, Defaults, etime, eloc)
    except:
        print("ERROR: "+Event.eventname+" could not complete getwaves") 