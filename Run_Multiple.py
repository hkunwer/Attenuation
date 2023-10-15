#obspy
from obspy import UTCDateTime
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")
#rtergpy (Andy's code)
from rtergpy.run import defaults, event, etime2name
from rtergpy.waveforms import getwaves, get_respinv
#attenuation (Hiba's code)
from Attenuation.Drafts.AttenuationFunctions import processANSS, filtering, maxamp_calc, update_and_save_dataframe, maxamp_plot

# Processing and Reading information about event stored in ANSS_data.txt
Defaults = defaults()
Event = event()
Defaults.src='RASPISHAKE'
Defaults.network='AM'
Defaults.chan= 'EHZ'
Defaults.stationrange=[0.5,10.0]
Event.ecount='00'
Event.newData = True   # False means use already downloaded data
edateold=""

ANSS = processANSS() 

for index, EQ in ANSS.iterrows(): #organizing data to use details
    rads = Defaults.stationrange
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
        
    if Defaults.src == 'RASPISHAKE':
        Event.iter = 'RS'

    else: 
        Event.iter = 'IRIS'
        
    edateold = EQ.Date
    Event.eventname=etime2name(etime,ecount=Event.ecount)+Event.iter
    Event.origin=[eloc,etime]
    
    print("\n\n"+Event.eventname+" ===============================")
    
    try:
        st, df = [], []
        st, df = getwaves(Defaults,Event)
        #inventory = get_respinv(Defaults.network,eloc,etime,rads,Defaults.chan,Defaults.src)
        stp = filtering(st) # filter stream for instrument response and taper
        df_new = maxamp_calc(stp, EQ, Defaults, etime, eloc) #calculate max amps and dist for each tr at freq bands
        update_and_save_dataframe(df_new)
    except:
        print("ERROR: "+Event.eventname+" could not complete getwaves") 
        
df_ALL = maxamp_plot()