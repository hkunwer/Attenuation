# %%
#attenuation
# from attenuation import processANSS, readANSS

#rtergpy
from rtergpy.run import defaults, event, etime2name, src2ergs

#obspy
from obspy import UTCDateTime
from obspy import read
from obspy import read_inventory

import pandas as pd 

import csv

def processANSS():

        # Read the data from the file
    with open('ANSS_data.txt', 'r') as file:
        lines = file.readlines()

    # Remove the header line
    header = lines[0].strip().split()[:-3]
    lines = lines[1:]

    # Modify the data
    modified_lines = []
    for line in lines:
        columns = line.split()
        modified_line = columns[0:6]
        modified_lines.append(modified_line)
    
    modified_lines.insert(0, header)

    # Write the modified data to a CSV file
    with open('ANSS_processed_data.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(modified_lines)
# %%
        
def readANSS(): # should be used after ANSS is processed
    
    Defaults=defaults()
    Event=event()
    Defaults.src='RASPISHAKE'
    Event.newData=False # True would collect new data. False if already downloaded.
    Event.ecount='00'
    edateold=""
    ANSS=pd.read_csv('ANSS_processed_data.csv', sep='\s+', comment="#")
    for index, EQ in ANSS.iterrows():
        eloc = [EQ.LAT,EQ.LONG,EQ.DEPTH] 
        year,mo,dd = EQ.DATE.split('/')
        hh,mn,sec = EQ.TIME.split(':')
        etime=(UTCDateTime(int(year),int(mo),int(dd),int(hh),int(mn),float(sec)))
        # iterate ecount
        if EQ.DATE == edateold:
            Event.ecount=str(int(Event.ecount)+1).zfill(2)
        else:
            Event.ecount='00'
        edateold=EQ.DATE
        Event.eventname=etime2name(etime,ecount=Event.ecount)
        Event.origin=[eloc,etime]
        Event.focmech=[EQ.STK, EQ.DP, EQ.RKE] # phi,delta,lmbda

        print("\n\n"+Event.eventname+" ===============================")
        try:
            src2ergs(Defaults=Defaults,Event=Event)
        except:
            print("ERROR: running on "+Event.eventname+" failed!!!!\n\n")


# %%
processANSS()
readANSS()

# %%



