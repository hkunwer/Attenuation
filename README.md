# READ ME

#edit an ANSS_data.txt to input data
#searchANSS -C-73.3710/18.4938/15 -T2019/01/01/2022/12/31 -M4
#Use searchANSS to collect data to put in the ANSS_data.txt file
#remember to run in rtergpy environment

# Processing ANSS
#done by function processANSS, 
#removes last two columns, prints with header, separates everything by commas
#creates a modified CSV file


# ATTENUATION

#we are calculating the peak amplitude for individual frequency bands using Velocity squared data
#We then obtain a near-field attenuation correction, t*, from the slope of a 3D plot of frequency and Velocity squared over distance
#This will then used, with corrections for S-wave dominance and majority horizontal crustal pathway, in energy calculations in the near-field for tsunamigenic early warning system, rtergpy. 


