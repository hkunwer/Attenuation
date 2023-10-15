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

#we are calculating the peak amplitude for individual frequency bands
#we will then use this max amplitude along with a geometric spreading estimate to get
#attenuation using the methods in Ismail et al 2023 and 2020 and McNamara et al 2014 as #guidance.
