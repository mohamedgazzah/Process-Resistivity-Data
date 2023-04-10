# Process-Resistivity-Data

This code takes in a dat file of Resistivity vs Magnetic field for n-different temperatures. Rather than creating n different files for eache temperature measurement, you can compile them all into one measurement in one file and process that data using this function. The other inputs are whether the measurement is MR or Hall(if it's a hall measurement we asymmetrize, if it's MR we symmetrize), Hysteresis temperature(if your measurement is 4 loop there's a different process than a two loop measurement), minimum difference between temperatures(when a machine is set to 300K for example, it's never actually at 300K but instead fluctuates, so this input gives a minimum temperature between one measuremnt at 300K and 305K to know that it's a different temperature rather than a fluctuation of temperature). 

This code will produce 4 plots if given 3 channels for measurment; assume you do a measurment for MR Hall and Hall with a 3 channel machine, then you will have 4 plots, first is the MR resistivity from channel 1, the second is Change in MR, third is Hall resistivity in 2 channel and fourth is Hall resistivity in channel 3(All temperatures are plotted on the same figure). If you have a measurment with less than 3 channels, put the other channel entries as "".

The file which processes the data is attached and works in the following way if interested in variable post processing procedures. 
[Data] = ProcessResData(filename,Ch1,Ch2,Ch3,Hysteresis,Diff)
 HOW TO USE
 filename is name of the file
 Ch1 is a string either Hall or MR for your first channel
 Ch2 is a string either Hall or MR for your second channel
 Ch3 is a string either Hall or MR for your third channel
 Hysteresis is the temperature at which hysteresis stops
