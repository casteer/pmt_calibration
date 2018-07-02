# PMT Calibration Library

This is a small python library for fitting of PMT charge histograms for calibration. 

In order to calibrate photomultipliers (PMTs) often the gain of the PMT needs to be measured before pressing ahead and building the latest and greatest radiation detection instrumentation. This is typically achieved using an LED, with an appropriate wavelength, which is flashed in front of the PMT; the LED is connected to a pulse generator with a short pulse width and its sync output connected to the data acquisition equipment. 

The fit is to the model described [here](https://www.iaea.org/inis/collection/NCLCollectionStore/_Public/25/044/25044984.pdf), which (I think) was published in Nuclear Instruments and Methods A. 

### How the example data was generated

In our case, for which this libary was written, our data acquisition is performed by a DRS4 evaluation board produced by PSI. The histogram data is obtained using the DRS4 Oscilloscope program as follows: First, the measure button is clicked and then the 'Gated Charge' button is switch on; second, the integration region is selected using the A and B cross-hairs; third, histograms are displayed (click the display button at the bottom), whilst selecting the number of events in the drop-down box; lastly, when the histogram has been fully populated, the (alternative histogram) save button is clicked (n.b. this is in the top-right of the black histogram window, and **not** the other save button on the far right-hand side). 

The example data has the format of a couple of headerlines, then column 1 is the left-hand edge of a bin, column 2 is the right-hand edge of the same bin, and column 3 is the number of events in the bin. 

### How do I use this for my data? 

Well, I guess the value of this libray is that it's freely available for you to study and understand. If you find any corrections or comments, please do get in touch and let me know. 

If you would like to modify this for your own use, please fork the repo and you'll need to change lines 128-130 in the python library. So change the line that loads the data, i.e. this one: 

data = np.loadtxt(filename,skiprows=2);

Then if you need to set the class variables self.charge and self.events, which in the library are given by

self.charge = -0.5*(data[:,0]+data[:,1])
self.events = data[:,2]/sum(data[:,2])

I hope you find it useful and wish you the best of luck with it. 
