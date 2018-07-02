#https://www.iaea.org/inis/collection/NCLCollectionStore/_Public/25/044/25044984.pdf

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.special as ss
from scipy.optimize import curve_fit
from pathlib import Path

class drs4_data_processor(object):
    def __init__(self,folder_name):

        self.data_folder = folder_name
        self.glob_text = '**/*charge*.txt'

        self.voltage_settings=[]
        self.params=dict()
        self.param_errors=dict()

        pathlist = Path(self.data_folder).glob(self.glob_text)

        self.data_files = []
        self.charge_histograms = dict()

        for path in pathlist:
            # because path is object not string
            path_in_str = str(path)
            self.data_files.append(path_in_str)
            print('Processing ', path_in_str)

            # Find voltage from filename
            str_num=''
            for i, c in enumerate(path_in_str):
                if c.isdigit():
                    str_num+=c

            if(len(str_num)==3):
                self.voltage_settings.append((1250.0/5.0)*float(str_num)/float(100.0))
            if(len(str_num)==2):
                self.voltage_settings.append((1250.0/5.0)*float(str_num)/float(10.0))

            # Process the file
            print("Voltage (kV): ", (self.voltage_settings[-1]))
            self.charge_histograms[(self.voltage_settings[-1])] = drs4_charge_histogram(path_in_str);


    def recalibrate(self):
        for i,ch in enumerate(sorted(self.charge_histograms)):
            print("Voltage (kV): ", (self.voltage_settings[i]))
            self.charge_histograms[ch].calibration.multires_calibrate()
            self.charge_histograms[ch].calibration.plot(str("Voltage = {}".format(self.voltage_settings[i])))
            self.charge_histograms[ch].calibration.print_fit()

    def set_consistent_bounds(self):

        kys = list(sorted(self.charge_histograms.keys()))
        for i,hv in enumerate(sorted(self.charge_histograms)):
            print("Filename: ", self.data_files[i])
            if(i!=0):
                self.charge_histograms[kys[i]].calibration.set_lower_bounds(self.charge_histograms[kys[i-1]].calibration)
            if(i!=len(self.charge_histograms)-1):
                self.charge_histograms[kys[i]].calibration.set_upper_bounds(self.charge_histograms[kys[i+1]].calibration)

    def plot_charge_resolution(self):
        voltage_scale=[]
        charge_res=[]
        error_charge_res=[]
        for i,hv in enumerate(sorted(self.charge_histograms)):
            voltage_scale.append(hv)
            charge_res.append(self.charge_histograms[hv].calibration.charge_resolution())
            error_charge_res.append(self.charge_histograms[hv].calibration.charge_resolution_error())

        plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k');
        plt.errorbar(voltage_scale,charge_res,yerr=error_charge_res,fmt='ro')
        plt.xlabel('Voltage (V)',fontsize=16)
        plt.ylabel('Charge Resolution',fontsize=16)
        plt.grid(b=True, which='minor')
        plt.grid(b=True, which='major')
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.tick_params(axis='both', which='minor', labelsize=10)

    def plot_gain(self):
        voltage_scale=[]
        gain = []
        error_gain=[]
        for i,hv in enumerate(sorted(self.charge_histograms)):
            voltage_scale.append(hv)
            gain.append(self.charge_histograms[hv].calibration.gain)
            error_gain.append(self.charge_histograms[hv].calibration.gain_error)

        plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k');
        plt.errorbar(voltage_scale,gain,yerr=error_gain,fmt='ro')
        plt.xlabel('Voltage (V)',fontsize=16)
        plt.ylabel('Gain',fontsize=16)
        ax=plt.gca()
        ax.set_yscale("log", nonposy='clip')
        plt.grid(b=True, which='minor')
        plt.grid(b=True, which='major')
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.tick_params(axis='y', which='major', labelsize=14)
        plt.tick_params(axis='both', which='minor', labelsize=12)

        popt,pcov = curve_fit(self.voltage_gain,voltage_scale,gain,p0=[1e-3,8],maxfev = 1800000)
        fit_gains = self.voltage_gain(voltage_scale,popt[0],popt[1])
        plt.plot(voltage_scale, fit_gains ,'b-.',label='fit')

        self.params['k'] = popt[1]
        errors = np.sqrt(np.diag(pcov))
        self.param_errors['k'] =  errors[1]
        print()
        print("==============")
        print("GAIN FIT PARAMETERS")
        print("==============")
        for var in self.params:
            print(var,' \t = \t{:8.6f} +/- {:8.6f}'.format( self.params[var], self.param_errors[var]))
        print("==============")
        print()


    def voltage_gain(self,volts,A,k):
        return A*(volts**k)


class drs4_charge_histogram():
    def __init__(self,filename):

        data = np.loadtxt(filename,skiprows=2);
        self.charge = -0.5*(data[:,0]+data[:,1])
        self.events = data[:,2]/sum(data[:,2])
        print(self.events.size)
        self.number_of_events = len(self.events)

        self.calibration = pmt_calibration(self.charge, self.events)
        self.calibration.calibrate()

class pmt_calibration():

    def __init__(self,charge,events):

        A= 0.0558348688201
        mu= 2.59419846474
        w= 0.102925765097
        alpha= 0.164162004611
        Q0= 1.7304425357
        sig0= 0.0770888074934
        Q1= 0.272786581019
        sig1= 0.226272555706

        # Initial default parameters for the fit
        self.params = dict();
        self.params['Q0'] = Q0
        self.params['sig0'] = sig0
        self.params['Q1'] = Q1
        self.params['sig1'] = sig1
        self.params['A'] = A
        self.params['alpha'] = alpha
        self.params['w'] = w
        self.params['mu'] = mu
        self.n_max = 2
        self.param_errors = dict()
        self.param_errors['Q0'] = 0.0
        self.param_errors['Q1'] = 0.0
        self.param_errors['sig0'] = 0.0
        self.param_errors['sig1'] = 0.0
        self.param_errors['w'] = 0.0
        self.param_errors['A'] = 0.0
        self.param_errors['mu'] = 0.0
        self.param_errors['alpha'] = 0.0
        self.gain = 0
        self.max_function_evals = int(1000000000)
        self.charge = charge
        self.events = events

        w_bounds  = (0.0,1.0)
        mu_bounds  = (0.1,10.0)
        alpha_bounds  = (0.0,10.0)
        mu_bounds  = (0.01,3.0)
        A_bounds  = (0.0001,10.0)

        # self.lower_param_bounds = (Q0*0.99, Q1*0.6, sig0*0.5,sig1*0.001, alpha_bounds[0],w_bounds[0],mu_bounds[0],A_bounds[0])
        # self.upper_param_bounds = (Q0*1.01, Q1*1.4, sig0*1.5,sig1*2.6, alpha_bounds[1],w_bounds[1],mu_bounds[1],A_bounds[1])
        self.lower_param_bounds = (0, 0,0,0,0 ,w_bounds[0],0,A_bounds[0])
        self.upper_param_bounds = (1e9,1e9,1e9,1e9,1e9,w_bounds[1],1e9,A_bounds[1])

    def set_lower_bounds(self,voltage_below):
        self.lower_param_bounds = (0, voltage_below.params['Q1'], 0,  0, 0, 0, voltage_below.params['mu'],0.001)

        if(self.params['Q1']<voltage_below.params['Q1']):
            self.params['Q1'] = voltage_below.params['Q1'] + 0.01

        if(self.params['mu']<voltage_below.params['mu']):
            self.params['mu'] = voltage_below.params['mu'] + 0.01

        # print('Lower bounds: ', self.lower_param_bounds)
        # print('Upper bounds: ', self.upper_param_bounds)

    def set_upper_bounds(self,voltage_above):
        self.upper_param_bounds = (1e9, voltage_above.params['Q1'], 1e9,  1e9, 1e9, 1, voltage_above.params['mu'], 10.0)

        if(self.params['Q1']>voltage_above.params['Q1']):
            self.params['Q1'] = voltage_above.params['Q1'] - 0.01

        if(self.params['mu']>voltage_above.params['mu']):
            self.params['mu'] = voltage_above.params['mu'] - 0.01

        # print('Lower bounds: ', self.lower_param_bounds)
        # print('Upper bounds: ', self.upper_param_bounds)

    def multires_calibrate(self):
        old_events = self.events
        old_charge = self.charge
        for i in range(4,0,-1):
            self.events = self.events[::int(2**i)]
            self.charge = self.charge[::int(2**i)]
            self.calibrate()
            print('Resolution Level : ', i, ' Chi2 : ', self.chi2())
            self.events = old_events
            self.charge = old_charge

    def calibrate(self):
        # Do the fit
        try:
            popt,pcov = curve_fit(      self.calibration_fit_function,
                                        self.charge,
                                        self.events,
                                        p0=[    self.params['Q0'],
                                                self.params['Q1'],
                                                self.params['sig0'],
                                            self.params['sig1'],
                                            self.params['alpha'],
                                            self.params['w'],
                                            self.params['mu'],
                                            self.params['A']
                                            ],
                                    bounds = (self.lower_param_bounds,self.upper_param_bounds),
                                    maxfev = self.max_function_evals )
        except:
            popt,pcov = curve_fit(      self.calibration_fit_function,
                                self.charge,
                                self.events,
                                p0=[    self.params['Q0'],
                                        self.params['Q1'],
                                        self.params['sig0'],
                                    self.params['sig1'],
                                    self.params['alpha'],
                                    self.params['w'],
                                    self.params['mu'],
                                    self.params['A']
                                    ],
                            maxfev = self.max_function_evals )





        # Save the parameters
        Q0, Q1, sig0, sig1, alpha, w, mu, A = popt;
        self.params['Q0'] = Q0
        self.params['sig0'] = sig0
        self.params['Q1'] = Q1
        self.params['sig1'] = sig1
        self.params['A'] = A
        self.params['alpha'] = alpha
        self.params['w'] = w
        self.params['mu'] = mu

        # Save the fit errors
        errors = np.sqrt(np.diag(pcov))
        self.param_errors['Q0'] =  errors[0]
        self.param_errors['Q1'] =  errors[1]
        self.param_errors['sig0'] =  errors[2]
        self.param_errors['sig1'] =  errors[3]
        self.param_errors['alpha'] =  errors[4]
        self.param_errors['w'] =  errors[5]
        self.param_errors['mu'] =  errors[6]
        self.param_errors['A'] =  errors[7]

        self.gain = self.params['Q1']*1e-12/1.6021765e-19
        self.gain_error = self.param_errors['Q1']*1e-12/1.6021765e-19

    def chi2(self):
        yfit = self.calibration_fit_function(self.charge,self.params['Q0'],self.params['Q1'],self.params['sig0'],self.params['sig1'],self.params['alpha'],self.params['w'],self.params['mu'],self.params['A'])
        return sum((self.events - yfit)**2)/float(self.events.size)

    def calibration_fit_function(self,x,Q0,Q1,sig0,sig1,alpha,w,mu,A):

        # Variable to return from this function
        calib = 0
        # Main sum over a number of photoelectrons
        for n in range(self.n_max):
            calib += self.calibration_photoelectron_spectrum(x,n,Q0,Q1,sig0,sig1,alpha,w,mu,A)

        return A*calib;

    def calibration_photoelectron_spectrum(self,x,n,Q0,Q1,sig0,sig1,alpha,w,mu,A):
        # Sigma for the funtion with n photoelectons
        sig_n = np.sqrt(sig0*sig0+n*sig1*sig1)
        Qn = Q0 + n*Q1

        # The Poisson function
        factorial = np.prod(range(1,n+1))
        poiss = (mu**n)*np.exp(-mu)/factorial

        # The first term
        term1 = self.gaus( x, Qn, sig_n )
        # The second term
        # - arguments of the functions to help write in more compactly
        Qx = (x - Qn - alpha*sig_n*sig_n)
        Q0n = np.abs(Q0 - Qn - sig_n*sig_n*alpha) / (sig_n*np.sqrt(2))
        term2 = 0.5*alpha*np.exp( - alpha*Qx )*(ss.erf( Q0n ) * np.ones_like(x) + np.sign(Qx) * ss.erf( Qx ) )

        # Add this to the running sum of photoelecton terms
        return poiss*((1-w)*term1 + w*term2)

    # Helper functions
    def gaus(self,x,x0,sigma):
        return (1.0/np.sqrt(2*np.pi*sigma*sigma))*np.exp(-((x-x0)**2)/(2*sigma**2))

    def print_fit(self):
        print()
        print("==============")
        print("FIT PARAMETERS")
        print("==============")
        for var in self.params:
            print(var,' \t = \t{:8.6f} +/- {:8.6f}'.format( self.params[var], self.param_errors[var]))
        print("==============")
        print('Gain \t = \t{:8.0f} +/- {:8.0f}'.format( self.gain, self.gain_error))
        print()


    def plot(self,in_label=''):
        plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k');
        plt.plot(self.charge,self.events,'ro')
        plt.plot(self.charge,self.calibration_fit_function(
            self.charge,
            self.params['Q0'],
            self.params['Q1'],
            self.params['sig0'],
            self.params['sig1'],
            self.params['alpha'],
            self.params['w'],
            self.params['mu'],
            self.params['A']),'b-',label=in_label,linewidth=2)
        plt.legend()
        # Plot the components
        # for n in range(2):
        #     plt.plot(self.charge,self.calibration_photoelectron_spectrum(
        #         n,
        #         self.charge,
        #         self.params['Q0'],
        #         self.params['Q1'],
        #         self.params['sig0'],
        #         self.params['sig1'],
        #         self.params['alpha'],
        #         self.params['w'],
        #         self.params['mu'],
        #         self.params['A']),'b-.',label='fit components')

        self.plot_photoelectron_component(0)
        self.plot_photoelectron_component(1)
        plt.xlabel('Charge (pC)',fontsize=16)
        plt.ylim([0.0,0.024])
        plt.ylabel('Events',fontsize=16)
        plt.grid(b=True, which='minor')
        plt.grid(b=True, which='major')
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.tick_params(axis='both', which='minor', labelsize=10)

        # v=0
        # for var in self.params:
        #     label_line = ' \t = \t{:8.6f} +/- {:8.6f}'.format( self.params[var], self.param_errors[var])
            # plt.text(.025,0.018-v,label_line)
            # v+=0.03

    def charge_resolution(self):
        return self.params['sig1']/self.params['Q1']
    def charge_resolution_error(self):
        return self.charge_resolution() * np.sqrt(((self.param_errors['sig1']/self.params['sig1'])**2)+((self.param_errors['Q1']/self.params['Q1'])**2))

    def plot_photoelectron_component(self,n):
        plt.plot(self.charge,self.params['A']*self.calibration_photoelectron_spectrum(self.charge,n,
            self.params['Q0'],
            self.params['Q1'],
            self.params['sig0'],
            self.params['sig1'],
            self.params['alpha'],
            self.params['w'],
            self.params['mu'],
            self.params['A']),'g-.',linewidth=2)
