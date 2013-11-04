# Finds the value of the weibull shape parameter that minimizes
# the L2 difference of a patient PSD and FPP simulation
#
# Kristian Weegink, Nov 2013

import numpy as np
from scipy import signal
from matplotlib import pylab
import FPPwrapper
import BrainSim as BS
import math
import scipy.io.wavfile as wav
from scipy.optimize import leastsq
import Logger as lg

log = lg.logger(1)
log.info('creating patient object')

patient = FPPwrapper.Patient(patient_file='/home/uqkweegi/Documents/Data/dbsdata/OLD MER/Patient32/32_S_L_E1_20.wav');

def fit_func(x0):
    return patient.fpp(x0[0],x0[1])

log.info('starting minimization')
x0 = [30, 0.5]
Cmin = leastsq(fit_func, x0,full_output=True)
