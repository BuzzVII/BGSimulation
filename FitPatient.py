#!/usr/bin/python2.7
# Finds the value of the weibull shape parameter that minimizes
# the residuals of a patient PSD and FPP simulation
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
import os
import sys
import fnmatch

log = lg.logger(1)
log.info('creating patient object')

fileList = []

for path, dirs, files in os.walk(os.path.abspath("/home/uqkweegi/Documents/Data/dbsdata/OLDMER/Patient38")):
    for filename in fnmatch.filter(files,"*.wav"):
        fileList.append(os.path.join(path, filename))


fid = open('/home/uqkweegi/Documents/Data/leastsq38.xml','w')

fid.write('<?xml version="1.0" encoding="UTF-8" ?>\n')

fid.write('<PATIENTFITS>\n')


for fileName in fileList:
    log.info(fileName)
    patient = FPPwrapper.Patient(patient_file=fileName);

    def fit_func(x0):
        return patient.fpp(x0[0],x0[1])

    log.info('starting minimization')
    x0 = [30, 0.5]
    Cmin = leastsq(fit_func, x0,full_output=True,xtol=1e-8,ftol=1e-8,epsfcn=0.1)
    log.info('Weibull parameters A = ' + str(Cmin[0][0]) + ', B = ' + str(Cmin[0][1]))
    fid.write("\t<RECORDING file='"+fileName+"' a='"+str(Cmin[0][0])+"' b='"+str(Cmin[0][1])+" />\n")
    fid.flush()

fid.write('</PATIENTFITS>')
fid.close()
