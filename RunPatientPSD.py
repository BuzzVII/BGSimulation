#
# This script finds the mean MER PSD for each patient then writes an XML file containing all the data
#
# Kristian Weegink - 07/2013

import numpy as np
from scipy import signal
from matplotlib import pylab
import FPP
import BrainSim as BS
import math
import scipy.io.wavfile as wav

# Files containing directory information of MERs
Fid=open('C:\Users\uqkweegi\Documents\Data\DAfit\E1R.txt')
E1R=Fid.read()
Fid.close()
Fid=open('C:\Users\uqkweegi\Documents\Data\DAfit\E1L.txt')
E1L=Fid.read()
Fid.close()
Fid=open('C:\Users\uqkweegi\Documents\Data\DAfit\patients.txt')
patients=Fid.read()
Fid.close()
PRList=E1R.split('\n')
PLList=E1L.split('\n')
PatientList=patients.split('\\')

#Determine how many patients
PatientNumbers=[]
for item in PatientList:
	if len(item) > 7:
		if item[0:7]=='Patient':
			PatientNumbers.append(item[7:-3])
sr,Vt=wav.read(PLList[0])
nfft=2**int(math.log(len(Vt),2))+1
Pxx,freqs=pylab.psd(x=Vt,Fs=sr,NFFT=nfft/10,window=pylab.window_none, noverlap=100)   

#find mean PSD and write to XML file
#
#Add in NFFT,Fs,Window,nooverlap to XML file
fid=open('C:\Users\uqkweegi\Documents\Data\DAfit\patientspsd.xml','w')
fid.write("<PATIENTDATA>\n")
for PNitem in PatientNumbers:
    Pxx-=Pxx; 
    fid.write("\t<PATIENT>\n")
    fid.write("\t\t<ID>"+str(PNitem)+"</ID>\n")
    fid.write("\t\t<SIDE>\n\t\t\t<ID>Left</ID>\n")
    count=1
    for item in PLList:
        if item[56:56+len(PNitem)] == PNitem: 
            sr,Vt=wav.read(item)
            Pxi,freqs=pylab.psd(x=Vt,Fs=sr,NFFT=nfft/10,window=pylab.window_none, noverlap=100)
            Pxx+=Pxi
            count+=1
    fid.write("\t\t\t<NORECORDINGS>"+str(count)+"</NORECORDINGS>\n")
    fid.write("\t\t\t<MEANPSD>\n\t\t\t\t")
    for i in Pxx:
        fid.write(str(i/count)+",")
    fid.write("\n\t\t\t</MEANPSD>\n")
    fid.write("\t\t\t<FREQ>\n\t\t\t\t")
    for i in freqs:
        fid.write(str(i)+",")
    fid.write("\n\t\t\t</FREQ>\n")
    fid.write("\t\t</SIDE>\n")
    fid.write("\t\t<SIDE>\n\t\t\t<ID>Right</ID>\n")
    count=1
    for item in PRList:
        if item[56:56+len(PNitem)] == PNitem: 
            sr,Vt=wav.read(item)
            Pxi,freqs=pylab.psd(x=Vt,Fs=sr,NFFT=nfft/10,window=pylab.window_none, noverlap=100)
            Pxx+=Pxi
            count+=1
    fid.write("\t\t\t<NORECORDINGS>"+str(count)+"</NORECORDINGS>\n")
    fid.write("\t\t\t<MEANPSD>\n\t\t\t\t")
    for i in Pxx:
        fid.write(str(i/count)+",")
    fid.write("\n\t\t\t</MEANPSD>\n")
    fid.write("\t\t\t<FREQ>\n\t\t\t\t")
    for i in freqs:
        fid.write(str(i)+",")
    fid.write("\n\t\t\t</FREQ>\n")
    fid.write("\t\t</SIDE>\n")
    fid.write("\t</PATIENT>\n")
fid.write("</PATIENTDATA>")
fid.close()
