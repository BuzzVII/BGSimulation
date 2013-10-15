import numpy as np
from scipy import signal
from matplotlib import pylab
import FPP
import BrainSim as BS
import math
import scipy.io.wavfile as wav


#----------------#
# Loop over different DA levels.
#----------------#

DAstep=0.01

DAlevels=[]

for i in range(100):
    DAlevels.append(DAstep*i)

#----------------#
# load patient data.
#----------------#
Vt,t=FPP.FPP(N=10000, dt=1./24000, STNrate='C:\\Users\\uqkweegi\\Documents\\Data\\DAfit\\STN',plotAll=False)
nfft=2**int(math.log(len(Vt),2))+1
sr=1/(t[2]-t[1])
Pxx,freqs=pylab.psd(x=Vt,Fs=sr,NFFT=nfft/10,window=pylab.window_none, noverlap=100)

fid=open('C:\Users\uqkweegi\Documents\Data\DAfit\simulationpsd24k.xml','w')
fid.write("<SIMULATIONDATA>\n")
for Dopamine in DAlevels:
#----------------#
# run simulation need to make AP sampled at dt.
#----------------#
    print "DA level ="+str(Dopamine)
    fid.write("\t<SIMULATION>\n")
    fid.write("\t\t<DA>"+str(Dopamine)+"</DA>\n")
    fid.write("\t\t<DT>1/24000</DT>\n")
    fid.write("\t\t<N>10000</N>\n")
    fid.write("\t\t<NUMBER>10</NUMBER>\n")
    Pxx-=Pxx; 
    BS.main(dt=1./24000,DA=Dopamine)
    for i in range(10):
        Vt,t=FPP.FPP(N=10000, dt=1./24000, STNrate='C:\\Users\\uqkweegi\\Documents\\Data\\DAfit\\STN',plotAll=False)
        sr=1/(t[2]-t[1])
        nfft=2**int(math.log(len(Vt),2))+1
        Pxi,freqs=pylab.psd(x=17*Vt,Fs=sr,NFFT=nfft/10,window=pylab.window_none, noverlap=100)
        Pxx+=Pxi
    fid.write("\t\t<MEANPSD>\n\t\t\t")
    for i in Pxx:
        fid.write(str(i/10)+",")
    fid.write("\n\t\t</MEANPSD>\n")
    fid.write("\t\t<FREQ>\n\t\t\t")
    for i in freqs:
        fid.write(str(i)+",")
    fid.write("\n\t\t</FREQ>\n")
    fid.write("\t</SIMULATION>\n")
fid.write("</SIMULATIONDATA>\n")    
fid.close()
#----------------#
# Measure difference.
#----------------#

#----------------#
# Write results to file.
#----------------#
