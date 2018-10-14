# Neural mass simulation of basal ganglia, thalamus and cortex to investigate how MER noise is produced and if it can be used to replace the STN. This is then extended to include DBS.
#
# Kristian Weegink: uqkweegi@uq.edu.au
# 12/12/2012 

import math
import numpy as np
import numpy.random


#----------------------------------------------#
class CNeuralStruct:
#----------------------------------------------#
# This class defines an object that represents a neural nucleus. The neural nucleus is represented by a neural mass that can then be stepped forward in time using an Euler method. 	
#----------------------------------------------#
	#construct object
	def __init__(self,name,parameters,path='C:\\Users\\Kristian\\Dropbox\\phd\\Data\\'):
		
		self.name=name							#set the name of the structure type
		
		self.paramInit(parameters)				#set all the parameters for coupling each structure
		
		self.Reset()							#set the structure to it's initial state
		try:
			self.file_name=open(path+name,'w') #'H:Data/BrainSim'+name,'w')	#prepare the log file for saving the simulation
			print('log file '+path+name+' opened')
		except:
			print('could not open '+path+name)
			exit(1)
	
#----------------------------------------------#	
	#destroy object
	def __del__(self):
		try:
			self.file_name.close()						#if the log file has been open for writing, close it
			print ('log file for '+self.name+' closed')
		except:											 
			print ('no log file open for '+self.name)	#inform user that the log file wasn't open for this object

#----------------------------------------------#
	#set all the parameters for coupling each structure object
	def paramInit(self, parameters):
		if (len(parameters)!=14):							#check that there are the correct amount of parameters
			print (self.name+' : Invalid paramater list')
			exit()
		self.weightGPi	=	parameters[0]					#assign weightings between structure XXX and this object
		self.weightGPe	=	parameters[1]
		self.weightD1	=	parameters[2]
		self.weightD2	=	parameters[3]
		self.weightSTN	=	parameters[4]
		self.weightFS	=	parameters[5]
		self.weightCtx	=	parameters[6]
		self.weightTha	=	parameters[7]
		self.weightTRN	=	parameters[8]
		self.tau		= 	parameters[9]					#set synaptic time constant
		self.maxrate	= 	parameters[10]					#configure the sigmoid function
		self.slope		=	parameters[11]
		self.maxpot		=	parameters[12]
		self.R			=	parameters[13]					#set synaptic resistance
		self.ex_in		=	1								#configure if structure is excitatory of inhibitory
		if (self.name=='GPi' or self.name=='D2' or self.name=='D1' or self.name=='GPe'):
			self.ex_in	=	-1

#----------------------------------------------#
	#set the structure to the initial resting state
	def Reset(self):
		numpy.random.seed(ord(self.name[0])*ord(self.name[-1]))
		self.history		=	numpy.random.random_sample(10)	#create random history
		
		self.GPicurrent_old	=	0
		self.GPecurrent_old	=	0
		self.D1current_old	=	0
		self.D2current_old	=	0
		self.STNcurrent_old	=	0
		self.FScurrent_old	=	0
		self.Ctxcurrent_old	=	0
		self.Thacurrent_old	=	0
		self.TRNcurrent_old =	0
		self.now			=	0
	
#----------------------------------------------#
	#increment the state of the structure by dt for N steps, set show to 'Y' to display result to user 
	def stepForward(self,inputs,dt=0.001,N=1,show='N'):
		
		#repeat increment N times
		for steps in range(N):
			
			#solve the new synaptic currents using: tau*dJ/dt=-J+v
			self.GPicurrent	=	self.GPicurrent_old+1/self.tau*(inputs[0]-self.GPicurrent_old)*dt
			self.GPecurrent	=	self.GPecurrent_old+1/self.tau*(inputs[1]-self.GPecurrent_old)*dt
			self.D1current	=	self.D1current_old+1/self.tau*(inputs[2]-self.D1current_old)*dt
			self.D2current	=	self.D2current_old+1/self.tau*(inputs[3]-self.D2current_old)*dt
			self.STNcurrent	=	self.STNcurrent_old+1/self.tau*(inputs[4]-self.STNcurrent_old)*dt
			self.FScurrent	=	self.FScurrent_old+1/self.tau*(inputs[5]-self.FScurrent_old)*dt
			self.Ctxcurrent	=	self.Ctxcurrent_old+1/self.tau*(inputs[6]-self.Ctxcurrent_old)*dt
			self.Thacurrent	=	self.Thacurrent_old+1/self.tau*(inputs[7]-self.Thacurrent_old)*dt
			self.TRNcurrent	=	self.TRNcurrent_old+1/self.tau*(inputs[8]-self.TRNcurrent_old)*dt

			#sum all the synaptic currents into the neuron weighted by a gain term G: I=sum(G*j)
			self.synaptic_current=(self.weightGPi*self.GPicurrent+self.weightGPe*self.GPecurrent+self.weightD1*self.D1current
													+self.weightD2*self.D2current+self.weightSTN*self.STNcurrent+self.weightFS*self.FScurrent
													+self.weightCtx*self.Ctxcurrent+self.weightTha*self.Thacurrent+self.weightTRN*self.TRNcurrent+inputs[9])
			
			#set the currents to be used as the n-1 input for the next step forward
			self.GPicurrent_old	=	self.GPicurrent
			self.GPecurrent_old	=	self.GPecurrent
			self.D1current_old	=	self.D1current
			self.D2current_old	=	self.D2current
			self.STNcurrent_old	=	self.STNcurrent
			self.FScurrent_old	=	self.FScurrent
			self.Ctxcurrent_old	=	self.Ctxcurrent
			self.Thacurrent_old	=	self.Thacurrent
			self.TRNcurrent_old	=	self.TRNcurrent
			
			#solve the firing rate for the structure using the synaptic inputs: dv/dt=-v+S[R*I+/-v]
			self.noise		=	numpy.random.random_sample(1)
			self.now	=	self.history[-1]+1/self.tau*(self.sigmoid((self.R*self.synaptic_current+self.ex_in*self.history[-1]))-self.history[-1]+self.noise*0.5)*dt
			
			if self.now>1: self.now=1
			
			if self.now<0: self.now=0
			
			#shift all the stored values back in the history vector
			item_n=0
			for hist_count in self.history:
				if (item_n+1 == len(self.history)):
					self.history[-1]=self.now
				else:
					self.history[item_n]=self.history[item_n+1]
					item_n+=1
			
			#display the current structural firing rate to the user
			if (show == 'Y'):
				print self.history[-1]
		
#----------------------------------------------#
	#allows the user to read out information about the current state of the structure
	def readStruct(self, read_what, delta_t=1):
		if (read_what == 'name'):
			print self.name						#displays the name of the structure to the user when requested
		elif (read_what == 'state'):
			if (delta_t > len(self.history)):	#informs the user they have tried to look to far back in the structures stored history
				print ('Neural mass history only held for '+ str(len(self.history))+' time steps\nlongest history item returned')
				return self.history[0]			#returns the longest history item because requested item is out of range
			else:
				return self.history[-delta_t]	#returns the firing rate delta_t steps back in history
		elif (read_what == 'parameters'):
			print self.variable					#prints the current parameters set for the strucutre 
		else:
			print (self.name+' does not have a property called: '+read_what)
												#informs the user that what they requested doesn't exist
												
#----------------------------------------------#
	#save the currently stored firing rate to file including the timestep used 
	def log(self,dt):
		self.file_name.writelines(str(dt)+' '+str(self.history[-1])+'\n')
	
#----------------------------------------------#
	#calculate the nonlinear (sigmoidal) response of a structure to synaptic inputs	
	def sigmoid(self,u):
		try:									#make sure the exponential doesn't blow up
			self.sigrate=2*self.maxrate/(1+math.exp(self.slope*(self.maxpot-u)))
		except:									#let the user know that the exponential was too large to compute so set sigmoid to 0
			self.sigrate=0
		#		print 'terminating simulation: math overflow error, sigmoid input=' 
		#		print (self.slope*(self.maxpot-u))
		#		exit()
		return self.sigrate						#returns the sigmoidal response for the given input
	
#==============================================#
def main(DA=1, dt=1./24000, sim_time=1.):
#==============================================#
	#this is the main program that prepares each structure, couples them and then simulates their time evolution
	print 'setting up simulation paramaters'
	
				#This matrix represents the coupling between each strucutre, and sets the structures synaptic response
				#	 GPi  GPe  D1  D2   STN FS  Ctx  	   Tha  TRN  tau   [Sigmoid]  R
	GPi_parama	=	[   0,-0.3, -1,   0,0.9,  0,   0	  ,   0,   0,0.020,0.5, 1,  0, 100]
	GPe_parama	=	[   0,   0,  0,  -1,0.9,  0,   0	  ,   0,   0,0.020,0.5, 1,  0, 100]
	D1_parama	=	[   0,1+DA,  0,   0,  0, -1,(1+DA)*0.5,   0,   0,0.020,0.5, 1,  0, 100]
	D2_parama	=	[   0,1-DA,  0,   0,  0, -1,(1-DA)*0.5,   0,   0,0.020,0.5, 1,  0, 100]
	STN_parama	=	[   0,  -1,  0,   0,  0,  0,   0.5	  ,   0,   0,0.005,0.5, 1,  0, 100]
	FS_parama	=	[   0,   1,  0,   0,  0,  0,   0	  ,   0,   0,0.005,0.5, 1,  0, 100]
	Ctx_parama	=	[   0,   0,  0,   0,  0,  0,   0 	  , 0.6,   0,0.080,0.5, 1,  0, 100]
	Tha_parama	=	[0.18,   0,  0,   0,  0,  0,   0.6	  ,   0,0.35,0.005,0.5, 1,  0, 100]
	TRN_parama	=	[   0,   0,  0,   0,  0,  0,   1	  ,0.35,   0,0.005,0.5, 1,  0, 100]
  
	print '\n'*100
	print 'constructing structure objects'
  
	GPi	= CNeuralStruct('GPi',GPi_parama)	#Creates each structure using the neural structure class, including setting their coupling to each other
	GPe	= CNeuralStruct('GPe',GPe_parama)
	D1	= CNeuralStruct('D1',D1_parama)
	D2	= CNeuralStruct('D2',D2_parama)
	STN	= CNeuralStruct('STN',STN_parama)
	FS	= CNeuralStruct('FS',FS_parama)
	Ctx	= CNeuralStruct('Ctx',Ctx_parama)
	Tha	= CNeuralStruct('Tha',Tha_parama)
	TRN	= CNeuralStruct('TRN',TRN_parama)
  
	print 'structure objects constructed'

	GPi.log(0.0)							#logs ths initial state of each structure to file
	GPe.log(0.0)
	D1.log(0.0)
	D2.log(0.0)
	STN.log(0.0)
	FS.log(0.0)
	Ctx.log(0.0)
	Tha.log(0.0)
	TRN.log(0.0)
	
	print 'running simulation for time:'
	print sim_time
	print 'with time step:'
	print dt
	
	total_time=0.0
	while (total_time<sim_time):			#continue incrementing each structure until total time reaches the simulation time
		
		#read the presynaptic current from each structure, including time delays, to be used for the step forward in time
		GPiinputs	=	[GPi.readStruct('state'),GPe.readStruct('state'),D1.readStruct('state'),D2.readStruct('state',5),
							STN.readStruct('state',2),FS.readStruct('state'),Ctx.readStruct('state'),Tha.readStruct('state'),TRN.readStruct('state'),0]
		
		GPeinputs	=	[GPi.readStruct('state'),GPe.readStruct('state'),D1.readStruct('state'),D2.readStruct('state',5),
							STN.readStruct('state',2),FS.readStruct('state'),Ctx.readStruct('state'),Tha.readStruct('state'),TRN.readStruct('state'),0]
		
		D1inputs	=	[GPi.readStruct('state'),GPe.readStruct('state'),D1.readStruct('state'),D2.readStruct('state'),
							STN.readStruct('state'),FS.readStruct('state'),Ctx.readStruct('state',10),Tha.readStruct('state'),TRN.readStruct('state'),0]
		
		D2inputs	=	[GPi.readStruct('state'),GPe.readStruct('state'),D1.readStruct('state'),D2.readStruct('state'),
							STN.readStruct('state'),FS.readStruct('state'),Ctx.readStruct('state',10),Tha.readStruct('state'),TRN.readStruct('state'),0]
		
		STNinputs	=	[GPi.readStruct('state',4),GPe.readStruct('state',4),D1.readStruct('state'),D2.readStruct('state'),
							STN.readStruct('state'),FS.readStruct('state'),Ctx.readStruct('state',3),Tha.readStruct('state'),TRN.readStruct('state'),0]
		
		FSinputs	=	[GPi.readStruct('state'),GPe.readStruct('state'),D1.readStruct('state'),D2.readStruct('state'),
							STN.readStruct('state'),FS.readStruct('state'),Ctx.readStruct('state'),Tha.readStruct('state'),TRN.readStruct('state'),0]
		
		Ctxinputs	=	[GPi.readStruct('state'),GPe.readStruct('state'),D1.readStruct('state'),D2.readStruct('state'),
							STN.readStruct('state'),FS.readStruct('state'),Ctx.readStruct('state'),Tha.readStruct('state'),TRN.readStruct('state'),0]
		
		Thainputs	=	[GPi.readStruct('state'),GPe.readStruct('state'),D1.readStruct('state'),D2.readStruct('state'),
							STN.readStruct('state'),FS.readStruct('state'),Ctx.readStruct('state'),Tha.readStruct('state'),TRN.readStruct('state'),0]
		
		TRNinputs	=	[GPi.readStruct('state'),GPe.readStruct('state'),D1.readStruct('state'),D2.readStruct('state'),
							STN.readStruct('state'),FS.readStruct('state'),Ctx.readStruct('state'),Tha.readStruct('state'),TRN.readStruct('state'),0]
		
		status='N'									#do not display current time step to user unless the current time is a whole second
		#if ((total_time*1000)%1000 < 1):
		#	status='Y'
		#	print (str(total_time)+'s simulated')
		
		#add input currents eg. DBS
		#if ((total_time*1000)%135 < 10):
		#	GPiinputs[9]	= -1	
		#	GPeinputs[9]	= -1
		#	D1inputs[9]		= -1
		#	D2inputs[9]		= -1
		#	STNinputs[9]	= -1
		#	FSinputs[9]		= -1
		#	Ctxinputs[9]	= -1
		#	Thainputs[9]	= -1
		#	TRNinputs[9]	= -1
		
		GPi.stepForward(GPiinputs,dt,1,status)		#increment each structure forward by dt once
		GPe.stepForward(GPeinputs,dt,1,status)
		D1.stepForward(D1inputs,dt,1,status)
		D2.stepForward(D2inputs,dt,1,status)
		STN.stepForward(STNinputs,dt,1,status)
		FS.stepForward(FSinputs,dt,1,status)
		Ctx.stepForward(Ctxinputs,dt,1,status)
		Tha.stepForward(Thainputs,dt,1,status)
		TRN.stepForward(TRNinputs,dt,1,status)
		
		GPi.log(dt)									#log the current firing rate from each structure to file
		GPe.log(dt)
		D1.log(dt)
		D2.log(dt)
		STN.log(dt)
		FS.log(dt)
		Ctx.log(dt)
		Tha.log(dt)
		TRN.log(dt)
		
		total_time+=dt								#increment time by timestep
	
	print 'simulation complete'

	#exit()

#----------------------------------------------#	
if (__name__ == '__main__'):						#only runs the simulation if this file is executed (allows the class CNeuralStruct to be called by an external program
  main()
