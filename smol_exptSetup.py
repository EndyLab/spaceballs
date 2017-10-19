"""
Created by Akshay Maheshwari
09/05/2017

Produces smoldyn parameter specification files that sweep parameters in a specified way.
"""

import os
import datetime
import numpy as np

_cellRAD_ = 0.1 # diameter 100 nm
phi = 0.25
_mol1RAD_ = 0.00165 #Typical volume of a protein assumed to be sphere w/ diameter 2.5nm (GFP=2.5nm, tRNA=2nm, ~8µm^2/s)
_mol1Num_ = 100
_Dc1_ = 1.2 #Diffusion coefficient = 10µm^2/s

_mol2RAD_ = 0.00125
_mol2Num_ = 1
_mol2Position_ = "0 0"
_Dc2_ = 0.12 

_bins_ = 40
_molposTS_ = 1e-6
 
_simtime_ = 0.01
_ts_ = 1e-8


#_simtimeTS_= _simtime_
exptList = open("expts/expt_list.txt","w")
outputSimtimeList=open("expts/outputSimtimeList.txt","w")
outputMolposList=open("expts/outputMolposList.txt","w")
expt_description = open("expts/expt_description.txt","w")

j_max=1
for i in range(0,1):
	#_mol1RAD_=0.00125
	#if(i==1):
	#	_mol2RAD_= 0.01
	#_mol1Num_=int(np.ceil((phi*_cellRAD_**2-_mol2RAD_**2)/_mol1RAD_**2))
	#_mol2Num_+=1


	for j in range(0, j_max):
		exptname = "expt_"+str(j_max*i+j)+".txt"
		_outputSimtime_ = "expt-"+str(j_max*i+j)+"-Simtime-"+datetime.date.today().strftime('%Y%m%d')+ ".csv"
		_outputMolpos_ = "expt-"+str(j_max*i+j)+"-Molpos-"+datetime.date.today().strftime('%Y%m%d')+ ".csv"

		expt = open("expts/" + exptname,'w')
		expt.write("variable _cellRAD_ " + str(_cellRAD_) + "\n")
		
		expt.write("variable _mol1RAD_ " + str(_mol1RAD_)+ "\n")
		expt.write("variable _mol1Num_ " + str(_mol1Num_)+ "\n")
		expt.write("variable _Dc1_ " + str(_Dc1_) + "\n")

		expt.write("variable _mol2RAD_ " + str(_mol2RAD_)+ "\n")
		expt.write("variable _mol2Num_ " + str(_mol2Num_)+ "\n")
		expt.write("define_global _mol2Position_ " + str(_mol2Position_) +"\n")
		expt.write("variable _Dc2_ " + str(_Dc2_) + "\n")

		expt.write("variable _ts_ " + str(_ts_)+ "\n")
		expt.write("variable _bins_ " + str(_bins_)+ "\n")

		expt.write("define_global _simtime_ " + str(_simtime_)+ "\n")
		#expt.write("define_global _simtimeTS_ " + str(_simtimeTS_)+ "\n")
		expt.write("define_global _outputSimtime_ " + "data/"+_outputSimtime_ +"\n")

		expt.write("define_global _molposTS_ " + str(_molposTS_) + "\n")
		expt.write("define_global _outputMolpos_ " + "data/"+_outputMolpos_)

		####### Variables to sweep########
		#_mol2RAD_ += 0.0005
		#_mol1Num_ = int(400 - np.ceil(((_mol2RAD_*1000)**2)/(1.25**2)))
		#_mol1RAD_=_mol1RAD_+0.001
		#_mol1Num_=int(np.ceil((phi*_cellRAD_**2-_mol2RAD_**2)/_mol1RAD_**2))
		_mol1Num_=_mol1Num_+100
		####### Write list of experiment and output
		exptList.write(exptname + "\n")
		outputSimtimeList.write(_outputSimtime_ + "\n")
		outputMolposList.write(_outputMolpos_ + "\n")


expt_description.write("Sampling of "+ str(_mol2Num_) +"molecule position (molpos) in" + str(_mol2Num_)+" molecules (monodisperse)"+" \n")
expt_description.write("Total simulation time: "+ str(_simtime_) + "s. ; Simulation time step: " + str(_ts_)+" s.")