"""
Created by Akshay Maheshwari
09/05/2017

Computes individual subexperiment analyses on a server (Map part of MapReduce)
"""

import glob
import os
import pickle as pkl
import numpy as np
import pandas as pd
import sys
from simanalysis_methods import *


# def molpos_1Dbin(data,bins,diameter):
#     """Creates a 1D histogram from X,Y location data of a single tracked molecule over time"""

#     pos=np.linspace(-diameter/2,diameter/2,bins+1)

#     speciesNum = int(data.shape[1]/2)
#     data_hist1D_total = np.zeros(bins)
#     print(speciesNum)
#     for i in range(speciesNum):
#         data_hist1D = np.histogram(data.loc[:,1],bins=pos)[0]
#         data_hist1D_total += data_hist1D
#     return data_hist1D_total
#     #for chunk in pd.read_csv(datapath,header=None,chunksize=10**6): #Chunk size 10^6 runtime: 72.68, chunk size 10^7 runtime:  72.29
#         #data_hist += np.histogram(chunk.loc[:,1],bins=pos)[0]

# def molpos_2Dbin(data, bins, diameter):
#     """Creates a 2D histogram from X,Y location data of a single tracked molecule over time"""

#     pos_x = np.linspace(-diameter/2,diameter/2,bins+1)
#     pos_y = np.linspace(-diameter/2,diameter/2,bins+1)
#     #for chunk in pd.read_csv(datapath,header=None,chunksize=10**6):
#     speciesNum = int(data.shape[1]/2)
#     data_hist2D_total = np.zeros((bins,bins))
#     for i in range(speciesNum):
#         data_hist2D = np.histogram2d(np.array(data.loc[:,1]),np.array(data.loc[:,2]),bins=[pos_x,pos_y])[0]
#         data_hist2D=np.flip(data_hist2D.T,0)
#         data_hist2D_total += data_hist2D
#     return data_hist2D_total

# def molpos2D_dispersion(data,diameter):
#     """Computes the total variation distance between the X,Y distribution of a tracked molecule over time &
#     the uniform distribution (expected over time in a dilute system). 
#     Returns: [0,1]. 0 means distribution is uniform. 1 is approached as tracked molecule stays in one place."""

#     data_hist2D=molpos_2Dbin(data,bins=10,diameter=diameter)
#     normalized_data=data_hist2D/sum(sum(data_hist2D))
#     newdata = np.abs(normalized_data-(sum(sum(normalized_data)))/88)
#     return 0.5*(sum(newdata[0][2:-2])+sum(newdata[-1][2:-2])+sum(newdata[2][1:-1])+sum(newdata[-2][1:-1])+sum(sum(newdata[2:-2][:])))

# def covertime(datapath,diameter,molposTS):
#     """Finds the earliest time point at which every 2D bin in the 2D area of the cell has been traversed
#     by the tracked molecule. Bins is hardcoded to 10x10 to account for circle edges in 2D array (need to
#     manually remove edges of 10x10 grid that don't fall in cell circle). MolposTS represents the data sampling
#     rate (independent from the underlying simulation time step) """

#     bins=10
#     data_hist = np.zeros((bins, bins))
#     pos_x = np.linspace(-diameter/2,diameter/2,bins+1)
#     pos_y = np.linspace(-diameter/2,diameter/2,bins+1)
#     coverTime = 10
#     foundCoverTime = False;
#     timesteps=0;
#     chunksize=10**3
#     for chunk in pd.read_csv(datapath,header=None,chunksize=chunksize):
#         timesteps+=1;
#         data_hist+=molpos_2Dbin(chunk,bins=bins,diameter=diameter) #np.histogram2d(np.array(chunk.loc[:,1]),np.array(chunk.loc[:,2]),bins=[pos_x,pos_y])[0]
#         if not 0 in data_hist[:,2:-2] and not 0 in data_hist[1:-1,1:-1] and not 0 in data_hist[2:-2,:] and not foundCoverTime:
#             coverTime = timesteps
#             foundCoverTime = True
#             return (molposTS*chunksize)*coverTime
#     return coverTime

            

bins=sys.argv[1]
radius=sys.argv[2]
molposTS=sys.argv[3]

molpospath=glob.glob('expts/data/*Molpos*')
molpos2D = pd.read_csv(molpospath[0], header=None).loc[:,1:3]

data_hist1D=molpos_1Dbin(molpos2D,bins=int(bins),diameter=float(radius)*2)
data_hist2D=molpos_2Dbin(molpos2D,bins=int(bins),diameter=float(radius)*2)
covtime=covertime(molpospath[0],diameter=float(radius)*2,molposTS=float(molposTS))
molpos2D_disperse = molpos2D_dispersion(molpos2D,diameter=float(radius)*2)

if not os.path.exists('expts/data/analysis'):
    os.mkdir('expts/data/analysis/')
path='expts/data/analysis/'+os.path.basename(molpospath[0]).split('.')[0]+'Histogram.pkl'
outputFile = open(path,"wb")
pkl.dump([data_hist1D,data_hist2D,covtime,molpos2D_disperse],outputFile)
