"""
Created by Akshay Maheshwari
09/05/2017

Produces analysis figures from experiment data
"""
from simanalysis_methods import *
import matplotlib.pyplot as plt
import time;

start_time=time.time()
expt_name = "171018_2219"
outputlist = loadOutputList(expt_name,'molpos')
histlistpklpath = combinePkls(expt_name,outputlist,covertime=True)
#histlistpklpath = saveHist(outputlist, expt_name,bins=10,diameter=0.1,molposTS=1e-7)
fig = plotHist(histlistpklpath,expt_name,diameter=0.1, graphs="all", logscale=False,step=1,start=1.25,simtime=1,x_label="R_crowder (nm)")
fig.suptitle("Effects of crowding molecule size on covertime, and dispersion of a single tracked molecule. \n[1s. sim] -- R_tracked=7.25nm -- R_crowder=[1.25nm,2.25nm,...9.25nm] -- $\phi$=0.25 -- time step=1e-7s.")
plt.savefig("data/"+expt_name+"/"+expt_name+"_analysis1.png")
print("--- %s seconds ---" % (time.time() - start_time))
