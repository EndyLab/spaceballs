"""
Created by Akshay Maheshwari
09/05/2017

Simulation analysis methods for 2D hard disk systems
"""

def loadOutputList(expt_name,outputType):
    """Loads a file containing all file names with a specified output data into a local list
    
    Args:
        expt_name (String): Name of experiment (which contains how many ever simulation output files)
        outputType (String): molpos or simtime. Determines which set of output filenames to load.
    
    Returns:
        TYPE: List
    """

    datalist = []
    if(outputType=='molpos'):
        path='data/'+expt_name+'/outputMolposList.txt'
    elif(outputType=='simtime'):
        path='data/'+expt_name+'/outputSimtimeList.txt'
    else:
        raise Exception('outputType required to be either \'molpos\' or \'simtime\'') 
    with open(path) as f:
        for outputfile in f:
            datalist.append("data/"+expt_name+"/"+outputfile.rstrip())
    return datalist 

def molpos_1Dbin(data,bins,diameter):
    """Creates a 1D histogram from X,Y location data of a single tracked molecule over time
    
    Args:
        data (pandas dataframe): time series 2D location data of a tracked molecule
        bins (int): # of rectangular bins to discretize cell with (bins are equidistant on x-axis)
        diameter (float): diameter of simulated cell
    
    Returns:
        TYPE: List
    """
    import numpy as np;
    pos=np.linspace(-diameter/2,diameter/2,bins+1)

    speciesNum = int(data.shape[1]/2)
    data_hist1D_total = np.zeros(bins)
    print(speciesNum)
    for i in range(speciesNum):
        data_hist1D = np.histogram(data.loc[:,1],bins=pos)[0]
        data_hist1D_total += data_hist1D
    return data_hist1D_total
    #for chunk in pd.read_csv(datapath,header=None,chunksize=10**6): #Chunk size 10^6 runtime: 72.68, chunk size 10^7 runtime:  72.29
        #data_hist += np.histogram(chunk.loc[:,1],bins=pos)[0]

def molpos_2Dbin(data, bins, diameter):
    """Creates a 2D histogram from X,Y location data of a single tracked molecule over time

    Args:
        data (pandas dataframe): time series 2D location data of a tracked molecule
        bins (int): # of bins to discretize cell with (bins x bins on x and y)
        diameter (float): diameter of simulated cell
    
    Returns:
        TYPE: List
    """

    import numpy as np;
    pos_x = np.linspace(-diameter/2,diameter/2,bins+1)
    pos_y = np.linspace(-diameter/2,diameter/2,bins+1)
    #for chunk in pd.read_csv(datapath,header=None,chunksize=10**6):
    speciesNum = int(data.shape[1]/2)
    data_hist2D_total = np.zeros((bins,bins))
    for i in range(speciesNum):
        data_hist2D = np.histogram2d(np.array(data.loc[:,1]),np.array(data.loc[:,2]),bins=[pos_x,pos_y])[0]
        data_hist2D=np.flip(data_hist2D.T,0)
        data_hist2D_total += data_hist2D
    return data_hist2D_total

def molpos2D_dispersion(data,diameter):
    """Computes the total variation distance between the X,Y distribution of a tracked molecule over time &
    the uniform distribution (expected over time in a dilute system). Bins is hardcoded to 10x10 to account 
    for circle edges in 2D array (need to manually remove edges of 10x10 grid that don't fall in cell circle).

    Args:
        data (pandas dataframe): time series 2D location data of a tracked molecule
        diameter (float): diameter of simulated cell
    
    Returns:[0,1]. 0 means distribution is uniform. 1 is approached as tracked molecule stays in one place.
        TYPE: float 
    """

    import numpy as np;
    bins=10
    circle_bins = 88 # of bins in a 10x10 grid that a circle would fall into.
    data_hist2D=molpos_2Dbin(data,bins=bins,diameter=diameter)
    normalized_data=data_hist2D/sum(sum(data_hist2D))
    newdata = np.abs(normalized_data-(sum(sum(normalized_data)))/circle_bins)
    return 0.5*(sum(newdata[0][2:-2])+sum(newdata[-1][2:-2])+sum(newdata[2][1:-1])+sum(newdata[-2][1:-1])+sum(sum(newdata[2:-2][:])))
    #return np.sqrt(sum(data.std(0)**2))
    #return sum(np.abs(data.skew()))

def covertime(datapath,diameter,molposTS):
    """Finds the earliest time point at which every 2D bin in the 2D area of the cell has been traversed
    by the tracked molecule. Bins is hardcoded to 10x10 to account for circle edges in 2D array (need to
    manually remove edges of 10x10 grid that don't fall in cell circle).

    Args:
        datapath (String): File path to datafile
        diameter (float): diameter of simulated cell        
        molposTS (float):  the data sampling rate (independent from the underlying simulation time step) 
    
    Returns:[0,molposTS*total timesteps].
        TYPE: float 
    """

    import numpy as np;
    import pandas as pd;

    bins=10
    data_hist = np.zeros((bins, bins))
    pos_x = np.linspace(-diameter/2,diameter/2,bins+1)
    pos_y = np.linspace(-diameter/2,diameter/2,bins+1)
    coverTime = 10
    foundCoverTime = False;
    timesteps=0;
    chunksize=10**3
    for chunk in pd.read_csv(datapath,header=None,chunksize=chunksize):
        timesteps+=1;
        data_hist+=molpos_2Dbin(chunk,bins=bins,diameter=diameter) #np.histogram2d(np.array(chunk.loc[:,1]),np.array(chunk.loc[:,2]),bins=[pos_x,pos_y])[0]
        if not 0 in data_hist[:,2:-2] and not 0 in data_hist[1:-1,1:-1] and not 0 in data_hist[2:-2,:] and not foundCoverTime:
            coverTime = timesteps
            foundCoverTime = True
            return (molposTS*chunksize)*coverTime
    return coverTime

def timeplot(ax, expt_name,step=1,scalefactor=1,logscale=False,start=0):
    """Creates a plot of runtime (y-axis) vs. whatever dependent variable is being experimentally swept (x-axis).
    
    Args:
        ax (TYPE): figure object to plot upon
        expt_name (String): Name of experiment file
        step (int, optional): Step size between sweep of dependent variable
        scalefactor (int, optional): Used to scale time data
        logscale (bool, optional): True if want plot to have a x logscale
        start (int, optional): 1st value of x-axis to begin plot with
    """
    import warnings
    import numpy as np

    simdata = []
    outputlist= loadOutputList(expt_name,'simtime')[:] #[a:b] = a to b sub-range of experiments that need to be plotted.

    for outputfile in outputlist:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simdata.append(np.loadtxt(open(outputfile), delimiter=",")/scalefactor)
    simdata = np.array(simdata)
    runtimes = simdata[:,1]
    simtime = simdata[-1,0]
    print("plotting ", runtimes)
    
    if(logscale):
        num_varx = np.logspace(start,start+len(runtimes)-1,num=len(runtimes))
        ax.set_xscale("log")
    else:
        num_varx = np.arange(start,start+len(runtimes)*step,step)

    ax.plot(num_varx,runtimes,'-o')
    ax.set_xlim(0)

    #Fits a quadratic curve to time data (note in order for poly1d to fit, need at least 2 data points for runtime):
    #ax.plot(np.unique(num_varx), np.poly1d(np.polyfit(num_varx, runtimes, 2))(np.unique(num_varx)))

def saveHist(outputlist,expt_name, bins, diameter,molposTS):
    """Creates a pkl dump of a list of each subexperiment filename with its associated data_hist1D, data_hist2D, covertime, and molpos2D_disperse
    
    Args:
        outputlist (list): List of all data file names for an experiment
        expt_name (string): Name of experiment
        bins (int): # of bins to compute 1D and 2D histograms with
        diameter (float): diameter of simulated cell        
        molposTS (float):  the data sampling rate (independent from the underlying simulation time step) 
    
    Returns: A oath to the pkl dump
        TYPE: string
    """
    import pickle as pkl
    import pandas as pd
    import os

    if not os.path.exists('data/'+expt_name+'/analysis/'):
        os.makedirs('data/'+expt_name+'/analysis/')
    path='data/'+expt_name+'/analysis/outputMolposHistogramList.pkl'
    outputFile = open(path,"wb")
    histlist = []
    for i in range(len(outputlist)):
        molpos2D = pd.read_csv(outputlist[i],header=None).loc[:,1:]
        data_hist1D=molpos_1Dbin(molpos2D,bins=bins,diameter=diameter)
        data_hist2D=molpos_2Dbin(molpos2D,bins=bins,diameter=diameter)
        covtime=covertime(outputlist[i],diameter=diameter,molposTS=molposTS)
        molpos2D_disperse = molpos2D_dispersion(molpos2D,diameter=diameter)
        print("saving",data_hist1D)
        histlist.append([outputlist[i], data_hist1D,data_hist2D,covtime,molpos2D_disperse])

    pkl.dump(histlist, outputFile)
    return path;

def plotHist(histlistpklpath, expt_name, diameter=0.1, graphs="both",logscale=True,step=1,start=0, simtime=1, x_label="# time samples"):
    """Generates a figure with (for each subexperiment): 
       plots of 2D molecule position heatmaps, simulation time, cover time, and total variation distance

    
    Args:
        histlistpklpath (string): Path to pkl with a list containing a list for each subexperiment containing all analysis (hist1D, hist2D, covertime, molpos2D)
        expt_name (string): Name of experiment
        diameter (float): diameter of simulated cell        
        graphs (str, optional): 'all', 'molpos', 'simtime', 'molpossim' -- determines which graphs included in analysis figure
        logscale (bool, optional): True if want x_axis to logscale for all data plots
        step (int, optional): Step size between sweep of dependent variable
        start (int, optional): 1st value of x-axis to begin plot with
        simtime (int, optional): Total simtime to set y-axis scale for covertime
        x_label (str, optional): X-label for all sub-graphs
    
    Returns:
        TYPE: matplotlib figure
    """

    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import pickle as pkl
    from matplotlib.ticker import ScalarFormatter
    import numpy as np
    pkl_file=open(histlistpklpath,'rb')
    histlist=pkl.load(pkl_file)
    plot_dim = int(np.ceil(len(histlist)/3))
    #fig = plt.figure(figsize=(12,4))
    fig = plt.figure(figsize=(10,15))
    outer=gridspec.GridSpec(4,1,height_ratios=[3,1,1,1])
    inner = gridspec.GridSpecFromSubplotSpec(plot_dim,3,subplot_spec=outer[0],wspace=0.2,hspace=0.25)
    xfmt=ScalarFormatter()
    xfmt.set_powerlimits((-1,1))

    if graphs=="all" or graphs=="molpos" or graphs=="molpossim":
        ###### Plot a heatmap of particle location for each subexperiment #######

       for i in range(plot_dim):
           for j in range(3):
               if len(histlist)-1>=i*3+j:
                   ax = plt.Subplot(fig,inner[i,j])
                   print("plotting", histlist[i*3+j][1])
                   data_hist=histlist[i*3+j][1]
                   pos=np.linspace(-diameter/2,diameter/2,len(data_hist)+1)
                   ax.bar(pos[0:-1], data_hist,width=diameter/(len(histlist[i*3+j][1])),align='edge')
                   ax.imshow(histlist[i*3+j][2].T)
                   #ax.yaxis.set_major_formatter(xfmt)
                   fig.add_subplot(ax)
#    np.set_printoptions(threshold=np.inf)

    if graphs=="all" or graphs=="simtime" or graphs=="molpossim":
        ###### Plot of simtime for each subexperiment #######

       ax2 = plt.Subplot(fig, outer[1])
       timeplot(ax2,expt_name, logscale=logscale,step=step,start=start)
       ax2.set_xlabel(x_label)
       ax2.set_ylabel('Runtime (s.)')
       fig.add_subplot(ax2)

    
    if graphs=="all":
        ####### Plot of covertime for each subexperiment #######
        ax3 = plt.Subplot(fig,outer[2]) 
        covertimearr = [item[3] for item in histlist]
        if(logscale):
            num_varx = np.logspace(start,start+len(covertimearr)-1,num=len(covertimearr))
            ax.set_xscale("log")
        else:
            num_varx = np.arange(start,start+len(covertimearr)*step,step)
        print("covertimearr", covertimearr)
        ax3.plot(num_varx,covertimearr,'-o')
        ax3.set_ylim(0,simtime)
        #ax3.set_yscale('log')
        fig.add_subplot(ax3)
        ax3.set_xlim(0)
        ax3.set_xlabel(x_label)
        ax3.set_ylabel('Cover time (s.)')

        ####### Plot of total variation distance for each subexperiment #######
        ax4 = plt.Subplot(fig,outer[3])
        stdarr = [item[4] for item in histlist]
        if(logscale):
            num_varx = np.logspace(start,start+len(stdarr)-1,num=len(stdarr))
            ax.set_xscale("log")
        else:
            num_varx = np.arange(start,start+len(stdarr)*step,step)
        print("stdarr, ", stdarr)
        ax4.plot(num_varx,stdarr,'-o')
        fig.add_subplot(ax4)
        ax4.set_xlim(0)
        ax4.set_xlabel(x_label)
        ax4.set_ylabel('Diffusive irregularity')

        outer.tight_layout(fig, rect=[0,0.03,1,0.90]) #rect args needed to leave space on top for title
    return fig

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    import re
    return [tryint(c) for c in re.split('([0-9]+)', s)]

def combinePkls(expt_name,outputlist,covertime):
    """Combines pkls for each subexperiment into a list of lists for each subexperiment with: 
    its associated filename, data_hist1D, data_hist2D, covertime, and molpos2D_disperse.data_hist2D
    (Reduce part of MapReduce)
    
    Args:
        expt_name (TYPE): folder name of experiment
        outputlist (TYPE): List with molpos output data subexperiment names
        covertime (boolean): If the experiment computed covertime, include it in the pkl [used for experiments performed before cover time was implemented]
    
    Returns: file path to combined pkl
        TYPE: String
    """
    import os
    import pickle as pkl

    histlist = []
    path='data/'+expt_name+'/analysis/'
    i=0;
    print("i'm here")
    for f in sorted(os.listdir(path),key=alphanum_key)[:]: #can put range here to only plot subset of experiments
        if not f.startswith('.') and f.startswith('expt'):
            data_hist_path = open(path+'/'+f,'rb')
            data_hist = pkl.load(data_hist_path)
            print("test",data_hist)
            if(covertime):
                histlist.append([outputlist[i],data_hist[0],data_hist[1],data_hist[2],data_hist[3]])
            else:
                histlist.append([outputlist[i],data_hist])
            i+=1;
    outputFile = open(path+'/outputMolposHistogramList.pkl',"wb")
    pkl.dump(histlist,outputFile)
    return 'data/'+expt_name+'/analysis/outputMolposHistogramList.pkl'

