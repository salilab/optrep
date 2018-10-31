
import IMP
import IMP.optrep
import os,sys,string,math
import argparse
import subprocess
import datetime
import time
import scipy.stats
import stats_helper
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib import gridspec

from mpl_toolkits.axes_grid1 import make_axes_locatable

def parse_args():
    
    parser = argparse.ArgumentParser(description="Plot bead maps! Usage: plot_beadmaps.py -s 2IDO -e 9. Flag -h for more details.")
    
    parser.add_argument("-s","--system",dest="system",type=str,help="the name of the system")
    
    parser.add_argument("-e","--experiment",dest="experiment",type=str,help="the experiment number")
    
    result = parser.parse_args()
 
    return result
    
def get_xlinked_residues(xlink_file,proteins_optimized):
    
    xlinked_res= []
    xlf = open(xlink_file,'r')
    for ln in xlf.readlines():
        
        if ln.startswith("res1"):
            continue
            
        [res1,prot1,res2,prot2]=ln.strip().split(',')
        
        if prot1 in proteins_optimized and not int(res1) in xlinked_res: #TODO for multiple proteins ideally xlinked_res is a dictionary keyed by protein name!
            xlinked_res.append(int(res1))
        if prot2 in proteins_optimized and not int(res2) in xlinked_res:
            xlinked_res.append(int(res2))
        
    return sorted(xlinked_res)
    
    xlf.close()
    
def get_beadwise_precisions(precisions_file,proteins_optimized):
    
    beadwise_precisions = {}
    pf=open(precisions_file,'r')
    for ln in pf.readlines():
        fields=ln.strip().split()
        
        if not fields[0] in proteins_optimized:
            continue
        
        beadwise_precisions[int(fields[2])]=float(fields[3])
       
    pf.close()
     
    return beadwise_precisions
    
def generate_precisions_per_residue(precision_file,beadmap_file,proteins_optimized):
    
    beadwise_precisions = get_beadwise_precisions(precision_file,proteins_optimized)
    
    precision_per_residue=[] # contiguous range
    
    residue_range=[] # sequence of residues
    
    bead_boxes=[]
    
    bmf=open(beadmap_file,'r')
    
    bead_count=0
    
    for ln in bmf.readlines():
        
        fields=ln.strip().split()
        if not fields[0] in proteins_optimized:
            continue
        
        bead_start = int(fields[2])
        bead_end = int(fields[3])
        
        for i in range(bead_start,bead_end+1):
            residue_range.append(i)
        
            precision_per_residue.append(beadwise_precisions[bead_count])    
        
        #bead_boxes.append(Rectangle((xrect,0.0),rect_width,1.0))
        bead_boxes.append(bead_end+1)
        
        bead_count+=1

    bmf.close()
    
    precision_per_residue=numpy.expand_dims(numpy.array(precision_per_residue),axis=0)
    
    
    return (precision_per_residue,residue_range,bead_boxes)
                
def plot_beadmaps():
    # Authors: Shruthi Viswanath
     
    # process input
    arg=parse_args()
  
    data_dir = os.path.join(os.path.expanduser('~'),"optrep","expts/expt"+arg.experiment,arg.system)
    
    config_file = os.path.join(os.path.expanduser('~'),"optrep","input",arg.system,arg.system+".config."+arg.experiment)
    
    config_params = stats_helper.parse_config_file(config_file)
    
    xlinked_residues = get_xlinked_residues(os.path.join(os.path.expanduser('~'),"optrep","input",arg.system,config_params["XLINKS_FILE"][0]),config_params["PROTEINS_TO_OPTIMIZE_LIST"])
   
    precision_residues = {}
    bead_boxes={}
    
    min_precision=10000
    max_precision=0
                    
    for ires,resolution in enumerate(config_params["RESOLUTIONS_LIST"]): 
        
        precision_file = os.path.join(data_dir,"r"+resolution,"bead_precisions_"+resolution+".txt")
        
        beadmap_file = os.path.join(data_dir,"r"+resolution,"bead_map_"+resolution+".txt")
        
        precision_residues[resolution],residue_range,bead_boxes[resolution]= generate_precisions_per_residue(precision_file,beadmap_file, config_params["PROTEINS_TO_OPTIMIZE_LIST"])

        min_precision=min(min_precision,min(precision_residues[resolution][0]))
        max_precision=max(max_precision,max(precision_residues[resolution][0]))
    
    # Get simple bar info from xlinks, to show along sequence
    xlink_bars = []
    for r in residue_range:
        if r in xlinked_residues:
            xlink_bars.append(0.5)
        else:
            xlink_bars.append(0.0)
         
         
    # changing the font sizes of axes and so on
    # need to put it in the front of the file
    matplotlib.rcParams["xtick.labelsize"] = 24 
     
    # for getting colorbar  
    #fig = plt.figure(figsize=(15, 4))
    #gs = gridspec.GridSpec(15, 1)

    fig = plt.figure(figsize=(14, 4))
    gs = gridspec.GridSpec(12, 1)

    # first add the xlinks to show where the data is along the sequence
    #ax1 = fig.add_subplot(7,1,1)
    ax1 = fig.add_subplot(gs[0,0])
    ax1.bar(residue_range,xlink_bars,width=1.0,align='edge',color='black')
    ax1.set_xlim(min(residue_range), max(residue_range)+1)
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    # get an empty subplot
    #ax2 = fig.add_subplot(7,1,2)
    ax2  = fig.add_subplot(gs[1,0])
    ax2.set_visible(False)
        
    #ax3 = fig.add_subplot(7,1,3)
    ax3 = fig.add_subplot(gs[2:4,0])
    ax3.imshow(precision_residues['1'],cmap='hot',interpolation='nearest',aspect='auto',extent=[min(residue_range),max(residue_range)+1,0,0.5],norm=matplotlib.colors.Normalize(vmin=min_precision,vmax=max_precision))
    ax3.vlines(bead_boxes['1'],ymin=0,ymax=0.5,colors='black',linestyles='dotted',linewidth=1.5)
    ax3.set_xlim(min(residue_range), max(residue_range)+1)
    ax3.set_xticks([])
    ax3.set_yticks([])
    
    # Make sure vlines, precisions and the x axes are all sync'ed
    
    #normalize precision over full heatmap. Can you compare precisions of different runs? You cannot but the numerical values need to be compared, so yes. Make it uniform. 
    
    #ax4 = fig.add_subplot(7,1,4)
    ax4 = fig.add_subplot(gs[4:6,0])
    ax4.imshow(precision_residues['5'],cmap='hot',interpolation='nearest',aspect='auto',extent=[min(residue_range),max(residue_range)+1,0,0.5],norm=matplotlib.colors.Normalize(vmin=min_precision,vmax=max_precision))
    ax4.vlines(bead_boxes['5'],ymin=0,ymax=0.5,colors='black',linestyles='dotted',linewidth=1.5)
    ax4.set_xlim(min(residue_range), max(residue_range)+1)
    ax4.set_xticks([])
    ax4.set_yticks([])
    
    #ax5 = fig.add_subplot(7,1,5)
    ax5 = fig.add_subplot(gs[6:8,0])
    ax5.imshow(precision_residues['10'],cmap='hot',interpolation='nearest',aspect='auto',extent=[min(residue_range),max(residue_range)+1,0,0.5],norm=matplotlib.colors.Normalize(vmin=min_precision,vmax=max_precision))
    ax5.vlines(bead_boxes['10'],ymin=0,ymax=0.5,colors='black',linestyles='dotted',linewidth=1.5)
    ax5.set_xlim(min(residue_range), max(residue_range)+1)
    ax5.set_xticks([])
    ax5.set_yticks([])
  
   
    #ax6 = fig.add_subplot(7,1,6)
    ax6 = fig.add_subplot(gs[8:10,0])
    ax6.imshow(precision_residues['20'],cmap='hot',interpolation='nearest',aspect='auto',extent=[min(residue_range),max(residue_range)+1,0,0.5],norm=matplotlib.colors.Normalize(vmin=min_precision,vmax=max_precision))
    ax6.vlines(bead_boxes['20'],ymin=0,ymax=0.5,colors='black',linestyles='dotted',linewidth=1.5)
    ax6.set_xlim(min(residue_range), max(residue_range)+1)
    ax6.set_xticks([])
    ax6.set_yticks([])
   
    #ax7 = fig.add_subplot(7,1,7)    
    ax7 = fig.add_subplot(gs[10:12,0])
    color_image = ax7.imshow(precision_residues['30'],cmap='hot',interpolation='nearest',aspect='auto',extent=[min(residue_range),max(residue_range)+1,0,0.5],norm=matplotlib.colors.Normalize(vmin=min_precision,vmax=max_precision))
    ax7.vlines(bead_boxes['30'],ymin=0,ymax=0.5,colors='black',linestyles='dotted',linewidth=1.5)
    ax7.set_xlim(min(residue_range), max(residue_range)+1)
    ax7.set_xticks(numpy.arange(min(residue_range), max(residue_range)+1, 20))
    ax7.set_yticks([])
    
    # this is for getting the color bar only
    ## empty subplot
    #ax8  = fig.add_subplot(gs[12:14,0])
    #ax8.set_visible(False)
    
    #ax9 = fig.add_subplot(gs[14:15,0])
    #cbar = fig.colorbar(color_image, cax=ax9, ticks=[min_precision, (min_precision+max_precision)/2.0, max_precision], orientation='horizontal',fraction=1.5,shrink=0.5,pad=1.0,aspect=100)
    #cbar.set_ticklabels(['High', 'Medium', 'Low'])  # horizontal colorbar
    
    plt.tight_layout()
   
    # remove the space between subplots
    plt.subplots_adjust(wspace=0.01,hspace=0.01) 
    # note that this must be called AFTER any call to tight_layout
    
    plt.savefig("beadheatmap_"+arg.system+".pdf",dpi=600)
      
if __name__ == "__main__" :
    plot_beadmaps()
