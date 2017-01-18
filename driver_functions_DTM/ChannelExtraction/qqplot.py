import numpy as np, matplotlib.pyplot as plt
from matplotlib import rcParams

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1

def read_q_q_file(file_name):
    f = open(file_name, 'r')
    lines = f.readlines()
    N_data = len(lines)-1
    
    # Initialise vectors
    quantiles = np.zeros(N_data)
    values = np.zeros(N_data)
    mn_values = np.zeros(N_data)
    # Load in data
    for i in range (0, N_data):        
       line = lines[i+1].strip().split(" ")
       quantiles[i]=float(line[0])
       values[i]=float(line[1])
       mn_values[i]=float(line[2])
    f.close()
    return quantiles,values,mn_values
    
def make_q_q_plot(snv,values,mn_values):

    flag = 0
    for i in range(0,len(snv)):
        if (snv[i] >= 0):
            if (values[i]*0.95>mn_values[i]):
                if (flag == 0):
                    flag = 1
                    x_thresh = snv[i]
            else:
                flag = 0
    plt.figure(1, facecolor='White',figsize=[4,4])
    ax1 = plt.subplot(1,1,1)            
    ax1.plot(snv,values,linewidth=2,color="blue",label="data")
    ax1.plot(snv,mn_values,"--",linewidth=2,color="red",label="normal distribution")
    #ax1.axvline(x=x_thresh,linestyle='--',color='black')
    xmin,xmax = ax1.get_xlim()       
    ax1.axvspan(x_thresh, xmax, alpha = 0.2, color='blue') 
    ax1.legend(loc = 2)
    ax1.set_xlabel("Standard Normal Variate", fontsize=rcParams['font.size']+2)
    ax1.set_ylabel("Curvature $(m^{-1})$", fontsize=rcParams['font.size']+2)
    ax1.annotate('Channels', xy=(0.9,0.1), xycoords='axes fraction',color='black',horizontalalignment='right', verticalalignment='center', fontsize=rcParams['font.size']+2) 
    ax1.annotate('Hillslopes', xy=(0.2,0.1), xycoords='axes fraction',color='black',horizontalalignment='left', verticalalignment='center', fontsize=rcParams['font.size']+2) 
    ax1.set_xlim(xmin,xmax)
    ax1.grid(True)
    plt.tight_layout()
    
if __name__ == "__main__":
    file = "K:/data/SurfaceRoughness/RayleighPeak/Dataset2/co_qq_.txt"
    x,y1,y2=read_q_q_file(file)
    plt.show()
    make_q_q_plot(x,y1,y2)
    plt.show()