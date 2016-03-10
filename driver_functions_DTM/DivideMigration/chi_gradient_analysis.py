import numpy as np, matplotlib.pyplot as plt
from matplotlib import rcParams
from itertools import cycle
from scipy import stats as stats

def read_cosmo_file(file_name):
    print "Reading comso file"    
    f = open(file_name, 'r')
    lines = f.readlines()
    N_data = len(lines)-1
    SampleID = []
    Area = np.zeros(N_data)
    Erate = np.zeros(N_data)
    Error = np.zeros(N_data)
    MeanChiGradient = np.zeros(N_data)
    SD = np.zeros(N_data)
    SErr = np.zeros(N_data)
    
    for i in range(0,N_data):
        data_elements = lines[i+1].strip().split(" ")
        SampleID.append(data_elements[0])
        Area[i] = data_elements[1]
        Erate[i] = data_elements[2]
        Error[i] = data_elements[3]
        MeanChiGradient[i] = data_elements[4]
        SD[i] = data_elements[5]
        SErr[i] = data_elements[6]
    f.close()
    return SampleID,Area,Erate,Error,MeanChiGradient,SD,SErr

def read_histogram_file(file_name):
    print "Reading histogram file"    
    data = np.loadtxt(file_name,skiprows=1)        
    mids = data[:,0]
    llim = data[:,1]
    ulim = data[:,2]
    PD = data[:,3]
    Count = data[:,4]
    return mids,llim,ulim,PD,Count
    

# This function plots validation parameters in three subplots   
def make_plots(cosmo_file,hist1_file,hist2_file):
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 10
    rcParams['legend.fontsize'] = 10.0
    rcParams['legend.numpoints'] = 1
    
    SampleID,Area,Erate,Error,MeanChiGradient,SD,SErr = read_cosmo_file(cosmo_file)
    mids,llim,ulim,PD,Count = read_histogram_file(hist1_file)
    temp,temp,temp,PD_ned,Count_ned = read_histogram_file(hist2_file)

    histogram_x = np.zeros(3*len(mids))    
    histogram_y = np.zeros(3*len(mids))    
    # make a histogram
    for i in range(0,len(mids)):
        histogram_x[3*i+0]=llim[i]
        histogram_x[3*i+1]=mids[i]
        histogram_x[3*i+2]=ulim[i]
        histogram_y[3*i+0]=PD[i]
        histogram_y[3*i+1]=PD[i]
        histogram_y[3*i+2]=PD[i]
        
    
    plt.figure(1, facecolor='white',figsize=(7,3.5))        
    ax1 = plt.subplot2grid((1,2),(0,0))   
    ax2 = plt.subplot2grid((1,2),(0,1))       
    
    ax1.errorbar(MeanChiGradient, Erate, xerr=SD, yerr=Error,fmt=".",color="black",linewidth=0.5)
    ax1.plot(MeanChiGradient, Erate, "o", markersize=5, color='blue', markeredgecolor="blue")

    #slope, intercept, r_value, p_value, std_err = stats.linregress(MeanChiGradient, Erate)
    #print "Linear fit statistics: slope=", slope, "; intercept=", intercept, "; p=", p_value, "; R-squared=", r_value**2, "; std error=", std_err
    exponent, constant, r_value, p_value, std_err = stats.linregress(np.log(MeanChiGradient), np.log(Erate))
    print "Log fit statistics: exponent=", exponent, "; constant=", constant, "; p=", p_value, "; R-squared=", r_value**2, "; std error=", std_err
        
    x_model = np.linspace(0,2.5,1000)
    y_model = np.exp(constant)*x_model**exponent
    
    ax1.plot(x_model,y_model,"-", color="black")

    
    ax2.plot(histogram_x,histogram_y, "-", markersize=5, color='blue')
    plt.fill_between(histogram_x, 0, histogram_y, where=None, color="blue", alpha = 0.2) 
        
    #ax2.plot(mids,PD_ned, "-", markersize=5, color='red')
    ax1.set_xlim(xmin=0)    
    ax2.set_xlim(0,6)
    ax2.set_ylim(0,0.15)
    #ax2.set_ylim(0,0.179)
    
    predicted_erate = np.exp(constant)*(MeanChiGradient**exponent)
    RMSE = np.sqrt(np.mean((Erate-predicted_erate)*(Erate-predicted_erate)))
    print "RMSE: ", RMSE    
    
    # Finishing off plots    
    ax1.annotate('a', xy=(0.05,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='bottom', fontsize=rcParams['font.size']+2) 
    ax1.set_ylabel("CRN Erosion Rate / mm ka$^{-1}$", fontsize=rcParams['font.size']+2)
    ax1.set_xlabel("Chi-gradient", fontsize=rcParams['font.size']+2)
    ax1.grid(True)
    
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.annotate('b', xy=(0.05,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='bottom', fontsize=rcParams['font.size']+2) 
    ax2.set_xlabel("Chi-gradient", fontsize=rcParams['font.size']+2)
    ax2.set_ylabel("Probability Density", fontsize=rcParams['font.size']+2)
    ax2.grid(True)
    
    #stats annotation
    output_stats = r"$R^2 =$ " + '%.3f' % (r_value*r_value) + "\n$RMSE =$ " + '%.1f' % (RMSE) + " mm ka$^{-1}$" + "\n$p-value =$ " + '%.3f' % (p_value)
    ax1.annotate(output_stats, xy=(0.05,0.67), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='center', fontsize=rcParams['font.size']) 
    
    
    print np.mean(SD/MeanChiGradient*100)    
    
    
    plt.tight_layout()    
    #plt.savefig("Cosmo_chi_analysis.pdf",format="pdf")    
    plt.show()    
    
    print "ciao"

if __name__ == "__main__":
    make_plots("fr_5m_cosmo_basin_chi_gradient.txt","fr_5m_chi_gradient_histogram.txt","fr_ned_chi_gradient_histogram.txt")