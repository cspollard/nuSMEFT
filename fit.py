
from  matplotlib import pyplot as plt
import numpy as np
import math
# Function to read data from the provided .dat file
def read_data(file_path):
    data = {"bins": [], "values": []}

    with open(file_path, 'r') as file:
        reading_data = False
        for line in file:
            if line.startswith("# BEGIN HISTO1D"):
                reading_data = True
            elif line.startswith("# END HISTO1D"):
                reading_data = False
            elif reading_data and not line.startswith("#") and not line.startswith("L", 0, 1)  and not line.startswith("T",0,1) and not line.startswith("S", 0,1) and not line.startswith("P",0,1):
                values = [float(val) for val in line.split()]
                data["bins"].append((values[0], values[1]))
                data["values"].append(values[2])

    #print(data)
    return data

# Function to plot histograms
def plot_histograms(data, title, xlabel, ylabel):
    for i, hist_data in enumerate(data):
        bins = [bin[0] for bin in hist_data["bins"]]
        #print(bins)
        values = hist_data["values"]
        MAX = max(values)
        print(bins[values.index(MAX)])

        for i in range(len(values)):
                values[i] = values[i]/MAX
        #plt.hist(bins, bins=bins + [bins[-1]], weights=values, alpha=0.5, label='Histogram {i + 1}')
        width = bins[1] - bins[0]
        #plt.bar(bins,height = values, align='center', width=width, alpha = 0.2)
        plt.step(bins,values, where = "post")

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend("Missing Momentum")

def Overall(file_path):
    # Specify the file path


    # Read data from the file
    hist_data_list = [read_data(file_path)]

    # Plot histograms
    plot_histograms(hist_data_list, 'Histograms', 'X-axis', 'Y-axis')

    make_breit_wigner(hist_data_list, 34.5)

    likelyhood(hist_data_list)

def make_breit_wigner(data,M):
    sigma = []
    lambda_full = 20

    for x, hist_data in enumerate(data):
        bins = [bin[0] for bin in hist_data["bins"]]
        
    
    for i in range(len(bins)):
        E = bins[i]

        sigma.append((lambda_full/2)**2/(((E)-(M))**2+(lambda_full/2)**2))

    return sigma

   

def likelyhood(data):
    L_List = {"Mass":[],"Likelyhood":[]}
    for M in np.arange(1,150,0.2):
       
        sigma=make_breit_wigner(data,M)
        
        for i, hist_data in enumerate(data):
            
            bins = [bin[0] for bin in hist_data["bins"]]
            values = hist_data["values"]
            MAX = max(values)
            for i in range(len(values)):
                values[i] = values[i]/MAX
        L = 0
        for x in range(len(values)):
            n = sigma[x]
            m = values[x]
           
            L = L + (np.exp(-m)*n**(m))/(math.gamma(n+1))
        L_List["Likelyhood"].append(L)
        L_List["Mass"].append(M)
    L_Max = max(L_List["Likelyhood"])
    index = L_List["Likelyhood"].index(L_Max)
    print(L_List["Mass"][index])

    sigma = make_breit_wigner(data,L_List["Mass"][index])

    plt.step(bins,sigma,where = "post")
        




#Overall("/home/barkerj/rivet/plots/new_10gev/MyAnalysis_hist_mT.dat")
Overall("/home/barkerj/rivet/plots/new_20gev/MyAnalysis_hist_mT.dat")
#Overall("/home/barkerj/rivet/plots/new_40gev/MyAnalysis_hist_mT.dat")
#Overall("/home/barkerj/rivet/plots/new_80gev/MyAnalysis_hist_mT.dat")

plt.show() 

