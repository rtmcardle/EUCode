import numpy as np
import os

#Averages over abundance data of entire halo for final data output

################################################################################

#COUNTS STARS/GRID POINTS

grid=len(os.listdir("./Results/Grid/"))-1
stars=len(np.loadtxt("./Results/Stars/stardata.txt"))
weight=[np.loadtxt("./Results/Grid/point_"+str(i).zfill(3)+"/point.data")[0] for i in range(grid)]


################################################################################

#AVERAGES REIONIZATION DATA FOR GRID POINTS

#grid point loop
for i in range(grid):
    k = 0
    nearby=len(os.listdir("./Results/Grid/point_"+str(i).zfill(3)))-1
    dists=np.loadtxt("./Results/Grid/point_"+str(i).zfill(3)+"/point.data")[1:]
    #star loop
    if (nearby!=0) and (weight[i]!=0): #if stars nearby, average abundances
        for j in range(stars):
            if dists[j]!=0:
                textdata = np.loadtxt("./Results/Grid/point_"+str(i).zfill(3)+"/abund."+str(j).zfill(3)+".txt")
                tdata = textdata[:,1:].copy()
                sweight = np.full(tdata.shape,1.0/nearby)
                sdata = np.multiply(tdata,sweight)
                if k == 0:
                   olddata = sdata
                   k+=1
                else:
                    data = [[olddata[m][n] + sdata[m][n] for n in range(len(olddata[0]))] for m in range(len(olddata))]
                    olddata=data
    else: #if no stars nearby, use background abundance data
        textdata = np.loadtxt("./Code/eucode_background/fort.16")
        data = textdata[:,1:].copy()
    datafile=open("./Results/Grid/point_"+str(i).zfill(3)+"/abund.avg", 'w+')
    np.savetxt(datafile, data, fmt='%.3e')
    datafile.close()
    gweight = np.full(np.array(data).shape, weight[i]/np.sum(weight))
    gdata = np.multiply(data,gweight)
    if i == 0:
       totdata = gdata
    else:
        totdata = totdata + gdata


################################################################################

#PREPARES DATA

#recombination data included for dimensional consistency of abundance files
recombdata = np.loadtxt("./Code/eucode_recomb/fort.16")

#adds z values back to reionization abundances data
zdata=textdata[:,0].copy()
z=np.reshape(zdata, (len(zdata),1))
iondata = np.append(z,totdata,1)


################################################################################

#OUTPUTS DATA

datafile = open("./Results/grid.avg.data", 'w+')
np.savetxt(datafile, recombdata, fmt='%.3e')
np.savetxt(datafile, iondata, fmt='%.3e')
datafile.close()


################################################################################
