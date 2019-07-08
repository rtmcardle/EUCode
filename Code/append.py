import numpy as np
import os

#Stitches together abundance data from before and after star death

################################################################################

#RECORDS WHICH STAR/POINT PAIR IS BEING EVALUATED
point=int(np.loadtxt("./Code/eucode_grid/sloop.data")[0])
star=int(np.loadtxt("./Code/eucode_grid/sloop.data")[1])

#LOADS DATA FOR BEFORE AND AFTER STAR DEATH
livedata=np.loadtxt("./Code/eucode_grid/abund."+str(star).zfill(3)+".txt")
deaddata=np.loadtxt("./Code/eucode_dead/"+str(star).zfill(3)+".txt")

#SAVES COMBINED DATA TO OUTPUT FILE
datafile=open("./Results/Grid/point_"+str(point).zfill(3)+"/abund."+str(star).zfill(3)+".txt", 'w+')
np.savetxt(datafile,livedata, fmt='%.3e')
np.savetxt(datafile,deaddata, fmt='%.3e')
datafile.close()


################################################################################
