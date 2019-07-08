import numpy as np
import math as m
import fileinput
import os

#Creates the grid points to be evaluated in the reionization code

################################################################################

#FUNCTIONS

#calculates the distance between two spherical coordinates
def distance(v,w):
    return m.sqrt(v[0]**2+w[0]**2-2*v[0]*w[0]*(m.sin(v[1]*m.sin(w[1])
        *m.cos(v[2]-w[2])+m.cos(v[1])*m.cos(w[1]))))

#calculates the spherical volume weight for a spherical element at coordinate
def sweight(v):
    r=v[0]
    t=v[1]
    p=v[2]
    dr=rmin
    dt=m.acos((2.0*((tmax-tmin)/tsteps))-1.0)
    dp=2*m.pi*((phmax-phmin)/phsteps)
    return np.array([r**2*m.sin(t)*dr*dt*dp])


################################################################################

#LOADS HALO INSTANCE DATA

data=np.loadtxt("./Results/Stars/stardata.txt")
strom=np.loadtxt("./Results/Stars/frontdata.txt")
try:
    stars=data[:,6:]
except:
    stars=data[0,6:]
halo=[(stars[x][0]+strom[x]) for x in range(len(stars))]


################################################################################

#DISTRIBUTES GRID POINTS

#radius range
rsteps=5
rmax=max(halo)
rmin=rmax/rsteps
rs=np.linspace(0,rmax,rsteps)

#theta range
tmin=0.0
tmax=1.0
tsteps=8
ts = [m.acos((2.0*t)-1.0) for t in np.linspace(tmin,tmax,tsteps)]
# ts=np.linspace(tmin,tmax,tsteps)
# for t in range(len(ts)):
#     ts[t]=m.acos((2.0*ts[t])-1.0)
    # ts[t]=ts[t]*m.pi

#phi range
phmin=0.0
phmax=1.0
phsteps=8
phs=np.linspace(phmin,phmax,phsteps)
for p in range(len(phs)):
    # phs[p]=2*m.pi*phs[p]
    phs[p]=p*2*m.pi/phsteps

#grid points
grid=[(x,y,z) for x in rs for y in ts for z in phs]


################################################################################

#CREATES DIRECTORIES FOR GRID DATA

if not os.path.exists("./Results/Grid/"):
    os.makedirs("./Results/Grid/")
datafile=open("./Results/Grid/gridpos.txt", "wb+")
np.savetxt(datafile, grid, fmt='%.3e')
datafile.close()


################################################################################

#SAVES DISTANCES TO STARS FOR EACH GRID POINT

for i in range(len(grid)):
    os.makedirs("./Results/Grid/point_"+str(i).zfill(3))
    dist=[]
    for n in range(len(stars)):
        if distance(stars[n],grid[i])<strom[n]:
            dist.append(distance(stars[n],grid[i]))
        else:
            dist.append(0)
    datafile = open("./Results/Grid/point_"+str(i).zfill(3)+"/point.data", "wb+")
    np.savetxt(datafile, sweight(grid[i]), fmt='%.3e')
    np.savetxt(datafile, dist, fmt='%.3e')
    datafile.close()

################################################################################

#UPDATES EUCODE_GRID XSTIFF FILE FOR HALO INSTANCE

zrd=np.loadtxt("./Code/eucode_grid/zrd_000.txt")
replace = {'var0':str(len(stars)), 'var1':str(len(zrd))}
f=fileinput.FileInput("./Code/xstiff.f", inplace=True)
for line in f:
    for i, j in replace.items():
        line = line.replace(i,j)
    print(line, end='')
f.close()

f=fileinput.FileInput("./Code/xstiff.d.f", inplace=True)
for line in f:
    line=line.replace('var0',str(len(stars)))
    print(line, end='')
f.close()


################################################################################
