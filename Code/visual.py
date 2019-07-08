import vpython as vi
import numpy as np
import os

#Visualizes the abundance data calculated by eugrid.sh

################################################################################

#CLASSES

class Star(vi.sphere):
    def __init__(self,*args,**kwargs):
        global fronts
        super(Star,self).__init__(*args,**kwargs)
        self.color = vi.vec(255,255,0)
        self.opacity = 0.0
        self.born = False
        self.dead = False
        self.strom = Front(pos=self.pos, radius=self.radius)
        fronts.append(self.strom)

    def birth(self):
        self.born = True
        self.opacity = 1.0

    def death(self):
        self.dead = True
        self.opacity = 0.2

class Front(vi.sphere):
    def __init__(self,*args,**kwargs):
        super(Front,self).__init__(*args,**kwargs)
        self.color = vi.vec(255,255,255)
        self.opacity = 0.0
        self.born = True
        self.dead = False
        self.zrd = []

    def death(self):
        self.dead = True
        self.opacity = 0

class Grid(vi.sphere):
    def __init__(self,*args,**kwargs):
        super(Grid,self).__init__(*args,**kwargs)
        self.color = vi.vec(255,255,255)
        self.opacity = 0.1
        self.colors = []
        self.opacities = []


################################################################################

#SETS VPYTHON SCREEN

scene=vi.canvas(title="EUpdate Visualization", x=0, y=0, height=720, width=1280,
                center=vi.vec(0,0,0))


################################################################################

#LOADS HALO INSTANCE DATA

stardata=np.loadtxt("./Results/Stars/stardata.txt")
gridpos=np.loadtxt("./Results/Grid/gridpos.txt")
position=stardata[:,6:]
radius=stardata[:,3]


################################################################################

#CREATES STARS, FRONTS, AND GRID POINTS

stars = []
fronts = []
for n,m in zip(position,radius):
    x=n[0]*np.sin(n[1])*np.cos(n[2])
    y=n[0]*np.sin(n[1])*np.sin(n[2])
    z=n[0]*np.cos(n[1])
    #r=m[3]
    r=m*10**8
    stars.append(Star(pos=vi.vec(x,y,z), radius = r))

grid = []
for n in range(len(gridpos)):
    x=gridpos[n][0]*np.sin(gridpos[n][1])*np.cos(gridpos[n][2])
    y=gridpos[n][0]*np.sin(gridpos[n][1])*np.sin(gridpos[n][2])
    z=gridpos[n][0]*np.cos(gridpos[n][1])
    r=10**19.5
    #r=radius[n][3]
    grid.append(Grid(pos=vi.vec(x,y,z), radius = r))


################################################################################

#PREPARES RADIUS AND OPACITY DATA

for n in range(len(stars)):
    data = np.loadtxt("./Code/eucode_grid/zrd_"+str(n).zfill(3)+".txt")
    fronts[n].zrd = data[:,1]

for g in range(len(grid)):
    data = np.loadtxt("./Results/Grid/point_"+str(g).zfill(3)+"/abund.avg")
    for n in range(len(data)):
        grid[g].colors.append(vi.vec(data[n][0]*1e0,data[n][1]*1e1,data[n][2]*1e1))
        grid[g].opacities.append(1-(data[n][0]*.9))


################################################################################

#VISUAL LOOP

for n in range(len(fronts[0].zrd)):
    vi.rate(10)
    for m in range(len(stars)):
        if fronts[m].zrd[n] != 0:
            stars[m].opacity = 1.0
            stars[m].born = True
            fronts[m].opacity = 0.25
        elif stars[m].born == True:
            stars[m].death()
        fronts[m].radius = fronts[m].zrd[n]
    if n%100 == 0:
        print(n)
exit()


################################################################################
