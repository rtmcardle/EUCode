import numpy as np
import math as m
import fileinput
import os

################################################################################

#FUNCTIONS

#returns dt/dz at redshift z
def dtdz(z):
    return 3.15e7*0.98e10/h100/(1.0+z)**2/np.sqrt(1.0+w0*z)

#returns probability of star being born at redshift z
def birthrate(z): #extracted from Xu et al. 2013 fig. 2
    return (3.28e-10*z**6-3.97e-8*z**5+1.99e-6*z**4-5.31e-5*z**3-0.01*z+0.02) \
            *pctocmcube*vhalo*dtdz(z)


################################################################################

#CLASSES

class Star():
    def __init__(self):
        #ASSIGNS MASS
        self.mass = np.random.choice(ms, p=mprobs)

        #ASSIGNS MASS DEPENDENT PROPERTIES (SCHAERER 2002)
        if self.mass < 5:
            self.lifetime = yrs * 6.190e7
            self.temp = 10.0 ** 4.440
            self.Q = 1.097e45
            self.lum = 10 ** 2.870
        elif self.mass < 9:
            self.lifetime = yrs * 2.022e7
            self.temp = 10.0 ** 4.622
            self.Q = 1.794e47
            self.lum = 10 ** 3.709
        elif self.mass < 15:
            self.lifetime = yrs * 1.040e7
            self.temp = 10.0 ** 4.759
            self.Q = 1.398e48
            self.lum = 10 ** 4.324
        elif self.mass < 25:
            self.lifetime = yrs * 6.459e6
            self.temp = 10.0 ** 4.850
            self.Q = 5.446e48
            self.lum = 10 ** 4.890
        elif self.mass < 40:
            self.lifetime = yrs * 3.864e6
            self.temp = 10.0 ** 4.900
            self.Q = 1.873e49
            self.lum = 10 ** 5.420
        elif self.mass < 60:
            self.lifetime = yrs * 3.464e6
            self.temp = 10.0 ** 4.943
            self.Q = 3.481e49
            self.lum = 10 ** 5.715
        elif self.mass < 80:
            self.lifetime = yrs * 3.012e6
            self.temp = 10.0 ** 4.970
            self.Q = 5.938e49
            self.lum = 10 ** 5.947
        elif self.mass < 120:
            self.lifetime = yrs * 2.521e6
            self.temp = 10.0 ** 4.981
            self.Q = 1.069e50
            self.lum = 10 ** 6.243
        elif self.mass < 200:
            self.lifetime = yrs * 2.204e6
            self.temp = 10.0 ** 4.999
            self.Q = 2.292e50
            self.lum = 10 ** 6.574
        elif self.mass < 300:
            self.lifetime = yrs * 2.047e6
            self.temp = 10.0 ** 5.007
            self.Q = 4.029e50
            self.lum = 10 ** 6.819
        elif self.mass < 400:
            self.lifetime = yrs * 1.974e6
            self.temp = 10.0 ** 5.028
            self.Q = 5.573e50
            self.lum = 10 ** 6.984
        else:
            self.lifetime = yrs * 1.899e6
            self.temp = 10.0 ** 5.029
            self.Q = 7.380e50
            self.lum = 10 ** 7.106


        #CALCULATED FROM L = R**2 * T**4
        self.rad = solrad*m.sqrt(self.lum)*(soltemp/self.temp)**2.0


        #ASSIGNES GALACTIC POSITION
        self.zbirth=np.random.choice(zbs, p=zbprobs)
        self.r=np.random.choice(rs, p=rprobs)
        self.theta=np.random.choice(ts, p=tprobs)
        self.phi=np.random.choice(phs, p=phprobs)


################################################################################

#PHYSICAL CONSTANTS

solarmass = 1.9891e30 #in kgs
hmass = 1.6737e-27 #in kgs
yrs=3.154e7 #in seconds
sollum=3.828e26 #in watts
solrad=6.957e10 #in cm
soltemp=5778.0 #in K
boltz=5.67e-8 #Stefan-Boltzman const
pctocmcube=2.94e-73 #conversion factor


################################################################################

#MODEL

#halo model constants
halo = 10.0**5.5 #halo mass in solar masses
sfe = 0.6e-3 #Hartwig et al. 2015
imf=1.35 #value for Salpeter IMF
nden = float(np.loadtxt("./Code/eucode_background/denhalo.txt"))

#calculates halo model properties
mhalo = halo*solarmass
mden = nden*hmass
vhalo=mhalo/mden
rhalo=(3.0*vhalo/(4.0*m.pi))**(1.0/3.0)
stellar_mass=halo*sfe

#imports universe model parameters (model.data)
replace = {'d':'e'}
f=fileinput.FileInput("./Code/model.grid.data", inplace=True)
for line in f:
    for i, j in replace.items():
        line = line.replace(i,j)
    print(line, end='')
model=np.genfromtxt("./Code/model.grid.data", dtype=float, delimiter=" ")
h100=model[0][1]
w0=model[0][2]


################################################################################

#STELLAR PROPERTY RANGES

#mass distribution
mmin=25.0
mmax=500.0
mstep=0.1
msteps=int((mmax-mmin)/mstep)

#birth distribution
zbmin=15.0
zbmax=30.0
zstep=model[2][1]
zbsteps=int((zbmax-zbmin)/zstep)

#position distribution
#radius
rsteps=1000
rmax=rhalo
rmin=rmax/rsteps
rmin1=rmax/4.0 #radius of core w/ uniform number density

#theta
tmin=0.0
tmax=1.0
tsteps=1000

#phi
phmin=0.0
phmax=1.0
phsteps=1000


################################################################################

#STELLAR PROPERTY DISTRIBUTIONS

#masses
ms=np.linspace(mmin,mmax,msteps)
ma = [x**-imf for x in ms]
tot = sum(ma)
mprobs = [x/tot for x in ma]


#births
zbs=np.linspace(zbmin,zbmax,zbsteps)
ztot=sum(birthrate(zbs))
zbprobs = [birthrate(x)/ztot for x in zbs]


#position
#radius
rs=np.linspace(rmin,rmax,rsteps)
r1=[]
for n in range(rsteps):
    if rs[n] <= rmin1:
        r1.append(rmin1**-2.0)
    else:
        r1.append(rs[n]**-2.0)
tot = sum(r1)
rprobs=[r1[n]/tot for n in range(rsteps)]

#theta
ts=np.linspace(tmin,tmax,tsteps)
for t in range(len(ts)):
    ts[t]=m.acos((2.0*ts[t])-1.0)
tprobs = [1.0/tsteps for i in range(tsteps)]

#phi
phs=np.linspace(phmin,phmax,phsteps)
for p in range(len(phs)):
    phs[p]=2*m.pi*phs[p]
phprobs = [1.0/phsteps for i in range(phsteps)]


###############################################################################

#CREATES STARS

star_list = []

#adds stars to list until stellar mass is reached
while sum(s.mass for s in star_list) < stellar_mass:
    star_list.append(Star())

#sorts stars by mass
star_list.sort(key=lambda s: s.mass)


################################################################################

#WRITES DATA TO DATA FILES

if not os.path.exists("./Results/Stars/"):
    os.makedirs("./Results/Stars/")
datafile = open("./Results/Stars/stardata.txt", 'w+')
data = []
for s in star_list:
    data.append(np.array([s.mass,s.lifetime,s.temp,s.rad,s.Q,s.zbirth,s.r,s.theta,s.phi]).T)
np.savetxt(datafile, data, fmt='%.3e')
datafile.close()


################################################################################
