FC=gfortran
FOPTS=-O
lirec: xstiff.o input.o rates.o indx.o fcn1.o fcnj.o ke.o \
	jacder.o tspln.o estate.o derivsp.o fcn1p.o kspln.o \
	kspln1.o kspln2.o hspln.o fcnjt.o \
	recfast.o stifbs.o simpr.o pzextr.o ludcmp.o \
	luksub.o odeint.o spline.o splint.o
	$(FC) -o lirec $(FOPTS) xstiff.o input.o \
	rates.o indx.o fcn1.o fcnj.o ke.o jacder.o \
	tspln.o estate.o derivsp.o fcn1p.o fcnjt.o \
	kspln.o kspln1.o kspln2.o hspln.o recfast.o \
		stifbs.o \
		simpr.o \
		pzextr.o \
		ludcmp.o \
		luksub.o \
		odeint.o \
		spline.o \
		splint.o

xstiff.o:
	$(FC) $(FOPTS) -c xstiff.f
input.o:
	$(FC) $(FOPTS) -c input.f
rates.o:
	$(FC) $(FOPTS) -c rates.f
indx.o:
	$(FC) $(FOPTS) -c indx.f
fcn1.o:
	$(FC) $(FOPTS) -c fcn1.f
fcnj.o:
	$(FC) $(FOPTS) -c fcnj.f
ke.o:
	$(FC) $(FOPTS) -c ke.f
jacder.o:
	$(FC) $(FOPTS) -c jacder.f
tspln.o:
	$(FC) $(FOPTS) -c tspln.f
estate.o:
	$(FC) $(FOPTS) -c estate.f
derivsp.o:
	$(FC) $(FOPTS) -c derivsp.f
fcn1p.o:
	$(FC) $(FOPTS) -c fcn1p.f
kspln.o:
	$(FC) $(FOPTS) -c kspln.f
kspln1.o:
	$(FC) $(FOPTS) -c kspln1.f
kspln2.o:
	$(FC) $(FOPTS) -c kspln2.f
fcnjt.o:
	$(FC) $(FOPTS) -c fcnjt.f
hspln.o:
	$(FC) $(FOPTS) -c hspln.f
recfast.o:
	$(FC) $(FOPTS) -c recfast.f
stifbs.o:
	$(FC) $(FOPTS) -c stifbs.f
simpr.o:
	$(FC) $(FOPTS) -c simpr.f
pzextr.o:
	$(FC) $(FOPTS) -c pzextr.f
ludcmp.o:
	$(FC) $(FOPTS) -c ludcmp.f
luksub.o:
	$(FC) $(FOPTS) -c luksub.f
odeint.o:
	$(FC) $(FOPTS) -c odeint.f
spline.o:
	$(FC) $(FOPTS) -c spline.f
splint.o:
	$(FC) $(FOPTS) -c splint.f

clean:
	/bin/rm  xstiff.o input.o core xstiff
