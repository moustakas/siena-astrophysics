EXEC   = MakeDiskGalaxy

OBJS   = main.o cooling.o read_param.o effmodel.o  globvars.o \
	halo.o haloset.o escape.o \
	disk.o diskset.o \
	bulge.o bulgeset.o \
	gas.o  gasset.o \
	save.o init.o structure.o \
	structure_withcontraction.o \
	misc.o toomre.o force.o \
	forcetree.o \
	nrsrc/bessi0.o   nrsrc/bessj1.o   \
	nrsrc/erff.o     nrsrc/polint.o   nrsrc/spline.o \
	nrsrc/bessi1.o   nrsrc/bessk0.o   nrsrc/gasdev.o  \
	nrsrc/qromb.o    nrsrc/splint.o  \
	nrsrc/bessj0.o   nrsrc/bessk1.o   \
	nrsrc/nrutil.o   nrsrc/ran1.o     nrsrc/trapzd.o \
	nrsrc/gammp.o    nrsrc/gcf.o      nrsrc/gser.o  \
        nrsrc/gammln.o   nrsrc/indexx.o   nrsrc/zriddr.o


INCL   = prototypes.h globvars.h nrsrc/nrutil.h nrsrc/nr.h

OPTIMIZE =  -O2  -Wall

OPT    =    -DREDSHIFT_SCALING  -DV_SCALING  # -DADIABATIC_CONTRACTION

CFLAGS = $(OPTIMIZE) $(OPT)

LIBS   =  -lm 

CC     =  gcc


$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC)  $(LIBS)

$(OBJS): $(INCL) 

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)












