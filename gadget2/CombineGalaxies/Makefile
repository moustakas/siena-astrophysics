
#OPT    +=  -DECC_NONZERO
OPT    +=  -DORBITCORRECTION   # this will correct the orbit if the halos already overlap
                               # the dark halos are assumed to follow a Hernquist profile in this

EXEC   = CombineGalaxies

OBJS   = main.o  globvars.o \
	 load.o save.o \
	 move.o turn.o \
	 nrsrc/nrutil.o



INCL   = globvars.h nrsrc/nrutil.h 


SYSTYPE="Grawp"
#SYSTYPE="MPA"
#SYSTYPE="Mako"
#SYSTYPE="Regatta"
#SYSTYPE="RZG_LinuxCluster"
#SYSTYPE="RZG_LinuxCluster-gcc"


CC       =    cc        # sets the C-compiler (default)
OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)


ifeq ($(SYSTYPE),"Grawp")
GSL_INCL =  -I/usr/local/gsl/gcc/include
GSL_LIBS =  -L/usr/local/gsl/gcc/lib
endif

ifeq ($(SYSTYPE),"MAKO")
GSL_INCL = 
GSL_LIBS = 
endif


ifeq ($(SYSTYPE),"MPA")
GSL_INCL = -I/usr/common/pdsoft/include
GSL_LIBS = -L/usr/common/pdsoft/lib
endif


ifeq ($(SYSTYPE),"Regatta")
CC       =   cc
OPTIMIZE =   -g -q64 #  -O5 -qstrict -qipa 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/include  -I/afs/rzg/@sys/lib/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/lib                
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster")
CC       =   icc
OPTIMIZE =   -O3 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_linux/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_linux/lib                -static
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster-gcc")
CC       =   gcc
OPTIMIZE =   -O3 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_linux_gcc3.2/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_linux_gcc3.2/lib
endif



LIBS   =    -lm  $(GSL_LIBS) -lgsl -lgslcblas 

CFLAGS =    $(OPTIONS) $(OPT)  $(OPTIMIZE)  $(GSL_INCL)


$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)  -o $(EXEC)  

$(OBJS): $(INCL) 

clean:
	 rm -f $(OBJS) $(EXEC)

