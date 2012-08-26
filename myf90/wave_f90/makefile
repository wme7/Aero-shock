.KEEP_STATE:
all: main
FOR=f90 
#mdir=-Imodule
#wmdir=-mdir module

msource= mod_wave.f90
mobj  = $(msource:%.f90=%.o)

obj=$(mobj) main.o bdryflux.o   errors.o     impvar.o\
legendre.o   ordr.o calc.o       euler.o      init.o \
rk.o cdt.o   phifun.o func.o       input.o exac.o\
superconv.o coefb.o      gauleg.o     intflux.o claim.o post.o\
 pordr.o eno.o  interpolation.o

main: $(obj)
	$(FOR) -o run $(mdir) $(obj)
	rm -f *.o; rm -f *.dat
	run

%.o:    %.f90
	$(FOR)  $(mdir) $(wmdir) -c $<
  
clear : 
	rm -f *.o ; rm -f *~

clean :
	rm -f *~

