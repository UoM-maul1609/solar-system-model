# LINUX or WIN64
#PLATFORM = WIN64   
PLATFORM = LINUX
DEBUG =#-gstabs -fbounds-check #-pg
MPI    =#-DMPI1
OPT    =-O3  -fopenmp-simd -fopenmp 

# these three lines should be edited for your system. On systems 
# that do not have separate fortran and c libraries, set NETCDF_FOR and NETCDF_C
# to the same, and set NETCDF_LIB to -lnetcdf (i.e. without the extra f)
NETCDF_FOR=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.4-mac/
NETCDF_C=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.1.1-mac/
NETCDF_LIB=-lnetcdff 


ifeq ($(strip $(PLATFORM)),$(strip LINUX))
NETCDFLIB=-L ${NETCDF_FOR}/lib/  \
          -L ${NETCDF_C}/lib/  #/usr/lib64/
NETCDFMOD= ${NETCDF_FOR}/include/ #/usr/lib64/gfortran/modules/

FOR = gfortran -c  #-fno-underscoring 
FOR2 = gfortran  #-fno-underscoring 
AR = ar 
RANLIB = ranlib 
SIXTY_FOUR_F= 
else
NETCDFLIB=../netcdf-4.1.3-mingw/lib/
NETCDFMOD=../netcdf-4.1.3-mingw/include/
FOR = x86_64-w64-mingw32-gfortran -c  #-fno-underscoring 
FOR2 = x86_64-w64-mingw32-gfortran  #-fno-underscoring 
AR = x86_64-w64-mingw32-ar
RANLIB = x86_64-w64-mingw32-ranlib
SIXTY_FOUR_F= -static -static-libgfortran #-m32 # note, uncomment the m32 for 32-bit model
endif



OBJ = o

FFLAGS = $(SIXTY_FOUR_F) $(OPT)  $(DEBUG) -Dheliocentric -o # helicentric means that we keep the sun at (0,0,0)
FFLAGS2 = $(SIXTY_FOUR_F) $(OPT) $(DEBUG) -o # 


main.exe	:  odelib.a  solar_system.$(OBJ) 
	$(FOR2) $(FFLAGS2)main.exe solar_system.$(OBJ)   \
	 -lm odelib.a ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG) 

odelib.a	:   d1mach.$(OBJ) vode.$(OBJ) dlinpk.$(OBJ) acdc.$(OBJ) colmod.$(OBJ) opkdmain.$(OBJ) opkda1.$(OBJ) opkda2.$(OBJ)
	$(AR) rc odelib.a d1mach.$(OBJ) vode.$(OBJ) dlinpk.$(OBJ) acdc.$(OBJ) colmod.$(OBJ) opkdmain.$(OBJ) opkda1.$(OBJ) opkda2.$(OBJ); $(RANLIB) odelib.a
solar_system.$(OBJ)   : solar_system.f90 
	$(FOR) -cpp solar_system.f90 -I ${NETCDFMOD}  -cpp $(FFLAGS)solar_system.$(OBJ)
d1mach.$(OBJ) 	: ./other/d1mach.f 
	$(FOR) 	./other/d1mach.f $(FFLAGS)d1mach.$(OBJ)  
opkdmain.$(OBJ) : ./other/opkdmain.f 
	$(FOR) 	./other/opkdmain.f $(FFLAGS)opkdmain.$(OBJ)
opkda1.$(OBJ) : ./other/opkda1.f 
	$(FOR) 	./other/opkda1.f -w $(FFLAGS)opkda1.$(OBJ)
opkda2.$(OBJ) : ./other/opkda2.f 
	$(FOR) ./other/opkda2.f -w $(FFLAGS)opkda2.$(OBJ)
vode.$(OBJ) : ./other/vode.f 
	$(FOR) 	./other/vode.f $(FFLAGS)vode.$(OBJ)
acdc.$(OBJ) : ./other/acdc.f 
	$(FOR) 	./other/acdc.f $(FFLAGS)acdc.$(OBJ)
dlinpk.$(OBJ) : ./other/dlinpk.f 
	$(FOR) 	./other/dlinpk.f -w $(FFLAGS)dlinpk.$(OBJ)
colmod.$(OBJ) : ./other/colmod.f 
	$(FOR) 	./other/colmod.f -Wno-all $(FFLAGS)colmod.$(OBJ)
clean :
	rm *.exe  *.o *.mod *~ \
	odelib.a

