OSNF_DIR = osnf

.PHONY: osnf_code cleanall
CLEANDIRS = $(OSNF_DIR) ./

# LINUX or WIN64
#PLATFORM = WIN64   
PLATFORM = LINUX
DEBUG =#-gstabs -fbounds-check #-pg
MPI    =#-DMPI1
OPT    =-O3  #-fopenmp-simd -fopenmp 

# these three lines should be edited for your system. On systems 
# that do not have separate fortran and c libraries, set NETCDF_FOR and NETCDF_C
# to the same, and set NETCDF_LIB to -lnetcdf (i.e. without the extra f)
#NETCDF_FOR=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.4-mac/
#NETCDF_C=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.1.1-mac/
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
VAR_TYPE = 1 # 0 single, 1 double


main.exe	:  odelib.a  solar_system.$(OBJ) 
	$(FOR2) $(FFLAGS2)main.exe solar_system.$(OBJ)   \
	 -lm odelib.a ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG) -L$(OSNF_DIR) 

odelib.a	:   osnf_code
	$(AR) rc odelib.a $(OSNF_DIR)/numerics.$(OBJ) $(OSNF_DIR)/zeroin.$(OBJ) $(OSNF_DIR)/sfmin.$(OBJ) \
				$(OSNF_DIR)/fmin.$(OBJ) $(OSNF_DIR)/r1mach.$(OBJ) \
                $(OSNF_DIR)/d1mach.$(OBJ) $(OSNF_DIR)/dfsid1.$(OBJ) \
                $(OSNF_DIR)/poly_int.$(OBJ) $(OSNF_DIR)/find_pos.$(OBJ) \
                $(OSNF_DIR)/svode.$(OBJ) \
                $(OSNF_DIR)/slinpk.$(OBJ) $(OSNF_DIR)/vode.$(OBJ) \
                $(OSNF_DIR)/dlinpk.$(OBJ) $(OSNF_DIR)/vode_integrate.$(OBJ) \
                $(OSNF_DIR)/erfinv.$(OBJ) $(OSNF_DIR)/tridiagonal.$(OBJ) \
                $(OSNF_DIR)/hygfx.$(OBJ) $(OSNF_DIR)/random.$(OBJ)			

solar_system.$(OBJ)   : solar_system.f90 osnf_code
	$(FOR) -cpp solar_system.f90 -I$(OSNF_DIR) -I${NETCDFMOD}  -cpp $(FFLAGS)solar_system.$(OBJ)
osnf_code:
	$(MAKE) -C $(OSNF_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	odelib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	
	

