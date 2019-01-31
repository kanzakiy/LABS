# Start of the makefile
# Defining variables
objects = Module_Globals.o Module_Gif_1.o Module_Gif_2.o \
 pfit.o  mod_o2_dif+adv.o mod_NS_MAC_2D.o LABS.o  
wrapper = umf4_f77wrapper.o
# switch = -g -fcheck=all
switch = -O2 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising \
  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=gnu  -pedantic  -fbacktrace
# switch = 
libs = -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -llapack -lopenblas
f90comp = gfortran
# Makefile
labs: $(objects)
	$(f90comp) -o $@ $(switch) $(objects) $(wrapper) $(libs)
globalvariables.mod: Module_Globals.o Module_Globals.f90
	$(f90comp) -c $(switch) Module_Globals.f90
bin_io.mod: Module_Gif_1.o Module_Gif_1.f90
	$(f90comp) -c $(switch) Module_Gif_1.f90
gif_util.mod: bin_io.mod Module_Gif_2.o Module_Gif_2.f90
	$(f90comp) -c $(switch) Module_Gif_2.f90
mod_pfit.mod: pfit.o pfit.f90
	$(f90comp) -c $(switch) pfit.f90
o2_diffusion.mod: globalvariables.mod mod_o2_dif+adv.o mod_o2_dif+adv.f90
	$(f90comp) -c $(switch) mod_o2_dif+adv.f90
ns_mac_2d.mod: globalvariables.mod mod_NS_MAC_2D.o mod_NS_MAC_2D.f90
	$(f90comp) -c $(switch) mod_NS_MAC_2D.f90
%.o: %.f90
	$(f90comp) -c $(switch) $<
# Cleaning everything
clean:
	rm $(objects)
# End of the makefile