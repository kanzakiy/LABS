# Start of the makefile
# Defining variables
objects = Module_Globals.o LABS.o Module_Gif_1.o Module_Gif_2.o \
 pfit.o  mod_o2_dif+adv.o mod_NS_MAC_2D.o 
wrapper = umf4_f77wrapper.o
switch = -g -fcheck=all
libs = -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -llapack -lopenblas
f90comp = gfortran
# Makefile
labs: $(objects)
	$(f90comp) -o labs $(switch) $(objects) $(wrapper) $(libs)
globalvariables.mod: Module_Globals_yk6.o Module_Globals_yk6.f90
	$(f90comp) -c $(switch) Module_Globals_yk6.f90
bin_io.mod: Module_Gif_1.o Module_Gif_1.f90
	$(f90comp) -c $(switch) Module_Gif_1.f90
gif_util.mod: bin_io.mod Module_Gif_2.o Module_Gif_2.f90
	$(f90comp) -c $(switch) Module_Gif_2.f90
mod_pfit.mod: pfit.o pfit.f90
	$(f90comp) -c $(switch) pfit.f90
o2_diffusion.mod: globalvariables.mod mod_o2_dif+adv_v7.o mod_o2_dif+adv_v7.f90
	$(f90comp) -c $(switch) mod_o2_dif+adv_v7.f90
ns_mac_2d.mod: globalvariables.mod mod_NS_MAC_2D_v3.o mod_NS_MAC_2D_v3.f90
	$(f90comp) -c $(switch) mod_NS_MAC_2D_v3.f90
%.o: %.f90
	$(f90comp) -c $(switch) $<
# Cleaning everything
clean:
	rm $(objects)
# End of the makefile