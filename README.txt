# LABS
Bioturbation model and its extension eLABS v0.1

Bioturbation code from Choi et al. (2002) Conmuters &amp; Geosciences 28, 213-222
https://www.iamg.org/index.php/publisher/articleview/frmArticleID/115

Extended to include calculation of water flow and oxygen and organic matter concentration fields. 

You need BLAS & UMFPACK libraries installed. 

Then, to make a simulation, following steps need be followed. 

(1) Specify parameters in LABS.f90, Module_Globals.f90 and Parameters_IN.txt.
    (a) Specify directory where output is stored in line 40 of LABS.f90. 
    (b) Turn switches on/off depending on simulations in lines 91-134 in LABS.f90.
    (c) Change reaction rate, shear velocity and oxygen concentration (if you want) in lines 203, 196-197; 200-201, and 193, respectively, in Module_Globals.f90.  
    (d) Change sedimentation rate, porosity and simulation duration (if you want) in lines 15, 17 and 21 in Parameters_IN.txt. 
(2) Compile codes by typing 'make'. If successful, labs.exe is created. If you previously created labs.exe, type 'make clear' before compiling. 
(3) Simulate an experiment by typing './labs your_experiment_name'.
(4) Analyze results in a directory of your_experiment_name in the directory you specified in LABS.f90.
    (a) Flux can be immediately plotted (with/without experiments finished) by double clicking biot-flx2.plt (you need gnuplot for this). 
    (b) Oxygen profile is stored in data directory in the your_experiment_name directory as O2data-xxxxx.txt files where xxxx records time steps in the experiments.
        To visualize oxygen profile, use any software. For example, you can use ImageJ.
    (c) Burrow geometry is stored in the same directory as O2 profile, with file name of txtimg-xxxx.txt where xxxx again records time steps in the experiments.
    (d) Stream function is stored in the same directory as O2 profile and burrow geometry, with file name of test2d(F)-Qcy-xxxx.txt (see above for xxx). 
    (e) You can examine biodiffusion coefficient after the experiment is finished by double clicking biot-Db.plt (again you need gnuplot for this).

Without specifying any parameters, you can run an original LABS simulation (Fig. 2a). 
To modify food control on organisms to change the burrow density in LABS (Fig. 2b), you need switch on 'i-shape' (i.e.,i-shape =.true.) in LABS.f90. 
To simulate oxygen reactive transport, you need switch on 'oxygen_ON' (i.e.,oxygen_ON =.true.) in LABS.f90.
To calculate advective water flow, you need switch on 'flow_ON' (i.e.,flow_ON =.true.) in LABS.f90.

More specifically, to simulate individual runs in Section 3, you need following changes in the codes.

~~~Section 3.1~~~
(a)  oxygen_ON = .true., flow_ON = true. (in LABS.f90)       
(b)  oxygen_ON = .true., flow_ON = true., oxFB_ON = .false. (in LABS.f90)        
(c)  oxygen_ON = .true., flow_ON = true. (in LABS.f90), kdcy = 1e-1/220.e-7, bio_fact = 1e3 (in Module_Globals.f90)      
(d)  oxygen_ON = .true., flow_ON = true. (in LABS.f90), 0.6 Porosity of sediment (in Parameters_IN.txt)           
(e)  oxygen_ON = .true., flow_ON = true. (in LABS.f90), shearfact = 10d0 (in Module_Globals.f90)
(f)  oxygen_ON = .true. (in LABS.f90)  
~~~Section 3.2~~~ 
(a)  oxygen_ON = .true.(in LABS.f90), 0.6 Porosity of sediment, 10 End of output (in Parameters_IN.txt)         
(b)  oxygen_ON = .true., oxFB_ON = .false. (in LABS.f90), 0.6 Porosity of sediment, 10 End of output (in Parameters_IN.txt)        
(c)  oxygen_ON = .true.(in LABS.f90), kdcy = 1e-1/220.e-7, bio_fact = 1e3 (in Module_Globals.f90),0.6 Porosity of sediment, 10 End of output (in Parameters_IN.txt)      
(d)  oxygen_ON = .true.(in LABS.f90), 0.6 Porosity of sediment, 10 End of output, 0.15 Sedimentation rate (in Parameters_IN.txt) 
~~~Section 3.3~~~ 
(a)  oxygen_ON = .true.(in LABS.f90), iox = 220.e-6 (in Module_Globals.f90), 0.6 Porosity of sediment, 10 End of output (in Parameters_IN.txt)         
(b)  oxygen_ON = .true., oxFB_ON = .false. (in LABS.f90), iox = 220.e-6 (in Module_Globals.f90), 0.6 Porosity of sediment, 10 End of output (in Parameters_IN.txt)        
(c)  oxygen_ON = .true.(in LABS.f90), kdcy = 1e-1/220.e-7, bio_fact = 1e3,iox = 220.e-6 (in Module_Globals.f90),0.6 Porosity of sediment, 10 End of output (in Parameters_IN.txt)      
(d)  oxygen_ON = .true.(in LABS.f90), iox = 220.e-6 (in Module_Globals.f90), 0.6 Porosity of sediment, 10 End of output, 0.15 Sedimentation rate (in Parameters_IN.txt) 