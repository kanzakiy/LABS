
   ! ***********************************************************************
   !         Program name  LABS  (Lattice-Automaton Bioturbation Simulator)
   !
   !         Project PI:   Bernard P. Boudreau
   !         Programmers:    Jae Choi and Frederique Fracois-Carcaillet
   !         Locale        Halifax, Nova Scotia, Canada
   !         Institution   Department of Oceanography, Dalhousie University
   !
   !         Language      Fortran 90
   !         Libraries      IMSL numerical libs  used in subroutine data_anaysis
   !
   !         Date modified 19 July 2000
   ! ***********************************************************************


   PROGRAM LABS

   Use GlobalVariables
   use O2_diffusion
   use ieee_arithmetic
   use NS_MAC_2D

   implicit none
   Logical               :: TimeToOutput
   Integer               :: i, DepthToScan, TimeNew, ParticlesToSediment
   Character*14          :: CurrentTimeString
   Real                  :: Area_Period
   Type(CoOrdinates)     :: Point
   integer :: xx, yy
   logical               :: IsWater, IsOrganism, IsParticle
   type(coordinates)     :: org_mid
   logical :: rec_flow
   
   write(Today,*) "_testing_flow_shear5_v5_"   ! this string goes into names of output files
   

       OPEN(unit = File_Profile,  file = 'C:/Users/YK/Desktop/biot-res/DepthProfiles'          &
                               //trim(adjustl(today))//'.OUT', status = 'unknown')
       OPEN(unit = File_Profile2,  file = 'C:/Users/YK/Desktop/biot-res/DepthProfiles2'        &
                               //trim(adjustl(today))//'.OUT', status = 'unknown')
       OPEN(unit = File_Profile3,  file = 'C:/Users/YK/Desktop/biot-res/DepthProfiles3'        &
                              //trim(adjustl(today))//'.OUT', status = 'unknown')
       OPEN(unit = File_Profile_Re,  file = 'C:/Users/YK/Desktop/biot-res/DepthProfiles_Re'        &
                             //trim(adjustl(today))//'.OUT', status = 'unknown')
       OPEN(unit = File_Scaling, file = 'C:/Users/YK/Desktop/biot-res/PorosityScale'        &
                              //trim(adjustl(today))//'.OUT', status = 'unknown')
       OPEN(unit = File_Activity, file = 'C:/Users/YK/Desktop/biot-res/ActivitySlope'        &
                               //trim(adjustl(today))//'.OUT', status = 'unknown')
       OPEN(unit = File_Displace, file = 'C:/Users/YK/Desktop/biot-res/MeanDisplacement'        &
                               //trim(adjustl(today))//'.OUT', status = 'unknown')
       OPEN(unit = File_Displace2, file = 'C:/Users/YK/Desktop/biot-res/MeanDisplacement2'        &
                               //trim(adjustl(today))//'.OUT', status = 'unknown')
       OPEN(unit = Poly_fit, file = 'C:/Users/YK/Desktop/biot-res/polyFit'                           &
                                  //trim(adjustl(today))//'.OUT', status = 'unknown')
       OPEN(unit = File_flux, file = 'C:/Users/YK/Desktop/biot-res/flux'                           &
                                  //trim(adjustl(today))//'.OUT', status = 'unknown')
       OPEN(unit = File_Log, file = 'C:/Users/YK/Desktop/biot-res/log'                           &
                                  //trim(adjustl(today))//'.txt', status = 'replace')
       OPEN(unit = File_flux_txt, file = 'C:/Users/YK/Desktop/biot-res/flux'                           &
                                  //trim(adjustl(today))//'.txt', status = 'unknown')

   ! Begin main program loop
   
   oxygen_ON = .true.
   ! oxygen_ON = .false.
   
   oxFB_ON = .true.
   ! oxFB_ON = .false.
   
   Resp_ON = .true.
   ! resp_ON = .false.
   
   errCHk = .true.
   errChk = .false.
   
   I_shape = .true.
   I_shape = .false.
   
   flow_ON = .true.
   ! flow_ON = .false.

   only_sed = .True.
   ! only_sed = .FALSE.
   
   Org_ID_ishape = 1



       CALL RANDOM_SEED    ! initialise random seed generator
       Call GetUserInput() ! user inputs
   
     call O2i_setup()
    call O2pre_setup()
    Call Output_O2txtImg()
    Call Output_txtImg()
   
   Savetime = 1 

       DO Time = 1, TimeMax

             Pr_Activity_year = 0.5 * (1. + cos(pi + 2.*pi*(Real(Time)/Real(Year)) ) )
             Pr_Activity_day = 0.5 * (1. + cos(pi + 2.*pi*(Real(Time)/Real(Day)) ) )
             CurrentTime = CurrentTimeString()
			 
         Do i = 1, N_Ind           ! for each organism

           Call Timetable(i)       ! obtain information about physiological/metabolic constraints

             IF (Org(i)%Can%Move) then

               Select Case (Org(i)%FeedingType)
                 Case ("D","d")                    ! Deposit feeders
                   CALL Rules_Of_Motion(i)         ! Locate the target direction for each organism
                     If (Org(i)%Can%Move) then
                       Call Mouth_locate(i, Point, Org(i)%Orientation) ! defines the block that represents the mouth region
                       Call Particles_ingest(i, Point, Org(i)%Orientation)
                       Call Particles_push(i, Point, Org(i)%Orientation)
                         If (Org(i)%Can%Move) then
                           Call Organism_move(i)
                           Call Particles_egest(i)
                         End if
                     End if
               End Select

             End if    ! organism can move

             if (resp_ON) call OrgResp(i)
			 
			 call org_midloc(i, 1, org(i)%headsize, org_mid)
			 flow_loc(1,i)%x = org_mid%x
			 flow_loc(1,i)%y = org_mid%y
			 call org_midloc(i,org(i)%bodysize-org(i)%headsize+1,org(i)%bodysize, org_mid)
			 flow_loc(2,i)%x = org_mid%x
			 flow_loc(2,i)%y = org_mid%y
         
         End Do        ! for each individual
         
         !  check for error 
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after org_move"; write(file_log,*)time, "err after org_move"; end if 
         end if 

       ! add new particles and sediment them down
           If (Mod(Time,Time_Sediment) .eq. 0) then
             If (Sedimentation_On) then
               ! determine number of particles to be input (the fraction of the area under the sinusoid) and do it
                   Area_Period = Real(Time_Sediment) * 0.5 *       &
                                 (1.+ cos(pi+ 2.*pi*(Real(Time - Time_Sediment*0.5)/Real(Year))))
                   ParticlesToSediment = ParticlesToSediment_y * ( Area_Period / Area_Total )
                   Call Particle_sediment(ParticlesToSediment)
                   Call LocateOrganisms_scan()  ! subroutine that collects info of the location of organisms ..
             End if

             ! Scan the water column and sediment down any free particles
                 DepthToScan = N_RowWater
                 Call WaterColumn_scan(DepthToScan)

             ! sneek an output of the current time to show how far the sim has progressed
                 write(*,*) CurrentTime, Time
                 
                 
         
         !  check for error 
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after sed"; write(file_log,*)time, "err after sed"; end if 
         end if 
                 

           End if

       ! check total number of particle is within tolerance .. if too large then shift observation window up one row
           If ((Total_N_Particles - Total_N_Particles0) .gt. ParticleTolerance) then
             If (Sedimentation_On) then
               Call Matrix_constrain()
             End if
             
             
         
         !  check for error 
         
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after constrain"; write(file_log,*)time, "err after constrain"; end if 
         end if 
             
           End if
		   
		   rec_flow = .false.
		   if (flow_ON) then 
		     if (any(abs(Vb)/=0).or.any(abs(Ub)/=0)) then 
			   call flow_calc()
			   rec_flow = .true.
			   ! if (time> 18700) call output_flow()
			   ! call output_flow()
			 end if 
			 
			 if (any(ieee_is_nan(Vo))) then
			   print *, 'NAN in Vo'
			   stop 1
			 end if 
			 if (any(ieee_is_nan(Uo))) then 
			   print *, "NAN in Uo"
			   stop 1
			 end if			 
		   end if 
   
   !!   working only when oxygen is switched ------------
         if (oxygen_ON) then
   
         ! call O2pre_setup()
         
         call OrgDecay()
         
         
         if (any(ieee_is_nan(O2%oxygen_use))) then
           print *, "NAN in oxygen_use"
           do yy = 1, n_row
             do xx = 1, n_col
               if (ieee_is_nan(O2(yy,xx)%oxygen_use)) then
                 print *, yy, xx
                 if (matrix(yy,xx)%class == p) then 
                   O2(yy,xx)%oxygen_use = 0.0
                 end if 
               end if 
             end do 
           end do 
		   stop 1
         end if 
         
		 if (.not.flow_ON) then
		   Vo = 0.0
		   Uo = 0.0
		 end if 
		 
         call oxygen_profile()
         
         ! Call Output_O2txtImg()  !  output every time step
         ! Call Output_txtImg()   !  output every time step
         
         call fluxes()
         
         If (Any(Time .eq. Time_Output)) then
           continue
         else 
           O2%oxygen_use = 0.0
         end if
         
         
         if (any(ieee_is_nan(O2%oxygen))) print *, "NAN in oxygen"
         
         do yy = 1, n_row
           do xx = 1, n_col
             if ((matrix(yy,xx)%class == p).and.(o2(yy,xx)%oxygen/=0.0)) then
               print *, "particle has oxygen --- permanent effect ---", yy, xx
             end if 
           end do
         end do
         
         end if
		 
		 
		 !!   working only when oxygen is switched ------------ END
         
         if (errDetect .and. errChk) then; call Output_txtImg(); call Output_txtImg_chk(); endif
		 
       ! check if time to output data
           If (Any(Time .eq. Time_Output)) then
             TimeToOutput = .true.
             ! write(*,*) time      !! Added to create movies -------YK 5/22/2017
             ! call OutputData()    !! Added to create movies -------YK 5/22/2017
             TimeNew = Time
           End If

           If (TimeTOOutput) then
             If (Movies_on) then                       ! movies ... output graphics along a linear timescale
               IF (Time .ge. TimeMin) THEN
                 TimeMin = TimeMin + TimeStep
                 Call OutputData()
                 if (oxygen_on) Call Output_O2txtImg()
			     if (rec_flow) call output_flow()
                 Call Output_txtImg()
                 call Output_txtImg_chk()
                 Savetime = savetime + 1       
               End If
               If ((Time - TimeNEW) .GT. 10*Day) then  ! continue output for 10 days
                 TimeToOutput = .false.                ! turn off output until the next output time
               End if
             Else              
                 call Make_matrix_chk()
                 call Matrix_err_chk()
                 if (errDetect) then; print *,time, "err when outputting"; write(file_log,*)time, "err when outputting"; end if                         ! snapshots at user-defined times
               
               Call OutputData()
               If (oxygen_ON) Call Output_O2txtImg()
			   if (rec_flow) call output_flow()
               Call Output_txtImg()
               call Output_txtImg_chk()
               O2%oxygen_use = 0.0
               Savetime = savetime + 1       
               TimeToOutput = .false.                  ! turn off output until the next output time
             End if
           End if
		   
		 
		 
		 Vo = 0.
		 Uo = 0.
		 flow_loc = coordinates(0,0)
		 Vb = 0.
		 Ub = 0.  
		   

       END Do    !  main time loop

   ! end main program

   ! clean up data allocations

       Deallocate (Org, Matrix, Particle, Org_Loc, Guts, MovementHistory, IngestionHistory, RespHistory, Dir_rec)
       
       deallocate (O2)
	   
	   deallocate (flow_loc, Ub, Vb)
	   deallocate (Ug,Vg, Pg, Dg) 
	   deallocate (Uo,Vo) 
	   deallocate (edif)

       Close(File_Displace)
       Close(File_Displace2)
       Close(File_Activity)
       Close(File_Scaling)
       Close(File_Profile)
       Close(File_Profile2)
       Close(File_Profile3)
       Close(File_Profile_Re)
       Close(File_oneD)
       Close(Poly_fit)
       Close(File_flux)
       Close(File_log)
       Close(File_flux_txt)

   End Program

   ! ********************************************************

   Subroutine GetUserInput()
     ! load data from files and allocate necessary memory
     ! called once by the main program, Automatons.f90

   Use GlobalVariables
   implicit none
   Character*1 :: FeedingType
   Integer     :: i, j, n, MaxOrgSize, MaxGutCapacity, MinHeadWidth, Lability_tmp
   Real        :: MatrixDepth, MatrixWidth, Depth_DeadZone, HalfLife
   Real        :: r_TimeMin, r_TimeMax, r_TimeStep
   Real        :: ParticleDensity, PArticleSize, Sum_Lability, MaxOrgSpeed, MinWidth, MaxWidth
   Real        :: tmp_MaxOrgSize, tmp_OrgSize, tmp_Width, tmp_Length
   Real        :: r_Width, r_Length, r_Speed
   Real        :: r_IngestRate, r_IngestSelectivity, r_Density
   
   ! input user defined variables ..
       OPEN(unit = File_Parameters, file = 'Parameters_IN.txt', status = 'old')
         READ(File_Parameters,*) Movies_on          ! OUTPUT IS (0) a movie or (1) user defined in the main program
         READ(File_Parameters,*) DepthAvoidance_on  ! Turn on/off DEPTH AVOIDANCE FUNCTION?
         READ(File_Parameters,*) Lability_ON        ! Turn on/off lability searching, selection, degradation functions
         READ(File_Parameters,*) Sedimentation_ON   ! Turn on or off the sedimentation of particles and the windowing function
         READ(File_Parameters,*) LocalMixing        ! completely turn off ingestion/egstion and particle loss through bottom boundary
         READ(File_Parameters,*) SaveData           ! save data to disk?
         READ(File_Parameters,*) PlotMatrix         ! plot a colour representation to screen ?
         READ(File_Parameters,*) PlotPorosity       ! plot to screen a locally averaged porosities
         READ(File_Parameters,*) PlotActivity       !   "  of the activity ...
         READ(File_Parameters,*) OutPutScaling      ! how porosity scales with the scale of resolution ...
         READ(File_Parameters,*) N_Ind              ! Enter the number of organisms
         READ(File_Parameters,*) MatrixDepth        ! Depth of the Matrix (Sediments and Water, inclusive, cm)
         READ(File_Parameters,*) MatrixWidth        ! Width of the Matrix (cm)
         READ(File_Parameters,*) ParticleSize       ! The size of particles .. grain size (length or diameter; cm)
         READ(File_Parameters,*) SedimentationRate  ! Sedimentation rate (e.g., w = 0.1, 0.01, 0.001 cm/yr)
         READ(File_Parameters,*) ParticleDensity    ! Average density of sedimenting particles (~ 2.5 g/cm3)
         READ(File_Parameters,*) Porosity0          ! % water at surface (0.8 = 80%)  .. the porosity
         READ(File_Parameters,*) Halflife           ! decay constant (22.3 yr for 210Pb)
         READ(File_Parameters,*) Depth_DeadZone     ! Depth at which organisms have ~0 Probability of being found (cm)
         READ(File_Parameters,*) r_TimeMin          ! OUTPUT begins at time = n
         READ(File_Parameters,*) r_TimeMax          ! Endtime (days)
         READ(File_Parameters,*) r_TimeStep         ! Time steps ... Frequency of data output (days)
       CLOSE(File_Parameters)

   ! preliminary reading of organism data to obtain scaling parameters
       tmp_MaxOrgSize = 0
       MaxOrgSpeed = 0
       MinWidth = ParticleSize
       MaxWidth = 0

       MissingValue= -1.
       Porosity_Type = 'UNI'     ! The kind of porosity function: UNIform, RANdom, LINear, EXPonential
       Porosity_DecayRate = -1   ! The extinction coefficient of the porosity change with depth (-1 if N/A)
       POROSITY_THRESHOLD = 0.9  ! THE CUTOFF AT WHICH ORGANISMS WILL TURN AWAY FROM THE WATER COLUMN
       WindowSize = 2            ! the number of pixels used for calculation of local porosity (function Porosity_local)

     OPEN(unit = File_Organisms, file = 'Organisms.IN', status = 'old')
       DO n = 1, N_Ind
         READ(File_Organisms,*) i
         Read(File_Organisms,*) r_Width
         Read(File_Organisms,*) r_Length
         Read(File_Organisms,*) r_Density
         Read(File_Organisms,*) FeedingType
         Read(File_Organisms,*) r_Speed   ! cm/day
         Read(File_Organisms,*) r_IngestRate   ! /day
         Read(File_Organisms,*) r_IngestSelectivity

         tmp_OrgSize = r_Width * r_Length
           If (tmp_OrgSize .gt. tmp_MaxOrgSize) then
             tmp_MaxOrgSize  = tmp_OrgSize
             tmp_Width       = r_Width
             tmp_Length      = r_Length
           End if
           If (r_Width .lt. MinWidth) then
             MinWidth = r_Width         
           End if
           If (r_Width .gt. MaxWidth) then
             MaxWidth = r_Width
           End if
           IF (r_Speed .gt. MaxOrgSpeed) then
             MaxOrgSpeed = r_Speed
           End if
       END Do
     CLOSE(File_Organisms)

   ! The ** critical scaling parameter for spatial dimensions **
   ! If the size of particles specified is larger or equal to the size of the head of the
   ! smallest organism, the simulation will continue but with particle size reduced to
   ! be (MinHeadWidth) times smaller than the smallest head width (PixelSize)

       If (MinWidth .lt. 2*ParticleSize) then  ! head is smaller than particles .. 
         MinHeadWidth = 2    ! defines the number of pixels of the smallest head
         PixelSize = MinWidth / Real(MinHeadWidth)   ! cm/pixel
       Else
         PixelSize = ParticleSize
       End IF

   ! redefine/rescale the spatial dimensions of the variables
       N_Row = CEILING(MatrixDepth/ PixelSize)    ! the number of rows of the matrix
       N_Col = CEILING(MatrixWidth/ PixelSize)    ! the number of columns of the matrix
       N_Cell =  N_Row * N_Col           ! the number of cells in the matrix
       Buffer_Zone = CEILING(MaxWidth/PixelSize)  ! depth above the sediment interface that organisms are seeded
       N_RowWater = 3 * N_Row / 10
       N_RowSed = N_Row - N_RowWater
       DepthToConstrain = N_RowWater    ! this is the depth to which the totalnumber of cells are constrained

       DepthDeadZone = Ceiling(Depth_DeadZone / PixelSize)  ! the inverse .. used for calculations elsewhere

   ! the ** critical scaling parameter for time **
   ! this gives the length of real time (days) for each unit of simulation time (RealTime/SimulTime)
       TimeScale = PixelSize / MaxOrgSpeed       !  day

   ! define some useful time constants

         Day   =   1 / TimeScale + Int(2.*mod(1.,Timescale))
         Year  = 365 / TimeScale + Int(2.*mod(365.,Timescale))
		 
      N_Outputs = 15
      If (oxygen_ON) N_Outputs = 73*3
      N_Outputs2 = 100000

    Allocate (Time_output(N_Outputs))
    Allocate (Time_output2(N_Outputs2))
   
   if (.not. oxygen_ON) then 
     Time_Output =(/1,25*Day,50*Day,100*Day,150*Day,200*Day,250*Day,1*Year,2*Year,3*Year,5*Year,10*year,20*year,30*Year,100*year/)
	 end if 
	 
     do i = 1, N_Outputs2
        Time_Output2(i) = i
     end do
     
	 If (oxygen_ON) then
     do i = 1, N_Outputs
        Time_Output(i) = 5000*i
     end do
	 end if 

   ! rescale time to units of simulation
       TimeMin =  r_TimeMin / TimeScale + Int(2.*mod(r_TimeMin,Timescale))
       TimeMax =  r_TimeMax / TimeScale + Int(2.*mod(r_TimeMax,Timescale))
       TimeStep = r_TimeStep / TimeScale + Int(2.*mod(r_TimeStep,Timescale))  ! only eefective when movies on 

   ! Rescale sedimentation rate:
   !   no of pixels in 1cm = 1/pixelSize
   !   no of pixels in 1cm2 = 1/(pixelSize^2)
   !   no particles in 1cm2 = (1-porosity) * no pixels in 1cm2 = ((1-Porosity) * 1/pixelsize^2)
   !   the n cm of sediments sedimented in 1 yr = sedrate * 1 yr
   !   this represents and area = Sedrate*MatrixWidth
   !   therefore, the no. of particles to sediment in 1 yr
   !             = no of particles in that area
   !             = Sedrate*MAtrixWidth * no particles in 1cm2
   !             = SedRate*MatrixWidth * ((1-Porosity) * 1/pixelsize^2)

       SedRate = ((1.-Porosity0) /(PixelsIZE*PixelSize)) * SedimentationRate * MatrixWidth
!                                   cm/pixel   cm/pixel     cm/yr               cm
!         this calculate the pixel number sedminentaed per year 
   ! The frequency at which the sedimentation of particles is to be made
       ParticlesToSediment_y = SedRate
       Area_Total            = Year / 2              ! area under the sinusoidal in one year
         If (SedimentationRate .le. 0.01) then       ! update yearly
           Time_Sediment   =  Year
         Else If (SedimentationRate .gt. 0.01) then  ! else update every quarter year
           Time_Sediment   =  Year / 4
         End if


     MaxOrgSize = (1 + (tmp_Width + PixelSize)/PixelSize) * (1 + (tmp_Length+PixelSize)/PixelSize)
!  YK  why adding 1? does not matter because it is for maximum?

     ALLOCATE (Matrix(N_Row, N_Col), MatrixOccupancy(N_Row))
     Allocate (Particle(N_Cell))
     Allocate (Org(N_Ind))
     Allocate (Org_Loc(MaxOrgSize, N_Ind))
     Allocate (IngestionHistory(N_Ind, Day), MovementHistory(N_Ind, Day))
     
     Allocate (RespHistory(N_Ind, Day))      
     allocate (O2(N_row,n_col))
     allocate (matrix_chk(N_row,n_col))
	 allocate (Dir_rec(Day,N_IND))
	 allocate (Ub(2,N_ind), Vb(2,N_ind)) !  in flow subroutines, coodinate is listed in order of x, y
	 allocate (flow_loc(2,N_ind)) !  in flow subroutines, coodinate is listed in order of x, y
	 allocate (Ug(N_col+1,N_row+2),Vg(N_col,N_row+2+1), Pg(N_col,N_row+2), Dg(N_col,N_row+2)) 
	 allocate (Uo(N_col,N_row),Vo(N_col,N_row))
	 allocate (edif(n_row, n_col)) 

   ! Read in again organism info and rescale organism traits to simulation time/space units
       OPEN(unit = File_Organisms, file = 'Organisms.IN', status = 'old')
       DO n = 1, N_Ind

         READ(File_Organisms,*) i
         Read(File_Organisms,*) r_Width
         Read(File_Organisms,*) r_Length
         Read(File_Organisms,*) r_Density
         Read(File_Organisms,*) FeedingType
         Read(File_Organisms,*) r_Speed
         Read(File_Organisms,*) r_IngestRate
         Read(File_Organisms,*) r_IngestSelectivity

         Org(i)%Density                = r_Density
         Org(i)%Pr_GoStraight          = 0.75
         Org(i)%FeedingType            = FeedingType
         Org(i)%Length                 = Ceiling((r_Length + PixelSize) / PixelSize)
         Org(i)%Width                  = Ceiling((r_Width + PixelSize) / PixelSize)
         Org(i)%HeadSize               = Org(i)%Width * Org(i)%Width
         Org(i)%BodySize               = Org(i)%Width * Org(i)%Length
         Org(i)%Orientation            = 0
         Org(i)%SensoryRange           = Org(i)%Width      ! this will eventually be made a user input
         Org(i)%FaecesWidth            = Ceiling(Real(Org(i)%Width) / 4.)
         Org(i)%Is%Reversing           = .false.
         Org(i)%Is%WrappingAround      = .false.
         Org(i)%Is%OnTheEdge           = .false.

         Org(i)%Can%Egest              = .true.
         Org(i)%Can%Ingest             = .false.
         Org(i)%Can%Move               = .true.

         Org(i)%Gut%Capacity           = Org(i)%BodySize / 4
         Org(i)%Gut%Content            = 0
         Org(i)%Gut%Pb_Activity        = 0

         Org(i)%Move%Rate              = r_Speed * TimeScale / PixelSize  ! number of pixels per unit simulation time
         Org(i)%Move%Distance          = 0
         Org(i)%Move%Bearing_Distance  = 0
         Org(i)%Move%Bearing           = 0.0
         Org(i)%Move%ReverseCount      = 0

         MovementHistory(i,:)          = 0

         Org(i)%Move%StoppedTime       = 0

         Org(i)%Ingest%Rate        = r_IngestRate*Org(i)%Bodysize*Org(i)%Density*Timescale/ParticleDensity ! n particles per day
         Org(i)%Ingest%Selectivity = r_IngestSelectivity
         Org(i)%Ingest%Amount      = 0

         IngestionHistory(i,:)     = 0
         
         respHistory(i,:)     = 0

       End do

       CLOSE(File_Organisms)

       MatrixOccupancy = 0   ! initial the frequency of occurence counts of organisms in the matrix

       MaxGutCapacity = 1
       do i=1, n_ind
         if (MaxGutCapacity .lt. Org(i)%Gut%Capacity) MaxGutCapacity = Org(i)%Gut%Capacity
       End do

       DecayConstant = -log(0.5) / HalfLife   ! the decay constant for 210Pb   !  /yr

   ! transition probabilities of lability change (1) == 1 to 0, (2) == 2 to 1, etc. .. as halflife in year**(-1)
       N_LabilityClasses = 10
       Allocate (Lability_decayConstant(N_LabilityClasses), Lability_proportion(N_LabilityClasses))
       
       do j = 1, N_labilityClasses
         lability_decayConstant(j) = 1e-1   ! /yr for labile OM cf., Canfield, 1994
		 ! lability_decayConstant(j) = lability_decayConstant(j)*1e1   ! when high rate is considered 
         lability_proportion(j) = j 
       end do

       Sum_Lability = Real(Sum(Lability_proportion))

      ! rescale and make them cumulative values
         Lability_tmp = 0
         Do j = 1, N_LabilityClasses
           Lability_tmp = Lability_tmp + Lability_proportion(j)
           Lability_proportion(j) = Real(Lability_tmp) / Sum_Lability
         End do
         
   ! Fill the matrix with particles and organisms
       Call Matrix_fill()

       Allocate (Guts(MaxGutCapacity, N_ind))

       Guts = CellContent(w,w)    ! initialise the guts
       Org_Loc = Coordinates(0,0)
	   
       Do i = 1, N_Ind
         Call Matrix_populate(i)
       End do

		 flow_loc = coordinates(0,0)
		 Vb = 0.0
		 Ub = 0.0

   End Subroutine GetUserInput


   ! ********************************************************************


   Subroutine Timetable(i)
     ! process information of the acitivity of organism i

   Use GlobalVariables
   implicit none
   Integer, Intent(in):: i
   REal IngestRAte, MovementRate, URand(3), Fullness, Pr_Activity
   integer :: j
   Character*21 numtemp
   
   real :: RespRate_ave
   integer :: k
   real :: ave_OM, ptcl_num
   
   ! shift elements right
       MovementHistory(i,:) = EOShift (MovementHistory(i,:), -1)
       IngestionHistory(i,:) = EOShift (IngestionHistory(i,:), -1)
       RespHistory(i,:) = EOShift (RespHistory(i,:), -1)
       
   ! calculate moving averages of the rates
       MovementRate = Real(Sum(MovementHistory(i,:))) / Day
       IngestRAte = Real(Sum(IngestionHistory(i,:))) / Day  !  n particle per single iteration 
       RespRate_ave = Real(Sum(RespHistory(i,:))) / Day
	   
   ! check ingestion rates
       If (MovementRate .eq. 0) then               ! not moving
         If (Org(i)%Move%StoppedTime .eq. 0) Org(i)%Move%StoppedTime = Time  ! check stomach (i.e., caloric) intake is within tolerance

         if((Time - Org(i)%Move%StoppedTime) .gt. (Year*.25)) then ! if not moving for a long time  .. let it die
 
! this is a fallback routine to catch any anomalies
           Call Organism_dead(i) 
! dead .. remove and start a new organism
           Org(i)%Move%StoppedTime = 0                              ! reset the counter
           Return
         End if

       End if

   ! if there are no problems of being stuck, etc. (above) ...
   ! calculate the activity of the organism = Pr_Activity
       Call Random_Number(URand)
	   
       Org(i)%Can%Move = .false.
       Org(i)%Can%Ingest = .false.
       Org(i)%Can%Egest = .false.

       ! assume that the amplitude of yearly variation is 10 time greater than the daily cycles
       ! an additional 0.2 is added to the probability to make the lower threshold atleast 20% activity
       ! other, organism-specific cycles can be added here ...
           Pr_Activity = Pr_activity_year * 0.9  + Pr_Activity_day * 0.1  + 0.2
		   
             If (URand(1) .lt. Pr_Activity)  then
             If (MovementRate .lt. Org(i)%Move%Rate) then      ! if within tolerance ...
               Org(i)%Can%Move = .true.
               If (IngestRate .lt. Org(i)%Ingest%Rate) then    ! if within tolerance ...
                 Fullness = Real(Org(i)%Gut%Content)/Real(Org(i)%Gut%Capacity) ! calculate gut fullness
				 if (resp_ON) then 
				   ave_OM = 0.
				   ptcl_num = 0.
				   do k = 1, org(i)%gut%capacity
				     if (guts(k,i)%class == p) then 
				     ptcl_num =  ptcl_num + 1.
					 ave_OM = ave_OM + particle(guts(k,i)%value)%OM%OMact
					 end if 
				   end do 
				   if (ptcl_num>0) then 
				     ave_OM = ave_OM/ptcl_num
				   end if 
				   
                   IF ( ave_OM < 0.1) then     ! if OM is low in guts
				     Org(i)%Can%Egest = .true.       
                     if (fullness < 1.) Org(i)%Can%Ingest = .true.
                   end if 					 
				 end if 
				 
                 if (.not.resp_ON) then 
                   IF (URand(2) .lt. Fullness) Org(i)%Can%Egest = .true.         ! too full .. defecate
                   If (Urand(3) .gt. Fullness) Org(i)%Can%Ingest = .true.        ! too empty .. eat
				 end if 
               End if
             End if
             End if
			 
   End Subroutine Timetable


   ! ******************************************************************


   SUBROUTINE Rules_of_Motion(i)
     ! the rules of motion for deposit feeder i ..

   Use GlobalVariables
   Implicit None
   Integer, Intent(In) :: i
   Integer :: Possible(4), Direction(4), DirectionToGo, j
   Real    :: RandomDir(4), ToGo(4), PREFERABLE(4), ToAvoid(4), Bearings(4), RandomWeight, URand(4)

   ! As the organism is about to move .. check for exceptions/boundaries
       Call Random_Number(URand)

       j = Int(URand(1) * Real(Org(i)%Width)) + 1   !  why using random nmber? to randomize calling organism_relocate?
       Org(i)%Is%OnTheEdge = .false.
         IF ((Org_loc(j,i)%X .eq. N_Col) .or. (Org_loc(j,i)%X .eq. 1)) THEN
           Org(i)%Is%OnTheEdge = .true.
           Call Organism_relocate(i)
           Org(i)%Can%Move = .false.
           Return
         End if

       If (Org(i)%Is%Reversing) then
         Call Organism_reverse(i)
         Org(i)%Can%Move = .false.
         Return
       End if

       IF (Org(i)%Is%WrappingAround) then          ! first time through
         Call Organism_relocate(i)
         Org(i)%Can%Move = .false.
         Return
       End if

   ! define the cardinal directions  ... 0=up, 90=right, 180=down, 270=left
       Direction = (/0, 90, 180, 270/)

   ! randomly weight the directions
       RandomWeight = 0.25
       RandomDir = URand * RandomWeight

   ! call rules of motions .. sumarised in the variable "ToGo"
       Call Rule_possible(i, Possible, Direction)    ! Return the possible directions
       Call Rule_avoid(i, Possible, ToAvoid)
       Call Rule_bearings(i, Bearings, Direction)  ! Randomise choice of directions a bit
         ToGo = Real(Possible) * RandomDir * ToAvoid * Bearings

       If (Lability_ON) then
         Call Rule_labilityCheck(i, ToGo, Preferable, Direction)
         ToGo = ToGo * Preferable
       End if

       If (Sum(ToGo) .eq. 0) then
         Call Organism_reverse(i)
         Return
       End if

       DirectionToGo = Direction(MaxLoc(ToGo, Dim=1))

       Call Head_rotate(i, DirectionToGo)

   End SUBROUTINE Rules_of_Motion

   ! ************************************************

   Subroutine Head_locate(i, Head_MinVal, Head_MaxVal)
     ! find the co-ordinates of the corners of the head: (miny, minx) and (maxy, maxx)

   Use GlobalVariables
   Implicit None
   Integer i
   Type(Coordinates) :: Head_MinVal, Head_MaxVal

   Select Case (Org(i)%Orientation)
     Case (0)
       Head_MinVal = Coordinates(Org_Loc(1,i)%Y, Org_Loc(1,i)%X)
       Head_MaxVal = Coordinates(Org_Loc(1,i)%Y + (Org(I)%Width-1), Org_Loc(1,i)%X + (Org(I)%Width-1)) 
     Case (90)
       Head_MinVal = Coordinates(Org_Loc(1,i)%Y, Org_Loc(1,i)%X - (Org(I)%Width-1))  
       Head_MaxVal = Coordinates(Org_Loc(1,i)%Y + (Org(I)%Width-1), Org_Loc(1,i)%X)
     Case (180)
       Head_MinVal = Coordinates(Org_Loc(1,i)%Y - (Org(I)%Width-1), Org_Loc(1,i)%X - (Org(I)%Width-1))
       Head_MaxVal = Coordinates(Org_Loc(1,i)%Y, Org_Loc(1,i)%X)
     Case (270)
       Head_MinVal = Coordinates(Org_Loc(1,i)%Y - (Org(I)%Width-1), Org_Loc(1,i)%X)
       Head_MaxVal = Coordinates(Org_Loc(1,i)%Y, Org_Loc(1,i)%X + (Org(I)%Width-1))    
   End Select
  
   End Subroutine Head_locate

   ! ************************************************

   Subroutine Rule_possible(i, Possible, Direction)
     ! Determine which directions are possible for motion

   Use GlobalVariables
   Implicit None
   Integer   :: i, j, N_Particle, N_Organism, N_Water
   Integer, Intent(In) :: direction(4)
   Integer, Intent(Out) :: Possible(4)
   Type(Coordinates) :: Point
   
   Possible = (/1, 1, 1, 1/)

   Do j = 1, 4
     Call Mouth_locate(i, Point, direction(j))
     Call Block_read(i, Point, direction(j), N_Particle, N_Organism, N_Water)
       If (Block_outsideMatrix)  Possible(j) = 0
       If (N_Organism .gt.0)     Possible(j) = 0
   End Do

   End Subroutine Rule_possible

   ! ****************************************************************************

   Subroutine Rule_avoid(i, Possible, ToAvoid)
     ! to avoid any direction due to being on the extremes of the matrix

   Use GlobalVariables
   Implicit None
   Integer, Intent(In) :: i, Possible(4)
   Real, Intent(Out)   :: ToAvoid(4)
   Type(CoOrdinates)   ::  Head_MinVal, Head_MaxVal
   Real                :: Porosity_local
   logical             :: Isfloating

     ToAvoid = (/1., 1., 1., 1./)  ! nothing to avoid yet
     Call Head_locate(i, Head_MinVal, Head_MaxVal)

     ! Can we go deeper ?
       If (Possible(3) .ne. 0) then
         If (Head_MaxVal%Y .eq. N_Row) ToAvoid = ToAvoid * (/1., 1., 0.0001, 1./)
         If (DepthAvoidance_on) then
           if (Head_MaxVal%Y .gt. (N_RowWAter+1)) ToAvoid = ToAvoid *        &
               (/1., 1., 1. - ( Real(Head_MaxVal%Y-N_RowWAter)/Real(DepthDeadZone)), 1./)
                ! linear decrease in Pr(going deeper)
         End if
       END IF

     ! Can we go higher ?
       If (Possible(1) .ne. 0) then
         If (Head_MinVal%Y .eq. 1) ToAvoid = ToAvoid * (/0.0001, 1., 1., 1./)
         If (Head_MinVal%Y .Lt. N_RowWater) then
           If ( Porosity_local ( CoOrdinates ( Head_MinVal%Y-Org(I)%Width-1, Head_MinVal%X+(Org(I)%Width/2) ), Org(I)%Width) &
             .gt.  POROSITY_THRESHOLD )   ToAvoid = ToAvoid * (/0.01, 1., 1., 1./)
         End if
		 
		 If (I_shape .and. i == Org_ID_ishape) then   
		   if (Isfloating(i)) ToAvoid = ToAvoid * (/0.0001, 0.0001, 1., 0.0001/)
		 end if 
		 
       End if
       
   End Subroutine Rule_avoid

   ! ****************************************************************************

   Subroutine Rule_bearings(i, Bearings, Direction)
     ! Shall we change our bearings ?

   Use GlobalVariables
   Implicit None
   Logical :: GetNewBearings
   Integer     :: j, minVal, DirectionToPoint
   Integer, Intent(In)  :: i, Direction(4)
   Real, Intent(Out) :: Bearings(4)
   Real  :: URand(2)
   
   Bearings = (/0.1, 0.1, 0.1, 0.1/)
   GetNewBearings = .false.
   Call Random_Number(URand)
   
     If (URand(1) .gt. Org(i)%Pr_GoStraight) GetNewBearings = .true.


     If (GetNewBearings) then
       Org(i)%Move%Bearing = Org(i)%Move%Bearing + Real(Int(URand(2)*360 + 1 - 180))
       Org(i)%Move%Bearing_Distance = 0
       IF (Org(i)%Move%Bearing .ge. 360) then
         Org(i)%Move%Bearing = Org(i)%Move%Bearing - 360.
       Else if (Org(i)%Move%Bearing .lt. 0) then
         Org(i)%Move%Bearing = Org(i)%Move%Bearing + 360.
       End if
     End if

     minVal = 360
     Do j = 0, 270, 90
       if (abs(j-int(Org(i)%Move%Bearing)) .le. minVal) then
         minVal = abs(Real(j)-Org(i)%Move%Bearing)
         DirectionToPoint = j
       End if
     End do

     WHERE (Direction .eq. DirectionToPoint) Bearings = 1

   End Subroutine Rule_bearings

   ! ***************************************************************

   Subroutine Rule_labilityCheck(i, ToGo, Preferable, Direction)
     ! Scan regions to determine which direction is preferable based upon particle lability

   Use GlobalVariables
   Implicit None
   Integer             :: j, Y, X, index
   Integer             :: N_Water(4), N_Particles(4), N_Organism(4)
   Integer             :: Lability_sum(4)
   Integer, Intent(In) :: i, Direction(4)
   Real, Intent(Out)   :: PREFERABLE(4)
   Real, Intent(In)    :: ToGo(4)
   Real                :: Lability_mean(4), O2_sum(4), O2_mean(4)
   Type(Coordinates)   :: Head_MinVal, Head_MaxVal
   
   real                ::  Fullness    
   real                ::  water_mean(4)   

   ! initialise the counts
       Lability_sum  = 0
       Lability_mean = 0
       N_Water   = 0
       N_Organism  = 0
       N_Particles = 0
       O2_sum = 0
       O2_mean = 0

     Call Head_locate(i, Head_MinVal, Head_MaxVal) ! obtain coordinates of (miny,minx) and (maxy, maxx)

     Do j = 1, 4
       If (ToGo(j) .eq. 0) Cycle

       index = -1                   ! index makes the search region an isoceles triangle

       Select Case (Direction(j))
         Case (0)
         Do Y = (Head_MinVal%Y-1), Head_MinVal%Y-Org(i)%SensoryRange, -1
           index = index+1
           Do X = (Head_MinVal%X-index), (Head_MaxVal%X+index), 1
             Call TabulateInfo(Y, X, j, Lability_sum, N_Particles, N_Water, N_Organism, O2_sum)
           End do
         End do

         Case (90)
         Do X = (Head_MaxVal%X+1), (Head_MaxVal%X+1) + Org(i)%SensoryRange, 1
           index = index+1
           Do Y = (Head_MinVal%Y-index), (Head_MaxVal%Y+index), 1
             Call TabulateInfo(Y, X, j, Lability_sum, N_Particles, N_Water, N_Organism, O2_sum)
           End do
         End do

         Case (180)
         Do Y = (Head_MaxVal%Y+1), (Head_MaxVal%Y+1) + Org(i)%SensoryRange, 1
           index = index+1
           Do X = (Head_MinVal%X-index), (Head_MaxVal%X+index), 1
             Call TabulateInfo(Y, X, j, Lability_sum, N_Particles, N_Water, N_Organism, O2_sum)
           End do
         End do

         Case (270)
         Do X = (Head_MinVal%X-1), (Head_MinVal%X-1) - Org(i)%SensoryRange, -1
           index = index+1
           Do Y = (Head_MinVal%Y-index), (Head_MaxVal%Y+index), 1
             Call TabulateInfo(Y, X, j, Lability_sum, N_Particles, N_Water, N_Organism, O2_sum)
           End do
         End do

       End Select

       Lability_mean(j) = Real(1 + Lability_sum(j)) / Real(N_Particles(j) + 1)
       O2_mean(j) = Real(1 + O2_sum(j)) / Real(N_water(j) + 1)
         
	   water_mean(j) = real(N_water(j)/(N_water(j)+N_particles(j)+N_organism(j)))

     End do        ! j index

     Preferable = Lability_mean
     if (oxFB_ON) Preferable = Lability_mean*O2_mean    !!  if considering oxygen_feedback on organism movement
     ! if (oxFB_ON) Preferable = O2_mean    !!   oxygen_feedback alone 
	 
	 If ((I_shape) .and. i == Org_ID_ishape) then 
	   Fullness = Real(Org(i)%Gut%Content)/Real(Org(i)%Gut%Capacity) ! calculate gut fullness
	   preferable = preferable*(water_mean*fullness + (1.0 - fullness)*(1.-water_mean))
     end if 
	 
   ! avoidance of other organisms
     Where (N_Organism .gt. 0) Preferable = Preferable * 0.0001 / Real(N_Organism)
     
   End Subroutine Rule_labilityCheck

   ! ****************************************************************************

   SUBROUTINE TabulateInfo(Y, X, j, Lability_sum, N_Particles, N_Water, N_Organism, O2_sum)
     ! low level routine that tabulates numbers of particles ,etc. for a given co-ordinate (Y,X)

   Use GlobalVariables
   Implicit None
   Logical :: InMatrix   
   Integer, intent(in) :: Y, X, j
   Integer, intent(inout) :: N_Particles(4), N_Water(4), N_Organism(4), Lability_sum(4)
   real, intent(inout) :: O2_sum(4)  
   real     :: Lab_real_New        
   
   If (InMAtrix(CoOrdinates(Y,X))) then
     Select Case (Matrix(Y,X)%Class)
       Case(p)
         N_Particles(j) = N_Particles(j) + 1
           IF (Lab_real_New(Matrix(Y,X)%Value) .le. 0) then
           Else
             Lability_sum(j) = Lability_sum(j) + int(Lab_real_New(Matrix(Y,X)%Value))
           End if
       Case(w)
         N_Water(j) = N_water(j) + 1
       Case(1:)
         N_Organism(j) = N_Organism(j) + 1
     End Select
     O2_sum(j) = O2_sum(j) + O2(Y,X)%oxygen
   End if
   
   End SUBROUTINE TabulateInfo

   ! ******************************************************************

   SUBROUTINE Head_rotate(i, DirectiontoGo)
     ! low-level routine that converts the absolute direction to a direction relative the orientation of the head

   Use GlobalVariables
   Implicit None
   Integer, Intent(In) :: i, DirectiontoGo
   
   Select Case (Org(i)%Orientation - DirectiontoGo)
     Case (0)
       Return
     Case (180, -180)
	   call rotate(i, 180)  
     Case (90, -270)
       Call Rotate (i, -90)
     Case (-90, 270)
       Call Rotate (i, 90)
   End Select
		 
   End Subroutine Head_rotate

   ! ***************

   Subroutine Rotate(i, Angle)
     ! low level routine to rotate the head a relative angle

   Use GlobalVariables
   Implicit None
   Integer :: j
   Integer, Intent(In) :: i, Angle
   Type(CellContent) :: Image_ID(Org(i)%Width * Org(i)%Width)
   Type(Coordinates) :: Image(Org(i)%Width * Org(i)%Width)
   Type(Coordinates) :: Head_MinVal, Head_MaxVal
   
   integer :: xnew, ynew, xnew_2, ynew_2
   
   Call Head_locate(i, Head_MinVal, Head_MaxVal)
    
	if (Head_MaxVal%X > n_col) then  !!  wrapping around
	
     do j = 1, Org(i)%HeadSize
       Image_ID(j) = Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X)
       Image(j) = Org_Loc(j,i)
	   if (image(j)%X < Org(i)%width ) image(j)%X = image(j)%X + n_col
     end do
	 
	else if (Head_MinVal%X < 1) then  !!  wrapping arround 
	
     do j = 1, Org(i)%HeadSize
       Image_ID(j) = Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X)
       Image(j) = Org_Loc(j,i)
	   if (image(j)%X > n_col - Org(i)%width) image(j)%X = image(j)%X - n_col
     end do
	 
	 else   !! mormal condition
	
     do j = 1, Org(i)%HeadSize
       Image_ID(j) = Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X)
       Image(j) = Org_Loc(j,i)
     end do
	 
	 end if 

   Select Case  (Angle)
     Case (-90)            ! rotate the head left (counter-clockwise) 90 degrees
       do j = 1, Org(i)%HeadSize
	     xnew = Image(j)%Y-Head_MinVal%Y+Head_MinVal%X
		 ynew = -(Image(j)%X-Head_MinVal%X)+Head_MinVal%Y+(Org(I)%Width-1)
		 if (xnew < 1) xnew = xnew + n_col
         if (xnew > n_col) xnew = xnew - n_col		 
         Org_Loc(j,i) = Coordinates(ynew, xnew)       
		 
         Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = Image_ID(j)
       end do

       Org(i)%Orientation = Org(i)%Orientation-90
       If (Org(i)%Orientation .eq. -90) Org(i)%Orientation = 270

     Case (90)             ! rotate the head right (clockwise) 90 degrees
       do j = 1, Org(i)%HeadSize
	     xnew = - (Image(j)%Y-Head_MinVal%Y)+Head_MinVal%X+(Org(I)%Width-1)
		 ynew = (Image(j)%X-Head_MinVal%X)+Head_MinVal%Y
		 if (xnew < 1 ) xnew = xnew + n_col
		 if (xnew > n_col) xnew = xnew - N_col
         Org_Loc(j,i) = Coordinates(Ynew,xnew)     
         Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = Image_ID(j)
       end do

       Org(i)%Orientation = Org(i)%Orientation+90
         If (Org(i)%Orientation .eq. 360) Org(i)%Orientation = 0
		 
	 Case(180)
	   do j = 1, Org(i)%HeadSize
	     xnew = - (Image(j)%Y-Head_MinVal%Y)+Head_MinVal%X+(Org(I)%Width-1)
		 ynew = (Image(j)%X-Head_MinVal%X)+Head_MinVal%Y
		 xnew_2 = - (ynew-Head_MinVal%Y)+Head_MinVal%X+(Org(I)%Width-1)
		 ynew_2 = (xnew-Head_MinVal%X)+Head_MinVal%Y
		 
		 if (xnew_2 < 1 ) xnew_2 = xnew_2 + n_col
		 if (xnew_2 > n_col) xnew_2 = xnew_2 - N_col
		 
         Org_Loc(j,i) = Coordinates(Ynew_2,xnew_2)     
         Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = Image_ID(j)
		 
	  end do
	  
       Org(i)%Orientation = Org(i)%Orientation+180
         If (Org(i)%Orientation >= 360) Org(i)%Orientation = Org(i)%Orientation - 360
		 
   End Select

   End SUBROUTINE Rotate

   ! *****************************************************************

   Subroutine Mouth_locate(i, Point, Orientation)
     ! returns the leftmost coordinates of the mouth, relative to the orientation of the head

   Use GlobalVariables
   Implicit None
   Integer i, Orientation
   Type(CoOrdinates) :: Point, Head_MinVal, Head_MaxVal
   
   Call Head_locate(i, Head_MinVal, Head_MaxVal)
   
   Select Case (Orientation)

   Case (0, 180, 360)
     Select Case (Orientation)
       Case (0, 360)
         Point = CoOrdinates(Head_MinVal%Y-1,Head_MinVal%X)
       Case (180)
         Point = CoOrdinates(Head_MaxVal%Y+1,Head_MaxVal%X)
   End Select

   Case (-90, 90, 270)
     Select Case (Orientation)

     Case (90)
       Point = CoOrdinates(Head_MinVal%Y,Head_MaxVal%X+1)
        IF (Org(i)%Is%OnTheEdge .and. (.not. Org(i)%Is%WrappingAround)) Point = CoOrdinates(Head_MinVal%Y,1)

     Case (-90, 270)
       Point = CoOrdinates(Head_MaxVal%Y,Head_MinVal%X-1)
       IF (Org(i)%Is%OnTheEdge .and. (.not. Org(i)%Is%WrappingAround)) Point = CoOrdinates(Head_MaxVal%Y,N_Col)

   End Select

   End Select

   End Subroutine Mouth_locate

   ! ********************************

   Subroutine Block_read(i, Point, Orientation, N_Particle, N_Organism, N_Water)
     ! given a point and a orientation, inventory what is found in a block of dimensions: 1 X BW, (BW=OrganismWidth)
     ! and read into a matrix called "block" and their coordinates in a matrix called "block_loc"
     ! if they are inside the sediment matrix

   Use GlobalVariables
   Implicit None
   Integer i, X, Y, k, Orientation
   Integer N_Particle, N_Organism, N_Water
   Type(Coordinates) :: Point
   
   k     = 0
   N_Particle  = 0
   N_Organism  = 0
   N_Water = 0
   Block_outsideMatrix = .false.

   Select Case (Orientation)
     Case (0, 180, 360)

       Select Case (Orientation)
         Case (0, 360)
              Y = Point%Y
           DO X = Point%X, Point%X+(Org(I)%Width-1)
             Call GatherInfo(Y, X, N_Particle, N_Organism, N_Water)
           End do

           If (.not. Block_outsideMatrix) then
             If (Allocated (Block)) Deallocate (Block)
             If (Allocated (Block_loc)) Deallocate (Block_loc)
             Allocate (Block(Org(i)%Width), Block_loc(Org(i)%Width))

               DO X = Point%X, Point%X+(Org(I)%Width-1)
                 k = k+1
                 Block(k) = Matrix(y,X)
                 Block_loc(k) = Coordinates(y,X)
               END Do
           End if

         Case (180)
              Y = Point%Y
           DO X = Point%X, Point%X-(Org(I)%Width-1), -1
             Call GatherInfo(Y, X, N_Particle, N_Organism, N_Water)
           End do
           If (.not. Block_outsideMatrix) then
             If (Allocated (Block)) Deallocate (Block)
             If (Allocated (Block_loc)) Deallocate (Block_loc)
             Allocate (Block(Org(i)%Width), Block_loc(Org(i)%Width))

             DO X = Point%X, Point%X-(Org(I)%Width-1), -1
                 k = k+1
                 Block(k) = Matrix(y,X)
                 Block_loc(k) = Coordinates(y,X)
             End do
           End if

       End Select

     Case (-90, 90, 270)

       Select Case (Orientation)
         Case (90)
              X = Point%X
           DO Y = Point%Y, Point%Y+(Org(I)%Width-1)
             Call GatherInfo(Y, X, N_Particle, N_Organism, N_Water)
           End do
           If (.not. Block_outsideMatrix) then
             If (Allocated (Block)) Deallocate (Block)
             If (Allocated (Block_loc)) Deallocate (Block_loc)
             Allocate (Block(Org(i)%Width), Block_loc(Org(i)%Width))

             DO Y = Point%Y, Point%Y+(Org(I)%Width-1)
                 k = k+1
                 Block(k) = Matrix(Y,x)
                 Block_loc(k) = Coordinates(Y,x)
             End do
           End if

         Case (-90, 270)
              X = Point%X
           DO Y = Point%Y, Point%Y-(Org(I)%Width-1), -1
             Call GatherInfo(Y, X, N_Particle, N_Organism, N_Water)
           End do
           If (.not. Block_outsideMatrix) then
             If (Allocated (Block)) Deallocate (Block)
             If (Allocated (Block_loc)) Deallocate (Block_loc)
             Allocate (Block(Org(i)%Width), Block_loc(Org(i)%Width))

             DO Y = Point%Y, Point%Y-(Org(I)%Width-1), -1
                 k = k+1
                 Block(k) = Matrix(Y,x)
                 Block_loc(k) = Coordinates(Y,x)
             End do
           End if

         End Select

   End Select
   
   End Subroutine Block_read

   ! ********************************

   Subroutine GatherInfo(Y, X, N_Particle, N_Organism, N_Water)
     ! low level routine to count matrix elements

   Use GlobalVariables
   Implicit None
   Logical Inmatrix  !YK  function
   Integer X, Y
   Integer N_Particle, N_Organism, N_Water
   
     if (InMatrix(CoOrdinates(Y,X))) then
       Select Case (Matrix(Y,X)%Class)
         Case (p)
           N_Particle  = N_Particle + 1
         Case (w)
           N_Water     = N_Water + 1
         Case (1:)
           N_Organism  = N_Organism + 1
       End Select
     Else
       Block_outsideMatrix = .true.
     End if
     
   End Subroutine GatherInfo

   ! ******************************

   Subroutine Organism_move(i)
     ! move the organism i by translation in the direction that it is oriented

   Use GlobalVariables
   use ieee_arithmetic
   implicit none
   Integer ::  i, j, dX, dY
   Type(Coordinates) :: Image(Org(i)%BodySize)
   Type(Coordinates) :: Image_O2(Org(i)%BodySize)
   
   If (Org(i)%Can%Move) then
   
     If (Allocated(tail)) Deallocate(tail)
     Allocate (Tail(Org(i)%Width))

     MovementHistory(i,1) = 1
     Org(i)%Move%Bearing_Distance = Org(i)%Move%Bearing_Distance + 1

     Call LookupDirection(Org(i)%Orientation, dy, dX)

       ! take an image of the body
         DO j =1, Org(i)%BodySize
           Image(j) = Org_Loc(j,i)
         END do
       ! move the head
         DO j =1, Org(i)%HeadSize                                    ! move the head
           Org_Loc(j,i) = Coordinates(Image(j)%Y+dY, Image(j)%X+dX)
           Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = CellContent(j,i)
           If (oxygen_ON) O2(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = O2(image(j)%Y, image(j)%X)
         END Do
       ! move the rest of the body
         DO j = Org(i)%HeadSize+1, Org(i)%BodySize
           Org_Loc(j,i) = Coordinates(Image(j-Org(I)%Width)%Y, Image(j-Org(I)%Width)%X)   !YK  translation of previous imange
           Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = CellContent(j,i)
           If (oxygen_ON) O2(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = O2(image(j)%Y, image(j)%X)
         end do
       ! finish the organism's tail-fill empty cells with water
         Do j = 1, Org(I)%Width
           Tail(j) = Image(Org(i)%BodySize-Org(I)%Width+j)
           Matrix(Image(Org(i)%BodySize-Org(I)%Width+j)%Y, Image(Org(i)%BodySize-Org(I)%Width+j)%X) = CellContent(w,w)
           If (oxygen_ON) O2(tail(j)%Y, tail(j)%X) = O2(image(j)%Y, image(j)%X)  
           !! exchange oxygen between occupied and emptied spaces by the movement of organism
         End do

     End if
	 
	 
       Dir_rec(:,i) = EOShift (Dir_rec(:,i), -1)
	   Dir_rec(1,i) = Org(i)%Orientation
	   
	   select case(Org(i)%Orientation)
	     case(0)
		   Vb(1,i) = Vb(1,i) + org(i)%width
	     case(180)
		   Vb(1,i) = Vb(1,i) - org(i)%width
	     case(90)
		   Ub(1,i) = Ub(1,i) - org(i)%width
	     case(270)
		   Ub(1,i) = Ub(1,i) + org(i)%width
		 case default
		   print *, 'error in setting flow const---moving',org(i)%orientation
		 end select
	 
     
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after org_move(sub)"; write(file_log,*)time, "err after org_move(sub)"; end if 
		 if ((any(Org_loc%X<0)) .or. (any(Org_loc%Y <0)) .or. (any(Org_loc%Y > N_row)) .or. (any(Org_loc%Y > N_row))) then
		   print *, 'errrrrrrr in org_loc'
		 end if 
         end if 
     
   End Subroutine Organism_move

! ****************************************************************

   Subroutine LookupDirection(MoveOrient, dy, dX)
     ! Covert direction to move into translational vectors (dY,dX)

   Integer, Intent (in)  :: MoveOrient
   Integer, Intent (out) :: dy, dx

   
     Select Case (MoveOrient)  ! (dX =+1 is right;  -1 is left);  (dY =+1  is down;-1 is up)
         Case(0, 360)
           dX = 0; dy = -1
         Case(180)
           dX = 0; dy = 1
         Case(90)
           dX = 1; dy = 0
         Case(270, -90)
           dX = -1; dy = 0
     End Select

   End Subroutine LookupDirection

  ! ********************************

   Subroutine Organism_relocate(i)
     ! the wraparound routine that moves organisms going through the lateral boundaries to the other side
     ! this step tests and prepares the mouth and tail region ..
     ! the final movement is completed with another routine: Organaism_wraparound(i)

   Use GlobalVariables
   Implicit None
   Logical :: OnEdge        ! YK function 
   Integer :: j, X_new, DirectionToGo, dX
   Integer, Intent(IN) :: i
   Type(CoOrdinates) :: Point

   If (.not. Org(i)%Is%WrappingAround) then     ! first time around  

       If (Org_loc(1,i)%X .gt. N_Col/2) then ! to the right .. force move left
         DirectionToGo = 90
         X_new = 1
       Else
         DirectionToGo = 270
         X_new = N_col
       End if

       Call Head_rotate(i, DirectionToGo)
       Call Mouth_locate(i, Point, Org(i)%Orientation)
       Call Particles_ingest(i, Point, Org(i)%Orientation)
       Call Particles_push(i, Point, Org(i)%Orientation)

       If (Org(i)%Can%Move) then
         Call Organism_wrapAround(i)        ! the head tip has been moved .. shift body along
         Call Particles_egest(i)
           Org(i)%Is%OnTheEdge = .false.    ! reset flags -- ready for the next part
           Org(i)%Is%WrappingAround = .true.
       End if

   End if
   
       Call Mouth_locate(i, Point, Org(i)%Orientation)
       Call Particles_ingest(i, Point, Org(i)%Orientation)
       Call Particles_push(i, Point, Org(i)%Orientation)

       If (Org(i)%Can%Move) then
         Call Organism_wrapAround(i)         ! the head tip has been moved .. shift body along
         Call Particles_egest(i)


       Select Case(LocalMixing)
         Case (.false.) ! i.e.  a normal run
           Do j = 1, 3*Org(i)%HeadSize/2         ! scan the head and "neck" to see if it clear of the edge
             If (OnEdge(Org_Loc(j,i))) RETURN    ! not clear .. continue wrapping
           End do

           Org(i)%Is%WrappingAround = .FALSE.    ! finished wraparound of the head .. 
                                                 ! reset the flag and continue with normal operations

           If (Org_loc(1,i)%X .gt. N_Col/2) then ! looped around from the left ..
             dX = -1   
           Else
             dX = 1
           End if

           do j = 1, Org(i)%Gut%Capacity         ! set location of particles ..
             IF (Guts(j,i)%Class .eq. p) Particle(Guts(j,i)%Value)%Plane = Particle(Guts(j,i)%Value)%Plane + dX
           End do

         Case (.true.) !  if local mixing, the whole body must be scanned to make sure rms particle displacements are correct.
           Do j = 1, Org(i)%BodySize
             If (OnEdge(Org_Loc(j,i))) RETURN    ! not clear .. continue wrapping
           End do

           Org(i)%Is%WrappingAround = .FALSE.    ! finished wraparound of the head ..

           Return

         End Select

       End if
	   
	   
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after relocate"; write(file_log,*)time, "err after relocate"; end if 
		 if ((any(Org_loc%X<0)) .or. (any(Org_loc%Y <0)) .or. (any(Org_loc%Y > N_row)) .or. (any(Org_loc%Y > N_row))) then
		   print *, 'errrrrrrr in org_loc'
		 end if 
         end if 

   End Subroutine Organism_relocate

   ! ***********************

   Subroutine Organism_reverse(i)
     ! move the organism backwards

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Logical :: InMatrix                               
   Integer :: j, neck, dx, dy, Reversedirection
   Integer, Intent(IN) :: i
   Type(Coordinates) :: Image(Org(i)%BodySize), Point
   
   Type(Coordinates) :: Image_2(Org(i)%BodySize)    
   integer           :: RevDir           
   integer :: xx, yy, k, o, oo , xmx, xmn, ymx, ymn
   integer :: DirChk(4), PreDir, PreOri, rec_image(Org(i)%length-Org(i)%width)
   
   real :: V_tmp(2), U_tmp(2)
   
   
   neck = Org(i)%HeadSize+1    ! index of the neck area

   ! check if we can reverse .. i.e. is any part of the body on the edge ..
       Do j = 1, Org(I)%Width
         IF ((Org_Loc(j,i)%X .eq. 1) .or. (Org_Loc(j,i)%X .eq. N_col)) then
           Org(i)%Is%Reversing = .false.
           Return
         End if
       End do
       Do j = Org(i)%BodySize-(Org(i)%HeadSize+2*Org(I)%Width), Org(i)%BodySize
         IF ((Org_Loc(j,i)%X .eq. 1) .or. (Org_Loc(j,i)%X .eq. N_col)) then
           Org(i)%Is%Reversing = .false.
           Return
         End if
       End do
       
   Org(i)%Is%Reversing = .true.
   ! find direction to reverse by comparing location of tail
       If (Org_Loc(Org(i)%BodySize,i)%X .eq. Org_Loc(Org(i)%BodySize-(Org(I)%Width-1),i)%X) then ! 90 or 270
         If (.NOT. InMatrix(CoOrdinates(Org_Loc(Org(i)%BodySize,i)%Y,Org_Loc(Org(i)%BodySize,i)%X+1))) Return

         If (Matrix(Org_Loc(Org(i)%BodySize,i)%Y,Org_Loc(Org(i)%BodySize,i)%X+1)%Class .eq. i) then
           Reversedirection = 270
         Else
           Reversedirection = 90
         End if

       Else if (Org_Loc(Org(i)%BodySize,i)%Y .eq. Org_Loc(Org(i)%BodySize-(Org(I)%Width-1),i)%Y) then  ! 0 or 180
         If (.NOT. InMatrix(CoOrdinates(Org_Loc(Org(i)%BodySize,i)%Y+1,Org_Loc(Org(i)%BodySize,i)%X))) Return

         If (Matrix(Org_Loc(Org(i)%BodySize,i)%Y+1,Org_Loc(Org(i)%BodySize,i)%X)%Class .eq. i) then
           Reversedirection = 0
         Else
           Reversedirection = 180
         End if

       End if
       
     Call LookupDirection (Reversedirection, dy, dx)

     Point = CoOrdinates(Org_loc(Org(i)%BodySize,i)%Y+dY, Org_loc(Org(i)%BodySize,i)%X+dX)
       ! the point and direction defines the region to be tested

       ! we are in this routine because move%possible is false ..
       ! set this temporarily to check if backward motion is possible

     Org(i)%Can%Move = .true.
     Call Particles_push(i, Point, Reversedirection)  ! clear away particles and test if move is possible
       If (Org(i)%Can%Move) then
         Org(i)%Can%Move = .false.                    ! return the flag to its original state

           ! find orientation of head by locating the "neck" and aligning the head to it
               If (Org_Loc(neck,i)%X .eq. Org_Loc(neck+(Org(I)%Width-1),i)%X) then ! it is vertically oriented facing 90 or 270
                 If (Org_Loc(neck,i)%Y .gt. Org_Loc(neck+(Org(I)%Width-1),i)%Y) then ! 270          
                   Call Head_rotate (i, 270)
				   RevDir = 90
				   PreDir = 4
                 Else
                   Call Head_rotate (i, 90)
				   RevDir = 270
				   PreDir = 2
                 End if
               Else if (Org_Loc(neck,i)%Y .eq. Org_Loc(neck+(Org(I)%Width-1),i)%Y) then  ! 0 or 180
                 If (Org_Loc(neck,i)%X .gt. Org_Loc(neck+(Org(I)%Width-1),i)%X) then ! 180
                   Call Head_rotate (i, 180)
				   RevDir = 0
				   PreDir = 3
                 Else
                   Call Head_rotate (i, 0)
				   RevDir = 180
				   PreDir = 1
                 End if
               End if

           DO j =1, Org(i)%BodySize            ! image the body
             Image(j) = Org_Loc(j,i)
           END Do
           DO j =Org(i)%BodySize, Org(i)%BodySize-(Org(I)%Width-1), -1 ! move the tip of the tail
             Org_Loc(j,i) = Coordinates(Image(j)%Y+dY, Image(j)%X+dX)
             Matrix(Org_loc(j,i)%Y, Org_loc(j,i)%X) = CellContent(j,i)       
             If (oxygen_ON) O2(Org_loc(j,i)%Y,Org_loc(j,i)%X) = O2(image(j)%Y,image(j)%X)   !! need be checked if working correctly
           END Do
           DO j = Org(i)%BodySize-Org(I)%Width, 1, -1     ! move the rest of the body
             Org_Loc(j,i) = Coordinates(Image(j+Org(I)%Width)%Y, Image(j+Org(I)%Width)%X)
             Matrix(Org_loc(j,i)%Y,Org_loc(j,i)%X) = cellContent(j,i)   
             If (oxygen_ON) O2(Org_loc(j,i)%Y,Org_loc(j,i)%X) = O2(image(j)%Y,image(j)%X)   !! need be checked if working correctly
           end do
           Do j = Org(I)%Width, 1, -1
             Matrix(Image(j)%Y, Image(j)%X) = CellContent(w,w)
             If (oxygen_ON) O2(image(j)%Y,image(j)%X) =               &
                       O2(image(Org(i)%BodySize- j + 1)%Y+dY,image(Org(i)%BodySize- j + 1)%X + dX)   !! need be checked if working correctly
           End do
		   
		   
       Dir_rec(:,i) = EOShift (Dir_rec(:,i), 1)
	   PreOri = Reversedirection - 180 
	   if (PreOri < 0) then 
	   PreOri = PreOri +360
	   end if 
	   
	   Dir_rec(Org(i)%length -  Org(i)%width,i) = PreOri

         Org(i)%Move%ReverseCount = Org(i)%Move%ReverseCount + 1
	   
	   select case(Reversedirection)
	     case(0)
		   Vb(2,i) = Vb(2,i) + org(i)%width
	     case(180)
		   Vb(2,i) = Vb(2,i) - org(i)%width
	     case(90)
		   Ub(2,i) = Ub(2,i) - org(i)%width
	     case(270)
		   Ub(2,i) = Ub(2,i) + org(i)%width
		 case default
		   print *, 'error in setting flow const---reversing',Reversedirection
		 end select
		 
           IF (Org(i)%Move%ReverseCount .gt. Org(I)%Width) then    ! reversed enough to be able to swap head and tail
             Org(i)%Move%ReverseCount = 0
             Org(i)%Is%Reversing = .false.
			   
			   do j = 1, Org(i)%BodySize
			   image_2(j) = Org_loc(j,i)
			   end do 
			   
			   do k = 1, Org(i)%length - Org(i)%width 
			   
			   RevDir = Dir_rec(k,i) - 180
			  if (RevDir < 0) RevDir = RevDir + 360 
			  
			  rec_image(Org(i)%length - Org(i)%width - k +1) = RevDir 
			   
			   call Head_rotate(i, RevDir)
			   
			   select case(RevDir)
			     case (0) 
				   dx = 0; dy = -1
				 case (90) 
				   dx = 1; dy = 0
				 case (180) 
				   dx = 0; dy = 1
				 case (270)
				   dx = -1; dy = 0
			  end select 
			  
			  do j = 1, Org(i)%BodySize
			  image(j) = Org_loc(j,i)
			  end do 
			  
			  do j = 1, Org(i)%HeadSize
			    Org_loc(j,i)%X = image(j)%x+dx 
			    Org_loc(j,i)%Y = image(j)%y+dy
				if (Org_loc(j,i)%X < 1) Org_loc(j,i)%X = Org_loc(j,i)%X + n_col
				if (Org_loc(j,i)%X > n_col) Org_loc(j,i)%X = Org_loc(j,i)%X - n_col
				if ((Org_loc(j,i)%Y > n_row).or.(Org_loc(j,i)%Y < 1)) print *,'error in head_move in Y direction'
              end do 	
			  
			  do j = Org(i)%HeadSize+ 1, Org(i)%HeadSize + Org(i)%width*k  !! width*(width + k)
			    Org_loc(j,i)%X = image(j-Org(i)%width)%x
			    Org_loc(j,i)%Y = image(j-Org(i)%width)%y
		      end do 
			  
			  end do
			  
			  do j = 1, Org(i)%length - Org(i)%width
			    Dir_rec(j,i) = rec_image(j)
		      end do
			  
			  Org(i)%Orientation = Reversedirection
			  
               DO j =1, Org(i)%BodySize      ! image the body
                 Image(j) = Org_Loc(j,i)
               END Do
               DO j =1, Org(i)%BodySize      ! do the swap of head and tail
                 Matrix(Org_loc(j,i)%Y, Org_loc(j,i)%X) = cellcontent(j,i)  
               END Do 			   
			   
			   U_tmp = Ub(:,i)
			   V_tmp = Vb(:,i)
			   Ub(2,i) = U_tmp(1)
			   Ub(1,i) = U_tmp(2)
			   Vb(2,i) = V_tmp(1)
			   Vb(1,i) = V_tmp(2)
			   
           End if
       Else ! stuck, try forward again
         Org(i)%Move%ReverseCount = 0
         Org(i)%Is%Reversing = .false.
       End if
       
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after reverse"; write(file_log,*)time, "err after reverse"; end if 
		 if ((any(Org_loc%X<0)) .or. (any(Org_loc%Y <0)) .or. (any(Org_loc%Y > N_row)) .or. (any(Org_loc%Y > N_row))) then
		   print *, 'errrrrrrr in org_loc'
		 end if 
         end if 
       
   End subroutine Organism_reverse

   ! ***************

   Subroutine Organism_dead(i)
     ! remove organism from matrix, releasing the particles in its guts and re-initialise at the water interface

   Use GlobalVariables
   Implicit None
   Logical :: IsWater
   Integer :: j, k
   Integer, Intent(IN) :: i
   Real :: URand
   
   integer :: xx, yy
   
     Do j = 1, Org(i)%BodySize
       Matrix(Org_loc(j,i)%Y,Org_loc(j,i)%X) = CellContent(w,w)
     End do

     Do j = 1, Org(i)%Gut%Capacity
       If (Guts(j,i)%Class .eq. p) then
         do
             Call Random_Number(URand)
             k = URand * Org(i)%BodySize + 1
             If (IsWater(Org_loc(k,i))) then
                 MAtrix(Org_loc(k,i)%Y, Org_loc(k,i)%X) = Guts(j,i)
				 particle(Guts(j,i)%value)%loc = Org_loc(k,i)   !! YK 
                 if (oxygen_ON) O2(Org_loc(k,i)%Y, Org_loc(k,i)%X)%oxygen = 0.0
                 Guts(j,i) = CellContent(w,w)
                 Exit
             End if
         End do
       End if
     End do

     Call Matrix_populate(i)
         
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after dead"; write(file_log,*)time, "err after dead"; end if 
		 if ((any(Org_loc%X<0)) .or. (any(Org_loc%Y <0)) .or. (any(Org_loc%Y > N_row)) .or. (any(Org_loc%Y > N_row))) then
		   print *, 'errrrrrrr in org_loc'
		 end if 
         end if 
   End subroutine Organism_dead

   ! *****************************************************************

   Subroutine Organism_wrapAround(i)
     ! low level routine that translates the organism i when they are wrapping around

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Integer :: j, k, dx, dy, X_new
   Integer, Intent(IN) :: i
   Type(Coordinates) :: Image(Org(i)%BodySize)
   
   integer :: xx, yy
   logical :: IsWater

   dy = 0

   If (Allocated (Tail)) Deallocate (tail)
   Allocate (Tail(Org(I)%Width))

   Select Case (Org(i)%Orientation)  ! (dX =+1 is right;  -1 is left); (dY =+1  is down;-1 is up)
     Case(90)
     dX = 1; X_new = 1
     Case(270)
     dX = -1; X_new = N_col
   End Select

   ! take an initial image of the body
     DO j =1, Org(i)%BodySize
       Image(j) = Org_Loc(j,i)
     END Do

       Select Case (Org(i)%Is%WrappingAround)
         Case (.false.)
           k = Org(i)%HeadSize
             ! move the tip of the head to the other side
                 DO j =1, Org(I)%Width
                   Org_Loc(j,i) = Coordinates(Image(j)%Y+dY, X_new)
                   Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = CellContent(j,i)
                   if (oxygen_ON) O2(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = O2(image(j)%Y, image(j)%X)
                 END Do
               ! move the rest of the head
                 DO j =Org(I)%Width+1, k
                   Org_Loc(j,i) = Coordinates(Image(j-Org(I)%Width)%Y, Image(j-Org(I)%Width)%X)
                   Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = CellContent(j,i)
                   if (oxygen_ON) O2(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = O2(image(j)%Y, image(j)%X)
                 END Do

         Case (.true.)
           k = Org(i)%Width
             ! move the head
               DO j =1, k                                    ! move the head
                 Org_Loc(j,i) = Coordinates(Image(j)%Y+dy, Image(j)%X+dX)
                 Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = CellContent(j,i)
                 if (oxygen_ON) O2(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = O2(image(j)%Y, image(j)%X)
               END Do
         End Select

   ! complete the rest of the body
       DO j = k+1, Org(i)%BodySize
         Org_Loc(j,i) = Coordinates(Image(j-Org(I)%Width)%Y, Image(j-Org(I)%Width)%X)
         Matrix(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = CellContent(j,i)
         if (oxygen_ON) O2(Org_Loc(j,i)%Y, Org_Loc(j,i)%X) = O2(image(j)%Y, image(j)%X)
       end do
   ! fix the tail
       Do j = 1, Org(I)%Width
         Tail(j) = Image(Org(i)%BodySize-Org(I)%Width+j)
         Matrix(Image(Org(i)%BodySize-Org(I)%Width+j)%Y, Image(Org(i)%BodySize-Org(I)%Width+j)%X) = CellContent(w,w)
         if (oxygen_ON) O2(tail(j)%Y, tail(j)%X) = O2(image(j)%Y, image(j)%X)
       End do

     MovementHistory(i,1) = 1
     Org(i)%Move%Bearing_Distance = Org(i)%Move%Bearing_Distance + 1
	 
	 
       Dir_rec(:,i) = EOShift (Dir_rec(:,i), -1)
	   Dir_rec(1,i) = Org(i)%Orientation
	   
	   
	   select case(Org(i)%Orientation)
	     case(0)
		   Vb(1,i) = Vb(1,i) + org(i)%width
	     case(180)
		   Vb(1,i) = Vb(1,i) - org(i)%width
	     case(90)
		   Ub(1,i) = Ub(1,i) - org(i)%width
	     case(270)
		   Ub(1,i) = Ub(1,i) + org(i)%width
		 case default
		   print *, 'error in setting flow const---wrappingArround',Org(i)%Orientation
		 end select
     
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after wrap"; write(file_log,*)time, "err after wrap"; end if 
		 if ((any(Org_loc%X<0)) .or. (any(Org_loc%Y <0)) .or. (any(Org_loc%Y > N_row)) .or. (any(Org_loc%Y > N_row))) then
		   print *, 'errrrrrrr in org_loc'
		 end if 
         end if 
     
   End Subroutine Organism_wrapAround


   ! *************************************

   Subroutine Matrix_fill()
     ! fill the matrix with particles

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Integer ::  X, Y, k, Particle_ID, Lability
   Real :: Porosity, DEpth, URand
   real :: lab_real
   
   ! initialise the particles and matrix
       Matrix = CellContent(w,w)  ! fill with water
       Particle%Loc = Coordinates(0,0)
       PArticle%Plane = 0
       Particle%Lability = 0
       Particle%Lability_Time0 = 0
       Particle%Pb%Activity0 = 0
       Particle%Pb%Activity = 0
       Particle%Pb%Time0 = 0
       
       Particle%OM%OMact = 0
       Particle%OM%OMact_0 = 0
       Particle%OM%OM_Time0 = 0

   ! fill the matrix :: based upon user input porosity, activity, lability, etc,  gradients/functional forms
       Particle_Id = 0
         DO X = 1, N_Col
         DO Y = N_RowWater + 1, N_Row
           Depth = Real(Y - N_RowWater-1) * PixelSize
           Call Random_Number(URand)
             IF (URand .gt. Porosity(Y)) THEN    ! first porosity
               Particle_ID = Particle_ID + 1
               Particle(Particle_ID)%Loc = Coordinates(Y,X)
               Particle(Particle_ID)%Loc_Init = Coordinates(Y,X)
               Particle(Particle_ID)%Lability = int(Lab_real(Y,X))
               Particle(Particle_ID)%Pb%Activity = Exp( - DecayConstant * Depth / SedimentationRate)
                 ! the activity of Pb210 assuming constant sed rate
               Particle(Particle_ID)%Pb%Activity0 = Particle(Particle_ID)%Pb%Activity
               Matrix(Y,X) = CellContent(Particle_ID, p)
               
               Particle(Particle_ID)%OM%OMact = real(Lab_real(Y,X))
               Particle(Particle_ID)%OM%OMact_0 = Particle(Particle_ID)%OM%OMact
               !!!  YK  treat lability as continuous value, not discrete ones
               
             END IF
         END Do
         END Do

     Total_N_Particles = Particle_ID
     Total_N_Particles0 = Total_N_Particles
       ! the initial number of particles that constrains the total allowable number in the simulation

   ! tolerance to deviations from the intial number of particles ..
       Tolerance = 0.5 ! tolerate n% deviations from the initial condition (total number of particles)
       ParticleTolerance = Int(Tolerance * Real(Total_N_Particles0) / 100.)

     Allocate (Particle_ID_free(2*Total_N_Particles0))

     PArticle_ID_Free = 0
     PointerToFreeParticle = 0
       Do k = Total_N_Particles+1, 2*Total_N_Particles  ! create the pointers to available particle numbers
         PointerToFreeParticle = PointerToFreeParticle + 1
         Particle_ID_free(PointerToFreeParticle) = k
       End do

   End Subroutine Matrix_fill

   !****************************************************************

   Subroutine Matrix_populate(i)
     ! populate the matrix with organisms

   Use GlobalVariables
   Implicit None
   Logical :: IsMoveable, InMatrix, Finished, IsParticle
   Integer :: i, j, k, X, Y, Orientation
   Real    :: URand(3)

   integer :: xx, yy

     j = 0
     Finished = .false.

   ! choose a co-ordinate/orientation; test for overlap with other organisms or matrix edge
       Do While (.not. Finished)

         Call Random_Number(URand)

     ! start in the water column
         Org_Loc(1,i)%X = INT(URand(1)*(N_col-2*Buffer_Zone)+1) + Buffer_Zone   !!   Buffer_zone is equivalent to N_col in default setting
         Org_Loc(1,i)%Y = INT(URand(2)*(2 * Buffer_Zone) + 1) + N_RowWater - 2*Buffer_Zone


         Orientation = Int(URand(3) * 4)
         Finished = .true.

         if (Orientation .le. 2) then
           Orientation = 1 ! force horizontal
         else
           Orientation = 3
         End if

         Select Case (Orientation)
         Case(0)
           Org(i)%Orientation = 0
main0:    DO Y = Org_Loc(1,i)%Y, Org_Loc(1,i)%Y+(Org(i)%Length-1), 1
             DO X = Org_Loc(1,i)%X, Org_Loc(1,i)%X+(Org(i)%Width-1), 1
                 IF (.not. InMatrix(CoOrdinates(Y,X))) then
                   Finished = .false.
                   Exit main0
                 Else
                   if (.not. IsMoveable(CoOrdinates(Y,X))) Then
                   Finished = .false.
                   End if
                 END IF
             END Do
           END Do main0

         Case(1)
           Org(i)%Orientation = 90
main90:   DO X = Org_Loc(1,i)%X, Org_Loc(1,i)%X-(Org(i)%Length-1), -1
             DO Y = Org_Loc(1,i)%Y, Org_Loc(1,i)%Y+(Org(i)%Width-1), 1
                 IF (.not. InMatrix(CoOrdinates(Y,X))) then
                   Finished = .false.
                   Exit main90
                 Else
                   if (.not. IsMoveable(CoOrdinates(Y,X))) Then
                   Finished = .false.
                   End if
                 END IF
             END Do
           END Do main90

           Case(2)
             Org(i)%Orientation = 180
main180:    DO Y = Org_Loc(1,i)%Y, Org_Loc(1,i)%Y - (Org(i)%Length-1), -1
               DO X = Org_Loc(1,i)%X, Org_Loc(1,i)%X - (Org(i)%Width-1), -1
                   IF (.not. InMatrix(CoOrdinates(Y,X))) then
                     Finished = .false.
                     Exit main180
                   Else
                     if (.not. IsMoveable(CoOrdinates(Y,X))) Then
                     Finished = .false.
                     End if
                   END IF
               END Do
             END Do main180

           Case(3)
             Org(i)%Orientation = 270
main270:    DO X = Org_Loc(1,i)%X, Org_Loc(1,i)%X+(Org(i)%Length-1), 1
               DO Y = Org_Loc(1,i)%Y, Org_Loc(1,i)%Y-(Org(i)%Width-1), -1
                   IF (.not. InMatrix(CoOrdinates(Y,X))) then
                     Finished = .false.
                     Exit main270
                   Else
                     if (.not. IsMoveable(CoOrdinates(Y,X))) Then
                     Finished = .false.
                     End if
                   END IF
               END Do
             END Do main270

           End Select
       End do

   ! now populate the matrix, overwritten particles are placed into the gut
    k = 0
    Select Case (Orientation)
     Case(0)
       DO Y = Org_Loc(1,i)%Y, Org_Loc(1,i)%Y+(Org(i)%Length-1), 1
       DO X = Org_Loc(1,i)%X, Org_Loc(1,i)%X+(Org(i)%Width-1), 1
       if (IsParticle(CoOrdinates(Y,X))) then
         if (k .lt. Org(i)%Gut%Capacity) then
           k = k + 1
           Guts(k,i) = Matrix(Y,X)
         Else
           Call Particle_remove(CoOrdinates(Y,X))
         End if
       End if
         j = j + 1
         Org_Loc(j,i) = Coordinates(Y,X)
         Matrix(Y,X)  = CellContent(j,i)
       END Do
       END Do

     Case(1)
       DO X = Org_Loc(1,i)%X, Org_Loc(1,i)%X-(Org(i)%Length-1), -1
       DO Y = Org_Loc(1,i)%Y, Org_Loc(1,i)%Y+(Org(i)%Width-1), 1
       if (IsParticle(CoOrdinates(Y,X))) then
         if (k .lt. Org(i)%Gut%Capacity) then
           k = k + 1
           Guts(k,i) = Matrix(Y,X)
         Else
           Call Particle_remove(CoOrdinates(Y,X))
         End if
       End if
         j = j + 1
         Org_Loc(j,i) = Coordinates(Y,X)
         Matrix(Y,X)  = CellContent(j,i)
       END Do
       END Do

     Case(2)
       DO Y = Org_Loc(1,i)%Y, Org_Loc(1,i)%Y - (Org(i)%Length-1), -1
       DO X = Org_Loc(1,i)%X, Org_Loc(1,i)%X - (Org(i)%Width-1), -1
       if (IsParticle(CoOrdinates(Y,X))) then
         if (k .lt. Org(i)%Gut%Capacity) then
           k = k + 1
           Guts(k,i) = Matrix(Y,X)
         Else
           Call Particle_remove(CoOrdinates(Y,X))
         End if
       End if
         j = j + 1
         Org_Loc(j,i) = Coordinates(Y,X)
         Matrix(Y,X)  = CellContent(j,i)
       END Do
       END Do

     Case(3)
       DO X = Org_Loc(1,i)%X, Org_Loc(1,i)%X+(Org(i)%Length-1), 1
       DO Y = Org_Loc(1,i)%Y, Org_Loc(1,i)%Y-(Org(i)%Width-1), -1
       if (IsParticle(CoOrdinates(Y,X))) then
         if (k .lt. Org(i)%Gut%Capacity) then
           k = k + 1
           Guts(k,i) = Matrix(Y,X)
         Else
           Call Particle_remove(CoOrdinates(Y,X))
         End if
       End if
         j = j + 1
         Org_Loc(j,i) = Coordinates(Y,X)
         Matrix(Y,X)  = CellContent(j,i)
       END Do
       END Do

   End Select
   
     
	   Dir_rec(1:Org(i)%length - Org(i)%width,i) = Org(i)%Orientation

         if (errChk) then 
         if (time > 0 ) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after populate"; write(file_log,*)time, "err after populate"; end if 
         end if 
         end if 
 
   End Subroutine Matrix_populate

   ! ***************************************

   Subroutine Particle_add(Point)
     ! add a particle at the co-ordinates (Point%Y, Point%X)

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Integer :: Particle_ID, Lability
   Type(Coordinates) :: Point
   real :: lab_real

   integer :: xx, yy
   
     Particle_ID                             = Particle_ID_free(PointerToFreeParticle)
     Particle_ID_free(PointerToFreeParticle) = 0
     PointerToFreeParticle                   = PointerToFreeParticle - 1
     Total_N_Particles                       = Total_N_Particles + 1
     Matrix(Point%Y,Point%X)                 = CellContent(PARTICLE_ID,P)

       Particle(Particle_ID)%Plane           = 0
       Particle(Particle_ID)%Loc             = Coordinates(Point%Y, Point%X)
       Particle(Particle_ID)%Loc_Init        = Coordinates(Point%Y, Point%X)
       Particle(Particle_ID)%Lability        = int(Lab_real(Point%Y, Point%X))
       Particle(Particle_ID)%Lability_Time0  = Time
       Particle(Particle_ID)%Pb%Activity     = 1
       Particle(Particle_ID)%Pb%Activity0    = 1
       Particle(Particle_ID)%Pb%Time0        = Time
       
       
       Particle(Particle_ID)%OM%OMact        = real(Lab_real(Point%Y, Point%X))
       Particle(Particle_ID)%OM%OMact_0      = Particle(Particle_ID)%OM%OMact
       Particle(Particle_ID)%OM%OM_Time0     = Time
       
       O2(point%Y, point%X)%oxygen           = 0.0 
       O2(point%Y, point%X)%oxygen_use       = 0.0 
       O2(point%Y, point%X)%oxygen_pre       = 0.0 
       O2(point%Y, point%X)%value_pre        = 0
       O2(point%Y, point%X)%mark             = 0
       
   End Subroutine Particle_add

   ! **************

   Subroutine Particle_remove(Point)
     ! remove the particle at the co-ordinates (Point%Y, Point%X)

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Integer :: Particle_ID
   Type(Coordinates) :: Point

   integer :: xx, yy
   logical :: IsWater
   
     Particle_ID                             = Matrix(Point%Y,Point%X)%Value
     PointerToFreeParticle                   = PointerToFreeParticle + 1
       ! add one more to the table of unused particles
     Particle_ID_free(PointerToFreeParticle) = Particle_ID
     Total_N_Particles                       = Total_N_Particles - 1
     Matrix(Point%Y,Point%X)                 = CellContent(0,0)

       Particle(Particle_ID)%Plane           = 0
       Particle(Particle_ID)%Loc             = Coordinates(0,0)
       Particle(Particle_ID)%Loc_Init        = Coordinates(0,0)
       Particle(Particle_ID)%Lability        = 0
       Particle(Particle_ID)%Lability_Time0  = 0
       Particle(Particle_ID)%Pb%Activity     = 0
       Particle(Particle_ID)%Pb%Activity0    = 0
       Particle(Particle_ID)%Pb%Time0        = 0
       
       
       Particle(Particle_ID)%OM%OMact        = 0
       Particle(Particle_ID)%OM%OMact_0      = 0
       Particle(Particle_ID)%OM%OM_Time0     = 0
   
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after add"; write(file_log,*)time, "err after add"; end if 
		 end if 
		 
   End Subroutine Particle_remove

   ! **********************

   Subroutine Particle_sediment(ParticlesToSediment)
     ! controlling loop to input new particles to the water column

   Use GlobalVariables
   Implicit None
   Logical :: IsWater
   Integer :: m, ParticlesToSediment
   Real    :: URand(ParticlesToSediment)
   Type(Coordinates) :: Point

     Call Random_Number (URand)
       DO m = 1, ParticlesToSediment   ! the number of particle to sediment in this time frame
           Point = CoOrdinates(1, (Int(URand(m)*Real(N_Col))+1))   ! randomly locate the point at the top
             if (IsWater(Point)) then
               Call Particle_add(Point)
               Call SedimentDown(Point)
               RainCount = RainCount + 1
             End if
       END Do
       
   End Subroutine Particle_sediment

! ******************************************************************

   Subroutine SedimentDown(Point)
     ! low level routine - move particle downwards/diagonally until the local density is stable

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Logical :: InMAtrix, IsWater
   Integer ::  TestWindow
   Real :: Porosity_local, URand(3)
   Type(Coordinates) :: Point, PointBelow, PointBelowBelow, PointBeside
   Type(oxygen_conc) :: O2_buff  !! YK

   integer :: xx, yy
   
   ! Initialisations:
       Call Random_Number (URand)
       TestWindow = 3 * Urand(3) + 1      ! i.e. a box of 3X3 cells wide

   ! Main loop
       Do

         PointBelow = CoOrdinates(Point%Y+1,Point%X)
         PointBelowBelow = CoOrdinates(PointBelow%Y+1,PointBelow%X)

         If (.not. IsWater(PointBelow)) then
           Exit                                      ! i.e., no need to move down
         Else                                        ! check if water below that point ..

           If (IsWater(PointBelowBelow)) then        ! all clear, move down
             O2_buff = O2(PointBelow%Y,PointBelow%X) 
             Matrix(PointBelow%Y,PointBelow%X) = Matrix(Point%Y,Point%X)
             if (oxygen_ON) O2(PointBelow%Y,PointBelow%X) = O2(Point%Y,Point%X)
             Particle(Matrix(PointBelow%Y,PointBelow%X)%Value)%loc = PointBelow
             
 
Particle(Matrix(PointBelow%Y,PointBelow%X)%Value)%loc_INIT = PointBelow
             Matrix(Point%Y,Point%X) = CellContent(w,w)
             if (oxygen_ON) O2(Point%Y,Point%X) = O2_buff
             Point = PointBelow
             Cycle                                   ! do the next particle
           Else                                      ! else, not all clear, check densities and move down
                                                     !       to fill last space if densities match ..
             If (Porosity_local(Point, TestWindow) .lt. Porosity_local(PointBelowBelow, TestWindow)) then
               O2_buff = O2(PointBelow%Y,PointBelow%X)  
               ! move down and fill the space
               Matrix(PointBelow%Y,PointBelow%X) = Matrix(Point%Y,Point%X)
               if (oxygen_ON) O2(PointBelow%Y,PointBelow%X) = O2(Point%Y,Point%X)
               Particle(Matrix(PointBelow%Y,PointBelow%X)%Value)%loc = PointBelow
 
Particle(Matrix(PointBELOW%Y,PointBeLOW%X)%Value)%loc_Init = PointBelow
                 ! THIS IS THE FINAL SETTLING POINT
               Matrix(Point%Y,Point%X) = CellContent(w,w)
               if (oxygen_ON) O2(Point%Y,Point%X) = O2_buff
               Point = PointBelow
             End if
             Exit
           End if
         End if
       End do

     If (URand(1) .gt. 0.3) then                     ! in 30% cases do not shift particles to the side
       If (URand(2) .gt. 0.5) then
         PointBeside = CoOrdinates(Point%Y+1, Point%X+1)
       Else
         PointBeside = CoOrdinates(Point%Y+1, Point%X-1)
       End if
         If (InMatrix(PointBeside)) then
           If (IsWater(PointBeside)) then
             O2_buff = O2(PointBeside%Y,PointBeside%X) 
             Matrix(PointBeside%Y,PointBeside%X) = Matrix(Point%Y,Point%X)              
             if (oxygen_ON) O2(PointBeside%Y,PointBeside%X) = O2(Point%Y,Point%X)
             Matrix(Point%Y,Point%X) = CellContent(w,w)
             if (oxygen_ON) O2(Point%Y,Point%X) = O2_buff
             Particle(Matrix(PointBeside%Y,PointBeside%X)%Value)%loc = PointBeside
 
Particle(Matrix(PointBeside%Y,PointBeside%X)%Value)%loc_Init = PointBeside
               ! THIS IS THE FINAL SETTLING POINT
           End if
         End if
     End if

         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after sedown"; write(file_log,*)time, "err after sedown"; end if 
		 end if 
     
   End Subroutine SedimentDown

   ! *******************************

   Subroutine WaterColumn_scan(DepthToScan)
     ! scan the water column and sediment down any particle found .. search intensity decreases with depth

   Use GlobalVariables
   Implicit none
   Integer :: Y,X, DepthToScan
   Logical :: IsPArticle
   Real    :: Pr
   Real    :: URand(DepthToScan)

   Call Random_Number(URand)

   Do Y = 1, DepthToScan
     Pr = Real(DepthToScan-Y+1) / Real(DepthToScan) + 0.25
       ! start with a 25 % chance at the interface and increase ..

     Do X = 1, N_Col
       If (URand(Y) .lt. Pr) then
         If (IsParticle(CoOrdinates(Y,X))) Call SedimentDown(CoOrdinates(Y,X))
       End if
     End do

   End do
   End Subroutine WaterColumn_scan

   ! ****************************

   Subroutine Matrix_constrain()
     ! constrain the total number of particles in the simulation to be ~constant within a given tolerance

   Use GlobalVariables
   use ieee_arithmetic
   Implicit none
   Logical :: IsOrganism, IsParticle
   Integer :: Y_new, Y, X, DummyVariable, PArticle_ID
   real :: Lab_real_New

   integer :: xx, yy
   integer :: i

   ! check bottom boundary to see if the shift is possible
       Y = N_Row

       Do X = 1, N_Col
         If (IsOrganism(CoOrdinates(Y,X))) Return  ! SHIFT NOT POSSIBLE
       End do

   ! As there are no boundary problems .. do the shift
       Do X = 1, N_Col
         if (IsParticle(CoOrdinates(Y,X))) Call Particle_remove(CoOrdinates(Y,X))
       End do

       Do X = 1, N_Col
         Do Y = N_Row - 1, 1, -1
           Y_new = Y + 1
             Select Case (Matrix(Y,X)%Class)
               Case(w)
                 Matrix(Y_new,X) = CellContent(w,w)
                 if (oxygen_ON) O2(Y_new, X) = O2(Y, X)
               Case(p)
                 Particle_ID = Matrix(Y,X)%Value
                 Particle(Particle_ID)%loc = Coordinates(Y_new, X)
                 Particle(Particle_ID)%loc_Init%Y = Particle(Particle_ID)%loc_Init%Y + 1 ! adjust initial position for shifting down of matrix
                 Matrix(Y_new,X) = Matrix(Y,X)
                 if (oxygen_ON) O2(Y_new, X) = O2(Y, X)
                 DummyVariable = Lab_real_New(Particle_ID)
               Case(1:)    ! organism
                 Org_Loc(Matrix(Y,X)%Value, Matrix(Y,X)%Class) = Coordinates(Y_new, X)
                 Matrix(Y_new,X) = Matrix(Y,X)
                 if (oxygen_ON) O2(Y_new, X) = O2(Y, X)
             End Select
           End do
       End do

       Y = 1
         Do X = 1, N_Col
           Matrix(Y,X) = CellContent(w,w)
         End do          

         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after constrain(sub)"; write(file_log,*)time, "err after constrain(sub)"; end if 
         end if 
		 
   End Subroutine Matrix_constrain

   ! ***********************************************

   SUBROUTINE Particles_push(i, Point, Orientation)
     ! push particles away in the block defined by the coordinates (Point) and Orientation

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Logical :: Finished, DoIt, LoopIt, Stuck
   Integer :: i, k, m, o, X, Y, Left, Forward, Right, Orientation
   Integer :: DirPossible(3), StuckCount
   Integer :: N_Particle, N_Organism, N_Water
   Integer :: X_start, X_end, X2_start, X2_end, dX, dY, Y_start, Y_end
   Real  :: DirRandom(3), DirWeights(3), URand(3)
   Type(CellContent) :: ToMove
   Type(Coordinates) :: ToMove_loc, Point
   Type(oxygen_conc) :: O2_buff   !!  YK
   integer :: X_buff, Y_buff

   integer :: xx, yy
   logical :: IsWater
   type(cellcontent) :: blcimg(org(i)%width)
   type(coordinates) :: blcimg_loc(org(i)%width)
   
   integer :: xxg, xxp

   Stuck = .false.
   StuckCount = 0
   m = 0
   blcimg_loc = coordinates(0,0)

Main: Do While (m .lt. Org(i)%Width)

       m = m + 1

       Call Block_read(i, Point, Orientation, N_Particle, N_Organism, N_Water)
         If (Block_outsideMatrix) Exit Main
         If (N_Organism .gt. 0) Exit Main
         if (N_Particle .eq. 0) Return       ! finished
		 
	   if (m == 1) then   ! YK make a copy of initial block
	      blcimg = block
		  do o = 1, org(i)%width
		    if (block(o)%class == p) then 
			   blcimg_loc(o)%x = particle(block(o)%value)%loc%x
			   blcimg_loc(o)%y = particle(block(o)%value)%loc%y
			end if 
		  end do 
	   end if 

       Call Random_Number(URand)   ! pick a particle
       if (URand(1) .ge. 0.5) then
         Do o = Org(i)%Width, 1, -1
           if (Block(o)%Class .eq. p) then
             k = o
             Exit
           End if
         End Do
       else
         Do o = 1, Org(i)%Width, 1
           if (Block(o)%Class .eq. p) then
             k = o
             Exit
           End if
         End Do
       End if

     ToMove      = Block(k)
     ToMove_loc  = Block_loc(k)

     Left = Orientation-90;  Forward = Orientation;  Right = Orientation+90
     DirPossible = (/Left, Forward, Right/)
     DirWeights  = (/1., 0.5, 1./) ! UnEQUAL WEIGHTING in the forward orientation
     Finished = .false.
     
ChooseDir:  Do While (Sum(DirWeights) .ne. 0) ! else return to the main do loop

     Call Random_Number(URand)
     DirRandom = URand * DirWeights
     o = (MaxLoc(DirRandom, Dim=1))
     DoIT = .false.
     LoopIt = .false.

   Select Case (DirPossible(o))

     Case(0, 360)
       X = ToMove_loc%X
       Y_start = ToMove_loc%Y-1; Y_end = 1;  dY = -1
         DO  Y = Y_start, Y_end, dY
           If      (Matrix(Y,X)%Class .eq. p) then
             Cycle
           Else if (Matrix(Y,X)%Class .eq. w) then       ! i.e. it is possible to compress .. so do it
             DoIt = .true.
             Exit
           Else
              DirWeights(o) = 0
              Cycle ChooseDir
           End if
         End DO
         If (DoIt) then
            DO Y = Y_start, Y_end, dY
             Call Particles_swap(X, Y, Finished, ToMove)
             If (Finished) then
               O2_buff = O2(Y_start+1, X)
               Matrix(ToMove_loc%Y, ToMove_loc%X) = CellContent(w,w)
               if (oxygen_ON) then
                  O2(Y_start+1, X) = O2(Y, X)
                  
                  O2(Y, X) = O2_buff 
               end if 
			   
               Finished = .false.
               Exit Choosedir
             End if
            End do
         End if
         DirWeights(o) = 0

     Case(180)
       X = ToMove_loc%X
       Y_start = ToMove_loc%Y+1; Y_end = N_Row; dY = 1

         DO  Y = Y_start, Y_end, dY
           If (Matrix(Y,X)%Class .eq. p) then
             If (Y .eq. N_row) then

               If (.not. LocalMixing) then               ! this turns off the pushing of particles out of the bottom boundary
                 Call PArticle_remove(CoOrdinates(Y,X))
                 DoIt = .true.
               End if

               Exit
             Else
               Cycle
             End if
           Else if (Matrix(Y,X)%Class .eq. w) then       ! i.e. it is possible to compress .. so do it
             DoIt = .true.
             Exit
           Else
              DirWeights(o) = 0
              Cycle ChooseDir
           End if
         End DO
         If (DoIt) then
            DO Y = Y_start, Y_end, dY 
             Call Particles_swap(X, Y, Finished, ToMove)
             If (Finished) then
               O2_buff = O2(Y_start -1, X)
               Matrix(ToMove_loc%Y, ToMove_loc%X) = CellContent(w,w)
               if (oxygen_ON) then 
                 O2(Y_start -1,X) = O2(Y,x)
                  O2(Y,x) = O2_buff
               end if 
			   
               Finished = .false.
               Exit Choosedir
             End if
            End do
         End if
         DirWeights(o) = 0

     Case(-90, 90, 270)
     Y = ToMove_loc%Y

       Select Case (DirPossible(o))
         Case (90)
           X_start = ToMove_loc%X+1; X_end = N_col; dX = 1
           X2_start = 1; X2_end = Point%X
         Case (-90, 270)
           X_start = ToMove_loc%X-1; X_end = 1; dX = -1
           X2_start = N_col; X2_end   = Point%X
       End Select
       
         X_buff = ToMove_loc%X

         DO  X = X_start, X_end, dX
           If      (Matrix(Y,X)%Class .eq. p) then
             Cycle
           Else if (Matrix(Y,X)%Class .eq. w) then
             DoIt = .true.
             Exit
           Else
             DirWeights(o) = 0
             Cycle ChooseDir
           End if
         End DO          ! finished to the edge .. loop around ..
           If (DoIt) then
             DO  X = X_start, X_end, dX
               Call Particles_swap(X, Y, Finished, ToMove)
               If (Finished) then
                 O2_buff = O2(Y, X_buff) 
                 Matrix(ToMove_loc%Y, ToMove_loc%X) = CellContent(w,w)
                 if (oxygen_ON) then 
                   O2(Y, X_buff) = O2(Y, X)
                   
                   O2(Y, X) = O2_buff
                 end if
                 
                 Finished = .false.
                 Exit Choosedir
               End if
             End do
           END IF

         DO  X = X2_start, X2_end, dX
           If      (Matrix(Y,X)%Class .eq. p) then
             Cycle
           Else if (Matrix(Y,X)%Class .eq. w) then
             LoopIt = .True.
             Exit
           Else
             DirWeights(o) = 0
             Cycle ChooseDir
           End if
         End DO
           If (loopIt) then
             DO  X = X_start, X_end, dX
               Call Particles_swap(X, Y, Finished, ToMove)
             End do

             Particle(ToMove%Value)%Plane = Particle(ToMove%Value)%Plane + dX

             DO  X = X2_start, X2_end, dX
               Call Particles_swap(X, Y, Finished, ToMove)
               If (Finished) then
                 O2_buff = O2(Y, X_buff)
                 Matrix(ToMove_loc%Y, ToMove_loc%X) = CellContent(w,w)
                 if (oxygen_ON) then 
                   O2(Y, X_buff) = O2(Y, X)
                   if (ieee_is_nan(O2(Y, X2_start)%oxygen)) then
                     print *,"NAN in oxygen",Y, X2_start
                     write(File_log, *) Time, "NAN in oxygen",Y, X2_start
                   endif 
                   if (ieee_is_nan(O2(Y, X2_start)%oxygen_use)) then
                     print *,"NAN in oxygen_use",Y, X2_start
                     write(File_log, *) Time, "NAN in oxygen",Y, X2_start
                   end if 
                   
                   O2(Y, X) = O2_buff
                   if (ieee_is_nan(O2(Y, X)%oxygen)) then 
                     print *,"NAN in oxygen",Y, X
                     write(File_log, *) Time, "NAN in oxygen",Y, X
                   end if 
                   if (ieee_is_nan(O2(Y, X)%oxygen_use)) then
                     print *,"NAN in oxygen_use",Y, X
                     write(File_log, *) Time, "NAN in oxygen",Y, X
                   end if 
                 end if 
                 
                 Finished = .false.
                 Exit Choosedir
               End if
             End do
           End if
         DirWeights(o) = 0
     End Select

   End Do ChooseDir

   End do Main

   Org(i)%Can%Move = .false.
   
   
   do o = 1, org(i)%width
     if (blcimg(o)%class/=p) cycle
     if (particle(blcimg(o)%value)%loc%x /= blcimg_loc(o)%x) then 
		if (particle(blcimg(o)%value)%loc%x  < blcimg_loc(o)%x) then 
		   Ub(1,i) = Ub(1,i) + 2.5
		else if (particle(blcimg(o)%value)%loc%x  > blcimg_loc(o)%x) then 
		   Ub(1,i) = Ub(1,i) - 2.5
		end if 
    end if 

    if (particle(blcimg(o)%value)%loc%y /= blcimg_loc(o)%y) then 
        if (particle(blcimg(o)%value)%loc%y  < blcimg_loc(o)%y) then 
		   Vb(1,i) = Vb(1,i) + 2.5
        else if (particle(blcimg(o)%value)%loc%y  > blcimg_loc(o)%y) then 
		   Vb(1,i) = Vb(1,i) - 2.5
		end if 
	end if 
   end do 
 
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after push"; write(file_log,*)time, "err after push"; end if 
         end if 

   End Subroutine Particles_push

   ! **************************************************************************

   SUBROUTINE Particles_swap(X, Y, Finished, ToMove)
     ! low level routine to swap particles used by the pushing routine

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Logical, Intent(Out) :: Finished
   Integer, Intent(in) :: X, Y
   Type(CellContent), Intent(InOut) :: ToMove
   Type(CellContent)  :: Buffer
   
   integer :: xx, yy
   logical :: IsWater

   Select Case (Matrix(Y,X)%Class)
     Case (p)
       Buffer = Matrix(Y,X)
       Matrix(Y,X) = ToMove
       Particle(Matrix(Y,X)%Value)%Loc = Coordinates(Y,X)   ! YK assign loc of particle id in tomove to x,y 
       ToMove = Buffer        !!  YK   tomove is replaced by Matrix(Y,X)
       Finished = .false.
     Case (w)
       Buffer = Matrix(Y,X)  !  YK added  buffer is now water 
       Matrix(Y,X) = ToMove   !  YK matrix(Y,X) is now particle(ToMove)
       Particle(Matrix(Y,X)%Value)%Loc = Coordinates(Y,X)      !!  particle location is moved to Y,X
       ToMove = Buffer    !  YK added          ToMove is now water   
       Finished = .true.
     Case Default
       Write (*,*) ' We have a problem in Particles_swap'
       write(file_log,*) " We have a problem in Particles_swap"
       Finished = .false.
   End Select
   
   
   End subroutine Particles_swap

   ! **************************************************************************

   Subroutine Particles_ingest(i, Point, Orientation)
     ! organism (i) ingests particles for the block defined by co-ordinates (Point%Y, Point%X) and the orientation input

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Logical :: EAt
   Integer :: i, k, m, n, o, Orientation
   Integer :: N_Particle, N_Organism, N_Water
   Real    :: URand(Org(i)%Width)
   Type(Coordinates) :: Point
   real :: n_real
   Type(oxygen_conc) :: O2_buff   !!  YK

   integer ::  xx, yy
   logical :: IsWater

   If (.not. Org(i)%Can%Ingest) then; Return; End if

   Call Block_read(i, Point, Orientation, N_Particle, N_Organism, N_Water)
     if (N_Particle .eq. 0) Return ! nothing to do

     if (N_Organism .gt. 0) then
       ! this is necessary as this routine is used in thw wraparound routine without any other tests.
       Org(i)%Can%Move = .false.
       Return
     End if

   If (LocalMixing) Return   ! if forcing local mixing, turn off the ingestion routine

   CALL Organism_peristalsis(i)
			   
   Call Random_Number(URand)

     DO k = 1, Org(I)%Width
       n_real = -1
       ! find maximum lability particles
         do m = 1, Org(i)%Width
           IF (Block(m)%Class .eq. p) then
             if (Particle(Block(m)%Value)%OM%OMact .gt. n_real) then
               n_real = Particle(Block(m)%Value)%OM%OMact
               o = m
             end if
           End if
         End do

         IF (n_real .eq. -1) Return

       EAT = .FALSE.
         If (Lability_ON) then                                 ! lability searching function is on
           If (Urand(k) .gt. Org(i)%Ingest%Selectivity) Eat = .true.
           If (n_real .ge. 3) Eat = .true.                          ! lability threshold: < class 2 particles are rejected

         Else                                                  ! if lability is not important .. make a random choice
           If (Urand(k) .gt. Org(i)%Ingest%Selectivity) Eat = .true.

         End if

         If (Eat) then
             If (Guts(k,i)%Class .eq. w) then
               O2_buff = O2(Block_loc(o)%Y, Block_loc(o)%X) 
               Guts(k,i) = Block(o)
               Matrix(Block_loc(o)%Y, Block_loc(o)%X) = CellContent(w,w)
               if (oxygen_ON) then 
                 O2(Block_loc(o)%Y, Block_loc(o)%X) = O2(org_loc(k,i)%Y,org_loc(k,i)%X) 
                 O2(org_loc(k,i)%Y,org_loc(k,i)%X)  = O2_buff
               end if 
               
               Block(o) = CellContent(w,w)
               IngestionHistory(i,1) = IngestionHistory(i,1) + 1   ! increment the running average
	           
			   
	   select case(Org(i)%Orientation)
	     case(0)
		   Vb(1,i) = Vb(1,i) + 2.5
	     case(180)
		   Vb(1,i) = Vb(1,i) - 2.5
	     case(90)
		   Ub(1,i) = Ub(1,i) - 2.5
	     case(270)
		   Ub(1,i) = Ub(1,i) + 2.5
		 case default
		   print *, 'error in setting flow const---ingest',Org(i)%Orientation
		 end select
	   
             End if
         End if

   ENd DO
   
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after ingest"; write(file_log,*)time, "err after ingest"; end if 
         end if 
		 
   End SUBROUTINE Particles_ingest

   ! **************************************************************************

   SUBROUTINE Particles_egest(i)
     ! organism i egests particles

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Logical  :: IsWater
   Integer :: i, j, k, loc, centre
   real :: biodecay               !!!YK 
   real :: timeIncrement, ox_buf
   type(coordinates) :: buff_loc
   
   integer :: xx, yy
   type(coordinates) :: org_mid
   
   biodecay = lability_decayConstant(N_LabilityClasses)*bio_fact*fact

   If (LocalMixing) return                            ! no need if forcing local mixing.

   If (Org(i)%Can%Egest) then

     Call Organism_peristalsis(i)      ! Peristalisis/compaction .. move down the gut contents

     k = -1
     centre = Org(I)%Width/2
       ! starting location of the placement of faeces (at the midpoint of the body width)

       Do j = Org(i)%Gut%Capacity - (Org(i)%FaecesWidth-1), Org(i)%Gut%Capacity
         ! egest the oldest particles (put into table Faeces which will be emptied with org movement)

         If (Guts(j,i)%Class .eq. p) then    ! if it is a particle
           k = k + 1
           loc = centre + (-1  ** k ) * k    ! distribute over the space

             If (IsWater(Tail(loc))) then
               IF (Particle(Guts(j,i)%Value)%OM%OMact .gt. 0) then
               Matrix(Tail(loc)%Y, Tail(loc)%X) = Guts(j,i)  !! fill the position of water in tail with a particle 
               Particle(Guts(j,i)%Value)%loc = Tail(loc)   !! location of particle is set at tail position 
               if (oxygen_ON) then 
                 ox_buf = O2(Tail(loc)%Y,Tail(loc)%X)%oxygen 
                 O2(Tail(loc)%Y,Tail(loc)%X)%oxygen = 0.0  !! because now tail is particle, oxygen conc. should be zero
			   else 
			     ox_buf = 1.
               end if 
			   
			   if (.not.resp_ON) then 
                 If (Lability_ON) then       ! lability functions are on ..
                     TimeIncrement = Real(Time - Particle(Guts(j,i)%Value)%OM%OM_Time0)*TimeScale
                     Particle(Guts(j,i)%Value)%OM%OMact = Particle(Guts(j,i)%Value)%OM%OMact   &
                              - biodecay*Particle(Guts(j,i)%Value)%OM%OMact*ox_buf*iox/OM_uni*TimeIncrement/365.
                     Particle(Guts(j,i)%Value)%OM%OM_Time0 = Time
                     
                     if (Particle(Guts(j,i)%Value)%OM%OMact < 0) Particle(Guts(j,i)%Value)%OM%OMact = 0.0
                 End if
			  end if 

               Guts(j,i) = CellContent(w,w)
	   
	   
	   select case(Dir_rec(org(i)%length,i))
	     case(0)
		   Vb(2,i) = Vb(2,i) + 2.5
	     case(180)
		   Vb(2,i) = Vb(2,i) - 2.5
	     case(90)
		   Ub(2,i) = Ub(2,i) - 2.5
	     case(270)
		   Ub(2,i) = Ub(2,i) + 2.5
		 case default
		   print *, 'error in setting flow const---egest',Dir_rec(org(i)%length,i)
		 end select
               
               end if 
             End if
         End if
       End do
       
   End if
   
         If (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after egest"; write(file_log,*)time, "err after egest"; end if 
         end if 
		 
   End Subroutine Particles_egest

   ! ********************************

   Subroutine Organism_peristalsis(i)
     ! shift gut contents down in the gut ..

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Integer :: i, j, k
   Type(CellContent) :: Image_Gut(Org(i)%Gut%Capacity)
   type(oxygen_conc) :: O2_buff(Org(i)%Gut%Capacity)          !!! YK

   integer :: xx, yy
   logical :: IsWater
   
   ! make an image of the gut contents
       DO j = 1, Org(i)%Gut%Capacity
         Image_Gut(j) = Guts(j,i)
         Guts(j,i) = CellContent(w,w)
       end do

   ! Pack the particles down the gut and count
       Org(i)%Gut%Content = 0
       k = Org(i)%Gut%Capacity + 1

       Do j = Org(i)%Gut%Capacity, 1, -1      ! do the resorting
         IF (Image_Gut(j)%Class .eq. p) then
           Org(i)%Gut%Content = Org(i)%Gut%Content + 1
           k = k - 1
           Guts(k,i) = Image_Gut(j)
         End if
       End do
       
         if (errChk) then 
         call Make_matrix_chk()
         call Matrix_err_chk()
         if (errDetect) then; print *,time, "err after peristalsis"; write(file_log,*)time, "err after peristalsis"; end if 
         end if 
		 
   End Subroutine Organism_peristalsis

   ! ***********************************************************************

   SubRoutine OutputData()
     ! main data output routine

   Use GlobalVariables
   use O2_diffusion
   Implicit None

   Call Output_stats()
   Call Output_Graphic()
     Call Output_Ascii() !  if ascii data desired
     Call Output_txtImg() 
     Call Output_O2txtImg() 

   End Subroutine OutputData

   ! ***********************************************************************

   SubRoutine Output_ascii()
     ! output the data in ascii format

   Use GlobalVariables
   Implicit None
   Logical :: IsWater, IsPArticle, IsOrganism
   Integer :: X, Y, Code
   Character*21 AsciiFile, AsciiFile2, numtemp
   integer :: tmppp
   integer :: txtimg(N_Col)
   
   tmppp = 1
    
   AsciiFile = CurrentTime //trim(adjustl(today))//'.ascii'
   write(numtemp,'(i10.1)') SaveTime
   OPEN(unit = File_ASCII, File = 'C:/Users/YK/Desktop/biot-res/data'//trim(adjustl(today))    &
                //trim(adjustl(numtemp))//'.ascii', status = 'unknown') 
   OPEN(unit = File_txtImg_2, File = 'C:/Users/YK/Desktop/biot-res/data_2'//trim(adjustl(today))    &
                //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO X = 1, N_col
     DO Y = 1, N_row
       If (IsParticle(CoOrdinates(Y,X))) then
         Code = 10
       Else if (IsWater(CoOrdinates(Y,X))) then
         Code = 0
         WRITE(File_txtImg_2,*) X, Y, 1 
       Else if (IsOrganism(CoOrdinates(Y,X))) then
            If (MAtrix(Y,X)%Value .le. Org(Matrix(Y,X)%Class)%HeadSize) then
            code = 30
            else 
            Code = 20
            end if
            WRITE(File_txtImg_2,*) X, Y, 1  
       End if
       if (tmppp == X) then 
       WRITE(File_ASCII,*) X, Y, code, Matrix(Y,X)
       else if (tmppp /= X) then 
       write(file_ascii, *) ""
       WRITE(File_ASCII,*) X, Y, code, Matrix(Y,X)
       tmppp = x
       end if
     END Do
     END Do
   Close(File_ASCII)

   End Subroutine Output_ascii

   ! *********************************************************

   SubRoutine Output_txtImg()
     ! output the data in ascii format

   Use GlobalVariables
   Implicit None
   Logical :: IsPArticle, IsOrganism
   Integer :: X, Y
   Character*21 numtemp
   integer :: txtimg(N_Col), particle_ID
   real :: Lab_real, txtimg_real(N_col)
   integer :: i , j
   ! logical :: rec_detail = .false.
   logical :: rec_detail = .true.
    
   write(numtemp,'(i10.1)') Time
   OPEN(unit = File_txtImg, File = 'C:/Users/YK/Desktop/biot-res/data'//trim(adjustl(today))          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO Y = 1, N_row
     txtimg = 0 
     DO X = 1, N_col
       If (IsParticle(CoOrdinates(Y,X))) then
         txtimg(X) = 1
       else if (IsOrganism(CoOrdinates(Y,X))) then
         if (matrix(y,x)%value <= Org(matrix(y,x)%class)%headSize) then
           txtimg(x) = 2
         else 
            txtimg(x) = 3
         end if 
       End if
     END Do
     write(File_txtImg, *) (txtimg(x), x = 1, n_col)
     END Do
   Close(File_txtimg)
   
   
   OPEN(unit = File_txtImg, File = 'C:/Users/YK/Desktop/biot-res/data_lability'//trim(adjustl(today))          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO Y = 1, N_row
     txtimg_real = 0 
     DO X = 1, N_col
       If (IsParticle(CoOrdinates(Y,X))) then
         Particle_ID = Matrix(Y,X)%Value
         Lab_real = particle(particle_ID)%OM%OMact
         txtimg_real(X) = 1 + Lab_real
       End if
     END Do
     write(File_txtImg, *) (txtimg_real(x), x = 1, n_col)
     END Do
   Close(File_txtimg)

   if (rec_detail) then 
      OPEN(unit = File_txtImg, File = 'C:/Users/YK/Desktop/biot-res/ptcl_list'//trim(adjustl(today))          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO Y = 1, N_row
     txtimg = 0 
     DO X = 1, N_col
         if (Matrix(Y,X)%Class==p) then 
           write(File_txtimg,*) Matrix(Y,X)%Class, Matrix(y,x)%value,Particle(Matrix(y,x)%value)%OM%OMact 
         else 
           write(File_txtimg,*) Matrix(Y,X)%Class, Matrix(Y,X)%Value 
         End if
     END Do
     END Do
   Close(File_txtimg)
   
      OPEN(unit = File_txtImg, File = 'C:/Users/YK/Desktop/biot-res/ptcl_list_2'//trim(adjustl(today))          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO i = 1, N_cell
         write(File_txtimg,*) i, Particle(i)%loc%Y,Particle(i)%loc%X, Particle(i)%plane
     END Do
   Close(File_txtimg)

   
   OPEN(unit = File_txtImg, File = 'C:/Users/YK/Desktop/biot-res/guts_list'//trim(adjustl(today))          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO i = 1, N_ind
     do j = 1,Org(i)%Gut%Capacity 
	    if (Guts(j,i)%Class == p) then 
           write(File_txtimg,*) Guts(j,i)%Class, guts(j,i)%Value, particle(guts(j,i)%value)%OM%OMact
		else 
		   write(File_txtimg,*) Guts(j,i)%Class, guts(j,i)%Value
		end if 
     END Do
     write(file_txtImg,*) ""
     end do 
   Close(File_txtimg)
   
   
   OPEN(unit = File_txtImg, File = 'C:/Users/YK/Desktop/biot-res/org_loc_list'//trim(adjustl(today))          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO i = 1, N_ind
     do j = 1, Org(i)%bodySize 
           write(File_txtimg,*) Org_loc(j,i)%Y, Org_loc(j,i)%X
     END Do
     write(file_txtimg, *) ""
     end do 
   Close(File_txtimg)
   
   end if 
   
   End Subroutine Output_txtimg
   
   ! *********************************************************

   SubRoutine Make_matrix_chk()
   
   use GlobalVariables
   implicit none 
   integer :: i, j, xx, yy
   integer :: txtimg(N_Col)
   
   
   matrix_chk = cellcontent(0,0)  ! all water
   errDetect = .false.
   
   do i = 1, N_cell
      if ((particle(i)%loc%x == 0) .and. (particle(i)%loc%y == 0)) cycle
      if (any(Guts(:,:)%value == i)) cycle
      matrix_chk(particle(i)%loc%y, particle(i)%loc%x) = cellcontent(i,p)
   end do 
   
   do i = 1, N_ind
      do j = 1, Org(i)%bodySize
        if (matrix_chk(Org_loc(j,i)%Y,Org_loc(j,i)%X)%class == p) then  !!  if particle (not in guts) exists
          print *, time, "error in (making mtx_chk): particle exist where organism",i," exist at", Org_loc(j,i)%Y,Org_loc(j,i)%X
          write(file_log,*) time, "error in (making mtx_chk): particle exist where organism",i," exist at",& 
		                    Org_loc(j,i)%Y,Org_loc(j,i)%X
          errDetect = .true.
        end if 
        matrix_chk(Org_loc(j,i)%Y,Org_loc(j,i)%X) = cellcontent(j,i)
      end do
   end do
   
   
   End Subroutine Make_matrix_chk
   
 ! ****************************************************************************************************

   SubRoutine Matrix_err_chk()
   
   use GlobalVariables
   implicit none 
   integer :: xx, yy
   
   errDetect = .false.
   
   do yy = 1, N_row
     do xx = 1, N_col
       if (matrix(yy,xx)%class /= matrix_chk(yy,xx)%class) then 
            write(file_log,*) time, 'error in water/organism assignment in matrix (err_chk_mtx)', yy, xx, & 
                matrix_chk(yy,xx)%class, matrix(yy,xx)%class, matrix_chk(yy,xx)%value, matrix(yy,xx)%value
            errDetect = .true.
       else 
         if (matrix(yy,xx)%value /= matrix_chk(yy,xx)%value) then
           write(file_log,*) time, '(err_chk_mtx) Id is different at', yy, xx, & 
                matrix_chk(yy,xx)%class, matrix(yy,xx)%class, matrix_chk(yy,xx)%value, matrix(yy,xx)%value
           errDetect = .true.
         end if 
       end if 
     end do
   end do 
   
   End Subroutine Matrix_err_chk
   
 ! ****************************************************************************************************

   SubRoutine Output_txtImg_chk()
   
   use GlobalVariables
   implicit none 
   integer :: i, j, xx, yy
   integer :: txtimg(N_Col)
   Character*21 numtemp
   
   
   write(numtemp,'(i10.1)') Time
   
   OPEN(unit = File_txtImg, File = 'C:/Users/YK/Desktop/biot-res/data_chk'//trim(adjustl(today))          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO YY = 1, N_row
     txtimg = 0 
     DO XX = 1, N_col
       If ((matrix_chk(yy,xx)%class == p).and.(matrix_chk(yy,xx)%value > 0)) then
         txtimg(XX) = 1
       else if ((matrix_chk(yy,xx)%class == p).and.(matrix_chk(yy,xx)%value == -1)) then 
         txtimg(XX) = 4
       else if (matrix_chk(yy,xx)%class >= 1) then
         if (matrix_chk(yy,xx)%value <= Org(matrix_chk(yy,xx)%class)%headSize) then
           txtimg(xx) = 2
         else 
            txtimg(xx) = 3
         end if 
       End if
     END Do
     write(File_txtImg, *) (txtimg(xx), xx = 1, n_col)
     END Do
   Close(File_txtimg)
   
   
   End Subroutine Output_txtimg_chk
   
 ! ****************************************************************************************************
   
   subroutine OrgDecay()
   
   use GlobalVariables
   implicit none 
   integer :: y, x
   real :: lab_real_new
   real :: aaa
   
   do y = 1, n_row
     do x = 1, n_col
     
     if (Matrix(y,x)%class == p) then 
       aaa= Lab_real_new(Matrix(y,x)%value)
     end if 
     
     end do 
   end do 
   
   end subroutine OrgDecay
   ! *********************************************************
   subroutine OrgResp(i)
   
   use GlobalVariables
   implicit none 
   integer :: k, o, i
   real :: resp_rate, lab_max, ox_resp, TimeIncrement
   
   lab_max = -1.0
   do k = Org(i)%Gut%Capacity, 1, -1
     if (guts(k,i)%class == p) then 
       
       if (particle(guts(k,i)%value)%OM%OMact > lab_max) then
         lab_max = particle(guts(k,i)%value)%OM%OMact 
         o = k
       end if 
       
     end if 
   end do 
   
   if (lab_max <= 0) return
   
   resp_rate = lab_max*Lability_decayConstant(N_LabilityClasses)*fact*bio_fact  ! calc wt% loss/yr once oxygen conc. is given
   
   ! resp_rate = resp_rate*1e-1   !! if decomposition rate is increased
   
   ox_resp = 0.0
   
   do k = 1, Org(i)%headSize
     O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen_use = &  
       O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen_use + & 
       resp_rate
     ox_resp = ox_resp + &
       resp_rate*O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen*TimeScale/365.*iox/OM_uni
   end do 
   
   
   particle(guts(o,i)%value)%OM%OMact = particle(guts(o,i)%value)%OM%OMact - ox_resp
   
   if (particle(guts(o,i)%value)%OM%OMact < 0) particle(guts(o,i)%value)%OM%OMact = 0.0
   
   RespHistory(i,1) = 1
   
   
   end subroutine OrgResp
   ! *********************************************************
   subroutine fluxes()
   
   use GlobalVariables
   implicit none 
   integer :: xx, yy, k, i
   integer :: xg, xp, yp, yg
   real :: tmpV, tmpU
   
   TotO2Dif = 0.0
   do xx = 1, n_col
     totO2Dif = totO2Dif +   &
       (O2(y_int,xx)%oxygen - O2(y_int+1,xx)%oxygen)*edif(y_int,xx)*iox*1e-3   &
          /pixelSize*(0.2*pixelSize)
   end do
   
   TotO2Dif = TotO2dif/0.2/(n_col*pixelSize)  !  mol /cm2 / yr flux
   
   TotO2Adv = 0.0
   do xx = 1, n_col
     TotO2Adv = TotO2Adv +  &
	   O2(y_int,xx)%oxygen*Vo(xx,y_int)*iox*1e-3*(0.2*pixelSize*pixelsize)
   end do 
   
   TotO2Adv = TotO2Adv/0.2/(n_col*pixelSize) 
   
   if (time == 1) then 
     pretotO2 = 0.0
   else 
     pretotO2 = totO2
   end if 
   
   totO2 = sum(O2(:,:)%oxygen)*iox*1e-3*(pixelSize)/(n_col)
   
   TotOrgDecay = 0.0
   
   do yy = 1, n_row
     do xx = 1, n_col
       TotOrgDecay = TotOrgDecay +   &
        O2(yy,xx)%Oxygen_use*O2(yy,xx)%oxygen*iox*1e-3*(0.2*pixelSize*pixelsize)
      end do 
    end do 
	
   TotResp = 0.0
   if (resp_ON) then 
   do i = 1, N_ind
   do k = 1, Org(i)%headSize
   if (matrix(Org_loc(k,i)%Y,Org_loc(k,i)%X)%class/=i) print *, "error in calc. flux"
   TotResp = TotResp +  &   
       O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen_use*O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen*iox   &
         *1e-3*0.2*(pixelSize)*(PixelSize)
   end do 
   end do
   end if 
   
   ToTAbio = TotOrgDecay - TotResp
   
   TotAbio = TotAbio/0.2/(n_col*pixelSize)
   
   TotResp = TotResp/0.2/(n_col*pixelSize) ! 
   
   TotOrgDecay = TotOrgDecay /0.2/(n_col*pixelSize) 
   
   write(File_flux,*) Time*Timescale, TotO2dif,TotO2Adv, TotAbio, TotResp, TotOrgDecay, (totO2-pretotO2)/timescale*365.0
   write(File_flux_txt,*) Time*Timescale, TotO2dif,TotO2Adv, TotAbio, TotResp, TotOrgDecay, (totO2-pretotO2)/timescale*365.0
   write(*,*) "FLUXES --- :", TotO2dif,TotO2Adv,TotAbio, TotResp, TotOrgDecay, (totO2-pretotO2)/timescale*365.0
   
   
   End Subroutine fluxes
   
   ! ***********************************************************   
   
   subroutine org_midloc(i, startj, endj,  org_mid)
   
   use globalvariables
   implicit none 
   type(coordinates) :: org_mid
   integer :: i, startj, endj, j
   integer :: maxx, maxy, minx, miny
   integer :: org_copy(endj - startj + 1)
   
   maxx = maxval(org_loc(startj:endj,i)%x)
   maxy = maxval(org_loc(startj:endj,i)%y)
   minx = minval(org_loc(startj:endj,i)%x)
   miny = minval(org_loc(startj:endj,i)%y)
   
   org_mid%y = int((maxy + miny )/2.0)
   
   if ((maxx - minx <= 2*org(i)%width)) then 
   org_mid%x = int((maxx + minx)/2.0)
   else if (maxx - minx > 2*org(i)%width) then 
      do j = 1, endj -startj + 1
	     if (org_loc(j,i)%x >=n_col/2) then 
		    org_copy(j) = org_loc(startj + j-1,i)%x - n_col
	     else 
		    org_copy(j) = org_loc(startj+j-1,i)%x
		 end if 
	  end do
   maxx = maxval(org_copy(1:endj-startj+1))
   minx = minval(org_copy(1:endj-startj+1))
   org_mid%x = int((maxx + minx )/2.0)
   else 
   print *, 'error in org_midloc'
   end if 
   
   if (org_mid%x < 1) org_mid%x = org_mid%x + n_col  
   if (org_mid%x > n_col) org_mid%x = org_mid%x - n_col  
   
   end subroutine org_midloc
   
   ! *************************************************************

   Subroutine Output_stats()

   Use GlobalVariables
   use ieee_arithmetic
   use mod_pfit

   Implicit None
   Logical :: IsParticle, InSediments, IfNAN

   Integer :: X, Y, Window, N_Samples, j, k, m, n
   Integer :: part_rowsum, Particle_ID

   Real :: Porosity_local, por_mean, por_stderr, por_rowsum
   Real :: URand(2), Depth, Db, minvalue
   Real :: Activity_New, Activity_mean, act_rowsum

   Character*8 :: Timevalue

   Integer :: lab_rowsum(N_LabilityClasses), Particle_Lability, P_count
   real :: Lab_real_New
   real :: Sum_Dy, dy, Sum_dx, dx, SUM_RMS, RMS, DEV_X, DEV_Y
   real :: Sum_Dy2(N_Row), Sum_dx2(N_Row), SUM_RMS2(N_Row), RMS2(N_Row), DEV_X2(N_Row), DEV_Y2(N_Row), P_count2(N_Row)

   Real :: Depth_vector(N_Row), Activity_vector(N_Row), weight_vector(N_Row)
   Real :: AMACH
   Real :: Db_lability, Lability_vector(N_Row, N_LabilityClasses)

   INTEGER NRMISS

   REAL AOV(15), CASE(N_Row,12), COEF(2,5), COVB(2,2), TESTLF(10), Profile(N_Row,3)
   REAL STAT(15,1)
   Real, allocatable :: PorosityData(:)
   
   real :: avx, avy, sumx2, sumy2, sumxy, slp, intcpt, crr     !!!  YK 5/22/2017
   
   real, allocatable ::  xsel(:), ysel(:), xlsel(:), ylsel(:)   !!!  YK 5/22/2017
   
   integer :: cnt, yint, yfin, cnt2, cnt3, row
   integer ::   list(3,n_row)
   
   double precision, allocatable :: amx(:,:), amx2(:,:), ymx(:), ymx2(:), cmx(:)
   
   real :: Db_list(n_row),  Activity_vector2(N_Row), Db_list2(n_row)
   
   logical :: matrixCheck, int_activity_1
   
   integer info
   double precision imbr
   integer, allocatable :: ipiv(:)
   
   real(8), allocatable :: a(:), sig(:), coeff(:,:), cov(:,:)
   
   real :: poly
   
   integer :: poly_fit_n = 6
   
   real :: lab_rowsum_real , lab_mean
   real :: lab_real_vect(N_row)

   ! IMSL Routine RONE is used to do the weighted linear regression:
   ! CALL RONE (NOBS(=N_row), NCOL, DataMatrix, ldx(=n_row), INTCEP(1=yes), IRSP(dep.var.=2), IND(indep.var.=1),
   !            IFRQ(freq option=0), IWT(weight.col=3), IPRED(prediction.option=0), CONPCM(conf.), CONPCP(conf2),
   !            IPRINT(print option=0), AOV(output), COEF(output), LDCOEF(leading dim of coef=2), COVB(output), LDCOVB(=2)
   !            TESTLF (lack of fit tests), CASE(ouput), LDCASE(dim=N_row), NRMISS(n missing values))

   !IMSL Routine UVSTA is used to obtain means and variances of porosities

   !! file names !! 
   
! File_Profile,  /biot-res/DepthProfiles          
! File_Profile2,  /biot-res/DepthProfiles2        
! File_Profile3,  /biot-res/DepthProfiles3        
! File_Profile_Re, /biot-res/DepthProfiles_Re        
! File_Scaling, /biot-res/PorosityScale        
! File_Activity, /biot-res/ActivitySlope        
! File_Displace, /biot-res/MeanDisplacement       
! File_Displace2, /biot-res/MeanDisplacement2        
! Poly_fit, /biot-res/polyFit
 
   If (.not. SaveData) Return

   ! DO SOME INITIAL DATA MANIPULTAIONS AND CALCULATIONS

     Timevalue = CurrentTime(1:8)
     InSediments = .false.
     MinValue = 1.

     sUM_DY = 0.
     SUM_DX = 0.
     SUM_RMS = 0.
     P_COUNT = 0

     DO Y = 1, n_row   ! INTEGRATE ACROSS
       
       part_rowsum = 0
       lab_rowsum = 0
       por_rowsum = 0.
       act_rowsum = 0.

       sUM_DY2(y) = 0.
       SUM_DX2(y) = 0.
       SUM_RMS2(y) = 0.
       P_COUNT2(y) = 0
       
       lab_rowsum_real = 0.

       DO X = 1, N_Col

         por_rowsum = por_rowsum + Porosity_local(CoOrdinates(Y,X), WindowSize)

         If (IsParticle(CoOrdinates(Y,X))) then

           Particle_ID = Matrix(Y,X)%Value
           part_rowsum = part_rowsum + 1
           act_rowsum = act_rowsum + Activity_New(Particle_ID)
           Particle_Lability = Lab_real_New(Particle_ID)
		   
             If (Particle_Lability .gt. 0) lab_rowsum_real = lab_rowsum_real + Particle_Lability 

             dY = Particle(Particle_ID)%loc%Y - Particle(Particle_ID)%loc_Init%Y
             ! here the %Plane variable is keeping track of the number of time the particle has
             ! been wrapped around, from its original location (plane/window)

               Select Case (Particle(Particle_ID)%Plane)
                 Case(0)   !  still the same plane
                   dX =  Particle(Particle_ID)%loc%X - Particle(Particle_ID)%loc_Init%x
                 Case(:-1) ! wrapped left
                   dX = - ( N_col - Particle(Particle_ID)%loc%x + Particle(Particle_ID)%loc_Init%X   &
                   + N_Col * (Particle(Particle_ID)%Plane + 1) )
                 Case(1:)  ! wrapped right
                   dX =   ( N_col - Particle(Particle_ID)%loc_Init%x + Particle(Particle_ID)%loc%X   &
                   + N_Col * (Particle(Particle_ID)%Plane - 1) )
               End Select

           Sum_Dy = Sum_dy +  Real(dy)
           Sum_Dx = Sum_dx +  Real(dx)
           sUM_rms = sUM_RMS + Real(dY*dY + Dx*Dx)  ! square of the distance displaced

           Sum_Dy2(y) = Sum_dy2(y) +  Real(dy)
           Sum_Dx2(y) = Sum_dx2(y) +  Real(dx)
           sUM_rms2(y) = sUM_RMS2(y) + Real(dY*dY + Dx*Dx)  ! square of the distance displaced
           
           P_count2(y) = P_count2(y) + part_rowsum

         End if
       END Do    ! for each column
       

         P_count = P_count + part_rowsum    ! P_count counts the total number of particles in the whole matrix
		 
       ! set up data vectors for linear regression analysis

         Depth = Real(Y) * PixelSize
         Depth_vector(y) = Depth
         por_mean = por_rowsum / Real(N_Col)
         lab_mean = lab_rowsum_real/ real(N_col)
           If (part_rowsum .gt. 0) then
             Activity_mean = act_rowsum / Real(part_rowsum)
               IF (Activity_mean .lt. 0.0001) Activity_mean = 0
           ELSE
             Activity_mean = 0
           END IF
		   
         activity_vector2(y) = activity_mean  ! YK 
         lab_real_vect(y) = lab_mean
		 
         WRITE(File_Profile,*) Timevalue, Depth, Activity_mean, por_mean, lab_mean, MatrixOccupancy(Y)
         
         if (y .eq. n_row) WRITE(File_Profile,*) ""       !      YK 5/22/2017

         If (Activity_mean .le. 0.0001) then   ! remove small numbers from analysis
           Activity_vector(y) = 0         !!        YK  5/22/2017
           Weight_vector(y) = 9999
         Else
           Activity_vector(y) = log(Activity_mean)
           Weight_vector(y) = Activity_vector(y)
           iF (Activity_vector(y) .LT. MinValue) MinValue = Activity_vector(y)
         End if

         WRITE(File_Profile2,*) Timevalue, Depth, Activity_vector(y), (Weight_vector(y)), &
         MatrixOccupancy(Y) !      YK 5/22/2017
         
         if (y .eq. n_row) WRITE(File_Profile2,*) ""       !      YK 5/22/2017

         Do k = 2, N_LabilityClasses
           If (lab_rowsum(k) .le. 1) then  ! remove small numbers and zero values from analysis
             Lability_vector(y,k) = 0  ! YK  5/22/2017
           Else
             Lability_vector(y,k) = log(Real(lab_rowsum(k)))
           End if
         
         End do

         WRITE(File_Profile3,*) Timevalue, depth,(Lability_vector(y,k),  k = 2, N_LabilityClasses), &  
         (Lability_vector(y,k)*Lability_vector(y,k),  k = 2, N_LabilityClasses)         !      YK 5/22/2017
         
         if (y .eq. n_row) WRITE(File_Profile3,*) ""       !      YK 5/22/2017

     END Do  ! for each row

    ! do the linear regression estimate of the slope FOR THE TRACER

       Profile(:,1)=depth_vector
       Profile(:,2)=activity_vector
       Profile(:,3)=weight_vector - MinValue

       m = N_Row
       do j = 1, N_Row
         if ( Profile(j,2).eq. 0 ) m = m-1
       End do
       
        Allocate (xsel(m),ysel(m))
       cnt = 1
       do j = 1, N_Row
         if ( Profile(j,2).eq. 0 ) cycle
         xsel(cnt) = Profile(j,1); ysel(cnt) = Profile(j,2); cnt = cnt + 1
       End do
       if (cnt-1/=m) then
         print *, "something", cnt, m
         write(File_log, *) Time, "subroutine stats", cnt, m
       end if 

         If (m .gt. 5) then     ! if enough data ..
           avx = sum(xsel(:))/m
           avy = sum(ysel(:))/m
           sumx2 = sum((xsel(:)-avx)**2.0)
           sumy2 = sum((ysel(:)-avy)**2.0)
           sumxy = sum((ysel(:)-avy)*(xsel(:)-avx))
           slp = sumxy/sumx2
           intcpt = avy - slp*avx              ! YK, 5/22/2017
           if ((sumx2/=0).and.(sumy2/=0)) then 
           crr = sumxy**2.0d0/sumx2/sumy2     !  YK, 5/22/2017
           else 
           crr = 0.
           end if
           if (slp .ne. 0) then 
           Db = DecayConstant/(slp*slp)   !  cm^2 yr^-1
                       
               WRITE (File_Activity,*) TimeValue, ' 210Pb', Db, slp, intcpt, crr
             End if
         End if
      
    ! do the linear regression estimate of the slope FOR THE VARIOUS LABILITY CLASSES
       Do n = 2, N_LabilityClasses
         Profile(:,2)=Lability_vector(:,n)
         Profile(:,3)=Lability_vector(:,n)

         m = N_Row
         do j = 1, N_Row
           if (Profile(j,1) .eq. 0) m = m-1
         End do
           If (m .gt. 5) then
           avx = sum(Profile(:,1))/m
           avy = sum(Profile(:,2))/m        
           sumx2 = sum((Profile(:,1)-avx)**2.0)
           sumy2 = sum((Profile(:,2)-avy)**2.0)
           sumxy = sum((Profile(:,2)-avy)*(Profile(:,1)-avx))     
           slp = sumxy/sumx2     
           intcpt = avy - slp*avx              ! YK, 5/22/2017   
           if ((sumx2/=0) .and. (sumy2/=0)) then 
           crr = sumxy**2.0d0/sumx2/sumy2     ! YK, 5/22/2017
           else 
           crr = 0.
           end if          
           
           if (slp .ne. 0) then 
           Db_lability = Lability_decayConstant(n)*365.0/(slp*slp)      !! cm^2/yr
           
                 WRITE (File_Activity,*) TimeValue, n, Db_Lability, slp, intcpt, crr
               End if
           End if
       End do
   ! OUTPUT THE PARTICLE DEVIATIONS

     RMS = PixelSize * SQRT( SUM_RMS / rEAL(P_COUNT) )   ! CALC OF rms .. == SQRT(MEAN(SQUARED DEVIATIONS))
     DEV_Y = PixelSize * SUM_DY / rEAL(P_COUNT)
     DEV_X = PixelSize * SUM_Dx / rEAL(P_COUNT)

       WRITE(File_Displace,40) Timevalue, DEV_Y, DEV_X, RMS
   
     do y = 1, n_row
       if (P_COUNT2(y) == 0) cycle
       RMS2(y) = PixelSize *PixelSize * ( SUM_RMS2(y) / rEAL(P_COUNT2(y)) )   ! CALC OF rms .. == SQRT(MEAN(SQUARED DEVIATIONS))
       DEV_Y2(y) = PixelSize * SUM_DY2(y) / rEAL(P_COUNT2(y))
       DEV_X2(y) = PixelSize * SUM_Dx2(y) / rEAL(P_COUNT2(y))
       
       if (time.ne.0) then 
       Db_list(y) = RMS2(y)*365./time/timescale/4.
       else 
       Db_list(y) = 0
       end if
	   
       write(File_Displace2,*) Timevalue, Depth_vector(y), DEV_Y2(y), DEV_X2(y), RMS2(y),   &
           Db_list(y), P_COUNT2(y)
       
       if (y.eq. n_row) write(file_displace2,*) ""
       
    end do   
     
     m = N_Row
     cnt = 0
       do j = 1, N_Row
         if (.not.ieee_is_nan(Db_list(j)) .and.(Db_list(j) > 0.)) cnt = cnt +1
       End do
        m = cnt
    
    if (m > max(5,poly_fit_n)) then
    
        DeAllocate (xsel,ysel)
        Allocate (xsel(m),ysel(m))
        Allocate (a(poly_fit_n),cov(poly_fit_n,poly_fit_n), coeff(m,poly_fit_n),sig(m))
     cnt = 1
     do j = 1, N_row
        if (cnt.eq. m) exit
        if (.not.ieee_is_nan(Db_list(j)) .and.(Db_list(j) > 0.)) then
        
        xsel(cnt) = depth_vector(j)
        ysel(cnt) = log10(Db_list(j))
        sig(cnt) = 1.0 
        
        cnt = cnt + 1
        end if
     end do
           
           call pfit(dble(xsel),dble(ysel),sig,a,coeff,cov)
     
     do j = 1, m-1
        
       write(poly_fit,*) xsel(j),ysel(j), poly(xsel(j),a,poly_fit_n)
     end do
     
     write(poly_fit,*) ""
     
     do y = 1, N_row
        Db_list2(y) = 10.0**poly(depth_vector(y),a,poly_fit_n)
     end do
     
     
     end if
     
     yfin = 0
     yint = 0
     cnt3 = 0
     list = 0
    
    do y = 1, n_row
       if (y .ge. yfin) then 
         cnt2 = 0
         j = y
         do while (((.not.ieee_is_nan(Db_list(j))).and. (j .le. n_row).and.(Db_list(j).ne.0))) 
           if (cnt2 .eq. 0) then 
             yint = j
           end if
           cnt2 = cnt2 + 1
           j = j + 1
           if (j.eq.n_row+1) then
             exit
           end if
         end do
       end if
       
       if (cnt2 /=0) then 
       yfin =yint + cnt2 - 1
       cnt3 = cnt3 + 1
       
       list(1,cnt3) = yint
       list(2,cnt3) = yfin
       list(3,cnt3) = cnt2
       
       end if
     
     end do
     
     matrixCheck = .false.
     int_activity_1 = .true.
     
     if (maxval(list(3,:)) > 50) then
        cnt2 = maxval(list(3,:)) 
        yint = maxval(list(1,maxloc(list(3,:))))
        yfin = maxval(list(2,maxloc(list(3,:))))
        allocate (amx(cnt2,cnt2))
        allocate (amx2(cnt2,cnt2))
        allocate (ymx(cnt2))
        allocate (ymx2(cnt2))
        allocate (cmx(cnt2))
        allocate (ipiv(cnt2))
        amx = 0.0
        ymx = 0.0
        ipiv = 0
        if (.not.ieee_is_nan(sum(Db_list2))) Db_list = Db_list2
        do y = yint , yfin 
          row = y -  yint + 1
          if (y ==   yint ) then
            amx(row,row) = (Db_list(y) - Db_list(y))/pixelSize/pixelSize + Db_list(y)*(-2.0)/pixelSize/pixelSize - decayconstant
            amx(row,row + 1) =  Db_list(y)*(1.0)/pixelSize/pixelSize 
            ymx(row) = (Db_list(y) - Db_list(y))/pixelSize/pixelSize* merge(Activity_vector2(y),1.      &
                       , (Activity_vector2(y).ne.0.0).and.(.not.ieee_is_nan(Activity_vector2(y)))     &
                       .and.(yint.ne.1).and.(.not.int_activity_1)) &
                     +  Db_list(y)*(-1.0)/pixelSize/pixelSize * merge(Activity_vector2(y),1.   &
                     ,  (Activity_vector2(y).ne.0.0).and.(.not.ieee_is_nan(Activity_vector2(y)))       &
                     .and.(yint.ne.1).and.(.not.int_activity_1))
          else if (y == yfin) then 
            amx(row,row) = (Db_list(y) - Db_list(y-1))/pixelSize/pixelSize + Db_list(y)*(-1.0)/pixelSize/pixelSize - decayconstant
            amx(row,row - 1) =  (Db_list(y) - Db_list(y-1))*(-1.0)/pixelSize/pixelSize +  Db_list(y)*(1.0)/pixelSize/pixelSize 
          else 
            amx(row,row) = (Db_list(y) - Db_list(y-1))/pixelSize/pixelSize + Db_list(y)*(-2.0)/pixelSize/pixelSize - decayconstant
            amx(row,row + 1) =  Db_list(y)*(1.0)/pixelSize/pixelSize 
            amx(row,row - 1) =  (Db_list(y) - Db_list(y-1))*(-1.0)/pixelSize/pixelSize +  Db_list(y)*(1.0)/pixelSize/pixelSize 
          end if
          
        end do
        
        if (matrixCheck) then
        
        open (unit=100, file="C:/Users/YK/Desktop/biot-res/check.txt",status = 'unknown')
        do y = yint, yfin
          row = y - yint + 1
          write(100,*) (amx(row, j), j = 1, cnt2)
        end do 
        close(100)
        
        open (unit=200, file="C:/Users/YK/Desktop/biot-res/check2.txt",status = 'unknown')
        do y = yint, yfin
          row = y - yint + 1
          write(200,*) ymx(row)
        end do 
        close(200)
        
        end if
        
        call DGESV(cnt2,1,amx,cnt2,IPIV,ymx,cnt2,INFO) 
        
        do j = 1, cnt2
        write(File_Profile_Re,*) depth_vector(j+yint-1),ymx(j),Activity_vector2(j+yint-1), ymx(j)/(Activity_vector2(j+yint-1))   &
              , merge(Activity_vector2(yint),1.,.not.int_activity_1)                     &
              *exp(-sqrt(decayconstant/Db)*(depth_vector(j+yint-1)-depth_vector(yint)))  &
              , merge(Activity_vector2(yint),1.,.not.int_activity_1)                         &
              *exp(-sqrt(decayconstant/Db)*(depth_vector(j+yint-1)-depth_vector(yint)))  &
                /(Activity_vector2(j+yint-1))
        end do
        
        write(File_Profile_Re,*) ""
          
     end if 
     
     do y = 1, cnt3
     write(File_oneD, *) list(1,y), list(2,y), list(3,y)
     end do
     
     write(File_oneD, *) ""

     ! OUTPUT SCALING FUNCTION

       If (OutputScaling) then
           N_Samples = 500
           Window = 1
           Allocate (PorosityData(N_Samples))

             Do While (Window .lt. (N_Row)/2)
                 Do j = 1, N_Samples
                  ! obtain a random value within the borders of the matrix,
                  ! adjusted for the size of the sampling window

                   Call Random_Number(URand)
                     X = URand(1) * Real(N_Col - 2 * Window) + 1 + Window
                     Y = URand(2) * Real(N_RowSed - 2 * Window) + 1 + Window + N_RowWater
                     PorosityData(j) = Porosity_local(CoOrdinates(Y,X), Window)
                 End do
				 
                 por_mean = sum(PorosityData(:))/N_samples
                 por_stderr   = (sum((PorosityData(:)-por_mean)**2.0)/(N_samples - 1.0))**0.5

                 WRITE(File_Scaling,*) Timevalue, Real(1+2*Window)*PixelSize, por_mean, por_stderr

                 Window = 2 * Window
             End do
       End if

10  FORMAT (A8, 3f9.4, 6I7)
20  FORMAT (A8, 3f10.6)
30  FORMAT (A8, A6, 8f10.4, f6.0)
31  FORMAT (A8, I6, 8f10.4, f6.0)
40  FORMAT (A8, 3f11.6)


   End Subroutine Output_stats

   ! ***********************

   Subroutine LocateOrganisms_scan()
     ! locate the organisms in the matrix and accumulate the verticle depth distribution in "MatrixOccupancy"

   Use GlobalVariables
   Implicit None
   Integer :: i, j

   j = 1 ! the index 1 of the body used to tag organisms ...
     Do i = 1, N_Ind
       MatrixOccupancy(Org_loc(j,i)%Y) = MatrixOccupancy(Org_loc(j,i)%Y) + 1
     End do

   End Subroutine LocateOrganisms_scan

   ! ***********************


   Subroutine Output_Graphic()
     ! Create images of the matrix

   Use GlobalVariables
   use gif_util

   implicit none
   Integer    :: X, Y
   integer, allocatable   :: image_data(:,:)
   Integer(2)   :: LookupPorosityColour, LookupActivityColour, LookupLabilityColour
   
   character*21 numtemp

   type(color), dimension(0:31)  :: palette

         palette(red)   = color(240,  0,  0)       ! red
         palette(green) = color(0  ,240,  0)       ! green
         palette(blue)  = color(0  ,  0,240)       ! blue
         palette(white) = color(240,240,240)       ! white
         palette(black) = color(  0,  0,  0)       ! black
         palette(yellow)= color(200,200,  0)       ! yellow

         palette(green1) = color(0,120,0)
         palette(green2) = color(0,160,0)
         palette(green3) = color(0,190,0)
         palette(green4) = color(0,220,0)
         palette(green5) = color(0,240,0)

         palette(blue1) = color(0,0, 90)
         palette(blue2) = color(0,0,120)
         palette(blue3) = color(0,0,150)
         palette(blue4) = color(0,0,180)
         palette(blue5) = color(0,0,200)

         palette(red1) = color(100,0,0)
         palette(red2) = color(128,0,0)
         palette(red3) = color(160,0,0)
         palette(red4) = color(192,0,0)
         palette(red5) = color(224,0,0)

         palette(turq1) = color(0,150,150)
         palette(turq2) = color(0,208,208)
         palette(turq3) = color(0,224,224)
         palette(turq4) = color(0,240,240)
         palette(turq5) = color(0,255,255)

         palette(grey1) = color(80 ,80 ,80)
         palette(grey2) = color(120,120,120)
         palette(grey3) = color(140 ,140 ,140)
         palette(grey4) = color(180,180,180)
         palette(grey5) = color(200,200,200)
         palette(grey6) = color(230,230,230)

   Allocate(image_data(n_col, n_row))
   
   write(numtemp,'(i10.1)') SaveTime

   If (PlotMatrix) then
         DO X = 1, N_Col
         DO Y = 1, N_Row
           image_data(X,Y) = LookupLabilityColour(Y,X)
         End Do
         End Do

       call writegif("C:/Users/YK/Desktop/biot-res/lab/lability"          &
                //trim(adjustl(today))//trim(adjustl(numtemp))//".gif", image_data, palette)

   End if

   If (.not. Movies_On) then
    If (PlotPorosity) then
         DO X = 1, N_Col
         DO Y = 1, N_Row
           image_data(X,Y) = LookupPorosityColour(Y,X) ! Write the stored data to the image array.
         End Do
         End Do
 
       call writegif("C:/Users/YK/Desktop/biot-res/por/porosity"           &
                 //trim(adjustl(today))//trim(adjustl(numtemp))//".gif", image_data, palette)

    End if

    If (PlotActivity) then
         DO X = 1, N_Col
         DO Y = 1, N_Row
           image_data(X,Y) = LookupActivityColour(Y,X)    ! Write the stored data to the image array.
         End Do
         End Do

       call writegif("C:/Users/YK/Desktop/biot-res/act/activity"//            &
                   trim(adjustl(today))//trim(adjustl(numtemp))//".gif", image_data, palette)

     End if
    End if

   End Subroutine Output_Graphic

   ! *******************************

   Function LookupLabilityColour(Y,X)
   ! return a colour value (lookup table)

   Use GlobalVariables
   Implicit none
   Logical   :: IsPArticle, IsWater, IsOrganism
   Integer(2) :: LookupLabilityColour
   Integer   :: Y,X

   If (IsPArticle(CoOrdinates(Y,X))) then
     Select Case (Particle(MAtrix(Y,X)%Value)%Lability)
       Case (-1)
         LookupLabilityColour = white  ! sand .. used when it is a plug
       Case (0)
         LookupLabilityColour = grey2  ! lowest lability
       Case (1)
         LookupLabilityColour = green1
       Case (2)
         LookupLabilityColour = green2
       Case (3)
         LookupLabilityColour = green3
       Case (4)
         LookupLabilityColour = green4
       Case (5)
         LookupLabilityColour = green5 ! highest lability
     End Select
   Else if (IsWater(CoOrdinates(Y,X))) then
     LookupLabilityColour = blue       ! water  .. blue
   Else if (IsOrganism(CoOrdinates(Y,X))) then
     If (MAtrix(Y,X)%Value .le. Org(Matrix(Y,X)%Class)%HeadSize) then
       LookupLabilityColour = red       ! head  ..  red
     Else
       LookupLabilityColour = yellow  ! other body parts .. yellow
     End If
   End if

   End Function LookupLabilityColour

   ! ***********************************************

   Function LookupPorosityColour(Y,X)
   ! return a colour value (lookup table)

   Use GlobalVariables
   Implicit none
   Integer(2) :: LookupPorosityColour
   Integer   :: Y,X
   Real :: Por, Delta_Porosity

   ! functions
   Real :: Porosity_local


   Por = Porosity_local(CoOrdinates(Y,X),WindowSize)
   Delta_Porosity = Porosity0 - Por

   If ( Por .gt. 0.98 ) then
     LookupPorosityColour = blue    ! water = 0
   Else
     If (Delta_Porosity .eq. 0) then
       LookupPorosityColour = black
     Else if (Delta_Porosity .gt. 0) then
       if (Delta_Porosity .lt. 0.05) then
         LookupPorosityColour = blue2
       Else if (Delta_Porosity .lt. 0.10) then
         LookupPorosityColour = blue3
       Else if (Delta_Porosity .lt. 0.15) then
         LookupPorosityColour = blue4
       Else
         LookupPorosityColour = blue5
       End if

     Else if (Delta_Porosity .lt. 0) then
       if (Delta_Porosity .gt. -0.05) then
         LookupPorosityColour = red1
       Else if (Delta_Porosity .gt. -0.10) then
         LookupPorosityColour = red2
       Else if (Delta_Porosity .gt. -0.15) then
         LookupPorosityColour = red3
       Else if (Delta_Porosity .gt. -0.20) then
         LookupPorosityColour = red4
       Else
         LookupPorosityColour = red5
       End if
     End if
   End if

   End Function LookupPorosityColour

   ! **********************

   Function LookupActivityColour(Y,X)
   ! return a colour value (lookup table)

   Use GlobalVariables
   Implicit none
   Integer(2) :: LookupActivityColour
   Real      :: Activity
   Integer   :: Y,X

       If (Matrix(Y,X)%Class .eq. p) then
         Activity = Particle(Matrix(Y,X)%Value)%Pb%Activity

           If (Activity .lt. 0.0001) then
             LookupActivityColour = turq1
           Else if (Activity .lt. 0.001) then
             LookupActivityColour = turq2
           Else if (Activity .lt. 0.01) then
             LookupActivityColour = turq3
           Else if (Activity .lt. 0.1) then
             LookupActivityColour = turq4
           Else
             LookupActivityColour = turq5
           End if
       Else
         LookupActivityColour = blue  ! water or organism = 0
       End if

   End Function LookupActivityColour

   ! *******************************************************************
   
   function Isfloating(i)
   
   use GlobalVariables
   implicit none
   logical Isfloating 
   logical :: IsParticle
   integer :: i,j, xx, yy, xxp, xxg, yyp, yyg, k , o
   Type(Coordinates) ::  Points(4)
   
   Isfloating = .false.
   k = 0
   do j = 1, Org(i)%BodySize
     xx = Org_loc(j,i)%X
	 yy = Org_loc(j,i)%Y
	 xxp = xx + 1
	 if (xxp > n_col) xxp = xxp - n_col
	 xxg = xx - 1
	 if (xxg < 1) xxg = xxg + n_col
	 yyp = yy + 1
	 yyg = yy - 1
	 if ((yyp > n_row) .or.(yyg < 1)) print *, "error in 'IsFloating'"
	 Points(1) = Coordinates(yy,xxp)
	 Points(2) = Coordinates(yy,xxg)
	 Points(3) = Coordinates(yyp,xx)
	 Points(4) = Coordinates(yyg,xx)
	 do o = 1, 4
	 if (IsParticle(Points(o))) k = k + 1
	 end do
   end do 
   
   if (k == 0) Isfloating = .true.
   
   End function Isfloating
   
   ! *******************************************************************

   Function IsWater(Point)
     ! tests whether the cell content of coordinate (Point%Y,Point%X) is water

   Use GlobalVariables
   Implicit None
   Logical IsWater
   Type (CoOrdinates) :: Point
   IsWater = .false.
   IF (Matrix(Point%Y,Point%X)%Class .eq. w) IsWater = .true.
   End Function Iswater

   !****************************************************************

   Function IsParticle(Point)
     ! tests whether the cell content of coordinate (Point%Y,Point%X) is a particle

   Use GlobalVariables
   Implicit None
   Logical IsParticle
   Type (CoOrdinates) :: Point
   
   IsParticle = .false.
   IF (Matrix(Point%Y,Point%X)%Class .eq. p) IsParticle = .true.
   
   End Function IsParticle

  !****************************************************************

   Function IsOrganism(Point)
     ! tests whether the cell content of coordinate (Point%Y,Point%X) is an organism

   Use GlobalVariables
   Implicit None
   Logical IsOrganism
   Type (CoOrdinates) :: Point
   
   IsOrganism = .false.
   IF (Matrix(Point%Y,Point%X)%Class .ge. 1) IsOrganism = .true.
   
   End Function IsOrganism

   !****************************************************************

   Function IsMoveable(Point)
     ! tests whether the cell contents are (particles or water) and so can be moved

   Use GlobalVariables
   Implicit None
   Logical IsMoveable, IsWater, IsPArticle
   Type (CoOrdinates) :: Point
   
   IsMoveable = .false.
   IF (IsParticle(Point) .or. IsWater(Point)) IsMoveable = .true.
   
   End Function IsMoveable

  !****************************************************************

   Function InMatrix(Point)
     ! tests whether the coordinates are within the matrix

   Use GlobalVariables
   Implicit None
   Logical InMatrix
   Type (CoOrdinates) :: Point

   InMatrix = .false.
   IF ((Point%Y .ge. 1) .and. (Point%Y .le. N_Row) .and. (Point%X .ge. 1) .and.  &
       (Point%X .le. N_Col)) InMatrix = .true.
       
   End Function InMatrix

  ! ************************8

   Function OnEdge(Point)
     ! tests whether the coordinates are on the lateral edges of the matrix

   Use GlobalVariables
   Implicit None
   Logical OnEdge
   Type (CoOrdinates) :: Point

   
   OnEdge = .false.
   IF ((Point%X .eq. 1) .or. (Point%X .eq. N_Col))   OnEdge = .true.

   End Function OnEdge

   !****************************************************************

   Function CanMoveTo(Point)
     ! tests whether the organism can move into a cell

   Use GlobalVariables
   Implicit None
   Logical CanMoveTo, IsMoveable, InMAtrix
   Type (CoOrdinates) :: Point

   CanMoveTo = .false.
   IF (InMAtrix(Point)) then
     If (IsMoveable(Point)) CanMoveTo = .true.
   End if

   End Function CanMoveTo

   !****************************************************************

   Function Porosity_local(Point,shift)
     ! returns the locally averaged porosity with window size of "shift*2+1"

   Use GlobalVariables
   Implicit None
   Logical :: InMAtrix, IsWater
   Integer :: Water_Count, Total_Count, X, Y, shift
   Real :: Porosity_local
   Type (CoOrdinates) :: Point

     Water_Count = 0
     Total_Count = 0

     Do X = Point%X-shift, Point%X+shift
     Do Y = Point%Y-shift, Point%Y+shift
       If (InMatrix(CoOrdinates(Y,X))) then
         Total_Count = Total_Count + 1
         If (IsWater(CoOrdinates(Y,X))) Water_count = Water_Count + 1
       End if
     End do
     End do

     Porosity_local = Real(Water_Count) / Real(Total_Count)

   End Function Porosity_local

   ! ******************************************************************

   Function Porosity(Y)
     ! returns the porosity of a given depth dependent upon the Porosity_Type random, uniform, linear, exponential

   Use GlobalVariables
   Implicit None
   Integer Y
   Real Porosity, Depth, DepthMAx, URand

     Call Random_Number(URand)

     Select Case (Porosity_Type)
       Case ("ran", "RAN", "Ran")
         Porosity = URand * Porosity0     ! Porosity0 is the porosity at the sediment/water interface
       Case ("uni", "UNI", "Uni")
         Porosity = Porosity0
       Case ("lin", "LIN", "Lin")
         Depth = Real(Y - N_RowWAter) * PixelSize
         DepthMAx = Real(N_Row - N_RowWAter) *  PixelSize
         Porosity = 1. - Real(Depth) / DepthMAx
       Case ("exp", "EXP", "Exp")
         Depth = Real(Y) * PixelSize
         Porosity = Porosity0 * EXP( - Porosity_DecayRate * Depth )
     End Select
     
   End Function Porosity

   ! ******************

   Function Lab_real(Y,x)
     ! returns the lability of a particle based upon the relative abundance of particles
     ! Y (the depth) is unused for the moment .. I plan to use it to vary the relative frequencies
     ! as a function of time or simply to layer the matrix

   Use GlobalVariables
   Implicit None
   Integer Y, x, j
   Real URand, depth, depthmax, lability_pr
   real :: lab_real

   Depth = Real(Y - N_RowWAter) * PixelSize
   DepthMAx = Real(N_Row - N_RowWAter) *  PixelSize
   Lability_pr = 1. - Real(Depth) / DepthMAx

       Call Random_Number(URand)
           If (Urand .lt. Lability_pr) then
               Do j = 2, N_LabilityClasses
                 If (Urand .lt. Lability_proportion(j)) then
                   Lab_real= real(j)/10.0
                   Return
                 End if
               End do
           END IF
   End Function Lab_real

  ! *********************

   Function Lab_real_New(ID)
     ! returns the Lability after probabilistic degradation

   Use GlobalVariables
   Implicit None
   Integer :: ID
   real :: Lab_real, Lab_real_New
   Real    :: URand, Pr_Cummulative, TimeIncrement
   real :: lab_decay
   
   logical :: inMatrix, IsParticle

     Lab_real = Particle(ID)%OM%OMact
     Lab_real_new = Particle(ID)%OM%OMact
     
       if (Lab_real .le. 0.01) Then
         Lab_real_New = 0
         Return
       END IF
	   
           TimeIncrement = Real(Time - Particle(ID)%OM%OM_Time0) * TimeScale / 365.   !  (yr)
		   if ((int(Lab_real*10)>=1).and.(int(Lab_real*10)<=N_LabilityClasses)) then 
           Pr_Cummulative = Lability_decayConstant(int(Lab_real*10))*TimeIncrement*Lab_real*fact  
		   else 
		   Pr_Cummulative = Lability_decayConstant(1)*TimeIncrement*Lab_real*fact
		   end if 
    
    if (oxygen_ON) then
    
    if (particle(ID)%loc%X == 1) then 
      if (particle(ID)%loc%Y == n_row) then 
    O2(particle(ID)%loc%Y,n_col)%Oxygen_use = &
            O2(particle(ID)%loc%Y,n_col)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
            O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
            O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    Lab_real_new = Lab_real -( &
             Pr_Cummulative*O2(particle(ID)%loc%Y,n_col)%oxygen  +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen    & 
             )*iox/OM_uni
      else if (particle(ID)%loc%Y == 1) then 
    O2(particle(ID)%loc%Y,n_col)%Oxygen_use = &
            O2(particle(ID)%loc%Y,n_col)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
            O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
            O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    Lab_real_new = Lab_real -( &
             Pr_Cummulative*O2(particle(ID)%loc%Y,n_col)%oxygen  +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen    & 
             )*iox/OM_uni
       else 
    O2(particle(ID)%loc%Y,n_col)%Oxygen_use = &
             O2(particle(ID)%loc%Y,n_col)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    Lab_real_new = Lab_real -( &
             Pr_Cummulative*O2(particle(ID)%loc%Y,n_col)%oxygen   +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen  + &
             Pr_Cummulative*O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen    & 
             ) *iox/OM_uni
        end if 
        
   else if (particle(ID)%loc%X == n_col) then 
      if (particle(ID)%loc%Y == n_row) then 
    O2(particle(ID)%loc%Y,1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    Lab_real_new = Lab_real - ( &
             Pr_Cummulative*O2(particle(ID)%loc%Y,1)%oxygen   +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen    & 
             )*iox/OM_uni
      else if (particle(ID)%loc%Y == 1) then 
    O2(particle(ID)%loc%Y,1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    Lab_real_new = Lab_real - ( &
             Pr_Cummulative*O2(particle(ID)%loc%Y,1)%oxygen   +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen    & 
             )*iox/OM_uni
       else 
    O2(particle(ID)%loc%Y,1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    Lab_real_new = Lab_real -( &
             Pr_Cummulative*O2(particle(ID)%loc%Y,1)%oxygen   +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen  +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen  +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen    & 
             )*iox/OM_uni
      end if
      
    else 
      if (particle(ID)%loc%Y == n_row) then 
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    Lab_real_new = Lab_real - ( &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen    & 
             )*iox/OM_uni
      else if (particle(ID)%loc%Y == 1) then 
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    Lab_real_new = Lab_real - ( &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen    & 
             )*iox/OM_uni
       else 
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)
             
    Lab_real_new = Lab_real -( &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen  +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen +  &
             Pr_Cummulative*O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen    & 
             )*iox/OM_uni
      end if
    end if
   
    ! in units of /yr
    
             Particle(ID)%OM%OM_Time0 = Time
             Particle(ID)%OM%OMact = Lab_real_New
    
     end if
     
     if (.not. oxygen_ON) then 
          
          Lab_real_new = Lab_real -( &
             Pr_Cummulative  & 
             )*iox/OM_uni
             
             Particle(ID)%OM%OM_Time0 = Time
             Particle(ID)%OM%OMact = Lab_real_New
     
     end if 
           
   End Function Lab_real_New

   !****************************************************************

   Function Activity_New(n)
     ! returns the new activity of particle n after decay

   Use GlobalVariables
   implicit none
   Integer n
   Real TimeIncrement, Activity_New

     TimeIncrement = Real (Time - Particle(n)%Pb%Time0) * TimeScale / 365.
     Activity_New = Particle(n)%Pb%Activity0 * (Exp( - DecayConstant * TimeIncrement ) )
     Particle(n)%Pb%Activity = Activity_New

   End Function Activity_New

   ! ************************************

   Function CurrentTimeString()
     ! returns the current time as a string character

   Use GlobalVariables
   implicit none
   Real :: YearValue, DayValue, HourValue, MinValue
   CHARACTER*3 Year_str, day_str
   Character*2 Hour, Minutes
   Character*14 CurrentTimeString

     YearValue = Real(Time) * TimeScale / 365.   ! no years elapsed
     DayValue = 365. * (YearValue - Real(INT(YearValue)))
     HourValue = 24. * (DayValue - Real(Int(DayValue)))
     MinValue = 60. * (HourValue - Real(Int(HourValue)))

     write(Year_str, '(I3.3)') Int(YearValue)
     write(day_str, '(I3.3)') Int(DayValue)
     write(Hour, '(I2.2)') Int(HourValue)
     write(Minutes, '(I2.2)') Int(MinValue)

     CurrentTimeString(1:1) = 'Y'
     CurrentTimeString(2:4) = Year_str
     CurrentTimeString(5:5) = 'D'
     CurrentTimeString(6:8) = day_str
     CurrentTimeString(9:9) = 'H'
     CurrentTimeString(10:11) = Hour
     CurrentTimeString(12:12) = 'M'
     CurrentTimeString(13:14) = Minutes

   End Function CurrentTimeString

  ! **********************

   Function P_Count()
     ! count the number of particles .. used to debug, not functional in the simulation

   Use GlobalVariables
   implicit none
   Logical IsParticle
   Integer  i, j, x, y, P_Count

     P_Count = 0

       Do Y = 1, N_row
       Do X = 1, N_Col
         If (IsParticle(CoOrdinates(Y,X))) P_count = P_count + 1
       End do
       End do

       Do i = 1, N_ind
       Do j = 1, Org(i)%Gut%Capacity
         if (Guts(j,i)%Class .eq. p)  P_Count = P_Count + 1
       End do
       End do

       If (P_Count .ne. Total_N_Particles) stop

   End Function P_Count

   ! ***************************
   
   function poly(x,a,n)
   
   implicit none
   
   integer :: n, i 
   double precision :: a(n)
   real x, clc, poly
   
   clc = 0.
   do i = 1, n
      clc = clc + a(i)*x**(i-1)
   end do 
   
   poly = clc
   
   end function poly
   
   
   
  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

!  ******************************************

  subroutine dot(a,c,b,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n), b(n)
! double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, k

do k=1, n
   coeff = 0.0
   do i=1,n
      coeff = coeff + a(k,i)*c(i)
   end do
   b(k) = coeff
end do

end subroutine dot
