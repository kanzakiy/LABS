 
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
   integer(kind=4)               :: i, DepthToScan, TimeNew, ParticlesToSediment
   Character*15          :: CurrentTimeString
   real(kind=8)                  :: Area_Period
   Type(CoOrdinates)     :: Point
   integer(kind=4)               :: xx, yy, j
   type(coordinates)     :: org_mid
   logical               :: rec_flow
   
   character*10 dumchr(3)
   integer(kind=4) dumint(8), dumPop, accpat
   real(kind=8) rnum 
   
   call date_and_time(dumchr(1),dumchr(2),dumchr(3),dumint)
   
   write(WorkDir,*) 'C:/Users/YK/Desktop/biot-res/'
   call getarg(1,WorkName)
   ! write(WorkName,*) 'Egest_test1'
   write(Today,*) trim(adjustl(WorkDir))//trim(adjustl(Workname))//'-'//trim(adjustl(dumchr(1)))
   
   call system ('mkdir -p '//trim(adjustl(Today)))
   call system ('mkdir -p '//trim(adjustl(Today))//'/data')
   
   call system ('cp labs.exe '//trim(adjustl(Today)))

       OPEN(unit = File_Profile,  file = trim(adjustl(Today))//          &
                               '/DepthProfiles.OUT', status = 'unknown')
       OPEN(unit = File_Profile2,  file = trim(adjustl(Today))//        &
                               '/DepthProfiles2.OUT', status = 'unknown')
       OPEN(unit = File_Profile3,  file = trim(adjustl(Today))//        &
                              '/DepthProfiles3.OUT', status = 'unknown')
       OPEN(unit = File_Profile_Re,  file = trim(adjustl(Today))//        &
                             '/DepthProfiles_Re.OUT', status = 'unknown')
       OPEN(unit = File_Scaling, file = trim(adjustl(Today))//        &
                              '/PorosityScale.OUT', status = 'unknown')
       OPEN(unit = File_Activity, file = trim(adjustl(Today))//        &
                               '/ActivitySlope.OUT', status = 'unknown')
       OPEN(unit = File_Displace, file = trim(adjustl(Today))//        &
                               '/MeanDisplacement.OUT', status = 'unknown')
       OPEN(unit = File_Displace2, file = trim(adjustl(Today))//        &
                               '/MeanDisplacement2.OUT', status = 'unknown')
       OPEN(unit = Poly_fit, file = trim(adjustl(Today))//                           &
                                  '/polyFit.OUT', status = 'unknown')
       OPEN(unit = File_flux, file = trim(adjustl(Today))//                           &
                                  '/flux.OUT', status = 'unknown')
       OPEN(unit = File_Log, file = trim(adjustl(Today))//                           &
                                  '/log.OUT', status = 'replace')
       OPEN(unit = File_flux2, file = trim(adjustl(Today))//                           &
                                  '/flux2.OUT', status = 'unknown')
       ! OPEN(unit = File_Diet, file = trim(adjustl(Today))//                           &
                                  ! '/Diet.OUT', status = 'unknown')
       OPEN(unit = File_Core, file = trim(adjustl(Today))//                           &
                                  '/Core.OUT', status = 'unknown')
       OPEN(unit = File_Core_M, file = trim(adjustl(Today))//                           &
                                  '/Core_M.txt', status = 'unknown')
       OPEN(unit = File_Core_L, file = trim(adjustl(Today))//                           &
                                  '/Core_L.txt', status = 'unknown')
       OPEN(unit = File_Core_A, file = trim(adjustl(Today))//                           &
                                  '/Core_A.txt', status = 'unknown')
       OPEN(unit = File_Pop, file = trim(adjustl(Today))//                           &
                                  '/Pop.OUT', status = 'unknown')
       OPEN(unit = File_sedrate, file = trim(adjustl(Today))//                           &
                                  '/sedrate.OUT', status = 'unknown')

   ! Begin main program loop
   
   oxygen_ON = .true.
   oxygen_ON = .false.
   
   oxFB_ON = .true.
   ! oxFB_ON = .false.
   
   Resp_ON = .true.
   ! resp_ON = .false.
   
   errCHk = .true.
   ! errChk = .false.
   
   I_shape = .true.
   I_shape = .false.
   
   flow_ON = .true.
   flow_ON = .false.

   only_sed = .True.   
   only_sed = .FALSE.
   
   Detail_Log = .True.
   Detail_Log = .False.

   Incl_ASH = .True. 
   Incl_ASH = .False. 
   
   O2ratelaw = 'linear'
   ! O2ratelaw = 'monod'
   ! O2ratelaw = 'zero'
	
   Org_ID_ishape = 1
   
   Mod_Pop = .True.
   Mod_Pop = .False.
   
   non_pop =  .true.
   non_pop = .false.
   
   trans_making = .true.
   ! trans_making = .false.
   
   Long_Run = .true.
   ! Long_Run = .false.



       CALL RANDOM_SEED    ! initialise random seed generator
       Call GetUserInput() ! user inputs
   
     call O2i_setup()
    call O2pre_setup()
    Call Output_O2txtImg()
    Call Output_txtImg()
    
    call gnuplot_flux(' ')
    call gnuplot_flux('2')
   
   Savetime = 1 
   
   DumPop = N_Ind
   accpat=0

       DO Time = 1, TimeMax

             Pr_Activity_year = 0.5 * (1. + cos(pi + 2.*pi*(Real(Time)/Real(Year)) ) )
             Pr_Activity_day = 0.5 * (1. + cos(pi + 2.*pi*(Real(Time)/Real(Day)) ) )
             CurrentTime = CurrentTimeString()
			 
			 print*,trim(adjustl(workname)),': ', currenttime
             
             
            call chk_org('tmtble')
			 
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
         call Make_matrix_chk('Tim')
         if (errDetect) then; print *,time, "err after org_move"; write(file_log,*)time, "err after org_move"; end if
         call Matrix_err_chk('Tim')
         if (errDetect) then; print *,time, "err after org_move"; write(file_log,*)time, "err after org_move"; end if 
         end if 
         
         ! YK: calculation of depth of seawater-sediment interface at each x
         do xx = 1, n_col
            do yy=1,n_row
                if (matrix(yy,xx)%class==w) cycle
                if (matrix(yy,xx)%class==p) then 
                    swi(xx)=yy
                    exit
                endif
            enddo
        enddo

       ! add new particles and sediment them down
           If (Mod(Time,Time_Sediment) .eq. 0) then
             If (Sedimentation_On) then
               ! determine number of particles to be input (the fraction of the area under the sinusoid) and do it
                   Area_Period = Real(Time_Sediment) * 0.5 *       &
                                 (1.+ cos(pi+ 2.*pi*(Real(Time - Time_Sediment*0.5)/Real(Year))))
                   ParticlesToSediment = ParticlesToSediment_y * ( Area_Period / Area_Total )
                   Call Particle_sediment(ParticlesToSediment)
                   accpat=accpat+ParticlesToSediment
                   write(file_sedrate,*) (time/real(year)),accpat*pixelsize/real(n_col)/(time/real(year))*1d3  ! cm/kyr
                   print*, (time/real(year)),accpat*pixelsize/real(n_col)/(time/real(year))*1d3
                   Call LocateOrganisms_scan()  ! subroutine that collects info of the location of organisms ..
             End if

             ! Scan the water column and sediment down any free particles
                 DepthToScan = N_RowWater
                 Call WaterColumn_scan(DepthToScan)

             ! sneek an output of the current time to show how far the sim has progressed
                 write(*,*) "sedimentation!!",CurrentTime, Time,ParticlesToSediment
                 
         call disperse()        ! subroutine to dispers in x direction
         
         
         !  check for error 
         ! if (errChk) then 
         call Make_matrix_chk('Sed')
         if (errDetect) then; print *,time, "err after sed"; write(file_log,*)time, "err after sed"; end if 
         call Matrix_err_chk('Sed')
         if (errDetect) then; print *,time, "err after sed"; write(file_log,*)time, "err after sed"; end if 
         ! end if 
                 

           End if

       ! Ashing  ew particles and sediment them down
           If (Incl_ASH) then
		     Time_Ash = 49.  ! days
			 Ash_ON = .false.
			 if (Time==int(Time_ASH/TimeScale)) then 
			 Ash_ON = .true.
			 Ash_porosity=porosity0
		   ! determine number of particles to be input (the fraction of the area under the sinusoid) and do it
			   Ash_thickness = 1.  ! cm
			   ParticlesToSediment = Ash_thickness/pixelsize * n_col*(1d0-ash_porosity)
			   Call Particle_sediment(ParticlesToSediment)
			   Call LocateOrganisms_scan()  ! subroutine that collects info of the location of organisms ..
			   write(*,*) "Ash!!!!",CurrentTime, Time, ParticlesToSediment
             End if

             ! Scan the water column and sediment down any free particles
                 DepthToScan = N_RowWater
                 Call WaterColumn_scan(DepthToScan)
         
         !  check for error 
         if (errChk) then 
         call Make_matrix_chk('Ash')
         if (errDetect) then; print *,time, "err after sed"; write(file_log,*)time, "err after sed"; end if 
         call Matrix_err_chk('Ash')
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
         call Make_matrix_chk('Bur')
         if (errDetect) then; print *,time, "err after constrain"; write(file_log,*)time, "err after constrain"; end if 
         call Matrix_err_chk('Bur')
         if (errDetect) then; print *,time, "err after constrain"; write(file_log,*)time, "err after constrain"; end if 
         end if 
             
           End if
		   
		   rec_flow = .false.
		   if (flow_ON) then 
		     if (any(abs(Vb)/=0).or.any(abs(Ub)/=0)) then 
               
               do i = 1,N_ind
                 do j = 1,2
                   if (matrix(flow_loc(j,i)%y,flow_loc(j,i)%x )%class /= i) then 
                     write(File_log,*) Time,"err:flow boundary at x,y=", flow_loc(j,i)%x &
                     ,   flow_loc(j,i)%y, 'class =',matrix(flow_loc(j,i)%y,flow_loc(j,i)%x )%class 
                   endif
                 enddo
               enddo
             
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
             
             if (errChk) then 
               do yy = 1, N_Row
                 do xx = 1, N_Col
                    if ((matrix(yy,xx)%class == p)) then
                      if ((Vo(xx,yy)==0.0).and.(Uo(xx,yy)==0.0)) then 
                        cycle
                      else 
                        print *, "particle has flow --- ERROR --- at (y,x) =", yy, xx &
                          , 'and Vo and Uo = ',Vo(xx,yy),Uo(xx,yy)
                        write(File_log, *) time, "particle has flow --- ERROR --- at (y,x) =", yy, xx  &
                          , 'and Vo and Uo = ',Vo(xx,yy),Uo(xx,yy)
                      endif 
                    end if 
                 end do
               end do
             endif

             
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
               write(File_log, *) "particle has oxygen --- permanent effect ---", yy, xx
             end if 
           end do
         end do
         
         end if
         
         If (.not. oxygen_ON) then
            o2(:,:)%oxygen =1.0
            call OrgDecay()
            call fluxes()
           O2%oxygen_use = 0.0
		 Endif
		 
		 !!   working only when oxygen is switched ------------ END
         
         if (errDetect .and. errChk) then; call Output_txtImg(); call Output_txtImg_chk(); endif
		 
       ! check if time to output data
           If (Any(Time .eq. Time_Output)) then
             TimeToOutput = .true.
             ! write(*,*) time      !! Added to create movies -------YK 5/22/2017
             ! call OutputData()    !! Added to create movies -------YK 5/22/2017
             TimeNew = Time
           End If
           
           ! if (time==25001) TimeToOutput = .true.

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
                 call Make_matrix_chk('Out')
                 if (errDetect) then; print *,time, "err when outputting"; write(file_log,*)time, "err when outputting"; end if  
                 call Matrix_err_chk('Out')
                 if (errDetect) then; print *,time, "err when outputting"; write(file_log,*)time, "err when outputting"; end if                         ! snapshots at user-defined times
               
               Call OutputData()
               If (oxygen_ON) Call Output_O2txtImg()
			   if (rec_flow) call output_flow()
               Call Output_txtImg()
               call Output_txtImg_chk()
               O2%oxygen_use = 0.0
          
          if (trans_making) then 
               call organism_dead(1,.true.)
               call make_trans()
               call addpop()
          endif
          
               Savetime = savetime + 1       
               TimeToOutput = .false.                  ! turn off output until the next output time
             End if
           End if
		   
           !  population modification 
           if (Mod_Pop) then 
           
           ! if (Time ==int(25./timescale) .or. Time ==int(50./timescale) .or. Time ==int(100./timescale)) then 
             ! call addpop()
           ! endif
           
           If (Mod(Time,Time_Sediment) .eq. 0) then
                 Call Random_Number(rnum)
                 print*,'just before addpop'
                 call addpop()
                 print*,'finished addpop'
                 call chk_org('addpop')
           End if
           
           ! if (Time ==int(325./timescale) ) then 
             ! do while(N_ind>0)  ! extinction 
             ! i = 1
             ! call organism_dead(i)
             ! enddo
           ! endif
           
           if (sum(DeathFlg)/=0) then 
               do while (sum(DeathFlg)/=0)
                 do i = 1,N_ind
                   if (DeathFlg(i)/=0) exit
                 enddo 
                 call organism_dead(i,.true.)  ! within mod_pop
               enddo
          endif
		 
         
          End if 
          
          ! testing no-org
          if (non_pop .and. time==1) then 
            call organism_dead(1,.true.)
            print*,'dead'
          endif
          
          ! if (trans_making) then 
            ! If (Mod(Time,Time_Sediment) .eq. 0) then
               ! call organism_dead(1)
               ! call make_trans()
               ! call addpop()
            ! endif
          ! endif
          
         
         if (Time==1 .or. DumPop/=N_Ind .or. Time==TimeMax) then 
            write(File_Pop,*) (Time-1)*Timescale/365.25, DumPop
            write(File_Pop,*) Time*Timescale/365.25, N_Ind
         endif
         
         DumPop = N_Ind
		 
		 Vo = 0.
		 Uo = 0.
		 flow_loc = coordinates(0,0)
		 Vb = 0.
		 Ub = 0.  
		   

       END Do    !  main time loop
       
       savetime = savetime-1
        call Gnuplot_Db(" ")
        call Gnuplot_Db("x")
        call Gnuplot_Db("y")

   ! end main program

   ! clean up data allocations

       Deallocate (Org, Matrix, Particle, Org_Loc, Guts, MovementHistory, IngestionHistory, RespHistory, Dir_rec)
       
       deallocate (O2)
       deallocate (particle2)
	   
	   deallocate (flow_loc, Ub, Vb)
	   deallocate (Ug,Vg, Pg, Dg) 
	   deallocate (Uo,Vo) 
	   deallocate (edif)
	   deallocate (RespCnst)
	   deallocate (DeathFlg)

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
       Close(File_flux2)
       ! Close(File_Diet)
       Close(File_Core)
       Close(File_Core_M)
       Close(File_Core_L)
       Close(File_Core_A)
       Close(File_Pop)
       Close(File_sedrate)
       
       do i=1,N_Ind
       close(1000+PopLogID(i))
       enddo
       
       ! do i =1,PopTot
         ! call gnuplot_diet(i)
       ! enddo

   End Program

   ! ********************************************************

   Subroutine GetUserInput()
     ! load data from files and allocate necessary memory
     ! called once by the main program, Automatons.f90

   Use GlobalVariables
   implicit none
   Character*1 :: FeedingType
   integer(kind=4)     :: i, j, n, MaxOrgSize, MaxGutCapacity, MinHeadWidth, Lability_tmp
   real(kind=8)        :: MatrixDepth, MatrixWidth, Depth_DeadZone, HalfLife
   real(kind=8)        :: r_TimeMin, r_TimeMax, r_TimeStep
   real(kind=8)        :: ParticleDensity, PArticleSize, Sum_Lability, MaxOrgSpeed, MinWidth, MaxWidth
   real(kind=8)        :: tmp_MaxOrgSize, tmp_OrgSize, tmp_Width, tmp_Length
   real(kind=8)        :: r_Width, r_Length, r_Speed
   real(kind=8)        :: r_IngestRate, r_IngestSelectivity, r_Density
   character*32 :: IDdum
   
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
         READ(File_Parameters,*) r_TimeMax          ! Endtime (days) changed to years
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
       ! N_RowWater = 3 * N_Row / 10
       N_RowWater = ceiling(3.5/pixelsize)
       N_RowSed = N_Row - N_RowWater
       DepthToConstrain = N_RowWater    ! this is the depth to which the totalnumber of cells are constrained

       DepthDeadZone = Ceiling(Depth_DeadZone / PixelSize)  ! the inverse .. used for calculations elsewhere

   ! the ** critical scaling parameter for time **
   ! this gives the length of real time (days) for each unit of simulation time (RealTime/SimulTime)
       TimeScale = PixelSize / MaxOrgSpeed       !  day
       
       
         print '("$$$$$$$$$$$ Time and space configuration $$$$$$$$$$$$$")'
         print'("PixelSize (cm/pixel): ",f6.4)',PixelSize
         print'("TimeScale (day/step): ",f6.4)',TimeScale
         print '("\\\\\\\\\ End of config. \\\\\\\\\\\\\")'
         print*,''

   ! define some useful time constants

         Day   =   1 / TimeScale + Int(2.*mod(1.,Timescale))
         Year  = 365 / TimeScale + Int(2.*mod(365.,Timescale))
		 
      ! N_Outputs = 15
      ! If (oxygen_ON) N_Outputs = 73*3
      N_Outputs = 10000
      N_Outputs2 = 100000

    Allocate (Time_output(N_Outputs))
    Allocate (Time_output2(N_Outputs2))
   
   ! if (.not. oxygen_ON) then 
     ! Time_Output =(/1,25*Day,50*Day,100*Day,150*Day,200*Day,250*Day,1*Year,2*Year,3*Year,5*Year,10*year,20*year,30*Year,100*year/)
	 ! end if 
	 
     do i = 1, N_Outputs2
        Time_Output2(i) = i
     end do
     
	 ! If (oxygen_ON) then
     do i = 1, N_Outputs
        Time_Output(i) = 5000*i
     end do
	 ! end if 
     
     if (Long_run) then
        ! Time_Output = Time_Output*10
        
         do i = 1, N_Outputs
            Time_Output(i) = r_TimeStep*year*i
         end do
     endif
     
     if (Detail_log) then
     deallocate(Time_output)
     allocate(Time_output(N_outputs2))
     time_output= time_output2
     endif

   ! rescale time to units of simulation
       TimeMin =  r_TimeMin / TimeScale + Int(2.*mod(r_TimeMin,Timescale))
       TimeMax =  (r_TimeMax*365 + 1) / TimeScale + Int(2.*mod(r_TimeMax,Timescale))
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
     MaxOrgSize = Ceiling(real(MaxOrgSize)*2.0*2.0)
!  YK  why adding 1? does not matter because it is for maximum?

     ALLOCATE (Matrix(N_Row, N_Col), MatrixOccupancy(N_Row))
     Allocate (Particle(N_Cell))
     Allocate (Particle2(N_Cell))
     Allocate (Org(N_Ind))
     Allocate (Org_Loc(MaxOrgSize, N_Ind))
     Allocate (IngestionHistory(N_Ind, Day), MovementHistory(N_Ind, Day))
     
     Allocate (RespHistory(N_Ind, Day),EgestHistory(N_Ind, Day))      
     Allocate (EnergyHistory(N_Ind, Day), CurrentEnergy(N_Ind))      
     allocate (O2(N_row,n_col))
     allocate (matrix_chk(N_row,n_col))
	 allocate (Dir_rec(Day,N_IND))
	 allocate (Ub(2,N_ind), Vb(2,N_ind)) !  in flow subroutines, coodinate is listed in order of x, y
	 allocate (flow_loc(2,N_ind)) !  in flow subroutines, coodinate is listed in order of x, y
	 allocate (Ug(N_col+1,N_row+2),Vg(N_col,N_row+2+1), Pg(N_col,N_row+2), Dg(N_col,N_row+2)) 
	 allocate (Uo(N_col,N_row),Vo(N_col,N_row))
	 allocate (edif(n_row, n_col)) 
	 allocate (DeathFlg(N_Ind)) 
	 allocate (RespCnst(N_Ind)) 
	 allocate (PopLogID(N_Ind)) 
     
     PopTot = 0

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
         EgestHistory(i,:)     = 0
         
         CurrentEnergy(i)   = 1.
         
         RespCnst(i) = kdcy*bio_fact   ! wt%-1 yr-1
         DeathFlg(i) = 0
         
         PopTot = PopTot + 1
         PopLogID(i) = PopTot
         
         Write(IDdum,'(i3.3)') PopLogID(i)
         
         OPEN(unit = 1000+PopLogID(i), file = trim(adjustl(Today))//                           &
                                  '/Diet-'//trim(adjustl(IDdum))//'.OUT', status = 'unknown')
         
         
         OPEN(unit = File_temp, file = trim(adjustl(Today))//                           &
                                  '/OrgInfo-'//trim(adjustl(IDdum))//'.txt', status = 'unknown')
         
         print '("=========== properties of organism #",i2," ================")',i
         print'("Density (g/cm3): ",f6.2)',Org(i)%Density
         print'("Length (pixels): ",i3.2)',Org(i)%Length
         print'("Width (pixels): ",i3.2)',Org(i)%Width
         print'("HeadSize (pixels): ",i3.2)',Org(i)%HeadSize
         print'("BodySize (pixels): ",i3.2)',Org(i)%BodySize
         print'("MoveRate (pixels/day): ",f6.2)',Org(i)%Move%Rate/TimeScale
         print'("IngestRate (particles/day): ",f6.2)',Org(i)%Ingest%Rate/TimeScale
         print '("~~~~~~~~~~~ End of properties of organism #",i2," ~~~~~~~~~~~~")',i
         print*,''
         
         write(FIle_temp, '("=========== properties of organism #",i2," ================")')i
         write(FIle_temp,'("Density (g/cm3): ",f6.2)')Org(i)%Density
         write(FIle_temp,'("Length (pixels): ",i3.2)')Org(i)%Length
         write(FIle_temp,'("Width (pixels): ",i3.2)')Org(i)%Width
         write(FIle_temp,'("HeadSize (pixels): ",i3.2)')Org(i)%HeadSize
         write(FIle_temp,'("BodySize (pixels): ",i3.2)')Org(i)%BodySize
         write(FIle_temp,'("MoveRate (pixels/day): ",f6.2)')Org(i)%Move%Rate/TimeScale
         write(FIle_temp,'("IngestRate (particles/day): ",f6.2)')Org(i)%Ingest%Rate/TimeScale
         write(FIle_temp, '("~~~~~~~~~~~ End of properties of organism #",i2," ~~~~~~~~~~~~")')i
         write(FIle_temp,*)''
         
         
         Close(File_temp)
         
         call gnuplot_diet(i)

       End do
       
       ! stop

       CLOSE(File_Organisms)

       MatrixOccupancy = 0   ! initial the frequency of occurence counts of organisms in the matrix

       MaxGutCapacity = 1
       do i=1, n_ind
         if (MaxGutCapacity .lt. Org(i)%Gut%Capacity) MaxGutCapacity = Org(i)%Gut%Capacity
       End do
       
       MaxGutCapacity = ceiling(MaxGutCapacity*2.0*2.0)

       DecayConstant = -log(0.5) / HalfLife   ! the decay constant for 210Pb   !  /yr

   ! transition probabilities of lability change (1) == 1 to 0, (2) == 2 to 1, etc. .. as halflife in year**(-1)
       N_LabilityClasses = 10
       Allocate (Lability_decayConstant(N_LabilityClasses), Lability_proportion(N_LabilityClasses))
       
       do j = 1, N_labilityClasses
         lability_decayConstant(j) = kdcy   ! /yr 
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
         
       allocate(swi(n_col))
       swi=0

   End Subroutine GetUserInput


   ! ********************************************************************

   Subroutine AddPop()
     ! increasing organism population (11/7/2018 YK added)

   Use GlobalVariables
   implicit none
   Character*1 :: FeedingType
   integer(kind=4)     :: i, j, n, MaxOrgSize, MaxGutCapacity, MinHeadWidth, Lability_tmp
   real(kind=8)        :: MatrixDepth, MatrixWidth, Depth_DeadZone, HalfLife
   real(kind=8)        :: r_TimeMin, r_TimeMax, r_TimeStep
   real(kind=8)        :: ParticleDensity, PArticleSize, Sum_Lability, MaxOrgSpeed, MinWidth, MaxWidth
   real(kind=8)        :: tmp_MaxOrgSize, tmp_OrgSize, tmp_Width, tmp_Length
   real(kind=8)        :: r_Width, r_Length, r_Speed
   real(kind=8)        :: r_IngestRate, r_IngestSelectivity, r_Density
   
   Type(Organism_Char),allocatable  :: Org_dum(:)
   Type(Coordinates), Allocatable :: Org_loc_dum(:,:)
   Type(Coordinates) :: flow_loc_dum(2,N_Ind)
   Type(CellContent), Allocatable :: Guts_dum(:,:)
   real(kind=8)        ::  Energy_dum(N_ind),  V_dum(2,N_IND), Resp_dum(N_Ind)
   integer(kind=4),allocatable     ::  History_dum(:,:),Dir_dum(:,:)
   integer(kind=4)     :: flg_dum(N_Ind), ID_dum(N_ind)
   character*36 :: name_dum
   real(kind=8) :: rnum, rfact 
   
   
   Call Random_Number(rnum)
   
   rfact = rnum + 0.5d0
   if (trans_making) rfact = 1d0

   MaxOrgSize = size(Org_loc(:,1))
   MaxGutCapacity = size(Guts(:,1))
   
   
   if (allocated(Org_loc_dum)) deallocate(Org_loc_dum)  
   if (allocated(Guts_dum)) deallocate(Guts_dum)
   if (allocated(Org_dum)) deallocate(Org_dum)
   if (allocated(Dir_dum)) deallocate(Dir_dum)
   if (allocated(History_dum)) deallocate(History_dum)
   
   allocate(Org_loc_dum(MaxOrgSize, N_Ind))   
   allocate (Guts_dum(MaxGutCapacity, N_ind))
   allocate (Org_dum(N_ind))
   allocate (History_dum(N_ind,day))
   allocate (Dir_dum(Day,N_ind))
   
   N_Ind = N_Ind + 1
   
   ! making copies, deallocate and save existing animals 
   Org_dum = Org 
   deallocate(Org)
   allocate(Org(N_ind))
   Org(1:N_ind-1) = Org_dum(:)
   
   Org_loc_dum = Org_loc 
   deallocate(Org_loc)
   Allocate (Org_Loc(MaxOrgSize, N_Ind))
   Org_loc(:,1:N_Ind-1) = Org_loc_dum(:,:)
   
   Guts_dum = guts
   deallocate(guts)
   allocate(guts(Maxgutcapacity,N_ind))
   guts(:,1:N_ind-1)=guts_dum
   
   History_dum = IngestionHistory 
   deallocate(IngestionHistory) 
   Allocate (IngestionHistory(N_Ind, Day))
   IngestionHistory(1:N_ind-1,:) = History_dum
   
   History_dum = MovementHistory
   deallocate(MovementHistory)
   allocate(MovementHistory(N_Ind,Day))
   MovementHistory(1:N_ind-1,:)=History_dum
   
   History_dum=RespHistory
   deallocate(RespHistory)
   allocate(RespHistory(N_Ind,Day))
   RespHistory(1:N_Ind-1,:) = History_dum
   
   History_dum = EgestHistory
   deallocate(EgestHistory)
   allocate(EgestHistory(N_ind,Day))
   EgestHistory(1:N_ind-1,:)=History_dum
   
   History_dum=EnergyHistory
   deallocate(EnergyHistory)
   allocate(EnergyHistory(N_ind,Day))
   EnergyHistory(1:N_ind-1,:)=History_dum
   
   Energy_dum = CurrentEnergy
   deallocate(CurrentEnergy)
   allocate(CurrentEnergy(N_ind))
   CurrentEnergy(1:N_ind-1)=Energy_dum
   
   Dir_dum = Dir_rec 
   deallocate(Dir_rec)
   allocate(Dir_rec(Day,N_Ind))
   Dir_rec(:,1:N_ind-1)=Dir_dum
   
   V_dum = Ub
   deallocate(Ub)
   allocate(Ub(2,N_ind))
   Ub(:,1:N_Ind-1)=V_dum
   
   V_dum = Vb
   deallocate(Vb)
   allocate(Vb(2,N_ind))
   Vb(:,1:N_ind-1)=V_dum
   
   flow_loc_dum = flow_loc
   deallocate(flow_loc)
   allocate(flow_loc(2,N_ind))
   flow_loc(:,1:N_ind-1)=flow_loc_dum
   
   Resp_dum = RespCnst
   deallocate(RespCnst)
   allocate(RespCnst(N_ind))
   RespCnst(1:N_ind-1)=Resp_dum
   
   flg_dum = DeathFlg
   deallocate(DeathFlg)
   allocate(DeathFlg(N_ind))
   DeathFlg(1:N_ind-1)=flg_dum
   
   ID_dum = PopLogID
   deallocate(PopLogID)
   allocate(PopLogID(N_ind))
   PopLogID(1:N_ind-1)=ID_dum

   ! Read in again organism info and rescale organism traits to simulation time/space units
       OPEN(unit = File_Organisms, file = 'Organisms.IN', status = 'old')

         READ(File_Organisms,*) i
         Read(File_Organisms,*) r_Width
         Read(File_Organisms,*) r_Length
         Read(File_Organisms,*) r_Density
         Read(File_Organisms,*) FeedingType
         Read(File_Organisms,*) r_Speed
         Read(File_Organisms,*) r_IngestRate
         Read(File_Organisms,*) r_IngestSelectivity
         
         i = N_Ind
         
         ParticleDensity = 2.5 

         Org(i)%Density                = r_Density
         Org(i)%Pr_GoStraight          = 0.75
         Org(i)%FeedingType            = FeedingType
         Org(i)%Length                 = max(Ceiling(((r_Length + PixelSize) / PixelSize)*rfact), 5)
         Org(i)%Width                  = max(Ceiling(((r_Width + PixelSize) / PixelSize)*rfact), 3)
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

         Org(i)%Gut%Capacity           = min(ceiling(real(Org(i)%BodySize) / 4.),MaxGutCapacity)
         Org(i)%Gut%Content            = 0
         Org(i)%Gut%Pb_Activity        = 0

         Org(i)%Move%Rate              = r_Speed * TimeScale / PixelSize*rfact  ! number of pixels per unit simulation time
         Org(i)%Move%Distance          = 0
         Org(i)%Move%Bearing_Distance  = 0
         Org(i)%Move%Bearing           = 0.0
         Org(i)%Move%ReverseCount      = 0

         MovementHistory(i,:)          = 0

         Org(i)%Move%StoppedTime       = 0

         Org(i)%Ingest%Rate        = r_IngestRate*Org(i)%Bodysize*Org(i)%Density*Timescale/ParticleDensity*rfact ! n particles per day
         Org(i)%Ingest%Selectivity = r_IngestSelectivity
         Org(i)%Ingest%Amount      = 0

         IngestionHistory(i,:)     = 0
         
         respHistory(i,:)     = 0
         EgestHistory(i,:)     = 0
         
         RespCnst(i) = kdcy*bio_fact/rfact   ! wt-1 yr-1
         DeathFlg(i) = 0
         
         PopTot = PopTot + 1
         PopLogID(i) = PopTot
         
         Write(name_dum,'(i3.3)') PopLogID(i)
         
         OPEN(unit = 1000+PopLogID(i), file = trim(adjustl(Today))//                           &
                                  '/Diet-'//trim(adjustl(name_dum))//'.OUT', status = 'unknown')
         
         
         OPEN(unit = File_temp, file = trim(adjustl(Today))//                           &
                                  '/OrgInfo-'//trim(adjustl(name_dum))//'.txt', status = 'unknown')
         
         print '("=========== properties of organism #",i2," ================")',i
         print'("Density (g/cm3): ",f6.2)',Org(i)%Density
         print'("Length (pixels): ",i3.2)',Org(i)%Length
         print'("Width (pixels): ",i3.2)',Org(i)%Width
         print'("HeadSize (pixels): ",i3.2)',Org(i)%HeadSize
         print'("BodySize (pixels): ",i3.2)',Org(i)%BodySize
         print'("MoveRate (pixels/day): ",f6.2)',Org(i)%Move%Rate/TimeScale
         print'("IngestRate (particles/day): ",f6.2)',Org(i)%Ingest%Rate/TimeScale
         print '("~~~~~~~~~~~ End of properties of organism #",i2," ~~~~~~~~~~~~")',i
         print*,''
         
         write(FIle_temp, '("=========== properties of organism #",i2," ================")')i
         write(FIle_temp,'("Density (g/cm3): ",f6.2)')Org(i)%Density
         write(FIle_temp,'("Length (pixels): ",i3.2)')Org(i)%Length
         write(FIle_temp,'("Width (pixels): ",i3.2)')Org(i)%Width
         write(FIle_temp,'("HeadSize (pixels): ",i3.2)')Org(i)%HeadSize
         write(FIle_temp,'("BodySize (pixels): ",i3.2)')Org(i)%BodySize
         write(FIle_temp,'("MoveRate (pixels/day): ",f6.2)')Org(i)%Move%Rate/TimeScale
         write(FIle_temp,'("IngestRate (particles/day): ",f6.2)')Org(i)%Ingest%Rate/TimeScale
         write(FIle_temp, '("~~~~~~~~~~~ End of properties of organism #",i2," ~~~~~~~~~~~~")')i
         write(FIle_temp,*)''
         
         
         Close(File_temp)
         
         call gnuplot_diet(i)

       
       ! stop

       CLOSE(File_Organisms)

       Guts(:,i) = CellContent(w,w)    ! initialise the guts
       Org_Loc(:,i) = Coordinates(0,0)
	   
       print*,'before populating'
       
       Call Matrix_populate(i)
       
		 flow_loc(:,i) = coordinates(0,0)
		 Vb(:,i) = 0.0
		 Ub(:,i) = 0.0
         
         CurrentEnergy(i) = 1.

   End Subroutine AddPop


   ! ********************************************************************

   Subroutine KillPop(i)
     ! decreasing organism population; i is going to be killed (11/7/2018 YK added)

   Use GlobalVariables
   implicit none
   integer(kind=4),Intent(in) :: i 
   Character*1 :: FeedingType
   integer(kind=4)     :: j, n, MaxOrgSize, MaxGutCapacity, MinHeadWidth, Lability_tmp
   real(kind=8)        :: MatrixDepth, MatrixWidth, Depth_DeadZone, HalfLife
   real(kind=8)        :: r_TimeMin, r_TimeMax, r_TimeStep
   real(kind=8)        :: ParticleDensity, PArticleSize, Sum_Lability, MaxOrgSpeed, MinWidth, MaxWidth
   real(kind=8)        :: tmp_MaxOrgSize, tmp_OrgSize, tmp_Width, tmp_Length
   real(kind=8)        :: r_Width, r_Length, r_Speed
   real(kind=8)        :: r_IngestRate, r_IngestSelectivity, r_Density
   
   Type(Organism_Char),allocatable  :: Org_dum(:)
   Type(Coordinates), Allocatable :: Org_loc_dum(:,:)
   Type(Coordinates) :: flow_loc_dum(2,N_Ind)
   Type(CellContent), Allocatable :: Guts_dum(:,:)
   integer(kind=4)     :: Dir_dum(Day,N_IND), History_dum(5,N_Ind,day), Flg_dum(N_Ind), ID_dum(N_ind)
   real(kind=8)        :: Energy_dum(N_ind),  V_dum(2,2,N_IND), Resp_dum(N_Ind) 
   integer(kind=4)     :: cnt1,cnt2, xx, yy

   close(1000+PopLogID(i))
   
   MaxOrgSize = size(Org_loc(:,1))
   MaxGutCapacity = size(Guts(:,1))
   
   
   if (allocated(Org_loc_dum)) deallocate(Org_loc_dum)  
   if (allocated(Guts_dum)) deallocate(Guts_dum)
   if (allocated(Org_dum)) deallocate(Org_dum)
   
   allocate(Org_loc_dum(MaxOrgSize, N_Ind))   
   allocate (Guts_dum(MaxGutCapacity, N_ind))
   allocate (Org_dum(N_ind))
   
   N_Ind = N_Ind - 1
   
   ! making copies, deallocate and save existing animals 
   Org_dum = Org 
   deallocate(Org)
   allocate(Org(N_ind))
   
   Org_loc_dum = Org_loc 
   deallocate(Org_loc)
   Allocate (Org_Loc(MaxOrgSize, N_Ind))
   
   Guts_dum = guts
   deallocate(guts)
   allocate(guts(Maxgutcapacity,N_ind))
   
   History_dum(1,:,:) = IngestionHistory 
   deallocate(IngestionHistory) 
   Allocate (IngestionHistory(N_Ind, Day))
   
   History_dum(2,:,:) = MovementHistory
   deallocate(MovementHistory)
   allocate(MovementHistory(N_Ind,Day))
   
   History_dum(3,:,:)=RespHistory
   deallocate(RespHistory)
   allocate(RespHistory(N_Ind,Day))
   
   History_dum(4,:,:) = EgestHistory
   deallocate(EgestHistory)
   allocate(EgestHistory(N_ind,Day))
   
   History_dum(5,:,:)=EnergyHistory
   deallocate(EnergyHistory)
   allocate(EnergyHistory(N_ind,Day))
   
   Energy_dum = CurrentEnergy
   deallocate(CurrentEnergy)
   allocate(CurrentEnergy(N_ind))
   
   Dir_dum = Dir_rec 
   deallocate(Dir_rec)
   allocate(Dir_rec(Day,N_Ind))
   
   V_dum(1,:,:) = Ub
   deallocate(Ub)
   allocate(Ub(2,N_ind))
   
   V_dum(2,:,:) = Vb
   deallocate(Vb)
   allocate(Vb(2,N_ind))
   
   flow_loc_dum = flow_loc
   deallocate(flow_loc)
   allocate(flow_loc(2,N_ind))
   
   Resp_dum = RespCnst
   deallocate(RespCnst)
   allocate(RespCnst(N_ind))
   
   Flg_dum = DeathFlg
   deallocate(DeathFlg)
   allocate(DeathFlg(N_ind))
   
   ID_dum = PopLogID
   deallocate(PopLogID)
   allocate(PopLogID(N_ind))
   
   
   cnt1=1
   cnt2=1
   do j=1,N_ind
   if (j==i) then 
      cnt2=cnt2+1
   endif
   Org(cnt1) = Org_dum(cnt2)
   Org_loc(:,cnt1) = Org_loc_dum(:,cnt2)
   guts(:,cnt1)=guts_dum(:,cnt2)
   IngestionHistory(cnt1,:) = History_dum(1,cnt2,:)
   MovementHistory(cnt1,:)=History_dum(2,cnt2,:)
   RespHistory(cnt1,:) = History_dum(3,cnt2,:)
   EgestHistory(cnt1,:)=History_dum(4,cnt2,:)
   EnergyHistory(cnt1,:)=History_dum(5,cnt2,:)
   CurrentEnergy(cnt1)=Energy_dum(cnt2)
   Dir_rec(:,cnt1)=Dir_dum(:,cnt2)
   Ub(:,cnt1)=V_dum(1,:,cnt2)
   Vb(:,cnt1)=V_dum(2,:,cnt2)
   flow_loc(:,cnt1)=flow_loc_dum(:,cnt2)
   RespCnst(cnt1) = Resp_dum(cnt2)
   DeathFlg(cnt1) = Flg_dum(cnt2)
   PopLogID(cnt1) = ID_dum(cnt2)
   ! if (j/=N_ind+1) then 
   do n = 1, size(Org_Loc(:,cnt1))
       xx = Org_Loc(n,cnt1)%X
       yy = Org_Loc(n,cnt1)%Y
       if (xx==0 .or. yy==0) cycle 
       if (Matrix(yy,xx)%class==cnt2) then 
            Matrix(yy,xx)  = CellContent(n,cnt1) 
       endif            
   enddo
   ! Endif
   cnt1=cnt1+1
   cnt2=cnt2+1
   enddo
   
   if (cnt1/=N_ind+1) print*,'Error in KillPop'
   if (cnt1/=N_ind+1) write(File_Log,*) Time,'Error in KillPop'
   
   ! if (j/=N_ind) then 
   ! do i = 1, size(Org_Loc(:,N_ind+1))
       ! yy = Org_Loc(i,N_ind+1)%Y
       ! xx = Org_Loc(i,N_ind+1)%X
       ! if (xx==0 .or. yy==0) cycle 
       ! if (Matrix(yy,xx)%class==N_ind+1) then 
       ! Matrix(yy,xx)  = CellContent(w,w)  
       ! Endif
   ! enddo
   ! endif

   End Subroutine KillPop


   ! ********************************************************************


   Subroutine Timetable(i)
     ! process information of the acitivity of organism i

   Use GlobalVariables
   implicit none
   integer(kind=4), Intent(in):: i
   real(kind=8) IngestRAte, MovementRate, URand(3), Fullness, Pr_Activity
   integer(kind=4) :: j
   Character*21 numtemp
   
   real(kind=8) :: RespRate_ave
   integer(kind=4) :: k
   real(kind=8) :: ave_OM, ptcl_num
   
   real(kind=8) :: EgestRate, Energy_ave 
   
   ! logical :: IsFloating
   
   ! shift elements right
       MovementHistory(i,:) = EOShift (MovementHistory(i,:), -1)
       IngestionHistory(i,:) = EOShift (IngestionHistory(i,:), -1)
       RespHistory(i,:) = EOShift (RespHistory(i,:), -1)
       EgestHistory(i,:) = EOShift (EgestHistory(i,:), -1)
       EnergyHistory(i,:) = EOShift (EnergyHistory(i,:), -1)
       
       EnergyHistory(i,1) = CurrentEnergy(i)
       
   ! calculate moving averages of the rates
       MovementRate = Real(Sum(MovementHistory(i,:))) / Day
       IngestRate = Real(Sum(IngestionHistory(i,:))) / Day  !  n particle per single iteration 
       EgestRate = Real(Sum(EgestHistory(i,:))) / Day  !  n particle per single iteration 
       RespRate_ave = Real(Sum(RespHistory(i,:))) / Day
       Energy_ave = Real(Sum(EnergyHistory(i,:))) / Day
	   
       write(1000+PopLogID(i),*) Time*Timescale,i,MovementRate, IngestRate, EgestRate, RespRate_ave,Org(i)%Gut%Content &
                        ,Org(i)%Move%Rate,Org(i)%Ingest%Rate,Org(i)%Gut%Capacity, currentEnergy(i)
       
   ! check ingestion rates
       If (MovementRate .eq. 0) then               ! not moving
         If (Org(i)%Move%StoppedTime .eq. 0) Org(i)%Move%StoppedTime = Time  ! check stomach (i.e., caloric) intake is within tolerance

         if((Time - Org(i)%Move%StoppedTime) .gt. (Year*.25)) then ! if not moving for a long time  .. let it die
 
! this is a fallback routine to catch any anomalies
           if (.not.Mod_Pop) Call Organism_dead(i,.false.)  ! false death, i.e., regeneration  
           DeathFlg(i) = 1 ! raise flg to die
! dead .. remove and start a new organism
           Org(i)%Move%StoppedTime = 0                              ! reset the counter
           ! Write(File_Pop,*) Time*Timescale
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
           
           CurrentEnergy(i) = 0.
           do j = 1, Org(i)%Gut%Capacity         ! set location of particles ..
             IF (Guts(j,i)%Class .eq. p) CurrentEnergy(i) = CurrentEnergy(i) + Particle(Guts(j,i)%Value)%OM%OMact 
           End do
           
           CurrentEnergy(i) = CurrentEnergy(i)/real(Org(i)%Gut%Capacity) 
           
           If (Mod_Pop) Pr_Activity = 0.05+ CurrentEnergy(i)
		   
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
				   
                   ! IF ( ave_OM < 0.1) then     ! if OM is low in guts
				     ! Org(i)%Can%Egest = .true.       
                     ! if (fullness < 1.) Org(i)%Can%Ingest = .true.
                   ! end if 

                   if (fullness < 1.) Org(i)%Can%Ingest = .true.
                   
                   ! if (currentenergy > 1.) Org(i)%Can%Ingest = .true.                   
                   ! if (currentenergy > 1.) Org(i)%Can%Ingest = .false.                   
				 end if 
				 
                 if (.not.resp_ON) then 
                   IF (URand(2) .lt. Fullness) Org(i)%Can%Egest = .true.         ! too full .. defecate
                   If (Urand(3) .gt. Fullness) Org(i)%Can%Ingest = .true.        ! too empty .. eat
				 end if 
               End if
               
               
               If (EgestRate .ge. IngestRate) Org(i)%Can%Egest = .true.    !! this make sure a low fullness and efficient ingestion/egestion cycles
               ! If (EgestRate .lt. Org(i)%Ingest%Rate) Org(i)%Can%Egest = .true.    !! this make sure a low fullness and efficient ingestion/egestion cycles
               
               ! If (I_shape .and. i == Org_ID_ishape) Org(i)%Can%Egest = .false.   ! let worm floating around
               
               ! if (currentenergy > 1.) Org(i)%Can%Egest = .True.  
               
               
             End if
             End if
			 
   End Subroutine Timetable


   ! ******************************************************************


   SUBROUTINE Rules_of_Motion(i)
     ! the rules of motion for deposit feeder i ..

   Use GlobalVariables
   Implicit None
   integer(kind=4), Intent(In) :: i
   integer(kind=4) :: Possible(4), Direction(4), DirectionToGo, j
   real(kind=8)    :: RandomDir(4), ToGo(4), PREFERABLE(4), ToAvoid(4), Bearings(4), RandomWeight, URand(4)

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
   integer(kind=4) i
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
   integer(kind=4)   :: i, j, N_Particle, N_Organism, N_Water
   integer(kind=4), Intent(In) :: direction(4)
   integer(kind=4), Intent(Out) :: Possible(4)
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
   integer(kind=4), Intent(In) :: i, Possible(4)
   real(kind=8), Intent(Out)   :: ToAvoid(4)
   Type(CoOrdinates)   ::  Head_MinVal, Head_MaxVal
   real(kind=8)                :: Porosity_local
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
   integer(kind=4)     :: j, minVal, DirectionToPoint
   integer(kind=4), Intent(In)  :: i, Direction(4)
   real(kind=8), Intent(Out) :: Bearings(4)
   real(kind=8)  :: URand(2)
   
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
   integer(kind=4)             :: j, Y, X, index
   integer(kind=4)             :: N_Water(4), N_Particles(4), N_Organism(4)
   integer(kind=4)             :: Lability_sum(4)
   integer(kind=4), Intent(In) :: i, Direction(4)
   real(kind=8), Intent(Out)   :: PREFERABLE(4)
   real(kind=8), Intent(In)    :: ToGo(4)
   real(kind=8)                :: Lability_mean(4), O2_sum(4), O2_mean(4)
   Type(Coordinates)   :: Head_MinVal, Head_MaxVal
   
   real(kind=8)                ::  Fullness    
   real(kind=8)                ::  water_mean(4)   

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
     ! if (oxFB_ON) Preferable = Lability_mean*O2_mean    !!  if considering oxygen_feedback on organism movement
     if (oxFB_ON) Preferable = Lability_mean*O2_mean*O2_mean    !!  stronger O2 feedback
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
   integer(kind=4), intent(in) :: Y, X, j
   integer(kind=4), intent(inout) :: N_Particles(4), N_Water(4), N_Organism(4), Lability_sum(4)
   real(kind=8), intent(inout) :: O2_sum(4)  
   real(kind=8)     :: Lab_real_New        
   
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
   integer(kind=4), Intent(In) :: i, DirectiontoGo
   
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
   integer(kind=4) :: j
   integer(kind=4), Intent(In) :: i, Angle
   Type(CellContent) :: Image_ID(Org(i)%Width * Org(i)%Width)
   Type(Coordinates) :: Image(Org(i)%Width * Org(i)%Width)
   Type(Coordinates) :: Head_MinVal, Head_MaxVal
   
   integer(kind=4) :: xnew, ynew, xnew_2, ynew_2
   
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
   integer(kind=4) i, Orientation
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
   integer(kind=4) i, X, Y, k, Orientation
   integer(kind=4) N_Particle, N_Organism, N_Water
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
   integer(kind=4) X, Y
   integer(kind=4) N_Particle, N_Organism, N_Water
   
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
   integer(kind=4) ::  i, j, dX, dY
   Type(Coordinates) :: Image(Org(i)%BodySize)
   Type(Coordinates) :: Image_O2(Org(i)%BodySize)
   
   ! if (CurrentEnergy < 1.) return
   
   
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
           if (Image(j)%X+dX> n_col) then 
           Org_Loc(j,i) = Coordinates(Image(j)%Y+dY, Image(j)%X+dX-n_col)
           else if (Image(j)%X+dX < 1) then
           Org_Loc(j,i) = Coordinates(Image(j)%Y+dY, Image(j)%X+dX+N_col)
           else 
           Org_Loc(j,i) = Coordinates(Image(j)%Y+dY, Image(j)%X+dX)
           endif
           ! print *, Org_Loc(j,i)%Y, Org_Loc(j,i)%X
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
               
       ! EnergyHistory(i,1) = EnergyHistory(i,1) - 0.01   ! increment the running average
       ! CurrentEnergy = currentEnergy - 0.01
	   
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
         call Make_matrix_chk('Mov')
         if (errDetect) then; print *,time, "err after org_move(sub)"; write(file_log,*)time, "err after org_move(sub)"; end if 
         call Matrix_err_chk('Mov')
         if (errDetect) then; print *,time, "err after org_move(sub)"; write(file_log,*)time, "err after org_move(sub)"; end if 
		 if ((any(Org_loc%X<0)) .or. (any(Org_loc%Y <0)) .or. (any(Org_loc%Y > N_row)) .or. (any(Org_loc%Y > N_row))) then
		   print *, 'errrrrrrr in org_loc'
		 end if 
         end if 
     
   End Subroutine Organism_move

! ****************************************************************

   Subroutine LookupDirection(MoveOrient, dy, dX)
     ! Covert direction to move into translational vectors (dY,dX)

   integer(kind=4), Intent (in)  :: MoveOrient
   integer(kind=4), Intent (out) :: dy, dx

   
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
   integer(kind=4) :: j, X_new, DirectionToGo, dX
   integer(kind=4), Intent(IN) :: i
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
         call Make_matrix_chk('Rel')
         if (errDetect) then; print *,time, "err after relocate"; write(file_log,*)time, "err after relocate"; end if
         call Matrix_err_chk('Rel')
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
   integer(kind=4) :: j, neck, dx, dy, Reversedirection
   integer(kind=4), Intent(IN) :: i
   Type(Coordinates) :: Image(Org(i)%BodySize), Point
   
   Type(Coordinates) :: Image_2(Org(i)%BodySize)    
   integer(kind=4)           :: RevDir           
   integer(kind=4) :: xx, yy, k, o, oo , xmx, xmn, ymx, ymn
   integer(kind=4) :: DirChk(4), PreDir, PreOri, rec_image(Org(i)%length-Org(i)%width)
   
   real(kind=8) :: V_tmp(2), U_tmp(2)
   
   
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
         call Make_matrix_chk('REv')
         if (errDetect) then; print *,time, "err after reverse"; write(file_log,*)time, "err after reverse"; end if 
         call Matrix_err_chk('Rev')
         if (errDetect) then; print *,time, "err after reverse"; write(file_log,*)time, "err after reverse"; end if 
		 if ((any(Org_loc%X<0)) .or. (any(Org_loc%Y <0)) .or. (any(Org_loc%Y > N_row)) .or. (any(Org_loc%Y > N_row))) then
		   print *, 'errrrrrrr in org_loc'
		 end if 
         end if 
       
   End subroutine Organism_reverse

   ! ***************

   Subroutine Organism_dead(i, kill)
     ! remove organism from matrix, releasing the particles in its guts and re-initialise at the water interface

   Use GlobalVariables
   Implicit None
   Logical :: IsWater, InMatrix
   integer(kind=4) :: j, k
   integer(kind=4), Intent(IN) :: i
   real(kind=8) :: URand
   logical, intent(IN) :: kill
   
   integer(kind=4) :: xx, yy
   
     Do j = 1, Org(i)%BodySize
       Matrix(Org_loc(j,i)%Y,Org_loc(j,i)%X) = CellContent(w,w)
     End do

     Do j = 1, Org(i)%Gut%Capacity
       If (Guts(j,i)%Class .eq. p) then
         do
             Call Random_Number(URand)
             k = URand * Org(i)%BodySize + 1  ! randomly choosing a location within previous body space
             if (.not.inmatrix(Org_loc(k,i))) then 
                print*,'error: in Organism_dead'
                write(file_log,*)'error: in Organism_dead'
             endif
             If (IsWater(Org_loc(k,i))) then  
                 MAtrix(Org_loc(k,i)%Y, Org_loc(k,i)%X) = Guts(j,i)
				 particle(Guts(j,i)%value)%loc = Org_loc(k,i)   !! YK 
				 particle2(Guts(j,i)%value)%loc = Org_loc(k,i)   !! YK 
                 if (oxygen_ON) O2(Org_loc(k,i)%Y, Org_loc(k,i)%X)%oxygen = 0.0
                 Guts(j,i) = CellContent(w,w)
                 Exit
             End if
         End do
       End if
     End do
     
     ! If (.not. non_pop) then 
        
        ! If (Mod_pop .or. trans_making) then 
           ! Call KillPop(i)
        ! else 
            ! call Matrix_populate(i)
        ! endif
     ! Else if (non_pop) then
        ! Call KillPop(i)
     ! Endif
     
     if (kill) then 
       call killpop(i)
     else 
       call matrix_populate(i)
     endif
         
         if (errChk) then 
         call Make_matrix_chk('Ded')
         if (errDetect) then; print *,time, "err after dead"; write(file_log,*)time, "err after dead"; end if 
         call Matrix_err_chk('Ded')
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
   integer(kind=4) :: j, k, dx, dy, X_new
   integer(kind=4), Intent(IN) :: i
   Type(Coordinates) :: Image(Org(i)%BodySize)
   
   
   ! if (CurrentEnergy < 1.) return

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
               
     ! EnergyHistory(i,1) = EnergyHistory(i,1) - 0.01   ! increment the running average
     ! CurrentEnergy = currentEnergy - 0.01
	 
	 
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
         call Make_matrix_chk('Wra')
         if (errDetect) then; print *,time, "err after wrap"; write(file_log,*)time, "err after wrap"; end if 
         call Matrix_err_chk('Wra')
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
   integer(kind=4) ::  X, Y, k, Particle_ID, Lability
   real(kind=8) :: Porosity, DEpth, URand
   real(kind=8) :: lab_real
   
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
	   
       Particle2%Ash%Ashact = 0
       Particle2%Ash%Ashact_0 = 0
       Particle2%Ash%Ash_Time0 = 0
	   
	   Particle2%Loc = Coordinates(0,0)

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
			   
               Particle2(Particle_ID)%Ash%Ashact = 0.
               Particle2(Particle_ID)%Ash%Ashact_0 = Particle2(Particle_ID)%Ash%Ashact
			   
               Particle2(Particle_ID)%Loc = Particle(Particle_ID)%Loc
               Particle2(Particle_ID)%Loc_Init = Particle(Particle_ID)%Loc_init
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
   integer(kind=4),intent(In) :: i
   integer(kind=4) :: j, k, X, Y, Orientation
   real(kind=8)    :: URand(3)

   integer(kind=4) :: xx, yy, p_cnt

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
         
         p_cnt = 0

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
                 
                 if (IsParticle(CoOrdinates(Y,X))) p_cnt = p_cnt + 1
                 if (p_cnt >= Org(i)%Gut%Capacity) Finished = .false.
                 
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
             
                 if (IsParticle(CoOrdinates(Y,X))) p_cnt = p_cnt + 1
                 if (p_cnt >= Org(i)%Gut%Capacity) Finished = .false.
                 
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
                 
                 if (IsParticle(CoOrdinates(Y,X))) p_cnt = p_cnt + 1
                 if (p_cnt >= Org(i)%Gut%Capacity) Finished = .false.
                 
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
                 
                 if (IsParticle(CoOrdinates(Y,X))) p_cnt = p_cnt + 1
                 if (p_cnt >= Org(i)%Gut%Capacity) Finished = .false.
                 
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
         if (X < 0 .or. X>N_col .or. Y<0 .or. Y>N_row) then 
            print *, 'Error in populating:','time,i,j,X,Y',time,i,j,X,Y
            write(File_log, *) 'Error in populating:','time,i,j,X,Y',time,i,j,X,Y
         Endif
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
         if (X < 0 .or. X>N_col .or. Y<0 .or. Y>N_row) then 
            print *, 'Error in populating:','time,i,j,X,Y',time,i,j,X,Y
            write(File_log, *) 'Error in populating:','time,i,j,X,Y',time,i,j,X,Y
         Endif
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
         if (X < 0 .or. X>N_col .or. Y<0 .or. Y>N_row) then 
            print *, 'Error in populating:','time,i,j,X,Y',time,i,j,X,Y
            write(File_log, *) 'Error in populating:','time,i,j,X,Y',time,i,j,X,Y
         Endif
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
         if (X < 0 .or. X>N_col .or. Y<0 .or. Y>N_row) then 
            print *, 'Error in populating:','time,i,j,X,Y',time,i,j,X,Y
            write(File_log, *) 'Error in populating:','time,i,j,X,Y',time,i,j,X,Y
         Endif
       END Do
       END Do

   End Select
   
   if (j /= Org(i)%BodySize) then
      print*,'error during population;','time,i, j, Org(i)%BodySize;',time,i, j, Org(i)%BodySize
      write(File_log,*)'error during population;','time,i, j, Org(i)%BodySize;',time,i, j, Org(i)%BodySize
   Endif
   
   call chk_org('pplate') 
     
	   Dir_rec(1:Org(i)%length - Org(i)%width,i) = Org(i)%Orientation

         if (errChk) then 
         if (time > 0 ) then 
         call Make_matrix_chk('Pop')
         if (errDetect) then; print *,time, "err after populate"; write(file_log,*)time, "err after populate"; end if 
         call Matrix_err_chk('Pop')
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
   integer(kind=4) :: Particle_ID, Lability
   Type(Coordinates) :: Point
   real(kind=8) :: lab_real

   integer(kind=4) :: xx, yy
   
     Particle_ID                             = Particle_ID_free(PointerToFreeParticle)
     Particle_ID_free(PointerToFreeParticle) = 0
     PointerToFreeParticle                   = PointerToFreeParticle - 1
     Total_N_Particles                       = Total_N_Particles + 1
     Matrix(Point%Y,Point%X)                 = CellContent(PARTICLE_ID,P)

       Particle(Particle_ID)%Plane           = 0
       Particle(Particle_ID)%Loc             = Coordinates(Point%Y, Point%X)
       Particle(Particle_ID)%Loc_Init        = Coordinates(Point%Y, Point%X)
       Particle(Particle_ID)%Lability        = int(Lab_real(Point%Y, Point%X))
	   if (ASH_ON) Particle(Particle_ID)%Lability = 0
       Particle(Particle_ID)%Lability_Time0  = Time
       Particle(Particle_ID)%Pb%Activity     = 1
       Particle(Particle_ID)%Pb%Activity0    = 1
       Particle(Particle_ID)%Pb%Time0        = Time
       
       
       Particle(Particle_ID)%OM%OMact        = real(Lab_real(Point%Y, Point%X))
       If (ASH_ON) Particle(Particle_ID)%OM%OMact = 0.
       Particle(Particle_ID)%OM%OMact_0      = Particle(Particle_ID)%OM%OMact
       Particle(Particle_ID)%OM%OM_Time0     = Time
       
       Particle2(Particle_ID)%Ash%Ashact        = 0.
       If (ASH_ON) Particle2(Particle_ID)%Ash%Ashact = 1.
       Particle2(Particle_ID)%Ash%Ashact_0      = Particle2(Particle_ID)%Ash%Ashact
       Particle2(Particle_ID)%Ash%Ash_Time0     = Time

       Particle2(Particle_ID)%Loc             = Particle(Particle_ID)%Loc  
       Particle2(Particle_ID)%Loc_Init        = Particle(Particle_ID)%Loc_INIT  
	   
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
   integer(kind=4) :: Particle_ID
   Type(Coordinates),intent(IN) :: Point
   
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
       
       Particle2(Particle_ID)%Ash%Ashact        = 0
       Particle2(Particle_ID)%Ash%Ashact_0      = 0
       Particle2(Particle_ID)%Ash%Ash_Time0     = 0
	   
       Particle2(Particle_ID)%Loc             = Particle(Particle_ID)%Loc
       Particle2(Particle_ID)%Loc_Init        = Particle(Particle_ID)%Loc_init
   
         if (errChk) then 
         call Make_matrix_chk('RmP')
         if (errDetect) then; print *,time, "err after add"; write(file_log,*)time, "err after add"; end if 
         call Matrix_err_chk('RmP')
         if (errDetect) then; print *,time, "err after add"; write(file_log,*)time, "err after add"; end if 
		 end if 
		 
   End Subroutine Particle_remove

   ! **********************

   Subroutine Particle_sediment(ParticlesToSediment)
     ! controlling loop to input new particles to the water column

   Use GlobalVariables
   Implicit None
   Logical :: IsWater, InMAtrix
   integer(kind=4) :: m, ParticlesToSediment
   real(kind=8)    :: URand(ParticlesToSediment)
   Type(Coordinates) :: Point
   integer(kind=4) :: xx,yy

     Call Random_Number (URand)
       DO m = 1, ParticlesToSediment   ! the number of particle to sediment in this time frame
           ! Point = CoOrdinates(1, (Int(URand(m)*Real(N_Col))+1))   ! randomly locate the point at the top
           Point = CoOrdinates(1, maxloc(swi,dim=1))   !
             if (.not.inmatrix(point)) then 
                print*,'error: in particle_sediment'
                write(file_log,*)'error: in particle_sediment'
             endif
             if (IsWater(Point)) then
               Call Particle_add(Point)
               Call SedimentDown(Point)
               RainCount = RainCount + 1
             End if
             
             do xx=1,n_col
                do yy=1,n_row
                    if (matrix(yy,xx)%class==w)cycle
                    if (matrix(yy,xx)%class==p) then 
                        swi(xx)=yy
                        exit
                    endif
                enddo
            enddo
             
       END Do
       
   End Subroutine Particle_sediment

! ******************************************************************

   Subroutine SedimentDown(Point)
     ! low level routine - move particle downwards/diagonally until the local density is stable

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Logical :: InMAtrix, IsWater
   integer(kind=4) ::  TestWindow
   real(kind=8) :: Porosity_local, URand(3)
   Type(Coordinates) :: Point, PointBelow, PointBelowBelow, PointBeside
   Type(oxygen_conc) :: O2_buff  !! YK

   integer(kind=4) :: xx, yy
   
   ! Initialisations:
       Call Random_Number (URand)
       TestWindow = 3 * Urand(3) + 1      ! i.e. a box of 3X3 cells wide

   ! Main loop
       Do

         PointBelow = CoOrdinates(Point%Y+1,Point%X)
         PointBelowBelow = CoOrdinates(PointBelow%Y+1,PointBelow%X)
         
         if (.not.inmatrix(pointbelow)) then 
            print*,'error: in sedimentDown, pointbelow'
            write(file_log,*)'error: in sedimentDown, pointbelow'
         endif

         If (.not. IsWater(PointBelow)) then
           Exit                                      ! i.e., no need to move down
         Else 
         if (.not.inmatrix(pointbelowbelow)) then 
            print*,'error: in sedimentDown, pointbelowbelow'
            write(file_log,*)'error: in sedimentDown, pointbelowbelow'
         endif                                       ! check if water below that point ..

           If (IsWater(PointBelowBelow)) then        ! all clear, move down
             O2_buff = O2(PointBelow%Y,PointBelow%X) 
             Matrix(PointBelow%Y,PointBelow%X) = Matrix(Point%Y,Point%X)
             if (oxygen_ON) O2(PointBelow%Y,PointBelow%X) = O2(Point%Y,Point%X)
             Particle(Matrix(PointBelow%Y,PointBelow%X)%Value)%loc = PointBelow
             Particle2(Matrix(PointBelow%Y,PointBelow%X)%Value)%loc = PointBelow
             
 
Particle(Matrix(PointBelow%Y,PointBelow%X)%Value)%loc_INIT = PointBelow
Particle2(Matrix(PointBelow%Y,PointBelow%X)%Value)%loc_INIT = PointBelow
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
               Particle2(Matrix(PointBelow%Y,PointBelow%X)%Value)%loc = PointBelow
 
Particle(Matrix(PointBELOW%Y,PointBeLOW%X)%Value)%loc_Init = PointBelow
Particle2(Matrix(PointBELOW%Y,PointBeLOW%X)%Value)%loc_Init = PointBelow
                 ! THIS IS THE FINAL SETTLING POINT
               Matrix(Point%Y,Point%X) = CellContent(w,w)
               if (oxygen_ON) O2(Point%Y,Point%X) = O2_buff
               Point = PointBelow
             End if
             Exit
           End if
         End if
       End do

     ! If (URand(1) .gt. 0.3) then                     ! in 30% cases do not shift particles to the side
     If (URand(1) .gt. 1.0) then                     ! in 100% cases do not shift particles to the side
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
             Particle2(Matrix(PointBeside%Y,PointBeside%X)%Value)%loc = PointBeside
 
Particle(Matrix(PointBeside%Y,PointBeside%X)%Value)%loc_Init = PointBeside
Particle2(Matrix(PointBeside%Y,PointBeside%X)%Value)%loc_Init = PointBeside
               ! THIS IS THE FINAL SETTLING POINT
           End if
         End if
     End if

         if (errChk) then 
         call Make_matrix_chk('Dow')
         if (errDetect) then; print *,time, "err after sedown"; write(file_log,*)time, "err after sedown"; end if 
         call Matrix_err_chk('Dow')
         if (errDetect) then; print *,time, "err after sedown"; write(file_log,*)time, "err after sedown"; end if 
		 end if 
     
   End Subroutine SedimentDown

   ! *******************************

   Subroutine WaterColumn_scan(DepthToScan)
     ! scan the water column and sediment down any particle found .. search intensity decreases with depth

   Use GlobalVariables
   Implicit none
   integer(kind=4) :: Y,X, DepthToScan
   Logical :: IsPArticle
   real(kind=8)    :: Pr
   real(kind=8)    :: URand(DepthToScan)

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
   integer(kind=4) :: Y_new, Y, X, DummyVariable, PArticle_ID
   real(kind=8) :: Lab_real_New

   integer(kind=4) :: xx, yy
   integer(kind=4) :: i
   
   real(kind=8) :: m_OM,m_poro, m_ash,tmp_lab(N_Col),tmp_ash(N_Col)
   integer(kind=4) :: tmp_img(N_Col)

   ! check bottom boundary to see if the shift is possible
       Y = N_Row

       Do X = 1, N_Col
         If (IsOrganism(CoOrdinates(Y,X))) Return  ! SHIFT NOT POSSIBLE
       End do
       
       m_OM = 0d0
       m_poro = 0d0
	   m_ash = 0d0
       tmp_img = 0
       tmp_lab = 0d0
       tmp_ash = 0d0

   ! As there are no boundary problems .. do the shift
       Do X = 1, N_Col
         if (IsParticle(CoOrdinates(Y,X))) then
           tmp_img(X) = 1
           tmp_lab(X) = Particle(Matrix(Y,X)%Value)%OM%OMact + 1d0
           tmp_ash(X) = Particle2(Matrix(Y,X)%Value)%Ash%Ashact + 1d0
           m_OM = m_OM + Particle(Matrix(Y,X)%Value)%OM%OMact 
           m_ash = m_ash + Particle2(Matrix(Y,X)%Value)%Ash%Ashact 
           m_poro = m_poro + 1d0
           Call Particle_remove(CoOrdinates(Y,X))
         endif
       End do
       
       m_OM = m_OM/N_col
       m_poro = m_poro/N_col
       m_ash = m_ash/N_col
       
       write(File_Core_M,*) (tmp_img(x),x=1,N_col)
       write(File_Core_L,*) (tmp_lab(x),x=1,N_col)
       write(File_Core_A,*) (tmp_ash(x),x=1,N_col)
       write(File_Core,*) Time*Timescale, m_poro, m_OM, m_ash

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
				 Particle2(Particle_ID)%loc=particle(particle_id)%loc
				 Particle2(Particle_ID)%loc_init=particle(particle_id)%loc_init
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
         call Make_matrix_chk('Cns')
         if (errDetect) then; print *,time, "err after constrain(sub)"; write(file_log,*)time, "err after constrain(sub)"; end if 
         call Matrix_err_chk('Cns')
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
   integer(kind=4),intent(IN) :: i, Orientation
   integer(kind=4) :: k, m, o, X, Y, Left, Forward, Right
   integer(kind=4) :: DirPossible(3), StuckCount
   integer(kind=4) :: N_Particle, N_Organism, N_Water
   integer(kind=4) :: X_start, X_end, X2_start, X2_end, dX, dY, Y_start, Y_end
   real(kind=8)  :: DirRandom(3), DirWeights(3), URand(3)
   Type(CellContent) :: ToMove
   Type(Coordinates) :: ToMove_loc
   Type(Coordinates),intent(IN) :: Point
   Type(oxygen_conc) :: O2_buff   !!  YK
   integer(kind=4) :: X_buff, Y_buff

   integer(kind=4) :: xx, yy, xxx,yyy
   type(cellcontent) :: blcimg(org(i)%width)
   type(coordinates) :: blcimg_loc(org(i)%width)
   integer(kind=4) :: blcimg_loc_x(org(i)%width),blcimg_loc_y(org(i)%width)
   
   ! if (CurrentEnergy < 1.) return

   Stuck = .false.
   StuckCount = 0
   m = 0
   blcimg_loc = coordinates(0,0)
   blcimg_loc_x = 0
   blcimg_loc_y = 0

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
			   blcimg_loc_x(o) = particle(block(o)%value)%loc%x
			   blcimg_loc(o)%y = particle(block(o)%value)%loc%y
			   blcimg_loc_y(o) = particle(block(o)%value)%loc%y
               ! if (blcimg_loc(o)%x==0 .or. blcimg_loc(o)%y==0)then 
               if (blcimg_loc_x(o)==0 .or. blcimg_loc_y(o)==0)then 
                    print*,'error in making blk img: ',time,block(o)%class,block(o)%value,blcimg_loc_x(o) ,blcimg_loc_y(o)
                    write(file_log,*)'error in making blk img: ',time,block(o)%class,block(o)%value,blcimg_loc_x(o) ,blcimg_loc_y(o)
                endif
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
               
   ! EnergyHistory(i,1) = EnergyHistory(i,1) - 0.01   ! increment the running average
   ! CurrentEnergy = currentEnergy - 0.01
   
   if (size(blcimg)/=Org(i)%width .or. size(block)/=Org(i)%width) then 
      print*,'error during push; block is wierd'
      write(File_log,*)'error during push; block is wierd'
   endif
   
   do o = 1, org(i)%width
     if (.not. (blcimg(o)%class==p .or. blcimg(o)%class==w)) then 
       if (.not.(N_Organism .gt. 0)) then 
      print*,'error during push; organisms in block ?',blcimg(o)%class, blcimg(o)%value
      write(File_log,*)'error during push;  organisms in block ?',blcimg(o)%class, blcimg(o)%value
      endif
    endif
   enddo
   
   do o = 1, org(i)%width
     if (blcimg(o)%class/=p) cycle
     if (blcimg(o)%value >= N_col*N_row) then 
        print *,'error during particle push',blcimg(o)%class,blcimg(o)%value
        write(File_log,*) 'error during particle push',blcimg(o)%class,blcimg(o)%value
     endif
     xx = particle(blcimg(o)%value)%loc%x
     yy = particle(blcimg(o)%value)%loc%y
     ! xxx = blcimg_loc(o)%x
     ! yyy = blcimg_loc(o)%y
     xxx = blcimg_loc_x(o)
     yyy = blcimg_loc_y(o)
     if (xx==0 .or. yy==0 .or. xxx==0 .or. yyy==0) then 
     if (.not. any( Particle_ID_free == blcimg(o)%value)) then 
     print*,'error in pushing; time, i, o, xx, yy, class, value ',time, i, o, xx, yy,xxx,yyy &
        ,blcimg(o)%class,blcimg(o)%value
    print*,Particle_ID_free(:)
     write(file_log,*)'error in pushing; time, i, o, xx, yy, class, value ',time, i, o, xx, yy,xxx,yyy &
        ,blcimg(o)%class,blcimg(o)%value
     stop
     endif
     endif
     if (particle(blcimg(o)%value)%loc%x /= blcimg_loc(o)%x) then  ! particle is pushed 
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
         call Make_matrix_chk('Psh')
         if (errDetect) then; print *,time, "err after push"; write(file_log,*)time, "err after push"; end if 
         call Matrix_err_chk('Psh')
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
   integer(kind=4), Intent(in) :: X, Y
   Type(CellContent), Intent(InOut) :: ToMove
   Type(CellContent)  :: Buffer

   Select Case (Matrix(Y,X)%Class)
     Case (p)
       Buffer = Matrix(Y,X)
       Matrix(Y,X) = ToMove
       Particle(Matrix(Y,X)%Value)%Loc = Coordinates(Y,X)   ! YK assign loc of particle id in tomove to x,y 
       Particle2(Matrix(Y,X)%Value)%Loc = Coordinates(Y,X)   ! YK assign loc of particle id in tomove to x,y 
       ToMove = Buffer        !!  YK   tomove is replaced by Matrix(Y,X)
       Finished = .false.
     Case (w)
       Buffer = Matrix(Y,X)  !  YK added  buffer is now water 
       Matrix(Y,X) = ToMove   !  YK matrix(Y,X) is now particle(ToMove)
       Particle(Matrix(Y,X)%Value)%Loc = Coordinates(Y,X)      !!  particle location is moved to Y,X
       Particle2(Matrix(Y,X)%Value)%Loc = Coordinates(Y,X)      !!  particle location is moved to Y,X
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
   integer(kind=4) :: i, k, m, n, o, Orientation
   integer(kind=4) :: N_Particle, N_Organism, N_Water
   real(kind=8)    :: URand(Org(i)%Width)
   Type(Coordinates) :: Point
   real(kind=8) :: n_real
   Type(oxygen_conc) :: O2_buff   !!  YK
   
   ! if (CurrentEnergy < 1.) return

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
               
               
               ! EnergyHistory(i,1) = EnergyHistory(i,1) - 0.01  ! increment the running average
               ! CurrentEnergy = currentEnergy - 0.01
	           
			   
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
         call Make_matrix_chk('Igt')
         if (errDetect) then; print *,time, "err after ingest"; write(file_log,*)time, "err after ingest"; end if 
         call Matrix_err_chk('Igt')
         if (errDetect) then; print *,time, "err after ingest"; write(file_log,*)time, "err after ingest"; end if 
         end if 
		 
   End SUBROUTINE Particles_ingest

   ! **************************************************************************

   SUBROUTINE Particles_egest(i)
     ! organism i egests particles

   Use GlobalVariables
   use ieee_arithmetic
   Implicit None
   Logical  :: IsWater, inmatrix
   integer(kind=4) :: i, j, k, loc, centre
   real(kind=8) :: biodecay               !!!YK 
   real(kind=8) :: timeIncrement, ox_buf
   type(coordinates) :: buff_loc
   
   
   ! if (CurrentEnergy < 1.) return
   
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
             
             
         if (.not.inmatrix(tail(loc))) then 
            print*,'error: in particles_egest'
            write(file_log,*)'error: in  particles_egest'
         endif
             
             If (IsWater(Tail(loc))) then
               ! IF (Particle(Guts(j,i)%Value)%OM%OMact .gt. 0) then
               Matrix(Tail(loc)%Y, Tail(loc)%X) = Guts(j,i)  !! fill the position of water in tail with a particle 
               Particle(Guts(j,i)%Value)%loc = Tail(loc)   !! location of particle is set at tail position 
               Particle2(Guts(j,i)%Value)%loc = Tail(loc)   !! location of particle is set at tail position 
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
              
              ! End if

               Guts(j,i) = CellContent(w,w)
               EgestHistory(i,1) = EgestHistory(i,1) + 1   ! increment the running average
               
               
               ! EnergyHistory(i,1) = EnergyHistory(i,1) - 0.01   ! increment the running average
               ! CurrentEnergy = currentEnergy - 0.01
	   
	   
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
               
               ! end if 
             End if
         End if
       End do
       
   End if
   
         If (errChk) then 
         call Make_matrix_chk('Egt')
         if (errDetect) then; print *,time, "err after egest"; write(file_log,*)time, "err after egest"; end if 
         call Matrix_err_chk('Egt')
         if (errDetect) then; print *,time, "err after egest"; write(file_log,*)time, "err after egest"; end if 
         end if 
		 
   End Subroutine Particles_egest

   ! ********************************

   Subroutine Organism_peristalsis(i)
     ! shift gut contents down in the gut ..

   Use GlobalVariables
   Implicit None
   integer(kind=4) :: i, j, k
   Type(CellContent) :: Image_Gut(Org(i)%Gut%Capacity)
   
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
         call Make_matrix_chk('Pss')
         if (errDetect) then; print *,time, "err after peristalsis"; write(file_log,*)time, "err after peristalsis"; end if 
         call Matrix_err_chk('Pss')
         if (errDetect) then; print *,time, "err after peristalsis"; write(file_log,*)time, "err after peristalsis"; end if 
         end if 
		 
   End Subroutine Organism_peristalsis

   ! ********************************

   Subroutine disperse()
     ! disperse sediment particles when local porosity approaches 0 ..

   Use GlobalVariables
   Implicit None
   integer(kind=4) :: xx,yy, xx_l,xx_r,xx_ll, xx_rr
   real(kind=8) :: porosity_local
   real(kind=8) :: poro_tmp
   integer(kind=4) :: poro_shift = 2  
   type (oxygen_conc) :: o2_buf
   
   do yy=1,n_row
   do xx=1,n_col
        poro_tmp = porosity_local(coordinates(yy,xx),poro_shift)  ! calc local porosity from x - shift to x + shift 
        if (poro_tmp == 0.05) then 
            xx_l = xx - poro_shift
            xx_ll = xx_l - 1
            xx_r = xx + poro_shift
            xx_rr = xx_r + 1
            if (xx_l < 1) xx_l = xx_l + n_col
            if (xx_r > n_col) xx_r = xx_r - n_col
            if (xx_ll < 1) xx_ll = xx_ll + n_col
            if (xx_rr > n_col) xx_rr = xx_rr - n_col
            
            if (matrix(yy,xx_ll)%class==w .and. &
                matrix(yy,xx_l)%class==p ) then
                ! recording a new x in particle 
                particle(matrix(yy,xx_l)%value)%loc%x &
                    = xx_ll
                particle2(matrix(yy,xx_l)%value)%loc%x &
                    = xx_ll
                ! matrix content is exchanged 
                matrix(yy,xx_ll)%class = p
                matrix(yy,xx_ll)%value &
                    = matrix(yy,xx_l)%value
                matrix(yy,xx_l)%class = w
                matrix(yy,xx_l)%value = 0
                ! oxygen conc. exchange
                o2_buf = o2(yy,xx_ll)
                o2(yy,xx_ll)  &
                    = o2(yy,xx_l)
                o2(yy,xx_l) = o2_buf
            endif
            if (matrix(yy,xx_rr)%class==w .and. &
                matrix(yy,xx_r)%class==p ) then
                ! recording a new x in particle 
                particle(matrix(yy,xx_r)%value)%loc%x &
                    = xx_rr
                particle2(matrix(yy,xx_r)%value)%loc%x &
                    = xx_rr
                ! matrix content is exchanged 
                matrix(yy,xx_rr)%class = p
                matrix(yy,xx_rr)%value &
                    = matrix(yy,xx_r)%value
                matrix(yy,xx_r)%class = w
                matrix(yy,xx_r)%value = 0
                ! oxygen conc. exchange
                o2_buf = o2(yy,xx_rr)
                o2(yy,xx_rr)  &
                    = o2(yy,xx_r)
                o2(yy,xx_r) = o2_buf
            endif
        endif
    enddo
    enddo
       
         if (errChk) then 
         call Make_matrix_chk('Dsp')
         if (errDetect) then; print *,time, "err after dispersion"; write(file_log,*)time, "err after dispersion"; end if 
         call Matrix_err_chk('Dsp')
         if (errDetect) then; print *,time, "err after dispersion"; write(file_log,*)time, "err after dispersion"; end if 
         end if 
		 
   End Subroutine disperse


   ! ***********************************************************************

   SubRoutine OutputData()
     ! main data output routine

   Use GlobalVariables
   use O2_diffusion
   Implicit None

   Call Output_stats()
   ! Call Output_Graphic()
     ! Call Output_Ascii() !  if ascii data desired
     Call Output_txtImg() 
     Call Output_O2txtImg() 

   End Subroutine OutputData

   ! ***********************************************************************

   SubRoutine Output_ascii()
     ! output the data in ascii format

   Use GlobalVariables
   Implicit None
   Logical :: IsWater, IsPArticle, IsOrganism
   integer(kind=4) :: X, Y, Code
   Character*21 AsciiFile, AsciiFile2, numtemp
   integer(kind=4) :: tmppp
   integer(kind=4) :: txtimg(N_Col)
   
   tmppp = 1
    
   AsciiFile = CurrentTime //trim(adjustl(today))//'.ascii'
   write(numtemp,'(i10.1)') SaveTime
   OPEN(unit = File_ASCII, File = trim(adjustl(today))//'/data/'    &
                //trim(adjustl(numtemp))//'.ascii', status = 'unknown') 
   OPEN(unit = File_txtImg_2, File = trim(adjustl(today))//'/data/'    &
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
   integer(kind=4) :: X, Y
   Character*21 numtemp
   integer(kind=4) :: txtimg(N_Col), particle_ID
   real(kind=8) :: Lab_real, txtimg_real(N_col)
   integer(kind=4) :: i , j
   ! logical :: rec_detail = .false.
   logical :: rec_detail = .true.
    
   write(numtemp,'(i10.1)') Time
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/txtimg-'          &
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
   
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/lability-'          &
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
   
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/ash-'          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO Y = 1, N_row
     txtimg_real = 0 
     DO X = 1, N_col
       If (IsParticle(CoOrdinates(Y,X))) then
         Particle_ID = Matrix(Y,X)%Value
         Lab_real = particle2(particle_ID)%Ash%Ashact
         txtimg_real(X) = 1 + Lab_real
       End if
     END Do
     write(File_txtImg, *) (txtimg_real(x), x = 1, n_col)
     END Do
   Close(File_txtimg)

   if (rec_detail) then 
      OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/ptcl_list-'          &
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
   
      OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/ptcl_list_2-'          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO i = 1, N_cell
         write(File_txtimg,*) i, Particle(i)%loc%Y,Particle(i)%loc%X, Particle(i)%plane
     END Do
   Close(File_txtimg)

   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/guts_list-'          &
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
   
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/org_loc_list-'         &
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

   SubRoutine Make_matrix_chk(chr)
   
   use GlobalVariables
   implicit none 
   integer(kind=4) :: i, j, xx, yy
   integer(kind=4) :: txtimg(N_Col)
   character*3, intent(IN) :: chr
   
   if (allocated(matrix_chk)) deallocate(matrix_chk)
   allocate(matrix_chk(N_row,N_col))
   matrix_chk = cellcontent(0,0)  ! all water
   errDetect = .false.
   
   do i = 1, N_cell
      if ((particle(i)%loc%x == 0) .and. (particle(i)%loc%y == 0)) cycle
      if (any(Guts(:,:)%value == i)) cycle
      matrix_chk(particle(i)%loc%y, particle(i)%loc%x) = cellcontent(i,p)
	  if (particle(i)%loc%x /= particle2(i)%loc%x .or.particle(i)%loc%y /= particle2(i)%loc%y &
	    .or. particle(i)%loc_init%x /= particle2(i)%loc_init%x  &
		.or. particle(i)%loc_init%y /= particle2(i)%loc_init%y) then 
	    print*, "error making matrix_chk from "//chr// &
          ": inconsistency between particle 1 & 2: time,cell#, loc1, loc2,locint1,locint2 ="  &
		  ,time,i,particle(i)%loc,particle(i)%loc%x,particle(i)%loc,particle(i)%loc%y,  &
		  particle2(i)%loc,particle(i)%loc%x,particle2(i)%loc,particle(i)%loc%y,  &
		  particle(i)%loc,particle(i)%loc_init%x,particle(i)%loc,particle(i)%loc_init%y,  &
		  particle2(i)%loc,particle(i)%loc_init%x,particle2(i)%loc,particle(i)%loc_init%y 
	    write(file_log,*) "error making matrix_chk from "//chr//  &
          ": inconsistency between particle 1 & 2: time,cell#, loc1, loc2,locint1,locint2 ="  &
		  ,time,i,particle(i)%loc,particle(i)%loc%x,particle(i)%loc,particle(i)%loc%y,  &
		  particle2(i)%loc,particle(i)%loc%x,particle2(i)%loc,particle(i)%loc%y,  &
		  particle(i)%loc,particle(i)%loc_init%x,particle(i)%loc,particle(i)%loc_init%y,  &
		  particle2(i)%loc,particle(i)%loc_init%x,particle2(i)%loc,particle(i)%loc_init%y 
      endif
   end do 
   
   do i = 1, N_ind
      do j = 1, Org(i)%bodySize
        if (time==0) print *, Org_loc(j,i)%Y,Org_loc(j,i)%X, j, Org(i)%bodysize
        xx = Org_loc(j,i)%X
        yy = Org_loc(j,i)%Y
        if (yy == 0 .or. xx ==0) then
          write(file_log,*) time, 'error making matrix_chk from '//chr//  &
             ': Org_loc records 0 in x or y within Org(i)%bodysize'  &
               ,'i, j, Org_loc(j,i)%Y, Org_loc(j,i)%X, Org(i)%bodySize:'  &
               ,i, j, yy, xx, Org(i)%bodySize
          write(file_log,*) Org_loc(:,i)%Y
          write(file_log,*) Org_loc(:,i)%X
          errDetect = .true.
        endif 
        if (Org_loc(j,i)%Y==0 .or. Org_loc(j,i)%X==0) print *,Org_loc(j,i)%Y,Org_loc(j,i)%X,j, Org(i)%bodysize, i,chr
        if (matrix_chk(Org_loc(j,i)%Y,Org_loc(j,i)%X)%class == p) then  !!  if particle (not in guts) exists
          print *, time, "error making matrix_chk from "//chr//  &
            ": particle exist where organism",i," exist at", Org_loc(j,i)%Y,Org_loc(j,i)%X
          write(file_log,*) time, "error making matrix_chk from "//chr//  &
            ": particle exist where organism",i," exist at",& 
		                    Org_loc(j,i)%Y,Org_loc(j,i)%X
          errDetect = .true.
        end if 
        matrix_chk(Org_loc(j,i)%Y,Org_loc(j,i)%X) = cellcontent(j,i)
      end do
   end do
   
   
   End Subroutine Make_matrix_chk
   
 ! ****************************************************************************************************

   SubRoutine Matrix_err_chk(chr)
   
   use GlobalVariables
   implicit none 
   integer(kind=4) :: xx, yy
   character*3, intent(IN) :: chr
   
   errDetect = .false.
   
   do yy = 1, N_row
     do xx = 1, N_col
       if (matrix(yy,xx)%class /= matrix_chk(yy,xx)%class) then 
            write(file_log,*) time, 'error in water/organism assignment in matrix (err_chk_mtx) from ' & 
            //trim(adjustl(chr))//':', yy, xx, & 
                matrix_chk(yy,xx)%class, matrix(yy,xx)%class, matrix_chk(yy,xx)%value, matrix(yy,xx)%value
            errDetect = .true.
       else 
         if (matrix(yy,xx)%value /= matrix_chk(yy,xx)%value) then
           write(file_log,*) time, '(err_chk_mtx) from ' & 
            //trim(adjustl(chr))//': Id is different at', yy, xx, & 
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
   integer(kind=4) :: i, j, xx, yy
   integer(kind=4) :: txtimg(N_Col)
   Character*21 numtemp
   
   
   write(numtemp,'(i10.1)') Time
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/chk-'          &
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
   integer(kind=4) :: y, x
   real(kind=8) :: lab_real_new
   real(kind=8) :: aaa
   
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
   integer(kind=4) :: k, o, i
   real(kind=8) :: resp_rate, lab_max, ox_resp, TimeIncrement
   real(kind=8) :: fact_law1, fact_law2
   real(kind=8) :: merge2
   
   
   select case (trim(adjustl(O2ratelaw)))
     case('linear','LINEAR','Linear')
        fact_law1 = 0d0
        fact_law2 = 1d0
     case('zero','ZERO','Zero')
        fact_law1 = 0d0
        fact_law2 = 0d0
     case('monod','MONOD','Monod')
        fact_law1 = 1d0
        fact_law2 = 0d0
   end select 
   
   
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
   
   ! resp_rate = lab_max*Lability_decayConstant(N_LabilityClasses)*fact*bio_fact  ! calc wt% loss/yr once oxygen conc. is given
   resp_rate = lab_max*RespCnst(i) ! calc wt% loss/yr once oxygen conc. is given
   
   ! resp_rate = resp_rate*1e-1   !! if decomposition rate is increased
   
   ox_resp = 0.0
   
   do k = 1, Org(i)%headSize
     O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen_use = &  
       O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen_use + & 
       resp_rate   & 
       *merge(1d0/mo2,1d0, fact_law1 == 1d0 .and.O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen <= mo2/iox)
     ox_resp = ox_resp + &
       resp_rate*merge2(1.d0,O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen   &
       ,(fact_law1 == 0d0 .and. fact_law2 == 0d0).or.   &
       (fact_law1==1d0 .and.O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen > mo2/iox))*TimeScale/365.*iox/OM_uni
   end do 
   
   
   particle(guts(o,i)%value)%OM%OMact = particle(guts(o,i)%value)%OM%OMact - ox_resp
   
   if (particle(guts(o,i)%value)%OM%OMact < 0) particle(guts(o,i)%value)%OM%OMact = 0.0
   
   RespHistory(i,1) = 1
   
   
   end subroutine OrgResp
   ! *********************************************************
   subroutine fluxes()
   
   use GlobalVariables
   implicit none 
   integer(kind=4) :: xx, yy, k, i
   integer(kind=4) :: xg, xp, yp, yg
   real(kind=8) :: tmpV, tmpU
   real(kind=8) :: fact_law1, fact_law2
   ! real(kind=8) :: width_3d
   real(kind=8) :: merge2

  ! width_3d = 2.5d0   
   
   if (.not.oxygen_ON) y_int = 1
   
   select case (trim(adjustl(O2ratelaw)))
     case('linear','LINEAR','Linear')
        fact_law1 = 0d0
        fact_law2 = 1d0
     case('zero','ZERO','Zero')
        fact_law1 = 0d0
        fact_law2 = 0d0
     case('monod','MONOD','Monod')
        fact_law1 = 1d0
        fact_law2 = 0d0
   end select 
   
   TotO2Dif = 0.0
   do xx = 1, n_col
     totO2Dif = totO2Dif +   &
       (O2(y_int,xx)%oxygen - O2(y_int+1,xx)%oxygen)*edif(y_int,xx)*iox*1e-3   &
          /pixelSize*(width_3d*pixelSize)
   end do
   
   TotO2Dif = TotO2dif/width_3d/(n_col*pixelSize)  !  mol /cm2 / yr flux
   
   TotO2Adv = 0.0
   do xx = 1, n_col
     TotO2Adv = TotO2Adv +  &
	   O2(y_int,xx)%oxygen*Vo(xx,y_int)*iox*1e-3*(width_3d*pixelSize*pixelsize)
   end do 
   
   TotO2Adv = TotO2Adv/width_3d/(n_col*pixelSize) 
   
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
        O2(yy,xx)%Oxygen_use*merge2(1.d0,O2(yy,xx)%oxygen   &
       ,(fact_law1 == 0d0 .and. fact_law2 == 0d0).or.   &
       (fact_law1==1d0 .and. O2(yy,xx)%oxygen > mo2/iox))  &
          *iox*1e-3*(width_3d*pixelSize*pixelsize)
      end do 
    end do 
	
   TotResp = 0.0
   if (resp_ON) then 
   do i = 1, N_ind
   do k = 1, Org(i)%headSize
   if (matrix(Org_loc(k,i)%Y,Org_loc(k,i)%X)%class/=i) print *, "error in calc. flux"
   TotResp = TotResp +  &   
       O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen_use*merge2(1.d0,O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen   &
       ,(fact_law1 == 0d0 .and. fact_law2 == 0d0).or.   &
       (fact_law1==1d0 .and. O2(Org_loc(k,i)%Y,Org_loc(k,i)%X)%oxygen > mo2/iox))  &
       *iox   &
         *1e-3*width_3d*(pixelSize)*(PixelSize)
   end do 
   
   ! EnergyHistory(i,1) = EnergyHistory(i,1) + TotResp*timescale/365.*5d6   ! increment the running average
   ! CurrentEnergy = currentEnergy + TotResp*timescale/365.*5d6
   
   end do
   end if 
   
   ToTAbio = TotOrgDecay - TotResp
   
   TotAbio = TotAbio/width_3d/(n_col*pixelSize)
   
   TotResp = TotResp/width_3d/(n_col*pixelSize) ! 
   
   TotOrgDecay = TotOrgDecay /width_3d/(n_col*pixelSize) 
   
   write(File_flux,*) Time*Timescale, TotO2dif,TotO2Adv, TotAbio, TotResp, TotOrgDecay, (totO2-pretotO2)/timescale*365.0
   write(*,*) "FLUXES_v1 --- :", TotO2dif,TotO2Adv,TotAbio, TotResp, TotOrgDecay, (totO2-pretotO2)/timescale*365.0
   
   
   End Subroutine fluxes
   
   !****************************************
   
   subroutine Gnuplot_flux(chr)
   
   use GlobalVariables
   implicit none 
   character*1, intent(IN) :: chr
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/biot-flx'//trim(adjustl(chr))//'.plt', status = 'unknown')
   write(File_txtImg,*)'set pm3d map  corners2color max'
   write(File_txtImg,*)'set size square'
   write(File_txtImg,*)'set border lw 3'
   write(File_txtImg,*)'set size ratio 0.24'
   write(File_txtImg,*)"set tics font 'arial, 35'"
   write(File_txtImg,*)'set palette rgbformula 22,13,-31'
   write(File_txtImg,*)"r1='"//trim(adjustl(today))//'/flux'//trim(adjustl(chr))//'.OUT'//"'"
   write(File_txtImg,*)'unset colorbox'
   write(File_txtImg,*)'set key outside horizontal bottom'
   write(File_txtImg,*)"plot r1 u ($1/365.0):($6*1e6) title 'TotOrgDecay' w l lw 1, \"
   write(File_txtImg,*)"  r1 u ($1/365.0):($5*1e6) title 'TotResp' w l lw 1, \"
   write(File_txtImg,*)"  r1 u ($1/365.0):($4*1e6) title 'TotAbio' w l lw 1, \"
   write(File_txtImg,*)"  r1 u ($1/365.0):($2*1e6) title 'TotO2dif' w l lw 1, \"
   write(File_txtImg,*)"  r1 u ($1/365.0):($7*1e3) title 'dO2/dt' w l lw 1"
   write(File_txtImg,*)'set terminal emf enhanced "Arial, 25"'
   write(File_txtImg,*)'set terminal emf enhanced "Arial, 22.4"'
   write(File_txtImg,*)'set output "biot-flux.emf"'
   write(File_txtImg,*)'replot'
   write(File_txtImg,*)'set output'
   write(File_txtImg,*)'set terminal wxt'
   Close(File_txtimg)
   
   
   End subroutine Gnuplot_flux
   
   !****************************************
   
   subroutine Gnuplot_diet(i)
   
   use GlobalVariables
   implicit none 
   character*36 dumID
   integer(kind=4), intent(in) :: i
   
   write(DumID,'(i3.3)') PopLogID(i)
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/biot-diet-'//trim(adjustl(dumID))//'.plt', status = 'unknown')
   write(File_txtImg,*)'set pm3d map  corners2color max'
   write(File_txtImg,*)'set size square'
   write(File_txtImg,*)'set border lw 3'
   write(File_txtImg,*)'set size ratio 0.24'
   write(File_txtImg,*)"set tics font 'arial, 35'"
   write(File_txtImg,*)'set palette rgbformula 22,13,-31'
   write(File_txtImg,*)"r1='"//trim(adjustl(today))//'/Diet-'//trim(adjustl(dumID))//'.OUT'//"'"
   write(File_txtImg,*)'unset colorbox'
   write(File_txtImg,*)'set key outside horizontal bottom'
   write(File_txtImg,*)'# Time*Timescale,i,MovementRate, IngestRate, EgestRate, RespRate_ave,Org(i)%Gut%Content' &
                        //',Org(i)%Move%Rate,Org(i)%Ingest%Rate,Org(i)%Gut%Capacity'
   write(File_txtImg,*)"plot r1 u ($1/365.0):($3/$8) title 'Rel. Move Rate' w l lw 1, \"
   write(File_txtImg,*)"  r1 u ($1/365.0):($4/$9) title 'Rel. Ingest Rate' w l lw 1, \"
   write(File_txtImg,*)"  r1 u ($1/365.0):($5/$9) title 'Rel. Egest Rate' w l lw 0.8 dt (10,5), \"
   write(File_txtImg,*)"  r1 u ($1/365.0):($6) title 'Resp. Rate' w l lw 1, \"
   write(File_txtImg,*)"  r1 u ($1/365.0):($7/$10) title 'Fullness' w l lw 1, \"
   write(File_txtImg,*)"  r1 u ($1/365.0):($11) title 'Energy' w l lw 1 "
   write(File_txtImg,*)'set terminal emf enhanced "Arial, 25"'
   write(File_txtImg,*)'set terminal emf enhanced "Arial, 22.4"'
   write(File_txtImg,*)'set output "biot-diet-'//trim(adjustl(dumID))//'.emf"'
   write(File_txtImg,*)'replot'
   write(File_txtImg,*)'set output'
   write(File_txtImg,*)'set terminal wxt'
   Close(File_txtimg)
   
   
   End subroutine Gnuplot_diet
   
   !****************************************
   
   subroutine Gnuplot_Db(dimchr)
   
   use GlobalVariables
   implicit none 
   character*1, intent(in) :: dimchr
   character*10 :: colchr
   
   select case (trim(adjustl(dimchr)))
     case('x','X')
       colchr = '8'
     case('y','Y')
       colchr = '7'
     case default
       colchr = '6'
   end select 
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/biot-Db'//trim(adjustl(dimchr))//'.plt', status = 'unknown')
   write(File_txtImg,*)'set border lw 3.5'
   write(File_txtImg,*)'set size ratio 1.25'
   write(File_txtImg,*)"set tics font 'arial, 35'"
   write(File_txtImg,*)'set format x2 "%.1f"'
   write(File_txtImg,*)'set palette rgbformula 22,13,-31'
   write(File_txtImg,*)"r1='"//trim(adjustl(today))//'/MeanDisplacement2.OUT'//"'"
   write(File_txtImg,*)'unset colorbox'
   write(File_txtImg,*)'set key outside vertical right'
   write(File_txtImg,*)'set yrange [] reverse'
   ! write(File_txtImg,*)'set yr[12:0]'
   ! write(File_txtImg,*)'set xr[-7:1]'
   ! write(File_txtImg,*)'set xtics 2'
   write(File_txtImg,*)'n = ',Savetime
   write(File_txtImg,*)'plot for [i=1:n] r1 u (log10($'//trim(adjustl(colchr))//')):2 every :::(i)::(i) title' &
     //' sprintf("%g",(i-1)*25) w l  lw 1  lc palette frac (i)/(n*1.)'
   write(File_txtImg,*)'set terminal emf enhanced "Arial, 25"'
   write(File_txtImg,*)'set output "biot-Db'//trim(adjustl(dimchr))//'.emf"'
   write(File_txtImg,*)'replot'
   write(File_txtImg,*)'set output'
   write(File_txtImg,*)'set terminal wxt'
   Close(File_txtimg)
   
   
   End subroutine Gnuplot_Db
   
   ! ***********************************************************   
   
   subroutine org_midloc(i, startj, endj,  org_mid)
   
   use globalvariables
   implicit none 
   type(coordinates) :: org_mid
   integer(kind=4) :: i, startj, endj, j
   integer(kind=4) :: maxx, maxy, minx, miny
   integer(kind=4) :: org_copy(endj - startj + 1)
   
   maxx = maxval(org_loc(startj:endj,i)%x)
   maxy = maxval(org_loc(startj:endj,i)%y)
   minx = minval(org_loc(startj:endj,i)%x)
   miny = minval(org_loc(startj:endj,i)%y)
   
   ! print*,maxx,maxy,minx,miny
   
   org_mid%y = int((maxy + miny )/2.0)
   
   if ((maxx - minx <= 2*org(i)%width)) then 
   org_mid%x = int((maxx + minx)/2.0)
   else if (maxx - minx > 2*org(i)%width) then 
      ! do j = 1, endj -startj + 1
	     ! if (org_loc(j,i)%x >=n_col/2) then 
		    ! org_copy(j) = org_loc(startj + j-1,i)%x - n_col
	     ! else 
		    ! org_copy(j) = org_loc(startj+j-1,i)%x
		 ! end if 
	  ! end do
   ! maxx = maxval(org_copy(1:endj-startj+1))
   ! minx = minval(org_copy(1:endj-startj+1))
   org_mid%x = int((maxx + minx + n_col )/2.0)
   else 
   print *, 'error in org_midloc'
   end if 
   
   if (org_mid%x < 1) org_mid%x = org_mid%x + n_col  
   if (org_mid%x > n_col) org_mid%x = org_mid%x - n_col  
   
   end subroutine org_midloc
   
   ! *************************************************************
   
   subroutine make_trans()
   
   use globalvariables
   implicit none 
   integer(kind=4) :: xx, yy, j, k, particle_id, y0
   real(kind=8),allocatable :: dz(:)    ! 1D z grid to be read
   real(kind=8) :: prob0, probEnd, probchk, dumreal, probchk2   
   real(kind=8),allocatable :: trans(:,:) ! transition matrix units in /yr (?)
   real(kind=8),allocatable :: part_num_t(:)
   integer(kind=4) :: layer0_max, layer0_min, layerEnd_max, layerEnd_min, part_num(n_row), nz
   Character*21 numtemp
   logical :: isparticle
   real(kind=8) :: tol = 1d-6
   integer(kind=4) :: transx(n_row, n_col)
   integer(kind=4) :: transy(n_row, n_col)
   character*256 gridfile
   real(kind=8) :: part_num_tmp
   
   write(numtemp,'(i10.1)') Time
   write(gridfile,*) '1dgrid_turbo.txt'
   ! write(gridfile,*) '1dgrid.txt'
   gridfile=trim(adjustl(gridfile))
   
   if (N_ind/=0) then 
        print*, 'Error: Org exists', N_ind
        write(File_log,*) 'Error: Org exists', N_ind
        stop
    endif
    
    transx = 0
    transy = 0
    
    if (.false.) then 
    do yy = 1, n_row
        do xx = 1, n_col
            if (isparticle(coordinates(yy,xx))) then 
                particle_id = matrix(yy,xx)%value 
                if (particle(particle_id)%loc_init%y<1 .or. &
                    particle(particle_id)%loc_init%y>N_row .or.  &
                    particle(particle_id)%loc_init%x<1 .or. &
                    particle(particle_id)%loc_init%x>N_col) then 
                    print*,'error in making_trans at ',particle(particle_id)%loc_init%y, &
                        particle(particle_id)%loc_init%x,particle_id
                    write(file_log,*)'error in making_trans at ',particle(particle_id)%loc_init%y, &
                        particle(particle_id)%loc_init%x,particle_id
                    stop
                endif
                transx(particle(particle_id)%loc_init%y, particle(particle_id)%loc_init%x) = xx  
                ! matrix denoting the x point to move at matrix point of original location
                transy(particle(particle_id)%loc_init%y, particle(particle_id)%loc_init%x) = yy  
                ! matrix denoting the y point to move at matrix point of original location
            endif
        enddo
    enddo
    
    ! print*,transx 
    ! print*,''
    ! print*,transy 
    ! stop
    
   OPEN(unit = File_temp, File = trim(adjustl(today))//'/data/trans_x-'         &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     
      do yy = 1,n_row
        write(file_temp,*) (transx(yy,xx),xx=1,n_col)
      enddo 
   close(File_temp)
    
   OPEN(unit = File_temp, File = trim(adjustl(today))//'/data/trans_y-'         &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     
      do yy = 1,n_row
        write(file_temp,*) (transy(yy,xx),xx=1,n_col)
      enddo 
   close(File_temp)
   
   endif
   
   if (allocated(dz)) deallocate(dz)
   if (allocated(trans)) deallocate(trans)
   nz = 0
   
   open(unit=File_temp, file=gridfile,action='read',status='old')
   do 
      read(file_temp,*,end=100) dumreal
      nz = nz+1
   enddo 
 100  close(File_temp)
   
   ! nz = 100
   
   print*,nz
   ! stop
   allocate(dz(nz),trans(nz,nz),part_num_t(nz))
   
   open(unit=File_temp, file=gridfile,action='read',status='old')
   do yy = 1,nz
      read(file_temp,*) dz(yy)
   enddo 
   close(File_temp)
   
   !                   -----------  layer0_min
   !          y_int ----------
   !                   -----------  layer0_max
   
   !                   -----------  layerEnd_min
   !          y ----------
   !                   -----------  layerEnd_max
   
   
   y0 = n_row  ! highest vertical point of original and current sediment particles
   part_num = 0 
   do yy = 1, n_row
     do xx = 1, n_col
       if (isparticle(coordinates(yy,xx))) then
          part_num(yy) = part_num(yy)+1
          particle_id = matrix(yy,xx)%value 
          ! if (particle(particle_id)%loc_init%y /=particle(particle_id)%loc%y) then 
              if (particle(particle_id)%loc_init%y < y0) y0 = particle(particle_id)%loc_init%y 
              if (particle(particle_id)%loc%y < y0) y0 = particle(particle_id)%loc%y 
          ! endif
          ! exit
       endif
     ! if (y0/=0) exit
     enddo
   enddo
   trans = 0d0
   part_num_t = 0d0
   do yy = 1, n_row
        do xx = 1, n_col
            if (isparticle(coordinates(yy,xx))) then
                particle_id = matrix(yy,xx)%value 
                if ((particle(particle_id)%loc_init%y - y0)==0) then 
                    layer0_min = 1
                    layer0_max = 1
                    do j = 1, nz-1
                        if ( ( real((particle(particle_id)%loc_init%y - y0 + 1),kind=8)*pixelsize - sum(dz(:j)) ) &
                            *( real((particle(particle_id)%loc_init%y - y0 + 1),kind=8)*pixelsize - sum(dz(:j+1)) ) <= 0d0 ) then 
                            layer0_max = j+1
                            exit
                        endif
                    enddo 
                else 
                    layer0_min = 1
                    layer0_max = 1
                    do j = 1, nz-1
                        if ( ( real((particle(particle_id)%loc_init%y - y0),kind=8)*pixelsize - sum(dz(:j)) ) &
                            *( real((particle(particle_id)%loc_init%y - y0),kind=8)*pixelsize - sum(dz(:j+1)) ) <= 0d0 ) then 
                            ! print*,y0,particle(particle_id)%loc_init%y,j &
                                ! ,real((particle(particle_id)%loc_init%y - y0),kind=8)*pixelsize ,sum(dz(:j)),sum(dz(:j+1))
                            layer0_min = j + 1
                            exit
                        endif
                    enddo 
                    do j = 1, nz-1
                        if ( ( real((particle(particle_id)%loc_init%y - y0+1),kind=8)*pixelsize - sum(dz(:j)) ) &
                            *( real((particle(particle_id)%loc_init%y - y0+1),kind=8)*pixelsize - sum(dz(:j+1)) ) <= 0d0 ) then 
                            layer0_max = j + 1
                            exit
                        endif
                    enddo 
                endif   
          
                if ((particle(particle_id)%loc%y - y0)==0) then 
                    layerEnd_min = 1
                    layerEnd_max = 1
                    do j = 1, nz-1
                        if ( ( real((particle(particle_id)%loc%y - y0 + 1),kind=8)*pixelsize - sum(dz(:j)) ) &
                            *( real((particle(particle_id)%loc%y - y0 + 1),kind=8)*pixelsize - sum(dz(:j+1)) ) <= 0d0 ) then 
                            layerEnd_max = j+1
                            exit
                        endif
                    enddo 
                else 
                    layerEnd_min = 1
                    layerEnd_max = 1
                    do j = 1, nz-1
                        if ( ( real((particle(particle_id)%loc%y - y0),kind=8)*pixelsize - sum(dz(:j)) ) &
                            *( real((particle(particle_id)%loc%y - y0),kind=8)*pixelsize - sum(dz(:j+1)) ) <= 0d0 ) then 
                            layerEnd_min = j + 1
                            exit
                        endif
                    enddo 
                    do j = 1, nz-1
                        if ( ( real((particle(particle_id)%loc%y - y0+1),kind=8)*pixelsize - sum(dz(:j)) ) &
                            *( real((particle(particle_id)%loc%y - y0+1),kind=8)*pixelsize - sum(dz(:j+1)) ) <= 0d0 ) then 
                            layerEnd_max = j + 1
                            exit
                        endif
                    enddo 
                endif   
                
                print*, y0, pixelsize,particle(particle_id)%loc_init%y, particle(particle_id)%loc%y, layer0_min, layer0_max  &
                    , layerEnd_min, layerEnd_max
          
                probchk = 0d0
                probchk2 = 0d0
          
                do j = layer0_min, layer0_max
                    if (j==1) then 
                        prob0 = min(sum(dz(:j)),real((particle(particle_id)%loc_init%y - y0+1),kind=8)*pixelsize )      &
                            -max(0d0, real((particle(particle_id)%loc_init%y - y0),kind=8)*pixelsize)
                    else 
                        prob0 = min(sum(dz(:j)),real((particle(particle_id)%loc_init%y - y0+1),kind=8)*pixelsize )      &
                            -max(sum(dz(:j-1)), real((particle(particle_id)%loc_init%y - y0),kind=8)*pixelsize)
                    endif
              
                    prob0 = prob0/pixelsize
                    
                    print*, j, prob0
              
                    probchk2 = probchk2 + prob0
                        
                    part_num_t(j) = part_num_t(j) + prob0
                                        
                    do k = layerEnd_min, layerEnd_max
              
                        if (k==1) then 
                            probEnd = min(sum(dz(:k)),real((particle(particle_id)%loc%y - y0+1),kind=8)*pixelsize )      &
                                -max(0d0, real((particle(particle_id)%loc%y - y0),kind=8)*pixelsize)
                        else 
                            probEnd = min(sum(dz(:k)),(particle(particle_id)%loc%y - y0+1)*pixelsize )      &
                                -max(sum(dz(:k-1)), real((particle(particle_id)%loc%y - y0),kind=8)*pixelsize)
                        endif
                        
                        probEnd = probEnd/pixelsize
              
                        probchk = probchk + prob0*probEnd
              
                        ! trans(j,k) = trans(j,k) + prob0*probEnd/part_num(yy)
                        ! trans(j,j) = trans(j,j) - prob0*probEnd/part_num(yy)
                        trans(j,k) = trans(j,k) + prob0*probEnd
                        trans(j,j) = trans(j,j) - prob0*probEnd
                    enddo
                enddo 
          
                ! if (nint(probchk)/=1) then 
                if (abs(probchk-1d0)> tol) then 
                    print*,'error in making transition matrix', probchk, xx, yy
                    write(File_log,*) 'error in making transition matrix', probchk, xx, yy
                    stop
                endif 
                if (abs(probchk2-1d0)> tol) then 
                    print*,'error2 in making transition matrix', probchk2, xx, yy
                    write(File_log,*) 'error2 in making transition matrix', probchk2, xx, yy
                    stop
                endif
            endif
        enddo 
    enddo 
    
   OPEN(unit = File_temp, File = trim(adjustl(today))//'/data/transmtx-'         &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     
      do yy = 1,nz
        if (part_num_t(yy)==0d0) then
            part_num_tmp = 1d10
        else 
            part_num_tmp = part_num_t(yy)
        endif
        ! part_num_tmp = part_num_t(yy)
        ! write(file_temp,*) (trans(yy,xx)/(Time*Timescale/365d0),xx=1,nz)
        ! write(file_temp,*) (trans(yy,xx)/(Time*Timescale/365d0)/part_num_tmp,xx=1,nz)
        write(file_temp,*) (trans(yy,xx)/part_num_tmp,xx=1,nz)
      enddo 
   close(File_temp)
   
   ! stop
   
   end subroutine make_trans
   
   ! *************************************************************

   Subroutine Output_stats()

   Use GlobalVariables
   use ieee_arithmetic
   use mod_pfit

   Implicit None
   Logical :: IsParticle, InSediments, IfNAN

   integer(kind=4) :: X, Y, Window, N_Samples, j, k, m, n
   integer(kind=4) :: part_rowsum, Particle_ID

   real(kind=8) :: Porosity_local, por_mean, por_stderr, por_rowsum
   real(kind=8) :: URand(2), Depth, Db, minvalue
   real(kind=8) :: Activity_New, Activity_mean, act_rowsum

   Character*8 :: Timevalue

   integer(kind=4) :: lab_rowsum(N_LabilityClasses), Particle_Lability, P_count
   real(kind=8) :: Lab_real_New
   real(kind=8) :: Sum_Dy, dy, Sum_dx, dx, SUM_RMS, RMS, DEV_X, DEV_Y
   real(kind=8) :: Sum_Dy2(N_Row), Sum_dx2(N_Row), SUM_RMS2(N_Row), RMS2(N_Row), DEV_X2(N_Row), DEV_Y2(N_Row), P_count2(N_Row)
   
   real(kind=8) :: Sum_Dy3(N_Row), Sum_dx3(N_Row),DEV_X3(N_Row), DEV_Y3(N_Row)     !! YK added
   
   real(kind=8) :: Depth_vector(N_Row), Activity_vector(N_Row), weight_vector(N_Row)
   real(kind=8) :: AMACH
   real(kind=8) :: Db_lability, Lability_vector(N_Row, N_LabilityClasses)

   integer(kind=4) NRMISS

   real(kind=8) AOV(15), CASE(N_Row,12), COEF(2,5), COVB(2,2), TESTLF(10), Profile(N_Row,3)
   real(kind=8) STAT(15,1)
   real(kind=8), allocatable :: PorosityData(:)
   
   real(kind=8) :: avx, avy, sumx2, sumy2, sumxy, slp, intcpt, crr     !!!  YK 5/22/2017
   
   real(kind=8), allocatable ::  xsel(:), ysel(:), xlsel(:), ylsel(:)   !!!  YK 5/22/2017
   
   integer(kind=4) :: cnt, yint, yfin, cnt2, cnt3, row
   integer(kind=4) ::   list(3,n_row)
   
   double precision, allocatable :: amx(:,:), amx2(:,:), ymx(:), ymx2(:), cmx(:)
   
   real(kind=8) :: Db_list(n_row),  Activity_vector2(N_Row), Db_list2(n_row), Dbx_list(n_row), Dby_list(n_row)
   
   logical :: matrixCheck, int_activity_1
   
   integer(kind=4) info
   double precision imbr
   integer(kind=4), allocatable :: ipiv(:)
   
   real(8), allocatable :: a(:), sig(:), coeff(:,:), cov(:,:)
   
   real(kind=8) :: poly
   
   integer(kind=4) :: poly_fit_n = 6
   
   real(kind=8) :: lab_rowsum_real , lab_mean
   real(kind=8) :: lab_real_vect(N_row)

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
       
       sum_DY3(y) = 0.
       sum_DX3(y) = 0.
       
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
           
           sum_dx3(y) = sum_dx3(y) + real(dx*dx)
           sum_dy3(y) = sum_dy3(y) + real(dy*dy)
           
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
       if (P_COUNT2(y) == 0) then 
       write(File_Displace2,*) Time, Depth_vector(y), 0., 0., 0.,   &
           0., 0., 0., 0.
       cycle 
       endif
       
       RMS2(y) = PixelSize *PixelSize * ( SUM_RMS2(y) / rEAL(P_COUNT2(y)) )   ! CALC OF rms .. == SQRT(MEAN(SQUARED DEVIATIONS))
       DEV_Y2(y) = PixelSize * SUM_DY2(y) / rEAL(P_COUNT2(y))
       DEV_X2(y) = PixelSize * SUM_Dx2(y) / rEAL(P_COUNT2(y))
       DEV_Y3(y) = PixelSize *PixelSize * SUM_DY3(y) / rEAL(P_COUNT2(y))
       DEV_X3(y) = PixelSize * PixelSize *SUM_Dx3(y) / rEAL(P_COUNT2(y))
       
       if (time.ne.0) then 
       Db_list(y) = RMS2(y)*365./time/timescale/4.
       Dby_list(y) = Dev_Y3(y)*365./time/timescale/4.
       Dbx_list(y) = Dev_X3(y)*365./time/timescale/4.
       else 
       Db_list(y) = 0
       Dby_list(y) = 0
       Dbx_list(y) = 0
       end if
	   
       write(File_Displace2,*) Time, Depth_vector(y), DEV_Y3(y), DEV_X3(y), RMS2(y),   &
           Db_list(y), Dby_list(y), Dbx_list(y), P_COUNT2(y)
       
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
            ymx(row) = (Db_list(y) - Db_list(y))/pixelSize/pixelSize* merge(Activity_vector2(y),1.d0      &
                       , (Activity_vector2(y).ne.0.0).and.(.not.ieee_is_nan(Activity_vector2(y)))     &
                       .and.(yint.ne.1).and.(.not.int_activity_1)) &
                     +  Db_list(y)*(-1.0)/pixelSize/pixelSize * merge(Activity_vector2(y),1.d0   &
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
              , merge(Activity_vector2(yint),1.d0,.not.int_activity_1)                     &
              *exp(-sqrt(decayconstant/Db)*(depth_vector(j+yint-1)-depth_vector(yint)))  &
              , merge(Activity_vector2(yint),1.d0,.not.int_activity_1)                         &
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
   integer(kind=4) :: i, j

   j = 1 ! the index 1 of the body used to tag organisms ...
     Do i = 1, N_Ind
       MatrixOccupancy(Org_loc(j,i)%Y) = MatrixOccupancy(Org_loc(j,i)%Y) + 1
     End do

   End Subroutine LocateOrganisms_scan

   ! ***********************

   Subroutine Chk_Org(chr)
     ! locate the organisms in the matrix and accumulate the verticle depth distribution in "MatrixOccupancy"

   Use GlobalVariables
   Implicit None
   integer(kind=4) :: i, j, xx, yy
   character*6, intent(IN) :: chr

       
     do i = 1, N_ind
         do j=1,Org(i)%BodySize
            xx= Org_loc(j,i)%X
            yy= Org_loc(j,i)%Y
            ! print*,j,i,xx,yy,Org(i)%BodySize
            if (.not.(matrix(yy,xx)%class==i .and. matrix(yy,xx)%value==j)) then 
                print*,'Error in organisms in '//chr//' ; time, i,j,xx,yy',time, i,j,xx,yy
                write(file_log,*)'Error in organisms in '//chr//' ; time, i,j,xx,yy',time, i,j,xx,yy
            endif
        enddo 
     enddo

   End Subroutine Chk_Org

   ! ***********************


   Subroutine Output_Graphic()
     ! Create images of the matrix

   Use GlobalVariables
   use gif_util

   implicit none
   integer(kind=4)    :: X, Y
   integer(kind=4), allocatable   :: image_data(:,:)
   Integer(kind=4)   :: LookupPorosityColour, LookupActivityColour, LookupLabilityColour
   
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

       call writegif(trim(adjustl(today))//"/data/lability-"          &
                //trim(adjustl(numtemp))//".gif", image_data, palette)

   End if

   If (.not. Movies_On) then
    If (PlotPorosity) then
         DO X = 1, N_Col
         DO Y = 1, N_Row
           image_data(X,Y) = LookupPorosityColour(Y,X) ! Write the stored data to the image array.
         End Do
         End Do
 
       call writegif(trim(adjustl(today))//"/data/porosity-"           &
                 //trim(adjustl(numtemp))//".gif", image_data, palette)

    End if

    If (PlotActivity) then
         DO X = 1, N_Col
         DO Y = 1, N_Row
           image_data(X,Y) = LookupActivityColour(Y,X)    ! Write the stored data to the image array.
         End Do
         End Do

       call writegif(trim(adjustl(today))//"/data/activity-"//            &
                   trim(adjustl(numtemp))//".gif", image_data, palette)

     End if
    End if

   End Subroutine Output_Graphic

   ! *******************************

   Function LookupLabilityColour(Y,X)
   ! return a colour value (lookup table)

   Use GlobalVariables
   Implicit none
   Logical   :: IsPArticle, IsWater, IsOrganism
   Integer(kind=4) :: LookupLabilityColour
   integer(kind=4)   :: Y,X

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
   Integer(kind=4) :: LookupPorosityColour
   integer(kind=4)   :: Y,X
   real(kind=8) :: Por, Delta_Porosity

   ! functions
   real(kind=8) :: Porosity_local


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
   Integer(kind=4) :: LookupActivityColour
   real(kind=8)      :: Activity
   integer(kind=4)   :: Y,X

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
   integer(kind=4) :: i,j, xx, yy, xxp, xxg, yyp, yyg, k , o
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
   Type (CoOrdinates), intent (IN) :: Point
   
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
   Logical IsMoveable, IsWater, IsPArticle, InMatrix
   Type (CoOrdinates) :: Point
   
   if (.not.inmatrix(point)) then
    print *, 'error: ismoveable function' 
    write(file_log, *) 'error: ismoveable function' 
   endif
   
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
   integer(kind=4) :: Water_Count, Total_Count, X, Y, shift
   real(kind=8) :: Porosity_local
   Type (CoOrdinates) :: Point
   integer(kind=4) :: xx,yy

     Water_Count = 0
     Total_Count = 0

     Do X = Point%X-shift, Point%X+shift
     Do Y = Point%Y-shift, Point%Y+shift
       xx = x
       yy= y
       if (xx<1 )xx =xx + n_col
       if (xx>n_col )xx =xx - n_col
       If (InMatrix(CoOrdinates(YY,XX))) then
         Total_Count = Total_Count + 1
         If (IsWater(CoOrdinates(YY,XX))) Water_count = Water_Count + 1
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
   integer(kind=4) Y
   real(kind=8) Porosity, Depth, DepthMAx, URand

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
   integer(kind=4) Y, x, j
   real(kind=8) URand, depth, depthmax, lability_pr
   real(kind=8) :: lab_real

   Depth = Real(Y - N_RowWAter) * PixelSize
   DepthMAx = Real(N_Row - N_RowWAter) *  PixelSize
   ! Lability_pr = 1. - Real(Depth) / DepthMAx
   Lability_pr = 1. - 0.9*Real(Depth) / DepthMAx  ! some organic matter can remain
   ! Lability_pr = 1.   ! random distribution

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
   integer(kind=4) :: ID
   real(kind=8) :: Lab_real, Lab_real_New
   real(kind=8)    :: URand, Pr_Cummulative, TimeIncrement
   real(kind=8) :: lab_decay
   real(kind=8) :: fact_law1, fact_law2
   real(kind=8) merge2
   
   logical :: inMatrix, IsParticle
   
   select case (trim(adjustl(O2ratelaw)))
     case('linear','LINEAR','Linear')
        fact_law1 = 0d0
        fact_law2 = 1d0
     case('zero','ZERO','Zero')
        fact_law1 = 0d0
        fact_law2 = 0d0
     case('monod','MONOD','Monod')
        fact_law1 = 1d0
        fact_law2 = 0d0
   end select 
   

     Lab_real = Particle(ID)%OM%OMact
     Lab_real_new = Particle(ID)%OM%OMact
     
       if (Lab_real .le. 0.01) Then
         Lab_real_New = 0
         Return
       END IF
	   
           TimeIncrement = Real(Time - Particle(ID)%OM%OM_Time0,kind=8) * TimeScale / 365.d0   !  (yr)
           
           Pr_Cummulative = kdcy*TimeIncrement*Lab_real 
    
    ! if (oxygen_ON) then
    
    if (particle(ID)%loc%X == 1) then 
      if (particle(ID)%loc%Y == n_row) then 
    O2(particle(ID)%loc%Y,n_col)%Oxygen_use = &
            O2(particle(ID)%loc%Y,n_col)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)   &
             * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,n_col)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
            O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)   &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
            O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.) &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    Lab_real_new = Lab_real -( &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,n_col)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,n_col)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))    &
             )*iox/OM_uni
      else if (particle(ID)%loc%Y == 1) then 
    O2(particle(ID)%loc%Y,n_col)%Oxygen_use = &
            O2(particle(ID)%loc%Y,n_col)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.) &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,n_col)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
            O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)  &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use = &
            O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)   &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    Lab_real_new = Lab_real -( &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,n_col)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,n_col)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))    &
             )*iox/OM_uni
       else 
    O2(particle(ID)%loc%Y,n_col)%Oxygen_use = &
             O2(particle(ID)%loc%Y,n_col)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)   &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,n_col)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.) &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)  &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.) &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    Lab_real_new = Lab_real -( &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,n_col)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,n_col)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))    &
             ) *iox/OM_uni
        end if 
        
   else if (particle(ID)%loc%X == n_col) then 
      if (particle(ID)%loc%Y == n_row) then 
    O2(particle(ID)%loc%Y,1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)  &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)  &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)  &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    Lab_real_new = Lab_real - ( &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))    &
             )*iox/OM_uni
      else if (particle(ID)%loc%Y == 1) then 
    O2(particle(ID)%loc%Y,1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)  &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.) &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)  &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    Lab_real_new = Lab_real - ( &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))    &
             )*iox/OM_uni
       else 
    O2(particle(ID)%loc%Y,1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.) &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.) &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.) &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.) &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    Lab_real_new = Lab_real -( &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))    &
             )*iox/OM_uni
      end if
      
    else 
      if (particle(ID)%loc%Y == n_row) then 
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)  &
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)&
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)&
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    Lab_real_new = Lab_real - ( &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))    &
             )*iox/OM_uni
      else if (particle(ID)%loc%Y == 1) then 
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)&
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)&
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)&
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    Lab_real_new = Lab_real - ( &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))    &
             )*iox/OM_uni
       else 
    O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)&
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use = &
             O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)&
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)&
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use = &
             O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen_use + &
             Pr_Cummulative/(TimeScale/365.)&
      * merge(1d0/mo2,1d0,fact_law1 == 1d0 .and. O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%Oxygen<=mo2/iox)
             
    Lab_real_new = Lab_real -( &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X+1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y,particle(ID)%loc%X-1)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y-1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))  +  &
             Pr_Cummulative*(merge2(O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen*(fact_law1/mo2+fact_law2)  &
                ,1.d0,O2(particle(ID)%loc%Y+1,particle(ID)%loc%X)%oxygen<= fact_Law1*mo2/iox + fact_law2*10d0))    &
             )*iox/OM_uni
      end if
    end if
   
    ! in units of /yr
    
             Particle(ID)%OM%OM_Time0 = Time
             Particle(ID)%OM%OMact = Lab_real_New
    
     ! end if
     
     ! if (.not. oxygen_ON) then 
          
          ! Lab_real_new = Lab_real -( &
             ! Pr_Cummulative  & 
             ! )*iox/OM_uni
             
             ! Particle(ID)%OM%OM_Time0 = Time
             ! Particle(ID)%OM%OMact = Lab_real_New
     
     ! end if 
           
   End Function Lab_real_New

   !****************************************************************

   Function Activity_New(n)
     ! returns the new activity of particle n after decay

   Use GlobalVariables
   implicit none
   integer(kind=4) n
   real(kind=8) TimeIncrement, Activity_New

     TimeIncrement = Real (Time - Particle(n)%Pb%Time0) * TimeScale / 365.
     Activity_New = Particle(n)%Pb%Activity0 * (Exp( - DecayConstant * TimeIncrement ) )
     Particle(n)%Pb%Activity = Activity_New

   End Function Activity_New

   ! ************************************

   Function CurrentTimeString()
     ! returns the current time as a string character

   Use GlobalVariables
   implicit none
   real(kind=8) :: YearValue, DayValue, HourValue, MinValue
   CHARACTER*4 Year_str
   CHARACTER*3 day_str
   Character*2 Hour, Minutes
   Character*15 CurrentTimeString

     YearValue = Real(Time) * TimeScale / 365.   ! no years elapsed
     DayValue = 365. * (YearValue - Real(INT(YearValue)))
     HourValue = 24. * (DayValue - Real(Int(DayValue)))
     MinValue = 60. * (HourValue - Real(Int(HourValue)))

     write(Year_str, '(I4.4)') Int(YearValue)
     write(day_str, '(I3.3)') Int(DayValue)
     write(Hour, '(I2.2)') Int(HourValue)
     write(Minutes, '(I2.2)') Int(MinValue)

     CurrentTimeString(1:1) = 'Y'
     CurrentTimeString(2:5) = Year_str
     CurrentTimeString(6:6) = 'D'
     CurrentTimeString(7:9) = day_str
     CurrentTimeString(10:10) = 'H'
     CurrentTimeString(11:12) = Hour
     CurrentTimeString(13:13) = 'M'
     CurrentTimeString(14:15) = Minutes

   End Function CurrentTimeString

  ! **********************

   Function P_Count()
     ! count the number of particles .. used to debug, not functional in the simulation

   Use GlobalVariables
   implicit none
   Logical IsParticle
   integer(kind=4)  i, j, x, y, P_Count

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
   
   integer(kind=4) :: n, i 
   double precision :: a(n)
   real(kind=8) x, clc, poly
   
   clc = 0.
   do i = 1, n
      clc = clc + a(i)*x**(i-1)
   end do 
   
   poly = clc
   
   end function poly

   ! ***************************
   
   function merge2(a,b,n)
   
   implicit none
   
   logical :: n 
   real(kind=8) :: a,b,merge2
   
   if (n) then
     merge2 = a
   else 
     merge2 = b
   endif
   
   end function merge2
   
   ! *************************
   
   
   
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
integer(kind=4) n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer(kind=4) i, j, k

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
integer(kind=4) n
double precision a(n,n), c(n), b(n)
! double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer(kind=4) i, k

do k=1, n
   coeff = 0.0
   do i=1,n
      coeff = coeff + a(k,i)*c(i)
   end do
   b(k) = coeff
end do

end subroutine dot


! ====================================================
subroutine heapsortR(n,array,turn)
!!!  from http://slpr.sakura.ne.jp/qp/sortf90/
  implicit none
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(out)::turn(1:n)
  real(kind=8),intent(inout)::array(1:n)
 
  integer(kind=4)::i,k,j,l,m
  real(kind=8)::t
 
  if(n.le.0)then
     write(6,*)"Error, at heapsort"; stop
  endif
  if(n.eq.1)return

  do i=1,N
     turn(i)=i
  enddo

  l=n/2+1
  k=n
  do while(k.ne.1)
     if(l.gt.1)then
        l=l-1
        t=array(l)
        m=turn(l)
     else
        t=array(k)
        m=turn(k)
        array(k)=array(1)
        turn(k)=turn(1)
        k=k-1
        if(k.eq.1) then
           array(1)=t
           turn(1)=m
           exit
        endif
     endif
     i=l
     j=l+l
     do while(j.le.k)
        if(j.lt.k)then
           if(array(j).lt.array(j+1))j=j+1
        endif
        if (t.lt.array(j))then
           array(i)=array(j)
           turn(i)=turn(j)
           i=j
           j=j+j
        else
           j=k+1
        endif
     enddo
     array(i)=t
     turn(i)=m
  enddo

  return
end subroutine heapsortR
! =====================================================