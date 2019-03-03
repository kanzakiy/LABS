   ! ***********************************************************************
   !         Program name  AUTOMATONS
   !
   !         Investigator  Bernard P. Boudreau
   !         Programmers   Frederique FRANCOIS, Jae S. Choi
   !         Locale        Halifax, Nova Scotia, Canada
   !         Institution   Department of Oceanography, Dalhousie University
   !
   !         Language      Fortran 90
   !         Libraries     IMSL numerical libs for data analysis
   !
   !         Date modified 19 July 2000
   ! ***********************************************************************

   Module GlobalVariables
     ! this module defines the global variables open to all subroutines

   Implicit None

   Logical, Save :: Lability_ON, Block_outsideMatrix, LocalMixing
   Logical, Save :: Sedimentation_ON, Movies_on, DepthAvoidance_on
   Logical, Save :: SaveData, PlotMatrix, PlotPorosity, PlotActivity,OutputScaling
   !  added parameters to eLABS
   logical, save :: oxygen_ON, oxFB_ON, Resp_ON , only_sed ,errDetect ,errChk , I_shape 
   logical, save :: flow_ON, Detail_Log ,Ash_ON  ,Incl_ASH, Mod_Pop,non_Pop,trans_making 
   logical, save :: Long_run,  O2lim_on 

   integer(kind=4), Save :: N_Ind, N_Row, N_Col, N_Cell, Buffer_Zone, N_RowWater, N_RowSed, N_LabilityClasses
   integer(kind=4), Save :: Total_N_Particles, Total_N_Particles0, ParticleTolerance
   integer(kind=4), Save :: StoppedTolerance, N_Outputs, DepthDeadZone, N_outputs2
   integer(kind=4), Save :: RainCount, PointerToFreeParticle, ParticlesToSediment_y, Time_Sediment
   integer(kind=4), Save :: TimeMin, TimeMax, TimeStep, RandSeed, DEPTHTOCONSTRAIN, MissingValue, WindowSize
   integer(kind=4), Save :: Time, Day, Year
   !  added 
   integer(kind=4), Save :: Savetime ,y_int, Org_ID_ishape , PopTot, Time_addpop, Time_dif_st, Time_dif_rec       

   real(kind=8), Save    :: TimeScale, PixelSize, Porosity0, POROSITY_THRESHOLD, Porosity_DecayRate, Tolerance
   real(kind=8), Save    :: SedimentationRate, SedRate, DecayConstant, AspectRatio
   real(kind=8), Save    :: Pr_Activity_day, Pr_Activity_year
   real(kind=8), Save    :: AREA_Total
   ! added 
   real(kind=8), Save    :: TotOrgDecay, TotResp, TotO2Dif, TotAbio, TotO2Adv   
   real(kind=8), save    :: totO2, pretotO2, reso2
   real(kind=8), save    :: Time_Ash , Ash_thickness , Ash_porosity    
   real(kind=8), save    :: pal ,iox , OM_uni ,fact , bio_fact , dif_0 ,mo2 ,kdcy ,width_3d ,shearfact , corg_0 

   Character*3   :: Porosity_Type
   Character*14  :: CurrentTime
   Character*255  :: WorkDir, Workname, Today, O2ratelaw

   integer(kind=4), Parameter :: w=0, p=-1   ! values of water (w) and particles (p) used in Particle%Class (see below)
   real(kind=8), Parameter :: pi = 3.14159

   Type CellContent
     integer(kind=4) :: Value, Class  ! the id value itself .. of the particle, organism, or water
   End Type CellContent       ! the class identifier - particle=p (=-1), water=w (=0), organism id

   Type Coordinates
     integer(kind=4) :: Y, X
   End Type Coordinates

   Type Tracer
     real(kind=8)    :: Activity, Activity0
     integer(kind=4) :: Time0
   End Type Tracer

   Type Tracer_OM
     real(kind=8)    :: OMact, OMact_0
     integer(kind=4) :: OM_Time0
   End Type Tracer_OM
   
   Type Tracer_Ash
     real(kind=8)    :: Ashact, Ashact_0
     integer(kind=4) :: Ash_Time0
   End Type Tracer_Ash

   Type Particle_Char
     integer(kind=4)           :: Lability, Lability_Time0, Plane
     Type(Tracer)      :: Pb
     Type(Tracer_OM)      :: OM
     ! Type(Tracer_Ash)      :: Ash
     Type(Coordinates) :: Loc, Loc_Init, Loc_prev
   End Type Particle_Char
   
   Type Particle_Char2
     Type(Tracer_Ash)      :: Ash
     Type(Coordinates) :: Loc, Loc_Init, Loc_prev
   End Type Particle_Char2

   Type Intake
     real(kind=8)    :: Rate, Selectivity
     integer(kind=4) :: Amount
   End Type Intake

   Type Gut_Char
     integer(kind=4) :: Capacity, Content
     real(kind=8)    :: Pb_Activity
   End type Gut_Char

   Type LocoMotory
     real(kind=8)        :: Rate
     integer(kind=4)     :: Distance, Bearing_Distance, rEVERSEcOUNT, StoppedTime
     real(kind=8)        :: Bearing
   End Type Locomotory

   Type CurrentState
     Logical     :: WrappingAround, OnTheEdge, Reversing
   End Type CurrentState

   Type PotentialState
     Logical     :: Egest, Ingest, Move
   End Type PotentialState

   Type Organism_Char
     real(kind=8)                :: Density, Pr_GoStraight
     Character (Len=1)   :: FeedingType
     integer(kind=4)             :: Length, Width, HeadSize, BodySize, Orientation, SensoryRange, FaecesWidth
     Type(PotentialState):: Can
     Type(CurrentState)  :: Is
     Type(Gut_Char)      :: Gut
     Type(Intake)        :: Ingest
     Type(Locomotory)    :: Move
   End Type Organism_Char
   
   type oxygen_conc
     real(kind=8) :: oxygen, oxygen_pre, oxygen_use
     integer(kind=4) :: value_pre, mark
   end type oxygen_conc

   Type(CellContent), Allocatable, Save    :: Matrix(:,:), Guts(:,:), Block(:)
   Type(Coordinates), Allocatable, Save    :: Org_Loc(:,:), Block_loc(:), Tail(:)
   Type(Particle_Char), Allocatable, Save  :: Particle(:)
   Type(Organism_Char), Allocatable, Save  :: Org(:)
   Type(Particle_Char2), Allocatable, Save  :: Particle2(:)
   ! added
   Type(oxygen_conc), Allocatable, Save    :: O2(:,:)  
   Type(CellContent), Allocatable, Save    :: Matrix_chk(:,:)  
   Type(coordinates), Allocatable, Save    :: flow_loc(:,:)  


   real(kind=8), Allocatable, Save                 :: Lability_decayConstant(:), Lability_proportion(:)
   integer(kind=4), Allocatable, Save              :: Particle_ID_free(:), MatrixOccupancy(:), Time_output(:), Time_output2(:)
   integer(kind=4), Allocatable, Save              :: IngestionHistory(:,:), MovementHistory(:,:)
   ! added 
   integer(kind=4), Allocatable, Save              :: RespHistory(:,:), EgestHistory(:,:), Dir_rec(:,:) , swi(:)      
   integer(kind=4), allocatable, save              :: DeathFlg(:), PopLogID(:)    
   real(kind=8), Allocatable, Save                 :: Ub(:,:),Vb(:,:), Ug(:,:), Vg(:,:), Pg(:,:), Dg(:,:)  
   real(kind=8), allocatable, save                 :: Uo(:,:), Vo(:,:), edif(:,:) 
   real(kind=8), allocatable, save                 :: EnergyHistory(:,:), CurrentEnergy(:), RespCnst(:)

   ! global colour definitions

   integer(kind=4), parameter :: black=0,  red=1,    green=2,  blue=3, white=4,  yellow=5,   &
                               green1=6, green2=7, green3=8, green4=9, green5=10,   &
                               blue1=11, blue2=12, blue3=13, blue4=14, blue5=15,     &
                               red1=16,  red2=17,  red3=18,  red4=19, red5=20,      &
                               turq1=21, turq2=22, turq3=23, turq4=24, turq5=25,     &
                               grey1=26, grey2=27, grey3=28, grey4=29, grey5=30, grey6=31

   integer(kind=4), Parameter :: File_Profile = 10  , File_Scaling = 20 , File_Activity=30, File_Displace=40
   integer(kind=4), Parameter :: File_Parameters = 1, File_Organisms = 2, File_ASCII = 3, File_txtImg = 4, File_txtImg_2 = 5
   
   integer(kind=4), Parameter :: File_Profile2 = 50, File_Profile3 = 60, File_Displace2=70, File_oneD=80,  File_Profile_Re = 90
   integer(kind=4), parameter :: File_temp = 950, File_flux = 230, File_Log = 340, file_flux2 = 220
   
   integer(kind=4), Parameter :: poly_fit = 100
   
   integer(kind=4), Parameter :: File_Diet = 19, File_Core = 29, File_Core_M = 39, File_Core_L = 49, File_Core_A = 59, File_pop = 55
   integer(kind=4),parameter  :: File_sedrate = 18, File_profile_st = 88, File_dbio_st = 77
   
   
   !  Summary of the structure of the User defined types
   !
   !  Cells / Sediment Matrix ***
   !  col-X or column index of of the matrix
   !  row-Y or row index of of the matrix
   !  Matrix(row,col)%Value -the unique ID number of each particle or Organism index of the bdy Part
   !          in the matrix (== ID, or "PNber") -- originally, Npar,
   !          & the unique %Organism

   !  Particles ***
   !  ID-a unique particle identification number for all particles (index) -- originally, Pnber
   !              - starting from the (total number of particles+1) any below are reserved for organism id's
   !  Particle(ID)%Lability-the lability class of that particle (0 to .. Lab_max) -- originally, Apar
   !  Particle(ID)%Pb%Activity0-the initial activity of Pb210 at the time of entry to the matrix (range?) --
   !  Particle(ID)%Pb%Time0-the initial time of entry into the sediment matrix (range?) --
   !  Particle(ID)%Loc%X-the x or column index of the particle (0 if not in the matrix) -- originally, CPar
   !  Particle(ID)%Loc%Y-the y or row index of the particle (0 if not in the matrix) -- originally, RPar

   !  Organisms ***
   !  i = organism ID number .. total (maximum) possible == N_Ind
   !  Org(i)%Can%Move-is the organism active in the matrix or on "holidays" (.t., .f.) -- originally, Energy
   !  Org(i)%FeedingType-the feeding mode of organism n (Deposit, Predator, Filter) -- originally OTyp
   !  Org(i)%Ingest%Rate-"rate" of Ingest, range (0 to 1?) -- originally, ORing
   !  Org(i)%Ingest%Selectivity-"Selectivity" of Ingest, "range (0 to 1) -- OSing
   !  Org(i)%Ingest%Period-time period between ingestions, "range (0 to 1) --
   !  Org(i)%Gut%Capacity- the max number of particles ingested, "range (?) -- Imax
   !  Org(i)%Move%Rate-"rate" of motion (0 to 1) -- OSped
   !  Org(i)%Width-(No. of Cells)
   !  Org(i)%Length-(No. of Cells)
   !  Org(i)%N_Segments-(No. of abdominal segments)
   !  Org(i)%Orientation-Orientation of the organism (0-top, 90-right, 180- bottom, 270-left)
   !    index J = 1 .. (r*c). and 1 is the head closest to the top right,
   !            relative to the orientation of the organism
   !  Org_loc(j,i)%X-the x or column index of the particle
   !  Org_loc(j,i)%Y-the y or row index of the particle
   !  etc ...

   End Module GlobalVariables