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
   logical, save :: oxygen_ON   !!! YK
   logical, save :: oxFB_ON   !!! YK
   logical, save :: Resp_ON   !!! YK
   logical, save :: only_sed  !!! YK
   logical, save :: errDetect !!!  YK
   logical, save :: errChk    !!!  YK
   logical, save :: I_shape    !!!  YK
   logical, save :: flow_ON    !!!  YK
   logical, save :: Detail_Log    !!!  YK
   logical, save :: Ash_ON    !!!  YK
   logical, save :: Incl_ASH    !!!  YK

   Integer, Save :: N_Ind, N_Row, N_Col, N_Cell, Buffer_Zone, N_RowWater, N_RowSed, N_LabilityClasses
   Integer, Save :: Total_N_Particles, Total_N_Particles0, ParticleTolerance
   Integer, Save :: StoppedTolerance, N_Outputs, DepthDeadZone, N_outputs2
   Integer, Save :: RainCount, PointerToFreeParticle, ParticlesToSediment_y, Time_Sediment
   Integer, Save :: TimeMin, TimeMax, TimeStep, RandSeed, DEPTHTOCONSTRAIN, MissingValue, WindowSize
   Integer, Save :: Time, Day, Year
   Integer, Save :: Savetime  !  added -------- YK 5/22/2017
   integer, save :: y_int       !! YK
   integer, save :: Org_ID_ishape       !! YK

   Real, Save    :: TimeScale, PixelSize, Porosity0, POROSITY_THRESHOLD, Porosity_DecayRate, Tolerance
   Real, Save    :: SedimentationRate, SedRate, DecayConstant, AspectRatio
   Real, Save    :: Pr_Activity_day, Pr_Activity_year
   Real, Save    :: AREA_Total
   Real, Save    :: TotOrgDecay, TotResp, TotO2Dif, TotAbio, TotO2Adv   !! YK
   real, save    :: totO2, pretotO2
   real, save    :: Time_Ash       !! YK
   real, save    :: Ash_thickness       !! YK
   real, save    :: Ash_porosity       !! YK

   Character*3   :: Porosity_Type
   Character*14  :: CurrentTime
   Character*255  :: WorkDir, Workname, Today, O2ratelaw

   Integer, Parameter :: w=0, p=-1   ! values of water (w) and particles (p) used in Particle%Class (see below)
   Real, Parameter :: pi = 3.14159

   Type CellContent
     Integer :: Value, Class  ! the id value itself .. of the particle, organism, or water
   End Type CellContent       ! the class identifier - particle=p (=-1), water=w (=0), organism id

   Type Coordinates
     Integer :: Y, X
   End Type Coordinates

   Type Tracer
     Real    :: Activity, Activity0
     Integer :: Time0
   End Type Tracer

   Type Tracer_OM
     Real    :: OMact, OMact_0
     Integer :: OM_Time0
   End Type Tracer_OM
   
   Type Tracer_Ash
     Real    :: Ashact, Ashact_0
     Integer :: Ash_Time0
   End Type Tracer_Ash

   Type Particle_Char
     Integer           :: Lability, Lability_Time0, Plane
     Type(Tracer)      :: Pb
     Type(Tracer_OM)      :: OM
     ! Type(Tracer_Ash)      :: Ash
     Type(Coordinates) :: Loc, Loc_Init
   End Type Particle_Char
   
   Type Particle_Char2
     Type(Tracer_Ash)      :: Ash
     Type(Coordinates) :: Loc, Loc_Init
   End Type Particle_Char2

   Type Intake
     Real    :: Rate, Selectivity
     Integer :: Amount
   End Type Intake

   Type Gut_Char
     Integer :: Capacity, Content
     Real    :: Pb_Activity
   End type Gut_Char

   Type LocoMotory
     Real        :: Rate
     Integer     :: Distance, Bearing_Distance, rEVERSEcOUNT, StoppedTime
     Real        :: Bearing
   End Type Locomotory

   Type CurrentState
     Logical     :: WrappingAround, OnTheEdge, Reversing
   End Type CurrentState

   Type PotentialState
     Logical     :: Egest, Ingest, Move
   End Type PotentialState

   Type Organism_Char
     Real                :: Density, Pr_GoStraight
     Character (Len=1)   :: FeedingType
     Integer             :: Length, Width, HeadSize, BodySize, Orientation, SensoryRange, FaecesWidth
     Type(PotentialState):: Can
     Type(CurrentState)  :: Is
     Type(Gut_Char)      :: Gut
     Type(Intake)        :: Ingest
     Type(Locomotory)    :: Move
   End Type Organism_Char
   
   type oxygen_conc
     real :: oxygen, oxygen_pre, oxygen_use
     integer :: value_pre, mark
   end type oxygen_conc

   Type(CellContent), Allocatable, Save    :: Matrix(:,:), Guts(:,:), Block(:)
   Type(Coordinates), Allocatable, Save    :: Org_Loc(:,:), Block_loc(:), Tail(:)
   Type(Particle_Char), Allocatable, Save  :: Particle(:)
   Type(Organism_Char), Allocatable, Save  :: Org(:)
   Type(Particle_Char2), Allocatable, Save  :: Particle2(:)
   
   Type(oxygen_conc), Allocatable, Save    :: O2(:,:)  !! YK
   Type(CellContent), Allocatable, Save    :: Matrix_chk(:,:)  !! YK
   Type(coordinates), Allocatable, Save    :: flow_loc(:,:)  !! YK


   Real, Allocatable, Save                 :: Lability_decayConstant(:), Lability_proportion(:)
   Integer, Allocatable, Save              :: Particle_ID_free(:), MatrixOccupancy(:), Time_output(:), Time_output2(:)
   Integer, Allocatable, Save              :: IngestionHistory(:,:), MovementHistory(:,:)
   Integer, Allocatable, Save              :: RespHistory(:,:)   !! YK 
   Integer, Allocatable, Save              :: EgestHistory(:,:)   !! YK 
   Integer, Allocatable, Save              :: Dir_rec(:,:)   !! YK 
   real, Allocatable, Save                 :: Ub(:,:),Vb(:,:)   !! YK 
   real, allocatable, save                 :: Ug(:,:), Vg(:,:), Pg(:,:), Dg(:,:)  !!  YK matrix recording velocity info
   real, allocatable, save                 :: Uo(:,:), Vo(:,:)  !!  YK matrix recording velocity info for oxygen concentration calc
   real, allocatable, save                 :: edif(:,:) 

   ! global colour definitions

   integer, parameter :: black=0,  red=1,    green=2,  blue=3, white=4,  yellow=5,   &
                               green1=6, green2=7, green3=8, green4=9, green5=10,   &
                               blue1=11, blue2=12, blue3=13, blue4=14, blue5=15,     &
                               red1=16,  red2=17,  red3=18,  red4=19, red5=20,      &
                               turq1=21, turq2=22, turq3=23, turq4=24, turq5=25,     &
                               grey1=26, grey2=27, grey3=28, grey4=29, grey5=30, grey6=31

   Integer, Parameter :: File_Profile = 10  , File_Scaling = 20 , File_Activity=30, File_Displace=40
   Integer, Parameter :: File_Parameters = 1, File_Organisms = 2, File_ASCII = 3, File_txtImg = 4, File_txtImg_2 = 5
   
   Integer, Parameter :: File_Profile2 = 50, File_Profile3 = 60, File_Displace2=70, File_oneD=80,  File_Profile_Re = 90
   integer, parameter :: File_temp = 950, File_flux = 230, File_Log = 340, file_flux_txt = 220
   
   Integer, Parameter :: poly_fit = 100
   
   Integer, Parameter :: File_Diet = 19, File_Core = 29, File_Core_M = 39, File_Core_L = 49, File_Core_A = 59
   
   real, parameter :: po_particle = 0.9
   real, parameter :: iox = 220.e-6  ! mol/L
   real, parameter :: OM_uni = 0.833 ! mol/L/wt%: conversion of wt% OM to mol/L OM
   real, parameter :: fact = 1e0/220.e-6   ! YK
   real, parameter :: bio_fact = 1e4  
   real, parameter :: dif_0 = 387.8928 ! cm^2/yr
   real, parameter :: mo2 = 8e-6 ! mol/L
   
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