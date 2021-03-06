module NS_MAC_2D
!  - - - - - - - - - - - - - - - - - 
!    based on NS_MAC_2D.f90  
!  - - - - - - - - - - - - - - - - -
use GlobalVariables
implicit none
contains 
!===============================================
subroutine flow_calc()

use GlobalVariables
implicit none
integer(kind=4) :: xx, yy, cnt, cnt2, xp, xg, yp, yg
integer(kind=4) :: nx , ny

integer(kind=4), allocatable  :: B(:,:)
integer(kind=4), allocatable  :: UP(:,:) 
integer(kind=4), allocatable  :: DOWN(:,:) 
integer(kind=4), allocatable  :: CONN(:,:) 
integer(kind=4), allocatable  :: B2(:,:) 
integer(kind=4), allocatable  :: Tp(:,:) 
integer(kind=4), allocatable  :: Tp2(:,:) 

integer(kind=4) :: np, pp, ppp
integer(kind=4) :: dirc
integer(kind=4) :: ntc  ! number of total cell 
integer(kind=4) :: step
integer(kind=4) :: nstep = 10

integer(kind=4), allocatable :: Cloc(:,:) ! location of cell associated with cell number
integer(kind=4), allocatable :: Cnum(:,:) ! cell number at (x,y)
integer(kind=4), allocatable :: CP(:) ! cell number where P is calculated (referred to P-cell)
integer(kind=4), allocatable :: CPi(:) ! P-cell number assigned to a cell
integer(kind=4), allocatable :: CDN(:,:) ! cell numbers of neighbors of a P-cell 
integer(kind=4), allocatable :: CDNt(:,:) ! cell type of neighbors of a P-cell

logical :: flow_dir_x 
logical :: choice_done

real(kind=8)  :: rdm(5)
real(kind=8)  :: Ptop = 1d0  ! in g/cm/s/s (*** note: Pa = 10 g/cm/s/s and 1 atm = 1e5 Pa, i.e., P = 1 (g/cm/s/s) corresponds 0.1 Pa and 1 uatm )
real(kind=8)  :: Pbot = 1d0
real(kind=8)  :: VelC          ! in cm/s
real(kind=8)  :: new = 0.015d0 ! viscosity in g/cm/s
real(kind=8)  :: delh   ! cm/pixel 
real(kind=8)  :: delt   ! s
real(kind=8)  :: delx   ! cm/pixel 
real(kind=8)  :: dely   ! cm/pixel 

real(kind=8) ,allocatable :: Pc(:,:)
real(kind=8) ,allocatable :: Uc(:,:)
real(kind=8) ,allocatable :: Vc(:,:)
real(kind=8) ,allocatable :: Um(:,:)
real(kind=8) ,allocatable :: Vm(:,:)
real(kind=8) ,allocatable :: Dc(:,:)
real(kind=8) ,allocatable :: Fc(:,:)
real(kind=8) ,allocatable :: Gc(:,:)
real(kind=8) ,allocatable :: Rc(:,:)
real(kind=8)  :: Ut_l, Ut_r, Vt_t, Vt_b

integer ( kind = 4 ) :: n, nnz
integer ( kind = 4 ), allocatable :: ai(:)  ! row number where /=0
integer ( kind = 4 ), allocatable  :: ap(:)   ! number of non-zero component at each column
real ( kind = 8 ), allocatable  :: ax(:)   !  conponent where /=0
real ( kind = 8 ), allocatable :: bx(:) 
real  ( kind = 8 ) :: control(20)
integer ( kind = 4 ) :: filenum
integer  ( kind = 4 ) :: i
real ( kind = 8 ) :: info(90)
integer ( kind = 8 ) :: numeric
integer ( kind = 4 ) :: status
integer ( kind = 8 ) :: symbolic
integer  ( kind = 4 ) :: sys
real ( kind = 8 ), allocatable  :: kai(:)

integer(kind=4) :: temp1(5)
real(kind=8)  :: temp2(5)
integer(kind=4) :: temp3(5)

logical :: chk_matrix = .false.
logical :: show_display = .true.
! logical :: show_display = .false.
logical :: simple_test = .false.
logical :: random_choice = .false.
logical :: p_create = .true.

real(kind=8)  :: poro = 0.80d0

logical :: const_bot = .false.
logical :: const_top = .false.
logical :: rec_binary = .false.
! logical :: rec_binary = .true.
logical :: chk_particles = .true.
logical :: flg_stop = .false.
logical :: show_display_details = .true.

character*21 numtemp, cat_tmp, adv_chr, flow_chr, loc_chr
integer(kind=4) xxx,yyy, x_tmp, y_tmp
real(kind=8) :: perm_u, perm_d
!---------------------------
write(numtemp,'(i10.1)') Time
!  allocate 
nx = n_col*flw_rsltn
ny = n_row*flw_rsltn
ntc = nx*(ny+2)
if (.not.allocated(B)) allocate(B(nx,ny))
if (.not.allocated(UP)) allocate(UP(nx,ny))
if (.not.allocated(DOWN)) allocate(DOWN(nx,ny))
if (.not.allocated(CONN)) allocate(CONN(nx,ny))

if (.not.allocated(B2)) allocate(B2(nx,ny+2))
if (.not.allocated(Tp)) allocate(Tp(nx,ny+2))
if (.not.allocated(Tp2)) allocate(Tp2(nx,ny+2))

if (.not.allocated(Cloc)) allocate(Cloc(ntc,2))
if (.not.allocated(Cnum)) allocate(Cnum(nx,ny+2))
if (.not.allocated(CPi)) allocate(CPi(ntc))

if (.not.allocated(Pc)) allocate(Pc(nx,ny+2))
if (.not.allocated(Uc)) allocate(Uc(nx,ny+2))
if (.not.allocated(Vc)) allocate(Vc(nx,ny+2))
if (.not.allocated(Um)) allocate(Um(nx+1,ny+2))
if (.not.allocated(Vm)) allocate(Vm(nx,ny+2+1))
if (.not.allocated(Dc)) allocate(Dc(nx,ny+2))
if (.not.allocated(Fc)) allocate(Fc(nx,ny+2))
if (.not.allocated(Gc)) allocate(Gc(nx,ny+2))
if (.not.allocated(Rc)) allocate(Rc(nx,ny+2))

delh = pixelsize/flw_rsltn
delx = delh
dely = delh
delt = delh/2.0d0

VelC = pixelsize/timescale/24.0d0/60.0d0/60.0d0  ! cm/s 

if (calc_perm) then 
    const_bot = .true.
    const_top = .true.
    nstep = 100
    Pbot = 0d0
    delt = delh/2.0d0/flw_rsltn
    cat_tmp = 'geo'
else 
    const_bot = .false.
    const_top = .false.
    nstep = 10
    cat_tmp = 'o2'
    if (max(maxval(abs(Vb)),maxval(abs(Ub)))*VelC > 1.) then  
        nstep =100
        delt = delt/10.
    endif         
endif

! reading file
if (.not.simple_test) then
	B = 0
	do xx = 1, nx/flw_rsltn
		do yy = 1, ny/flw_rsltn
			if (matrix(yy,xx)%class == p) then 
                if (flw_rsltn/=1) then 
                    B((xx-1)*flw_rsltn+1:(xx-1)*flw_rsltn +1+ flw_rsltn-1,(yy-1)*flw_rsltn+1:(yy-1)*flw_rsltn+1+ flw_rsltn-1) = 1
                elseif(flw_rsltn ==1) then 
                    B(xx,yy) = 1
                endif 
            endif 
		end do 
	end do 
end if 

call random_seed
call random_number(rdm)
if (simple_test) then 
	B = 0
	if (p_create) then 
		do pp = 1, int((1.0d0-poro)*nx*ny)
			xx = int(1+rdm(1)*nx)
			yy = int(1+rdm(2)*ny)
			if (B(xx,yy) ==0) B(xx,yy) = 1
	    end do
	end if 
end if 

! x100 bindary file (unnecessary?); now, 100 = obstacle, 0 = water
B = 100*B
if (calc_perm) then ! ignoring organisms
	where (B == 200)
		B = 0
	else where (B == 300)
		B = 0
	end where
end if 

!!!! check the read file
if ((rec_binary).and.(time>18700)) then 
	open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-BI-'  & 
		//trim(adjustl(numtemp))//'.txt', status = 'replace')
	do yy = 1,ny
		write(100,*) (B(xx,yy), xx = 1,nx)
	end do 
	close(100)
end if 

! check connections (the following choice may not necessary if B is transformed instead?)
flow_dir_x = .true.
flow_dir_x = .false.  
UP = 0; DOWN = 0
cnt = 100
cnt2 = 0
! in case flow direction is y
if (.not.flow_dir_x) then 
	! check top to bottom connection
    do while (cnt > 0) 
        cnt2 = cnt2 + 1
        cnt = 0
        do yy = 1, ny
            do xx = 1, nx
                if (B(xx,yy) /= 0) cycle ! obstacles
                if (yy == 1) UP(xx,yy) = 1
                if (yy /= 1) then
                    if (UP(xx,yy) == 1) cycle 
                    if (UP(xx,yy-1) == 1) UP(xx,yy) = 1
                    if (xx == 1) then
                        if (UP(2,yy) == 1) UP(xx,yy) = 1
                        if (UP(nx,yy) == 1) UP(xx,yy) = 1
                    else if (xx == nx) then
                        if (UP(nx-1,yy) == 1) UP(xx,yy) = 1
                        if (UP(1,yy) == 1) UP(xx,yy) = 1
                    else
                        if (UP(xx-1,yy) == 1) UP(xx,yy) = 1
                        if (UP(xx+1,yy) == 1) UP(xx,yy) = 1
                    end if 
                    if (UP(xx,yy) == 1) cnt = cnt + 1
                end if 
            end do
            do xx = nx, 1, -1
                if (B(xx,yy) /= 0) cycle
                if (yy == 1) UP(xx,yy) = 1
                if (yy /= 1) then
                    if (UP(xx,yy) == 1) cycle 
                    if (UP(xx,yy-1) == 1) UP(xx,yy) = 1
                    if (xx == 1) then
                        if (UP(2,yy) == 1) UP(xx,yy) = 1
                        if (UP(nx,yy) == 1) UP(xx,yy) = 1
                    else if (xx == nx) then
                        if (UP(nx-1,yy) == 1) UP(xx,yy) = 1
                        if (UP(1,yy) == 1) UP(xx,yy) = 1
                    else
                        if (UP(xx-1,yy) == 1) UP(xx,yy) = 1
                        if (UP(xx+1,yy) == 1) UP(xx,yy) = 1
                    end if 
                    if (UP(xx,yy) == 1) cnt = cnt + 1 
                end if 
            end do 
        end do 
        ! check bottom to top connection
        do yy = ny, 1, -1
            do xx = 1, nx
                if (B(xx,yy) /= 0) cycle
                if (yy == ny) cycle
                if (yy /= ny) then
                    if (UP(xx,yy) == 1) cycle 
                    if (UP(xx,yy+1) == 1) UP(xx,yy) = 1
                    if (xx == 1) then
                        if (UP(2,yy) == 1) UP(xx,yy) = 1
                        if (UP(nx,yy) == 1) UP(xx,yy) = 1
                    else if (xx == nx) then
                        if (UP(nx-1,yy) == 1) UP(xx,yy) = 1
                        if (UP(1,yy) == 1) UP(xx,yy) = 1
                    else
                        if (UP(xx-1,yy) == 1) UP(xx,yy) = 1
                        if (UP(xx+1,yy) == 1) UP(xx,yy) = 1
                    end if 
                    if (UP(xx,yy) == 1) cnt = cnt + 1 
                end if 
            end do 
            do xx = nx, 1, -1
                if (B(xx,yy) /= 0) cycle
                if (yy == ny) cycle
                if (yy /= ny) then
                    if (UP(xx,yy) == 1) cycle 
                    if (UP(xx,yy+1) == 1) UP(xx,yy) = 1
                    if (xx == 1) then
                        if (UP(2,yy) == 1) UP(xx,yy) = 1
                        if (UP(nx,yy) == 1) UP(xx,yy) = 1
                    else if (xx == nx) then
                        if (UP(nx-1,yy) == 1) UP(xx,yy) = 1
                        if (UP(1,yy) == 1) UP(xx,yy) = 1
                    else
                        if (UP(xx-1,yy) == 1) UP(xx,yy) = 1
                        if (UP(xx+1,yy) == 1) UP(xx,yy) = 1
                    end if 
                    if (UP(xx,yy) == 1) cnt = cnt + 1 
                end if 
            end do 
        end do
        ! check bottom to top connection
        do yy = ny, 1, -1
            do xx = 1, nx
                if (B(xx,yy) /= 0) cycle
                if (yy == ny) DOWN(xx,yy) = 1
                if (yy /= ny) then
                    if (DOWN(xx,yy) == 1) cycle 
                    if (DOWN(xx,yy+1) == 1) DOWN(xx,yy) = 1
                    if (xx == 1) then
                        if (DOWN(2,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(nx,yy) == 1) DOWN(xx,yy) = 1
                    else if (xx == nx) then
                        if (DOWN(nx-1,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(1,yy) == 1) DOWN(xx,yy) = 1
                    else
                        if (DOWN(xx-1,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(xx+1,yy) == 1) DOWN(xx,yy) = 1
                    end if 
                    if (DOWN(xx,yy) == 1) cnt = cnt + 1 
                end if 
            end do 
            do xx = nx, 1, -1
                if (B(xx,yy) /= 0) cycle
                if (yy == ny) DOWN(xx,yy) = 1
                if (yy /= ny) then
                    if (DOWN(xx,yy) == 1) cycle 
                    if (DOWN(xx,yy+1) == 1) DOWN(xx,yy) = 1
                    if (xx == 1) then
                        if (DOWN(2,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(nx,yy) == 1) DOWN(xx,yy) = 1
                    else if (xx == nx) then
                        if (DOWN(nx-1,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(1,yy) == 1) DOWN(xx,yy) = 1
                    else
                        if (DOWN(xx-1,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(xx+1,yy) == 1) DOWN(xx,yy) = 1
                    end if 
                    if (DOWN(xx,yy) == 1) cnt = cnt + 1 
                end if 
            end do 
        end do
        ! check top to bottom 
        do yy = 1, ny
            do xx = 1, nx
                if (B(xx,yy) /= 0) cycle ! obstacles
                if (yy == 1) cycle
                if (yy /= 1) then
                    if (DOWN(xx,yy) == 1) cycle 
                    if (DOWN(xx,yy-1) == 1) DOWN(xx,yy) = 1
                    if (xx == 1) then
                        if (DOWN(2,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(nx,yy) == 1) DOWN(xx,yy) = 1
                    else if (xx == nx) then
                        if (DOWN(nx-1,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(1,yy) == 1) DOWN(xx,yy) = 1
                    else
                        if (DOWN(xx-1,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(xx+1,yy) == 1) DOWN(xx,yy) = 1
                    end if 
                    if (DOWN(xx,yy) == 1) cnt = cnt + 1
                end if 
            end do
            do xx = nx, 1, -1
                if (B(xx,yy) /= 0) cycle
                if (yy == 1) cycle
                if (yy /= 1) then
                    if (DOWN(xx,yy) == 1) cycle 
                    if (DOWN(xx,yy-1) == 1) DOWN(xx,yy) = 1
                    if (xx == 1) then
                        if (DOWN(2,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(nx,yy) == 1) DOWN(xx,yy) = 1
                    else if (xx == nx) then
                        if (DOWN(nx-1,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(1,yy) == 1) DOWN(xx,yy) = 1
                    else
                        if (DOWN(xx-1,yy) == 1) DOWN(xx,yy) = 1
                        if (DOWN(xx+1,yy) == 1) DOWN(xx,yy) = 1
                    end if 
                    if (DOWN(xx,yy) == 1) cnt = cnt + 1 
                end if 
            end do 
        end do 
    enddo
    
    ! print *, cnt2

! in case flow direction is x
else if (flow_dir_x) then 
	! check left to right connection
	do xx = 1, nx
		do yy = 1, ny
			if (B(xx,yy) /= 0) cycle
			if (xx == 1) UP(xx,yy) = 1
			if (xx /= 1) then
				if (UP(xx-1,yy) == 1) UP(xx,yy) = 1
				if (yy == 1) then
					if (UP(xx,2) == 1) UP(xx,yy) = 1
					if (UP(xx,ny) == 1) UP(xx,yy) = 1
				else if (yy == ny) then
					if (UP(xx,ny-1) == 1) UP(xx,yy) = 1
					if (UP(xx,1) == 1) UP(xx,yy) = 1
				else
					if (UP(xx,yy-1) == 1) UP(xx,yy) = 1
					if (UP(xx,yy+1) == 1) UP(xx,yy) = 1
				end if 
			end if 
		end do
		do yy = ny, 1, -1
			if (B(xx,yy) /= 0) cycle
			if (xx == 1) UP(xx,yy) = 1
			if (xx /= 1) then
				if (UP(xx-1,yy) == 1) UP(xx,yy) = 1
				if (yy == 1) then
					if (UP(xx,2) == 1) UP(xx,yy) = 1
					if (UP(xx,ny) == 1) UP(xx,yy) = 1
				else if (yy == ny) then
					if (UP(xx,ny-1) == 1) UP(xx,yy) = 1
					if (UP(xx,1) == 1) UP(xx,yy) = 1
				else
					if (UP(xx,yy-1) == 1) UP(xx,yy) = 1
					if (UP(xx,yy+1) == 1) UP(xx,yy) = 1
				end if 
			end if 
		end do 
	end do 
	! check right to left connection
	do xx = nx, 1, -1
		do yy = 1, ny
			if (B(xx,yy) /= 0) cycle
			if (xx == nx) DOWN(xx,yy) = 1
			if (xx /= nx) then
				if (DOWN(xx+1,yy) == 1) DOWN(xx,yy) = 1
				if (yy == 1) then
					if (DOWN(xx,2) == 1) DOWN(xx,yy) = 1
					if (DOWN(xx,ny) == 1) DOWN(xx,yy) = 1
				else if (yy == ny) then
					if (DOWN(xx,ny-1) == 1) DOWN(xx,yy) = 1
					if (DOWN(xx,1) == 1) DOWN(xx,yy) = 1
				else
					if (DOWN(xx,yy-1) == 1) DOWN(xx,yy) = 1
					if (DOWN(xx,yy+1) == 1) DOWN(xx,yy) = 1
				end if 
			end if 
		end do 
		do yy = ny, 1, -1
			if (B(xx,yy) /= 0) cycle
			if (xx == nx) DOWN(xx,yy) = 1
			if (xx /= nx) then
				if (DOWN(xx+1,yy) == 1) DOWN(xx,yy) = 1
				if (yy == 1) then
					if (DOWN(xx,2) == 1) DOWN(xx,yy) = 1
					if (DOWN(xx,ny) == 1) DOWN(xx,yy) = 1
				else if (yy == ny) then
					if (DOWN(xx,ny-1) == 1) DOWN(xx,yy) = 1
					if (DOWN(xx,1) == 1) DOWN(xx,yy) = 1
				else
					if (DOWN(xx,yy-1) == 1) DOWN(xx,yy) = 1
					if (DOWN(xx,yy+1) == 1) DOWN(xx,yy) = 1
				end if 
			end if 
		end do 
	end do
end if 

! making a connected binary 
CONN = 100
where ((UP == 1).or.(DOWN == 1))
	CONN = 0   ! connected water = 0; obstacles and non-connected waters = 100
end where 

!!!! check the connected bindaries
if ((rec_binary)) then 
	open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-UP-'  & 
		//trim(adjustl(numtemp))//'.txt', status = 'replace')
	do yy = 1,ny
		write(100,*) (UP(xx,yy), xx = 1,nx)
	end do 
	close(100)

	open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-DOWN-'  & 
		//trim(adjustl(numtemp))//'.txt', status = 'replace')
	do yy = 1,ny
		write(100,*) (DOWN(xx,yy), xx = 1,nx)
	end do 
	close(100)

	open(unit = 100, file = trim(adjustl(today))//''//trim(adjustl(cat_tmp))//'/test2d(F)-CONN-'  & 
		//trim(adjustl(numtemp))//'.txt', status = 'replace')
	do yy = 1,ny
		write(100,*) (CONN(xx,yy), xx = 1,nx)
	end do 
	close(100)
end if 

! add outer boundaries
B2 = 0
B2(:,2:ny+1) = CONN(:,:)

! type of each cell
Tp = 0
Tp2 = 0
Tp = B2
Tp(:,1) = 300  ! upper boundary
Tp(:,ny+2) = 200  ! lower boudary

! modifying lower boundary to help calculation
yy = ny+1
do xx = 1, nx
	if (B2(xx,yy) == 100) then
		B2(xx, ny+2) = 100
		Tp(xx, ny+2) = 100
	end if 
end do 

! classify the type of obstacle depending on surrounding water
do xx = 1, nx
	do yy = 1, ny+2
		if (B2(xx,yy) /= 0) then              
			cnt = B2(xx,yy)  ! i.e., 100
			xp = xx + 1
			xg = xx - 1
			yp = yy + 1
			yg = yy - 1
			if (xp == nx+1) xp = 1
			if (xg == 0) xg = nx
			if (yp == ny+3) yp = yy
			if (yg == 0) yg = yy
			
			if (B2(xp,yy) == 0) cnt = cnt + 2   ! numbering clockwise
			if (B2(xg,yy) == 0) cnt = cnt + 8   ! 
			if (B2(xx,yp) == 0) cnt = cnt + 4   ! 
			if (B2(xx,yg) == 0) cnt = cnt + 1   ! 
			Tp(xx,yy) = cnt
		end if 
	end do
end do

! classify the type of water depedning on surrounding obstacles
do xx = 1, nx
	do yy = 1, ny+2
		if (B2(xx,yy) == 0) then 
			cnt = 0
			xp = xx + 1
			xg = xx - 1
			yp = yy + 1
			yg = yy - 1
			if (xp == nx+1) xp = 1
			if (xg == 0) xg = nx
			if (yp == ny+3) yp = yy
			if (yg == 0) yg = yy
			
			if (B2(xp,yy) == 100) cnt = cnt + 2  ! numbering clockwise
			if (B2(xg,yy) == 100) cnt = cnt + 8  ! 
			if (B2(xx,yp) == 100) cnt = cnt + 4  ! 
			if (B2(xx,yg) == 100) cnt = cnt + 1  ! 
			Tp2(xx,yy) = cnt
		end if 
	end do
end do 

! choose randomly a location for constant flow boundary 
if (.not.calc_perm) then 
    if (random_choice) then 
        choice_done = .false.
        call random_seed
        do while (.not.choice_done)
            call random_number(rdm)
            xx = 3 + int(rdm(1)*(nx-3))
            yy = 3 + int(rdm(2)*(ny-3))
            xp = xx + 1
            xg = xx - 1
            yp = yy + 1
            yg = yy - 1
            if (.not. simple_test) then 
                if ((B2(xx,yy) == 0).and.(yy > int(3*ny/10))) then 
                    dirc = 1 + int(rdm(3)*4)
                    if ((dirc == 1).and.(B2(xp,yy) == 0)) Tp(xx,yy) = 500; choice_done = .true.  ! rightward flow
                    if ((dirc == 2).and.(B2(xg,yy) == 0)) Tp(xx,yy) = 600; choice_done = .true.  ! leftward flow
                    if ((dirc == 3).and.(B2(xx,yp) == 0)) Tp(xx,yy) = 700; choice_done = .true.  ! downward flow
                    if ((dirc == 4).and.(B2(xx,yg) == 0)) Tp(xx,yy) = 800; choice_done = .true.  ! upward flow
                end if 
            else if (simple_test) then 
                if ((B2(xx,yy) == 0)) then 
                    dirc = 1 + int(rdm(4)*4)
                    if ((dirc == 1).and.(B2(xp,yy) == 0)) Tp(xx,yy) = 500; choice_done = .true.  ! rightward flow
                    if ((dirc == 2).and.(B2(xg,yy) == 0)) Tp(xx,yy) = 600; choice_done = .true.  ! leftward flow
                    if ((dirc == 3).and.(B2(xx,yp) == 0)) Tp(xx,yy) = 700; choice_done = .true.  ! downward flow
                    if ((dirc == 4).and.(B2(xx,yg) == 0)) Tp(xx,yy) = 800; choice_done = .true.  ! upward flow
                end if 
            end if 
        end do 
    end if 
    if (.not. random_choice) then 
        do i = 1, n_ind
            do pp = 1,2
                if ((Vb(pp,i) /= 0).or.(Ub(pp,i) /= 0)) then 
                    xx = flow_loc(pp,i)%X
                    yy = flow_loc(pp,i)%Y
                    if (flw_rsltn == 1) then 
                        Tp(xx,yy+1) = 900
                    elseif (flw_rsltn/=1) then 
                        Tp((xx-1)*flw_rsltn+1,(yy-1)*flw_rsltn+1+1) = 900
                    endif 
                    ! print*, xx,yy,matrix(yy,xx)%class
                    if (matrix(yy,xx)%class /= i) then 
                      print *, '+++ error in flow boundary +++'
                      write(File_log,*) '+++ error in flow boundary +++ (x,y)=',xx,yy,'class = ',matrix(yy,xx)%class 
                      flg_stop=.true.
                      ! stop
                    endif
                end if 
            end do
        end do
    end if 
endif 
				
! check the type-classified map
if (rec_binary) then 
	open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-Tp-'  & 
		//trim(adjustl(numtemp))//'.txt', status = 'replace')
	do yy = 1,ny+2
		write(100,*) (Tp(xx,yy), xx = 1,nx)
	end do 
	close(100)

	open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-Tp2-'  & 
		//trim(adjustl(numtemp))//'.txt', status = 'replace')
	do yy = 1,ny+2
		write(100,*) (Tp2(xx,yy), xx = 1,nx)
	end do 
	close(100)
end if 

! determine the cell number where P is calculated
np = 0
do yy = 1, ny+2
	do xx = 1, nx
		if (Tp(xx,yy) == 0) np = np + 1
	end do
end do

! assign numbers on cell grids and associate the numbers with locations 
if (allocated(CP)) deallocate(CP)
allocate(CP(np))
Cloc = 0
Cnum = 0
CP = 0
CPi = 0
cnt = 0 
cnt2 = 0
do xx = 1, nx
	do yy = 1, ny+2
		Cnum(xx,yy) = cnt + 1
		Cloc(cnt+1,1) = xx
		Cloc(cnt+1,2) = yy
		if (Tp(xx,yy) == 0) then
			CP(cnt2+1) = cnt + 1
			CPi(cnt+1) = cnt2 + 1
			cnt2 = cnt2 + 1
		end if 
		cnt = cnt + 1
	end do
end do 
if (cnt /= ntc) write(File_log, *) Time, 'error cnt', cnt, ntc
if (cnt2 /= np) write(File_log, *) Time, 'error cnt2', cnt2, np

! information on neighbors of P-cells
if (allocated(CDN)) deallocate(CDN)
if (allocated(CDNt)) deallocate(CDNt)
allocate(CDN(np,4), CDNt(np,4))
CDN = 0
CDNt = 0
do pp = 1, np
	cnt = CP(pp)
	xx = Cloc(cnt,1)
	yy = Cloc(cnt,2)
	
	xp = xx + 1
	xg = xx - 1
	yp = yy + 1
	yg = yy - 1
	if (xp == nx+1) xp = 1
	if (xg == 0) xg = nx
	if (yp == ny+3) yp = yy
	if (yg == 0) yg = yy
	
	CDN(pp,1) = Cnum(xp,yy)
	CDN(pp,2) = Cnum(xg,yy)
	CDN(pp,3) = Cnum(xx,yp)
	CDN(pp,4) = Cnum(xx,yg)
	
	CDNt(pp,1) = Tp(xp,yy)
	CDNt(pp,2) = Tp(xg,yy)
	CDNt(pp,3) = Tp(xx,yp)
	CDNt(pp,4) = Tp(xx,yg)
end do

! couting the non-zero component in matrix for P calculation 
cnt = 0
do pp = 1, np
	cnt = cnt + 1
	do ppp = 1,4
		if (CDNt(pp,ppp)==0) cnt = cnt + 1
	end do
end do 
nnz = cnt ! number of non-zero component
n = np ! number of row 
allocate (ai(nnz),ap(n+1),ax(nnz),bx(n),kai(n))

! set all initial values at zero
Pc = 0.0d0
Uc = 0.0d0
Vc = 0.0d0
Um = 0.0d0
Vm = 0.0d0
Dc = 0.0d0
Fc = 0.0d0
Gc = 0.0d0
Rc = 0.0d0 

! start calculation 
step = 0
do while (step < nstep)
	! calculate u,v from previous iteration  
	do xx = 1, nx
		do yy = 1, ny + 2
			xp = xx + 1
			xg = xx - 1
			yp = yy + 1
			yg = yy - 1
			if (xp == nx+1) xp = 1
			if (xg == 0) xg = nx
			if (yp == ny+3) yp = yy
			if (yg == 0) yg = yy
			
			if (B2(xx,yy)==0) then
				Um(xx,yy) = Um(xx,yy) + delt*(-Fc(xx,yy) - (Pc(xx,yy) - Pc(xg,yy))/delx)
				Vm(xx,yy) = Vm(xx,yy) + delt*(-Gc(xx,yy) - (Pc(xx,yy) - Pc(xx,yg))/dely)
			end if 
		end do
	end do 
	! u, v boundary
    if (show_display_details .and. step == 0 .and. (.not.calc_perm)) then 
        print '(A8,A)', '','Boundary conditions'
        print '(A8,A)', '', 'p     =   particle               w     =   water'
        print '(A8,A)', '', 'U     =   upper boundary         L     =   lower boundary'
        print '(A8,A)', '', 'Arrow =   imposed water flow     ?     =   ?????'
        print *
    endif 
	do xx = 1, nx
		do yy = 1, ny+2
            xp = xx + 1
            xg = xx - 1
            yp = yy + 1
            yg = yy - 1
            if (xp == nx+1) xp = 1
            if (xg == 0) xg = nx
            if (yp == ny+3) yp = yy
            if (yg == 0) yg = yy
			
			if (Tp(xx,yy) == 900) then ! const flow in any direction
                ! print *, xx, yy
                ! print *, flow_loc
                ! print *, Ub
                ! print *, Vb
				ppp = 0
                if (flw_rsltn==1) then 
                    xxx = xx
                    yyy = yy - 1
                elseif (flw_rsltn/=1) then 
                    xxx = (xx -1)/flw_rsltn + 1
                    yyy = (yy- 1 - 1)/flw_rsltn + 1
                endif 
				do pp = 1,2
					if ((Vb(pp,matrix(yyy,xxx)%class)==0).and.(Ub(pp,matrix(yyy,xxx)%class)==0)) cycle
					if (flow_loc(pp, matrix(yyy,xxx)%class )%X == xxx   &
                        .and. flow_loc(pp, matrix(yyy,xxx)%class )%Y == yyy) then
                        ppp = pp
                    endif 
                    ! print *, '----', pp, Vb(pp,matrix(yyy,xxx)%class), Ub(pp,matrix(yyy,xxx)%class), ppp
                    ! print *, 'xxxx', flow_loc(pp, matrix(yyy,xxx)%class ) 
				end do 
				if (ppp == 0) then 
                    print *, time,'error to set boundary 900',xx,yy
                    write(File_log, *) time, 'error to set boundary 900',xx,yy
                    flg_stop = .true.
                endif
                if (B2(xp,yy) == 0) Um(xp,yy) = Ub(ppp,matrix(yyy,xxx)%class)*VelC
                if (B2(xg,yy) == 0) Um(xg+1,yy) = Ub(ppp,matrix(yyy,xxx)%class)*VelC
                if (B2(xx,yp) == 0) Vm(xx,yp) = Vb(ppp,matrix(yyy,xxx)%class)*VelC
                if (B2(xx,yg) == 0) Vm(xx,yg+1) = Vb(ppp,matrix(yyy,xxx)%class)*VelC
                
                if (Ub(ppp,matrix(yyy,xxx)%class) > 0d0 .and. Vb(ppp,matrix(yyy,xxx)%class) > 0d0) then 
                    flow_chr = '\v'
                elseif (Ub(ppp,matrix(yyy,xxx)%class) < 0d0 .and. Vb(ppp,matrix(yyy,xxx)%class) < 0d0) then 
                    flow_chr = '^\'
                elseif (Ub(ppp,matrix(yyy,xxx)%class) < 0d0 .and. Vb(ppp,matrix(yyy,xxx)%class) > 0d0) then 
                    flow_chr = 'v/'
                elseif (Ub(ppp,matrix(yyy,xxx)%class) > 0d0 .and. Vb(ppp,matrix(yyy,xxx)%class) < 0d0) then 
                    flow_chr = '/^'
                elseif (Ub(ppp,matrix(yyy,xxx)%class) > 0d0 .and. Vb(ppp,matrix(yyy,xxx)%class) == 0d0) then 
                    flow_chr = '->'
                elseif (Ub(ppp,matrix(yyy,xxx)%class) < 0d0 .and. Vb(ppp,matrix(yyy,xxx)%class) == 0d0) then 
                    flow_chr = '<-'
                elseif (Ub(ppp,matrix(yyy,xxx)%class) == 0d0 .and. Vb(ppp,matrix(yyy,xxx)%class) > 0d0) then 
                    flow_chr = 'v'
                elseif (Ub(ppp,matrix(yyy,xxx)%class) == 0d0 .and. Vb(ppp,matrix(yyy,xxx)%class) < 0d0) then 
                    flow_chr = '^'
                else
                    flow_chr = 'F?'
                endif 
                
                if (ppp == 1) loc_chr = 'head'
                if (ppp == 2) loc_chr = 'tail'
                
                if (show_display_details .and. step == 0 .and. (.not.calc_perm)) then 
                    print '(A8,A,I0,A,A)', '','Org # = ',matrix(yyy,xxx)%class,' around ',trim(adjustl(loc_chr))
                    ! print *, ppp, matrix(yyy,xxx)%class
                    print *
                    if (trim(adjustl(loc_chr)) == 'tail' .and. ppp /=2) pause
                    if (trim(adjustl(loc_chr)) == 'head' .and. ppp /=1) pause
                    do y_tmp = 1,21
                        yyy = yy + (y_tmp - 11) 
                        if (yyy > ny+2) yyy = ny+2
                        if (yyy < 1) yyy = 1
                        if (y_tmp==1 ) then 
                            do x_tmp = 1,21
                                xxx = xx + (x_tmp - 11)
                                if (xxx > nx) xxx = xxx - nx
                                if (xxx < 1) xxx = xxx + nx
                                adv_chr = 'no'
                                if (x_tmp==21) adv_chr = 'yes'
                                if (x_tmp==1) write (*,fmt='(a11)', advance=trim(adjustl(adv_chr))) '' 
                                write (*,fmt='(I3)', advance=trim(adjustl(adv_chr))) xxx 
                            enddo
                        endif 
                        do x_tmp = 1,21
                            xxx = xx + (x_tmp - 11)
                            if (xxx > nx) xxx = xxx - nx
                            if (xxx < 1) xxx = xxx + nx
                            adv_chr = 'no'
                            if (x_tmp==21) adv_chr = 'yes'
                            if (x_tmp==1) write (*,fmt='(i11)', advance=trim(adjustl(adv_chr))) yyy
                            if (Tp(xxx,yyy) == 900) then 
                                ! write (*,fmt='(a3)', advance=trim(adjustl(adv_chr))) ' F '
                                write (*,fmt='(a3)', advance=trim(adjustl(adv_chr))) trim(adjustl(flow_chr))
                            elseif (Tp(xxx,yyy) == 0) then 
                                write (*,fmt='(a3)', advance=trim(adjustl(adv_chr))) ' w '
                            elseif (Tp(xxx,yyy) == 200) then 
                                write (*,fmt='(a3)', advance=trim(adjustl(adv_chr))) ' L '
                            elseif (Tp(xxx,yyy) == 300) then 
                                write (*,fmt='(a3)', advance=trim(adjustl(adv_chr))) ' U '
                            elseif (Tp(xxx,yyy) >= 100) then 
                                write (*,fmt='(a3)', advance=trim(adjustl(adv_chr))) ' p '
                            else 
                                write (*,fmt='(a3)', advance=trim(adjustl(adv_chr))) ' ? '
                            endif 
                        enddo 
                    enddo 
                    print *
                endif 
                loc_chr = ' '
			end if 
            
			if (B2(xx,yy)==100) then 
				if (B2(xp,yy) == 0) Um(xp,yy) = 0.0d0
				if (B2(xg,yy) == 0) Um(xg+1,yy) = 0.0d0
				if (B2(xx,yp) == 0) Vm(xx,yp) = 0.0d0
				if (B2(xx,yg) == 0) Vm(xx,yg+1) = 0.0d0
			end if 
				
			if ((.not.const_bot).and.(Tp(xx,yy) == 200)) then ! lower boundary; no pressure dradient in y
				Vm(xx,yy) = 0.0d0
				Vm(xx,yy+1) = 0.0d0
			end if 
			if ((.not.const_top).and.(Tp(xx,yy) == 300)) then ! upper boundary; no pressure dradient in y 
				Vm(xx,yy) = 0.0d0
				Vm(xx,yy+1) = 0.0d0
			end if 
			
			if (xx == nx) Um(xx+1,yy) = Um(1,yy) ! continuous boundary
		end do
	end do 
	
	do xx = 1, nx
		do yy = 1, ny+2
			xp = xx + 1
			xg = xx - 1
			yp = yy + 1
			yg = yy - 1
			if (xp == nx+1) xp = 1
			if (xg == 0) xg = nx
			if (yp == ny+3) yp = yy
			if (yg == 0) yg = yy
			
			Dc(xx,yy) = (Um(xx+1,yy)-Um(xx,yy))/delx + (Vm(xx,yy+1)-Vm(xx,yy))/dely
			
			Fc(xx,yy) = ((Um(xx+1,yy)+Um(xx,yy))**2 - (Um(xg,yy)+Um(xx,yy))**2)/(4*delx) &
				+((Um(xx,yp)+Um(xx,yy))*(Vm(xx,yy+1)+Vm(xg,yy+1)) - (Um(xx,yg)+Um(xx,yy))*(Vm(xx,yy)+Vm(xg,yy)))/(4*dely) &
				-new*((Um(xx+1,yy)+Um(xg,yy)-2*Um(xx,yy))/(delx**2) + (Um(xx,yp)+Um(xx,yg)-2*Um(xx,yy))/(dely**2))
				
			Gc(xx,yy) = ((Vm(xx,yy)+Vm(xx,yy+1))**2 - (Vm(xx,yy)+Vm(xx,yg))**2)/(4*dely) &
				+((Um(xx+1,yy)+Um(xx+1,yg))*(Vm(xp,yy)+Vm(xx,yy)) - (Um(xx,yy)+Um(xx,yg))*(Vm(xx,yy)+Vm(xg,yy)))/(4*delx) &
				-new*((Vm(xp,yy)+Vm(xg,yy)-2*Vm(xx,yy))/(delx**2) + (Vm(xx,yy+1)+Vm(xx,yg)-2*Vm(xx,yy))/(dely**2))
			
			if (B2(xx,yy) == 0) then 
				if (Tp2(xx,yy) == 2) then  ! right wall
					Fc(xp,yy) = -new*(2*Um(xx,yy)/(delx**2))
					
					Gc(xx,yy) = ((Vm(xx,yy)+Vm(xx,yy+1))**2 - (Vm(xx,yy)+Vm(xx,yg))**2)/(4*dely) &
						+(0.0d0 - (Um(xx,yy)+Um(xx,yg))*(Vm(xx,yy)+Vm(xg,yy)))/(4*delx) &
						-new*((Vm(xg,yy)-3*Vm(xx,yy))/(delx**2) + (Vm(xx,yy+1)+Vm(xx,yg)-2*Vm(xx,yy))/(dely**2))
					
					Gc(xx,yp) = ((Vm(xx,yp)+Vm(xx,yp+1))**2 - (Vm(xx,yp)+Vm(xx,yy))**2)/(4*dely) &
						+(0.0d0 - (Um(xx,yp)+Um(xx,yy))*(Vm(xx,yp)+Vm(xg,yp)))/(4*delx) &
						-new*((Vm(xg,yp)-3*Vm(xx,yp))/(delx**2) + (Vm(xx,yp+1)+Vm(xx,yy)-2*Vm(xx,yp))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 8) then  ! left wall 
					Fc(xx,yy) = -new*((2*Um(xx+1,yy))/(delx**2))
			
					Gc(xx,yy) = ((Vm(xx,yy)+Vm(xx,yy+1))**2 - (Vm(xx,yy)+Vm(xx,yg))**2)/(4*dely) &
						+((Um(xx+1,yy)+Um(xx+1,yg))*(Vm(xp,yy)+Vm(xx,yy)) - 0.0d0)/(4*delx) &
						-new*((Vm(xp,yy)-3*Vm(xx,yy))/(delx**2) + (Vm(xx,yy+1)+Vm(xx,yg)-2*Vm(xx,yy))/(dely**2))
						
					Gc(xx,yp) = ((Vm(xx,yp)+Vm(xx,yp+1))**2 - (Vm(xx,yp)+Vm(xx,yy))**2)/(4*dely) &
						+((Um(xx+1,yp)+Um(xx+1,yy))*(Vm(xp,yp)+Vm(xx,yp)) - 0.0d0)/(4*delx) &
						-new*((Vm(xp,yp)-3*Vm(xx,yp))/(delx**2) + (Vm(xx,yp+1)+Vm(xx,yy)-2*Vm(xx,yp))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 4) then ! down wall 
					Fc(xx,yy) = ((Um(xx+1,yy)+Um(xx,yy))**2 - (Um(xg,yy)+Um(xx,yy))**2)/(4*delx) &
						+(0.0d0 - (Um(xx,yg)+Um(xx,yy))*(Vm(xx,yy)+Vm(xg,yy)))/(4*dely) &
						-new*((Um(xx+1,yy)+Um(xg,yy)-2*Um(xx,yy))/(delx**2) + (Um(xx,yg)-3*Um(xx,yy))/(dely**2))
						
					Fc(xp,yy) = ((Um(xp+1,yy)+Um(xp,yy))**2 - (Um(xx,yy)+Um(xp,yy))**2)/(4*delx) &
						+(0.0d0 - (Um(xp,yg)+Um(xp,yy))*(Vm(xp,yy)+Vm(xx,yy)))/(4*dely) &
						-new*((Um(xp+1,yy)+Um(xx,yy)-2*Um(xp,yy))/(delx**2) + (Um(xp,yg)-3*Um(xp,yy))/(dely**2))
						
					Gc(xx,yp) = -new*((2*Vm(xx,yy))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 1) then ! up wall 
					Fc(xx,yy) = ((Um(xx+1,yy)+Um(xx,yy))**2 - (Um(xg,yy)+Um(xx,yy))**2)/(4*delx) &
						+((Um(xx,yp)+Um(xx,yy))*(Vm(xx,yy+1)+Vm(xg,yy+1)) - 0.0d0)/(4*dely) &
						-new*((Um(xx+1,yy)+Um(xg,yy)-2*Um(xx,yy))/(delx**2) + (Um(xx,yp)-3*Um(xx,yy))/(dely**2))
						
					Fc(xp,yy) = ((Um(xp+1,yy)+Um(xp,yy))**2 - (Um(xx,yy)+Um(xp,yy))**2)/(4*delx) &
						+((Um(xp,yp)+Um(xp,yy))*(Vm(xp,yy+1)+Vm(xx,yy+1)) - 0.0d0)/(4*dely) &
						-new*((Um(xp+1,yy)+Um(xx,yy)-2*Um(xp,yy))/(delx**2) + (Um(xp,yp)-3*Um(xp,yy))/(dely**2))
						
					Gc(xx,yy) = -new*((2*Vm(xx,yy+1))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 3) then ! up and right walls 
					Fc(xx,yy) = ((Um(xx+1,yy)+Um(xx,yy))**2 - (Um(xg,yy)+Um(xx,yy))**2)/(4*delx) &
						+((Um(xx,yp)+Um(xx,yy))*(Vm(xx,yy+1)+Vm(xg,yy+1)) - 0.0d0)/(4*dely) &
						-new*((Um(xx+1,yy)+Um(xg,yy)-2*Um(xx,yy))/(delx**2) + (Um(xx,yp)-3*Um(xx,yy))/(dely**2))
						
					Fc(xp,yy) = -new*((2*Um(xx,yy))/(delx**2))
						
					Gc(xx,yy) = -new*((2*Vm(xx,yy+1))/(dely**2))
						
					Gc(xx,yp) = ((Vm(xx,yp)+Vm(xx,yp+1))**2 - (Vm(xx,yp)+Vm(xx,yy))**2)/(4*dely) &
						+(0.0d0 - (Um(xx,yp)+Um(xx,yy))*(Vm(xx,yp)+Vm(xg,yp)))/(4*delx) &
						-new*((Vm(xg,yp)-3*Vm(xx,yp))/(delx**2) + (Vm(xx,yp+1)+Vm(xx,yy)-2*Vm(xx,yp))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 5) then ! up and down walls 
					Fc(xx,yy) = ((Um(xx+1,yy)+Um(xx,yy))**2 - (Um(xg,yy)+Um(xx,yy))**2)/(4*delx) &
						-new*((Um(xx+1,yy)+Um(xg,yy)-2*Um(xx,yy))/(delx**2) + (-4*Um(xx,yy))/(dely**2))
						
					Fc(xp,yy) = ((Um(xp+1,yy)+Um(xp,yy))**2 - (Um(xx,yy)+Um(xp,yy))**2)/(4*delx) &
						-new*((Um(xp+1,yy)+Um(xx,yy)-2*Um(xp,yy))/(delx**2) + (-4*Um(xp,yy))/(dely**2))
						
					Gc(xx,yy) = 0.0d0
						
					Gc(xx,yp) = 0.0d0
				end if 
				
				if (Tp2(xx,yy) == 9) then ! up and left walls 
					Fc(xx,yy) = -new*((2*Um(xx+1,yy))/(delx**2))
						
					Fc(xp,yy) = ((Um(xp+1,yy)+Um(xp,yy))**2 - (Um(xx,yy)+Um(xp,yy))**2)/(4*delx) &
						+((Um(xp,yp)+Um(xp,yy))*(Vm(xp,yy+1)+Vm(xx,yy+1)) - 0.0d0)/(4*dely) &
						-new*((Um(xp+1,yy)+Um(xx,yy)-2*Um(xp,yy))/(delx**2) + (Um(xp,yp)-3*Um(xp,yy))/(dely**2))
						
					Gc(xx,yy) = -new*((2*Vm(xx,yy+1))/(dely**2))
						
					Gc(xx,yp) = ((Vm(xx,yp)+Vm(xx,yp+1))**2 - (Vm(xx,yp)+Vm(xx,yy))**2)/(4*dely) &
						+((Um(xx+1,yp)+Um(xx+1,yy))*(Vm(xp,yp)+Vm(xx,yp)) - 0.0d0)/(4*delx) &
						-new*((Vm(xp,yp)-3*Vm(xx,yp))/(delx**2) + (Vm(xx,yp+1)+Vm(xx,yy)-2*Vm(xx,yp))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 6) then ! right and down walls 
					Fc(xx,yy) = ((Um(xx+1,yy)+Um(xx,yy))**2 - (Um(xg,yy)+Um(xx,yy))**2)/(4*delx) &
						+(0.0d0 - (Um(xx,yg)+Um(xx,yy))*(Vm(xx,yy)+Vm(xg,yy)))/(4*dely) &
						-new*((Um(xx+1,yy)+Um(xg,yy)-2*Um(xx,yy))/(delx**2) + (Um(xx,yg)-3*Um(xx,yy))/(dely**2))
						
					Fc(xp,yy) = -new*((2*Um(xx,yy))/(delx**2))
						
					Gc(xx,yy) = ((Vm(xx,yy)+Vm(xx,yy+1))**2 - (Vm(xx,yy)+Vm(xx,yg))**2)/(4*dely) &
						+(0.0d0 - (Um(xx,yy)+Um(xx,yg))*(Vm(xx,yy)+Vm(xg,yy)))/(4*delx) &
						-new*((Vm(xg,yy)-3*Vm(xx,yy))/(delx**2) + (Vm(xx,yy+1)+Vm(xx,yg)-2*Vm(xx,yy))/(dely**2))
						
					Gc(xx,yp) = -new*((2*Vm(xx,yy))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 10) then ! right and left walls 
					Fc(xx,yy) = 0.0d0
						
					Fc(xp,yy) = 0.0d0
						
					Gc(xx,yy) = ((Vm(xx,yy)+Vm(xx,yy+1))**2 - (Vm(xx,yy)+Vm(xx,yg))**2)/(4*dely) &
						-new*((-4*Vm(xx,yy))/(delx**2) + (Vm(xx,yy+1)+Vm(xx,yg)-2*Vm(xx,yy))/(dely**2))
						
					Gc(xx,yp) = ((Vm(xx,yp)+Vm(xx,yp+1))**2 - (Vm(xx,yp)+Vm(xx,yy))**2)/(4*dely) &
						-new*((-4*Vm(xx,yp))/(delx**2) + (Vm(xx,yp+1)+Vm(xx,yy)-2*Vm(xx,yp))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 12) then ! down and left walls 
					Fc(xx,yy) = -new*((2*Um(xx+1,yy))/(delx**2))
						
					Fc(xp,yy) = ((Um(xp+1,yy)+Um(xp,yy))**2 - (Um(xx,yy)+Um(xp,yy))**2)/(4*delx) &
						+(0.0d0 - (Um(xp,yg)+Um(xp,yy))*(Vm(xp,yy)+Vm(xx,yy)))/(4*dely) &
						-new*((Um(xp+1,yy)+Um(xx,yy)-2*Um(xp,yy))/(delx**2) + (Um(xp,yg)-3*Um(xp,yy))/(dely**2))
						
					Gc(xx,yy) = ((Vm(xx,yy)+Vm(xx,yy+1))**2 - (Vm(xx,yy)+Vm(xx,yg))**2)/(4*dely) &
						+((Um(xx+1,yy)+Um(xx+1,yg))*(Vm(xp,yy)+Vm(xx,yy)) - 0.0d0)/(4*delx) &
						-new*((Vm(xp,yy)-3*Vm(xx,yy))/(delx**2) + (Vm(xx,yy+1)+Vm(xx,yg)-2*Vm(xx,yy))/(dely**2))
						
					Gc(xx,yp) = -new*((2*Vm(xx,yy))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 14) then ! right, down and left walls 
					Fc(xx,yy) = 0.0d0
						
					Fc(xp,yy) = 0.0d0
						
					Gc(xx,yy) = ((Vm(xx,yy)+Vm(xx,yy+1))**2 - (Vm(xx,yy)+Vm(xx,yg))**2)/(4*dely) &
						-new*((-4*Vm(xx,yy))/(delx**2) + (Vm(xx,yy+1)+Vm(xx,yg)-2*Vm(xx,yy))/(dely**2))
						
					Gc(xx,yp) = -new*((2*Vm(xx,yy))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 13) then ! up, down and left walls 
					Fc(xx,yy) = -new*((2*Um(xx+1,yy))/(delx**2))
						
					Fc(xp,yy) = ((Um(xp+1,yy)+Um(xp,yy))**2 - (Um(xx,yy)+Um(xp,yy))**2)/(4*delx) &
						-new*((Um(xp+1,yy)+Um(xx,yy)-2*Um(xp,yy))/(delx**2) + (-4*Um(xp,yy))/(dely**2))
						
					Gc(xx,yy) = 0.0d0
						
					Gc(xx,yp) = 0.0d0
				end if 
				
				if (Tp2(xx,yy) == 11) then ! up, right and left walls 
					Fc(xx,yy) = 0.0d0
						
					Fc(xp,yy) = 0.0d0
						
					Gc(xx,yy) = -new*((2*Vm(xx,yy+1))/(dely**2))
						
					Gc(xx,yp) = ((Vm(xx,yp)+Vm(xx,yp+1))**2 - (Vm(xx,yp)+Vm(xx,yy))**2)/(4*dely) &
						-new*((-4*Vm(xx,yp))/(delx**2) + (Vm(xx,yp+1)+Vm(xx,yy)-2*Vm(xx,yp))/(dely**2))
				end if 
				
				if (Tp2(xx,yy) == 7) then ! up, right and down walls 
					Fc(xx,yy) = ((Um(xx+1,yy)+Um(xx,yy))**2 - (Um(xg,yy)+Um(xx,yy))**2)/(4*delx) &
						-new*((Um(xx+1,yy)+Um(xg,yy)-2*Um(xx,yy))/(delx**2) + (-4*Um(xx,yy))/(dely**2))
						
					Fc(xp,yy) = -new*((2*Um(xx,yy))/(delx**2))
						
					Gc(xx,yy) = 0.0d0
						
					Gc(xx,yp) = 0.0d0
				end if 
				
				if (Tp2(xx,yy) == 15) then ! up, right, down and left walls 
					Fc(xx,yy) = 0.0d0
						
					Fc(xp,yy) = 0.0d0
						
					Gc(xx,yy) = 0.0d0
						
					Gc(xx,yp) = 0.0d0
				end if 
			end if 
		end do
	end do 
	
	! calculate P using sparse solver
	do xx = 1, nx
		do yy = 1, ny+2
			xp = xx + 1
			xg = xx - 1
			yp = yy + 1
			yg = yy - 1
			if (xp == nx+1) xp = 1
			if (xg == 0) xg = nx
			if (yp == ny+3) yp = yy
			if (yg == 0) yg = yy
			Rc(xx,yy) = -(Fc(xp,yy)-Fc(xx,yy))/delx - (Gc(xx,yp)-Gc(xx,yy))/dely + Dc(xx,yy)/delt
		end do
	end do
	! making matrix
	ai = 0
	ap = 0
	ax = 0.0d0
	bx = 0.0d0
	! filling matrix
	ap(1) = 0
	do pp = 1, np
		temp1 = 0
		temp2 = 0.0d0
		temp3 = 0
		temp1(1) = pp
		temp2(1) = temp2(1) -4.0d0/(delh**2)
		xx = Cloc(CP(pp),1)
		yy = Cloc(CP(pp),2)
        xp = xx + 1  ! (20180912 added YK)
        ! yp = yy + 1  ! (20180912 added YK)
        if (xp == nx+1) xp = 1  ! (20180912 added YK)
        ! if (yp == ny+3) yp = yy  ! (20180912 added YK)
		bx(pp) = bx(pp) + Rc(xx,yy)
		cnt = 1
		do ppp = 1,4
			if (CDN(pp,ppp) == CP(pp)) print *, "error: pointing itself"
			if (CDNt(pp,ppp) == 0) then
				cnt = cnt + 1
				temp1(cnt) = CPi(CDN(pp,ppp))
				temp2(cnt) = temp2(cnt) + 1.0d0/(delh**2)
			else if (CDNt(pp,ppp) == 200) then  ! lower boundary 
				if (.not.const_bot) temp2(1) = temp2(1) + 1.0d0/(delh**2)  ! no pressure gradient in y 
				if (const_bot) bx(pp) = bx(pp) - Pbot/(delh**2)   ! constant pressure
			else if (CDNt(pp,ppp) == 300) then   !  upper boundary
				if (.not.const_top) temp2(1) = temp2(1) + 1.0d0/(delh**2)  ! no pressure gradient in y 
				if (const_top) bx(pp) = bx(pp) - Ptop/(delh**2)    ! constant pressure
			else if (CDNt(pp,ppp) >=500) then    ! const. flow boundary
				temp2(1) = temp2(1) + 1.0d0/(delh**2)
				if (ppp == 1) bx(pp) = bx(pp) - (-delx*Fc(xp,yy)/(delh**2))  !  xx+1 --> xp  20180912
				if (ppp == 2) bx(pp) = bx(pp) - (delx*Fc(xx,yy)/(delh**2))
				if (ppp == 3) bx(pp) = bx(pp) - (-dely*Gc(xx,yy+1)/(delh**2))  !  yy+1 --> yp  20180912
				if (ppp == 4) bx(pp) = bx(pp) - (dely*Gc(xx,yy)/(delh**2))
			else ! CDNt(pp,ppp) is 100-115, i.e., obstacles
				temp2(1) = temp2(1) + 1.0d0/(delh**2)
				if (ppp == 1) then 
					if (xx == nx) then
						bx(pp) = bx(pp) - (-delx*Fc(1,yy)/(delh**2))
					else 
						bx(pp) = bx(pp) - (-delx*Fc(xx+1,yy)/(delh**2))  !  xx+1 --> xp  20180912
					end if 
				else if (ppp == 2) then 
					bx(pp) = bx(pp) - (delx*Fc(xx,yy)/(delh**2))
				else if (ppp == 3) then 
					bx(pp) = bx(pp) - (-dely*Gc(xx,yy+1)/(delh**2))  !  yy+1 --> yp  20180912
				else if (ppp == 4) then 
					bx(pp) = bx(pp) - (dely*Gc(xx,yy)/(delh**2))
				end if 
			end if 
		end do 
		ap(pp+1) = ap(pp) + cnt
		call heapsort2(cnt,temp1(1:cnt),temp3(1:cnt))
        if (cnt==1) temp3 = 1  ! (20180912 added YK)
		do ppp = 1, cnt
			ai(ap(pp)+ppp) = temp1(ppp) - 1  !! matrix index must start with 0 in UMFPACK
			ax(ap(pp)+ppp) = temp2(temp3(ppp))
		end do 
	end do
	! chk matrix 
	if (chk_matrix) then 
		open(500, file="chck_ap.txt", status = 'replace')
		do pp = 1, np+1
			write(500,*) ap(pp)
		end do 
		close(500)

		open(500, file="chck_ai.txt", status = 'replace')
		do pp = 1, nnz
			write(500,*) ai(pp)
		end do 
		close(500)
	
		open(500, file="chck_ax.txt", status = 'replace')
		do pp = 1, nnz
			write(500,*) ax(pp)
		end do 
		close(500)
	
		open(500, file="chck_bx.txt", status = 'replace')
		do pp = 1, n
			write(500,*) bx(pp)
		end do 
		close(500)
	end if
	! solving matrix with UMFPACK (following is pasted from umfpack_simple.f90)
	! Set the default control parameters.
	call umf4def( control )
	! From the matrix data, create the symbolic factorization information.
	call umf4sym ( n, n, ap, ai, ax, symbolic, control, info )
	if ( info(1) < 0.0D+00 ) then
		write ( *, * ) ''
		write ( *, *) 'UMFPACK_SIMPLE - Fatal error!'
		write ( *, * ) '  UMF4SYM returns INFO(1) = ', info(1)
		stop 1
	end if
	! From the symbolic factorization information, carry out the numeric factorization.
	call umf4num ( ap, ai, ax, symbolic, numeric, control, info )
	if ( info(1) < 0.0D+00 ) then
		write ( *, '(a)' ) ''
		write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
		write ( *, '(a,g14.6)' ) '  UMF4NUM returns INFO(1) = ', info(1)
		stop 1
	end if
	!  Free the memory associated with the symbolic factorization.
	call umf4fsym ( symbolic )
	! Solve the linear system.
	sys = 0
	call umf4sol ( sys, kai, bx, numeric, control, info )
	if ( info(1) < 0.0D+00 ) then
		write ( *, '(a)' ) ''
		write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
		write ( *, '(a,g14.6)' ) '  UMF4SOL returns INFO(1) = ', info(1)
		stop 1
	end if
	! Free the memory associated with the numeric factorization.
	call umf4fnum ( numeric )
	! chk the solution 
	if (chk_matrix) then 
		open(500, file="chck_kai.txt", status = 'replace')
		do pp = 1, np
			write(500,*) kai(pp)
		end do 
		close(500)
	end if
	! solution put into P
	do pp = 1, np
		xx = Cloc(CP(pp),1)
		yy = Cloc(CP(pp),2)
		Pc(xx,yy) = kai(pp)
        if (chk_particles) then 
          if (yy==1 .or. yy==ny+2) then
            cycle
          else 
            if (flw_rsltn==1) then 
                if (matrix(yy-1,xx)%class==p) then 
                   write(*,*) time,"error: P is calculated on particles x,y=",xx,yy
                   write(File_log,*) time,"error: P is calculated on particles x,y=",xx,yy
                   flg_stop=.true.
                   ! stop
                endif 
            elseif (flw_rsltn/=1) then 
                do x_tmp = 1,flw_rsltn
                    if (mod(xx-x_tmp,flw_rsltn)==0) exit 
                enddo 
                do y_tmp = 1,flw_rsltn
                    if (mod(yy-y_tmp-1,flw_rsltn)==0) exit 
                enddo 
                xxx = (xx-x_tmp)/flw_rsltn + 1
                yyy = (yy-y_tmp-1)/flw_rsltn + 1
                if (matrix(yyy,xxx)%class==p) then 
                   write(*,*) time,"error: P is calculated on particles x,y=",xxx,yyy
                   write(File_log,*) time,"error: P is calculated on particles x,y=",xxx,yyy
                   flg_stop=.true.
                   ! stop
                endif 
            endif 
          endif
        endif
	end do 
	! P boundary 
	do xx = 1,nx
		do yy = 1,ny+2
			xp = xx + 1
			xg = xx - 1
			yp = yy + 1
			yg = yy - 1
			if (xp == nx+1) xp = 1
			if (xg == 0) xg = nx
			if (yp == ny+3) yp = yy
			if (yg == 0) yg = yy
			! lower boundary 
			if (Tp(xx,yy) == 200) then 
				if (.not.const_bot) then
					Pc(xx,yy) = Pc(xx,yg)  ! no pressure gradient in y 
				else if (const_bot) then 
					Pc(xx,yy) = Pbot  ! const. pressure
				end if 
			end if 
			! upper boundary
			if (Tp(xx,yy) == 300) then
				if (const_top) then 
					Pc(xx,yy) = Ptop    ! const. pressure    
				else if (.not.const_top) then
					Pc(xx,yy) = Pc(xx,yp)    ! no pressure gradient in y 
				end if 
			end if 
			
			if (Tp(xx,yy) == 500) Pc(xx,yy) = Pc(xp,yy) + delx*Fc(xp,yy)  ! rightward flow
			if (Tp(xx,yy) == 600) Pc(xx,yy) = Pc(xg,yy) - delx*Fc(xx,yy)  !  left
			if (Tp(xx,yy) == 700) Pc(xx,yy) = Pc(xx,yp) + dely*Gc(xx,yp)  !  down
			if (Tp(xx,yy) == 800) Pc(xx,yy) = Pc(xx,yg) - dely*Gc(xx,yy)  !  up
			
			if (Tp(xx,yy) == 900) then 
				if (Tp(xp,yy) == 0) Pc(xx,yy) = Pc(xp,yy) + delx*Fc(xp,yy)
				if (Tp(xg,yy) == 0) Pc(xx,yy) = Pc(xg,yy) - delx*Fc(xx,yy)
				if (Tp(xx,yp) == 0) Pc(xx,yy) = Pc(xx,yp) + dely*Gc(xx,yp)
				if (Tp(xx,yg) == 0) Pc(xx,yy) = Pc(xx,yg) - dely*Gc(xx,yy)
			end if 
		end do
	end do
	! show results on display
	if (show_display) then 
        print '(A,I3,A)','        Flow time-step: ',step,' ---------------------------'
        print '(5A11)', '','U[cm/s]','V[cm/s]','P[x0.1Pa]','D[/s]'
		print '(A11,4E11.3)', 'MAX',maxval(Um(:,2:ny+1)),maxval(Vm(:,2:ny+1)),maxval(Pc(:,2:ny+1)),maxval(Dc(:,2:ny+1))
		print '(A11,4E11.3)', 'MIN',minval(Um(:,2:ny+1)),minval(Vm(:,2:ny+1)),minval(Pc(:,2:ny+1)),minval(Dc(:,2:ny+1))
		print '(4A11,E11.3)', 'SUM','-','-','-',sum(Dc(:,2:ny+1))
		print '(A)','        -----------------------------------------------'
        print *
	end if 
    if (calc_perm) then 
		Vt_t = sum((Vm(:,2) + Vm(:,3))/2.0d0)
		Vt_b = sum((Vm(:,ny+1) + Vm(:,ny+2))/2.0d0 )
		perm_u = (Vt_t)/nx*new*ny*delh/(Ptop-Pbot)*1e-4   
		perm_d = (Vt_b)/nx*new*ny*delh/(Ptop-Pbot)*1e-4  
        if (show_display) then 
            print '(A)','        *****************************'
            print '(A15,2A11)','','UP','DOWN'
            print '(A8,A7,2E11.3)', '','K[m2]: ',perm_u,perm_d
            print '(A)','        *****************************'
            print *
        endif 
    endif 
	! end calculation 1 step
	step = step + 1
end do
Ug = Um
Vg = Vm 
Pg = Pc
Dg = Dc
Uc(:,:) = (Um(1:nx,:) + Um(2:nx+1,:))/2.0d0
Vc(:,:) = (Vm(:,1:ny+2) + Vm(:,2:ny+3))/2.0d0 
if (flw_rsltn == 1) then 
    Uo(:,:) = Uc(:,2:ny+1)*60.d0*60.d0*24.0d0*365.0d0  !! cm/yr
    Vo(:,:) = Vc(:,2:ny+1)*60.d0*60.d0*24.0d0*365.0d0  !! cm/yr
    Uom(:,:) = Um(:,2:ny+1)*60.d0*60.d0*24.0d0*365.0d0  !! cm/yr
    Vom(:,:) = Vm(:,2:ny+2)*60.d0*60.d0*24.0d0*365.0d0  !! cm/yr
elseif (flw_rsltn/=1) then 
    do yy=1,ny/flw_rsltn
        do xx=1,nx/flw_rsltn
            Uo(xx,yy) = 0
            Vo(xx,yy) = 0
            Uom(xx,yy) = 0
            Vom(xx,yy) = 0
            do xxx=1,flw_rsltn
                do yyy = 1,flw_rsltn
                    Uo(xx,yy) = Uo(xx,yy) + Uc((xx-1)*flw_rsltn+1+xxx-1,(yy-1)*flw_rsltn+1+1+yyy-1)*60.d0*60.d0*24.0d0*365.0d0  !! cm/yr
                    Vo(xx,yy) = Vo(xx,yy) + Vc((xx-1)*flw_rsltn+1+xxx-1,(yy-1)*flw_rsltn+1+1+yyy-1)*60.d0*60.d0*24.0d0*365.0d0  !! cm/yr
                    Uom(xx,yy) = Uom(xx,yy) + Um((xx-1)*flw_rsltn+1+xxx-1,(yy-1)*flw_rsltn+1+1+yyy-1)*60.d0*60.d0*24.0d0*365.0d0  !! cm/yr
                    Vom(xx,yy) = Vom(xx,yy) + Vm((xx-1)*flw_rsltn+1+xxx-1,(yy-1)*flw_rsltn+1+1+yyy-1)*60.d0*60.d0*24.0d0*365.0d0  !! cm/yr
                enddo
            enddo
            Uo(xx,yy) = Uo(xx,yy)/flw_rsltn/flw_rsltn
            Vo(xx,yy) = Vo(xx,yy)/flw_rsltn/flw_rsltn
            Uom(xx,yy) = Uom(xx,yy)/flw_rsltn/flw_rsltn
            Vom(xx,yy) = Vom(xx,yy)/flw_rsltn/flw_rsltn
        enddo
    enddo
    yy = ny/flw_rsltn + 1
    do xx=1,nx/flw_rsltn
        Vom(xx,yy) = 0
            do xxx=1,flw_rsltn
                yyy = ny+1 + 1
                Vom(xx,yy) = Vom(xx,yy) + Vm((xx-1)*flw_rsltn+1+xxx-1,yyy)*60.d0*60.d0*24.0d0*365.0d0  !! cm/yr
            enddo
        Vom(xx,yy) = Vom(xx,yy)/flw_rsltn
    enddo
    xx = nx/flw_rsltn + 1
    do yy=1,ny/flw_rsltn
        Uom(xx,yy) = 0
            do yyy=1,flw_rsltn
                xxx = nx+1 
                Uom(xx,yy) = Uom(xx,yy) + Um(xxx,(yy-1)*flw_rsltn+1+1+yyy-1)*60.d0*60.d0*24.0d0*365.0d0  !! cm/yr
            enddo
        Uom(xx,yy) = Uom(xx,yy)/flw_rsltn
    enddo
endif
if (chk_particles) then 
  do yy = 1,ny/flw_rsltn
    do xx=1,nx/flw_rsltn
      if (matrix(yy,xx)%class ==p) then 
        if (Uo(xx,yy)==0 .and. Vo(xx,yy)==0) then
          cycle
        else 
          print *, time,'error in NS module check for flow in particles at x,y = ',xx,yy     &
                  , 'uo,vo=', Uo(xx,yy),Vo(xx,yy)!, 'B2, Tp, Tp2 =',B2(xx,yy+1),Tp(xx,yy+1),Tp2(xx,yy+1)
          write(File_log, *) time,'error in NS module check for flow in particles at x,y = ',xx,yy     &
                  , 'uo,vo=', Uo(xx,yy),Vo(xx,yy)!, 'B2, Tp, Tp2 =',B2(xx,yy+1),Tp(xx,yy+1),Tp2(xx,yy+1)
          flg_stop=.true.
          ! stop
       endif 
     endif 
   enddo
 enddo       
endif            

! check matrix file vs type/binary files
if (chk_particles .and. flg_stop) then 
	open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/matrix-class-'  & 
		//trim(adjustl(numtemp))//'.txt', status = 'replace')
	do yy = 1,ny
		write(100,*) (matrix(yy,xx)%class, xx = 1,nx)
	end do 
	close(100)
    
	open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/Tp(nxxny)-'  & 
		//trim(adjustl(numtemp))//'.txt', status = 'replace')
	do yy = 2,ny+1
		write(100,*) (Tp(xx,yy), xx = 1,nx)
	end do 
	close(100)

	open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/Tp2(nxxny)-'  & 
		//trim(adjustl(numtemp))//'.txt', status = 'replace')
	do yy = 2,ny+1
		write(100,*) (Tp2(xx,yy), xx = 1,nx)
	end do 
	close(100)

	open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/B2(nxxny)-'  & 
		//trim(adjustl(numtemp))//'.txt', status = 'replace')
	do yy = 2,ny+1
		write(100,*) (B2(xx,yy), xx = 1,nx)
	end do 
	close(100)
    
    STOP !!!!!!!!!!!
    
end if 

if (calc_perm) then 
    write(File_Perm,*) Time*Timescale, perm_u,perm_d
endif 

end subroutine flow_calc
! =================================================================
subroutine output_flow()
use GlobalVariables
implicit none 
integer(kind=4) :: nx , ny , xx, yy
real(kind=8)  :: delx, dely
Character*21 numtemp, cat_tmp 

real(kind=8) ,allocatable :: Pc(:,:)
real(kind=8) ,allocatable :: Uc(:,:)
real(kind=8) ,allocatable :: Vc(:,:)
real(kind=8) ,allocatable :: Um(:,:)
real(kind=8) ,allocatable :: Vm(:,:)
real(kind=8) ,allocatable :: Dc(:,:)
real(kind=8) ,allocatable :: Qcx(:,:)
real(kind=8) ,allocatable :: Qcy(:,:)
real(kind=8) ,allocatable :: vabs(:,:)
real(kind=8)  :: Ut_l, Ut_r, Vt_t, Vt_b

logical :: show_display = .false.

write(numtemp,'(i10.1)') Time

nx = n_col*flw_rsltn
ny = n_row*flw_rsltn  

if (.not.allocated(Pc)) allocate(Pc(nx,ny+2))
if (.not.allocated(Uc)) allocate(Uc(nx,ny+2))
if (.not.allocated(Vc)) allocate(Vc(nx,ny+2))
if (.not.allocated(vabs)) allocate(vabs(nx,ny+2))
if (.not.allocated(Um)) allocate(Um(nx+1,ny+2))
if (.not.allocated(Vm)) allocate(Vm(nx,ny+2+1))
if (.not.allocated(Dc)) allocate(Dc(nx,ny+2))
if (.not.allocated(Qcx)) allocate(Qcx(nx,ny+2))
if (.not.allocated(Qcy)) allocate(Qcy(nx,ny+2))

delx = pixelsize/flw_rsltn
dely = delx

Um = Ug
Vm = Vg
Pc = Pg
Dc = Dg

if (calc_perm) then
    cat_tmp = 'geo'
else 
    cat_tmp = 'o2'
endif

open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-Um-'  & 
	//trim(adjustl(numtemp))//'.txt', status = 'replace')
do yy = 1,ny+2
	write(100,*) (Um(xx,yy), xx = 1,nx+1)
end do 
close(100)
open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-Vm-'  & 
	//trim(adjustl(numtemp))//'.txt', status = 'replace')
do yy = 1,ny+2+1
	write(100,*) (Vm(xx,yy), xx = 1,nx)
end do 
close(100)
Uc(:,:) = (Um(1:nx,:) + Um(2:nx+1,:))/2.0d0
Vc(:,:) = (Vm(:,1:ny+2) + Vm(:,2:ny+3))/2.0d0 
vabs = sqrt(Uc*Uc + Vc*Vc)
Ut_l = sum(Uc(1,:))
Ut_r = sum(Uc(nx,:))
Vt_t = sum(Vc(:,2))
Vt_b = sum(Vc(:,ny+2))
Qcx = 0.0d0
do yy = 1, ny+2
	if (yy == 1) cycle
	Qcx(1,yy) = Qcx(1,yy-1) + dely*(Um(1,yy-1)+Um(2,yy)+Um(1,yy)+Um(2,yy-1))/4.0d0
	do xx = 1, nx
		if (xx == 1) cycle
		Qcx(xx,yy) = Qcx(xx-1,yy) - delx*(Vm(xx-1,yy) + Vm(xx-1, yy+1) &
                 +Vm(xx,yy) + Vm(xx, yy+1))/4.0d0
	end do
end do 
Qcy = 0.0d0
do xx = 1, nx
	if (xx == 1) cycle
	Qcy(xx,1) = Qcy(xx-1,1) - delx*(Vm(xx-1,1) + Vm(xx, 2)+Vm(xx,1) + Vm(xx-1, 2))/4.0d0
	do yy = 1, ny + 2
		if (yy == 1) cycle
		Qcy(xx,yy) = Qcy(xx,yy-1) + dely*(Um(xx,yy-1) + Um(xx+1,yy-1)  &
                     +Um(xx,yy) + Um(xx+1,yy))/4.0d0
	end do 
end do

open(unit = 100, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-Uc-'  & 
	//trim(adjustl(numtemp))//'.txt', status = 'replace')
open(unit = 200, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-Vc-'  & 
	//trim(adjustl(numtemp))//'.txt', status = 'replace')
open(unit = 300, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-Pc-'  & 
	//trim(adjustl(numtemp))//'.txt', status = 'replace')
open(unit = 400, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-Dc-'  & 
	//trim(adjustl(numtemp))//'.txt', status = 'replace')
open(unit = 500, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-vabs-'  & 
	//trim(adjustl(numtemp))//'.txt', status = 'replace')
open(unit = 600, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-Qcx-'  & 
	//trim(adjustl(numtemp))//'.txt', status = 'replace')
open(unit = 700, file = trim(adjustl(today))//'/'//trim(adjustl(cat_tmp))//'/test2d(F)-Qcy-'  & 
	//trim(adjustl(numtemp))//'.txt', status = 'replace')
do yy = 2,ny+1
	write(100,*) (Uc(xx,yy), xx = 1,nx)
	write(200,*) (Vc(xx,yy), xx = 1,nx)
	write(300,*) (Pc(xx,yy), xx = 1,nx)
	write(400,*) (Dc(xx,yy), xx = 1,nx)
	write(500,*) (vabs(xx,yy), xx = 1,nx)
	write(600,*) (Qcx(xx,yy), xx = 1,nx)
	write(700,*) (Qcy(xx,yy), xx = 1,nx)
end do 
close(100)
close(200)
close(300)
close(400)
close(500)
close(600)
close(700)
close(800)
if (show_display ) then
	print *,'totU:',Ut_l,'/',Ut_r,'/',abs((Ut_l-Ut_r)/Ut_r*100),'/'
	print *,'totV:',Vt_t,'/',Vt_b,'/',abs((Vt_t-Vt_b)/Vt_t*100),'/'
	print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
end if 

end subroutine output_flow
! ====================================================
subroutine heapsort2(n,array,turn)
!!!  from http://slpr.sakura.ne.jp/qp/sortf90/
  implicit none
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(out)::turn(1:n)
  integer(kind=4),intent(inout)::array(1:n)
 
  integer(kind=4)::i,k,j,l,m
  integer(kind=4)::t
 
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
end subroutine heapsort2
! =====================================================
end module