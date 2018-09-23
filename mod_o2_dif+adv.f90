module o2_diffusion


use Globalvariables

implicit none
 
contains

   ! *********************************************************

   SubRoutine Output_O2txtImg()

   Use GlobalVariables
   Implicit None
   Integer :: X, Y
   Character*21 numtemp
   real :: txtimg(N_Col)
    
   write(numtemp,'(i10.1)') Time
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/O2data-'          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO Y = 1, N_row
     txtimg = 0 
     DO X = 1, N_col
         txtimg(X) = O2(y,x)%oxygen
     END Do
     write(File_txtImg, *) (txtimg(x), x = 1, n_col)
     END Do
   Close(File_txtimg)
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/O2_use_data-'          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO Y = 1, N_row
     txtimg = 0 
     DO X = 1, N_col
         txtimg(X) = O2(y,x)%oxygen_use
     END Do
     write(File_txtImg, *) (txtimg(x), x = 1, n_col)
     END Do
   Close(File_txtimg)
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/O2_chan_data-'          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO Y = 1, N_row
     txtimg = 0 
     DO X = 1, N_col
         txtimg(X) = O2(y,x)%value_pre
     END Do
     write(File_txtImg, *) (txtimg(x), x = 1, n_col)
     END Do
   Close(File_txtimg)
   
   OPEN(unit = File_txtImg, File = trim(adjustl(today))//'/data/O2_dif-'          &
                        //trim(adjustl(numtemp))//'.txt', status = 'unknown')
     DO Y = 1, N_row
     txtimg = 0 
     DO X = 1, N_col
         txtimg(X) = edif(y,x)
     END Do
     write(File_txtImg, *) (txtimg(x), x = 1, n_col)
     END Do
   Close(File_txtimg)

   End Subroutine Output_O2txtimg
   
   ! *********************************************************
   
   subroutine O2i_setup()
   
   use globalvariables
   
   implicit none
   
   integer :: x, y
   real :: porosity 
   logical :: O2_sat_initial
   
     O2(:,:)%Oxygen = 0.
     O2(:,:)%Oxygen_pre = 0.
     O2(:,:)%Oxygen_use = 0.
	 
   O2(1,:)%Oxygen = 1.
   O2(1,:)%Oxygen_pre = 1.
   
   O2_sat_initial = .true.
   ! O2_sat_initial = .false.
   
   if (O2_sat_initial) then 
   
   do y = 1, n_row
   
     do x = 1, N_Col
	 
	 if (matrix(y,x)%class==p) cycle
     O2(y,x)%Oxygen = 1.
     O2(y,x)%Oxygen_pre = 1.
     O2(y,x)%Oxygen_use = 0.
	 
     
     end do
   
   end do
   
   end if
   
   end subroutine O2i_setup
   
   ! *********************************************************
   
   subroutine O2pre_setup()
   
   use globalvariables
   
   implicit none
   
   integer :: x, y, i
   
   integer :: up(n_row,n_col), down(n_row,n_col)
   
   up = 0
   down = 0
   
   O2%value_pre = 0
   
   up(1,:) = 1
   down(n_row,:) = 1
   
   do y = 2, n_row
   
     do x = 1, N_Col
          
       if (matrix(y,x)%class == p ) cycle
	   
       if (up(y-1,x)/=0) up(y,x) = 1
	   
	   if (x == 1) then
	     if (up(y,2)/=0) up(y,x) = 1
	     if (up(y,n_col)/=0) up(y,x) = 1
	   else if (x==n_col) then  
	     if (up(y,1)/=0) up(y,x) = 1
	     if (up(y,n_col-1)/=0) up(y,x) = 1
	   else
	     if (up(y,x+1)/=0) up(y,x) = 1
	     if (up(y,x-1)/=0) up(y,x) = 1
	   end if 
		 
     end do
	 
     do x = n_col, 1, -1
          
       if (matrix(y,x)%class == p ) cycle
	   
	   if (x == 1) then
	     if (up(y,2)/=0) up(y,x) = 1
	     if (up(y,n_col)/=0) up(y,x) = 1
	   else if (x==n_col) then  
	     if (up(y,1)/=0) up(y,x) = 1
	     if (up(y,n_col-1)/=0) up(y,x) = 1
	   else
	     if (up(y,x+1)/=0) up(y,x) = 1
	     if (up(y,x-1)/=0) up(y,x) = 1
	   end if 
		 
     end do
   
   end do
   
   do y = n_row-1, 1, -1
   
     do x = 1, N_Col
       
       if (matrix(y,x)%class == p ) cycle
	   
       if (down(y+1,x)/=0) down(y,x) = 1
	   
	   if (x == 1) then
	     if (down(y,2)/=0) down(y,x) = 1
	     if (down(y,n_col)/=0) down(y,x) = 1
	   else if (x==n_col) then  
	     if (down(y,1)/=0) down(y,x) = 1
	     if (down(y,n_col-1)/=0) down(y,x) = 1
	   else
	     if (down(y,x+1)/=0) down(y,x) = 1
	     if (down(y,x-1)/=0) down(y,x) = 1
	   end if 
		 
     end do
	 
     do x = n_col, 1, -1
          
       if (matrix(y,x)%class == p ) cycle
	   
	   if (x == 1) then
	     if (down(y,2)/=0) down(y,x) = 1
	     if (down(y,n_col)/=0) down(y,x) = 1
	   else if (x==n_col) then  
	     if (down(y,1)/=0) down(y,x) = 1
	     if (down(y,n_col-1)/=0) down(y,x) = 1
	   else
	     if (down(y,x+1)/=0) down(y,x) = 1
	     if (down(y,x-1)/=0) down(y,x) = 1
	   end if 
		 
     end do
   
   end do
   
   where ((up/=0).or.(down/=0)) O2%value_pre = 1
   
   end subroutine O2pre_setup

   ! *********************************************************
   
   subroutine oxygen_profile()
   
   !!!!
   !! solving   - [c(x,t)-c(x,t-1)]/dt + <D[c(x+1,t)-c(x,t)]/dx - D[c(x,t)-c(x-1,t)]/dx>/dx - v[c(x,t)-c(x-1,t)]/dx (if v >0) - R = 0
   !!  
   !!  (using 1st order upwind scheme for advection, 2nd order central differenciation for diffusion term. Time is explicit, 1st order)
   
   use GlobalVariables
   use ieee_arithmetic
   
   implicit none
    
   integer ( kind = 4 ) :: n, nnz
   
   integer ( kind = 4 ), allocatable :: ai( : )  ! row number where /=0
   integer ( kind = 4 ), allocatable  :: ap(:)   ! number of non-zero component at each column
   real ( kind = 8 ), allocatable  :: ax( : )   !  component where /=0
   real ( kind = 8 ), allocatable :: b( :) 
   real  ( kind = 8 ) :: control(20)
   integer ( kind = 4 ) :: filenum
   integer  ( kind = 4 ) :: i
   real ( kind = 8 ) :: info(90)
   integer ( kind = 8 ) :: numeric
   integer ( kind = 4 ) :: status
   integer ( kind = 8 ) :: symbolic
   integer  ( kind = 4 ) :: sys
   real ( kind = 8 ), allocatable  :: x(:)
   
   integer :: n_row_calc       !  number of row used for the calculation
   integer :: y_fin       !  bottom y layer
   
   real :: org_dif_fact = 1.  !  porosity of particle
   integer :: j, k, cnt, cnt2, xx, yy, cnt3, cnt4, cnt5
   logical :: mtx_chk
   
   real :: ox_l, ox_r, ox_u, ox_d, ox_c
   real :: po_l, po_r, po_u, po_d, po_c, po_pre
   real :: dif_l, dif_r, dif_u, dif_d, dif_c
   real :: rxn
   
   real :: dt, dx
   
   integer :: cnt_rec(n_row,n_col)
   
   real, parameter :: imp = 1.0
   real, parameter :: epl = 1- imp
   real, parameter :: cntgrad = 1.0
   
   integer :: step, totstep
   
   real :: tmpU,tmpV
   
   logical :: divide = .false.
   ! logical :: divide = .true.
   
   ! logical :: only_sed
   
   real, parameter :: rx_ex = 0.0  ! 1.0 if reaction is explicitly reflected
   real, parameter :: rx_imp = 1.0 - rx_ex  ! 1.0 if reaction is implicitly reflected
   
   real ::  shear, visc
   integer :: sedloc
   
   
   y_fin = n_row   
   if (only_sed) then 
   do yy = 1, n_row
      if (any(matrix(yy,:)%class==p)) exit
   end do 
   y_int = yy - 3
   
   else 
   y_int  = 1       !  when upper water column is boudary 
   end if 
   
   y_int = 1
   
   
   n_row_calc = y_fin - (y_int - 1)
   
   if (n_row_calc /= y_fin - y_int + 1) write(*,*) "fatal error"
   
   cnt = 0    ! count the number of water 
   cnt2 = 0    ! count necessary non-zero matrix components
   cnt_rec = 0
   do yy = y_int  + 1, y_fin
   
     do xx = 1, n_col
	 
	 if (matrix(yy,xx)%class /=p) then
       
       cnt = cnt + 1	
       cnt2 = cnt2 + 1	  
       
       cnt_rec(yy,xx) = cnt	   
	   
	   if (yy == y_int+1) then 
	   
	     if (xx == 1) then 
		 
		   if (matrix(yy,n_col)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy,2)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy+1,1)%class /= p) cnt2 = cnt2 + 1
		 
		 else if (xx == n_col) then 
		 
		   if (matrix(yy,1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy,n_col-1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy+1,n_col)%class /= p) cnt2 = cnt2 + 1
		 
		 else 
		 
		   if (matrix(yy,xx-1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy,xx+1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy+1,xx)%class /= p) cnt2 = cnt2 + 1
		 
		 end if 
	   
	   else if (yy == y_fin) then 
	   
	     if (xx == 1) then 
		 
		   if (matrix(yy,n_col)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy,2)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy-1,1)%class /= p) cnt2 = cnt2 + 1
		 
		 else if (xx == n_col) then 
		 
		   if (matrix(yy,1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy,n_col-1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy-1,n_col)%class /= p) cnt2 = cnt2 + 1
		 
		 else 
		 
		   if (matrix(yy,xx-1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy,xx+1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy-1,xx)%class /= p) cnt2 = cnt2 + 1
		 
		 end if 
	   
	   else 
	   
	     if (xx == 1) then 
		 
		   if (matrix(yy,n_col)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy,2)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy+1,1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy-1,1)%class /= p) cnt2 = cnt2 + 1
		 
		 else if (xx == n_col) then 
		 
		   if (matrix(yy,1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy,n_col-1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy+1,n_col)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy-1,n_col)%class /= p) cnt2 = cnt2 + 1
		 
		 else 
		 
		   if (matrix(yy,xx-1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy,xx+1)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy+1,xx)%class /= p) cnt2 = cnt2 + 1
		   if (matrix(yy-1,xx)%class /= p) cnt2 = cnt2 + 1
		 
		 end if 
	   
	   end if
	   
	 end if
	 
	 end do 

  end do
  
  
  edif = dif_0
  shear = 3.171e-2*60.*60.*24.*365. ! cm/yr (compare 1e-5 m/s; Volkenborn et al. 2012)
  ! shear = 3.171e-1*60.*60.*24.*365. ! cm/yr (compare 1e-5 m/s; Volkenborn et al. 2012)
  visc =  4.79e5   ! kinematic viscosity in cm^2/yr
  
  !! lateral change in diffusitivity allowed 
  
  do xx = 1, n_col  
     sedloc = 1
     do yy = 1, n_row
	   if (matrix(yy,xx)%class == p) then 
	     sedloc = yy
		 exit 
	   end if 
	 end do 
	 
	 do yy = sedloc, 1, -1
	   edif(yy,xx) = edif(yy,xx) + 0.4*visc/363.*((sedloc - yy)*pixelsize*shear/visc)**3
	 end do 
  end do
  
	
  
  n = cnt
  nnz = cnt2
	     
		 
   
   Allocate (ai(nnz),ap(n+1),ax(nnz),b(n),x(n))
   
   totstep = 1
   
   if (divide) then 
   do i = 1,N_ind
      totstep = max(totstep, 5*ceiling(sqrt(maxval(Vb(:,i))**2+maxval(Vb(:,i))**2)/org(i)%width))
	  print *, totstep, org(i)%width, ceiling(sqrt(maxval(Vb(:,i))**2+maxval(Vb(:,i))**2))
   end do
   end if
   
   do step = 1, totstep
   
   ai = 0
   ap = 0
   ax = 0.
   b = 0.
   
   dt = TImescale/365./real(totstep)  ! yr
   
   
   dx = PixelSIze  ! cm
   
   cnt2 = 1
   
   do yy = y_int  + 1, y_fin
   
     do xx = 1, n_col
       
	   if (matrix(yy,xx)%class == p) cycle
	   
       j = cnt_rec(yy,xx)    ! matrix number
       cnt = 0
       
       ox_c = O2(yy,xx)%Oxygen
       
	   dif_c = dif_0
       if (matrix(yy,xx)%Class >= 1) dif_c = dif_c*org_dif_fact 
		   
	   if (Uo(xx,yy) == 0.) tmpU = 1e20
	   if (Uo(xx,yy) /= 0.) tmpU = Uo(xx,yy)
	   
	   if (Vo(xx,yy) == 0.) tmpV = 1e20
	   if (Vo(xx,yy) /= 0.) tmpV = Vo(xx,yy)
       
       rxn = rx_ex*O2(yy,xx)%oxygen_use*O2(yy,xx)%oxygen
       
       b(j) = b(j) - ox_c/dt + rxn
     
       if (yy == y_int  + 1) then
           
         b(j) = b(j) -edif(yy,xx)*1./dx/dx  + imp*(Vo(xx,yy)+abs(Vo(xx,yy)))*(-1.)/dx/2.   &	 !  upper boundary condition
              + epl*((Vo(xx,yy)+abs(Vo(xx,yy)))*(O2(yy,xx)%oxygen-1.) )/dx/2.   &
              + epl*cntgrad*(Vo(xx,yy)+abs(Vo(xx,yy)))/tmpV*O2(yy,xx)%oxygen*(Vo(xx,yy)-Vo(xx,yy-1))/dx/2.   &
			  -(edif(yy,xx)-edif(yy-1,xx))*(-1.)/dx/dx
         
		 if (xx == 1) then
		   
		   k = cnt_rec(yy,xx)
		   ax(cnt2) = ax(cnt2) -edif(yy,xx)*1./dx/dx - 1./dt - imp*(Vo(xx,yy)+abs(Vo(xx,yy)))*(1.)/dx/2.   &
		               - imp*cntgrad*(Vo(xx,yy)+abs(Vo(xx,yy)))/tmpV/2.*(Vo(xx,yy)-Vo(xx,yy-1))/dx   &
					   - rx_imp*o2(yy,xx)%oxygen_use  &
					   + (edif(yy,xx)-edif(yy-1,xx))*(1.)/dx/dx
           ai(cnt2) = k 
           cnt2 = cnt2 + 1
           cnt = cnt + 1
		   
		   cnt3 = 0
		   
		   if (matrix(yy,2)%class /=p) then 
		     
			 k = cnt_rec(yy,2) 
			 ax(cnt2) = ax(cnt2) + edif(yy,2)/dx/dx - imp*(Uo(2,yy)+abs(Uo(2,yy)))*(-1./dx)/2.  &
			               +(edif(yy,2) - edif(yy,1))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Uo(xx,yy)-abs(Uo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Uo(xx,yy)-abs(Uo(xx,yy)))*(Uo(2,yy) - Uo(xx,yy))/dx/2./tmpU  &
							   + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)-abs(Uo(xx,yy)))*(O2(yy,2)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(2,yy) - Uo(xx,yy)))/dx
		   end if
		   
		   if (matrix(yy,n_col)%class /=p) then 
		     
			 k = cnt_rec(yy,n_col) 
			 ax(cnt2) = ax(cnt2) + edif(yy,n_col)/dx/dx - imp*(Uo(n_col,yy)-abs(Uo(n_col,yy)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Uo(xx,yy)+abs(Uo(xx,yy)))*(1./dx)/2. & 
			                  - imp*cntgrad*(Uo(xx,yy)+abs(Uo(xx,yy)))*(Uo(xx,yy) - Uo(n_col,yy))/dx/2./tmpU   &
							  + edif(yy,xx)*(-1.)/dx/dx + (edif(yy,xx) - edif(yy, n_col))*(1.)/dx/dx
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)+abs(Uo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy,n_col)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(n_col,yy)))/dx
		   end if
			 
		   if (matrix(yy+1,1)%class /=p) then 
		     
			 k = cnt_rec(yy+1,1) 
			 ax(cnt2) = ax(cnt2) + edif(yy+1,xx)/dx/dx - imp*(Vo(1,yy+1)+abs(Vo(1,yy+1)))*(-1./dx)/2.   &
			                 + (edif(yy+1,xx) - edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Vo(xx,yy)-abs(Vo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Vo(xx,yy)-abs(Vo(xx,yy)))*(Vo(xx,yy+1) - Vo(xx,yy))/dx/2./tmpV   &
							  + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)-abs(Vo(xx,yy)))*(O2(yy+1,xx)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)-abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy+1) - Vo(xx,yy)))/dx
		   end if
		   
         else if (xx == n_col) then
		   
		   cnt3 = 0
		   cnt4 = 0
		   
		   if (matrix(yy,1)%class /=p) then 
		     
			 k = cnt_rec(yy,1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,1)/dx/dx - imp*(Uo(1,yy)+abs(Uo(1,yy)))*(-1./dx)/2.  &
			                  + (edif(yy,1)-edif(yy,n_col))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 cnt4 = 1
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)-abs(Uo(xx,yy)))*(O2(yy,1)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(1,yy) - Uo(xx,yy)))/dx
		   end if
		   
		   if (matrix(yy,n_col - 1)%class /=p) then 
		     
			 k = cnt_rec(yy,n_col - 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,n_col-1)/dx/dx - imp*(Uo(xx-1,yy)-abs(Uo(xx-1,yy)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
		     ax(cnt2) = ax(cnt2) - imp*(Uo(xx,yy)+abs(Uo(xx,yy)))*(1./dx)/2. & 
			               - imp*cntgrad*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy)-Uo(n_col-1,yy))/dx  &
						   + edif(yy,xx)*(-1.)/dx/dx + (edif(yy,xx)- edif(yy,xx-1))*(1.)/dx/dx
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)+abs(Uo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy,n_col-1)%oxygen)/2. & 
			              + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(n_col-1,yy)))/dx
		   end if
		   
		   k = cnt_rec(yy,xx)
		   ax(cnt2) = ax(cnt2) -edif(yy,xx)*(1.)/dx/dx - 1./dt - imp*(Vo(xx,yy)+abs(Vo(xx,yy)))*(1.)/dx/2.  &
		               - imp*cntgrad*(Vo(xx,yy)+abs(Vo(xx,yy)))/tmpV/2.*(Vo(xx,yy)-Vo(xx,yy-1))/dx   &
					   - rx_imp*o2(yy,xx)%oxygen_use    &
					   + (edif(yy,xx) - edif(yy-1,xx))*(1.)/dx/dx              &
					   + cnt4*(- imp*(Uo(xx,yy)-abs(Uo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Uo(xx,yy)-abs(Uo(xx,yy)))*(Uo(1,yy) - Uo(xx,yy))/dx/2./tmpU  &
							  + edif(yy,xx)*(-1.)/dx/dx )
           ai(cnt2) = k 
           cnt2 = cnt2 + 1
           cnt = cnt + 1
		   
		   cnt3 = 0
			 
		   if (matrix(yy+1,xx)%class /=p) then 
		     
			 k = cnt_rec(yy+1,xx) 
			 ax(cnt2) = ax(cnt2) + edif(yy+1,xx)/dx/dx - imp*(Vo(xx,yy+1)+abs(Vo(xx,yy+1)))*(-1./dx)/2.   &
			                     + (edif(yy+1,xx)-edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Vo(xx,yy)-abs(Vo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Vo(xx,yy)-abs(Vo(xx,yy)))*(Vo(xx,yy+1) - Vo(xx,yy))/dx/2./tmpV   &
							  + edif(yy,xx)*(-1.)/dx/dx
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)-abs(Vo(xx,yy)))*(O2(yy+1,xx)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)-abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy+1) - Vo(xx,yy)))/dx
		   end if
		   
         else 
		   
		   cnt3 = 0
		   
		   if (matrix(yy,xx-1)%class /=p) then 
		     
			 k = cnt_rec(yy,xx - 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,xx-1)/dx/dx - imp*(Uo(xx-1,yy)-abs(Uo(xx-1,yy)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 ax(cnt2) = ax(cnt2) - imp*(Uo(xx,yy)+abs(Uo(xx,yy)))*(1./dx)/2.  &
             			 - imp*cntgrad*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy)-Uo(xx-1,yy))/dx   &
						 + edif(yy,xx)*(-1.)/dx/dx + (edif(yy,xx) - edif(yy,xx-1))*(1.)/dx/dx
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)+abs(Uo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy,xx-1)%oxygen)/2. &
 			           + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(xx-1,yy)))/dx
		   end if
		   
		   k = cnt_rec(yy,xx)
		   ax(cnt2) = ax(cnt2) -edif(yy,xx)*(1.)/dx/dx - 1./dt - imp*(Vo(xx,yy)+abs(Vo(xx,yy)))*(1.)/dx/2. &
		               - imp*cntgrad*(Vo(xx,yy)+abs(Vo(xx,yy)))/tmpV/2.*(Vo(xx,yy)-Vo(xx,yy-1))/dx   &
					   - rx_imp*o2(yy,xx)%oxygen_use   &
					    + (edif(yy,xx) - edif(yy-1,xx))*(1.)/dx/dx
           ai(cnt2) = k 
           cnt2 = cnt2 + 1
           cnt = cnt + 1
		   
		   cnt3 = 0
		   
		   if (matrix(yy,xx + 1)%class /=p) then 
		     
			 k = cnt_rec(yy,xx + 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,xx+1)/dx/dx - imp*(Uo(xx+1,yy)+abs(Uo(xx+1,yy)))*(-1./dx)/2.   &
			            + (edif(yy,xx+1) - edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Uo(xx,yy)-abs(Uo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Uo(xx,yy)-abs(Uo(xx,yy)))*(Uo(xx+1,yy) - Uo(xx,yy))/dx/2./tmpU  &
							  + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)-abs(Uo(xx,yy)))*(O2(yy,xx+1)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(xx+1,yy) - Uo(xx,yy)))/dx
		   end if
			 
		   if (matrix(yy+1,xx)%class /=p) then 
		     
			 k = cnt_rec(yy+1,xx) 
			 ax(cnt2) = ax(cnt2) + edif(yy+1,xx)/dx/dx - imp*(Vo(xx,yy+1)+abs(Vo(xx,yy+1)))*(-1./dx)/2.   &
			                     + (edif(yy+1,xx) - edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Vo(xx,yy)-abs(Vo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Vo(xx,yy)-abs(Vo(xx,yy)))*(Vo(xx,yy+1) - Vo(xx,yy))/dx/2./tmpV   &
							  + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)-abs(Vo(xx,yy)))*(O2(yy+1,xx)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)-abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy+1) - Vo(xx,yy)))/dx
		   end if
		   
		 end if
     
       ELSE IF (yy == y_fin) THEN
         
         if (xx == 1) then
		   
		   cnt3 = 0
		   
		   if (matrix(yy - 1,xx )%class /=p) then 
		     
			 k = cnt_rec(yy - 1,xx) 
			 ax(cnt2) = ax(cnt2) + edif(yy-1,xx)/dx/dx - imp*(Vo(xx,yy-1)-abs(Vo(xx,yy-1)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 ax(cnt2) = ax(cnt2) - imp*(Vo(xx,yy)+abs(Vo(xx,yy)))*(1./dx)/2. &
			              - cntgrad*imp*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy)-Vo(xx,yy-1))/dx   &
						   + edif(yy,xx)*(-1.)/dx/dx + (edif(yy,xx) -edif(yy-1,xx))*(1.)/dx/dx
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)+abs(Vo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy-1,xx)%oxygen)/2. &
			          + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy) - Vo(xx,yy-1)))/dx
		   end if
		   
		   k = cnt_rec(yy,xx)
		   ax(cnt2) = ax(cnt2)  - 1./dt   &
					   - rx_imp*o2(yy,xx)%oxygen_use
           ai(cnt2) = k 
           cnt2 = cnt2 + 1
           cnt = cnt + 1
		   
		   cnt3 = 0
		   
		   if (matrix(yy,xx + 1)%class /=p) then 
		     
			 k = cnt_rec(yy,xx + 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,xx+1)/dx/dx - imp*(Uo(xx+1,yy)+abs(Uo(xx+1,yy)))*(-1./dx)/2.  &
			             + (edif(yy,xx+1) - edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Uo(xx,yy)-abs(Uo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Uo(xx,yy)-abs(Uo(xx,yy)))*(Uo(xx+1,yy) - Uo(xx,yy))/dx/2./tmpU   &
							  + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)-abs(Uo(xx,yy)))*(O2(yy,xx+1)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(xx+1,yy) - Uo(xx,yy)))/dx
		   end if
		   
		   if (matrix(yy,n_col)%class /=p) then 
		     
			 k = cnt_rec(yy,n_col) 
			 ax(cnt2) = ax(cnt2) + edif(yy,n_col)/dx/dx - imp*(Uo(n_col,yy)-abs(Uo(n_col,yy)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Uo(xx,yy)+abs(Uo(xx,yy)))*(1./dx)/2. &
			             - cntgrad*imp*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(n_col,yy))/dx  &
						 + edif(yy,xx)*(-1.)/dx/dx + (edif(yy,xx) - edif(yy,n_col))*(1.)/dx/dx
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)+abs(Uo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy,n_col)%oxygen)/2. &
			          + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(n_col,yy)))/dx
		   end if
		   
         else if (xx == n_col) then
		   
		   cnt3 = 0
		   cnt4 = 0
		   cnt5 = 0
		   
		   if (matrix(yy - 1,xx)%class /=p) then 
		     
			 k = cnt_rec(yy - 1,xx) 
			 ax(cnt2) = ax(cnt2) + edif(yy-1,xx)/dx/dx - imp*(Vo(xx,yy-1)-abs(Vo(xx,yy-1)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 cnt4 = 1
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)+abs(Vo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy-1,xx)%oxygen)/2. &
			       + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy) - Vo(xx,yy-1)))/dx
		   end if
		   
		   if (matrix(yy,1)%class /=p) then 
		     
			 k = cnt_rec(yy,1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,1)/dx/dx - imp*(Uo(1,yy)+abs(Uo(1,yy)))*(-1./dx)/2.   &
			              + (edif(yy,1) - edif(yy, n_col))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 cnt5 = 1
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)-abs(Uo(xx,yy)))*(O2(yy,1)%oxygen-O2(yy,xx)%oxygen)/2. &
			     + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(1,yy) - Uo(xx,yy)))/dx
		   end if
		   
		   if (matrix(yy,xx - 1)%class /=p) then 
		     
			 k = cnt_rec(yy,xx - 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,xx-1)/dx/dx - imp*(Uo(xx-1,yy)-abs(Uo(xx-1,yy)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 ax(cnt2) = ax(cnt2) - imp*(Uo(xx,yy)+abs(Uo(xx,yy)))*(1./dx)/2. &
			      - cntgrad*imp*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(xx-1,yy))/dx  &
				   + edif(yy,xx)*(-1.)/dx/dx + (edif(yy,xx) - edif(yy,xx-1))*(1.)/dx/dx
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)+abs(Uo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy,xx-1)%oxygen)/2. &
			     + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(xx-1,yy)))/dx
		   end if
		   
		   k = cnt_rec(yy,xx)
		   ax(cnt2) = ax(cnt2)  - 1./dt & 
		      - imp*cnt4*(Vo(xx,yy)+abs(Vo(xx,yy)))*(1./dx)/2. &
			  -imp*cnt4*cntgrad*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy) - Vo(xx,yy-1))/dx  &
			  - imp*cnt5*(Uo(xx,yy)-abs(Uo(xx,yy)))*(-1./dx)/2. &
			  - imp*cnt5*cntgrad*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(1,yy) - Uo(xx,yy))/dx   &
					   - rx_imp*o2(yy,xx)%oxygen_use    &
			  + cnt4*edif(yy,xx)*(-1.)/dx/dx + cnt4*(edif(yy,xx)-edif(yy-1,xx))*(1.)/dx/dx   &
			  + cnt5*edif(yy,xx)*(-1.)/dx/dx 
           ai(cnt2) = k 
           cnt2 = cnt2 + 1
           cnt = cnt + 1
		   
		   
         else 
		   
		   cnt3 = 0
		   cnt4 = 0
		   
		   if (matrix(yy - 1,xx)%class /=p) then 
		     
			 k = cnt_rec(yy - 1,xx) 
			 ax(cnt2) = ax(cnt2) + edif(yy-1,xx)/dx/dx - imp*(Vo(xx,yy-1)-abs(Vo(xx,yy-1)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 cnt4 = 1
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)+abs(Vo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy-1,xx)%oxygen)/2. &
			     + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy) - Vo(xx,yy-1)))/dx
		   end if
		   
		   if (matrix(yy,xx - 1)%class /=p) then 
		     
			 k = cnt_rec(yy,xx - 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,xx-1)/dx/dx - imp*(Uo(xx-1,yy)-abs(Uo(xx-1,yy)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 ax(cnt2) = ax(cnt2) - imp*(Uo(xx,yy)+abs(Uo(xx,yy)))*(1./dx)/2. &
			     - cntgrad*imp*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(xx-1,yy))/dx   & 
				 + edif(yy,xx)*(-1.)/dx/dx +(edif(yy,xx) - edif(yy, xx-1))*(1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)+abs(Uo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy,xx-1)%oxygen)/2. &
			     + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(xx-1,yy)))/dx
		   end if
		   
		   k = cnt_rec(yy,xx)
		   ax(cnt2) = ax(cnt2)  - 1./dt &
		     - imp*cnt4*(Vo(xx,yy)+abs(Vo(xx,yy)))*(1./dx)/2. &
			 -imp*cnt4*cntgrad*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy)-Vo(xx,yy-1))/dx   &
					   - rx_imp*o2(yy,xx)%oxygen_use   &
			 +cnt4*edif(yy,xx)*(-1.)/dx/dx + cnt4*(edif(yy,xx) - edif(yy-1,xx))*(1.)/dx/dx  
           ai(cnt2) = k 
           cnt2 = cnt2 + 1
           cnt = cnt + 1
		   
		   cnt3 = 0
			 
		   if (matrix(yy,xx + 1)%class /=p) then 
		     
			 k = cnt_rec(yy,xx + 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,xx+1)/dx/dx - imp*(Uo(xx+1,yy)+abs(Uo(xx+1,yy)))*(-1./dx)/2.  & 
			            + (edif(yy,xx+1) - edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Uo(xx,yy)-abs(Uo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Uo(xx,yy)-abs(Uo(xx,yy)))*(Uo(xx+1,yy) - Uo(xx,yy))/dx/2./tmpU   &
							  + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)-abs(Uo(xx,yy)))*(O2(yy,xx+1)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(xx+1,yy) - Uo(xx,yy)))/dx
		   end if
         
		 end if 
		 
		 
       ELSE 
         
         if (xx == 1) then
		   
		   cnt3 = 0
		   
		   if (matrix(yy - 1,1)%class /=p) then 
		     
			 k = cnt_rec(yy - 1,1) 
			 ax(cnt2) = ax(cnt2) + edif(yy-1,xx)/dx/dx - imp*(Vo(xx,yy-1) - abs(Vo(xx,yy-1)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 ax(cnt2) = ax(cnt2) - imp*(Vo(xx,yy)+abs(Vo(xx,yy)))*(1./dx)/2. &
			   - cntgrad*imp*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy) - Vo(xx,yy-1))/dx  & 
			   + edif(yy,xx)*(-1.)/dx/dx + (edif(yy,xx) - edif(yy-1,xx))*(1.0)/dx/dx
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)+abs(Vo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy-1,xx)%oxygen)/2. &
			    + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy) - Vo(xx,yy-1)))/dx
		   end if
		   
		   k = cnt_rec(yy,xx)
		   ax(cnt2) = ax(cnt2) - 1./dt   &
					   - rx_imp*o2(yy,xx)%oxygen_use
           ai(cnt2) = k 
           cnt2 = cnt2 + 1
           cnt = cnt + 1
		   
		   cnt3 = 0
		   
		   if (matrix(yy,2)%class /=p) then 
		     
			 k = cnt_rec(yy,2) 
			 ax(cnt2) = ax(cnt2) + edif(yy,2)/dx/dx - imp*(Uo(2,xx)+abs(Uo(2,yy)))*(-1./dx)/2.  &
			             + (edif(yy,2) - edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Uo(xx,yy)-abs(Uo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Uo(xx,yy)-abs(Uo(xx,yy)))*(Uo(xx+1,yy) - Uo(xx,yy))/dx/2./tmpU   &
							  + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)-abs(Uo(xx,yy)))*(O2(yy,xx+1)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(xx+1,yy) - Uo(xx,yy)))/dx
		   end if
			 
		   if (matrix(yy,n_col)%class /=p) then 
		     
			 k = cnt_rec(yy,n_col) 
			 ax(cnt2) = ax(cnt2) + edif(yy,n_col)/dx/dx - imp*(Uo(n_col,yy)-abs(Uo(n_col,yy)))*(1./dx)/2.  
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Uo(xx,yy)+abs(Uo(xx,yy)))*(1./dx)/2. &
			    - cntgrad*imp*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy)-Uo(n_col,yy))/dx   & 
				+ edif(yy,xx)*(-1.)/dx/dx + (edif(yy,xx) - edif(yy,n_col))*(1.)/dx/dx
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)+abs(Uo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy,n_col)%oxygen)/2. &
			    + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(n_col,yy)))/dx
		   end if
			 
		   if (matrix(yy + 1,xx)%class /=p) then 
		     
			 k = cnt_rec(yy + 1,xx) 
			 ax(cnt2) = ax(cnt2) + edif(yy+1,xx)/dx/dx - imp*(Vo(xx,yy+1)+abs(Vo(xx,yy+1)))*(-1./dx)/2.  &
			            + (edif(yy+1,xx) - edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Vo(xx,yy)-abs(Vo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Vo(xx,yy)-abs(Vo(xx,yy)))*(Vo(xx,yy+1) - Vo(xx,yy))/dx/2./tmpV   &
							  + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)-abs(Vo(xx,yy)))*(O2(yy+1,xx)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)-abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy+1) - Vo(xx,yy)))/dx
		   end if
		   
         else if (xx == n_col) then
		   
		   cnt3 = 0
		   cnt4 = 0
		   cnt5 = 0
		   
		   if (matrix(yy - 1,xx)%class /=p) then 
		     
			 k = cnt_rec(yy - 1,xx) 
			 ax(cnt2) = ax(cnt2) + edif(yy-1,xx)/dx/dx - imp*(Vo(xx,yy-1)-abs(Vo(xx,yy-1)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 cnt4 = 1
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)+abs(Vo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy-1,xx)%oxygen)/2. &
			    + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy) - Vo(xx,yy-1)))/dx
		   end if
		   
		   if (matrix(yy,  1)%class /=p) then 
		     
			 k = cnt_rec(yy, 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,1)/dx/dx - imp*(Uo(1,yy)+abs(Uo(1,yy)))*(-1./dx)/2.  &
			             + (edif(yy,1) - edif(yy,n_col))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 cnt5 = 1
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)-abs(Uo(xx,yy)))*(O2(yy,1)%oxygen-O2(yy,xx)%oxygen)/2. &
			    + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(1,yy) - Uo(xx,yy)))/dx
		   end if
		   
		   if (matrix(yy, n_col - 1)%class /=p) then 
		     
			 k = cnt_rec(yy, n_col - 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy, n_col-1)/dx/dx - imp*(Uo(xx-1,yy)-abs(Uo(xx-1,yy)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 ax(cnt2) = ax(cnt2) - imp*(Uo(xx,yy)+abs(Uo(xx,yy)))*(1./dx)/2. &
			    - cntgrad*imp*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy)- Uo(xx-1,yy))/dx   & 
				+ edif(yy,xx)*(-1.)/dx/dx + (edif(yy,xx) - edif(yy, n_col-1))*(1.)/dx/dx
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)+abs(Uo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy,n_col-1)%oxygen)/2. &
			    + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(n_col-1,yy)))/dx
		   end if
		   
		   k = cnt_rec(yy,xx)
		   ax(cnt2) = ax(cnt2)  - 1./dt  &
		     - imp*cnt4*(Vo(xx,yy)+abs(Vo(xx,yy)))*(1./dx)/2.  & 
			 - imp*cnt4*cntgrad*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy)-Vo(xx,yy-1))/dx  &
		     - imp*cnt5*(Uo(xx,yy)-abs(Uo(xx,yy)))*(-1./dx)/2.   &
			 - imp*cnt5*cntgrad*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(1,yy)-Uo(xx,yy))/dx   &
					   - rx_imp*o2(yy,xx)%oxygen_use   &
			 + cnt4*edif(yy,xx)*(-1.)/dx/dx + cnt4*(edif(yy,xx)- edif(yy-1,xx))*(1.)/dx/dx   &
			 + cnt5*edif(yy,xx)*(-1.)/dx/dx
           ai(cnt2) = k 
           cnt2 = cnt2 + 1
           cnt = cnt + 1
		   
		   cnt3 = 0
			 
		   if (matrix(yy+1,xx)%class /=p) then 
		     
			 k = cnt_rec(yy+1,xx) 
			 ax(cnt2) = ax(cnt2) + edif(yy+1,xx)/dx/dx - imp*(Vo(xx,yy+1)+abs(Vo(xx,yy+1)))*(-1./dx)/2.   &
			                 + (edif(yy+1,xx) - edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
		     ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Vo(xx,yy)-abs(Vo(xx,yy)))*(-1./dx)/2. & 
			                  - imp*cntgrad*(Vo(xx,yy)-abs(Vo(xx,yy)))*(Vo(xx,yy+1) - Vo(xx,yy))/dx/2./tmpV  &
							  + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)-abs(Vo(xx,yy)))*(O2(yy+1,xx)%oxygen-O2(yy,xx)%oxygen)/2. &
			   + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)-abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy+1) - Vo(xx,yy)))/dx
		   end if
		   
         else 
		   
		   cnt3 = 0
		   cnt4 = 0
		   
		   if (matrix(yy - 1,xx)%class /=p) then 
		     
			 k = cnt_rec(yy - 1,xx) 
			 ax(cnt2) = ax(cnt2) + edif(yy-1,xx)/dx/dx - imp*(Vo(xx,yy-1)-abs(Vo(xx,yy-1)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 cnt4 = 1
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)+abs(Vo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy-1,xx)%oxygen)/2. &
			    + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy) - Vo(xx,yy-1)))/dx
		   end if
		   
		   if (matrix(yy ,xx - 1)%class /=p) then 
		     
			 k = cnt_rec(yy ,xx - 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,xx-1)/dx/dx - imp*(Uo(xx-1,yy)-abs(Uo(xx-1,yy)))*(1./dx)/2.
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 ax(cnt2) = ax(cnt2) - imp*(Uo(xx,yy)+abs(Uo(xx,yy)))*(1./dx)/2. &
			     - cntgrad*imp*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy)-Uo(xx-1,yy))/dx   &
				 + edif(yy,xx)*(-1.)/dx/dx + (edif(yy,xx) - edif(yy,xx-1))*(1.)/dx/dx
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)+abs(Uo(xx,yy)))*(O2(yy,xx)%oxygen-O2(yy,xx-1)%oxygen)/2. &
			    + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)+abs(Uo(xx,yy)))/2./tmpU*(Uo(xx,yy) - Uo(xx-1,yy)))/dx
		   end if
		   
		   k = cnt_rec(yy,xx)
		   ax(cnt2) = ax(cnt2) - 1./dt &
		      - imp*cnt4*(Vo(xx,yy)+abs(Vo(xx,yy)))*(1./dx)/2. &
			  - imp*cnt4*cntgrad*(Vo(xx,yy)+abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy)-Vo(xx,yy-1))/dx   &
					   - rx_imp*o2(yy,xx)%oxygen_use   &
					   + cnt4*edif(yy,xx)*(-1.)/dx/dx  + cnt4*(edif(yy,xx) - edif(yy-1,xx))*(1.)/dx/dx
           ai(cnt2) = k 
           cnt2 = cnt2 + 1
           cnt = cnt + 1
		   
		   cnt3 = 0
		   
		   if (matrix(yy,xx + 1)%class /=p) then 
		     
			 k = cnt_rec(yy,xx + 1) 
			 ax(cnt2) = ax(cnt2) + edif(yy,xx+1)/dx/dx - imp*(Uo(xx+1,yy)+abs(Uo(xx+1,yy)))*(-1./dx)/2.  &
			     + (edif(yy,xx+1) - edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
			 ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Uo(xx,yy)-abs(Uo(xx,yy)))*(-1./dx)/2. &
			     - cntgrad*imp*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(xx+1,yy)-Uo(xx,yy))/dx   & 
				 + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Uo(xx,yy)-abs(Uo(xx,yy)))*(O2(yy,xx+1)%oxygen-O2(yy,xx)%oxygen)/2. &
			    + cntgrad*O2(yy,xx)%oxygen*(Uo(xx,yy)-abs(Uo(xx,yy)))/2./tmpU*(Uo(xx+1,yy) - Uo(xx,yy)))/dx
		   end if
			 
		   if (matrix(yy+1,xx)%class /=p) then 
		     
			 k = cnt_rec(yy+1,xx) 
			 ax(cnt2) = ax(cnt2) + edif(yy+1,xx)/dx/dx - imp*(Vo(xx,yy+1)+abs(Vo(xx,yy+1)))*(-1./dx)/2.   &
			                   + (edif(yy+1,xx) - edif(yy,xx))*(-1.)/dx/dx
             ai(cnt2) = k 
			 cnt3 = cnt3 + 1
             cnt2 = cnt2 + 1
             cnt = cnt + 1
			 
			 ax(cnt2 - 1 - cnt3) = ax(cnt2 - 1 - cnt3) - imp*(Vo(xx,yy)-abs(Vo(xx,yy)))*(-1./dx)/2. &
			     - cntgrad*imp*(Vo(xx,yy)-abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy+1)-Vo(xx,yy))/dx    & 
				 + edif(yy,xx)*(-1.)/dx/dx 
			 
			 b(j) = b(j) + epl*((Vo(xx,yy)-abs(Vo(xx,yy)))*(O2(yy+1,xx)%oxygen-O2(yy,xx)%oxygen)/2. &
			    + cntgrad*O2(yy,xx)%oxygen*(Vo(xx,yy)-abs(Vo(xx,yy)))/2./tmpV*(Vo(xx,yy+1) - Vo(xx,yy)))/dx
		   end if
		 
		 end if
		 
	  END IF
       
	  ap(j + 1) =   ap(j) + cnt 
	  
	end do 
  
  end do
	   
   
   ai = ai - 1
   
   mtx_chk = .FALSE.
   ! mtx_chk = .TRUE.
   
   if (mtx_chk) then
   
   open(500, file="chck_ap.txt", status = 'unknown')
   do j = 1, n+1
   write(500,*) ap(j)
   end do 
   close(500)
   
   open(500, file="chck_ai.txt", status = 'unknown')
   do j = 1, nnz
   write(500,*) ai(j)
   end do 
   close(500)
   
   open(500, file="chck_ax.txt", status = 'unknown')
   do j = 1, nnz
   write(500,*) ax(j)
   end do 
   close(500)
   
   open(500, file="chck_b.txt", status = 'unknown')
   do j = 1, n
   write(500,*) b(j)
   end do 
   close(500)
   
   end if
   
   if (mtx_chk) stop
   !
!  Set the default control parameters.
!
  call umf4def( control )
!
!  From the matrix data, create the symbolic factorization information.
!
  call umf4sym ( n, n, ap, ai, ax, symbolic, control, info )

  if ( info(1) < 0.0D+00 ) then
    write ( *, * ) ''
    write ( *, *) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, * ) '  UMF4SYM returns INFO(1) = ', info(1)
    stop 1
  end if
!
!  From the symbolic factorization information, carry out the numeric factorization.
!
  call umf4num ( ap, ai, ax, symbolic, numeric, control, info )

  if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4NUM returns INFO(1) = ', info(1)
    stop 1
  end if
!
!  Free the memory associated with the symbolic factorization.
!
  call umf4fsym ( symbolic )
!
!  Solve the linear system.
!
  sys = 0
  call umf4sol ( sys, x, b, numeric, control, info )
  
  if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4SOL returns INFO(1) = ', info(1)
    stop 1
  end if
!
!  Free the memory associated with the numeric factorization.
!
  call umf4fnum ( numeric )
!
!  Print the solution.
!
  
  write ( *, * ) ''
  write ( *, * ) '  Computed solution:', TIme
  write ( *, * ) ''
  
  call O2pre_setup()
  
  do yy = y_int  + 1, y_fin
   do xx = 1, n_col
     if (matrix(yy,xx)%class == p ) cycle
	 
     j = cnt_rec(yy,xx)
     O2(yy,xx)%oxygen = x(j)
	   
	 if (O2(yy,xx)%oxygen < 0) then 
	   if (O2(yy,xx)%value_pre /=1) then 
	     O2(yy,xx)%oxygen = 0.0
	   else 
		 print *, "================="
		 print *, "=negative oxygen=" , O2(yy,xx)%oxygen, xx,yy, O2(yy,xx)%value_pre
		 print *, "================="
	     O2(yy,xx)%oxygen = 0.0
	   end if
     end if 

     if (O2(yy,xx)%oxygen > 1.) then
       	 
	   print *, "%%%%%%%%%%%%%%%%%%%%"
	   print *, "=oxygen production?=" , O2(yy,xx)%oxygen, xx,yy, O2(yy,xx)%value_pre
	   print *, "%%%%%%%%%%%%%%%%%%%%"
		 
      end if 
	  
	if (ieee_is_nan(O2(yy,xx)%oxygen)) then
	
	  print *,"//////////////////"
	  print *,"//  NAN at",yy,xx,"///"
	  print *,"//////////////////"
	  
	  if (matrix(yy,xx)%class == p) O2(yy,xx)%oxygen = 0.0
	  
	  
	  
	end if
		 
   end do
  end do 
  
  end do
  
   end subroutine oxygen_profile

   ! *********************************************************
   
   end module
 