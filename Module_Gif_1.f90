
   module bin_io

   !  Buffered reading and writing data to binary files. This module can use
   !  the Fortran file-units MinUnit .. MaxUnit for this purpose.
   !  Version 1.0, February 1998
   !  Written by Jos Bergervoet

   implicit none                  ! Check all declarations

   public  :: open_for_read, open_for_write, close_block_io,  &
              integer_block_read, integer_block_write,        &
              char_block_write, char_read

   private :: init_file_block, read_last_rec, flush_last_rec

   integer, private, parameter    :: MinUnit = 90, MaxUnit = 99, BufLen = 256

   type, public :: FileHandle     ! Can be used outside this unit (is public)
     private                      ! but all fields are unknown outside (private)
     character(len=BufLen), pointer  :: Buf
     character(len=255)              :: Name
     integer                         :: FilePos, BufPos
     logical                         :: ForWriting
   end type FileHandle

   type(FileHandle), private, dimension(MinUnit:MaxUnit), save  :: Info

   contains

   subroutine init_file_block(Funit, Fname, OpenForWrite)
     integer, intent(in)           :: Funit
     character(len=*), intent(in)  :: Fname
     logical, intent(in)           :: OpenForWrite

     if (Funit<MinUnit .or. Funit>MaxUnit) then
       write(unit=*,fmt=*) "Error: binary file-unit",Funit, "out of range."
       stop
     end if
     allocate ( Info(Funit)%Buf )
     Info(Funit) % FilePos = 0
     Info(Funit) % Name = Fname
     Info(Funit) % BufPos = 0
     Info(Funit) % ForWriting = OpenForWrite
     return
   end subroutine init_file_block

   subroutine open_for_read(Funit, Fname)
     integer, intent(in)           :: Funit
     character(len=*), intent(in)  :: Fname

     call init_file_block(Funit, Fname, OpenForWrite=.false.)
     open (unit=Funit, file=Fname, form="unformatted", access="direct", &
           recl=BufLen, status="old", action="read")
     return
   end subroutine open_for_read

   subroutine open_for_write(Funit, Fname)
     integer, intent(in)           :: Funit
     character(len=*), intent(in)  :: Fname

     call init_file_block(Funit, Fname, OpenForWrite=.true.)
     open (unit=Funit, file=Fname, form="unformatted", access="direct", &
           recl=BufLen, status="replace", action="write")
     return
   end subroutine open_for_write

   subroutine flush_last_rec(Funit)           ! Write (last) record < BufLen
     integer, intent(in)           :: Funit
     integer                       :: fp, i
                                              ! close, then re-open with recl=1
     close (unit=Funit)
     open (unit=Funit,file=Info(Funit)%Name,form="unformatted",access="direct", &
           recl=1, status="old", action="write")

     fp = Info(Funit) % FilePos * BufLen
     do i = 1,Info(Funit) % BufPos
       write(unit=Funit,rec=fp+i) Info(Funit) % buf(i:i)
     end do
     return
   end subroutine flush_last_rec

   subroutine close_block_io(Funit)
     integer, intent(in)           :: Funit
     if (Info(Funit)%ForWriting .and. Info(Funit)%BufPos > 0) then
       call flush_last_rec(Funit)
     end if
     close(unit=Funit)
     deallocate ( Info(Funit)%Buf )
     return
   end subroutine close_block_io

   subroutine char_block_write(Funit, s)
     integer, intent(in)           :: Funit
     character(len=*), intent(in)  :: s
     integer                       :: ls, ps, fp, bp, slice

     fp = Info(Funit) % FilePos
     bp = Info(Funit) % BufPos
     ls = len(s)
     ps = 0
     do  ! write the string in slices of our buffersize BufLen
       if (ls <= ps) then
         exit
       else if (ls-ps > BufLen-bp) then
         slice = BufLen-bp
       else
         slice = ls-ps
       end if
       Info(Funit) % buf(bp+1 : bp+slice) = s(ps+1 : ps+slice)
       bp = bp+slice
       ps = ps+slice
       if (bp == BufLen) then
         write(unit=Funit,rec=fp+1) Info(Funit) % buf
         fp = fp+1
         bp = 0
       end if
     end do
     Info(Funit) % FilePos = fp
     Info(Funit) % BufPos = bp
     return
   end subroutine char_block_write

   subroutine read_last_rec(Funit, NewPtr, ErrCod)  ! Last record < BufLen
     integer, intent(in)  :: Funit
     integer, intent(out) :: NewPtr, ErrCod
     integer              :: fp, bp, Nfound

     ! We'll close and re-open with recl=1 to read the last few bytes one by one.
     close (unit=Funit)
     open (unit=Funit,file=Info(Funit)%Name,form="unformatted",access="direct", &
           recl=1, status="old", action="read")
     fp = Info(Funit) % FilePos * BufLen

     Nfound = 0
     do bp = 1,BufLen
       read(unit=Funit,rec=fp+bp,iostat=ErrCod) Info(Funit) % buf(bp:bp)
       if (ErrCod /= 0) then
         Nfound = bp-1
         exit
       end if
       Nfound = bp
     end do
     if (Nfound > 0) then      ! at least we found some valid bytes
       ErrCod = 0
     end if
                               ! Now we shift the bytes to high side of buffer
     NewPtr = BufLen-Nfound+1
     Info(Funit)%Buf (NewPtr:BufLen) = Info(Funit)%Buf (1:Nfound)
     return
   end subroutine read_last_rec

   subroutine char_read(Funit, c, ReadErr)
     integer, intent(in)            :: Funit
     character(len=*), intent(out)  :: c        ! only c(1:1) will be written
     integer, intent(out), optional :: ReadErr
     integer                        :: fp, bp, ErrCod

     ErrCod = 0

     bp = Info(Funit) % BufPos
     if (bp == 0) then
       fp = Info(Funit) % FilePos
       read(unit=Funit,rec=fp+1,iostat=ErrCod) Info(Funit) % buf
       if (ErrCod == 0) then
         bp = 1
       else
         call read_last_rec(Funit, bp, ErrCod)
       end if
       Info(Funit) % FilePos = fp+1
     end if

     if (present(ReadErr)) then
       ReadErr = ErrCod
       if (ErrCod /= 0) then
         return
       end if
     else if (ErrCod /= 0) then
       write(unit=*,fmt=*) "Error reading binary file ", trim(Info(Funit)%Name)
       stop
     end if

     c = Info(Funit) % buf(bp:bp)
     Info(Funit) % BufPos = modulo(bp+1, BufLen+1)
     return
   end subroutine char_read

   subroutine integer_block_read(Funit, iread, length)
     integer, intent(in)   :: Funit, length
     integer, intent(out)  :: iread
     integer               :: i
     character(len=1)      :: c

     iread = 0
     do i=1,length
       call char_read(Funit, c)
       iread = iread + ichar(c)*256**(i-1)
     end do
     return
   end subroutine integer_block_read

   subroutine integer_block_write(Funit, iwrite, length)
     integer, intent(in)   :: Funit, length, iwrite
     integer               :: i
     character(len=length) :: c

     do i=1,length
       c(i:i) = char(modulo( iwrite/256**(i-1), 256))
     end do
     call char_block_write(Funit, c)
     return
   end subroutine integer_block_write

   end module bin_io
