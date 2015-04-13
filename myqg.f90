
program main
  use myqg_module
  implicit none

  integer :: i,j,k,l,n
  integer :: is,ie,js,je
  character(len=128) :: fname
  character(len=128) :: fprfx
  !character(len=128) :: path_data
  integer, dimension(3) :: dims
  real*8, allocatable, dimension(:,:,:) :: A
  integer :: nt
  logical :: doio
  real*8 :: tmpreal

  call initialize_setup
  write(*,*) "Finished: Initialize variables"

  ! write grid data 
  call write_3d("XC        ", xt, (/ nx,  1,  1 /), -1)
  call write_3d("YC        ", yt, (/  1, ny,  1 /), -1)
  call write_3d("RC        ", zt, (/  1,  1, nz /), -1)
  call write_3d("XG        ", xu, (/ nx,  1,  1 /), -1)
  call write_3d("YG        ", yu, (/  1, ny,  1 /), -1)
  !call write_3d("RF        ", zu, (/  1,  1, nz /), -1)
  call write_3d("RF        ", (/ zu(:), 0.d0 /), (/  1,  1, nz+1 /), -1)

  !call cyclic_exchange(pv)
  !call cyclic_exchange(forc)
  !call write_3d("psi       ", psi(1:nx,1:ny,1:nz),  (/ nx,  ny,  nz /), 0)
  !call write_3d("pv        ", pv (1:nx,1:ny,1:nz),  (/ nx,  ny,  nz /), 0)
  !call write_3d("u         ", u  (1:nx,1:ny,1:nz),  (/ nx,  ny,  nz /), 0)
  !call write_3d("v         ", v  (1:nx,1:ny,1:nz),  (/ nx,  ny,  nz /), 0)
  !call write_3d("hpr       ", hpr(1:nx,1:ny,1:nz),  (/ nx,  ny,  nz /), 0)
  !call write_3d("wek       ", wek(1:nx,1:ny),       (/ nx,  ny,   1 /), 0)
  !call write_3d("pvr       ", pvr(1:nx,1:ny,1:nz),  (/ nx,  ny,  nz /), 0)
  !call write_3d("pvp       ", pvp(1:nx,1:ny,1:nz),  (/ nx,  ny,  nz /), 0)
  call write_snapshot

! ================================================================================ 
! main time stepping loop
! ================================================================================ 
  write(*,*) "Start: Time stepping"
  do while ( time < t_end )
    tstep = tstep + 1
    time  = tstep*dt

    if ( floor(time/time_monitor) == time/time_monitor ) then
      doio = .true.
    else
      doio = .false.
    endif

    if ( doio ) then
      write(*,*) "================================================================================"
      write(*,*) " tstep = ", tstep
      write(*,*) "================================================================================"
    endif
 
  ! assume psi is known already (initial psi or from last iteration) 

    ! calculate interface dissplacement
    do j=1-ox,ny+ox
      do i=1-ox,nx+ox
        do k=1,nz
          if ( k==1 ) then
            hpr(i,j,k) = f0/gred(k)*( psi(i,j,k) )
          else
            hpr(i,j,k) = f0/gred(k)*( psi(i,j,k) - psi(i,j,k-1) )
          endif
        enddo
      enddo
    enddo

    ! calculate layer width and reciprocal volume
    do j=1-ox,ny+ox
      do i=1-ox,nx+ox
        do k=1,nz
          dz(i,j,k)       = Hk(k) + hpr(i,k,k)
          recepvol(i,j,k) = 1.0 / (dx*dy*dz(i,j,k))
        enddo
      enddo
    enddo

    ! inhomogenity of poisson equation
    pvr = pv - pvp
    forc = pvr
    !do j=1-ox,ny+ox
    !  do i=1-ox,nx+ox
    !    do k=1,nz
    !      forc(i,j,k) = pv(i,j,k) - beta * yu(j) 
    !      !forc(i,j,k) = pv(i,j,k) - fCoru(j)
    !    enddo
    !  enddo
    !enddo

    ! calculate streamfunction
    call absmax_element(1-ox,nx+ox,1-ox,ny+ox,nz,forc,tmpreal)
    !write(*,*) "tmpreal = ", tmpreal 
    !write(*,*) "forc = ", forc(32,:,1)
    !if ( tstep == 2 ) then
    !  stop
    !endif
    call solve_poisson_cg(1-ox, nx+ox, 1-ox, ny+ox, nz, forc, psi, max_itt, crit, est_error, doio)
 
    ! calculate velocity
    call calc_curl_psi
    !call cyclic_exchange(u)
    !call cyclic_exchange(v)

    ! calculate tendencies
    call calc_Gpv

    ! pv forward step
    if ( tstep==1 ) then
      Gpvm1 = Gpv
    endif
    pv(1:nx,1:ny,1:nz) = pv(1:nx,1:ny,1:nz) + ( (1.5+epsab)*Gpv - (0.5+epsab)*Gpvm1 ) * dt

    ! reset tendencies
    Gpvm1 = Gpv   ! save old tendency for Adams-Bashforth
    Gpv   = 0.0   ! reset new tendency
   
    ! do cyclic_exchange (Does this occur at correct place???)
    !call cyclic_exchange(psi)
    !call cyclic_exchange(pv)

    ! do model I/O
    if ( floor(time/timeio) == time/timeio ) then
      call write_snapshot
      !write(*,*) "================================================================================"
      !write(*,*) "Model I/O at tstep ", tstep
      !call write_3d('pv        ', pv (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
      !call write_3d('pvp       ', pvp(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
      !call write_3d('pvr       ', pvr(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
      !call write_3d('psi       ', psi(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
      !call write_3d('u         ', u  (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
      !call write_3d('v         ', v  (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
      !call write_3d('hpr       ', hpr(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
      !write(*,*) "================================================================================"
    endif

  enddo 
! ================================================================================ 
! end main time stepping loop
! ================================================================================ 

  write(*,*) "nx = ", nx, " ny = ", ny, " nz = ", nz
  write(*,*) "All done"
end program main

subroutine calc_Gpv
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  !Gpv = -pv/t_end
  call calc_pvadv
  call calc_pvdiff
  call calc_ekman_pumping
end subroutine calc_Gpv

subroutine calc_curl_psi
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  do i=1-ox,nx+ox
    do j=1-ox,ny+ox-1
      do k=1,nz
        u(i,j,k) = - (psi(i,j+1,k)-psi(i,j,k)) / dy
      enddo
    enddo
  enddo
  do i=1-ox,nx+ox-1
    do j=1-ox,ny+ox
      do k=1,nz
        v(i,j,k) =   (psi(i+1,j,k)-psi(i,j,k)) / dx
      enddo
    enddo
  enddo
end subroutine calc_curl_psi

subroutine calc_pvdiff
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  ! prepare outer domain pv values for no-flux boundary conditions
  pvr(0,:,:)    = pvr(1,:,:)
  pvr(nx+1,:,:) = pvr(nx,:,:)
  pvr(:,0,:)    = pvr(:,1,:)
  pvr(:,ny+1,:) = pvr(:,ny,:)

  Gpvdif = 0.0
  do i=1,nx
    do j=1,ny
      do k=1,nz
        Gpvdif(i,j,k) =   1.0/dx**2 * ( pvr(i+1,j,k) - 2.0*pvr(i,j,k) + pvr(i-1,j,k) ) & 
                      & + 1.0/dy**2 * ( pvr(i,j+1,k) - 2.0*pvr(i,j,k) + pvr(i,j-1,k) )  
      enddo
    enddo
  enddo
  Gpv = Gpv + Gpvdif
end subroutine calc_pvdiff
  
subroutine calc_pvdiff_old
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox,nz) :: fZon
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox,nz) :: fMer

  fZon = 0.0
  fMer = 0.0

  do i=1-ox,nx+ox-1
    do j=1-ox,ny+ox
      do k=1,nz
        fZon(i,j,k) = -diffPVh * (pvr(i+1,j,k)-pvr(i,j,k))/dx * dy*dz(i,j,k)
      enddo
    enddo
  enddo

  if ( .not. cyclic_x ) then
    fZon(0, :,:) = 0.0
    fZon(nx,:,:) = 0.0
  endif

  do i=1-ox,nx+ox
    do j=1-ox,ny+ox-1
      do k=1,nz
        fMer(i,j,k) = -diffPVh * (pvr(i,j+1,k)-pvr(i,j,k))/dy * dx*dz(i,j,k)
      enddo
    enddo
  enddo

  if ( .not. cyclic_y ) then
    fMer(:,0, :) = 0.0
    fMer(:,ny,:) = 0.0
  endif

  ! diffusive tendency
  Gpvdif = 0.0
  do i=1,nx
    do j=1,ny
      do k=1,nz
        Gpvdif(i,j,k)  = Gpvdif(i,j,k) -  recepvol(i,j,k) * ( &
                    & fZon(i,j,k) - fZon(i-1,j,k) + &
                    & fMer(i,j,k) - fMer(i,j-1,k)   &
                    & )
      enddo
    enddo
  enddo
  Gpv = Gpv + Gpvdif
end subroutine calc_pvdiff_old

subroutine calc_pvadv
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  real*8  :: J1, J2, J3
  Gpvadv = 0.0
  do j=1,ny
    do i=1,nx
      do k=1,nz
        J1 = 1.0/(4.0*dx*dy) * ( &
           &   (psi(i+1,j,k)-psi(i-1,j,k)) * (pv (i,j+1,k)-pv (i,j-1,k)) &
           & - (psi(i,j+1,k)-psi(i,j-1,k)) * (pv (i+1,j,k)-pv (i-1,j,k)) &
           & )
        J2 = 1.0/(4.0*dx*dy) * ( &
           &   psi(i+1,j,k)*(pv (i+1,j+1,k)-pv (i+1,j-1,k)) - psi(i-1,j,k)*(pv (i-1,j+1,k)-pv (i-1,j-1,k) ) &
           & - psi(i,j+1,k)*(pv (i+1,j+1,k)-pv (i-1,j+1,k)) - psi(i,j-1,k)*(pv (i+1,j-1,k)-pv (i-1,j-1,k) ) &
           & )
        J3 = 1.0/(4.0*dx*dy) * ( &
           &   pv (i+1,j,k)*(psi(i+1,j+1,k)-psi(i+1,j-1,k)) - pv (i-1,j,k)*(psi(i-1,j+1,k)-psi(i-1,j-1,k) ) &
           & - pv (i,j+1,k)*(psi(i+1,j+1,k)-psi(i-1,j+1,k)) - pv (i,j-1,k)*(psi(i+1,j-1,k)-psi(i-1,j-1,k) ) &
           & )
        Gpvadv(i,j,k) = Gpvadv(i,j,k) + (J1 + J2 + J3) / 3.0
      enddo
    enddo
  enddo
  Gpv = Gpv + Gpvadv
end subroutine calc_pvadv

subroutine calc_pvadv_old
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox,nz) :: uTrans
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox,nz) :: vTrans
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox,nz) :: fZon
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox,nz) :: fMer

  uTrans = 0.0
  vTrans = 0.0
  fZon   = 0.0
  fMer   = 0.0

  ! zonal advective flux
  do i=1-ox,nx+ox
    do j=1-ox,ny+ox
      do k=1,nz
        uTrans(i,j,k) = 0.25 * ( u(i,j,k) + u(i+1,j,k) + u(i,j+1,k) + u(i+1,j+1,k) )*dy*dz(i,j,k)
      enddo
    enddo
  enddo

  do i=1-ox,nx+ox
    do j=1-ox,ny+ox
      do k=1,nz
        fZon(i,j,k) = uTrans(i,j,k) * 0.5 * ( pv(i,j,k)+pv(i+1,j,k) )
      enddo
    enddo
  enddo

  ! meridional advective flux
  do i=1-ox,nx+ox
    do j=1-ox,ny+ox
      do k=1,nz
        vTrans(i,j,k) = 0.25 * ( v(i,j,k) + v(i+1,j,k) + v(i,j+1,k) + v(i+1,j+1,k) )*dx*dz(i,j,k)
      enddo
    enddo
  enddo

  do i=1-ox,nx+ox
    do j=1-ox,ny+ox
      do k=1,nz
        fMer(i,j,k) = vTrans(i,j,k) * 0.5 * ( pv(i,j,k)+pv(i,j+1,k) )
      enddo
    enddo
  enddo

  ! advective tendency
  do i=1,nx
    do j=1,ny
      do k=1,nz
        Gpv(i,j,k)  = Gpv(i,j,k) -  recepvol(i,j,k) * ( &
                    & fZon(i,j,k) - fZon(i-1,j,k) + &
                    & fMer(i,j,k) - fMer(i,j-1,k)   &
                    & )
      enddo
    enddo
  enddo
end subroutine calc_pvadv_old

subroutine calc_ekman_pumping
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  Gpvfor = 0.0 
  do i=1,nx
    do j=1,ny
      Gpvfor(i,j,1)  = Gpvfor(i,j,1) + ( &
                  & f0/Hk(1) * wek(i,j) & 
                  & )
    enddo
  enddo
  Gpv = Gpv + Gpvfor
end subroutine calc_ekman_pumping


subroutine cyclic_exchange(var)
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox,nz), intent(inout) :: var

  call myqg_error("cyclic_exchange should not been called!")

  if ( cyclic_x .and. cyclic_y ) then
    do i=1,ox
      do j=1,ny
        var(1-i,j,:)   = var(nx+1-i,j,:)
        var(nx+i,j,:)  = var(i,j,:)
      enddo
    enddo
    do i=1,nx
      do j=1,ox
        var(i,1-j,:)   = var(i,ny+1-j,:)
        var(i,ny+j,:)  = var(i,j,:)
      enddo
    enddo
    do i=1,ox
      do j=1,ox
        var(1-i, 1-j, :) = var(nx+1-i,ny+1-j,:)
        var(1-i, ny+j,:) = var(nx+1-i,j,:)
        var(nx+i,1-j, :) = var(i,ny+1-j,:)
        var(nx+i,ny+j,:) = var(i,j,:)
      enddo
    enddo
  else
    call myqg_error("Both, cyclic_x and cyclic_y, have to be enabled so far")
  endif
end subroutine cyclic_exchange

subroutine cyclic_exchange_2d(var)
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox), intent(inout) :: var
  if ( cyclic_x .and. cyclic_y ) then
    do i=1,ox
      do j=1,ny
        var(1-i,j)   = var(nx+1-i,j)
        var(nx+i,j)  = var(i,j)
      enddo
    enddo
    do i=1,nx
      do j=1,ox
        var(i,1-j)   = var(i,ny+1-j)
        var(i,ny+j)  = var(i,j)
      enddo
    enddo
    do i=1,ox
      do j=1,ox
        var(1-i, 1-j)  = var(nx+1-i,ny+1-j)
        var(1-i, ny+j) = var(nx+1-i,j)
        var(nx+i,1-j)  = var(i,ny+1-j)
        var(nx+i,ny+j) = var(i,j)
      enddo
    enddo
  else
    call myqg_error("Both, cyclic_x and cyclic_y, have to be enabled so far")
  endif
end subroutine cyclic_exchange_2d

subroutine write_3d(fprfx, outdata, dims, tstepout)
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  character(len=10),            intent(in)  :: fprfx
  integer, dimension(3),        intent(in)  :: dims
  real*8, dimension(dims(1),dims(2),dims(3)),  intent(in)  :: outdata
  integer,                      intent(in)  :: tstepout
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  character(len=128)    :: fname
  character(len=10)     :: tstepstr 

  !if (.not. present(path_name)) then
  !  path_name = "./"
  !endif

! write .data file
  !if (present(tstepout)) then
  if (tstepout>=0) then
    write(tstepstr,"(I10.10)") tstepout
    fname = trim(path_data) // trim(fprfx) // "." // tstepstr // ".data"
  else
    fname = trim(path_data) // trim(fprfx) // ".data"
  endif
  open( unit=fid, file=fname, form='unformatted', status='replace', &
      & access='direct', recl=bytes*product(dims), convert=endian )
  write(fid, rec=1) sngl(outdata)
  close(fid)
  write(*,*) "Write: ", fname

! write .meta file
  !if (present(tstepout)) then
  if (tstepout>=0) then
    write(tstepstr,"(I10.10)") tstepout
    fname = trim(path_data) // trim(fprfx) // "." // tstepstr // ".meta"
  else
    fname = trim(path_data) // trim(fprfx) // ".meta"
  endif
  open( unit=fid, file=fname, form='formatted', status='replace', &
      & access='sequential')
  write(fid,*) "nDims = [", 3, "];"
  write(fid,*) "dimList = [" 
  write(fid,*) dims(1),",", 1, ",", dims(1), ","
  write(fid,*) dims(2),",", 1, ",", dims(2), ","
  write(fid,*) dims(3),",", 1, ",", dims(3)
  write(fid,*) "];"
  if ( bytes==4) then
    write(fid,*) "dataprec = [ 'float32' ];"
  elseif (bytes==8) then
    write(fid,*) "dataprec = [ 'float64' ];"
  endif
  write(fid,*) "nrecords = [ 1 ];"
  write(fid,*) "timeStepNumber = [", tstepout, " ];"
  write(fid,*) "nFlds = [    1 ];"
  write(fid,*) "fldList = {"
  write(fid,*) "'", fprfx, "'"
  write(fid,*) "};"
  close(fid)
end subroutine write_3d

! ---
subroutine write_2d(fprfx, outdata, dims, tstepout)
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  character(len=10),            intent(in)  :: fprfx
  integer, dimension(3),        intent(in)  :: dims
  !!
  real*8, dimension(dims(1),dims(2)),  intent(in)  :: outdata
  integer,                      intent(in)  :: tstepout
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  character(len=128)    :: fname
  character(len=10)     :: tstepstr 

  !if (.not. present(path_name)) then
  !  path_name = "./"
  !endif

! write .data file
  !if (present(tstepout)) then
  if (tstepout>=0) then
    write(tstepstr,"(I10.10)") tstepout
    fname = trim(path_data) // trim(fprfx) // "." // tstepstr // ".data"
  else
    fname = trim(path_data) // trim(fprfx) // ".data"
  endif
  open( unit=fid, file=fname, form='unformatted', status='replace', &
      & access='direct', recl=bytes*product(dims), convert=endian )
  write(fid, rec=1) sngl(outdata)
  close(fid)
  write(*,*) "Write: ", fname

! write .meta file
  !if (present(tstepout)) then
  if (tstepout>=0) then
    write(tstepstr,"(I10.10)") tstepout
    fname = trim(path_data) // trim(fprfx) // "." // tstepstr // ".meta"
  else
    fname = trim(path_data) // trim(fprfx) // ".meta"
  endif
  open( unit=fid, file=fname, form='formatted', status='replace', &
      & access='sequential')
  write(fid,*) "nDims = [", 3, "];"
  write(fid,*) "dimList = [" 
  write(fid,*) dims(1),",", 1, ",", dims(1), ","
  write(fid,*) dims(2),",", 1, ",", dims(2), ","
  !!write(fid,*) dims(3),",", 1, ",", dims(3)
  write(fid,*) "];"
  if ( bytes==4) then
    write(fid,*) "dataprec = [ 'float32' ];"
  elseif (bytes==8) then
    write(fid,*) "dataprec = [ 'float64' ];"
  endif
  write(fid,*) "nrecords = [ 1 ];"
  write(fid,*) "timeStepNumber = [", tstepout, " ];"
  write(fid,*) "nFlds = [    1 ];"
  write(fid,*) "fldList = {"
  write(fid,*) "'", fprfx, "'"
  write(fid,*) "};"
  close(fid)
end subroutine write_2d
!---

subroutine myqg_error(errormsg)
  implicit none
  character(len=*), intent(in) :: errormsg
  write(*,*) "::: Error: " // trim(errormsg) // " :::" 
  stop 1
end subroutine myqg_error

subroutine write_snapshot
  use myqg_module
  implicit none
  write(*,*) "================================================================================"
  write(*,*) "Model I/O at tstep ", tstep
  call write_3d('pv        ', pv (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('pvp       ', pvp(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('pvr       ', pvr(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('psi       ', psi(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('u         ', u  (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('v         ', v  (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('hpr       ', hpr(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('Gpvadv    ', Gpvadv(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('Gpvdif    ', Gpvdif(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('Gpvfor    ', Gpvfor(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  write(*,*) "================================================================================"
end subroutine write_snapshot
