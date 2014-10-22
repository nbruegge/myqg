
program main
  use myqg_module
  implicit none

  integer :: i,j,k,l,n
  character(len=128) :: fname
  character(len=128) :: fprfx
  character(len=128) :: path_data
  integer, dimension(3) :: dims
  integer :: nt
  real*8 :: d,r

  nx = 50
  ny = 50
  nz = 2
  path_data = "/scratch/uni/ifmto/u241161/myqg/test/"

  call allocate_myqg_module

  dx = 1.0
  dy = 1.0
  dz = (/ 50.d0, 50.d0 /)

  nt = 1000
  dt = 1.e3
  tstep = 0
  t_start = dt * tstep
  time    = t_start
  t_end   = dt * nt
  timeio  = t_end / 40. 

  diffPVh = 0.1 * dx**2/dt/4.0

  do i=1,nx
    xu(i) = dx*(i-1)
    xt(i) = dx*(i-0.5)
  enddo 

  do j=1,ny
    yu(j)=dy*(j-1)
    yt(j)=dy*(j-0.5)
  enddo 

  do k=2,nz
    zu(k)=zu(k-1) + dz(k)
  enddo
  do k=1,nz
    zt(k)=zu(k) + 0.5*dz(k)
  enddo

  do i=1,nx
    do j=1,ny
      do k=1,nz
        recepvol(i,j,k) = 1.0 / (dx*dy*dz(k))
      enddo
    enddo
  enddo

  Lx = nx*dx
  Ly = ny*dy

  do k=1,nz
    do i=1-ox,nx+ox
      do j=1-ox,ny+ox
        !u(i,j,k) = sin(4*xu(i)/Lx*2*pi)*cos(yu(j)/Ly*2*pi)
        !u(i,j,k) = (xu(i)/Lx)**2 + (yu(j)/Ly)**2 
        !u(i,j,k) = i*j*k
        !u(i,j,k) = 1.
        !v(i,j,k) = 1.
        psi(i,j,k) = - dy*j + dx*i
      enddo
    enddo
  enddo
  !fprfx = "test"
  !call write_3d(fprfx, u)

  f0 = 1e-7
  beta = 0.0
  do j=1,ny
    fCoru = f0 + beta*yu(j)
    fCort = f0 + beta*yt(j)
  enddo
  fCoru=0.0

  do i=1,nx
    do j=1,ny
      do k=1,nz
        !pv(i,j,k) = exp(-( ((xt(i)-Lx/2.)/(0.1*Lx))**2 + ((yt(j)-Ly/2.)/(0.1*Ly))**2 ))
        !pv(i,j,k) = (2*pi)**2*(-1/Lx**2-1/Ly**2) * sin(xu(i)/Lx*2*pi)*cos(yu(j)/Ly*2*pi)
        d = 0.1 * (Lx**2+Ly**2)**0.5
        r = ( (xt(i)-Lx/2.)**2 + (yt(j)-Ly/2.)**2 )**0.5
        !pv(i,j,k) = 1.0 / ( (r/d) + 0.1 )
        pv(i,j,k) = 4.0/d**2* ( r**2/d**2 - 1) * exp(-  r**2/d**2  )
        ! pv(i,j,k) = exp(- r/d)
        forc(i,j,k) = pv(i,j,k)
        !pv(i,j,k) = j*nx + i
      enddo
    enddo
  enddo

!  write(*,*) sngl(pv(:,1,1))
!  call cyclic_exchange(pv)
!  write(*,*) "start"
!  do j=1-ox,ny+ox
!    write(*,*) sngl(pv(:,j,1))
!  enddo
!  stop
  write(*,*) "Finished: Initialize variables"

  ! write grid data out
  call write_3d("XC        ", xt, (/ nx,  1,  1 /), -1, path_data)
  call write_3d("YC        ", yt, (/  1, ny,  1 /), -1, path_data)
  call write_3d("RC        ", zt, (/  1,  1, nz /), -1, path_data)
  call write_3d("XG        ", xu, (/ nx,  1,  1 /), -1, path_data)
  call write_3d("YG        ", yu, (/  1, ny,  1 /), -1, path_data)
  !call write_3d("RF        ", zu, (/  1,  1, nz /), -1, path_data)
  call write_3d("RF        ", (/ zu(:), 0.d0 /), (/  1,  1, nz+1 /), -1, path_data)

  !call cyclic_exchange(pv)
  !call cyclic_exchange(forc)
  do k=1,1!nz
    call solve_poisson_cg(1-ox,nx+ox,1-ox,ny+ox,dx,dy,forc(:,:,k),psi(:,:,k),max_itt,crit, est_error)
  enddo
  call write_3d("psi       ", psi(1:nx,1:ny,1:nz),  (/ nx,  ny,  nz /), 10, path_data)
  call write_3d("pv        ", pv (1:nx,1:ny,1:nz),  (/ nx,  ny,  nz /), 10, path_data)
  stop

  ! main time stepping loop
  write(*,*) "Start: Time stepping"
  do while ( time < t_end )
    tstep = tstep + 1
    time  = tstep*dt

    ! calculate velocity
    call calc_curl_psi
    call cyclic_exchange(u)
    call cyclic_exchange(v)

    ! calculate tendencies
    call calc_Gpv

    ! pv forward step
    if ( tstep==1 ) then
      Gpvm1 = Gpv
    endif
    pv(1:nx,1:ny,1:nz) = pv(1:nx,1:ny,1:nz) + ( (1.5+epsab)*Gpv - (0.5+epsab)*Gpvm1 ) * dt

    ! reset tendencies
    Gpvm1 = Gpv
    Gpv = 0.0

    ! do cyclic_exchange (Does this occur at correct place???)
    !call cyclic_exchange(psi)
    call cyclic_exchange(pv)

    ! do model I/O
    if ( floor(time/timeio) == time/timeio ) then
      write(*,*) "Model I/O at tstep ", tstep
      call write_3d('pv        ', pv (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep, path_data)
      call write_3d('psi       ', psi(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep, path_data)
      call write_3d('u         ', u  (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep, path_data)
      call write_3d('v         ', v  (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep, path_data)
    endif

  enddo

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
  !call calc_pvdiff
end subroutine calc_Gpv

subroutine calc_curl_psi
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  do i=1-ox,nx+ox-1
    do j=1-ox,ny+ox
      do k=1,nz
        u(i,j,k) = - (psi(i,j+1,k)-psi(i,j,k)) / dy
      enddo
    enddo
  enddo
  do i=1-ox,nx+ox
    do j=1-ox,ny+ox-1
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
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox,nz) :: fZon
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox,nz) :: fMer

  do i=1-ox,nx+ox-1
    do j=1-ox,ny+ox
      do k=1,nz
        fZon(i,j,k) = -diffPVh * (pv(i+1,j,k)-pv(i,j,k))/dx * dy*dz(k)
      enddo
    enddo
  enddo

  do i=1-ox,nx+ox
    do j=1-ox,ny+ox-1
      do k=1,nz
        fMer(i,j,k) = -diffPVh * (pv(i,j+1,k)-pv(i,j,k))/dy * dx*dz(k)
      enddo
    enddo
  enddo

  ! diffusive tendency
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
end subroutine calc_pvdiff

subroutine calc_pvadv
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
        uTrans(i,j,k) = 0.25 * ( u(i,j,k) + u(i+1,j,k) + u(i,j+1,k) + u(i+1,j+1,k) )*dy*dz(k)
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
        vTrans(i,j,k) = 0.25 * ( v(i,j,k) + v(i+1,j,k) + v(i,j+1,k) + v(i+1,j+1,k) )*dx*dz(k)
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
end subroutine calc_pvadv

subroutine cyclic_exchange(var)
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  real*8, dimension(1-ox:nx+ox,1-ox:ny+ox,nz), intent(inout) :: var
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

subroutine write_3d(fprfx, outdata, dims, tstepout, path_data)
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  character(len=10),            intent(in)  :: fprfx
  integer, dimension(3),        intent(in)  :: dims
  real*8, dimension(dims(1),dims(2),dims(3)),  intent(in)  :: outdata
  integer,                      intent(in)  :: tstepout
  character(len=128),           intent(in)  :: path_data
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  character(len=128)    :: fname
  character(len=10)     :: tstepstr 

  !if (.not. present(path_name)) then
  !  path_name = "./"
  !endif

! write .data file
  !if (present(tstepout)) then
  if (tstepout>0) then
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
  if (tstepout>0) then
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

subroutine myqg_error(errormsg)
  implicit none
  character(len=*), intent(in) :: errormsg
  write(*,*) "::: Error: " // trim(errormsg) // " :::" 
  stop 1
end subroutine myqg_error
