
program main
  use myqg_module
  implicit none

  integer :: i,j,k,l,n
  character(len=128) :: fname
  character(len=48) :: fprfx
  integer :: nt

  nx = 120
  ny = 60
  nz = 2

  call allocate_myqg_module

  dx = 10.e3
  dy = 10.e3
  dz = 5.

  do i=1,nx
    xu(i) = dx*(i-1)
    xt(i) = dx*(i-0.5)
  enddo 

  do j=1,ny
    yu(j)=dy*(j-1)
    yt(j)=dy*(j-0.5)
  enddo 

  Lx = nx*dx
  Ly = ny*dy

  do i=1,nx
    do j=1,ny
      u(i,j,1) = sin(4*xu(i)/Lx*2*pi)*cos(yu(j)/Ly*2*pi)
      !u(i,j,1) = (xu(i)/Lx)**2 + (yu(j)/Ly)**2 
      !u(i,j,1) = i*j
      if (i == 0) then
        write(*,*) u(i,j,1)
      endif
    enddo
  enddo
  fprfx = "test"
  call write_3d(fprfx, u)
 
  nt = 1000
  dt = 10.
  tstep = 0
  t_start = dt * tstep
  time    = t_start
  t_end   = dt * nt

  ! main time stepping loop
  do while ( time < t_end )
    write(*,*) tstep
    tstep = tstep + 1
    time  = tstep*dt

    call calc_Gpv
    pv = pv + ( (1.5+epsab)*Gpv - (0.5+epsab)*Gpvm1 ) *dt
    Gpvm1 = Gpv

    ! do model I/O
    if ( tstep/tstepio == tstep/tstepio ) then
      fprfx = "pv." // tstepform
      call write_3d(fprfx, pv)
    endif

  enddo

  write(*,*) "nx = ", nx, " ny = ", ny, " nz = ", nz

  write(*,*) "All done"
end program main

subroutine calc_Gpv
  implicit none
!INPUT PARAMETERS:  ======================================== 
  real*8, dimension(nx,ny,nz)   :: dummy
  character(len=10)             :: dummy 
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(nx,ny,nz)   :: dummy
!LOCAL VARIABLES:   ======================================== 

  Gpv = pv*time/t_end

end subroutine calc_Gpv


subroutine write_3d(fprfx, outdata)
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  real*8, dimension(nx,ny,nz)   :: outdata
  character(len=48)             :: fprfx
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 

  open( unit=fid, file=trim(fprfx)//".data", form='unformatted', status='replace', &
      & access='direct', recl=bytes*nx*ny*nz, convert=endian )
  write(fid, rec=1) sngl(outdata)
  close(fid)

  open( unit=fid, file=trim(fprfx)//".meta", form='formatted', status='replace', &
      & access='sequential')
  write(fid,*) "dims = [", nz, ",", ny, ",", nx, "]"
  close(fid)
end subroutine write_3d
