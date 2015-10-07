program main
  use myqg_module
  implicit none

  integer :: i,j,k,l,n
  integer :: is,ie,js,je
  character(len=128) :: fname
  character(len=128) :: fprfx
  !character(len=128) :: path_data
  integer, dimension(3) :: dims
  !real*8, allocatable, dimension(:,:,:) :: A
  integer :: nt
  logical :: doio
  real*8 :: tmpreal

  ! start and end indices (no parallelization yet)
  !is =  1 - ox
  !ie = nx + ox
  !js =  1 - ox
  !je = ny + ox

  call initialize_setup
  write(*,*) "Finished: Initialize variables"

  ! write grid data 
  call write_3d("XC        ", xt, (/ nx,  1,  1 /), -1)
  call write_3d("YC        ", yt, (/  1, ny,  1 /), -1)
  call write_3d("RC        ", zt, (/  1,  1, nz /), -1)
  call write_3d("XG        ", xu, (/ nx,  1,  1 /), -1)
  call write_3d("YG        ", yu, (/  1, ny,  1 /), -1)
  call write_3d("RF        ", (/ zu(:), zu(nz)+Hk(nz) /), (/  1,  1, nz+1 /), -1)

  ! write initial snapshot
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


    ! inhomogenity of poisson equation
    pvr = pv - pvp
    forc = pvr

    ! calculate streamfunction
    call solve_poisson_cg(1-ox, nx+ox, 1-ox, ny+ox, nz, forc, psi, max_itt, crit, est_error, doio, tstep, matA, C)
    
    ! apply boundary condition for psi
    !psi(1-ox:1,:,:)   = 0.0
    !psi(nx:nx+ox,:,:) = 0.0
    !psi(:,1-ox:1,:)   = 0.0
    !psi(:,ny:ny+ox,:) = 0.0

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
      call diagnostics
      call write_snapshot
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

subroutine calc_pvdiff
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  ! prepare outer domain pv values for no-flux boundary conditions
  !pvr(0,:,:)    = pvr(1,:,:)
  !pvr(nx+1,:,:) = pvr(nx,:,:)
  !pvr(:,0,:)    = pvr(:,1,:)
  !pvr(:,ny+1,:) = pvr(:,ny,:)
  ! prepare outer domain pv values for no-slip boundary conditions
  pvr(0,:,:)    = 0.0 
  pvr(nx+1,:,:) = 0.0 
  pvr(:,0,:)    = 0.0 
  pvr(:,ny+1,:) = 0.0 

  Gpvdif = 0.0
  do i=1,nx
    do j=1,ny
      do k=1,nz
        Gpvdif(i,j,k) =   diffPVh * (  1.0/dx**2 * ( pvr(i+1,j,k) - 2.0*pvr(i,j,k) + pvr(i-1,j,k) ) & 
                      &              + 1.0/dy**2 * ( pvr(i,j+1,k) - 2.0*pvr(i,j,k) + pvr(i,j-1,k) ) ) 
      enddo
    enddo
  enddo
  Gpv = Gpv + Gpvdif
end subroutine calc_pvdiff
  
subroutine calc_pvadv
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!OUTPUT PARAMETERS: ======================================== 
!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k
  real*8  :: J1, J2, J3
  ! Arakawa-Jacobian operator
  Gpvadv = 0.0
  do j=1,ny
    do i=1,nx
      do k=1,nz
        J1 = 1.0/(4.0*dx*dy) * ( &
           &   (psi(i+1,j,k)-psi(i-1,j,k)) * (pv (i,j+1,k)-pv (i,j-1,k)) &
           & - (psi(i,j+1,k)-psi(i,j-1,k)) * (pv (i+1,j,k)-pv (i-1,j,k)) &
           & )
        J2 = 1.0/(4.0*dx*dy) * ( &
           & + psi(i+1,j,k)*(pv (i+1,j+1,k)-pv (i+1,j-1,k)) - psi(i-1,j,k)*(pv (i-1,j+1,k)-pv (i-1,j-1,k) ) &
          !& - psi(i,j+1,k)*(pv (i+1,j+1,k)-pv (i-1,j+1,k)) - psi(i,j-1,k)*(pv (i+1,j-1,k)-pv (i-1,j-1,k) ) &
           & - psi(i,j+1,k)*(pv (i+1,j+1,k)-pv (i-1,j+1,k)) + psi(i,j-1,k)*(pv (i+1,j-1,k)-pv (i-1,j-1,k) ) &
           & )
        J3 = 1.0/(4.0*dx*dy) * ( &
          !&   pv (i+1,j,k)*(psi(i+1,j+1,k)-psi(i+1,j-1,k)) - pv (i-1,j,k)*(psi(i-1,j+1,k)-psi(i-1,j-1,k) ) &
           & - pv (i+1,j,k)*(psi(i+1,j+1,k)-psi(i+1,j-1,k)) + pv (i-1,j,k)*(psi(i-1,j+1,k)-psi(i-1,j-1,k) ) &
           & + pv (i,j+1,k)*(psi(i+1,j+1,k)-psi(i-1,j+1,k)) - pv (i,j-1,k)*(psi(i+1,j-1,k)-psi(i-1,j-1,k) ) &
           & )
        Gpvadv(i,j,k) = Gpvadv(i,j,k) - (J1 + J2 + J3) / 3.0
      enddo
    enddo
  enddo
  Gpv = Gpv + Gpvadv
end subroutine calc_pvadv

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

subroutine myqg_error(errormsg)
  implicit none
  character(len=*), intent(in) :: errormsg
  write(*,*) "::: Error: " // trim(errormsg) // " :::" 
  stop 1
end subroutine myqg_error

