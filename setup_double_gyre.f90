subroutine initialize_setup
  use myqg_module
  implicit none

!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k,l,n
  real*8, allocatable, dimension(:,:,:) :: A
  real*8  :: d, r
  integer :: nt 
  real*8  :: Lb, Lek
  real*8  :: R2, R3
  real*8  :: tau0

! ================================================================================ 
! initialize model
! ================================================================================ 
  nx = 64
  ny = 64
  nz = 3
  path_data = "/scratch/uni/ifmto/u241161/myqg/test/"
  !path_data = "../test_myqg/"

  call allocate_myqg_module
  allocate( A(3,3,nz) ); A=0.0

  Lb = 1920.0e3
  Lx = 2.0*Lb
  Ly = 2.0*Lb

  dx  = Lx/nx
  dy  = Ly/ny 

  f0    = 1e-4
  beta  = 2.e-11

  Hk(1)   =  250.0
  Hk(2)   =  750.0
  Hk(3)   = 1000.0

  R2      = 40.0e3
  R3      = 23.0e3

  ! Lr = N H / f
  ! R1 = sqrt( -g/rho_0*(rho(1)-rho(2))/(zt(1)-zt(2)) ) * Hk(1) / f0
  ! R1 = ( gred(k)*Hk(k) )**0.5 / f0

  rho(1)  = 1000.
  gred(1) = 10.
  gred(2) = (R2*f0)**2 / Hk(2)
  gred(3) = (R3*f0)**2 / Hk(3)
 ! do k=2,nz
 !   rho(k) = rho(k-1) + 10.0
 !   gred(k) = gred(1)/rho(1) * (rho(k) - rho(k-1))
 ! enddo
  
  ! time stepping
  nt      = 10
  !nt      = 1
  !dt      = 1200.
  !dt      = 86400.
  dt      = 1200.
  tstep   = 0
  t_start = dt * tstep
  time    = t_start
  t_end   = dt * nt
  timeio        = t_end / 10. 
  time_monitor  = t_end / 10.

  diffPVh = dx**2/dt/4.0
  ! beta v = - Ah v_xxx
  ! beta   = Ah / delata**3
  !diffPVh = beta * ( 4*dx )**3

  write(*,*) "dx = ", dx, "dy = ", dy
  write(*,*) "gred = ", gred
  write(*,*) "diffPVh = ", diffPVh

  ! grid
  do i=1-ox,nx+ox
    xu(i) = dx*(i-1)
    xt(i) = dx*(i-0.5)
  enddo 

  do j=1-ox,ny+ox
    yu(j)=dy*(j-1)
    yt(j)=dy*(j-0.5)
  enddo 

  do k=2,nz
    zu(k)=zu(k-1) + Hk(k)
  enddo
  do k=1,nz
    zt(k)=zu(k) + 0.5*Hk(k)
  enddo

  tau0 = 0.8
  Lek  = Lb-dy/2.0
  do j=1-ox,ny+ox
    do i=1-ox,nx+ox
      wek(i,j) = - pi*tau0/(rho(1)*f0*Lek) * sin(pi*(Lek+yu(j))/Lek)
    enddo
  enddo

  ! initial pv
!  do i=1,nx
!    do j=1,ny
!      do k=1,nz
!        d = 0.1 * (Lx**2+Ly**2)**0.5
!        r = ( (xt(i)-Lx/2.)**2 + (yt(j)-Ly/2.)**2 )**0.5
!        !pv (i,j,k) = 1.e2 * 4.0/d**2* ( r**2/d**2 - 1) *  exp(- r**2/d**2 )
!        psi(i,j,k) = 1.e2 *                               exp(- r**2/d**2 )
!      enddo
!    enddo
!  enddo

!  do j=1-ox,ny+ox
!    do i=1-ox,nx+ox
!      do k=1,nz
!        !u(i,j,k) = sin(4*xu(i)/Lx*2*pi)*cos(yu(j)/Ly*2*pi)
!        !u(i,j,k) = (xu(i)/Lx)**2 + (yu(j)/Ly)**2 
!        !u(i,j,k) = i*j*k
!        !u(i,j,k) = 1.
!        !v(i,j,k) = 1.
!        psi(i,j,k) = - dy*j + dx*i
!      enddo
!    enddo
!  enddo

  do j=1,ny
    fCoru(j) = f0 + beta*( yu(j) - Ly/2.0 )
    fCort(j) = f0 + beta*( yt(j) - Ly/2.0 )
  enddo

  ! inhomogenity of poisson equation
  !call solve_poisson_cg(1-ox,nx+ox,1-ox,ny+ox,nz,forc,psi,max_itt,crit, est_error, doio)

  ! calculate initial relative vorticity from streamfunction
  call make_matrix(1-ox, nx+ox, 1-ox, ny+ox, A)
  call matrix_prod(1-ox, nx+ox, 1-ox, ny+ox, nz, A, psi, pv)
  !call solve_poisson_cg(1-ox, nx+ox, 1-ox, ny+ox, nz, pv,   hpr, max_itt, crit, est_error, doio)
  
  ! add planetary vorticity
  !write(*,*) "rel pv = ", pv(32,:,1)
  do j=1-ox,ny+ox
    do i=1-ox,nx+ox
      do k=1,nz
        pv(i,j,k) = pv(i,j,k) + beta * yu(j) 
        pvp(i,j,k) = pvp(i,j,k) + beta * yu(j)
        !pv(i,j,k) = pv(i,j,k) + fCoru(j) 
      enddo
    enddo
  enddo
  !write(*,*) "abs pv = ", pv(32,:,1)


end subroutine initialize_setup
