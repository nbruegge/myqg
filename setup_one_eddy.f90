subroutine initialize_setup
  use myqg_module
  implicit none

!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k,l,n
  real*8, allocatable, dimension(:,:,:) :: A
  real*8  :: d, r
  integer :: nt 

! ================================================================================ 
! initialize model
! ================================================================================ 
  nx = 50
  ny = 50
  nz = 3
  !path_data = "/scratch/uni/ifmto/u241161/myqg/test/"
  path_data = "../test_myqg/"

  call allocate_myqg_module
  allocate( A(3,3,nz) ); A=0.0

  dx  = 1.0e3
  dy  = 1.0e3

  f0    = 1e-4
  beta  = 0.0

  ! this results in N2 = 1e-3
  Hk      = 100.0
  rho(1)  = 1000.
  gred(1) = 10.
  do k=2,nz
    rho(k) = rho(k-1) + 10.0
    gred(k) = gred(1)/rho(1) * (rho(k) - rho(k-1))
  enddo

  ! time stepping
  nt      = 100
  !nt      = 1
  dt      = 1.e1
  tstep   = 0
  t_start = dt * tstep
  time    = t_start
  t_end   = dt * nt
  timeio        = t_end / 10. 
  time_monitor  = t_end / 10.

  diffPVh = 0.1 * dx**2/dt/4.0
  
  ! grid
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

  do k=2,nz
    zu(k)=zu(k-1) + Hk(k)
  enddo
  do k=1,nz
    zt(k)=zu(k) + 0.5*Hk(k)
  enddo

  ! initial pv
  do i=1,nx
    do j=1,ny
      do k=1,nz
        d = 0.1 * (Lx**2+Ly**2)**0.5
        r = ( (xt(i)-Lx/2.)**2 + (yt(j)-Ly/2.)**2 )**0.5
        !pv (i,j,k) = 1.e2 * 4.0/d**2* ( r**2/d**2 - 1) *  exp(- r**2/d**2 )
        psi(i,j,k) = 1.e2 *                               exp(- r**2/d**2 )
      enddo
    enddo
  enddo

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
    fCoru = f0 + beta*yu(j)
    fCort = f0 + beta*yt(j)
  enddo

  ! inhomogenity of poisson equation
  do j=1-ox,ny+ox
    do i=1-ox,nx+ox
      do k=1,nz
        forc(i,j,k) = pv(i,j,k) - beta * yu(j) 
      enddo
    enddo
  enddo

  !call solve_poisson_cg(1-ox,nx+ox,1-ox,ny+ox,nz,forc,psi,max_itt,crit, est_error, doio)

  call make_matrix(1-ox, nx+ox, 1-ox, ny+ox, A)
  call matrix_prod(1-ox, nx+ox, 1-ox, ny+ox, nz, A, psi, pv)
  !call solve_poisson_cg(1-ox, nx+ox, 1-ox, ny+ox, nz, pv,   hpr, max_itt, crit, est_error, doio)

end subroutine initialize_setup
