subroutine initialize_setup
  use myqg_module
  implicit none

!LOCAL VARIABLES:   ======================================== 
  integer :: i,j,k,l,n
  !real*8, allocatable, dimension(:,:,:,:) :: matA
  real*8  :: d, r
  integer :: nt 
  real*8  :: Lb, Lek
  real*8  :: R2, R3
  real*8  :: tau0

  character(len=20) :: fmtrea
  character(len=20) :: fmtint

  real*8 :: tmpreal

! ================================================================================ 
! initialize model
! ================================================================================ 
  nx = 64
  ny = 64
  nz = 2
  !path_data = "/scratch/uni/ifmto/u241161/myqg/test/"
  path_data = "/Users/nbruggemann/work/test_myqg/"

  call allocate_myqg_module
  !allocate( matA(3,3,nz,3) ); matA=0.0

  Lb = 1920.0e3
  Lx = 2.0*Lb
  Ly = 2.0*Lb

  dx  = Lx/nx
  dy  = Ly/ny 

  f0    = 1e-4
  beta  = 2.e-11

  ! three layer example (use nz=3)
!  Hk(1)   =  250.0
!  Hk(2)   =  750.0
!  Hk(3)   = 1000.0
!  R2      = 40.0e3
!  R3      = 23.0e3
!  ! Lr = N H / f
!  ! R1 = sqrt( -g/rho_0*(rho(1)-rho(2))/(zt(1)-zt(2)) ) * Hk(1) / f0
!  ! R1 = ( gred(k)*Hk(k) )**0.5 / f0
!  rho(1)  = 1000.
!  gred(1) = 10.
!  gred(2) = (R2*f0)**2 / Hk(2)
!  gred(3) = (R3*f0)**2 / Hk(3)
! ! do k=2,nz
! !   rho(k) = rho(k-1) + 10.0
! !   gred(k) = gred(1)/rho(1) * (rho(k) - rho(k-1))
! ! enddo

! two layer example (use nz=2)
  tmpreal = 0.
  Hk(1)   = 500.
  Hk(2)   = 500.+tmpreal
  gred(1) = 10.
  gred(2) = 0.02 
  !gred(2) = 10. 
  rho(1)  = 1000.
  
  ! time stepping
  nt      = 500
  !nt      = 1
  !dt      = 1200.
  dt      = 86400./10.
  !dt      = 3*1200.
  tstep   = 0
  t_start = dt * tstep
  time    = t_start
  t_end   = dt * nt
  timeio        = t_end / 20. 
  time_monitor  = t_end / 20.

  !diffPVh = dx**2/dt/4.0
  ! beta v = - Ah v_xxx
  ! beta   = Ah / delata**3
  diffPVh = beta * ( 2*dx )**3
  !diffPVh = 10

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
    zu(k)=zu(k-1) + Hk(k-1)
  enddo
  do k=1,nz
    zt(k)=zu(k) + 0.5*Hk(k)
  enddo

  !tau0 = 0.8
  !Lek  = Lb-dy/2.0
  tau0 = 0.1
  do j=1-ox,ny+ox
    do i=1-ox,nx+ox
      !wek(i,j) = - pi*tau0/(rho(1)*f0*Lek) * sin(pi*(Lek+yu(j))/Lek)
      wek(i,j) = - tau0/(rho(1)*f0) * 2*pi/(Ly-dy) * sin(2*pi*yu(j)/(Ly-dy))
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
    !fCoru(j) = f0 + beta*( yu(j) - Ly/2.0 )
    !fCort(j) = f0 + beta*( yt(j) - Ly/2.0 )

    fCoru(j) = f0 + beta * yu(j)
    fCort(j) = f0 + beta * yt(j)
  enddo

  ! inhomogenity of poisson equation
  !call solve_poisson_cg(1-ox,nx+ox,1-ox,ny+ox,nz,forc,psi,max_itt,crit, est_error, doio)

  ! calculate initial relative vorticity from streamfunction
  tstep = 0
  call make_matrix
  call matrix_prod(1-ox, nx+ox, 1-ox, ny+ox, nz, matA, matRelx, matRely, matStr, psi, pv)
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

  fmtint="(A,I14,A)"
  fmtrea="(A,ES14.8,A)"
  !fmtvec="(A,2ES14.8,A)"
  open(20, file="parameter.m", status="unknown")
  write(20,fmtint) "nx = ",nx,";"
  write(20,fmtint) "ny = ",ny,";"
  write(20,fmtint) "nz = ",nz,";"
  write(20,fmtint) "nt = ",nt,";"
  write(20,fmtint) ""
  write(20,fmtrea) "dx = ",dx,";"
  write(20,fmtrea) "dy = ",dy,";"
  write(20,fmtrea) "dt = ",dt,";" 
  write(20,fmtint) ""
  write(20,fmtrea) "Lx = ",Lx,";"
  write(20,fmtrea) "Ly = ",Ly,";"
  write(20,fmtint) ""
  write(20,fmtrea) "f0      = ",f0,";" 
  write(20,fmtrea) "beta    = ",beta,";" 
  write(20,fmtrea) "tau0    = ",tau0,";" 
  write(20,fmtrea) "diffPVh = ",diffPVh,";" 
  write(20,fmtint) ""
  write(20,"(A,2ES14.8,A)") "gred = ",gred,";" 
  write(20,"(A,2ES14.8,A)") "Hk = ",Hk,";" 
  close(20)


end subroutine initialize_setup
