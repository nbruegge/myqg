
module myqg_module
  implicit none

  real*8, parameter  :: pi      = 3.14159265358979323846264338327950588

! grid parameters
  integer :: nx
  integer :: ny
  integer :: nz
  integer :: ox=1

  real*8  :: dx
  real*8  :: dy
  real*8, allocatable, dimension(:,:,:)  :: dz

  real*8  :: Lx
  real*8  :: Ly
  real*8  :: Lz

  logical :: cyclic_x=.false.
  logical :: cyclic_y=.false.

  real*8, allocatable, dimension(:)       :: xt, xu
  real*8, allocatable, dimension(:)       :: yt, yu
  real*8, allocatable, dimension(:)       :: zt, zu

  real*8, allocatable, dimension(:,:,:)   :: recepvol 

! time stepping parameters
  real*8  :: dt
  integer :: tstep
  real*8  :: t_start
  real*8  :: t_end
  real*8  :: time
  real*8  :: timeio
  real*8  :: epsab=0.01
  real*8  :: time_monitor

! parameters for Poisson solver
  !integer :: max_itt=100
  integer :: max_itt=1000     ! pyOM value
  !real*8  :: crit=0.1
  real*8  :: crit=1e-12       ! pyOM value
  !real*8  :: crit=1e-4       ! Carsten's qg_model value
  real*8  :: est_error
  real*8, allocatable, dimension(:,:,:,:)   :: matA
  real*8, allocatable, dimension(:,:,:,:)   :: matRelx
  real*8, allocatable, dimension(:,:,:,:)   :: matRely
  real*8, allocatable, dimension(:,:,:,:)   :: matStr
  real*8, allocatable, dimension(:,:,:)     :: C 

! mpi parameters
  integer :: numpr
  integer :: npr
  !integer :: is, ie
  !integer :: js, je
  
! model variables
  real*8, allocatable, dimension(:,:,:)   :: u, v, w
  real*8, allocatable, dimension(:,:,:)   :: b, p
  real*8, allocatable, dimension(:,:,:)   :: psi, pv
  real*8, allocatable, dimension(:,:,:)   :: pvp, pvr
  real*8, allocatable, dimension(:,:,:)   :: pt
  real*8, allocatable, dimension(:,:,:)   :: hpr

! tendencies for diagnostic variables
  real*8, allocatable, dimension(:,:,:)   :: Gpv
  real*8, allocatable, dimension(:,:,:)   :: Gpvadv, Gpvdif, Gpvfor
  real*8, allocatable, dimension(:,:,:)   :: Gpvm1
  real*8, allocatable, dimension(:,:,:)   :: Gpt
  real*8, allocatable, dimension(:,:,:)   :: Gptm1

  real*8, allocatable, dimension(:)       :: rho, gred, Hk

! forcing for poisson equation (forc = pv - beta*y)
  real*8, allocatable, dimension(:,:,:)   :: forc 

! Ekman velocity
  real*8, allocatable, dimension(:,:)     :: wek
 
! parameters for I/O
  character(len=24) :: endian = "big_endian"
  integer           :: bytes  = 4
  integer           :: fid    = 25
  character(len=128) :: path_data="./"

! model parameters
  real*8, allocatable, dimension(:)       :: fCort, fCoru
  real*8 :: f0
  real*8 :: beta
  real*8 :: diffPVh 

end module myqg_module

subroutine allocate_myqg_module
!=======================================================================
! allocate all arrays within main module
!=======================================================================
  use myqg_module
  implicit none
  
  allocate( xt(1-ox:nx+ox), xu(1-ox:nx+ox) ); xt=0; xu=0
  allocate( yt(1-ox:ny+ox), yu(1-ox:ny+ox) ); yt=0; yu=0
  allocate( zt(nz), zu(nz) ); zt=0; zu=0
  allocate( fCort(ny), fCoru(ny) ); fCort=0; fCort=0

  allocate( dz      (1-ox:nx+ox,1-ox:ny+ox,nz) );                               dz=0
  allocate( recepvol(1-ox:nx+ox,1-ox:ny+ox,nz) );                               recepvol=0

  allocate( u    (1-ox:nx+ox,1-ox:ny+ox,nz), v(1-ox:nx+ox,1-ox:ny+ox,nz) );     u=0; v=0
  allocate( w    (1-ox:nx+ox,1-ox:ny+ox,nz), p(1-ox:nx+ox,1-ox:ny+ox,nz) );     w=0; p=0
  allocate( b    (1-ox:nx+ox,1-ox:ny+ox,nz) );                                  b=0
  allocate( psi  (1-ox:nx+ox,1-ox:ny+ox,nz), pv(1-ox:nx+ox,1-ox:ny+ox,nz) );    psi=0; pv=0
  allocate( pvp(1-ox:nx+ox,1-ox:ny+ox,nz), pvr(1-ox:nx+ox,1-ox:ny+ox,nz) );    pvp=0; pvr=0
  allocate( pt   (1-ox:nx+ox,1-ox:ny+ox,nz) );                                  pt=0
  allocate( hpr  (1-ox:nx+ox,1-ox:ny+ox,nz) );                                  hpr=0
  !allocate( Gpv(1-ox:nx+ox,1-ox:ny+ox,nz), Gpvm1(1-ox:nx+ox,1-ox:ny+ox,nz) ); Gpv=0; Gpvm1=0
  allocate( Gpv (nx,ny,nz), Gpvm1(nx,ny,nz) );                                  Gpv=0; Gpvm1=0
  allocate( Gpvadv(nx,ny,nz), Gpvdif(nx,ny,nz) );                                  Gpvadv=0; Gpvdif=0
  allocate( Gpvfor(nx,ny,nz)                   );                                  Gpvfor=0
  allocate( Gpt (1-ox:nx+ox,1-ox:ny+ox,nz), Gptm1(1-ox:nx+ox,1-ox:ny+ox,nz) );  Gpt=0; Gptm1=0
  allocate( forc(1-ox:nx+ox,1-ox:ny+ox,nz) );                                   forc=0
  allocate( wek(1-ox:nx+ox,1-ox:ny+ox) );                                       wek=0

  allocate( rho(nz), gred(nz), Hk(nz) );                                        rho=0; gred=0; Hk=0;

  allocate( matA(3,3,nz,3) );                                        matA=0;
  allocate( matRelx(3,3,nz,3) );                                     matRelx=0;
  allocate( matRely(3,3,nz,3) );                                     matRely=0;
  allocate( matStr(3,3,nz,3) );                                      matStr=0;
  allocate( C(1-ox:nx+ox,1-ox:ny+ox,nz) );                           C=0;

end subroutine allocate_myqg_module
