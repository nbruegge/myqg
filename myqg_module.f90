
module myqg_module
  implicit none

  real*8, parameter  :: pi      = 3.14159265358979323846264338327950588

! grid parameters
  integer :: nx
  integer :: ny
  integer :: nz

  real*8  :: dx
  real*8  :: dy
  real*8  :: dz

  real*8  :: Lx
  real*8  :: Ly
  real*8  :: Lz

! time stepping parameters
  real*8  :: dt
  integer :: tstep
  real*8  :: t_start
  real*8  :: t_end
  real*8  :: time

  real*8, allocatable, dimension(:)       :: xt, xu
  real*8, allocatable, dimension(:)       :: yt, yu
  real*8, allocatable, dimension(:)       :: zt, zu

  logical :: testlog

! mpi parameters
  integer :: numpr
  integer :: npr
  integer :: is, ie
  integer :: js, je
  
! model variables
  real*8, allocatable, dimension(:,:,:)   :: u, v, w
  real*8, allocatable, dimension(:,:,:)   :: b, p
  real*8, allocatable, dimension(:,:,:)   :: psi, pv
  real*8, allocatable, dimension(:,:,:)   :: pt

! tendencies for diagnostic variables
  real*8, allocatable, dimension(:,:,:)   :: Gpv  
  real*8, allocatable, dimension(:,:,:)   :: Gpt
 

! parameters for I/O
  character(len=24) :: endian = "big_endian"
  integer           :: bytes  = 4
  integer           :: fid    = 25

end module myqg_module

subroutine allocate_myqg_module
!=======================================================================
! allocate all arrays within main module
!=======================================================================
  use myqg_module
  implicit none
  
  allocate( xt(nx), xu(nx) ); xt=0; xu=0
  allocate( yt(ny), yu(ny) ); yt=0; yu=0
  allocate( zt(nz), zu(nz) ); zt=0; zu=0

  allocate( u(nx,ny,nz), v(ny,ny,nz) ); u=0; v=0
  allocate( w(nx,ny,nz), p(ny,ny,nz) ); w=0; p=0
  allocate( b(nx,ny,nz) ); b=0
  allocate( psi(nx,ny,nz), pv(ny,ny,nz) ); psi=0; pv=0
  allocate( pt(nx,ny,nz) ); pt=0
  allocate( Gpv(nx,ny,nz) ); Gpv=0
  allocate( Gpt(nx,ny,nz) ); Gpt=0

end subroutine allocate_myqg_module
