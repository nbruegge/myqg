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
