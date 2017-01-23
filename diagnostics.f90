
subroutine diagnostics 
  use myqg_module
  implicit none
  integer :: i,j,k

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

    !!! calculate layer width and reciprocal volume
    !!do j=1-ox,ny+ox
    !!  do i=1-ox,nx+ox
    !!    do k=1,nz
    !!      dz(i,j,k)       = Hk(k) + hpr(i,k,k)
    !!      recepvol(i,j,k) = 1.0 / (dx*dy*dz(i,j,k))
    !!    enddo
    !!  enddo
    !!enddo

    ! calculate velocity
    call calc_curl_psi
    call calc_pv_components
    !call cyclic_exchange(u)
    !call cyclic_exchange(v)
end subroutine diagnostics

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

subroutine calc_pv_components
  !!! this should be the same function as in solve_posson_cg.f90/matrix_prod
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
!!!  real*8, dimension(is:ie,js:je,nz), intent(in)  :: d
  !real*8, dimension(is:ie,js:je,3,3), intent(out) :: matA
!!!  real*8, dimension(3,3,nz,3), intent(in) :: matA
!!!  real*8, dimension(3,3,nz,3), intent(in) :: matRelx
!!!  real*8, dimension(3,3,nz,3), intent(in) :: matRely
!!!  real*8, dimension(3,3,nz,3), intent(in) :: matStr
!OUTPUT PARAMETERS: ======================================== 
!  real*8, dimension(is:ie,js:je,nz), intent(out) :: pv_r, pv_s
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: is,ie,js,je
  integer                                     :: i, j, k, ii, jj, kk
  integer                                     :: k1, k2
  pv_r = 0.0
  pv_s = 0.0
  is = 1-ox
  ie = nx+ox
  js = 1-ox
  je = ny+ox

  do k=1,nz
    if ( k==1 ) then
      k1 =  0
      k2 =  1
    elseif ( k==nz ) then
      k1 = -1
      k2 =  0
    else
      k1 = -1
      k2 =  1
    endif
    do j=js+1,je-1
      do i=is+1,ie-1
        do ii = -1,1
          pv_r(i,j,k) = pv_r(i,j,k) + matRelx(ii+2,2,k,2)*psi(i+ii,j,k)
        enddo
        do jj = -1,1
          pv_r(i,j,k) = pv_r(i,j,k) + matRely(ii,jj+2,k,2)*psi(i,j+jj,k)
        enddo

        do kk = k1,k2
          pv_s(i,j,k) = pv_s(i,j,k) + matStr(2,2,k,kk+2)*psi(i,j,k+kk)
        enddo
      enddo
    enddo
  enddo
end subroutine calc_pv_components

subroutine write_snapshot
  use myqg_module
  implicit none
  write(*,*) "================================================================================"
  write(*,*) "Model I/O at tstep ", tstep
  write(*,*) "  snap/all_snaps = ", int(time/timeio),'/', int(t_end/timeio)
  call write_3d('pv        ', pv (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('pvp       ', pvp(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('pvr       ', pvr(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('pv_r      ', pv_r(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('pv_s      ', pv_s(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('psi       ', psi(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  !call write_3d('u         ', u  (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  !call write_3d('v         ', v  (1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  call write_3d('hpr       ', hpr(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  !call write_3d('Gpvadv    ', Gpvadv(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  !call write_3d('Gpvdif    ', Gpvdif(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  !call write_3d('Gpvfor    ', Gpvfor(1:nx,1:ny,1:nz), (/ nx, ny, nz /), tstep)
  write(*,*) "================================================================================"
end subroutine write_snapshot

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
