

subroutine solve_poisson_cg(is,ie,js,je,dx,dy,forc,sol,max_itt,crit)
  !use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  real*8, intent(in)                          :: dx, dy
  real*8, dimension(is:ie,js:je), intent(in)  :: forc
  integer, intent(in)                         :: max_itt       
  real*8, intent(in)                          :: crit
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,js:je), intent(out) :: sol
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k,n
  !real*8, dimension(is:ie,js:je,3,3)          :: A
  real*8, dimension(3,3)                      :: A
  real*8, dimension(is:ie,js:je)              :: res, d, z
  logical, save                               :: first=.true.
  real*8                                      :: alpha, beta
  real*8                                      :: resTres, resnTresn, dTz
  real*8, dimension(is:ie,js:je)              :: appinvA

! sol:  best guess of solution
! res:  residual vector res = forc - A*sol 
! d:    search direction sol_new = sol + beta*d
!       sol_new = sol + alpha*d
!       d_new   = res + beta*d
! initialize variables
  A   = 0.0
  res = 0.0
  sol = 0.0
  d   = 0.0
  z   = 0.0

! make coefficient matrix A
  if ( first) then
    call make_matrix(is,ie,js,je,dx,dy,A)
    first = .false.
  endif

! initial guess
  call approx_inv_mat(is,ie,js,je,A,appinvA)
  call prod_invA_forc(is,ie,js,je,appinvA,forc,sol)
  call cyclic_exchange_2d(sol)

! r0 = b - A*x
  call matrix_prod(is,ie,js,je,A,sol,res)
  call cyclic_exchange_2d(res)
  res = forc - res
! d0 = r0
  d = res

  do n=1,max_itt
    ! z = A * d
    call matrix_prod(is,ie,js,je,A,d,z)
    call cyclic_exchange_2d(z)

    ! r^T*r and d^T*z
    call scalar_prod(is,ie,js,je,res,res,resTres)
    call scalar_prod(is,ie,js,je,d,z,dTz)

    ! alpha = r^T*r / d^T*z
    alpha = resTres / dTz

    ! x_n+1 = x_n + alpha * d
    sol = sol + alpha * d

    ! res_n+1 = res_n - alpha*z
    res = res - alpha * z

    ! res_n+1^T*res_n+1
    call scalar_prod(is,ie,js,je,res,res,resnTresn)

    ! beta = res_n+1^T*res_n+1 / res^T*res  
    beta = resnTresn / resTres

    ! d_n+1 = r_n+1 + beta*d_n
    d = res + beta * d

    write(*,*) abs((beta-1.0))
    if ( abs((beta-1.0)) < crit ) then
    ! stop iteration
      goto 100
    endif

  enddo
100 write(*,*) "Converged after n = ", n 

end subroutine solve_poisson_cg

!subroutine make_matrix2(A)
!  implicit none
!  real*8, dimension(3,3), intent(out):: A
!  integer :: i,j
!  A = 0.0
!  A(1,2) = -1.0
!  A(2,1) = 2.0
!  A(3,3) = 3.0
!end subroutine make_matrix2
!
!subroutine matrix_prod2(A,v,mprod)
!  real*8, dimension(3,3), intent(in):: A
!  real*8, dimension(3), intent(in):: v
!  real*8, dimension(3), intent(out):: mprod
!  mprod = 0.0
!  do j=1,3
!    mprod(i) = mprod(i) + A(i,j)*v(j) 
!  enddo
!end subroutine matrix_prod2

subroutine make_matrix(is,ie,js,je,dx,dy,A)
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  real*8, intent(in)                          :: dx, dy
!OUTPUT PARAMETERS: ======================================== 
  !real*8, dimension(is:ie,js:je,3,3), intent(out):: A
  real*8, dimension(3,3), intent(out):: A
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k
  A = 0.0
  do j=js,je
    do i=is,ie
      !A(i,j,1,2) =  1.0/dx**2
      !A(i,j,3,2) =  1.0/dx**2
      !A(i,j,2,2) = -2.0/dx**2-2.0/dy**2
      !A(i,j,2,1) =  1.0/dy**2
      !A(i,j,2,3) =  1.0/dy**2

      A(1,2) =  1.0/dx**2
      A(3,2) =  1.0/dx**2
      A(2,2) = -2.0/dx**2-2.0/dy**2
      A(2,1) =  1.0/dy**2
      A(2,3) =  1.0/dy**2
    enddo
  enddo
end subroutine make_matrix

subroutine scalar_prod(is,ie,js,je,p1,p2,sprod)
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  real*8, dimension(is:ie,js:je), intent(in)  :: p1, p2
!OUTPUT PARAMETERS: ======================================== 
  real*8, intent(out)                         :: sprod
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j
  sprod = 0.0
  do i=is,ie
    do j=js,je
      sprod = sprod  + p1(i,j)*p2(i,j)
    enddo
  enddo
end subroutine scalar_prod

subroutine matrix_prod(is,ie,js,je,A,v,mprod)
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  real*8, dimension(is:ie,js:je), intent(in)  :: v
  !real*8, dimension(is:ie,js:je,3,3), intent(out) :: A
  real*8, dimension(3,3), intent(out) :: A
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,js:je), intent(out) :: mprod
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i, j, ii, jj
  mprod = 0.0
  do j=js+1,je-1
    do i=is+1,ie-1
      do ii = -1,1
        do jj = -1,1
          !mprod(i,j) = mprod(i,j) + A(i,j,ii+2,jj+2)*v(i+ii,j+jj)
          mprod(i,j) = mprod(i,j) + A(ii+2,jj+2)*v(i+ii,j+jj)
        enddo
      enddo
    enddo
  enddo
end subroutine matrix_prod

subroutine approx_inv_mat(is,ie,js,je,A,appinvA)
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  !real*8, dimension(is:ie,js:je,3,3), intent(in):: A
  real*8, dimension(3,3), intent(in):: A
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,js:je), intent(out):: appinvA
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k
  appinvA=0.0
  ! copy diagonal elements
  do j=js,je
    do i=is,ie
      !appinvA(i,j) = A(i,j,2,2)
      appinvA(i,j) = A(2,2)
    enddo
  enddo
  ! invert A
  do j=js,je
    do i=is,ie
      if ( appinvA(i,j) .ne. 0.0 ) then
        appinvA(i,j) = 1/appinvA(i,j)
      else
        appinvA(i,j) = 0.0
      endif
    enddo
  enddo
end subroutine approx_inv_mat

subroutine prod_invA_forc(is,ie,js,je,appinvA,forc,sol)
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  real*8, dimension(is:ie,js:je), intent(in)  :: appinvA
  real*8, dimension(is:ie,js:je), intent(in)  :: forc 
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,js:je), intent(out) :: sol
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k
  sol=0.0
  do j=js,je
    do i=is,ie
      sol(i,j) = appinvA(i,j)*forc(i,j)
    enddo
  enddo
end subroutine prod_invA_forc
