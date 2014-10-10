

subroutine solve_poisson_qg(is,ie,js,je,dx,dy,forc,sol)
  !use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  real*8, intent(in)                          :: dx, dy
  real*8, dimension(is:ie,js:je), intent(in)  :: forc
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,,js:je), intent(out):: sol
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k,n
  integer                                     :: max_itt       

! initialize

! r0 = b - A*x
! d0 = r0
  res = forc - A*sol
  d = res

! make A
  call make_matrix(is,ie,js,je,dx,dy,A)


  do n=1,max_itt
    ! z = A * d
    call matrix_prod(is,ie,js,je,A,d,z)

    ! r^T*r and d^T*z
    call scalar_prod(is,ie,js,je,res,res,resTres)
    call scalar_prod(is,ie,js,je,d,z,dTz)

    ! alpha = r^T*r / d^T*z
    alpha = resTres / dTz

    ! x_n+1 = x_n + alpha * d
    sol = sol + alpha * d

    ! res_n+1 = res_n - alpha*z
    res = res + alpha * z

    ! res_n+1^T*res_n+1
    call scalar_prod(is,ie,js,je,res,res,resnTresn)

    ! beta = res_n+1^T*res_n+1 / res^T*res  
    beta = resnTresn / resTres

    ! d_n+1 = r_n+1 + beta*d_n
    d = res + beta * d

    if ( resnTresn < crit ) then
    ! stop iteration
      goto 100
    endif

  enddo
100 

end subroutine solve_poisson_qg



subroutine make_matrix(is,ie,js,je,dx,dy,A)
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  real*8, intent(in)                          :: dx, dy
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,,js:je), intent(out):: A
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k

end subroutine make_matrix


subroutine scalar_prod(is,ie,js,je,p1,p2,sprod)
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  real*8, intent(in)                          :: dx, dy
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
  real*8, intent(in)                          :: dx, dy
  real*8, dimension(is:ie,js:je), intent(in)  :: v
  real*8, dimension(is:ie,js:je), intent(out) :: A
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,js:je), intent(out) :: mprod
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j
  mprod = 0.0
  do i=is,ie
    do j=js,je
      mprod(i,j) = mprod  + p1(ii,jj)*p2(ii,jj)
    enddo
  enddo
end subroutine scalar_prod


