

subroutine solve_poisson_cg(is,ie,js,je,nz,forc,sol,max_itt,crit, est_error)
  !use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je,nz
  !real*8, intent(in)                          :: dx, dy
  real*8, dimension(is:ie,js:je,nz), intent(in)  :: forc
  integer, intent(in)                         :: max_itt       
  real*8, intent(in)                          :: crit
!UPDATE PARAMETERS: ======================================== 
  ! INPUT   sol is used as initial guess
  ! OUTPUT: sol is calculated by inversion of forc
  real*8, dimension(is:ie,js:je,nz), intent(inout) :: sol
!OUTPUT PARAMETERS: ======================================== 
  real*8, intent(out)                         :: est_error
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k,n
  !real*8, dimension(is:ie,js:je,3,3)          :: A
  real*8, dimension(3,3,nz)                   :: A
  real*8, dimension(is:ie,js:je,nz)           :: res, d, z, h
  logical, save                               :: first=.true.
  real*8                                      :: alpha, beta
  real*8                                      :: resTh, resp1Thp1, dTz, resp1Tresp1
  real*8, dimension(is:ie,js:je,nz)           :: C
  real*8                                      :: dmax, hmax, beta_min
  real*8                                      :: convergence_rate
  real*8                                      :: maximprovement, maximprovement1
  real*8                                      :: solmax

! This routine solves the Poisson equation with the Preconditioned Conjugate Gradient method (PCG)
!
! -------------------------------------------------------------------------------
! Numenclatura:
! -------------------------------------------------------------------------------
! forc = A * sol    : Poisson equation
!
! sol:    iterative solution
! A:      Laplace operator
! forc:   inhomogenity in Poisson equation
! res:    residual of solution: res = forc - A*sol 
!
! -------------------------------------------------------------------------------
! The Algorithm (after Wikipedia "CG-Verfahren"):
! -------------------------------------------------------------------------------
! Calculated once per iteration and than saved to save computational time:
!   z       = A*d 
! Find new solution:
!   alpha   = res*h / d*z
!   sol_new = sol + alpha * d
! Update residuum:
!   res_new = res - alpha * z
! Update h:
!   h_new   = C * res_new
! Update search direction:
!   beta    = r_new*h_new / r*h
!   d_new   = h_new + beta * d

! ----------------------------------------
! make coefficient matrix A
! ----------------------------------------
  !if ( first) then
  !  call make_matrix(is,ie,js,je,A)
  !  first = .false.
  !endif
  call make_matrix(is,ie,js,je,A)

! ----------------------------------------
! make preconditioned matrix C
! ----------------------------------------
  call precond_mat(is,ie,js,je,nz,A,C)

! ----------------------------------------
! initial guess
! ----------------------------------------
  !call matrixC_prod(is,ie,js,je,nz,C,forc,sol)
  !call cyclic_exchange_2d(sol)
  sol = 0.0

! ----------------------------------------
! r0 = b - A*x0
! ----------------------------------------
  call matrix_prod(is,ie,js,je,nz,A,sol,res)
  !call cyclic_exchange_2d(res)  !?? no cyclic
  res = forc - res

! ----------------------------------------
! h0 = C*r0
! ----------------------------------------
  call matrixC_prod(is,ie,js,je,nz,C,res,h)

! security check (if h=0, we divide by 0 when calculating beta) 
  n = 0
  call absmax_element(is,ie,js,je,nz, h, hmax)
  if ( 100.0*hmax .lt. crit ) then
    est_error = 100.0*hmax
    write(*,*) "Converged before iteration"
    !goto 101    ! success
  endif

! ----------------------------------------
! d0 = h0
! ----------------------------------------
  d = h

! ----------------------------------------
! enter main loop
! ----------------------------------------
  do n=1,max_itt
    ! ----------------------------------------
    ! z = A * d_k
    ! ----------------------------------------
    call matrix_prod(is,ie,js,je,nz,A,d,z)
    !call cyclic_exchange_2d(z) !?? not necessary

    ! ----------------------------------------
    ! r^T*h and d^T*z
    ! ----------------------------------------
    call scalar_prod(is,ie,js,je,nz,res,h,resTh)
    call scalar_prod(is,ie,js,je,nz,d,z,dTz)

    ! security check (if d=0, we divide by 0 when calculating alpha)
    if ( abs(dTz) .lt. abs(beta)*1e-10 ) then
      call absmax_element(is,ie,js,je,nz, d, dmax)
      est_error = 100.0 * dmax
      !goto 101    ! success
    endif

    ! ----------------------------------------
    ! alpha_k = r_k^T*h_k / d_k^T*z
    ! ----------------------------------------
    alpha = resTh / dTz

    ! ----------------------------------------
    ! x_k+1 = x_k + alpha_k * d_k
    ! ----------------------------------------
    sol = sol + alpha * d

    ! ----------------------------------------
    ! res_k+1 = res_k - alpha_k * z
    ! ----------------------------------------
    res = res - alpha * z

    ! ----------------------------------------
    ! h_k+1 = C * r_k+1
    ! ----------------------------------------
    call matrixC_prod(is,ie,js,je,nz,C,res,h)
    !call cyclic_exchange_2d(h)

    ! ----------------------------------------
    ! res_k+1^T*h_h+1
    ! ----------------------------------------
    call scalar_prod(is,ie,js,je,nz,res,h,resp1Thp1)

    ! ----------------------------------------
    ! beta = res_k+1^T*h_k+1 / res_k^T*h_k  
    ! ----------------------------------------
    beta = resp1Thp1 / resTh

    !if ( n==10 ) then
    !  write(*,*) "alpha = ", alpha
    !  write(*,*) "beta  = ", beta
    !  write(*,*) sngl(res)
    !  stop
    !endif

    ! security check
    if ( n==1 ) then
      beta_min = abs(beta)
    else
      beta_min = min(beta_min,abs(beta))
      if ( abs(beta) .gt. 100.0*beta_min ) then
        !write(*,*) "Warning: Solver is diverging at itt = ??"
      endif
    endif

    ! ----------------------------------------
    ! d_k+1 = h_k+1 + beta_k*d_k
    ! ----------------------------------------
    d = h + beta * d
    !?? call cyclic_exchange_2d(d)

    ! ----------------------------------------
    ! test for convergence
    ! ----------------------------------------
    call absmax_element(is,ie,js,je,nz, d, dmax)
    maximprovement = abs(alpha) * dmax
    if ( n == 1 ) then
      maximprovement1 = maximprovement
      est_error = maximprovement
      write(*,*) "converged at n=1"
      if ( est_error .lt. crit )  goto 101    !success
    else
      convergence_rate = exp(log(maximprovement/maximprovement1)/(n-1))
      est_error = maximprovement*convergence_rate/(1.0-convergence_rate)
      call scalar_prod(is,ie,js,je,nz,res,res,resp1Tresp1)
      write(*,*) "n = ", n, "; convergence_rate = ", sngl(convergence_rate), &
               & "; est_error = ", sngl(est_error), "; abs(resp1) = ", resp1Tresp1
      call absmax_element(is,ie,js,je,nz, sol, solmax)
      write(*,*) "absmax_element(sol) = ", solmax
      !if ( est_error .lt. crit )  goto 101    !success
    endif
      
  enddo ! end of iteration

  ! not converged after max_itt steps 
  !goto 99    ! error
  goto 101

! error
99  write(*,*) ":: There is no convergence! ::"
  stop

! success
101 write(*,*) "Converged after n = ", n 

end subroutine solve_poisson_cg

subroutine make_matrix(is,ie,js,je,A)
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  !real*8, intent(in)                          :: dx, dy
!OUTPUT PARAMETERS: ======================================== 
  !real*8, dimension(is:ie,js:je,3,3), intent(out):: A
  real*8, dimension(3,3,nz), intent(out):: A
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k
  A = 0.0
  !do j=js,je
  !  do i=is,ie
      !A(i,j,1,2) =  1.0/dx**2
      !A(i,j,3,2) =  1.0/dx**2
      !A(i,j,2,2) = -2.0/dx**2-2.0/dy**2
      !A(i,j,2,1) =  1.0/dy**2
      !A(i,j,2,3) =  1.0/dy**2
  !  enddo
  !enddo

  do k=1,nz
    A(1,2,k)   = A(1,2,k)   + 1.0/dx**2
    A(2,2,k)   = A(2,2,k)   - 2.0/dx**2
    A(3,2,k)   = A(3,2,k)   + 1.0/dx**2

    A(2,1,k)   = A(2,1,k)   + 1.0/dy**2
    A(2,2,k)   = A(2,2,k)   - 2.0/dy**2
    A(2,3,k)   = A(2,3,k)   + 1.0/dy**2

    if ( k > 1 .and. k < nz ) then
      A(2,2,k-1) = A(2,2,k-1) + f0**2/Hk(k)*(1.0/gred(k))
      A(2,2,k)   = A(2,2,k)   - f0**2/Hk(k)*(1.0/gred(k)+1.0/gred(k+1)) 
      A(2,2,k+1) = A(2,2,k+1) + f0**2/Hk(k)*(1.0/gred(k+1))
    elseif ( k == 1 ) then
      !A(2,2,k-1) = A(2,2,k-1) + f0**2/Hk(k)*(1.0/gred(k))
      A(2,2,k)   = A(2,2,k)   - f0**2/Hk(k)*(1.0/gred(k)+1.0/gred(k+1)) 
      A(2,2,k+1) = A(2,2,k+1) + f0**2/Hk(k)*(1.0/gred(k+1))
    elseif ( k == nz ) then
      A(2,2,k-1) = A(2,2,k-1) + f0**2/Hk(k)*(1.0/gred(k))
      A(2,2,k)   = A(2,2,k)   - f0**2/Hk(k)*(1.0/gred(k)) 
      !A(2,2,k+1) = A(2,2,k+1) + f0**2/Hk(k)*(1.0/gred(k+1))
    endif
  enddo

end subroutine make_matrix

subroutine scalar_prod(is,ie,js,je,nz,p1,p2,sprod)
  !use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je,nz
  real*8, dimension(is:ie,js:je,nz), intent(in)  :: p1, p2
!OUTPUT PARAMETERS: ======================================== 
  real*8, intent(out)                         :: sprod
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k
  sprod = 0.0
  do i=is,ie
    do j=js,je
      do k=1,nz
        sprod = sprod  + p1(i,j,k)*p2(i,j,k)
      enddo
    enddo
  enddo
end subroutine scalar_prod

subroutine matrix_prod(is,ie,js,je,nz,A,d,z)
  !use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je,nz
  real*8, dimension(is:ie,js:je,nz), intent(in)  :: d
  !real*8, dimension(is:ie,js:je,3,3), intent(out) :: A
  real*8, dimension(3,3,nz), intent(out) :: A
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,js:je,nz), intent(out) :: z
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i, j, k, ii, jj
  z = 0.0
  do j=js+1,je-1
    do i=is+1,ie-1
      do k=1,nz
        do ii = -1,1
          do jj = -1,1
            !z(i,j) = z(i,j) + A(i,j,ii+2,jj+2)*d(i+ii,j+jj)
            z(i,j,k) = z(i,j,k) + A(ii+2,jj+2,k)*d(i+ii,j+jj,k)
          enddo
        enddo
      enddo
    enddo
  enddo
end subroutine matrix_prod

subroutine precond_mat(is,ie,js,je,nz,A,C)
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je,nz
  !real*8, dimension(is:ie,js:je,3,3), intent(in):: A
  real*8, dimension(3,3,nz), intent(in):: A
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,js:je,nz), intent(out):: C
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k
  C=0.0
  ! copy diagonal elements
  do j=js,je
    do i=is,ie
      do k=1,nz
        !C(i,j) = A(i,j,2,2)
        C(i,j,k) = A(2,2,k)
      enddo
    enddo
  enddo
  ! invert A
  do j=js,je
    do i=is,ie
      do k=1,nz
        if ( C(i,j,k) .ne. 0.0 ) then
          C(i,j,k) = 1/C(i,j,k)
        else
          C(i,j,k) = 0.0
        endif
      enddo
    enddo
  enddo
end subroutine precond_mat

subroutine matrixC_prod(is,ie,js,je,nz,C,res,h)
  !use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                             :: is,ie,js,je,nz
  real*8, dimension(is:ie,js:je,nz), intent(in)   :: C
  real*8, dimension(is:ie,js:je,nz), intent(in)   :: res
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,js:je,nz), intent(out)  :: h
!LOCAL VARIABLES:   ======================================== 
  integer                                         :: i,j,k
  h=0.0
  do j=js,je
    do i=is,ie
      do k=1,nz
        h(i,j,k) = C(i,j,k)*res(i,j,k)
      enddo
    enddo
  enddo
end subroutine matrixC_prod

subroutine absmax_element(is,ie,js,je,nz, d, dmax)
  !use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                             :: is,ie,js,je,nz
  real*8, dimension(is:ie,js:je,nz), intent(in)   :: d 
!OUTPUT PARAMETERS: ======================================== 
  real*8, intent(out)                             :: dmax
!LOCAL VARIABLES:   ======================================== 
  integer                                         :: i,j,k
  dmax = 0.0
  do j=js,je
    do i=is,ie
      do k=1,nz
        dmax = max(abs(d(i,j,k)),dmax)
      enddo
    enddo
  enddo
end subroutine absmax_element
