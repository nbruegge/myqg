

subroutine solve_poisson_cg(is,ie,js,je,nz,forc,sol,max_itt,crit, est_error, doio, tstep, matA, C)
  !use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je,nz
  !real*8, intent(in)                          :: dx, dy
  real*8, dimension(is:ie,js:je,nz), intent(in)  :: forc
  integer, intent(in)                         :: max_itt       
  real*8, intent(in)                          :: crit
  logical, intent(in)                         :: doio
  integer, intent(in)                         :: tstep
  real*8, dimension(3,3,nz,3), intent(in)     :: matA
  real*8, dimension(is:ie,js:je,nz,3), intent(in) :: C 
!UPDATE PARAMETERS: ======================================== 
  ! INPUT   sol is used as initial guess
  ! OUTPUT: sol is calculated by inversion of forc
  real*8, dimension(is:ie,js:je,nz), intent(inout) :: sol
!OUTPUT PARAMETERS: ======================================== 
  real*8, intent(out)                         :: est_error
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k,n
  real*8, dimension(is:ie,js:je,nz)           :: res, d, z, h
  logical, save                               :: first=.true.
  real*8                                      :: alpha, beta
  real*8                                      :: resTh, resp1Thp1, dTz, resp1Tresp1
  real*8                                      :: forcTforc 
  !real*8, dimension(is:ie,js:je,nz)           :: C
  real*8                                      :: dmax, hmax, beta_min
  real*8                                      :: convergence_rate
  real*8                                      :: maximprovement, maximprovement1
  real*8                                      :: solmax

! This routine solves the Poisson equation with the Preconditioned Conjugate Gradient method (PCG)
!
! -------------------------------------------------------------------------------
! Numenclatura:
! -------------------------------------------------------------------------------
! forc = matA * sol    : Poisson equation
!
! sol:    iterative solution
! matA:   Laplace operator
! forc:   inhomogenity in Poisson equation
! res:    residual of solution: res = forc - matA*sol 
!
! -------------------------------------------------------------------------------
! The Algorithm (after Wikipedia "CG-Verfahren"):
! -------------------------------------------------------------------------------
! Calculated once per iteration and than saved to save computational time:
!   z       = matA*d 
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
! only needed once at the beginning 
! ----------------------------------------
  if ( first) then
! ----------------------------------------
! make coefficient matrix A
! ----------------------------------------
    call make_matrix
! ----------------------------------------
! make preconditioned matrix C
! ----------------------------------------
    call precond_mat(is,ie,js,je)
    first = .false.
  endif
  !call make_matrix(matA)

! ----------------------------------------
! initial guess
! ----------------------------------------
  !call matrixC_prod(is,ie,js,je,nz,C,forc,sol)
  !call cyclic_exchange_2d(sol)
  sol = 0.0

! ----------------------------------------
! r0 = b - A*x0
! ----------------------------------------
  call matrix_prod(is,ie,js,je,nz,matA,sol,res)
  !call cyclic_exchange_2d(res)  !?? no cyclic
  res = forc - res
  !write(*,*) res(:,:,2)

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
    goto 101    ! success
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
    call matrix_prod(is,ie,js,je,nz,matA,d,z)
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
        write(*,*) "Warning: Solver is diverging at n = ", n
        !goto 99
      endif
    endif

    ! ----------------------------------------
    ! d_k+1 = h_k+1 + beta_k*d_k
    ! ----------------------------------------
    d = h + beta * d
    !?? call cyclic_exchange_2d(d)

    ! ----------------------------------------
    ! calculate convergence
    ! ----------------------------------------
    call absmax_element(is,ie,js,je,nz, d, dmax)
    maximprovement = abs(alpha) * dmax

    if ( n == 1 ) then
      maximprovement1 = maximprovement
      est_error = maximprovement
      if ( est_error .lt. crit )  goto 101    !success
    else
      convergence_rate = exp(log(maximprovement/maximprovement1)/(n-1))
      est_error = maximprovement*convergence_rate/(1.0-convergence_rate)
      call scalar_prod(is,ie,js,je,nz,res,res,resp1Tresp1)
      call scalar_prod(is,ie,js,je,nz,forc,forc,forcTforc)
      ! ----------------------------------------
      ! monitor solver
      ! ----------------------------------------
      if ( doio ) then
        write(*,*) "n = ", n, "; convergence_rate = ", sngl(convergence_rate), &
                 & "; est_error = ", sngl(est_error), "; abs(resTres/forcTforc) = ", resp1Tresp1/forcTforc
      endif

      ! ----------------------------------------
      ! if converged break out
      ! ----------------------------------------
      !if ( est_error .lt. crit )  goto 101    !success
      if ( resp1Tresp1/forcTforc .lt. 1.e-10 )  goto 101    !success

      !! temp solution to survive first step
      !if ( tstep == 2 .and. abs(sngl(est_error)) < 10. ) goto 101
    endif
      
  enddo ! end of iteration

  ! not converged after max_itt steps 
  !goto 99    ! error
  goto 101

! error
99  write(*,*) ":: There is no convergence! ::"
  stop

! success
!101 write(*,*) "Converged after n = ", n 
101 return 

end subroutine solve_poisson_cg

subroutine make_matrix
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  !integer, intent(in)                         :: is,ie,js,je
  !real*8, intent(in)                          :: dx, dy
!OUTPUT PARAMETERS: ======================================== 
  !real*8, dimension(is:ie,js:je,3,3), intent(out):: matA
  !real*8, dimension(3,3,nz,3), intent(out):: matA
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k
  integer                                     :: ii,jj,kk

  !do j=js,je
  !  do i=is,ie
      !matA(i,j,1,2) =  1.0/dx**2
      !matA(i,j,3,2) =  1.0/dx**2
      !matA(i,j,2,2) = -2.0/dx**2-2.0/dy**2
      !matA(i,j,2,1) =  1.0/dy**2
      !matA(i,j,2,3) =  1.0/dy**2
  !  enddo
  !enddo
  
  matA=0.0
  do k=1,nz
    matA(1,2,k,2)   = matA(1,2,k,2)   + 1.0/dx**2
    matA(2,2,k,2)   = matA(2,2,k,2)   - 2.0/dx**2
    matA(3,2,k,2)   = matA(3,2,k,2)   + 1.0/dx**2

    matA(2,1,k,2)   = matA(2,1,k,2)   + 1.0/dy**2
    matA(2,2,k,2)   = matA(2,2,k,2)   - 2.0/dy**2
    matA(2,3,k,2)   = matA(2,3,k,2)   + 1.0/dy**2

    if ( k > 1 .and. k < nz ) then
      matA(2,2,k,1) = matA(2,2,k,1) + f0**2/Hk(k)*(1.0/gred(k))
      matA(2,2,k,2) = matA(2,2,k,2) - f0**2/Hk(k)*(1.0/gred(k)+1.0/gred(k+1)) 
      matA(2,2,k,3) = matA(2,2,k,3) + f0**2/Hk(k)*(1.0/gred(k+1))
    elseif ( k == 1 ) then
      matA(2,2,k,2) = matA(2,2,k,2) - f0**2/Hk(k)*(1.0/gred(k)+1.0/gred(k+1)) 
      matA(2,2,k,3) = matA(2,2,k,3) + f0**2/Hk(k)*(1.0/gred(k+1))
      !matA(2,2,k,2) = matA(2,2,k,2) - f0**2/500.*(1.0/gred(k)+1.0/gred(k+1)) 
      !matA(2,2,k,3) = matA(2,2,k,3) + f0**2/500.*(1.0/gred(k+1))
    elseif ( k == nz ) then
      matA(2,2,k,1) = matA(2,2,k,1) + f0**2/Hk(k)*(1.0/gred(k))
      matA(2,2,k,2) = matA(2,2,k,2) - f0**2/Hk(k)*(1.0/gred(k)) 
      !matA(2,2,k,1) = matA(2,2,k,1) + f0**2/500.*(1.0/gred(k))
      !matA(2,2,k,2) = matA(2,2,k,2) - f0**2/500.*(1.0/gred(k)) 
    endif
    
    !!if ( k > 1 .and. k < nz ) then
    !!  matA(2,2,k-1,0) = matA(2,2,k-1,0) + f0**2/Hk(k)*(1.0/gred(k))
    !!  matA(2,2,k,0)   = matA(2,2,k,0)   - f0**2/Hk(k)*(1.0/gred(k)+1.0/gred(k+1)) 
    !!  matA(2,2,k+1,0) = matA(2,2,k+1,0) + f0**2/Hk(k)*(1.0/gred(k+1))
    !!elseif ( k == 1 ) then
    !!  !matA(2,2,k-1,0) = matA(2,2,k-1,0) + f0**2/Hk(k)*(1.0/gred(k))
    !!  matA(2,2,k,0)   = matA(2,2,k,0)   - f0**2/Hk(k)*(1.0/gred(k)+1.0/gred(k+1)) 
    !!  matA(2,2,k+1,0) = matA(2,2,k+1,0) + f0**2/Hk(k)*(1.0/gred(k+1))
    !!elseif ( k == nz ) then
    !!  matA(2,2,k-1,0) = matA(2,2,k-1,0) + f0**2/Hk(k)*(1.0/gred(k))
    !!  matA(2,2,k,0)   = matA(2,2,k,0)   - f0**2/Hk(k)*(1.0/gred(k)) 
    !!  !matA(2,2,k+1,0) = matA(2,2,k+1,0) + f0**2/Hk(k)*(1.0/gred(k+1))
    !!endif
    !if ( k > 1 .and. k < nz ) then
    !  matA(2,2,k-1) = matA(2,2,k-1) + f0**2/Hk(k)*(1.0/gred(k-1))
    !  matA(2,2,k)   = matA(2,2,k)   - f0**2/Hk(k)*(1.0/gred(k-1)+1.0/gred(k)) 
    !  matA(2,2,k+1) = matA(2,2,k+1) + f0**2/Hk(k)*(1.0/gred(k))
    !elseif ( k == 1 ) then
    !  !matA(2,2,k-1) = matA(2,2,k-1) + f0**2/Hk(k)*(1.0/gred(k-1))
    !  matA(2,2,k)   = matA(2,2,k)   - f0**2/Hk(k)*(1.0/grav+1.0/gred(k)) 
    !  matA(2,2,k+1) = matA(2,2,k+1) + f0**2/Hk(k)*(1.0/gred(k))
    !elseif ( k == nz ) then
    !  matA(2,2,k-1) = matA(2,2,k-1) + f0**2/Hk(k)*(1.0/gred(k-1))
    !  matA(2,2,k)   = matA(2,2,k)   - f0**2/Hk(k)*(1.0/gred(k-1)) 
    !  !matA(2,2,k+1) = matA(2,2,k+1) + f0**2/Hk(k)*(1.0/gred(k+1))
    !endif
  enddo

  !!!!!do k=1,nz
  !!!!!  do ii=1,3
  !!!!!    do jj=1,3
  !!!!!      do kk=1,3
  !!!!!        write(*,*) 'ii, jj, k, kk ', ii, jj, k, kk
  !!!!!        write(*,*) matA(ii,jj,k,kk)
  !!!!!      enddo
  !!!!!    enddo
  !!!!!  enddo
  !!!!!enddo
  !!!!!  write(*,*) '_____________________________________________________'
  !!!!!  write(*,*) 'tstep = ', tstep
  !!!!!  !write(*,*) pvr(:,:,1)
  !!!!!  write(*,*) '_____________________________________________________'
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

subroutine matrix_prod(is,ie,js,je,nz,matA,d,z)
  !use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je,nz
  real*8, dimension(is:ie,js:je,nz), intent(in)  :: d
  !real*8, dimension(is:ie,js:je,3,3), intent(out) :: matA
  real*8, dimension(3,3,nz,3), intent(in) :: matA
!OUTPUT PARAMETERS: ======================================== 
  real*8, dimension(is:ie,js:je,nz), intent(out) :: z
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i, j, k, ii, jj, kk
  integer                                     :: k1, k2
  z = 0.0
  !do j=js+1,je-1
  !  do i=is+1,ie-1
  !    do k=1,nz
  !      do ii = -1,1
  !        do jj = -1,1
  !          !z(i,j) = z(i,j) + matA(i,j,ii+2,jj+2)*d(i+ii,j+jj)
  !          z(i,j,k) = z(i,j,k) + matA(ii+2,jj+2,k)*d(i+ii,j+jj,k)
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

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
          do jj = -1,1
            do kk = k1,k2
              !if ( i==10 .and. j==10 ) then
              !  write(*,*) 'i,j,k,ii,jj,kk', i,j,k,ii,jj,kk
              !  write(*,*) matA(ii+2,jj+2,k,kk+2)*d(i+ii,j+jj,k+kk)
              !endif
              z(i,j,k) = z(i,j,k) + matA(ii+2,jj+2,k,kk+2)*d(i+ii,j+jj,k+kk)
              !z(i,j,k) = z(i,j,k) + matA(ii+2,jj+2,k,2)*d(i+ii,j+jj,k)
            enddo
          enddo
        enddo
        !do kk = k1,k2
        !  z(i,j,k) = z(i,j,k) + matA(2,2,k,kk+2)*d(i,j,k+kk)
        !enddo
      enddo
    enddo
  enddo
end subroutine matrix_prod

subroutine precond_mat(is,ie,js,je)
  use myqg_module
  implicit none
!INPUT PARAMETERS:  ======================================== 
  integer, intent(in)                         :: is,ie,js,je
  !integer, intent(in)                         :: nz
  !real*8, dimension(is:ie,js:je,3,3), intent(in):: matA
  !real*8, dimension(3,3,nz,3), intent(in):: matA
!OUTPUT PARAMETERS: ======================================== 
  !real*8, dimension(is:ie,js:je,nz), intent(out):: C
!LOCAL VARIABLES:   ======================================== 
  integer                                     :: i,j,k
  C=0.0
  ! copy diagonal elements
  do j=js,je
    do i=is,ie
      do k=1,nz
        !C(i,j) = matA(i,j,2,2)
        C(i,j,k) = matA(2,2,k,2)
      enddo
    enddo
  enddo
  ! invert matA
  do j=js,je
    do i=is,ie
      do k=1,nz
        if ( C(i,j,k) .ne. 0.0 ) then
          C(i,j,k) = 1./C(i,j,k)
        !else
        !  C(i,j,k) = 0.0
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
