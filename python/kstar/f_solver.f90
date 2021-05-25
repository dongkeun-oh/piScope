!---------------------------------------------------------------------------
!  f_solver.f90 !
!  Solver of elliptic equation with irregular boundary.
!
!  created by Insik Choi
!  2015. 3
!----------------------------------------------------------------------------

	  SUBROUTINE f_solver(Rbound,Zbound,N,idxmap,ij2dof,SR,PSI)
	   use input_params

	    implicit none
            real(kind=8), intent(in) :: RBound(0:NR,0:NZ),ZBound(0:NR,0:NZ)
	    integer, intent(in) :: N
            integer, intent(in) :: idxmap(1:NR*NZ,1:3), ij2dof(1:NR-1,1:NZ-1)
            real(kind=8), intent(in) :: SR(0:NR,0:NZ)
	    real(kind=8), intent(inout) :: PSI(0:NR,0:NZ)

	    real(kind=8) :: DR2, DR2DZ2, temp

	    integer :: i, j, k, flag

! Linear algebraic equation Variables  
!A system of linear equations with a real general matrix
!(Crout's method)
            real(kind=8),dimension(1:N,1:N) :: M_A
            real(kind=8),dimension(1:N) :: M_B, VW
            integer, dimension(1:N) :: IP
            real(kind=8) :: EPSZ
            integer :: ISW, IS, ICON

	    DR2 = DR * DR
            DR2DZ2 = DR2/Dz/Dz
 
            M_A = 0.0

            do k = 1, N
               i = idxmap(k, 1)
               j = idxmap(k, 2)
               flag = idxmap(k, 3)
               if (flag == 0) then
                  temp = 1./(1. + 0.5d0*DR/R(i))
                  M_A(k,k) = - temp
                  M_A(k,ij2dof(i+1,j)) = temp
                  temp = 1./(1. - 0.5d0*DR/R(i))
                  M_A(k,ij2dof(i-1,j)) = temp
                  M_A(k,k) = M_A(k,k) - temp - 2.*DR2DZ2
                  M_A(k,ij2dof(i,j-1)) = DR2DZ2
                  M_A(k,ij2dof(i,j+1)) = DR2DZ2

                  M_B(k) = DR2*SR(i,j)
               else
                  select case(flag)
                     case (1)
                        temp = Rbound(i-1,j)
                        M_A(k,k) = temp
                        M_A(k,ij2dof(i-1,j)) = 1. - temp
                     case (2)
                        temp = Rbound(i,j)
                        M_A(k,k) = temp
                        M_A(k,ij2dof(i+1,j)) = 1. - temp
                     case (3)
                        temp = Zbound(i,j-1)
                        M_A(k,k) = temp
                        M_A(k,ij2dof(i,j-1)) = 1. - temp
                     case (4)
                        temp = Zbound(i,j)
                        M_A(k,k) = temp
                        M_A(k,ij2dof(i,j+1)) = 1. - temp
                  endselect

                  M_B(k) = 0.0
               endif
            enddo

            EPSZ = 1.0E-5
            ISW = 1
!            CALL DLAX(M_A(1:N,1:N), N, N, M_B(1:N), EPSZ, ISW, IS, VW, IP, ICON)
            CALL DGESV(N,1,M_A,N,IP,M_B,N,ICON)

            WRITE(6,620) ICON
            if (ICON /=0) stop

            PSI = 0.0
            do k = 1, N
               i = idxmap(k, 1)
               j = idxmap(k, 2)
               PSI(i,j) = M_B(k)
            enddo

620 FORMAT('0',10X,'CONDITION CODE=',I5)
630 FORMAT((2X,45(F6.2)))
	 
	  RETURN
	  END SUBROUTINE f_solver

	  SUBROUTINE set_coef(RBound,ZBound,N,idxmap,ij2dof)
	   use input_params

            implicit none
            real(kind=8), intent(in) :: RBound(0:NR,0:NZ),ZBound(0:NR,0:NZ)
            integer, intent(out) :: N
            integer, intent(out) :: idxmap(1:NR*NZ,1:3)  
            integer, intent(out) :: ij2dof(1:NR-1,1:NZ-1)
            real(kind=8) :: temp
            integer :: i, j, k, ki, nnn, flag

            N = 0
            ij2dof = -1
        
            do i = 1, NR-1
               do j = 1, NZ-1
                  nnn = 0
                  temp = 0.0
                  if (Rbound(i-1,j) > 0.) then
                     nnn = nnn + 1
                     temp = Rbound(i-1,j)
                     flag = 1
                  endif
                  if (Rbound(i,j) > 0.) then
                     nnn = nnn + 1
                     if (temp < Rbound(i,j)) then
                        temp = Rbound(i,j)
                        flag = 2
                     endif
                  endif
                  if (Zbound(i,j-1) > 0.) then
                     nnn = nnn + 1
                     if (temp < Zbound(i,j-1)) then
                        temp = Rbound(i,j-1)
                        flag = 3
                     endif
                  endif
                  if (Zbound(i,j) > 0.) then
                     nnn = nnn + 1
                     if (temp < Zbound(i,j)) then
                        temp = Zbound(i,j)
                        flag = 4
                     endif
                  endif
                  if (nnn > 0) then
                     N = N + 1
                     idxmap(N, 1) = i
                     idxmap(N, 2) = j
                     if (nnn == 4) then 
                        idxmap(N, 3) = 0
                     else
                        idxmap(N, 3) = flag
                     endif
                     ij2dof(i,j) = N
                  else
                     ij2dof(i,j) = -1
                  endif
               enddo
            end do

	  RETURN
	  END SUBROUTINE set_coef
