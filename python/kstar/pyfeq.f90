      SUBROUTINE run
      USE params
      USE output

      IMPLICIT NONE

      TYPE extremum_info
          REAL(kind=8) :: X, Y, PSI
      END TYPE extremum_info

! Constant parameters
      REAL(kind=8), PARAMETER :: pi=3.1415926535897932385,  &
          twopi=2*pi,pisq=pi*pi

! MKS Units
      REAL(kind=8), PARAMETER :: mu0=4.0d0*pi*1.0d-7
      REAL(kind=8), PARAMETER :: epsilon0=8.8541878d-12
      REAL(kind=8), PARAMETER :: lc=2.99792458d8

      real(kind=8),dimension(:,:),allocatable :: J_region, SR_new, SR_old
      real(kind=8), dimension(:,:), allocatable :: PSI_old
      real(kind=8), dimension(:,:), allocatable :: PSIbar, PSIbar_old
      real(kind=8),dimension(:,:),allocatable :: RBoundary, ZBoundary
      integer, dimension(:,:), allocatable :: idxmap, ij2dof
! set coefs of Grad-Shfranov eq.
      integer :: N
      
      real(kind=8) :: alpha_g
      real(kind=8) :: grid_area, pl_area, chk_current, RI2
      real(kind=8) :: sumpp, sumgp, sumgp1, delta_psi
      real(kind=8) :: chk_val, chk_val1, diff_R, diff_Z
      real(kind=8) :: rb, zb
      
      type(extremum_info) :: PSI_l, PSI_ext, Xpoint
      
      integer :: i,j, k, i1, nfound, nbfound
      logical :: is_boundary

      logical :: is_done = .false.

      INTERFACE
         REAL(kind=8) FUNCTION flux_R(R0, R_shift, s_radius,delta, theta)
           IMPLICIT NONE
           real(kind=8), intent(in) ::  R0, R_shift, s_radius,delta, theta
         END FUNCTION flux_R
         REAL(kind=8) FUNCTION flux_Z(Z0, Z_shift, ellipticity, s_radius, theta)
           IMPLICIT NONE
           real(kind=8), intent(in) ::  Z0, Z_shift, ellipticity, s_radius, theta
         END FUNCTION flux_Z
         INTEGER FUNCTION InsidePlasma(x1, y1,ellip_R,ellip_Z, nfound)
           use input_params, only : pi, N_divide
           IMPLICIT NONE
           real(kind=8), intent(in) ::  x1,y1,ellip_R(N_divide),ellip_Z(N_divide)
           integer, intent(in) ::  nfound
         END FUNCTION InsidePlasma
      END INTERFACE

      if (is_set) then
         if (allocated(PSI_new)) deallocate(PSI_new)
         
         allocate(RBoundary(0:NR,0:NZ),ZBoundary(0:NR,0:NZ),J_region(0:NR,0:NZ), &
              SR_new(0:NR,0:NZ), SR_old(0:NR,0:NZ), PSI_new(0:NR,0:NZ), &
              PSI_old(0:NR,0:NZ),PSIbar(0:NR,0:NZ), PSIbar_old(0:NR,0:NZ), &
              idxmap(1:NR*NZ,1:3), ij2dof(1:NR-1,1:NZ-1))
         
         PSI_l%PSI = 0.0
!         print*,'Calling Set Boundary ...'
         CALL set_boundary(npoints, pl_R, pl_Z,RBoundary, ZBoundary, nfound)
!         print*,'Calling Set Coeffient ...'
         CALL set_coef(RBoundary,ZBoundary,N,idxmap,ij2dof)
!         
         CALL init()
         
         PSI_old = 0.0
         PSIbar=0.0
         PSIbar_old=0.0
         is_done = .false.
         is_good = .false.

         do i1 = 1, max_loop
            print*,' '
            print*, 'Check ii iteration = ',i1
! ------------------------------------------------------------------------------------------------------
            print*,'Calling elliptic_solver ...'
           CALL f_solver(Rboundary,Zboundary,N,idxmap,ij2dof,SR_new,PSI_new)
! ----------------------------------------------------------------------------------------
            chk_val = 0.0d0
            do i = 1, NR-1
               do j = 1, NZ-1
!  Compare old and new values of PSI(i,j) in the plasma region
                  chk_val1 = DABS(PSI_new(i,j) - PSI_old(i,j))
                  if (chk_val1 > chk_val) then
                     chk_val = chk_val1
                     diff_R = R(i)
                     diff_Z = Z(j)
                  end if
               end do
            end do
! ----------------------------------------------------------------------------------------
            print*, ' '
            print*,'i1 = ',i1, ' chk_value =  ',chk_val,' crit2 = ',crit2
            if (chk_val < crit2 ) then
               print*,'Calling Critical values ...'
               CALL critical(N,idxmap,ij2dof,R0,Z0,PSI_new,PSI_ext,Xpoint)
        
               nbfound = 0
               do i = 0, NR
                  do j = 0, NZ
                     is_boundary = .false.
                     if ((Rboundary(i,j) < 1.0) .and. (Rboundary(i+1,j) == 1.0) ) then
                        rb = R(i) + (1. - Rboundary(i,j))*DR
                        zb = Z(j)
                        is_boundary = .true.
                     else if (i > 0 .and. (Rboundary(i,j) == 0.0) .and. (Rboundary(i-1,j) > 0.)) then
                        rb = R(i) - (1. - Rboundary(i-1,j))*DR
                        zb = Z(j)
                        is_boundary = .true.
                     else if ((Zboundary(i,j) < 1.0) .and. (Zboundary(i,j+1) == 1.0)) then
                        rb = R(i)
                        zb = Z(j) + (1. - Zboundary(i,j))*DZ
                        is_boundary = .true.
                     else if (j > 0 .and. (Zboundary(i,j) == 0.0) .and. (Zboundary(i,j-1) > 0.0)) then
                        rb = R(i)
                        zb = Z(j) - (1. - Zboundary(i,j-1))*DZ
                        is_boundary = .true.
                     endif
                     if (is_boundary) then
                        nbfound = nbfound + 1
                        RBout(nbfound) = rb
                        ZBout(nbfound) = zb
                     endif
                     print*, R(i), Z(j), PSI_new(i,j)
                  end do
               end do
               
              print*,' '
              print*,'ADJUSTED BOUNDARY'
               call sort_order(RBout, ZBout, nbfound, 1.8d0, 0.d0)
               do i = 1, nbfound
                  print*,RBout(i),ZBout(i)
               end do
               is_done = .true.
               is_good = .true.
            end if
! set old as new values
            PSI_old = PSI_new
!---------------------------------------------------------------------
            if (i1 == max_loop) then
               print*, 'Fail to Simulate Grad-Shafranov Eq. within i1 iteration ',max_loop
               print*, 'chk_value =  ',chk_val
               ! Check PSI values at the Simulation Fail
               is_done = .true.
               is_good = .false.
            end if
!---------------------------------------------------------------------
            
            if (is_done) exit

            ! Mixing old and new psi as next psi for eliminating odd-even limit cycle

            if (i1 /= 1) then
               do i = 1, NR-1
                  do j = 1, NZ-1
                     PSI_new(i,j) = (1.0 - mixing)*PSI_new(i,j) +  mixing*PSI_old(i,j)
                  end do
               end do
            end if

!---------------- TEST ----------------------
!              do i = 0, NR
!                do j = 0, NZ
!                 print*,R(i),Z(j),PSI_new(i,j)
!                end do
!              end do
!---------------- TEST ----------------------

            print*,'Calling Critical values ...'
            CALL critical(N,idxmap,ij2dof,R0,Z0,PSI_new,PSI_ext,Xpoint)

            delta_psi = PSI_l%PSI - PSI_ext%PSI

            if (pl_flag == 1) then
               sumpp = 0.0
               sumgp = 0.0
               sumgp1 = 0.0
               do i = 1, NR-1
                  do j = 1, NZ-1
                     SR_old(i,j) = SR_new(i,j)
                     PSIbar(i,j) =  (PSI_l%PSI - PSI_new(i,j))/delta_psi 
                     do k = 1, NP
                        sumpp = sumpp + alpha(k)*R(i)*PSIbar(i,j)**k
                     end do
                     sumgp = sumgp + PSIbar(i,j)/R(i)
                     if (NF > 1) then
                        do k = 2, NF
                           sumgp1 = sumgp1 + beta(k)*(PSIbar(i,j)**k)/R(i)
                        end do
                     end if
                  end do
               end do
               alpha_g = mu0*pl_cur/grid_area + sumpp + sumgp1
               alpha_g = -alpha_g/sumgp
!    Source term=mu0*R*J
               chk_current = 0.0
               do i = 1, NR-1
                  RI2 = R(i)*R(i)
                  do j = 1, NZ-1
                     SR_new(i,j) = 0.0 
                     do k = 1, NP
                        SR_new(i,j) = SR_new(i,j) - RI2*alpha(k)*PSIbar(i,j)**k
                     end do
                     SR_new(i,j) = SR_new(i,j) - alpha_g*PSIbar(i,j)
                     if (NF > 1) then
                        do k = 2, NF
                           SR_new(i,j) = SR_new(i,j) - beta(k)*PSIbar(i,j)**k
                        end do
                     end if
                     chk_current  = chk_current + SR_new(i,j)/R(i)
                  end do
               end do
               print*,' beta(1) = ',beta(1),' alpha_g = ',alpha_g
               
            else
!         set alpha_g
               sumpp = 0.0
               sumgp = 0.0
               do i = 1, NR-1
                  do j = 1, NZ-1
                     SR_old(i,j) = SR_new(i,j)
                     PSIbar(i,j) =  (PSI_l%PSI - PSI_new(i,j))/delta_psi 
                     sumpp = sumpp + R(i)*NP*PSIbar(i,j)**(NP-1)
                     sumgp = sumgp + NF*PSIbar(i,j)**(NF-1)/R(i)
                  end do
               end do
               sumpp = -p0*sumpp
               sumgp = 0.5*g02*sumgp
               alpha_g = sumpp + pl_cur*delta_psi/grid_area
               alpha_g = mu0*alpha_g/sumgp

!    Source term=mu0*R*J
               chk_current = 0.0
               do i = 1, NR-1
                  RI2 = R(i)*R(i)
                  do j = 1, NZ-1
                     SR_new(i,j) = 0.0 
                     SR_new(i,j) = mu0*RI2*p0*NP*PSIbar(i,j)**(NP-1)
                     SR_new(i,j) = SR_new(i,j) + 0.5d0*g02*alpha_g*NF*PSIbar(i,j)**(NF-1)
                     SR_new(i,j) = SR_new(i,j)/delta_psi
                     chk_current  = chk_current + SR_new(i,j)/R(i)
                  end do
               end do
            end if
            
            chk_current  = grid_area*chk_current/mu0
!            print*,'chk_current = ',chk_current,' plasma_current = ',pl_cur
         end do   ! end loop for i1

         
         deallocate(RBoundary, ZBoundary, J_region, &
              SR_new, SR_old,  &
              PSI_old, PSIbar, PSIbar_old, &
              idxmap, ij2dof)         
      endif
      
    contains
      
      subroutine init()
        implicit none
        real(kind=8) :: sumi, cur_density
        integer :: i, j

        sumi = 0.0
        pl_area = 0.0
        grid_area = DR*DZ
        do i=0,NR
           do j=0,NZ
              J_region(i,j) = 0.0d0
              if (InsidePlasma(R(i),Z(j),pl_R,pl_Z,nfound) == 0) cycle
              J_region(i,j) = 1.0d0
              pl_area = pl_area + 1.0d0
              sumi = sumi + 1.0/R(i)
           end do
        end do
        pl_area = pl_area*grid_area
        cur_density = pl_cur/pl_area
             print*,'sumi = ',sumi,' pl_area = ',pl_area
        chk_current  = 0.0
        do i=0,NR
           do j=0,NZ
              SR_new(i,j)= mu0*pl_cur*J_region(i,j)/sumi/grid_area
              chk_current  = chk_current + pl_cur*J_region(i,j)/R(i)               
              !print*,i,j,'  SR_new(i,j) = ',SR_new(i,j)
           end do
        end do
        chk_current  = chk_current/sumi
!             print*,'chk_current = ',chk_current,' plasma_current = ',pl_cur             
        return
      end subroutine init

      SUBROUTINE sort_order(Rout, Zout, nfound, R0, Z0)

        IMPLICIT NONE
        integer, intent(in) ::  nfound
        real(kind=8), intent(in) ::  R0, Z0
        real(kind=8), intent(inout) ::  Rout(nfound), Zout(nfound)
        
        real(kind=8), parameter ::  twopi = 2*3.1415926535897932385
        real(kind=8) ::  angle(nfound), min, R_min, Z_min
        integer :: i, j, minindex
        
        do i = 1, nfound
           angle(i) = datan2(Zout(i) - Z0,Rout(i) - R0)
           if (angle(i) < 0.0 ) angle(i) = twopi + angle(i)
        end do

! sort points with key angle(i) by selection sort algorithm

        do i = 1, nfound-1
           minindex = i
           min = angle(i)
           do j = i+1, nfound
              if (min > angle(j)) then
                 minindex = j
                 min = angle(j)
              end if
           end do
           
           R_min = Rout(minindex)
           Z_min = Zout(minindex)
           
           angle(minindex) = angle(i)
           Rout(minindex) = Rout(i)
           Zout(minindex) = Zout(i)
           
           angle(i) = min
           Rout(i) = R_min
           Zout(i) = Z_min
        end do
        
        RETURN
      END SUBROUTINE sort_order
      
    END SUBROUTINE run

    subroutine set_input()
      
      use output, only: R, Z, NR, NZ, DR, DZ, PSI_new
      use params
      
      implicit none
      
      real(kind=8) :: RMIN, RMAX
      real(kind=8) :: ZMIN, ZMAX
      integer :: i, j
      real(kind=8) :: val1
      
      if (allocated(pl_R).and.allocated(pl_Z)) then
         if (size(pl_R,1).ge.min_npl.and.size(pl_Z,1).ge.min_npl) then
            is_set = .true.
            npoints = min(size(pl_R,1), size(pl_Z,1))
         else
            npoints = 0
            is_set = .false.
         endif
      endif
      !        
      if (NF.le.0.or.NP.le.0) then
         is_set = .false.
      elseif (.not.allocated(alpha).or..not.allocated(beta)) then
         is_set = .false.
      elseif (size(alpha,1).ne.NP.or.size(beta,1).ne.NF) then
         is_set = .false.
      else
         is_set = .true.
      endif
      !
      if (is_set) then
         if (allocated(R)) deallocate(R) 
         if (allocated(Z)) deallocate(Z)            
         if (NR.lt.minN) then
            is_set = .false.
         elseif (NR.lt.minN) then
            is_set = .false.
         elseif (NZ.lt.minN) then
            is_set = .false.
         else
            val1 = dlog(1.D0*NZ)/dlog(2.D0)
            if (val1.ne.floor(val1)) then ! check NZ = 2^n
               is_set =.false.
            else
               allocate (R(0:NR), Z(0:NZ))
               is_set = .true.
            endif
         endif
      endif
      
      if (is_set) then
         g0 = R0*BT0
         g02 = g0*g0
         RMIN = pl_R(1)
         RMAX = pl_R(1)
         ZMIN = pl_Z(1)
         ZMAX = pl_Z(1)
         do i = 1, npoints
            RMIN = MIN(RMIN, pl_R(i))
            RMAX = MAX(RMAX, pl_R(i))
            ZMIN = MIN(ZMIN, pl_Z(i))
            ZMAX = MAX(ZMAX, pl_Z(i))
         end do
         ! Trick for InsidePlasma routine       
         DR = (RMAX - RMIN)/NR
         DZ = (ZMAX - ZMIN)/NZ
         RMIN = RMIN - 2*DR
         RMAX = RMAX + 2*DR
         ZMIN = ZMIN - 2*DZ
         ZMAX = ZMAX + 2*DZ
         
         DR = (RMAX - RMIN)/NR
         DZ = (ZMAX - ZMIN)/NZ
         do i = 0, NR
            R(i) = RMIN + i * DR
         end do
         do i= 0, NZ
            Z(i) = ZMIN + i * DZ
         end do
         if (allocated(PSI_new)) deallocate(PSI_new)
      endif
      
      return
    end subroutine set_input

