!---------------------------------------------------------------------------
!  set_boundary.f90
!
!  created by Insik Choi
!  2015. 3.
!----------------------------------------------------------------------------

       SUBROUTINE set_boundary(Npl, pl_Rin,pl_Zin,RBoundary,ZBoundary, nfound)
        use input_params

        IMPLICIT NONE
        integer :: Npl
        real(kind=8),intent(inout) :: pl_Rin(Npl), pl_Zin(Npl)
        real(kind=8),intent(out),dimension(0:NR,0:NZ) :: RBoundary, ZBoundary
        integer,intent(out) :: nfound

        real(kind=8) :: theta, dtheta, ptmp, ptmp1
        real(kind=8) :: a_theta, temp, rb, zb
!        real(kind=8)  :: VW(28), EPS, pl1_R(7),pl1_Z(7)
        real(kind=8)  :: VW(12), EPS, pl1_R(3),pl1_Z(3)
        real(kind=8) :: pzmin, prmin
        integer :: i, j, k, kk, mindex
        character(len=25) :: line(3)
        integer :: open_status, stat
        integer :: ii1, ii2, ni, ICON
        character(len=20)::data_fn

  
        INTERFACE
          REAL(kind=8) FUNCTION flux_R(R0, R_shift, s_radius, delta, theta)
            IMPLICIT NONE
            real(kind=8), intent(in) ::  R0, R_shift, s_radius, delta, theta
          END FUNCTION
          REAL(kind=8) FUNCTION flux_Z(Z0, Z_shift, ellipticity, s_radius, theta)
           IMPLICIT NONE
           real(kind=8), intent(in) ::  Z0, Z_shift, ellipticity, s_radius, theta
          END FUNCTION
          INTEGER FUNCTION InsidePlasma(x1, y1,ellip_R,ellip_Z,nfound)
            use input_params, only : pi, N_divide
            IMPLICIT NONE
            real(kind=8), intent(in) ::  x1,y1,ellip_R(N_divide),ellip_Z(N_divide)
            integer, intent(in) ::  nfound
           END FUNCTION
        END INTERFACE


! Define Plasma Shape
        print*,'******', size(pl_Rin), size(pl_Zin)
        data_fn = fn(1:len_trim(fn))
        if (pl_flag == 1) then
           open(5, file=data_fn, status='old', action='read', iostat=open_status)
           if (open_status > 0 ) stop "=============== Error:: Plasma Shape Data File is not opened properly =============="
             nfound =0
             do 
                read(5, *, iostat=stat) line
                print*, line
              if (stat < 0) then
               print*,'stat = ',stat
               print*,'END OF FILE DETECTED'
               close(5)
               exit
              else if (stat == 0) then
      
                if (line(1)(1:8) == "alpha(1)") read(line(3)(1:25),'(f20.12)') alpha(1) 
                if (line(1)(1:8) == "alpha(2)") read(line(3)(1:25),'(f20.12)') alpha(2) 
                if (line(1)(1:8) == "alpha(3)") read(line(3)(1:25),'(f20.12)') alpha(3) 
                if (line(1)(1:7) == "beta(1)") read(line(3)(1:25),'(e20.12)') beta(1) 

                print*,' alpha = ',alpha,'  beta = ', beta


               if (line(1)(1:6) == "nfound") read(line(3)(1:3),'(i3)') nfound 
                print*, nfound
               if (nfound > 0) then 
                  read(5, *) 
                  do i = 1, nfound
                    read(5, *) pl_Rin(i), pl_Zin(i)
                    print*,pl_Rin(i), pl_Zin(i)
                  end do
               end if
             end if
           end do
        elseif (pl_flag == -1) then
           nfound = Npl
        else
          dtheta = 2.0d0*pi/N_divide
          do i = 1, N_divide
            theta = (i-1)*dtheta
            pl_Rin(i) = flux_R(R0,R_shift,p_radius,delta,theta)
            pl_Zin(i) = flux_Z(Z0,Z_shift,kappa,p_radius,theta)
!            pl_Rin(i) = R0 + p_radius*dcos(theta - delta*dsin(theta))
!            pl_Zin(i) = Z0 + kappa*p_radius*dsin(theta)
!           print*,pl_Rin(i), pl_Zin(i)
          end do
        end if

print*,' '
print*,'DZ = ',DZ,' DR = ',DR
print*,' '


! Initialize RBoundary and ZBoundary = 1.0
       RBoundary = 0.0
       ZBoundary = 0.0

        do i = 1, NR-1
         do j = 1, NZ-1

           if (InsidePlasma(R(i),Z(j),pl_Rin,pl_Zin, nfound)==1) then
              RBoundary(i,j) = 1.0
              ZBoundary(i,j) = 1.0

             if (InsidePlasma(R(i+1),Z(j),pl_Rin,pl_Zin, nfound)==0) then

              if (pl_flag == 1 .or. pl_flag == -1) then
                   ptmp = DSQRT((pl_Rin(1) - R(i+1))**2 + (pl_Zin(1) - Z(j))**2)
                   ni = 1
                 do k = 2, nfound
                   ptmp1 = DSQRT((pl_Rin(k) - R(i+1))**2 + (pl_Zin(k) - Z(j))**2)
                   if (ptmp1 < ptmp) then
                      ptmp = ptmp1
                      ni = k
                   end if
                 end do

!print*,ni,' ptmp = ',ptmp,' pl_Rin+(ni) = ',pl_Rin(ni),'  pl_Zin(ni) =',pl_Zin(ni),' R(i)=',R(i) 

                 ii1 = 3
                 EPS = 1.0E-6

                 if (ni == 1) then
                    pl1_Z(1) = pl_Zin(nfound)
                    pl1_R(1) = pl_Rin(nfound)
                    pl1_Z(2:3) = pl_Zin(1:2)
                    pl1_R(2:3) = pl_Rin(1:2)
                 else if (ni == nfound) then
                    pl1_Z(1:2) = pl_Zin(nfound-1:nfound)
                    pl1_R(1:2) = pl_Rin(nfound-1:nfound)
                    pl1_Z(3) = pl_Zin(1)
                    pl1_R(3) = pl_Rin(1)
                 else
                    pl1_Z = pl_Zin(ni-1:ni+1)
                    pl1_R = pl_Rin(ni-1:ni+1)
                 end if



!                 ii1 = 7
!                 EPS = 1.0E-6
!                 if (ni == 3) then
!                    pl1_Z(1) = pl_Zin(nfound)
!                    pl1_R(1) = pl_Rin(nfound)
!                    pl1_Z(2:7) = pl_Zin(1:6)
!                    pl1_R(2:7) = pl_Rin(1:6)
!                 else if (ni == 2) then
!                    pl1_Z(1:2) = pl_Zin(nfound-1:nfound)
!                    pl1_R(1:2) = pl_Rin(nfound-1:nfound)
!                    pl1_Z(3:7) = pl_Zin(1:5)
!                    pl1_R(3:7) = pl_Rin(1:5)
!                 else if (ni == 1) then
!                    pl1_Z(1:3) = pl_Zin(nfound-2:nfound)
!                    pl1_R(1:3) = pl_Rin(nfound-2:nfound)
!                    pl1_Z(4:7) = pl_Zin(1:4)
!                    pl1_R(4:7) = pl_Rin(1:4)
!                 else if (ni == nfound-2) then
!                    pl1_Z(1:6) = pl_Zin(nfound-5:nfound)
!                    pl1_R(1:6) = pl_Rin(nfound-5:nfound)
!                    pl1_Z(7) = pl_Zin(1)
!                    pl1_R(7) = pl_Rin(1)
!                 else if (ni == nfound-1) then
!                    pl1_Z(1:5) = pl_Zin(nfound-4:nfound)
!                    pl1_R(1:5) = pl_Rin(nfound-4:nfound)
!                    pl1_Z(6:7) = pl_Zin(1:2)
!                    pl1_R(6:7) = pl_Rin(1:2)
!                 else if (ni == nfound) then
!                    pl1_Z(1:4) = pl_Zin(nfound-3:nfound)
!                    pl1_R(1:4) = pl_Rin(nfound-3:nfound)
!                    pl1_Z(5:7) = pl_Zin(1:3)
!                    pl1_R(5:7) = pl_Rin(1:3)
!                 else
!                    pl1_Z = pl_Zin(ni-3:ni+3)
!                    pl1_R = pl_Rin(ni-3:ni+3)
!                 end if

! Sort the pl_Zin as pl_Zin(k+1) > pl_Zin(k) for Aitken-Lagrange interpolation
                  do k = 1, ii1-1
                    mindex = k
                    pzmin = pl1_Z(k)
                    do kk = k+1, ii1
                       if (pzmin > pl1_Z(kk))  then
                          mindex = kk
                          pzmin = pl1_Z(kk)
                       end if
                     end do
                     prmin = pl1_R(mindex)
! Exchange k and mindex                      
                     pl1_Z(mindex) = pl1_Z(k)
                     pl1_R(mindex) = pl1_R(k)
                     
                     pl1_Z(k) = pzmin
                     pl1_R(k) = prmin
                   end do

                 CALL DAKLAG(pl1_Z, pl1_R, ii1, Z(j),ii1, EPS, rb, VW, ICON)
!                 WRITE(6,620) ICON

                 RBoundary(i,j) = (rb - R(i))/DR
print*,i,j,'R(i+1)=',R(i+1),' Rb = ',Rb,' R(i) = ',R(i),' RBoundary+ =', RBoundary(i,j)
                 if (ICON ==30000) stop
              else

                temp = p_radius*kappa
                temp = (Z(j)-Z0+Z_shift)/temp
                a_theta = dasin(temp) 
                rb = flux_R(R0,R_shift,p_radius,delta,a_theta)
                RBoundary(i,j) = (rb - R(i))/DR

!print*,i,j,' R(i+1) = ',R(i+1),' R(i) = ',R(i),' DR = ',R(i+1)- R(i),' dif = ',rb - R(i)
print*,i,j,' arcsine = ',a_theta,' Rb = ',Rb,' R(i) = ',R(i),' RBoundary+ =', RBoundary(i,j)
               end if

                 if (RBoundary(i,j) < 0.0) then
                    print*,'ADJUSTED: Can not have RBoundary+(i,j) < 0.0'
                    RBoundary(i,j) = 0.0
                    RBoundary(i-1,j) = (rb - R(i-1))/DR
                 end if

                 if ( RBoundary(i,j) > 1.0 ) then
                    print*,'ADJUSTED: Can not have RBoundary+(i,j) > 1.0'
                    RBoundary(i,j) = 1.0
                    RBoundary(i+1,j) = (rb - R(i+1))/DR
                 end if

             end if

             if (InsidePlasma(R(i-1),Z(j),pl_Rin,pl_Zin,nfound)==0) then

              if (pl_flag == 1 .or. pl_flag == -1) then
                   ptmp = DSQRT((pl_Rin(1) - R(i-1))**2 + (pl_Zin(1) - Z(j))**2)
                   ni = 1
                 do k = 2, nfound
                   ptmp1 = DSQRT((pl_Rin(k) - R(i-1))**2 + (pl_Zin(k) - Z(j))**2)
                   if (ptmp1 < ptmp) then
                      ptmp = ptmp1
                      ni = k
                   end if
                 end do

                 ii1 = 3
                 EPS = 1.0E-6

                 if (ni == 1) then
                    pl1_Z(1) = pl_Zin(nfound)
                    pl1_R(1) = pl_Rin(nfound)
                    pl1_Z(2:3) = pl_Zin(1:2)
                    pl1_R(2:3) = pl_Rin(1:2)
                 else if (ni == nfound) then
                    pl1_Z(1:2) = pl_Zin(nfound-1:nfound)
                    pl1_R(1:2) = pl_Rin(nfound-1:nfound)
                    pl1_Z(3) = pl_Zin(1)
                    pl1_R(3) = pl_Rin(1)
                 else
                    pl1_Z = pl_Zin(ni-1:ni+1)
                    pl1_R = pl_Rin(ni-1:ni+1)
                 end if


!print*,'Z(j)=',Z(j),' pl1_Z = ',pl1_Z
! Sort the pl_Zin as pl_Zin(k+1) > pl_Zin(k) for Aitken-Lagrange interpolation
                  do k = 1, ii1-1
                    mindex = k
                    pzmin = pl1_Z(k)
                    do kk = k+1, ii1
                       if (pzmin > pl1_Z(kk))  then
                          mindex = kk
                          pzmin = pl1_Z(kk)
                       end if
                     end do
                     prmin = pl1_R(mindex)
! Exchange k and mindex                      
                     pl1_Z(mindex) = pl1_Z(k)
                     pl1_R(mindex) = pl1_R(k)
                     
                     pl1_Z(k) = pzmin
                     pl1_R(k) = prmin
                   end do

!print*,ni,' ptmp = ',ptmp,' pl_Rin-(ni) = ',pl_Rin(ni),'  pl_Zin(ni) =',pl_Zin(ni),' R(i-1)=',R(i-1),' R(i)=',R(i) 
!print*,'Z(j)=',Z(j),' pl1_Z = ',pl1_Z
!print*,pl1_R

                 CALL DAKLAG(pl1_Z, pl1_R, ii1, Z(j),ii1,EPS, rb, VW, ICON)

!                 WRITE(6,620) ICON

                 RBoundary(i-1,j) = (R(i)-rb)/DR
print*,i,j,'R(i-1)=',R(i-1),' Rb = ',Rb,' R(i) = ',R(i),' RBoundary- =', i-1, RBoundary(i-1,j)
                 if (ICON ==30000) stop

              else

                temp = p_radius*kappa
                temp = (Z(j)-Z0+Z_shift)/temp
                a_theta = -dasin(temp) + pi
                rb = flux_R(R0,R_shift,p_radius,delta,a_theta)
                RBoundary(i-1,j) = (R(i) - rb)/DR
print*,i,j,' arcsine = ',a_theta,'R(i-1) = ',R(i-1),' Rb = ',Rb,' R(i) = ',R(i),' RBoundary- =', RBoundary(i-1,j)
              end if

                 if (RBoundary(i-1,j) < 0.0) then
                    print*,'ADJUSTED: Can not have RBoundary-(i-1,j) < 0.0'
                    RBoundary(i-1,j) = 0.0
                    RBoundary(i,j) = (R(i+1)-rb)/DR
                 end if

                 if ( RBoundary(i-1,j) > 1.0 ) then
                    print*,'ADJUSTED: Can not have RBoundary-(i-1,j) > 1.0'
                    RBoundary(i-1,j) = 1.0
                    RBoundary(i-2,j) = (R(i-1)-rb)/DR
                 end if

             end if

             if (InsidePlasma(R(i),Z(j+1),pl_Rin,pl_Zin,nfound)==0) then

              if (pl_flag == 1 .or. pl_flag == -1) then
                   ptmp = DSQRT((pl_Rin(1) - R(i))**2 + (pl_Zin(1) - Z(j+1))**2)
                   ni = 1
                 do k = 2, nfound
                   ptmp1 = DSQRT((pl_Rin(k) - R(i))**2 + (pl_Zin(k) - Z(j+1))**2)
                   if (ptmp1 < ptmp) then
                      ptmp = ptmp1
                      ni = k
                   end if
                 end do

                 ii1 = 3
                 EPS = 1.0E-6

                 if (ni == 1) then
                    pl1_Z(1) = pl_Zin(nfound)
                    pl1_R(1) = pl_Rin(nfound)
                    pl1_Z(2:3) = pl_Zin(1:2)
                    pl1_R(2:3) = pl_Rin(1:2)
                 else if (ni == nfound) then
                    pl1_Z(1:2) = pl_Zin(nfound-1:nfound)
                    pl1_R(1:2) = pl_Rin(nfound-1:nfound)
                    pl1_Z(3) = pl_Zin(1)
                    pl1_R(3) = pl_Rin(1)
                 else
                    pl1_Z = pl_Zin(ni-1:ni+1)
                    pl1_R = pl_Rin(ni-1:ni+1)
                 end if

!print*,'ni=',ni,'R(i)=',R(i),' pl1_R = ',pl1_R
!print*,'pl1_Z = ',pl1_Z
! Sort the pl1_R as pl1_R(k+1) > pl1_R(k) for Aitken-Lagrange interpolation
                  do k = 1, ii1-1
                    mindex = k
                    prmin = pl1_R(k)
                    do kk = k+1, ii1
                       if (prmin > pl1_R(kk))  then
                          mindex = kk
                          prmin = pl1_R(kk)
                       end if
                     end do
                     pzmin = pl1_Z(mindex)
! Exchange k and mindex                      
                     pl1_Z(mindex) = pl1_Z(k)
                     pl1_R(mindex) = pl1_R(k)
                     
                     pl1_Z(k) = pzmin
                     pl1_R(k) = prmin
                   end do

!print*,'R(i)=',R(i),' pl1_R = ',pl1_R
!print*,'pl1_Z = ',pl1_Z

                 CALL DAKLAG(pl1_R(1:3), pl1_Z(1:3), ii1, R(i),ii1,EPS, zb, VW, ICON)
!                 WRITE(6,620) ICON

                 ZBoundary(i,j) = (zb - Z(j))/DZ
print*,i,j,'Z(j+1) = ',Z(j+1),' zb = ',zb,' Z(j) = ',Z(j),' ZBoundary+ =', ZBoundary(i,j)
                 if (ICON ==30000) stop

                else

                 temp = (R(i) - R0 + R_shift)/p_radius
                 if (delta ==  0.0) then
                    a_theta = dacos(temp)
                    zb = flux_Z(Z0,Z_shift,kappa,p_radius,a_theta)
                    ZBoundary(i,j) = (zb - Z(j))/DZ
                  else 
                    CALL get_atheta(0.d0,pi,temp,delta,a_theta)
                    zb = flux_Z(Z0,Z_shift,kappa, p_radius,a_theta)
                    ZBoundary(i,j) = (zb - Z(j))/DZ
                end if

!print*,i,j,' Z(j+1) = ',Z(j+1),' Z(j) = ',Z(j),' DZ = ',Z(j+1)-Z(j),' dif = ',zb - Z(j)
print*,i,j,' arcsine = ',a_theta,'Z(j+1) = ',Z(j+1),' zb = ',zb,' Z(j) = ',Z(j),' ZBoundary+ =', ZBoundary(i,j)
               end if

                 if (ZBoundary(i,j) < 0.0) then
                    print*,'ADJUSTED: Can not have ZBoundary+(i,j+1) < 0.0'
                    ZBoundary(i,j) = 0.0
                    ZBoundary(i,j-1) = (zb - Z(j-1))/DZ  
               end if

                 if ( ZBoundary(i,j) > 1.0 ) then
                    print*,'ADJUSTED: Can not have ZBoundary+(i,j+1) > 1.0'
                     ZBoundary(i,j)=1.0
                     ZBoundary(i,j+1) = (zb - Z(j+1))/DZ
                 end if

             end if

             if (InsidePlasma(R(i),Z(j-1),pl_Rin,pl_Zin,nfound)==0) then
              if (pl_flag == 1 .or. pl_flag == -1) then
                   ptmp = DSQRT((pl_Rin(1) - R(i))**2 + (pl_Zin(1) - Z(j-1))**2)
                   ni = 1
                 do k = 2, nfound
                   ptmp1 = DSQRT((pl_Rin(k) - R(i))**2 + (pl_Zin(k) - Z(j-1))**2)
                   if (ptmp1 < ptmp) then
                      ptmp = ptmp1
                      ni = k
                   end if
                 end do

                 ii1 = 3
                 EPS = 1.0E-6

                 if (ni == 1) then
                    pl1_Z(1) = pl_Zin(nfound)
                    pl1_R(1) = pl_Rin(nfound)
                    pl1_Z(2:3) = pl_Zin(1:2)
                    pl1_R(2:3) = pl_Rin(1:2)
                 else if (ni == nfound) then
                    pl1_Z(1:2) = pl_Zin(nfound-1:nfound)
                    pl1_R(1:2) = pl_Rin(nfound-1:nfound)
                    pl1_Z(3) = pl_Zin(1)
                    pl1_R(3) = pl_Rin(1)
                 else
                    pl1_Z = pl_Zin(ni-1:ni+1)
                    pl1_R = pl_Rin(ni-1:ni+1)
                 end if

! Sort the pl1_R as pl1_R(k+1) > pl1_R(k) for Aitken-Lagrange interpolation
                  do k = 1, ii1-1
                    mindex = k
                    prmin = pl1_R(k)
                    do kk = k+1, ii1
                       if (prmin > pl1_R(kk))  then
                          mindex = kk
                          prmin = pl1_R(kk)
                       end if
                     end do
                     pzmin = pl1_Z(mindex)
! Exchange k and mindex                      
                     pl1_Z(mindex) = pl1_Z(k)
                     pl1_R(mindex) = pl1_R(k)
                     
                     pl1_Z(k) = pzmin
                     pl1_R(k) = prmin
                   end do

!print*,'R(i)=',R(i),' pl1_R = ',pl1_R
!print*,'pl1_Z = ',pl1_Z

                 CALL DAKLAG(pl1_R, pl1_Z, ii1, R(i),ii1,EPS, zb, VW, ICON)
!                 WRITE(6,620) ICON

                 ZBoundary(i,j-1) = (Z(j) - zb)/DZ
print*,i,j,' Z(j-1)=',Z(j-1),' zb = ',zb,' Z(j) = ',Z(j),' ZBoundary- =', ZBoundary(i,j-1)

                 if (ICON ==30000) stop

                 else

                 temp = (R(i) - R0 + R_shift)/p_radius
                 if (delta ==  0.0) then
                    a_theta = -dacos(temp)
                    zb = flux_Z(Z0,Z_shift,kappa,p_radius,a_theta)
                    ZBoundary(i,j-1) = (Z(j) - zb)/DZ
                  else 
                    CALL get_atheta(0.d0,pi,temp,delta,a_theta)
                    a_theta = -a_theta
                    zb = flux_Z(Z0,Z_shift,kappa,p_radius,a_theta)
                    ZBoundary(i,j-1) = (Z(j) - zb)/DZ
                end if


!print*,i,j,' Z(j-1) = ',Z(j-1),' Z(j) = ',Z(j),' DZ = ',Z(j)-Z(j-1),' dif = ',Z(j) - zb
print*,i,j,' arcsine = ',a_theta,' Z(j-1)=',Z(j-1),' zb = ',zb,' Z(j) = ',Z(j),' ZBoundary- =', ZBoundary(i,j-1)
                end if

                 if (ZBoundary(i,j-1) < 0.0) then
                    print*,'ADJUSTED: Can not have ZBoundary-(i,j-1) < 0.0'
                    ZBoundary(i,j-1) = 0.0
                    ZBoundary(i,j) = (Z(j+1) - zb)/DZ
                 end if

                if ( ZBoundary(i,j-1) > 1.0 ) then
                   print*,'ADJUSTED: Can not have ZBoundary-(i,j-1) > 1.0'
                     ZBoundary(i,j-1)=1.0
                     ZBoundary(i,j-2) = (Z(j-1) - zb)/DZ
                end if

             end if

           end if ! end if of  (InsidePlasma(R(i),Z(j),pl_Rin,pl_Zin, nfound)==1) 

         end do
        end do
! Print contents of Rbound and ZBound for debug purpose
!        do i = 0, NR
!         do j = 0, NZ
!          print*,i,j,' RBound=',RBoundary(i,j),' ZBoundary = ',ZBoundary(i,j)
!         end do
!        end do
620 FORMAT('0',10X,'CONDITION CODE=',I5)
        RETURN
       END SUBROUTINE set_boundary

       SUBROUTINE get_atheta(AI,BI,cval1, delta1, x)
        use input_params, only : pi
        IMPLICIT NONE
        real(kind=8), intent(in) :: AI, BI
        real(kind=8), intent(in) :: cval1,delta1
        real(kind=8), intent(out) :: x
        external rfun
        real(kind=8) :: cval, delta
        real(kind=8) :: EPST=0.0
        integer :: ICON
        common /TSD1V/ cval, delta

        cval = cval1
        delta = delta1
 
!        AI = pi/2.0
!        BI = 3.0*pi/2.0
        
! SSL II subroutine 
         CALL DTSD1(AI,BI,rfun,EPST,x,ICON)

       END SUBROUTINE get_atheta

       REAL(kind=8) FUNCTION rfun(x)
        IMPLICIT NONE
        real(kind=8), intent(in) :: x
        real(kind=8) :: cval, delta
        common /TSD1V/ cval, delta

        rfun = cval - dcos(x + delta*dsin(x))

       RETURN
       END FUNCTION rfun

       REAL(kind=8) FUNCTION flux_R(R0, R_shift, s_radius, delta, theta)
        IMPLICIT NONE
        real(kind=8), intent(in) ::  R0, R_shift, s_radius, delta,theta

        flux_R = R0 - R_shift + s_radius*dcos(theta + delta*dsin(theta))

        return
       END FUNCTION flux_R

       REAL(kind=8) FUNCTION flux_Z(Z0, Z_shift, ellipticity, s_radius, theta)
        IMPLICIT NONE
        real(kind=8), intent(in) ::  Z0, Z_shift, ellipticity, s_radius, theta

        flux_Z = Z0 - Z_shift + s_radius*ellipticity*dsin(theta)
 
        return
       END FUNCTION flux_Z

       INTEGER FUNCTION InsidePlasma(x1, y1,ellip_R,ellip_Z, nfound)
         use input_params, only : pi, N_divide, pl_flag

          IMPLICIT NONE
          real(kind=8), intent(in) ::  x1,y1,ellip_R(N_divide),ellip_Z(N_divide)
          integer, intent(in) ::  nfound
          real(kind=8) ::  ellip1_R, ellip1_Z, ellip2_R, ellip2_Z
          real(kind=8) ::  angle
          real(kind=8) ::  dtheta, theta
          integer :: i

          INTERFACE
            REAL(kind=8) FUNCTION Angle2D(x1,y1,x2,y2)
            REAL(kind=8), intent(in) ::  x1,y1,x2,y2
            END FUNCTION

          END INTERFACE

          angle = 0.0d0

          if (pl_flag == 1 .or. pl_flag == -1) then
           do i = 1, nfound-1
             ellip1_R = ellip_R(i) - x1
             ellip1_Z = ellip_Z(i) - y1
             ellip2_R = ellip_R(i+1) - x1
             ellip2_Z = ellip_Z(i+1) - y1

             angle = angle + Angle2D(ellip1_R, ellip1_Z, ellip2_R, ellip2_Z)

           end do
         
          else
           do i = 1, N_divide-1
             ellip1_R = ellip_R(i) - x1
             ellip1_Z = ellip_Z(i) - y1
             ellip2_R = ellip_R(i+1) - x1
             ellip2_Z = ellip_Z(i+1) - y1

             angle = angle + Angle2D(ellip1_R, ellip1_Z, ellip2_R, ellip2_Z)

           end do
          end if

          if (DABS(angle) < pi) then
               InsidePlasma = 0    ! FLALSE
          else
               InsidePlasma = 1    ! TRUE
          end if

          RETURN
       END FUNCTION InsidePlasma

       REAL(kind=8) FUNCTION Angle2D(x1,y1,x2,y2)
        use input_params, only : pi, twopi

          IMPLICIT NONE
          real(kind=8), intent(in) ::  x1,y1,x2,y2
          real(kind=8) ::  theta1, theta2

          theta1 = datan2(y1,x1)
          theta2 = datan2(y2,x2)

          Angle2D = theta1 - theta2

          if (Angle2D > pi)  Angle2D =  Angle2D - twopi
          if (Angle2D < -pi)  Angle2D = Angle2D + twopi

          RETURN
       END FUNCTION Angle2D

!----- ad hoc implementations of SSL II routines (by L. Terzolo) --------- 

       SUBROUTINE DAKLAG(XI, FI, N, X,DUM, EPS, F, VW, ICON)
         implicit none
         INTEGER, PARAMETER :: NMAX=21
         INTEGER, INTENT (IN) :: N,DUM,ICON
         real(kind=8)  :: VW(12)
         INTEGER :: I,J
         REAL(kind=8), INTENT (IN) :: X,EPS
         REAL(kind=8), INTENT (OUT) :: F
         REAL(kind=8) :: X1, X2, F1, F2
         REAL(kind=8), INTENT (IN), DIMENSION (N):: XI, FI
         REAL(kind=8), DIMENSION (NMAX):: FT
         !
         IF (N.GT.NMAX) STOP 'Dimension of the data is too large.'
         DO I = 1, N
            FT(I) = FI(I)
         END DO
         !
         DO I = 1, N-1  
            DO J = 1, N-I
               X1 = XI(J)
               X2 = XI(J+I)
               F1 = FT(J)
               F2 = FT(J+1)
               FT(J) = (X-X1)/(X2-X1)*F2+(X-X2)/(X1-X2)*F1
            END DO
         END DO
         F = FT(1) 
!         DF = (ABS(F-F1)+ABS(F-F2))/2.0
         
       end SUBROUTINE DAKLAG

       SUBROUTINE DTSD1(A,B,f,EPS,sol,IERR)
        real(kind=8), intent(in) :: A, B,EPS
        real(kind=8), intent(out) :: sol
        integer :: IERR
        integer :: i, MAX
        parameter (MAX = 10000)
        real(kind=8) :: gauche, droite, fg,fd, c, f
        real(kind=8) :: step
        logical :: s, boucle

        
        gauche = a
        droite = b
        boucle=.true.
        do while(boucle)
           fg = f(gauche)
           fd = f(droite)
           if(fg*fd.lt. 0.d0) then
              boucle=.false.
!              print*,fg*fd
           else
!              print*,A,B,gauche,fg,droite,fd
              droite=droite-(B-A)*1.d-2
              if(droite .le. gauche) then
                 boucle=.false.
              end if
           end if
        end do

        if(fg*fd.lt.0) then

           do i = 1,MAX
              c = (gauche + droite)/2.d0
              fc = f(c)
              if (fg*fc .lt. 0.d0) then
                 droite = c
              else
                 gauche = c
                 fg = fc
              endif
           enddo
!        if((droite - gauche).lt. EPS) then
           sol=    (gauche + droite)/2.d0
           if(f(sol).lt. EPS) then
!          sol=    (gauche + droite)/2.d0
              IERR=0
!          print*, (droite-gauche), eps
!          pause
           else
              IERR=1
!              print*,'DTSD1 Error: Convergence: ',sol,f(sol),EPS
!              print*,A,B
           end if
        else
           IERR=2
! probably two roots => other method to find the first root
! too long !
!           step=eps!*1d-3
!           c=A
!           s=(f(c) > 0)
!           boucle=.true.
!           do while(boucle)
!              fc=f(c)
!              if(abs(fc) < EPS) then
!                 boucle=.false.
!              end if
!              c=c+step
!              if(c >B) boucle=.false.
!              if((f(c) > 0) .neqv. s) boucle=.false.
!              print*, c,B
!           end do
!           sol=c
           print*,'DTSD1 Error: Convergence: f(a)&f(b) same sign'
!           print*,A,B
!           pause

        end if

        
       END SUBROUTINE DTSD1

