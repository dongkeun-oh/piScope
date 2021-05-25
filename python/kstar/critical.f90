           SUBROUTINE critical(N,idxmap,ij2dof,M_AXISR,M_AXISZ,PSI_new,PSI_ext, Xpoint)
           
	   use input_params, only: extremum_info,NR,NZ,R,Z,DR,DZ,max_loop

	   implicit none
           integer, intent(in) :: N
            integer, intent(in) :: idxmap(1:NR*NZ,1:3), ij2dof(1:NR-1,1:NZ-1)
           real(kind=8), intent(in) :: M_AXISR, M_AXISZ, PSI_new(0:NR,0:NZ)
           real(kind=8) :: HDR, HDZ, DR2, DZ2, D, delta_r, delta_z
           real(kind=8) :: R_diff, Z_diff
           !integer, intent(in) :: N_nod
           integer :: i,j,k,Saddle_flag,Extremum_flag
           type(extremum_info) :: Saddle_P(30), Extremum_P(30)
           type(extremum_info) :: PSI_ext, Xpoint

! Bivariate interpolation variables for limiter positions
           integer (kind = 4 ), parameter :: ncp=4
           integer (kind=4) ::  md
           integer (kind=4) ::  ndp, ndp1

           integer(kind=4), dimension(:), allocatable :: ipc
           real (kind = 8), dimension(:), allocatable :: XD, YD, ZD, pd

           real (kind = 8) :: R_0, Z_0, P_0


! Bivariate interpolation variables for partial differentiation at given position
           integer (kind = 4 ), parameter :: nxid =3, nyid =3

           integer(kind=4),dimension(:),allocatable :: iwkd
           real (kind = 8), dimension(:), allocatable:: wkd

           real ( kind = 8 ) :: xp, yp
           real ( kind = 8 ) :: xid(nxid), yid(nyid), zid(nxid,nyid)
           real ( kind = 8 ) :: ZX, ZY, ZXX, ZXY, ZYY


           ndp = N
           ndp1 = max(31,27+ncp)*N+nxid+nyid
          
           allocate (ipc(ncp*ndp), XD(ndp), YD(ndp), ZD(ndp), pd(5*ndp), &
                    iwkd(ndp1), wkd(8*ndp) )
    
           HDR = DR/2.0d0
           HDZ = DZ/2.0d0
           DR2 = DR*DR
           DZ2 = DZ*DZ

! Reshape arrays for bivariate interpolation and smooth surface fitting
! NDP = i*(NZ-1)+j+1  for i=0,1,...,NR,   j=0,1,2,...,NZ

           do k=1,N
              i = idxmap(k, 1)
              j = idxmap(k, 2)
              XD(k)=R(i)
              YD(k)=Z(j)
              ZD(k)=PSI_new(i,j)
           end do

! Look for critical values at first glance
          call idcldp(ndp, XD, YD, ncp, ipc)
          call idpdrv(ndp, XD, YD, zd, ncp, ipc, pd)

          Extremum_flag = 0
          Saddle_flag = 0

          do i = 1, NDP
            D=pd(5*i)*pd(5*i-2) - pd(5*i-1)*pd(5*i-1)

            delta_r = - pd(5*i)*pd(5*i-4) + pd(5*i-1)*pd(5*i-3)
            delta_r = delta_r/D
            delta_z =  pd(5*i-1)*pd(5*i-4) - pd(5*i-2)*pd(5*i-3)
            delta_z =  delta_z/D

            if ((DABS(delta_r) <= HDR) .and. (DABS(delta_z) <= HDZ)) then
               if ((D > 0.0d0) .and. (DSQRT((1.80 - XD(i))*(1.80 - XD(i)) +  YD(i)*YD(i)) < 0.5)) then
                     Extremum_flag = Extremum_flag + 1
                     print*, 'Exist Extremum Point ', Extremum_flag
               print*,'D = ',D
               print*,'delta_r = ',delta_r, '   delta_z = ',delta_z
               print*,'R = ',XD(i), '   Z = ',YD(i), '    PSI = ',ZD(i)
               print*,'Extremum R = ',XD(i)+delta_r, ' Extremum  Z = ',YD(i)+delta_z, '    PSI = ',ZD(i)
               print*,' '
                  if (Extremum_flag < 30) then
                     Extremum_P(Extremum_flag)%X=XD(i)
                     Extremum_P(Extremum_flag)%Y=YD(i)
                     Extremum_P(Extremum_flag)%PSI=ZD(i)
                   else
                     print*,'CHECK:: There are more than 30 minimum points'
                   end if
                else
                     Saddle_flag = Saddle_flag + 1
!               print*, 'Exist Saddle Point ', Saddle_flag
!               print*,'D = ',D
!               print*,'delta_r = ',delta_r, '   delta_z = ',delta_z
!               print*,'R = ',XD(i), '   Z = ',YD(i), '    PSI = ',ZD(i)
!               print*,' '
                     if (Saddle_flag < 30) then
                       Saddle_P(Saddle_flag)%X=XD(i)
                       Saddle_P(Saddle_flag)%Y=YD(i)
                       Saddle_P(Saddle_flag)%PSI=ZD(i)
                     else
                       print*,'CHECK:: There are more than 30 saddle points'
                     end if
                 end if
             end if
          end do

!Find extremum in PSI_new(i,j) or ZD
         if (Extremum_flag /= 0) then
                  PSI_ext%X = Extremum_P(1)%X
                  PSI_ext%Y = Extremum_P(1)%Y
                  PSI_ext%PSI = Extremum_P(1)%PSI
                  R_diff = (Extremum_P(1)%X - M_AXISR)**2
                  R_diff = R_diff + (Extremum_P(1)%Y - M_AXISZ)**2
            do i = 2, Extremum_flag
                  Z_diff = (Extremum_P(i)%X - M_AXISR)**2
                  Z_diff = Z_diff + (Extremum_P(i)%Y - M_AXISZ)**2
              if (Z_diff < R_diff) then
                 R_diff = Z_diff
                 PSI_ext%X = Extremum_P(i)%X
                 PSI_ext%Y = Extremum_P(i)%Y
                 PSI_ext%PSI = Extremum_P(i)%PSI
              end if
            end do
          else 
! Find extreme values near axis by using Newton iteration method
! Guess initial value about psi extrimum
            R_0 = M_AXISR
            Z_0 = M_AXISZ

            md = 1
              
            xp=R_0
            yp=Z_0
            do k = 1, max_loop

              xid(1) = xp - DR
              yid(1) = yp - DZ
              xid(2) = xp
              yid(2) = yp
              xid(3) = xp + DR
              yid(3) = yp + DZ


              call idsfft ( md, ncp, ndp, XD, YD, ZD, nxid, nyid, xid, yid, zid, iwkd, wkd)

              if (k==1) P_0 = zid(2,2)
   
              print*,' '
              print*,' AKIM_DIFF xp = ',xp,'  yp = ',yp,'   PSI = ',zid(2,2)
   
              ZX = zid(3,2) - zid(1,2)
              ZX = 0.5d0 * ZX/DR
   
              ZY = zid(2,3) - zid(2,1)
              ZY = 0.5d0 * ZY/DZ
   
              ZXX = zid(3,2) - 2.0d0 * zid(2,2) + zid(1,2)
              ZXX = ZXX/DR2
   
              ZXY = zid(3,3) - zid(3,1) - zid(1,3) + zid(1,1)
              ZXY = 0.25d0 * ZXY/DR
              ZXY = ZXY/DZ

              ZYY = zid(2,3) - 2.0d0 * zid(2,2) + zid(2,1)
              ZYY = ZYY/DZ2

              D = ZXX * ZYY - ZXY * ZXY

              delta_r = (ZXY*ZY - ZYY*ZX)/D
              delta_z = (ZXY*ZX - ZXX*ZY)/D

              xp = xp+delta_r
              yp = yp+delta_z

              md = 2

               print*,'D = ',D
               print*,'AKIM delta_r = ',delta_r, '   delta_z = ',delta_z
               print*,' '

              if (DABS(delta_r) <HDR .and. DABS(delta_z) < HDZ) exit

          end do

          if (D > 0.0d0) then
              Extremum_flag = Extremum_flag + 1
              PSI_ext%X = xp - delta_r
              PSI_ext%Y = yp - delta_z
              PSI_ext%PSI = zid(2,2)
           else
              print*,'WARNING::  NOT FOUND EXTREMUM VALUES '
! set extremum as center
              PSI_ext%X = R_0
              PSI_ext%Y = Z_0
              PSI_ext%PSI = P_0

              Saddle_flag = Saddle_flag + 1
              Saddle_P(Saddle_flag)%X = xp - delta_r
              Saddle_P(Saddle_flag)%Y = yp - delta_z
              Saddle_P(Saddle_flag)%PSI = zid(2,2)
           end if
          print*,' '
          print*,'  PSI_ext.R = ',PSI_ext%X,'  PSI_ext.Z = ',PSI_ext%Y,'   PSI_ext.PSI = ',PSI_ext%PSI

          end if


! Find X point
          if (Saddle_flag /= 0) then
             Xpoint%X = Saddle_P(1)%X
             Xpoint%Y = Saddle_P(1)%Y
             Xpoint%PSI = Saddle_P(1)%PSI
            do i = 1, Saddle_flag
              if (Saddle_P(i)%Y < Xpoint%Y) then
                 Xpoint%X = Saddle_P(i)%X
                 Xpoint%Y = Saddle_P(i)%Y
                 Xpoint%PSI = Saddle_P(i)%PSI
              end if
            end do
          print*,'  Xpoint.R = ',Xpoint%X,'  Xpoint.Z = ',Xpoint%Y,' Xpoint.PSI = ',Xpoint%PSI
          end if

          return
          end subroutine

