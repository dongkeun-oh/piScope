!---------------------------------------------------------------------------
!  input_params.f90
!
!  module containing defintions of kstar parameter
!  created by Insik Choi
!  2015. 3. 10
!----------------------------------------------------------------------------
      MODULE input_params

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

! Define Plasma Shape and Boundary
        REAL(kind=8) :: R0=1.8, R_shift=0.0
        REAL(kind=8) :: Z0=0.0, Z_shift=0.0
        REAL(kind=8) :: p_radius = 0.5, kappa=1.0, delta=0.0

!Define number of  closed poloidal points
      INTEGER, PARAMETER :: N_divide =1000 

!Define free function parameters
!            p =  p0*p(psi_bar(i,j))**NP    
!            g**2/2 = g0**2(1+alpha*psi_bar(i,j)**NF)/2   
        REAL(kind=8) :: p0, g02, alpha_g
        REAL(kind=8) :: g0, BT0
        INTEGER :: NP=2, NF=2
        real(kind=8), dimension(:), allocatable :: alpha, beta 

! Given Plasma Current 
        REAL(kind=8) ::  pl_cur

! Define Mesh where NZ must be 2^l
        INTEGER :: NR=64, NZ=64
        REAL(kind=8) :: RMIN, RMAX
        REAL(kind=8) :: ZMIN, ZMAX

!Define max. number of  closed flux points
        integer :: npoints = 300 

!Define max. number of q value points
        integer :: N_q = 15 

! criterions & mixing values
        integer :: max_loop = 50  
        real(kind=8) ::  crit2 = 1.0d-2, mixing = 0.5d0

! negative current flag
        integer :: neg_cur = 1  ! negative current  = 1 

! plasma shape data file name and flag
        character(len=20) :: fn
        integer :: pl_flag = 1

! DEBUG FLAG
        integer :: debug = 1

! R, Z Coordinates values
        real(kind=8), dimension(:), allocatable :: R, Z 
        real(kind=8) :: DR, DZ


        namelist /input_values/R0,R_shift,Z0,Z_shift, &
                  p_radius, kappa, delta,p0,BT0,NP,NF, &
                  pl_cur,NR, NZ, npoints,N_q, & 
                  max_loop, crit2, mixing, neg_cur, &
                  fn, pl_flag, debug

        contains 

        subroutine read_input()
        implicit none
        integer :: i, j
        open(unit=5, file='input', status='old')
        read(5, nml=input_values)
        close (unit=5)

        allocate (R(0:NR), Z(0:NZ), alpha(NP), beta(NF))

        g0 = R0*BT0
        g02 = g0*g0

        RMIN = R0 - p_radius*1.2
        RMAX = R0 + p_radius*1.1
        ZMIN = Z0 - p_radius*kappa*1.2
        ZMAX = Z0 + p_radius*kappa*1.2
! Trick for InsidePlasma routine       
        RMIN = RMIN - 0.0001
        RMAX = RMAX + 0.0001
        ZMIN = ZMIN - 0.0001
        ZMAX = ZMAX + 0.0001

        DR = (RMAX - RMIN)/NR
        DZ = (ZMAX - ZMIN)/NZ

          do i = 0, NR
             R(i) = RMIN + i * DR
          end do

          do j = 0, NZ
             Z(j) = ZMIN + j * DZ
          end do


        return
        end subroutine read_input
        
      END MODULE input_params
