      MODULE output
        IMPLICIT NONE

        INTEGER, PARAMETER :: N_divide = 2000
! Define Mesh where NZ must be 2^l
        INTEGER :: NR=64, NZ=64
! R, Z Coordinates values + PSI_new (magnetic flux)
        real(kind=8), dimension(:), allocatable :: R, Z
        real(kind=8), dimension(:,:), allocatable :: PSI_new
        real(kind=8) :: DR, DZ
!        real(kind=8), dimension(:), allocatable :: RBout, ZBout
        real(kind=8) :: RBout(N_divide), ZBout(N_divide)   
        logical :: is_good
      END MODULE output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE params

      IMPLICIT NONE

!Define number of  closed poloidal points
      INTEGER, PARAMETER :: min_npl =10, minN = 16

        INTEGER :: NP=2, NF=2
        real(kind=8), dimension(:), allocatable :: alpha, beta 

! Given Plasma Current 
        REAL(kind=8) ::  pl_cur
! Free Parameters
        REAL(kind=8) :: p0=1.0D0, BT0=1.9994D0
!
        REAL(kind=8) :: g0, g02 ! g0 = R0*BT0

!Define max. number of  closed flux points
        integer :: npoints = 300
        REAL(kind=8), dimension(:), allocatable :: pl_R(:), pl_Z(:)
        REAL(kind=8) :: R0=1.8D0, Z0=0.0D0

!!Define max. number of q value points
!        integer :: N_q = 15 

! criterions & mixing values
        integer :: max_loop = 50  
        real(kind=8) ::  crit2 = 1.0d-2, mixing = 0.5d0

! negative current flag
        integer :: neg_cur = 1  ! negative current  = 1 

! plasma shape data flag
        integer :: pl_flag = -1

! DEBUG FLAG
        integer :: debug = 1
        logical :: is_set = .false.
        
      END MODULE params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE input_params
        use params
        use output

        implicit none

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
        REAL(kind=8) :: R_shift=0.0
        REAL(kind=8) :: Z_shift=0.0
        REAL(kind=8) :: p_radius = 0.5, kappa=1.0, delta=0.0

!Define free function parameters
!            p =  p0*p(psi_bar(i,j))**NP    
!            g**2/2 = g0**2(1+alpha*psi_bar(i,j)**NF)/2   
        REAL(kind=8) :: alpha_g

!Define max. number of q value points
        integer :: N_q = 15 
! plasma shape data file name and flag
        character(len=20) :: fn
        
        namelist /input_values/R0,R_shift,Z0,Z_shift, &
                  p_radius, kappa, delta,p0,BT0,NP,NF, &
                  pl_cur,NR, NZ, npoints,N_q, & 
                  max_loop, crit2, mixing, neg_cur, &
                  fn, pl_flag, debug

        contains 

          subroutine read_input()
            implicit none
            integer :: i, j
            real(kind=8) :: RMIN,RMAX,ZMIN,ZMAX

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
