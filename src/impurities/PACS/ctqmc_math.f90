module ctqmc_math
  
  use GlobalVariables, only: dp, pi, xi, Zero, Half, One, Two

  implicit none

  interface
     ! compute an LU factorization of a general M-by-N matrix A using
     ! partial pivoting with row interchanges
     
     subroutine ZGETRF( M, N, A, LDA, IPIV, INFO )
       integer, intent(in)    :: LDA, M, N
       integer, intent(out)   :: IPIV(*), INFO
       complex(Selected_Real_Kind(8)), intent(inout) :: A(LDA, N)
     end subroutine ZGETRF
     
     ! compute  the  inverse of a matrix using the LU factorization
     ! computed by ZGETRF

     subroutine ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO)
       integer, intent(in)    :: N, LDA, LWORK
       integer, intent(in)    :: IPIV(*)
       integer, intent(out)   :: INFO
       complex(Selected_Real_Kind(8)), intent(inout)  :: A(LDA, N)
       complex(Selected_Real_Kind(8)), intent(out) :: WORK(LWORK)
     end subroutine ZGETRI

     subroutine DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
       character, intent(in) :: JOBZ
       integer, intent(in)   :: LDZ, N
       integer, intent(out)  :: INFO
       real(selected_real_kind(8)), intent(inout) :: D( * ), E( * )
       real(selected_real_kind(8)), intent(out) :: WORK( * ), Z( LDZ, * )
     end subroutine DSTEV
     
     SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       
       INTEGER, intent(in)   ::    INFO, LDA, LDB, N, NRHS
       INTEGER, intent(out)  ::    IPIV( * )       
       real(selected_real_kind(8)), intent(in)  :: A( LDA, * )
       real(selected_real_kind(8)), intent(out) :: B( LDB, * )
     END SUBROUTINE DGESV

     SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,  &
                       WORK, LWORK, RWORK, INFO)
       CHARACTER(LEN = 1), intent(in) :: JOBVL, JOBVR
       INTEGER, intent(in)  :: N, LDA, LDVL, LDVR, LWORK
       INTEGER, intent(out) :: INFO
       REAL(selected_real_kind(8)), intent(out) :: RWORK(2*N)
       COMPLEX(selected_real_kind(8)), intent(inout) :: A(LDA, N)
       COMPLEX(selected_real_kind(8)), intent(out) :: W(N), VL(LDVL, N), VR(N, LDVR)
       COMPLEX(selected_real_kind(8)), intent(out) :: WORK(LWORK)
     END SUBROUTINE ZGEEV
    end interface

    interface Simpson
       module procedure Simpson_Dble, Simpson_Cmplx
    end interface
     
    interface inverse 
       module procedure inverse_cmplx, inverse_dble
    end interface 

contains
  !------------------------------------------------------------------
  subroutine RanGen_Ini(clock)
    !
    ! Purpose
    ! =======
    !   Initialize the Random Number Generator. Here there are two Gernerators for choosing.
    ! One is the Fortran90 default Random_Number, the other is the MT19937.
    !
    integer, intent(in) :: clock

    ! ... local vars ... 
    integer :: i, n
    integer, allocatable :: seed(:)

    call RANDOM_SEED(size = n)
    ALLOCATE(SEED(n))
    do i = 1, n
       SEED(i) = clock + 37*(i-1)
    end do
    call RANDOM_SEED(PUT = SEED)
    DEALLOCATE(SEED)

  end subroutine RanGen_Ini

  !------------------------------------------------------------------
  real(dp) function ranw()
    !
    ! Purpose
    ! =======
    !   random number generator, use the f90 intrinsic routin, random numbers
    !   distuributes in (0, 1)

    call random_number(ranw)
      
  end function ranw

  !-------------------------------------------------------------------
  subroutine diag_ZHEEV(n, A, W)
    !
    ! Purpose
    ! =======
    !
    INTEGER, intent(in)  :: n
    COMPLEX(dp), intent(inout) :: A(n, n)
    REAL(dp), intent(out)      :: W(n)      

    ! ... local vars ...
    INTEGER      :: INFO
    COMPLEX(dp)  :: WORK(4*n)
    REAL(dp)     :: RWORK(3*n-2)

    CALL ZHEEV('V', 'U', n, A, n, W, WORK, 4*n, RWORK, INFO)

    IF (INFO /= 0) THEN
       write(*, *) 'diag_Zheev is not correctly returned, aborting ...'
    END IF

  end subroutine diag_ZHEEV

  !-------------------------------------------------------------------
  subroutine diag_ZGEEV(n, A, W, VL, VR)

     IMPLICIT NONE
     INTEGER,     INTENT(IN)    :: n
     COMPLEX(DP), INTENT(IN)    :: A(n, n)
     COMPLEX(DP), INTENT(OUT)   :: VL(n, n), VR(n, n), W(n) 
 
     ! ... local vars ...
     INTEGER     :: INFO
     REAL(DP)    :: RWORK(2*N)
     COMPLEX(DP) :: WORK(3*N), B(n, n)

     ! ... executable ...
     B = A
     CALL ZGEEV('V', 'V', N, B, N, W, VL, N, VR, N,  WORK, 3*N, RWORK, INFO)
     IF (INFO /= 0) THEN
        write(*, *) 'diag_ZGEEV does not return correctly, aborting ...'
     END IF

  end subroutine diag_ZGEEV

  !--------------------------------------------------------------------
  SUBROUTINE diag_DSYEV(n, A, W)
    !
    ! Purpose
    ! =======
    !
    INTEGER, intent(in)  :: n
    REAL(KIND=8), intent(inout) :: A(n, n)
    REAL(KIND=8), intent(out)   :: W(n)
    
    ! ... local vars ...
    INTEGER      :: INFO
    REAL(KIND=8) :: WORK(4*n)
     
    CALL DSYEV('V', 'U', n, A, n, W, WORK, 4*n, INFO)
    
    IF (INFO /= 0) THEN
       write(*, *) 'diag_DSYEV is not correctly returned, aborting ...'
    END IF
    
  END SUBROUTINE diag_DSYEV

  !-------------------------------------------------------------------
  subroutine timer( current_time )
    real(dp), intent(out) :: current_time

    integer :: t1, clock_rate, clock_max

    call system_clock ( t1, clock_rate, clock_max )
    current_time = dble(t1)/dble(clock_rate)
    
  end subroutine

  !-------------------------------------------------------------------
  subroutine Trapezoid_dble(a, n, h, sum_value)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: a(n), h
    real(dp), intent(out) :: sum_value

    ! ... local vars ...
    integer :: i
    if (n < 2) then
       write(*, '(a)') 'n is smaller than 2 for Trapezoid rule integration!'
       stop
    end if
    sum_value = Half*h*(a(1) + a(n))
    do i = 2, n-1
       sum_value = sum_Value + a(i)*h
    end do
  end subroutine Trapezoid_dble

  !-------------------------------------------------------------------
  subroutine Trapezoid_cmplx(a, n, h, sum_value)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: h
    complex(dp), intent(in)  :: a(n)
    complex(dp), intent(out) :: sum_value

    ! ... local vars ...
    integer :: i
    if (n < 2) then
       write(*, '(a)') 'n is smaller than 2 for Trapezoid rule integration!'
       stop
    end if
    sum_value = Half*h*(a(1) + a(n))
    do i = 2, n-1
       sum_value = sum_value + a(i)*h
    end do
  end subroutine Trapezoid_cmplx

  !-------------------------------------------------------------------
  subroutine Simpson_Cmplx(Func, nL, h, solution)
    
    implicit none
    integer, intent(in) :: nL
    complex(dp), intent(in) ::  Func(nL)
    complex(dp), intent(out) :: solution
    real(dp), intent(in) :: h
    
    integer :: i
    
    if (Mod(nL, 2) == 0) then
       print*, "nL must be odd!"
       stop
    endif
    
    solution = (Func(1) + Func(nL))
    do i = 2, nL-1, 2
       solution = solution + 4.d0 * Func(i)
    enddo
    do i = 3, nL-2, 2
       solution = solution + 2.d0 * Func(i)
    enddo
    solution = solution * h/3.d0
    
  end subroutine Simpson_Cmplx
  
  !-----------------------------------------------------!
  subroutine Simpson_Dble(Func, nL, h, solution)
    
    implicit none
    integer, intent(in)   :: nL
    real(dp), intent(in)  :: Func(nL), h
    real(dp), intent(out) :: solution
    
    integer :: i
    
    if (Mod(nL, 2) == 0) then
       print*, "nL must be odd!"
       print*, nL
       stop
    endif
    
    solution = (Func(1) + Func(nL))
    do i = 2, nL-1, 2
       solution = solution + 4.d0 * Func(i)
    end do
    
    do i = 3, nL-2, 2
       solution = solution + 2.d0 * Func(i)
    enddo
    solution = solution * h/3.d0
    
  end subroutine Simpson_Dble
  
  !-------------------------------------------------------------------
  subroutine inverse_cmplx(matrix_a, matrix_b, n)
    !
    ! Purpose
    ! =======
    !   complex matrix inversion operation
    !
    ! Arguments
    ! =========
    !
    implicit none
    integer, intent(in) :: n
    complex(dp), intent(inout)  :: matrix_a(n, n)
    complex(dp), intent(inout) :: matrix_b(n, n)

    integer :: ipiv(n), info
    complex(dp) :: y(n, n), work(n)

    y = matrix_a
    call zgetrf(n,n,y,n,ipiv,info)
    call zgetri(n,y,n,ipiv,work,n,info)
    matrix_b = y

  end subroutine inverse_cmplx

  subroutine inverse_dble(matrix_a, matrix_b, n)
    !
    ! Purpose
    ! =======
    !   complex matrix inversion operation
    !
    ! Arguments
    ! =========
    !
    implicit none
    integer, intent(in) :: n
    real(dp), intent(inout)  :: matrix_a(n, n)
    real(dp), intent(inout) :: matrix_b(n, n)

    integer :: ipiv(n), info
    real(dp) :: y(n, n), work(n)

    y = matrix_a
    call dgetrf(n,n,y,n,ipiv,info)
    call dgetri(n,y,n,ipiv,work,n,info)
    matrix_b = y

  end subroutine inverse_dble

  !-------------------------------------------------------------------------------------------------------
  subroutine asymptotic_moment(wmax, GRe, GIm, a, b)

    ! determine the asymptotic by fitting the Green's function at the last
    ! Matsubara frequency.
    ! Note that this route is designed only for Green's function like functions,
    ! whose leading
    ! frequency dependence is 1/iw. Other functions with different scaling of
    ! 1/iw
    ! is not
    ! compatible with this routine.

    real(8), intent(in)  :: wmax      ! the largest Matsubara frequency
    real(8), intent(in)  :: GRe, GIm  ! real and imaginary part of the Green's function at wmax
    real(8), intent(out) :: a, b

    ! ... local vars ...
    real(8) :: square, eps

    square = GRe*GRe + GIm*GIm

    eps = 1.d-12 - wmax*(GIm+wmax*square)*(1.d0+4.d0*wmax*(GIm+wmax*square))
    if (eps < 0.d0) then
       a = 0.d0
       b = 0.d0
    else
       a = (-wmax*GRe + sqrt(eps))/(GIm+2.d0*wmax*square)
       b = (-wmax*GRe - sqrt(eps))/(GIm+2.d0*wmax*square)
    end if

  end subroutine asymptotic_moment

  !----------------------------------------------------------------------------
  subroutine Legendre_Polynomials(nLength, nOrder, a, x, LegenPoly)
    !
    ! Purpose
    ! =======
    !   Generate the Lengendre Polynomials for given Order = nOrder, from the recurrence formula.
    !
    integer,  intent(in)  :: nLength, nOrder
    real(dp), intent(in)  :: x(0:nLength), a
    real(dp), intent(out) :: LegenPoly(0:nOrder, 0:nLength)

    ! ... local vars ...
    integer  :: i, j
    real(dp) :: y(0:nLength)

    LegenPoly = Zero
    do j = 0, nLength
       y(j) = Two/a*X(j) - One
       LegenPoly(0, j) = One
       LegenPoly(1, j) = y(j)
       do i = 1, nOrder-1
          LegenPoly(i+1, j) = dble(2*i+1)/dble(i+1)*y(j)*LegenPoly(i, j) - dble(i)/dble(i+1)*LegenPoly(i-1, j)
       end do
    end do
    
  end subroutine

  !-----------------------------------------------------------------------------------------------------------
  subroutine Legendre_Polynomials_Matsubara(n, nOrder, a, LegenPolyMatsubara)
    !
    ! Purpose
    ! =======
    !   Generate the Legendre Polynomials in the Matsubara frquency space.
    !
    integer,     intent(in)   :: n, nOrder
    real(dp),    intent(in)   :: a
 !   complex(dp), intent(in)   :: W(n)
    complex(dp), intent(out)  :: LegenPolyMatsubara(0:nOrder, n)

    ! ... local vars ...
    integer  :: i, i1, k, IFAIL
    real(8) :: dummy, coeff_temp, w
    real(dp) :: J(0:nOrder), Y(0:nOrder), Jp(0:nOrder), Yp(0:nOrder)

    complex(kind=16) :: cdummy, tlp

    coeff_temp = a
    do i = 1, n
       dummy = Half*(Two*i - One)*Pi
       call SBESJY(dummy, nOrder, J(0:nOrder), Y(0:nOrder), JP(0:nOrder), YP(0:nOrder), IFAIL)
       if (IFAIL /= 0) then
          print*, "No convergence reached for the spherical Bessl function."
          stop
       end if

       do k = 0, nOrder
          LegenPolyMatsubara(k, i) = xi**(k+1)*J(k)*coeff_temp
       end do
       coeff_temp = -coeff_temp
    end do

  end subroutine Legendre_Polynomials_Matsubara

  !-----------------------------------------------------------------------------------------------------------
  subroutine Coefficient_LegendrePolynomials(nLength, nOrder, a, Func, LegenPoly, CoeffLegenPoly)
    !
    ! Purposes
    ! ========
    !   calculate the coefficient of the Legendre Polynomials for given function in [0, a].
    !
    ! Arguments
    ! =========
    integer,  intent(in) :: nLength             ! the length of the fitted function.
    integer,  intent(in) :: nOrder              ! order of the Legendre Polynomials
    real(dp), intent(in) :: a                   ! [0, a] is the range in which the Legendre Polynomials are orthogonal.
    real(dp), intent(in) :: LegenPoly(0:nOrder, 0:nLength)   ! Legendre Polynomials
    real(dp), intent(in) :: Func(0:nLength)       ! function wanted to fit
    real(dp), intent(out) :: CoeffLegenPoly(0:nOrder)  ! the coefficient value of the Lgendre Polynomials.

    ! ... local vars ...
    integer  :: i, j
    real(dp) :: dummy, dx, D1(0:nLength), q

    ! --- executable ---
    dx = a/DBLE(nLength)
    q  = Pi/(nOrder + 1)
    do j = 0, nOrder
       CoeffLegenPoly(j) = Zero
       dummy = Zero
       do i = 0, nLength
          D1(i) = Func(i)*LegenPoly(j, i)
       end do
       call Simpson(D1(0:nLength), nLength+1, dx, CoeffLegenPoly(j))
       CoeffLegenPoly(j) = CoeffLegenPoly(j)*(2*j+1)/a 
    end do

  end subroutine Coefficient_LegendrePolynomials

  !----------------------------------------------------------------------------
  subroutine LegendrePolynomials_Moments(nOrder, a, CoeffLegenPoly, C1, C2, C3)
    ! 
    ! Purposes
    ! ========
    !   Determine the high frequency tail of Gw from moment expansion: Gw_i = C1/(iw_i) + C2/(iw_i)**2 + C3/(iw_i)**3
    !
    integer,  intent(in)  :: nOrder
    real(dp), intent(in)  :: a, CoeffLegenPoly(0:nOrder)
    real(dp), intent(out) :: C1, C2, C3

    ! ... local vars ...
    integer :: i
    
    C1 = Zero
    C2 = Zero
    C3 = Zero
    do i = 0, nOrder, 2   ! even term
       C1 = C1 - CoeffLegenPoly(i)*Two*(Two*i + One)/a 
       C3 = C3 - CoeffLegenPoly(i)*(i+2)*i*(i+1)*(i-1)*(Two*i + One)/a/a/a 
    end do

    do i = 1, nOrder, 2   ! odd term
       C2 = C2 + CoeffLegenPoly(i)*i*(i+1)*Two*(Two*i+One)/a/a 
    end do

  end subroutine LegendrePolynomials_Moments

  !-----------------------------------------------------------------------------------
  subroutine MessageBox( Term, a )
    ! 
    ! Purpose
    ! =======
    !  display message on screen to report error or debuging information.
    !
    integer, intent(in) :: Term ! terminal 
    character(*) :: a
    write(Term, *)
    write(Term, *) TRIM(a)

  end subroutine MessageBox
  
  !-----------------------------------------------------------------------------------
  SUBROUTINE SBESJY(X, LMAX, J,Y, JP, YP, IFAIL )
    !   REAL SPHERICAL BESSEL FUNCTIONS AND X DERIVATIVES
    !            j , y , j', y'                    FROM L=0 TO L=LMAX
    !        FOR REAL X > SQRT(ACCUR) (E.G. 1D-7)    AND INTEGER LMAX
    !  J (L)  =      j/L/(X) STORES   REGULAR SPHERICAL BESSEL FUNCTION:
    !  JP(L)  = D/DX j/L/(X)            j(0) =  SIN(X)/X
    !  Y (L)  =      y/L/(X) STORES IRREGULAR SPHERICAL BESSEL FUNCTION:
    !  YP(L)  = D/DX y/L/(X)            y(0) = -COS(X)/X
    !                                                
    !    IFAIL = -1 FOR ARGUMENTS OUT OF RANGE
    !          =  0 FOR ALL RESULTS SATISFACTORY
    ! 
    !   USING LENTZ-THOMPSON EVALUATION OF CONTINUED FRACTION CF1,
    !   AND TRIGONOMETRIC FORMS FOR L = 0 SOLUTIONS.
    !   LMAX IS LARGEST L NEEDED AND MUST BE <= MAXL, THE ARRAY INDEX.
    !   MAXL CAN BE DELETED AND ALL THE ARRAYS DIMENSIONED (0:*)
    !   SMALL IS MACHINE DEPENDENT, ABOUT SQRT(MINIMUM REAL NUMBER),
    !         SO 1D-150 FOR DOUBLE PRECISION ON VAX, PCS ETC.
    !   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
    !   IN OSCILLATING REGION X .GE.  [ SQRT{LMAX*(LMAX+1)} ]
    !   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8
    !   IS THE SMALLEST NUMBER WITH 1+ACC8.NE.1 FOR OUR WORKING PRECISION
    !   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
    !   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
    !   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADILY WORSE:
    !   THE VARIABLE ERR IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
    !   
    !   NOTE: FOR X=1 AND L=100  J = 7.4 E-190     Y = -6.7+E186    1.4.94
    !
    !   AUTHOR :   A.R.BARNETT       MANCHESTER    12 MARCH 1990/95
    !                                AUCKLAND      12 MARCH 1991
    !---------------------------------------------------------------------  
    IMPLICIT    NONE
    INTEGER,  intent(in)  ::  LMAX
    INTEGER,  intent(out) ::  IFAIL
    real(dp), intent(in)  ::  X
    real(dp), intent(out) ::  J(0:LMAX), Y(0:LMAX), JP(0:LMAX), YP(0:LMAX)

    integer, PARAMETER :: LIMIT = 20000
    integer  :: NFP, L
    real(dp) ::  ERR, ACCUR, TK,SL
    real(dp) ::  XINV, CF1,DCF1, DEN, C,D, OMEGA, TWOXI
    real(dp), parameter :: SMALL = 1.D-150, THREE = 3.0D0

    ACCUR = 1.D-14                  ! suitable for double precision
    IFAIL = -1                      ! user to check on exit
    IF (X .LT. SQRT(ACCUR) )   GOTO 50

    !-------TAKES CARE OF NEGATIVE X ... USE REFLECTION FORMULA
    !-------BEGIN CALCULATION OF CF1 UNLESS LMAX = 0, WHEN SOLUTIONS BELOW

    XINV  = ONE / X
    DEN   = ONE
    IF (LMAX .GT. 0) THEN
       TWOXI =     XINV + XINV
       SL  =  REAL(LMAX)* XINV     ! used also in do loop 3
       TK  =  TWO * SL  + XINV * THREE
       CF1 =  SL                   ! initial value of CF1
       DEN =  ONE                  ! unnormalised j(Lmax,x)
       IF ( ABS(CF1) .LT. SMALL ) CF1 = SMALL
       C   = CF1                   ! inverse ratio of A convergents
       D   = ZERO                  ! direct  ratio of B convergents   
       DO  L = 1, LIMIT
          C   = TK - ONE / C
          D   = TK - D
          IF ( ABS(C) .LT. SMALL ) C = SMALL
          IF ( ABS(D) .LT. SMALL ) D = SMALL
          D   = ONE / D
          DCF1= D   * C
          CF1 = CF1 * DCF1
          IF ( D .LT. ZERO ) DEN = - DEN
          IF ( ABS(DCF1 - ONE) .LE. ACCUR ) GOTO 20
          TK   = TK + TWOXI
          NFP  = L                 ! ie number in loop
       end DO
       GOTO 50                     ! error exit, no convergence
20     CONTINUE

       ERR = ACCUR * SQRT(DBLE(NFP))    ! error estimate
       J (LMAX) = DEN                   ! lower-case j's  really
       JP(LMAX) = CF1 * DEN
       !------ DOWNWARD RECURSION TO L=0  AS SPHERICAL BESSEL FUNCTIONS
       DO L =  LMAX , 1, -1
          J (L-1)  = (SL + XINV) * J(L)   + JP(L)
          SL  =  SL - XINV
          JP(L-1)  =  SL * J(L-1)          - J(L)
       END DO
       DEN = J(0)
    ENDIF                           ! end loop for Lmax GT 0

    !------ CALCULATE THE L=0 SPHERICAL BESSEL FUNCTIONS DIRECTLY
    J (0)   =  XINV * SIN(X)
    Y (0)   = -XINV * COS(X)
    JP(0)   = -Y(0) - XINV * J(0)
    YP(0)   =  J(0) - XINV * Y(0)
    IF (LMAX .GT. 0) THEN
       OMEGA  =  J(0) / DEN
       SL  = ZERO
       DO L = 1 , LMAX
          J (L) = OMEGA * J (L)
          JP(L) = OMEGA * JP(L)
          Y (L) = SL * Y(L-1)   -   YP(L-1)
          SL  = SL + XINV
          YP(L) = Y(L-1)  -  (SL + XINV) * Y(L)
       END DO
    ENDIF

    IFAIL = 0                       ! calculations successful

    RETURN

!---------------------------------------------------------------------
!       ERROR TRAPS
!---------------------------------------------------------------------

50  IF (X .LT. ZERO) THEN
       WRITE(6,1000) X
    ELSEIF (ABS(X) < 1.d-8) THEN
       IFAIL = 0
       J(0) = ONE
       DO L = 1, LMAX
          J(L) = ZERO     ! remaining arrays untouched
       END DO
    ELSE                          ! x .le. sqrt(accur), e.g. 1D-7
       WRITE(6,1001) X
    ENDIF

1000 FORMAT(' X NEGATIVE !',1PE15.5,'    USE REFLECTION FORMULA'/)
1001 FORMAT(' WITH X = ',1PE15.5,'    TRY SMALL-X SOLUTIONS',/, &
          '    j/L/(X)  ->   X**L / (2L+1)!!          AND',/,   &
          '    y/L/(X)  ->  -(2L-1)!! / X**(L+1)'/)
    RETURN
  END SUBROUTINE SBESJY

  !------------------------------------------------------------------------------------------------
  SUBROUTINE set_sym_bl(at, nrot, s)
    !---------------------------------------------------------------------
    !! Provides symmetry operations for all bravais lattices.
    !! Tests the 24 proper rotations for the cubic lattice first, then
    !! the 8 rotations specific for the hexagonal axis (special axis c),
    !! then inversion is added.
    !

    IMPLICIT NONE
    INTEGER,  INTENT(OUT) :: nrot
    INTEGER,  INTENT(OUT) :: s(3,3,48)
    !! symmetry matrices, in crystal axis
    REAL(DP), INTENT(IN)  :: at(3, 3)

    !
    ! sin3 = sin(pi/3), cos3 = cos(pi/3), msin3 = -sin(pi/3), mcos3 = -cos(pi/3)
    !
    REAL(DP), PARAMETER :: eps1 = 1.0d-6
    REAL(DP), PARAMETER :: sin3 = 0.866025403784438597d0,  cos3 =  0.5d0, &
         msin3 =-0.866025403784438597d0, mcos3 = -0.5d0
    !
    REAL(DP) :: s0(3,3,32), overlap(3,3), rat(3), rot(3,3), value
    ! s0: the s matrices in cartesian axis
    ! overlap: inverse overlap matrix between direct lattice
    ! rat: the rotated of a direct vector ( cartesian )
    ! rot: the rotated of a direct vector ( crystal axis )
    ! value: component of the s matrix in axis basis
    INTEGER :: jpol, kpol, mpol, irot, imat(32)
    ! counters over the polarizations and the rotations
    !
    DATA s0/ 1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
            -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
            -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
             1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
             0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
             0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
             0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
             0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
             0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
             0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
             0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
             0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
            -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, &
            -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, &
             1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0, &
             1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0, &
             0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
             0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
             0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
             0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
             0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0, &
             0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0, &
             0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0, &
             0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0, &
             cos3,  sin3, 0.d0, msin3,  cos3, 0.d0, 0.d0, 0.d0,  1.d0, &
             cos3, msin3, 0.d0,  sin3,  cos3, 0.d0, 0.d0, 0.d0,  1.d0, &
             mcos3,  sin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0,  1.d0, &
             mcos3, msin3, 0.d0,  sin3, mcos3, 0.d0, 0.d0, 0.d0,  1.d0, &
             cos3, msin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
             cos3,  sin3, 0.d0,  sin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
             mcos3, msin3, 0.d0, msin3,  cos3, 0.d0, 0.d0, 0.d0, -1.d0, &
             mcos3,  sin3, 0.d0,  sin3,  cos3, 0.d0, 0.d0, 0.d0, -1.d0 /
    !
    ! ... compute the overlap matrix for crystal axis
    DO jpol = 1, 3
       DO kpol = 1, 3
          rot(kpol,jpol) = at(1,kpol)*at(1,jpol) + at(2,kpol)*at(2,jpol) + at(3,kpol)*at(3,jpol)
       ENDDO
    ENDDO
    !
    ! ... then its inverse (rot is used as work space)
    CALL inverse(rot, overlap, 3)
    !
    nrot = 1
    !
    DO irot = 1, 32
       !
       ! ... for each possible symmetry
       DO jpol = 1, 3
          DO mpol = 1, 3
             !
             ! ... compute, in cartesian coordinates the rotated vector
             rat(mpol) = s0(mpol,1,irot)*at(1,jpol) + s0(mpol,2,irot)*at(2,jpol) + &
                  s0(mpol,3,irot)*at(3,jpol)
          ENDDO

          DO kpol = 1, 3
             !
             ! ... the rotated vector is projected on the direct lattice
             rot(kpol,jpol) = at(1,kpol)*rat(1) + at(2,kpol)*rat(2) + at(3,kpol)*rat(3)
          ENDDO
       ENDDO
       !
       ! ... and the inverse of the overlap matrix is applied
       DO jpol = 1,3
          DO kpol = 1,3
             value = overlap(jpol,1)*rot(1,kpol) + overlap(jpol,2)*rot(2,kpol) + &
                  overlap(jpol,3)*rot(3,kpol)
             !
             IF ( ABS(DBLE(NINT(value))-value) > eps1 ) THEN
                !
                ! ... if a noninteger is obtained, this implies that this
                ! operation
                ! is not a symmetry operation for the given lattice
                !
                GOTO 10
             ENDIF
             !
             s(kpol,jpol,nrot) = NINT(value)
          ENDDO
       ENDDO
       !
       imat(nrot) = irot
       nrot = nrot+1
       !
10     CONTINUE
       !
    ENDDO
    !
    ! ... set the inversion symmetry (Bravais lattices have always inversion
    ! symmetry)
    DO irot = 1, nrot
       DO kpol = 1, 3
          DO jpol = 1, 3
             s(kpol,jpol,irot+nrot) = -s(kpol,jpol,irot)
          ENDDO
       ENDDO
    ENDDO
    !
    nrot = 2*nrot
    !
    RETURN
    !
  END SUBROUTINE set_sym_bl

  !------------------------------------------------------------------------------------------
  SUBROUTINE special_kpoints(at, s)
    !
    REAL(DP), INTENT(IN)  :: at(3, 3)
    INTEGER,  INTENT(OUT) :: s(3,3,48)
    !! symmetry matrices, in crystal axis    
    
    integer, parameter :: nptx=20000
    INTEGER  :: i, j, l, nk, n, n1, n2, n3
    INTEGER  :: nrot, nmax(3), nshift(3), nstart(3)
    INTEGER  :: k(3,nptx), kw(nptx), ieq(nptx)
    REAL(DP) :: bg(3, 3), xk(3,nptx), xkw(nptx)

    ! ... executable ...
    call recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
    !
    call set_sym_bl ( at, nrot, s )
    !
    do i=1,3
       nshift(i)=2
       nmax(i)=nshift(i)*nmax(i)
       nstart(i)=1
    enddo
    !
    n=0
    do n3=nstart(3),nmax(3)-1,nshift(3)
       do n2=nstart(2),nmax(2)-1,nshift(2)
          do n1=nstart(1),nmax(1)-1,nshift(1)
             n=n+1
             k(1,n)=n1
             k(2,n)=n2
             k(3,n)=n3
             kw(n)=1
             ieq(n)=0
             call check(n,k,kw,ieq,s,nrot,nmax)
          enddo
       enddo
    enddo
    !
    nk=0
    do j=1,n
       if(kw(j).gt.0) then
          nk=nk+1
          xkw(nk)=kw(j)
          do l=1,3
             xk(l,nk)=0.d0
             do i=1,3
                xk(l,nk)=xk(l,nk)+k(i,j)*bg(l,i)/nmax(i)
             enddo
          end do
       endif
    enddo
    
  END SUBROUTINE special_kpoints
  
  !-----------------------------------------------------------------------------
  subroutine check(n,k,kw,ieq,s,nrot,nmax)
    !-----------------------------------------------------------------------
    !
    integer :: n, nrot
    integer :: k(3,n),kw(n), s(3,3,nrot),kr(3),ieq(n),nmax(3)
    logical :: flag
    integer :: irot, j, naux, np
    !
    irot=1
    flag=.true.
    do while(irot.le.nrot.and.flag)
       kr(1)=0
       kr(2)=0
       kr(3)=0
       call ruotaijk ( s(1,1,irot),k(1,n),k(2,n),k(3,n),kr(1),kr(2),kr(3) )
       do j=1,3
          do while(kr(j).ge.nmax(j))
             kr(j)=kr(j)-nmax(j)
          enddo
          do while(kr(j).le.-1)
             kr(j)=kr(j)+nmax(j)
          enddo
       enddo
       np=1
       do while(flag.and.np.le.n-1)
          if( kr(1).eq.k(1,np) .and. &
               kr(2).eq.k(2,np) .and. &
               kr(3).eq.k(3,np) ) then
             kw(n)=0
             naux =np
             do while(kw(naux).eq.0)
                naux=ieq(naux)
             enddo
             ieq(n)=naux
             kw(naux)=kw(naux)+1
             flag=.false.
          endif
          np=np+1
       enddo
       irot=irot+1
    enddo
    !
    return
  end subroutine check
  !
  !-----------------------------------------------------------------------
  subroutine ruotaijk(s,i,j,k,ri,rj,rk)
    !-----------------------------------------------------------------------
    !
    integer, intent(in)  :: s(3,3), i, j, k
    integer, intent(out) :: ri,rj,rk
    !
    ri=s(1,1)*i+s(1,2)*j+s(1,3)*k
    rj=s(2,1)*i+s(2,2)*j+s(2,3)*k
    rk=s(3,1)*i+s(3,2)*j+s(3,3)*k
    !
    return
  end subroutine ruotaijk

  !--------------------------------------------------------------------
  subroutine recips (a1, a2, a3, b1, b2, b3)
    !---------------------------------------------------------------------
    !
    !   This routine generates the reciprocal lattice vectors b1,b2,b3
    !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
    !
    !     first the input variables
    !
    implicit none
    real(DP) :: a1 (3), a2 (3), a3 (3), b1 (3), b2 (3), b3 (3)
    ! input: first direct lattice vector
    ! input: second direct lattice vector
    ! input: third direct lattice vector
    ! output: first reciprocal lattice vector
    ! output: second reciprocal lattice vector
    ! output: third reciprocal lattice vector
    !
    !   then the local variables
    !
    real(DP) :: den, s
    ! the denominator
    ! the sign of the permutations
    integer :: iperm, i, j, k, l, ipol
    ! counter on the permutations
    !\
    !  Auxiliary variables
    !/
    !
    ! Counter on the polarizations
    !
    !    first we compute the denominator
    !
    den = 0
    i = 1
    j = 2
    k = 3
    s = 1.d0
100 do iperm = 1, 3
       den = den + s * a1 (i) * a2 (j) * a3 (k)
       l = i
       i = j
       j = k
       k = l
    enddo
    i = 2
    j = 1
    k = 3
    s = - s
    if (s.lt.0.d0) goto 100
    !
    !    here we compute the reciprocal vectors
    !
    i = 1
    j = 2
    k = 3
    do ipol = 1, 3
       b1 (ipol) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) ) / den
       b2 (ipol) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) ) / den
       b3 (ipol) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) ) / den
       l = i
       i = j
       j = k
       k = l
    enddo
    return
  end subroutine recips

  !----------------------------------------------------------------------------------------------------
  SUBROUTINE QDIAG(N, K, C0_in, C_in, W)
    !
    ! Purpose
    ! =======
    !   >> QDIAG joint matrix diagonalization <<
    ! W = qdiag( C0, C) finds a matrix W so that W'*C0*W has diagonal elements
    ! equal to 1 and W'*C(:,:,i)*W has smallest possible off-diagonal
    ! elements. C0 is a positive definite NxN matrix, and C is a NxNxK array
    ! of K correlation matrices.
    !
    
    INTEGER,      INTENT(IN)    :: N, K
    REAL(KIND=8), INTENT(IN)    :: C0_in(N, N), C_in(N, N, K)
    REAL(KIND=8), INTENT(OUT)   :: W(N, N)

    ! ... local vars ...
    INTEGER      :: i, j, ite, ite_max, row, col, index
    REAL(KIND=8) :: A(N, N), B(N, N), D(N, N), E(N), P(N, N), P_not(N, N), C(N, N, K)
    REAL(KIND=8) :: MM1(N, N), MM2(N, N), weight(K), C0(N, N)
    REAl(KIND=8) :: W_row(N, 1), W_new(N, 1), m1(N, 1), m2(N, 1), delta_w, Emin
    
    ! ... executable ...
    ! initialize W matrix with random number. Here the simple RAND() is used.
    ! For better random number
    ! generator, one should consider random_number.

    C0 = C0_in
    C  = C_in
   
    weight(1:K) = 1.d0/dble(K)

    ! normalize the length
    A = matmul(C0, W)
    A = matmul(Transpose(W), A) ! A now becomes a real symmetric matrix
    DO i = 1, N
       DO j = 1, N
          W(i, j) = W(i, j)/SQRT(A(j, j))
       END DO
    END DO
    
    ! do the sphering of the whole proble wrt. C0
    A = C0
    call diag_DSYEV(N, A, E)

    DO i = 1, N
       DO j = 1, N
          P(i, j)     = A(j, i)/SQRT(E(i))
          P_not(i, j) = A(i, j)*SQRT(E(j))
       END DO
    END DO
    
    do i = 1, k
       A = C(:,:, i)
       A = matmul(transpose(A), transpose(P))
       C(:, :, i) = matmul(transpose(A), transpose(P))
    end do

    C0 = matmul(P, C0)
    C0 = matmul(C0, transpose(P))
    
    W  = matmul(transpose(P_not), W)

    D = 0.d0
    DO i = 1, K
       MM1 = matmul(C(1:N, 1:N, i), W)
       MM2 = matmul(Transpose(C(1:N, 1:N, i)), W)
       D  = D + weight(i)*(matmul(MM1, transpose(MM1)) + matmul(MM2, transpose(MM2)))
    END DO
    
    ! the iteration of main loop
    ite_max = 100
    
    do ite = 1, ite_max
       delta_w = 0.d0

       do i = 1, N
          W_row(1:N, 1) = W(1:N, i)
          
          do j = 1, K
             m1(1:N, 1) = matmul(C(1:N, 1:N, j), W_row(1:N, 1))
             m2(1:N, 1) = matmul(transpose(C(1:N, 1:N, j)), W_row(1:N, 1))
             do row = 1, N
                do col = 1, N                 
                   D(row, col) = D(row, col) - weight(j)*(m1(row, 1)*m1(col, 1) + m2(row, 1)*m2(col, 1))
                end do
             end do
          end do
          
          B = D
          call diag_DSYEV(N, B, E)
          
          Emin  = E(1)
          index = 1
          do row = 1, N
             if (Emin > E(row)) then
                Emin  = E(row)
                index = row
             end if
          end do
          w_new(1:N, 1) = B(1:N, index)
          
          delta_w = max(delta_w, min(norm2(w_row-w_new), norm2(w_row+w_new)))
          
          do j = 1, K
             m1(1:N, 1)  = matmul(C(1:N, 1:N, j), W_new(1:N, 1))
             m2(1:N, 1) = matmul(transpose(C(1:N, 1:N, j)), W_new(1:N, 1))
             
             do row = 1, N
                do col = 1, N
                   D(row, col) = D(row, col) + weight(j)*(m1(row, 1)*m1(col, 1) + m2(row, 1)*m2(col, 1))
                end do
             end do
          end do  
          W(1:N, i) = w_new(1:N, 1)

       end do ! end of innier loop
       
       if (delta_w < 1.d-4) then
          write(*, *) 'convergence is achieved, stop'
          exit
       end if
       
    end do ! end of main loop

    ! revert the sphering
    W = matmul(transpose(P), W)

  END SUBROUTINE QDIAG
    
end module ctqmc_math
