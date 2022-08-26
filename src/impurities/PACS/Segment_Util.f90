!******************************COPYRIGHT************************************************
! (c) Crown copyright, Gang Li < 2020 >. All rights reserved.
!
! This routine has been licensed to the valid partners for use and distribution
! under the collaboration agreement, subject to the terms and conditions set out
! therein.
!
! [ligang at Shanghaitech.edu.cn]
!******************************COPYRIGHT************************************************

MODULE Segment_Util
  
  USE GlobalVariables
  USE MPI_MOD
  
  IMPLICIT NONE
  
  !---------------------------------------------------------------------------------------
  ! Description:
  !
  !   Module Segment_Util contains the supporting routines for CT-HYB. The central task is
  !   to generate hybridization function in imaginary-time space.
  !
  ! Routines list:
  !
  !    ReadIn_Parameters
  !    Segment_Initialize
  !    Segment_Finalize
  !
  ! Code hisotry: 2021.1.5 adapted from code 2007.10
  !
  ! Current Code Owner: Gang Li
  !
  ! Code Description:
  !   Language: Fortran 90.
  !---------------------------------------------------------------------------------------

  ! NAMELIST_Model_Parameter
  INTEGER     ::         nSpin = 2                   ! two different spin species
  INTEGER     ::        nOrbit = 3                   ! number of orbitals
 
  REAL(DP)    ::          Beta = 50.d0               ! inverse temperatures
  REAL(DP)    ::            xU = 0.d0                ! interaction parameters
  REAL(DP)    ::            xJ = 0.0d0               ! interaction parameters
  INTEGER     ::        nOmega = 200                 ! number of Matsubara frequencies
  INTEGER     ::          nTau = 2000                ! mesh of imaginary-time domain
  
  ! NAMELIST_MC_Parameter
  INTEGER     ::      nDMFT    = 1                   ! number of DMFT iterations
  INTEGER     ::      nWarm    = 50000               ! number of MC warm-up steps
  INTEGER     ::      nMeasure = 1000000             ! number of MC measure steps
  INTEGER     ::      nBin     = 20                  ! binning analysis

  INTEGER, PARAMETER ::   nLegenPoly = 150           ! maximal order of Legendre polynomial, it is adversized to not exceed 150
  ! otherwise, numerically unstable in evaluating Bessel function.

  LOGICAL     ::      DEBUG    = .False.     ! internal switch for coding debug, not open to users
 
  COMPLEX(DP), ALLOCATABLE ::    Omega(:)
  COMPLEX(DP), ALLOCATABLE ::      G0w(:, :, :)
  COMPLEX(DP), ALLOCATABLE ::       Gw(:, :, :)
  COMPLEX(DP), ALLOCATABLE ::    Sigma(:, :, :)
  COMPLEX(DP), ALLOCATABLE ::      H0w(:, :, :)
  REAL(DP),    ALLOCATABLE ::     Eimp(:, :)           ! impurity energy level, containing chemical potential  
  REAl(DP),    ALLOCATABLE ::      H0t(:, :, :)
  REAL(DP),    ALLOCATABLE ::  ImpurityOccupancy(:, :)  ! DMFT particle number, i.e. -G(tau=beta) 
  REAL(DP),    ALLOCATABLE ::       Gt(:, :, :)
  REAL(DP),    ALLOCATABLE ::      Gst(:, :, :)            ! improved estimator, G*Sigma
  REAL(DP),    ALLOCATABLE ::     Umat(:, :, :, :)

  REAL(DP),    ALLOCATABLE :: LegenPoly(:, :)           ! Legendre Polynormial basis in imaginary-time space
  REAL(DP),    ALLOCATABLE :: CoeffLegenPoly(:)         ! Coefficient of Legnedre Polynormial expansion of fw
  COMPLEX(DP), ALLOCATABLE :: LegenPolyMatsubara(:, :)  ! Legendre Polynormial basis in Matsubara frequency space
  
  INTERFACE Matmul_GEMM
     ! double real precision 
     REAL(KIND=8) FUNCTION Matmul_DGEMM(Matrix_A, Matrix_B, nDim)
       INTEGER,  INTENT(IN) :: nDim
       REAL(KIND=8), INTENT(IN) :: Matrix_A(nDim, nDim), Matrix_B(nDim, nDim)
     END FUNCTION Matmul_DGEMM
     
     ! double complex
     COMPLEX(KIND=8) FUNCTION Matmul_ZGEMM(Matrix_A, Matrix_B, nDim)
       INTEGER,     INTENT(IN) :: nDim
       COMPLEX(KIND=8), INTENT(IN) :: Matrix_A(nDim, nDim), Matrix_B(nDim, nDim)
     END FUNCTION Matmul_ZGEMM
  END INTERFACE Matmul_GEMM
  
CONTAINS
  !---------------------------------------------------------------------------------------
  SUBROUTINE Segment_Initialize

    IMPLICIT NONE

    !
    ! Purposes
    ! ========
    !    (1) allocate memory 
    !    (1) Determine the interaction matrix from the input xU and xJ parameters
    !
    ! Code developer: Gang Li
    ! Code history:   2021.1.5
    !

    ! ... local vars ...
    INTEGER  :: iOmega, iOrbit, iSpin, jOrbit, jSpin, iTau
    REAL(dp) :: xTau(0:nTau)
    
    ALLOCATE(Omega(nOmega))
    ALLOCATE(G0w(nOmega, nOrbit, nSpin))
    ALLOCATE(Gw(nOmega, nOrbit, nSpin))
    ALLOCATE(Sigma(nomega, nOrbit, nSpin))
    ALLOCATE(H0w(nOmega, nOrbit, nSpin))
    ALLOCATE(H0t(0:nTau, nOrbit, nSpin))
    ALLOCATE(Umat(nOrbit, nSpin, nOrbit, nSpin))
    ALLOCATE(Gt(0:nTau, nOrbit, nSpin))
    ALLOCATE(Gst(0:nTau, nOrbit, nSpin))
    ALLOCATE(Eimp(nOrbit, nSpin))
    ALLOCATE(ImpurityOccupancy(nOrbit, nSpin))

    ! Legendre Polynomial in imaginary-time and frequency spaces
    ALLOCATE(LegenPoly(0:nLegenPoly, 0:nTau))
    ALLOCATE(CoeffLegenPoly(0:nLegenPoly))
    ALLOCATE(LegenPolyMatsubara(0:nLegenPoly, nOmega))
    
    !   --- generate interaction matrix ---
    Umat = Zero
    DO iSpin = 1, nSpin
       DO jSpin = 1, nSpin
          DO iOrbit = 1, nOrbit
             DO jOrbit = 1, nOrbit
                
                IF (iSpin == jSpin .AND. iOrbit /= jOrbit) THEN
                   IF (iOrbit /= jOrbit) UMat(iOrbit, iSpin, jOrbit, jSpin) = xU - Three*xJ
                ELSE
                   IF (iOrbit == jOrbit) THEN
                      Umat(iOrbit, iSpin, jOrbit, jSpin) = xU
                   ELSE
                      Umat(iOrbit, iSpin, jOrbit, jSpin) = xU - Two*xJ
                   END IF
                END IF
                
             END DO
          END DO
       END DO
    END DO
    
    !   --- generate Matsubrar frequencies ----
    DO iOmega = 1, nOmega
       Omega(iOmega) = xi*Pi/Beta*(Two*iOmega - One)
    END DO

    !   ---- generate Legendre Polynomial basis in imaginary-time and frequency space
    DO iTau = 0, nTau
       xTau(iTau) = Beta/DBLE(nTau)*DBLE(iTau)
    END DO
    CALL Legendre_Polynomials(nTau, nLegenPoly, Beta, xTau(0:nTau), LegenPoly(0:nLegenPoly, 0:nTau))    
    CALL Legendre_Polynomials_Matsubara(nOmega, nLegenPoly, Beta, LegenPolyMatsubara(0:nLegenPoly, 1:nOmega))

  END SUBROUTINE SEGMENT_INITIALIZE

  !--------------------------------------------------------------------------------------------------
  SUBROUTINE Segment_Finalize
    
    DEALLOCATE(Omega)
    DEALLOCATE(G0w)
    DEALLOCATE(Gw)
    DEALLOCATE(Sigma)
    DEALLOCATE(H0w)
    DEALLOCATE(H0t)
    DEALLOCATE(Gt)
    DEALLOCATE(Gst)
    DEALLOCATE(Umat)
    DEALLOCATE(Eimp)
    DEALLOCATE(ImpurityOccupancy)

    DEALLOCATE(LegenPoly)
    DEALLOCATE(CoeffLegenPoly)
    DEALLOCATE(LegenPolyMatsubara)
    
  END SUBROUTINE SEGMENT_FINALIZE  
  
  !----------------------------------------------------------------------------------------------------------
  SUBROUTINE Segment_HybridizationFunc
    USE ctqmc_math, only : asymptotic_moment

    IMPLICIT NONE
    
    !
    ! Purposes
    ! ========
    !    (1) Calculate the hybridization function in imaginary-time domain
    !    (2) Output Weiss Green's fuction and the hybridization function to files
    !
    !  NOTE: This is the only routine one needs to modify for a given Hamiltonian
    !
    ! Code developer: Gang Li
    ! Code history:   2021.1.5
    !    
    ! ... executable ...
    
    ! ... local vars ...
    INTEGER           :: iOmega, iOrbit, iSpin, iTau
    REAL(DP)          :: Tau, wmax, GRe, GIm, a, b
    CHARACTER(Len=30) :: FLE

    ! from H0w to G0w
    wmax = imag(Omega(nOmega))
    DO iSpin = 1, nSpin
       DO iOrbit = 1, nOrbit
          GRe = dble(G0w(nOmega, iOrbit, iSpin))
          GIm = imag(G0w(nOmega, iOrbit, iSpin))
          CALL Asymptotic_moment(wmax, GRe, GIm, a, b)
          
          ! from G0w to G0t
          DO iTau = 0, nTau
             Gt(iTau, iOrbit, iSpin) = Zero
             tau = Beta/DBLE(nTau)*iTau
             DO iOmega = 1, nOmega
                Gt(iTau, iOrbit, iSpin) = Gt(iTau, iOrbit, iSpin) + Two/Beta*( exp(-Omega(iOmega)*tau) * &
                     ( G0w(iOmega, iOrbit, iSpin) - Half/(Omega(iOmega)-a) - Half/(Omega(iOmega)-b) ) )
             END DO
             Gt(iTau, iOrbit, iSpin) = Gt(iTau, iOrbit, iSpin) - Half*exp((Beta-tau)*a)/(exp(Beta*a) + One) - Half*exp((Beta-tau)*b)/(exp(Beta*b) + One)
          END DO

          ! from H0w to H0t
          CALL FFT_Matsubara2Tau(H0w(1:nOmega, iOrbit, iSpin), nOmega, nTau, H0t(0:nTau, iOrbit, iSpin))
!          DO iTau = 0, nTau
!             H0t(iTau, iOrbit, iSpin) = Zero
!             Tau = Beta/DBLE(nTau)*iTau
!             DO iOmega = 1, nOmega
!                H0t(iTau, iOrbit, iSpin) = H0t(iTau, iOrbit, iSpin) + Two/Beta*( exp(-Omega(iOmega)*tau) * & 
!                     ( H0w(iOmega, iOrbit, iSpin) - Eimp(iOrbit, iSpin) - Half*(a+b) - (a-b)**2/4.d0/(Omega(iOmega) - Half*(a+b) )) )
!             END DO
!             H0t(iTau, iOrbit, iSpin) = H0t(iTau, iOrbit, iSpin) - (a-b)**2/4.d0*exp((Beta-tau)*(a+b)*Half)/(exp(Beta*(a+b)*Half) + One)
!          END DO
          
       END DO
    END DO

    IF (ID == master) THEN
       ! output weiss Green's function in Matsubara frequency domain
       WRITE(FLE, '(A)') 'G0w.dat'
       OPEN(UNIT=1, FILE=TRIM(FLE), STATUS='UNKNOWN')
       DO iOmega = 1, nOmega
          WRITE(1, '(100f21.15)') IMAG(Omega(iOmega)), (G0w(iOmega, 1:nOrbit, iSpin), iSpin=1, 2)
       END DO
       CLOSE(1)
       
       ! output weiss Green's function in imaginary-time space
       WRITE(FLE, '(A)') 'G0t.dat'
       OPEN(UNIT=1, FILE=TRIM(FLE), STATUS='UNKNOWN')
        DO iTau = 0, nTau
           Tau  = Beta/DBLE(nTau)*DBLE(iTau)
           WRITE(1, '(100f21.15)') Beta/DBLE(nTau)*DBLE(iTau), (Gt(iOmega, 1:nOrbit, iSpin), iSpin=1, 2)
        END DO
       CLOSE(1)

       ! output hybridization function in Matsubara frequency domain
       WRITE(FLE, '(A)') 'H0w.dat'
       OPEN(UNIT=1, FILE=TRIM(FLE), STATUS='UNKNOWN')
       DO iOmega = 1, nOmega
          WRITE(1, '(100f21.15)') IMAG(Omega(iOmega)), (H0w(iOmega, 1:nOrbit, iSpin), iSpin=1, 2)
       END DO
       CLOSE(1)
       
       ! output hybridization function in imaginary-time domain
       WRITE(FLE, '(a)') 'H0t.dat'
       OPEN(UNIT=1, FILE=TRIM(FLE), STATUS='UNKNOWN')
        DO iTau = 0, nTau
           Tau  = Beta/DBLE(nTau)*DBLE(iTau)
           WRITE(1, '(f21.15\)') Beta/DBLE(nTau)*DBLE(iTau)
           DO iSpin = 1, 2
             DO iOrbit = 1, nOrbit
               WRITE(1, '(f21.15\)') HybridizationFunc(Tau, iOrbit, iSpin)
             ENDDO
           ENDDO
           WRITE(1,*)
        END DO
       CLOSE(1)
    END IF
    
    Sigma = Zero
    Gw    = Zero
    Gt    = Zero
    Gst   = Zero

  END SUBROUTINE Segment_HybridizationFunc

  !--------------------------------------------------------------------------------------------------------------------------------------------
  REAL(KIND=8) FUNCTION HybridizationFunc(Tau, Idx_Orbit, Idx_Spin)

    IMPLICIT NONE
    
    !
    ! Purpose
    ! =======
    !    We calculate the hybridization for an arbitary given Tau by linear interpolating
    ! the H0t array determined at uniform tau-mesh.
    !
    ! Code developer: Gang Li
    ! Code history: 2021.1.5
    !
    
    REAL(DP), INTENT(IN) :: Tau
    INTEGER,  INTENT(IN) :: Idx_Orbit, Idx_Spin
    
    ! ... local vars ...
    INTEGER  :: Idx_Tau
    REAL(DP) :: Tau_Local

    ! ... executable ...
    IF (Tau >= Zero) THEN
       Tau_local = -Tau + Beta 
    ELSE
       Tau_local = -Tau 
    END IF

    Tau_local = Tau_local/Beta*DBLE(nTau)
    Idx_Tau = INT(Tau_local)
    
    IF (Idx_Tau == nTau) THEN
       HybridizationFunc = H0t(Idx_Tau, Idx_Orbit, Idx_Spin)
    ELSE
       HybridizationFunc = H0t(Idx_Tau, Idx_Orbit, Idx_Spin) + &
            (Tau_local - Idx_Tau)*(H0t(Idx_Tau+1, Idx_Orbit, Idx_Spin) - H0t(Idx_Tau, Idx_Orbit, Idx_Spin))
    END IF

    IF (Tau >= Zero) HybridizationFunc = -HybridizationFunc

  END FUNCTION HybridizationFunc
  
  !-------------------------------------------------------------------------------------
  SUBROUTINE FFT_Matsubara2Tau(fw, n, m, ftau)
    
    IMPLICIT NONE

    !
    ! Purposes
    ! ========
    !   Calculate the imaginary-time hybridization function Delta(tau) using the
    !   discrete
    !   Matsubara frequency summation augmented with high-frequency tail:
    !   f(iw_n) = C1/iw_n + C2/(iw_n)^2 + c3/(iw_n)^3
    !    f(tau) = -C1/2 + C2/4*(-beta + 2tau) + C3/4(beta tau - tau^2)
    !
    !   NOTE: C1, C2, C3 are averaged over the last 20 Matsubara frequency
    !   points, so the user
    !    has to enusre the asymptotic tail to be recovered by the 20 Matsubara
    !    frequencies.
    !
    ! Code developer: Gang Li
    ! Code history: 2021.1.5
    !
    
    INTEGER,     INTENT(IN)  :: n             ! number of Matsubara frequency
    INTEGER,     INTENT(IN)  :: m             ! number of tau mesh
    COMPLEX(DP), INTENT(IN)  :: fw(n)
    REAL(DP),    INTENT(OUT) :: ftau(0:m)
    
    ! ... local vars ...
    INTEGER     :: i, j
    !    REAL(DP)    :: G_asp_m1, G_asp_m2, re1, im1, wmax, sq

    REAL(DP)  :: tau, C0, C1, C2, C0_new, error

    C0 = Zero
    C1 = Zero
    C2 = Zero
   
    IF (n <= 30) STOP ' >>> ERROR: Not enough Matsubara frequency for FFT.'

    ! calculating C0, C1, C3 only needs two fw. We average their values over the last 30 fw. 
    ! A stable soluion of C0 is more reliably determined by iteration
    C0    = DBLE(fw(n))  ! initial value
    error = 1.d0
    j     = 0
    do while(error > 1.d-5)  
       C0_new = 0.d0
       j      = j + 1
       do i = n-29, n
          C0_new = C0_new + ( DBLE(fw(n)) - DBLE(Omega(i)/Omega(n))*DIMAG(fw(n))/DIMAG(fw(i))*(DBLE(fw(i)) - C0) )/30.d0
       end do
       error = ABS(C0_new - C0)
       C0    = C0_new
       if (j > 100) then 
          write(*, *) ' ERROR: convergence is not achieved for C0'
          STOP
       end if
    end do

    DO i = n-29, n
       C2  = C2 + DBLE(fw(i) - C0)*DBLE(Omega(i)/fw(i))/30.d0
       C1  = C1 - (C2*C2 - DBLE(Omega(i)*Omega(i)))*DBLE(fw(i)/Omega(i))/30.d0
    END DO

    ! fourier transformation by substracting the high-frequency tail
    ftau = Zero
    DO i = 0, m
       tau = Beta/DBLE(m)*i
       DO j = 1, n-30
          ftau(i) = ftau(i) + Two/Beta*DBLE( exp(-Omega(j)*tau) * (fw(j) - C0 - C1/(Omega(j) - C2)) )
       END DO
       ftau(i) = ftau(i) - C1*exp((Beta-tau)*C2)/(exp(Beta*C2) + One)
       if (ftau(i) > -1.d-12) ftau(i) = -1.d-12
    END DO

  END SUBROUTINE FFT_Matsubara2Tau

END MODULE Segment_Util

