!******************************COPYRIGHT************************************************
! (c) Crown copyright, Gang Li < 2020 >. All rights reserved.
!
! This routine has been licensed to the valid partners for use and distribution
! under the collaboration agreement, subject to the terms and conditions set out
! therein.
!
! [ligang at Shanghaitech.edu.cn]
!******************************COPYRIGHT************************************************

MODULE Segment_Phys

  USE Segment_MonteCarlo
  IMPLICIT NONE

  !---------------------------------------------------------------------------------------
  ! Description:
  !
  !   Here, we determine some interesting physical observables from statistitical average.
  !   These include
  !   (1) Gt --- imaginary-time Green's function from both raw average and Legendre Polynomail fit
  !   (2) Gw --- Matsubara Green's function from Legendre Polynomails
  !   (3) Statistical information on Monte Carlo update
  !   (4) 
  !
  ! Code hisotry: 2021.1.5 adapted from code June 2007
  !
  ! Current Code Owner: Gang Li
  !
  ! Code Description:
  !   Language: Fortran 90.
  !-------------------------------------------------------------------------------------------------------------------------------------

  ! ... variables used to store the reduced sum
  INTEGER,     ALLOCATABLE, DIMENSION(:, :)  :: nInsertSegment, nInsertAntiSegment, nRemoveSegment, nRemoveAntiSegment, nShiftSegment, nEmptyFull

  REAL(DP),    ALLOCATABLE :: Gt_Bin_Reduce(:, :, :, :), Gt_Bin(:, :, :, :), MCsign_Bin(:, :), MCsign_Reduce(:, :)
  REAL(DP),    ALLOCATABLE :: Gst_Bin(:, :, :, :), Gst_Bin_Reduce(:, :, :, :)
  REAL(DP),    ALLOCATABLE :: nk_Probability(:, :, :)

CONTAINS
  !-------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Segment_Update
    
    IMPLICIT NONE
    
    !
    ! Purposes
    ! ========
    !
    !
    ! Code developer: Gang Li
    ! Code history:  2021.1.6 adapted from code June 2006
    !
    
    ! ... local vars ...
    INTEGER  :: iMeasure, i, iBin
    INTEGER  :: iOrbit, iSpin
    
    ! ... executable ...
    CALL Segment_MC_Init
    Gt  = Zero
    Gst = Zero
    GW  = Zero
    
    DO iMeasure = 1, nWarm
       ! >>> warm up to thermal equilibrium
       DO iOrbit = 1, nOrbit
          DO iSpin = 1, nSpin
             CALL MonteCarlo_WorkFlow(iOrbit, iSpin)
          END DO
       END DO
    END DO
    
    IF (id == master) WRITE(*, *) ' >>> warm up is done! starting to accumulate results'
    
    ! >>> initialize the variables and prepare to reduce from all nodes
    CALL Segment_Phys_Initialize
    
    ! >>> start to measure imaginary-time Green's function and equal-time correlators
    DO iBin = 1, nBin
       Gt         = Zero
       Gst        = Zero
       MCsign_Bin(:, iBin) = Zero
       DO iMeasure = 1, nMeasure
          
          DO iOrbit = 1, nOrbit
             DO iSpin = 1, nSpin
                
                CALL MonteCarlo_WorkFlow(iOrbit, iSpin)
                
                ! >>> measure the imaginary-time Green's function
                CALL Segment_MeasureGt(iOrbit, iSpin)
                ! >>> determine the average sign of configurations
                MCsign_Bin(1:2, iBin) = MCSign_Bin(1:2, iBin) + Segment_Statistic%MCsign(1:2)/DBLE(nOrbit*nSpin*nMeasure)
                
             END DO  ! do iSpin = 1, nSpin
          END DO  ! do iOrbit = 1, nOrbit
          
       END DO
       
       ! monitoring G(tau = beta) on standard output, which should be -1/2 for half-filled Hubbard model
       DO iOrbit = 1, nOrbit
          IF (id == master) WRITE(*, '(5x, a, i3, a, i3, a, 2f12.6)') 'iBin =', iBin, ' Orbital=', iOrbit, ', G(Beta, iOrbit)=', Gt(nTau, iOrbit, 1:nSpin)
       END DO

       ! special treatement of Gt(0) and Gt(beta)
       DO iOrbit = 1, nOrbit
          DO iSpin = 1, nSpin
             Gt(0, iOrbit, iSpin)                 = -One - Gt(nTau, iOrbit, iSpin)
             Gt_Bin(0:nTau, iOrbit, iSpin, iBin)  = Gt(0:nTau, iOrbit, iSpin)
             Gst_Bin(0:ntau, iOrbit, iSpin, iBin) = Gst(0:ntau, iOrbit, iSpin)
          END DO
       END DO
       
    END DO

    ! reduce statistical message to the master node
    Gt_Bin_Reduce  = Zero    
    Gst_Bin_Reduce = Zero
    CALL MPI_REDUCE(Gt_Bin, Gt_Bin_Reduce, (nTau+1)*nOrbit*nSpin*nBin,        MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, rc)
    CALL MPI_REDUCE(Gst_Bin, Gst_Bin_Reduce, (nTau+1)*nOrbit*nSpin*nBin,      MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, rc)
    CALL MPI_REDUCE(Segment_Statistic%nInsertSegment,     nInsertSegment,     nOrbit*nSpin, MPI_INTEGER, MPI_SUM, master, MPI_COMM_WORLD, rc)
    CALL MPI_REDUCE(Segment_Statistic%nInsertAntiSegment, nInsertAntiSegment, nOrbit*nSpin, MPI_INTEGER, MPI_SUM, master, MPI_COMM_WORLD, rc)
    CALL MPI_REDUCE(Segment_Statistic%nRemoveSegment,     nRemoveSegment,     nOrbit*nSpin, MPI_INTEGER, MPI_SUM, master, MPI_COMM_WORLD, rc)
    CALL MPI_REDUCE(Segment_Statistic%nRemoveAntiSegment, nRemoveAntiSegment, nOrbit*nSpin, MPI_INTEGER, MPI_SUM, master, MPI_COMM_WORLD, rc)
    CALL MPI_REDUCE(Segment_Statistic%nShiftSegment,      nShiftSegment,      nOrbit*nSpin, MPI_INTEGER, MPI_SUM, master, MPI_COMM_WORLD, rc)
    CALL MPI_REDUCE(Segment_Statistic%nEmptyFull,         nEmptyFull,         nOrbit*nSpin, MPI_INTEGER, MPI_SUM, master, MPI_COMM_WORLD, rc)
    CALL MPI_REDUCE(Segment_Statistic%nk_Probability,     nk_Probability,    (nMax+1)*nOrbit*nSpin, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, rc)
    CALL MPI_REDUCE(MCsign_Bin,                           MCsign_Reduce,      2*nBin,               MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, rc)
    
    IF (id == master) THEN

       nk_Probability = nk_Probability/DBLE(ntasks)
       MCsign_Reduce  = MCsign_Reduce/DBLE(ntasks)       
       OPEN(UNIT=1, FILE='MC_Status.dat', STATUS='unknown')
       
       WRITE(1, '(3x, a50, 100(f8.2, "%"))') ' # of successful inserting segment',      (DBLE(nInsertSegment(iOrbit, 1:nSpin)*100)/DBLE(nMeasure*nBin*ntasks), iOrbit=1, nOrbit)
       WRITE(1, '(3x, a50, 100(f8.2, "%"))') ' # of successful inserting antisegment',  (DBLE(nInsertAntiSegment(iOrbit, 1:nSpin)*100)/DBLE(nMeasure*nBin*ntasks), iOrbit=1, nOrbit)
       WRITE(1, '(3x, a50, 100(f8.2, "%"))') ' # of successful removing segment',       (DBLE(nRemoveSegment(iOrbit, 1:nSpin)*100)/DBLE(nMeasure*nBin*ntasks), iOrbit=1, nOrbit)
       WRITE(1, '(3x, a50, 100(f8.2, "%"))') ' # of successful removing antisegment',   (DBLE(nRemoveAntiSegment(iOrbit, 1:nSpin)*100)/DBLE(nMeasure*nBin*ntasks), iOrbit=1, nOrbit)
       WRITE(1, '(3x, a50, 100(f8.2, "%"))') ' # of successful shifting segment',       (DBLE(nShiftSegment(iOrbit, 1:nSpin)*100)/DBLE(nMeasure*nBin*ntasks), iOrbit=1, nOrbit)
       WRITE(1, '(3x, a50, 100(f8.2, "%"))') ' # of successful switching empty-full segment',  (DBLE(nEmptyFull(iOrbit, 1:nSpin)*100)/DBLE(nMeasure*nBin*ntasks), iOrbit=1, nOrbit)
       WRITE(1, '(3x, a50, 2f12.6)') ' average sign for all and accepted configurations', SUM(MCsign_Reduce(1, 1:nBin))/DBLE(nBin), SUM(MCsign_Reduce(2, 1:nBin))/DBLE(nBin)
       WRITE(1, *)
       WRITE(1, *) '# histogram'
       DO i = 0, nMax
          WRITE(1, '(i5, 100f12.6)') i, (nk_Probability(i, iOrbit, 1:nSpin), iOrbit = 1, nOrbit)
       END DO
       CLOSE(1)

       CALL Segment_CalGtw
       
    END IF
    
    CALL MPI_BCAST(Sigma(1:nOmega, 1:nOrbit, 1:nSpin), nOmega*nOrbit*nSpin, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    CALL MPI_BCAST(ImpurityOccupancy(1:nOrbit, 1:nSpin), nOrbit*nSpin, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, rc)

    ! >>> clear memory
    CALL Segment_Phys_Finalize
    
    CALL Segment_MC_Finalize
    
  END SUBROUTINE Segment_Update

  !--------------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Segment_Phys_Initialize

    IMPLICIT NONE

    ! >>> initialize the accumulators for all measurements
    
    Segment_Statistic%nInsertSegment     = 0
    Segment_Statistic%nInsertAntiSegment = 0
    Segment_Statistic%nRemoveSegment     = 0
    Segment_Statistic%nRemoveAntiSegment = 0
    Segment_Statistic%nShiftSegment      = 0
    Segment_Statistic%nEmptyFull         = 0
    Segment_Statistic%MCsign               = Zero

    ALLOCATE(nInsertSegment(nOrbit, nSpin))
    ALLOCATE(nInsertAntiSegment(nOrbit, nSpin))
    ALLOCATE(nRemoveSegment(nOrbit, nSpin))
    ALLOCATE(nRemoveAntiSegment(nOrbit, nSpin))
    ALLOCATE(nShiftSegment(nOrbit, nSpin))
    ALLOCATE(nEmptyFull(nOrbit, nSpin))
    nInsertSegment     = 0
    nInsertAntiSegment = 0
    nRemoveSegment     = 0
    nRemoveAntiSegment = 0
    nShiftSegment      = 0
    nEmptyFull         = 0
    
    ALLOCATE(Gt_Bin(0:nTau, nOrbit, nSpin, nBin))
    ALLOCATE(Gt_Bin_Reduce(0:nTau, nOrbit, nSpin, nBin))
    ALLOCATE(Gst_Bin(0:nTau, nOrbit, nSpin, nBin))
    ALLOCaTE(Gst_Bin_Reduce(0:nTau, nOrbit, nSpin, nBin))
    
    ALLOCATE(nk_Probability(0:nMax, nOrbit, nSpin))
    nk_Probability = Zero
    
    ALLOCATE(MCsign_Bin(2, nBin))
    ALLOCATE(MCsign_Reduce(2, nBin))
    
  END SUBROUTINE Segment_Phys_Initialize
  
  !--------------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Segment_Phys_Finalize

    DEALLOCATE(nInsertSegment)
    DEALLOCATE(nInsertAntiSegment)
    DEALLOCATE(nRemoveSegment)
    DEALLOCATE(nRemoveAntiSegment)
    DEALLOCATE(nShiftSegment)
    DEALLOCATE(nEmptyFull)

    DEALLOCATE(Gt_Bin)
    DEALLOCATE(Gt_Bin_Reduce)
    DEALLOCATE(Gst_Bin)
    DEALLOCATE(Gst_Bin_Reduce)
    
    DEALLOCATE(nk_Probability)

    DEALLOCATE(MCsign_Bin)
    DEALLOCATE(MCsign_Reduce)

  END SUBROUTINE Segment_Phys_Finalize

  !--------------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Segment_MeasureGt(iOrbit, iSpin)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iOrbit, iSpin

    ! ... local vars ...
    INTEGER  :: i, j, Idx_Tau, i1, j1 
    REAL(DP) :: Tau !, Length
    
    ! loop over all accepted segment and calculate the imaginary-time Green's function
    DO i = 1, Kink(iOrbit, iSpin)%nk
       DO j = 1, Kink(iOrbit, iSpin)%nk
          
          Tau = Kink(iOrbit, iSpin)%Tau_END(i) - Kink(iOrbit, iSpin)%Tau_Start(j)  
          IF (Tau >= Zero) THEN
             Idx_Tau = NINT(Tau/Beta*nTau)
             IF (Idx_Tau /= nTau .AND. Idx_Tau /= 0) THEN
                ! >>> standard estimator
                Gt(Idx_Tau, iOrbit, iSpin)  =  Gt(Idx_Tau, iOrbit, iSpin) &
                     - Kink(iOrbit, iSpin)%Det(j, i)/Beta*DBLE(nTau)/Beta/DBLE(nMeasure)

                ! >>> improved estimator
                DO i1 = 1, nOrbit
                   DO j1 = 1, nSpin
                      
                      IF (i1 /= iOrbit .OR. j1 /= iSpin) THEN
                         Gst(Idx_Tau, iOrbit, iSpin) = Gst(Idx_Tau, iOrbit, iSpin) - Half*( Umat(iOrbit, iSpin, i1, j1) + Umat(i1, j1, iOrbit, iSpin) ) &
                              * Kink(iOrbit, iSpin)%n_Start(i, i1, j1) * Kink(iOrbit, iSpin)%Det(j, i)*Segment_Statistic%MCSign(2)/Beta*DBLE(nTau)/Beta/DBLE(nMeasure)
                      END IF
                      
                   END DO
                END DO
            END IF
             
          ELSE
             Tau = Tau + Beta
             Idx_Tau = NINT(Tau/Beta*nTau)
             IF (Idx_Tau /= nTau .AND. Idx_Tau /= 0) THEN
                ! >>> standard estimator
                Gt(Idx_Tau, iOrbit, iSpin)  =  Gt(Idx_Tau, iOrbit, iSpin) &
                     + Kink(iOrbit, iSpin)%Det(j, i)/Beta*DBLE(nTau)/Beta/DBLE(nMeasure)
             
                ! >>> improved estimator
                DO i1 = 1, nOrbit
                   DO j1 = 1, nSpin
                      IF (i1 /= iOrbit .OR. j1 /= iSpin) THEN
                         Gst(Idx_Tau, iOrbit, iSpin) = Gst(Idx_Tau, iOrbit, iSpin) + Half*( Umat(iOrbit, iSpin, i1, j1) + Umat(i1, j1, iOrbit, iSpin) ) &
                              * Kink(iOrbit, iSpin)%n_Start(i, i1, j1) * Kink(iOrbit, iSpin)%Det(j, i)*Segment_Statistic%MCSign(2)/Beta*DBLE(nTau)/Beta/DBLE(nMeasure)
                      END IF
                  END DO
                END DO
                
             END IF

          END IF
          
       END DO
    END DO

    ! >>> G(beta) corresponds to the average occupancy, which is simply measured from the total length of segments
    Gt(nTau, iOrbit, iSpin) = Gt(nTau, iOrbit, iSpin) - Kink(iOrbit, iSpin)%Total_Length*Segment_Statistic%MCSign(2)/Beta/DBLE(nMeasure)

    ! >>> Gst(tau=0) and Gst(tau=nTau) are not directly measured, we simply linearly extrapolate from neighboring points
    Gst(0, iOrbit, iSpin)    = Two*Gst(1, iOrbit, iSpin) - Gst(2, iOrbit, iSpin)
    Gst(nTau, iOrbit, iSpin) = Two*Gst(nTau-1, iOrbit, iSpin) - Gst(nTau-2, iOrbit, iSpin) 
    
    ! probability of configurations with the same number of segments p(k)
    Segment_Statistic%nk_Probability(Kink(iOrbit, iSpin)%nk, iOrbit, iSpin) =  Segment_Statistic%nk_Probability(Kink(iOrbit, iSpin)%nk, iOrbit, iSpin) + One/DBLE(nMeasure*nBin)
    
  END SUBROUTINE Segment_MeasureGt

  !---------------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Segment_CalGtw

    IMPLICIT NONE

    ! ... local vars ...
    INTEGER               :: iOrbit, iSpin, iTau, iBin
    INTEGER               :: i, j, nLegenPolyOptimal(nOrbit, nSpin), n_smooth
    REAL(DP)              :: Gt_Error(0:nTau, nOrbit, nSpin), Gt_Temp, Chi2, Chi2_min, Temp(0:nTau)
    CHARACTER(Len=30)     :: FLE
    
    ! >>> don't forget to average over nodes
    Gt_Bin_Reduce  = Gt_Bin_Reduce/DBLE(ntasks)
    Gst_Bin_Reduce = Gst_Bin_Reduce/DBLE(ntasks)
    
    ! >>> statistical error analysis
    Gt       = Zero
    Gst      = Zero
    Gt_Error = Zero
    
    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          DO iTau = 0, nTau
             DO iBin = 1, nBin
                Gt(iTau, iOrbit, iSpin)  = Gt(iTau, iOrbit, iSpin)  + Gt_Bin_Reduce(iTau, iOrbit, iSpin, iBin)/DBLE(nBin)/MCsign_Bin(2, iBin)
                Gst(iTau, iOrbit, iSpin) = Gst(iTau, iOrbit, iSpin) + Gst_Bin_Reduce(iTau, iOrbit, iSpin, iBin)/DBLE(nBin)/MCsign_Bin(2, iBin) 
                Gt_Error(iTau, iOrbit, iSpin) = Gt_Error(iTau, iOrbit, iSpin) + Gt_Bin_Reduce(iTau, iOrbit, iSpin, iBin)**2/DBLE(nBin)/MCsign_Bin(2, iBin)
             END DO
             Gt_Error(iTau, iOrbit, iSpin) = sqrt(ABS(Gt_Error(iTau, iOrbit, iSpin) - Gt(iTau, iOrbit, iSpin)**2)/DBLE(nBin-1))
          END DO
       END DO
    END DO

    ! >>> smooth the raw data by averaging over neighboring 12 points
    n_smooth=30
    do iOrbit = 1, nOrbit
       do iSpin = 1, nSpin
          Temp = Zero
          do iTau = nTau/3, 2*nTau/3
             do i = 1, n_smooth
                Temp(iTau) = Temp(iTau) + (Gt(iTau-i, iOrbit, iSpin) + Gt(iTau+i, iOrbit, iSpin))/Two
             end do   
             Temp(iTau)=Temp(iTau)/DBLE(n_smooth)
          end do
          Gt(nTau/3:2*nTau/3, iOrbit, iSpin) = Temp(nTau/3:2*nTau/3)
          Temp = Zero
          do iTau = nTau/3, 2*nTau/3
             do i = 1, n_smooth
                Temp(iTau) = Temp(iTau) + (Gst(iTau-i, iOrbit, iSpin) + Gst(iTau+i, iOrbit, iSpin))/Two
             end do
             Temp(iTau)=Temp(iTau)/DBLE(n_smooth)
          end do
          Gst(nTau/3:2*nTau/3, iOrbit, iSpin) = Temp(nTau/3:2*nTau/3)
       end do
    end do

    ! >>> store the impurity occupancy from Gt(nTau)
    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          ImpurityOccupancy(iOrbit, iSpin) = -Gt(nTau, iOrbit, iSpin)
       END DO
    END DO 

    ! >>> write raw imaginary-time Green's function to file
    WRITE(FLE, '(a)') 'GtRaw.dat'
    OPEN(UNIT=1, FILE=TRIM(FLE), STATUS='UNKNOWN')
    WRITE(1, '(a)') '# Tau   Gt(tau, iOrbit, iSpin=1)   Gt(tau, iOrbit, iSpin=2)  Gt_Error(tau, iOrbit, iSpin=1)  Gt_Error(tau, iOrbit, iSpin=2)'
    DO i = 0, nTau
       WRITE(1, '(100f21.15)') Beta/DBLE(nTau)*DBLE(i), (Gt(i, 1:nOrbit, iSpin), iSpin=1, 2)
    END DO
    CLOSE(1)

    WRITE(FLE, '(a)') 'GtRaw_error.dat'
    OPEN(UNIT=1, FILE=TRIM(FLE), STATUS='UNKNOWN')
    WRITE(1, '(a)') '# Tau   Gt(tau, iOrbit, iSpin=1)   Gt(tau, iOrbit, iSpin=2)  Gt_Error(tau, iOrbit, iSpin=1)  Gt_Error(tau, iOrbit, iSpin=2)'
    DO i = 0, nTau
       WRITE(1, '(100f21.15)') Beta/DBLE(nTau)*DBLE(i), (Gt_Error(i, 1:nOrbit, iSpin), iSpin=1, 2)
    END DO
    CLOSE(1)
    
    ! >>> Legendre Polynomail filtering
    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          CALL Coefficient_LegendrePolynomials(nTau, nLegenPoly, Beta, Gt(0:nTau, iOrbit, iSpin), LegenPoly(0:nLegenPoly, 0:nTau), CoeffLegenPoly(0:nLegenPoly))

          nLegenPolyOptimal(iOrbit, iSpin) = nLegenPoly
          ! determine the optimal order of Legendre Polynomial
          Chi2_min = 100.d0
          DO j = 20, nLegenPoly
             chi2 = Zero
             DO iTau = 0, nTau
                Gt_Temp = Zero
                DO i = 0, j
                   Gt_Temp = Gt_Temp + CoeffLegenPoly(i)*LegenPoly(i, iTau)
                END DO
                chi2 = chi2 + ABS(Gt(iTau, iOrbit, iSpin) - Gt_Temp)/DBLE(nTau+1)
             END DO
             IF (chi2 < Chi2_min) THEN
                Chi2_min = Chi2
                nLegenPolyOptimal(iOrbit, iSpin) = j
             END IF
          END DO

          ! >>> calculate the imaginary-time Green's function
          DO iTau = 0, nTau             
             Gt(iTau, iOrbit, iSpin) =  Zero
             DO j = 0, nLegenPolyOptimal(iOrbit, iSpin)
                Gt(iTau, iOrbit, iSpin) = Gt(iTau, iOrbit, iSpin) + CoeffLegenPoly(j)*LegenPoly(j, iTau)
             END DO
          END DO

          ! >>> calculate the Matsubara Green's function
          DO i = 1, nOmega
             Gw(i, iOrbit, iSpin) = Zero
             DO j = 0, nLegenPolyOptimal(iOrbit, iSpin)
                Gw(i, iOrbit, iSpin) = Gw(i, iOrbit, iSpin) + CoeffLegenPoly(j)*LegenPolyMatsubara(j, i)
             END DO
          END DO

       END DO
    END DO
    CLOSE(1)
     
    ! >>> output the Matsubara Green's function calculated from improved estimatorLegendre polynomial fit
    WRITE(FLE, '(a)') 'Gw.dat'
    OPEN(unit=1, file=TRIM(FLE), STATUS='UNKNOWN')
    WRITE(1, '(a,100i5)') '# Optimal value of Legendre Polynomial', (nLegenPolyOptimal(1:nOrbit, iSpin), iSpin=1, 2)
    DO i = 1, nOmega
       WRITE(1, '(100f21.15)') IMAG(Omega(i)), (Gw(i, 1:nOrbit, iSpin), iSpin=1, 2)
    END DO
    CLOSE(1)
    
    ! >>> write the smoothed imaginary-time Green's function to file
    WRITE(FLE, '(a)') 'GtLegenPoly.dat'
    OPEN(UNIT=1, FILE=TRIM(FLE), STATUS='UNKNOWN')
    WRITE(1, '(a, 100i5)') '# Optimal value of Legendre Polynomial', (nLegenPolyOptimal(1:nOrbit, iSpin), iSpin=1, 2)
    DO i = 0, nTau
       WRITE(1, '(100f21.15)') Beta/DBLE(nTau)*DBLE(i), (Gt(i, 1:nOrbit, iSpin), iSpin=1, 2)
    END DO
    CLOSE(1)

    ! >>> improved self-energy estmator
    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          CALL Coefficient_LegendrePolynomials(nTau, nLegenPoly, Beta, Gst(0:nTau, iOrbit, iSpin), LegenPoly(0:nLegenPoly, 0:nTau), CoeffLegenPoly(0:nLegenPoly))

          nLegenPolyOptimal(iOrbit, iSpin) = nLegenPoly
          ! determine the optimal order of Legendre Polynomial
          Chi2_min = 100.d0
          DO j = 20, nLegenPoly
             chi2 = Zero
             DO iTau = 0, nTau
                Gt_Temp = Zero
                DO i = 0, j
                   Gt_Temp = Gt_Temp + CoeffLegenPoly(i)*LegenPoly(i, iTau)
                END DO
                chi2 = chi2 + ABS(Gst(iTau, iOrbit, iSpin) - Gt_Temp)/DBLE(nTau+1)
             END DO
             IF (chi2 < Chi2_min) THEN
                Chi2_min = Chi2
                nLegenPolyOptimal(iOrbit, iSpin) = j
             END IF
          END DO

          ! >>> calculate the Matsubara self-energy function, substracted with Hartreen energy
          CALL Legendre_Polynomials_Matsubara(nOmega, nLegenPolyOptimal(iOrbit, iSpin), Beta, LegenPolyMatsubara(0:nLegenPolyOptimal(iOrbit, iSpin), 1:nOmega))
          ! only the first 1/2 is directly measured, the last 1/2 will be replaced by the second-order Self-energy
          DO i = 1, nOmega
             Sigma(i, iOrbit, iSpin) = Zero   ! here we use Sigma to temporarily store F(iW_n)
             DO j = 0, nLegenPolyOptimal(iOrbit, iSpin)
                Sigma(i, iOrbit, iSpin) = Sigma(i, iOrbit, iSpin) + CoeffLegenPoly(j)*LegenPolyMatsubara(j, i)
             END DO
          END DO
       END DO
    END DO

    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          DO i = 1, nOmega
             Gw(i, iOrbit, iSpin) = G0w(i, iOrbit, iSpin) * (One + Sigma(i, iOrbit, iSpin))  ! &
!                   / (One + G0w(i, iOrbit, iSpin)*Half*(xU + (nOrbit-1)*(xU - 2*xJ) + (nOrbit-1)*(xU - 3*xJ)))
             Sigma(i, iOrbit, iSpin) = Sigma(i, iOrbit, iSpin)/Gw(i, iOrbit, iSpin) !- Half*(xU + (nOrbit-1)*(xU - 2*xJ) + (nOrbit-1)*(xU - 3*xJ))
          END DO
       END DO
    END DO
    
    ! >>> output the self-energy function from improved estimator
    WRITE(FLE, '(a, i1, a)') 'Sigma.dat'
    
    OPEN(unit=1, file=TRIM(FLE), STATUS='UNKNOWN')
    WRITE(1, '(a, 100i5)') '# Optimal value of Legendre Polynomial', (nLegenPolyOptimal(1:nOrbit, iSpin), iSpin=1, 2)
    DO i = 1, nOmega
       WRITE(1, '(100f21.15)') IMAG(Omega(i)), (Sigma(i, 1:nOrbit, iSpin), iSpin=1, 2)
    END DO
    CLOSE(1)        

    ! >>> output the Matsubara Green's function calculated from improved estimatorLegendre polynomial fit
    WRITE(FLE, '(a)') 'GwP.dat'
    OPEN(unit=1, file=TRIM(FLE), STATUS='UNKNOWN')
    WRITE(1, '(a, 100i5)') '# Optimal value of Legendre Polynomial', (nLegenPolyOptimal(1:nOrbit, iSpin), iSpin=1, 2)
    DO i = 1, nOmega
       WRITE(1, '(100f21.15)') IMAG(Omega(i)), (Gw(i, 1:nOrbit, iSpin), iSpin=1, 2)
    END DO
    CLOSE(1)
    
  END SUBROUTINE Segment_CalGtw
  
END MODULE Segment_Phys
