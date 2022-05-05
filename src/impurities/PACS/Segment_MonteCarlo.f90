!******************************COPYRIGHT************************************************
! (c) Crown copyright, Gang Li < 2020 >. All rights reserved.
!
! This routine has been licensed to the valid partners for use and distribution
! under the collaboration agreement, subject to the terms and conditions set out
! therein.
!
! [ligang at Shanghaitech.edu.cn]
!******************************COPYRIGHT************************************************
  
MODULE Segment_MonteCarlo

  USE Segment_Util

  IMPLICIT NONE

  !---------------------------------------------------------------------------------------
  ! Description:
  !
  !   Module Segment_MonteCarlo performs the Monte Carlo simulation by manipulating the segment,
  !   which is the hybridization function in imaginary-time space.
  !
  ! Code hisotry: 2021.1.5 adapted from code June 2007
  !
  ! Current Code Owner: Gang Li
  !
  ! Code Description:
  !   Language: Fortran 90.
  !---------------------------------------------------------------------------------------

  INTEGER, PARAMETER   :: nMax = 256       ! maximum number of kinks possibly added into [0, beta)

  TYPE Segment
     INTEGER  :: nk
     REAL(DP) :: Tau_Start(nMax)
     REAL(DP) :: Tau_End(nMax)
     REAL(DP) :: Det(nMax, nMax)
     REAL(DP) :: Total_Length
     LOGICAL  :: EmptySegment         ! only make sense when nk == 0. .True. for empty segment,
                                      !                               .False. for full segment
     INTEGER, ALLOCATABLE :: n_Start(:, :, :)   ! if Tau_Start is occupied in other flavors. 
                                                ! 1 for occupied case
                                                ! 0 for unoccupied case
  END type Segment
    
  TYPE(Segment), DIMENSION(:, :), ALLOCATABLE :: Kink

  TYPE MC_Statistic
     INTEGER, ALLOCATABLE  :: nInsertSegment(:, :)
     INTEGER, ALLOCATABLE  :: nInsertAntiSegment(:, :)
     INTEGER, ALLOCATABLE  :: nRemoveSegment(:, :)
     INTEGER, ALLOCATABLE  :: nRemoveAntiSegment(:, :)
     INTEGER, ALLOCATABLE  :: nShiftSegment(:, :)
     INTEGER, ALLOCATABLE  :: nEmptyFull(:, :)
     REAL(DP), ALLOCATABLE :: nk_Probability(:, :, :)
     
     ! We will measure two signs, Sign(1) records all negative sign configurations; Sign(2) records only
     ! those accepted updates with ranw() < ABS(Update_Ratio) 
     REAL(DP) :: MCsign(2)                  
                                                        
  END type MC_Statistic

  TYPE(MC_Statistic) :: Segment_Statistic

  PRIVATE :: Insert_Kink, Remove_Kink, Shift_Kink, Empty_Full, Segment_Overlap
  
CONTAINS
  !---------------------------------------------------------------------------------------
  SUBROUTINE Segment_MC_Init

    IMPLICIT NONE

    !
    ! Purposes
    ! ========
    !
    !
    ! Code developer: Gang Li
    ! Code history:  2021.1.6 adapted from code June 2006
    !
    INTEGER :: i, j
    
    ALLOCATE(Kink(nOrbit, nSpin))
    Kink(1:nOrbit, 1:nSpin)%nk               = 0
    DO i = 1, nMax
       Kink(1:nOrbit, 1:nSpin)%Tau_Start(i)  = Zero
       Kink(1:nOrbit, 1:nSpin)%Tau_End(i)    = Zero
       DO j = 1, nMax
          Kink(1:nOrbit, 1:nSpin)%Det(i, j)  = Zero          
       END DO
    END DO

    Kink(1:nOrbit, 1:nSpin)%EmptySegment     = .True.
    Kink(1:nOrbit, 1:nSpin)%Total_Length     = Zero

    DO i = 1, nOrbit
       DO j = 1, nSpin
          ALLOCATE(Kink(i, j)%n_Start(nMax, nOrbit, nSpin))
          Kink(i, j)%n_Start(1:nMax, 1:nOrbit, 1:nSpin) = 0
       END DO
    END DO  

    ALLOCATE(Segment_Statistic%nInsertSegment(nOrbit, nSpin))
    ALLOCATE(Segment_Statistic%nInsertAntiSegment(nOrbit, nSpin))
    ALLOCATE(Segment_Statistic%nRemoveSegment(nOrbit, nSpin))
    ALLOCATE(Segment_Statistic%nRemoveAntiSegment(nOrbit, nSpin))
    ALLOCATE(Segment_Statistic%nShiftSegment(nOrbit, nSpin))
    ALLOCATE(Segment_Statistic%nEmptyFull(nOrbit, nSpin))
    Segment_Statistic%nInsertSegment     = 0
    Segment_Statistic%nInsertAntiSegment = 0
    Segment_Statistic%nRemoveSegment     = 0
    Segment_Statistic%nRemoveAntiSegment = 0
    Segment_Statistic%nShiftSegment      = 0
    Segment_Statistic%nEmptyFull         = 0
    Segment_Statistic%MCsign             = Zero

    ALLOCATE(Segment_Statistic%nk_Probability(0:nMax, nOrbit, nSpin))
    Segment_Statistic%nk_Probability     = Zero
    
  END SUBROUTINE Segment_MC_Init

  !---------------------------------------------------------------------------------------
  SUBROUTINE Segment_MC_Finalize

    IMPLICIT NONE

    !
    ! Purposes
    ! ========
    !
    !
    ! Code developer: Gang Li
    ! Code history:  2021.1.6 adapted from code June 2006
    !

    DEALLOCATE(Kink)
    DEALLOCATE(Segment_Statistic%nInsertSegment)
    DEALLOCATE(Segment_Statistic%nInsertAntiSegment)
    DEALLOCATE(Segment_Statistic%nRemoveSegment)
    DEALLOCATE(Segment_Statistic%nRemoveAntiSegment)
    DEALLOCATE(Segment_Statistic%nShiftSegment)
    DEALLOCATE(Segment_Statistic%nEmptyFull)
    DEALLOCATE(Segment_Statistic%nk_Probability)
    
  END SUBROUTINE Segment_MC_Finalize

  !---------------------------------------------------------------------------------------
  SUBROUTINE MonteCarlo_WorkFlow(Idx_Orbit, Idx_Spin)

    INTEGER, INTENT(IN) :: Idx_Orbit, Idx_Spin

    ! ... local vars ...
    REAL(DP) :: rand
    
    ! --- Monte Carlo sweep ---
    rand = ranw()
    IF (rand > 0.75d0) THEN
       CALL Insert_Kink(Idx_Orbit, Idx_Spin)
    ELSEIF (rand > 0.5d0) THEN
       CALL Remove_Kink(Idx_Orbit, Idx_Spin)
    ELSE IF (rand > 0.25d0) THEN
      CALL Shift_Kink(Idx_Orbit, Idx_Spin)
    ELSE
      CALL Empty_Full(Idx_Orbit, Idx_Spin)
    END IF
    ! --- End Monte Carlo sweep ---
    
  END SUBROUTINE MonteCarlo_WorkFlow
  
  !---------------------------------------------------------------------------------------
  SUBROUTINE Insert_Kink(Idx_Orbit, Idx_Spin)

    IMPLICIT NONE

    !
    ! Purposes
    ! ========
    !
    INTEGER, INTENT(IN)  :: Idx_Orbit
    INTEGER, INTENT(IN)  :: Idx_Spin
    
    ! ... local vars ...
    INTEGER  :: Idx_Ts, Idx_Te                         ! the location of the inserted kink in the 'nk' kink list
    INTEGER  :: Idx_Start(2), Idx_End(2)
    REAL(DP) :: Ts, Te                                 ! the starting and ending point of the inserted kink
    REAL(DP) :: Tau_Max, Tau_Length, Tau
    REAL(DP) :: Update_Ratio, Det_Ratio, Temp
    REAL(DP) :: Overlap(nOrbit, nSpin)
    REAL(DP), ALLOCATABLE :: Gi(:), Gj(:), XL(:), XR(:)  ! We will only use part of these matrices to avoid frequently allocating memory
    
    LOGICAL  :: LSegment = .TRUE.
    INTEGER  :: nk
    INTEGER  :: i, j
    INTEGER  :: iOrbit, iSpin, Occupancy_temp(nOrbit, nSpin)
     
    ! ... exectuable ...
    nk = Kink(Idx_Orbit, Idx_Spin)%nk           ! get the current order

    ! randomly chose a starting point
    Ts         = ranw()*Beta
    Te         = Zero
    Idx_Ts     = 0
    Idx_Start  = 0
    Idx_End    = 0
    Tau_Max    = Zero
    Tau_Length = Zero

    ! STEP[1]: randomly add a segment or an antisegment.
    IF (ranw() > Half) THEN
       
       ! >>> try to insert new segment
       LSegment = .TRUE.
       
       IF (nk == 0) THEN  ! when nk = 0, be causion with the full segment
          ! --- start block ---
          IF (Kink(Idx_Orbit, Idx_Spin)%EmptySegment) THEN
             Te              = ranw()*Beta
             Tau_Length      = Te - Ts
             Tau_Max         = Beta
             IF (Tau_Length < Zero) Tau_Length = Tau_Length + Beta          
          END IF
          Idx_Ts    = 1
          Idx_Start = (/nk+1, 1/)
          Idx_End   = (/nk+1, 1/)
          ! --- end block ---
          
       ELSE  ! -> IF (nk /= 0)
          
          Idx_Ts = 0
          
          DO i = 2, nk
             IF (Ts < Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i) .AND. Ts > Kink(Idx_Orbit, Idx_Spin)%Tau_End(i-1)) THEN
                Idx_Ts     = i
                Idx_Te     = i
                Tau_Max    = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i) - Ts
                Te         = ranw()*Tau_Max + Ts
                Tau_Length = Te - Ts
             END IF
          END DO
          
          IF (Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) < Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1)) THEN
             IF (Ts < Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1) .AND. Ts > Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)) THEN
                Idx_Ts     = 1
                Idx_Te     = 1
                Tau_Max    = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1) - Ts
                Te         = ranw()*Tau_Max + Ts
                Tau_Length = Te - Ts
             END IF
          ELSE
             IF (Ts < Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1)) THEN
                Idx_Ts     = 1
                Idx_Te     = 1
                Tau_Max    = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1) - Ts
                Te         = ranw()*Tau_Max + Ts
                Tau_Length = Te - Ts
             END IF
             IF (Ts > Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) ) THEN
                Idx_Ts     = nk+1
                Idx_Te     = nk+1
                Tau_Max    = Beta + Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1) - Ts
                Te         = ranw()*Tau_Max + Ts
                Tau_Length = Te - Ts
                IF (Te > Beta) Te = Te - Beta
             END IF
          END IF
          
          Idx_Start = (/nk+1, Idx_Ts/)
          Idx_End   = (/nk+1, Idx_Te/)
       END IF  ! -> IF (nk == 0)

       IF (DEBUG) WRITE(*, '(3x, a, 2i5, 2f12.6)')  'Runtime Info: insert segment', Idx_Orbit, Idx_Spin, Ts, Te
       
    ELSE
              
       ! try to insert new antisegment
       LSegment = .FALSE.
       Idx_Ts = 0
       
       IF (nk == 0) THEN
          IF (.NOT. Kink(Idx_Orbit, Idx_Spin)%EmptySegment) THEN
             Tau_Max    = Beta
             Te         = ranw()*Beta
             Tau_Length = Te - Ts
             IF (Tau_Length < Zero) Tau_Length = Beta + Tau_Length
             Idx_Ts    = 1
             Idx_Start = (/nk+1, 1/)
             Idx_End   = (/nk+1, 1/)
          END IF
       ELSE
          
          Idx_Ts = 0
          
          IF (Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) < Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1)) THEN
             IF (Ts > Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk)) THEN
                ! >>> It will be an anti-segment breaking the last segment in the kink list
                Tau_Max      = Beta + Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) - Ts
                Te           = ranw()*Tau_Max + Ts
                Tau_Length   = Te - Ts
                Idx_Ts       = nk
                IF (Te > Beta) THEN
                   Te        = Te - Beta
                   Idx_Start = (/nk+1, 1/)
                   Idx_End   = (/nk, 1/)
                ELSE
                   Idx_Start = (/nk+1, nk+1/)
                   Idx_End   = (/nk+1, nk/)
                END IF
             END IF
             IF (Ts < Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)) THEN
                ! >>> It will be an anti-segment breaking the last segment in the kink list, and adding a new segment in the begining
                ! >>> of the kink list
                Tau_Max     = Kink(Idx_Orbit, Idx_Spin)%Tau_end(nk) - Ts
                Te          = ranw()*Tau_Max + Ts
                Tau_Length  = Te - Ts
                Idx_Ts      = 1
                Idx_Start   = (/nk+1, 1/)
                Idx_End     = (/nk,   1/)
             END IF
          ELSE
             IF (Ts > Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk) .AND. Ts < Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)) THEN
                ! >>> It will be an anti-segment breaking the last segment in the kink list
                Tau_Max    = Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) - Ts
                Te         = ranw()*Tau_Max + Ts
                Tau_Length = Te - Ts
                Idx_Te     = nk
                Idx_Ts     = nk+1
                Idx_Start  = (/nk+1, nk+1/)
                Idx_End    = (/nk+1, nk/)
             END IF
          END IF

          DO i = 1, nk-1
             IF (Ts > Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i) .AND. Ts < Kink(Idx_Orbit, Idx_Spin)%Tau_End(i)) THEN
                Tau_Max    = Kink(Idx_Orbit, Idx_Spin)%Tau_End(i) - Ts
                Te         = ranw()*Tau_Max + Ts
                Tau_Length = Te - Ts
                Idx_Te     = i
                Idx_Ts     = i+1
                Idx_Start  = (/nk+1, Idx_Ts/)
                Idx_End    = (/nk+1, Idx_Te/)
             END IF
          END DO
          
       END IF
       IF (DEBUG) WRITE(*, '(3x, a, 2i5, 2f12.6)')  'Runtime Info: insert antisegment', Idx_Orbit, Idx_Spin, Ts, Te
       
    END IF

    IF (Idx_Ts == 0) RETURN   ! Ts locates on a segment when inserting segment, or on an anit-segment when inserting antisegment
    IF (abs(Tau_Max) < 1.d-4) RETURN
    IF (abs(Tau_Length) < 1.d-3) RETURN

    IF (LSegment) THEN
       Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk+1) = Ts
       Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk+1)   = Te
    ELSE
       Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk+1) = Te
       Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk+1)   = Ts
    END IF
    ! 0 for not occupied. To be determined. 
    Kink(Idx_Orbit, Idx_Spin)%n_Start(nk+1, 1:nOrbit, 1:nSpin) = 0       
    
    ! STEP[2] : Calculate the overlap of the newly inserted kink with all other exisiting segments
    CALL Segment_Overlap(Ts, Te, Idx_Orbit, Idx_Spin, Overlap)

    ALLOCATE(Gi(nk+1))
    ALLOCATE(Gj(nk+1))
    ALLOCATE(XL(nk))
    Gi = Zero
    Gj = Zero
    DO i = 1, nk+1
       Tau = Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk+1) - Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i)
       Gi(i) =  HybridizationFunc(Tau, Idx_Orbit, Idx_Spin)
       
       Tau = Kink(Idx_Orbit, Idx_Spin)%Tau_End(i) - Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk+1)
       Gj(i) = HybridizationFunc(Tau, Idx_Orbit, Idx_Spin)
    END DO
    
    ! >>> calculate Update_Ratio = Update_Ratio + Gi(i)*Kink(Idx_Orbit, Idx_Spin)%Det(i, j)*Gj(j)    
!    CALL DGEMV('N', nk, nk, One, Kink(Idx_Orbit, Idx_Spin)%Det(1:nk, 1:nk), nk, Gj(1:nk), 1, Zero, XL(1:nk), 1)
    XL = Zero
    DO i = 1, nk
       DO j = 1, nk
          XL(i) = XL(i) + Kink(Idx_Orbit, Idx_Spin)%Det(i, j)*Gj(j)
       END DO
    END DO
    Update_Ratio = Zero    
    DO i = 1, nk
       Update_Ratio = Update_Ratio + Gi(i)*XL(i)
    END DO    
    Det_Ratio = Gi(nk+1) - Update_Ratio
!    print*, Gi(nk+1), Update_Ratio
    
    IF (LSegment) THEN
       Temp = Eimp(Idx_Orbit, Idx_Spin)*Tau_Length
       Update_Ratio =  Beta*Tau_Max/DBLE(nk+1)*Det_Ratio 
    ELSE
       Temp = -Eimp(Idx_Orbit, Idx_Spin)*Tau_Length
       Update_Ratio = -Beta*Tau_Max/DBLE(nk+1)*Det_Ratio
    END IF
    IF (Ts > Te) Update_Ratio = -Update_Ratio

    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          IF (LSegment) THEN
             Temp = Temp - Umat(Idx_Orbit, Idx_Spin, iOrbit, iSpin)*Overlap(iOrbit, iSpin)
          ELSE
             Temp = Temp + Umat(Idx_Orbit, Idx_Spin, iOrbit, iSpin)*Overlap(iOrbit, iSpin)
          END IF
       END DO
    END DO 
    Update_Ratio = Update_Ratio*DEXP(Temp)

    ! >>> Record the sign of the current configuration into Statistic%MCsign(1)  
    IF (Update_Ratio < Zero) THEN
       Segment_Statistic%MCsign(1) = -One
     ELSE
       Segment_Statistic%MCsign(1) = One
    END IF
   
    IF (DEBUG) WRITE(*, *) 'Insert segment ', Update_Ratio
    IF (ABS(Update_Ratio) < 1.d-3) RETURN 
 
    IF (ranw() < (Update_Ratio)) THEN

       ! >>> Record Sign(2) 
       IF (Update_Ratio < Zero) THEN
          Segment_Statistic%MCsign(2) = -One
        ELSE
          Segment_Statistic%MCsign(2) = One
       END IF
       
       IF (LSegment) THEN
          Segment_Statistic%nInsertSegment(Idx_Orbit, Idx_Spin)     = Segment_Statistic%nInsertSegment(Idx_Orbit, Idx_Spin)     + 1
       ELSE
          Segment_Statistic%nInsertAntiSegment(Idx_Orbit, Idx_Spin) = Segment_Statistic%nInsertAntiSegment(Idx_Orbit, Idx_Spin) + 1
       END IF

       ! >>> update the determinant matrix 
       ALLOCATE(XR(nk))
 !      CALL DGEMV('T', nk, nk, One, Kink(Idx_Orbit, idx_Spin)%Det(1:nk, 1:nk), nk, Gi(1:nk), 1, Zero, XR(1:nk), 1)
       XR = Zero
       DO j = 1, nk
          DO i = 1, nk
             XR(j) = XR(j) + Gi(i)*Kink(Idx_Orbit, Idx_Spin)%Det(i, j)
          END DO
       END DO
       
       DO i = 1, nk
          DO j = 1, nk
             Kink(Idx_Orbit, Idx_Spin)%Det(i, j) =  Kink(Idx_Orbit, Idx_Spin)%Det(i, j)  + XL(i)*XR(j)/Det_Ratio
          END DO
          Kink(Idx_Orbit, Idx_Spin)%Det(i, nk+1) = -XL(i)/Det_Ratio
          Kink(Idx_Orbit, Idx_Spin)%Det(nk+1, i) = -XR(i)/Det_Ratio
       END DO
       Kink(Idx_Orbit, Idx_Spin)%Det(nk+1, nk+1) = One/Det_Ratio

       ! >>> check if the interted Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk+1) is occupied in other flavor or not.
       CALL Check_Occupancy( LSegment, Ts, Te, Idx_Orbit, Idx_Spin)
       
       
       ! >>> Sorting Tau_Start, Tau_End and update them to Kink%Tau_Start and Kink%Tau_End
       Temp  = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_Start(1))
       DO i = Idx_Start(1)-1, Idx_Start(2), -1
          Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i+1)                    = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i)
       END DO
       Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_Start(2)) = Temp

       Temp = Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_End(1))
       Occupancy_Temp(1:nOrbit, 1:nSpin) = Kink(Idx_Orbit, Idx_Spin)%n_Start(Idx_End(1), 1:nOrbit, 1:nSpin)
       DO i = Idx_End(1)-1, Idx_End(2), -1
          Kink(Idx_Orbit, Idx_Spin)%Tau_End(i+1)                      =  Kink(Idx_Orbit, Idx_Spin)%Tau_End(i)
          Kink(Idx_Orbit, Idx_Spin)%n_Start(i+1, 1:nOrbit, 1:nSpin)   = Kink(Idx_Orbit, Idx_Spin)%n_Start(i, 1:nOrbit, 1:nSpin)         
       END DO
       Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_End(2)) = Temp
       Kink(Idx_Orbit, Idx_Spin)%n_Start(Idx_End(2), 1:nOrbit, 1:nSpin) = Occupancy_Temp(1:nOrbit, 1:nSpin) 

       ! >>> rearrange Det matrix according to the order of the sorted Tau_Start and Tau_End
       DO j = 1, nk+1
          Temp = Kink(Idx_Orbit, Idx_Spin)%Det(Idx_Start(1), j)
          DO i = Idx_Start(1)-1, Idx_Start(2), -1
             Kink(Idx_Orbit, Idx_Spin)%Det(i+1, j) = Kink(Idx_Orbit, Idx_Spin)%Det(i, j)
          END DO
          Kink(Idx_Orbit, Idx_Spin)%Det(Idx_Start(2), j) = Temp
       END DO

       DO i = 1, nk+1
          Temp = Kink(Idx_Orbit, Idx_Spin)%Det(i, Idx_End(1))
          DO j = Idx_End(1)-1, Idx_End(2), -1
             Kink(Idx_Orbit, Idx_Spin)%Det(i, j+1) = Kink(Idx_Orbit, Idx_Spin)%Det(i, j)
          END DO
          Kink(Idx_Orbit, Idx_Spin)%Det(i, Idx_End(2)) = Temp
       END DO

       IF (LSegment) THEN
          Kink(Idx_Orbit, Idx_Spin)%Total_Length  = Kink(Idx_Orbit, Idx_Spin)%Total_Length + Tau_Length
       ELSE
          Kink(Idx_Orbit, Idx_Spin)%Total_Length  = Kink(Idx_Orbit, Idx_Spin)%Total_Length - Tau_Length
       END IF
       Kink(Idx_Orbit, Idx_Spin)%nk = nk + 1
       DEALLOCATE(XR)
       
    ELSE
       Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk+1) = Zero
       Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk+1)   = Zero
    END IF

    IF (DEBUG) THEN
       WRITE(*, *) 'Updated Tau_Start and Tau_End list'
       DO i = 1, Kink(Idx_Orbit, Idx_Spin)%nk
          WRITE(*, '(i5, 2f12.6)') i, Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i), Kink(Idx_Orbit, Idx_Spin)%Tau_End(i)
       END DO
      
       CALL Segment_Debug(Idx_Orbit, Idx_Spin)
       
    END IF

    DEALLOCATE(Gi)
    DEALLOCATE(Gj)
    DEALLOCATE(XL)
    
  END SUBROUTINE Insert_Kink

  !---------------------------------------------------------------------------------------
  SUBROUTINE Remove_Kink(Idx_Orbit, Idx_Spin)

    IMPLICIT NONE
    
    !
    ! Purposes
    ! =======
    !

    INTEGER, INTENT(IN)  :: Idx_Orbit, Idx_Spin
    
    ! ... local vars ...
    INTEGER  :: Idx_Kink
    INTEGER  :: Idx_Start(2), Idx_End(2)
    INTEGER  :: Is, Ie
    REAL(DP) :: Ts, Te                                 ! the starting and ending point of the inserted kink
    REAL(DP) :: Tau_Max, Tau_Length
    REAL(DP) :: Update_Ratio, Det_Ratio, Temp
    REAL(DP) :: Overlap(nOrbit, nSpin)
    
    INTEGER  :: nk
    INTEGER  :: i, j
    INTEGER  :: iOrbit, iSpin
    INTEGER  :: Occupancy_Temp(nOrbit, nSpin)
    LOGICAL  :: LSegment = .True.                      ! .True. for adding segment, .False. for adding antisegment
    LOGICAL  :: EmptySegment = .True.                  ! when nk=0, .True. for empty segment, .False. for full segment
    
    ! ... exectuable ...
    nk = Kink(Idx_Orbit, Idx_Spin)%nk           ! get the current order
    IF (nk == 0) RETURN
    
    ! STEP[1]: randomly select one kink to remove, can be either segment or anti-segment
    Idx_Kink = CEILING( ranw() * nk )

    IF (ranw() > Half) THEN
       
       IF (DEBUG) WRITE(*, '(3x, a, 2i5)')  'Runtime Info: remove segment', Idx_Orbit, Idx_Spin
              
       !  >>> segment is always taken as #Idx_Kink
       LSegment = .True.
       Ts       = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_Kink)
       Te       = Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_Kink)
       IF (Idx_Kink == nk) THEN
          Tau_Max =  Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1) - Ts + Beta
       ELSE
          Tau_Max =  Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_Kink+1) -  Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_kink)
       END IF
       Tau_Length = Te - Ts
       IF (Tau_Length < -epsilon) Tau_Length = Tau_Length + Beta
       
       Idx_Start = (/Idx_Kink, nk/)
       Idx_End   = (/Idx_Kink, nk/)
       IF (nk == 1) EmptySegment = .True.

       Is = Idx_Kink
       Ie = Idx_Kink
    ELSE
       
       IF (DEBUG) WRITE(*, '(3x, a, 2i5)')  'Runtime Info: remove antisegment', Idx_Orbit, Idx_Spin
              
       !  >>> anti-segment is always taken as the one before Tau_Start(Idx_Kink)
       LSegment = .False.
       IF (Idx_Kink == 1) THEN
          Ts            = Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)
          Te            = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1)
          IF ( Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) <  Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk)) THEN
             Tau_Length = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1) - Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)
             Tau_Max    = Kink(Idx_Orbit, Idx_Spin)%Tau_End(1) - Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)
             IF (nk == 1) Tau_Max = Beta
          ELSE
             Tau_Length = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1) + Beta - Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)
             Tau_Max    = Kink(Idx_Orbit, Idx_Spin)%Tau_End(1)   + Beta - Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)
          END IF

          IF (nk > 1) THEN
             Idx_Start = (/1, nk/)
             Idx_End   = (/1, nk-1/)
          ELSE
             Idx_Start = 1
             Idx_End   = 1
          END IF
          Is = 1
          Ie = nk
       ELSE
          Ts         = Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_Kink-1)
          Te         = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_Kink)
          Tau_Length = Te - Ts
          Tau_Max    = Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_Kink) - Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_kink-1)
          IF (Tau_Max < -epsilon) Tau_Max = Tau_Max + Beta  
          
          Idx_Start = (/Idx_Kink, nk/)
          Idx_End   = (/Idx_Kink-1, nk/)
          Is = Idx_Kink
          Ie = Idx_Kink-1
       END IF
       IF (nk == 1) EmptySegment = .False.   ! This will be a full segment
       
    END IF

    IF (Te > Ts) THEN
       Det_Ratio    =  Kink(Idx_Orbit, Idx_Spin)%Det(Is, Ie)
    ELSE
       Det_Ratio    = -Kink(Idx_Orbit, Idx_Spin)%Det(Is, Ie)
    END IF
    
    IF (LSegment) THEN
       Temp = -Eimp(Idx_Orbit, Idx_Spin)*Tau_Length
       Update_Ratio =  Det_Ratio*nk/Beta/Tau_Max
    ELSE
       Temp = Eimp(Idx_Orbit, Idx_Spin)*Tau_Length
       Update_Ratio = -Det_Ratio*nk/Beta/Tau_Max
    END IF

    ! STEP[2] : calculate the overlap change after removing the proposed kink
    CALL Segment_Overlap(Ts, Te, Idx_Orbit, Idx_Spin, Overlap)

    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          IF (LSegment) THEN
             ! >>> remove a segment
             Temp = Temp + Umat(Idx_Orbit, Idx_Spin, iOrbit, iSpin)*Overlap(iOrbit, iSpin)
          ELSE
             ! >>> remove an antisegment
             Temp = Temp - Umat(Idx_Orbit, Idx_Spin, iOrbit, iSpin)*Overlap(iOrbit, iSpin)
          END IF
       END DO
    END DO
    Update_Ratio = Update_Ratio*DEXP(Temp)

    IF (DEBUG) WRITE(*, *) 'Update Ratio is: ', Update_Ratio
    IF (ABS(Update_Ratio) < 1.d-3) RETURN        
 
    ! >>> Record the sign of the current configuration into Statistic%MCsign(1)  
    IF (Update_Ratio < Zero) THEN
       Segment_Statistic%MCsign(1) = -One
     ELSE
       Segment_Statistic%MCsign(1) = One
    END IF
    
    IF (DEBUG) WRITE(*, *) 'Remove Segment ', Update_Ratio

    IF (ranw() < (Update_Ratio)) THEN
       
       ! >>> Record Sign(2) 
       IF (Update_Ratio < Zero) THEN
          Segment_Statistic%MCsign(2) = -One
       ELSE
          Segment_Statistic%MCsign(2) = One
       END IF
  
       ! >>> record the number of elementary operations
       IF (LSegment) THEN
          Segment_Statistic%nRemoveSegment(Idx_Orbit, Idx_Spin)     = Segment_Statistic%nRemoveSegment(Idx_Orbit, Idx_Spin) + 1
       ELSE
          Segment_Statistic%nRemoveAntiSegment(Idx_Orbit, Idx_Spin) = Segment_Statistic%nRemoveAntiSegment(Idx_Orbit, Idx_Spin) + 1
       END IF

       !  >>> Check the occupancy change of Tau_Start in other flavors after the remove operation
       DO iOrbit = 1, nOrbit
          DO iSpin = 1, nSpin
             IF (iOrbit /= Idx_Orbit .OR. iSpin /= Idx_Spin) THEN
                IF (LSegment) THEN
                   ! >>> remove a segment
                   DO i = 1, Kink(iOrbit, iSpin)%nk
                      IF (Te > Ts) THEN
                         IF (Kink(iOrbit, iSpin)%Tau_End(i) > Ts .AND. Kink(iOrbit, iSpin)%Tau_End(i) < Te) Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 0
                      ELSE
                         IF (Kink(iOrbit, iSpin)%Tau_End(i) < Te .OR. Kink(iOrbit, iSpin)%Tau_End(i) > Ts)  Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 0
                      END IF
                   END DO
                ELSE
                   ! >>> remove an antisegment
                   DO i = 1, Kink(iOrbit, iSpin)%nk
                      IF (Te > Ts) THEN
                         IF (Kink(iOrbit, iSpin)%Tau_End(i) > Ts .AND. Kink(iOrbit, iSpin)%Tau_End(i) < Te) Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 1
                      ELSE
                         IF (Kink(iOrbit, iSpin)%Tau_End(i) < Te .OR. Kink(iOrbit, iSpin)%Tau_End(i) > Ts)  Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 1
                      END IF
                   END DO
                END IF
             END IF
          END DO
       END DO
              
       ! rearrange the determinental matrix
       DO j = 1, nk
          Temp = Kink(Idx_Orbit, Idx_Spin)%Det(Idx_Start(1), j) 
          DO i = Idx_Start(1), Idx_Start(2)-1
             Kink(Idx_Orbit, Idx_Spin)%Det(i, j) = Kink(Idx_Orbit, Idx_Spin)%Det(i+1, j)
          END DO
          Kink(Idx_Orbit, Idx_Spin)%Det(Idx_Start(2), j) = Temp
       END DO
       DO i = 1, nk
          Temp = Kink(Idx_Orbit, Idx_Spin)%Det(i, Idx_End(1))
          DO j = Idx_End(1), Idx_End(2)-1
             Kink(Idx_Orbit, Idx_Spin)%Det(i, j) = Kink(Idx_Orbit, Idx_Spin)%Det(i, j+1)
          END DO
          Kink(Idx_Orbit, Idx_Spin)%Det(i, Idx_End(2)) = Temp
       END DO
       
       ! update the determinal matrix
       DO i = 1, nk-1
          DO j = 1, nk-1
             Kink(Idx_Orbit, Idx_Spin)%Det(i, j) = Kink(Idx_Orbit, Idx_Spin)%Det(i, j) &
                  - Kink(Idx_Orbit, Idx_Spin)%Det(i, nk)*Kink(Idx_Orbit, Idx_Spin)%Det(nk, j)/Kink(Idx_Orbit, Idx_Spin)%Det(nk, nk)
          END DO
       END DO


       Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk+1) = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_Start(1))
       DO i = Idx_Start(1), Idx_Start(2)-1
          Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i) = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i+1)
       END DO
       Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_Start(2)) = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk+1)
       Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk+1)         = Zero

       Occupancy_Temp(1:nOrbit, 1:nSpin)       = Kink(Idx_Orbit, Idx_Spin)%n_Start(Idx_End(1), 1:nOrbit, 1:nSpin)
       Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk+1) = Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_End(1))
       DO i = Idx_End(1), Idx_End(2)-1
          Kink(Idx_Orbit, Idx_Spin)%Tau_End(i) = Kink(Idx_Orbit, Idx_Spin)%Tau_End(i+1)
          Kink(Idx_Orbit, Idx_Spin)%n_Start(i, 1:nOrbit, 1:nSpin)   = Kink(Idx_Orbit, Idx_Spin)%n_Start(i+1, 1:nOrbit, 1:nSpin) 
       END DO
       Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_End(2)) = Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk+1)
       Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk+1)       = Zero
       Kink(Idx_Orbit, Idx_Spin)%n_Start(Idx_End(2), 1:nOrbit, 1:nSpin) = Occupancy_Temp(1:nOrbit, 1:nSpin)
       
       IF (nk > 1) THEN
          Kink(Idx_Orbit, Idx_Spin)%nk = nk - 1
          Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk)     = Zero
          Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)       = Zero
       ELSE
          Kink(Idx_Orbit, Idx_Spin)%nk = 0
          Kink(Idx_Orbit, Idx_Spin)%Tau_Start = Zero
          Kink(Idx_Orbit, Idx_Spin)%Tau_End   = Zero
          Kink(Idx_Orbit, Idx_Spin)%Det       = Zero
       END IF

       IF (LSegment) THEN
          Kink(Idx_Orbit, Idx_Spin)%Total_Length = Kink(Idx_Orbit, Idx_Spin)%Total_Length - Tau_Length
       ELSE
          Kink(Idx_Orbit, Idx_Spin)%Total_Length = Kink(Idx_Orbit, Idx_Spin)%Total_Length + Tau_Length          
       END IF
       IF (Kink(Idx_Orbit, Idx_Spin)%nk == 0) THEN
          Kink(Idx_Orbit, Idx_Spin)%EmptySegment = EmptySegment
          IF (EmptySegment) THEN
             Kink(Idx_Orbit, Idx_Spin)%Total_Length = Zero
          ELSE
             Kink(Idx_Orbit, Idx_Spin)%Total_Length = Beta
          END IF
       END IF
       
    END IF

    IF (DEBUG) THEN
       WRITE(*, *) 'Updated Tau_Start and Tau_End list'
       DO i = 1, Kink(Idx_Orbit, Idx_Spin)%nk
          WRITE(*, '(i5, 2f12.6)') i, Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i), Kink(Idx_Orbit, Idx_Spin)%Tau_End(i)
       END DO
       
       CALL Segment_Debug(Idx_Orbit, Idx_Spin)
    END IF
    
  END SUBROUTINE Remove_Kink

  !--------------------------------------------------------------------------------------------------------
  SUBROUTINE Shift_Kink(Idx_Orbit, Idx_Spin)

    IMPLICIT NONE

    !
    ! Purpose
    ! =======
    !   shift the end-point of a randomly selected segment
    !

    INTEGER, INTENT(IN)  :: Idx_Orbit, Idx_Spin
    
    ! ... local vars ...
    INTEGER  :: Idx_Kink
    REAL(DP) :: Ts, Te                                 ! the starting and ending point of the inserted kink
    REAL(DP) :: Tau_Max, Tau_Length
    REAL(DP) :: Update_Ratio, Tau
    REAL(DP) :: Overlap(nOrbit, nSpin), Temp
    REAL(DP), ALLOCATABLE :: AB(:), Binv(:, :), Gi(:), Det(:, :), Tau_Start(:), Tau_End(:)  

    INTEGER  :: i, j, nk
    INTEGER  :: iOrbit, iSpin
    
    
    nk = Kink(Idx_Orbit, Idx_Spin)%nk           ! get the current order
    ALLOCATE(Tau_Start(nk))
    ALLOCATE(Tau_End(nk))
    Tau_Start(1:nk) = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1:nk)
    Tau_End(1:nk)   = Kink(Idx_Orbit, Idx_Spin)%Tau_End(1:nk)
    
    IF (DEBUG) WRITE(*, '(3x, a, 2i5)')  'Runtime Info: try to shift', Idx_Orbit, Idx_Spin
           
    IF (nk <= 0) RETURN

    ! Randomly pich a segment to shift
    Idx_Kink   = CEILING( ranw()*nk )
    
    IF (Idx_Kink /= nk) THEN
       Tau_Max    = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_Kink+1) -  Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_Kink)
       Ts         = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(Idx_Kink)
       Te         = ranw()*Tau_Max + Ts
       Ts         = Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_Kink)
       Tau_Length = Te - Ts   ! This is the "change" of segment length, when >0 it is a segment; when <0 it is an antisegment 
    ELSE
       Tau_Max    = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1) - Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk) + Beta
       Ts         = Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk)
       Te         = ranw()*Tau_Max + Ts
       Ts         = Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)

       IF (Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) < Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1)) THEN
          Tau_Length = Te - Ts - Beta
       ELSE
          Tau_Length = Te - Ts  ! This is the "change" of segment length
       END IF
       IF (Te > Beta) Te = Te - Beta
    END IF
    
    Tau_End(Idx_Kink) = Te
   
    IF (ABS(Tau_Max) < 1.d-4) RETURN
    IF (ABS(Tau_Length) < 1.d-3) RETURN
 
    IF (Tau_Length > Zero) THEN
       CALL Segment_Overlap(Ts, Te, Idx_Orbit, Idx_Spin, Overlap)
    ELSE
       CALL Segment_Overlap(Te, Ts, Idx_Orbit, Idx_Spin, Overlap)
    END IF
    
    ! Calculate the update ratio
    ALLOCATE(Gi(nk))
    ALLOCATE(AB(nk))
    
    DO i = 1, nk
       Tau = Te - Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i)
       Gi(i) = HybridizationFunc(Tau, Idx_Orbit, Idx_Spin)
    END DO

    AB   = Zero
    DO i = 1, nk
       DO j = 1, nk
          AB(i) = AB(i) + Gi(j)*Kink(Idx_Orbit, Idx_Spin)%Det(j, i)
       END DO
    END DO
    Update_Ratio = AB(Idx_Kink)

    ! the correct sign of the exponent is ensured by the sign of Tau_Length. Conceptually, it can be always viewd as a segment 
    Temp = Eimp(Idx_Orbit, Idx_Spin)*Tau_Length

    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          IF (Tau_Length > Zero) THEN
             Temp = Temp - Umat(Idx_Orbit, Idx_Spin, iOrbit, iSpin)*Overlap(iOrbit, iSpin)
          ELSE
             Temp = Temp + Umat(Idx_Orbit, Idx_Spin, iOrbit, iSpin)*Overlap(iOrbit, iSpin)
          END IF
       END DO
    END DO
    Update_Ratio = Update_Ratio*DEXP(Temp)  
 
    IF ( (Tau_End(nk)-Tau_Start(nk))*(Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) - Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk)) < 0 ) THEN
       Update_Ratio = -Update_Ratio
    END IF

    ! >>> Record the sign of the current configuration into Statistic%MCsign(1)  
    IF (Update_Ratio < Zero) THEN
       Segment_Statistic%MCsign(1) = -One
    ELSE
       Segment_Statistic%MCsign(1) = One
    END IF
    
    IF (DEBUG) WRITE(*, *) "Shift Segment ", Update_Ratio
    IF (ABS(Update_Ratio) < 1.d-3) RETURN     

    IF (ranw() < (Update_Ratio)) THEN
       
       ! >>> Record Sign(2) 
       IF (Update_Ratio < Zero) THEN
          Segment_Statistic%MCsign(2) = -One
        ELSE
          Segment_Statistic%MCsign(2) = One
       END IF

       ! >>> record the number of elementary operations
       Segment_Statistic%nShiftSegment(Idx_Orbit, Idx_Spin) = Segment_Statistic%nShiftSegment(Idx_Orbit, Idx_Spin) + 1

       ! >>> check the change of Occupancy after the shift
       DO iOrbit = 1, nOrbit
          DO iSpin = 1, nSpin

             IF (iOrbit /= Idx_Orbit .OR. iSpin /= Idx_Spin) THEN

                !  >>> occupancy change of the shifted Kink
                Kink(Idx_Orbit, Idx_Spin)%n_Start(Idx_Kink, iOrbit, iSpin) = 0
                IF ( Kink(iOrbit, iSpin)%nk == 0 ) THEN
                   IF ( .NOT. Kink(iOrbit, iSpin)%EmptySegment ) THEN
                      Kink(Idx_Orbit, Idx_Spin)%n_Start(Idx_Kink, iOrbit, iSpin) = 1
                   ELSE
                      Kink(Idx_Orbit, Idx_Spin)%n_Start(Idx_Kink, iOrbit, iSpin) = 0
                   END IF
                ELSE
                   IF ( Kink(iOrbit, iSpin)%Tau_End(Kink(iOrbit, iSpin)%nk) < Kink(iOrbit, iSpin)%Tau_Start(Kink(iOrbit, iSpin)%nk) ) THEN
                      IF (Te < Kink(iOrbit, iSpin)%Tau_End(Kink(iOrbit, iSpin)%nk) .OR. Te > Kink(iOrbit, iSpin)%Tau_Start(Kink(iOrbit, iSpin)%nk)) THEN
                         Kink(Idx_Orbit, Idx_Spin)%n_Start(Idx_Kink, iOrbit, iSpin) = 1
                      ELSE
                         DO i = 1, Kink(iOrbit, iSpin)%nk-1
                            IF (Te > Kink(iOrbit, iSpin)%Tau_Start(i) .AND. Te < Kink(iOrbit, iSpin)%Tau_End(i)) THEN
                               Kink(Idx_Orbit, Idx_Spin)%n_Start(Idx_Kink, iOrbit, iSpin) = 1
                            END IF
                         END DO
                      END IF                      
                   ELSE
                      DO i = 1, Kink(iOrbit, iSpin)%nk
                         IF (Te < Kink(iOrbit, iSpin)%Tau_End(i) .AND. Te > Kink(iOrbit, iSpin)%Tau_Start(i)) THEN
                            Kink(Idx_Orbit, Idx_Spin)%n_Start(Idx_Kink, iOrbit, iSpin) = 1
                         END IF
                      END DO
                   END IF
                END IF
                
                ! >>> changes in other flavors
                IF ( Idx_Kink /= nk ) THEN
                   IF ( Te > Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_Kink) ) THEN
                      ! segment is enlongated. 
                      DO i = 1, Kink(iOrbit, iSpin)%nk
                         IF (Kink(iOrbit, iSpin)%Tau_End(i) < Te .AND. Kink(iOrbit, iSpin)%Tau_End(i) > Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_Kink)) THEN
                            Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 1 
                         END IF
                      END DO
                   ELSE
                      ! segment is shorten
                      DO i = 1, Kink(iOrbit, iSpin)%nk
                         IF (Kink(iOrbit, iSpin)%Tau_End(i) < Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_Kink) .AND. Kink(iOrbit, iSpin)%Tau_End(i) > Te) THEN
                            Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 0
                         END IF
                      END DO
                   END IF
                ELSE
                   IF (( Te - Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk))*(Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) - Kink(Idx_Orbit, Idx_Spin)%Tau_Start(nk)) > Zero ) THEN
                      IF (Te > Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)) THEN
                         DO i = 1, Kink(iOrbit, iSpin)%nk
                            IF ( Kink(iOrbit, iSpin)%Tau_End(i) < Te .AND. Kink(iOrbit, iSpin)%Tau_End(i) > Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) ) THEN
                               Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 1 
                            END IF
                         END DO
                      ELSE
                         DO i = 1, Kink(iOrbit, iSpin)%nk
                            IF ( Kink(iOrbit, iSpin)%Tau_End(i) > Te .AND. Kink(iOrbit, iSpin)%Tau_End(i) < Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) ) THEN
                               Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 0
                            END IF
                         END DO
                      END IF
                   ELSE
                      ! New Te leads to the wraping abound beta.
                      IF (Te < Kink(Idx_Orbit, Idx_Spin)%Tau_Start(1)) THEN
                         ! Te is the smallest kink point, Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) is the largest kink point
                         DO i = 1, Kink(iOrbit, iSpin)%nk
                            IF (Kink(iOrbit, iSpin)%Tau_End(i) < Te .OR. Kink(iOrbit, iSpin)%Tau_End(i) > Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)) THEN
                               Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 1
                            END IF
                         END DO
                      ELSE
                         ! Te is the largest Kink point, Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk) is the smallest kink point.
                         DO i = 1, Kink(iOrbit, iSpin)%nk
                            IF (Kink(iOrbit, iSpin)%Tau_End(i) > Te .OR. Kink(iOrbit, iSpin)%Tau_End(i) < Kink(Idx_Orbit, Idx_Spin)%Tau_End(nk)) THEN
                               Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 0
                            END IF
                         END DO 
                      END IF
                   END IF    
                END IF
                
             END IF
          END DO
       END DO

       ! >>> update the determinant matrix
       ALLOCATE(Binv(nk, nk))
       ALLOCATE(Det(nk, nk))
       DO i = 1, nk
          DO j = 1, nk
             Binv(i, j) = Zero
             IF (i == Idx_Kink) THEN
                Binv(i, j) = -AB(j)/AB(Idx_Kink)
             END IF
          END DO
          BinV(i, i) = One
       END DO
       Binv(Idx_Kink, Idx_Kink) = One/AB(Idx_Kink)
       Temp = Kink(Idx_Orbit, Idx_Spin)%Det(Idx_Kink, Idx_Kink)
       
       Det      = Zero
       CALL DGEMM('N', 'N', nk, nk, nk, One, Kink(Idx_Orbit, Idx_Spin)%Det(1:nk, 1:nk), nk, Binv(1:nk, 1:nk), nk, Zero, Det(1:nk, 1:nk), nk)
       Kink(Idx_Orbit, Idx_Spin)%Det(1:nk, 1:nk)         =  Det(1:nk, 1:nk)
       Kink(Idx_Orbit, Idx_Spin)%Det(Idx_Kink, Idx_Kink) =  Temp/AB(Idx_Kink) 
       Kink(Idx_Orbit, Idx_Spin)%Tau_End(Idx_Kink)       =  Tau_End(Idx_Kink)
       
       Kink(Idx_Orbit, Idx_Spin)%Total_Length = Kink(Idx_Orbit, Idx_Spin)%Total_Length + Tau_Length
       
       DEALLOCATE(Binv)
       DEALLOCATE(Det)
    END IF

    DEALLOCATE(Gi)
    DEALLOCATE(AB)
    DEALLOCATE(Tau_Start)
    DEALLOCATE(Tau_End)
    
        
  END SUBROUTINE SHIFT_KINK
  
  !-------------------------------------------------------------------------------------------------------
  SUBROUTINE Empty_Full(Idx_Orbit, Idx_Spin)

    IMPLICIT NONE
    !

    INTEGER, INTENT(IN)  :: Idx_Orbit
    INTEGER, INTENT(IN)  :: Idx_Spin
    
    ! ... local vars ...
    INTEGER  :: nk, iOrbit, iSpin, i
    REAL(DP) :: Overlap(nOrbit, nSpin), Tau_Length, Update_Ratio
    
    nk = Kink(Idx_Orbit, Idx_Spin)%nk           ! get the current order

    IF (nk /= 0) RETURN

    Overlap = Zero
    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          IF (iOrbit /= Idx_Orbit .OR. iSpin /= Idx_Spin) THEN
             Tau_Length = Zero
             DO i = 1, Kink(iOrbit, iSpin)%nk
                Tau_Length = Kink(iOrbit, iSpin)%Tau_End(i) - Kink(iOrbit, iSpin)%Tau_Start(i)
                IF (Tau_Length < Zero) Tau_Length = Tau_Length + Beta
                Overlap(iOrbit, iSpin) = Overlap(iOrbit, iSpin) + Tau_Length
             END DO
             IF (Kink(iOrbit, iSpin)%nk == 0 .AND. .NOT. Kink(iOrbit, iSpin)%EmptySegment) Overlap(iOrbit, iSpin) = Beta
          END IF
       END DO
    END DO

    Update_Ratio = Zero
    IF (ranw() >= Half) THEN
       IF ( .NOT. Kink(Idx_Orbit, Idx_Spin)%EmptySegment )  THEN
          Segment_Statistic%MCsign(1) = One ! Sign is always possitive
          DO iOrbit = 1, nOrbit
             DO iSpin = 1, nSpin
                IF (iOrbit /= Idx_Orbit .OR. iSpin /= Idx_Spin) THEN
                   Update_Ratio = DEXP( Umat(Idx_Orbit, Idx_Spin, iOrbit, iSpin)*Overlap(iOrbit, iSpin) - Beta*Eimp(Idx_Orbit, Idx_Spin))
                END IF
             END DO
          END DO
          IF (ranw() < Update_Ratio) THEN
             Segment_Statistic%MCsign(2) = One   ! sign is always positive
             Kink(Idx_Orbit, Idx_Spin)%EmptySegment = .False.

             ! >>> record the number of elementary operations
             Segment_Statistic%nEmptyFull(Idx_Orbit, Idx_Spin) = Segment_Statistic%nEmptyFull(Idx_Orbit, Idx_Spin) + 1

             DO iOrbit = 1, nOrbit
                DO iSpin = 1, nSpin
                   IF (iOrbit /= Idx_Orbit .OR. iSpin /= Idx_Spin) Kink(iOrbit, iSpin)%n_Start(1:Kink(iOrbit, iSpin)%nk, Idx_Orbit, Idx_Spin) = 0
                END DO
             END DO
          END IF
       END IF
    ELSE    
       IF (Kink(Idx_Orbit, Idx_Spin)%EmptySegment) THEN
          Segment_Statistic%MCsign(1) = One   ! Sign is always positive

          DO iOrbit = 1, nOrbit
             DO iSpin = 1, nSpin
                IF (iOrbit /= Idx_Orbit .OR. iSpin /= Idx_Spin) THEN
                   Update_Ratio = Update_Ratio * DEXP( -Umat(Idx_Orbit, Idx_Spin, iOrbit, iSpin)*Overlap(iOrbit, iSpin) + Beta*Eimp(Idx_Orbit, Idx_Spin))
                END IF
             END DO
          END DO
          IF (ranw() < Update_Ratio) THEN
             Segment_Statistic%MCsign(2) = One   ! sign is always positive
             Kink(Idx_Orbit, Idx_Spin)%EmptySegment = .True.

             ! >>> record the number of elementary operations
             Segment_Statistic%nEmptyFull(Idx_Orbit, Idx_Spin) = Segment_Statistic%nEmptyFull(Idx_Orbit, Idx_Spin) + 1

             DO iOrbit = 1, nOrbit
                DO iSpin = 1, nSpin
                   IF (iOrbit /= Idx_Orbit .OR. iSpin /= Idx_Spin) Kink(iOrbit, iSpin)%n_Start(1:Kink(iOrbit, iSpin)%nk, Idx_Orbit, Idx_Spin) = 1
                END DO
             END DO

          END IF
       END IF
    END IF

  END SUBROUTINE Empty_Full
  
  !-------------------------------------------------------------------------------------------------------
  SUBROUTINE Segment_Overlap(Ts, Te, Idx_Orbit, Idx_Spin, Overlap)

    IMPLICIT NONE

    !
    ! Purposes
    ! =======
    !   Determine the overlap of the to-be inserted/removed kink with other segments. Here we always use the
    ! segment notation to calculate the overlap. The overlap associated with antisegment is simply calculated
    ! by subtracting the segment overlap from the total length
    !

    
    INTEGER,  INTENT(IN)  :: Idx_Orbit
    INTEGER,  INTENT(IN)  :: Idx_Spin
    REAL(DP), INTENT(IN)  :: Ts
    REAL(DP), INTENT(IN)  :: Te
    
    REAL(DP), INTENT(OUT) :: Overlap(nOrbit, nSpin)

    ! ... local vars ...
    INTEGER  :: iOrbit, iSpin
    INTEGER  :: i, nk
    REAL(DP) :: T1, T2             ! temperarily store Ts and Te in incremetal order
    REAL(DP) :: Overlap_check, Overlap_Sum    
 
    ! ... executable ...

    T1 = MINVAL( (/Ts, Te/) )
    T2 = MAXVAL( (/Ts, Te/) )
    Overlap      = Zero

    DO iSpin = 1, nSpin
       DO iOrbit = 1, nOrbit
          
          Overlap_Sum  = Zero
          IF (iSpin /= Idx_Spin .OR. iOrbit /= Idx_Orbit) THEN
             
             nk = Kink(iOrbit, iSpin)%nk
             
             IF (nk == 0 .AND. .NOT. Kink(iOrbit, iSpin)%EmptySegment)  THEN
                IF (Te < Ts) THEN
                   Overlap_Sum = Overlap_Sum + Beta + Te - Ts
                ELSE
                   Overlap_Sum = Overlap_Sum + Te - Ts
                END IF
             ELSE
                DO i = 1, nk
                   CALL Check_Overlap(Ts, Te, Kink(iOrbit, iSpin)%Tau_Start(i), Kink(iOrbit, iSpin)%Tau_End(i), Overlap_Check) 
                   Overlap_Sum = Overlap_Sum + Overlap_Check
                END DO
             END IF
             Overlap(iOrbit, iSpin) = Overlap_Sum
             
          END IF
       END DO
    END DO

    IF (DEBUG) THEN
       WRITE(*, *) '---- End of Debug Block ----'  
       WRITE(*, *)
    END IF
    
  END SUBROUTINE Segment_Overlap

  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Check_Overlap(Tau_S1, Tau_E1, Tau_S2, Tau_E2, Overlap)

    IMPLICIT NONE

    !

    REAL(DP), INTENT(IN)  :: Tau_S1, Tau_E1, Tau_S2, Tau_E2
    REAL(DP), INTENT(OUT) :: Overlap
    
    ! ... local vars ...
    REAL(DP) :: TauE1, TauE2, TauE3, TauE4, TauS1, TauS2, TauS3, TauS4, TauE_max, TauS_min
    REAL(DP) :: t
    
    Overlap = Zero
    IF (Tau_E1 < Tau_S1) THEN
       TauS1 = Zero  
       TauE1 = Tau_E1
       TauS2 = Tau_S1 
       TauE2 = Beta
       IF (Tau_E2 < Tau_S2) THEN
          TauS3 = Zero
          TauE3 = Tau_E2
          TauS4 = Tau_S2
          TauE4 = Beta

          TauS_min = MINVAL( (/TauS1, TauS3/) )
          TauE_max = MAXVAL( (/TauE1, TauE3/) )
          Overlap  = (TauE3 - TauS3) + (TauE1 - TauS1) - (TauE_max - TauS_min) 
          IF (Overlap < Zero) Overlap = Zero

          TauS_min = MINVAL( (/TauS1, TauS4/) )
          TauE_max = MAXVAL( (/TauE1, TauE4/) )
          t = (TauE4 - TauS4) + (TauE1 - TauS1) - (TauE_max - TauS_min)
          IF (t < Zero) t = Zero
          Overlap = Overlap + t

          TauS_min = MINVAL( (/TauS2, TauS3/) )
          TauE_max = MAXVAL( (/TauE2, TauE3/) )
          t = (TauE3 - TauS3) + (TauE2 - TauS2) - (TauE_max - TauS_min)
          IF (t < Zero) t = Zero
          Overlap = Overlap + t

          TauS_min = MINVAL( (/TauS2, TauS4/) )
          TauE_max = MAXVAL( (/TauE2, TauE4/) )
          t = (TauE4 - TauS4) + (TauE2 - TauS2) - (TauE_max - TauS_min)
          IF (t < Zero) t = Zero
          Overlap = Overlap + t

       ELSE
          TauS_min = MINVAL( (/TauS1, Tau_S2/) )
          TauE_max = MAXVAL( (/TauE1, Tau_E2/) )
          t = (Tau_E2 - Tau_S2) + (TauE1 - TauS1) - (TauE_max - TauS_min)
          IF (t < Zero) t = Zero
          Overlap = Overlap + t

          TauS_min = MINVAL( (/TauS2, Tau_S2/) )
          TauE_max = MAXVAL( (/TauE2, Tau_E2/) )
          t = (Tau_E2 - Tau_S2) + (TauE2 - TauS2) - (TauE_max - TauS_min)
          IF (t < Zero) t = Zero
          Overlap = Overlap + t

       END IF  
    ELSE
       IF (Tau_E2 < Tau_S2) THEN
          TauS3 = Zero 
          TauE3 = Tau_E2
          TauS4 = Tau_S2
          TauE4 = Beta

          TauS_min = MINVAL( (/Tau_S1, TauS3/) )
          TauE_max = MAXVAL( (/Tau_E1, TauE3/) )
          t = (TauE3 - TauS3) + (Tau_E1 - Tau_S1) - (TauE_max - TauS_min)
          IF (t < Zero) t = Zero
          Overlap = Overlap + t

          TauS_min = MINVAL( (/Tau_S1, TauS4/) )
          TauE_max = MAXVAL( (/Tau_E1, TauE4/) )
          t = (TauE4 - TauS4) + (Tau_E1 - Tau_S1) - (TauE_max - TauS_min)
          IF (t < Zero) t = Zero
          Overlap = Overlap + t

       ELSE
          TauS_min = MINVAL( (/Tau_S1, Tau_S2/) )
          TauE_max = MAXVAL( (/Tau_E1, Tau_E2/) )
          Overlap  = (Tau_E2 - Tau_S2) + (Tau_E1 - Tau_S1) - (TauE_max - TauS_min)         
          IF (Overlap < Zero) Overlap = Zero
       END IF
    END IF

  END SUBROUTINE Check_Overlap
  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Check_Occupancy(LSegment_, Ts, Te, Idx_Orbit, Idx_Spin)

    IMPLICIT NONE

    !
    ! Purpose
    ! =======
    !   Supplementary routine. It determines a given Tau is occupied or not in a given flavor.
    !
    LOGICAL,  INTENT(IN)  :: LSegment_
    REAL(DP), INTENT(IN)  :: Ts, Te
    INTEGER,  INTENT(IN)  :: Idx_Orbit
    INTEGER,  INTENT(IN)  :: Idx_Spin

    ! ... local vars ...
    INTEGER  :: nk
    INTEGER  :: i
    INTEGER  :: iOrbit, iSpin
    REAL(DP) :: Tau
    
    ! ... exectuable ...
    IF (LSegment_) THEN
       Tau = Te
    ELSE
       Tau = Ts
    END IF

    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          
          Kink(Idx_Orbit, Idx_Spin)%n_Start(Kink(Idx_Orbit, Idx_Spin)%nk+1, iOrbit, iSpin) = 0
          
          IF (iOrbit /= Idx_Orbit .OR. iSpin /= Idx_Spin) THEN
             
             ! >>> first, determine the occupancy of the inserted Tau_Start in other flavors             
             nk = Kink(iOrbit, iSpin)%nk
             IF (nk == 0) THEN
                IF ( .NOT. Kink(iOrbit, iSpin)%EmptySegment ) THEN
                   Kink(Idx_Orbit, Idx_Spin)%n_Start(Kink(Idx_Orbit, Idx_Spin)%nk+1, iOrbit, iSpin) = 1
                END IF
             ELSE
                IF ( Kink(iOrbit, iSpin)%Tau_End(nk) < Kink(iOrbit, iSpin)%Tau_Start(nk) ) THEN
                   IF ( Tau < Kink(iOrbit, iSpin)%Tau_End(nk) .OR. Tau > Kink(iOrbit, iSpin)%Tau_Start(nk) ) THEN
                      Kink(Idx_Orbit, Idx_Spin)%n_Start(Kink(Idx_Orbit, Idx_Spin)%nk+1, iOrbit, iSpin) = 1
                   ELSE
                      DO i = 1, nk-1
                         IF ( Tau < Kink(iOrbit, iSpin)%Tau_End(i) .AND. Tau > Kink(iOrbit, iSpin)%Tau_Start(i) ) THEN
                            Kink(Idx_Orbit, Idx_Spin)%n_Start(Kink(Idx_Orbit, Idx_Spin)%nk+1, iOrbit, iSpin) = 1
                            EXIT
                         END IF
                      END DO
                   END IF
                ELSE
                   DO i = 1, nk
                      IF ( Tau < Kink(iOrbit, iSpin)%Tau_End(i) .AND. Tau > Kink(iOrbit, iSpin)%Tau_Start(i) ) THEN
                         Kink(Idx_Orbit, Idx_Spin)%n_Start(Kink(Idx_Orbit, Idx_Spin)%nk+1, iOrbit, iSpin) = 1
                         EXIT
                      END IF
                   END DO
                END IF
             END IF
             
             ! >>> second, update the occupancy change in other flavors
             DO i = 1, nk
                IF (Te > Ts) THEN
                   IF (Kink(iOrbit, iSpin)%Tau_End(i) < Te .AND. Kink(iOrbit, iSpin)%Tau_End(i) > Ts) THEN
                      IF (LSegment_) THEN
                         Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 1
                      ELSE
                         Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 0
                      END IF
                   END IF
                ELSE
                   IF (Kink(iOrbit, iSpin)%Tau_End(i) > Ts .OR. Kink(iOrbit, iSpin)%Tau_End(i) < Te) THEN
                      IF (LSegment_) THEN
                         Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 1
                      ELSE
                         Kink(iOrbit, iSpin)%n_Start(i, Idx_Orbit, Idx_Spin) = 0
                      END IF
                   END IF
                END IF
             END DO
             
          END IF
       END DO
    END DO
    
    
  END SUBROUTINE Check_Occupancy
  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Segment_Debug(Idx_Orbit, Idx_Spin)
    
    IMPLICIT NONE

    !
    ! Purpose
    ! =======
    !    supplementary routine used to debug the update of determinal matrix. We calcualte the hybridization function matrix
    ! directly from list of segment. The determianl matrix is obtained by inversing the hybridization matrix
    !

    INTEGER, INTENT(IN) :: Idx_Orbit
    INTEGER, INTENT(IN) :: Idx_Spin
    
    ! ... local vars ...
    INTEGER  :: i, j
    INTEGER  :: iOrbit, iSpin
    INTEGER  :: nk
    REAL(DP) :: Tau
    REAL(DP), ALLOCATABLE :: XG(:, :), XM(:, :)
    
    nk = Kink(Idx_Orbit, Idx_Spin)%nk

    ALLOCATE(XG(nk, nk))
    ALLOCATE(XM(nk, nk))

    WRITE(*, *) '    ===== Kink list ===='
    DO i = 1, nk
       WRITE(*, '(3x, 2f20.12, 100i3)') Kink(Idx_Orbit, Idx_Spin)%Tau_Start(i), Kink(Idx_Orbit, Idx_Spin)%Tau_End(i), Kink(Idx_Orbit, Idx_Spin)%n_Start(i, 1:nOrbit, 1:nSpin)
       DO j = 1, nk
          Tau = Kink(Idx_Orbit, Idx_Spin)%Tau_End(i) -  Kink(Idx_Orbit, Idx_Spin)%Tau_Start(j)
          XG(i, j) = HybridizationFunc(Tau, Idx_Orbit, Idx_Spin)
       END DO
    END DO
    WRITE(*, *)

    DO iOrbit = 1, nOrbit
       DO iSpin = 1, nSpin
          WRITE(*, '("iOrbit=", i3, " iSpin=", i3, L)') iOrbit, iSpin, Kink(iOrbit, iSpin)%EmptySegment
          DO i = 1, Kink(iOrbit, iSpin)%nk
             WRITE(*, '(2f20.12, 100i3)') Kink(iOrbit, iSpin)%Tau_Start(i), Kink(iOrbit, iSpin)%Tau_End(i), Kink(iOrbit, iSpin)%n_Start(i, 1:nOrbit, 1:nSpin)
          END DO
       END DO
    END DO
    WRITE(*, *)
    
    IF (nk > 0) THEN
       CALL INVERSE(XG, XM, nk)

       WRITE(*, *)
       WRITE(*, *) '    ====== Hybridization matrix ====='
       DO i = 1, nk
          WRITE(*, '(30f20.8)') (XG(i, j), j = 1, nk)
       END DO
       
       WRITE(*, *)
       WRITE(*, *) '    ====== From fast update ====='
       DO i = 1, nk
          WRITE(*, '(30f20.8)') (Kink(Idx_Orbit, Idx_Spin)%Det(i, j), j = 1, nk)
       END DO
       
       WRITE(*, *)
       WRITE(*, *) '    ====== From matrix inversion ===== '
       DO i = 1, nk
          WRITE(*, '(30f20.8)') (XM(i, j), j = 1, nk)
       END DO
    END IF

    
    DEALLOCATE(XG)
    DEALLOCATE(XM)
    
  END SUBROUTINE Segment_Debug
  
END MODULE Segment_MonteCarlo
