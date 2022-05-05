module MPI_mod
  
  USE MPI
  USE CTQMC_Math
  
  implicit none
  !--- mpi ---
  integer, parameter :: master = 0
  integer :: id, ntasks, rc
  integer, allocatable :: status(:)

contains
  !-------------------------------------------------------------------------------------------
  subroutine parallel_init
    IMPLICIT NONE
    LOGICAL :: MPI_On

    call MPI_Initialized(MPI_On, rc)
    if (rc /= MPI_SUCCESS) then
       print*, "MPI initialization failed"
       stop
    end if

    IF (.NOT. MPI_On) THEN   
      call MPI_INIT(rc)
      if (rc /= MPI_SUCCESS) then
         print*, "MPI initialization failed"
         stop
      end if
    END IF

  end subroutine parallel_init

  subroutine parallel_start

    IMPLICIT NONE
    INTEGER :: clock, i
    ! LOGICAL :: MPI_On

    ! call MPI_Initialized(MPI_On, rc)
    ! if (rc /= MPI_SUCCESS) then
    !    print*, "MPI initialization failed"
    !    stop
    ! end if
    
    ! IF (.NOT. MPI_On) THEN   

    !    call MPI_INIT(rc)
    !    if (rc /= MPI_SUCCESS) then
    !       print*, "MPI initialization failed"
    !       stop
    !    end if
       call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
       call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc)
       allocate(status(MPI_STATUS_SIZE))
       
       ! Initialize the Random Number Generator
       IF (id == master) then
          call SYSTEM_CLOCK(COUNT = clock)
          CALL RanGen_Ini(clock)
       END IF
       
       ! --- distribute the random seed ---
       if (id == master) then
          do i = 1, ntasks-1
             clock = clock + Int(ranw()*1000000)
             call MPI_SEND(clock, 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, rc)
          end do
       else
          call MPI_RECV(clock, 1, MPI_INTEGER, master, id, MPI_COMM_WORLD, status, rc)
          call RanGen_Ini(clock)
       end if

    ! END IF
  end subroutine parallel_start
  
  !-------------------------------------------------------------------------------------------
  subroutine parallel_end
    call MPI_FINALIZE(rc)
  end subroutine parallel_end

  !-------------------------------------------------------------------------------------------
  SUBROUTINE Start_Time
    IMPLICIT NONE

    INTEGER             :: Time_Start(8)
    
    CALL DATE_AND_TIME(VALUES=Time_Start)
    IF (id == master) THEN
       WRITE(*, '(A)') " --- >>> CTQMC calculation <<< --- "
       WRITE(*, *)
       WRITE(*, '(5x, "Starts at ", i4, 5(a, i2.2), " CST")') Time_Start(1), '/', Time_Start(2), '/', Time_Start(3), ' ', &
            Time_Start(5), ':', Time_Start(6), ':', Time_Start(7)
       WRITE(*, *)
    END IF
  END SUBROUTINE Start_Time

  SUBROUTINE End_Time
    IMPLICIT NONE
    INTEGER      :: Time_End(8)
    
    CALL DATE_AND_TIME(VALUES=Time_End)
    IF (id == master) THEN
       WRITE(*, *)
       WRITE(*, '(5x, "Ends at ", i4, 5(a, i2.2), " CST")') Time_End(1), '/', Time_End(2), '/', Time_End(3), ' ', &
            Time_End(5), ':', Time_End(6), ':', Time_End(7)
       WRITE(*, *)
    END IF
    
  END SUBROUTINE End_Time
  
end module MPI_mod
