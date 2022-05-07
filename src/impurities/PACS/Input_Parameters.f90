! This module is used only to read in input parameters given by DFT+DMFT codes. 
! It's only called in program main. Do not confuse with MODULE Segment_Util

MODULE Inputs
  USE MPI_MOD
  USE GlobalVariables, only: DP

  IMPLICIT NONE

  INTEGER     ::       nOrbit = 3                   ! number of orbitals
  REAL(DP)    ::         Beta = 50.d0               ! inverse temperatures
  REAL(DP)    ::           xU = 0.d0                ! interaction parameters
  REAL(DP)    ::           xJ = 0.0d0               ! interaction parameters
  
  INTEGER     ::       nOmega = 200                 ! number of Matsubara frequencies
  INTEGER     ::         nTau = 2000                ! mesh of imaginary-time domain
  ! INTEGER     ::   nLegenPoly = 40                  ! optimal order of Legendre polynomial
  
  INTEGER     ::     nWarm    = 50000               ! number of MC warm-up steps
  INTEGER     ::     nMeasure = 1000000             ! number of MC measure steps
  INTEGER     ::     nBin     = 20                  ! binning analysis

  REAL(DP),    ALLOCATABLE :: Mu_vec(:)             ! Eloc_Imp(nspin*norbit), mu_vector
  COMPLEX(DP), ALLOCATABLE :: H0w(:, :)             ! H0w(nOmega, nspin*nOrbit_), hybridization function

  NAMELIST/Model_Parameter/nOrbit, Beta, xU, xJ, nOmega, nTau
  NAMELIST/MC_Parameter/nWarm, nMeasure, nBin

  CONTAINS

  SUBROUTINE Read_inputs
    call Read_Parameters
    call Read_H0w
    call Read_mu_vec
  END SUBROUTINE Read_inputs

  SUBROUTINE Read_Parameters

    IMPLICIT NONE
    
    !
    ! Purposes
    ! ========
    !    A user-defined parameter files can be provided to overwrite the defaul parameters
    !    defined and packed in this module.
    !
    ! Code developer: Gang Li
    ! Code history:   2021.1.5 adapted from code June, 2006
    !
    
    IF (ID == master) THEN
       OPEN(UNIT=1, FILE='parameters.in', STATUS='OLD')
       READ(1, nml=Model_Parameter)
       READ(1, nml=MC_Parameter)
       CLOSE(1)
    END IF

    CALL MPI_BCAST(nOrbit,       1, MPI_INTEGER, master, MPI_COMM_WORLD, rc)
    CALL MPI_BCAST(nOmega,       1, MPI_INTEGER, master, MPI_COMM_WORLD, rc)
    CALL MPI_BCAST(nTau,         1, MPI_INTEGER, master, MPI_COMM_WORLD, rc)
    ! CALL MPI_BCAST(nLegenPoly,   1, MPI_INTEGER, master, MPI_COMM_WORLD, rc)
    CALL MPI_BCAST(Beta,         1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, rc)
    CALL MPI_BCAST(xU,           1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, rc)
    CALL MPI_BCAST(xJ,           1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, rc)
    CALL MPI_BCAST(nWarm,        1, MPI_INTEGER, master, MPI_COMM_WORLD, rc)
    CALL MPI_BCAST(nMeasure,     1, MPI_INTEGER, master, MPI_COMM_WORLD, rc)
    CALL MPI_BCAST(nBin,         1, MPI_INTEGER, master, MPI_COMM_WORLD, rc)

  END SUBROUTINE READ_PARAMETERS

  SUBROUTINE Read_H0w
    IMPLICIT NONE

    INTEGER  :: iOmega, iOrbit
    REAL(DP) :: omegan
    REAL(DP) :: data(4*nOrbit)
    CHARACTER*10000 :: line

    IF(ALLOCATED(H0w)) DEALLOCATE(H0w)
    ALLOCATE(H0w(nOmega, 2*nOrbit))

    IF (ID == master) THEN
      ! Read in the imaginary frequency hybridization function      
      OPEN(UNIT=1, FILE='delta.dat', STATUS='OLD')
      DO iOmega = 1, nOmega
        READ(1,'(A)') line
        READ(LINE,*) omegan, (data(:))
        DO iOrbit = 1, 2*nOrbit
          H0w(iOmega, iOrbit) = DCMPLX(data(2*iOrbit-1), data(2*iOrbit))
        END DO
      END DO
      CLOSE(1)
    END IF

    DO iOrbit = 1, 2*nOrbit
      CALL MPI_BCAST(H0w(:,iorbit), nOmega, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    ENDDO

  END SUBROUTINE Read_H0w

  SUBROUTINE Read_mu_vec
    IMPLICIT NONE

    CHARACTER*1000 :: line

    IF(ALLOCATED(Mu_vec)) DEALLOCATE(Mu_vec) 
    ALLOCATE(Mu_vec(2*nOrbit))

    IF (ID == master) THEN
      OPEN(UNIT=1, FILE='mu_vector.dat', STATUS='OLD')
      READ(1, '(A)') line
      READ(line, *) (Mu_vec(:))
      CLOSE(1) 
    END IF

    CALL MPI_BCAST(Mu_vec, 2*nOrbit, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, rc)

  END SUBROUTINE Read_mu_vec

END MODULE Inputs