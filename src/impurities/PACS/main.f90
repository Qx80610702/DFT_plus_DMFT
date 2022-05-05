! program main
!   use MPI_mod
!   use Inputs

!   IMPLICIT NONE

!   integer :: iDMFT, nOrbit_, nOmega_, nTau_, nWarmUp_, nMeasure_, nBin_
!   real(8) :: Beta_, xU_, xJ_
!   real(8),    allocatable :: Eloc_Imp(:), ImpurityOccupancy_(:)
!   complex(8), allocatable :: H0w_Imp(:, :), Sigma_Imp(:, :)


!   CALL parallel_init

!   ! CALL PACS_Run

!   CALL parallel_start

!   CALL Read_inputs

!   ! iDMFT     = 1
!   ! nOrbit_   = 3
!   ! nOmega_   = 300
!   ! nTau_     = 2000
!   ! nWarmUp_  = 100000
!   ! nMeasure_ = 1000000
!   ! nBin_     = 10
  
!   ALLOCATE(Eloc_Imp(2*nOrbit))
!   ALLOCATE(ImpurityOccupancy_(2*nOrbit))
!   ALLOCATE(H0w_Imp(nOmega, 2*nOrbit))
!   ALLOCATE(Sigma_Imp(nOmega, 2*nOrbit))

!   CALL CTSEG_Run(nOrbit, nOmega, nTau, nWarm, nMeasure, nBin, Beta, xU, xJ, mu_vec, H0w, Sigma_Imp, ImpurityOccupancy_) 

!   DEALLOCATE(Eloc_Imp, ImpurityOccupancy_, H0w_Imp, Sigma_Imp)

!   CALL parallel_end

! end program 

subroutine PACS_Run
  use MPI_mod
  use Inputs

  IMPLICIT NONE
  ! Arguments
  ! INTEGER,     INTENT(IN)  :: myrank
  ! INTEGER,     INTENT(IN)  :: nprocess

  ! Variables
  real(8),    allocatable :: Eloc_Imp(:), ImpurityOccupancy_(:)
  complex(8), allocatable :: H0w_Imp(:, :), Sigma_Imp(:, :)

  CALL parallel_start

  CALL Read_inputs

  ! iDMFT     = 1
  ! nOrbit_   = 3
  ! nOmega_   = 300
  ! nTau_     = 2000
  ! nWarmUp_  = 100000
  ! nMeasure_ = 1000000
  ! nBin_     = 10
  
  ALLOCATE(Eloc_Imp(2*nOrbit))
  ALLOCATE(ImpurityOccupancy_(2*nOrbit))
  ALLOCATE(H0w_Imp(nOmega, 2*nOrbit))
  ALLOCATE(Sigma_Imp(nOmega, 2*nOrbit))

  CALL CTSEG_Run(nOrbit, nOmega, nTau, nWarm, nMeasure, nBin, Beta, xU, xJ, mu_vec, H0w, Sigma_Imp, ImpurityOccupancy_) 

  DEALLOCATE(Eloc_Imp, ImpurityOccupancy_, H0w_Imp, Sigma_Imp)

end subroutine
