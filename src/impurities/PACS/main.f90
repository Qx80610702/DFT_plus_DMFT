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

  CALL CLEAN_MEMORY

end subroutine

SUBROUTINE CLEAN_MEMORY
  use MPI_mod
  use Segment_Phys
  use Segment_Util

  IMPLICIT NONE

  ! Moudle MPI_mod
  IF(ALLOCATED(status)) DEALLOCATE(status)

  ! Moudle Segment_Phys
  IF(ALLOCATED(nInsertSegment)) DEALLOCATE(nInsertSegment)
  IF(ALLOCATED(nInsertAntiSegment)) DEALLOCATE(nInsertAntiSegment)
  IF(ALLOCATED(nRemoveSegment)) DEALLOCATE(nRemoveSegment)
  IF(ALLOCATED(nRemoveAntiSegment)) DEALLOCATE(nRemoveAntiSegment)
  IF(ALLOCATED(nShiftSegment)) DEALLOCATE(nShiftSegment)
  IF(ALLOCATED(nEmptyFull)) DEALLOCATE(nEmptyFull)
  IF(ALLOCATED(Gt_Bin_Reduce)) DEALLOCATE(Gt_Bin_Reduce)
  IF(ALLOCATED(Gt_Bin)) DEALLOCATE(Gt_Bin)
  IF(ALLOCATED(MCsign_Bin)) DEALLOCATE(MCsign_Bin)
  IF(ALLOCATED(MCsign_Reduce)) DEALLOCATE(MCsign_Reduce)
  IF(ALLOCATED(Gst_Bin)) DEALLOCATE(Gst_Bin)
  IF(ALLOCATED(Gst_Bin_Reduce)) DEALLOCATE(Gst_Bin_Reduce)
  IF(ALLOCATED(nk_Probability)) DEALLOCATE(nk_Probability)

  ! Module Segment_Util
  IF(ALLOCATED(Omega)) DEALLOCATE(Omega)
  IF(ALLOCATED(G0w)) DEALLOCATE(G0w)
  IF(ALLOCATED(Gw)) DEALLOCATE(Gw)
  IF(ALLOCATED(Sigma)) DEALLOCATE(Sigma)
  IF(ALLOCATED(H0w)) DEALLOCATE(H0w)
  IF(ALLOCATED(Eimp)) DEALLOCATE(Eimp)
  IF(ALLOCATED(H0t)) DEALLOCATE(H0t)
  IF(ALLOCATED(ImpurityOccupancy)) DEALLOCATE(ImpurityOccupancy)
  IF(ALLOCATED(Gt)) DEALLOCATE(Gt)
  IF(ALLOCATED(Gst)) DEALLOCATE(Gst)
  IF(ALLOCATED(Umat)) DEALLOCATE(Umat)
  IF(ALLOCATED(LegenPoly)) DEALLOCATE(LegenPoly)
  IF(ALLOCATED(CoeffLegenPoly)) DEALLOCATE(CoeffLegenPoly)
  IF(ALLOCATED(LegenPolyMatsubara)) DEALLOCATE(LegenPolyMatsubara)

END SUBROUTINE CLEAN_MEMORY


