!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE CTSEG_Run(nOrbit_, nOmega_, nTau_, nWarmUp_, nMeasure_, nBin_, Beta_, xU_, xJ_, Eloc_Imp, H0w_Imp, Sigma_Imp, ImpurityOccupancy_)

  ! Arguments
  ! =========

  !  iDMFT   : the current iteration number of DMFT loop, it is used to label the output files for the current loop.
  !  nOrbit_ : number of Orbital flavors. The current implementation assumes a orbital-diagonal form of hybridization function.
  !  nSpin_  : number of spin flavors. It is either 1 or 2.
  !  nOmega_ : number of Matsubara frequencies. We only determine and store Matsubara quantities on the positive axis.
  !  Beta_   : inverse temperature
  !  xU_     : Coulomb interaction parameter
  !  xJ_     : Hund's coupling constant. We always assume the SU(2) symmetric form of the interactions, i.e., U' = U - 2J  
  !  Eloc_Imp: local enery level of the impurity problem. 
  !  H0w_Imp : H0w_Imp(nOmega_, 2*nOrbit_), the hybridization function "H". The first (nOmega_, nOrbit_) enties store spin-up component, 
  !            the second half store spin-down component
  !            The Weiss Green's function G, hybridization function H, and the impurity level are related by the following relation
  !            G^-1 = iw_n + Eloc_Imp - H
  !  Sigma_Imp  : Impurity self-energy function calculated by the impurity solver.
  !  ImpurityOccupancy_(2*nOrbit_) : DMFT impurity orbital occupancy 

  USE MPI_MOD
  USE Segment_Util
  USE Segment_Phys, ONLY : Segment_Update

  IMPLICIT NONE
  
  !
  ! Purposes
  ! ========
  !   This is the top most level of routines working as the impurity
  !   solver.
  !   It is single purposed routine and only solve the impurity problem
  !   defined by chemical potential Xmu and hybridization function H0w.
  !   DMFT
  !   iterationn and bisection of chemical potential etc. are not
  !   contained
  !   here, which need to be supplmented by users.
  !
  ! Code developer: Gang Li
  ! Code history:  2021.7.25, brand new library mode
  !
  INTEGER,     INTENT(IN)  :: nOrbit_
  INTEGER,     INTENT(IN)  :: nOmega_
  INTEGER,     INTENT(IN)  :: nTau_
  INTEGER,     INTENT(IN)  :: nWarmUp_
  INTEGER,     INTENT(IN)  :: nMeasure_
  INTEGER,     INTENT(IN)  :: nBin_
  REAL(DP),    INTENT(IN)  :: Beta_
  REAL(DP),    INTENT(IN)  :: xU_
  REAL(DP),    INTENT(IN)  :: xJ_
  REAL(DP),    INTENT(IN)  :: Eloc_Imp(2*nOrbit_)
  COMPLEX(DP), INTENT(IN)  :: H0w_Imp(nOmega_, 2*nOrbit_)
  COMPLEX(DP), INTENT(OUT) :: Sigma_Imp(nOmega_, 2*nOrbit_)

  REAl(DP),    INTENT(OUT) :: ImpurityOccupancy_(2*nOrbit_)

  ! ... local vars ...
  INTEGER     :: iOrbit, iOmega
 
  ! pass the arguments to global variables used in CT_SEG
  Beta     = Beta_
  nOrbit   = nOrbit_
  nOmega   = nOmega_
  xU       = xU_
  xJ       = xJ_
  nTau     = nTau_
  nWarm    = nWarmUp_
  nMeasure = nMeasure_
  nBin     = nBin_
  
  ! --- allocate memory  ---
  CALL Segment_Initialize

  ! --- determine the hybridization function from the
  DO iOrbit = 1, nOrbit
     DO iOmega = 1, nOmega
        G0w(iOmega, iOrbit, 1) = One/(Omega(iOmega) + Eloc_Imp(iOrbit) - H0w_Imp(iOmega, iOrbit)) 
        G0w(iOmega, iOrbit, 2) = One/(Omega(iOmega) + Eloc_Imp(iOrbit+nOrbit) - H0w_Imp(iOmega, iOrbit+nOrbit))
     END DO
     Eimp(iOrbit, 1) = Eloc_Imp(iOrbit)  
     Eimp(iOrbit, 2) = Eloc_Imp(iOrbit+nOrbit)
     DO iOmega = 1, nOmega
        H0w(iOmega, iOrbit, 1) = H0w_Imp(iOmega, iOrbit)  
        H0w(iOmega, iOrbit, 2) = H0w_Imp(iOmega, iOrbit+nOrbit)
     END DO
  END DO
 
  IF (id == master) CALL Start_Time
  
  CALL Segment_HybridizationFunc
  CALL Segment_Update

  Sigma_Imp(1:nOmega, 1:nOrbit)          = Sigma(1:nOmega, 1:nOrbit, 1)
  Sigma_Imp(1:nOmega, nOrbit+1:2*nOrbit) = Sigma(1:nOmega, 1:nOrbit, 2)
  ImpurityOccupancy_(1:nOrbit)           = ImpurityOccupancy(1:nOrbit, 1)
  ImpurityOccupancy_(nOrbit+1:2*nOrbit)  = ImpurityOccupancy(1:nOrbit, 2)

  ! --- calculate the runtime cost for this iteration ---
  IF (id == master) CALL End_Time
  
  ! --- clear and deallocate memory ---
  CALL Segment_Finalize
  
END SUBROUTINE CTSEG_Run

