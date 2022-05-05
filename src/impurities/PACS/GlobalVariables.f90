!******************************COPYRIGHT************************************************
! (c) Crown copyright, Gang Li < 2020 >. All rights reserved.
!
! This routine has been licensed to the valid partners for use and distribution
! under the collaboration agreement, subject to the terms and conditions set out therein.
!
! [ligang at Shanghaitech.edu.cn]
!******************************COPYRIGHT************************************************

MODULE GlobalVariables
  
  IMPLICIT NONE

!---------------------------------------------------------------------------------------
! Description:
!
!   Module GlobalVariables defines the top-level variables used throughout the program
!   and all derived functions, subroutines. Such as the presision, constants, etc. Every
!   variable is public, and all other routines should respect the definition here.
!
! Code hisotry: 2020.11.25
!
! Current Code Owner: Gang Li
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!---------------------------------------------------------------------------------------

  INTEGER,     PARAMETER :: SP = SELECTED_REAL_KIND(4)
  INTEGER,     PARAMETER :: DP = SELECTED_REAL_KIND(8)

  INTEGER,     PARAMETER :: STDOUT = 6

  REAL(DP),    PARAMETER :: ZERO = 0.d0, HALF = 0.5d0, ONE = 1.d0, TWO = 2.d0, THREE = 3.d0
  REAL(DP),    PARAMETER :: PI = ACOS(-ONE)
  REAL(DP),    PARAMETER :: epsilon = 1.d-8
  
  COMPLEX(DP), PARAMETER :: XI = DCMPLX(ZERO, ONE)


END MODULE GlobalVariables
