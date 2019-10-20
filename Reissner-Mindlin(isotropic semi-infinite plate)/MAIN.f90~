PROGRAM MAIN

     USE NEWTYPE
     USE GLBVAR
     USE DOMAIN
     USE PLOT
	 USE PACKAGE
     USE ASSEMBLE
     USE BOUNDARY
     USE BICGSTAB
     USE ERRORESTIMATE

     IMPLICIT NONE

	TYPE(SUBMAT) :: K
    REAL*8 :: F(3*DOF), X(3*DOF), DNR, ERROR
    REAL*8, ALLOCATABLE :: FULL_K(:)
    INTEGER, ALLOCATABLE :: FULL_R(:), FULL_C(:)
	INTEGER :: NZ

!    Get domain information and plot the domain
     CALL GETDOMAIN(OMEGA)
!     CALL PLOTDOMAIN()

!    Get background mesh and plot the mesh
     CALL GETBGMESH(OMEGA,BGMESH)
     CALL PLOTBGMESH()

!	 Check pu and special (interpolation) functions
	 CALL CHECKBASIS()     
!    Assemble stiffness matrix
     CALL STIF(K)
!     ALLOCATE (FULLA(2*NZ), FULLROW(2*NZ), FULLCOL(2*NZ))
!    Assemble load vector
     CALL LOADVEC(F)
     
     ALLOCATE(FULL_K(9*DOF**2), FULL_R(9*DOF**2), FULL_C(9*DOF**2))
     
!    Impose the boundary (Essential Boundary everywhere)
	 CALL IMPOSEBD(FULL_R, FULL_C, FULL_K, NZ, K, F)
!    Solve the linearsystem Ax=B using BiConjugate Gradient Stabilized method
     CALL SOLVER(X, DNR, FULL_R, FULL_C, FULL_K, F, 3*DOF, NZ)
     PRINT*, ""
     PRINT*, "RESIDUAL NORM ERROR : ", DNR

	 DEALLOCATE(FULL_K, FULL_R, FULL_C)
	 
!    Estimate the Relative Energy Norm Error
!     CALL ENRGY(ROW, COL, K, X, DOF, NZ_COPY)
!    Estimate the max. Norm Error
     CALL MAXNORM(X)

END PROGRAM MAIN
