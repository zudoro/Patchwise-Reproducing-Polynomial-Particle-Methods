PROGRAM MAIN

     USE NEWTYPE
     USE GLBVAR
     USE DOMAIN
     USE PLOT
	 USE PACKAGE
     USE ASSEMBLE
     USE BOUNDARY
     USE LUDECOMPOSITION
     USE ERRORESTIMATE

     IMPLICIT NONE

	TYPE(SUBMAT) :: K
    REAL*8 :: FAT_K(3*DOF,3*DOF), F(3*DOF), DD, FCOPY(3*DOF), FAT_KCOPY(3*DOF, 3*DOF)
    INTEGER :: INDX(3*DOF)
    REAL :: CPTS, CPTE
	
	 CALL CPU_TIME(CPTS)
!    Get domain information and plot the domain
     CALL GETDOMAIN(OMEGA)
!     CALL PLOTDOMAIN()

!    Get background mesh and plot the mesh
     CALL GETBGMESH(OMEGA,BGMESH)
!      CALL PLOTBGMESH()

!	 Check pu and special (interpolation) functions
! 	 CALL CHECKBASIS()     

!    Assemble stiffness matrix
     CALL STIF(K)
	 
!	 ALLOCATE(F(3*DOF))
	 
!    Assemble load vector
     CALL LOADVEC(F)

!	 ALLOCATE(FAT_K(3*DOF, 3*DOF))
	 
!    Impose the boundary (Essential Boundary everywhere)
	 CALL IMPOSEBD(FAT_K, K, F)
	 
!	 ALLOCATE(INDX(3*DOF))
	 FAT_KCOPY(:,:) = FAT_K
	 FCOPY(:) = F
!    Solve the linearsystem Ax=B using LU decomposition.
     CALL LUDCMP(FAT_K,3*DOF,3*DOF,INDX,DD)
     CALL LUBKSB(FAT_K,3*DOF,3*DOF,INDX,F)
	 CALL RESIDUALNORM(FAT_KCOPY, F, FCOPY)
	 
	 CALL CPU_TIME(CPTE)
	 PRINT*, ""
	 PRINT*, "ELAPSED CPU TIME : ", CPTE - CPTS
	 PRINT*, ""

!    Estimate the Relative Energy Norm Error
!     CALL ENRGY(ROW, COL, K, X, DOF, NZ_COPY)
!    Estimate the max. Norm Error
     CALL MAXNORM(F)

END PROGRAM MAIN
